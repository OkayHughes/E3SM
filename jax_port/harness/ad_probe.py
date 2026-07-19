#!/usr/bin/env python3
"""Differentiability audit of the scream_jax process steps (host venv).

For each process: attempt reverse-mode grad and forward-mode jvp of a
scalar reduction of the outputs w.r.t. a key input on golden states;
report success/failure, NaN/Inf counts in the gradient, and (for
processes where both modes work) a finite-difference directional-
derivative check at several step sizes.

Usage: cd jax_port && .venv/bin/python harness/ad_probe.py [--ncol 8]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / "physics_suite_218x72_dt1800_2steps.npz"
DATA = REPO.parent / "e3sm-inputdata" / "atm" / "scream"


def load_state(ncol):
    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    s = {name: np.asarray(z[f"{name}__step0"], dtype=np.float64)[:ncol]
         for name in meta["fields"]}
    import netCDF4
    ds = netCDF4.Dataset(DATA / "init" /
                         "screami_unit_tests_ne2np4L72_20220822.nc")
    geo = {k: np.array(ds[k][:]) for k in ("hyam", "hybm", "lat", "lon",
                                           "area")}
    ds.close()
    return s, geo, meta


def probe(name, f, x0, out_names=None):
    """Try value, reverse grad, forward jvp of scalar(f). Returns dict."""
    res = {"name": name}
    v = np.random.default_rng(0).normal(size=x0.shape)
    v /= np.linalg.norm(v)

    def scalar(x):
        out = f(x)
        return jnp.sum(sum(jnp.sum(o) for o in out)) if isinstance(
            out, (tuple, list)) else jnp.sum(out)

    try:
        t0 = time.time()
        val = float(scalar(jnp.asarray(x0)))
        res["primal"] = f"{val:.6e} ({time.time()-t0:.1f}s)"
    except Exception as e:  # noqa: BLE001
        res["primal"] = f"FAIL: {type(e).__name__}: {str(e)[:120]}"
        return res

    try:
        t0 = time.time()
        _, tang = jax.jvp(scalar, (jnp.asarray(x0),), (jnp.asarray(v),))
        tang = float(tang)
        res["jvp"] = f"{tang:.6e} ({time.time()-t0:.1f}s)"
        res["_tang"] = tang
    except Exception as e:  # noqa: BLE001
        res["jvp"] = f"FAIL: {type(e).__name__}: {str(e)[:120]}"

    try:
        t0 = time.time()
        g = jax.grad(scalar)(jnp.asarray(x0))
        g = np.asarray(g)
        nan, inf = int(np.isnan(g).sum()), int(np.isinf(g).sum())
        res["grad"] = (f"|g|={np.linalg.norm(g[np.isfinite(g)]):.3e} "
                       f"NaN={nan}/{g.size} Inf={inf} "
                       f"({time.time()-t0:.1f}s)")
        res["_grad_dot_v"] = float(np.where(np.isfinite(g), g, 0.0)
                                   .ravel() @ v.ravel())
        res["_nan"] = nan
    except Exception as e:  # noqa: BLE001
        res["grad"] = f"FAIL: {type(e).__name__}: {str(e)[:160]}"

    # FD directional derivative at several relative steps
    if "_tang" in res:
        scale = max(np.abs(x0).mean(), 1e-30)
        fd = {}
        for eps_rel in (1e-3, 1e-5, 1e-7):
            eps = eps_rel * scale
            fp = float(scalar(jnp.asarray(x0 + eps * v)))
            fm = float(scalar(jnp.asarray(x0 - eps * v)))
            fd[eps_rel] = (fp - fm) / (2 * eps)
        res["fd"] = {k: f"{val:.6e}" for k, val in fd.items()}
        best = min(abs(fd[k] - res["_tang"]) /
                   max(abs(res["_tang"]), 1e-30) for k in fd)
        res["fd_vs_jvp_best_rel"] = f"{best:.2e}"
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncol", type=int, default=8)
    args = ap.parse_args()
    s, geo, meta = load_state(args.ncol)
    n = args.ncol

    from scream_jax.driver import (CLDFRAC_ICE_4OUT_THRESHOLD,
                                   CLDFRAC_ICE_THRESHOLD, P0, SHOC_PARAMS,
                                   ScreamPhysics)
    from scream_jax.cld_fraction.main import cld_fraction_main
    from scream_jax.shoc.process import shoc_process_step
    from scream_jax.p3.process import p3_process_step
    from scream_jax.rrtmgp.process import rrtmgp_process_step
    from scream_jax.foundation.thermo import calculate_dx_from_area
    from scream_jax import shoc

    phys = ScreamPhysics(DATA, geo["hyam"], geo["hybm"], geo["lat"][:n],
                         geo["lon"][:n],
                         np.asarray(calculate_dx_from_area(
                             geo["area"][:n], geo["lat"][:n])))
    dt = 300.0
    results = []

    # ---- cld_fraction: d(sum cldfrac_tot)/d(qi) ----
    results.append(probe(
        "cld_fraction  d(sum cldfrac)/d(qi)",
        lambda qi: cld_fraction_main(CLDFRAC_ICE_THRESHOLD,
                                     CLDFRAC_ICE_4OUT_THRESHOLD,
                                     qi, jnp.asarray(s["cldfrac_liq"]))[1],
        s["qi"]))

    # ---- shoc: d(sum T_out)/d(T_mid) ----
    def shoc_f(T):
        out = shoc_process_step(
            dt, phys.npbl, jnp.asarray(phys.cell_length),
            SHOC_PARAMS["lambda_low"], SHOC_PARAMS["lambda_high"],
            SHOC_PARAMS["lambda_slope"], SHOC_PARAMS["lambda_thresh"],
            SHOC_PARAMS["thl2tune"], SHOC_PARAMS["qw2tune"],
            SHOC_PARAMS["qwthl2tune"], SHOC_PARAMS["w2tune"],
            SHOC_PARAMS["length_fac"], SHOC_PARAMS["c_diag_3rd_mom"],
            SHOC_PARAMS["ckh"], SHOC_PARAMS["ckm"], False, False,
            T, jnp.asarray(s["p_mid"]), jnp.asarray(s["p_int"]),
            jnp.asarray(s["pseudo_density"]), jnp.asarray(s["omega"]),
            jnp.asarray(s["phis"]), jnp.asarray(s["surf_sens_flux"]),
            jnp.asarray(s["surf_evap"]),
            jnp.asarray(s["surf_mom_flux"][:, 0]),
            jnp.asarray(s["surf_mom_flux"][:, 1]),
            jnp.asarray(s["qv"]), jnp.asarray(s["qc"]),
            jnp.asarray(s["tke"]), jnp.asarray(s["horiz_winds"][:, 0, :]),
            jnp.asarray(s["horiz_winds"][:, 1, :]),
            jnp.asarray(s["cldfrac_liq"]),
            jnp.asarray(s["sgs_buoy_flux"]),
            jnp.asarray(s["eddy_diff_mom"]))
        return out["T_mid"]
    results.append(probe("shoc          d(sum T_out)/d(T_mid)", shoc_f,
                         s["T_mid"]))

    # ---- p3: d(sum qc_out + qr_out)/d(qc) ----
    def p3_f(qc):
        out = p3_process_step(
            dt, True, True, True, False, False, False, False, False,
            jnp.asarray(s["T_mid"]), jnp.asarray(s["p_mid"]),
            jnp.asarray(s["p_dry_mid"]), jnp.asarray(s["pseudo_density"]),
            jnp.asarray(s["pseudo_density_dry"]),
            jnp.asarray(s["cldfrac_tot"]),
            jnp.asarray(s["qv"]), qc, jnp.asarray(s["nc"]),
            jnp.asarray(s["qr"]), jnp.asarray(s["nr"]),
            jnp.asarray(s["qi"]), jnp.asarray(s["qm"]),
            jnp.asarray(s["ni"]), jnp.asarray(s["bm"]),
            jnp.asarray(s["qv_prev_micro_step"]),
            jnp.asarray(s["T_prev_micro_step"]),
            jnp.asarray(s["nc_nuceat_tend"]), jnp.asarray(s["nccn"]),
            jnp.asarray(s["ni_activated"]),
            jnp.asarray(s["inv_qc_relvar"]),
            jnp.asarray(s["precip_liq_surf_mass"]),
            jnp.asarray(s["precip_ice_surf_mass"]),
            phys.p3_tables, phys.p3_opts,
            sed_use_while_loop=False)
        return out["qc"] + out["qr"]
    results.append(probe("p3            d(sum qc+qr out)/d(qc)", p3_f,
                         s["qc"]))

    # ---- rrtmgp: d(sum LW_flux_up)/d(T_mid) ----
    def rad_f(T):
        out = rrtmgp_process_step(
            phys.kd_sw, phys.kd_lw, phys.co_sw, phys.co_lw,
            phys.rrtmgp_params, 1800.0, 0, 2021, 285.5,
            jnp.asarray(phys.lat), jnp.asarray(phys.lon),
            T, jnp.asarray(s["p_mid"]), jnp.asarray(s["p_int"]),
            jnp.asarray(s["pseudo_density"]),
            jnp.asarray(s["sfc_alb_dir_vis"]),
            jnp.asarray(s["sfc_alb_dir_nir"]),
            jnp.asarray(s["sfc_alb_dif_vis"]),
            jnp.asarray(s["sfc_alb_dif_nir"]),
            jnp.asarray(s["qv"]), jnp.asarray(s["qc"]),
            jnp.asarray(s["nc"]), jnp.asarray(s["qi"]),
            jnp.asarray(s["cldfrac_tot"]),
            jnp.asarray(s["eff_radius_qc"]),
            jnp.asarray(s["eff_radius_qi"]),
            jnp.asarray(s["surf_lw_flux_up"]),
            jnp.asarray(s["o3_volume_mix_ratio"]),
            jnp.asarray(s["rad_heating_pdel"]),
            aero_tau_sw=jnp.asarray(s["aero_tau_sw"]),
            aero_ssa_sw=jnp.asarray(s["aero_ssa_sw"]),
            aero_g_sw=jnp.asarray(s["aero_g_sw"]),
            aero_tau_lw=jnp.asarray(s["aero_tau_lw"]))
        return out["LW_flux_up"]
    results.append(probe("rrtmgp        d(sum LWup)/d(T_mid)", rad_f,
                         s["T_mid"]))

    print("\n" + "=" * 78)
    for r in results:
        print(f"\n{r['name']}")
        for k in ("primal", "jvp", "grad", "fd", "fd_vs_jvp_best_rel"):
            if k in r:
                print(f"  {k:20s} {r[k]}")
    print("=" * 78)


if __name__ == "__main__":
    main()
