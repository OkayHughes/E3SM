#!/usr/bin/env python3
"""Generate CLUBB slice-1 goldens (run in the scream-dev container).

Drives the REAL CLUBB Fortran (eam_clubb_f, the full unmodified stack,
see build_clubb.py) exactly as EAM does and records golden archives:

  golden/clubb_grid.npz         setup_grid arrays + zt2zm/zm2zt/ddzt/
                                ddzm on regime-sweeping test fields,
                                on 3 grids (EAM-like 73-level, coarse
                                41-level, irregular).
  golden/clubb_tridag.npz       tridag_solve (LAPACK dgtsv) systems:
                                diagonally-dominant CLUBB-like,
                                pivot-exercising non-dominant, and
                                multi-RHS cases.
  golden/clubb_sat.npz          sat_mixrat_liq/ice over a (p,T) grid
                                (150-330 K incl. both polynomial-clip
                                branches and the esat~p guard).
  golden/clubb_pdf_closure.npz  pdf_closure (ADG1) over 40 regime-
                                sweeping columns x 73 levels: stable /
                                convective-cloudy / extreme-skewness /
                                degenerate-wp2 / degenerate-xp2 /
                                saturated / cold-cirrus / randomized;
                                plus the Skx and sigma_sqd_w helper
                                inputs+outputs used to build them.

EAMv3 tunable parameters: CLUBB compiled-in defaults from the real
read_parameters(-99) with the namelist_defaults_eam.xml phys="default"
clubb_param_nl overrides applied, INCLUDING the coupling rules
read_parameters applies when a value is namelist-set (verbatim from
parameters_tunable.F90): clubb_C2rt set -> C2thl = C2rt and
C2rtthl = C2rt * 1.3 (clubb_c2thl/clubb_c2rtthl are NOT in the EAMv3
namelist); clubb_C6rt NOT set -> C6rt/C6thl keep their defaults.
The full packed params vector is recorded in the golden metadata.
"""
import json
import subprocess
from pathlib import Path

import numpy as np

from fbuild import REPO, import_extension

GOLDEN = Path(__file__).resolve().parent.parent / "golden"
GOLDEN.mkdir(exist_ok=True)

f = import_extension("eam_clubb_f")
d = f.clubb_driver

IDX_NAMES = ["C1", "C1b", "C1c", "C2rt", "C2thl", "C2rtthl", "C6rt",
             "C6rtb", "C6rtc", "C6thl", "C6thlb", "C6thlc", "C7", "C7b",
             "C8", "C11", "C11b", "C11c", "C14", "beta", "gamma_coef",
             "gamma_coefb", "gamma_coefc", "mu", "nu1", "c_K10",
             "c_K10h", "wpxp_L_thresh", "altitude_threshold",
             "Skw_denom_coef", "Skw_max_mag"]

# EAMv3 phys="default" clubb_param_nl values
# (bld/namelist_files/namelist_defaults_eam.xml; clubb_wpxp_L_thresh
# is the general default at line 871).  clubb_C2rt=1.75 triggers the
# read_parameters coupling C2thl=C2rt, C2rtthl=1.3*C2rt.
EAMV3_OVERRIDES = {
    "C1": 2.4, "C1b": 2.8, "C1c": 0.75,
    "C2rt": 1.75, "C2thl": 1.75, "C2rtthl": 1.75 * 1.3,
    "C6rtb": 7.50, "C6rtc": 0.50, "C6thlb": 7.50, "C6thlc": 0.50,
    "C8": 5.2, "C11": 0.70, "C11b": 0.20, "C11c": 0.85,
    "gamma_coef": 0.12, "gamma_coefb": 0.28, "gamma_coefc": 1.2,
    "mu": 0.0005, "wpxp_L_thresh": 60.0,
}


def eamv3_params():
    nparams = d.drv_get_nparams()
    params = np.array(d.drv_default_params(nparams))
    idx = {n: i - 1 for n, i in zip(IDX_NAMES, d.drv_param_indices())}
    for name, val in EAMV3_OVERRIDES.items():
        params[idx[name]] = val
    return params, idx


def eam_like_grid(nz=73, ztop=41000.0, seed=None):
    """CLUBB-oriented zi/zt heights (index 1 = surface): EAM-like
    hyperbolic-stretched interfaces (dz ~ 25 m at the surface to
    ~2.5 km at the model top), plus the below-surface ghost point."""
    s = np.linspace(0.0, 1.0, nz)
    zi = ztop * (np.expm1(3.2 * s) / np.expm1(3.2))
    zi[0] = 0.0
    if seed is not None:
        rng = np.random.default_rng(seed)
        mid = zi[1:-1]
        gaps = np.diff(zi)
        jitter = rng.uniform(-0.25, 0.25, mid.size) \
            * np.minimum(gaps[:-1], gaps[1:])
        zi[1:-1] = mid + jitter
    zt = np.empty(nz)
    zt[1:] = 0.5 * (zi[1:] + zi[:-1])
    zt[0] = -zt[1]
    return zi, zt


def gen_grid():
    # nz is fixed at 73 across grids: setup_grid (which sets gr%nz and
    # allocates) runs once inside setup_clubb_core_api, exactly as in
    # EAM where pverp never changes; subsequent drv_setup calls take
    # the setup_grid_heights_api refresh path.
    grids = {"a": eam_like_grid(73), "b": eam_like_grid(73, 28000.0),
             "c": eam_like_grid(73, 41000.0, seed=7)}
    params, _ = eamv3_params()
    out = {}
    rng = np.random.default_rng(42)
    for tag, (zi, zt) in grids.items():
        nz = zi.size
        err = d.drv_setup(29, params, zi, zt)
        assert err == 0
        zm_o, zt_o, dzm, dzt, idzm, idzt, wt2m, wm2t = d.drv_grid_arrays(nz)
        out[f"{tag}_zi_in"] = zi
        out[f"{tag}_zt_in"] = zt
        out[f"{tag}_zm"] = zm_o
        out[f"{tag}_zt"] = zt_o
        out[f"{tag}_dzm"] = dzm
        out[f"{tag}_dzt"] = dzt
        out[f"{tag}_invrs_dzm"] = idzm
        out[f"{tag}_invrs_dzt"] = idzt
        out[f"{tag}_w_zt2zm"] = wt2m
        out[f"{tag}_w_zm2zt"] = wm2t
        # operator sweeps: smooth, oscillatory, random, constant fields
        fields = np.stack([
            300.0 + 0.01 * zt, np.sin(zt / 2000.0),
            rng.normal(size=nz), np.full(nz, 5.0),
            np.exp(-zt / 8000.0)])
        fields_m = np.stack([
            1.0 + zm_o / 1e4, np.cos(zm_o / 3000.0),
            rng.normal(size=nz), np.full(nz, -2.0),
            np.exp(-zm_o / 7000.0)])
        z2m = np.empty_like(fields)
        m2z = np.empty_like(fields)
        dzt_o = np.empty_like(fields)
        dzm_o = np.empty_like(fields)
        for i in range(fields.shape[0]):
            z2m[i], m2z[i], dzt_o[i], dzm_o[i] = \
                d.drv_grid_ops(fields[i], fields_m[i])
        out[f"{tag}_fields_zt"] = fields
        out[f"{tag}_fields_zm"] = fields_m
        out[f"{tag}_zt2zm"] = z2m
        out[f"{tag}_zm2zt"] = m2z
        out[f"{tag}_ddzt"] = dzt_o
        out[f"{tag}_ddzm"] = dzm_o
    return out


def gen_tridag():
    rng = np.random.default_rng(1)
    cases = []
    for n, nrhs in [(73, 1), (73, 2), (73, 5), (41, 1), (5, 3), (2, 1)]:
        # CLUBB-like diagonally dominant (implicit diffusion matrix)
        sub = -rng.uniform(0.1, 1.0, n)
        sup = -rng.uniform(0.1, 1.0, n)
        dia = 1.0 + rng.uniform(0.5, 2.0, n) + np.abs(sub) + np.abs(sup)
        rhs = rng.normal(size=(n, nrhs))
        cases.append((sup, dia, sub, rhs))
        # non-dominant, pivoting exercised
        sub = rng.normal(size=n) * 10.0
        sup = rng.normal(size=n)
        dia = rng.normal(size=n) * 0.1
        rhs = rng.normal(size=(n, nrhs))
        cases.append((sup, dia, sub, rhs))
        # mixed magnitudes
        sub = rng.normal(size=n) * 10.0 ** rng.integers(-6, 6, n)
        sup = rng.normal(size=n) * 10.0 ** rng.integers(-6, 6, n)
        dia = rng.normal(size=n) * 10.0 ** rng.integers(-6, 6, n)
        rhs = rng.normal(size=(n, nrhs))
        cases.append((sup, dia, sub, rhs))
    out = {"ncases": np.array(len(cases))}
    for i, (sup, dia, sub, rhs) in enumerate(cases):
        sol, err = d.drv_tridag_solve(sup, dia, sub, rhs)
        assert err == 0, (i, err)
        out[f"{i}_supd"] = sup
        out[f"{i}_diag"] = dia
        out[f"{i}_subd"] = sub
        out[f"{i}_rhs"] = rhs
        out[f"{i}_solution"] = sol
    return out


def gen_sat():
    t = np.concatenate([np.linspace(150.0, 330.0, 181),
                        [273.15, 273.15 - 85.0, 273.15 - 90.0,
                         273.15 - 84.999, 273.15 - 89.999, 200.0]])
    p = np.concatenate([np.geomspace(50.0, 1.05e5, 181),
                        [1.0e5, 500.0, 120.0, 3000.0, 2.0, 700.0]])
    # all (p, T) pairs, flattened
    pp, tt = np.meshgrid(p, t, indexing="ij")
    pp, tt = pp.ravel(), tt.ravel()
    rsl, rsi = d.drv_sat(pp, tt)
    return {"p": pp, "t": tt, "rsl": rsl, "rsi": rsi}


def _build_columns(nz, zt, params, idx, rng):
    """40 regime-sweeping pdf_closure input columns (CLUBB
    orientation, level 1 = below-surface ghost)."""
    gamma_coef = params[idx["gamma_coef"]]
    gamma_coefb = params[idx["gamma_coefb"]]
    gamma_coefc = params[idx["gamma_coefc"]]
    skw_denom_coef = params[idx["Skw_denom_coef"]]

    ncol = 40
    z = zt.copy()
    cols = {k: np.zeros((ncol, nz)) for k in
            ["p", "exner", "thv_ds", "wm", "wp2", "wp3", "skw", "skthl",
             "skrt", "rtm", "rtp2", "wprtp", "thlm", "thlp2", "wpthlp",
             "um", "up2", "upwp", "vm", "vp2", "vpwp", "rtpthlp",
             "gamma_skw_fnc", "sigma_sqd_w", "wp2_in3", "rtp3", "thlp3"]}

    for c in range(ncol):
        fam = c % 8
        h = rng.uniform(7000.0, 8500.0)
        p = 1.0e5 * np.exp(-np.maximum(z, 0.0) / h)
        p[0] = p[1]
        exner = (p / 1.0e5) ** (287.042 / 1004.64)
        t_sfc = rng.uniform(275.0, 302.0)
        lapse = rng.uniform(0.005, 0.0075)
        T = np.maximum(t_sfc - lapse * np.maximum(z, 0.0), 195.0)
        thlm = T / exner
        rsl, _ = d.drv_sat(p, T)

        zpbl = rng.uniform(800.0, 2500.0)
        bl = np.exp(-np.maximum(z, 0.0) / zpbl)

        if fam == 0:    # stable BL, negative skewness
            rh = rng.uniform(0.3, 0.9)
            wp2 = 0.02 * bl + 1e-3
            skw_t = -rng.uniform(0.2, 1.5) * bl
        elif fam == 1:  # convective, cloudy layer
            rh = rng.uniform(0.85, 1.0)
            wp2 = rng.uniform(0.3, 1.0) * bl + 1e-3
            skw_t = rng.uniform(0.5, 2.5) * bl
        elif fam == 2:  # extreme skewness, both signs
            rh = rng.uniform(0.5, 1.0)
            wp2 = rng.uniform(0.2, 0.8) * bl + 1e-3
            skw_t = np.sign(rng.normal()) * rng.uniform(3.5, 4.4) * bl
        elif fam == 3:  # degenerate wp2 (below w_tol_sqd = 4e-4)
            rh = rng.uniform(0.4, 1.0)
            wp2 = np.full(nz, rng.uniform(1e-6, 3e-4))
            skw_t = np.zeros(nz)
        elif fam == 4:  # degenerate scalar variances in layers
            rh = rng.uniform(0.5, 0.98)
            wp2 = 0.3 * bl + 1e-3
            skw_t = rng.uniform(0.2, 1.0) * bl
        elif fam == 5:  # saturated column (cf ~ 1)
            rh = rng.uniform(1.01, 1.08)
            wp2 = 0.4 * bl + 1e-3
            skw_t = rng.uniform(0.3, 1.5) * bl
        elif fam == 6:  # cold/cirrus emphasis (ice supersat branch)
            rh = rng.uniform(0.7, 1.05)
            wp2 = 0.05 + 0.2 * np.exp(-((z - 9000.0) / 2500.0) ** 2)
            skw_t = 0.8 * np.exp(-((z - 9000.0) / 2500.0) ** 2)
        else:           # randomized
            rh = rng.uniform(0.2, 1.05)
            wp2 = rng.uniform(1e-4, 1.2, nz) * bl + 5e-4
            skw_t = rng.normal(0.0, 1.5, nz) * bl

        rtm = np.clip(rh * rsl, 1e-7, 0.025)
        # moments
        rtp2 = (rng.uniform(0.02, 0.25) * rtm) ** 2 * bl + 1e-14
        thlp2 = rng.uniform(0.05, 1.5) ** 2 * bl + 1e-6
        if fam == 4:
            rtp2[nz // 2:] = 1e-17          # below rt_tol^2
            thlp2[: nz // 4] = 1e-5         # below thl_tol^2
        corr_wrt = rng.uniform(-0.7, 0.9, nz)
        corr_wthl = rng.uniform(-0.9, 0.7, nz)
        corr_rtthl = rng.uniform(-0.95, 0.95, nz)
        wprtp = corr_wrt * np.sqrt(wp2 * rtp2)
        wpthlp = corr_wthl * np.sqrt(wp2 * thlp2)
        rtpthlp = corr_rtthl * np.sqrt(rtp2 * thlp2)
        wp3 = skw_t * wp2 ** 1.5
        rtp3 = rng.normal(0.0, 0.5, nz) * rtp2 ** 1.5
        thlp3 = rng.normal(0.0, 0.5, nz) * thlp2 ** 1.5

        um = 5.0 + 10.0 * np.tanh(z / 3000.0) + rng.normal(0, 2)
        vm = -3.0 + rng.normal(0, 2) * bl
        up2 = rng.uniform(0.05, 0.6) * bl + 1e-4
        vp2 = rng.uniform(0.05, 0.6) * bl + 1e-4
        upwp = rng.uniform(-0.6, 0.6, nz) * np.sqrt(wp2 * up2)
        vpwp = rng.uniform(-0.6, 0.6, nz) * np.sqrt(wp2 * vp2)
        wm = rng.normal(0.0, 0.02, nz) * bl

        thv_ds = thlm.copy()
        # Skewness/sigma exactly as pdf_closure_driver builds them
        skw = d.drv_skx(wp2, wp3, 2.0e-2)
        skthl = d.drv_skx(thlp2, thlp3, 1.0e-2)
        skrt = d.drv_skx(rtp2, rtp3, 1.0e-8)
        gam = gamma_coefb + (gamma_coef - gamma_coefb) * np.exp(
            -0.5 * (skw / gamma_coefc) ** 2)
        sig = d.drv_sigma_sqd_w(gam, wp2, thlp2, rtp2, up2, vp2,
                                wpthlp, wprtp, upwp, vpwp)

        for k, v in [("p", p), ("exner", exner), ("thv_ds", thv_ds),
                     ("wm", wm), ("wp2", wp2), ("wp3", wp3),
                     ("skw", skw), ("skthl", skthl), ("skrt", skrt),
                     ("rtm", rtm), ("rtp2", rtp2), ("wprtp", wprtp),
                     ("thlm", thlm), ("thlp2", thlp2),
                     ("wpthlp", wpthlp), ("um", um), ("up2", up2),
                     ("upwp", upwp), ("vm", vm), ("vp2", vp2),
                     ("vpwp", vpwp), ("rtpthlp", rtpthlp),
                     ("gamma_skw_fnc", gam), ("sigma_sqd_w", sig),
                     ("wp2_in3", wp2), ("rtp3", rtp3),
                     ("thlp3", thlp3)]:
            cols[k][c] = v
    del cols["wp2_in3"]
    assert skw_denom_coef == 0.0  # CLUBB_CAM compile-time default
    return cols


def gen_pdf_closure(params, idx):
    zi, zt = eam_like_grid(73)
    nz = zi.size
    err = d.drv_setup(29, params, zi, zt)
    assert err == 0
    rng = np.random.default_rng(2024)
    cols = _build_columns(nz, zt, params, idx, rng)
    ncol = cols["p"].shape[0]

    moments = np.zeros((ncol, nz, 17))
    pdfp = np.zeros((ncol, nz, 47))
    sig_out = np.zeros((ncol, nz))
    for c in range(ncol):
        m, pp, so, perr = d.drv_pdf_closure(
            cols["p"][c], cols["exner"][c], cols["thv_ds"][c],
            cols["wm"][c], cols["wp2"][c], cols["wp3"][c],
            cols["sigma_sqd_w"][c], cols["skw"][c], cols["skthl"][c],
            cols["skrt"][c], cols["rtm"][c], cols["rtp2"][c],
            cols["wprtp"][c], cols["thlm"][c], cols["thlp2"][c],
            cols["wpthlp"][c], cols["um"][c], cols["up2"][c],
            cols["upwp"][c], cols["vm"][c], cols["vp2"][c],
            cols["vpwp"][c], cols["rtpthlp"][c])
        assert perr == 0, (c, perr)
        moments[c] = m
        pdfp[c] = pp
        sig_out[c] = so
    out = {f"in_{k}": v for k, v in cols.items()}
    out["zi"] = zi
    out["zt"] = zt
    out["moments"] = moments
    out["pdf_params"] = pdfp
    out["sigma_sqd_w_out"] = sig_out
    return out


def main():
    params, idx = eamv3_params()
    sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    meta = dict(
        source_sha=sha,
        fflags="-O2 -fPIC -ffree-line-length-none "
               "-fallow-argument-mismatch -std=legacy -ffp-contract=off "
               "-DCLUBB_CAM -DCLUBB_SGS -DCLUBB_REAL_TYPE=dp",
        lapack="Ubuntu reference liblapack (container)",
        eamv3_overrides=EAMV3_OVERRIDES,
        param_indices={k: int(v) for k, v in idx.items()},
        config=dict(grid_type=3, l_implemented=True, sclr_dim=0,
                    hydromet_dim=0, edsclr_dim=29, theta0=300.0,
                    ts_nudge=86400.0, saturation_formula="flatau",
                    iiPDF_type="ADG1", l_stats=False, debug_level=0),
    )
    meta_s = json.dumps(meta)

    np.savez_compressed(GOLDEN / "clubb_grid.npz", meta=meta_s,
                        params=params, **gen_grid())
    print("wrote clubb_grid.npz")
    np.savez_compressed(GOLDEN / "clubb_tridag.npz", meta=meta_s,
                        **gen_tridag())
    print("wrote clubb_tridag.npz")
    np.savez_compressed(GOLDEN / "clubb_sat.npz", meta=meta_s,
                        **gen_sat())
    print("wrote clubb_sat.npz")
    np.savez_compressed(GOLDEN / "clubb_pdf_closure.npz", meta=meta_s,
                        params=params, **gen_pdf_closure(params, idx))
    print("wrote clubb_pdf_closure.npz")


if __name__ == "__main__":
    main()
