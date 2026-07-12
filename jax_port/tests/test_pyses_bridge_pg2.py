"""End-to-end smoke of the pg2 bridge: element-shaped pySEs-style state
through dyn_to_fv_phys -> full physics suite on FV columns ->
fv_phys_to_dyn, for two coupled steps.

The horizontal grid is synthetic (96 IC columns reshaped into 6 fake
np4 elements, identity map tensors, IC areas as the metric), which is
fine for a smoke: the physics is column-local and the remap operators
only need a positive metric. Checks: finiteness of all forcings, the
limiter's positivity guarantee on updated tracers, physical T range
after applying the forcing, and tracer-mass bookkeeping between the FV
physics increment and the GLL update.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

PYSES_ROOT = Path(os.environ.get(
    "PYSES_PATH", "/Users/ostensiblyowen/development/python/pyses_07_04_26"))
if not (PYSES_ROOT / "pyses").is_dir():
    pytest.skip("pySEs source tree not available", allow_module_level=True)
sys.path.insert(0, str(PYSES_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

REPO = Path(__file__).resolve().parents[2]
DATA = REPO.parent / "e3sm-inputdata" / "atm" / "scream"

from pyses.dynamical_cores.finite_volume_grid import init_fv_grid  # noqa: E402

from pyses_ext import finite_volume_grid_operational as opfv  # noqa: E402
from scream_jax.pyses_bridge import (PysesScreamCouplerPg2,  # noqa: E402
                                     columnize)

NELEM, NPT, NF = 6, 4, 2
NCOL_GLL = NELEM * NPT * NPT     # 96 IC columns reshaped into elements
NCOL_FV = NELEM * NF * NF


def elem(a):
    """(96, ...) IC columns -> (6, 4, 4, ...) fake elements."""
    return np.asarray(a)[:NCOL_GLL].reshape((NELEM, NPT, NPT)
                                            + np.asarray(a).shape[1:])


@pytest.mark.slow
def test_pg2_coupler_end_to_end():
    if not (DATA / "init" / "rrtmgp-data-sw-g112-210809.nc").exists():
        pytest.skip("scream data files not available")
    nc4 = pytest.importorskip("netCDF4")
    ic = DATA / "init" / "screami_unit_tests_ne2np4L72_20220822.nc"
    ds = nc4.Dataset(ic)
    hyam, hybm = np.array(ds["hyam"][:]), np.array(ds["hybm"][:])
    lat, lon, area = (np.array(ds[k][:]) for k in ("lat", "lon", "area"))
    T = elem(np.array(ds["T_mid"][0]))
    qv_wet = elem(np.array(ds["qv"][0]))
    qc_wet = elem(np.array(ds["qc"][0]))
    ps = elem(np.array(ds["ps"][0]))
    o3 = elem(np.array(ds["o3_volume_mix_ratio"][0]))
    alb = {k: elem(np.array(ds[k][0])) for k in
           ("sfc_alb_dir_vis", "sfc_alb_dir_nir",
            "sfc_alb_dif_vis", "sfc_alb_dif_nir")}
    slw = elem(np.array(ds["surf_lw_flux_up"][0]))
    ds.close()

    # --- synthetic element grid over the reshaped IC columns ---
    eye = np.broadcast_to(np.eye(2), (NELEM, NPT, NPT, 2, 2)).copy()
    h_grid = {
        "metric_determinant": elem(area),
        "physical_coords": np.stack([elem(lat), elem(lon)], axis=-1),
        "contra_to_physical": eye,
        "physical_to_contra": eye.copy(),
    }
    fv_grid = init_fv_grid(h_grid, {"npt": NPT, "num_elem": NELEM}, nf=NF)
    fv_grid = opfv.extend_fv_grid_operational(fv_grid, h_grid)

    p0 = 100000.0
    ai = np.concatenate([[hyam[0] * 0.5], 0.5 * (hyam[:-1] + hyam[1:]),
                         [0.0]])
    bi = np.concatenate([[0.0], 0.5 * (hybm[:-1] + hybm[1:]), [1.0]])
    v_grid = {
        "hybrid_a_i": ai, "hybrid_b_i": bi,
        "hybrid_a_m": hyam, "hybrid_b_m": hybm,
        "reference_surface_mass": p0,
    }

    # dry state from the wet IC (same conversion as the np4 smoke test)
    p_int = p0 * ai[None, None, None] + ps[..., None] * bi[None, None, None]
    dp_wet = np.diff(p_int, axis=-1)
    qsum = qv_wet + qc_wet
    dp_dry = dp_wet * (1 - qsum)
    q_dry = {"qv": qv_wet * dp_wet / dp_dry,
             "qc": qc_wet * dp_wet / dp_dry,
             "qi": np.zeros_like(dp_dry), "qr": np.zeros_like(dp_dry)}
    # remap uses ptop + sum(dp_dry) for ps; consistent by construction of ai

    # FV-column geometry/forcing data via the remap operators
    from pyses.dynamical_cores.finite_volume_grid import gll_to_fv
    from scream_jax.foundation.thermo import calculate_dx_from_area
    lat_fv = columnize(np.asarray(gll_to_fv(elem(lat), fv_grid)))
    lon_fv = columnize(np.asarray(gll_to_fv(elem(lon), fv_grid)))
    area_fv = np.repeat(elem(area).sum(axis=(1, 2)) / (NF * NF), NF * NF)
    o3_fv = columnize(np.asarray(gll_to_fv(o3, fv_grid)))
    surface = dict(
        surf_evap=np.zeros(NCOL_FV), surf_sens_flux=np.zeros(NCOL_FV),
        surf_mom_flux=np.zeros((NCOL_FV, 2)),
        surf_lw_flux_up=columnize(np.asarray(gll_to_fv(slw, fv_grid))),
        **{k: columnize(np.asarray(gll_to_fv(v, fv_grid)))
           for k, v in alb.items()})

    # SPA data is per-column on the np4 grid; give each FV cell a
    # representative source column (its lower-left GLL point)
    spa_idx = np.array([e * NPT * NPT + (2 * c) * NPT + 2 * d
                        for e in range(NELEM)
                        for c in range(NF) for d in range(NF)])
    coupler = PysesScreamCouplerPg2(
        DATA, h_grid, fv_grid, v_grid, lat_fv, lon_fv,
        np.asarray(calculate_dx_from_area(area_fv, lat_fv)),
        surface, o3_fv, spa_col_indices=spa_idx)

    hw = np.zeros(T.shape + (2,))
    hw[..., 0] = 5.0
    hw[..., 1] = 1.0
    omega = np.zeros_like(T)
    dt = 1800.0
    doy0 = 284 + 45000.0 / 86400.0

    spheremp = np.asarray(fv_grid["spheremp"])
    spheremp_fv = np.asarray(fv_grid["spheremp_fv"])
    for nstep in range(2):
        forcing, diags = coupler.step(T, hw, omega, dp_dry, q_dry, dt,
                                      nstep, doy0 + nstep * dt / 86400.0)
        for k in ("FT", "FU", "FV"):
            assert np.isfinite(forcing[k]).all(), k
            assert forcing[k].shape == T.shape, k
        for k, v in forcing["FQ"].items():
            assert np.isfinite(v).all(), k

        # tracer-mass bookkeeping: the GLL update carries exactly the
        # FV-grid physics mass increment (the limiter conserves it)
        dp_fv = np.asarray(opfv.calc_dp_fv(
            gll_to_fv(coupler.ptop + dp_dry.sum(-1), fv_grid), v_grid))
        q0_fv = np.asarray(opfv.gll_to_fv_mixing_ratio(
            q_dry["qv"], dp_dry, dp_fv, fv_grid))
        # post-physics DRY qv on FV, reconstructed the way the bridge
        # does (diags qv is wet; the wet dp is implied by the new water)
        tot_q_wet = sum(diags[s] for s in ("qv", "qc", "qi", "qr", "qm")
                        if s in diags)
        dp_wet_new = diags["pseudo_density_dry"] / (1.0 - tot_q_wet)
        qv1_dry_fv = (diags["qv"] * dp_wet_new
                      / diags["pseudo_density_dry"]).reshape(
                          NELEM, NF, NF, -1)
        dq_fv_mass = np.sum(dp_fv * (qv1_dry_fv - q0_fv)
                            * spheremp_fv[..., None])
        dq_gll_mass = np.sum(dp_dry * dt * forcing["FQ"]["qv"]
                             * spheremp[..., None])
        mass_scale = np.sum(dp_dry * q_dry["qv"] * spheremp[..., None])
        assert abs(dq_gll_mass - dq_fv_mass) < 1e-11 * mass_scale, (
            dq_gll_mass, dq_fv_mass)

        # apply the forcing as pySEs' lump coupling would
        T = T + dt * forcing["FT"]
        for k in q_dry:
            q_dry[k] = q_dry[k] + dt * forcing["FQ"][k]
            # limiter guarantee: updated state within [min(q0, q1_fv), ...]
            # and in particular non-negative
            assert q_dry[k].min() >= -1e-12, k
            q_dry[k] = np.maximum(q_dry[k], 0.0)
        assert np.isfinite(T).all() and T.min() > 150 and T.max() < 350

    assert coupler.inner.state["tke"].max() > 0.0004
    print("pg2 coupled 2 steps: dT max", np.abs(dt * forcing["FT"]).max(),
          "K; qv range", q_dry["qv"].min(), q_dry["qv"].max())
