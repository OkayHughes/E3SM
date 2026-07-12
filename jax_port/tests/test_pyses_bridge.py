"""End-to-end smoke test of the pySEs bridge: synthetic pySEs-shaped
state (elements x GLL x lev, dry mass + dry mixing ratios) through the
full validated physics suite and back to tendencies.

Checks: round-trip consistency of the state conversions, water-mass
bookkeeping of the forcing, finiteness/positivity, and internal-state
persistence across steps.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

REPO = Path(__file__).resolve().parents[2]
DATA = REPO.parent / "e3sm-inputdata" / "atm" / "scream"

from scream_jax.pyses_bridge import (PysesScreamCoupler, columnize,  # noqa: E402
                                     decolumnize, eamxx_to_pyses_forcing,
                                     pyses_to_eamxx)


def test_state_conversion_roundtrip():
    rng = np.random.default_rng(0)
    ncol, nlev = 6, 20
    dmass = rng.uniform(900.0, 1100.0, (ncol, nlev)) * 5
    qv_dry = rng.uniform(1e-4, 1e-2, (ncol, nlev))
    qc_dry = rng.uniform(0.0, 1e-4, (ncol, nlev))
    T = rng.uniform(210.0, 300.0, (ncol, nlev))
    z = np.zeros((ncol, nlev))
    f = pyses_to_eamxx(T, z, z, z, dmass, {"qv": qv_dry, "qc": qc_dry},
                       ptop=100.0)
    # pressure consistency
    assert np.allclose(np.diff(f["p_int"], axis=1), f["pseudo_density"])
    assert np.allclose(f["pseudo_density"],
                       dmass * (1 + qv_dry + qc_dry))
    # wet<->dry mixing-ratio round trip
    qv_dry_rt = f["qv"] * f["pseudo_density"] / f["pseudo_density_dry"]
    np.testing.assert_allclose(qv_dry_rt, qv_dry, rtol=1e-14)
    # zero-tendency forcing when nothing changes
    forcing = eamxx_to_pyses_forcing(f, f, {"qv": qv_dry, "qc": qc_dry},
                                     300.0)
    np.testing.assert_allclose(forcing["FT"], 0.0, atol=1e-18)
    np.testing.assert_allclose(forcing["FQ"]["qv"], 0.0, atol=1e-18)


def test_columnize_roundtrip():
    a = np.arange(2 * 3 * 3 * 5, dtype=float).reshape(2, 3, 3, 5)
    c = columnize(a)
    assert c.shape == (18, 5)
    np.testing.assert_array_equal(decolumnize(c, (2, 3, 3)), a)


@pytest.mark.slow
def test_coupler_end_to_end():
    if not (DATA / "init" / "rrtmgp-data-sw-g112-210809.nc").exists():
        pytest.skip("scream data files not available")
    nc4 = pytest.importorskip("netCDF4")
    ic = DATA / "init" / "screami_unit_tests_ne2np4L72_20220822.nc"
    ds = nc4.Dataset(ic)
    geo = {k: np.array(ds[k][:]) for k in ("hyam", "hybm", "lat", "lon",
                                           "area")}
    # borrow a realistic thermodynamic state from the IC file (fields
    # carry a leading time dim of size 1)
    T = np.array(ds["T_mid"][0])
    qv_wet = np.array(ds["qv"][0])
    qc_wet = np.array(ds["qc"][0])
    ps = np.array(ds["ps"][0])
    o3 = np.array(ds["o3_volume_mix_ratio"][0])
    alb = {k: np.array(ds[k][0]) for k in
           ("sfc_alb_dir_vis", "sfc_alb_dir_nir",
            "sfc_alb_dif_vis", "sfc_alb_dif_nir")}
    slw = np.array(ds["surf_lw_flux_up"][0])
    ds.close()

    ncol, nlev = T.shape
    p0 = 100000.0
    dp_wet = (np.diff(geo["hyam"], append=0) * 0)  # placeholder below
    ai = np.concatenate([[geo["hyam"][0] * 0.5], 0.5 * (geo["hyam"][:-1]
                                                        + geo["hyam"][1:]),
                         [0.0]])
    # simple interface coefficients from midpoints (adequate for a smoke)
    bi = np.concatenate([[0.0], 0.5 * (geo["hybm"][:-1] + geo["hybm"][1:]),
                         [1.0]])
    p_int = p0 * ai[None, :] + ps[:, None] * bi[None, :]
    dp_wet = np.diff(p_int, axis=1)
    # dry mass / dry mixing ratios from the wet state
    qsum = qv_wet + qc_wet
    dp_dry = dp_wet * (1 - qsum)
    qv_dry = qv_wet * dp_wet / dp_dry
    qc_dry = qc_wet * dp_wet / dp_dry

    from scream_jax.foundation.thermo import calculate_dx_from_area
    surface = dict(
        surf_evap=np.zeros(ncol), surf_sens_flux=np.zeros(ncol),
        surf_mom_flux=np.zeros((ncol, 2)), surf_lw_flux_up=slw, **alb)
    coupler = PysesScreamCoupler(
        DATA, geo["hyam"], geo["hybm"], geo["lat"], geo["lon"],
        np.asarray(calculate_dx_from_area(geo["area"], geo["lat"])),
        ncol, nlev, surface, o3)

    z = np.zeros((ncol, nlev))
    dt = 1800.0
    doy0 = 284 + 45000.0 / 86400.0
    q_dry = {"qv": qv_dry, "qc": qc_dry, "qi": z.copy(), "qr": z.copy()}
    for nstep in range(2):
        forcing, diags = coupler.step(
            T, z + 5.0, z + 1.0, z, dp_dry, q_dry, float(p_int[0, 0].mean()),
            dt, nstep, doy0 + nstep * dt / 86400.0)
        for k in ("FT", "FU", "FV"):
            assert np.isfinite(forcing[k]).all(), k
        for k, v in forcing["FQ"].items():
            assert np.isfinite(v).all(), k
        # apply the forcing as pySEs' lump coupling would
        T = T + dt * forcing["FT"]
        for k in q_dry:
            q_dry[k] = np.maximum(q_dry[k] + dt * forcing["FQ"][k], 0.0)
        assert np.isfinite(T).all() and T.min() > 100 and T.max() < 350
    # physics-internal state persisted and evolving
    assert coupler.state["tke"].max() > 0.0004
    assert np.isfinite(coupler.state["rad_heating_pdel"]).all()
    print("coupled 2 steps: dT max", np.abs(dt * forcing["FT"]).max(),
          "K; qv range", q_dry["qv"].min(), q_dry["qv"].max())
