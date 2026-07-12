"""Tier-0 tests for P3 sedimentation and homogeneous freezing.

Property tests: band finding, exact column mass closure (integrated mass
loss == accumulated surface precipitation), no-op behavior without
condensate, downward-only transport, and CFL substepping robustness.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.p3 import tables  # noqa: E402
from scream_jax.p3.sedimentation import (  # noqa: E402
    cloud_sedimentation,
    find_top_bottom,
    homogeneous_freezing,
    ice_sedimentation,
    rain_sedimentation,
)

REPO = Path(__file__).resolve().parents[2]
TDIR = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "tables"

NLEV = 72
OPTS = {
    "constant_mu_rain": 1.0,
    "min_rime_rho": 50.0,
    "max_rime_rho": 900.0,
    "ice_sedimentation_factor": 1.0,
}


@pytest.fixture(scope="module")
def tbl():
    if not (TDIR / f"p3_lookup_table_1.dat-v{tables.P3_VERSION}").exists():
        pytest.skip("P3 tables not available")
    return tables.p3_init(str(TDIR))


def _atmosphere(ncol=6, nlev=NLEV, seed=0):
    """Simple hydrostatic-ish test columns (k=0 top)."""
    rng = np.random.default_rng(seed)
    p = np.linspace(3000.0, 101325.0, nlev)[None, :] * np.ones((ncol, 1))
    T = 300.0 - 60.0 * (p[:, -1:] - p) / (p[:, -1:] - p[:, :1])
    rho = p / (287.0 * T)
    dz = np.diff(np.linspace(30000.0, 0.0, nlev + 1))[None, :] * np.ones((ncol, 1))
    dz = -dz  # decreasing z: positive thickness
    inv_dz = 1.0 / dz
    cld = np.clip(rng.uniform(0.2, 1.0, (ncol, nlev)), 1e-4, 1.0)
    return p, T, rho, 1.0 / rho, dz, inv_dz, cld


def _col_mass(q, rho, dz):
    return np.sum(np.asarray(q) * rho * dz, axis=-1)


def test_find_top_bottom():
    q = np.zeros((3, 10))
    q[0, 4] = 1e-8
    q[0, 7] = 1e-8
    q[2, 0] = 1e-8
    k_top, k_bot, present = find_top_bottom(q, c.QSMALL)
    assert present.tolist() == [True, False, True]
    assert k_top[0] == 4 and k_bot[0] == 7
    assert k_top[2] == 0 and k_bot[2] == 0


def test_cloud_sed_mass_closure_and_noop():
    p, T, rho, inv_rho, dz, inv_dz, cld = _atmosphere()
    ncol = rho.shape[0]
    rng = np.random.default_rng(1)

    qc = np.zeros((ncol, NLEV))
    qc[:, 40:60] = rng.uniform(1e-5, 5e-4, (ncol, 20))
    qc[0, :] = 0.0  # column 0: no condensate at all
    nc = np.where(qc > 0, 1e8 * inv_rho, 0.0)
    qc_incld = qc / cld
    nc_incld = nc / cld
    mu = 1.496e-6 * T ** 1.5 / (T + 120.0)
    acn = c.gravit * c.RHO_H2O / (18.0 * mu)

    dt, inv_dt = 300.0, 1.0 / 300.0
    zeros = np.zeros((ncol, NLEV))
    prt_in = np.full(ncol, 7.0)  # sentinel: must survive in empty columns

    out = cloud_sedimentation(qc_incld, rho, inv_rho, cld, acn, inv_dz,
                              dt, inv_dt, True, qc, nc, nc_incld,
                              zeros, zeros, qc, nc, prt_in)

    qc_new = np.asarray(out["qc"])
    prt = np.asarray(out["precip_liq_surf"])
    assert np.all(np.isfinite(qc_new)) and np.all(qc_new >= 0)

    # no-op column keeps state and the sentinel precip value
    assert np.array_equal(qc_new[0], qc[0])
    assert prt[0] == 7.0

    # exact mass closure: column mass loss == surface precip
    lost = _col_mass(qc, rho, dz)[1:] - _col_mass(qc_new, rho, dz)[1:]
    surf = prt[1:] * c.RHO_H2O * dt
    np.testing.assert_allclose(lost, surf, rtol=1e-12, atol=1e-14)
    assert np.all(prt[1:] >= 0)

    # transport is downward only: nothing appears above the initial top
    k_top = 40
    assert np.array_equal(qc_new[1:, :k_top], qc[1:, :k_top])

    # tendency definition
    np.testing.assert_allclose(np.asarray(out["qc_tend"]),
                               (qc_new - qc) * inv_dt, rtol=1e-12)


def test_cloud_sed_cfl_substepping():
    # thin layers force Co > 1 so multiple substeps run; closure must hold
    p, T, rho, inv_rho, dz, inv_dz, cld = _atmosphere(ncol=2)
    dz = dz * 0.02
    inv_dz = 1.0 / dz
    qc = np.zeros((2, NLEV))
    qc[:, 30:50] = 4e-4
    nc = 1e8 * inv_rho
    mu = 1.496e-6 * T ** 1.5 / (T + 120.0)
    acn = c.gravit * c.RHO_H2O / (18.0 * mu)
    dt, inv_dt = 1800.0, 1.0 / 1800.0
    zeros = np.zeros_like(qc)

    out = cloud_sedimentation(qc / cld, rho, inv_rho, cld, acn, inv_dz,
                              dt, inv_dt, True, qc, nc, nc / cld,
                              zeros, zeros, qc, nc, np.zeros(2))
    qc_new = np.asarray(out["qc"])
    lost = _col_mass(qc, rho, dz) - _col_mass(qc_new, rho, dz)
    surf = np.asarray(out["precip_liq_surf"]) * c.RHO_H2O * dt
    np.testing.assert_allclose(lost, surf, rtol=1e-12, atol=1e-16)
    assert np.all(qc_new >= 0) and np.all(np.isfinite(qc_new))


def test_rain_sed_mass_closure_and_flux(tbl):
    p, T, rho, inv_rho, dz, inv_dz, cld = _atmosphere(seed=2)
    ncol = rho.shape[0]
    rng = np.random.default_rng(3)

    qr = np.zeros((ncol, NLEV))
    qr[:, 45:71] = rng.uniform(1e-5, 2e-3, (ncol, 26))
    qr[0, :] = 0.0
    nr = np.where(qr > 0, 2e5 * inv_rho, 0.0)
    rhofacr = (c.RHO_1000MB * inv_rho) ** 0.54
    dt, inv_dt = 300.0, 1.0 / 300.0
    zeros = np.zeros((ncol, NLEV))
    flux_in = np.zeros((ncol, NLEV + 1))
    prt_in = np.full(ncol, 0.5)

    out = rain_sedimentation(rho, inv_rho, rhofacr, cld, inv_dz, qr / cld,
                             tbl["vn_table_vals"], tbl["vm_table_vals"],
                             dt, inv_dt, qr, nr, nr / cld, zeros, zeros,
                             flux_in, qr, nr, prt_in, OPTS)

    qr_new = np.asarray(out["qr"])
    prt = np.asarray(out["precip_liq_surf"])
    flux = np.asarray(out["precip_liq_flux"])
    assert np.all(np.isfinite(qr_new)) and np.all(qr_new >= 0)

    # empty column: state untouched, precip accumulates nothing
    assert np.array_equal(qr_new[0], qr[0])
    assert prt[0] == 0.5
    assert np.all(flux[0] == 0.0)

    # accumulation semantics: contribution = prt - prt_in
    # closure to FP-accumulation level (rho*inv_rho, dz*inv_dz are not
    # exactly 1, and many CFL substeps accumulate the rounding)
    lost = _col_mass(qr, rho, dz)[1:] - _col_mass(qr_new, rho, dz)[1:]
    surf = (prt[1:] - 0.5) * c.RHO_H2O * dt
    np.testing.assert_allclose(lost, surf, rtol=1e-9, atol=1e-12)
    assert np.all(flux >= 0.0)
    # flux only below the rain top
    assert np.all(flux[1:, :45] == 0.0)


def test_ice_sed_mass_closure(tbl):
    p, T, rho, inv_rho, dz, inv_dz, cld = _atmosphere(seed=4)
    ncol = rho.shape[0]
    rng = np.random.default_rng(5)

    qi = np.zeros((ncol, NLEV))
    qi[:, 20:55] = rng.uniform(1e-6, 5e-4, (ncol, 35))
    qi[0, :] = 0.0
    ni = np.where(qi > 0, 1e5 * inv_rho, 0.0)
    qm = 0.3 * qi
    bm = qm / 400.0
    rhofaci = (c.RHO_600MB * inv_rho) ** 0.54
    dt, inv_dt = 300.0, 1.0 / 300.0
    zeros = np.zeros((ncol, NLEV))

    out = ice_sedimentation(rho, inv_rho, rhofaci, cld, inv_dz, dt, inv_dt,
                            qi, qi / cld, ni, ni / cld, qm, qm / cld,
                            bm, bm / cld, tbl["ice_table_vals"],
                            qi, ni, np.zeros(ncol), OPTS)

    qi_new = np.asarray(out["qi"])
    prt = np.asarray(out["precip_ice_surf"])
    assert np.all(np.isfinite(qi_new)) and np.all(qi_new >= 0)
    assert np.array_equal(qi_new[0], qi[0]) and prt[0] == 0.0

    lost = _col_mass(qi, rho, dz)[1:] - _col_mass(qi_new, rho, dz)[1:]
    surf = prt[1:] * c.RHO_H2O * dt
    np.testing.assert_allclose(lost, surf, rtol=1e-9, atol=1e-12)

    # rime mass stays bounded by total ice mass where ice remains
    qm_new = np.asarray(out["qm"])
    assert np.all(qm_new <= np.asarray(out["qi"]) + 1e-12)


def test_homogeneous_freezing():
    ncol, nlev = 3, 8
    T = np.full((ncol, nlev), 250.0)
    T[0, :] = 220.0  # below T_homogfrz = 233.15
    inv_exner = np.full((ncol, nlev), 1.2)
    qc = np.full((ncol, nlev), 1e-4)
    nc = np.full((ncol, nlev), 1e7)
    qr = np.full((ncol, nlev), 5e-5)
    nr = np.full((ncol, nlev), 1e5)
    qi = np.full((ncol, nlev), 1e-5)
    ni = np.full((ncol, nlev), 1e4)
    qm = np.full((ncol, nlev), 2e-6)
    bm = qm / 400.0
    th = np.full((ncol, nlev), 300.0)

    out = homogeneous_freezing(T, inv_exner, qc, nc, qr, nr, qi, ni, qm, bm, th)

    # warm columns untouched
    for k in ("qc", "nc", "qr", "nr", "qi", "ni", "qm", "bm", "th_atm"):
        assert np.array_equal(np.asarray(out[k])[1:], locals()[
            {"th_atm": "th"}.get(k, k)][1:])

    # cold column: all liquid converts to ice, ni gains the drop numbers
    assert np.all(np.asarray(out["qc"])[0] == 0.0)
    assert np.all(np.asarray(out["qr"])[0] == 0.0)
    np.testing.assert_allclose(np.asarray(out["qi"])[0],
                               qi[0] + qc[0] + qr[0], rtol=1e-15)
    np.testing.assert_allclose(np.asarray(out["qm"])[0],
                               qm[0] + qc[0] + qr[0], rtol=1e-15)
    np.testing.assert_allclose(np.asarray(out["ni"])[0],
                               ni[0] + nc[0] + nr[0], rtol=1e-15)
    # latent heating warms theta
    dth = np.asarray(out["th_atm"])[0] - th[0]
    # two sequential updates (qc then qr) vs one combined product: FP assoc.
    np.testing.assert_allclose(
        dth, inv_exner[0] * (qc[0] + qr[0]) * c.LatIce * c.INV_CP, rtol=1e-12)
