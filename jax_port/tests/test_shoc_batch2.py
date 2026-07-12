"""Tier-0 tests for SHOC batch 2: energy, obklen, vertflux/varorcovar."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.foundation.thermo import exner_function  # noqa: E402
from scream_jax.shoc import (  # noqa: E402
    calc_shoc_varorcovar,
    calc_shoc_vertflux,
    constants as sc,
    shoc_diag_obklen,
    shoc_energy_fixer,
    shoc_energy_integrals,
    update_host_dse,
)


def test_energy_integrals_analytic():
    ncol, nlev = 3, 8
    pdel = np.full((ncol, nlev), c.gravit)  # weight w = 1 per level
    ones = np.ones((ncol, nlev))
    se, ke, wv, wl = (np.asarray(a) for a in shoc_energy_integrals(
        2.0 * ones, pdel, 0.01 * ones, 0.004 * ones, 3.0 * ones, 4.0 * ones))
    np.testing.assert_allclose(se, 2.0 * nlev, rtol=1e-14)
    np.testing.assert_allclose(ke, 0.5 * 25.0 * nlev, rtol=1e-14)
    np.testing.assert_allclose(wv, 0.006 * nlev, rtol=1e-12)
    np.testing.assert_allclose(wl, 0.004 * nlev, rtol=1e-12)


def test_update_host_dse_formula():
    rng = np.random.default_rng(0)
    shape = (4, 16)
    thlm = rng.uniform(250.0, 320.0, shape)
    ql = rng.uniform(0.0, 2e-3, shape)
    inv_exner = rng.uniform(1.0, 2.5, shape)
    z = rng.uniform(10.0, 30000.0, shape)
    phis = rng.uniform(0.0, 3000.0, 4)
    dse = np.asarray(update_host_dse(thlm, ql, inv_exner, z, phis))
    T = thlm / inv_exner + (c.LatVap / c.CP) * ql
    np.testing.assert_allclose(dse, c.CP * T + c.gravit * z + phis[:, None], rtol=1e-14)


def _fixer_inputs(ncol=2, nlev=10, seed=1):
    rng = np.random.default_rng(seed)
    zi = np.linspace(20000.0, 0.0, nlev + 1)[None, :] * np.ones((ncol, 1))
    zt = 0.5 * (zi[:, :-1] + zi[:, 1:])
    rho = rng.uniform(0.3, 1.2, (ncol, nlev))
    pint = np.linspace(20000.0, 101325.0, nlev + 1)[None, :] * np.ones((ncol, 1))
    host_dse = rng.uniform(2.9e5, 3.2e5, (ncol, nlev))
    return zt, zi, rho, pint, host_dse


def test_energy_fixer_noop_when_balanced():
    zt, zi, rho, pint, host_dse = _fixer_inputs()
    ncol, nlev = host_dse.shape
    tke = np.full((ncol, nlev), 5.0)  # SHOC active everywhere -> shoctop=0
    z = np.zeros(ncol)
    out = np.asarray(shoc_energy_fixer(
        300.0, 6, zt, zi,
        1.0, 2.0, 3.0, 4.0,   # before
        1.0, 2.0, 3.0, 4.0,   # after (identical)
        z, z,                 # no surface fluxes
        rho, tke, pint, host_dse))
    np.testing.assert_allclose(out, host_dse, rtol=0, atol=1e-10)


def test_energy_fixer_spreads_disbalance_below_shoctop():
    zt, zi, rho, pint, host_dse = _fixer_inputs(seed=2)
    ncol, nlev = host_dse.shape
    ktop = 4
    tke = np.full((ncol, nlev), sc.mintke)
    tke[:, ktop:] = 1.0   # active only from k=ktop down
    z = np.zeros(ncol)
    se_a = 7.5e3  # after has extra energy -> se_dis>0 -> dse decreases
    out = np.asarray(shoc_energy_fixer(
        300.0, 6, zt, zi,
        0.0, 0.0, 0.0, 0.0,
        se_a, 0.0, 0.0, 0.0,
        z, z, rho, tke, pint, host_dse))
    # Untouched above shoctop
    np.testing.assert_array_equal(out[:, :ktop], host_dse[:, :ktop])
    # Uniform shift below
    se_dis = se_a / (pint[:, -1] - pint[:, ktop])
    np.testing.assert_allclose(out[:, ktop:],
                               host_dse[:, ktop:] - (se_dis * c.gravit)[:, None],
                               rtol=1e-12)


def test_energy_fixer_surface_fluxes_enter_te_b():
    # With wthl_sfc/wqw_sfc supplying exactly the after-minus-before energy,
    # the fixer must again be a no-op.
    zt, zi, rho, pint, host_dse = _fixer_inputs(seed=3)
    ncol, nlev = host_dse.shape
    tke = np.full((ncol, nlev), 5.0)
    dtime, nadv = 300.0, 6

    # rho_zi at surface equals rho_zt extrapolated; compute via the same interp
    from scream_jax.shoc import linear_interp
    rho_zi = np.asarray(linear_interp(zt, zi, rho, nlev, nlev + 1, 0.0))
    exner_sfc = np.asarray(exner_function(pint[:, -1]))
    wthl = np.full(ncol, 0.01)
    shf = wthl * c.CP * rho_zi[:, -1] * exner_sfc
    se_a = shf * dtime * nadv  # exactly the flux contribution

    out = np.asarray(shoc_energy_fixer(
        dtime, nadv, zt, zi,
        0.0, 0.0, 0.0, 0.0,
        se_a, 0.0, 0.0, 0.0,
        wthl, np.zeros(ncol), rho, tke, pint, host_dse))
    np.testing.assert_allclose(out, host_dse, rtol=0, atol=1e-9)


def test_diag_obklen_signs_and_floor():
    # ustar floor
    ustar, kbfs, obklen = (np.asarray(a) for a in shoc_diag_obklen(
        np.array(0.0), np.array(0.0), np.array(0.02), np.array(0.0),
        np.array(300.0), np.array(0.0), np.array(0.01)))
    assert ustar == sc.ustar_min
    # Unstable surface (kbfs > 0) => negative Obukhov length
    assert kbfs > 0 and obklen < 0
    # Stable surface (kbfs < 0) => positive Obukhov length
    _, kbfs2, obklen2 = (np.asarray(a) for a in shoc_diag_obklen(
        np.array(0.1), np.array(0.1), np.array(-0.02), np.array(0.0),
        np.array(300.0), np.array(0.0), np.array(0.01)))
    assert kbfs2 < 0 and obklen2 > 0


def test_vertflux_downgradient_and_band():
    ncol, nlev = 2, 6
    rng = np.random.default_rng(4)
    tkh_zi = rng.uniform(1.0, 50.0, (ncol, nlev + 1))
    dz_zi = rng.uniform(50.0, 500.0, (ncol, nlev + 1))
    sentinel = np.full((ncol, nlev + 1), -777.0)

    # Uniform scalar -> zero flux (interior); boundaries untouched.
    const = np.full((ncol, nlev), 3.5)
    out = np.asarray(calc_shoc_vertflux(tkh_zi, dz_zi, const, sentinel))
    np.testing.assert_array_equal(out[:, 0], -777.0)
    np.testing.assert_array_equal(out[:, -1], -777.0)
    np.testing.assert_allclose(out[:, 1:-1], 0.0, atol=1e-15)

    # invar increasing downward (k) -> invar(k-1)-invar(k) < 0 -> flux > 0
    inc = np.tile(np.arange(nlev, dtype=float), (ncol, 1))
    out2 = np.asarray(calc_shoc_vertflux(tkh_zi, dz_zi, inc, sentinel))
    assert np.all(out2[:, 1:-1] > 0)


def test_varorcovar_variance_nonnegative():
    ncol, nlev = 2, 6
    rng = np.random.default_rng(5)
    iso = rng.uniform(1.0, 1000.0, (ncol, nlev + 1))
    tkh = rng.uniform(1.0, 50.0, (ncol, nlev + 1))
    dz = rng.uniform(50.0, 500.0, (ncol, nlev + 1))
    x = rng.uniform(250.0, 320.0, (ncol, nlev))
    sentinel = np.full((ncol, nlev + 1), -777.0)
    var = np.asarray(calc_shoc_varorcovar(1.0, iso, tkh, dz, x, x, sentinel))
    assert np.all(var[:, 1:-1] >= 0)      # variance of x with itself
    np.testing.assert_array_equal(var[:, 0], -777.0)
    np.testing.assert_array_equal(var[:, -1], -777.0)
