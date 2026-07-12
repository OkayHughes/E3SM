"""Tier-0 tests for SHOC batch 4: implicit vertical diffusion."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.shoc import (  # noqa: E402
    update_prognostics_implicit,
    vd_shoc_decomp,
    vd_shoc_solve,
)


def _setup(ncol=3, nlev=24, nq=5, seed=0):
    rng = np.random.default_rng(seed)
    zi = np.linspace(12000.0, 0.0, nlev + 1)[None, :] * np.ones((ncol, 1))
    zt = 0.5 * (zi[:, :-1] + zi[:, 1:])
    dz_zt = zi[:, :-1] - zi[:, 1:]
    dz_zi = np.concatenate(
        [np.full((ncol, 1), 1.0), zt[:, :-1] - zt[:, 1:], zt[:, -1:]], axis=-1)
    rho = np.linspace(0.3, 1.2, nlev)[None, :] * np.ones((ncol, 1))
    tk = rng.uniform(0.1, 20.0, (ncol, nlev))
    tkh = rng.uniform(0.1, 20.0, (ncol, nlev))
    return rng, zt, zi, dz_zt, dz_zi, rho, tk, tkh


def test_vd_decomp_structure():
    rng, zt, zi, dz_zt, dz_zi, rho, tk, tkh = _setup()
    ncol, nlev = zt.shape
    kv = rng.uniform(0.1, 20.0, (ncol, nlev + 1))
    tmpi = rng.uniform(0.0, 5.0, (ncol, nlev + 1))
    rdp = rng.uniform(1e-4, 1e-2, (ncol, nlev))
    dl, d, du = (np.asarray(a) for a in vd_shoc_decomp(kv, tmpi, rdp, 300.0, 0.0))
    np.testing.assert_array_equal(dl[:, 0], 0.0)
    np.testing.assert_array_equal(du[:, -1], 0.0)
    assert np.all(dl <= 0) and np.all(du <= 0)
    np.testing.assert_allclose(d, 1.0 - du - dl, rtol=1e-14)  # flux=0 case
    # Rows sum to 1 => uniform profiles are steady states of the solve.
    np.testing.assert_allclose(dl + d + du, 1.0, rtol=1e-13)
    # Diagonal dominance -> well-posed Thomas solve
    assert np.all(d >= 1.0)


def test_vd_solve_preserves_uniform_profile():
    rng, zt, zi, dz_zt, dz_zi, rho, tk, tkh = _setup(seed=1)
    ncol, nlev = zt.shape
    kv = rng.uniform(0.1, 20.0, (ncol, nlev + 1))
    tmpi = rng.uniform(0.0, 5.0, (ncol, nlev + 1))
    rdp = rng.uniform(1e-4, 1e-2, (ncol, nlev))
    dl, d, du = vd_shoc_decomp(kv, tmpi, rdp, 300.0, 0.0)
    x = np.full((ncol, nlev), 7.25)
    sol = np.asarray(vd_shoc_solve(dl, d, du, x))
    np.testing.assert_allclose(sol, 7.25, rtol=1e-12)


def test_update_prognostics_implicit_conservation_and_smoothing():
    rng, zt, zi, dz_zt, dz_zi, rho, tk, tkh = _setup(seed=2)
    ncol, nlev = zt.shape
    nq = 4
    zeros = np.zeros(ncol)

    thetal = 300.0 + rng.uniform(-5, 5, (ncol, nlev))
    qw = rng.uniform(1e-3, 2e-2, (ncol, nlev))
    qtr = rng.uniform(0.0, 1e-3, (ncol, nq, nlev))
    tke = rng.uniform(0.01, 0.5, (ncol, nlev))
    u = rng.uniform(-20, 20, (ncol, nlev))
    v = rng.uniform(-20, 20, (ncol, nlev))

    # No surface fluxes and no drag floor effects on scalars: with
    # wthl=wqw=wtracer=0 the thermo solve conserves the mass-weighted
    # integral of thetal/qw/tracers (pure diffusion with zero-flux BCs),
    # apart from the explicit tke surface source.
    out = update_prognostics_implicit(
        300.0, dz_zt, dz_zi, rho, zt, zi, tk, tkh,
        zeros, zeros, zeros, zeros, np.zeros((ncol, nq)),
        thetal.copy(), qw.copy(), qtr.copy(), tke.copy(), u.copy(), v.copy())
    thetal2, qw2, qtr2, tke2, u2, v2 = (np.asarray(a) for a in out)

    pdel = c.gravit * rho * dz_zt  # dp = g rho dz
    for before, after in ((thetal, thetal2), (qw, qw2)):
        np.testing.assert_allclose((after * pdel).sum(-1), (before * pdel).sum(-1),
                                   rtol=1e-11)
    np.testing.assert_allclose((qtr2 * pdel[:, None, :]).sum(-1),
                               (qtr * pdel[:, None, :]).sum(-1), rtol=1e-11)

    # Diffusion must not create new extrema in the scalar profiles.
    assert thetal2.max() <= thetal.max() + 1e-9
    assert thetal2.min() >= thetal.min() - 1e-9

    # Winds: surface drag (ksrf floor) removes momentum; KE must not grow.
    ke_before = (0.5 * (u ** 2 + v ** 2) * pdel).sum(-1)
    ke_after = (0.5 * (u2 ** 2 + v2 ** 2) * pdel).sum(-1)
    assert np.all(ke_after <= ke_before + 1e-9)


def test_update_prognostics_surface_fluxes_add_mass():
    rng, zt, zi, dz_zt, dz_zi, rho, tk, tkh = _setup(seed=3)
    ncol, nlev = zt.shape
    zeros = np.zeros(ncol)
    qw = np.full((ncol, nlev), 1e-2)
    wqw = np.full(ncol, 1e-4)

    out = update_prognostics_implicit(
        300.0, dz_zt, dz_zi, rho, zt, zi, tk, tkh,
        zeros, zeros, zeros, wqw, np.zeros((ncol, 0)),
        np.full((ncol, nlev), 300.0), qw.copy(),
        np.zeros((ncol, 0, nlev)), np.full((ncol, nlev), 0.01),
        np.zeros((ncol, nlev)), np.zeros((ncol, nlev)))
    qw2 = np.asarray(out[1])

    # Explicit surface flux adds exactly cmnfac*wqw at the lowest level,
    # then diffusion conserves it: column gain = dtime * g * rho_sfc * wqw
    # / (g rho dz)|sfc * pdel|sfc  =  dtime * wqw * g * rho_zi_sfc.
    pdel = c.gravit * rho * dz_zt
    gain = (qw2 * pdel).sum(-1) - (qw * pdel).sum(-1)
    from scream_jax.shoc import linear_interp
    rho_zi = np.asarray(linear_interp(zt, zi, rho, nlev, nlev + 1, 0.0))
    expected = 300.0 * c.gravit * rho_zi[:, -1] * wqw
    np.testing.assert_allclose(gain, expected, rtol=1e-10)
