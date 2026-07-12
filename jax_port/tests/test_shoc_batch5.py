"""Tier-0 tests for SHOC batch 5: second-moment driver chain."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.shoc import (  # noqa: E402
    diag_second_shoc_moments,
    shoc_diag_second_moments_lbycond,
    shoc_diag_second_moments_srf,
)


def _inputs(ncol=3, nlev=24, seed=0):
    rng = np.random.default_rng(seed)
    zi = np.linspace(12000.0, 0.0, nlev + 1)[None, :] * np.ones((ncol, 1))
    zt = 0.5 * (zi[:, :-1] + zi[:, 1:])
    dz_zi = np.concatenate(
        [np.full((ncol, 1), 1.0), zt[:, :-1] - zt[:, 1:], zt[:, -1:]], axis=-1)
    return dict(
        rng=rng, zt=zt, zi=zi, dz_zi=dz_zi,
        thetal=300.0 + rng.uniform(-5, 5, (ncol, nlev)),
        qw=rng.uniform(1e-3, 2e-2, (ncol, nlev)),
        u=rng.uniform(-20, 20, (ncol, nlev)),
        v=rng.uniform(-20, 20, (ncol, nlev)),
        tke=rng.uniform(0.01, 0.5, (ncol, nlev)),
        iso=rng.uniform(1.0, 1000.0, (ncol, nlev)),
        tkh=rng.uniform(0.1, 20.0, (ncol, nlev)),
        tk=rng.uniform(0.1, 20.0, (ncol, nlev)),
        mix=rng.uniform(20.0, 1000.0, (ncol, nlev)))


def test_srf_and_lbycond():
    ustar2, wstar = (np.asarray(a) for a in shoc_diag_second_moments_srf(
        np.array([0.02, -0.02]), np.array([0.1, 0.0]), np.array([0.0, 0.0])))
    np.testing.assert_allclose(ustar2, [0.1, 0.0], rtol=1e-14)
    # wstar zero for downward heat flux; positive and cube-rooted otherwise
    assert wstar[1] == 0.0
    np.testing.assert_allclose(
        wstar[0], (c.gravit * 0.02 / c.basetemp) ** (1.0 / 3.0), rtol=1e-13)

    out = [np.asarray(a) for a in shoc_diag_second_moments_lbycond(
        np.array(0.02), np.array(1e-4), np.array(0.1), np.array(0.05),
        np.asarray(ustar2[0]), np.asarray(wstar[0]))]
    wthl_b, wqw_b, uw_b, vw_b, wtke_b, thl_b, qw_b, qwthl_b = out
    # Fluxes pass through
    assert wthl_b == 0.02 and wqw_b == 1e-4 and uw_b == 0.1 and vw_b == 0.05
    # Variances non-negative; covariance sign = product of flux signs
    assert thl_b >= 0 and qw_b >= 0 and qwthl_b >= 0 and wtke_b > 0


def test_diag_second_shoc_moments_structure():
    d = _inputs()
    ncol, nlev = d["thetal"].shape
    wthl = np.full(ncol, 0.02)
    wqw = np.full(ncol, 1e-4)
    uw = np.full(ncol, 0.1)
    vw = np.full(ncol, -0.05)

    (thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec,
     wtke_sec, w_sec, ustar2, wstar) = (np.asarray(a) for a in
        diag_second_shoc_moments(
            1.0, 1.0, 1.0, 1.0, False,
            d["thetal"], d["qw"], d["u"], d["v"], d["tke"], d["iso"],
            d["tkh"], d["tk"], d["dz_zi"], d["zt"], d["zi"], d["mix"],
            wthl, wqw, uw, vw))

    # Shapes: interfaces for *_sec, midpoints for w_sec
    assert thl_sec.shape == (ncol, nlev + 1) and w_sec.shape == (ncol, nlev)
    # Upper bc: all zero at k=0
    for a in (thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec):
        np.testing.assert_array_equal(a[:, 0], 0.0)
    # Lower bc: fluxes equal surface fluxes
    np.testing.assert_allclose(wthl_sec[:, -1], wthl, rtol=1e-14)
    np.testing.assert_allclose(wqw_sec[:, -1], wqw, rtol=1e-14)
    np.testing.assert_allclose(uw_sec[:, -1], uw, rtol=1e-14)
    np.testing.assert_allclose(vw_sec[:, -1], vw, rtol=1e-14)
    # Variances non-negative everywhere
    assert np.all(thl_sec >= 0) and np.all(qw_sec >= 0)
    # w_sec = (2/3) tke with w2tune=1
    np.testing.assert_allclose(w_sec, (2.0 / 3.0) * d["tke"], rtol=1e-14)
    # Cauchy-Schwarz-like consistency: |qwthl| <= sqrt(thl_sec qw_sec)
    interior = slice(1, nlev)
    cs = np.sqrt(thl_sec[:, interior] * qw_sec[:, interior])
    assert np.all(np.abs(qwthl_sec[:, interior]) <= cs + 1e-12)


def test_diag_second_shoc_moments_1p5tke_zeroes_variances():
    d = _inputs(seed=1)
    ncol, nlev = d["thetal"].shape
    args = (1.0, 1.0, 1.0, 1.0, True,
            d["thetal"], d["qw"], d["u"], d["v"], d["tke"], d["iso"],
            d["tkh"], d["tk"], d["dz_zi"], d["zt"], d["zi"], d["mix"],
            np.full(ncol, 0.02), np.full(ncol, 1e-4),
            np.full(ncol, 0.1), np.full(ncol, 0.0))
    out = [np.asarray(a) for a in diag_second_shoc_moments(*args)]
    thl_sec, qw_sec, _, _, qwthl_sec = out[0], out[1], out[2], out[3], out[4]
    w_sec = out[8]
    np.testing.assert_array_equal(thl_sec, 0.0)
    np.testing.assert_array_equal(qw_sec, 0.0)
    np.testing.assert_array_equal(qwthl_sec, 0.0)
    np.testing.assert_array_equal(w_sec, 0.0)
