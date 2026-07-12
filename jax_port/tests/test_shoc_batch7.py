"""Tier-0 tests for SHOC batch 7: pblintd chain."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.shoc import (  # noqa: E402
    pblintd,
    shoc_pblintd_cldcheck,
    shoc_pblintd_init_pot,
)


def _column(ncol=3, nlev=32):
    zi = np.linspace(12000.0, 0.0, nlev + 1)[None, :] * np.ones((ncol, 1))
    zt = 0.5 * (zi[:, :-1] + zi[:, 1:])
    return zt, zi


def test_init_pot_matches_formula():
    rng = np.random.default_rng(0)
    thl = rng.uniform(280.0, 320.0, (2, 8))
    ql = rng.uniform(0.0, 2e-3, (2, 8))
    q = rng.uniform(1e-3, 2e-2, (2, 8))
    thv = np.asarray(shoc_pblintd_init_pot(thl, ql, q))
    th = thl + (c.LatVap / c.Cpair) * ql
    np.testing.assert_allclose(thv, th * (1 + c.ZVIR * q - ql), rtol=1e-14)


def test_cldcheck_bumps_pblh():
    pblh = np.array([100.0, 900.0])
    zi_low = np.array([300.0, 300.0])
    out = np.asarray(shoc_pblintd_cldcheck(zi_low, np.array([0.2, 0.2]), pblh))
    np.testing.assert_allclose(out, [350.0, 900.0])


def test_pblintd_neutral_vs_stable_columns():
    ncol, nlev = 3, 32
    zt, zi = _column(ncol, nlev)
    npbl = nlev  # allow the whole column

    # Column 0: strongly stable stratification -> Ri crosses 0.3 very low
    # Column 1: well-mixed (neutral thv) with strong wind shear -> Ri stays
    #           below critical in the mixed layer, deeper PBL
    thl = np.full((ncol, nlev), 300.0)
    thl[0] += np.linspace(30.0, 0.0, nlev)
    thl[1] = 300.0
    thl[2] += np.linspace(15.0, 0.0, nlev)
    ql = np.zeros((ncol, nlev))
    q = np.full((ncol, nlev), 5e-3)
    u = np.tile(np.linspace(20.0, 0.0, nlev), (ncol, 1))
    v = np.zeros((ncol, nlev))
    ustar = np.full(ncol, 0.3)
    obklen = np.full(ncol, -50.0)
    kbfs = np.full(ncol, -0.01)   # stable surface: no convective correction
    cldn = np.full((ncol, nlev), -1.0)  # disable cldcheck (cldn < 0)

    pblh = np.asarray(pblintd(zt, zi, thl, ql, q, u, v,
                              ustar, obklen, kbfs, cldn, npbl))
    # All bounded below by the mechanical minimum and above by column top
    assert np.all(pblh >= 700.0 * 0.3 - 1e-9)
    assert np.all(pblh <= zt[:, 0])
    # Stable column shallower than neutral sheared column
    assert pblh[0] < pblh[1]
    # Intermediate stratification: intermediate depth
    assert pblh[0] <= pblh[2] <= pblh[1] + 1e-9


def test_pblintd_convective_correction_deepens():
    ncol, nlev = 2, 32
    zt, zi = _column(ncol, nlev)
    npbl = nlev
    # Weakly stable profile, calm winds
    thl = 300.0 + np.linspace(3.0, 0.0, nlev)[None, :] * np.ones((ncol, 1))
    ql = np.zeros((ncol, nlev))
    q = np.full((ncol, nlev), 5e-3)
    u = np.zeros((ncol, nlev))
    v = np.zeros((ncol, nlev))
    ustar = np.full(ncol, 0.2)
    cldn = np.full((ncol, nlev), -1.0)

    # Column 0 stable surface (kbfs<0), column 1 convective (kbfs>0)
    kbfs = np.array([-0.01, 0.05])
    obklen = np.array([100.0, -30.0])
    pblh = np.asarray(pblintd(zt, zi, thl, ql, q, u, v,
                              ustar, obklen, kbfs, cldn, npbl))
    assert pblh[1] > pblh[0]  # convective correction deepens the PBL
