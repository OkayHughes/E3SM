"""Tier-0 property tests for scream_jax.tms."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.tms import compute_tms  # noqa: E402


def _column_inputs(ncol=6, nlev=72, seed=0):
    rng = np.random.default_rng(seed)
    u = rng.uniform(-20, 20, (ncol, nlev))
    v = rng.uniform(-20, 20, (ncol, nlev))
    t = rng.uniform(230, 300, (ncol, nlev))
    p = np.linspace(2000, 101000, nlev)[None, :] * np.ones((ncol, 1))
    exner = (p / 1e5) ** (c.RD / c.CP)
    z = np.linspace(35000, 10.0, nlev)[None, :] * np.ones((ncol, 1))
    sgh = rng.uniform(50, 500, ncol)
    landfrac = rng.uniform(0, 1, ncol)
    return u, v, t, p, exner, z, sgh, landfrac


def test_zero_below_horomin_and_over_ocean():
    u, v, t, p, ex, z, sgh, lf = _column_inputs()
    sgh[0] = 0.5   # horo = orocnst*sgh < 1 m -> inactive
    lf[1] = 0.0    # ocean -> ksrf = 0 via landfrac factor
    ksrf, tau = (np.asarray(a) for a in compute_tms(u, v, t, p, ex, z, sgh, lf))
    assert ksrf[0] == 0.0 and np.all(tau[0] == 0.0)
    assert ksrf[1] == 0.0 and np.all(tau[1] == 0.0)
    assert np.all(ksrf >= 0.0)
    assert np.all(np.isfinite(ksrf)) and np.all(np.isfinite(tau))


def test_stress_opposes_bottom_wind():
    u, v, t, p, ex, z, sgh, lf = _column_inputs(seed=1)
    ksrf, tau = (np.asarray(a) for a in compute_tms(u, v, t, p, ex, z, sgh, lf))
    np.testing.assert_allclose(tau[:, 0], -ksrf * u[:, -1], rtol=1e-14)
    np.testing.assert_allclose(tau[:, 1], -ksrf * v[:, -1], rtol=1e-14)


def test_single_column_hand_value():
    # One neutral-ish column, checked against the formula chain by hand.
    nlev = 4
    u = np.array([[0.0, 0.0, 5.0, 10.0]])
    v = np.zeros((1, nlev))
    t = np.array([[280.0, 280.0, 280.0, 280.0]])
    p = np.array([[50000.0, 70000.0, 90000.0, 100000.0]])
    exner = np.ones((1, nlev))  # theta == T -> ri = 0 -> stabfri = 1
    z = np.array([[8000.0, 3000.0, 500.0, 50.0]])
    sgh = np.array([200.0])
    landfrac = np.array([1.0])

    ksrf, tau = (np.asarray(a) for a in
                 compute_tms(u, v, t, p, exner, z, sgh, landfrac))

    z0oro = min(c.z0fac * 200.0, 100.0)               # 15 m
    cd = (c.Karman / np.log((50.0 + z0oro) / z0oro)) ** 2
    rho = 100000.0 / (c.Rair * 280.0)
    expected_ksrf = rho * cd * 10.0
    np.testing.assert_allclose(ksrf[0], expected_ksrf, rtol=1e-13)
    np.testing.assert_allclose(tau[0, 0], -expected_ksrf * 10.0, rtol=1e-13)
    assert tau[0, 1] == 0.0


def test_stable_ri_kills_stress():
    # Strongly stable lowest layers (theta increasing sharply with height and
    # tiny shear) -> ri > 1 -> stabfri = 0 -> no stress.
    nlev = 3
    u = np.array([[1.0, 1.0, 1.001]])   # shear^2 << dv2min floor
    v = np.zeros((1, nlev))
    t = np.array([[320.0, 320.0, 250.0]])
    p = np.full((1, nlev), 90000.0)
    exner = np.ones((1, nlev))
    z = np.array([[2000.0, 800.0, 50.0]])
    ksrf, tau = (np.asarray(a) for a in
                 compute_tms(u, v, t, p, exner, z, np.array([300.0]), np.array([1.0])))
    assert ksrf[0] == 0.0
    assert np.all(tau == 0.0)
