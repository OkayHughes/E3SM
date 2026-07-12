"""Tier-0 property tests for scream_jax.foundation.

These are sanity/invariant checks only — the authoritative validation is the
golden-data comparison against EAMxx (not yet wired; see STATUS.md). Values
in the spot checks below are physical ballparks, not EAMxx golden data.

Run:  pytest jax_port/tests/test_foundation.py   (needs: pip install jax pytest)
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import (  # noqa: E402
    SaturationFcn,
    constants as c,
    murphy_koop_svp,
    polysvp1,
    qv_sat_dry,
    qv_sat_wet,
)


def test_constants_derived_values():
    # Derived constants must reproduce the C++ arithmetic exactly.
    assert c.ep_2 == pytest.approx(18.016 / 28.966, rel=0, abs=0)
    assert c.RHOSUR == pytest.approx(100000.0 / (287.042 * 273.15), rel=0, abs=0)
    assert c.ZVIR == pytest.approx((c.Rgas / c.MWWV) / 287.042 - 1.0, rel=0, abs=0)
    # Expressions, not decimal literals: the C++ computes Tmelt - 40 in
    # binary FP (233.14999999999998), and we must match that bit pattern.
    assert c.T_homogfrz == 273.15 - 40.0
    assert c.T_rainfrz == 273.15 - 4.0


@pytest.mark.parametrize("svp", [polysvp1, murphy_koop_svp])
def test_svp_at_triple_point(svp):
    # Both formulations give ~611 Pa at 273.15/273.16 K, ice and liquid alike.
    for ice in (False, True):
        es = float(svp(np.array(273.15), ice))
        assert 605.0 < es < 615.0


@pytest.mark.parametrize("svp", [polysvp1, murphy_koop_svp])
@pytest.mark.parametrize("ice", [False, True])
def test_svp_monotonic_in_t(svp, ice):
    # Start above 193.15 K: polysvp1 clamps dt = max(T-273.15, -80), so es is
    # deliberately constant below that (matching the C++).
    t = np.linspace(195.0, 330.0, 271)
    es = np.asarray(svp(t, ice))
    assert np.all(np.diff(es) > 0)
    assert np.all(np.isfinite(es))


def test_polysvp1_clamped_below_193K():
    # The Flatau polynomial is evaluated at dt = -80 for all T < 193.15 K.
    lo = np.asarray(polysvp1(np.array([150.0, 180.0, 193.0]), True))
    np.testing.assert_allclose(lo, lo[-1], rtol=0, atol=0)


def test_ice_supersaturation_sign():
    # Below freezing, es over ice must be lower than es over liquid.
    t = np.linspace(200.0, 272.0, 73)
    for svp in (polysvp1, murphy_koop_svp):
        assert np.all(np.asarray(svp(t, True)) < np.asarray(svp(t, False)))
    # At/above freezing the ice flag must be inert.
    t_warm = np.linspace(273.15, 320.0, 47)
    for svp in (polysvp1, murphy_koop_svp):
        np.testing.assert_array_equal(np.asarray(svp(t_warm, True)),
                                      np.asarray(svp(t_warm, False)))


@pytest.mark.parametrize("func", [SaturationFcn.POLYSVP1, SaturationFcn.MURPHY_KOOP])
def test_qv_sat_dry_definition(func):
    # qv_sat_dry must equal ep_2 * es / p_dry for ordinary inputs.
    t, p = np.array(285.0), np.array(90000.0)
    es = polysvp1(t, False) if func == SaturationFcn.POLYSVP1 else murphy_koop_svp(t, False)
    expected = c.ep_2 * float(es) / 90000.0
    assert float(qv_sat_dry(t, p, False, func)) == pytest.approx(expected, rel=1e-15)


def test_qv_sat_wet_scaling():
    # qsat_wet = qsat_dry * dp_dry/dp_wet, elementwise.
    t = np.array([250.0, 285.0, 300.0])
    p = np.array([50000.0, 85000.0, 100000.0])
    dp_wet = np.array([800.0, 900.0, 1000.0])
    dp_dry = np.array([790.0, 880.0, 985.0])
    qd = np.asarray(qv_sat_dry(t, p, True))
    qw = np.asarray(qv_sat_wet(t, p, True, dp_wet, dp_dry))
    np.testing.assert_allclose(qw, qd * dp_dry / dp_wet, rtol=1e-15)


def test_float64_default():
    out = np.asarray(polysvp1(np.array(273.15), False))
    assert out.dtype == np.float64


# --- thermo + column_ops ------------------------------------------------------

from scream_jax.foundation import column_ops, thermo  # noqa: E402


def test_theta_T_round_trip():
    rng = np.random.default_rng(0)
    T = rng.uniform(180.0, 320.0, (4, 72))
    p = rng.uniform(2000.0, 103000.0, (4, 72))
    theta = np.asarray(thermo.calculate_theta_from_T(T, p))
    np.testing.assert_allclose(np.asarray(thermo.calculate_T_from_theta(theta, p)),
                               T, rtol=1e-14)


def test_virtual_temperature_round_trip_and_sign():
    rng = np.random.default_rng(1)
    T = rng.uniform(180.0, 320.0, 100)
    qv = rng.uniform(0.0, 0.02, 100)
    Tv = np.asarray(thermo.calculate_virtual_temperature(T, qv))
    assert np.all(Tv >= T)  # moist air is less dense -> Tv >= T
    np.testing.assert_allclose(
        np.asarray(thermo.calculate_temperature_from_virtual_temperature(Tv, qv)),
        T, rtol=1e-14)


def test_dse_round_trip():
    rng = np.random.default_rng(2)
    T = rng.uniform(180.0, 320.0, 72)
    z = rng.uniform(0.0, 40000.0, 72)
    phis = 123.4
    dse = np.asarray(thermo.calculate_dse(T, z, phis))
    np.testing.assert_allclose(
        np.asarray(thermo.calculate_temperature_from_dse(dse, z, phis)), T, rtol=1e-12)


def test_mmr_vmr_round_trip():
    rng = np.random.default_rng(3)
    qv = rng.uniform(0.0, 0.02, 50)
    mmr = rng.uniform(1e-9, 1e-3, 50)
    w = 44.0095  # co2
    vmr = np.asarray(thermo.calculate_vmr_from_mmr(w, qv, mmr))
    np.testing.assert_allclose(
        np.asarray(thermo.calculate_mmr_from_vmr(w, qv, vmr)), mmr, rtol=1e-14)


def test_wet_dry_mmr_round_trip():
    rng = np.random.default_rng(4)
    qv_wet = rng.uniform(0.0, 0.02, 50)
    wet = rng.uniform(1e-9, 1e-2, 50)
    dry = np.asarray(thermo.calculate_drymmr_from_wetmmr(wet, qv_wet))
    qv_dry = np.asarray(thermo.calculate_drymmr_from_wetmmr(qv_wet, qv_wet))
    np.testing.assert_allclose(
        np.asarray(thermo.calculate_wetmmr_from_drymmr(dry, qv_dry)), wet, rtol=1e-13)


def test_z_int_z_mid():
    # dz all ones, z_surf=0: interfaces count down from nlev (k=0 is model top).
    dz = np.ones((3, 5))
    z_int = np.asarray(thermo.calculate_z_int(dz, 0.0))
    np.testing.assert_allclose(z_int, np.tile(np.arange(5, -1, -1.0), (3, 1)))
    z_mid = np.asarray(thermo.calculate_z_mid(z_int))
    np.testing.assert_allclose(z_mid, np.tile(np.arange(4.5, 0.0, -1.0), (3, 1)))


def test_column_scan_directions():
    dx = np.array([1.0, 2.0, 3.0])
    np.testing.assert_allclose(
        np.asarray(column_ops.column_scan(dx, 10.0, from_top=True)),
        [10.0, 11.0, 13.0, 16.0])
    np.testing.assert_allclose(
        np.asarray(column_ops.column_scan(dx, 10.0, from_top=False)),
        [16.0, 15.0, 13.0, 10.0])


def test_midpoint_delta_inverts_scan():
    rng = np.random.default_rng(5)
    dx = rng.standard_normal((2, 7))
    # Forward delta inverts the from-top scan; the from-bottom scan
    # accumulates downward-decreasing values, so the delta is -dx.
    x_i = column_ops.column_scan(dx, 3.0, from_top=True)
    np.testing.assert_allclose(
        np.asarray(column_ops.compute_midpoint_delta(x_i)), dx, atol=1e-12)
    x_i = column_ops.column_scan(dx, 3.0, from_top=False)
    np.testing.assert_allclose(
        np.asarray(column_ops.compute_midpoint_delta(x_i)), -dx, atol=1e-12)


def test_interface_values_linear_constant_field():
    # A constant midpoint field with matching bcs must stay constant.
    x_m = np.full((4,), 7.5)
    dz = np.array([1.0, 2.0, 0.5, 3.0])
    x_i = np.asarray(column_ops.compute_interface_values_linear(x_m, dz, 7.5, 7.5))
    np.testing.assert_allclose(x_i, 7.5)


def test_interface_values_compatible_round_trip():
    # Documented property: midpoints of the result reproduce x_m exactly.
    rng = np.random.default_rng(6)
    x_m = rng.standard_normal(6)
    for fix_top in (True, False):
        x_i = column_ops.compute_interface_values_compatible(x_m, 0.25, fix_top=fix_top)
        np.testing.assert_allclose(
            np.asarray(column_ops.compute_midpoint_values(x_i)), x_m, atol=1e-12)
        bc_idx = 0 if fix_top else -1
        assert np.asarray(x_i)[bc_idx] == pytest.approx(0.25)


def test_psl_branches():
    # Near sea level: psl == p_ground exactly.
    assert float(thermo.calculate_psl(288.0, 101325.0, 0.0)) == 101325.0
    # Elevated, ordinary temperature: psl > p_ground.
    psl = float(thermo.calculate_psl(280.0, 70000.0, 3000.0 * 9.80616))
    assert psl > 70000.0
    # Sanity: roughly 1000 hPa slot for a ~3 km mountain at 70 kPa.
    assert 90000.0 < psl < 110000.0


def test_rayleigh_friction_qualitative():
    # The scheme (transcribed exactly from the C++) damps winds isotropically
    # and heats. NOTE: it is not exactly KE-conserving in this direct-update
    # form; the authoritative check is Tier-0 golden data vs EAMxx.
    dt, otau = 100.0, np.array([1e-3])
    u, v, T = np.array([30.0]), np.array([-20.0]), np.array([250.0])
    u2, v2, T2 = (np.asarray(a).item() for a in
                  thermo.apply_rayleigh_friction(dt, otau, u, v, T))
    assert 0.0 < u2 < 30.0
    assert -20.0 < v2 < 0.0
    np.testing.assert_allclose(u2 / v2, 30.0 / -20.0, rtol=1e-13)
    assert T2 > 250.0
    # otau = 0 must be an exact no-op for winds and temperature.
    u3, v3, T3 = (np.asarray(a).item() for a in
                  thermo.apply_rayleigh_friction(dt, np.array([0.0]), u, v, T))
    assert u3 == 30.0 and v3 == -20.0 and T3 == 250.0


# --- tridiag ------------------------------------------------------------------

from scream_jax.foundation import tridiag  # noqa: E402


def _dense_from_diags(dl, d, du):
    n = d.shape[-1]
    T = np.zeros((n, n))
    T[np.arange(n), np.arange(n)] = d
    T[np.arange(1, n), np.arange(n - 1)] = dl[1:]
    T[np.arange(n - 1), np.arange(1, n)] = du[:-1]
    return T


def test_thomas_matches_dense_solve():
    rng = np.random.default_rng(7)
    n = 72
    d = rng.uniform(4.0, 6.0, n)          # diagonally dominant
    dl = rng.uniform(-1.0, 1.0, n)
    du = rng.uniform(-1.0, 1.0, n)
    x = rng.standard_normal(n)
    sol = np.asarray(tridiag.thomas(dl, d, du, x))
    expected = np.linalg.solve(_dense_from_diags(dl, d, du), x)
    np.testing.assert_allclose(sol, expected, rtol=1e-12)


def test_thomas_batched():
    rng = np.random.default_rng(8)
    ncol, n = 5, 33
    d = rng.uniform(4.0, 6.0, (ncol, n))
    dl = rng.uniform(-1.0, 1.0, (ncol, n))
    du = rng.uniform(-1.0, 1.0, (ncol, n))
    x = rng.standard_normal((ncol, n))
    sol = np.asarray(tridiag.thomas(dl, d, du, x))
    for i in range(ncol):
        expected = np.linalg.solve(_dense_from_diags(dl[i], d[i], du[i]), x[i])
        np.testing.assert_allclose(sol[i], expected, rtol=1e-11)


def test_thomas_unused_corner_entries_ignored():
    rng = np.random.default_rng(9)
    n = 8
    d = rng.uniform(4.0, 6.0, n)
    dl = rng.uniform(-1.0, 1.0, n)
    du = rng.uniform(-1.0, 1.0, n)
    x = rng.standard_normal(n)
    base = np.asarray(tridiag.thomas(dl, d, du, x))
    dl2, du2 = dl.copy(), du.copy()
    dl2[0], du2[-1] = 1e9, -1e9   # must have no effect
    np.testing.assert_array_equal(np.asarray(tridiag.thomas(dl2, d, du2, x)), base)
