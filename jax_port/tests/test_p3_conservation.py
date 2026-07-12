"""Tier-0 tests for P3 conservation limiters (batch 3)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.p3 import conservation as pc  # noqa: E402


def test_cloud_water_conservation_scales_sinks():
    ctx = np.ones(2, dtype=bool)
    dt = 300.0
    qc = np.array([1e-5, 1e-2])   # first: sinks deplete qc; second: plenty
    t = np.full(2, 1e-7)          # each sink 1e-7*300 = 3e-5 > 1e-5
    out = pc.cloud_water_conservation(
        qc, dt, t, t, t, t, t, t, np.full(2, 1e-8), np.full(2, 1e-8),
        np.zeros(2), np.zeros(2), False, ctx, np.ones(2), np.ones(2))
    qcaut, qcacc, qccol, qchetc, qcshd, qcberg, qisub, qidep = \
        (np.asarray(a) for a in out[:8])
    sinks_after = (qcaut + qcacc + qccol + qchetc + qcshd + qcberg)[0] * dt
    np.testing.assert_allclose(sinks_after, qc[0], rtol=1e-12)  # exactly limited
    np.testing.assert_allclose(qcaut[1], 1e-7)                  # untouched
    # ratio=1 lanes zero vapdep/sublim wherever qc > qtendsmall (C++ quirk)
    assert qidep[1] == 0.0 and qisub[1] == 0.0
    assert qidep[0] > 0.0                                       # (1-ratio) > 0


def test_rain_and_ice_water_conservation():
    ctx = np.ones(1, dtype=bool)
    dt = 300.0
    qrevp, qrcol, qrheti = pc.rain_water_conservation(
        np.array([1e-6]), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1),
        dt, np.array([1e-7]), np.array([1e-7]), np.array([1e-7]), ctx)
    total = (np.asarray(qrevp) + np.asarray(qrcol) + np.asarray(qrheti)) * dt
    np.testing.assert_allclose(total, 1e-6, rtol=1e-12)

    qisub, qimlt = pc.ice_water_conservation(
        np.array([1e-6]), *(np.zeros(1),) * 7, dt,
        np.zeros(1), np.zeros(1), np.zeros(1),
        np.array([1e-7]), np.array([1e-7]), False, ctx)
    np.testing.assert_allclose((np.asarray(qisub) + np.asarray(qimlt)) * dt,
                               1e-6, rtol=1e-12)


def test_number_conservations():
    ctx = np.ones(1, dtype=bool)
    dt = 300.0
    out = pc.nc_conservation(np.array([1e3]), np.zeros(1), dt,
                             np.array([10.0]), np.array([10.0]),
                             np.array([10.0]), np.array([10.0]),
                             np.zeros(1), np.zeros(1), False, ctx)
    total = sum(np.asarray(a).item() for a in (out[0], out[1], out[2], out[3])) * dt
    np.testing.assert_allclose(total, 1e3, rtol=1e-12)

    nrcol, nrheti, nrslf, nrevp = pc.nr_conservation(
        np.array([1e3]), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1),
        dt, 1.0, np.array([10.0]), np.array([10.0]), np.array([10.0]),
        np.array([10.0]), ctx)
    total = sum(np.asarray(a).item() for a in (nrcol, nrheti, nrslf, nrevp)) * dt
    np.testing.assert_allclose(total, 1e3, rtol=1e-12)

    nimlt, nisub, nislf = pc.ni_conservation(
        np.array([1e3]), np.zeros(1), np.zeros(1), np.zeros(1),
        np.zeros(1), np.zeros(1), np.zeros(1), dt,
        np.array([10.0]), np.array([10.0]), np.array([10.0]), False, ctx)
    total = sum(np.asarray(a).item() for a in (nimlt, nisub, nislf)) * dt
    np.testing.assert_allclose(total, 1e3, rtol=1e-12)


def test_ice_supersat_conservation_limits_to_available_vapor():
    ctx = np.ones(1, dtype=bool)
    qidep, qinuc, _ = pc.ice_supersat_conservation(
        np.array([1e-5]), np.array([1e-5]), np.zeros(1), np.array([0.5]),
        np.array([5e-3]), np.array([4.9e-3]), np.array([250.0]), 300.0,
        np.zeros(1), np.zeros(1), False, ctx)
    # sink was 2e-5/s * 300s = 6e-3 >> available ~1e-4: must be scaled down
    assert np.asarray(qidep).item() < 1e-5
    assert np.asarray(qinuc).item() < 1e-5


def test_impose_max_total_ni():
    ctx = np.ones(2, dtype=bool)
    ni = np.array([1e6, 1e2])
    out = np.asarray(pc.impose_max_total_ni(ni, 1e5, np.ones(2), ctx))
    np.testing.assert_allclose(out, [1e5, 1e2])


def test_incloud_mixingratios_limits():
    ctx = np.ones(2, dtype=bool)
    qc = np.array([1e-3, 1e-15])
    out = pc.calculate_incloud_mixingratios(
        qc, np.full(2, 1e-3), np.full(2, 1e-3), np.full(2, 5e-4),
        np.full(2, 1e8), np.full(2, 1e5), np.full(2, 1e4), np.full(2, 1e-6),
        np.full(2, 10.0), np.full(2, 10.0), np.full(2, 10.0), ctx)
    qc_in, qr_in, qi_in = (np.asarray(a) for a in out[:3])
    assert qc_in[0] <= 5.1e-3 + 1e-15    # incloud_limit
    assert qr_in[0] <= 1.0e-2 + 1e-15    # precip_limit
    assert qc_in[1] == 0.0               # below qsmall
