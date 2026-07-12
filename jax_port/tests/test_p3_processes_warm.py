"""Tier-0 sanity tests for P3 warm/freezing process rates (batch 1)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.p3 import processes_warm as pw  # noqa: E402

OPTS = dict(autoconversion_prefactor=1350.0, autoconversion_qc_exponent=2.47,
            autoconversion_nc_exponent=1.79, autoconversion_radius=25e-6,
            accretion_prefactor=67.0, accretion_qc_exponent=1.15,
            accretion_qr_exponent=1.15,
            rain_selfcollection_breakup_diameter=0.00028,
            rain_selfcollection_prefactor=5.78,
            deposition_nucleation_exponent=0.304,
            immersion_freezing_exponent=0.65)


def test_autoconversion_masks_and_signs():
    ctx = np.ones(4, dtype=bool)
    rho = np.full(4, 1.0)
    qc = np.array([0.0, 1e-9, 1e-4, 2e-3])
    nc = np.full(4, 1e8)
    qcaut, ncautc, ncautr = (np.asarray(a) for a in
        pw.cloud_water_autoconversion(rho, qc, nc, np.ones(4), OPTS, ctx))
    assert qcaut[0] == 0 and qcaut[1] == 0        # below 1e-8 threshold
    assert qcaut[2] > 0 and qcaut[3] > qcaut[2]   # increasing in qc
    assert np.all(ncautc >= 0) and np.all(ncautr >= 0)


def test_accretion_requires_both_species():
    ctx = np.ones(3, dtype=bool)
    args = (np.ones(3), np.ones(3),
            np.array([1e-4, 0.0, 1e-4]), np.full(3, 1e8),
            np.array([1e-4, 1e-4, 0.0]), np.ones(3))
    qccol, ncacc = (np.asarray(a) for a in
                    pw.cloud_rain_accretion(*args, OPTS, ctx))
    assert qccol[0] > 0 and qccol[1] == 0 and qccol[2] == 0


def test_rain_self_collection_breakup_reduces_rate():
    ctx = np.ones(2, dtype=bool)
    qr = np.full(2, 1e-3)
    nr = np.array([1e2, 1e6])
    tend = np.asarray(pw.rain_self_collection(np.ones(2), qr, nr, OPTS, ctx))
    d = np.cbrt(qr / (c.Pi * c.RHO_H2O * nr))
    assert d[0] > OPTS["rain_selfcollection_breakup_diameter"] > d[1]
    assert tend[0] / nr[0] < tend[1] / nr[1]


def test_ice_nucleation_branches():
    # do_log = (!predict_nc || prescribed): True -> Cooper-type formula,
    # False (predict_nc alone) -> ni_activated relaxation. (The C++ mask
    # names `any_if_log`/`any_if_not_log` are inverted relative to do_log.)
    ctx = np.ones(3, dtype=bool)
    cold = np.array([240.0, 240.0, 280.0])
    supersat = np.array([0.1, 0.01, 0.1])
    # predict_nc=False -> Cooper formula: needs cold + supersaturated,
    # independent of ni_activated
    q, n = (np.asarray(a) for a in pw.ice_nucleation(
        cold, np.ones(3), np.zeros(3), np.zeros(3), supersat,
        1.0 / 300.0, False, False, OPTS, ctx))
    assert n[0] > 0 and q[0] > 0
    assert n[1] == 0 and n[2] == 0
    # predict_nc=True -> relax ni toward ni_activated
    q2, n2 = (np.asarray(a) for a in pw.ice_nucleation(
        cold, np.ones(3), np.zeros(3), np.full(3, 5e4), supersat,
        1.0 / 300.0, True, False, OPTS, ctx))
    np.testing.assert_allclose(n2[0], 5e4 / 300.0, rtol=1e-14)
    assert q2[0] > 0 and n2[1] == 0 and n2[2] == 0
    # with nothing to activate, the predict_nc branch nucleates nothing
    q3, n3 = (np.asarray(a) for a in pw.ice_nucleation(
        cold, np.ones(3), np.zeros(3), np.zeros(3), supersat,
        1.0 / 300.0, True, False, OPTS, ctx))
    assert np.all(n3 == 0) and np.all(q3 == 0)


def test_immersion_freezing_temperature_gate():
    ctx = np.ones(3, dtype=bool)
    T = np.array([260.0, 268.0, 271.0])  # T_rainfrz = 269.15
    lamc = np.full(3, 5e4)
    qchetc, nchetc = (np.asarray(a) for a in pw.cldliq_immersion_freezing(
        T, lamc, np.full(3, 5.0), np.full(3, 1e8), np.full(3, 1e-4),
        np.ones(3), OPTS, ctx))
    assert qchetc[0] > qchetc[1] > 0
    assert qchetc[2] == 0
    qrcol, nrcol = (np.asarray(a) for a in pw.rain_immersion_freezing(
        T, np.full(3, 5e3), np.ones(3), np.full(3, 1e6), np.full(3, 1e-4),
        OPTS, ctx))
    assert qrcol[0] > qrcol[1] > 0 and qrcol[2] == 0


def test_calc_rime_density_bounds():
    ctx = np.ones(4, dtype=bool)
    T = np.array([250.0, 260.0, 272.0, 274.0])
    vtrmi1, rho_qm = (np.asarray(a) for a in pw.calc_rime_density(
        T, np.ones(4), np.full(4, 1.0), np.full(4, 1e5), np.full(4, 5e4),
        np.full(4, 5.0), np.full(4, 1e-4), np.full(4, 1e-6), ctx))
    assert np.all((rho_qm >= c.RHO_RIMEMIN - 1) | (rho_qm == 400.0))
    assert rho_qm[3] == 400.0 and vtrmi1[3] == 0.0
