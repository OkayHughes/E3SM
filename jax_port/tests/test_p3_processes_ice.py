"""Tier-0 sanity tests for P3 ice process rates (batch 2)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.p3 import processes_ice as pi  # noqa: E402

OPTS = dict(cldliq_to_ice_collection_factor=1.0,
            rain_to_ice_collection_factor=1.0)


def test_cldliq_collection_temperature_split():
    ctx = np.ones(3, dtype=bool)
    T = np.array([260.0, 275.0, 260.0])
    qi = np.array([1e-4, 1e-4, 0.0])   # third: no ice
    qccol, nccol, qcshd, ncshdc = (np.asarray(a) for a in pi.ice_cldliq_collection(
        np.ones(3), T, np.ones(3), np.full(3, 1e-3), qi, np.full(3, 1e-4),
        np.full(3, 1e4), np.full(3, 1e8), OPTS, ctx))
    assert qccol[0] > 0 and qcshd[0] == 0          # freezing: riming
    assert qccol[1] == 0 and qcshd[1] > 0          # warm: shedding
    assert ncshdc[1] > 0 and ncshdc[0] == 0
    assert qccol[2] == 0 and qcshd[2] == 0         # no ice: nothing
    assert nccol[0] > 0 and nccol[1] > 0


def test_ice_self_collection_rime_efficiency():
    ctx = np.ones(3, dtype=bool)
    qi = np.full(3, 1e-4)
    qm = np.array([0.0, 0.5e-4, 0.95e-4])          # rime fractions 0, .5, .95
    nis = (np.asarray(pi.ice_self_collection(
        np.ones(3), np.ones(3), np.full(3, 1e-3), np.full(3, 0.1),
        qm, qi, np.full(3, 1e4), ctx)))
    assert nis[0] > 0 and nis[1] > 0
    assert nis[2] == 0                              # rime frac >= 0.9: sticking 0


def test_ice_melting_only_above_freezing():
    ctx = np.ones(2, dtype=bool)
    T = np.array([275.0, 270.0])
    qimlt, nimlt = (np.asarray(a) for a in pi.ice_melting(
        np.ones(2), T, np.full(2, 9e4), np.ones(2), np.full(2, 1e-2),
        np.full(2, 1e-2), np.full(2, 2e-5), np.full(2, 0.6), np.full(2, 1.7e-5),
        np.full(2, 0.024), np.full(2, 5e-3), np.full(2, 1e-4), np.full(2, 1e4), ctx))
    assert qimlt[0] > 0 and nimlt[0] > 0
    assert qimlt[1] == 0 and nimlt[1] == 0


def test_deposition_sublimation_branches():
    ctx = np.ones(3, dtype=bool)
    qv = np.array([6e-3, 3e-3, 6e-3])
    qv_sat_i = np.full(3, 4e-3)
    qv_sat_l = np.full(3, 5e-3)
    T = np.array([260.0, 260.0, 275.0])
    dep, sub, nsub, berg = (np.asarray(a) for a in pi.ice_deposition_sublimation(
        np.full(3, 1e-4), np.full(3, 1e4), T, qv_sat_l, qv_sat_i,
        np.full(3, 1e-2), np.ones(3), qv, 1.0 / 300, ctx))
    assert dep[0] > 0 and sub[0] == 0 and berg[0] > 0   # supersat, cold
    assert dep[1] == 0 and sub[1] > 0 and nsub[1] > 0   # subsat: sublimation
    assert dep[2] == 0 and berg[2] == 0                 # warm: no deposition


def test_wet_growth_sheds_and_reduces_collection():
    ctx = np.ones(1, dtype=bool)
    qccol = np.array([1e-5])
    qrcol = np.array([2e-5])
    # tiny growth rate forces shedding of nearly all collected water
    log_wg, qrcol2, qccol2, growth, nrshdr, qcshd = (np.asarray(a) for a in
        pi.ice_cldliq_wet_growth(
            np.ones(1), np.array([271.0]), np.array([9e4]), np.ones(1),
            np.array([1e-10]), np.array([1e-10]), np.array([2e-5]),
            np.array([0.024]), np.array([1.7e-5]), np.array([0.6]),
            np.array([9e-3]),  # very moist -> small growth
            np.array([1e-4]), np.array([1e-4]), np.array([1e4]), np.array([1e-4]),
            qrcol, qccol, np.zeros(1), np.zeros(1), ctx))
    assert bool(log_wg[0])
    assert qccol2[0] < qccol[0] and qrcol2[0] < qrcol[0]
    assert nrshdr[0] > 0 and qcshd[0] > 0
    # shed + remaining collection ~ original collection + growth balance
    assert qccol2[0] >= 0 and qrcol2[0] >= 0


def test_evaporate_rain_subsaturated_only():
    ctx = np.ones(3, dtype=bool)
    qv = np.array([4e-3, 6e-3, 4e-3])
    qv_sat_l = np.full(3, 5e-3)
    cld_frac_r = np.array([0.8, 0.8, 0.2])
    cld_frac_l = np.array([0.1, 0.1, 0.5])
    # column 3 needs real cloud water so cld_frac_l is honored (the kernel
    # zeroes the effective cloud fraction when qc+qi < 1e-6, as in the C++)
    qc_incld = np.array([0.0, 0.0, 1e-5])
    tend, ntend = (np.asarray(a) for a in pi.evaporate_rain(
        np.full(3, 1e-4), qc_incld, np.full(3, 1e3), np.zeros(3),
        cld_frac_l, cld_frac_r, qv, qv, qv_sat_l, np.full(3, 4.5e-3),
        np.ones(3), np.ones(3), np.full(3, 1e-2), np.zeros(3),
        np.full(3, 280.0), np.full(3, 280.0), np.full(3, 3e-4), 300.0, ctx))
    assert tend[0] > 0 and ntend[0] > 0     # subsaturated, clear-sky rain
    assert tend[1] == 0                     # supersaturated: no evap
    assert tend[2] == 0                     # rain area <= cloudy area
    # bounded by removing all rain in one step
    assert tend[0] <= 1e-4 / 300.0 + 1e-20
