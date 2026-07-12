"""Tier-0 tests for P3 cell-average scaling and prognostic updates."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.p3.cell_average import (  # noqa: E402
    _SCALE_MAP,
    back_to_cell_average,
    get_time_space_phys_variables,
)
from scream_jax.p3.update import (  # noqa: E402
    update_prognostic_ice,
    update_prognostic_liquid,
)


def test_back_to_cell_average_fractions():
    ctx = np.ones(1, dtype=bool)
    tends = {k: np.ones(1) for k in _SCALE_MAP}
    out = back_to_cell_average(np.array([0.5]), np.array([0.3]),
                               np.array([0.7]), tends, ctx)
    assert np.asarray(out["qc2qr_autoconv_tend"]).item() == 0.5   # l
    assert np.asarray(out["qr2qv_evap_tend"]).item() == 0.3       # r
    assert np.asarray(out["qi2qr_melt_tend"]).item() == 0.7       # i
    assert np.asarray(out["qc2qr_accret_tend"]).item() == 0.3     # min(l,r)
    assert np.asarray(out["qc2qi_collect_tend"]).item() == 0.5    # min(i,l)
    assert np.asarray(out["qr2qi_collect_tend"]).item() == 0.3    # min(i,r)


def test_phys_variables_sane():
    ctx = np.ones(3, dtype=bool)
    T = np.array([230.0, 263.15, 290.0])
    mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii = (np.asarray(a) for a in
        get_time_space_phys_variables(T, np.full(3, 8e4), np.ones(3),
                                      np.full(3, 5e-3), np.full(3, 4e-3), ctx))
    assert np.all(mu > 0) and np.all(dv > 0) and np.all(sc > 0)
    assert np.all(ab > 1) and np.all(abi > 1)
    np.testing.assert_allclose(eii[0], 0.001)
    np.testing.assert_allclose(eii[1], 0.001 + 10 * (0.3 - 0.001) / 20)
    np.testing.assert_allclose(eii[2], 0.3)


def _zero_tends():
    keys = ["qc2qi_hetero_freeze_tend", "qc2qi_collect_tend",
            "qc2qr_ice_shed_tend", "nc_collect_tend",
            "nc2ni_immers_freeze_tend", "ncshdc", "qr2qi_collect_tend",
            "nr_collect_tend", "qr2qi_immers_freeze_tend",
            "nr2ni_immers_freeze_tend", "nr_ice_shed_tend", "qi2qr_melt_tend",
            "ni2nr_melt_tend", "qi2qv_sublim_tend", "qv2qi_vapdep_tend",
            "qv2qi_nucleat_tend", "ni_nucleat_tend", "ni_selfcollect_tend",
            "ni_sublim_tend", "qc2qi_berg_tend", "ncheti_cnt", "nicnt",
            "ninuc_cnt", "qcheti_cnt", "qicnt", "qinuc_cnt", "inv_exner"]
    return {k: np.zeros(1) if k != "inv_exner" else np.ones(1) for k in keys}


def test_update_prognostic_ice_water_and_energy_budget():
    ctx = np.ones(1, dtype=bool)
    dt = 300.0
    t = _zero_tends()
    t["qv2qi_vapdep_tend"] = np.array([1e-7])   # vapor -> ice
    t["qi2qr_melt_tend"] = np.array([2e-8])     # ice -> rain
    state = dict(th_atm=np.array([300.0]), qv=np.array([5e-3]),
                 qi=np.array([1e-4]), ni=np.array([1e4]),
                 qm=np.array([5e-5]), bm=np.array([1e-7]),
                 qc=np.array([1e-4]), nc=np.array([1e8]),
                 qr=np.array([1e-5]), nr=np.array([1e3]))
    tot_before = (state["qv"] + state["qi"] + state["qc"] + state["qr"]).item()
    out = update_prognostic_ice(t, True, np.zeros(1, dtype=bool), dt, 1.0,
                                np.full(1, 400.0), state, False, ctx)
    tot_after = (np.asarray(out["qv"]) + np.asarray(out["qi"])
                 + np.asarray(out["qc"]) + np.asarray(out["qr"])).item()
    np.testing.assert_allclose(tot_after, tot_before, rtol=1e-12)  # water closed
    # deposition heats (th up), melting cools; net here: heating dominates
    assert np.asarray(out["th_atm"]).item() > 300.0
    assert np.asarray(out["qv"]).item() < 5e-3
    assert np.asarray(out["qi"]).item() > 1e-4


def test_update_prognostic_liquid_budget():
    ctx = np.ones(1, dtype=bool)
    dt = 300.0
    th, qv, qc, nc, qr, nr = update_prognostic_liquid(
        np.array([1e-7]), np.array([10.0]), np.array([5e-8]), np.array([5.0]),
        np.array([2.0]), np.array([1.0]), np.array([2e-8]), np.array([0.5]),
        np.array([0.3]), True, False, np.ones(1), np.ones(1), dt,
        np.array([300.0]), np.array([5e-3]), np.array([1e-4]),
        np.array([1e8]), np.array([1e-5]), np.array([1e3]), ctx)
    # water conservation: qc + qr + qv unchanged
    tot = np.asarray(qv).item() + np.asarray(qc).item() + np.asarray(qr).item()
    np.testing.assert_allclose(tot, 5e-3 + 1e-4 + 1e-5, rtol=1e-12)
    # evaporation cools
    assert np.asarray(th).item() < 300.0
    # prescribed-CCN branch resets nc
    _, _, _, nc2, _, _ = update_prognostic_liquid(
        *(np.zeros(1),) * 9, False, False, np.full(1, 2.0), np.ones(1), dt,
        np.array([300.0]), np.array([5e-3]), np.array([1e-4]),
        np.array([1e8]), np.array([1e-5]), np.array([1e3]), ctx)
    np.testing.assert_allclose(np.asarray(nc2).item(), c.NCCNST * 2.0)
