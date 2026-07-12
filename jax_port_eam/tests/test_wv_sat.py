"""Tier-1 golden replay + Tier-0 property tests for eam_jax.wv_sat."""


import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "wv_sat_golden.npz"

from eam_jax import wv_sat  # noqa: E402

RTOL = 1e-12


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


# ---------------------------------------------------------------------------
# Tier-1: golden replay
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name,idx", [("oldgoffgratch", 0),
                                      ("goffgratch", 1),
                                      ("murphykoop", 2), ("bolton", 3)])
def test_svp_schemes(gold, name, idx):
    t = gold["t_sweep"]
    np.testing.assert_allclose(np.asarray(wv_sat.svp_water(t, idx)),
                               gold[f"svp_water_{name}"], rtol=RTOL)
    np.testing.assert_allclose(np.asarray(wv_sat.svp_ice(t, idx)),
                               gold[f"svp_ice_{name}"], rtol=RTOL)
    np.testing.assert_allclose(np.asarray(wv_sat.svp_trans(t, idx)),
                               gold[f"svp_trans_{name}"], rtol=RTOL)


def test_estblf(gold):
    np.testing.assert_allclose(np.asarray(wv_sat.estblf(gold["t_sweep"])),
                               gold["estblf"], rtol=RTOL)


@pytest.mark.parametrize("fam,fn", [
    ("qsat", lambda t, p: wv_sat.qsat(t, p)),
    ("qsat_water", lambda t, p: wv_sat.qsat_water(t, p)),
    ("qsat_ice", lambda t, p: wv_sat.qsat_ice(t, p))])
def test_qsat_family(gold, fam, fn):
    t, p = gold["qsat_t"], gold["qsat_p"]
    out = fn(t, p)
    for key in ("es", "qs", "gam", "dqsdt", "enthalpy"):
        np.testing.assert_allclose(np.asarray(out[key]),
                                   gold[f"{fam}_{key}"], rtol=RTOL,
                                   atol=1e-300, err_msg=f"{fam}.{key}")


@pytest.mark.parametrize("tag,use_ice", [("ice", True), ("noice", False)])
def test_findsp(gold, tag, use_ice):
    q, t, p = gold["findsp_q"], gold["findsp_t"], gold["findsp_p"]
    tsp, qsp, status = wv_sat.findsp(q, t, p, use_ice)
    np.testing.assert_array_equal(np.asarray(status),
                                  gold[f"findsp_status_{tag}"])
    np.testing.assert_allclose(np.asarray(tsp), gold[f"findsp_tsp_{tag}"],
                               rtol=1e-10)
    np.testing.assert_allclose(np.asarray(qsp), gold[f"findsp_qsp_{tag}"],
                               rtol=1e-9, atol=1e-14)


# ---------------------------------------------------------------------------
# Tier-0: properties
# ---------------------------------------------------------------------------
def test_svp_monotone_and_positive():
    t = np.linspace(150.0, 370.0, 500)
    for idx in range(4):
        es = np.asarray(wv_sat.svp_water(t, idx))
        assert (es > 0).all()
        assert (np.diff(es) > 0).all(), f"scheme {idx} not monotone"


def test_ice_below_water_near_freezing():
    # es_ice < es_water below 0 C (thermodynamic requirement)
    t = np.linspace(230.0, 273.0, 100)
    esw = np.asarray(wv_sat.svp_water(t, 1))
    esi = np.asarray(wv_sat.svp_ice(t, 1))
    assert (esi < esw).all()


def test_qs_limiter():
    # very low pressure: qs must clamp to 1, es to p
    out = wv_sat.qsat(np.array([350.0]), np.array([50.0]))
    assert float(out["qs"][0]) == 1.0
    assert float(out["es"][0]) == 50.0
    assert float(out["dqsdt"][0]) == 0.0


def test_findsp_conserves_enthalpy():
    rng = np.random.default_rng(7)
    t = rng.uniform(250.0, 310.0, 200)
    p = rng.uniform(5.0e4, 1.0e5, 200)
    q = rng.uniform(1e-4, 2e-2, 200)
    tsp, qsp, status = wv_sat.findsp(q, t, p, True)
    ok = np.asarray(status) == 0
    assert ok.mean() > 0.9
    hlt0 = np.asarray(wv_sat.calc_hltalt(np.asarray(t))[0])
    hlt1 = np.asarray(wv_sat.calc_hltalt(np.asarray(tsp))[0])
    en_in = wv_sat.CPAIR * t + hlt0 * q
    en_out = wv_sat.CPAIR * np.asarray(tsp) + hlt1 * np.asarray(qsp)
    rel = np.abs(en_in - en_out) / en_in
    assert rel[ok].max() < 2e-4
