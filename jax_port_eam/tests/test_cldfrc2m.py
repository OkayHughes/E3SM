"""Tier-1 golden replay + Tier-0 properties for the cldfrc2m
cloud-fraction port (eam_jax.cldfrc2m).

Golden archive: harness/gen_cldfrc2m_golden.py (run in the scream-dev
container) — a 40x72 flattened profile set through astG_PDF, astG_RHU
and aist for every iceopt 1-7 (plus an alternate rhmini/rhmaxi pair
for the default iceopt=5), and dense 1-D single-point sweeps of the
_single routines (U x pressure x surface type, recording a/G/orhmin;
qv/qi/T sweeps for aist_single at iceopt=5).

The port is like-for-like float64; every comparison holds at
rtol=1e-12 (libm last-ulp differences only). Sentinel G values
(1e10/1e20) and hard zeros replay exactly.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "cldfrc2m_golden.npz"

from eam_jax import cldfrc2m  # noqa: E402

ICEOPTS = [1, 2, 3, 4, 5, 6, 7]


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


@pytest.fixture(scope="module")
def prm(gold):
    p = json.loads(str(gold["__metadata__"]))["params"]
    return cldfrc2m.make_params(
        rhminl=p["rhminl"], rhminl_adj_land=p["rhminl_adj_land"],
        rhminh=p["rhminh"], premit=p["premit"], premib=p["premib"],
        iceopt=p["iceopt_default"], icecrit=p["icecrit"],
        minice=p["minice"], rhmini=p["rhmini"], rhmaxi=p["rhmaxi"])


# ---------------------------------------------------------------------------
# Tier-1 golden replay — vector profile set
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name,fn", [("pdf", cldfrc2m.astG_PDF),
                                     ("rhu", cldfrc2m.astG_RHU)])
def test_astg_profiles(gold, prm, name, fn):
    a, ga, _ = fn(gold["u"], gold["p"], gold["qv"], gold["landfrac"],
                  gold["snowh"], prm)
    np.testing.assert_allclose(np.asarray(a), gold[f"a_{name}"],
                               rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(np.asarray(ga), gold[f"ga_{name}"],
                               rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("opt", ICEOPTS)
def test_aist_profiles(gold, prm, opt):
    p = dict(prm, iceopt=opt)
    ist = cldfrc2m.aist(gold["qv"], gold["t"], gold["p"], gold["qi"],
                        gold["ni"], gold["landfrac"], gold["snowh"], p)
    np.testing.assert_allclose(np.asarray(ist), gold[f"aist_opt{opt}"],
                               rtol=1e-12, atol=0.0)


def test_aist_alt_rh_bounds(gold, prm):
    meta = json.loads(str(gold["__metadata__"]))["params"]
    p = dict(prm, rhmini=meta["rhmini_alt"], rhmaxi=meta["rhmaxi_alt"])
    ist = cldfrc2m.aist(gold["qv"], gold["t"], gold["p"], gold["qi"],
                        gold["ni"], gold["landfrac"], gold["snowh"], p)
    np.testing.assert_allclose(np.asarray(ist), gold["aist_opt5_alt"],
                               rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# Tier-1 golden replay — dense single-point sweeps
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name,fn", [("pdf", cldfrc2m.astG_PDF),
                                     ("rhu", cldfrc2m.astG_RHU)])
def test_astg_single_sweeps(gold, prm, name, fn):
    u = gold["sweep_u"]
    for ip, p in enumerate(gold["sweep_p"]):
        for isf, (lf, sh) in enumerate(zip(gold["sweep_landfrac"],
                                           gold["sweep_snowh"])):
            a, ga, rhmin = fn(u, p, 5.0e-3, lf, sh, prm)
            np.testing.assert_allclose(
                np.asarray(a), gold[f"sweep_a_{name}"][ip, isf],
                rtol=1e-12, atol=0.0)
            np.testing.assert_allclose(
                np.asarray(ga), gold[f"sweep_ga_{name}"][ip, isf],
                rtol=1e-12, atol=0.0)
            np.testing.assert_allclose(
                np.broadcast_to(np.asarray(rhmin), u.shape),
                gold[f"sweep_rhmin_{name}"][ip, isf],
                rtol=1e-12, atol=0.0)


def test_aist_single_sweeps(gold, prm):
    qv = gold["aist_qv_sweep_in"]
    ist = cldfrc2m.aist(qv, 230.0, 30000.0, 1.0e-5, 0.0, 0.0, 0.0, prm)
    np.testing.assert_allclose(np.asarray(ist), gold["aist_qv_sweep"],
                               rtol=1e-12, atol=0.0)
    qi = gold["aist_qi_sweep_in"]
    ist = cldfrc2m.aist(2.0e-4, 230.0, 30000.0, qi, 0.0, 0.0, 0.0, prm)
    np.testing.assert_allclose(np.asarray(ist), gold["aist_qi_sweep"],
                               rtol=1e-12, atol=0.0)
    t = gold["aist_t_sweep_in"]
    ist = cldfrc2m.aist(1.0e-4, t, 30000.0, 1.0e-5, 0.0, 0.0, 0.0, prm)
    np.testing.assert_allclose(np.asarray(ist), gold["aist_t_sweep"],
                               rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# Tier-0 properties
# ---------------------------------------------------------------------------
def _rand_state(n=4096, seed=7):
    rng = np.random.default_rng(seed)
    return dict(u=rng.uniform(-0.1, 1.3, n),
                p=rng.uniform(2000.0, 1.05e5, n),
                qv=rng.uniform(1e-9, 2e-2, n),
                t=rng.uniform(180.0, 310.0, n),
                qi=10.0 ** rng.uniform(-14.0, -2.0, n),
                ni=10.0 ** rng.uniform(0.0, 7.0, n),
                lf=rng.choice([0.0, 0.5, 1.0], n),
                sh=rng.choice([0.0, 1e-6, 0.1], n))


def test_bounds():
    """a in [0,1] for both liquid schemes, aist in [0, 0.999] for
    every iceopt, on broad random inputs."""
    s = _rand_state()
    prm = cldfrc2m.make_params()
    for fn in (cldfrc2m.astG_PDF, cldfrc2m.astG_RHU):
        a, ga, rhmin = fn(s["u"], s["p"], s["qv"], s["lf"], s["sh"],
                          prm)
        a = np.asarray(a)
        assert a.min() >= 0.0 and a.max() <= 1.0
        assert np.asarray(ga).min() > 0.0
        assert np.all((np.asarray(rhmin) >= 0.80)
                      & (np.asarray(rhmin) <= 0.95))
    for opt in ICEOPTS:
        ist = np.asarray(cldfrc2m.aist(
            s["qv"], s["t"], s["p"], s["qi"], s["ni"], s["lf"],
            s["sh"], dict(prm, iceopt=opt)))
        assert ist.min() >= 0.0 and ist.max() <= 0.999


def test_monotone_in_rh():
    """Stratus fraction is nondecreasing in RH in every pressure band
    and for every surface type, and aist is nondecreasing in qv."""
    prm = cldfrc2m.make_params()
    u = np.linspace(0.0, 1.2, 2401)
    for p in (95000.0, 50000.0, 10000.0):
        for lf, sh in ((0.0, 0.0), (1.0, 0.0), (1.0, 0.1)):
            for fn in (cldfrc2m.astG_PDF, cldfrc2m.astG_RHU):
                a = np.asarray(fn(u, p, 5e-3, lf, sh, prm)[0])
                assert np.all(np.diff(a) >= 0.0)
    qv = np.linspace(0.0, 6e-4, 1201)
    ist = np.asarray(cldfrc2m.aist(qv, 230.0, 30000.0, 1e-5, 0.0,
                                   0.0, 0.0, prm))
    assert np.all(np.diff(ist) >= -1e-15)


def test_zero_below_critical_rh_and_one_at_saturation():
    """a == 0 exactly for U <= rhmin and a == 1 exactly for U >= 1;
    aist == 0 exactly below rhmini (with negligible qi effect) and
    for qi < minice regardless of RH."""
    prm = cldfrc2m.make_params()
    s = _rand_state(seed=11)
    for fn in (cldfrc2m.astG_PDF, cldfrc2m.astG_RHU):
        a, _, rhmin = fn(s["u"], s["p"], s["qv"], s["lf"], s["sh"],
                         prm)
        a, rhmin = np.asarray(a), np.asarray(rhmin)
        assert np.all(a[s["u"] <= rhmin] == 0.0)
        assert np.all(a[s["u"] >= 1.0] == 1.0)
    # iceopt 5: dry + qi below minice -> exactly zero
    ist = np.asarray(cldfrc2m.aist(
        s["qv"], s["t"], s["p"], np.full_like(s["qi"], 1e-13),
        s["ni"], s["lf"], s["sh"], prm))
    assert np.all(ist == 0.0)


def test_land_lowers_critical_rh_only_when_snowfree():
    """Snow-free land lowers rhmin by rhminl_adj_land in the lowest
    band only; snow-covered land and ocean share rhminl; high band is
    land-independent."""
    prm = cldfrc2m.make_params()
    _, _, rh_ocean = cldfrc2m.astG_PDF(0.9, 9e4, 5e-3, 0.0, 0.0, prm)
    _, _, rh_land = cldfrc2m.astG_PDF(0.9, 9e4, 5e-3, 1.0, 0.0, prm)
    _, _, rh_snow = cldfrc2m.astG_PDF(0.9, 9e4, 5e-3, 1.0, 0.1, prm)
    assert float(rh_land) == prm["rhminl"] - prm["rhminl_adj_land"]
    assert float(rh_ocean) == float(rh_snow) == prm["rhminl"]
    _, _, rh_hi_l = cldfrc2m.astG_PDF(0.9, 1e4, 5e-3, 1.0, 0.0, prm)
    assert float(rh_hi_l) == prm["rhminh"]
