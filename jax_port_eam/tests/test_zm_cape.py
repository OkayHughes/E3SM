"""Tier-1 golden replay + Tier-0 properties for the ZM dilute CAPE
core (eam_jax.zm_cape).

Golden archive: harness/gen_zm_cape_golden.py (run in the scream-dev
container) — 42 columns x 72 levels, 6 regime families, 4 flag
configurations (a: EAMv3 defaults / ULL launch search, b: dcape
second call with frozen launch level, c: legacy pblt search with
num_cin=5, d: use_input_tq_mx). Level indices in the archive are
1-based (Fortran); the port is 0-based, hence the -1 shifts.

Tolerances: the Brent entropy inversion stops at a ~5e-4 K tolerance;
identical trajectories replay bit-for-bit, but a last-ulp libm
difference (container glibc vs host XLA log/exp) can flip a
convergence test one iteration early/late in rare lanes, shifting a
temperature by up to that tolerance. Observed agreement is recorded
per-field below; anything looser than 1e-10 relative is documented at
the assert.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "zm_cape_golden.npz"

from eam_jax import zm_cape  # noqa: E402

CFGS = ["a", "b", "c", "d"]
NPC = 7  # columns per regime family (see generator)


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


@pytest.fixture(scope="module")
def results(gold):
    """Run the port once per golden flag configuration."""
    meta = json.loads(str(gold["__metadata__"]))
    p = meta["params"]
    zc = zm_cape.make_zm_const()
    num_msg = int(gold["num_msg"])
    out = {}
    for cfg in CFGS:
        zp = zm_cape.make_zm_param(
            dmpdz=p["dmpdz"], tiedke_add=p["tiedke_add"],
            tpert_fac=p["tpert_fac"],
            mx_bot_lyr_adj=p["mx_bot_lyr_adj"],
            num_cin=p["num_cin"][cfg],
            tpert_fix=bool(p["tpert_fix"][cfg]),
            trig_ull=bool(p["trig_ull"][cfg]),
            trig_dcape=bool(p["trig_dcape"][cfg]))
        if p["state"][cfg] == "t/q":
            q, t = gold["q"], gold["t"]
        else:
            q, t = gold["q_star"], gold["t_star"]
        kw = {}
        if cfg == "b":
            kw = dict(calc_msemax_klev=False,
                      prev_msemax_klev=gold["msemax_klev_a"] - 1)
        elif cfg == "d":
            kw = dict(use_input_tq_mx=True,
                      prev_msemax_klev=gold["msemax_klev_a"] - 1,
                      q_mx=gold["q_mx_a"], t_mx=gold["t_mx_a"])
        out[cfg] = zm_cape.compute_dilute_cape(
            q, t, gold["zmid"], gold["pmid"], gold["pint"],
            gold["pblt"] - 1, gold["tpert"], num_msg, zc, zp, **kw)
    return out


# ---------------------------------------------------------------------------
# Tier-1 golden replay
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("cfg", CFGS)
def test_launch_and_level_indices(gold, results, cfg):
    r = results[cfg]
    np.testing.assert_array_equal(np.asarray(r["msemax_klev"]),
                                  gold[f"msemax_klev_{cfg}"] - 1)
    np.testing.assert_array_equal(np.asarray(r["lcl_klev"]),
                                  gold[f"lcl_klev_{cfg}"] - 1)
    np.testing.assert_array_equal(np.asarray(r["eql_klev"]),
                                  gold[f"eql_klev_{cfg}"] - 1)


@pytest.mark.parametrize("cfg", CFGS)
def test_parcel_profiles(gold, results, cfg):
    r = results[cfg]
    np.testing.assert_allclose(np.asarray(r["parcel_temp"]),
                               gold[f"parcel_temp_{cfg}"], rtol=1e-12)
    np.testing.assert_allclose(np.asarray(r["parcel_qsat"]),
                               gold[f"parcel_qsat_{cfg}"], rtol=1e-12,
                               atol=1e-18)
    np.testing.assert_allclose(np.asarray(r["lcl_temperature"]),
                               gold[f"lcl_temperature_{cfg}"],
                               rtol=1e-12)


@pytest.mark.parametrize("cfg", CFGS)
def test_cape(gold, results, cfg):
    r = results[cfg]
    # CAPE is an integral of near-cancelling buoyancies; 1e-10
    # relative with a 1e-8 J/kg floor for the zero-CAPE columns
    np.testing.assert_allclose(np.asarray(r["cape"]),
                               gold[f"cape_{cfg}"], rtol=1e-10,
                               atol=1e-8)


@pytest.mark.parametrize("cfg", CFGS)
def test_launch_tq_saved(gold, results, cfg):
    r = results[cfg]
    np.testing.assert_allclose(np.asarray(r["q_mx"]),
                               gold[f"q_mx_{cfg}"], rtol=1e-13)
    np.testing.assert_allclose(np.asarray(r["t_mx"]),
                               gold[f"t_mx_{cfg}"], rtol=1e-13)


# ---------------------------------------------------------------------------
# Tier-0 properties
# ---------------------------------------------------------------------------
def test_cape_nonnegative(results):
    for cfg in CFGS:
        assert float(np.min(np.asarray(results[cfg]["cape"]))) >= 0.0


def test_desert_high_lcl_gives_zero_cape(results):
    # dry desert family (columns 14-20): the LCL sits above the
    # 600 hPa threshold, which must suppress CAPE entirely
    r = results["a"]
    rows = slice(2 * NPC, 3 * NPC)
    assert (np.asarray(r["lcl_pmid"])[rows]
            < zm_cape.LCL_PRESSURE_THRESHOLD).all()
    np.testing.assert_array_equal(np.asarray(r["cape"])[rows], 0.0)


def test_parcel_equals_environment_below_launch(gold, results):
    # below the launch level the parcel arrays are the environment,
    # exactly (initialization copies, never overwritten)
    r = results["a"]
    pt = np.asarray(r["parcel_temp"])
    pq = np.asarray(r["parcel_qsat"])
    launch = np.asarray(r["msemax_klev"])
    ncol, pver = pt.shape
    for i in range(ncol):
        js = np.arange(launch[i] + 1, pver)
        np.testing.assert_array_equal(pt[i, js], gold["t"][i, js])
        np.testing.assert_array_equal(pq[i, js], gold["q"][i, js])


def test_parcel_temp_continuous_at_launch(gold, results):
    # at the launch level the parcel is the (entropy-inverted)
    # environment: continuity within the Brent tolerance (~5e-4 K)
    r = results["a"]
    pt = np.asarray(r["parcel_temp"])
    launch = np.asarray(r["msemax_klev"])
    rows = np.arange(pt.shape[0])
    dt = np.abs(pt[rows, launch] - gold["t"][rows, launch])
    assert dt.max() < 2.0e-3


def test_stable_isothermal_column_zero_cape():
    # an isothermal (absolutely stable) moist column: the parcel is
    # negatively buoyant everywhere above the LCL, so no equilibrium
    # level exists and CAPE is exactly zero
    nlev = 72
    ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
    pint = (225.5 + ai * (1.0e5 - 225.5)) * 1e-2
    pmid = 0.5 * (pint[:-1] + pint[1:])
    t = np.full(nlev, 280.0)
    es = 6.112 * np.exp(17.67 * (t - 273.15) / (t - 29.65))
    q = np.maximum(0.3 * 0.622 * es / (pmid - es), 1e-9)
    tvv = t * (1 + 0.608 * q)
    zmid = np.zeros(nlev)
    zint = np.zeros(nlev + 1)
    for k in range(nlev - 1, -1, -1):
        rog = 287.042 * tvv[k] / 9.80616
        zmid[k] = zint[k + 1] + rog * np.log(pint[k + 1] / pmid[k])
        zint[k] = zint[k + 1] + rog * np.log(pint[k + 1] / pint[k])
    zc = zm_cape.make_zm_const()
    zp = zm_cape.make_zm_param()
    r = zm_cape.compute_dilute_cape(
        q[None, :], t[None, :], zmid[None, :], pmid[None, :],
        pint[None, :], np.array([65]), np.array([0.0]), 11, zc, zp)
    assert float(r["cape"][0]) == 0.0
