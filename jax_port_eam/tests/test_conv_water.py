"""Tier-1 golden replay + Tier-0 properties for conv_water_4rad
(eam_jax.conv_water).

Golden archive: harness/gen_conv_water_golden.py (run in the
scream-dev container) - 16 columns x 72 levels of stratified-random
per-point states, 5 configurations (a: EAMv3 production default
mode=1/zm_microp=T/P3, b: mode=1 non-microp, c: mode=2 emissivity
average, d: pergro repartition, e: RK kabsi in mode 2).

Tolerances: mode 1 and the microp branch replay at <= 1e-13 relative
(exact branch arithmetic). The mode-2 configs (c, e) evaluate
log(sum(f*exp(alpha*w)))/alpha; when the in-cloud water sits just
above ic_limit the log argument is 1 + alpha*w with |alpha*w| ~ 1e-8,
and a 1-ulp container-glibc vs host-libm/XLA difference in exp/log
(measured: host libm == JAX != container Fortran on 2/1152 points)
is amplified by 1/|alpha*w| - measured max 2.25e-8 relative on values
of ~7e-14 kg/kg. The absolute error is bounded by ulp(1)/|alpha| <
4e-19 kg/kg, so mode-2 fields assert rtol=1e-12 with atol=1e-17
(orders below any physically significant condensate).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "conv_water_golden.npz"

from eam_jax import conv_water  # noqa: E402

CFGS = ["a", "b", "c", "d", "e"]
NAMES = ("totg_liq", "totg_ice", "sh_cldliq", "sh_cldice")


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


def _run(gold, cfg_params):
    return conv_water.conv_water_4rad(
        gold["t"], gold["pdel"], gold["q_cldliq"], gold["q_cldice"],
        gold["sh_icwmr"], gold["dp_icwmr"], gold["dp_icimr"],
        gold["fice"], gold["sh_frac"], gold["dp_frac"], gold["ast"],
        gold["rei"],
        conv_water_mode=cfg_params["conv_water_mode"],
        zm_microp=bool(cfg_params["zm_microp"]),
        microp_scheme=cfg_params["microp_scheme"],
        pergro_mods=cfg_params["pergro_mods"])


@pytest.fixture(scope="module")
def results(gold):
    meta = json.loads(str(gold["__metadata__"]))
    return {cfg: _run(gold, meta["params"][cfg]) for cfg in CFGS}


# ---------------------------------------------------------------------------
# Tier-1 golden replay
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("cfg", CFGS)
@pytest.mark.parametrize("field", NAMES)
def test_replay(gold, results, cfg, field):
    got = np.asarray(results[cfg][field])
    want = gold[f"{field}_{cfg}"]
    if cfg in ("c", "e") and field.startswith("totg"):
        # mode 2: exp/log 1-ulp cross-libm noise amplified by the
        # log(1+alpha*w)/alpha cancellation (see docstring)
        rtol, atol = 1e-12, 1e-17
    else:
        rtol, atol = 1e-13, 0.0
    np.testing.assert_allclose(got, want, rtol=rtol, atol=atol)


# ---------------------------------------------------------------------------
# Tier-0 properties
# ---------------------------------------------------------------------------
def test_outputs_nonnegative(results):
    for cfg in CFGS:
        for f in NAMES:
            v = np.asarray(results[cfg][f])
            assert np.isfinite(v).all(), (cfg, f)
            assert v.min() >= 0.0, (cfg, f)


def test_microp_branch_ignores_mode(gold):
    """The zm_microp partition never reads conv_water_mode."""
    meta = json.loads(str(gold["__metadata__"]))
    p = dict(meta["params"]["a"])
    r1 = _run(gold, {**p, "conv_water_mode": 1})
    r2 = _run(gold, {**p, "conv_water_mode": 2})
    for f in NAMES:
        np.testing.assert_array_equal(np.asarray(r1[f]),
                                      np.asarray(r2[f]))


@pytest.mark.parametrize("mode", [1, 2])
def test_no_convection_reduces_to_stratiform(gold, mode):
    """With zero convective cloud everywhere the non-microp result is
    the pure stratiform in-cloud value repartitioned by wrk1:
    totg_liq+totg_ice = (ql+qi)*ast/max(0.01, ast) for ast >= 0.01,
    exactly zero below the fraction threshold."""
    z = np.zeros_like(gold["t"])
    r = conv_water.conv_water_4rad(
        gold["t"], gold["pdel"], gold["q_cldliq"], gold["q_cldice"],
        z, z, z, gold["fice"], z, z, gold["ast"], gold["rei"],
        conv_water_mode=mode, zm_microp=False)
    tot = np.asarray(r["totg_liq"]) + np.asarray(r["totg_ice"])
    ql, qi, ast = gold["q_cldliq"], gold["q_cldice"], gold["ast"]
    want = np.where(ast >= 0.01,
                    (ql + qi) / np.maximum(0.01, ast) * ast, 0.0)
    np.testing.assert_allclose(tot, want, rtol=1e-13, atol=1e-30)
    assert np.all(tot[ast < 0.01] == 0.0)
    assert np.all(np.asarray(r["sh_cldliq"]) == 0.0)
    assert np.all(np.asarray(r["sh_cldice"]) == 0.0)


def test_cosp_outputs_identity_and_nan_guard(gold, results):
    """sh_cldliq + sh_cldice == sh_icwmr*sh_frac wherever FICE is
    defined; NaN FICE points are exactly zero (all configs share the
    same COSP outputs)."""
    nanm = np.isnan(gold["fice"])
    for cfg in CFGS:
        liq = np.asarray(results[cfg]["sh_cldliq"])
        ice = np.asarray(results[cfg]["sh_cldice"])
        assert np.all(liq[nanm] == 0.0) and np.all(ice[nanm] == 0.0)
        np.testing.assert_allclose(
            (liq + ice)[~nanm],
            (gold["sh_icwmr"] * gold["sh_frac"])[~nanm],
            rtol=1e-13, atol=1e-30)


def test_pergro_repartitions_only_condensate_starved_points(gold):
    """pergro_mods may only change points whose stratiform in-cloud
    condensate is below 100*ic_limit, and there it partitions by
    temperature: all ice below 243 K, all liquid above 263 K."""
    meta = json.loads(str(gold["__metadata__"]))
    p = meta["params"]["b"]
    rb = _run(gold, p)
    rd = _run(gold, {**p, "pergro_mods": True})
    dliq = np.asarray(rd["totg_liq"]) - np.asarray(rb["totg_liq"])
    changed = dliq != 0.0
    assert changed.any()
    ast = gold["ast"]
    ls_icwmr = np.where(ast < 0.01, 0.0,
                        (gold["q_cldliq"] + gold["q_cldice"])
                        / np.maximum(0.01, ast))
    # points with significant stratiform condensate are untouched
    assert not changed[ls_icwmr >= 1e-10].any()
    t = gold["t"]
    cold = changed & (t < 243.0)
    warm = changed & (t >= 263.0)
    np.testing.assert_array_equal(np.asarray(rd["totg_liq"])[cold], 0.0)
    np.testing.assert_array_equal(np.asarray(rd["totg_ice"])[warm], 0.0)
