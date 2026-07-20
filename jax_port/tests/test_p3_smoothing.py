"""Tests for the opt-in P3 jump-site smoothing (smooth_width kwarg).

Guards, on golden-derived states:
  1. default bitwise: p3_process_step(smooth_width=0.0) is array_equal
     to a call without the kwarg (the wider suite + strict p3 goldens
     guard the default path itself);
  2. an empty smooth_families tuple with a nonzero width is also exactly
     the hard path (family threading resolves to width 0 everywhere);
  3. smooth_width=0.3 runs, outputs are finite, and the reverse-mode
     grad through the full process step (scan sedimentation) is
     finite/NaN-free and differs from the hard grad (the surrogate
     actually re-routes gradient through the jump sites).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402

import scream_jax  # noqa: E402,F401  (enables x64)
from scream_jax.p3 import DEFAULT_OPTS, tables  # noqa: E402
from scream_jax.p3.process import p3_process_step  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / \
    "physics_suite_218x72_dt1800_2steps.npz"
TDIR = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "tables"

NCOL = 8
DT = 300.0


@pytest.fixture(scope="module")
def state():
    if not GOLDEN.exists():
        pytest.skip(f"golden archive not found: {GOLDEN}")
    if not (TDIR / f"p3_lookup_table_1.dat-v{tables.P3_VERSION}").exists():
        pytest.skip("P3 tables not available")
    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    s = {name: np.asarray(z[f"{name}__step0"], dtype=np.float64)[:NCOL]
         for name in meta["fields"]}
    tbl = tables.p3_init(str(TDIR))
    return s, tbl, dict(DEFAULT_OPTS)


def _run(s, tbl, opts, qc=None, **kw):
    return p3_process_step(
        DT, True, True, True, False, False, False, False, False,
        s["T_mid"], s["p_mid"], s["p_dry_mid"], s["pseudo_density"],
        s["pseudo_density_dry"], s["cldfrac_tot"],
        s["qv"], s["qc"] if qc is None else qc, s["nc"],
        s["qr"], s["nr"], s["qi"], s["qm"], s["ni"], s["bm"],
        s["qv_prev_micro_step"], s["T_prev_micro_step"],
        s["nc_nuceat_tend"], s["nccn"], s["ni_activated"],
        s["inv_qc_relvar"], s["precip_liq_surf_mass"],
        s["precip_ice_surf_mass"], tbl, opts, **kw)


def _assert_all_equal(a, b, ctx):
    assert set(a) == set(b)
    for k in a:
        av, bv = np.asarray(a[k]), np.asarray(b[k])
        assert np.array_equal(av, bv, equal_nan=True), (
            f"{ctx}: output '{k}' differs "
            f"(max abs diff {np.nanmax(np.abs(av - bv)):.3e})")


def test_smooth_width_zero_matches_default(state):
    s, tbl, opts = state
    a = _run(s, tbl, opts)
    b = _run(s, tbl, opts, smooth_width=0.0)
    _assert_all_equal(a, b, "smooth_width=0.0 vs no kwarg")


def test_empty_family_tuple_is_hard_path(state):
    s, tbl, opts = state
    a = _run(s, tbl, opts)
    b = _run(s, tbl, opts, smooth_width=0.3, smooth_families=())
    _assert_all_equal(a, b, "smooth_width=0.3, smooth_families=()")


@pytest.mark.parametrize("width", [0.02, 0.3])
def test_smooth_width_runs_finite_and_grad_nan_free(state, width):
    # width=0.02 exercises the deep sigmoid tails (values ~1e-150) that
    # can overflow guarded divisions in the VJP; width=0.3 the wide
    # transition regime.
    s, tbl, opts = state
    out = _run(s, tbl, opts, smooth_width=width)
    for k, val in out.items():
        arr = np.asarray(val)
        assert np.isfinite(arr).all(), (
            f"smooth_width={width}: non-finite values in output '{k}'")

    def scalar(qc, w):
        o = _run(s, tbl, opts, qc=qc, sed_use_while_loop=False,
                 smooth_width=w)
        return jnp.sum(o["qc"] + o["qr"]) + jnp.sum(o["T_mid"])

    x0 = jnp.asarray(s["qc"])
    g_smooth = np.asarray(jax.grad(scalar)(x0, width))
    assert int(np.isnan(g_smooth).sum()) == 0, (
        f"NaN in smoothed grad at width {width}")
    assert int(np.isinf(g_smooth).sum()) == 0, (
        f"Inf in smoothed grad at width {width}")
    assert np.linalg.norm(g_smooth) > 0.0

    # the surrogate must actually change the gradient signal
    g_hard = np.asarray(jax.grad(scalar)(x0, 0.0))
    assert not np.array_equal(g_smooth, g_hard), (
        f"smooth_width={width} grad identical to hard grad — smoothing "
        "not reached by the qc-dependent paths")
