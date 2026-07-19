"""Tests for the approximation-by-identity smoothing operators and
their first adoption site (cld_fraction).

Contract under test:
1. width = 0 is BITWISE the exact hard operator (no surrogate math).
2. width > 0 converges to the hard operator as width -> 0 (away from
   the switch point), with error O(width * scale).
3. width > 0 gives finite, nonzero gradients where the hard operator's
   gradient is zero (jumps) or undefined (kinks).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

from scream_jax.foundation import smoothing  # noqa: E402
from scream_jax.cld_fraction.main import cld_fraction_main  # noqa: E402


RNG = np.random.default_rng(42)
X = np.concatenate([RNG.normal(size=200), [0.0, -0.0, 1e-300, -1e-300]])
A = RNG.normal(size=204)
B = np.concatenate([RNG.normal(size=200), A[200:204]])  # incl. exact ties


# ---------------------------------------------------------------------------
# 1. width = 0 is bitwise exact
# ---------------------------------------------------------------------------
def test_width_zero_bitwise():
    np.testing.assert_array_equal(
        np.asarray(smoothing.step(X, 0.0)), (X > 0).astype(np.float64))
    np.testing.assert_array_equal(
        np.asarray(smoothing.blend(X, A, B, 0.0)), np.where(X > 0, A, B))
    np.testing.assert_array_equal(
        np.asarray(smoothing.smooth_abs(X, 0.0)), np.abs(X))
    np.testing.assert_array_equal(
        np.asarray(smoothing.smooth_max(A, B, 0.0)), np.maximum(A, B))
    np.testing.assert_array_equal(
        np.asarray(smoothing.smooth_min(A, B, 0.0)), np.minimum(A, B))
    np.testing.assert_array_equal(
        np.asarray(smoothing.smooth_clip(X, -0.5, 0.5, 0.0)),
        np.clip(X, -0.5, 0.5))


def test_cld_fraction_default_is_original():
    qi = RNG.uniform(0, 3e-12, (16, 72))
    liq = RNG.uniform(0, 1, (16, 72))
    thresh, thresh4 = 1e-12, 1e-5
    got = cld_fraction_main(thresh, thresh4, qi, liq)
    one, zero = 1.0, 0.0
    ice = np.where(qi > thresh, one, zero)
    ice4 = np.where(qi > thresh4, one, zero)
    np.testing.assert_array_equal(np.asarray(got[0]), ice)
    np.testing.assert_array_equal(np.asarray(got[1]), np.maximum(ice, liq))
    np.testing.assert_array_equal(np.asarray(got[2]), ice4)
    np.testing.assert_array_equal(np.asarray(got[3]),
                                  np.maximum(ice4, liq))


# ---------------------------------------------------------------------------
# 2. convergence as width -> 0
# ---------------------------------------------------------------------------
def test_convergence_away_from_switch():
    # points at least one width away from the switch: error bounded and
    # shrinking with width
    x = np.linspace(-2.0, 2.0, 401)
    x = x[np.abs(x) > 0.2]
    hard = (x > 0).astype(float)
    for w_big, w_small in [(0.1, 0.025), (0.025, 0.00625)]:
        e_big = np.abs(np.asarray(smoothing.step(x, w_big)) - hard).max()
        e_small = np.abs(np.asarray(smoothing.step(x, w_small))
                         - hard).max()
        assert e_small < e_big < 0.15  # sigmoid(-|x|min/w) at the edge

    for w in (0.1, 0.01):
        err = np.abs(np.asarray(smoothing.smooth_abs(x, w)) - np.abs(x))
        assert err.max() <= w  # sqrt(x^2+w^2) - |x| <= w


def test_smooth_max_bias_bound():
    a = RNG.normal(size=500)
    b = RNG.normal(size=500)
    w = 0.05
    sm = np.asarray(smoothing.smooth_max(a, b, w))
    hard = np.maximum(a, b)
    assert ((sm - hard) >= -1e-15).all()          # smooth_max >= max
    assert (sm - hard).max() <= 0.5 * w + 1e-15   # bias bound at ties


# ---------------------------------------------------------------------------
# 3. gradients: finite and nonzero where the hard ops are flat/kinked
# ---------------------------------------------------------------------------
def test_step_gradient_at_switch():
    for w, scale in [(0.1, 1.0), (0.05, 2.0)]:
        g = float(jax.grad(
            lambda s: smoothing.step(s, w, scale=scale))(0.0))
        np.testing.assert_allclose(g, 1.0 / (4.0 * w * scale), rtol=1e-12)


def test_smooth_ops_c1_at_kinks():
    # gradient finite (no NaN) exactly at the kink/tie points
    g_abs = float(jax.grad(lambda x: smoothing.smooth_abs(x, 0.1))(0.0))
    assert np.isfinite(g_abs) and abs(g_abs) < 1e-12  # symmetric at 0
    g_max = jax.grad(
        lambda a: smoothing.smooth_max(a, 1.0, 0.1))(1.0)
    assert np.isfinite(float(g_max))
    np.testing.assert_allclose(float(g_max), 0.5, rtol=1e-12)  # symmetry


def test_cld_fraction_smooth_gradient():
    qi = RNG.uniform(0, 3e-12, (8, 72))
    liq = RNG.uniform(0, 1, (8, 72))
    thresh, thresh4 = 1e-12, 1e-5

    def total(q, w):
        return jnp.sum(cld_fraction_main(thresh, thresh4, q, liq,
                                         smooth_width=w)[1])

    g_hard = np.asarray(jax.grad(lambda q: total(q, 0.0))(jnp.asarray(qi)))
    assert (g_hard == 0.0).all()          # the baseline pathology

    g_smooth = np.asarray(jax.grad(lambda q: total(q, 0.1))(
        jnp.asarray(qi)))
    assert np.isfinite(g_smooth).all()
    assert (np.abs(g_smooth) > 0).mean() > 0.5   # gradient actually flows
    assert g_smooth.min() >= 0.0   # more ice can only increase cloud

    # surrogate primal converges to the hard scheme as width -> 0
    hard = float(total(jnp.asarray(qi), 0.0))
    for w, tol in [(0.1, 0.35), (0.01, 0.05)]:
        rel = abs(float(total(jnp.asarray(qi), w)) - hard) / hard
        assert rel < tol, (w, rel)


def test_smooth_cld_fraction_range_note():
    # documented surrogate bias: smooth_max can exceed 1 by <= width/2
    qi = np.full((4, 4), 1e-10)
    liq = np.ones((4, 4))
    out = np.asarray(cld_fraction_main(1e-12, 1e-5, qi, liq,
                                       smooth_width=0.1)[1])
    assert out.max() <= 1.0 + 0.05 + 1e-12
