"""Approximation-by-identity smoothing operators.

Families of smooth surrogates for the discontinuous/non-smooth
primitives in the physics (Heaviside branches, min/max/clip, |x|),
parameterized by a width so that **width = 0 recovers the exact hard
operator bitwise** (an explicit Python-level branch returns the
original jnp op — no sigmoid evaluated at all). With width > 0 the
operators are C-infinity, giving automatic differentiation a nonzero,
finite gradient through what is otherwise a jump (where AD reports 0)
or a kink (where AD picks a subgradient).

Conventions
-----------
- Every op takes `width` (dimensionless) and `scale` (the
  characteristic magnitude of the switching variable at the call
  site); the effective transition half-width is ``width * scale``.
  Callers pick `scale` once per site (e.g. the threshold itself for a
  ``q > q_thresh`` gate); `width` is the user's single knob.
- ``width`` must be a static Python float when used inside ``jax.jit``
  (jit with ``static_argnames`` at adoption sites) — the width-0 exact
  branch is resolved at trace time.
- Surrogate error is O(width * scale) away from the switch point;
  ``smooth_max(a, b)`` has a +width*scale/2 bias exactly at a == b.

Adoption pattern (keeps the default path bitwise-original)::

    if smooth_width == 0.0:
        f = jnp.where(q > thresh, on, off)           # original line
    else:
        f = smoothing.blend(q - thresh, on, off,
                            scale=thresh, width=smooth_width)
"""

import jax
import jax.numpy as jnp


def step(s, width, scale=1.0):
    """Heaviside family: sigmoid(s / (width*scale)); exact ``s > 0``
    indicator (as float) at width 0."""
    if width == 0.0:
        return (s > 0).astype(jnp.result_type(s, float))
    return jax.nn.sigmoid(s / (width * scale))


def blend(s, on, off, width, scale=1.0):
    """``where(s > 0, on, off)`` family, smoothly mixing via `step`.

    NaN-safe under reverse mode by construction only if `on`/`off` are
    finite everywhere they are evaluated — sanitize operands first
    (this replaces the *selection*, not the double-where idiom)."""
    if width == 0.0:
        return jnp.where(s > 0, on, off)
    h = step(s, width, scale)
    return off + h * (on - off)


def smooth_abs(x, width, scale=1.0):
    """|x| family: sqrt(x^2 + (width*scale)^2); exact ``jnp.abs`` at
    width 0. Bias +width*scale at x == 0."""
    if width == 0.0:
        return jnp.abs(x)
    w = width * scale
    return jnp.sqrt(x * x + w * w)


def smooth_max(a, b, width, scale=1.0):
    """max family: (a + b + |a-b|_smooth)/2; exact ``jnp.maximum`` at
    width 0."""
    if width == 0.0:
        return jnp.maximum(a, b)
    return 0.5 * (a + b + smooth_abs(a - b, width, scale))


def smooth_min(a, b, width, scale=1.0):
    """min family: (a + b - |a-b|_smooth)/2; exact ``jnp.minimum`` at
    width 0."""
    if width == 0.0:
        return jnp.minimum(a, b)
    return 0.5 * (a + b - smooth_abs(a - b, width, scale))


def smooth_relu(x, width, scale=1.0):
    """Positive-part family: exact ``maximum(x, 0)`` at width 0."""
    return smooth_max(x, jnp.zeros_like(x), width, scale)


def smooth_clip(x, lo, hi, width, scale=1.0):
    """clip family: exact ``jnp.clip`` at width 0 (clip = min(max(x,
    lo), hi), applied in that order as jnp.clip does)."""
    if width == 0.0:
        return jnp.clip(x, lo, hi)
    return smooth_min(smooth_max(x, lo, width, scale), hi, width, scale)
