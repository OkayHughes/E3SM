"""Linear vertical interpolation between the zt and zi grids.

Source: components/eamxx/src/physics/shoc/impl/shoc_linear_interp_impl.hpp
"""

import functools

import jax
import jax.numpy as jnp


@functools.partial(jax.jit, static_argnames=("km1", "km2"))
def linear_interp(x1, x2, y1, km1: int, km2: int, minthresh):
    """Interpolate y1(x1) onto x2, then floor at minthresh.

    Two supported shapes (matching the C++, which hard-errors otherwise):
    - km1 == km2+1 (e.g. interfaces -> midpoints): target k2 interpolates
      between source points (k2, k2+1).
    - km2 == km1+1 (e.g. midpoints -> interfaces): interior target k2 uses
      source points (k2-1, k2); the k2=0 and k2=km2-1 ends *extrapolate*
      from the first/last source pair.

    x1, y1 have km1 level entries; x2 (and the result) km2. Leading axes
    broadcast. minthresh floors the result (use shoc constant `largeneg`
    for "no floor").
    """
    x1 = jnp.asarray(x1)
    x2 = jnp.asarray(x2)
    y1 = jnp.asarray(y1)

    if km1 == km2 + 1:
        # idx = k2+1 for every target level: pairs (k2, k2+1)
        lo = slice(0, km2)
        hi = slice(1, km2 + 1)
        x_lo, x_hi = x1[..., lo], x1[..., hi]
        y_lo, y_hi = y1[..., lo], y1[..., hi]
    elif km2 == km1 + 1:
        # idx = clamp(k2, 1, km1-1): pairs (idx-1, idx), extrapolating ends
        idx = jnp.clip(jnp.arange(km2), 1, km1 - 1)
        x_lo, x_hi = x1[..., idx - 1], x1[..., idx]
        y_lo, y_hi = y1[..., idx - 1], y1[..., idx]
    else:
        raise ValueError(f"Unsupported dimensions for linear interp: {km1=}, {km2=}")

    y2 = y_lo + (y_hi - y_lo) * (x2 - x_lo) / (x_hi - x_lo)
    return jnp.maximum(y2, minthresh)
