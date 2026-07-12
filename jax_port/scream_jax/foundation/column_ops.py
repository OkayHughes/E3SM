"""Vertical column operators (midpoint/interface conversions, scans).

Transcribed from components/eamxx/src/share/util/eamxx_column_ops.hpp
(struct scream::ColumnOps). The pack_size==1 code paths are the semantic
spec; all EKAT pack machinery collapses in JAX.

Conventions (identical to EAMxx):
- The level axis is the LAST axis. Midpoint arrays have nlev entries,
  interface arrays nlev+1.
- k=0 is the MODEL TOP; k=nlev is the surface interface.

The C++ CombineMode output-blending argument (y = beta*y + alpha*f(x)) is a
Kokkos in-place-output detail with no JAX analogue: compose instead.
"""

import jax.numpy as jnp


def compute_midpoint_values(x_i):
    """Midpoint values from interface values: x_m(k) = (x_i(k)+x_i(k+1))/2."""
    x_i = jnp.asarray(x_i)
    return 0.5 * (x_i[..., :-1] + x_i[..., 1:])


def compute_midpoint_delta(x_i):
    """Forward difference onto midpoints: dx_m(k) = x_i(k+1) - x_i(k)."""
    x_i = jnp.asarray(x_i)
    return x_i[..., 1:] - x_i[..., :-1]


def column_scan(dx_m, s0=0.0, *, from_top: bool):
    """Scan sum of a midpoint quantity, yielding its interface integral.

    from_top=True:  x_i(0) = s0,      x_i(k+1) = s0 + sum_{n=0..k}   dx_m(n)
    from_top=False: x_i(nlev) = s0,   x_i(k)   = s0 + sum_{n=k..N-1} dx_m(n)

    (logical inverse of compute_midpoint_delta; C++ column_scan<FromTop>).
    s0 may be a scalar or an array broadcastable to the column dimensions.
    """
    dx_m = jnp.asarray(dx_m)
    s0 = jnp.broadcast_to(jnp.asarray(s0, dtype=dx_m.dtype), dx_m.shape[:-1])[..., None]
    if from_top:
        return jnp.concatenate([s0, s0 + jnp.cumsum(dx_m, axis=-1)], axis=-1)
    suffix = jnp.cumsum(dx_m[..., ::-1], axis=-1)[..., ::-1]
    return jnp.concatenate([s0 + suffix, s0], axis=-1)


def compute_interface_values_linear(x_m, dz, bc_top, bc_bot):
    """Interface values from midpoints via dz-weighted linear interpolation.

    Interior: x_i(k) = (x_m(k)*dz(k-1) + x_m(k-1)*dz(k)) / (dz(k-1) + dz(k));
    boundaries imposed. Matches ColumnOps::compute_interface_values_linear.
    """
    x_m, dz = jnp.asarray(x_m), jnp.asarray(dz)
    interior = ((x_m[..., 1:] * dz[..., :-1] + x_m[..., :-1] * dz[..., 1:])
                / (dz[..., :-1] + dz[..., 1:]))
    top = jnp.broadcast_to(jnp.asarray(bc_top, dtype=x_m.dtype), x_m.shape[:-1])[..., None]
    bot = jnp.broadcast_to(jnp.asarray(bc_bot, dtype=x_m.dtype), x_m.shape[:-1])[..., None]
    return jnp.concatenate([top, interior, bot], axis=-1)


def compute_interface_values_compatible(x_m, bc, *, fix_top: bool):
    """Interface values from midpoints such that midpoints of the result
    reproduce x_m exactly (x_m(x_i(x_m)) == x_m).

    One boundary value bc is imposed (top if fix_top else bottom); the rest
    follow from the alternating-sign scan (C++
    compute_interface_values_compatible<FixTop>):

        fix_top:  x_i(k+1) = (-1)^k [ -bc + 2*sum_{n=0..k}   (-1)^n x_m(n) ]
        fix_bot:  x_i(k)   = (-1)^k [ (-1)^N bc + 2*sum_{n=k..N-1} (-1)^n x_m(n) ]

    CAVEAT (as in the C++): not monotonic — can create interface
    maxima/minima larger than any input midpoint value.
    """
    x_m = jnp.asarray(x_m)
    n = x_m.shape[-1]
    sign_m = 1 - 2 * (jnp.arange(n) % 2)          # (-1)^k at midpoints
    sign_i = 1 - 2 * (jnp.arange(n + 1) % 2)      # (-1)^k at interfaces
    scanned = column_scan(2 * sign_m * x_m, 0.0, from_top=fix_top)
    bc = jnp.asarray(bc, dtype=x_m.dtype)
    if fix_top:
        return (scanned - bc[..., None] if bc.ndim else scanned - bc) * (-sign_i)
    sign_n = 1 - 2 * (n % 2)                       # (-1)^N
    return (scanned + sign_n * (bc[..., None] if bc.ndim else bc)) * sign_i
