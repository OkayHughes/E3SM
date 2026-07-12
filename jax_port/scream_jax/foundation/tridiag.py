"""Batched tridiagonal solver (Thomas algorithm).

Transcribed from externals/ekat/src/algorithm/ekat_tridiag.hpp, serial
implementation `impl::thomas_a1x1` (the team-parallel `thomas`, `cr`
(cyclic reduction) and `bfb` variants perform the same arithmetic with
different parallel decompositions — the serial ordering is the reference).

Used by SHOC's implicit vertical diffusion (vd_shoc_decomp/vd_shoc_solve in
components/eamxx/src/physics/shoc/impl/shoc_tridiag_solver_impl.hpp).

Conventions:
- The system axis (vertical levels) is the LAST axis, matching the rest of
  scream_jax. Leading axes are batch (columns, RHS, ...).
- dl, d, du, x all have shape (..., n); dl[..., 0] and du[..., n-1] are
  unused, exactly as in EKAT.
- Pure function: returns the solution; inputs are not modified.
"""

import jax.numpy as jnp
from jax import lax


def thomas(dl, d, du, x):
    """Solve T @ sol = x for a tridiagonal T, batched over leading axes.

    T has subdiagonal dl[..., 1:], diagonal d, superdiagonal du[..., :-1].
    No pivoting (as in EKAT) — caller guarantees the system is well
    conditioned (SHOC's diffusion matrices are diagonally dominant).
    """
    dl, d, du, x = jnp.broadcast_arrays(*map(jnp.asarray, (dl, d, du, x)))

    # Scan along the system axis: move it to the front.
    dl_f = jnp.moveaxis(dl, -1, 0)
    d_f = jnp.moveaxis(d, -1, 0)
    du_f = jnp.moveaxis(du, -1, 0)
    x_f = jnp.moveaxis(x, -1, 0)

    # Forward elimination:
    #   dli   = dl[i] / d'[i-1]
    #   d'[i] = d[i] - dli * du[i-1]
    #   x'[i] = x[i] - dli * x'[i-1]
    def fwd(carry, row):
        d_prev, x_prev = carry
        dl_i, d_i, du_im1, x_i = row
        dli = dl_i / d_prev
        d_new = d_i - dli * du_im1
        x_new = x_i - dli * x_prev
        return (d_new, x_new), (d_new, x_new)

    (_, _), (d_rest, x_rest) = lax.scan(
        fwd, (d_f[0], x_f[0]), (dl_f[1:], d_f[1:], du_f[:-1], x_f[1:]))
    d_fac = jnp.concatenate([d_f[:1], d_rest], axis=0)
    x_fac = jnp.concatenate([x_f[:1], x_rest], axis=0)

    # Back substitution:
    #   sol[n-1] = x'[n-1] / d'[n-1]
    #   sol[j]   = (x'[j] - du[j] * sol[j+1]) / d'[j]   for j = n-2 .. 0
    sol_last = x_fac[-1] / d_fac[-1]

    def bwd(sol_next, row):
        x_j, du_j, d_j = row
        sol_j = (x_j - du_j * sol_next) / d_j
        return sol_j, sol_j

    _, sol_rest = lax.scan(
        bwd, sol_last, (x_fac[:-1], du_f[:-1], d_fac[:-1]), reverse=True)
    sol = jnp.concatenate([sol_rest, sol_last[None]], axis=0)

    return jnp.moveaxis(sol, 0, -1)
