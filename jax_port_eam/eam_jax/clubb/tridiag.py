"""CLUBB tridiagonal solver — exact port of LAPACK dgtsv.

PORT_NOTES
----------
CLUBB's lapack_wrap.F90 tridag_solve is a thin wrapper around LAPACK
DGTSV (reference netlib LAPACK; the scream-dev container links Ubuntu's
reference liblapack 3.12, which the goldens were generated against).
DGTSV is Gaussian elimination with partial pivoting on the three
diagonals, non-blocked, followed by back substitution — deterministic
and portable, so a faithful re-implementation replays it bitwise for
finite, nonsingular systems.

Algorithm (dgtsv.f, NRHS>=2 branch; the NRHS=1 branch differs only in
that a zero subdiagonal element skips the trivially-zero update, which
changes nothing for finite data — we always perform the arithmetic,
i.e. the NRHS=1 code structure):

  for i = 1 .. n-1:
    if |d(i)| >= |dl(i)|:            # no row interchange
      fact    = dl(i)/d(i)
      d(i+1) -= fact*du(i)
      b(i+1,:) -= fact*b(i,:)
      du2(i)  = 0
    else:                            # interchange rows i and i+1
      fact    = d(i)/dl(i)
      d(i)    = dl(i)
      temp    = d(i+1)
      d(i+1)  = du(i) - fact*temp
      if i < n-1:
        du2(i)   = du(i+1)
        du(i+1)  = -fact*du2(i)
      du(i)   = temp
      b(i,:), b(i+1,:) = b(i+1,:), b(i,:) - fact*b(i+1,:)

  back substitution:
    x(n)   = b(n)/d(n)
    x(n-1) = (b(n-1) - du(n-1)*x(n))/d(n-1)
    x(i)   = (b(i) - du(i)*x(i+1) - du2(i)*x(i+2))/d(i),  i = n-2..1

Interface mirrors tridag_solve: `supd`/`diag`/`subd` are the length-n
CLUBB arrays where supd(1:n-1) is DU and subd(2:n) is DL (supd(n) and
subd(1) are ignored, exactly as lapack_wrap slices them).

Singular systems: LAPACK returns info>0 and CLUBB sets a fatal error
and solution = -999.  The port returns (solution, singular_flag);
callers/tests treat singular_flag=True as the error path (the -999
fill is applied here too for parity).
"""

import jax
import jax.numpy as jnp


def _dgtsv(dl, d, du, b):
    """dl, d, du: (n,), (n,), (n,) with dl[0] and du[-1] unused;
    b: (n, nrhs). Returns (x, singular)."""
    n = d.shape[0]

    def fwd(carry, inp):
        d_cur, du_cur, b_cur, sing = carry
        dl_i, d_next, du_next, b_next, has_next = inp

        no_pivot = jnp.abs(d_cur) >= jnp.abs(dl_i)

        # --- no interchange ---
        fact_np = dl_i / d_cur
        d_next_np = d_next - fact_np * du_cur
        b_next_np = b_next - fact_np * b_cur
        sing_np = sing | (d_cur == 0.0)

        # --- interchange rows ---
        fact_p = d_cur / dl_i
        d_i_p = dl_i
        temp = d_next
        d_next_p = du_cur - fact_p * temp
        du2_p = jnp.where(has_next, du_next, 0.0)
        du_next_p = jnp.where(has_next, -fact_p * du_next, du_next)
        du_i_p = temp
        b_i_p = b_next
        b_next_p = b_cur - fact_p * b_next

        d_i = jnp.where(no_pivot, d_cur, d_i_p)
        du_i = jnp.where(no_pivot, du_cur, du_i_p)
        du2_i = jnp.where(no_pivot, 0.0, du2_p)
        b_i = jnp.where(no_pivot, b_cur, b_i_p)
        d_nx = jnp.where(no_pivot, d_next_np, d_next_p)
        du_nx = jnp.where(no_pivot, du_next, du_next_p)
        b_nx = jnp.where(no_pivot, b_next_np, b_next_p)
        sing = jnp.where(no_pivot, sing_np, sing)

        return (d_nx, du_nx, b_nx, sing), (d_i, du_i, du2_i, b_i)

    has_next = jnp.arange(n - 1) < n - 2
    du_next_in = jnp.where(has_next, du[1:], 0.0)
    carry0 = (d[0], du[0], b[0], jnp.asarray(False))
    (d_n, _du_n, b_n, singular), (d_out, du_out, du2_out, b_out) = \
        jax.lax.scan(fwd, carry0,
                     (dl[1:], d[1:], du_next_in, b[1:], has_next))
    singular = singular | (d_n == 0.0)

    d_f = jnp.concatenate([d_out, d_n[None]])
    b_f = jnp.concatenate([b_out, b_n[None]])
    # du_out[i] is the final du(i+1) for i=0..n-2; du2_out likewise.

    x_n = b_f[n - 1] / d_f[n - 1]
    x_nm1 = (b_f[n - 2] - du_out[n - 2] * x_n) / d_f[n - 2]

    def back(carry, inp):
        x_ip1, x_ip2 = carry
        b_i, du_i, du2_i, d_i = inp
        x_i = (b_i - du_i * x_ip1 - du2_i * x_ip2) / d_i
        return (x_i, x_ip1), x_i

    if n > 2:
        _, xs = jax.lax.scan(
            back, (x_nm1, x_n),
            (b_f[n - 3::-1], du_out[n - 3::-1], du2_out[n - 3::-1],
             d_f[n - 3::-1]))
        x = jnp.concatenate([xs[::-1], x_nm1[None], x_n[None]])
    else:
        x = jnp.concatenate([x_nm1[None], x_n[None]])

    return x, singular


def tridag_solve(supd, diag, subd, rhs):
    """CLUBB lapack_wrap tridag_solve.

    supd/diag/subd: (n,) CLUBB-layout diagonals (supd[n-1], subd[0]
    unused); rhs: (n,) or (n, nrhs).  Returns (solution, singular);
    solution is -999 everywhere when singular (CLUBB's error fill).
    """
    supd = jnp.asarray(supd, dtype=jnp.float64)
    diag = jnp.asarray(diag, dtype=jnp.float64)
    subd = jnp.asarray(subd, dtype=jnp.float64)
    rhs = jnp.asarray(rhs, dtype=jnp.float64)
    squeeze = rhs.ndim == 1
    if squeeze:
        rhs = rhs[:, None]
    x, singular = _dgtsv(subd, diag, supd, rhs)
    x = jnp.where(singular, -999.0, x)
    if squeeze:
        x = x[:, 0]
    return x, singular
