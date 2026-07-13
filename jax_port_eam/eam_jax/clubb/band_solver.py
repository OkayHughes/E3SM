"""CLUBB band-diagonal solver — exact port of lapack_wrap band_solve
(LAPACK dgbsv) for the advance_xm_wpxp 2x2-coupled system.

PORT_NOTES
----------
Source: components/eam/src/physics/clubb/lapack_wrap.F90 band_solve.

advance_xm_wpxp's xm_wpxp_solve calls band_solve (-> LAPACK DGBSV)
under the EAMv3 configuration: the band_solvex/DGBSVX
condition-estimate + equilibration + iterative-refinement path is only
taken when l_stats_samp is on AND the irtm/ithlm_matrix_condt_num
stats are registered — EAM runs l_stats=.false., so DGBSVX (and with
it dgbtrf/dgbcon/dgbtrs-with-refinement) is DEAD here.  Only DGBSV is
ported.

band_solve first restructures CLUBB's row-oriented band storage
  lhs(1..5, i):  equation i;  row 3 = coefficient of variable i,
                 rows 2/1 = variables i+1/i+2 (super-diagonals),
                 rows 4/5 = variables i-1/i-2 (sub-diagonals)
into LAPACK GB storage with kl=2 extra fill rows,
  AB(kl+ku+1+i-j, j) = A(i,j)   (1-based; LDAB = 2*kl+ku+1 = 7),
via the three verbatim index loops (equivalent to
AB0[r, c] = lhs0[r-2, c+r-4] where in range, else 0, in 0-based
indices), then calls DGBSV.

Reference DGBSV = DGBTRF + DGBTRS.  For kl = 2, DGBTRF's blocked path
is skipped (NB = ILAENV(1,'DGBTRF') = 32 > kl), so it reduces to the
unblocked DGBTF2:

  kv = kl + ku
  for j = 1..n:
     zero the fill-in column j+kv (rows 1..kl)   [already zero here]
     km = min(kl, n-j)
     jp = IDAMAX(km+1, AB(kv+1,j))               [first-of-ties argmax]
     ipiv(j) = jp + j - 1
     if AB(kv+jp,j) /= 0:
        ju = max(ju, min(j+ku+jp-1, n))
        if jp /= 1: banded row swap over columns j..ju
        if km > 0:
           DSCAL: multipliers *= 1/pivot        [reciprocal-MULTIPLY]
           if ju > j: DGER rank-1 update over columns j+1..ju,
                      SKIPPING columns whose pivot-row entry is 0
                      (the reference DGER y(j)==0 short-circuit)
     else: info = j (first zero pivot)

and DGBTRS ('No transpose'):

  L-solve: for j = 1..n-1: swap rhs rows (ipiv), then per rhs column
           temp = -b(j); b(j+1..j+lm) += multipliers * temp,
           SKIPPING columns with b(j) == 0 (DGER short-circuit)
  U-solve: DTBSV upper/no-transpose/non-unit with k = kl+ku,
           SKIPPING rows with x(j) == 0 (including the division).

Every skip-gate above is replicated with jnp.where so that signed
zeros and no-op columns replay exactly.  The scale step multiplies by
a precomputed reciprocal (NOT a division), exactly like DSCAL is
called.  Pivot search ties resolve to the FIRST index (IDAMAX uses a
strict > update), which is jnp.argmax's behavior.

Implementation: the column sweeps are lax.fori_loop's over a
column-padded AB / row-padded RHS so every dynamic_slice window is in
range; out-of-band touches are masked off exactly like the Fortran
loop bounds (masking selects, it never changes arithmetic).

Singular systems: LAPACK returns info > 0; CLUBB's band_solve then
sets a fatal error and never assigns `solution`.  The port returns
(solution, singular) with solution undefined (NaN-laden possible)
when singular; callers treat singular=True as the error path.

The container's Ubuntu liblapack is reference netlib LAPACK compiled
with default gcc options; like the dgtsv multi-RHS case in slice A,
FMA contraction inside its DGER/DTBSV loops makes the replay a few
ulps rather than bitwise — measured in tests/test_clubb_xm_wpxp.py.
"""

import jax.numpy as jnp
from jax import lax

KL = 2
KU = 2
KV = KL + KU            # 4: 0-based row of the diagonal in AB
LDAB = 2 * KL + KU + 1  # 7
NDIAG = KL + KU + 1     # 5: window width of a column's band coupling


def _assemble_lulhs(lhs):
    """CLUBB (5, n) row-oriented band storage -> LAPACK GB storage
    AB (7, n) with kl fill rows zeroed (band_solve's reorder loops)."""
    n = lhs.shape[1]
    cols = jnp.arange(n)
    rows = []
    for r in range(LDAB):
        if r < KL:
            rows.append(jnp.zeros(n, dtype=lhs.dtype))
        else:
            shift = r - KV           # -2..2
            src_col = cols + shift
            valid = (src_col >= 0) & (src_col < n)
            vals = lhs[r - KL, jnp.clip(src_col, 0, n - 1)]
            rows.append(jnp.where(valid, vals, 0.0))
    return jnp.stack(rows)


def _dgbtf2(ab):
    """Unblocked band LU with partial pivoting (reference DGBTF2),
    kl = ku = 2.  Returns (ab, ipiv (0-based row indices), info)."""
    n = ab.shape[1]
    abp = jnp.pad(ab, ((0, 0), (0, KL + KU)))   # in-range windows

    def body(j, carry):
        abp, ipiv, ju, info = carry
        block = lax.dynamic_slice(abp, (0, j), (LDAB, NDIAG))
        # IDAMAX over rows KV..KV+km of column j (first-of-ties);
        # candidates beyond the matrix are masked BELOW any |value|
        ivalid = j + jnp.arange(KL + 1) <= n - 1
        cand = jnp.where(ivalid, jnp.abs(block[KV:KV + KL + 1, 0]),
                         -1.0)
        jp = jnp.argmax(cand).astype(jnp.int32)
        ipiv = ipiv.at[j].set(j + jp)
        pivot = block[KV + jp, 0]
        nonzero = pivot != 0.0
        info = jnp.where((~nonzero) & (info == 0), j + 1, info)
        ju = jnp.where(
            nonzero,
            jnp.maximum(ju, jnp.minimum(j + KU + jp,
                                        jnp.int32(n - 1))), ju)
        w = ju - j                               # 0..KL+KU

        # banded row interchange over columns j..ju (jp != 0 only)
        do_swap = nonzero & (jp != 0)
        for m in range(NDIAG):
            hi = block[KV + jp - m, m]
            lo = block[KV - m, m]
            mask = do_swap & (m <= w)
            block = block.at[KV - m, m].set(jnp.where(mask, hi, lo))
            block = block.at[KV + jp - m, m].set(
                jnp.where(mask, lo, hi))

        # DSCAL: multipliers *= 1/pivot (reciprocal-multiply)
        recip = 1.0 / jnp.where(nonzero, pivot, 1.0)
        for i in range(1, KL + 1):
            valid = nonzero & (j + i <= n - 1)
            block = block.at[KV + i, 0].set(
                jnp.where(valid, block[KV + i, 0] * recip,
                          block[KV + i, 0]))

        # DGER rank-1 update over columns j+1..ju
        for mcol in range(1, NDIAG):
            y = block[KV - mcol, mcol]
            temp = -y
            colmask = nonzero & (mcol <= w) & (y != 0.0)
            for i in range(1, KL + 1):
                mask = colmask & (j + i <= n - 1)
                block = block.at[KV - mcol + i, mcol].set(
                    jnp.where(mask,
                              block[KV - mcol + i, mcol]
                              + block[KV + i, 0] * temp,
                              block[KV - mcol + i, mcol]))

        abp = lax.dynamic_update_slice(abp, block, (0, j))
        return abp, ipiv, ju, info

    ipiv0 = jnp.zeros(n, dtype=jnp.int32)
    abp, ipiv, _ju, info = lax.fori_loop(
        0, n, body,
        (abp, ipiv0, jnp.zeros((), jnp.int32),
         jnp.zeros((), jnp.int32)))
    return abp[:, :n], ipiv, info


def _dgbtrs(ab, ipiv, b):
    """Reference DGBTRS 'No transpose' on the DGBTF2 factorization.
    b: (n, nrhs)."""
    n, nrhs = b.shape

    # Solve L*X = B (row interchanges + rank-1 updates); the
    # multipliers of column j live in AB rows KV+1..KV+KL (0-based).
    bp = jnp.pad(b, ((0, KL), (0, 0)))

    def lbody(j, bp):
        block = lax.dynamic_slice(bp, (j, 0), (KL + 1, nrhs))
        l_off = ipiv[j] - j                      # 0..KL
        row0 = block[0]
        rowl = block[l_off]
        need = l_off != 0
        block = block.at[0].set(jnp.where(need, rowl, row0))
        block = block.at[l_off].set(jnp.where(need, row0, rowl))
        temp = -block[0]
        colmask = block[0] != 0.0                # DGER per-rhs skip
        mult = lax.dynamic_slice(ab, (KV + 1, j), (KL, 1))[:, 0]
        for i in range(1, KL + 1):
            mask = colmask & (j + i <= n - 1)
            block = block.at[i].set(
                jnp.where(mask, block[i] + mult[i - 1] * temp,
                          block[i]))
        return lax.dynamic_update_slice(bp, block, (j, 0))

    bp = lax.fori_loop(0, n - 1, lbody, bp)

    # Solve U*X = B (DTBSV upper, no transpose, non-unit, k = kl+ku);
    # leading pad rows absorb the (masked-off) below-bottom touches.
    bq = jnp.pad(bp[:n], ((KL + KU, 0), (0, 0)))

    def ubody(t, bq):
        j = n - 1 - t
        block = lax.dynamic_slice(bq, (j, 0), (KV + 1, nrhs))
        abcol = lax.dynamic_slice(ab, (0, j), (KV + 1, 1))[:, 0]
        xj = block[KV]
        mask = xj != 0.0                         # DTBSV zero-skip
        xj_new = xj / abcol[KV]
        block = block.at[KV].set(jnp.where(mask, xj_new, xj))
        temp = jnp.where(mask, xj_new, 0.0)
        for i in range(1, KV + 1):
            valid = mask & (j - i >= 0)
            block = block.at[KV - i].set(
                jnp.where(valid, block[KV - i] - temp * abcol[KV - i],
                          block[KV - i]))
        return lax.dynamic_update_slice(bq, block, (j, 0))

    bq = lax.fori_loop(0, n, ubody, bq)
    return bq[KL + KU:]


def band_solve(lhs, rhs):
    """lapack_wrap band_solve (nsup = nsub = 2): CLUBB band storage
    lhs (5, n), rhs (n,) or (n, nrhs).  Returns (solution, singular).
    solution is undefined where singular=True (CLUBB never assigns it
    and raises a fatal error)."""
    lhs = jnp.asarray(lhs, dtype=jnp.float64)
    rhs = jnp.asarray(rhs, dtype=jnp.float64)
    squeeze = rhs.ndim == 1
    if squeeze:
        rhs = rhs[:, None]
    ab = _assemble_lulhs(lhs)
    ab, ipiv, info = _dgbtf2(ab)
    # DGBSV only back-substitutes when the factorization succeeded
    sol = _dgbtrs(ab, ipiv, rhs)
    singular = info > 0
    if squeeze:
        sol = sol[:, 0]
    return sol, singular
