"""Port of eam/src/physics/cam/dadadj.F90 (GFDL-style dry adiabatic
adjustment).

PORT_NOTES
----------
- Operates on the TOP `nlvdry` layer pairs (level 1 = model top in
  CAM index order); requires nlvdry < pver.
- A column is adjusted only if any pair is superadiabatic by more than
  zeps (dodad detection uses the same formula as the sweep).
- The sweep is order-dependent: pair k updates t[k+1] then t[k], and
  pair k+1 sees the updated t[k+1] — ported as a sequential
  lax.scan over pairs inside a lax.while_loop over sweeps.
- Convergence: up to 15 sweeps at zeps=2e-5; on failure zeps doubles
  (restart sweeps) until zeps > 1e-4, which the Fortran treats as a
  fatal error — the port returns `converged=False` for such columns
  instead of aborting (caller decides; the golden cases all converge).
- q is mass-averaged over the adjusted pair (enthalpy/moisture
  conserving within the pair).
"""

import jax
import jax.numpy as jnp

from .constants import CAPPA

NITER = 15
ZEPS0 = 2.0e-5
ZEPS_MAX = 1.0e-4


def _sweep(t, q, pmid, pint, pdel, c1, c2, c3, c4, zeps, nlvdry):
    """One adjustment sweep over pairs k = 0..nlvdry-1 (sequential).
    Arrays are per-column 1-D (lev,). Returns (t, q, any_adjusted)."""

    def body(carry, k):
        t, q, adjusted = carry
        zepsdp = zeps * (pmid[k + 1] - pmid[k])
        zgamma = c1[k] * (t[k] + t[k + 1])
        do = (t[k + 1] - t[k]) >= (zgamma + zepsdp)
        t_kp1 = t[k] * c3[k] + t[k + 1] * c4[k]
        t_k = c2[k] * t_kp1
        qave = ((pdel[k + 1] * q[k + 1] + pdel[k] * q[k])
                / (pdel[k + 1] + pdel[k]))
        t = t.at[k + 1].set(jnp.where(do, t_kp1, t[k + 1]))
        t = t.at[k].set(jnp.where(do, t_k, t[k]))
        q = q.at[k + 1].set(jnp.where(do, qave, q[k + 1]))
        q = q.at[k].set(jnp.where(do, qave, q[k]))
        return (t, q, adjusted | do), None

    (t, q, adjusted), _ = jax.lax.scan(
        body, (t, q, jnp.asarray(False)), jnp.arange(nlvdry))
    return t, q, adjusted


def _adjust_column(t, q, pmid, pint, pdel, nlvdry):
    """Full per-column adjustment with sweep iteration and zeps
    doubling. Returns (t, q, converged)."""
    k = jnp.arange(nlvdry)
    c1 = CAPPA * 0.5 * (pmid[k + 1] - pmid[k]) / pint[k + 1]
    c2 = (1.0 - c1) / (1.0 + c1)
    rdenom = 1.0 / (pdel[k] * c2 + pdel[k + 1])
    c3 = rdenom * pdel[k]
    c4 = rdenom * pdel[k + 1]

    def cond(state):
        t, q, zeps, sweeps_left, done, failed = state[:6]
        return (~done) & (~failed)

    def body(state):
        t, q, zeps, sweeps_left, done, failed = state
        t, q, adjusted = _sweep(t, q, pmid, pint, pdel,
                                c1, c2, c3, c4, zeps, nlvdry)
        done = ~adjusted
        sweeps_left = sweeps_left - 1
        # out of sweeps at this zeps: double and restart the counter
        out = (sweeps_left == 0) & ~done
        zeps = jnp.where(out, zeps + zeps, zeps)
        failed = out & (zeps > ZEPS_MAX)
        sweeps_left = jnp.where(out, NITER, sweeps_left)
        return (t, q, zeps, sweeps_left, done, failed)

    state = (t, q, jnp.asarray(ZEPS0), jnp.asarray(NITER),
             jnp.asarray(False), jnp.asarray(False))
    t, q, zeps, _, done, failed = jax.lax.while_loop(cond, body, state)
    return t, q, ~failed


def dadadj(pmid, pint, pdel, t, q, nlvdry=3):
    """Dry adiabatic adjustment of the top `nlvdry` layer pairs.

    pmid/pdel/t/q: (ncol, pver); pint: (ncol, pver+1). Returns
    (t, q, converged) with converged (ncol,) bool — the Fortran calls
    endrun where this is False."""
    t = jnp.asarray(t, dtype=jnp.float64)
    q = jnp.asarray(q, dtype=jnp.float64)
    pmid = jnp.asarray(pmid, dtype=jnp.float64)
    pint = jnp.asarray(pint, dtype=jnp.float64)
    pdel = jnp.asarray(pdel, dtype=jnp.float64)

    # dodad detection (pairs 0..nlvdry-1)
    k = jnp.arange(nlvdry)
    gammad = CAPPA * 0.5 * (t[:, k + 1] + t[:, k]) / pint[:, k + 1]
    dtdp = (t[:, k + 1] - t[:, k]) / (pmid[:, k + 1] - pmid[:, k])
    dodad = ((dtdp + ZEPS0) > gammad).any(axis=1)

    t_adj, q_adj, conv = jax.vmap(
        lambda tc, qc, pm, pi, pd: _adjust_column(tc, qc, pm, pi, pd,
                                                  nlvdry)
    )(t, q, pmid, pint, pdel)

    sel = dodad[:, None]
    return (jnp.where(sel, t_adj, t), jnp.where(sel, q_adj, q),
            jnp.where(dodad, conv, True))
