"""EAM MCICA subcolumn generator: KISS RNG + maximum-random overlap.

FRESH PORT (the SCREAM/EAMxx MCICA is a different algorithm: JSF64
PRNG seeded per (col,lay,gpt); EAM uses the vectorized KISS generator
with a per-column stream). Sources, ported verbatim:

  share/RandNum/src/kissvec/kissvec.c  (kiss_rng: Marsaglia KISS —
    congruential + 3-shift register + two 16-bit multiply-with-carry;
    ran = (int32)(s1+s2+(s3<<16)+s4) * 2.328306e-10 + 0.5)
  share/RandNum/src/shr_RandNum_mod.F90 (ShrKissRandGen: random(array)
    draws array(:,i) column-vectors sequentially, one kiss step each)
  components/eam/src/physics/rrtmgp/mcica_subcol_gen.F90
    (mcica_subcol_mask): seeds from the fractional parts of the
    bottom-four layer pressures, `changeseed` warm-up draws, ngpt
    draws of (ncol, pver) random fields, top-down maximum-random
    overlap, iscloudy = cdf >= 1 - cldf.
  components/eam/src/physics/rrtmgp/cam_optics.F90
    (sample_cloud_optics_sw/lw): combined cloud+snow fraction,
    changeseed = 1 for BOTH SW and LW (the streams differ only through
    ngpt), band->gpt subsampling with tau=0/ssa=1/asm=0 outside the
    cloudy mask.

The KISS state is exactly reproduced with uint32 arithmetic
(jax_enable_x64 required for the final int32 -> float64 conversion).
"""

import jax.numpy as jnp
from jax import lax

_KISS_SCALE = 2.328306e-10  # C double literal in kissvec.c


def _kiss_step(state):
    """One kiss_rng step for a vector of states. state: 4 uint32
    arrays. Returns (new_state, ran) with ran float64 in (0,1)."""
    s1, s2, s3, s4 = state
    s1 = jnp.uint32(69069) * s1 + jnp.uint32(1327217885)
    s2 = s2 ^ (s2 << jnp.uint32(13))
    s2 = s2 ^ (s2 >> jnp.uint32(17))
    s2 = s2 ^ (s2 << jnp.uint32(5))
    s3 = jnp.uint32(18000) * (s3 & jnp.uint32(65535)) + (s3 >> jnp.uint32(16))
    s4 = jnp.uint32(30903) * (s4 & jnp.uint32(65535)) + (s4 >> jnp.uint32(16))
    total = s1 + s2 + (s3 << jnp.uint32(16)) + s4
    ran = total.astype(jnp.int32).astype(jnp.float64) * _KISS_SCALE + 0.5
    return (s1, s2, s3, s4), ran


def kiss_seeds_from_pmid(pmid):
    """mcica_subcol_mask seeding: kiss_seed(i,j) from the fractional
    part of pmid at the bottom four layers (j=1 -> pver, j=4 ->
    pver-3), truncated toward zero after *1e9."""
    pmid = jnp.asarray(pmid)
    seeds = []
    for j in range(4):
        p = pmid[:, -(j + 1)]
        seeds.append(jnp.trunc((p - jnp.trunc(p)) * 1000000000.0)
                     .astype(jnp.int64).astype(jnp.uint32))
    return tuple(seeds)


def mcica_subcol_mask(ngpt, pmid, cldfrac, changeseed=1):
    """mcica_subcol_gen.F90 mcica_subcol_mask. pmid/cldfrac
    (ncol, pver), top-to-surface. Returns a bool mask
    (ngpt, ncol, pver) matching the Fortran iscloudy."""
    pmid = jnp.asarray(pmid)
    cldf = jnp.asarray(cldfrac)
    ncol, pver = cldf.shape
    cldf = jnp.where(cldf < 1.0e-80, 0.0, cldf)

    state = kiss_seeds_from_pmid(pmid)
    # changeseed warm-up draws (each is one kiss step per column)
    for _ in range(changeseed):
        state, _ = _kiss_step(state)

    # ngpt fields of (ncol, pver), drawn level-by-level: draw order is
    # (isubcol, ilev) with one vectorized step over columns each
    def step(state, _):
        state, ran = _kiss_step(state)
        return state, ran

    _, draws = lax.scan(step, state, None, length=ngpt * pver)
    # draws: (ngpt*pver, ncol) -> cdf (ngpt, ncol, pver)
    cdf = jnp.moveaxis(draws.reshape(ngpt, pver, ncol), 2, 1)

    # maximum-random overlap, marching top -> surface
    def overlap(prev, args):
        cdf_k, cldf_km1 = args           # (ngpt,ncol), (ncol,)
        cloudy_above = prev > 1.0 - cldf_km1[None, :]
        new = jnp.where(cloudy_above, prev, cdf_k * (1.0 - cldf_km1[None, :]))
        return new, new

    cdf0 = cdf[:, :, 0]
    _, rest = lax.scan(overlap, cdf0,
                       (jnp.moveaxis(cdf[:, :, 1:], 2, 0), cldf[:, :-1].T))
    cdf = jnp.concatenate([cdf0[:, :, None], jnp.moveaxis(rest, 0, 2)],
                          axis=2)

    return cdf >= 1.0 - cldf[None, :, :]


def sample_cloud_optics_sw(pmid, cld, cldfsnow, tau_bnd, ssa_bnd, asm_bnd,
                           gpt2bnd):
    """cam_optics sample_cloud_optics_sw (changeseed=1). Band arrays
    (ncol, nlev, nbnd) in RRTMGP band order; returns gpt arrays
    (ncol, nlev, ngpt)."""
    combined = jnp.maximum(jnp.asarray(cld), jnp.asarray(cldfsnow))
    gpt2bnd = jnp.asarray(gpt2bnd)
    ngpt = gpt2bnd.shape[0]
    iscloudy = mcica_subcol_mask(ngpt, pmid, combined, changeseed=1)
    keep = (jnp.moveaxis(iscloudy, 0, -1)          # (ncol, nlev, ngpt)
            & (combined[:, :, None] > 0.0))
    tau_gpt = jnp.where(keep, jnp.asarray(tau_bnd)[:, :, gpt2bnd], 0.0)
    ssa_gpt = jnp.where(keep, jnp.asarray(ssa_bnd)[:, :, gpt2bnd], 1.0)
    asm_gpt = jnp.where(keep, jnp.asarray(asm_bnd)[:, :, gpt2bnd], 0.0)
    return tau_gpt, ssa_gpt, asm_gpt


def sample_cloud_optics_lw(pmid, cld, cldfsnow, tau_bnd, gpt2bnd):
    """cam_optics sample_cloud_optics_lw (changeseed=1)."""
    combined = jnp.maximum(jnp.asarray(cld), jnp.asarray(cldfsnow))
    gpt2bnd = jnp.asarray(gpt2bnd)
    ngpt = gpt2bnd.shape[0]
    iscloudy = mcica_subcol_mask(ngpt, pmid, combined, changeseed=1)
    keep = (jnp.moveaxis(iscloudy, 0, -1)
            & (combined[:, :, None] > 0.0))
    return jnp.where(keep, jnp.asarray(tau_bnd)[:, :, gpt2bnd], 0.0)
