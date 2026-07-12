"""Turbulent mixing-length chain.

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_compute_brunt_shoc_length_impl.hpp,
  shoc_compute_l_inf_shoc_length_impl.hpp,
  shoc_compute_shoc_mix_shoc_length_impl.hpp,
  shoc_check_length_scale_shoc_length_impl.hpp,
  shoc_length_impl.hpp (driver)
"""

import jax
import jax.numpy as jnp

from ..foundation import constants as c
from . import constants as sc
from .interp import linear_interp


@jax.jit
def compute_brunt_shoc_length(dz_zt, thv, thv_zi):
    """Brunt-Vaisala frequency (squared) on midpoints:
    brunt(k) = (g/thv(k)) * (thv_zi(k) - thv_zi(k+1)) / dz_zt(k).
    thv_zi has nlev+1 entries; result nlev."""
    thv = jnp.asarray(thv)
    thv_zi = jnp.asarray(thv_zi)
    return (c.gravit / thv) * (thv_zi[..., :-1] - thv_zi[..., 1:]) / jnp.asarray(dz_zt)


@jax.jit
def compute_l_inf_shoc_length(zt_grid, dz_zt, tke):
    """Asymptotic length scale: 0.1 * tke^1/2-weighted mean height."""
    w = jnp.sqrt(jnp.asarray(tke)) * jnp.asarray(dz_zt)
    numer = jnp.sum(w * jnp.asarray(zt_grid), axis=-1)
    denom = jnp.sum(w, axis=-1)
    return 0.1 * (numer / denom)


@jax.jit
def compute_shoc_mix_shoc_length(length_fac, tke, brunt, zt_grid, l_inf):
    """SHOC mixing length (Functions::compute_shoc_mix_shoc_length).

    tscale = 400 s eddy turnover; blends near-wall, asymptotic, and
    stability-limited scales; capped at maxlen. l_inf broadcasts over
    columns (shape (...,)).
    """
    tke = jnp.asarray(tke)
    zt_grid = jnp.asarray(zt_grid)
    l_inf = jnp.asarray(l_inf)
    if l_inf.ndim:
        l_inf = l_inf[..., None]

    tscale = 400.0
    tkes = jnp.sqrt(tke)
    brunt2 = jnp.maximum(0.0, jnp.asarray(brunt))

    inv_sum = (1.0 / (tscale * tkes * c.Karman * zt_grid)
               + 1.0 / (tscale * tkes * l_inf)
               + 0.01 * (brunt2 / tke))
    return jnp.minimum(sc.maxlen, 2.8284 * jnp.sqrt(1.0 / inv_sum) / length_fac)


@jax.jit
def check_length_scale_shoc_length(dx, dy, shoc_mix):
    """Clip mixing length to [minlen, sqrt(dx*dy)]; dx/dy broadcast over
    columns."""
    dx = jnp.asarray(dx)
    dy = jnp.asarray(dy)
    upper = jnp.sqrt(dx * dy)
    if upper.ndim:
        upper = upper[..., None]
    return jnp.minimum(upper, jnp.maximum(sc.minlen, jnp.asarray(shoc_mix)))


@jax.jit
def shoc_length(length_fac, dx, dy, zt_grid, zi_grid, dz_zt, tke, thv):
    """Mixing-length driver (Functions::shoc_length).

    Returns (brunt, shoc_mix), both on midpoints.
    """
    zt_grid = jnp.asarray(zt_grid)
    nlev = zt_grid.shape[-1]

    thv_zi = linear_interp(zt_grid, jnp.asarray(zi_grid), jnp.asarray(thv),
                           nlev, nlev + 1, 0.0)
    brunt = compute_brunt_shoc_length(dz_zt, thv, thv_zi)
    l_inf = compute_l_inf_shoc_length(zt_grid, dz_zt, tke)
    shoc_mix = compute_shoc_mix_shoc_length(length_fac, tke, brunt, zt_grid, l_inf)
    shoc_mix = check_length_scale_shoc_length(dx, dy, shoc_mix)
    return brunt, shoc_mix
