"""Column energy diagnostics and the SHOC total-energy fixer.

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_energy_integrals_impl.hpp, shoc_energy_fixer_impl.hpp,
  shoc_update_host_dse_impl.hpp
"""

import jax
import jax.numpy as jnp

from ..foundation import constants as c
from ..foundation.thermo import exner_function
from . import constants as sc
from .interp import linear_interp


@jax.jit
def shoc_energy_integrals(host_dse, pdel, rtm, rcm, u_wind, v_wind):
    """Mass-weighted column integrals (Functions::shoc_energy_integrals).

    Returns (se_int, ke_int, wv_int, wl_int), each of shape (...,):
    static energy, kinetic energy, water vapor, and liquid water integrals,
    all weighted by pdel/g.
    """
    w = jnp.asarray(pdel) / c.gravit
    se_int = jnp.sum(jnp.asarray(host_dse) * w, axis=-1)
    ke_int = jnp.sum(0.5 * (jnp.asarray(u_wind) ** 2 + jnp.asarray(v_wind) ** 2) * w,
                     axis=-1)
    wv_int = jnp.sum((jnp.asarray(rtm) - jnp.asarray(rcm)) * w, axis=-1)
    wl_int = jnp.sum(jnp.asarray(rcm) * w, axis=-1)
    return se_int, ke_int, wv_int, wl_int


@jax.jit
def update_host_dse(thlm, shoc_ql, inv_exner, zt_grid, phis):
    """Host-model dry static energy from SHOC prognostics
    (Functions::update_host_dse):

    T = thlm/inv_exner + (Lv/cp) ql;  dse = cp T + g z + phis.
    phis broadcasts over columns (shape (...,) or scalar).
    """
    temp = (jnp.asarray(thlm) / jnp.asarray(inv_exner)
            + (c.LatVap / c.CP) * jnp.asarray(shoc_ql))
    phis = jnp.asarray(phis)
    if phis.ndim:
        phis = phis[..., None]
    return c.CP * temp + c.gravit * jnp.asarray(zt_grid) + phis


@jax.jit
def shoc_energy_fixer(dtime, nadv, zt_grid, zi_grid,
                      se_b, ke_b, wv_b, wl_b,
                      se_a, ke_a, wv_a, wl_a,
                      wthl_sfc, wqw_sfc, rho_zt, tke, pint, host_dse):
    """Total-energy conservation fixer (Functions::shoc_energy_fixer).

    Compares before/after energy integrals (plus surface fluxes over the
    host step dtime*nadv), spreads the disbalance uniformly (per unit
    pressure) over the SHOC-active depth, and subtracts it from host_dse.
    Returns the updated host_dse.

    The *_b/*_a scalars are the shoc_energy_integrals outputs before/after
    the SHOC update; tke is post-update; pint is on interfaces.
    """
    zt_grid = jnp.asarray(zt_grid)
    zi_grid = jnp.asarray(zi_grid)
    rho_zt = jnp.asarray(rho_zt)
    tke = jnp.asarray(tke)
    pint = jnp.asarray(pint)
    host_dse = jnp.asarray(host_dse)

    nlev = zt_grid.shape[-1]
    nlevi = nlev + 1

    # Interface density (recomputed as in the C++, floor at 0)
    rho_zi = linear_interp(zt_grid, zi_grid, rho_zt, nlev, nlevi, 0.0)

    hdtime = dtime * nadv
    exner_int = exner_function(pint[..., -1])

    shf = wthl_sfc * c.CP * rho_zi[..., -1] * exner_int
    lhf = wqw_sfc * rho_zi[..., -1]
    te_a = se_a + ke_a + (c.LatVap + c.LatIce) * wv_a + c.LatIce * wl_a
    te_b = (se_b + ke_b + (c.LatVap + c.LatIce) * wv_b + c.LatIce * wl_b
            + (shf + lhf * (c.LatVap + c.LatIce)) * hdtime)

    # Highest level where SHOC is active: first k (from the model top) with
    # tke != mintke, capped at nlev-2 (mirrors the C++ masked min-reduction).
    k_idx = jnp.arange(nlev)
    cand = jnp.where(tke != sc.mintke, k_idx, nlev - 2)
    shoctop = jnp.minimum(jnp.min(cand, axis=-1), nlev - 2)

    p_top = jnp.take_along_axis(pint, shoctop[..., None], axis=-1)[..., 0]
    se_dis = (te_a - te_b) / (pint[..., -1] - p_top)

    active = k_idx >= shoctop[..., None]
    return jnp.where(active, host_dse - se_dis[..., None] * c.gravit, host_dse)
