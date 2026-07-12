"""Vertical grid quantities for SHOC.

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_grid_impl.hpp, shoc_dp_inverse_impl.hpp, shoc_compute_tmpi_impl.hpp
"""

import jax
import jax.numpy as jnp

from ..foundation import constants as c


@jax.jit
def shoc_grid(zt_grid, zi_grid, pdel):
    """Grid thicknesses and air density (Functions::shoc_grid).

    Args:
        zt_grid: (..., nlev) midpoint heights [m] (k=0 top).
        zi_grid: (..., nlev+1) interface heights [m]; zi(nlev) is the surface.
        pdel: (..., nlev) pressure thickness [Pa].

    Returns:
        dz_zt: (..., nlev) thermo-layer thickness zi(k) - zi(k+1).
        dz_zi: (..., nlev+1) interface-layer thickness:
            0 at k=0; zt(k-1) - zt(k) for 0<k<nlev; zt(nlev-1) at k=nlev.
        rho_zt: (..., nlev) air density pdel/(g*dz_zt) [kg/m3].
    """
    zt_grid = jnp.asarray(zt_grid)
    zi_grid = jnp.asarray(zi_grid)
    pdel = jnp.asarray(pdel)

    dz_zt = zi_grid[..., :-1] - zi_grid[..., 1:]

    zeros = jnp.zeros_like(zt_grid[..., :1])
    dz_zi = jnp.concatenate(
        [zeros,                                   # k = 0
         zt_grid[..., :-1] - zt_grid[..., 1:],    # 0 < k < nlev
         zt_grid[..., -1:]],                      # k = nlev (surface)
        axis=-1)

    rho_zt = (1.0 / c.gravit) * (pdel / dz_zt)
    return dz_zt, dz_zi, rho_zt


@jax.jit
def dp_inverse(rho_zt, dz_zt):
    """Inverse pressure thickness rdp_zt = 1/(g*rho*dz) (Functions::dp_inverse)."""
    return 1.0 / (c.gravit * jnp.asarray(rho_zt) * jnp.asarray(dz_zt))


@jax.jit
def compute_tmpi(dtime, rho_zi, dz_zi):
    """Implicit-solver coefficient on interfaces (Functions::compute_tmpi).

    tmpi(0) = 0; tmpi(k) = dtime * g * rho_zi(k) / dz_zi(k) for k > 0.
    All arrays are on the interface grid (..., nlev+1).
    """
    rho_zi = jnp.asarray(rho_zi)
    dz_zi = jnp.asarray(dz_zi)
    interior = dtime * (c.gravit * rho_zi[..., 1:]) / dz_zi[..., 1:]
    return jnp.concatenate([jnp.zeros_like(rho_zi[..., :1]), interior], axis=-1)
