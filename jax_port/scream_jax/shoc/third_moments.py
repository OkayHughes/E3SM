"""Third moment of vertical velocity (w'3), Canuto et al. closure.

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_compute_diag_third_shoc_moment_impl.hpp,
  shoc_clipping_diag_third_shoc_moments_impl.hpp,
  shoc_diag_third_shoc_moments_impl.hpp (driver)
"""

import functools

import jax
import jax.numpy as jnp

from ..foundation import constants as c
from . import constants as sc
from .interp import linear_interp


@functools.partial(jax.jit, static_argnames=("shoc_1p5tke",))
def compute_diag_third_shoc_moment(c_diag_3rd_mom, shoc_1p5tke: bool,
                                   w_sec, thl_sec, wthl_sec, tke,
                                   dz_zt, dz_zi, isotropy_zi, brunt_zi,
                                   w_sec_zi, thetal_zi):
    """w'3 on interfaces (Functions::compute_diag_third_shoc_moment).

    Midpoint inputs: w_sec, tke, dz_zt (nlev). Interface inputs: thl_sec,
    wthl_sec, dz_zi, isotropy_zi, brunt_zi, w_sec_zi, thetal_zi (nlev+1).
    Returns w3 on interfaces with w3(0) = w3(nlev) = 0.
    """
    w_sec = jnp.asarray(w_sec)
    thl_sec = jnp.asarray(thl_sec)
    wthl_sec = jnp.asarray(wthl_sec)
    tke = jnp.asarray(tke)
    dz_zt = jnp.asarray(dz_zt)
    dz_zi = jnp.asarray(dz_zi)
    isotropy_zi = jnp.asarray(isotropy_zi)
    brunt_zi = jnp.asarray(brunt_zi)
    w_sec_zi = jnp.asarray(w_sec_zi)
    thetal_zi = jnp.asarray(thetal_zi)

    zeros = jnp.zeros_like(thl_sec[..., :1])

    if shoc_1p5tke:
        return jnp.zeros_like(thl_sec)

    cd = c_diag_3rd_mom
    a0 = (0.52 * (1.0 / (cd * cd))) / (cd - 2.0)
    a1 = 0.87 / (cd * cd)
    a2 = 0.5 / cd
    a3 = 0.6 / (cd * (cd - 2.0))
    a4 = 2.4 / (3.0 * cd + 5.0)
    a5 = 0.6 / (cd * (3.0 + 5.0 * cd))

    # Interior interfaces k = 1..nlev-1. Slice shorthands:
    #   interface arrays: [k] -> [..., 1:-1], [k-1] -> [..., :-2], [k+1] -> [..., 2:]
    #   midpoint arrays:  [k] -> [..., 1:],   [k-1] -> [..., :-1]
    thedz = 1.0 / dz_zi[..., 1:-1]
    thedz2 = 1.0 / (dz_zt[..., 1:] + dz_zt[..., :-1])

    iso = isotropy_zi[..., 1:-1]
    isosqrd = iso ** 2
    buoy_sgs2 = isosqrd * brunt_zi[..., 1:-1]
    bet2 = c.gravit / thetal_zi[..., 1:-1]

    wthl_k = wthl_sec[..., 1:-1]
    wsz = w_sec_zi[..., 1:-1]

    thl_sec_diff = thl_sec[..., :-2] - thl_sec[..., 2:]
    wthl_sec_diff = wthl_sec[..., :-2] - wthl_sec[..., 2:]
    wsec_diff = w_sec[..., :-1] - w_sec[..., 1:]
    tke_diff = tke[..., :-1] - tke[..., 1:]

    f0 = thedz2 * bet2 ** 3 * (iso * iso) * (iso * iso) * wthl_k * thl_sec_diff
    f1 = (thedz2 * bet2 ** 2 * iso ** 3
          * (wthl_k * wthl_sec_diff + 0.5 * wsz * thl_sec_diff))
    f2 = (thedz * bet2 * isosqrd * wthl_k * wsec_diff
          + 2.0 * thedz2 * bet2 * isosqrd * wsz * wthl_sec_diff)
    f3 = (thedz2 * bet2 * isosqrd * wsz * wthl_sec_diff
          + thedz * bet2 * isosqrd * (wthl_k * tke_diff))
    f4 = thedz * iso * wsz * (wsec_diff + tke_diff)
    f5 = thedz * iso * wsz * wsec_diff

    omega0 = a4 / (1.0 - a5 * buoy_sgs2)
    omega1 = omega0 / (2.0 * cd)
    omega2 = omega1 * f3 + (5.0 / 4.0) * omega0 * f4

    x0 = (a2 * buoy_sgs2 * (1.0 - a3 * buoy_sgs2)) / (1.0 - (a1 + a3) * buoy_sgs2)
    y0 = (2.0 * a2 * buoy_sgs2 * x0) / (1.0 - a3 * buoy_sgs2)
    x1 = ((a0 * f0 + a1 * f1 + a2 * (1.0 - a3 * buoy_sgs2) * f2)
          / (1.0 - (a1 + a3) * buoy_sgs2))
    y1 = (2.0 * a2 * (buoy_sgs2 * x1 + (a0 / a1) * f0 + f1)) / (1.0 - a3 * buoy_sgs2)

    aa0 = omega0 * x0 + omega1 * y0
    aa1 = omega0 * x1 + omega1 * y1 + omega2

    interior = (aa1 - 1.2 * x1 - 1.5 * f5) / (cd - 1.2 * x0 + aa0)
    return jnp.concatenate([zeros, interior, zeros], axis=-1)


@jax.jit
def clipping_diag_third_shoc_moments(w_sec_zi, w3):
    """Clip unrealistically large |w3| to the (small, constant) default
    (Functions::clipping_diag_third_shoc_moments):

        where |w3(k)| > w3clip * sqrt(2 * w_sec_zi(k)^3):  w3(k) = 0.02
    """
    w3 = jnp.asarray(w3)
    w_sec_zi = jnp.asarray(w_sec_zi)
    w3clipdef = 0.02
    cond = sc.w3clip * jnp.sqrt(2.0 * w_sec_zi ** 3)
    return jnp.where(jnp.abs(w3) > cond, w3clipdef, w3)


@functools.partial(jax.jit, static_argnames=("shoc_1p5tke",))
def diag_third_shoc_moments(c_diag_3rd_mom, shoc_1p5tke: bool,
                            w_sec, thl_sec, wthl_sec, isotropy, brunt,
                            thetal, tke, dz_zt, dz_zi, zt_grid, zi_grid):
    """w'3 driver (Functions::diag_third_shoc_moments). Returns w3 on
    interfaces. Note the interpolation floors, transcribed from the C++:
    brunt floored at `largeneg` (i.e. effectively unfloored), w_sec at
    (2/3) mintke, isotropy and thetal at 0.
    """
    zt_grid = jnp.asarray(zt_grid)
    nlev = zt_grid.shape[-1]

    isotropy_zi = linear_interp(zt_grid, jnp.asarray(zi_grid),
                                jnp.asarray(isotropy), nlev, nlev + 1, 0.0)
    brunt_zi = linear_interp(zt_grid, jnp.asarray(zi_grid), jnp.asarray(brunt),
                             nlev, nlev + 1, sc.largeneg)
    w_sec_zi = linear_interp(zt_grid, jnp.asarray(zi_grid), jnp.asarray(w_sec),
                             nlev, nlev + 1, (2.0 / 3.0) * sc.mintke)
    thetal_zi = linear_interp(zt_grid, jnp.asarray(zi_grid), jnp.asarray(thetal),
                              nlev, nlev + 1, 0.0)

    w3 = compute_diag_third_shoc_moment(
        c_diag_3rd_mom, shoc_1p5tke, w_sec, thl_sec, wthl_sec, tke,
        dz_zt, dz_zi, isotropy_zi, brunt_zi, w_sec_zi, thetal_zi)
    return clipping_diag_third_shoc_moments(w_sec_zi, w3)
