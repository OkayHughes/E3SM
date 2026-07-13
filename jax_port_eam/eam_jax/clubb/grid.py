"""CLUBB vertical grid — port of grid_class.F90 (EAM path).

PORT_NOTES
----------
Source: components/eam/src/physics/clubb/grid_class.F90.

EAM runs CLUBB with l_implemented=.true. and grid_type=3, and passes
BOTH the momentum heights (zi_g, index 1 = surface) and thermodynamic
heights (zt_g, with the below-surface ghost zt_g(1) = -zt_g(2)) every
timestep (clubb_intr.F90 -> setup_grid_heights_api).  On that path
setup_grid_heights simply copies both arrays (no halfway
reconstruction), then derives:

  dzm(k)       = zt(k+1) - zt(k)   (k<nz;  dzm(nz)   = dzm(nz-1))
  dzt(k)       = zm(k) - zm(k-1)   (k>1;   dzt(1)    = dzt(2))
  invrs_dzm(k) = 1/(zt(k+1)-zt(k)) (k<nz;  [nz] = [nz-1])
  invrs_dzt(k) = 1/(zm(k)-zm(k-1)) (k>1;   [1]  = [2])

NB invrs_dzm is a fresh division, NOT 1/dzm elementwise of the stored
dzm (same value bitwise anyway since both are 1/(same difference)).

Interpolation weights (linear; l_cubic_interp = .false., the
model_flags default EAM uses):
  weights_zt2zm(1,k) = factor_m(k), weights_zt2zm(2,k) = 1-factor_m(k)
    factor_m(k) = (zm(k)-zt(k))/(zt(k+1)-zt(k))          (k<nz)
    factor_m(nz) = (zm(nz)-zt(nz-1))/(zt(nz)-zt(nz-1))   (extension >1)
  weights_zm2zt(1,k) = factor_t(k), weights_zm2zt(2,k) = 1-factor_t(k)
    factor_t(k) = (zt(k)-zm(k-1))/(zm(k)-zm(k-1))        (k>1)
    factor_t(1) = (zt(1)-zm(1))/(zm(2)-zm(1))            (extension <0)

Operators (linear versions; all return arrays on the OTHER grid):
  zt2zm: azm(k) = w(1,k)*(azt(k+1)-azt(k)) + azt(k) for k<nz;
         azm(nz) = linear extension from azt(nz-1:nz).
         (linear_interp_factor form: factor*(hi-lo)+lo.)
  zm2zt: azt(k) = w(1,k)*(azm(k)-azm(k-1)) + azm(k-1) for k>1;
         azt(1) = linear extension from azm(1:2).
  ddzt : gradzt — derivative of a zt-field, on zm levels:
         out(k) = (azt(k+1)-azt(k))*invrs_dzm(k) for k<nz;
         out(nz) = (azt(nz)-azt(nz-1))*invrs_dzm(nz-1).
  ddzm : gradzm — derivative of a zm-field, on zt levels:
         out(k) = (azm(k)-azm(k-1))*invrs_dzt(k) for k>1;
         out(1) = (azm(2)-azm(1))*invrs_dzt(2).

Not ported (unused in the EAM configuration): grid_type 1/2 standalone
setup, cubic (Steffen 1990) interpolation, read_grid_heights.
"""

from typing import NamedTuple

import jax.numpy as jnp


class Grid(NamedTuple):
    """CLUBB grid arrays (all shape (nz,) except the weight factors)."""
    zm: jnp.ndarray
    zt: jnp.ndarray
    dzm: jnp.ndarray
    dzt: jnp.ndarray
    invrs_dzm: jnp.ndarray
    invrs_dzt: jnp.ndarray
    factor_zt2zm: jnp.ndarray  # weights_zt2zm(1,:); weights(2,:) = 1 - this
    factor_zm2zt: jnp.ndarray  # weights_zm2zt(1,:); weights(2,:) = 1 - this


def setup_grid(momentum_heights, thermodynamic_heights):
    """setup_grid_heights for l_implemented=.true. (EAM host path)."""
    zm = jnp.asarray(momentum_heights, dtype=jnp.float64)
    zt = jnp.asarray(thermodynamic_heights, dtype=jnp.float64)

    dzm = jnp.concatenate([zt[1:] - zt[:-1], (zt[-1] - zt[-2])[None]])
    dzt = jnp.concatenate([(zm[1] - zm[0])[None], zm[1:] - zm[:-1]])
    invrs_dzm = jnp.concatenate([1.0 / (zt[1:] - zt[:-1]),
                                 (1.0 / (zt[-1] - zt[-2]))[None]])
    invrs_dzt = jnp.concatenate([(1.0 / (zm[1] - zm[0]))[None],
                                 1.0 / (zm[1:] - zm[:-1])])

    factor_zt2zm = jnp.concatenate(
        [(zm[:-1] - zt[:-1]) / (zt[1:] - zt[:-1]),
         ((zm[-1] - zt[-2]) / (zt[-1] - zt[-2]))[None]])
    factor_zm2zt = jnp.concatenate(
        [((zt[0] - zm[0]) / (zm[1] - zm[0]))[None],
         (zt[1:] - zm[:-1]) / (zm[1:] - zm[:-1])])

    return Grid(zm, zt, dzm, dzt, invrs_dzm, invrs_dzt,
                factor_zt2zm, factor_zm2zt)


def zt2zm(gr: Grid, azt):
    """linear_interpolated_azm: zt-field -> zm levels."""
    azt = jnp.asarray(azt)
    interior = gr.factor_zt2zm[:-1] * (azt[1:] - azt[:-1]) + azt[:-1]
    top = ((azt[-1] - azt[-2]) / (gr.zt[-1] - gr.zt[-2])) \
        * (gr.zm[-1] - gr.zt[-1]) + azt[-1]
    return jnp.concatenate([interior, top[None]])


def zm2zt(gr: Grid, azm):
    """linear_interpolated_azt: zm-field -> zt levels."""
    azm = jnp.asarray(azm)
    interior = gr.factor_zm2zt[1:] * (azm[1:] - azm[:-1]) + azm[:-1]
    bottom = ((azm[1] - azm[0]) / (gr.zm[1] - gr.zm[0])) \
        * (gr.zt[0] - gr.zm[0]) + azm[0]
    return jnp.concatenate([bottom[None], interior])


def ddzt(gr: Grid, azt):
    """gradzt: vertical derivative of a zt-field, on zm levels."""
    azt = jnp.asarray(azt)
    interior = (azt[1:] - azt[:-1]) * gr.invrs_dzm[:-1]
    top = (azt[-1] - azt[-2]) * gr.invrs_dzm[-2]
    return jnp.concatenate([interior, top[None]])


def ddzm(gr: Grid, azm):
    """gradzm: vertical derivative of a zm-field, on zt levels."""
    azm = jnp.asarray(azm)
    interior = (azm[1:] - azm[:-1]) * gr.invrs_dzt[1:]
    bottom = (azm[1] - azm[0]) * gr.invrs_dzt[1]
    return jnp.concatenate([bottom[None], interior])
