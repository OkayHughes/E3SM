"""SHOC diagnostic thermodynamics.

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_compute_shoc_vapor_impl.hpp, shoc_compute_shoc_temperature_impl.hpp
"""

import jax
import jax.numpy as jnp

from ..foundation import constants as c


@jax.jit
def compute_shoc_vapor(qw, ql):
    """Water vapor from total water and cloud liquid: qv = qw - ql."""
    return jnp.asarray(qw) - jnp.asarray(ql)


@jax.jit
def compute_shoc_temperature(thetal, ql, inv_exner):
    """Absolute temperature from liquid-water potential temperature:

    tabs = thetal/inv_exner + (Lv/cp)*ql
    (Functions::compute_shoc_temperature; inv_exner = 1/exner as provided
    by the SHOC process interface.)
    """
    return (jnp.asarray(thetal) / jnp.asarray(inv_exner)
            + (c.LatVap / c.CP) * jnp.asarray(ql))
