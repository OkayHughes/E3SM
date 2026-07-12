"""TKE-related kernels.

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_check_tke_impl.hpp
  (shoc_tke, adv_sgs_tke, compute_shr_prod, integ_column_stability,
   isotropic_ts, eddy_diffusivities to follow in this module.)
"""

import jax
import jax.numpy as jnp

from . import constants as sc


@jax.jit
def check_tke(tke):
    """Clip TKE to its minimum allowed value (Functions::check_tke)."""
    return jnp.maximum(jnp.asarray(tke), sc.mintke)
