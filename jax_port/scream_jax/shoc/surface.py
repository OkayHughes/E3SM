"""Surface-layer diagnostics.

Source: components/eamxx/src/physics/shoc/impl/shoc_diag_obklen_impl.hpp
"""

import jax
import jax.numpy as jnp

from ..foundation import constants as c
from . import constants as sc


@jax.jit
def shoc_diag_obklen(uw_sfc, vw_sfc, wthl_sfc, wqw_sfc,
                     thl_sfc, cldliq_sfc, qv_sfc):
    """Surface friction velocity, kinematic buoyancy flux, Obukhov length
    (Functions::shoc_diag_obklen). All arguments are per-column scalars
    (any broadcastable shape); returns (ustar, kbfs, obklen).
    """
    thl_sfc = jnp.asarray(thl_sfc)
    cldliq_sfc = jnp.asarray(cldliq_sfc)
    qv_sfc = jnp.asarray(qv_sfc)

    eps = c.ZVIR
    th_sfc = thl_sfc + (c.LatVap / c.CP) * cldliq_sfc
    thv_sfc = th_sfc * (1 + eps * qv_sfc - cldliq_sfc)

    ustar = jnp.maximum(sc.ustar_min,
                        jnp.sqrt(jnp.asarray(uw_sfc) ** 2 + jnp.asarray(vw_sfc) ** 2))
    kbfs = jnp.asarray(wthl_sfc) + eps * th_sfc * jnp.asarray(wqw_sfc)
    # Tiny signed offset avoids division by zero at neutral stability.
    sign_val = jnp.where(kbfs >= 0, 1e-10, -1e-10)
    obklen = -thv_sfc * ustar ** 3 / (c.gravit * c.Karman * (kbfs + sign_val))
    return ustar, kbfs, obklen
