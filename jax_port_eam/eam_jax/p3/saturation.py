"""Saturation vapor pressure / mixing ratio for the EAM P3 port.

Source: components/eam/src/physics/cam/wv_sat_scream.F90. EAM P3 calls
qv_sat (the "legacy" form): MurphyKoop svp with the wet denominator
    qv_sat = ep_2 * es / max(1e-3, p - es)
(scream_jax uses qv_sat_dry with denominator p — the sole adaptation
from the copied foundation/saturation.py, whose murphy_koop_svp is
retained verbatim; polysvp1 is dead code in EAM P3 and not ported).
"""

import jax.numpy as jnp

from . import constants as c


def murphy_koop_svp(t_atm, ice: bool):
    """Saturation vapor pressure [Pa], Murphy & Koop (2005)
    (wv_sat_scream.F90 MurphyKoop_svp). With ice=True the ice formula
    applies where T < T_zerodegc, the liquid formula elsewhere."""
    t = jnp.asarray(t_atm)
    logt = jnp.log(t)

    ic = (9.550426, 5723.265, 3.53068, 0.00728332)
    ice_result = jnp.exp(ic[0] - (ic[1] / t) + (ic[2] * logt) - (ic[3] * t))

    lq = (54.842763, 6763.22, 4.210, 0.000367, 0.0415, 218.8, 53.878,
          1331.22, 9.44523, 0.014025)
    liq_result = jnp.exp(
        lq[0] - (lq[1] / t) - (lq[2] * logt) + (lq[3] * t)
        + (jnp.tanh(lq[4] * (t - lq[5]))
           * (lq[6] - (lq[7] / t) - (lq[8] * logt) + lq[9] * t))
    )

    if not ice:
        return liq_result
    return jnp.where(t < c.T_zerodegc, ice_result, liq_result)


def qv_sat(t_atm, p_atm, ice: bool):
    """EAM qv_sat: ep_2*es/max(1e-3, p - es) (wv_sat_scream.F90)."""
    e_pres = murphy_koop_svp(t_atm, ice)
    return c.ep_2 * e_pres / jnp.maximum(jnp.asarray(p_atm) - e_pres, 1.0e-3)
