"""CLUBB saturation — port of saturation.F90, "flatau" formula only.

PORT_NOTES
----------
Source: components/eam/src/physics/clubb/saturation.F90.  EAM sets
saturation_equation = "flatau" (clubb_intr.F90); the bolton/gfdl/
lookup variants are not ported.

sat_vapor_press_liq_flatau: the FACTORED form of the 8th-order Flatau
polynomial (clubb issue 834 rewrite) with its verbatim root constants;
T clipped below at -85 C.

sat_vapor_press_ice_flatau: the Horner form with the Table-4 ice
coefficients (each multiplied by 100 at compile time); T clipped below
at -90 C.

sat_mixrat_*: rs = ep*esat/(p-esat), with the esat ~ p guard
(p - esat < 1 Pa -> rs = ep).  NB the liquid version is written
`ep * esatv / (p - esatv)` (left-to-right: (ep*esatv)/(p-esatv))
while the ice version is `ep * ( esat_ice / (p - esat_ice) )` —
different rounding order, reproduced verbatim.
"""

import jax.numpy as jnp

from .constants import EP, T_FREEZE_K

_ICE_A = tuple(100.0 * c for c in
               (6.09868993, 0.499320233, 0.184672631e-01,
                0.402737184e-03, 0.565392987e-05, 0.521693933e-07,
                0.307839583e-09, 0.105785160e-11, 0.161444444e-14))


def sat_vapor_press_liq_flatau(t_in_k):
    t_in_c = jnp.maximum(t_in_k - T_FREEZE_K, -85.0)
    t_sqd = t_in_c ** 2
    return (-3.21582393e-14
            * (t_in_c - 646.5835252598777)
            * (t_in_c + 90.72381630364440)
            * (t_sqd + 111.0976961559954 * t_in_c + 6459.629194243118)
            * (t_sqd + 152.3131930092453 * t_in_c + 6499.774954705265)
            * (t_sqd + 174.4279584934021 * t_in_c + 7721.679732114084))


def sat_vapor_press_ice_flatau(t_in_k):
    t = jnp.maximum(t_in_k - T_FREEZE_K, -90.0)
    a = _ICE_A
    return a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * (
        a[4] + t * (a[5] + t * (a[6] + t * (a[7] + t * a[8])))))))


def sat_mixrat_liq(p_in_pa, t_in_k):
    esatv = sat_vapor_press_liq_flatau(t_in_k)
    return jnp.where(p_in_pa - esatv < 1.0, EP,
                     (EP * esatv) / (p_in_pa - esatv))


def sat_mixrat_ice(p_in_pa, t_in_k):
    esati = sat_vapor_press_ice_flatau(t_in_k)
    return jnp.where(p_in_pa - esati < 1.0, EP,
                     EP * (esati / (p_in_pa - esati)))
