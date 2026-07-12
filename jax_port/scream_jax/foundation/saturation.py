"""Saturation vapor pressure and saturation mixing ratio.

Transcribed from
components/eamxx/src/share/physics/physics_saturation_impl.hpp
(declarations in physics_functions.hpp). The coefficient arrays are copied
digit-for-digit from the C++ — they are load-bearing; do not reformat them
against the original papers.

Differences from the C++, all deliberate:
- EKAT Pack/Mask vectorization is replaced by plain jnp arrays and
  jnp.where; the `range_mask` argument (which only masks pack padding)
  therefore has no equivalent and is dropped.
- check_temperature (an EKAT_KERNEL_REQUIRE abort on T<=0 or NaN) is not
  reproducible inside jit; use `check_temperature` below outside jit, or
  jax.debug callbacks, when validating inputs.
- The C++ wraps some literals in sp() (cast-to-float in single-precision
  builds). This port targets double precision only, where sp() is identity.
"""

import enum

import jax.numpy as jnp

from . import constants as c


class SaturationFcn(enum.IntEnum):
    """Which saturation vapor pressure formulation to use.

    Matches enum SaturationFcn in physics_functions.hpp: the runtime option
    plumbed from input.yaml. P3 currently defaults to Polysvp1 via its
    process interface; MurphyKoop is the qv_sat_* default argument.
    """

    POLYSVP1 = 0     # Flatau et al. 1992
    MURPHY_KOOP = 1  # Murphy & Koop 2005


def check_temperature(t_atm, caller: str = "check_temperature") -> None:
    """Host-side equivalent of Functions::check_temperature (abort on bad T).

    Only callable on concrete (non-traced) arrays — use in tests/adapters,
    not inside jitted code.
    """
    t = jnp.asarray(t_atm)
    if bool(jnp.any(t <= 0)):
        raise ValueError(f"{caller}: non-positive temperature encountered")
    if bool(jnp.any(jnp.isnan(t))):
        raise ValueError(f"{caller}: NaN temperature encountered")


def murphy_koop_svp(t_atm, ice: bool):
    """Saturation vapor pressure [Pa], Murphy & Koop (2005).

    Ice branch: eq. (7) (valid down to 110 K); liquid branch: eq. (10)
    (valid 123 K < T < 332 K). With ice=True the ice formula is used where
    T < Tmelt, the liquid formula elsewhere; ice=False is all-liquid
    (supersaturated with respect to liquid below freezing).
    """
    t = jnp.asarray(t_atm)
    logt = jnp.log(t)

    # Equation (7) — ice
    ic = (9.550426, 5723.265, 3.53068, 0.00728332)
    ice_result = jnp.exp(ic[0] - (ic[1] / t) + (ic[2] * logt) - (ic[3] * t))

    # Equation (10) — liquid
    lq = (54.842763, 6763.22, 4.210, 0.000367, 0.0415, 218.8, 53.878,
          1331.22, 9.44523, 0.014025)
    liq_result = jnp.exp(
        lq[0] - (lq[1] / t) - (lq[2] * logt) + (lq[3] * t)
        + (jnp.tanh(lq[4] * (t - lq[5]))
           * (lq[6] - (lq[7] / t) - (lq[8] * logt) + lq[9] * t))
    )

    if not ice:
        return liq_result
    return jnp.where(t < c.Tmelt, ice_result, liq_result)


def polysvp1(t_atm, ice: bool):
    """Saturation vapor pressure [Pa], Flatau et al. (1992), Table 4 RHS.

    Eighth-order polynomial in dt = max(T - 273.15, -80); coefficients give
    hPa, hence the trailing *100. Branch selection as in murphy_koop_svp.
    """
    t = jnp.asarray(t_atm)
    dt = jnp.maximum(t - 273.15, -80.0)

    # ice
    ai = (6.11147274, 0.503160820, 0.188439774e-1,
          0.420895665e-3, 0.615021634e-5, 0.602588177e-7,
          0.385852041e-9, 0.146898966e-11, 0.252751365e-14)
    ice_result = (ai[0] + dt*(ai[1] + dt*(ai[2] + dt*(ai[3] + dt*(ai[4]
                  + dt*(ai[5] + dt*(ai[6] + dt*(ai[7] + ai[8]*dt)))))))) * 100

    # liquid, V1.7
    a = (6.11239921, 0.443987641, 0.142986287e-1,
         0.264847430e-3, 0.302950461e-5, 0.206739458e-7,
         0.640689451e-10, -0.952447341e-13, -0.976195544e-15)
    liq_result = (a[0] + dt*(a[1] + dt*(a[2] + dt*(a[3] + dt*(a[4]
                  + dt*(a[5] + dt*(a[6] + dt*(a[7] + a[8]*dt)))))))) * 100

    if not ice:
        return liq_result
    return jnp.where(t < c.Tmelt, ice_result, liq_result)


def qv_sat_dry(t_atm, p_atm_dry, ice: bool,
               func: SaturationFcn = SaturationFcn.MURPHY_KOOP):
    """Saturation mixing ratio w.r.t. DRY air mass [kg/kg-dry].

    qsat_dry = ep_2 * es / max(p_dry, 1e-3); the pressure floor avoids
    division blow-up in padded/degenerate cells (as in the C++).
    """
    if func == SaturationFcn.POLYSVP1:
        e_pres = polysvp1(t_atm, ice)
    elif func == SaturationFcn.MURPHY_KOOP:
        e_pres = murphy_koop_svp(t_atm, ice)
    else:
        raise ValueError(f"Invalid saturation function: {func!r}")

    return c.ep_2 * e_pres / jnp.maximum(jnp.asarray(p_atm_dry), 1.0e-3)


def qv_sat_wet(t_atm, p_atm_dry, ice: bool, dp_wet, dp_dry,
               func: SaturationFcn = SaturationFcn.MURPHY_KOOP):
    """Saturation mixing ratio w.r.t. WET (total) air mass [kg/kg-wet].

    qsat_wet = qsat_dry * dp_dry / dp_wet, with dp_wet = pseudo_density and
    dp_dry = pseudo_density_dry.
    """
    qsatdry = qv_sat_dry(t_atm, p_atm_dry, ice, func)
    return qsatdry * jnp.asarray(dp_dry) / jnp.asarray(dp_wet)
