"""Turbulent mountain stress kernel (tms::Functions::compute_tms).

Source: components/eamxx/src/physics/tms/impl/compute_tms_impl.hpp
(single kernel; Fortran lineage components/eam/src/physics/cam/trb_mtn_stress.F90).

Computes a surface drag coefficient and wind stress from subgrid orography,
using only the bottom two model levels (k = nlev-2, nlev-1; k=0 is the model
top so these are the layers nearest the surface).

NOTE (transcription fidelity): the Richardson number uses t_mid*exner
products exactly as in the C++. What the `exner` argument contains is
defined by the process interface that calls this kernel — mirror it when
assembling, do not reinterpret.
"""

import jax
import jax.numpy as jnp

from ..foundation import constants as c

# Kernel-local constants (values from compute_tms_impl.hpp)
_HOROMIN = 1.0    # minimum subgrid orographic height for mountain stress [m]
_Z0MAX = 100.0    # maximum orographic z_0 [m]
_DV2MIN = 0.01    # minimum shear squared [m2/s2]


@jax.jit
def compute_tms(u_wind, v_wind, t_mid, p_mid, exner, z_mid, sgh, landfrac):
    """Compute the TMS drag coefficient and surface stress.

    Args:
        u_wind, v_wind: (ncol, nlev) horizontal wind components [m/s].
        t_mid, p_mid, exner, z_mid: (ncol, nlev) midpoint temperature [K],
            pressure [Pa], exner values as provided by the process interface,
            and geometric height above the surface [m].
        sgh: (ncol,) standard deviation of subgrid orography [m].
        landfrac: (ncol,) land fraction [-].

    Returns:
        ksrf: (ncol,) surface drag coefficient [kg/s/m2]
        tau_tms: (ncol, 2) surface stress (x, y components) [N/m2]
    """
    horo = c.orocnst * jnp.asarray(sgh)
    active = horo >= _HOROMIN

    # Guard inactive columns so the log/divisions below stay finite there
    # (their results are discarded by the where at the end).
    z0oro = jnp.where(active, jnp.minimum(c.z0fac * horo, _Z0MAX), 1.0)

    # Neutral drag coefficient from the lowest midpoint height
    tmp = c.Karman / jnp.log((z_mid[..., -1] + z0oro) / z0oro)
    cd = tmp * tmp

    # Richardson number over the lowest two layers (kb = nlev-1, kt = nlev-2)
    theta_kt = t_mid[..., -2] * exner[..., -2]
    theta_kb = t_mid[..., -1] * exner[..., -1]
    du = u_wind[..., -2] - u_wind[..., -1]
    dv = v_wind[..., -2] - v_wind[..., -1]
    dv2 = jnp.maximum(du * du + dv * dv, _DV2MIN)
    ri = (2.0 * c.gravit * (theta_kt - theta_kb) * (z_mid[..., -2] - z_mid[..., -1])
          / ((theta_kt + theta_kb) * dv2))

    # Crude stability function: 1 for ri<0, 0 for ri>1, linear ramp between.
    stabfri = jnp.maximum(0.0, jnp.minimum(1.0, 1.0 - ri))
    cd = cd * stabfri

    # Stress from bottom-level properties
    rho = p_mid[..., -1] / (c.Rair * t_mid[..., -1])
    vmag = jnp.sqrt(u_wind[..., -1] ** 2 + v_wind[..., -1] ** 2)

    ksrf = jnp.where(active, rho * cd * vmag * jnp.asarray(landfrac), 0.0)
    tau_x = -ksrf * u_wind[..., -1]
    tau_y = -ksrf * v_wind[..., -1]
    return ksrf, jnp.stack([tau_x, tau_y], axis=-1)
