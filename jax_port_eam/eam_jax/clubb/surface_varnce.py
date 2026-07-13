"""CLUBB surface variances — port of surface_varnce_module.F90
(calc_surface_varnce) for the EAM configuration.

PORT_NOTES
----------
Source: components/eam/src/physics/clubb/surface_varnce_module.F90.

Diagnoses the surface (zm level 1) second-order moments from the
surface fluxes.  Only the l_andre_1978 = .false. branch is ported —
l_andre_1978 is a COMPILE-TIME parameter, so the Andre et al. (1978)
alternative (with its Monin-Obukhov zeta branches) is dead in every
EAM build; likewise sclr_dim = 0 removes the passive-scalar block
(iisclr_rt/iisclr_thl are -1 in EAM).

Active formulas (Andre et al. 1976/1978 velocity scales):
    ustar2 = sqrt(upwp^2 + vpwp^2)
    wstar  = (1/T0 * g * wpthlp * 1m)^(1/3) if wpthlp > 0 else 0
    uf     = max(ufmin, sqrt(ustar2 + 0.3*wstar^2)),  ufmin = 0.01
    wp2    = a * uf^2                       (a = a_const = 1.8)
    up2 = vp2 = up2_vp2_factor * a * uf^2   (tunable, default 2.0,
                                             no EAMv3 override)
    thlp2  = max(thl_tol^2, 0.4*a*(wpthlp/uf)^2)
    rtp2   = max(rt_tol^2, 0.4*a*(wprtp/uf)^2)
    rtpthlp = 0.2*a*(wpthlp/uf)*(wprtp/uf)
then the wp2 splatting correction with a correlation guard:
    min_wp2 = max(w_tol_sqd, wprtp^2/(rtp2*0.99^2),
                  wpthlp^2/(thlp2*0.99^2))
    if wp2 + tau_zm(1)*wp2_splat(1) < min_wp2:
        correction = -wp2 + min_wp2;  wp2 = min_wp2
    else:
        correction = tau_zm(1)*wp2_splat(1);  wp2 += correction
    up2 -= 0.5*correction;  vp2 -= 0.5*correction

The um_sfc/vm_sfc/Lscale_up_sfc arguments are only read by the dead
Andre-1978 branch but remain part of the interface (advance_clubb_core
passes um(2), vm(2), Lscale_up(2)); they are accepted and ignored
here, mirroring the Fortran call.

advance_clubb_core only calls this when the lowest momentum level is
at the surface (|zm(1)-sfc_elevation| <= |zm(1)+sfc_elevation|*eps/2,
always true in EAM where zi_g(1)=0=sfc); the else branch (tolerance
values) lives in lscale_tau.py.
"""

import jax.numpy as jnp

from .constants import (GRAV, MAX_MAG_CORRELATION_FLUX, RT_TOL,
                        THL_TOL, W_TOL_SQD)

# max_mag_correlation_flux**2 — gfortran evaluates the integer power
# as a multiply; keep the identical double.
_MMCF_SQD = MAX_MAG_CORRELATION_FLUX * MAX_MAG_CORRELATION_FLUX

A_CONST = 1.8
Z_CONST = 1.0        # defined height of 1 m
UFMIN = 0.01         # minimum friction velocity [m/s]
ONE_THIRD = 1.0 / 3.0


def calc_surface_varnce(upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc,
                        um_sfc, vm_sfc, lscale_up_sfc, wp2_splat_sfc,
                        tau_zm_sfc, t0, up2_vp2_factor):
    """Returns (wp2, up2, vp2, thlp2, rtp2, rtpthlp) at the surface."""
    del um_sfc, vm_sfc, lscale_up_sfc   # dead-branch-only arguments

    ustar2 = jnp.sqrt(upwp_sfc * upwp_sfc + vpwp_sfc * vpwp_sfc)

    # abs() only guards the untaken (wpthlp <= 0) lane of the where
    # against NaN; it is the identity on the taken lane.
    wstar = jnp.where(
        wpthlp_sfc > 0.0,
        jnp.power(jnp.abs(1.0 / t0 * GRAV * wpthlp_sfc * Z_CONST),
                  ONE_THIRD),
        0.0)

    uf = jnp.sqrt(ustar2 + 0.3 * wstar * wstar)
    uf = jnp.maximum(UFMIN, uf)

    wp2_sfc = A_CONST * uf ** 2
    up2_sfc = up2_vp2_factor * A_CONST * uf ** 2
    vp2_sfc = up2_vp2_factor * A_CONST * uf ** 2

    thlp2_sfc = 0.4 * A_CONST * (wpthlp_sfc / uf) ** 2
    thlp2_sfc = jnp.maximum(THL_TOL ** 2, thlp2_sfc)

    rtp2_sfc = 0.4 * A_CONST * (wprtp_sfc / uf) ** 2
    rtp2_sfc = jnp.maximum(RT_TOL ** 2, rtp2_sfc)

    rtpthlp_sfc = 0.2 * A_CONST \
        * (wpthlp_sfc / uf) * (wprtp_sfc / uf)

    # splatting correction with the flux-correlation guard
    min_wp2_sfc_val = jnp.maximum(
        W_TOL_SQD,
        jnp.maximum(
            wprtp_sfc ** 2 / (rtp2_sfc * _MMCF_SQD),
            wpthlp_sfc ** 2 / (thlp2_sfc * _MMCF_SQD)))
    guard = wp2_sfc + tau_zm_sfc * wp2_splat_sfc < min_wp2_sfc_val
    wp2_splat_sfc_correction = jnp.where(
        guard, -wp2_sfc + min_wp2_sfc_val, tau_zm_sfc * wp2_splat_sfc)
    wp2_sfc = jnp.where(guard, min_wp2_sfc_val,
                        wp2_sfc + tau_zm_sfc * wp2_splat_sfc)
    up2_sfc = up2_sfc - 0.5 * wp2_splat_sfc_correction
    vp2_sfc = vp2_sfc - 0.5 * wp2_splat_sfc_correction

    return wp2_sfc, up2_sfc, vp2_sfc, thlp2_sfc, rtp2_sfc, rtpthlp_sfc
