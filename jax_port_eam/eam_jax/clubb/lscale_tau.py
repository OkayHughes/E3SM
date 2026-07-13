"""The advance_clubb_core Lscale/tau segment — the inline code between
the skewness assembly and advance_xm_wpxp, exactly as EAMv3 executes
it (slice D).

PORT_NOTES
----------
Source: components/eam/src/physics/clubb/advance_clubb_core_module.F90
(the inline segment at lines ~1006/1083-1516 under the EAMv3
configuration), plus set_Lscale_max (private one-liner,
l_implemented=T).  Validated against a verbatim Fortran transcription
of the same segment (harness drv_lscale_tau_segment) that is itself
validated BITWISE against the real advance_clubb_core through its
khzm/khzt diagnostics.

Sequence (EAMv3 flags; ipdf_call_placement=2 so the pre-advance pdf
block is skipped and the entry prognostic state feeds this segment
directly):
  Lscale_max = 0.25 * min(host_dx, host_dy)          [set_Lscale_max]
  wp2_zt   = max(zm2zt(wp2), w_tol_sqd)
  thvm     = thlm + ep1*thv_ds_zt*rtm
             + (Lv/(Cp*exner) - ep2*thv_ds_zt)*rcm
  em       = 0.5*(wp2 + vp2 + up2)                   [l_tke_aniso=T]
  sqrt_em_zt = sqrt(max(em_min, zm2zt(em)))          [em_min = 6e-4]
  Lscale   = compute_mixing_length(...)  — called ONCE
             (l_avg_Lscale is a compile-time .false. parameter: the
             perturbed-Lscale calls and 1/3-averaging are DEAD;
             l_diag_Lscale_from_tau=F so the invrs_tau branch is dead)
  tau_zt   = min(Lscale/sqrt_em_zt, taumax)          [taumax = 3600]
  tau_zm   = min(max(zt2zm(Lscale), 0)/sqrt(max(em_min, em)), taumax)
  Kh_zt    = c_K*Lscale*sqrt_em_zt;  Kh_zm likewise  [c_K = 0.2]
  wp2_splat/wp3_splat = term_wp2/wp3_splat           [C_wp2_splat = 0]
  surface variances:
    if |zm(1)-sfc_elevation| <= |zm(1)+sfc_elevation|*eps/2:
        calc_surface_varnce(fluxes, um(2), vm(2), Lscale_up(2),
                            wp2_splat(1), tau_zm(1))
    else: tolerance values (w_tol_sqd/thl_tol^2/rt_tol^2/0)
  stability_correction = calc_stability_correction(...)   [ACTIVE:
        l_stability_correct_tau_zm=T]
  tau_N2_zm = tau_zm / stability_correction   (= tau_C1_zm = tau_C6_zm)
  Cx_fnc_Richardson = 0.0   [l_use_C7_Richardson =
        l_use_C11_Richardson = l_use_wp3_pr3 = F: the
        compute_Cx_fnc_Richardson call is DEAD — documented in
        stability.py, not ported]

Not part of this segment (later slices): everything from
advance_xm_wpxp on.  wp2_zt is returned because advance_clubb_core
computes it at line ~1006 and slice C's advance_xp2_xpyp consumes it.

Tunables consumed (EAMv3 values recorded in the golden metadata): mu
(5e-4), c_K (0.2), taumax (3600), C_wp2_splat (0), lmin (derived,
lmin_coef*40 m = 4), lambda0_stability_coef (0.03), up2_vp2_factor
(2.0), T0 (theta0 = 300).
"""

from typing import NamedTuple

import jax.numpy as jnp

from .constants import EPS, THL_TOL, RT_TOL, W_TOL_SQD, CP, EP1, EP2, LV
from .grid import Grid, zm2zt, zt2zm
from .mixing_length import compute_mixing_length
from .stability import (calc_stability_correction, term_wp2_splat,
                        term_wp3_splat)
from .surface_varnce import calc_surface_varnce

EM_MIN = 1.5 * W_TOL_SQD    # constants_clubb: em >= (3/2)*w_tol_sqd


class LscaleTau(NamedTuple):
    em: jnp.ndarray
    thvm: jnp.ndarray
    sqrt_em_zt: jnp.ndarray
    lscale: jnp.ndarray
    lscale_up: jnp.ndarray
    lscale_down: jnp.ndarray
    tau_zt: jnp.ndarray
    tau_zm: jnp.ndarray
    kh_zt: jnp.ndarray
    kh_zm: jnp.ndarray
    wp2_splat: jnp.ndarray
    wp3_splat: jnp.ndarray
    stability_correction: jnp.ndarray
    tau_n2_zm: jnp.ndarray
    wp2_zt: jnp.ndarray
    cx_fnc_richardson: jnp.ndarray
    # surface (zm level 1) moments after calc_surface_varnce
    wp2_sfc: jnp.ndarray
    up2_sfc: jnp.ndarray
    vp2_sfc: jnp.ndarray
    thlp2_sfc: jnp.ndarray
    rtp2_sfc: jnp.ndarray
    rtpthlp_sfc: jnp.ndarray


def lscale_tau_segment(gr: Grid, dt, sfc_elevation, host_dx, host_dy,
                       wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
                       thlm, rtm, rcm, wp2, wp3, up2, vp2, um, vm,
                       p_in_pa, exner, thv_ds_zt, tunables):
    """Replays the advance_clubb_core Lscale/tau segment on an entry
    prognostic state.  `tunables` needs mu, c_K, taumax, C_wp2_splat,
    lmin, lambda0_stability_coef, up2_vp2_factor, T0."""
    mu = jnp.asarray(tunables["mu"], dtype=jnp.float64)
    c_k = tunables["c_K"]
    taumax = tunables["taumax"]
    c_wp2_splat = tunables["C_wp2_splat"]
    lmin = tunables["lmin"]
    lambda0 = tunables["lambda0_stability_coef"]
    t0 = tunables["T0"]

    lscale_max = 0.25 * jnp.minimum(host_dx, host_dy)

    wp2_zt = jnp.maximum(zm2zt(gr, wp2), W_TOL_SQD)

    thvm = thlm + EP1 * thv_ds_zt * rtm \
        + (LV / (CP * exner) - EP2 * thv_ds_zt) * rcm

    em = 0.5 * (wp2 + vp2 + up2)
    sqrt_em_zt = jnp.sqrt(jnp.maximum(EM_MIN, zm2zt(gr, em)))

    lscale, lscale_up, lscale_down = compute_mixing_length(
        gr, thvm, thlm, rtm, em, lscale_max, p_in_pa, exner,
        thv_ds_zt, mu, lmin)

    tau_zt = jnp.minimum(lscale / sqrt_em_zt, taumax)
    tau_zm = jnp.minimum(jnp.maximum(zt2zm(gr, lscale), 0.0)
                         / jnp.sqrt(jnp.maximum(EM_MIN, em)), taumax)

    kh_zt = c_k * lscale * sqrt_em_zt
    kh_zm = c_k * jnp.maximum(zt2zm(gr, lscale), 0.0) \
        * jnp.sqrt(jnp.maximum(em, EM_MIN))

    wp2_splat = term_wp2_splat(gr, c_wp2_splat, dt, wp2, wp2_zt,
                               tau_zm)
    wp3_splat = term_wp3_splat(gr, c_wp2_splat, dt, wp2, wp3, tau_zt)

    # surface-variance branch (verbatim comparison)
    at_sfc = jnp.abs(gr.zm[0] - sfc_elevation) \
        <= jnp.abs(gr.zm[0] + sfc_elevation) * EPS / 2
    sv = calc_surface_varnce(
        upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, um[1], vm[1],
        lscale_up[1], wp2_splat[0], tau_zm[0], t0,
        tunables["up2_vp2_factor"])
    off = (W_TOL_SQD, W_TOL_SQD, W_TOL_SQD, THL_TOL ** 2, RT_TOL ** 2,
           0.0)
    sfc_out = tuple(jnp.where(at_sfc, a, b) for a, b in zip(sv, off))

    stability_correction = calc_stability_correction(
        gr, thlm, lscale, em, exner, rtm, rcm, p_in_pa, thvm,
        lambda0, t0)
    # l_stability_correct_tau_zm = .true.
    tau_n2_zm = tau_zm / stability_correction

    cx_fnc_richardson = jnp.zeros_like(tau_zm)

    return LscaleTau(em, thvm, sqrt_em_zt, lscale, lscale_up,
                     lscale_down, tau_zt, tau_zm, kh_zt, kh_zm,
                     wp2_splat, wp3_splat, stability_correction,
                     tau_n2_zm, wp2_zt, cx_fnc_richardson, *sfc_out)
