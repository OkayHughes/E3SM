"""CLUBB PDF closure — port of pdf_closure_module.F90 (ADG1 path) and
its helpers (adg1_adg2_3d_luhar_pdf.F90, pdf_utilities.F90,
Skx_module.F90, sigma_sqd_w_module.F90).

PORT_NOTES
----------
Scope is the exact configuration EAM compiles and runs (see the
package docstring): iiPDF_type = iiPDF_ADG1 (a compile-time constant),
sclr_dim = 0, hydromet_dim = 0, l_explicit_turbulent_adv_wp3 =
l_explicit_turbulent_adv_xpyp = .false., all stats indices 0
(l_stats=.false.), debug_level 0.  Consequences reproduced here:

- Only ADG1_pdf_driver runs; corr_w_rt/thl/chi/eta and corr_u_w/
  corr_v_w PDF-component correlations are identically zero.
- wp4, wprtp2, wpthlp2, wprtpthlp are NEVER computed (their compute
  blocks are gated on the explicit-adv flags / stats indices) — they
  are intent(out) garbage in EAM and are not returned by this port.
- wp2rtp and wp2thlp ARE always computed (calc_wp2xp_pdf with zero
  component correlations).
- rcp2 is always computed under #ifdef CLUBB_CAM.
- pdf_implicit_coefs_terms / F_w etc. are zero-filled for ADG1 and not
  returned.
- The GFDL rsatl blend and RH_crit paths compile out;
  l_calc_ice_supersat_frac = .true. always.
- The l_liq_ice_loading_test hydrometeor loading block is dead
  (local parameter .false. + hydromet_dim = 0).

Functions
---------
skx_func            Skx_module Skx_func (Skw_denom_coef = 0.0 under
                    CLUBB_CAM, but taken from params for generality;
                    l_clipping_kluge=.false. so no Skw_max_mag clip).
gamma_skw_fnc       The gamma_Skw_fnc expression from
                    advance_clubb_core (l_gamma_Skw = .true. branch).
compute_sigma_sqd_w sigma_sqd_w_module (l_predict_upwp_vpwp=.false.,
                    so only the thl and rt correlations enter; u/v
                    moments are accepted for signature parity).
pdf_closure         The main closure.  Column arrays (nz,); vmap for
                    batches.  Returns (moments_dict, pdf_params_dict).

Numerical-faithfulness notes (lessons from PORTING_PLAN row 9/10):
- Divisions by scalar constants go through _divs (XLA's reciprocal
  rewrite is not correctly rounded); array/array divisions are left
  as true divisions.
- erf comes from jax.scipy.special.erf vs gfortran's libm erf:
  cross-libm 1-ulp differences are possible on the partly-cloudy
  branch (measured at the 1e-15 relative level in the goldens' cf/rc).
- x**2/x**3/x**4 integer powers match gfortran's expansion by
  repeated multiplication.
- The `where (wp2 > w_tol_sqd)` / `where (xp2 > x_tol**2)` nested
  masking of the Fortran is reproduced with jnp.where trees whose
  guarded denominators are made safe before dividing.

Tunables used: beta, mixt_frac_max_mag (derived in
setup_parameters_model from Skw_max_mag with the hard-coded
sigma_sqd_w=0.4: mixt_frac_max_mag = 1 - 0.5*(1 - Skw_max_mag/
sqrt(4*(1-0.4)^3 + Skw_max_mag^2))), Skw_denom_coef, gamma_coef/b/c.
The EAMv3 values are recorded in the goldens' params vector.
"""

import math

import jax.numpy as jnp
from jax.scipy.special import erf

from .constants import (CHI_TOL, CP, EP, EP1, EP2, EPS, ETA_TOL, LV,
                        MAX_MAG_CORRELATION, MAX_NUM_STDEVS, RD, RT_TOL,
                        SQRT_2, SQRT_2PI, T_FREEZE_K, THL_TOL, W_TOL,
                        W_TOL_SQD, ZERO_THRESHOLD)
from .saturation import sat_mixrat_ice, sat_mixrat_liq

TWO_THIRDS = 2.0 / 3.0


def _divs(a, s):
    """Fortran-faithful division of an array by a scalar constant."""
    return a / jnp.full_like(a, s)


def mixt_frac_max_mag(skw_max_mag):
    """setup_parameters_model derivation (verbatim magic 0.4)."""
    return 1.0 - (0.5 * (1.0 - skw_max_mag /
                         math.sqrt(4.0 * (1.0 - 0.4) ** 3
                                   + skw_max_mag ** 2)))


def skx_func(xp2, xp3, x_tol, skw_denom_coef):
    """Skx_module Skx_func (l_clipping_kluge = .false.)."""
    skx_denom_tol = skw_denom_coef * x_tol ** 2
    return xp3 / ((xp2 + skx_denom_tol) * jnp.sqrt(xp2 + skx_denom_tol))


def gamma_skw_fnc(skw_zm, gamma_coef, gamma_coefb, gamma_coefc):
    """advance_clubb_core gamma_Skw_fnc (l_gamma_Skw = .true. and
    gamma_coef /= gamma_coefb branch; else it is just gamma_coef)."""
    if gamma_coef != gamma_coefb:
        return gamma_coefb + (gamma_coef - gamma_coefb) * jnp.exp(
            -(1.0 / 2.0) * (skw_zm / jnp.full_like(skw_zm, gamma_coefc)) ** 2)
    return jnp.full_like(skw_zm, gamma_coef)


def compute_sigma_sqd_w(gamma_skw_fnc_in, wp2, thlp2, rtp2,
                        up2, vp2, wpthlp, wprtp, upwp, vpwp):
    """sigma_sqd_w_module compute_sigma_sqd_w
    (l_predict_upwp_vpwp = .false.: u/v terms not included)."""
    del up2, vp2, upwp, vpwp  # l_predict_upwp_vpwp = .false.
    max_corr_w_x_sqd = jnp.maximum(
        (wpthlp / (jnp.sqrt(wp2 * thlp2) + 0.01 * W_TOL * THL_TOL)) ** 2,
        (wprtp / (jnp.sqrt(wp2 * rtp2) + 0.01 * W_TOL * RT_TOL)) ** 2)
    return gamma_skw_fnc_in * (1.0 - jnp.minimum(max_corr_w_x_sqd, 1.0))


def _adg1_w_closure(wm, wp2, skw, sigma_sqd_w, sqrt_wp2, mf_max_mag):
    """ADG1_w_closure (adg1_adg2_3d_luhar_pdf.F90)."""
    wide = wp2 > W_TOL_SQD

    mf_formula = 0.5 * (1.0 - skw / jnp.sqrt(
        4.0 * (1.0 - sigma_sqd_w) ** 3 + skw ** 2))
    mixt_frac_w = jnp.where(jnp.abs(skw) <= 1.0e-5, 0.5, mf_formula)
    mixt_frac_w = jnp.minimum(
        jnp.maximum(mixt_frac_w, 1.0 - mf_max_mag), mf_max_mag)

    w_1_n_w = jnp.sqrt(((1.0 - mixt_frac_w) / mixt_frac_w)
                       * (1.0 - sigma_sqd_w))
    w_2_n_w = -jnp.sqrt((mixt_frac_w / (1.0 - mixt_frac_w))
                        * (1.0 - sigma_sqd_w))

    mixt_frac = jnp.where(wide, mixt_frac_w, 0.5)
    w_1_n = jnp.where(wide, w_1_n_w, jnp.sqrt(1.0 - sigma_sqd_w))
    w_2_n = jnp.where(wide, w_2_n_w, -jnp.sqrt(1.0 - sigma_sqd_w))
    w_1 = jnp.where(wide, wm + sqrt_wp2 * w_1_n, wm)
    w_2 = jnp.where(wide, wm + sqrt_wp2 * w_2_n, wm)
    varnce_w = jnp.where(wide, sigma_sqd_w * wp2, 0.0)
    return w_1, w_2, w_1_n, w_2_n, varnce_w, varnce_w, mixt_frac


def _adg1_responder(xm, xp2, wp2, sqrt_wp2, wpxp, w_1_n, w_2_n,
                    mixt_frac, sigma_sqd_w, x_tol, beta):
    """ADG1_ADG2_responder_params."""
    w_wide = wp2 > W_TOL_SQD
    x_wide = xp2 > x_tol ** 2

    width_factor_1 = TWO_THIRDS * beta \
        + 2.0 * mixt_frac * (1.0 - TWO_THIRDS * beta)
    width_factor_2 = 2.0 - width_factor_1

    # guarded divisions for the (w_wide & x_wide) branch
    safe_w2n = jnp.where(w_wide, w_2_n, 1.0)
    safe_w1n = jnp.where(w_wide, w_1_n, 1.0)
    x_1_c = xm - (wpxp / sqrt_wp2) / safe_w2n
    x_2_c = xm - (wpxp / sqrt_wp2) / safe_w1n
    alpha_c = 0.5 * (1.0 - wpxp * wpxp
                     / ((1.0 - sigma_sqd_w) * wp2 * xp2))
    alpha_c = jnp.maximum(jnp.minimum(alpha_c, 1.0), ZERO_THRESHOLD)
    varnce_x_1_c = (alpha_c / mixt_frac * xp2) * width_factor_1
    varnce_x_2_c = (alpha_c / (1.0 - mixt_frac) * xp2) * width_factor_2

    both = w_wide & x_wide
    x_1 = jnp.where(both, x_1_c, xm)
    x_2 = jnp.where(both, x_2_c, xm)
    varnce_x_1 = jnp.where(both, varnce_x_1_c,
                           jnp.where(w_wide, 0.0, xp2))
    varnce_x_2 = jnp.where(both, varnce_x_2_c,
                           jnp.where(w_wide, 0.0, xp2))
    alpha_x = jnp.where(both, alpha_c, 0.5)
    return x_1, x_2, varnce_x_1, varnce_x_2, alpha_x


def _calc_comp_corrs_binormal(xpyp, xm, ym, mu_x_1, mu_x_2, mu_y_1,
                              mu_y_2, sig_x_1_sqd, sig_x_2_sqd,
                              sig_y_1_sqd, sig_y_2_sqd, mixt_frac):
    """pdf_utilities calc_comp_corrs_binormal (corr_1 == corr_2)."""
    defined = (sig_x_1_sqd * sig_y_1_sqd > 0.0) \
        | (sig_x_2_sqd * sig_y_2_sqd > 0.0)
    denom = mixt_frac * jnp.sqrt(sig_x_1_sqd * sig_y_1_sqd) \
        + (1.0 - mixt_frac) * jnp.sqrt(sig_x_2_sqd * sig_y_2_sqd)
    corr = (xpyp
            - mixt_frac * (mu_x_1 - xm) * (mu_y_1 - ym)
            - (1.0 - mixt_frac) * (mu_x_2 - xm) * (mu_y_2 - ym)) \
        / jnp.where(defined, denom, 1.0)
    corr = jnp.maximum(-MAX_MAG_CORRELATION,
                       jnp.minimum(MAX_MAG_CORRELATION, corr))
    corr = jnp.where(defined, corr, 0.0)
    return corr, corr


def _calc_wp2xp_pdf(wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, varnce_w_2,
                    varnce_x_1, varnce_x_2, corr_w_x_1, corr_w_x_2,
                    mixt_frac):
    """pdf_closure_module calc_wp2xp_pdf."""
    return (mixt_frac
            * (((w_1 - wm) ** 2 + varnce_w_1) * (x_1 - xm)
               + 2.0 * corr_w_x_1 * jnp.sqrt(varnce_w_1 * varnce_x_1)
               * (w_1 - wm))
            + (1.0 - mixt_frac)
            * (((w_2 - wm) ** 2 + varnce_w_2) * (x_2 - xm)
               + 2.0 * corr_w_x_2 * jnp.sqrt(varnce_w_2 * varnce_x_2)
               * (w_2 - wm)))


def _transform_pdf_chi_eta_component(tl, rsatl, rt, exner, varnce_thl,
                                     varnce_rt, corr_rt_thl):
    """transform_pdf_chi_eta_component."""
    beta_sd = EP * (LV / (RD * tl)) * (LV / (CP * tl))
    chi = (rt - rsatl) / (1.0 + beta_sd * rsatl)
    crt = 1.0 / (1.0 + beta_sd * rsatl)
    cthl = (1.0 + beta_sd * rt) / (1.0 + beta_sd * rsatl) ** 2 \
        * (CP / LV) * beta_sd * rsatl * exner

    varnce_rt_term = crt ** 2 * varnce_rt
    varnce_thl_term = cthl ** 2 * varnce_thl
    covar_chi_eta = varnce_rt_term - varnce_thl_term
    corr_rt_thl_term = 2.0 * corr_rt_thl * crt * cthl \
        * jnp.sqrt(varnce_rt * varnce_thl)
    stdev_chi = jnp.sqrt(varnce_rt_term - corr_rt_thl_term
                         + varnce_thl_term)
    stdev_eta = jnp.sqrt(varnce_rt_term + corr_rt_thl_term
                         + varnce_thl_term)

    degenerate = (stdev_chi < CHI_TOL) | (stdev_eta < ETA_TOL)
    stdev_chi = jnp.where(stdev_chi < CHI_TOL, 0.0, stdev_chi)
    stdev_eta = jnp.where(stdev_eta < ETA_TOL, 0.0, stdev_eta)
    corr_c = covar_chi_eta / jnp.where(degenerate, 1.0,
                                       stdev_chi * stdev_eta)
    corr_c = jnp.minimum(MAX_MAG_CORRELATION,
                         jnp.maximum(-MAX_MAG_CORRELATION, corr_c))
    corr_chi_eta = jnp.where(degenerate, 0.0, corr_c)
    return chi, crt, cthl, stdev_chi, stdev_eta, covar_chi_eta, \
        corr_chi_eta


def _calc_cloud_frac_component(mean_chi_i, stdev_chi_i, chi_at_sat):
    """calc_cloud_frac_component (elemental)."""
    dchi = mean_chi_i - chi_at_sat
    all_clear = ((jnp.abs(dchi) <= EPS) & (stdev_chi_i <= CHI_TOL)) \
        | (dchi < -MAX_NUM_STDEVS * stdev_chi_i)
    all_cloud = dchi > MAX_NUM_STDEVS * stdev_chi_i

    safe_stdev = jnp.where(all_clear | all_cloud, 1.0, stdev_chi_i)
    zeta_i = dchi / safe_stdev
    cf_partly = 0.5 * (1.0 + erf(_divs(zeta_i, SQRT_2)))
    rc_partly = dchi * cf_partly \
        + _divs(stdev_chi_i * jnp.exp(-0.5 * zeta_i ** 2), SQRT_2PI)

    cloud_frac_i = jnp.where(all_clear, 0.0,
                             jnp.where(all_cloud, 1.0, cf_partly))
    rc_i = jnp.where(all_clear, 0.0,
                     jnp.where(all_cloud, dchi, rc_partly))
    return cloud_frac_i, rc_i


def _calc_xprcp_component(wm, rtm, thlm, um, vm, rcm, w_i, rt_i, thl_i,
                          u_i, v_i, varnce_w_i, chi_i, stdev_chi_i,
                          stdev_eta_i, corr_chi_eta_i, crt_i, cthl_i,
                          rc_i, cloud_frac_i):
    """calc_xprcp_component (ADG1: the corr_w_chi_i where-block is
    skipped because iiPDF_type is ADG1)."""
    del chi_i  # only used in the non-ADG1 where-block
    wprcp_c = (w_i - wm) * (rc_i - rcm)
    wp2rcp_c = ((w_i - wm) ** 2 + varnce_w_i) * (rc_i - rcm)
    rtprcp_c = (rt_i - rtm) * (rc_i - rcm) \
        + (corr_chi_eta_i * stdev_eta_i + stdev_chi_i) \
        / (2.0 * crt_i) * stdev_chi_i * cloud_frac_i
    thlprcp_c = (thl_i - thlm) * (rc_i - rcm) \
        + (corr_chi_eta_i * stdev_eta_i - stdev_chi_i) \
        / (2.0 * cthl_i) * stdev_chi_i * cloud_frac_i
    uprcp_c = (u_i - um) * (rc_i - rcm)
    vprcp_c = (v_i - vm) * (rc_i - rcm)
    return wprcp_c, wp2rcp_c, rtprcp_c, thlprcp_c, uprcp_c, vprcp_c


def pdf_closure(p_in_pa, exner, thv_ds, wm, wp2, wp3, sigma_sqd_w,
                skw, skthl, skrt, rtm, rtp2, wprtp, thlm, thlp2,
                wpthlp, um, up2, upwp, vm, vp2, vpwp, rtpthlp,
                beta, mf_max_mag):
    """pdf_closure (iiPDF_ADG1, sclr_dim=0, hydromet_dim=0).

    All array arguments are (nz,) on the same level set (EAM calls it
    with zt-level fields and again with zm-interpolated fields).
    `beta` and `mf_max_mag` are the tunables (params[ibeta],
    mixt_frac_max_mag).  skthl/skrt are accepted for interface parity
    with the Fortran but are unused on the ADG1 path.

    Returns (moments, pdf_params) dicts. moments keys: wp2rtp,
    wp2thlp, cloud_frac, ice_supersat_frac, rcm, wpthvp, wp2thvp,
    rtpthvp, thlpthvp, wprcp, wp2rcp, rtprcp, thlprcp, rcp2, uprcp,
    vprcp, rc_coef.  pdf_params keys: the 47 pdf_parameter fields.
    """
    del wp3, skthl, skrt  # wp3 enters only through Skw for ADG1

    sqrt_wp2 = jnp.sqrt(wp2)

    # --- ADG1 PDF driver ---
    w_1, w_2, w_1_n, w_2_n, varnce_w_1, varnce_w_2, mixt_frac = \
        _adg1_w_closure(wm, wp2, skw, sigma_sqd_w, sqrt_wp2, mf_max_mag)

    rt_1, rt_2, varnce_rt_1, varnce_rt_2, alpha_rt = _adg1_responder(
        rtm, rtp2, wp2, sqrt_wp2, wprtp, w_1_n, w_2_n, mixt_frac,
        sigma_sqd_w, RT_TOL, beta)
    thl_1, thl_2, varnce_thl_1, varnce_thl_2, alpha_thl = \
        _adg1_responder(thlm, thlp2, wp2, sqrt_wp2, wpthlp, w_1_n,
                        w_2_n, mixt_frac, sigma_sqd_w, THL_TOL, beta)
    # NB the Fortran passes thl_tol for u and v as well.
    u_1, u_2, varnce_u_1, varnce_u_2, alpha_u = _adg1_responder(
        um, up2, wp2, sqrt_wp2, upwp, w_1_n, w_2_n, mixt_frac,
        sigma_sqd_w, THL_TOL, beta)
    v_1, v_2, varnce_v_1, varnce_v_2, alpha_v = _adg1_responder(
        vm, vp2, wp2, sqrt_wp2, vpwp, w_1_n, w_2_n, mixt_frac,
        sigma_sqd_w, THL_TOL, beta)
    del varnce_u_1, varnce_u_2, alpha_u, varnce_v_1, varnce_v_2, alpha_v

    # --- component correlations ---
    corr_rt_thl_1, corr_rt_thl_2 = _calc_comp_corrs_binormal(
        rtpthlp, rtm, thlm, rt_1, rt_2, thl_1, thl_2, varnce_rt_1,
        varnce_rt_2, varnce_thl_1, varnce_thl_2, mixt_frac)

    zero = jnp.zeros_like(wm)
    corr_w_rt_1 = corr_w_rt_2 = zero
    corr_w_thl_1 = corr_w_thl_2 = zero

    # --- higher-order moments (always computed) ---
    wp2rtp = _calc_wp2xp_pdf(wm, rtm, w_1, w_2, rt_1, rt_2, varnce_w_1,
                             varnce_w_2, varnce_rt_1, varnce_rt_2,
                             corr_w_rt_1, corr_w_rt_2, mixt_frac)
    wp2thlp = _calc_wp2xp_pdf(wm, thlm, w_1, w_2, thl_1, thl_2,
                              varnce_w_1, varnce_w_2, varnce_thl_1,
                              varnce_thl_2, corr_w_thl_1, corr_w_thl_2,
                              mixt_frac)

    # --- chi/eta transform and cloud fraction ---
    tl1 = thl_1 * exner
    tl2 = thl_2 * exner
    rsatl_1 = sat_mixrat_liq(p_in_pa, tl1)
    rsatl_2 = sat_mixrat_liq(p_in_pa, tl2)

    chi_1, crt_1, cthl_1, stdev_chi_1, stdev_eta_1, covar_chi_eta_1, \
        corr_chi_eta_1 = _transform_pdf_chi_eta_component(
            tl1, rsatl_1, rt_1, exner, varnce_thl_1, varnce_rt_1,
            corr_rt_thl_1)
    cloud_frac_1, rc_1 = _calc_cloud_frac_component(chi_1, stdev_chi_1,
                                                    0.0)
    chi_2, crt_2, cthl_2, stdev_chi_2, stdev_eta_2, covar_chi_eta_2, \
        corr_chi_eta_2 = _transform_pdf_chi_eta_component(
            tl2, rsatl_2, rt_2, exner, varnce_thl_2, varnce_rt_2,
            corr_rt_thl_2)
    cloud_frac_2, rc_2 = _calc_cloud_frac_component(chi_2, stdev_chi_2,
                                                    0.0)

    # --- ice supersaturation fraction (l_calc_ice_supersat_frac) ---
    chi_at_ice_sat1 = (sat_mixrat_ice(p_in_pa, tl1) - rsatl_1) * crt_1
    isf1_c, rc1_ice_c = _calc_cloud_frac_component(
        chi_1, stdev_chi_1, chi_at_ice_sat1)
    freezing1 = tl1 <= T_FREEZE_K
    ice_supersat_frac_1 = jnp.where(freezing1, isf1_c, cloud_frac_1)
    rc_1_ice = jnp.where(freezing1, rc1_ice_c, rc_1)
    del rc_1_ice  # computed but unused downstream (as in the Fortran)

    chi_at_ice_sat2 = (sat_mixrat_ice(p_in_pa, tl2) - rsatl_2) * crt_2
    isf2_c, rc2_ice_c = _calc_cloud_frac_component(
        chi_2, stdev_chi_2, chi_at_ice_sat2)
    freezing2 = tl2 <= T_FREEZE_K
    ice_supersat_frac_2 = jnp.where(freezing2, isf2_c, cloud_frac_2)
    rc_2_ice = jnp.where(freezing2, rc2_ice_c, rc_2)
    del rc_2_ice

    ice_supersat_frac = mixt_frac * ice_supersat_frac_1 \
        + (1.0 - mixt_frac) * ice_supersat_frac_2

    # --- cloud fraction and mean liquid water ---
    cloud_frac = mixt_frac * cloud_frac_1 \
        + (1.0 - mixt_frac) * cloud_frac_2
    rcm = mixt_frac * rc_1 + (1.0 - mixt_frac) * rc_2
    rcm = jnp.maximum(ZERO_THRESHOLD, rcm)

    # ADG1: w-chi/eta and u/v-w correlations are zero
    corr_w_chi_1 = corr_w_chi_2 = zero
    corr_w_eta_1 = corr_w_eta_2 = zero
    corr_u_w_1 = corr_u_w_2 = zero
    corr_v_w_1 = corr_v_w_2 = zero
    del corr_u_w_1, corr_u_w_2, corr_v_w_1, corr_v_w_2

    # --- x'rc' contributions and th_v moments ---
    wprcp_c1, wp2rcp_c1, rtprcp_c1, thlprcp_c1, uprcp_c1, vprcp_c1 = \
        _calc_xprcp_component(wm, rtm, thlm, um, vm, rcm, w_1, rt_1,
                              thl_1, u_1, v_1, varnce_w_1, chi_1,
                              stdev_chi_1, stdev_eta_1, corr_chi_eta_1,
                              crt_1, cthl_1, rc_1, cloud_frac_1)
    wprcp_c2, wp2rcp_c2, rtprcp_c2, thlprcp_c2, uprcp_c2, vprcp_c2 = \
        _calc_xprcp_component(wm, rtm, thlm, um, vm, rcm, w_2, rt_2,
                              thl_2, u_2, v_2, varnce_w_2, chi_2,
                              stdev_chi_2, stdev_eta_2, corr_chi_eta_2,
                              crt_2, cthl_2, rc_2, cloud_frac_2)

    rc_coef = LV / (exner * CP) - EP2 * thv_ds

    wprcp = mixt_frac * wprcp_c1 + (1.0 - mixt_frac) * wprcp_c2
    wp2rcp = mixt_frac * wp2rcp_c1 + (1.0 - mixt_frac) * wp2rcp_c2
    rtprcp = mixt_frac * rtprcp_c1 + (1.0 - mixt_frac) * rtprcp_c2
    thlprcp = mixt_frac * thlprcp_c1 + (1.0 - mixt_frac) * thlprcp_c2
    uprcp = mixt_frac * uprcp_c1 + (1.0 - mixt_frac) * uprcp_c2
    vprcp = mixt_frac * vprcp_c1 + (1.0 - mixt_frac) * vprcp_c2

    wpthvp = wpthlp + EP1 * thv_ds * wprtp + rc_coef * wprcp
    wp2thvp = wp2thlp + EP1 * thv_ds * wp2rtp + rc_coef * wp2rcp
    rtpthvp = rtpthlp + EP1 * thv_ds * rtp2 + rc_coef * rtprcp
    thlpthvp = thlp2 + EP1 * thv_ds * rtpthlp + rc_coef * thlprcp

    # --- rc'^2 (always computed under CLUBB_CAM) ---
    rcp2 = mixt_frac * (chi_1 * rc_1 + cloud_frac_1 * stdev_chi_1 ** 2) \
        + (1.0 - mixt_frac) * (chi_2 * rc_2
                               + cloud_frac_2 * stdev_chi_2 ** 2) \
        - rcm ** 2
    rcp2 = jnp.maximum(ZERO_THRESHOLD, rcp2)

    moments = dict(
        wp2rtp=wp2rtp, wp2thlp=wp2thlp, cloud_frac=cloud_frac,
        ice_supersat_frac=ice_supersat_frac, rcm=rcm, wpthvp=wpthvp,
        wp2thvp=wp2thvp, rtpthvp=rtpthvp, thlpthvp=thlpthvp,
        wprcp=wprcp, wp2rcp=wp2rcp, rtprcp=rtprcp, thlprcp=thlprcp,
        rcp2=rcp2, uprcp=uprcp, vprcp=vprcp, rc_coef=rc_coef)

    pdf_params = dict(
        w_1=w_1, w_2=w_2, varnce_w_1=varnce_w_1, varnce_w_2=varnce_w_2,
        rt_1=rt_1, rt_2=rt_2, varnce_rt_1=varnce_rt_1,
        varnce_rt_2=varnce_rt_2, thl_1=thl_1, thl_2=thl_2,
        varnce_thl_1=varnce_thl_1, varnce_thl_2=varnce_thl_2,
        corr_w_rt_1=corr_w_rt_1, corr_w_rt_2=corr_w_rt_2,
        corr_w_thl_1=corr_w_thl_1, corr_w_thl_2=corr_w_thl_2,
        corr_rt_thl_1=corr_rt_thl_1, corr_rt_thl_2=corr_rt_thl_2,
        alpha_thl=alpha_thl, alpha_rt=alpha_rt, crt_1=crt_1,
        crt_2=crt_2, cthl_1=cthl_1, cthl_2=cthl_2, chi_1=chi_1,
        chi_2=chi_2, stdev_chi_1=stdev_chi_1, stdev_chi_2=stdev_chi_2,
        stdev_eta_1=stdev_eta_1, stdev_eta_2=stdev_eta_2,
        covar_chi_eta_1=covar_chi_eta_1, covar_chi_eta_2=covar_chi_eta_2,
        corr_w_chi_1=corr_w_chi_1, corr_w_chi_2=corr_w_chi_2,
        corr_w_eta_1=corr_w_eta_1, corr_w_eta_2=corr_w_eta_2,
        corr_chi_eta_1=corr_chi_eta_1, corr_chi_eta_2=corr_chi_eta_2,
        rsatl_1=rsatl_1, rsatl_2=rsatl_2, rc_1=rc_1, rc_2=rc_2,
        cloud_frac_1=cloud_frac_1, cloud_frac_2=cloud_frac_2,
        mixt_frac=mixt_frac, ice_supersat_frac_1=ice_supersat_frac_1,
        ice_supersat_frac_2=ice_supersat_frac_2)

    return moments, pdf_params
