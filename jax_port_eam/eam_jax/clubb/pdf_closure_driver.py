"""CLUBB pdf_closure_driver — port of the private subroutine
pdf_closure_driver (+ its private helpers trapezoidal_rule_zt/zm,
trapezoid_zt/zm, compute_cloud_cover, clip_rcm) from
components/eam/src/physics/clubb/advance_clubb_core_module.F90, as
advance_clubb_core invokes it under the EAMv3 configuration.

PORT_NOTES
----------
Scope is exactly EAMv3's flag set (goldened by
harness/gen_clubb_golden.py against a verbatim extraction of the real
private routine, itself validated BITWISE end-to-end against the
public advance_clubb_core):

- l_call_pdf_closure_twice = .true. (clubb_vert_avg_closure=T): the
  zm-level variables come from a SECOND pdf_closure call on
  momentum-level-interpolated inputs, not from interpolating the first
  call's outputs; the entire `else` interpolation branch is dead.
- l_trapezoidal_rule_zt = .true.: cloud_frac, ice_supersat_frac, rcm
  and wp2thvp are recomputed on thermo levels by the trapezoidal rule
  from their zt+zm values (l_stats=F so the wprtp2/wpthlp2/wprtpthlp
  trapezoids are skipped; l_apply_rule_to_pdf_params is a
  compile-time .false.).
- l_trapezoidal_rule_zm = .true.: wpthvp, thlpthvp, rtpthvp are
  recomputed on momentum levels from their zm+zt values.
- l_use_cloud_cover = .true.: after clip_rcm and compute_cloud_cover,
  cloud_frac := cloud_cover and rcm := rcm_in_layer, then both cloud
  fractions are clipped to <= 1.
- l_use_ice_latent = .false. (the whole _frz third/fourth pdf_closure
  block is dead; rtm_frz/thlm_frz are intent(out) garbage in EAM and
  not returned here), l_rtm_nudge = .false. (rtm passes through
  unchanged), l_rcm_supersat_adj = .false. (rcm_supersat_adj = 0; the
  rsat/rel_humidity computed just before it feed nothing else and are
  skipped), l_refine_grid_in_cloud = .false. (compile-time),
  ipdf_call_placement is irrelevant to this routine itself.
- l_gamma_Skw = .true. with EAMv3 gamma_coef=0.12 /= gamma_coefb=0.28:
  the Skw-dependent gamma branch is always taken (the eps-scaled
  runtime inequality reduces to gamma_coef /= gamma_coefb, handled in
  pdf_closure.gamma_skw_fnc).
- l_stats = .false.: Skw_velocity is gated on iSkw_velocity > 0 and is
  intent(out) garbage in EAM — not returned.  wp4, wprtp2, wpthlp2,
  wprtpthlp (and their _zm) are never computed by pdf_closure under
  the EAM configuration (see pdf_closure.py) — not returned.
- sclr_dim = 0, hydromet_dim = 0: all scalar/hydrometeor plumbing is
  empty.

Level placement of inputs (CLUBB orientation, index 0 = surface
ghost/surface): momentum-level: wprtp, wpthlp, rtp2, thlp2, rtpthlp,
wp2, up2, upwp, vp2, vpwp, wm_zm, thv_ds_zm; thermo-level: thlm, rtm,
rtp3, thlp3, wp3, wm_zt, um, vm, p_in_pa, exner, thv_ds_zt.

compute_cloud_cover's final `else` error branch is unreachable (its
guard is the exact complement of the three taken branches) and is not
represented.  clip_rcm uses the Fortran epsilon(x) for double
(EPSILON_R8 = 2**-52) as an absolute decrement.

Numerical-faithfulness: divisions by the scalar constant p0 go
through _divs (see pdf_closure.py); exner_zm passes through libm pow
(cross-libm 1-ulp tails possible); everything else is composition of
slice-A kernels (pdf_closure bitwise except erf/exp tails) and linear
grid operations.
"""

import jax.numpy as jnp

from . import grid as cgrid
from .constants import (EPSILON_R8, KAPPA, P0, RC_TOL, RT_TOL, THL_TOL,
                        W_TOL, W_TOL_SQD, ZERO_THRESHOLD)
from .pdf_closure import (_divs, compute_sigma_sqd_w, gamma_skw_fnc,
                          pdf_closure, skx_func)


def trapezoid_zt(gr, variable_zt, variable_zm):
    """trapezoid_zt: recompute a thermo-level variable with the
    trapezoidal rule (not valid at zt level 1)."""
    interior = (0.5 * (variable_zm[1:] + variable_zt[1:])
                * (gr.zm[1:] - gr.zt[1:]) * gr.invrs_dzt[1:]
                + 0.5 * (variable_zt[1:] + variable_zm[:-1])
                * (gr.zt[1:] - gr.zm[:-1]) * gr.invrs_dzt[1:])
    return jnp.concatenate([variable_zt[:1], interior])


def trapezoid_zm(gr, variable_zm, variable_zt):
    """trapezoid_zm: recompute a momentum-level variable with the
    trapezoidal rule (not valid at zm levels 1 and nz)."""
    interior = (0.5 * (variable_zt[2:] + variable_zm[1:-1])
                * (gr.zt[2:] - gr.zm[1:-1]) * gr.invrs_dzm[1:-1]
                + 0.5 * (variable_zm[1:-1] + variable_zt[1:-1])
                * (gr.zm[1:-1] - gr.zt[1:-1]) * gr.invrs_dzm[1:-1])
    return jnp.concatenate([variable_zm[:1], interior,
                            variable_zm[-1:]])


def clip_rcm(rtm, rcm):
    """clip_rcm: reduce rcm wherever it exceeds rtm (prevents
    rvm < 0); rcm := max(0, rtm - epsilon(rtm))."""
    return jnp.where(rtm < rcm,
                     jnp.maximum(ZERO_THRESHOLD, rtm - EPSILON_R8), rcm)


def compute_cloud_cover(gr, mixt_frac, chi_1, chi_2, cloud_frac, rcm):
    """compute_cloud_cover: cloud cover and in-layer liquid water.

    Interior levels (k = 2..nz-1 in Fortran) take one of three
    branches: no cloud at k / cloud above and below / cloud fails to
    reach gridbox top or base (the partial-fill branch, with its
    gradual-transition blending).  Boundary levels copy cloud_frac and
    rcm.  The Fortran's final `else` error branch is logically
    unreachable."""
    chi_mean = mixt_frac * chi_1 + (1.0 - mixt_frac) * chi_2

    r_k = rcm[1:-1]
    cf_k = cloud_frac[1:-1]
    r_up = rcm[2:]
    r_dn = rcm[:-2]
    chi_up = chi_mean[2:]
    chi_dn = chi_mean[:-2]

    no_cloud = r_k < RC_TOL
    full = (r_up >= RC_TOL) & (r_dn >= RC_TOL)
    partial = ~no_cloud & ~full

    # cloud top: rcm(k+1) < rc_tol
    denom_up = jnp.where(partial, r_k + jnp.abs(chi_up), 1.0)
    vcf_up = ((0.5 / gr.invrs_dzm[1:-1]) / (gr.zm[1:-1] - gr.zt[1:-1])
              * (r_k / denom_up))
    vcf_up = jnp.minimum(0.5, vcf_up)
    vcf_up = vcf_up + _divs(r_up, RC_TOL) * (0.5 - vcf_up)
    vcf_up = jnp.where(r_up < RC_TOL, vcf_up, 0.5)

    # cloud base: rcm(k-1) < rc_tol
    denom_dn = jnp.where(partial, r_k + jnp.abs(chi_dn), 1.0)
    vcf_dn = ((0.5 / gr.invrs_dzm[:-2]) / (gr.zt[1:-1] - gr.zm[:-2])
              * (r_k / denom_dn))
    vcf_dn = jnp.minimum(0.5, vcf_dn)
    vcf_dn = vcf_dn + _divs(r_dn, RC_TOL) * (0.5 - vcf_dn)
    vcf_dn = jnp.where(r_dn < RC_TOL, vcf_dn, 0.5)

    vcf = vcf_up + vcf_dn
    vcf = jnp.maximum(cf_k, jnp.minimum(1.0, vcf))
    vcf_safe = jnp.where(partial, vcf, 1.0)

    cc_k = jnp.where(partial, cf_k / vcf_safe, cf_k)
    ril_k = jnp.where(partial, r_k / vcf_safe, r_k)

    cloud_cover = jnp.concatenate([cloud_frac[:1], cc_k,
                                   cloud_frac[-1:]])
    rcm_in_layer = jnp.concatenate([rcm[:1], ril_k, rcm[-1:]])
    return cloud_cover, rcm_in_layer


def pdf_closure_driver(gr, wprtp, thlm, wpthlp, rtp2, rtp3, thlp2,
                       thlp3, rtpthlp, wp2, wp3, wm_zm, wm_zt, um, up2,
                       upwp, vm, vp2, vpwp, p_in_pa, exner, thv_ds_zm,
                       thv_ds_zt, rtm, beta, mf_max_mag, gamma_coef,
                       gamma_coefb, gamma_coefc):
    """pdf_closure_driver under the EAMv3 flag set (see module
    docstring).  All arrays (nz,) on the grid `gr` (a grid.Grid).
    Tunables: beta, mf_max_mag (= mixt_frac_max_mag(Skw_max_mag)),
    gamma_coef/b/c.  dt is not an argument: it only feeds the
    l_rtm_nudge branch, which is off (rtm passes through unchanged and
    is not returned).

    Returns (outs, pdf_params, pdf_params_zm): outs maps the 29
    defined output fields (the golden's outs_slots order), pdf_params
    / pdf_params_zm are the 47-field dicts of the zt / zm pdf_closure
    calls."""
    # ---- skewness on both grids -------------------------------------
    wp2_zt = jnp.maximum(cgrid.zm2zt(gr, wp2), W_TOL_SQD)
    wp3_zm = cgrid.zt2zm(gr, wp3)
    thlp2_zt = jnp.maximum(cgrid.zm2zt(gr, thlp2), THL_TOL ** 2)
    thlp3_zm = cgrid.zt2zm(gr, thlp3)
    rtp2_zt = jnp.maximum(cgrid.zm2zt(gr, rtp2), RT_TOL ** 2)
    rtp3_zm = cgrid.zt2zm(gr, rtp3)

    skw_zt = skx_func(wp2_zt, wp3, W_TOL, 0.0)
    skw_zm = skx_func(wp2, wp3_zm, W_TOL, 0.0)
    skthl_zt = skx_func(thlp2_zt, thlp3, THL_TOL, 0.0)
    skthl_zm = skx_func(thlp2, thlp3_zm, THL_TOL, 0.0)
    skrt_zt = skx_func(rtp2_zt, rtp3, RT_TOL, 0.0)
    skrt_zm = skx_func(rtp2, rtp3_zm, RT_TOL, 0.0)

    # ---- sigma_sqd_w (+ vertical smoothing) -------------------------
    gamma = gamma_skw_fnc(skw_zm, gamma_coef, gamma_coefb, gamma_coefc)
    sigma_sqd_w = compute_sigma_sqd_w(gamma, wp2, thlp2, rtp2, up2,
                                      vp2, wpthlp, wprtp, upwp, vpwp)
    sigma_sqd_w = cgrid.zt2zm(gr, cgrid.zm2zt(gr, sigma_sqd_w))
    sigma_sqd_w = jnp.maximum(ZERO_THRESHOLD, sigma_sqd_w)
    sigma_sqd_w_zt = jnp.maximum(cgrid.zm2zt(gr, sigma_sqd_w),
                                 ZERO_THRESHOLD)

    rtpthlp_zt = cgrid.zm2zt(gr, rtpthlp)

    # ---- first pdf_closure call (thermodynamic levels) --------------
    mom_zt, pdfp_zt = pdf_closure(
        p_in_pa, exner, thv_ds_zt, wm_zt, wp2_zt, wp3, sigma_sqd_w_zt,
        skw_zt, skthl_zt, skrt_zt, rtm, rtp2_zt,
        cgrid.zm2zt(gr, wprtp), thlm, thlp2_zt,
        cgrid.zm2zt(gr, wpthlp), um, cgrid.zm2zt(gr, up2),
        cgrid.zm2zt(gr, upwp), vm, cgrid.zm2zt(gr, vp2),
        cgrid.zm2zt(gr, vpwp), rtpthlp_zt, beta, mf_max_mag)

    # ---- second pdf_closure call (momentum levels) -------------------
    p_zm = cgrid.zt2zm(gr, p_in_pa)
    p_zm = p_zm.at[0].set(p_in_pa[0])
    p_zm = p_zm.at[-1].set(jnp.maximum(p_zm[-1], 0.5 * p_in_pa[-1]))
    exner_zm = _divs(p_zm, P0) ** KAPPA
    rtm_zm = cgrid.zt2zm(gr, rtm)
    rtm_zm = rtm_zm.at[-1].set(jnp.maximum(rtm_zm[-1], RT_TOL))
    thlm_zm = cgrid.zt2zm(gr, thlm)
    thlm_zm = thlm_zm.at[-1].set(jnp.maximum(thlm_zm[-1], THL_TOL))

    mom_zm, pdfp_zm = pdf_closure(
        p_zm, exner_zm, thv_ds_zm, wm_zm, wp2, wp3_zm, sigma_sqd_w,
        skw_zm, skthl_zm, skrt_zm, rtm_zm, rtp2, wprtp, thlm_zm,
        thlp2, wpthlp, cgrid.zt2zm(gr, um), up2, upwp,
        cgrid.zt2zm(gr, vm), vp2, vpwp, rtpthlp, beta, mf_max_mag)

    # ---- trapezoidal rules ------------------------------------------
    cloud_frac = trapezoid_zt(gr, mom_zt["cloud_frac"],
                              mom_zm["cloud_frac"])
    ice_supersat_frac = trapezoid_zt(gr, mom_zt["ice_supersat_frac"],
                                     mom_zm["ice_supersat_frac"])
    rcm = trapezoid_zt(gr, mom_zt["rcm"], mom_zm["rcm"])
    wp2thvp = trapezoid_zt(gr, mom_zt["wp2thvp"], mom_zm["wp2thvp"])

    wpthvp = trapezoid_zm(gr, mom_zm["wpthvp"], mom_zt["wpthvp"])
    thlpthvp = trapezoid_zm(gr, mom_zm["thlpthvp"], mom_zt["thlpthvp"])
    rtpthvp = trapezoid_zm(gr, mom_zm["rtpthvp"], mom_zt["rtpthvp"])

    # ---- clip_rcm, cloud cover, l_use_cloud_cover -------------------
    rcm = clip_rcm(rtm, rcm)

    cloud_cover, rcm_in_layer = compute_cloud_cover(
        gr, pdfp_zt["mixt_frac"], pdfp_zt["chi_1"], pdfp_zt["chi_2"],
        cloud_frac, rcm)

    cloud_frac = cloud_cover
    rcm = rcm_in_layer

    cloud_frac = jnp.minimum(1.0, cloud_frac)
    ice_supersat_frac = jnp.minimum(1.0, ice_supersat_frac)

    rcm_supersat_adj = jnp.zeros_like(rcm)

    outs = dict(
        rcm=rcm, cloud_frac=cloud_frac,
        ice_supersat_frac=ice_supersat_frac, wprcp=mom_zm["wprcp"],
        sigma_sqd_w=sigma_sqd_w, wpthvp=wpthvp, wp2thvp=wp2thvp,
        rtpthvp=rtpthvp, thlpthvp=thlpthvp, rc_coef=mom_zt["rc_coef"],
        rcm_in_layer=rcm_in_layer, cloud_cover=cloud_cover,
        rcp2_zt=mom_zt["rcp2"], thlprcp=mom_zm["thlprcp"],
        rc_coef_zm=mom_zm["rc_coef"], wp2rtp=mom_zt["wp2rtp"],
        wp2thlp=mom_zt["wp2thlp"], wp2rcp=mom_zt["wp2rcp"],
        rtprcp=mom_zm["rtprcp"], rcp2=mom_zm["rcp2"],
        uprcp=mom_zm["uprcp"], vprcp=mom_zm["vprcp"],
        cloud_frac_zm=mom_zm["cloud_frac"],
        ice_supersat_frac_zm=mom_zm["ice_supersat_frac"],
        rtm_zm=rtm_zm, thlm_zm=thlm_zm, rcm_zm=mom_zm["rcm"],
        rcm_supersat_adj=rcm_supersat_adj,
        sigma_sqd_w_zt=sigma_sqd_w_zt)

    return outs, pdfp_zt, pdfp_zm
