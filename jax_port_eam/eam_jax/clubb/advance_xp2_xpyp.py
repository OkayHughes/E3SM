"""CLUBB slice C: advance_xp2_xpyp + clip_covars_denom.

PORT_NOTES
----------
Sources (components/eam/src/physics/clubb/):
  advance_xp2_xpyp_module.F90  advance_xp2_xpyp, xp2_xpyp_lhs,
                               xp2_xpyp_rhs, xp2_xpyp_uv_rhs,
                               xp2_xpyp_solve, term_tp, term_dp1_lhs,
                               term_dp1_rhs, term_pr1, term_pr2,
                               pos_definite_variances
  turbulent_adv_pdf.F90        xpyp_term_ta_pdf_lhs_all/rhs_all
                               ("upwind" branch)
  diffusion.F90                diffusion_zm_lhs_all
  mean_adv.F90                 term_ma_zm_lhs_all
  clip_explicit.F90            clip_variance, clip_covar,
                               clip_covars_denom
  fill_holes.F90               fill_holes_vertical,
                               fill_holes_multiplicative, vertical_avg
                               (the "zm" path pos_definite_variances
                               uses)

EAMv3-active configuration (compile-time model_flags defaults +
clubb_intr settings; asserted against the Fortran module state in
harness/gen_clubb_golden.py gen_xp2 via drv_xp2_config):
  l_iter_xp2_xpyp    = T  (compile-time parameter in advance_clubb_core)
  l_single_C2_Skw    = F  -> C2rt/C2thl/C2rtthl are the constant
                            tunables (no Skw_zm dependence)
  l_C2_cloud_frac    = F  -> no cloud_frac modulation of C2
  l_explicit_turbulent_adv_xpyp = F and iiPDF_type = ADG1 (compile
                            time) -> the ADG1 semi-implicit turbulent
                            advection coefficients
  l_upwind_xpyp_ta   = T  -> "upwind" ta discretization; the turbulent
                            velocity sign array is wp3_on_wp2 itself
                            (sign tested as > 0)
  l_min_xp2_from_corr_wx = F -> plain tolerance clipping of rtp2/thlp2
  l_hole_fill        = T  -> pos_definite_variances (fill_holes
                            num_draw_pts=2, "zm") on all 4 variances
  l_clip_large_rtp2  = T  (compile-time parameter), rtp2_clip_coef=0.5
  up2_vp2_sponge_damp_settings%l_sponge_damping = F (never assigned in
                            an EAM build; static zero-init)
  l_tke_aniso        = T  (clip_covars_denom clips upwp/vpwp against
                            up2/vp2, not wp2)
  l_predict_upwp_vpwp = F, sclr_dim = 0, l_stats = F,
  l_clip_large_neg_mc = F (compile-time local in xp2_xpyp_rhs)

DEAD arguments dropped from the port signature (real arguments of the
Fortran routine, unused under the flags above; the goldens record them
so the scope is checkable): Skw_zm and cloud_frac (only the
l_single_C2_Skw / l_C2_cloud_frac branches), wprtp2/wpthlp2/wprtpthlp
(only l_explicit_turbulent_adv_xpyp=T), wp3 (debug print only), Lscale
and the wp2_zt(k+1)/wp2_zt(k) args of term_pr2 (only the compile-time
dead l_use_experimental_term_pr2 branch), pdf_implicit_coefs_terms
(only iiPDF_new).  Furthermore, under l_upwind_xpyp_ta=T the entire
zt-side ADG1 coefficient/term set (coef_*_implicit,
term_*_explicit, a1_zt, wprtp_zt/wpthlp_zt/upwp_zt/vpwp_zt, and their
wp2_zt / wp3_on_wp2_zt inputs) is computed by the Fortran but NEVER
READ — the upwind discretization uses only the *_zm arrays — so
wp2_zt and wp3_on_wp2_zt are dropped too (verified: every
xpyp_term_ta_pdf_{lhs,rhs}_all upwind branch reads only coef/term_zm).
rho_ds_zt is likewise dead on this path (centered-ta discretization
and the fill_holes "zt" branch only) and is dropped.

Numerical-fidelity notes:
  * Fortran objects are built -ffp-contract=off; expressions here
    reproduce the exact multiplication/addition order of the source
    (e.g. term_pr1 divides by tau_zm LAST; uv term_tp is (-t) - t with
    two identical products; vp2's term_pr2 subtracts the vpwp shear
    term BEFORE the upwp one — opposite order from up2's).
  * Divisions by broadcast scalars go through _divs (XLA otherwise
    rewrites x/const into x*(1/const), which is not correctly rounded).
  * The vertical averages in fill_holes accumulate STRICTLY
    left-to-right (lax.scan), matching the Fortran do-loop; jnp.sum's
    pairwise order would break bitwise agreement of mass_fraction.
  * fill_holes' window loop is sequential (window k reads values that
    window k-1 may have modified); it is a Python loop over the static
    nz with traced jnp.where gates, so an untriggered window is an
    exact no-op (matching the Fortran `if any(...)` gates).
  * The tridiagonal solves reuse slice A's exact dgtsv port
    (tridiag.tridag_solve); up2/vp2 share one matrix with 2 RHS,
    exactly like the Fortran xp2_xpyp_solve( xp2_xpyp_up2_vp2, 2, ...).

Tunables read (from the EAMv3 params vector; see
gen_clubb_golden.py): C2rt, C2thl, C2rtthl, C4, C5, C14, beta, c_K2,
c_K9, plus the setup_parameters-derived nu2/nu9_vert_res_dep profiles
(constant profiles under l_adj_low_res_nu=T with scalar mult_factor;
taken as inputs, captured from the Fortran module).
"""

import jax
import jax.numpy as jnp

from .constants import (EPS, GAMMA_OVER_IMPLICIT_TS, GRAV,
                        MAX_MAG_CORRELATION, MAX_MAG_CORRELATION_FLUX,
                        RT_TOL, THL_TOL, W_TOL_SQD, ZERO_THRESHOLD)
from .tridiag import tridag_solve

ONE_THIRD = 1.0 / 3.0
TWO_THIRDS = 2.0 / 3.0
RTP2_CLIP_COEF = 0.5    # l_clip_large_rtp2 coefficient (compile-time)


# ---------------------------------------------------------------------------
# LHS building blocks (interior levels k = 1..nz-2, 0-based; the
# Fortran *_all routines zero the boundary rows and
# set_boundary_conditions_lhs then overwrites them with identity rows)
# ---------------------------------------------------------------------------

def _diffusion_zm_lhs_interior(gr, k_zt, nu):
    """diffusion_zm_lhs_all, interior rows (super, main, sub)."""
    idzm = gr.invrs_dzm[1:-1]
    a_up = (k_zt[2:] + nu[2:]) * gr.invrs_dzt[2:]      # (K+nu)(k+1)*idzt(k+1)
    a_dn = (k_zt[1:-1] + nu[1:-1]) * gr.invrs_dzt[1:-1]
    sup = -idzm * (k_zt[2:] + nu[2:]) * gr.invrs_dzt[2:]
    dia = idzm * (a_up + a_dn)
    sub = -idzm * (k_zt[1:-1] + nu[1:-1]) * gr.invrs_dzt[1:-1]
    return sup, dia, sub


def _ma_zm_lhs_interior(gr, wm_zm):
    """term_ma_zm_lhs_all, interior rows.  weights_zm2zt(1,:) is
    grid.factor_zm2zt; weights_zm2zt(2,:) = 1 - factor."""
    w1 = gr.factor_zm2zt
    w2 = 1.0 - gr.factor_zm2zt
    wk = wm_zm[1:-1] * gr.invrs_dzm[1:-1]
    sup = wk * w1[2:]
    dia = wk * (w2[2:] - w1[1:-1])
    sub = -(wk * w2[1:-1])
    return sup, dia, sub


def _ta_pdf_lhs_interior(gr, sgn, coef_zm, rho_ds_zm, invrs_rho_ds_zm):
    """xpyp_term_ta_pdf_lhs_all, "upwind" branch, interior rows."""
    inv = invrs_rho_ds_zm[1:-1]
    idzt = gr.invrs_dzt[1:-1]
    idztp1 = gr.invrs_dzt[2:]
    up = sgn[1:-1] > 0.0
    # wind blowing upward
    sup_u = jnp.zeros_like(inv)
    dia_u = inv * idzt * rho_ds_zm[1:-1] * coef_zm[1:-1]
    sub_u = -(inv * idzt * rho_ds_zm[:-2] * coef_zm[:-2])
    # wind blowing downward
    sup_d = inv * idztp1 * rho_ds_zm[2:] * coef_zm[2:]
    dia_d = -(inv * idztp1 * rho_ds_zm[1:-1] * coef_zm[1:-1])
    sub_d = jnp.zeros_like(inv)
    return (jnp.where(up, sup_u, sup_d), jnp.where(up, dia_u, dia_d),
            jnp.where(up, sub_u, sub_d))


def _ta_pdf_rhs_interior(gr, sgn, term_zm, rho_ds_zm, invrs_rho_ds_zm):
    """xpyp_term_ta_pdf_rhs_all, "upwind" branch, interior values."""
    inv = invrs_rho_ds_zm[1:-1]
    up = sgn[1:-1] > 0.0
    rhs_u = -(inv * gr.invrs_dzt[1:-1]
              * (rho_ds_zm[1:-1] * term_zm[1:-1]
                 - rho_ds_zm[:-2] * term_zm[:-2]))
    rhs_d = -(inv * gr.invrs_dzt[2:]
              * (rho_ds_zm[2:] * term_zm[2:]
                 - rho_ds_zm[1:-1] * term_zm[1:-1]))
    return jnp.where(up, rhs_u, rhs_d)


def _xp2_xpyp_lhs(gr, dt, coef_zm, sgn, tau_zm, wm_zm, kw, rho_ds_zm,
                  invrs_rho_ds_zm, cn, nu, l_iter):
    """xp2_xpyp_lhs (l_upwind_xpyp_ta=T): returns (supd, diag, subd)
    full-length rows with the identity boundary rows applied."""
    nz = tau_zm.shape[0]
    d_sup, d_dia, d_sub = _diffusion_zm_lhs_interior(gr, kw, nu)
    m_sup, m_dia, m_sub = _ma_zm_lhs_interior(gr, wm_zm)
    t_sup, t_dia, t_sub = _ta_pdf_lhs_interior(gr, sgn, coef_zm,
                                               rho_ds_zm, invrs_rho_ds_zm)
    lhs_dp1 = (cn[1:-1] / tau_zm[1:-1]) * GAMMA_OVER_IMPLICIT_TS

    sup = d_sup + m_sup + t_sup * GAMMA_OVER_IMPLICIT_TS
    dia = (d_dia + m_dia + t_dia * GAMMA_OVER_IMPLICIT_TS) + lhs_dp1
    sub = d_sub + m_sub + t_sub * GAMMA_OVER_IMPLICIT_TS
    if l_iter:
        dia = dia + (1.0 / dt)

    zero = jnp.zeros((1,), dtype=tau_zm.dtype)
    one = jnp.ones((1,), dtype=tau_zm.dtype)
    supd = jnp.concatenate([zero, sup, zero])
    diag = jnp.concatenate([one, dia, one])
    subd = jnp.concatenate([zero, sub, zero])
    del nz
    return supd, diag, subd


# ---------------------------------------------------------------------------
# RHS
# ---------------------------------------------------------------------------

def _term_tp(xam, xbm, wpxap, wpxbp, invrs_dzm_i):
    """term_tp on interior momentum levels: xam/xbm are FULL zt
    arrays; returns interior values.  Fortran:
      - wpxbp*invrs_dzm*(xam(k+1)-xam(k)) - wpxap*invrs_dzm*(xbm(k+1)-xbm(k))
    """
    return (-(wpxbp[1:-1] * invrs_dzm_i * (xam[2:] - xam[1:-1]))
            - wpxap[1:-1] * invrs_dzm_i * (xbm[2:] - xbm[1:-1]))


def _xp2_xpyp_rhs(gr, dt, coef_zm, term_zm, sgn, wpxap, wpxbp, xam,
                  xbm, xapxbp, xpyp_forcing, rho_ds_zm,
                  invrs_rho_ds_zm, cn, tau_zm, threshold, l_iter):
    """xp2_xpyp_rhs (l_upwind_xpyp_ta=T, l_clip_large_neg_mc=F)."""
    idzm = gr.invrs_dzm[1:-1]
    rhs_ta = _ta_pdf_rhs_interior(gr, sgn, term_zm, rho_ds_zm,
                                  invrs_rho_ds_zm)
    t_sup, t_dia, t_sub = _ta_pdf_lhs_interior(gr, sgn, coef_zm,
                                               rho_ds_zm, invrs_rho_ds_zm)
    x_kp1, x_k, x_km1 = xapxbp[2:], xapxbp[1:-1], xapxbp[:-2]

    rhs = rhs_ta + (1.0 - GAMMA_OVER_IMPLICIT_TS) \
        * ((-(t_sup * x_kp1) - t_dia * x_k) - t_sub * x_km1)
    rhs = rhs + _term_tp(xam, xbm, wpxap, wpxbp, idzm)
    rhs = rhs + (cn[1:-1] / tau_zm[1:-1]) * threshold
    rhs = rhs + (1.0 - GAMMA_OVER_IMPLICIT_TS) \
        * (-((cn[1:-1] / tau_zm[1:-1]) * x_k))
    rhs = rhs + xpyp_forcing[1:-1]
    if l_iter:
        # Fortran: rhs(k) + one/dt * xapxbp(k)  (reciprocal-multiply)
        rhs = rhs + (1.0 / dt) * x_k

    return jnp.concatenate([xapxbp[:1], rhs,
                            jnp.full((1,), threshold,
                                     dtype=xapxbp.dtype)])


def _xp2_xpyp_uv_rhs(gr, dt, coef_zm, term_zm, sgn, wp2, wpthvp, tau_zm,
                     xam, xbm, wpxap, wpxbp, xap2, xbp2, rho_ds_zm,
                     invrs_rho_ds_zm, thv_ds_zm, c4, c5, c14,
                     c4_c14_1d, wp2_splat, l_iter):
    """xp2_xpyp_uv_rhs (original term_pr2, l_upwind_xpyp_ta=T).
    Argument roles mirror the Fortran call site: for up2, xam=um,
    xbm=vm, wpxap=upwp, wpxbp=vpwp, xap2=up2, xbp2=vp2; for vp2 they
    are all swapped (which also swaps the fp ORDER of the two shear
    terms inside term_pr2)."""
    idzm = gr.invrs_dzm[1:-1]
    rhs_ta = _ta_pdf_rhs_interior(gr, sgn, term_zm, rho_ds_zm,
                                  invrs_rho_ds_zm)
    t_sup, t_dia, t_sub = _ta_pdf_lhs_interior(gr, sgn, coef_zm,
                                               rho_ds_zm, invrs_rho_ds_zm)
    x_kp1, x_k, x_km1 = xap2[2:], xap2[1:-1], xap2[:-2]

    rhs = rhs_ta - 0.5 * wp2_splat[1:-1]
    rhs = rhs + (1.0 - GAMMA_OVER_IMPLICIT_TS) \
        * ((-(t_sup * x_kp1) - t_dia * x_k) - t_sub * x_km1)

    # turbulent production: term_tp(xam(k+1), xam(k), xam(k+1), xam(k),
    #                               wpxap(k), wpxap(k), invrs_dzm(k))
    # = (-t) - t with two IDENTICAL products (kept as such).
    t = wpxap[1:-1] * idzm * (xam[2:] - xam[1:-1])
    rhs = rhs + (1.0 - c5) * ((-t) - t)

    # term_pr1 (combined dp1/pr1 explicit part)
    rhs = rhs + (ONE_THIRD * (c4 - c14) * (xbp2[1:-1] + wp2[1:-1])
                 / tau_zm[1:-1]
                 + (c14 / tau_zm[1:-1]) * W_TOL_SQD)

    # balance for the over-implicit combined dp1/pr1 LHS term
    rhs = rhs + (1.0 - GAMMA_OVER_IMPLICIT_TS) \
        * (-((c4_c14_1d[1:-1] / tau_zm[1:-1]) * x_k))

    # term_pr2 (original version), floored at zero_threshold
    pr2 = TWO_THIRDS * c5 * (
        ((GRAV / thv_ds_zm[1:-1]) * wpthvp[1:-1]
         - wpxap[1:-1] * idzm * (xam[2:] - xam[1:-1]))
        - wpxbp[1:-1] * idzm * (xbm[2:] - xbm[1:-1]))
    rhs = rhs + jnp.maximum(pr2, ZERO_THRESHOLD)

    if l_iter:
        # Fortran: rhs(k) + one/dt * xap2(k)  (reciprocal-multiply)
        rhs = rhs + (1.0 / dt) * x_k

    return jnp.concatenate([xap2[:1], rhs,
                            jnp.full((1,), W_TOL_SQD, dtype=xap2.dtype)])


# ---------------------------------------------------------------------------
# fill_holes ("zm" path used by pos_definite_variances)
# ---------------------------------------------------------------------------

def _seq_sum(x):
    """Strict left-to-right accumulation (Fortran do-loop order)."""
    total, _ = jax.lax.scan(lambda c, v: (c + v, None), jnp.zeros((), x.dtype), x)
    return total


def _fill_holes_multiplicative(threshold, rho, dz, field):
    """fill_holes_multiplicative on one window; returns the filled
    field (caller gates on `any(field < threshold)`)."""
    w = rho * dz
    field_avg = _seq_sum(w * field) / _seq_sum(w)
    clipped = jnp.where(field_avg >= threshold,
                        jnp.maximum(threshold, field),
                        jnp.minimum(threshold, field))
    clipped_avg = _seq_sum(w * clipped) / _seq_sum(w)
    # round-off guard: return unchanged (Fortran early return)
    degenerate = (jnp.abs(clipped_avg - threshold)
                  <= jnp.abs(clipped_avg + threshold) * EPS / 2.0)
    mass_fraction = (field_avg - threshold) \
        / jnp.where(degenerate, 1.0, clipped_avg - threshold)
    filled = mass_fraction * (clipped - threshold) + threshold
    return jnp.where(degenerate, field, filled)


def fill_holes_vertical_zm(gr, num_draw_pts, threshold, rho_ds_zm, field):
    """fill_holes_vertical for a zm-grid field (the only path
    pos_definite_variances uses; num_draw_pts=2 there).
    upper_hf_level = nz-1 (Fortran) — the top zm level is never
    touched, nor is the surface level."""
    nz = field.shape[0]
    n = num_draw_pts
    # local pass: Fortran k = 2+n .. (nz-1)-n (1-based) — sequential,
    # each window sees the previous window's updates
    for kk in range(1 + n, nz - 1 - n):     # 0-based center
        lo, hi = kk - n, kk + n + 1
        window = field[lo:hi]
        gate = jnp.any(window < threshold)
        filled = _fill_holes_multiplicative(
            threshold, rho_ds_zm[lo:hi], gr.dzm[lo:hi], window)
        field = field.at[lo:hi].set(jnp.where(gate, filled, window))
    # global pass over Fortran 2..nz-1 (0-based 1..nz-2)
    window = field[1:nz - 1]
    gate = jnp.any(window < threshold)
    filled = _fill_holes_multiplicative(
        threshold, rho_ds_zm[1:nz - 1], gr.dzm[1:nz - 1], window)
    field = field.at[1:nz - 1].set(jnp.where(gate, filled, window))
    return field


def pos_definite_variances(gr, tolerance, rho_ds_zm, xp2):
    """pos_definite_variances: hole filling, gated (exact no-op when
    no level is below tolerance)."""
    gate = jnp.any(xp2 < tolerance)
    filled = fill_holes_vertical_zm(gr, 2, tolerance, rho_ds_zm, xp2)
    return jnp.where(gate, filled, xp2)


# ---------------------------------------------------------------------------
# clipping
# ---------------------------------------------------------------------------

def clip_variance(threshold, xp2):
    """clip_explicit clip_variance: floor at threshold on Fortran
    levels 1..nz-1 (the top level is NOT clipped; the surface is)."""
    clipped = jnp.where(xp2 < threshold, threshold, xp2)
    return jnp.concatenate([clipped[:-1], xp2[-1:]])


def clip_covar(max_mag_corr, xp2, yp2, xpyp):
    """clip_explicit clip_covar: bound |xpyp| by
    max_mag_corr*sqrt(xp2*yp2) on interior levels only."""
    bound = max_mag_corr * jnp.sqrt(xp2 * yp2)
    clipped = jnp.where(xpyp > bound, bound,
                        jnp.where(xpyp < -bound, -bound, xpyp))
    return jnp.concatenate([xpyp[:1], clipped[1:-1], xpyp[-1:]])


def clip_covars_denom(wp2, rtp2, thlp2, up2, vp2, wprtp, wpthlp, upwp,
                      vpwp):
    """clip_explicit clip_covars_denom under EAMv3: sclr_dim=0,
    l_tke_aniso=T (upwp/vpwp clipped against up2/vp2), pert pointers
    unassociated.  The *_cl_num arguments of the Fortran routine only
    gate statistics (l_stats=F) and are dropped."""
    wprtp = clip_covar(MAX_MAG_CORRELATION_FLUX, wp2, rtp2, wprtp)
    wpthlp = clip_covar(MAX_MAG_CORRELATION_FLUX, wp2, thlp2, wpthlp)
    upwp = clip_covar(MAX_MAG_CORRELATION, wp2, up2, upwp)
    vpwp = clip_covar(MAX_MAG_CORRELATION, wp2, vp2, vpwp)
    return wprtp, wpthlp, upwp, vpwp


# ---------------------------------------------------------------------------
# main routine
# ---------------------------------------------------------------------------

def advance_xp2_xpyp(gr, dt, tau_zm, wm_zm, rtm, wprtp, thlm, wpthlp,
                     wpthvp, um, vm, wp2, upwp, vpwp,
                     sigma_sqd_w, kh_zt, rtp2_forcing, thlp2_forcing,
                     rtpthlp_forcing, rho_ds_zm,
                     invrs_rho_ds_zm, thv_ds_zm, wp3_on_wp2,
                     wp2_splat, rtp2, thlp2, rtpthlp,
                     up2, vp2, nu2_vert_res_dep, nu9_vert_res_dep,
                     tunables, l_iter=True, l_hole_fill=True):
    """advance_xp2_xpyp under the EAMv3 configuration (see module
    docstring for the flag set and the dropped dead arguments).

    tunables: mapping with C2rt, C2thl, C2rtthl, C4, C5, C14, beta,
    c_K2, c_K9 (floats from the EAMv3 params vector).
    l_hole_fill=False exists only for Tier-0 branch-coverage tests.

    Returns (rtp2, thlp2, rtpthlp, up2, vp2) after the solve, hole
    filling, and all clipping steps.
    """
    c2rt, c2thl, c2rtthl = (tunables["C2rt"], tunables["C2thl"],
                            tunables["C2rtthl"])
    c4, c5, c14 = tunables["C4"], tunables["C5"], tunables["C14"]
    beta = tunables["beta"]

    nz = tau_zm.shape[0]
    dtype = tau_zm.dtype

    # constant C2 profiles (l_single_C2_Skw=F, l_C2_cloud_frac=F)
    c2rt_1d = jnp.full(nz, c2rt, dtype=dtype)
    c2thl_1d = jnp.full(nz, c2thl, dtype=dtype)
    c2rtthl_1d = jnp.full(nz, c2rtthl, dtype=dtype)
    c4_c14_1d = jnp.full(nz, (TWO_THIRDS * c4) + (ONE_THIRD * c14),
                         dtype=dtype)

    # a_1 = 1 / (1 - sigma_sqd_w) on momentum levels.  (The Fortran
    # also builds a1_zt and the whole zt-side coefficient set; all of
    # it is unread under l_upwind_xpyp_ta=T — see module docstring.)
    a1 = 1.0 / (1.0 - sigma_sqd_w)

    # ADG1 semi-implicit turbulent-advection coefficients on momentum
    # levels (l_explicit_turbulent_adv_xpyp=F, iiPDF_ADG1).  The same
    # coef expression serves rtp2/thlp2/rtpthlp and up2/vp2 (the
    # Fortran recomputes it identically for coef_wpup2_wpvp2_..._zm).
    coef_impl_zm = ONE_THIRD * beta * a1 * wp3_on_wp2
    term_wprtp2_expl_zm = (1.0 - ONE_THIRD * beta) * a1 ** 2 \
        * wprtp ** 2 * wp3_on_wp2 / wp2
    term_wpthlp2_expl_zm = (1.0 - ONE_THIRD * beta) * a1 ** 2 \
        * wpthlp ** 2 * wp3_on_wp2 / wp2
    term_wprtpthlp_expl_zm = (1.0 - ONE_THIRD * beta) * a1 ** 2 \
        * wprtp * wpthlp * wp3_on_wp2 / wp2
    term_wpup2_expl_zm = (1.0 - ONE_THIRD * beta) * a1 ** 2 \
        * upwp ** 2 * wp3_on_wp2 / wp2
    term_wpvp2_expl_zm = (1.0 - ONE_THIRD * beta) * a1 ** 2 \
        * vpwp ** 2 * wp3_on_wp2 / wp2
    # sign of the turbulent velocity is wp3_on_wp2 itself for ADG1
    sgn = wp3_on_wp2

    kw2 = tunables["c_K2"] * kh_zt
    kw9 = tunables["c_K9"] * kh_zt

    # ----- rt'^2 -----
    lhs = _xp2_xpyp_lhs(gr, dt, coef_impl_zm, sgn, tau_zm, wm_zm, kw2,
                        rho_ds_zm, invrs_rho_ds_zm, c2rt_1d,
                        nu2_vert_res_dep, l_iter)
    rhs = _xp2_xpyp_rhs(gr, dt, coef_impl_zm, term_wprtp2_expl_zm, sgn,
                        wprtp, wprtp, rtm, rtm, rtp2, rtp2_forcing,
                        rho_ds_zm, invrs_rho_ds_zm, c2rt_1d, tau_zm,
                        RT_TOL ** 2, l_iter)
    rtp2, _ = tridag_solve(*lhs, rhs)

    # ----- th_l'^2 -----
    lhs = _xp2_xpyp_lhs(gr, dt, coef_impl_zm, sgn, tau_zm, wm_zm, kw2,
                        rho_ds_zm, invrs_rho_ds_zm, c2thl_1d,
                        nu2_vert_res_dep, l_iter)
    rhs = _xp2_xpyp_rhs(gr, dt, coef_impl_zm, term_wpthlp2_expl_zm, sgn,
                        wpthlp, wpthlp, thlm, thlm, thlp2,
                        thlp2_forcing, rho_ds_zm, invrs_rho_ds_zm,
                        c2thl_1d, tau_zm, THL_TOL ** 2, l_iter)
    thlp2, _ = tridag_solve(*lhs, rhs)

    # ----- r_t'th_l' -----
    lhs = _xp2_xpyp_lhs(gr, dt, coef_impl_zm, sgn, tau_zm, wm_zm, kw2,
                        rho_ds_zm, invrs_rho_ds_zm, c2rtthl_1d,
                        nu2_vert_res_dep, l_iter)
    rhs = _xp2_xpyp_rhs(gr, dt, coef_impl_zm, term_wprtpthlp_expl_zm,
                        sgn, wprtp, wpthlp, rtm, thlm, rtpthlp,
                        rtpthlp_forcing, rho_ds_zm, invrs_rho_ds_zm,
                        c2rtthl_1d, tau_zm, ZERO_THRESHOLD, l_iter)
    rtpthlp, _ = tridag_solve(*lhs, rhs)

    # ----- u'^2 / v'^2 (one matrix, two RHS) -----
    lhs = _xp2_xpyp_lhs(gr, dt, coef_impl_zm, sgn, tau_zm, wm_zm, kw9,
                        rho_ds_zm, invrs_rho_ds_zm, c4_c14_1d,
                        nu9_vert_res_dep, l_iter)
    rhs_u = _xp2_xpyp_uv_rhs(gr, dt, coef_impl_zm, term_wpup2_expl_zm,
                             sgn, wp2, wpthvp, tau_zm, um, vm, upwp,
                             vpwp, up2, vp2, rho_ds_zm,
                             invrs_rho_ds_zm, thv_ds_zm, c4, c5, c14,
                             c4_c14_1d, wp2_splat, l_iter)
    rhs_v = _xp2_xpyp_uv_rhs(gr, dt, coef_impl_zm, term_wpvp2_expl_zm,
                             sgn, wp2, wpthvp, tau_zm, vm, um, vpwp,
                             upwp, vp2, up2, rho_ds_zm,
                             invrs_rho_ds_zm, thv_ds_zm, c4, c5, c14,
                             c4_c14_1d, wp2_splat, l_iter)
    uv_solution, _ = tridag_solve(*lhs, jnp.stack([rhs_u, rhs_v],
                                                  axis=1))
    up2 = uv_solution[:, 0]
    vp2 = uv_solution[:, 1]

    # ----- positive-definite hole filling (l_hole_fill=T) -----
    if l_hole_fill:
        rtp2 = pos_definite_variances(gr, RT_TOL ** 2, rho_ds_zm, rtp2)
        thlp2 = pos_definite_variances(gr, THL_TOL ** 2, rho_ds_zm,
                                       thlp2)
        up2 = pos_definite_variances(gr, W_TOL_SQD, rho_ds_zm, up2)
        vp2 = pos_definite_variances(gr, W_TOL_SQD, rho_ds_zm, vp2)

    # ----- clipping (l_min_xp2_from_corr_wx=F) -----
    rtp2 = clip_variance(RT_TOL ** 2, rtp2)

    # l_clip_large_rtp2 (all levels, incl. top)
    rtp2 = jnp.where(rtp2 > RTP2_CLIP_COEF * rtm ** 2,
                     RTP2_CLIP_COEF * rtm ** 2, rtp2)

    thlp2 = clip_variance(THL_TOL ** 2, thlp2)

    up2 = clip_variance(W_TOL_SQD, up2)
    up2 = jnp.where(up2 > 1000.0, 1000.0, up2)

    vp2 = clip_variance(W_TOL_SQD, vp2)
    vp2 = jnp.where(vp2 > 1000.0, 1000.0, vp2)

    # up2/vp2 sponge damping: OFF in EAM (settings never assigned)

    # r_t'th_l' correlation clipping
    rtpthlp = clip_covar(MAX_MAG_CORRELATION, rtp2, thlp2, rtpthlp)

    return rtp2, thlp2, rtpthlp, up2, vp2
