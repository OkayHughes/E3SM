"""CLUBB slice E: advance_xm_wpxp — the coupled semi-implicit advance
of (rtm, wprtp) and (thlm, wpthlp).

PORT_NOTES
----------
Sources (components/eam/src/physics/clubb/):
  advance_xm_wpxp_module.F90  advance_xm_wpxp, xm_wpxp_lhs,
                              xm_wpxp_rhs, xm_wpxp_solve,
                              xm_wpxp_clipping_and_stats,
                              xm_term_ta_lhs_all, wpxp_term_tp_lhs_all,
                              wpxp_terms_ac_pr2_lhs_all,
                              wpxp_term_pr1_lhs_all,
                              wpxp_terms_bp_pr3_rhs_all,
                              damp_coefficient
  turbulent_adv_pdf.F90       xpyp_term_ta_pdf_lhs_all/rhs_all
                              (CENTERED branch — l_upwind_wpxp_ta=F,
                              unlike the xpyp upwind branch of slice C)
  mean_adv.F90                term_ma_zm_lhs_all      (slice C helper)
  diffusion.F90               diffusion_zm_lhs_all    (slice C helper)
  lapack_wrap.F90             band_solve -> DGBSV     (band_solver.py)
  mono_flux_limiter.F90       monotonic flux limiter  (mono_flux_limiter.py)
  fill_holes.F90              fill_holes_vertical "zt" path
  clip_explicit.F90           clip_covar (slice C helper; wprtp/wpthlp
                              use max_mag_correlation_flux)

EAMv3-active configuration (asserted against the Fortran module state
in harness/gen_clubb_golden.py gen_xm_wpxp via drv_xm_wpxp_config):
  l_iter = T (compile-time), l_clip_semi_implicit = F,
  l_explicit_turbulent_adv_wpxp = F and iiPDF_type = ADG1 (compile
  time) -> the turbulent advection of <w'x'> is ENTIRELY implicit:
  <w'^2 x'> = a1_zt * wp3_on_wp2_zt * <w'x'>, explicit term = 0;
  l_upwind_wpxp_ta = F -> CENTERED ta discretization (the *_zm coef /
  sgn_turbulent_vel arrays are never assigned or read);
  l_predict_upwp_vpwp = F -> um/upwp/vm/vpwp and all their inputs
  (um_forcing, ug/vg, wpthvp, fcor, uprcp/vprcp/rc_coef, up2/vp2,
  um/vm sponge + l_uv_nudge nudging) are DEAD;
  sclr_dim = 0 -> the wpsclrp branch bodies never execute (loops run
  zero times; NOT ported);
  l_pos_def = F -> pos_definite_adj is dead (pos_definite_module NOT
  ported); l_hole_fill = T; l_clip_turb_adv = F -> the
  xm_correction_wpxp_cl adjuster is dead (NOT ported);
  l_diffuse_rtm_and_thlm = F and l_stability_correct_Kh_N2_zm = F ->
  Kh_zm, em, Lscale (in the lhs), thlm/exner/rtm/rcm/p/thvm arguments
  of xm_wpxp_lhs are dead; rtm/thlm sponge damping settings are never
  assigned in an EAM build (static zero-init -> OFF);
  l_stats = F -> xm_wpxp_solve is called WITHOUT rcond, i.e.
  band_solve -> LAPACK DGBSV.  The dgbsvx condition-estimate /
  equilibration / iterative-refinement path (band_solvex) is dead.

Under this flag set the Fortran takes the "simple case" branch: ONE
LHS band matrix is built (with wpxp := zeros_vector, unused), the rtm
and thlm right-hand sides are solved together as a 2-RHS DGBSV call,
and xm_wpxp_clipping_and_stats runs per field:
  solution unpack -> xm(1)=xm(2) -> monotonic flux limiter ->
  [pos_def skipped] -> fill_holes_vertical(2, xm_threshold, "zt") when
  any(xm < xm_threshold) -> clip_covar with RELAXED variance
  (max(1e-7, rtp2) / max(0.01, thlp2), l_enable_relaxed_clipping=T,
  max_mag_correlation_flux) -> [xm correction skipped].

C6/C7 skewness functions: C6x_Skw_fnc = C6xb + (C6x-C6xb) *
exp(-Skw_zm^2/(2*C6xc^2)) when |C6x-C6xb| > |C6x+C6xb|*eps/2, else the
constant C6xb (under EAMv3 C7=C7b=0.5, so C7_Skw_fnc starts as the
constant 0.5 — the branch is decided on concrete params at trace
time, exactly like the Fortran's scalar if).  All three are then
linearly damped toward C6rt/C6thl/C7_Lscale0 where
Lscale < wpxp_L_thresh (= 100.0, the EAMv3 phys="default" namelist
value, see the slice-C correction) AND gr%zt > altitude_threshold
(= 100.0): damped = Lscale0 + ((C - Lscale0)/thresh) * Lscale.  The
debug-level-0 assertion 0 <= C7_Skw_fnc <= 1 cannot trip for the
EAMv3 parameter values (0.5..0.85 by construction) and is not ported.

Numerical-fidelity notes:
  * Fortran objects are -ffp-contract=off; every term assembly here
    reproduces the source's association order (see the helpers).
  * xm_wpxp_rhs multiplies by a precomputed invrs_dt = 1/dt
    (reciprocal-multiply), while mono_flux_limiter's mfl_xm_rhs
    divides by dt — both replicated (see _divs there).
  * Divisions by scalar tunables (Skw_zm/C6rtc) go through _divs to
    force true IEEE division under XLA.
  * The band solve is the exact reference-LAPACK DGBSV port
    (band_solver.py); the container liblapack's FMA contraction makes
    it a few ulps rather than bitwise (measured in the tests).

Returns (rtm, wprtp, thlm, wpthlp, diags) where diags carries the
mono-flux-limiter adjustment profiles and hole-fill gates (real
Fortran locals, exposed for Tier-0 and coverage tests only).
"""

import jax.numpy as jnp
import numpy as np
from jax import lax

from .advance_xp2_xpyp import (_diffusion_zm_lhs_interior,
                               _fill_holes_multiplicative,
                               _ma_zm_lhs_interior, clip_covar)
from .band_solver import band_solve
from .constants import (EPS, GAMMA_OVER_IMPLICIT_TS, GRAV,
                        MAX_MAG_CORRELATION_FLUX, RT_TOL, THL_TOL,
                        ZERO_THRESHOLD)
from .grid import zm2zt
from .mono_flux_limiter import (MONO_FLUX_RTM, MONO_FLUX_THLM,
                                calc_turb_adv_range,
                                monotonic_turbulent_flux_limit)

RT_TOL_MFL = 1.0e-4     # constants_clubb rt_tol_mfl
THL_TOL_MFL = 1.0e-2    # constants_clubb thl_tol_mfl
RELAXED_XP2_FLOOR = {MONO_FLUX_RTM: 1.0e-7, MONO_FLUX_THLM: 0.01}


def _divs(a, s):
    """Fortran-faithful division of an array by a scalar (XLA rewrites
    x / const into x * (1/const), which is not correctly rounded)."""
    return a / jnp.full_like(a, s)


# ---------------------------------------------------------------------------
# coefficient functions
# ---------------------------------------------------------------------------

def skw_coefficient(c, cb, cc, skw_zm):
    """C6x/C7 skewness function: cb + (c-cb)*exp(-(Skw/cc)^2/2), or the
    constant cb when c ~= cb (Fortran compute-time-saving branch on
    concrete scalars)."""
    if abs(c - cb) > abs(c + cb) * EPS / 2.0:
        return cb + (c - cb) * jnp.exp(-0.5 * _divs(skw_zm, cc) ** 2)
    return jnp.full_like(skw_zm, cb)


def damp_coefficient(coefficient, cx_skw_fnc, max_coeff_value,
                     threshold, lscale, zt, altitude_threshold):
    """Linear Lscale damping of C6/C7 in stably stratified regions
    (CLUBB ticket #431); active where Lscale < threshold
    (wpxp_L_thresh) and zt > altitude_threshold."""
    slope = (coefficient - max_coeff_value) / threshold  # scalars
    return jnp.where((lscale < threshold) & (zt > altitude_threshold),
                     max_coeff_value + slope * lscale, cx_skw_fnc)


# ---------------------------------------------------------------------------
# LHS building blocks (interior rows k = 1..nz-2 0-based unless noted)
# ---------------------------------------------------------------------------

def _ta_wpxp_centered_interior(gr, coef_zt, rho_ds_zt,
                               invrs_rho_ds_zm):
    """xpyp_term_ta_pdf_lhs_all, CENTERED branch (super, main, sub)."""
    w1 = gr.factor_zm2zt          # weights_zm2zt(1,:)
    w2 = 1.0 - gr.factor_zm2zt    # weights_zm2zt(2,:)
    inv = invrs_rho_ds_zm[1:-1]
    idzm = gr.invrs_dzm[1:-1]
    sup = inv * idzm * rho_ds_zt[2:] * coef_zt[2:] * w1[2:]
    dia = inv * idzm * (rho_ds_zt[2:] * coef_zt[2:] * w2[2:]
                        - rho_ds_zt[1:-1] * coef_zt[1:-1] * w1[1:-1])
    sub = -(inv * idzm * rho_ds_zt[1:-1] * coef_zt[1:-1] * w2[1:-1])
    return sup, dia, sub


def _ta_rhs_centered_interior(gr, term_zt, rho_ds_zt,
                              invrs_rho_ds_zm):
    """xpyp_term_ta_pdf_rhs_all, CENTERED branch (term_zt = 0 under
    ADG1 implicit ta; computed anyway for -0.0 parity)."""
    return -(invrs_rho_ds_zm[1:-1] * gr.invrs_dzm[1:-1]
             * (rho_ds_zt[2:] * term_zt[2:]
                - rho_ds_zt[1:-1] * term_zt[1:-1]))


def _build_lhs_rhs(gr, dt, wp2, coef_zt, kw6, tau_c6_zm,
                   c6rt_skw_fnc, c6thl_skw_fnc, c7_skw_fnc, wm_zm,
                   wm_zt, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,
                   invrs_rho_ds_zt, thv_ds_zm, nu6, rtm, wprtp,
                   rtm_forcing, wprtp_forcing, rtpthvp, thlm, wpthlp,
                   thlm_forcing, wpthlp_forcing, thlpthvp):
    """xm_wpxp_lhs (once) + xm_wpxp_rhs for rtm (RHS 1) and thlm
    (RHS 2), interleaved band storage: CLUBB lhs (5, 2*nz), 0-based
    column 2k = xm(k), 2k+1 = wpxp(k)."""
    nz = wp2.shape[0]
    n2 = 2 * nz
    invrs_dt = 1.0 / dt
    dtype = wp2.dtype

    # --- shared interior term rows ---
    d_sup, d_dia, d_sub = _diffusion_zm_lhs_interior(gr, kw6, nu6)
    m_sup, m_dia, m_sub = _ma_zm_lhs_interior(gr, wm_zm)
    # ta of xm, Fortran k = 2..nz (0-based 1..nz-1)
    ta_xm1 = invrs_rho_ds_zt[1:] * gr.invrs_dzt[1:] * rho_ds_zm[1:]
    ta_xm2 = -(invrs_rho_ds_zt[1:] * gr.invrs_dzt[1:]
               * rho_ds_zm[:-1])
    # tp / pr1 / ac_pr2, Fortran k = 2..nz-1
    tp1 = wp2[1:-1] * gr.invrs_dzm[1:-1]
    tp2 = -(wp2[1:-1] * gr.invrs_dzm[1:-1])
    pr1_rt = c6rt_skw_fnc[1:-1] / tau_c6_zm[1:-1]
    pr1_thl = c6thl_skw_fnc[1:-1] / tau_c6_zm[1:-1]
    ac_pr2 = (1.0 - c7_skw_fnc[1:-1]) * gr.invrs_dzm[1:-1] \
        * (wm_zt[2:] - wm_zt[1:-1])
    t_sup, t_dia, t_sub = _ta_wpxp_centered_interior(
        gr, coef_zt, rho_ds_zt, invrs_rho_ds_zm)

    # --- LHS band matrix (single matrix; pr1 uses C6rt as in the
    # Fortran simple-case call with C6rt_Skw_fnc) ---
    cols_xm = 2 * np.arange(1, nz)          # xm(k), k = 1..nz-1
    cols_wp = 2 * np.arange(1, nz - 1) + 1  # wpxp(k), k = 1..nz-2
    lhs = jnp.zeros((5, n2), dtype=dtype)
    lhs = lhs.at[1, cols_xm].set(ta_xm1)
    lhs = lhs.at[2, cols_xm].set(invrs_dt)
    lhs = lhs.at[3, cols_xm].set(ta_xm2)
    lhs = lhs.at[0, cols_wp].set(
        m_sup + d_sup + GAMMA_OVER_IMPLICIT_TS * t_sup)
    lhs = lhs.at[1, cols_wp].set(tp1)
    lhs = lhs.at[2, cols_wp].set(
        (m_dia + d_dia + ac_pr2)
        + GAMMA_OVER_IMPLICIT_TS * (t_dia + pr1_rt) + invrs_dt)
    lhs = lhs.at[3, cols_wp].set(tp2)
    lhs = lhs.at[4, cols_wp].set(
        (m_sub + d_sub) + GAMMA_OVER_IMPLICIT_TS * t_sub)
    # boundary identity rows: xm(1), wpxp(1), wpxp(nz)
    lhs = lhs.at[:, 0].set(0.0).at[2, 0].set(1.0)
    lhs = lhs.at[:, 1].set(0.0).at[2, 1].set(1.0)
    lhs = lhs.at[:, n2 - 1].set(0.0).at[2, n2 - 1].set(1.0)

    # --- RHS vectors ---
    zeros_t = jnp.zeros(nz, dtype=dtype)
    rhs_ta = _ta_rhs_centered_interior(gr, zeros_t, rho_ds_zt,
                                       invrs_rho_ds_zm)

    def one_rhs(xm, wpxp, xm_forcing, wpxp_forcing, xpthvp, pr1):
        bp_pr3 = (GRAV / thv_ds_zm[1:-1]) \
            * (1.0 - c7_skw_fnc[1:-1]) * xpthvp[1:-1]
        inner = (((-(t_sup * wpxp[2:]) - t_dia * wpxp[1:-1])
                  - t_sub * wpxp[:-2]) - pr1 * wpxp[1:-1])
        rhs_wp = (bp_pr3 + wpxp_forcing[1:-1] + rhs_ta
                  + (1.0 - GAMMA_OVER_IMPLICIT_TS) * inner) \
            + wpxp[1:-1] * invrs_dt
        rhs_xm = xm[1:-1] * invrs_dt + xm_forcing[1:-1]
        rhs = jnp.zeros(n2, dtype=dtype)
        rhs = rhs.at[0].set(xm[0]).at[1].set(wpxp[0])
        rhs = rhs.at[jnp.asarray(cols_xm[:-1])].set(rhs_xm)
        rhs = rhs.at[jnp.asarray(cols_wp)].set(rhs_wp)
        rhs = rhs.at[n2 - 2].set(xm[nz - 1] * invrs_dt
                                 + xm_forcing[nz - 1])
        return rhs   # rhs[n2-1] stays 0 (upper wpxp boundary)

    rhs_rt = one_rhs(rtm, wprtp, rtm_forcing, wprtp_forcing, rtpthvp,
                     pr1_rt)
    rhs_thl = one_rhs(thlm, wpthlp, thlm_forcing, wpthlp_forcing,
                      thlpthvp, pr1_thl)
    return lhs, jnp.stack([rhs_rt, rhs_thl], axis=1)


# ---------------------------------------------------------------------------
# fill_holes_vertical, "zt" path
# ---------------------------------------------------------------------------

def fill_holes_vertical_zt(gr, num_draw_pts, threshold, rho_ds_zt,
                           field):
    """fill_holes_vertical for a zt-grid field: upper_hf_level = nz
    (the top thermo level IS fillable, unlike zm), density rho_ds_zt,
    thickness gr.dzt.  Sequential windows then a global pass, each
    exactly gated (an untriggered window is a bitwise no-op)."""
    nz = field.shape[0]
    n = num_draw_pts
    width = 2 * n + 1
    rho_ds_zt = jnp.asarray(rho_ds_zt, dtype=field.dtype)
    dzt = gr.dzt

    def window_pass(kk, field):
        lo = kk - n
        window = lax.dynamic_slice(field, (lo,), (width,))
        gate = jnp.any(window < threshold)
        filled = _fill_holes_multiplicative(
            threshold, lax.dynamic_slice(rho_ds_zt, (lo,), (width,)),
            lax.dynamic_slice(dzt, (lo,), (width,)), window)
        return lax.dynamic_update_slice(
            field, jnp.where(gate, filled, window), (lo,))

    # Fortran k = 2+n..nz-n (0-based centers n+1..nz-n-2), sequential
    field = lax.fori_loop(n + 1, nz - n - 1, window_pass, field)
    window = field[1:nz]
    gate = jnp.any(window < threshold)
    filled = _fill_holes_multiplicative(
        threshold, rho_ds_zt[1:nz], dzt[1:nz], window)
    field = field.at[1:nz].set(jnp.where(gate, filled, window))
    return field


# ---------------------------------------------------------------------------
# per-field clipping-and-stats sequence
# ---------------------------------------------------------------------------

def _clipping_sequence(gr, solve_type, dt, wp2, xp2, xm_forcing,
                       rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,
                       invrs_rho_ds_zt, xp2_threshold, xm_threshold,
                       xm_tol, low_lev, high_lev, xm_old, sol_xm,
                       sol_wpxp):
    """xm_wpxp_clipping_and_stats under EAMv3 (see module docstring):
    mono flux limiter -> hole filling -> relaxed covariance clip."""
    xm = sol_xm.at[0].set(sol_xm[1])
    wpxp = sol_wpxp

    xm, wpxp, wpxp_net_adjust = monotonic_turbulent_flux_limit(
        gr, solve_type, dt, xm_old, xp2, xm_forcing, rho_ds_zm,
        rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, xp2_threshold,
        low_lev, high_lev, xm, xm_tol, wpxp)

    # l_pos_def = F: no pos_definite_adj

    hole_gate = jnp.any(xm < xm_threshold)   # call-site gate
    filled = fill_holes_vertical_zt(gr, 2, xm_threshold, rho_ds_zt,
                                    xm)
    xm = jnp.where(hole_gate, filled, xm)

    # relaxed covariance clipping (l_enable_relaxed_clipping = T)
    xp2_relaxed = jnp.maximum(RELAXED_XP2_FLOOR[solve_type], xp2)
    wpxp = clip_covar(MAX_MAG_CORRELATION_FLUX, wp2, xp2_relaxed,
                      wpxp)

    # l_clip_turb_adv = F: no xm_correction_wpxp_cl
    return xm, wpxp, wpxp_net_adjust, hole_gate


# ---------------------------------------------------------------------------
# main routine
# ---------------------------------------------------------------------------

def advance_xm_wpxp(gr, dt, sigma_sqd_w, wm_zm, wm_zt, wp2, lscale,
                    wp3_on_wp2_zt, kh_zt, tau_c6_zm, skw_zm, rtpthvp,
                    rtm_forcing, wprtp_forcing, thlpthvp,
                    thlm_forcing, wpthlp_forcing, rho_ds_zm,
                    rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
                    thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm,
                    varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm, rtm,
                    wprtp, thlm, wpthlp, nu6_vert_res_dep, tunables):
    """advance_xm_wpxp under the EAMv3 configuration.

    DEAD Fortran arguments dropped from this signature (see module
    docstring): em, wp3_on_wp2 (zm), Kh_zm, wp2rtp, wp2thlp, exner,
    rcm, p_in_Pa, thvm, Cx_fnc_Richardson, rtm_ref/thlm_ref, all
    um/vm-prediction inputs, all sclr arrays,
    pdf_implicit_coefs_terms (iiPDF_new only).

    tunables: mapping with C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
    C7, C7b, C7c, c_K6, C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0,
    wpxp_L_thresh, altitude_threshold (floats from the EAMv3 params
    vector; wpxp_L_thresh = 100.0).

    Returns (rtm, wprtp, thlm, wpthlp, diags)."""
    # ----- C6/C7 skewness functions + Lscale damping -----
    c6rt_skw = skw_coefficient(tunables["C6rt"], tunables["C6rtb"],
                               tunables["C6rtc"], skw_zm)
    c6thl_skw = skw_coefficient(tunables["C6thl"], tunables["C6thlb"],
                                tunables["C6thlc"], skw_zm)
    c7_skw = skw_coefficient(tunables["C7"], tunables["C7b"],
                             tunables["C7c"], skw_zm)
    c7_skw = damp_coefficient(
        tunables["C7"], c7_skw, tunables["C7_Lscale0"],
        tunables["wpxp_L_thresh"], lscale, gr.zt,
        tunables["altitude_threshold"])
    c6rt_skw = damp_coefficient(
        tunables["C6rt"], c6rt_skw, tunables["C6rt_Lscale0"],
        tunables["wpxp_L_thresh"], lscale, gr.zt,
        tunables["altitude_threshold"])
    c6thl_skw = damp_coefficient(
        tunables["C6thl"], c6thl_skw, tunables["C6thl_Lscale0"],
        tunables["wpxp_L_thresh"], lscale, gr.zt,
        tunables["altitude_threshold"])

    # ----- eddy diffusivity for w'x' -----
    kw6 = tunables["c_K6"] * kh_zt

    # ----- monotonic-flux-limiter level range -----
    low_lev, high_lev = calc_turb_adv_range(
        gr, dt, w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm,
        mixt_frac_zm)

    # ----- ADG1 implicit turbulent-advection coefficient -----
    a1 = 1.0 / (1.0 - sigma_sqd_w)
    a1_zt = jnp.maximum(zm2zt(gr, a1), ZERO_THRESHOLD)
    coef_zt = a1_zt * wp3_on_wp2_zt

    # ----- one LHS, two RHS, one DGBSV solve -----
    lhs, rhs = _build_lhs_rhs(
        gr, dt, wp2, coef_zt, kw6, tau_c6_zm, c6rt_skw, c6thl_skw,
        c7_skw, wm_zm, wm_zt, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,
        invrs_rho_ds_zt, thv_ds_zm, nu6_vert_res_dep, rtm, wprtp,
        rtm_forcing, wprtp_forcing, rtpthvp, thlm, wpthlp,
        thlm_forcing, wpthlp_forcing, thlpthvp)
    solution, singular = band_solve(lhs, rhs)

    # ----- per-field clipping-and-stats -----
    rtm_new, wprtp_new, adj_rt, hole_rt = _clipping_sequence(
        gr, MONO_FLUX_RTM, dt, wp2, rtp2, rtm_forcing, rho_ds_zm,
        rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, RT_TOL ** 2,
        RT_TOL, RT_TOL_MFL, low_lev, high_lev, rtm,
        solution[0::2, 0], solution[1::2, 0])
    thlm_new, wpthlp_new, adj_thl, hole_thl = _clipping_sequence(
        gr, MONO_FLUX_THLM, dt, wp2, thlp2, thlm_forcing, rho_ds_zm,
        rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, THL_TOL ** 2,
        THL_TOL, THL_TOL_MFL, low_lev, high_lev, thlm,
        solution[0::2, 1], solution[1::2, 1])

    diags = dict(singular=singular, wprtp_mfl_adjust=adj_rt,
                 wpthlp_mfl_adjust=adj_thl, rtm_hole_fill=hole_rt,
                 thlm_hole_fill=hole_thl, low_lev_effect=low_lev,
                 high_lev_effect=high_lev,
                 c6rt_skw_fnc=c6rt_skw, c6thl_skw_fnc=c6thl_skw,
                 c7_skw_fnc=c7_skw)
    return rtm_new, wprtp_new, thlm_new, wpthlp_new, diags
