"""CLUBB monotonic turbulent flux limiter — port of
mono_flux_limiter.F90 (calc_turb_adv_range, mean_vert_vel_up_down,
monotonic_turbulent_flux_limit + mfl_xm_lhs/rhs/solve).

PORT_NOTES
----------
Source: components/eam/src/physics/clubb/mono_flux_limiter.F90.

EAMv3-active configuration (asserted against the Fortran module state
in harness/gen_clubb_golden.py gen_xm_wpxp):
  l_mfl_xm_imp_adj     = T (compile-time): when the limiter engages,
                         xm is re-solved semi-implicitly with a
                         tridiagonal system (slice A's exact dgtsv).
  l_constant_thickness = F (compile-time in calc_turb_adv_range): the
                         level range is found by walking grid levels
                         until the travel time (dz / mean one-signed
                         w) accumulates past dt.
  l_implemented        = T: the mfl_xm_lhs mean-advection term is
                         zero; the LHS is 1/dt on the diagonal.
  solve_type um/vm (l_predict_upwp_vpwp=T only) NOT ported: EAMv3 only
  reaches the rtm (max_xp2=5e-6, positive-definite floor) and thlm
  (max_xp2=5.0, positive-definite floor) cases.

Structure notes:
  * calc_turb_adv_range / mean_vert_vel_up_down are implemented in
    NUMPY with Python control flow — an exact transcription of the
    sequential per-level walks (data-dependent early exits).  They
    return Fortran 1-based integer level indices, exactly as recorded
    in the golden (the limiter converts internally).  This mirrors the
    Fortran: the outputs are integers, used only for index windows.
    mean_vert_vel_up_down leaves levels 1 and nz unset in the Fortran
    (never read by the walks); the port zeroes them.
  * The Fortran walk accumulates dt_all_grid_levs sequentially from
    the starting level outward; the port replays the identical
    accumulation order (no prefix-sum reassociation).
  * monotonic_turbulent_flux_limit's main loop is SEQUENTIAL in k: the
    bounds at level k use wpxp(k-1) *after* any adjustment made at
    level k-1.  Ported as a Python loop over the static nz with
    jnp.where gates.
  * The engagement gate is any(|wpxp_net_adjust| > eps) with
    eps = 1e-10 — wpxp itself is already reset to the bound even when
    the adjustment is below eps (no xm re-solve then).  Replicated.
  * The re-solved xm comes from tridag_solve (LAPACK dgtsv, slice A
    port) on a diagonal-only system (l_implemented=T), then
    xm(1) = xm(2).
  * The domain-top spike check |xm(nz)-xm_enter(nz)| > 10*xm_tol uses
    vertical_integral = SUM(field*rho*dz) over Fortran levels
    2..nz-1; the port accumulates strictly left-to-right (gfortran SUM
    loop order).  The xm_vert_integral < eps warning path (integral
    ~ 0) keeps xm conservation off exactly as the Fortran does.

Returns (xm, wpxp, wpxp_net_adjust); the adjustment profile is a real
Fortran local exposed for Tier-0/coverage tests only.
"""

import numpy as np
import jax
import jax.numpy as jnp
from jax.scipy.special import erf as _jerf

from .constants import EPS, SQRT_2, SQRT_2PI, ZERO_THRESHOLD
from .grid import zm2zt
from .tridiag import tridag_solve

# solve-type constants (values match advance_xm_wpxp_module /
# mono_flux_limiter named constants)
MONO_FLUX_THLM = 1
MONO_FLUX_RTM = 2

MAX_XP2 = {MONO_FLUX_RTM: 5.0e-6, MONO_FLUX_THLM: 5.0}


def _divs(a, s):
    """Fortran-faithful division of an array by a scalar (XLA rewrites
    x / const into x * (1/const), which is not correctly rounded)."""
    return a / jnp.full_like(a, s)


def mean_vert_vel_up_down(w_1_zm, w_2_zm, varnce_w_1_zm,
                          varnce_w_2_zm, mixt_frac_zm, w_ref=0.0):
    """Average one-signed vertical velocities of the two-component w
    PDF (numpy; interior levels 2..nz-1 Fortran, boundaries zeroed)."""
    w1 = np.asarray(w_1_zm, dtype=np.float64)
    w2 = np.asarray(w_2_zm, dtype=np.float64)
    s1 = np.sqrt(np.asarray(varnce_w_1_zm, dtype=np.float64))
    s2 = np.sqrt(np.asarray(varnce_w_2_zm, dtype=np.float64))
    mf = np.asarray(mixt_frac_zm, dtype=np.float64)
    nz = w1.shape[0]

    def comp(w_np, s_np):
        # exp/erf evaluated through jnp (XLA libm), like every other
        # transcendental in this port; results pulled back to numpy
        w = jnp.asarray(w_np)
        s = jnp.asarray(s_np)
        all_down = w + 3.0 * s <= w_ref
        all_up = w - 3.0 * s >= w_ref
        # straddling branch (guard the divisions; masked out elsewhere)
        s_safe = jnp.where(s > 0.0, s, 1.0)
        exp_cache = jnp.exp(-(w_ref - w) ** 2 / (2.0 * s_safe ** 2))
        erf_cache = _jerf((w_ref - w) / (SQRT_2 * s_safe))
        down_c = -(s / SQRT_2PI) * exp_cache \
            + w * 0.5 * (1.0 + erf_cache)
        up_c = +(s / SQRT_2PI) * exp_cache \
            + w * 0.5 * (1.0 - erf_cache)
        down = jnp.where(all_down, w, jnp.where(all_up, 0.0, down_c))
        up = jnp.where(all_down, 0.0, jnp.where(all_up, w, up_c))
        return np.asarray(down), np.asarray(up)

    d1, u1 = comp(w1, s1)
    d2, u2 = comp(w2, s2)
    mean_w_down = mf * d1 + (1.0 - mf) * d2
    mean_w_up = mf * u1 + (1.0 - mf) * u2
    # Fortran computes only k = 2..nz-1; levels 1 and nz are unset
    # (and never read by calc_turb_adv_range)
    mean_w_down[0] = mean_w_down[-1] = 0.0
    mean_w_up[0] = mean_w_up[-1] = 0.0
    return mean_w_down, mean_w_up


def calc_turb_adv_range(gr, dt, w_1_zm, w_2_zm, varnce_w_1_zm,
                        varnce_w_2_zm, mixt_frac_zm):
    """Lowest/highest thermodynamic level able to affect level k by
    turbulent advection within one dt (l_constant_thickness=F path).
    Returns FORTRAN 1-based integer index arrays (low, high)."""
    invrs_dzm = np.asarray(gr.invrs_dzm, dtype=np.float64)
    nz = invrs_dzm.shape[0]
    vert_vel_down, vert_vel_up = mean_vert_vel_up_down(
        w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm,
        0.0)
    dt = float(dt)

    low = np.zeros(nz, dtype=np.int64)
    high = np.zeros(nz, dtype=np.int64)
    for k in range(3, nz - 1):          # Fortran k = 3..nz-2
        # downward walk (upward velocities), Fortran j = k-1, ...
        j = k - 1
        acc = 0.0
        while True:
            if vert_vel_up[j - 1] > 0.0:
                acc = acc + (1.0 / invrs_dzm[j - 1]) / vert_vel_up[j - 1]
                if acc >= dt:
                    low[k - 1] = j
                    break
                elif j == 2:
                    low[k - 1] = j
                    break
                else:
                    j -= 1
            else:
                low[k - 1] = j + 1
                break
        # upward walk (downward velocities), Fortran j = k+1, ...
        j = k + 1
        acc = 0.0
        while True:
            if vert_vel_down[j - 2] < 0.0:
                acc = acc \
                    + (-(1.0 / invrs_dzm[j - 2])) / vert_vel_down[j - 2]
                if acc >= dt:
                    high[k - 1] = j
                    break
                elif j == nz:
                    high[k - 1] = j
                    break
                else:
                    j += 1
            else:
                high[k - 1] = j - 1
                break

    # boundary bookkeeping (Fortran end of calc_turb_adv_range)
    low[0], high[0] = 1, 1
    low[1], high[1] = 2, 2
    low[nz - 2], high[nz - 2] = nz - 1, nz
    low[nz - 1], high[nz - 1] = nz, nz
    return low, high


def _seq_sum(x):
    """Strict left-to-right accumulation (gfortran SUM loop order)."""
    total, _ = jax.lax.scan(lambda c, v: (c + v, None),
                            jnp.zeros((), x.dtype), x)
    return total


def monotonic_turbulent_flux_limit(gr, solve_type, dt, xm_old, xp2,
                                   xm_forcing, rho_ds_zm, rho_ds_zt,
                                   invrs_rho_ds_zm, invrs_rho_ds_zt,
                                   xp2_threshold, low_lev_effect,
                                   high_lev_effect, xm, xm_tol, wpxp):
    """monotonic_turbulent_flux_limit under EAMv3 (l_implemented=T,
    l_mfl_xm_imp_adj=T).  low/high_lev_effect are the Fortran 1-based
    arrays from calc_turb_adv_range.  wm_zt is a real argument of the
    Fortran routine but its mean-advection contribution is hard-coded
    OFF (m_adv_term = 0), so it is dropped here.

    Returns (xm, wpxp, wpxp_net_adjust)."""
    xm = jnp.asarray(xm, dtype=jnp.float64)
    wpxp = jnp.asarray(wpxp, dtype=jnp.float64)
    xm_old = jnp.asarray(xm_old, dtype=jnp.float64)
    nz = xm.shape[0]
    max_xp2 = MAX_XP2[solve_type]
    l_positive = True   # rtm/thlm/scalars (um/vm not ported)

    xm_enter_mfl = xm

    xp2_zt = jnp.maximum(zm2zt(gr, xp2), xp2_threshold)
    xp2_zt = jnp.minimum(xp2_zt, max_xp2)

    # per-level usable min/max (vectorized over Fortran k = 2..nz)
    stnd_dev_x = jnp.sqrt(xp2_zt)
    max_dev = jnp.maximum(2.0 * stnd_dev_x, xm_tol)
    xm_without_ta = xm_old + dt * xm_forcing        # m_adv_term = 0
    if l_positive:
        min_lev = jnp.maximum(xm_without_ta - max_dev, ZERO_THRESHOLD)
    else:
        min_lev = xm_without_ta - max_dev
    max_lev = xm_without_ta + max_dev
    # boundary level 1 (Fortran)
    xm_without_ta = xm_without_ta.at[0].set(xm[0])
    min_lev = min_lev.at[0].set(min_lev[1])
    max_lev = max_lev.at[0].set(max_lev[1])

    # sequential limiter sweep, Fortran k = 2..nz-1 (wpxp(k-1) is read
    # AFTER any adjustment made at the previous level) — lax.scan with
    # the full wpxp array as carry
    idx = jnp.arange(nz)
    low_j = jnp.asarray(low_lev_effect, dtype=jnp.int32)
    high_j = jnp.asarray(high_lev_effect, dtype=jnp.int32)
    invrs_rho_ds_zm = jnp.asarray(invrs_rho_ds_zm, dtype=xm.dtype)
    rho_ds_zm = jnp.asarray(rho_ds_zm, dtype=xm.dtype)
    rho_ds_zt = jnp.asarray(rho_ds_zt, dtype=xm.dtype)
    invrs_rho_ds_zt = jnp.asarray(invrs_rho_ds_zt, dtype=xm.dtype)
    xm_forcing = jnp.asarray(xm_forcing, dtype=xm.dtype)

    def sweep(carry, k):
        wpxp_c, adjust_c = carry
        lo0 = jnp.maximum(low_j[k], 2) - 1
        hi0 = jnp.minimum(high_j[k], nz) - 1
        sel = (idx >= lo0) & (idx <= hi0)
        min_x = jnp.min(jnp.where(sel, min_lev, jnp.inf))
        max_x = jnp.max(jnp.where(sel, max_lev, -jnp.inf))

        wpxp_mfl_max = invrs_rho_ds_zm[k] \
            * ((rho_ds_zt[k] / (dt * gr.invrs_dzt[k]))
               * (xm_without_ta[k] - min_x)
               + rho_ds_zm[k - 1] * wpxp_c[k - 1])
        wpxp_mfl_min = invrs_rho_ds_zm[k] \
            * ((rho_ds_zt[k] / (dt * gr.invrs_dzt[k]))
               * (xm_without_ta[k] - max_x)
               + rho_ds_zm[k - 1] * wpxp_c[k - 1])

        too_big = wpxp_c[k] > wpxp_mfl_max
        too_small = (~too_big) & (wpxp_c[k] < wpxp_mfl_min)
        adj = jnp.where(too_big, wpxp_mfl_max - wpxp_c[k],
                        jnp.where(too_small,
                                  wpxp_mfl_min - wpxp_c[k], 0.0))
        new = jnp.where(too_big, wpxp_mfl_max,
                        jnp.where(too_small, wpxp_mfl_min,
                                  wpxp_c[k]))
        return (wpxp_c.at[k].set(new), adjust_c.at[k].set(adj)), None

    (wpxp, wpxp_net_adjust), _ = jax.lax.scan(
        sweep, (wpxp, jnp.zeros(nz, dtype=xm.dtype)),
        jnp.arange(1, nz - 1))

    engaged = jnp.any(jnp.abs(wpxp_net_adjust) > EPS)

    # ---- semi-implicit xm re-solve (l_mfl_xm_imp_adj = T) ----------
    # mfl_xm_lhs: l_implemented=T -> no mean advection; diagonal 1/dt
    # on Fortran rows 2..nz, identity row 1.
    one_over_dt = jnp.zeros(nz, dtype=xm.dtype) + (1.0 / dt)
    diag = jnp.concatenate([jnp.ones((1,), dtype=xm.dtype),
                            one_over_dt[1:]])
    supd = jnp.zeros(nz, dtype=xm.dtype)
    subd = jnp.zeros(nz, dtype=xm.dtype)
    # mfl_xm_rhs (note: xm_old(k) / dt is a TRUE division here, unlike
    # the reciprocal-multiplies in xm_wpxp_rhs -> _divs)
    rhs_int = (_divs(xm_old[1:], dt)
               - invrs_rho_ds_zt[1:] * gr.invrs_dzt[1:]
               * (rho_ds_zm[1:] * wpxp[1:]
                  - rho_ds_zm[:-1] * wpxp[:-1])
               + xm_forcing[1:])
    rhs = jnp.concatenate([xm_old[:1], rhs_int])
    xm_solved, _sing = tridag_solve(supd, diag, subd, rhs)
    xm_solved = xm_solved.at[0].set(xm_solved[1])

    # ---- domain-top spike removal ----------------------------------
    spike = jnp.abs(xm_solved[nz - 1] - xm_enter_mfl[nz - 1]) \
        > 10.0 * xm_tol
    dz_top = gr.zm[nz - 1] - gr.zm[nz - 2]
    xm_density_weighted = rho_ds_zt[nz - 1] \
        * (xm_solved[nz - 1] - xm_enter_mfl[nz - 1]) * dz_top
    xm_vert_integral = _seq_sum(
        xm_solved[1:nz - 1] * rho_ds_zt[1:nz - 1] * gr.dzt[1:nz - 1])
    degenerate = xm_vert_integral < EPS
    xm_adj_coef = xm_density_weighted \
        / jnp.where(degenerate, 1.0, xm_vert_integral)
    xm_adj_coef = jnp.where(xm_adj_coef < -0.99, -0.99, xm_adj_coef)
    xm_scaled = xm_solved * (1.0 + xm_adj_coef)
    xm_spiked = jnp.where(degenerate, xm_solved, xm_scaled)
    xm_spiked = xm_spiked.at[nz - 1].set(xm_enter_mfl[nz - 1])
    xm_after = jnp.where(spike, xm_spiked, xm_solved)

    xm = jnp.where(engaged, xm_after, xm)

    return xm, wpxp, wpxp_net_adjust
