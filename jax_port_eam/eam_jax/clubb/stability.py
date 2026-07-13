"""CLUBB stability helpers — port of advance_helper_module.F90:
calc_brunt_vaisala_freq_sqd, calc_stability_correction,
term_wp2_splat, term_wp3_splat.

PORT_NOTES
----------
Source: components/eam/src/physics/clubb/advance_helper_module.F90.

calc_brunt_vaisala_freq_sqd — all three formula variants are ported
(the flags are module variables, toggled in the harness to golden
each):
  - l_brunt_vaisala_freq_moist=F, l_use_thvm_in_bv_freq=F (the EAMv3
    model_flags defaults; clubb_intr never touches these):
        N^2 = (grav / T0) * ddzt(thlm)          [T0 = theta0 = 300 K]
  - l_use_thvm_in_bv_freq=T: N^2 = (grav / zt2zm(thvm)) * ddzt(thvm)
  - l_brunt_vaisala_freq_moist=T: the Durran & Klemp (1982) Eq. 36
    saturated formula, applied at EVERY level (the routine has no
    in-cloud/clear conditional; cloud enters through rcm via
    T = thlm*exner + Lv*rcm/Cp and thm = thlm + Lv/(Cp*exner)*rcm).

calc_stability_correction (ACTIVE in EAMv3: l_stability_correct_tau_zm
= .true., the model_flags default; EAMv3's clubb_stabcorrect=F
namelist flag only controls the UNRELATED l_stability_correct_Kh_N2_zm
/ l_diffuse_rtm_and_thlm):
    lambda0 = lambda0_stability_coef where N^2 > 0 else 0
    correction = 1 + min(lambda0 * N^2 * zt2zm(Lscale)^2 / em, 3)
advance_clubb_core then uses tau_N2_zm = tau_zm / correction as
tau_C1_zm and tau_C6_zm.

term_wp2_splat / term_wp3_splat — the eddy "splatting" tendencies:
    wp2_splat = -wp2 * min(5/dt, C_wp2_splat * tau_zm *
                            ddzt(sqrt(wp2_zt))^2)
    wp3_splat = -wp3 * min(5/dt, 3 * C_wp2_splat * tau_zt *
                            ddzm(sqrt(wp2))^2)
EAMv3 C_wp2_splat = 0.0 (CLUBB default, no namelist override), which
makes both tendencies exactly (-0.0); ported verbatim anyway (the
goldens also cover C_wp2_splat = 2 incl. the 5/dt clip arm).

DEAD in EAM (documented, NOT ported): compute_Cx_fnc_Richardson.  It
is only called by advance_clubb_core when l_use_C7_Richardson,
l_use_C11_Richardson or l_use_wp3_pr3 is true; all three are false in
the EAMv3 configuration (model_flags defaults, never overridden), so
advance_clubb_core sets Cx_fnc_Richardson = 0.0 (the lscale_tau
segment port does the same).  Its helper Lscale_width_vert_avg has no
other caller.
"""

import jax.numpy as jnp

from .constants import CP, EP, GRAV, LV, RD
from .grid import Grid, ddzm, ddzt, zt2zm
from .saturation import sat_mixrat_liq


def calc_brunt_vaisala_freq_sqd(gr: Grid, thlm, exner, rtm, rcm,
                                p_in_pa, thvm, t0,
                                l_brunt_vaisala_freq_moist=False,
                                l_use_thvm_in_bv_freq=False):
    """Brunt-Vaisala frequency squared N^2 on momentum levels."""
    if not l_brunt_vaisala_freq_moist:
        if l_use_thvm_in_bv_freq:
            thvm_zm = zt2zm(gr, thvm)
            ddzt_thvm = ddzt(gr, thvm)
            return (GRAV / thvm_zm) * ddzt_thvm
        ddzt_thlm = ddzt(gr, thlm)
        return (GRAV / t0) * ddzt_thlm

    # moist (Durran & Klemp 1982, Eq. 36)
    t_in_k = thlm * exner + LV * rcm / CP          # thlm2T_in_K
    t_in_k_zm = zt2zm(gr, t_in_k)
    rsat = sat_mixrat_liq(p_in_pa, t_in_k)
    rsat_zm = zt2zm(gr, rsat)
    ddzt_rsat = ddzt(gr, rsat)
    thm = thlm + LV / (CP * exner) * rcm
    thm_zm = zt2zm(gr, thm)
    ddzt_thm = ddzt(gr, thm)
    ddzt_rtm = ddzt(gr, rtm)
    return GRAV * (((1.0 + LV * rsat_zm / (RD * t_in_k_zm))
                    / (1.0 + EP * (LV ** 2) * rsat_zm
                       / (CP * RD * t_in_k_zm ** 2)))
                   * ((1.0 / thm_zm * ddzt_thm)
                      + (LV / (CP * t_in_k_zm)) * ddzt_rsat)
                   - ddzt_rtm)


def calc_stability_correction(gr: Grid, thlm, lscale, em, exner, rtm,
                              rcm, p_in_pa, thvm,
                              lambda0_stability_coef, t0,
                              l_brunt_vaisala_freq_moist=False,
                              l_use_thvm_in_bv_freq=False):
    """Stability correction factor (momentum levels), >= 1."""
    bv_sqd = calc_brunt_vaisala_freq_sqd(
        gr, thlm, exner, rtm, rcm, p_in_pa, thvm, t0,
        l_brunt_vaisala_freq_moist, l_use_thvm_in_bv_freq)
    lambda0_stability = jnp.where(bv_sqd > 0.0,
                                  lambda0_stability_coef, 0.0)
    return 1.0 + jnp.minimum(
        lambda0_stability * bv_sqd * zt2zm(gr, lscale) ** 2 / em, 3.0)


def term_wp2_splat(gr: Grid, c_wp2_splat, dt, wp2, wp2_zt, tau_zm):
    """Negative wp2 tendency from eddy splatting (momentum levels)."""
    d_sqrt_wp2_dz = ddzt(gr, jnp.sqrt(wp2_zt))
    return -wp2 * jnp.minimum(
        5.0 / dt, c_wp2_splat * tau_zm * d_sqrt_wp2_dz ** 2)


def term_wp3_splat(gr: Grid, c_wp2_splat, dt, wp2, wp3, tau_zt):
    """wp3 damping from eddy splatting (thermodynamic levels)."""
    d_sqrt_wp2_dz = ddzm(gr, jnp.sqrt(wp2))
    return -wp3 * jnp.minimum(
        5.0 / dt, 3.0 * c_wp2_splat * tau_zt * d_sqrt_wp2_dz ** 2)
