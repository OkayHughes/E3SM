"""p3_main_part3: final consistency clipping and diagnostic outputs.

Source: micro_p3.F90 p3_main_part3 (eam variant). COPIED from
scream_jax/p3/main_part3.py; adaptations:

  * p3_mincdnc floor on nc_incld after the cloud dsd call (affects the
    nc output at cloudy levels).
  * diag_ze_rain / diag_ze_ice outputs (10*log10(ze*1e18) in the
    condensate branches; callers keep -99 elsewhere).
  * diag_eff_radius_qc is zeroed in the no-cloud branch and
    diag_eff_radius_qi / diag_vm_qi (plus diag_diam_qi) are zeroed in
    the no-ice branch (scream leaves the init values).
  * no diag_eff_radius_qr in EAM.
  * mu_c/lamc/mu_r/lamr merge with their incoming values where the
    dsd routines are not called (Fortran intent(out) scalars stay
    untouched in the no-condensate branches).
  * ni floor at nsmall and impose_max_total_ni only inside the
    qi >= qsmall branch (identical net effect, kept masked).
"""

import jax.numpy as jnp

from . import constants as c
from .conservation import impose_max_total_ni
from .dsd import calc_bulk_rho_rime, get_cloud_dsd2, get_rain_dsd2
from .table_lookups import apply_table_ice, lookup_ice


def p3_main_part3(p3_max_mean_rain_size, mincdnc,
                  exner, cld_frac_l, cld_frac_r, cld_frac_i,
                  rho, inv_rho, rhofaci,
                  qv, th, qc, nc, qr, nr, qi, ni, qm, bm,
                  mu_c_in, lamc_in, mu_r_in, lamr_in,
                  vap_liq_exchange, ze_rain, ze_ice, ice_table_vals):
    """Returns a dict with the updated state and diagnostics."""
    qv, th = jnp.asarray(qv), jnp.asarray(th)
    qc, nc, qr, nr = (jnp.asarray(a) for a in (qc, nc, qr, nr))
    qi, ni, qm, bm = (jnp.asarray(a) for a in (qi, ni, qm, bm))
    cld_frac_l = jnp.asarray(cld_frac_l)
    cld_frac_r = jnp.asarray(cld_frac_r)
    cld_frac_i = jnp.asarray(cld_frac_i)
    exner = jnp.asarray(exner)
    vap_liq_exchange = jnp.asarray(vap_liq_exchange)
    ze_rain = jnp.asarray(ze_rain)
    ze_ice = jnp.asarray(ze_ice)

    diag_vm_qi = jnp.zeros_like(qc)
    diag_diam_qi = jnp.zeros_like(qc)
    rho_qi = jnp.zeros_like(qc)
    diag_ze_rain = jnp.full_like(qc, -99.0)
    diag_ze_ice = jnp.full_like(qc, -99.0)

    # ---- cloud ----
    qc_gt = qc >= c.qsmall
    qc_incld = qc / cld_frac_l
    nc_incld = nc / cld_frac_l
    nc_incld2, mu_c2, nu, lamc2, _, _ = get_cloud_dsd2(qc_incld, nc_incld,
                                                       rho, qc_gt)
    mu_c = jnp.where(qc_gt, mu_c2, jnp.asarray(mu_c_in))
    lamc = jnp.where(qc_gt, lamc2, jnp.asarray(lamc_in))
    nc_incld = jnp.where(qc_gt, nc_incld2, nc_incld)
    # mincdnc floor (opts value may be traced under jit: use a mask)
    mincdnc = jnp.asarray(mincdnc)
    nc_incld = jnp.where(qc_gt & (mincdnc > 0.0),
                         jnp.maximum(nc_incld, mincdnc / rho), nc_incld)
    lamc_safe = jnp.where(qc_gt & (lamc > 0), lamc, 1.0)
    diag_eff_radius_qc = jnp.where(qc_gt, 0.5 * (mu_c + 3.0) / lamc_safe,
                                   0.0)
    nc = jnp.where(qc_gt, nc_incld * cld_frac_l, nc)
    qc_small = ~qc_gt
    qv = jnp.where(qc_small, qv + qc, qv)
    th = jnp.where(qc_small, th - exner * qc * c.latvap * c.inv_cp, th)
    vap_liq_exchange = jnp.where(qc_small, vap_liq_exchange - qc,
                                 vap_liq_exchange)
    qc = jnp.where(qc_small, 0.0, qc)
    nc = jnp.where(qc_small, 0.0, nc)

    # ---- rain ----
    qr_gt = qr >= c.qsmall
    qr_incld = qr / cld_frac_r
    nr_incld = nr / cld_frac_r
    nr_incld2, mu_r2, lamr2, _, _ = get_rain_dsd2(qr_incld, nr_incld,
                                                  p3_max_mean_rain_size,
                                                  qr_gt)
    mu_r = jnp.where(qr_gt, mu_r2, jnp.asarray(mu_r_in))
    lamr = jnp.where(qr_gt, lamr2, jnp.asarray(lamr_in))
    nr = jnp.where(qr_gt, nr_incld2 * cld_frac_r, nr)
    lamr_safe = jnp.where(qr_gt & (lamr > 0), lamr, 1.0)
    ze_rain = jnp.where(qr_gt,
                        nr * (mu_r + 6.0) * (mu_r + 5.0) * (mu_r + 4.0)
                        * (mu_r + 3.0) * (mu_r + 2.0) * (mu_r + 1.0)
                        / lamr_safe ** 6.0,
                        ze_rain)
    ze_rain = jnp.where(qr_gt, jnp.maximum(ze_rain, 1.0e-22), ze_rain)
    diag_ze_rain = jnp.where(qr_gt, 10.0 * jnp.log10(ze_rain * 1.0e18),
                             diag_ze_rain)
    qr_small = ~qr_gt
    qv = jnp.where(qr_small, qv + qr, qv)
    th = jnp.where(qr_small, th - exner * qr * c.latvap * c.inv_cp, th)
    vap_liq_exchange = jnp.where(qr_small, vap_liq_exchange - qr,
                                 vap_liq_exchange)
    qr = jnp.where(qr_small, 0.0, qr)
    nr = jnp.where(qr_small, 0.0, nr)

    # ---- ice ----
    qi_gt = qi >= c.qsmall
    ni = jnp.where(qi_gt, jnp.maximum(ni, c.nsmall), ni)
    qi_incld = qi / cld_frac_i
    ni_incld = ni / cld_frac_i
    qm_incld = qm / cld_frac_i
    bm_incld = bm / cld_frac_i

    rhop, qm_incld, bm_incld = calc_bulk_rho_rime(qi_incld, qm_incld,
                                                  bm_incld, qi_gt)
    qm = jnp.where(qi_gt, qm_incld * cld_frac_i, qm)
    bm = jnp.where(qi_gt, bm_incld * cld_frac_i, bm)

    ni_incld = impose_max_total_ni(ni_incld, c.max_total_ni, inv_rho, qi_gt)

    ti = lookup_ice(qi_incld, ni_incld, qm_incld, rhop, qi_gt)
    t_fallspd = jnp.where(qi_gt, apply_table_ice(1, ice_table_vals, ti), 0.0)
    t_eff_rad = jnp.where(qi_gt, apply_table_ice(5, ice_table_vals, ti), 0.0)
    t_lammax = jnp.where(qi_gt, apply_table_ice(6, ice_table_vals, ti), 0.0)
    t_lammin = jnp.where(qi_gt, apply_table_ice(7, ice_table_vals, ti), 0.0)
    t_refl = jnp.where(qi_gt, apply_table_ice(8, ice_table_vals, ti), 0.0)
    t_diam = jnp.where(qi_gt, apply_table_ice(10, ice_table_vals, ti), 0.0)
    t_dens = jnp.where(qi_gt, apply_table_ice(11, ice_table_vals, ti), 0.0)

    ni_incld = jnp.where(qi_gt, jnp.minimum(ni_incld, t_lammax * ni_incld),
                         ni_incld)
    ni_incld = jnp.where(qi_gt, jnp.maximum(ni_incld, t_lammin * ni_incld),
                         ni_incld)
    ni = jnp.where(qi_gt, ni_incld * cld_frac_i, ni)

    qm_small = (qm < c.qsmall) & qi_gt
    qm = jnp.where(qm_small, 0.0, qm)
    bm = jnp.where(qm_small, 0.0, bm)

    diag_vm_qi = jnp.where(qi_gt, t_fallspd * jnp.asarray(rhofaci), 0.0)
    diag_eff_radius_qi = jnp.where(qi_gt, t_eff_rad, 0.0)
    diag_diam_qi = jnp.where(qi_gt, t_diam, 0.0)
    rho_qi = jnp.where(qi_gt, t_dens, rho_qi)

    ze_ice = jnp.where(qi_gt,
                       ze_ice + 0.1892 * t_refl * ni_incld
                       * jnp.asarray(rho), ze_ice)
    ze_ice = jnp.where(qi_gt, jnp.maximum(ze_ice, 1.0e-22), ze_ice)
    ze_ice = jnp.where(qi_gt, ze_ice * cld_frac_i, ze_ice)
    diag_ze_ice = jnp.where(qi_gt, 10.0 * jnp.log10(ze_ice * 1.0e18),
                            diag_ze_ice)

    qi_small = ~qi_gt
    qv = jnp.where(qi_small, qv + qi, qv)
    th = jnp.where(qi_small, th - exner * qi * c.latsub * c.inv_cp, th)
    qi = jnp.where(qi_small, 0.0, qi)
    ni = jnp.where(qi_small, 0.0, ni)
    qm = jnp.where(qi_small, 0.0, qm)
    bm = jnp.where(qi_small, 0.0, bm)

    diag_equiv_reflectivity = 10.0 * jnp.log10((ze_rain + ze_ice) * 1.0e18)
    nr = jnp.where(qr < c.qsmall, 0.0, nr)

    return {
        "qv": qv, "th": th, "qc": qc, "nc": nc, "qr": qr, "nr": nr,
        "qi": qi, "ni": ni, "qm": qm, "bm": bm,
        "mu_c": mu_c, "lamc": lamc, "mu_r": mu_r, "lamr": lamr,
        "vap_liq_exchange": vap_liq_exchange,
        "ze_rain": ze_rain, "ze_ice": ze_ice,
        "diag_vm_qi": diag_vm_qi, "diag_eff_radius_qi": diag_eff_radius_qi,
        "diag_diam_qi": diag_diam_qi, "rho_qi": rho_qi,
        "diag_equiv_reflectivity": diag_equiv_reflectivity,
        "diag_eff_radius_qc": diag_eff_radius_qc,
        "diag_ze_rain": diag_ze_rain, "diag_ze_ice": diag_ze_ice,
    }
