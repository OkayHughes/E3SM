"""p3_main_part3: final consistency clipping and diagnostic outputs.

Sources: components/eamxx/src/physics/p3/impl/p3_main_impl_part3.hpp and
calc_bulk_rho_rime (p3_calc_bulk_rho_rime in p3_functions.hpp impl).
"""

import jax.numpy as jnp

from ..foundation import constants as c
from ..foundation import smoothing
from .dsd import get_cloud_dsd2, get_rain_dsd2
from .conservation import impose_max_total_ni
from .table_lookups import apply_table_ice, lookup_ice


def calc_bulk_rho_rime(qi_tot, qi_rim, bi_rim, opts, context,
                       smooth_width=0.0):
    """Bulk rime density with limiters (Functions::calc_bulk_rho_rime).
    Returns (rho_rime, qi_rim, bi_rim)."""
    qi_tot = jnp.asarray(qi_tot)
    qi_rim = jnp.asarray(qi_rim)
    bi_rim = jnp.asarray(bi_rim)

    if smooth_width == 0.0:
        gt = (bi_rim >= c.BSMALL) & context
        lt = (bi_rim < c.BSMALL) & context

        bi_safe = jnp.where(gt, bi_rim, 1.0)
        rho_rime = jnp.where(gt, qi_rim / bi_safe, 0.0)

        lo = rho_rime < opts["min_rime_rho"]
        hi = rho_rime > opts["max_rime_rho"]
        rho_rime = jnp.where(gt & lo, opts["min_rime_rho"], rho_rime)
        rho_rime = jnp.where(gt & hi, opts["max_rime_rho"], rho_rime)
        adjust = gt & (lo | hi)
        rho_safe = jnp.where(rho_rime > 0, rho_rime, 1.0)
        bi_rim = jnp.where(adjust, qi_rim / rho_safe, bi_rim)

        qi_rim = jnp.where(lt, 0.0, qi_rim)
        bi_rim = jnp.where(lt, 0.0, bi_rim)
        rho_rime = jnp.where(lt, 0.0, rho_rime)

        over = (qi_rim > qi_tot) & (rho_rime > 0) & context
        qi_rim = jnp.where(over, qi_tot, qi_rim)
        bi_rim = jnp.where(over, qi_rim / rho_safe, bi_rim)

        tiny = (qi_rim < c.QSMALL) & context
        qi_rim = jnp.where(tiny, 0.0, qi_rim)
        bi_rim = jnp.where(tiny, 0.0, bi_rim)
        return rho_rime, qi_rim, bi_rim

    # PORT_NOTES (smoothing, family "rime"): JUMP — s = bi_rim - BSMALL,
    # scale = BSMALL (1e-15). Where bi_rim falls below BSMALL the C++
    # snaps qi_rim (finite, NOT -> 0 at the threshold) and bi_rim to
    # zero; the on/off branches are blended with h = step(s). KINKS
    # left hard: the min/max rho clamps (bi_rim = qi_rim/rho at the
    # clamp boundary equals the unadjusted bi_rim -> continuous), the
    # `over` cap (qi_rim = qi_tot continuous at the crossing) and the
    # final qi_rim < QSMALL snap (jump size QSMALL, negligible).
    h = jnp.where(context,
                  smoothing.step(bi_rim - c.BSMALL, smooth_width,
                                 scale=c.BSMALL),
                  0.0)
    # NaN-safety: floor the denominator so 1/bi^2 in the VJP cannot
    # overflow to inf (0*inf = NaN through the blend). For
    # bi_rim < 1e-100 the blend weight h is exactly 0 (sigmoid argument
    # << -745), so outputs are unchanged.
    bi_safe = jnp.where(bi_rim > 0, jnp.maximum(bi_rim, 1.0e-100), 1.0)
    rho_on = qi_rim / bi_safe

    lo = rho_on < opts["min_rime_rho"]
    hi = rho_on > opts["max_rime_rho"]
    rho_on = jnp.where(lo, opts["min_rime_rho"], rho_on)
    rho_on = jnp.where(hi, opts["max_rime_rho"], rho_on)
    rho_safe = jnp.where(rho_on > 0, rho_on, 1.0)
    bi_on = jnp.where(lo | hi, qi_rim / rho_safe, bi_rim)

    # blend the >=BSMALL branch (qi_rim, bi_on, rho_on) with the snap
    # branch (0, 0, 0); leave non-context entries untouched
    qi_b = h * qi_rim
    bi_b = h * bi_on
    rho_b = h * rho_on

    rho_bsafe = jnp.where(rho_b > 0, rho_b, 1.0)
    over = (qi_b > qi_tot) & (rho_b > 0) & context
    qi_b = jnp.where(over, qi_tot, qi_b)
    bi_b = jnp.where(over, qi_b / rho_bsafe, bi_b)

    tiny = (qi_b < c.QSMALL) & context
    qi_b = jnp.where(tiny, 0.0, qi_b)
    bi_b = jnp.where(tiny, 0.0, bi_b)

    qi_rim = jnp.where(context, qi_b, qi_rim)
    bi_rim = jnp.where(context, bi_b, bi_rim)
    rho_rime = jnp.where(context, rho_b, 0.0)
    return rho_rime, qi_rim, bi_rim


def p3_main_part3(max_total_ni, dnu, ice_table_vals,
                  inv_exner, cld_frac_l, cld_frac_r, cld_frac_i,
                  rho, inv_rho, rhofaci,
                  qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm,
                  vap_liq_exchange, ze_rain, ze_ice, opts,
                  diag_eff_radius_qc_in=None, diag_eff_radius_qr_in=None,
                  diag_eff_radius_qi_in=None):
    """Returns a dict with the updated state and diagnostics
    (diag_eff_radius_qc/qr/qi, diag_vm_qi, diag_diam_qi, rho_qi,
    diag_equiv_reflectivity, mu_c, lamc, mu_r, lamr, updated ze arrays).

    The eff-radius arrays are masked-set over their prior contents (the
    p3_main_init values 10/25/500 um) — pass those via the *_in args;
    they default to zero. diag_vm_qi/diag_diam_qi/rho_qi are zero-init
    in the C++ (zero_init list), so zeros are correct for them."""
    qv, th_atm = jnp.asarray(qv), jnp.asarray(th_atm)
    qc, nc, qr, nr = (jnp.asarray(a) for a in (qc, nc, qr, nr))
    qi, ni, qm, bm = (jnp.asarray(a) for a in (qi, ni, qm, bm))
    cld_frac_l = jnp.asarray(cld_frac_l)
    cld_frac_r = jnp.asarray(cld_frac_r)
    cld_frac_i = jnp.asarray(cld_frac_i)
    inv_exner = jnp.asarray(inv_exner)
    vap_liq_exchange = jnp.asarray(vap_liq_exchange)
    ze_rain = jnp.asarray(ze_rain)
    ze_ice = jnp.asarray(ze_ice)

    diag_eff_radius_qc = (jnp.zeros_like(qc) if diag_eff_radius_qc_in is None
                          else jnp.asarray(diag_eff_radius_qc_in))
    diag_eff_radius_qr = (jnp.zeros_like(qc) if diag_eff_radius_qr_in is None
                          else jnp.asarray(diag_eff_radius_qr_in))
    diag_eff_radius_qi = (jnp.zeros_like(qc) if diag_eff_radius_qi_in is None
                          else jnp.asarray(diag_eff_radius_qi_in))
    diag_vm_qi = jnp.zeros_like(qc)
    diag_diam_qi = jnp.zeros_like(qc)
    rho_qi = jnp.zeros_like(qc)

    # ---- Cloud ----
    qc_gt = qc >= c.QSMALL
    qc_incld = qc / cld_frac_l
    nc_incld = nc / cld_frac_l
    nc_incld2, mu_c, nu, lamc, _, _ = get_cloud_dsd2(qc_incld, nc_incld, rho, qc_gt)
    nc = jnp.where(qc_gt, nc_incld2 * cld_frac_l, nc)
    lamc_safe = jnp.where(qc_gt & (lamc > 0), lamc, 1.0)
    diag_eff_radius_qc = jnp.where(qc_gt, 0.5 * (mu_c + 3.0) / lamc_safe,
                                   diag_eff_radius_qc)
    qc_small = ~qc_gt
    qv = jnp.where(qc_small, qv + qc, qv)
    th_atm = jnp.where(qc_small, th_atm - inv_exner * qc * c.LatVap * c.INV_CP,
                       th_atm)
    vap_liq_exchange = jnp.where(qc_small, vap_liq_exchange - qc,
                                 vap_liq_exchange)
    qc = jnp.where(qc_small, 0.0, qc)
    nc = jnp.where(qc_small, 0.0, nc)

    # ---- Rain ----
    qr_gt = qr >= c.QSMALL
    qr_incld = qr / cld_frac_r
    nr_incld = nr / cld_frac_r
    nr_incld2, mu_r, lamr = get_rain_dsd2(qr_incld, nr_incld,
                                          opts["constant_mu_rain"], qr_gt)
    nr = jnp.where(qr_gt, nr_incld2 * cld_frac_r, nr)
    lamr_safe = jnp.where(qr_gt & (lamr > 0), lamr, 1.0)
    ze_rain = jnp.where(qr_gt,
                        nr * (mu_r + 6) * (mu_r + 5) * (mu_r + 4)
                        * (mu_r + 3) * (mu_r + 2) * (mu_r + 1) / lamr_safe ** 6.0,
                        ze_rain)
    ze_rain = jnp.where(qr_gt, jnp.maximum(ze_rain, 1.0e-22), ze_rain)
    diag_eff_radius_qr = jnp.where(qr_gt, 1.5 / lamr_safe, diag_eff_radius_qr)
    qr_small = ~qr_gt
    qv = jnp.where(qr_small, qv + qr, qv)
    th_atm = jnp.where(qr_small, th_atm - inv_exner * qr * c.LatVap * c.INV_CP,
                       th_atm)
    vap_liq_exchange = jnp.where(qr_small, vap_liq_exchange - qr,
                                 vap_liq_exchange)
    qr = jnp.where(qr_small, 0.0, qr)
    nr = jnp.where(qr_small, 0.0, nr)

    # ---- Ice ----
    qi_gt = qi >= c.QSMALL
    ni = jnp.maximum(ni, c.NSMALL)
    qi_incld = qi / cld_frac_i
    ni_incld = ni / cld_frac_i
    qm_incld = qm / cld_frac_i
    bm_incld = bm / cld_frac_i

    rhop, qm_incld, bm_incld = calc_bulk_rho_rime(qi_incld, qm_incld, bm_incld,
                                                  opts, qi_gt)
    qm = jnp.where(qi_gt, qm_incld * cld_frac_i, qm)
    bm = jnp.where(qi_gt, bm_incld * cld_frac_i, bm)

    ni_incld = impose_max_total_ni(ni_incld, max_total_ni, inv_rho,
                                   jnp.ones_like(qi, dtype=bool))

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
    ni = ni_incld * cld_frac_i

    qm_small = (qm < c.QSMALL) & qi_gt
    qm = jnp.where(qm_small, 0.0, qm)
    bm = jnp.where(qm_small, 0.0, bm)

    diag_vm_qi = jnp.where(qi_gt, t_fallspd * jnp.asarray(rhofaci), diag_vm_qi)
    diag_eff_radius_qi = jnp.where(qi_gt, t_eff_rad, diag_eff_radius_qi)
    diag_diam_qi = jnp.where(qi_gt, t_diam, diag_diam_qi)
    rho_qi = jnp.where(qi_gt, t_dens, rho_qi)

    ze_ice = jnp.where(qi_gt,
                       ze_ice + 0.1892 * t_refl * ni_incld * jnp.asarray(rho),
                       ze_ice)
    ze_ice = jnp.where(qi_gt, jnp.maximum(ze_ice, 1.0e-22), ze_ice)
    ze_ice = jnp.where(qi_gt, ze_ice * cld_frac_i, ze_ice)

    qi_small = ~qi_gt
    qv = jnp.where(qi_small, qv + qi, qv)
    th_atm = jnp.where(qi_small,
                       th_atm - inv_exner * qi * (c.LatVap + c.LatIce) * c.INV_CP,
                       th_atm)
    qi = jnp.where(qi_small, 0.0, qi)
    ni = jnp.where(qi_small, 0.0, ni)
    qm = jnp.where(qi_small, 0.0, qm)
    bm = jnp.where(qi_small, 0.0, bm)
    diag_diam_qi = jnp.where(qi_small, 0.0, diag_diam_qi)

    diag_equiv_reflectivity = 10.0 * jnp.log10((ze_rain + ze_ice) * 1.0e18)
    nr = jnp.where(qr < c.QSMALL, 0.0, nr)

    return {
        "qv": qv, "th_atm": th_atm, "qc": qc, "nc": nc, "qr": qr, "nr": nr,
        "qi": qi, "ni": ni, "qm": qm, "bm": bm,
        "mu_c": mu_c, "lamc": lamc, "mu_r": mu_r, "lamr": lamr,
        "vap_liq_exchange": vap_liq_exchange,
        "ze_rain": ze_rain, "ze_ice": ze_ice,
        "diag_vm_qi": diag_vm_qi, "diag_eff_radius_qi": diag_eff_radius_qi,
        "diag_diam_qi": diag_diam_qi, "rho_qi": rho_qi,
        "diag_equiv_reflectivity": diag_equiv_reflectivity,
        "diag_eff_radius_qc": diag_eff_radius_qc,
        "diag_eff_radius_qr": diag_eff_radius_qr,
    }
