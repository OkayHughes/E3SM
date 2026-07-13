"""P3 warm-rain and freezing process rates.

Source: micro_p3.F90 (eam variant): cloud_water_autoconversion,
cloud_rain_accretion, droplet_self_collection, rain_self_collection,
subgrid_variance_scaling, ice_nucleation, CNT_couple,
cldliq_immersion_freezing, rain_immersion_freezing, calc_rime_density.
COPIED from scream_jax/p3/processes_warm.py; adaptations:

  * subgrid_variance_scaling is ACTIVE in autoconversion (expon =
    p3_qc_autocon_expon), accretion (p3_qc_accret_expon) and cloud
    immersion freezing (2.0) — scream hard-sets the coefficient to 1.
  * autoconversion: coefficient/exponents from opts, the nc exponent is
    used with its namelist SIGN (p3_nc_autocon_expon = -1.10),
    ncautr = qcaut*cons3*(1/p3_embryonic_rain_size**3), do_precip_off
    zeroing.
  * accretion: qc2qr = sgs * p3_accret_coeff * (qc*qr)**expon (the
    PRODUCT is exponentiated, as in the Fortran bfb_pow(qc*qr, e)).
  * ice_nucleation (EAM): Cooper branch needs N_nuc >= 1e-20; the
    do_Cooper_inP3 add-on applies Cooper for T > 236.18 K on top of the
    aerosol (ni_activated) term when predicting nc.
  * CNT_couple: qc gate is qsmall = 1e-14 (scream: 1e-18); mi0l_min =
    4/3*pi*rho_h2o*(4e-6)**3.0 (pow with real exponent, kept).
  * rain immersion freezing / rime density: unchanged from the copy
    (identical Fortran).
"""

import jax.numpy as jnp
from jax.scipy.special import gammaln

from . import constants as c

def _cbrt(x):
    """bfb_cbrt: x**(1/3) via pow (gfortran x**(1._rtype/3._rtype)) —
    NOT libm cbrt, which rounds differently in the last ulp."""
    return x ** (1.0 / 3.0)



def _tgamma(x):
    return jnp.exp(gammaln(x))


def subgrid_variance_scaling(relvar, expon):
    """gamma(relvar+expon)/(gamma(relvar)*relvar^expon)."""
    relvar = jnp.asarray(relvar)
    return _tgamma(relvar + expon) / (_tgamma(relvar) * relvar ** expon)


def cloud_water_autoconversion(rho, qc_incld, nc_incld, inv_qc_relvar,
                               opts, do_precip_off: bool, context):
    """KK2000 autoconversion with subgrid variance scaling
    (cloud_water_autoconversion). Returns
    (qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr)."""
    qc_incld = jnp.asarray(qc_incld)
    nc_incld = jnp.asarray(nc_incld)
    active = (qc_incld >= 1.0e-8) & context

    qc_safe = jnp.where(active, qc_incld, 1.0)
    nc_safe = jnp.where(active, nc_incld, 1.0)
    relvar_safe = jnp.where(active, jnp.asarray(inv_qc_relvar), 1.0)

    sgs = subgrid_variance_scaling(relvar_safe, opts["p3_qc_autocon_expon"])
    qcaut = jnp.where(
        active,
        sgs * opts["p3_autocon_coeff"]
        * qc_safe ** opts["p3_qc_autocon_expon"]
        * (nc_safe * 1.0e-6 * jnp.asarray(rho)) ** opts["p3_nc_autocon_expon"],
        0.0)
    ncautr = jnp.where(active,
                       qcaut * c.cons3
                       * (1.0 / opts["p3_embryonic_rain_size"] ** 3.0), 0.0)
    ncautc = jnp.where(active, qcaut * nc_safe / qc_safe, 0.0)

    if do_precip_off:
        ncautc = jnp.where(context, 0.0, ncautc)
        qcaut = jnp.where(context, 0.0, qcaut)
        ncautr = jnp.where(context, 0.0, ncautr)
    else:
        ncautc = jnp.where((qcaut == 0.0) & context, 0.0, ncautc)
        qcaut = jnp.where((ncautc == 0.0) & context, 0.0, qcaut)
    return qcaut, ncautc, ncautr


def cloud_rain_accretion(rho, inv_rho, qc_incld, nc_incld, qr_incld,
                         inv_qc_relvar, opts, context):
    """Accretion of cloud by rain (cloud_rain_accretion, iparam=3).
    Returns (qc2qr_accret_tend, nc_accret_tend)."""
    qc_incld = jnp.asarray(qc_incld)
    qr_incld = jnp.asarray(qr_incld)
    nc_incld = jnp.asarray(nc_incld)
    active = (qr_incld >= c.qsmall) & (qc_incld >= c.qsmall) & context

    qc_safe = jnp.where(active, qc_incld, 1.0)
    qr_safe = jnp.where(active, qr_incld, 1.0)
    relvar_safe = jnp.where(active, jnp.asarray(inv_qc_relvar), 1.0)

    sgs = subgrid_variance_scaling(relvar_safe, opts["p3_qc_accret_expon"])
    qcacc = jnp.where(
        active,
        sgs * opts["p3_accret_coeff"]
        * (qc_safe * qr_safe) ** opts["p3_qc_accret_expon"], 0.0)
    ncacc = jnp.where(active, qcacc * nc_incld / qc_safe, 0.0)
    ncacc = jnp.where((qcacc == 0.0) & context, 0.0, ncacc)
    qcacc = jnp.where((ncacc == 0.0) & context, 0.0, qcacc)
    return qcacc, ncacc


def droplet_self_collection(qc_incld, context):
    """Cloud droplet self-collection: zero for iparam = 3
    (droplet_self_collection)."""
    return jnp.zeros_like(jnp.asarray(qc_incld))


def rain_self_collection(rho, qr_incld, nr_incld, context):
    """Rain self-collection/breakup (rain_self_collection, iparam=3).
    Returns nr_selfcollect_tend."""
    qr_incld = jnp.asarray(qr_incld)
    nr_incld = jnp.asarray(nr_incld)
    active = (qr_incld >= c.qsmall) & context

    dum1 = 280.0e-6
    nr_safe = jnp.where(active & (nr_incld > 0), nr_incld, 1.0)
    dum2 = _cbrt(jnp.where(active, qr_incld, 1.0)
                    / (c.pi * c.rho_h2o * nr_safe))
    dum = jnp.where(dum2 < dum1, 1.0,
                    2.0 - jnp.exp(2300.0 * (dum2 - dum1)))
    return jnp.where(active,
                     dum * 5.78 * nr_incld * qr_incld * jnp.asarray(rho), 0.0)


def ice_nucleation(t_atm, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt,
                   do_predict_nc: bool, do_prescribed_ccn: bool,
                   do_cooper: bool, context):
    """Deposition/condensation-freezing ice nucleation (ice_nucleation,
    EAM form incl. the do_Cooper_inP3 branch).
    Returns (qinuc, ni_nucleat_tend)."""
    t_atm = jnp.asarray(t_atm)
    ni = jnp.asarray(ni)
    inv_rho = jnp.asarray(inv_rho)

    base = (t_atm < c.T_icenuc) & (jnp.asarray(qv_supersat_i) >= 0.05) & context

    def cooper(gate):
        dum = 0.005 * jnp.exp(0.304 * (c.T_zerodegc - t_atm)) * 1000.0 * inv_rho
        dum = jnp.minimum(dum, 100.0e3 * inv_rho)
        n_nuc = jnp.maximum(0.0, (dum - ni) * inv_dt)
        ok = gate & (n_nuc >= 1.0e-20)
        q_nuc = jnp.maximum(0.0, (dum - ni) * c.mi0 * inv_dt)
        return jnp.where(ok, q_nuc, 0.0), jnp.where(ok, n_nuc, 0.0)

    if (not do_predict_nc) or do_prescribed_ccn:
        return cooper(base)

    # aerosol-predicted branch
    n_aer = jnp.where(base,
                      jnp.maximum(0.0, (jnp.asarray(ni_activated) - ni)
                                  * inv_dt), 0.0)
    if do_cooper:
        _q_c, n_c = cooper(base & (t_atm > 236.18))
        n_tot = n_c + n_aer
        return jnp.where(base, n_tot * c.mi0, 0.0), n_tot
    return n_aer * c.mi0, n_aer


def cnt_couple(frzimm, frzcnt, frzdep, rho, qc_incld, nc_incld,
               iflag: int, context):
    """Classical-nucleation-theory coupling from the hetfrz_classnuc
    rates (CNT_couple). iflag=1 -> (ncheti_cnt, qcheti_cnt);
    iflag=2 -> (nicnt, qicnt, ninuc_cnt, qinuc_cnt).
    EAM qc gate: qsmall = 1e-14."""
    rho = jnp.asarray(rho)
    qc_incld = jnp.asarray(qc_incld)

    mi0l_min = 4.0 / 3.0 * c.pi * c.rho_h2o * (4.0e-6) ** 3.0
    nc_safe = jnp.maximum(jnp.asarray(nc_incld), 1.0e6 / rho)
    mi0l = jnp.maximum(mi0l_min, qc_incld / nc_safe)
    mask = (qc_incld > c.qsmall) & context

    if iflag == 1:
        ncheti = jnp.where(mask, jnp.asarray(frzimm) * 1.0e6 / rho, 0.0)
        return ncheti, jnp.where(mask, ncheti * mi0l, 0.0)
    elif iflag == 2:
        nicnt = jnp.where(mask, jnp.asarray(frzcnt) * 1.0e6 / rho, 0.0)
        ninuc = jnp.where(mask, jnp.asarray(frzdep) * 1.0e6 / rho, 0.0)
        return (nicnt, jnp.where(mask, nicnt * mi0l, 0.0),
                ninuc, jnp.where(mask, ninuc * c.mi0, 0.0))
    raise ValueError(f"Unhandled iflag {iflag}")


def cldliq_immersion_freezing(t_atm, lamc, mu_c, cdist1, qc_incld,
                              inv_qc_relvar, context):
    """Immersion freezing of cloud droplets (cldliq_immersion_freezing;
    subgrid variance scaling with exponent 2 ACTIVE in EAM). Returns
    (qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend)."""
    t_atm = jnp.asarray(t_atm)
    qc_incld = jnp.asarray(qc_incld)
    active = (qc_incld >= c.qsmall) & (t_atm <= c.T_rainfrz) & context

    dum1 = jnp.exp(c.aimm * (c.T_zerodegc - t_atm))
    lamc_safe = jnp.where(active & (jnp.asarray(lamc) > 0),
                          jnp.asarray(lamc), 1.0)
    dum2 = (1.0 / lamc_safe) * (1.0 / lamc_safe) * (1.0 / lamc_safe)
    relvar_safe = jnp.where(active, jnp.asarray(inv_qc_relvar), 1.0)
    sgs = subgrid_variance_scaling(relvar_safe, 2.0)

    q_nuc = sgs * c.cons6 * jnp.asarray(cdist1) \
        * _tgamma(7.0 + jnp.asarray(mu_c)) * dum1 * (dum2 * dum2)
    n_nuc = c.cons5 * jnp.asarray(cdist1) \
        * _tgamma(jnp.asarray(mu_c) + 4.0) * dum1 * dum2
    return (jnp.where(active, q_nuc, 0.0), jnp.where(active, n_nuc, 0.0))


def rain_immersion_freezing(t_atm, lamr, mu_r, cdistr, qr_incld, context):
    """Immersion freezing of rain (rain_immersion_freezing). Returns
    (qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend)."""
    t_atm = jnp.asarray(t_atm)
    qr_incld = jnp.asarray(qr_incld)
    active = (qr_incld >= c.qsmall) & (t_atm <= c.T_rainfrz) & context

    cdistr_safe = jnp.where(active & (jnp.asarray(cdistr) > 0),
                            jnp.asarray(cdistr), 1.0)
    lamr_safe = jnp.where(active & (jnp.asarray(lamr) > 0),
                          jnp.asarray(lamr), 1.0)
    exp_aimm = jnp.exp(c.aimm * (c.T_zerodegc - t_atm))
    mu_r = jnp.asarray(mu_r)

    q_nuc = c.cons6 * jnp.exp(jnp.log(cdistr_safe)
                              + jnp.log(_tgamma(7.0 + mu_r))
                              - 6.0 * jnp.log(lamr_safe)) * exp_aimm
    n_nuc = c.cons5 * jnp.exp(jnp.log(cdistr_safe)
                              + jnp.log(_tgamma(mu_r + 4.0))
                              - 3.0 * jnp.log(lamr_safe)) * exp_aimm
    return (jnp.where(active, q_nuc, 0.0), jnp.where(active, n_nuc, 0.0))


def calc_rime_density(t_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c,
                      qc_incld, qccol, context):
    """Rime density from riming rate and temperature (calc_rime_density).
    Returns (vtrmi1, rho_qm_cloud)."""
    t_atm = jnp.asarray(t_atm)
    qc_incld = jnp.asarray(qc_incld)
    qccol = jnp.asarray(qccol)
    mu_c = jnp.asarray(mu_c)

    active = (qccol >= c.qsmall) & (t_atm < c.T_zerodegc) & context
    vtrmi1 = jnp.where(active,
                       jnp.asarray(table_val_qi_fallspd) * jnp.asarray(rhofaci),
                       0.0)

    with_qc = active & (qc_incld >= c.qsmall)
    lamc_safe = jnp.where(with_qc & (jnp.asarray(lamc) > 0),
                          jnp.asarray(lamc), 1.0)
    vt_qc = (jnp.asarray(acn) * _tgamma(4.0 + c.bcn + mu_c)
             / (lamc_safe ** c.bcn * _tgamma(mu_c + 4.0)))
    d_c = (mu_c + 4.0) / lamc_safe
    v_impact = jnp.abs(vtrmi1 - vt_qc)
    inv_tc = 1.0 / jnp.minimum(-0.001, t_atm - c.T_zerodegc)
    ri = jnp.maximum(1.0, jnp.minimum(-0.5e6 * d_c * v_impact * inv_tc, 12.0))

    rho_ri = jnp.where(ri <= 8.0,
                       (0.051 + 0.114 * ri - 0.0055 * (ri * ri)) * 1000.0,
                       611.0 + 72.25 * (ri - 8.0))
    rho_qm_cloud = jnp.where(with_qc, rho_ri, 400.0)
    return vtrmi1, rho_qm_cloud
