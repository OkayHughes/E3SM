"""P3 warm-rain and freezing process rates (batch 1).

Sources (components/eamxx/src/physics/p3/impl/):
  p3_autoconversion_impl.hpp, p3_cloud_rain_acc_impl.hpp,
  p3_droplet_self_coll_impl.hpp, p3_rain_self_collection_impl.hpp,
  p3_subgrid_variance_scaling_impl.hpp, p3_ice_nucleation_impl.hpp,
  p3_ice_classical_nucleation_impl.hpp, p3_cldliq_imm_freezing_impl.hpp,
  p3_rain_imm_freezing_impl.hpp, p3_calc_rime_density_impl.hpp

Tendency outputs follow the C++ masked-set convention: values are written
only where the kernel's activity mask holds; elsewhere the incoming value
(here: an explicit `out`-style default of 0 provided by the caller/driver)
is preserved. These functions return the masked results directly with 0
defaults, matching the zero-initialized workspace arrays in p3_main.
"""

import jax.numpy as jnp
from jax.scipy.special import gammaln

from ..foundation import constants as c


def _tgamma(x):
    return jnp.exp(gammaln(x))


def subgrid_variance_scaling(relvar, expon):
    """Morrison & Gettelman (2008) eq. 9 scaling:
    gamma(relvar+expon)/(gamma(relvar)*relvar^expon)."""
    relvar = jnp.asarray(relvar)
    return _tgamma(relvar + expon) / (_tgamma(relvar) * relvar ** expon)


def cloud_water_autoconversion(rho, qc_incld, nc_incld, inv_qc_relvar,
                               opts, context):
    """KK2000-style autoconversion (Functions::cloud_water_autoconversion).
    Returns (qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr).
    NB: sgs_var_coef is hard-set to 1 in this branch of the C++."""
    qc_incld = jnp.asarray(qc_incld)
    nc_incld = jnp.asarray(nc_incld)
    active = (qc_incld >= 1e-8) & context

    cons3 = 1.0 / (c.CONS2 * opts["autoconversion_radius"] ** 3)
    qc_safe = jnp.where(active, qc_incld, 1.0)
    nc_safe = jnp.where(active, nc_incld, 1.0)

    qcaut = jnp.where(
        active,
        opts["autoconversion_prefactor"]
        * qc_safe ** opts["autoconversion_qc_exponent"]
        * (nc_safe * 1e-6 * jnp.asarray(rho)) ** (-opts["autoconversion_nc_exponent"]),
        0.0)
    ncautr = jnp.where(active, qcaut * cons3, 0.0)
    ncautc = jnp.where(active, qcaut * nc_safe / qc_safe, 0.0)

    ncautc = jnp.where((qcaut == 0.0) & context, 0.0, ncautc)
    qcaut = jnp.where((ncautc == 0.0) & context, 0.0, qcaut)
    return qcaut, ncautc, ncautr


def cloud_rain_accretion(rho, inv_rho, qc_incld, nc_incld, qr_incld,
                         inv_qc_relvar, opts, context):
    """Accretion of cloud by rain (Functions::cloud_rain_accretion).
    Returns (qc2qr_accret_tend, nc_accret_tend)."""
    qc_incld = jnp.asarray(qc_incld)
    qr_incld = jnp.asarray(qr_incld)
    nc_incld = jnp.asarray(nc_incld)
    active = (qr_incld >= c.QSMALL) & (qc_incld >= c.QSMALL) & context

    qc_safe = jnp.where(active, qc_incld, 1.0)
    qccol = jnp.where(
        active,
        opts["accretion_prefactor"]
        * qc_safe ** opts["accretion_qc_exponent"]
        * jnp.where(active, qr_incld, 1.0) ** opts["accretion_qr_exponent"],
        0.0)
    ncacc = jnp.where(active, qccol * nc_incld / qc_safe, 0.0)
    qccol = jnp.where((ncacc == 0.0) & context, 0.0, qccol)
    ncacc = jnp.where((qccol == 0.0) & context, 0.0, ncacc)
    return qccol, ncacc


def droplet_self_collection(qc_incld, context):
    """Cloud droplet self-collection: zero in this P3 version
    (Functions::droplet_self_collection sets the tendency to 0)."""
    return jnp.zeros_like(jnp.asarray(qc_incld))


def rain_self_collection(rho, qr_incld, nr_incld, opts, context):
    """Rain self-collection/breakup (Functions::rain_self_collection).
    Returns nr_selfcollect_tend."""
    qr_incld = jnp.asarray(qr_incld)
    nr_incld = jnp.asarray(nr_incld)
    active = (qr_incld >= c.QSMALL) & context

    d_break = opts["rain_selfcollection_breakup_diameter"]
    nr_safe = jnp.where(active & (nr_incld > 0), nr_incld, 1.0)
    dum2 = jnp.cbrt(jnp.where(active, qr_incld, 1.0)
                    / (c.Pi * c.RHO_H2O * nr_safe))
    dum = jnp.where(dum2 < d_break, 1.0,
                    2.0 - jnp.exp(2300.0 * (dum2 - d_break)))
    return jnp.where(active,
                     dum * opts["rain_selfcollection_prefactor"]
                     * nr_incld * qr_incld * jnp.asarray(rho), 0.0)


def ice_nucleation(temp, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt,
                   do_predict_nc: bool, do_prescribed_ccn: bool, opts, context):
    """Deposition/condensation-freezing ice nucleation
    (Functions::ice_nucleation). Returns (qv2qi_nucleat_tend, ni_nucleat_tend)."""
    temp = jnp.asarray(temp)
    ni = jnp.asarray(ni)
    t_icenuc = c.Tmelt - 15.0
    mi0 = 4.0 * c.PIOV3 * 900.0 * 1e-18

    base = (temp < t_icenuc) & (jnp.asarray(qv_supersat_i) >= 0.05) & context
    do_log = (not do_predict_nc) or do_prescribed_ccn

    if do_log:
        n_nuc = jnp.where(base, jnp.maximum(0.0, (jnp.asarray(ni_activated) - ni) * inv_dt), 0.0)
        q_nuc = n_nuc * mi0
        return q_nuc, n_nuc

    dum = 0.005 * jnp.exp(opts["deposition_nucleation_exponent"]
                          * (c.Tmelt - temp)) * 1.0e3 * jnp.asarray(inv_rho)
    dum = jnp.minimum(dum, 1.0e5 * jnp.asarray(inv_rho))
    n_nuc = jnp.maximum(0.0, (dum - ni) * inv_dt)
    ok = base & (n_nuc >= c.NSMALL)
    q_nuc = jnp.maximum(0.0, (dum - ni) * mi0 * inv_dt)
    return jnp.where(ok, q_nuc, 0.0), jnp.where(ok, n_nuc, 0.0)


def ice_classical_nucleation(frzimm, frzcnt, frzdep, rho, qc_incld, nc_incld,
                             iflag: int):
    """Heterogeneous (classical) ice nucleation from prescribed rates
    (Functions::ice_classical_nucleation).

    iflag=1: immersion freezing -> (ncheti_cnt, qcheti_cnt);
    iflag=2: contact+deposition  -> (nicnt, qicnt, ninuc_cnt, qinuc_cnt).
    """
    rho = jnp.asarray(rho)
    qc_incld = jnp.asarray(qc_incld)
    qsmall = 1.0e-18
    mi0 = 4.0 * (c.Pi / 3.0) * 900.0 * 1e-18
    mi0l_min = (4.0 / 3.0) * c.Pi * c.RHO_H2O * 4.0e-6 ** 3

    mi0l = qc_incld / jnp.maximum(jnp.asarray(nc_incld), 1.0e6 / rho)
    mi0l = jnp.maximum(mi0l_min, mi0l)
    mask = qc_incld > qsmall

    if iflag == 1:
        ncheti = jnp.where(mask, jnp.asarray(frzimm) * 1.0e6 / rho, 0.0)
        return ncheti, jnp.where(mask, ncheti * mi0l, 0.0)
    elif iflag == 2:
        nicnt = jnp.where(mask, jnp.asarray(frzcnt) * 1.0e6 / rho, 0.0)
        ninuc = jnp.where(mask, jnp.asarray(frzdep) * 1.0e6 / rho, 0.0)
        return (nicnt, jnp.where(mask, nicnt * mi0l, 0.0),
                ninuc, jnp.where(mask, ninuc * mi0, 0.0))
    raise ValueError(f"Unhandled iflag {iflag}")


def cldliq_immersion_freezing(T_atm, lamc, mu_c, cdist1, qc_incld,
                              inv_qc_relvar, opts, context):
    """Immersion freezing of cloud droplets (Bigg 1953 style)
    (Functions::cldliq_immersion_freezing). Returns
    (qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend).
    sgs_var_coef = 1 in this branch."""
    T_atm = jnp.asarray(T_atm)
    qc_incld = jnp.asarray(qc_incld)
    active = (qc_incld >= c.QSMALL) & (T_atm <= c.T_rainfrz) & context

    exp_aimm = jnp.exp(opts["immersion_freezing_exponent"] * (c.T_zerodegc - T_atm))
    lamc_safe = jnp.where(active & (jnp.asarray(lamc) > 0), jnp.asarray(lamc), 1.0)
    inv_lamc3 = (1.0 / lamc_safe) ** 3

    qchetc = jnp.where(active,
                       c.CONS6 * jnp.asarray(cdist1) * _tgamma(7.0 + jnp.asarray(mu_c))
                       * exp_aimm * inv_lamc3 ** 2, 0.0)
    nchetc = jnp.where(active,
                       c.CONS5 * jnp.asarray(cdist1) * _tgamma(4.0 + jnp.asarray(mu_c))
                       * exp_aimm * inv_lamc3, 0.0)
    return qchetc, nchetc


def rain_immersion_freezing(T_atm, lamr, mu_r, cdistr, qr_incld, opts, context):
    """Immersion freezing of rain (Functions::rain_immersion_freezing).
    Returns (qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend)."""
    T_atm = jnp.asarray(T_atm)
    qr_incld = jnp.asarray(qr_incld)
    active = (qr_incld >= c.QSMALL) & (T_atm <= c.T_rainfrz) & context

    cdistr_safe = jnp.where(active & (jnp.asarray(cdistr) > 0), jnp.asarray(cdistr), 1.0)
    lamr_safe = jnp.where(active & (jnp.asarray(lamr) > 0), jnp.asarray(lamr), 1.0)
    exp_aimm = jnp.exp(opts["immersion_freezing_exponent"] * (c.T_zerodegc - T_atm))
    mu_r = jnp.asarray(mu_r)

    qrcol = jnp.where(active,
                      c.CONS6 * jnp.exp(jnp.log(cdistr_safe)
                                        + jnp.log(_tgamma(7.0 + mu_r))
                                        - 6.0 * jnp.log(lamr_safe)) * exp_aimm, 0.0)
    nrcol = jnp.where(active,
                      c.CONS5 * jnp.exp(jnp.log(cdistr_safe)
                                        + jnp.log(_tgamma(4.0 + mu_r))
                                        - 3.0 * jnp.log(lamr_safe)) * exp_aimm, 0.0)
    return qrcol, nrcol


def calc_rime_density(T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c,
                      qc_incld, qc2qi_collect_tend, context):
    """Rime density from riming rate and temperature
    (Functions::calc_rime_density). Returns (vtrmi1, rho_qm_cloud)."""
    T_atm = jnp.asarray(T_atm)
    qc_incld = jnp.asarray(qc_incld)
    qccol = jnp.asarray(qc2qi_collect_tend)
    mu_c = jnp.asarray(mu_c)

    active = (qccol >= c.QSMALL) & (T_atm < c.T_zerodegc) & context
    vtrmi1 = jnp.where(active, jnp.asarray(table_val_qi_fallspd) * jnp.asarray(rhofaci), 0.0)

    with_qc = active & (qc_incld >= c.QSMALL)
    lamc_safe = jnp.where(with_qc & (jnp.asarray(lamc) > 0), jnp.asarray(lamc), 1.0)
    vt_qc = (jnp.asarray(acn) * _tgamma(4.0 + c.bcn + mu_c)
             / (lamc_safe ** c.bcn * _tgamma(4.0 + mu_c)))
    d_c = (4.0 + mu_c) / lamc_safe
    v_impact = jnp.abs(vtrmi1 - vt_qc)
    inv_tc = 1.0 / jnp.minimum(-0.001, T_atm - c.T_zerodegc)
    ri = jnp.maximum(1.0, jnp.minimum(-0.5e6 * d_c * v_impact * inv_tc, 12.0))

    rho_ri = jnp.where(ri <= 8.0,
                       (0.051 + 0.114 * ri - 0.0055 * ri ** 2) * 1000.0,
                       611.0 + 72.25 * (ri - 8.0))
    rho_qm_cloud = jnp.where(with_qc, rho_ri,
                             jnp.where(active, 400.0, 400.0))
    rho_qm_cloud = jnp.where(active, rho_qm_cloud, 400.0)
    return vtrmi1, rho_qm_cloud
