"""P3 ice-phase process rates and rain evaporation.

Source: micro_p3.F90 (eam variant): ice_cldliq_collection,
ice_rain_collection, ice_self_collection, ice_melting,
ice_cldliq_wet_growth, ice_deposition_sublimation,
calc_ice_relaxation_timescale, calc_liq_relaxation_timescale,
evaporate_rain (+ its 3 helpers). COPIED from
scream_jax/p3/processes_ice.py; adaptations:

  * qsat0 uses the EAM qv_sat (MurphyKoop with p - es denominator).
  * collection efficiencies are the module constants eci/eri (no opts).
  * wet growth does NOT clamp qccol/qrcol at zero after shedding
    (scream added max(0, ..) — "collection shouldn't work backwards").
  * ice_deposition_sublimation is the EAM (older) form: qidep =
    epsi/abi*(qv - qv_sat_i) with NO min(epsi/abi, inv_dt) limiter;
    deposition only for T < 0C (else it becomes sublimation);
    sublimation at any T; Bergeron rate scaled by p3_wbf_coeff;
    ni_sublim from the (possibly zero) sublimation at any T.
  * evaporate_rain: identical algorithm (the C++ inherited it); kept,
    with constants from this package.
"""

import jax.numpy as jnp
from jax.scipy.special import gammaln

from . import constants as c
from .saturation import qv_sat
from .table_lookups import apply_table3, lookup_table3

def _cbrt(x):
    """bfb_cbrt: x**(1/3) via pow (gfortran x**(1._rtype/3._rtype)) —
    NOT libm cbrt, which rounds differently in the last ulp."""
    return x ** (1.0 / 3.0)



def _tgamma(x):
    return jnp.exp(gammaln(x))


def ice_cldliq_collection(rho, t_atm, rhofaci, table_val_qc2qi_collect,
                          qi_incld, qc_incld, ni_incld, nc_incld, context):
    """Collection of cloud droplets by ice; sheds to rain above freezing
    (ice_cldliq_collection). Returns
    (qccol, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc)."""
    t_atm = jnp.asarray(t_atm)
    qi_incld = jnp.asarray(qi_incld)
    qc_incld = jnp.asarray(qc_incld)

    freezing = t_atm <= c.T_zerodegc
    both = (qi_incld >= c.qsmall) & (qc_incld >= c.qsmall) & context

    rate_q = (jnp.asarray(rhofaci) * jnp.asarray(table_val_qc2qi_collect)
              * qc_incld * c.eci * jnp.asarray(rho) * jnp.asarray(ni_incld))
    rate_n = (jnp.asarray(rhofaci) * jnp.asarray(table_val_qc2qi_collect)
              * jnp.asarray(nc_incld) * c.eci * jnp.asarray(rho)
              * jnp.asarray(ni_incld))

    qccol = jnp.where(both & freezing, rate_q, 0.0)
    qcshd = jnp.where(both & ~freezing, rate_q, 0.0)
    nccol = jnp.where(both, rate_n, 0.0)
    ncshdc = jnp.where(both & ~freezing, qcshd * c.inv_dropmass, 0.0)
    return qccol, nccol, qcshd, ncshdc


def ice_rain_collection(rho, t_atm, rhofaci, logn0r, table_val_nr_collect,
                        table_val_qr2qi_collect, qi_incld, ni_incld,
                        qr_incld, context):
    """Collection of rain by ice (ice_rain_collection).
    Returns (qrcol, nr_collect_tend)."""
    t_atm = jnp.asarray(t_atm)
    both = ((jnp.asarray(qi_incld) >= c.qsmall)
            & (jnp.asarray(qr_incld) >= c.qsmall) & context)
    freezing = t_atm <= c.T_zerodegc

    qrcol = jnp.where(both & freezing,
                      10.0 ** (jnp.asarray(table_val_qr2qi_collect)
                               + jnp.asarray(logn0r))
                      * jnp.asarray(rho) * jnp.asarray(rhofaci) * c.eri
                      * jnp.asarray(ni_incld), 0.0)
    nrcol = jnp.where(both,
                      10.0 ** (jnp.asarray(table_val_nr_collect)
                               + jnp.asarray(logn0r))
                      * jnp.asarray(rho) * jnp.asarray(rhofaci) * c.eri
                      * jnp.asarray(ni_incld), 0.0)
    return qrcol, nrcol


def ice_self_collection(rho, rhofaci, table_val_ni_self_collect, eii,
                        qm_incld, qi_incld, ni_incld, context):
    """Ice-ice aggregation with rime-fraction sticking efficiency
    (ice_self_collection). Returns ni_selfcollect_tend."""
    qi_incld = jnp.asarray(qi_incld)
    qm_incld = jnp.asarray(qm_incld)
    qi_ok = (qi_incld >= c.qsmall) & context
    qm_pos = qm_incld > 0.0

    qi_safe = jnp.where(qi_ok, qi_incld, 1.0)
    tmp1 = jnp.where(qi_ok & qm_pos, qm_incld / qi_safe, 0.0)
    eii_fact = jnp.where(
        qm_pos,
        jnp.where(tmp1 < 0.6, 1.0,
                  jnp.where(tmp1 < 0.9, 1.0 - (tmp1 - 0.6) / 0.3, 0.0)),
        1.0)
    return jnp.where(qi_ok,
                     jnp.asarray(table_val_ni_self_collect) * jnp.asarray(rho)
                     * jnp.asarray(eii) * eii_fact * jnp.asarray(rhofaci)
                     * jnp.asarray(ni_incld) * jnp.asarray(ni_incld), 0.0)


def ice_melting(rho, t_atm, pres, rhofaci, table_val_qi2qr_melting,
                table_val_qi2qr_vent_melt, latent_heat_vapor,
                latent_heat_fusion, dv, sc, mu, kap, qv,
                qi_incld, ni_incld, context):
    """Melting of ice to rain (ice_melting).
    Returns (qi2qr_melt_tend, ni2nr_melt_tend)."""
    t_atm = jnp.asarray(t_atm)
    qi_incld = jnp.asarray(qi_incld)
    active = (qi_incld >= c.qsmall) & (t_atm > c.T_zerodegc) & context

    qsat0 = qv_sat(jnp.full_like(t_atm, c.T_zerodegc), jnp.asarray(pres),
                   False)
    melt = ((jnp.asarray(table_val_qi2qr_melting)
             + jnp.asarray(table_val_qi2qr_vent_melt)
             * _cbrt(jnp.asarray(sc))
             * jnp.sqrt(jnp.asarray(rhofaci) * jnp.asarray(rho)
                        / jnp.asarray(mu)))
            * ((t_atm - c.T_zerodegc) * jnp.asarray(kap)
               - jnp.asarray(rho) * latent_heat_vapor * jnp.asarray(dv)
               * (qsat0 - jnp.asarray(qv)))
            * 2.0 * c.pi / latent_heat_fusion) * jnp.asarray(ni_incld)
    qimlt = jnp.maximum(jnp.where(active, melt, 0.0), 0.0)
    qi_safe = jnp.where(active, qi_incld, 1.0)
    nimlt = jnp.where(active, qimlt * (jnp.asarray(ni_incld) / qi_safe), 0.0)
    return qimlt, nimlt


def ice_cldliq_wet_growth(rho, t_atm, pres, rhofaci, table_val_qi2qr_melting,
                          table_val_qi2qr_vent_melt, latent_heat_vapor,
                          latent_heat_fusion, dv, kap, mu, sc, qv,
                          qc_incld, qi_incld, ni_incld, qr_incld,
                          qrcol, qccol, nr_ice_shed_tend,
                          qc2qr_ice_shed_tend, context):
    """Wet growth of rimed ice; excess collected water sheds
    (ice_cldliq_wet_growth; NO max(0,..) clamp on qccol/qrcol in EAM).
    Returns (log_wetgrowth, qrcol, qccol, qwgrth, nr_ice_shed_tend,
    qc2qr_ice_shed_tend)."""
    t_atm = jnp.asarray(t_atm)
    qi_incld = jnp.asarray(qi_incld)
    qccol = jnp.asarray(qccol)
    qrcol = jnp.asarray(qrcol)
    nrshdr = jnp.asarray(nr_ice_shed_tend)
    qcshd = jnp.asarray(qc2qr_ice_shed_tend)

    any_if = ((qi_incld >= c.qsmall)
              & ((jnp.asarray(qc_incld) + jnp.asarray(qr_incld)) >= 1.0e-6)
              & (t_atm < c.T_zerodegc) & context)

    qsat0 = qv_sat(jnp.full_like(t_atm, c.T_zerodegc), jnp.asarray(pres),
                   False)
    growth = ((jnp.asarray(table_val_qi2qr_melting)
               + jnp.asarray(table_val_qi2qr_vent_melt)
               * _cbrt(jnp.asarray(sc))
               * jnp.sqrt(jnp.asarray(rhofaci) * jnp.asarray(rho)
                          / jnp.asarray(mu)))
              * 2.0 * c.pi
              * (jnp.asarray(rho) * latent_heat_vapor * jnp.asarray(dv)
                 * (qsat0 - jnp.asarray(qv))
                 - (t_atm - c.T_zerodegc) * jnp.asarray(kap))
              / (latent_heat_fusion + c.cpw * (t_atm - c.T_zerodegc))) \
        * jnp.asarray(ni_incld)
    qwgrth = jnp.where(any_if, jnp.maximum(growth, 0.0), 0.0)

    dum = jnp.maximum(0.0, (qccol + qrcol) - qwgrth)
    dum_ok = any_if & (dum >= 1.0e-10)
    col_ok = dum_ok & ((qccol + qrcol) >= 1.0e-10)
    total_safe = jnp.where((qccol + qrcol) > 0, qccol + qrcol, 1.0)
    dum1 = 1.0 / total_safe

    nrshdr = jnp.where(dum_ok, nrshdr + dum * 1.923e6, nrshdr)
    qcshd = jnp.where(col_ok, qcshd + dum * qccol * dum1, qcshd)
    qccol_new = jnp.where(col_ok, qccol - dum * qccol * dum1, qccol)
    qrcol_new = jnp.where(col_ok, qrcol - dum * qrcol * dum1, qrcol)
    log_wetgrowth = dum_ok
    return log_wetgrowth, qrcol_new, qccol_new, qwgrth, nrshdr, qcshd


def ice_deposition_sublimation(qi_incld, ni_incld, t_atm, qv_sat_l, qv_sat_i,
                               epsi, abi, qv, p3_wbf_coeff, context):
    """Vapor deposition / sublimation of ice and Bergeron rate
    (ice_deposition_sublimation, EAM form). Returns
    (qidep, qi2qv_sublim_tend, ni_sublim_tend, qiberg)."""
    qi_incld = jnp.asarray(qi_incld)
    t_atm = jnp.asarray(t_atm)
    qv_sat_i = jnp.asarray(qv_sat_i)
    active = (qi_incld >= c.qsmall) & context

    abi_safe = jnp.where(jnp.asarray(abi) != 0, jnp.asarray(abi), 1.0)
    oabi = 1.0 / abi_safe
    dep0 = jnp.asarray(epsi) * oabi * (jnp.asarray(qv) - qv_sat_i)
    freezing = t_atm < c.T_zerodegc

    keep_dep = freezing & (dep0 > 0.0)
    qidep = jnp.where(active & keep_dep, dep0, 0.0)
    qisub = jnp.where(active & ~keep_dep, -jnp.minimum(dep0, 0.0), 0.0)

    qiberg = jnp.where(
        active & freezing,
        jnp.maximum(jnp.asarray(epsi) * oabi
                    * (jnp.asarray(qv_sat_l) - qv_sat_i), 0.0)
        * p3_wbf_coeff, 0.0)

    qi_safe = jnp.where(active, qi_incld, 1.0)
    ni_sublim = jnp.where(active,
                          qisub * (jnp.asarray(ni_incld) / qi_safe), 0.0)
    return qidep, qisub, ni_sublim, qiberg


def calc_ice_relaxation_timescale(rho, t_atm, rhofaci,
                                  table_val_qi2qr_melting,
                                  table_val_qi2qr_vent_melt, dv, mu, sc,
                                  qi_incld, ni_incld, epsi_tot, context):
    """Inverse ice supersaturation relaxation timescale
    (calc_ice_relaxation_timescale). Returns (epsi, updated epsi_tot)."""
    t_atm = jnp.asarray(t_atm)
    active = ((jnp.asarray(qi_incld) >= c.qsmall)
              & (t_atm < c.T_zerodegc) & context)

    epsi = jnp.where(
        active,
        ((jnp.asarray(table_val_qi2qr_melting)
          + jnp.asarray(table_val_qi2qr_vent_melt)
          * _cbrt(jnp.asarray(sc))
          * jnp.sqrt(jnp.asarray(rhofaci) * jnp.asarray(rho)
                     / jnp.asarray(mu)))
         * 2.0 * c.pi * jnp.asarray(rho) * jnp.asarray(dv))
        * jnp.asarray(ni_incld),
        0.0)
    return epsi, jnp.where(active, jnp.asarray(epsi_tot) + epsi,
                           jnp.asarray(epsi_tot))


def calc_liq_relaxation_timescale(revap_table_vals, rho, dv, mu, sc,
                                  mu_r, lamr, cdistr, cdist,
                                  qr_incld, qc_incld, context):
    """Inverse liquid supersaturation relaxation timescales
    (calc_liq_relaxation_timescale). Returns (epsr, epsc)."""
    qr_incld = jnp.asarray(qr_incld)
    qc_incld = jnp.asarray(qc_incld)
    rho = jnp.asarray(rho)
    mu_r = jnp.asarray(mu_r)
    lamr = jnp.asarray(lamr)

    qr_ok = (qr_incld >= c.qsmall) & context
    tab = lookup_table3(mu_r, lamr, qr_ok)
    lamr_safe = jnp.where(qr_ok & (lamr > 0), lamr, 1.0)
    epsr = jnp.where(
        qr_ok,
        2.0 * c.pi * jnp.asarray(cdistr) * rho * jnp.asarray(dv)
        * (c.f1r * _tgamma(mu_r + 2.0) / lamr_safe
           + c.f2r * jnp.sqrt(rho / jnp.asarray(mu))
           * _cbrt(jnp.asarray(sc))
           * apply_table3(revap_table_vals, tab)),
        0.0)

    qc_ok = (qc_incld >= c.qsmall) & context
    epsc = jnp.where(qc_ok,
                     2.0 * c.pi * rho * jnp.asarray(dv) * jnp.asarray(cdist),
                     0.0)
    return epsr, epsc


def evaporate_rain(qr_incld, qc_incld, nr_incld, qi_incld, cld_frac_l,
                   cld_frac_r, qv, qv_prev, qv_sat_l, qv_sat_i, ab, abi,
                   epsr, epsi_tot, t_atm, t_prev, latent_heat_sublim,
                   dqsdt, dt, context):
    """Timestep-averaged rain evaporation, Morrison & Milbrandt (2015)
    (evaporate_rain). Returns (qr2qv_evap_tend, nr_evap_tend)."""
    qr_incld = jnp.asarray(qr_incld)
    qv = jnp.asarray(qv)
    qv_sat_l = jnp.asarray(qv_sat_l)
    t_atm = jnp.asarray(t_atm)
    cld_frac_r = jnp.asarray(cld_frac_r)
    epsr = jnp.asarray(epsr)
    ab = jnp.asarray(ab)

    inv_dt = 1.0 / dt
    ssat_r = qv - qv_sat_l
    cld_frac = jnp.where((jnp.asarray(qc_incld) + jnp.asarray(qi_incld))
                         < 1.0e-6,
                         0.0, jnp.asarray(cld_frac_l))

    is_evap = ((cld_frac_r > cld_frac) & (ssat_r < 0.0)
               & (qr_incld >= c.qsmall) & context)

    epsr_safe = jnp.where(is_evap & (epsr > 0), epsr, 1.0)
    tau_r = 1.0 / epsr_safe

    freezing = t_atm < 273.15
    lat_fac = 1.0 + latent_heat_sublim * c.inv_cp * jnp.asarray(dqsdt)
    eps_eff = jnp.where(freezing,
                        epsr + jnp.asarray(epsi_tot) * lat_fac
                        / jnp.asarray(abi),
                        epsr)
    a_c_common = ((qv - jnp.asarray(qv_prev)) * inv_dt
                  - jnp.asarray(dqsdt) * (t_atm - jnp.asarray(t_prev))
                  * inv_dt)
    a_c = jnp.where(freezing,
                    a_c_common - (qv_sat_l - jnp.asarray(qv_sat_i))
                    * lat_fac / jnp.asarray(abi) * jnp.asarray(epsi_tot),
                    a_c_common)

    eps_eff = jnp.maximum(1.0e-20, eps_eff)
    tau_eff = 1.0 / eps_eff

    qr_tiny = (qr_incld < 1e-12) & (qv / qv_sat_l < 0.999)

    # blended instantaneous/equilibrium tendency (helper subroutines
    # rain_evap_tscale_weight / _equilib_tend / _instant_tend inlined)
    dt_over_tau = dt / tau_eff
    tscale_weight = -jnp.expm1(-dt_over_tau) / dt_over_tau
    equilib_tend = -a_c / ab * tau_eff / tau_r
    instant_tend = -ssat_r / (ab * tau_r)
    blended = instant_tend * tscale_weight \
        + equilib_tend * (1.0 - tscale_weight)

    tend = jnp.where(qr_tiny, qr_incld * inv_dt, blended)
    tend = jnp.where(is_evap, tend, 0.0)

    # limiter sequence, in the Fortran order
    tend = jnp.where(is_evap, jnp.minimum(tend, -ssat_r * inv_dt / ab), tend)
    tend = jnp.where(is_evap, jnp.maximum(0.0, tend), tend)
    tend = jnp.where(is_evap, jnp.minimum(tend, qr_incld * inv_dt), tend)
    cfr_safe = jnp.where(cld_frac_r != 0, cld_frac_r, 1.0)
    tend = jnp.where(is_evap, tend * (cld_frac_r - cld_frac) / cfr_safe, tend)

    qr_safe = jnp.where(is_evap, qr_incld, 1.0)
    nr_tend = jnp.where(is_evap, tend * (jnp.asarray(nr_incld) / qr_safe), 0.0)
    return tend, nr_tend
