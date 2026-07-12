"""P3 ice-phase process rates and rain evaporation (batch 2).

Sources (components/eamxx/src/physics/p3/impl/):
  p3_ice_collection_impl.hpp (3 kernels), p3_ice_melting_impl.hpp,
  p3_ice_cldliq_wet_growth_impl.hpp, p3_ice_deposition_sublimation_impl.hpp,
  p3_ice_relaxation_timescale_impl.hpp,
  p3_calc_liq_relaxation_timescale_impl.hpp, p3_evaporate_rain_impl.hpp

Masked-set conventions as in processes_warm.py (zero defaults). Kernels
that accumulate into their outputs in the C++ (wet growth's shed terms,
epsi_tot) take the incoming array and return the updated value.
"""

import jax.numpy as jnp
from jax.scipy.special import gammaln

from ..foundation import constants as c
from ..foundation.saturation import SaturationFcn, qv_sat_dry
from .table_lookups import apply_table3, lookup_table3


def _tgamma(x):
    return jnp.exp(gammaln(x))


def ice_cldliq_collection(rho, temp, rhofaci, table_val_qc2qi_collect,
                          qi_incld, qc_incld, ni_incld, nc_incld,
                          opts, context):
    """Collection of cloud droplets by ice (riming); above freezing the
    collected water sheds to rain (Functions::ice_cldliq_collection).
    Returns (qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc).
    """
    temp = jnp.asarray(temp)
    qi_incld = jnp.asarray(qi_incld)
    qc_incld = jnp.asarray(qc_incld)

    freezing = temp <= c.Tmelt
    both = (qi_incld >= c.QSMALL) & (qc_incld >= c.QSMALL) & context

    rate_q = (jnp.asarray(rhofaci) * jnp.asarray(table_val_qc2qi_collect)
              * qc_incld * opts["cldliq_to_ice_collection_factor"]
              * jnp.asarray(rho) * jnp.asarray(ni_incld))
    rate_n = (jnp.asarray(rhofaci) * jnp.asarray(table_val_qc2qi_collect)
              * jnp.asarray(nc_incld) * opts["cldliq_to_ice_collection_factor"]
              * jnp.asarray(rho) * jnp.asarray(ni_incld))

    qccol = jnp.where(both & freezing, rate_q, 0.0)
    qcshd = jnp.where(both & ~freezing, rate_q, 0.0)
    nccol = jnp.where(both, rate_n, 0.0)
    ncshdc = jnp.where(both & ~freezing, qcshd * (c.ONE / c.dropmass), 0.0)
    return qccol, nccol, qcshd, ncshdc


def ice_rain_collection(rho, temp, rhofaci, logn0r, table_val_nr_collect,
                        table_val_qr2qi_collect, qi_incld, ni_incld, qr_incld,
                        opts, context):
    """Collection of rain by ice (Functions::ice_rain_collection).
    Returns (qr2qi_collect_tend, nr_collect_tend)."""
    temp = jnp.asarray(temp)
    both = (jnp.asarray(qi_incld) >= c.QSMALL) & (jnp.asarray(qr_incld) >= c.QSMALL) & context
    freezing = temp <= c.Tmelt
    f = opts["rain_to_ice_collection_factor"]

    qrcol = jnp.where(both & freezing,
                      10.0 ** (jnp.asarray(table_val_qr2qi_collect) + jnp.asarray(logn0r))
                      * jnp.asarray(rho) * jnp.asarray(rhofaci) * f
                      * jnp.asarray(ni_incld), 0.0)
    nrcol = jnp.where(both,
                      10.0 ** (jnp.asarray(table_val_nr_collect) + jnp.asarray(logn0r))
                      * jnp.asarray(rho) * jnp.asarray(rhofaci) * f
                      * jnp.asarray(ni_incld), 0.0)
    return qrcol, nrcol


def ice_self_collection(rho, rhofaci, table_val_ni_self_collect, eii,
                        qm_incld, qi_incld, ni_incld, context):
    """Ice-ice aggregation with rime-fraction-dependent sticking efficiency
    (Functions::ice_self_collection). Returns ni_selfcollect_tend."""
    qi_incld = jnp.asarray(qi_incld)
    qm_incld = jnp.asarray(qm_incld)
    qi_ok = (qi_incld >= c.QSMALL) & context
    qm_pos = (qm_incld > 0.0) & context

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
                     * jnp.asarray(ni_incld) ** 2, 0.0)


def ice_melting(rho, T_atm, pres, rhofaci, table_val_qi2qr_melting,
                table_val_qi2qr_vent_melt, dv, sc, mu, kap, qv,
                qi_incld, ni_incld, context):
    """Melting of ice to rain (Functions::ice_melting).
    Returns (qi2qr_melt_tend, ni2nr_melt_tend)."""
    T_atm = jnp.asarray(T_atm)
    qi_incld = jnp.asarray(qi_incld)
    active = (qi_incld >= c.QSMALL) & (T_atm > c.Tmelt) & context

    qsat0 = qv_sat_dry(jnp.full_like(T_atm, c.Tmelt), jnp.asarray(pres),
                       False, SaturationFcn.MURPHY_KOOP)
    melt = ((jnp.asarray(table_val_qi2qr_melting)
             + jnp.asarray(table_val_qi2qr_vent_melt) * jnp.cbrt(jnp.asarray(sc))
             * jnp.sqrt(jnp.asarray(rhofaci) * jnp.asarray(rho) / jnp.asarray(mu)))
            * ((T_atm - c.Tmelt) * jnp.asarray(kap)
               - jnp.asarray(rho) * c.LatVap * jnp.asarray(dv) * (qsat0 - jnp.asarray(qv)))
            * 2.0 * c.Pi / c.LatIce) * jnp.asarray(ni_incld)
    qimlt = jnp.maximum(jnp.where(active, melt, 0.0), 0.0)
    qi_safe = jnp.where(active, qi_incld, 1.0)
    nimlt = jnp.where(active, qimlt * (jnp.asarray(ni_incld) / qi_safe), 0.0)
    return qimlt, nimlt


def ice_cldliq_wet_growth(rho, temp, pres, rhofaci, table_val_qi2qr_melting,
                          table_val_qi2qr_vent_melt, dv, kap, mu, sc, qv,
                          qc_incld, qi_incld, ni_incld, qr_incld,
                          qr2qi_collect_tend, qc2qi_collect_tend,
                          nr_ice_shed_tend, qc2qr_ice_shed_tend, context):
    """Wet growth of rimed ice; excess collected water sheds
    (Functions::ice_cldliq_wet_growth). Modifies the collection/shed
    tendencies; returns (log_wetgrowth, qr2qi_collect_tend,
    qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend,
    qc2qr_ice_shed_tend)."""
    temp = jnp.asarray(temp)
    qi_incld = jnp.asarray(qi_incld)
    qccol = jnp.asarray(qc2qi_collect_tend)
    qrcol = jnp.asarray(qr2qi_collect_tend)
    nrshdr = jnp.asarray(nr_ice_shed_tend)
    qcshd = jnp.asarray(qc2qr_ice_shed_tend)

    any_if = ((qi_incld >= c.QSMALL)
              & ((jnp.asarray(qc_incld) + jnp.asarray(qr_incld)) >= 1.0e-6)
              & (temp < c.Tmelt) & context)
    any_if_col = any_if & ((qccol + qrcol) >= 1.0e-10)

    qsat0 = qv_sat_dry(jnp.full_like(temp, c.Tmelt), jnp.asarray(pres),
                       False, SaturationFcn.MURPHY_KOOP)
    growth = ((jnp.asarray(table_val_qi2qr_melting)
               + jnp.asarray(table_val_qi2qr_vent_melt) * jnp.cbrt(jnp.asarray(sc))
               * jnp.sqrt(jnp.asarray(rhofaci) * jnp.asarray(rho) / jnp.asarray(mu)))
              * 2.0 * c.Pi
              * (jnp.asarray(rho) * c.LatVap * jnp.asarray(dv) * (qsat0 - jnp.asarray(qv))
                 - (temp - c.Tmelt) * jnp.asarray(kap))
              / (c.LatIce + c.CpLiq * (temp - c.Tmelt))) * jnp.asarray(ni_incld)
    qc_growth_rate = jnp.where(any_if, jnp.maximum(growth, 0.0), 0.0)

    dum = jnp.maximum(0.0, (qccol + qrcol) - qc_growth_rate)
    dum_ok = (dum >= 1.0e-10) & context
    total_safe = jnp.where((qccol + qrcol) > 0, qccol + qrcol, 1.0)
    dum1 = 1.0 / total_safe

    nrshdr = jnp.where(any_if & dum_ok, nrshdr + dum * 1.923e6, nrshdr)
    qcshd = jnp.where(any_if_col & dum_ok, qcshd + dum * qccol * dum1, qcshd)
    qccol_new = jnp.where(any_if_col & dum_ok,
                          jnp.maximum(0.0, qccol - dum * qccol * dum1), qccol)
    qrcol_new = jnp.where(any_if_col & dum_ok,
                          jnp.maximum(0.0, qrcol - dum * qrcol * dum1), qrcol)
    log_wetgrowth = any_if & dum_ok
    return log_wetgrowth, qrcol_new, qccol_new, qc_growth_rate, nrshdr, qcshd


def ice_deposition_sublimation(qi_incld, ni_incld, T_atm, qv_sat_l, qv_sat_i,
                               epsi, abi, qv, inv_dt, context):
    """Vapor deposition / sublimation of ice and Bergeron rate
    (Functions::ice_deposition_sublimation). Returns
    (qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend)."""
    qi_incld = jnp.asarray(qi_incld)
    T_atm = jnp.asarray(T_atm)
    qv_sat_i = jnp.asarray(qv_sat_i)
    active = (qi_incld > c.QSMALL) & context

    qi_tend = jnp.minimum(jnp.asarray(epsi) / jnp.asarray(abi), inv_dt) \
        * (jnp.asarray(qv) - qv_sat_i)
    neg = qi_tend < 0.0
    freezing = T_atm < c.T_zerodegc

    qi2qv_sublim = jnp.where(active & neg, -qi_tend, 0.0)
    qi_safe = jnp.where(active, qi_incld, 1.0)
    ni_sublim = jnp.where(active & neg,
                          qi2qv_sublim * (jnp.asarray(ni_incld) / qi_safe), 0.0)
    qv2qi_vapdep = jnp.where(active & freezing & ~neg, qi_tend, 0.0)
    qc2qi_berg = jnp.where(active & freezing,
                           jnp.maximum(jnp.asarray(epsi) / jnp.asarray(abi)
                                       * (jnp.asarray(qv_sat_l) - qv_sat_i), 0.0),
                           0.0)
    return qv2qi_vapdep, qi2qv_sublim, ni_sublim, qc2qi_berg


def ice_relaxation_timescale(rho, temp, rhofaci, table_val_qi2qr_melting,
                             table_val_qi2qr_vent_melt, dv, mu, sc,
                             qi_incld, ni_incld, epsi_tot, context):
    """Inverse ice supersaturation relaxation timescale
    (Functions::ice_relaxation_timescale). Returns (epsi, epsi_tot updated)."""
    temp = jnp.asarray(temp)
    active = (jnp.asarray(qi_incld) >= c.QSMALL) & (temp < c.Tmelt) & context

    epsi = jnp.where(
        active,
        ((jnp.asarray(table_val_qi2qr_melting)
          + jnp.asarray(table_val_qi2qr_vent_melt) * jnp.cbrt(jnp.asarray(sc))
          * jnp.sqrt(jnp.asarray(rhofaci) * jnp.asarray(rho) / jnp.asarray(mu)))
         * 2.0 * c.Pi * jnp.asarray(rho) * jnp.asarray(dv)) * jnp.asarray(ni_incld),
        0.0)
    return epsi, jnp.where(active, jnp.asarray(epsi_tot) + epsi,
                           jnp.asarray(epsi_tot))


def calc_liq_relaxation_timescale(revap_table_vals, rho, dv, mu, sc,
                                  mu_r, lamr, cdistr, cdist,
                                  qr_incld, qc_incld, context):
    """Inverse liquid (rain/cloud) supersaturation relaxation timescales
    (Functions::calc_liq_relaxation_timescale). Returns (epsr, epsc)."""
    qr_incld = jnp.asarray(qr_incld)
    qc_incld = jnp.asarray(qc_incld)
    rho = jnp.asarray(rho)
    mu_r = jnp.asarray(mu_r)
    lamr = jnp.asarray(lamr)

    qr_ok = (qr_incld >= c.QSMALL) & context
    tab = lookup_table3(mu_r, lamr, qr_ok)
    lamr_safe = jnp.where(qr_ok & (lamr > 0), lamr, 1.0)
    epsr = jnp.where(
        qr_ok,
        2.0 * c.Pi * jnp.asarray(cdistr) * rho * jnp.asarray(dv)
        * (c.f1r * _tgamma(mu_r + 2.0) / lamr_safe
           + c.f2r * jnp.sqrt(rho / jnp.asarray(mu)) * jnp.cbrt(jnp.asarray(sc))
           * apply_table3(revap_table_vals, tab)),
        0.0)

    qc_ok = (qc_incld >= c.QSMALL) & context
    epsc = jnp.where(qc_ok, 2.0 * c.Pi * rho * jnp.asarray(dv) * jnp.asarray(cdist), 0.0)
    return epsr, epsc


def evaporate_rain(qr_incld, qc_incld, nr_incld, qi_incld, cld_frac_l,
                   cld_frac_r, qv, qv_prev, qv_sat_l, qv_sat_i, ab, abi,
                   epsr, epsi_tot, t_atm, t_atm_prev, dqsdt, dt, context):
    """Timestep-averaged rain evaporation, Morrison & Milbrandt (2015)
    (Functions::evaporate_rain). Returns (qr2qv_evap_tend, nr_evap_tend)."""
    qr_incld = jnp.asarray(qr_incld)
    qv = jnp.asarray(qv)
    qv_sat_l = jnp.asarray(qv_sat_l)
    t_atm = jnp.asarray(t_atm)
    cld_frac_r = jnp.asarray(cld_frac_r)
    epsr = jnp.asarray(epsr)
    ab = jnp.asarray(ab)

    inv_dt = 1.0 / dt
    inv_cp = 1.0 / c.Cpair

    ssat_r = qv - qv_sat_l
    cld_frac = jnp.where((jnp.asarray(qc_incld) + jnp.asarray(qi_incld)) < 1.0e-6,
                         0.0, jnp.asarray(cld_frac_l))

    is_evap = ((qr_incld >= c.QSMALL) & (ssat_r < 0.0)
               & (cld_frac_r > cld_frac) & context)

    epsr_safe = jnp.where(is_evap & (epsr > 0), epsr, 1.0)
    tau_r = 1.0 / epsr_safe

    freezing = t_atm < c.Tmelt
    lat_fac = (1.0 + (c.LatVap + c.LatIce) * inv_cp * jnp.asarray(dqsdt))
    eps_eff = jnp.where(freezing,
                        epsr + jnp.asarray(epsi_tot) * lat_fac / jnp.asarray(abi),
                        epsr)
    A_c_common = ((qv - jnp.asarray(qv_prev)) * inv_dt
                  - jnp.asarray(dqsdt) * (t_atm - jnp.asarray(t_atm_prev)) * inv_dt)
    A_c = jnp.where(freezing,
                    A_c_common - (qv_sat_l - jnp.asarray(qv_sat_i))
                    * lat_fac / jnp.asarray(abi) * jnp.asarray(epsi_tot),
                    A_c_common)

    eps_eff = jnp.maximum(eps_eff, 1e-20)
    tau_eff = 1.0 / eps_eff

    qr_tiny = (qr_incld < 1e-12) & (qv / qv_sat_l < 0.999)

    # Blended instantaneous/equilibrium tendency (helpers inlined)
    dt_over_tau = dt / tau_eff
    tscale_weight = -jnp.expm1(-dt_over_tau) / dt_over_tau
    equilib_tend = -A_c / ab * tau_eff / tau_r
    instant_tend = -ssat_r / (ab * tau_r)
    blended = instant_tend * tscale_weight + equilib_tend * (1.0 - tscale_weight)

    tend = jnp.where(qr_tiny, qr_incld * inv_dt, blended)
    tend = jnp.where(is_evap, tend, 0.0)

    # Clips, in the C++ order
    tend = jnp.where(is_evap & (tend > -ssat_r * inv_dt / ab),
                     -ssat_r * inv_dt / ab, tend)
    tend = jnp.where(is_evap & (tend < 0.0), 0.0, tend)
    tend = jnp.where(is_evap & (tend > qr_incld * inv_dt), qr_incld * inv_dt, tend)
    cfr_safe = jnp.where(cld_frac_r != 0, cld_frac_r, 1.0)
    tend = jnp.where(is_evap, tend * (cld_frac_r - cld_frac) / cfr_safe, tend)

    qr_safe = jnp.where(is_evap, qr_incld, 1.0)
    nr_tend = jnp.where(is_evap, tend * (jnp.asarray(nr_incld) / qr_safe), 0.0)
    return tend, nr_tend
