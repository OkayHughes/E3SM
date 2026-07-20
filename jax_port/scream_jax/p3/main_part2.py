"""p3_main_part2: the main microphysical process k-loop.

Source: components/eamxx/src/physics/p3/impl/p3_main_impl_part2.hpp

Orchestrates the ported process-rate kernels level-by-level (vectorized
over all levels/columns at once): DSD parameters -> ice-table lookups ->
collection/melting/wet-growth/freezing/evaporation/deposition rates ->
cell-average scaling -> conservation limiters -> supersaturation guards ->
prognostic updates -> final clipping and in-cloud recomputation.

The C++ per-pack skip_all/skip_micro early exits become masks; levels that
would be skipped keep their incoming values (all tendency writes are
masked by not_skip_all exactly as in the C++).
"""

import functools

import jax
import jax.numpy as jnp
from jax.scipy.special import gammaln

from ..foundation import constants as c
from . import family_width
from . import cell_average as ca
from . import conservation as cons
from . import processes_ice as pice
from . import processes_warm as pwarm
from .dsd import get_cloud_dsd2, get_rain_dsd2
from .main_part3 import calc_bulk_rho_rime
from .table_lookups import (
    apply_table_coll,
    apply_table_ice,
    lookup_ice,
    lookup_rain,
)
from .update import update_prognostic_ice, update_prognostic_liquid


def _tgamma(x):
    return jnp.exp(gammaln(x))


def get_cdistr_logn0r(qr, nr, mu_r, lamr, context):
    """Rain distribution intercept helpers (Functions::get_cdistr_logn0r)."""
    qr = jnp.asarray(qr)
    gt = (qr >= c.QSMALL) & context
    cdistr = jnp.where(gt, jnp.asarray(nr) / _tgamma(jnp.asarray(mu_r) + 1.0), 0.0)
    cdistr_safe = jnp.where(gt & (cdistr > 0), cdistr, 1.0)
    lamr_safe = jnp.where(gt & (jnp.asarray(lamr) > 0), jnp.asarray(lamr), 1.0)
    logn0r = jnp.where(gt, jnp.log10(cdistr_safe)
                       + (jnp.asarray(mu_r) + 1.0) * jnp.log10(lamr_safe), 0.0)
    return cdistr, logn0r


@functools.partial(jax.jit, static_argnames=(
    "predict_nc", "do_prescribed_ccn", "do_ice_production",
    "use_hetfrz_classnuc", "use_separate_ice_liq_frac",
    "smooth_width", "smooth_families"))
def p3_main_part2(predict_nc: bool, do_prescribed_ccn: bool,
                  do_ice_production: bool, use_hetfrz_classnuc: bool,
                  use_separate_ice_liq_frac: bool,
                  dt, max_total_ni,
                  hetfrz_immersion_nucleation_tend,
                  hetfrz_contact_nucleation_tend,
                  hetfrz_deposition_nucleation_tend,
                  tables, pres, dpres, dz, nc_nuceat_tend, inv_exner, exner,
                  inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
                  ni_activated, inv_qc_relvar,
                  cld_frac_i, cld_frac_l, cld_frac_r, qv_prev, t_prev,
                  st, opts, smooth_width=0.0, smooth_families=None):
    """One pass of the process k-loop.

    st is the state dict from p3_main_part1 (updated in place semantically;
    a new dict is returned) — must contain T_atm, qv, th_atm, the 9
    prognostics, rho/inv_rho/qv_sat_l/qv_sat_i/qv_supersat_i/rhofacr/
    rhofaci/acn and the *_incld arrays. Returns (new_state_dict,
    diagnostics_dict, hydrometeors_present).
    """
    inv_dt = 1.0 / dt
    g = dict(st)  # shallow copy; values replaced functionally

    # static (trace-time) per-family smoothing widths; 0.0 = exact path
    w_tmelt = family_width(smooth_width, smooth_families, "tmelt")
    w_frz = family_width(smooth_width, smooth_families, "frz")
    w_evap = family_width(smooth_width, smooth_families, "evap")
    w_rime = family_width(smooth_width, smooth_families, "rime")
    w_nucl = family_width(smooth_width, smooth_families, "nucl")

    T_atm = g["T_atm"]
    qv_supersat_i = g["qv_supersat_i"]
    qc, qr, qi = g["qc"], g["qr"], g["qi"]

    not_skip_all = ~((qc < c.QSMALL) & (qr < c.QSMALL) & (qi < c.QSMALL)
                     & (T_atm < c.T_zerodegc) & (qv_supersat_i < -0.05))
    not_skip_micro = not_skip_all & ((g["qc_incld"] >= c.QSMALL)
                                     | (g["qr_incld"] >= c.QSMALL)
                                     | (g["qi_incld"] >= c.QSMALL))

    zero = jnp.zeros_like(qc)
    tend_names = [
        "qc2qr_accret_tend", "qr2qv_evap_tend", "qc2qr_autoconv_tend",
        "nc_accret_tend", "nc_selfcollect_tend", "nc2nr_autoconv_tend",
        "nr_selfcollect_tend", "nr_evap_tend", "ncautr",
        "qi2qv_sublim_tend", "nr_ice_shed_tend", "qc2qi_hetero_freeze_tend",
        "qr2qi_collect_tend", "qc2qr_ice_shed_tend", "qi2qr_melt_tend",
        "qc2qi_collect_tend", "qr2qi_immers_freeze_tend", "qv2qi_nucleat_tend",
        "ni2nr_melt_tend", "nc_collect_tend", "ncshdc",
        "nc2ni_immers_freeze_tend", "nr_collect_tend", "ni_selfcollect_tend",
        "ni_nucleat_tend", "qv2qi_vapdep_tend", "qc2qi_berg_tend",
        "nr2ni_immers_freeze_tend", "ni_sublim_tend",
        "ncheti_cnt", "qcheti_cnt", "nicnt", "qicnt", "ninuc_cnt", "qinuc_cnt",
    ]
    t = {k: zero for k in tend_names}
    rho_qm_cloud = jnp.full_like(qc, 400.0)
    epsi_tot = zero
    wetgrowth = jnp.zeros_like(qc, dtype=bool)

    # ---------- micro-active section ----------
    mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii = \
        ca.get_time_space_phys_variables(
            T_atm, pres, g["rho"], g["qv_sat_l"], g["qv_sat_i"], not_skip_micro)

    nc_incld, mu_c, nu, lamc, cdist, cdist1 = get_cloud_dsd2(
        g["qc_incld"], g["nc_incld"], g["rho"], not_skip_micro)
    g["nc_incld"] = nc_incld
    g["nc"] = jnp.where(not_skip_micro, nc_incld * cld_frac_l, g["nc"])

    nr_incld, mu_r, lamr = get_rain_dsd2(
        g["qr_incld"], g["nr_incld"], opts["constant_mu_rain"], not_skip_micro)
    g["nr_incld"] = nr_incld
    cdistr, logn0r = get_cdistr_logn0r(g["qr_incld"], nr_incld, mu_r, lamr,
                                       not_skip_micro)
    g["nr"] = jnp.where(not_skip_micro, nr_incld * cld_frac_r, g["nr"])

    g["ni_incld"] = cons.impose_max_total_ni(g["ni_incld"], max_total_ni,
                                             g["inv_rho"], not_skip_micro)

    qi_gt = (g["qi_incld"] >= c.QSMALL) & not_skip_micro
    g["ni_incld"] = jnp.where(qi_gt, jnp.maximum(g["ni_incld"], c.NSMALL),
                              g["ni_incld"])
    g["nr_incld"] = jnp.where(qi_gt, jnp.maximum(g["nr_incld"], c.NSMALL),
                              g["nr_incld"])
    rhop, qm_incld, bm_incld = calc_bulk_rho_rime(
        g["qi_incld"], g["qm_incld"], g["bm_incld"], opts, qi_gt,
        smooth_width=w_rime)
    g["qm_incld"] = qm_incld
    g["bm_incld"] = bm_incld
    g["qm"] = jnp.where(qi_gt, qm_incld * cld_frac_i, g["qm"])
    g["bm"] = jnp.where(qi_gt, bm_incld * cld_frac_i, g["bm"])

    ti = lookup_ice(g["qi_incld"], g["ni_incld"], qm_incld, rhop, qi_gt)
    tr = lookup_rain(g["qr_incld"], g["nr_incld"], qi_gt)
    ice_tbl = tables["ice_table_vals"]
    coll_tbl = tables["collect_table_vals"]

    def tv(idx):
        return jnp.where(qi_gt, apply_table_ice(idx, ice_tbl, ti), 0.0)

    t_fallspd = tv(1)
    t_ni_selfc = tv(2)
    t_qc2qi = tv(3)
    t_melt = tv(4)
    t_lammax = tv(6)
    t_lammin = tv(7)
    t_vent = tv(9)
    qr_gt = (g["qr_incld"] >= c.QSMALL) & qi_gt
    t_nr_coll = jnp.where(qr_gt, apply_table_coll(0, coll_tbl, ti, tr), 0.0)
    t_qr2qi = jnp.where(qr_gt, apply_table_coll(1, coll_tbl, ti, tr), 0.0)

    g["ni_incld"] = jnp.where(qi_gt,
                              jnp.minimum(g["ni_incld"], t_lammax * g["ni_incld"]),
                              g["ni_incld"])
    g["ni_incld"] = jnp.where(qi_gt,
                              jnp.maximum(g["ni_incld"], t_lammin * g["ni_incld"]),
                              g["ni_incld"])

    if do_ice_production:
        (t["qc2qi_collect_tend"], t["nc_collect_tend"],
         t["qc2qr_ice_shed_tend"], t["ncshdc"]) = pice.ice_cldliq_collection(
            g["rho"], T_atm, g["rhofaci"], t_qc2qi, g["qi_incld"],
            g["qc_incld"], g["ni_incld"], g["nc_incld"], opts, not_skip_micro,
            smooth_width=w_tmelt)
        t["qr2qi_collect_tend"], t["nr_collect_tend"] = pice.ice_rain_collection(
            g["rho"], T_atm, g["rhofaci"], logn0r, t_nr_coll, t_qr2qi,
            g["qi_incld"], g["ni_incld"], g["qr_incld"], opts, not_skip_micro,
            smooth_width=w_tmelt)
        t["ni_selfcollect_tend"] = pice.ice_self_collection(
            g["rho"], g["rhofaci"], t_ni_selfc, eii, g["qm_incld"],
            g["qi_incld"], g["ni_incld"], not_skip_micro)

    t["qi2qr_melt_tend"], t["ni2nr_melt_tend"] = pice.ice_melting(
        g["rho"], T_atm, pres, g["rhofaci"], t_melt, t_vent, dv, sc, mu, kap,
        g["qv"], g["qi_incld"], g["ni_incld"], not_skip_micro,
        smooth_width=w_tmelt)

    if do_ice_production:
        (wetgrowth, t["qr2qi_collect_tend"], t["qc2qi_collect_tend"],
         _qc_growth, t["nr_ice_shed_tend"], t["qc2qr_ice_shed_tend"]) = \
            pice.ice_cldliq_wet_growth(
                g["rho"], T_atm, pres, g["rhofaci"], t_melt, t_vent, dv, kap,
                mu, sc, g["qv"], g["qc_incld"], g["qi_incld"], g["ni_incld"],
                g["qr_incld"], t["qr2qi_collect_tend"],
                t["qc2qi_collect_tend"], t["nr_ice_shed_tend"],
                t["qc2qr_ice_shed_tend"], not_skip_micro,
                smooth_width=w_tmelt)

    epsi, epsi_tot = pice.ice_relaxation_timescale(
        g["rho"], T_atm, g["rhofaci"], t_melt, t_vent, dv, mu, sc,
        g["qi_incld"], g["ni_incld"], epsi_tot, not_skip_micro)

    _vtrmi1, rho_qm_cloud = pwarm.calc_rime_density(
        T_atm, g["rhofaci"], t_fallspd, g["acn"], lamc, mu_c,
        g["qc_incld"], t["qc2qi_collect_tend"], not_skip_micro)

    if do_ice_production:
        if use_hetfrz_classnuc:
            t["ncheti_cnt"], t["qcheti_cnt"] = pwarm.ice_classical_nucleation(
                hetfrz_immersion_nucleation_tend, hetfrz_contact_nucleation_tend,
                hetfrz_deposition_nucleation_tend, g["rho"], g["qc_incld"],
                g["nc_incld"], 1)
            (t["nicnt"], t["qicnt"], t["ninuc_cnt"], t["qinuc_cnt"]) = \
                pwarm.ice_classical_nucleation(
                    hetfrz_immersion_nucleation_tend,
                    hetfrz_contact_nucleation_tend,
                    hetfrz_deposition_nucleation_tend, g["rho"],
                    g["qc_incld"], g["nc_incld"], 2)
        else:
            (t["qc2qi_hetero_freeze_tend"], t["nc2ni_immers_freeze_tend"]) = \
                pwarm.cldliq_immersion_freezing(
                    T_atm, lamc, mu_c, cdist1, g["qc_incld"], inv_qc_relvar,
                    opts, not_skip_micro, smooth_width=w_frz)
        (t["qr2qi_immers_freeze_tend"], t["nr2ni_immers_freeze_tend"]) = \
            pwarm.rain_immersion_freezing(
                T_atm, lamr, mu_r, cdistr, g["qr_incld"], opts, not_skip_micro,
                smooth_width=w_frz)

    epsr, epsc = pice.calc_liq_relaxation_timescale(
        tables["revap_table_vals"], g["rho"], dv, mu, sc, mu_r, lamr, cdistr,
        cdist, g["qr_incld"], g["qc_incld"], not_skip_micro)

    t["qr2qv_evap_tend"], t["nr_evap_tend"] = pice.evaporate_rain(
        g["qr_incld"], g["qc_incld"], g["nr_incld"], g["qi_incld"],
        cld_frac_l, cld_frac_r, g["qv"], qv_prev, g["qv_sat_l"],
        g["qv_sat_i"], ab, abi, epsr, epsi_tot, T_atm, t_prev, dqsdt, dt,
        not_skip_micro, smooth_width=w_evap)

    if do_ice_production:
        (t["qv2qi_vapdep_tend"], t["qi2qv_sublim_tend"], t["ni_sublim_tend"],
         t["qc2qi_berg_tend"]) = pice.ice_deposition_sublimation(
            g["qi_incld"], g["ni_incld"], T_atm, g["qv_sat_l"], g["qv_sat_i"],
            epsi, abi, g["qv"], inv_dt, not_skip_micro,
            smooth_width=w_tmelt)

    # ---------- level-active (not_skip_all) section ----------
    if do_ice_production:
        t["qv2qi_nucleat_tend"], t["ni_nucleat_tend"] = pwarm.ice_nucleation(
            T_atm, g["inv_rho"], g["ni"], ni_activated, qv_supersat_i, inv_dt,
            predict_nc, do_prescribed_ccn, opts, not_skip_all,
            smooth_width=w_nucl)

    (t["qc2qr_autoconv_tend"], t["nc2nr_autoconv_tend"], t["ncautr"]) = \
        pwarm.cloud_water_autoconversion(
            g["rho"], g["qc_incld"], g["nc_incld"], inv_qc_relvar, opts,
            not_skip_all)
    t["nc_selfcollect_tend"] = pwarm.droplet_self_collection(
        g["qc_incld"], not_skip_all)
    t["qc2qr_accret_tend"], t["nc_accret_tend"] = pwarm.cloud_rain_accretion(
        g["rho"], g["inv_rho"], g["qc_incld"], g["nc_incld"], g["qr_incld"],
        inv_qc_relvar, opts, not_skip_all)
    t["nr_selfcollect_tend"] = pwarm.rain_self_collection(
        g["rho"], g["qr_incld"], g["nr_incld"], opts, not_skip_all)

    t = ca.back_to_cell_average(cld_frac_l, cld_frac_r, cld_frac_i, t,
                                not_skip_all, use_separate_ice_liq_frac)

    # ---------- conservation ----------
    (t["qc2qr_autoconv_tend"], t["qc2qr_accret_tend"],
     t["qc2qi_collect_tend"], t["qc2qi_hetero_freeze_tend"],
     t["qc2qr_ice_shed_tend"], t["qc2qi_berg_tend"], t["qi2qv_sublim_tend"],
     t["qv2qi_vapdep_tend"], t["qcheti_cnt"], t["qicnt"]) = \
        cons.cloud_water_conservation(
            g["qc"], dt, t["qc2qr_autoconv_tend"], t["qc2qr_accret_tend"],
            t["qc2qi_collect_tend"], t["qc2qi_hetero_freeze_tend"],
            t["qc2qr_ice_shed_tend"], t["qc2qi_berg_tend"],
            t["qi2qv_sublim_tend"], t["qv2qi_vapdep_tend"],
            t["qcheti_cnt"], t["qicnt"], use_hetfrz_classnuc, not_skip_all,
            cld_frac_l, cld_frac_i, use_separate_ice_liq_frac)

    (t["qr2qv_evap_tend"], t["qr2qi_collect_tend"],
     t["qr2qi_immers_freeze_tend"]) = cons.rain_water_conservation(
        g["qr"], t["qc2qr_autoconv_tend"], t["qc2qr_accret_tend"],
        t["qi2qr_melt_tend"], t["qc2qr_ice_shed_tend"], dt,
        t["qr2qv_evap_tend"], t["qr2qi_collect_tend"],
        t["qr2qi_immers_freeze_tend"], not_skip_all)

    t["qi2qv_sublim_tend"], t["qi2qr_melt_tend"] = cons.ice_water_conservation(
        g["qi"], t["qv2qi_vapdep_tend"], t["qv2qi_nucleat_tend"],
        t["qc2qi_berg_tend"], t["qr2qi_collect_tend"],
        t["qc2qi_collect_tend"], t["qr2qi_immers_freeze_tend"],
        t["qc2qi_hetero_freeze_tend"], dt, t["qinuc_cnt"], t["qcheti_cnt"],
        t["qicnt"], t["qi2qv_sublim_tend"], t["qi2qr_melt_tend"],
        use_hetfrz_classnuc, not_skip_all)

    (t["nc_collect_tend"], t["nc2ni_immers_freeze_tend"],
     t["nc_accret_tend"], t["nc2nr_autoconv_tend"], t["ncheti_cnt"],
     t["nicnt"]) = cons.nc_conservation(
        g["nc"], t["nc_selfcollect_tend"], dt, t["nc_collect_tend"],
        t["nc2ni_immers_freeze_tend"], t["nc_accret_tend"],
        t["nc2nr_autoconv_tend"], t["ncheti_cnt"], t["nicnt"],
        use_hetfrz_classnuc, not_skip_all)

    (t["nr_collect_tend"], t["nr2ni_immers_freeze_tend"],
     t["nr_selfcollect_tend"], t["nr_evap_tend"]) = cons.nr_conservation(
        g["nr"], t["ni2nr_melt_tend"], t["nr_ice_shed_tend"], t["ncshdc"],
        t["nc2nr_autoconv_tend"], dt, c.nmltratio, t["nr_collect_tend"],
        t["nr2ni_immers_freeze_tend"], t["nr_selfcollect_tend"],
        t["nr_evap_tend"], not_skip_all)

    (t["ni2nr_melt_tend"], t["ni_sublim_tend"], t["ni_selfcollect_tend"]) = \
        cons.ni_conservation(
            g["ni"], t["ni_nucleat_tend"], t["nr2ni_immers_freeze_tend"],
            t["nc2ni_immers_freeze_tend"], t["ncheti_cnt"], t["nicnt"],
            t["ninuc_cnt"], dt, t["ni2nr_melt_tend"], t["ni_sublim_tend"],
            t["ni_selfcollect_tend"], use_hetfrz_classnuc, not_skip_all)

    (t["qv2qi_vapdep_tend"], t["qv2qi_nucleat_tend"], t["qinuc_cnt"]) = \
        cons.ice_supersat_conservation(
            t["qv2qi_vapdep_tend"], t["qv2qi_nucleat_tend"], t["qinuc_cnt"],
            cld_frac_i, g["qv"], g["qv_sat_i"], g["th_atm"] / inv_exner, dt,
            t["qi2qv_sublim_tend"], t["qr2qv_evap_tend"],
            use_hetfrz_classnuc, not_skip_all)

    t["qi2qv_sublim_tend"], t["qr2qv_evap_tend"] = \
        cons.prevent_liq_supersaturation(
            pres, T_atm, g["qv"], dt, t["qv2qi_vapdep_tend"],
            t["qv2qi_nucleat_tend"], t["qi2qv_sublim_tend"],
            t["qr2qv_evap_tend"], not_skip_all)

    # ---------- prognostic updates ----------
    t["inv_exner"] = inv_exner
    ice_state = {k: g[k] for k in ("th_atm", "qv", "qi", "ni", "qm", "bm",
                                   "qc", "nc", "qr", "nr")}
    ice_out = update_prognostic_ice(t, predict_nc, wetgrowth, dt, c.nmltratio,
                                    rho_qm_cloud, ice_state,
                                    use_hetfrz_classnuc, not_skip_all)
    for k, v in ice_out.items():
        g[k] = v

    (g["th_atm"], g["qv"], g["qc"], g["nc"], g["qr"], g["nr"]) = \
        update_prognostic_liquid(
            t["qc2qr_accret_tend"], t["nc_accret_tend"],
            t["qc2qr_autoconv_tend"], t["nc2nr_autoconv_tend"], t["ncautr"],
            t["nc_selfcollect_tend"], t["qr2qv_evap_tend"], t["nr_evap_tend"],
            t["nr_selfcollect_tend"], predict_nc, do_prescribed_ccn,
            g["inv_rho"], inv_exner, dt, g["th_atm"], g["qv"], g["qc"],
            g["nc"], g["qr"], g["nr"], not_skip_all)

    # ---------- exchange diagnostics ----------
    dep_extra = t["qinuc_cnt"] if use_hetfrz_classnuc else 0.0
    diags = {
        "qv2qi_depos_tend": jnp.where(
            not_skip_all,
            t["qv2qi_vapdep_tend"] - t["qi2qv_sublim_tend"]
            + t["qv2qi_nucleat_tend"] + dep_extra, 0.0),
        "precip_total_tend": jnp.where(
            not_skip_all,
            t["qc2qr_accret_tend"] + t["qc2qr_autoconv_tend"]
            + t["qc2qr_ice_shed_tend"] + t["qc2qi_collect_tend"], 0.0),
        "nevapr": jnp.where(not_skip_all,
                            t["qi2qv_sublim_tend"] + t["qr2qv_evap_tend"], 0.0),
        "qr_evap_tend": jnp.where(not_skip_all, t["qr2qv_evap_tend"], 0.0),
        "vap_ice_exchange": jnp.where(
            not_skip_all,
            t["qv2qi_vapdep_tend"] - t["qi2qv_sublim_tend"]
            + t["qv2qi_nucleat_tend"], 0.0),
        "vap_liq_exchange": jnp.where(not_skip_all,
                                      -t["qr2qv_evap_tend"], 0.0),
        "liq_ice_exchange": jnp.where(
            not_skip_all,
            t["qc2qi_hetero_freeze_tend"] + t["qr2qi_immers_freeze_tend"]
            - t["qi2qr_melt_tend"] + t["qc2qi_berg_tend"]
            + t["qc2qi_collect_tend"] + t["qr2qi_collect_tend"], 0.0),
        "pratot": jnp.where(not_skip_all, t["qc2qr_accret_tend"], 0.0),
        "prctot": jnp.where(not_skip_all, t["qc2qr_autoconv_tend"], 0.0),
    }

    # ---------- final clipping (as in part1, but masked by activity) ----------
    for spec, latent in (("qc", c.LatVap), ("qr", c.LatVap)):
        small = (g[spec] < c.QSMALL) & not_skip_all
        g["qv"] = jnp.where(small, g["qv"] + g[spec], g["qv"])
        g["th_atm"] = jnp.where(
            small, g["th_atm"] - inv_exner * g[spec] * latent * c.INV_CP,
            g["th_atm"])
        n_name = "nc" if spec == "qc" else "nr"
        g[n_name] = jnp.where(small, 0.0, g[n_name])
        g[spec] = jnp.where(small, 0.0, g[spec])

    small = (g["qi"] < c.QSMALL) & not_skip_all
    g["qv"] = jnp.where(small, g["qv"] + g["qi"], g["qv"])
    g["th_atm"] = jnp.where(
        small, g["th_atm"] - inv_exner * g["qi"] * (c.LatVap + c.LatIce) * c.INV_CP,
        g["th_atm"])
    g["qi"] = jnp.where(small, 0.0, g["qi"])
    g["ni"] = jnp.where(small, 0.0, g["ni"])
    g["qm"] = jnp.where(small, 0.0, g["qm"])
    g["bm"] = jnp.where(small, 0.0, g["bm"])

    hydrometeors_present = jnp.any(
        (((g["qc"] >= c.QSMALL) | (g["qr"] >= c.QSMALL) | (g["qi"] >= c.QSMALL))
         & not_skip_all), axis=-1)

    ni_incld = jnp.where(not_skip_all, g["ni"] / cld_frac_i, g["ni_incld"])
    ni_incld = cons.impose_max_total_ni(ni_incld, max_total_ni, g["inv_rho"],
                                        not_skip_all)
    g["ni"] = jnp.where(not_skip_all, ni_incld * cld_frac_i, g["ni"])
    g["ni_incld"] = ni_incld

    incld = cons.calculate_incloud_mixingratios(
        g["qc"], g["qr"], g["qi"], g["qm"], g["nc"], g["nr"], g["ni"], g["bm"],
        inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, not_skip_all)
    for name, val in zip(("qc_incld", "qr_incld", "qi_incld", "qm_incld",
                          "nc_incld", "nr_incld", "ni_incld", "bm_incld"),
                         incld):
        g[name] = jnp.where(not_skip_all, val, g[name])

    # DSD params needed downstream (part3/sedimentation reuse workspaces)
    g["mu_c"], g["nu"], g["lamc"] = mu_c, nu, lamc
    g["cdist"], g["cdist1"], g["cdistr"] = cdist, cdist1, cdistr
    g["mu_r"], g["lamr"], g["logn0r"] = mu_r, lamr, logn0r

    return g, diags, hydrometeors_present
