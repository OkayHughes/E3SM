"""p3_main_part2: the main microphysical process k-loop.

Source: micro_p3.F90 p3_main_part2 (eam variant). COPIED from
scream_jax/p3/main_part2.py; adaptations (see package docstring):
EAM tuning parameters via opts, subgrid variance scaling active,
CNT_couple threshold, EAM dsd/dep-sublim/nucleation forms, the EAM
conservation sequence ending in prevent_ice_overdepletion ->
ice_supersat_conservation, qi_wetDepos-based precip_total_tend, and
the p3_tend_out slot recording.

The Fortran per-level goto 555 (skip_all) / goto 444 (skip in-cloud
micro, keep nucleation+autoconv section) become masks; masked levels
keep their incoming values exactly as in the Fortran.
"""

import functools

import jax
import jax.numpy as jnp

from . import cell_average as ca
from . import conservation as cons
from . import constants as c
from . import processes_ice as pice
from . import processes_warm as pwarm
from .dsd import calc_bulk_rho_rime, get_cloud_dsd2, get_rain_dsd2
from .table_lookups import (
    apply_table_coll,
    apply_table_ice,
    lookup_ice,
    lookup_rain,
)
from .update import update_prognostic_ice, update_prognostic_liquid

# p3_tend_out slot -> tendency name (1-based Fortran slots; slots 1,
# 8, 9, 10, 12, 34 stay zero; 36-49 are filled by p3_main)
TEND_SLOTS = {
    2: "qc2qr_accret_tend", 3: "qc2qr_autoconv_tend",
    4: "nc_accret_tend", 5: "nc2nr_autoconv_tend",
    6: "nc_selfcollect_tend", 7: "nr_selfcollect_tend",
    11: "qr2qv_evap_tend", 13: "nr_evap_tend", 14: "ncautr",
    15: "qccol", 16: "qwgrth", 17: "qidep", 18: "qrcol", 19: "qinuc",
    20: "nc_collect_tend", 21: "nr_collect_tend",
    22: "ni_nucleat_tend", 23: "qi2qv_sublim_tend",
    24: "qi2qr_melt_tend", 25: "ni2nr_melt_tend", 26: "ni_sublim_tend",
    27: "ni_selfcollect_tend", 28: "qc2qi_hetero_freeze_tend",
    29: "qr2qi_immers_freeze_tend", 30: "nc2ni_immers_freeze_tend",
    31: "nr2ni_immers_freeze_tend", 32: "nr_ice_shed_tend",
    33: "qc2qr_ice_shed_tend", 35: "ncshdc",
}


@functools.partial(jax.jit, static_argnames=(
    "do_predict_nc", "do_prescribed_ccn", "use_hetfrz_classnuc",
    "do_cooper", "do_precip_off"))
def p3_main_part2(do_predict_nc: bool, do_prescribed_ccn: bool,
                  use_hetfrz_classnuc: bool, do_cooper: bool,
                  do_precip_off: bool, dt,
                  tables, pres, dpres, dz, nc_nuceat_tend, exner, inv_exner,
                  inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
                  ni_activated, inv_qc_relvar,
                  cld_frac_i, cld_frac_l, cld_frac_r, qv_prev, t_prev,
                  frzimm, frzcnt, frzdep, st, opts):
    """One pass of the process k-loop. st is the state dict from
    p3_main_part1. Returns (new_state_dict, diagnostics_dict,
    tendencies_dict, hydrometeors_present)."""
    inv_dt = 1.0 / dt
    g = dict(st)  # shallow copy; values replaced functionally

    t_atm = g["t_atm"]
    qv_supersat_i = g["qv_supersat_i"]
    qc, qr, qi = g["qc"], g["qr"], g["qi"]

    not_skip_all = ~((qc < c.qsmall) & (qr < c.qsmall) & (qi < c.qsmall)
                     & (t_atm < c.T_zerodegc) & (qv_supersat_i < -0.05))
    not_skip_micro = not_skip_all & ((g["qc_incld"] >= c.qsmall)
                                     | (g["qr_incld"] >= c.qsmall)
                                     | (g["qi_incld"] >= c.qsmall))

    zero = jnp.zeros_like(qc)
    tend_names = [
        "qc2qr_accret_tend", "qr2qv_evap_tend", "qc2qr_autoconv_tend",
        "nc_accret_tend", "nc_selfcollect_tend", "nc2nr_autoconv_tend",
        "nr_selfcollect_tend", "nr_evap_tend", "ncautr",
        "qi2qv_sublim_tend", "nr_ice_shed_tend", "qc2qi_hetero_freeze_tend",
        "qrcol", "qc2qr_ice_shed_tend", "qi2qr_melt_tend",
        "qccol", "qr2qi_immers_freeze_tend", "qinuc",
        "ni2nr_melt_tend", "nc_collect_tend", "ncshdc",
        "nc2ni_immers_freeze_tend", "nr_collect_tend",
        "ni_selfcollect_tend", "ni_nucleat_tend", "qidep", "qiberg",
        "nr2ni_immers_freeze_tend", "ni_sublim_tend", "qwgrth",
        "ncheti_cnt", "qcheti_cnt", "nicnt", "qicnt", "ninuc_cnt",
        "qinuc_cnt",
    ]
    t = {k: zero for k in tend_names}
    epsi_tot = zero
    wetgrowth = jnp.zeros_like(qc, dtype=bool)

    # ---------- in-cloud micro section (skipped by goto 444) ----------
    mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii = \
        ca.get_time_space_phys_variables(
            t_atm, pres, g["rho"], g["qv_sat_l"], g["qv_sat_i"],
            not_skip_micro)

    nc_incld, mu_c, nu, lamc, cdist, cdist1 = get_cloud_dsd2(
        g["qc_incld"], g["nc_incld"], g["rho"], not_skip_micro)
    g["nc_incld"] = jnp.where(not_skip_micro, nc_incld, g["nc_incld"])
    g["nc"] = jnp.where(not_skip_micro, nc_incld * cld_frac_l, g["nc"])

    nr_incld, mu_r, lamr, cdistr, logn0r = get_rain_dsd2(
        g["qr_incld"], g["nr_incld"], opts["p3_max_mean_rain_size"],
        not_skip_micro)
    g["nr_incld"] = jnp.where(not_skip_micro, nr_incld, g["nr_incld"])
    g["nr"] = jnp.where(not_skip_micro, nr_incld * cld_frac_r, g["nr"])

    g["ni_incld"] = cons.impose_max_total_ni(
        g["ni_incld"], c.max_total_ni, g["inv_rho"], not_skip_micro)

    qi_gt = (g["qi_incld"] >= c.qsmall) & not_skip_micro
    g["ni_incld"] = jnp.where(qi_gt, jnp.maximum(g["ni_incld"], c.nsmall),
                              g["ni_incld"])
    g["nr_incld"] = jnp.where(qi_gt, jnp.maximum(g["nr_incld"], c.nsmall),
                              g["nr_incld"])
    rhop, qm_incld, bm_incld = calc_bulk_rho_rime(
        g["qi_incld"], g["qm_incld"], g["bm_incld"], qi_gt)
    g["qm_incld"] = jnp.where(qi_gt, qm_incld, g["qm_incld"])
    g["bm_incld"] = jnp.where(qi_gt, bm_incld, g["bm_incld"])
    g["qm"] = jnp.where(qi_gt, qm_incld * cld_frac_i, g["qm"])
    g["bm"] = jnp.where(qi_gt, bm_incld * cld_frac_i, g["bm"])

    ti = lookup_ice(g["qi_incld"], g["ni_incld"], g["qm_incld"], rhop, qi_gt)
    tr = lookup_rain(g["qr_incld"], g["nr_incld"], qi_gt)
    ice_tbl = tables["ice_table_vals"]
    coll_tbl = tables["collect_table_vals"]

    def tv(idx):
        return jnp.where(qi_gt, apply_table_ice(idx, ice_tbl, ti), 0.0)

    t_fallspd = tv(1)      # Fortran index 2
    t_ni_selfc = tv(2)     # 3
    t_qc2qi = tv(3)        # 4
    t_melt = tv(4)         # 5
    t_lammax = tv(6)       # 7
    t_lammin = tv(7)       # 8
    t_vent = tv(9)         # 10
    qr_gt = (g["qr_incld"] >= c.qsmall) & qi_gt
    t_nr_coll = jnp.where(qr_gt, apply_table_coll(0, coll_tbl, ti, tr), 0.0)
    t_qr2qi = jnp.where(qr_gt, apply_table_coll(1, coll_tbl, ti, tr), 0.0)

    g["ni_incld"] = jnp.where(
        qi_gt, jnp.minimum(g["ni_incld"], t_lammax * g["ni_incld"]),
        g["ni_incld"])
    g["ni_incld"] = jnp.where(
        qi_gt, jnp.maximum(g["ni_incld"], t_lammin * g["ni_incld"]),
        g["ni_incld"])

    (t["qccol"], t["nc_collect_tend"],
     t["qc2qr_ice_shed_tend"], t["ncshdc"]) = pice.ice_cldliq_collection(
        g["rho"], t_atm, g["rhofaci"], t_qc2qi, g["qi_incld"],
        g["qc_incld"], g["ni_incld"], g["nc_incld"], not_skip_micro)

    t["qrcol"], t["nr_collect_tend"] = pice.ice_rain_collection(
        g["rho"], t_atm, g["rhofaci"], logn0r, t_nr_coll, t_qr2qi,
        g["qi_incld"], g["ni_incld"], g["qr_incld"], not_skip_micro)

    t["ni_selfcollect_tend"] = pice.ice_self_collection(
        g["rho"], g["rhofaci"], t_ni_selfc, eii, g["qm_incld"],
        g["qi_incld"], g["ni_incld"], not_skip_micro)

    t["qi2qr_melt_tend"], t["ni2nr_melt_tend"] = pice.ice_melting(
        g["rho"], t_atm, pres, g["rhofaci"], t_melt, t_vent,
        c.latvap, c.latice, dv, sc, mu, kap,
        g["qv"], g["qi_incld"], g["ni_incld"], not_skip_micro)

    (wetgrowth, t["qrcol"], t["qccol"], t["qwgrth"],
     t["nr_ice_shed_tend"], t["qc2qr_ice_shed_tend"]) = \
        pice.ice_cldliq_wet_growth(
            g["rho"], t_atm, pres, g["rhofaci"], t_melt, t_vent,
            c.latvap, c.latice, dv, kap, mu, sc, g["qv"], g["qc_incld"],
            g["qi_incld"], g["ni_incld"], g["qr_incld"], t["qrcol"],
            t["qccol"], t["nr_ice_shed_tend"], t["qc2qr_ice_shed_tend"],
            not_skip_micro)

    epsi, epsi_tot = pice.calc_ice_relaxation_timescale(
        g["rho"], t_atm, g["rhofaci"], t_melt, t_vent, dv, mu, sc,
        g["qi_incld"], g["ni_incld"], epsi_tot, not_skip_micro)

    _vtrmi1, rho_qm_cloud = pwarm.calc_rime_density(
        t_atm, g["rhofaci"], t_fallspd, g["acn"], lamc, mu_c,
        g["qc_incld"], t["qccol"], not_skip_micro)

    if use_hetfrz_classnuc:
        t["ncheti_cnt"], t["qcheti_cnt"] = pwarm.cnt_couple(
            frzimm, frzcnt, frzdep, g["rho"], g["qc_incld"],
            g["nc_incld"], 1, not_skip_micro)
        (t["nicnt"], t["qicnt"], t["ninuc_cnt"], t["qinuc_cnt"]) = \
            pwarm.cnt_couple(frzimm, frzcnt, frzdep, g["rho"],
                             g["qc_incld"], g["nc_incld"], 2,
                             not_skip_micro)
    else:
        (t["qc2qi_hetero_freeze_tend"], t["nc2ni_immers_freeze_tend"]) = \
            pwarm.cldliq_immersion_freezing(
                t_atm, lamc, mu_c, cdist1, g["qc_incld"], inv_qc_relvar,
                not_skip_micro)

    (t["qr2qi_immers_freeze_tend"], t["nr2ni_immers_freeze_tend"]) = \
        pwarm.rain_immersion_freezing(
            t_atm, lamr, mu_r, cdistr, g["qr_incld"], not_skip_micro)

    epsr, epsc = pice.calc_liq_relaxation_timescale(
        tables["revap_table_vals"], g["rho"], dv, mu, sc, mu_r, lamr,
        cdistr, cdist, g["qr_incld"], g["qc_incld"], not_skip_micro)

    t["qr2qv_evap_tend"], t["nr_evap_tend"] = pice.evaporate_rain(
        g["qr_incld"], g["qc_incld"], g["nr_incld"], g["qi_incld"],
        cld_frac_l, cld_frac_r, g["qv"], qv_prev, g["qv_sat_l"],
        g["qv_sat_i"], ab, abi, epsr, epsi_tot, t_atm, t_prev,
        c.latsub, dqsdt, dt, not_skip_micro)

    (t["qidep"], t["qi2qv_sublim_tend"], t["ni_sublim_tend"],
     t["qiberg"]) = pice.ice_deposition_sublimation(
        g["qi_incld"], g["ni_incld"], t_atm, g["qv_sat_l"], g["qv_sat_i"],
        epsi, abi, g["qv"], opts["p3_wbf_coeff"], not_skip_micro)

    # ---------- level-active (not_skip_all) section (label 444) ----------
    t["qinuc"], t["ni_nucleat_tend"] = pwarm.ice_nucleation(
        t_atm, g["inv_rho"], g["ni"], ni_activated, qv_supersat_i, inv_dt,
        do_predict_nc, do_prescribed_ccn, do_cooper, not_skip_all)

    (t["qc2qr_autoconv_tend"], t["nc2nr_autoconv_tend"], t["ncautr"]) = \
        pwarm.cloud_water_autoconversion(
            g["rho"], g["qc_incld"], g["nc_incld"], inv_qc_relvar, opts,
            do_precip_off, not_skip_all)
    t["nc_selfcollect_tend"] = pwarm.droplet_self_collection(
        g["qc_incld"], not_skip_all)
    t["qc2qr_accret_tend"], t["nc_accret_tend"] = pwarm.cloud_rain_accretion(
        g["rho"], g["inv_rho"], g["qc_incld"], g["nc_incld"], g["qr_incld"],
        inv_qc_relvar, opts, not_skip_all)
    t["nr_selfcollect_tend"] = pwarm.rain_self_collection(
        g["rho"], g["qr_incld"], g["nr_incld"], not_skip_all)

    qwgrth_saved = t.pop("qwgrth")  # diagnostic only; not cell-averaged
    t = ca.back_to_cell_average(cld_frac_l, cld_frac_r, cld_frac_i, t,
                                not_skip_all)
    t["qwgrth"] = qwgrth_saved

    # ---------- conservation of water (EAM sequence) ----------
    (t["qc2qr_autoconv_tend"], t["qc2qr_accret_tend"], t["qccol"],
     t["qc2qi_hetero_freeze_tend"], t["qc2qr_ice_shed_tend"], t["qiberg"],
     t["qi2qv_sublim_tend"], t["qidep"], t["qcheti_cnt"], t["qicnt"]) = \
        cons.cloud_water_conservation(
            g["qc"], dt, t["qc2qr_autoconv_tend"], t["qc2qr_accret_tend"],
            t["qccol"], t["qc2qi_hetero_freeze_tend"],
            t["qc2qr_ice_shed_tend"], t["qiberg"], t["qi2qv_sublim_tend"],
            t["qidep"], t["qcheti_cnt"], t["qicnt"],
            use_hetfrz_classnuc, not_skip_all)

    (t["qr2qv_evap_tend"], t["qrcol"], t["qr2qi_immers_freeze_tend"]) = \
        cons.rain_water_conservation(
            g["qr"], t["qc2qr_autoconv_tend"], t["qc2qr_accret_tend"],
            t["qi2qr_melt_tend"], t["qc2qr_ice_shed_tend"], dt,
            t["qr2qv_evap_tend"], t["qrcol"],
            t["qr2qi_immers_freeze_tend"], not_skip_all)

    t["qi2qv_sublim_tend"], t["qi2qr_melt_tend"] = \
        cons.ice_water_conservation(
            g["qi"], t["qidep"], t["qinuc"], t["qiberg"], t["qrcol"],
            t["qccol"], t["qr2qi_immers_freeze_tend"],
            t["qc2qi_hetero_freeze_tend"], dt, t["qinuc_cnt"],
            t["qcheti_cnt"], t["qicnt"], t["qi2qv_sublim_tend"],
            t["qi2qr_melt_tend"], use_hetfrz_classnuc, not_skip_all)

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

    t["qi2qv_sublim_tend"], t["qr2qv_evap_tend"] = \
        cons.prevent_ice_overdepletion(
            pres, t_atm, g["qv"], c.latvap, c.latsub, inv_dt, dt,
            t["qidep"], t["qinuc"], t["qinuc_cnt"],
            t["qi2qv_sublim_tend"], t["qr2qv_evap_tend"], not_skip_all)

    (t["qidep"], t["qinuc"], t["qinuc_cnt"]) = \
        cons.ice_supersat_conservation(
            t["qidep"], t["qinuc"], t["qi2qv_sublim_tend"],
            t["qr2qv_evap_tend"], t["qinuc_cnt"], cld_frac_i, g["qv"],
            g["qv_sat_i"], c.latsub, g["th"] / exner, dt,
            use_hetfrz_classnuc, not_skip_all)

    # ---------- prognostic updates ----------
    t["exner"] = exner
    ice_state = {k: g[k] for k in ("th", "qv", "qi", "ni", "qm", "bm",
                                   "qc", "nc", "qr", "nr")}
    ice_out, qi_wet_depos = update_prognostic_ice(
        t, do_predict_nc, wetgrowth, dt, c.nmltratio, rho_qm_cloud,
        ice_state, use_hetfrz_classnuc, not_skip_all)
    for k, v in ice_out.items():
        g[k] = v

    (g["th"], g["qv"], g["qc"], g["nc"], g["qr"], g["nr"]) = \
        update_prognostic_liquid(
            t["qc2qr_accret_tend"], t["nc_accret_tend"],
            t["qc2qr_autoconv_tend"], t["nc2nr_autoconv_tend"], t["ncautr"],
            t["nc_selfcollect_tend"], t["qr2qv_evap_tend"],
            t["nr_evap_tend"], t["nr_selfcollect_tend"], do_predict_nc,
            opts["nccnst"], do_prescribed_ccn, g["inv_rho"], exner,
            c.latvap, dt, g["th"], g["qv"], g["qc"], g["nc"], g["qr"],
            g["nr"], not_skip_all)

    # ---------- exchange diagnostics ----------
    dep_extra = t["qinuc_cnt"] if use_hetfrz_classnuc else 0.0
    diags = {
        "qv2qi_depos_tend": jnp.where(
            not_skip_all,
            t["qidep"] - t["qi2qv_sublim_tend"] + t["qinuc"] + dep_extra,
            0.0),
        "precip_total_tend": jnp.where(
            not_skip_all,
            t["qc2qr_accret_tend"] + t["qc2qr_autoconv_tend"]
            + t["qc2qr_ice_shed_tend"] + qi_wet_depos, 0.0),
        "nevapr": jnp.where(not_skip_all,
                            t["qi2qv_sublim_tend"] + t["qr2qv_evap_tend"],
                            0.0),
        "qr_evap_tend": jnp.where(not_skip_all, t["qr2qv_evap_tend"], 0.0),
        "vap_ice_exchange": jnp.where(
            not_skip_all,
            t["qidep"] - t["qi2qv_sublim_tend"] + t["qinuc"], 0.0),
        "vap_liq_exchange": jnp.where(not_skip_all,
                                      -t["qr2qv_evap_tend"], 0.0),
        "liq_ice_exchange": jnp.where(
            not_skip_all,
            t["qc2qi_hetero_freeze_tend"] + t["qr2qi_immers_freeze_tend"]
            - t["qi2qr_melt_tend"] + t["qiberg"] + t["qccol"] + t["qrcol"],
            0.0),
        "pratot": jnp.where(not_skip_all, t["qc2qr_accret_tend"], 0.0),
        "prctot": jnp.where(not_skip_all, t["qc2qr_autoconv_tend"], 0.0),
    }

    # ---------- final clipping (masked by activity) ----------
    for spec, n_name in (("qc", "nc"), ("qr", "nr")):
        small = (g[spec] < c.qsmall) & not_skip_all
        g["qv"] = jnp.where(small, g["qv"] + g[spec], g["qv"])
        g["th"] = jnp.where(
            small, g["th"] - exner * g[spec] * c.latvap * c.inv_cp, g["th"])
        g[n_name] = jnp.where(small, 0.0, g[n_name])
        g[spec] = jnp.where(small, 0.0, g[spec])

    small = (g["qi"] < c.qsmall) & not_skip_all
    g["qv"] = jnp.where(small, g["qv"] + g["qi"], g["qv"])
    g["th"] = jnp.where(small,
                        g["th"] - exner * g["qi"] * c.latsub * c.inv_cp,
                        g["th"])
    g["qi"] = jnp.where(small, 0.0, g["qi"])
    g["ni"] = jnp.where(small, 0.0, g["ni"])
    g["qm"] = jnp.where(small, 0.0, g["qm"])
    g["bm"] = jnp.where(small, 0.0, g["bm"])

    hydrometeors_present = jnp.any(
        (((g["qc"] >= c.qsmall) | (g["qr"] >= c.qsmall)
          | (g["qi"] >= c.qsmall)) & not_skip_all), axis=-1)

    ni_incld = jnp.where(not_skip_all, g["ni"] / cld_frac_i, g["ni_incld"])
    ni_incld = cons.impose_max_total_ni(ni_incld, c.max_total_ni,
                                        g["inv_rho"], not_skip_all)
    g["ni"] = jnp.where(not_skip_all, ni_incld * cld_frac_i, g["ni"])
    g["ni_incld"] = ni_incld

    incld = cons.calculate_incloud_mixingratios(
        g["qc"], g["qr"], g["qi"], g["qm"], g["nc"], g["nr"], g["ni"],
        g["bm"], inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
        not_skip_all)
    for name, val in zip(("qc_incld", "qr_incld", "qi_incld", "qm_incld",
                          "nc_incld", "nr_incld", "ni_incld", "bm_incld"),
                         incld):
        g[name] = jnp.where(not_skip_all, val, g[name])

    # DSD params needed downstream (sedimentation/part3 reuse them)
    g["mu_c"], g["nu"], g["lamc"] = mu_c, nu, lamc
    g["cdist"], g["cdist1"], g["cdistr"] = cdist, cdist1, cdistr
    g["mu_r"], g["lamr"], g["logn0r"] = mu_r, lamr, logn0r

    # p3_tend_out slots recorded by part2 (masked like the Fortran:
    # skipped levels keep zero)
    tend_out = {}
    for slot, name in TEND_SLOTS.items():
        tend_out[slot] = jnp.where(not_skip_all, t[name], 0.0)
    tend_out[16] = jnp.where(not_skip_all, t["qwgrth"], 0.0)

    return g, diags, tend_out, hydrometeors_present
