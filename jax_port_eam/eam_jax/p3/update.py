"""Prognostic-state updates from process rates.

Source: micro_p3.F90 update_prognostic_ice / update_prognostic_liquid
(eam variant). COPIED from scream_jax/p3/update.py; adaptations:

  * EAM naming: the theta factor is `exner` (multiplies T).
  * update_prognostic_ice additionally computes qi_wetDepos (the total
    ice process rate used by EAM's precip_total_tend).
  * update_prognostic_liquid takes nccnst (runtime) for the
    specified-Nc branch.
"""

import jax.numpy as jnp

from . import constants as c


def update_prognostic_ice(t, do_predict_nc: bool, log_wetgrowth, dt,
                          nmltratio, rho_qm_cloud, state,
                          use_hetfrz_classnuc: bool, context):
    """Apply ice-phase tendencies (update_prognostic_ice).

    t: dict of tendencies (Fortran names + 'exner'); state: dict with
    th, qv, qi, ni, qm, bm, qc, nc, qr, nr. Returns (updated state
    dict, qi_wetDepos)."""
    s = {k: jnp.asarray(v) for k, v in state.items()}
    g = lambda k: jnp.asarray(t[k])  # noqa: E731

    if use_hetfrz_classnuc:
        s["qc"] = jnp.where(context, s["qc"] + (
            -g("qcheti_cnt") - g("qicnt") - g("qccol")
            - g("qc2qr_ice_shed_tend") - g("qiberg")) * dt, s["qc"])
    else:
        s["qc"] = jnp.where(context, s["qc"] + (
            -g("qc2qi_hetero_freeze_tend") - g("qccol")
            - g("qc2qr_ice_shed_tend") - g("qiberg")) * dt, s["qc"])

    if do_predict_nc:
        if use_hetfrz_classnuc:
            s["nc"] = jnp.where(context, s["nc"] + (
                -g("nc_collect_tend") - g("ncheti_cnt") - g("nicnt")) * dt,
                s["nc"])
        else:
            s["nc"] = jnp.where(context, s["nc"] + (
                -g("nc_collect_tend") - g("nc2ni_immers_freeze_tend")) * dt,
                s["nc"])

    s["qr"] = jnp.where(context, s["qr"] + (
        -g("qrcol") + g("qi2qr_melt_tend")
        - g("qr2qi_immers_freeze_tend") + g("qc2qr_ice_shed_tend")) * dt,
        s["qr"])
    s["nr"] = jnp.where(context, s["nr"] + (
        -g("nr_collect_tend") - g("nr2ni_immers_freeze_tend")
        + nmltratio * g("ni2nr_melt_tend") + g("nr_ice_shed_tend")
        + g("ncshdc")) * dt, s["nr"])

    # sublimation/melting deplete bm/qm proportionally, then qi
    qi_ok = (s["qi"] >= c.qsmall) & context
    qi_safe = jnp.where(qi_ok, s["qi"], 1.0)
    subl_melt = g("qi2qv_sublim_tend") + g("qi2qr_melt_tend")
    s["bm"] = jnp.where(qi_ok, s["bm"] - (subl_melt / qi_safe) * dt * s["bm"],
                        s["bm"])
    s["qm"] = jnp.where(qi_ok, s["qm"] - (subl_melt * s["qm"] / qi_safe) * dt,
                        s["qm"])
    s["qi"] = jnp.where(qi_ok, s["qi"] - subl_melt * dt, s["qi"])

    if use_hetfrz_classnuc:
        dum = (g("qrcol") + g("qccol") + g("qr2qi_immers_freeze_tend")
               + g("qcheti_cnt") + g("qicnt")) * dt
        qi_src = (g("qidep") + g("qinuc") + g("qiberg")
                  + g("qinuc_cnt")) * dt
        qi_wet_depos = ((g("qidep") + g("qinuc") + g("qiberg")
                         + g("qinuc_cnt"))
                        + (g("qrcol") + g("qccol")
                           + g("qr2qi_immers_freeze_tend")
                           + g("qcheti_cnt") + g("qicnt")))
        bm_src = (g("qrcol") * c.inv_rho_rimeMax
                  + g("qccol") / jnp.asarray(rho_qm_cloud)
                  + (g("qr2qi_immers_freeze_tend") + g("qcheti_cnt")
                     + g("qicnt")) * c.inv_rho_rimeMax) * dt
        ni_src = (g("ni_nucleat_tend") - g("ni2nr_melt_tend")
                  - g("ni_sublim_tend") - g("ni_selfcollect_tend")
                  + g("nr2ni_immers_freeze_tend") + g("ncheti_cnt")
                  + g("nicnt") + g("ninuc_cnt")) * dt
    else:
        dum = (g("qrcol") + g("qccol") + g("qr2qi_immers_freeze_tend")
               + g("qc2qi_hetero_freeze_tend")) * dt
        qi_src = (g("qidep") + g("qinuc") + g("qiberg")) * dt
        qi_wet_depos = ((g("qidep") + g("qinuc") + g("qiberg"))
                        + (g("qrcol") + g("qccol")
                           + g("qr2qi_immers_freeze_tend")
                           + g("qc2qi_hetero_freeze_tend")))
        bm_src = (g("qrcol") * c.inv_rho_rimeMax
                  + g("qccol") / jnp.asarray(rho_qm_cloud)
                  + (g("qr2qi_immers_freeze_tend")
                     + g("qc2qi_hetero_freeze_tend"))
                  * c.inv_rho_rimeMax) * dt
        ni_src = (g("ni_nucleat_tend") - g("ni2nr_melt_tend")
                  - g("ni_sublim_tend") - g("ni_selfcollect_tend")
                  + g("nr2ni_immers_freeze_tend")
                  + g("nc2ni_immers_freeze_tend")) * dt

    s["qi"] = jnp.where(context, s["qi"] + qi_src + dum, s["qi"])
    qi_wet_depos = jnp.where(context, qi_wet_depos, 0.0)
    s["qm"] = jnp.where(context, s["qm"] + dum, s["qm"])
    s["bm"] = jnp.where(context, s["bm"] + bm_src, s["bm"])
    s["ni"] = jnp.where(context, s["ni"] + ni_src, s["ni"])

    neg_qm = (s["qm"] < 0.0) & context
    s["qm"] = jnp.where(neg_qm, 0.0, s["qm"])
    s["bm"] = jnp.where(neg_qm, 0.0, s["bm"])

    # densify under wet growth
    wet = jnp.asarray(log_wetgrowth) & context
    s["qm"] = jnp.where(wet, s["qi"], s["qm"])
    s["bm"] = jnp.where(wet, s["qm"] * c.inv_rho_rimeMax, s["bm"])

    if use_hetfrz_classnuc:
        qv_tend = (-g("qidep") + g("qi2qv_sublim_tend") - g("qinuc")
                   - g("qinuc_cnt"))
        th_sub = ((g("qidep") - g("qi2qv_sublim_tend") + g("qinuc")
                   + g("qinuc_cnt")) * c.latsub * c.inv_cp
                  + (g("qrcol") + g("qccol") + g("qcheti_cnt") + g("qicnt")
                     + g("qr2qi_immers_freeze_tend") - g("qi2qr_melt_tend")
                     + g("qiberg")) * c.latice * c.inv_cp)
    else:
        qv_tend = -g("qidep") + g("qi2qv_sublim_tend") - g("qinuc")
        th_sub = ((g("qidep") - g("qi2qv_sublim_tend") + g("qinuc"))
                  * c.latsub * c.inv_cp
                  + (g("qrcol") + g("qccol")
                     + g("qc2qi_hetero_freeze_tend")
                     + g("qr2qi_immers_freeze_tend") - g("qi2qr_melt_tend")
                     + g("qiberg")) * c.latice * c.inv_cp)

    s["qv"] = jnp.where(context, s["qv"] + qv_tend * dt, s["qv"])
    s["th"] = jnp.where(context,
                        s["th"] + jnp.asarray(t["exner"]) * th_sub * dt,
                        s["th"])
    return s, qi_wet_depos


def update_prognostic_liquid(qc2qr_accret_tend, nc_accret_tend,
                             qc2qr_autoconv_tend, nc2nr_autoconv_tend,
                             ncautr, nc_selfcollect_tend, qr2qv_evap_tend,
                             nr_evap_tend, nr_selfcollect_tend,
                             do_predict_nc: bool, nccnst,
                             do_prescribed_ccn: bool,
                             inv_rho, exner, latent_heat_vapor, dt,
                             th, qv, qc, nc, qr, nr, context):
    """Apply liquid-phase tendencies (update_prognostic_liquid,
    iparam=3). Returns (th, qv, qc, nc, qr, nr)."""
    th, qv, qc, nc, qr, nr = (jnp.asarray(a) for a in
                              (th, qv, qc, nc, qr, nr))
    qcacc = jnp.asarray(qc2qr_accret_tend)
    qcaut = jnp.asarray(qc2qr_autoconv_tend)
    qrevp = jnp.asarray(qr2qv_evap_tend)

    qc = jnp.where(context, qc + (-qcacc - qcaut) * dt, qc)
    qr = jnp.where(context, qr + (qcacc + qcaut - qrevp) * dt, qr)

    if do_predict_nc or do_prescribed_ccn:
        nc = jnp.where(context,
                       nc + (-jnp.asarray(nc_accret_tend)
                             - jnp.asarray(nc2nr_autoconv_tend)
                             + jnp.asarray(nc_selfcollect_tend)) * dt, nc)
    else:
        nc = jnp.where(context, nccnst * jnp.asarray(inv_rho), nc)

    # iparam == 3
    nr = jnp.where(context,
                   nr + (jnp.asarray(ncautr)
                         - jnp.asarray(nr_selfcollect_tend)
                         - jnp.asarray(nr_evap_tend)) * dt, nr)

    qv = jnp.where(context, qv + qrevp * dt, qv)
    th = jnp.where(context,
                   th + jnp.asarray(exner)
                   * (-qrevp * latent_heat_vapor * c.inv_cp) * dt, th)
    return th, qv, qc, nc, qr, nr
