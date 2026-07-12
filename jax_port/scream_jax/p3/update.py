"""Prognostic-state updates from process rates (batch 4b).

Source: components/eamxx/src/physics/p3/impl/p3_update_prognostics_impl.hpp
(update_prognostic_ice, update_prognostic_liquid). Pure-functional: state
arrays in, updated state out. iparam = 3 (as compiled in the C++), so the
rain-number source uses ncautr.
"""

import jax.numpy as jnp

from ..foundation import constants as c


def update_prognostic_ice(t, do_predict_nc: bool, log_wetgrowth, dt,
                          nmltratio, rho_qm_cloud, state,
                          use_hetfrz_classnuc: bool, context):
    """Apply ice-phase tendencies (Functions::update_prognostic_ice).

    t: dict of tendencies (C++ names); state: dict with th_atm, qv, qi, ni,
    qm, bm, qc, nc, qr, nr. Returns the updated state dict.
    """
    s = {k: jnp.asarray(v) for k, v in state.items()}
    g = lambda k: jnp.asarray(t[k])

    if use_hetfrz_classnuc:
        s["qc"] = jnp.where(context, s["qc"] + (
            -g("qcheti_cnt") - g("qicnt") - g("qc2qi_collect_tend")
            - g("qc2qr_ice_shed_tend") - g("qc2qi_berg_tend")) * dt, s["qc"])
    else:
        s["qc"] = jnp.where(context, s["qc"] + (
            -g("qc2qi_hetero_freeze_tend") - g("qc2qi_collect_tend")
            - g("qc2qr_ice_shed_tend") - g("qc2qi_berg_tend")) * dt, s["qc"])

    if do_predict_nc:
        if use_hetfrz_classnuc:
            s["nc"] = jnp.where(context, s["nc"] + (
                -g("nc_collect_tend") - g("ncheti_cnt") - g("nicnt")) * dt, s["nc"])
        else:
            s["nc"] = jnp.where(context, s["nc"] + (
                -g("nc_collect_tend") - g("nc2ni_immers_freeze_tend")) * dt, s["nc"])

    s["qr"] = jnp.where(context, s["qr"] + (
        -g("qr2qi_collect_tend") + g("qi2qr_melt_tend")
        - g("qr2qi_immers_freeze_tend") + g("qc2qr_ice_shed_tend")) * dt, s["qr"])
    s["nr"] = jnp.where(context, s["nr"] + (
        -g("nr_collect_tend") - g("nr2ni_immers_freeze_tend")
        + nmltratio * g("ni2nr_melt_tend") + g("nr_ice_shed_tend")
        + g("ncshdc")) * dt, s["nr"])

    # Sublimation/melting deplete qm/bm proportionally, then qi
    qi_ok = (s["qi"] >= c.QSMALL) & context
    qi_safe = jnp.where(qi_ok, s["qi"], 1.0)
    depl = (g("qi2qv_sublim_tend") + g("qi2qr_melt_tend")) / qi_safe * dt
    s["bm"] = jnp.where(qi_ok, s["bm"] - depl * s["bm"], s["bm"])
    s["qm"] = jnp.where(qi_ok, s["qm"] - depl * s["qm"], s["qm"])
    s["qi"] = jnp.where(qi_ok,
                        s["qi"] - (g("qi2qv_sublim_tend") + g("qi2qr_melt_tend")) * dt,
                        s["qi"])

    if use_hetfrz_classnuc:
        dum = (g("qr2qi_collect_tend") + g("qc2qi_collect_tend")
               + g("qr2qi_immers_freeze_tend") + g("qcheti_cnt") + g("qicnt")) * dt
        qi_src = (g("qv2qi_vapdep_tend") + g("qv2qi_nucleat_tend")
                  + g("qc2qi_berg_tend") + g("qinuc_cnt")) * dt
        bm_src = (g("qr2qi_collect_tend") * c.INV_RHO_RIMEMAX
                  + g("qc2qi_collect_tend") / jnp.asarray(rho_qm_cloud)
                  + (g("qr2qi_immers_freeze_tend") + g("qcheti_cnt")
                     + g("qicnt")) * c.INV_RHO_RIMEMAX) * dt
        ni_src = (g("ni_nucleat_tend") - g("ni2nr_melt_tend") - g("ni_sublim_tend")
                  - g("ni_selfcollect_tend") + g("nr2ni_immers_freeze_tend")
                  + g("ncheti_cnt") + g("nicnt") + g("ninuc_cnt")) * dt
    else:
        dum = (g("qr2qi_collect_tend") + g("qc2qi_collect_tend")
               + g("qr2qi_immers_freeze_tend") + g("qc2qi_hetero_freeze_tend")) * dt
        qi_src = (g("qv2qi_vapdep_tend") + g("qv2qi_nucleat_tend")
                  + g("qc2qi_berg_tend")) * dt
        bm_src = (g("qr2qi_collect_tend") * c.INV_RHO_RIMEMAX
                  + g("qc2qi_collect_tend") / jnp.asarray(rho_qm_cloud)
                  + (g("qr2qi_immers_freeze_tend")
                     + g("qc2qi_hetero_freeze_tend")) * c.INV_RHO_RIMEMAX) * dt
        ni_src = (g("ni_nucleat_tend") - g("ni2nr_melt_tend") - g("ni_sublim_tend")
                  - g("ni_selfcollect_tend") + g("nr2ni_immers_freeze_tend")
                  + g("nc2ni_immers_freeze_tend")) * dt

    s["qi"] = jnp.where(context, s["qi"] + qi_src + dum, s["qi"])
    s["qm"] = jnp.where(context, s["qm"] + dum, s["qm"])
    s["bm"] = jnp.where(context, s["bm"] + bm_src, s["bm"])
    s["ni"] = jnp.where(context, s["ni"] + ni_src, s["ni"])

    neg_qm = (s["qm"] < 0.0) & context
    s["qm"] = jnp.where(neg_qm, 0.0, s["qm"])
    s["bm"] = jnp.where(neg_qm, 0.0, s["bm"])

    wet = jnp.asarray(log_wetgrowth) & context
    s["qm"] = jnp.where(wet, s["qi"], s["qm"])
    s["bm"] = jnp.where(wet, s["qm"] * c.INV_RHO_RIMEMAX, s["bm"])

    if use_hetfrz_classnuc:
        qv_tend = (-g("qv2qi_vapdep_tend") + g("qi2qv_sublim_tend")
                   - g("qv2qi_nucleat_tend") - g("qinuc_cnt"))
        th_sub = ((g("qv2qi_vapdep_tend") - g("qi2qv_sublim_tend")
                   + g("qv2qi_nucleat_tend") + g("qinuc_cnt"))
                  * (c.LatVap + c.LatIce) * c.INV_CP
                  + (g("qr2qi_collect_tend") + g("qc2qi_collect_tend")
                     + g("qcheti_cnt") + g("qicnt")
                     + g("qr2qi_immers_freeze_tend") - g("qi2qr_melt_tend")
                     + g("qc2qi_berg_tend")) * c.LatIce * c.INV_CP)
    else:
        qv_tend = (-g("qv2qi_vapdep_tend") + g("qi2qv_sublim_tend")
                   - g("qv2qi_nucleat_tend"))
        th_sub = ((g("qv2qi_vapdep_tend") - g("qi2qv_sublim_tend")
                   + g("qv2qi_nucleat_tend")) * (c.LatVap + c.LatIce) * c.INV_CP
                  + (g("qr2qi_collect_tend") + g("qc2qi_collect_tend")
                     + g("qc2qi_hetero_freeze_tend")
                     + g("qr2qi_immers_freeze_tend") - g("qi2qr_melt_tend")
                     + g("qc2qi_berg_tend")) * c.LatIce * c.INV_CP)

    s["qv"] = jnp.where(context, s["qv"] + qv_tend * dt, s["qv"])
    s["th_atm"] = jnp.where(context,
                            s["th_atm"] + jnp.asarray(t["inv_exner"]) * th_sub * dt,
                            s["th_atm"])
    return s


def update_prognostic_liquid(qc2qr_accret_tend, nc_accret_tend,
                             qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr,
                             nc_selfcollect_tend, qr2qv_evap_tend,
                             nr_evap_tend, nr_selfcollect_tend,
                             do_predict_nc: bool, do_prescribed_ccn: bool,
                             inv_rho, inv_exner, dt,
                             th_atm, qv, qc, nc, qr, nr, context):
    """Apply liquid-phase tendencies (Functions::update_prognostic_liquid).
    Returns (th_atm, qv, qc, nc, qr, nr)."""
    th_atm, qv, qc, nc, qr, nr = (jnp.asarray(a) for a in
                                  (th_atm, qv, qc, nc, qr, nr))
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
        nc = jnp.where(context, c.NCCNST * jnp.asarray(inv_rho), nc)

    # IPARAM == 3
    nr = jnp.where(context,
                   nr + (jnp.asarray(ncautr) - jnp.asarray(nr_selfcollect_tend)
                         - jnp.asarray(nr_evap_tend)) * dt, nr)

    qv = jnp.where(context, qv + qrevp * dt, qv)
    th_atm = jnp.where(context,
                       th_atm + jnp.asarray(inv_exner)
                       * (-qrevp * c.LatVap * c.INV_CP) * dt, th_atm)
    return th_atm, qv, qc, nc, qr, nr
