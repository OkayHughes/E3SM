"""P3 conservation limiters and in-cloud mixing ratios.

Source: micro_p3.F90 / micro_p3_utils.F90 (eam variant):
cloud/rain/ice_water_conservation, nc/nr/ni_conservation,
prevent_ice_overdepletion, ice_supersat_conservation,
impose_max_total_ni, calculate_incloud_mixingratios. COPIED from
scream_jax/p3/conservation.py; adaptations:

  * prevent_ice_overdepletion is the EAM-only limiter (scream's
    prevent_liq_supersaturation replaced it upstream): if end-of-step
    vapor would exceed LIQUID saturation, rain evaporation and ice
    sublimation are each rescaled against the sink-side saturation
    deficit (with liquid resp. sublimation latent heats).
  * threshold for the water conservation limiters is 1e-20 (as in the
    Fortran; scream_jax's QTENDSMALL has the same value).
  * ice_water_conservation hetfrz branch omits qc2qi_hetero (EAM).
  * cloud_water_conservation: unchanged semantics from the copy
    (qidep/qisub *(1-ratio) whenever qc > 1e-20).
"""

import jax.numpy as jnp

from . import constants as c
from .saturation import qv_sat

QTENDSMALL = 1.0e-20


def cloud_water_conservation(qc, dt, qcaut, qcacc, qccol, qchetc, qcshd,
                             qcberg, qisub, qidep, qcheti_cnt, qicnt,
                             use_hetfrz_classnuc: bool, context):
    """cloud_water_conservation. Returns the rescaled
    (qcaut, qcacc, qccol, qchetc, qcshd, qcberg, qisub, qidep,
    qcheti_cnt, qicnt)."""
    qc = jnp.asarray(qc)
    qcaut, qcacc, qccol, qchetc, qcshd, qcberg, qisub, qidep = (
        jnp.asarray(a) for a in (qcaut, qcacc, qccol, qchetc, qcshd,
                                 qcberg, qisub, qidep))
    qcheti_cnt = jnp.asarray(qcheti_cnt)
    qicnt = jnp.asarray(qicnt)

    if use_hetfrz_classnuc:
        sinks = (qcaut + qcacc + qccol + qcheti_cnt + qcshd + qcberg
                 + qicnt) * dt
    else:
        sinks = (qcaut + qcacc + qccol + qchetc + qcshd + qcberg) * dt

    enforce = (sinks > qc) & (sinks >= QTENDSMALL) & context
    sinks_safe = jnp.where(enforce, sinks, 1.0)
    ratio = jnp.where(enforce, qc / sinks_safe, 1.0)

    qcaut = jnp.where(enforce, qcaut * ratio, qcaut)
    qcacc = jnp.where(enforce, qcacc * ratio, qcacc)
    qccol = jnp.where(enforce, qccol * ratio, qccol)
    qcshd = jnp.where(enforce, qcshd * ratio, qcshd)
    qcberg = jnp.where(enforce, qcberg * ratio, qcberg)
    if use_hetfrz_classnuc:
        qcheti_cnt = jnp.where(enforce, qcheti_cnt * ratio, qcheti_cnt)
        qicnt = jnp.where(enforce, qicnt * ratio, qicnt)
    else:
        qchetc = jnp.where(enforce, qchetc * ratio, qchetc)

    # ratio is also the fraction of the step with liquid: Bergeron for
    # "ratio" of the step, deposition/sublimation for the rest
    enforce2 = (qc > 1.0e-20) & context
    qidep = jnp.where(enforce2, qidep * (1.0 - ratio), qidep)
    qisub = jnp.where(enforce2, qisub * (1.0 - ratio), qisub)

    return (qcaut, qcacc, qccol, qchetc, qcshd, qcberg, qisub, qidep,
            qcheti_cnt, qicnt)


def rain_water_conservation(qr, qcaut, qcacc, qimlt, qcshd, dt,
                            qrevp, qrcol, qrheti, context):
    """rain_water_conservation. Returns (qrevp, qrcol, qrheti)."""
    qrevp, qrcol, qrheti = (jnp.asarray(a) for a in (qrevp, qrcol, qrheti))
    sinks = (qrevp + qrcol + qrheti) * dt
    sources = jnp.asarray(qr) + (jnp.asarray(qcaut) + jnp.asarray(qcacc)
                                 + jnp.asarray(qimlt)
                                 + jnp.asarray(qcshd)) * dt
    enforce = (sinks > sources) & (sinks >= QTENDSMALL) & context
    ratio = jnp.where(enforce, sources / jnp.where(enforce, sinks, 1.0), 1.0)
    return (jnp.where(enforce, qrevp * ratio, qrevp),
            jnp.where(enforce, qrcol * ratio, qrcol),
            jnp.where(enforce, qrheti * ratio, qrheti))


def ice_water_conservation(qi, qidep, qinuc, qiberg, qrcol, qccol, qrheti,
                           qchetc, dt, qinuc_cnt, qcheti_cnt, qicnt,
                           qisub, qimlt, use_hetfrz_classnuc: bool, context):
    """ice_water_conservation. Returns (qisub, qimlt). NB the EAM
    hetfrz branch omits qc2qi_hetero_freeze_tend from the sources."""
    qisub, qimlt = jnp.asarray(qisub), jnp.asarray(qimlt)
    sinks = (qisub + qimlt) * dt
    if use_hetfrz_classnuc:
        sources = jnp.asarray(qi) + (
            jnp.asarray(qidep) + jnp.asarray(qinuc) + jnp.asarray(qrcol)
            + jnp.asarray(qccol) + jnp.asarray(qrheti) + jnp.asarray(qiberg)
            + jnp.asarray(qinuc_cnt) + jnp.asarray(qcheti_cnt)
            + jnp.asarray(qicnt)) * dt
    else:
        sources = jnp.asarray(qi) + (
            jnp.asarray(qidep) + jnp.asarray(qinuc) + jnp.asarray(qrcol)
            + jnp.asarray(qccol) + jnp.asarray(qrheti) + jnp.asarray(qchetc)
            + jnp.asarray(qiberg)) * dt
    enforce = (sinks > sources) & (sinks >= QTENDSMALL) & context
    ratio = jnp.where(enforce, sources / jnp.where(enforce, sinks, 1.0), 1.0)
    return (jnp.where(enforce, qisub * ratio, qisub),
            jnp.where(enforce, qimlt * ratio, qimlt))


def nc_conservation(nc, ncslf, dt, nccol, nchetc, ncacc, ncautc,
                    ncheti_cnt, nicnt, use_hetfrz_classnuc: bool, context):
    """nc_conservation. Returns
    (nccol, nchetc, ncacc, ncautc, ncheti_cnt, nicnt)."""
    nccol, nchetc, ncacc, ncautc = (jnp.asarray(a) for a in
                                    (nccol, nchetc, ncacc, ncautc))
    ncheti_cnt, nicnt = jnp.asarray(ncheti_cnt), jnp.asarray(nicnt)
    if use_hetfrz_classnuc:
        sink = (nccol + ncheti_cnt + ncacc + ncautc + nicnt) * dt
    else:
        sink = (nccol + nchetc + ncacc + ncautc) * dt
    source = jnp.asarray(nc) + jnp.asarray(ncslf) * dt
    mask = (sink > source) & context
    ratio = jnp.where(mask, source / jnp.where(mask, sink, 1.0), 1.0)
    nccol = jnp.where(mask, nccol * ratio, nccol)
    ncacc = jnp.where(mask, ncacc * ratio, ncacc)
    ncautc = jnp.where(mask, ncautc * ratio, ncautc)
    if use_hetfrz_classnuc:
        ncheti_cnt = jnp.where(mask, ncheti_cnt * ratio, ncheti_cnt)
        nicnt = jnp.where(mask, nicnt * ratio, nicnt)
    else:
        nchetc = jnp.where(mask, nchetc * ratio, nchetc)
    return nccol, nchetc, ncacc, ncautc, ncheti_cnt, nicnt


def nr_conservation(nr, nimlt, nrshdr, ncshdc, ncautc, dt, nmltratio,
                    nrcol, nrheti, nrslf, nrevp, context):
    """nr_conservation. Returns (nrcol, nrheti, nrslf, nrevp)."""
    nrcol, nrheti, nrslf, nrevp = (jnp.asarray(a) for a in
                                   (nrcol, nrheti, nrslf, nrevp))
    sink = (nrcol + nrheti + nrslf + nrevp) * dt
    source = jnp.asarray(nr) + (jnp.asarray(nimlt) * nmltratio
                                + jnp.asarray(nrshdr) + jnp.asarray(ncshdc)
                                + jnp.asarray(ncautc)) * dt
    mask = (sink > source) & context
    ratio = jnp.where(mask, source / jnp.where(mask, sink, 1.0), 1.0)
    return tuple(jnp.where(mask, a * ratio, a)
                 for a in (nrcol, nrheti, nrslf, nrevp))


def ni_conservation(ni, ninuc, nrheti, ncheti, ncheti_cnt, nicnt, ninuc_cnt,
                    dt, nimlt, nisub, nislf, use_hetfrz_classnuc: bool,
                    context):
    """ni_conservation. Returns (nimlt, nisub, nislf)."""
    nimlt, nisub, nislf = (jnp.asarray(a) for a in (nimlt, nisub, nislf))
    sink = (nimlt + nisub + nislf) * dt
    if use_hetfrz_classnuc:
        source = jnp.asarray(ni) + (jnp.asarray(ninuc) + jnp.asarray(nrheti)
                                    + jnp.asarray(ncheti_cnt)
                                    + jnp.asarray(nicnt)
                                    + jnp.asarray(ninuc_cnt)) * dt
    else:
        source = jnp.asarray(ni) + (jnp.asarray(ninuc) + jnp.asarray(nrheti)
                                    + jnp.asarray(ncheti)) * dt
    mask = (sink > source) & context
    ratio = jnp.where(mask, source / jnp.where(mask, sink, 1.0), 1.0)
    return tuple(jnp.where(mask, a * ratio, a) for a in (nimlt, nisub, nislf))


def prevent_ice_overdepletion(pres, t_atm, qv, latent_heat_vapor,
                              latent_heat_sublim, inv_dt, dt, qidep, qinuc,
                              qinuc_cnt, qisub, qrevp, context):
    """prevent_ice_overdepletion (EAM-only). Returns (qisub, qrevp)."""
    pres = jnp.asarray(pres)
    t_atm = jnp.asarray(t_atm)
    qv = jnp.asarray(qv)
    qidep, qinuc = jnp.asarray(qidep), jnp.asarray(qinuc)
    qinuc_cnt = jnp.asarray(qinuc_cnt)
    qisub, qrevp = jnp.asarray(qisub), jnp.asarray(qrevp)

    qtmp_all = qv - (qidep + qinuc + qinuc_cnt) * dt + (qisub + qrevp) * dt
    ttmp_all = t_atm + ((qidep - qisub + qinuc + qinuc_cnt)
                        * latent_heat_sublim * c.inv_cp
                        + (-qrevp * latent_heat_vapor * c.inv_cp)) * dt
    ttmp_safe = jnp.where(ttmp_all > 0, ttmp_all, 273.0)
    qv_sat_l = qv_sat(ttmp_safe, pres, False)

    limit = (qtmp_all > qv_sat_l) & context

    q_sink = qv - (qidep + qinuc + qinuc_cnt) * dt
    t_sink = t_atm + ((qidep + qinuc + qinuc_cnt) * latent_heat_sublim
                      * c.inv_cp) * dt
    t_sink_safe = jnp.where(t_sink > 0, t_sink, 273.0)
    dumqv_sat_l = qv_sat(t_sink_safe, pres, False)
    qv_source_evp = qrevp + qisub
    src_safe = jnp.maximum(qv_source_evp, 1.0e-20)

    qrevp_satadj = (q_sink - dumqv_sat_l) \
        / (1.0 + (latent_heat_vapor * latent_heat_vapor) * dumqv_sat_l
           / (c.cp * c.rv * (t_sink_safe * t_sink_safe))) * inv_dt
    qrevp_new = qrevp * jnp.minimum(
        1.0, jnp.maximum(0.0, -qrevp_satadj) / src_safe)

    dumqv_sat_i = qv_sat(t_sink_safe, pres, True)
    qisub_satadj = (q_sink - dumqv_sat_i) \
        / (1.0 + (latent_heat_sublim * latent_heat_sublim) * dumqv_sat_i
           / (c.cp * c.rv * (t_sink_safe * t_sink_safe))) * inv_dt
    qisub_new = qisub * jnp.minimum(
        1.0, jnp.maximum(0.0, -qisub_satadj) / src_safe)

    return (jnp.where(limit, qisub_new, qisub),
            jnp.where(limit, qrevp_new, qrevp))


def ice_supersat_conservation(qidep, qinuc, qisub, qrevp, qinuc_cnt,
                              cld_frac_i, qv, qv_sat_i, latent_heat_sublim,
                              t_atm, dt, use_hetfrz_classnuc: bool, context):
    """ice_supersat_conservation. Returns (qidep, qinuc, qinuc_cnt)."""
    qidep, qinuc = jnp.asarray(qidep), jnp.asarray(qinuc)
    qinuc_cnt = jnp.asarray(qinuc_cnt)
    qv_sat_i = jnp.asarray(qv_sat_i)
    t_atm = jnp.asarray(t_atm)

    if use_hetfrz_classnuc:
        qv_sink = qidep + qinuc + qinuc_cnt
    else:
        qv_sink = qidep + qinuc
    mask = (qv_sink > c.qsmall) & (jnp.asarray(cld_frac_i) > 1.0e-20) & context

    qv_avail = (jnp.asarray(qv) + (jnp.asarray(qisub) + jnp.asarray(qrevp))
                * dt - qv_sat_i) \
        / (1.0 + (latent_heat_sublim * latent_heat_sublim) * qv_sat_i
           / (c.cp * c.rv * (t_atm * t_atm))) / dt
    qv_avail = jnp.maximum(qv_avail, 0.0)

    limited = (qv_sink > qv_avail) & mask
    sink_safe = jnp.where(limited, qv_sink, 1.0)
    fract = jnp.where(limited, qv_avail / sink_safe, 1.0)
    if use_hetfrz_classnuc:
        qinuc_cnt = jnp.where(limited, qinuc_cnt * fract, qinuc_cnt)
    return (jnp.where(limited, qidep * fract, qidep),
            jnp.where(limited, qinuc * fract, qinuc), qinuc_cnt)


def impose_max_total_ni(ni_local, max_total_ni, inv_rho_local, context):
    """impose_max_total_ni. Returns clipped ni_local."""
    ni_local = jnp.asarray(ni_local)
    mask = (ni_local >= 1.0e-20) & context
    ni_safe = jnp.where(mask, ni_local, 1.0)
    return jnp.where(mask,
                     ni_local * jnp.minimum(
                         max_total_ni * jnp.asarray(inv_rho_local) / ni_safe,
                         1.0),
                     ni_local)


def calculate_incloud_mixingratios(qc, qr, qi, qm, nc, nr, ni, bm,
                                   inv_cld_frac_l, inv_cld_frac_i,
                                   inv_cld_frac_r, context):
    """calculate_incloud_mixingratios (micro_p3_utils.F90). Returns
    (qc_incld, qr_incld, qi_incld, qm_incld,
     nc_incld, nr_incld, ni_incld, bm_incld)."""
    qc, qr, qi, qm = (jnp.asarray(a) for a in (qc, qr, qi, qm))
    nc, nr, ni, bm = (jnp.asarray(a) for a in (nc, nr, ni, bm))
    il, ii, ir = (jnp.asarray(a) for a in
                  (inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r))

    qc_ok = (qc >= c.qsmall) & context
    qi_ok = (qi >= c.qsmall) & context
    qim_ok = (qm >= c.qsmall) & (qi >= c.qsmall) & context
    qr_ok = (qr >= c.qsmall) & context

    qc_incld = jnp.where(qc_ok, qc * il, 0.0)
    nc_incld = jnp.where(qc_ok, jnp.maximum(nc * il, 0.0), 0.0)
    qi_incld = jnp.where(qi_ok, qi * ii, 0.0)
    ni_incld = jnp.where(qi_ok, jnp.maximum(ni * ii, 0.0), 0.0)
    qm_incld = jnp.where(qim_ok, qm * ii, 0.0)
    # NB: bm uses the LIQUID inverse cloud fraction in the Fortran too
    bm_incld = jnp.where(qim_ok, jnp.maximum(bm * il, 0.0), 0.0)
    qr_incld = jnp.where(qr_ok, qr * ir, 0.0)
    nr_incld = jnp.where(qr_ok, jnp.maximum(nr * ir, 0.0), 0.0)

    over = ((qc_incld > c.incloud_limit) | (qi_incld > c.incloud_limit)
            | (qr_incld > c.precip_limit) | (bm_incld > c.incloud_limit)) \
        & context
    qc_incld = jnp.where(over, jnp.minimum(qc_incld, c.incloud_limit),
                         qc_incld)
    qi_incld = jnp.where(over, jnp.minimum(qi_incld, c.incloud_limit),
                         qi_incld)
    bm_incld = jnp.where(over, jnp.minimum(bm_incld, c.incloud_limit),
                         bm_incld)
    qr_incld = jnp.where(over, jnp.minimum(qr_incld, c.precip_limit),
                         qr_incld)

    return (qc_incld, qr_incld, qi_incld, qm_incld,
            nc_incld, nr_incld, ni_incld, bm_incld)
