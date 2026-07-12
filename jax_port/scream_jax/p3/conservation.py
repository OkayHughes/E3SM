"""P3 conservation limiters and in-cloud mixing ratios (batch 3).

Sources (components/eamxx/src/physics/p3/impl/):
  p3_conservation_impl.hpp (cloud/rain/ice water),
  p3_nc_conservation_impl.hpp, p3_nr_conservation_impl.hpp,
  p3_ni_conservation_impl.hpp, p3_ice_supersat_conservation_impl.hpp,
  p3_prevent_liq_supersaturation_impl.hpp, p3_impose_max_total_ni_impl.hpp,
  p3_incloud_mixingratios_impl.hpp

All limiters rescale groups of tendencies in place in the C++; here they
take and return the tendencies. use_hetfrz_classnuc is a static bool
switching which tendency set participates. NOTE (transcribed exactly):
cloud_water_conservation zeroes qv2qi_vapdep/qi2qv_sublim wherever cloud
water is present but unlimited (ratio == 1) — within liquid cloud, ice
growth proceeds via the Bergeron term instead.
"""

import jax.numpy as jnp

from ..foundation import constants as c
from ..foundation.saturation import SaturationFcn, qv_sat_dry


def cloud_water_conservation(qc, dt, qcaut, qcacc, qccol, qchetc, qcshd,
                             qcberg, qisub, qidep, qcheti_cnt, qicnt,
                             use_hetfrz_classnuc: bool, context,
                             cld_frac_l, cld_frac_i,
                             use_separate_ice_liq_frac: bool = False):
    """Functions::cloud_water_conservation. Returns the rescaled
    (qcaut, qcacc, qccol, qchetc, qcshd, qcberg, qisub, qidep,
    qcheti_cnt, qicnt)."""
    qc = jnp.asarray(qc)
    qcaut, qcacc, qccol, qchetc, qcshd, qcberg, qisub, qidep = (
        jnp.asarray(a) for a in (qcaut, qcacc, qccol, qchetc, qcshd,
                                 qcberg, qisub, qidep))
    qcheti_cnt = jnp.asarray(qcheti_cnt)
    qicnt = jnp.asarray(qicnt)

    if use_hetfrz_classnuc:
        sinks = (qcaut + qcacc + qccol + qcheti_cnt + qcshd + qcberg) * dt
    else:
        sinks = (qcaut + qcacc + qccol + qchetc + qcshd + qcberg) * dt
    sources = qc

    if use_separate_ice_liq_frac:
        il_cldm = jnp.minimum(jnp.asarray(cld_frac_i), jnp.asarray(cld_frac_l))
        cld_frac_glac = jnp.maximum(jnp.asarray(cld_frac_i) - il_cldm, 0.0001)
    else:
        il_cldm = 1.0
        cld_frac_glac = 1.0

    enforce = (sinks > sources) & (sinks >= c.QTENDSMALL) & context
    sinks_safe = jnp.where(enforce, sinks, 1.0)
    ratio = jnp.where(enforce, sources / sinks_safe, 1.0)

    qcaut = jnp.where(enforce, qcaut * ratio, qcaut)
    qcacc = jnp.where(enforce, qcacc * ratio, qcacc)
    qccol = jnp.where(enforce, qccol * ratio, qccol)
    if use_hetfrz_classnuc:
        qcheti_cnt = jnp.where(enforce, qcheti_cnt * ratio, qcheti_cnt)
        qicnt = jnp.where(enforce, qicnt * ratio, qicnt)
    else:
        qchetc = jnp.where(enforce, qchetc * ratio, qchetc)
    qcshd = jnp.where(enforce, qcshd * ratio, qcshd)
    qcberg = jnp.where(enforce, qcberg * ratio, qcberg)

    enforce2 = (sources > c.QTENDSMALL) & context
    if use_separate_ice_liq_frac:
        qidep = jnp.where(enforce2,
                          qidep + qidep * (1.0 - ratio) * (il_cldm / cld_frac_glac),
                          qidep)
        qisub = jnp.where(enforce2,
                          qisub + qisub * (1.0 - ratio) * (il_cldm / cld_frac_glac),
                          qisub)
    else:
        qidep = jnp.where(enforce2, qidep * (1.0 - ratio), qidep)
        qisub = jnp.where(enforce2, qisub * (1.0 - ratio), qisub)

    return (qcaut, qcacc, qccol, qchetc, qcshd, qcberg, qisub, qidep,
            qcheti_cnt, qicnt)


def rain_water_conservation(qr, qcaut, qcacc, qimlt, qcshd, dt,
                            qrevp, qrcol, qrheti, context):
    """Functions::rain_water_conservation. Returns (qrevp, qrcol, qrheti)."""
    qrevp, qrcol, qrheti = (jnp.asarray(a) for a in (qrevp, qrcol, qrheti))
    sinks = (qrevp + qrcol + qrheti) * dt
    sources = jnp.asarray(qr) + (jnp.asarray(qcaut) + jnp.asarray(qcacc)
                                 + jnp.asarray(qimlt) + jnp.asarray(qcshd)) * dt
    enforce = (sinks > sources) & (sinks >= c.QTENDSMALL) & context
    ratio = jnp.where(enforce, sources / jnp.where(enforce, sinks, 1.0), 1.0)
    return (jnp.where(enforce, qrevp * ratio, qrevp),
            jnp.where(enforce, qrcol * ratio, qrcol),
            jnp.where(enforce, qrheti * ratio, qrheti))


def ice_water_conservation(qi, qidep, qinuc, qcberg, qrcol, qccol, qrheti,
                           qchetc, dt, qinuc_cnt, qcheti_cnt, qicnt,
                           qisub, qimlt, use_hetfrz_classnuc: bool, context):
    """Functions::ice_water_conservation. Returns (qisub, qimlt)."""
    qisub, qimlt = jnp.asarray(qisub), jnp.asarray(qimlt)
    sinks = (qisub + qimlt) * dt
    common = (jnp.asarray(qidep) + jnp.asarray(qinuc) + jnp.asarray(qrcol)
              + jnp.asarray(qccol) + jnp.asarray(qrheti) + jnp.asarray(qcberg))
    if use_hetfrz_classnuc:
        sources = jnp.asarray(qi) + (common + jnp.asarray(qinuc_cnt)
                                     + jnp.asarray(qcheti_cnt) + jnp.asarray(qicnt)) * dt
    else:
        sources = jnp.asarray(qi) + (common + jnp.asarray(qchetc)) * dt
    enforce = (sinks > sources) & (sinks >= c.QTENDSMALL) & context
    ratio = jnp.where(enforce, sources / jnp.where(enforce, sinks, 1.0), 1.0)
    return (jnp.where(enforce, qisub * ratio, qisub),
            jnp.where(enforce, qimlt * ratio, qimlt))


def nc_conservation(nc, ncslf, dt, nccol, nchetc, ncacc, ncautc,
                    ncheti_cnt, nicnt, use_hetfrz_classnuc: bool, context):
    """Functions::nc_conservation. Returns
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
    if use_hetfrz_classnuc:
        ncheti_cnt = jnp.where(mask, ncheti_cnt * ratio, ncheti_cnt)
        nicnt = jnp.where(mask, nicnt * ratio, nicnt)
    else:
        nchetc = jnp.where(mask, nchetc * ratio, nchetc)
    return (nccol, nchetc, jnp.where(mask, ncacc * ratio, ncacc),
            jnp.where(mask, ncautc * ratio, ncautc), ncheti_cnt, nicnt)


def nr_conservation(nr, nimlt, nrshdr, ncshdc, ncautr, dt, nmltratio,
                    nrcol, nrheti, nrslf, nrevp, context):
    """Functions::nr_conservation. Returns (nrcol, nrheti, nrslf, nrevp)."""
    nrcol, nrheti, nrslf, nrevp = (jnp.asarray(a) for a in
                                   (nrcol, nrheti, nrslf, nrevp))
    sink = (nrcol + nrheti + nrslf + nrevp) * dt
    source = jnp.asarray(nr) + (jnp.asarray(nimlt) * nmltratio
                                + jnp.asarray(nrshdr) + jnp.asarray(ncshdc)
                                + jnp.asarray(ncautr)) * dt
    mask = (sink > source) & context
    ratio = jnp.where(mask, source / jnp.where(mask, sink, 1.0), 1.0)
    return tuple(jnp.where(mask, a * ratio, a)
                 for a in (nrcol, nrheti, nrslf, nrevp))


def ni_conservation(ni, ninuc, nrheti, ncheti, ncheti_cnt, nicnt, ninuc_cnt,
                    dt, nimlt, nisub, nislf, use_hetfrz_classnuc: bool, context):
    """Functions::ni_conservation. Returns (nimlt, nisub, nislf)."""
    nimlt, nisub, nislf = (jnp.asarray(a) for a in (nimlt, nisub, nislf))
    sink = (nimlt + nisub + nislf) * dt
    if use_hetfrz_classnuc:
        source = jnp.asarray(ni) + (jnp.asarray(ninuc) + jnp.asarray(nrheti)
                                    + jnp.asarray(ncheti_cnt) + jnp.asarray(nicnt)
                                    + jnp.asarray(ninuc_cnt)) * dt
    else:
        source = jnp.asarray(ni) + (jnp.asarray(ninuc) + jnp.asarray(nrheti)
                                    + jnp.asarray(ncheti)) * dt
    mask = (sink > source) & context
    ratio = jnp.where(mask, source / jnp.where(mask, sink, 1.0), 1.0)
    return tuple(jnp.where(mask, a * ratio, a) for a in (nimlt, nisub, nislf))


def ice_supersat_conservation(qidep, qinuc, qinuc_cnt, cld_frac_i, qv,
                              qv_sat_i, t_atm, dt, qi2qv_sublim_tend,
                              qr2qv_evap_tend, use_hetfrz_classnuc: bool,
                              context):
    """Functions::ice_supersat_conservation. Returns
    (qidep, qinuc, qinuc_cnt)."""
    qidep, qinuc = jnp.asarray(qidep), jnp.asarray(qinuc)
    qinuc_cnt = jnp.asarray(qinuc_cnt)
    qv_sat_i = jnp.asarray(qv_sat_i)
    t_atm = jnp.asarray(t_atm)

    latsublim2 = (c.LatVap + c.LatIce) ** 2
    qv_sink = qidep + qinuc + (qinuc_cnt if use_hetfrz_classnuc else 0.0)
    mask = (qv_sink > c.QSMALL) & (jnp.asarray(cld_frac_i) > 1e-20) & context

    qv_avail = (jnp.asarray(qv)
                + (jnp.asarray(qi2qv_sublim_tend) + jnp.asarray(qr2qv_evap_tend)) * dt
                - qv_sat_i) \
        / (1.0 + latsublim2 * qv_sat_i / (c.CP * c.RH2O * t_atm ** 2)) / dt
    qv_avail = jnp.maximum(qv_avail, 0.0)

    limited = (qv_sink > qv_avail) & mask
    sink_safe = jnp.where(limited, qv_sink, 1.0)
    fract = jnp.where(limited, qv_avail / sink_safe, 1.0)
    if use_hetfrz_classnuc:
        qinuc_cnt = jnp.where(limited, qinuc_cnt * fract, qinuc_cnt)
    return (jnp.where(limited, qidep * fract, qidep),
            jnp.where(limited, qinuc * fract, qinuc), qinuc_cnt)


def prevent_liq_supersaturation(pres, t_atm, qv, dt, qv2qi_vapdep_tend,
                                qinuc, qi2qv_sublim_tend, qr2qv_evap_tend,
                                context):
    """Functions::prevent_liq_supersaturation. Returns
    (qi2qv_sublim_tend, qr2qv_evap_tend)."""
    qisub = jnp.asarray(qi2qv_sublim_tend)
    qrevp = jnp.asarray(qr2qv_evap_tend)
    qv = jnp.asarray(qv)
    t_atm = jnp.asarray(t_atm)

    qv_sources = qisub + qrevp
    has = (qv_sources >= c.QSMALL) & context

    qv_sinks = jnp.asarray(qv2qi_vapdep_tend) + jnp.asarray(qinuc)
    T_end = t_atm + ((qv_sinks - qisub) * (c.LatVap + c.LatIce) * c.INV_CP
                     - qrevp * c.LatVap * c.INV_CP) * dt
    T_end_safe = jnp.where(has, T_end, 273.0)
    qsl = qv_sat_dry(T_end_safe, jnp.asarray(pres), False,
                     SaturationFcn.MURPHY_KOOP)
    A = (c.LatVap * qsl * dt * c.INV_CP / (c.RV * T_end_safe ** 2)
         * ((c.LatVap + c.LatIce) * qisub + c.LatVap * qrevp))
    src_safe = jnp.where(has, qv_sources * dt + A, 1.0)
    frac = jnp.clip((qsl - qv + qv_sinks * dt + A) / src_safe, 0.0, 1.0)
    return (jnp.where(has, frac * qisub, qisub),
            jnp.where(has, frac * qrevp, qrevp))


def impose_max_total_ni(ni_local, max_total_ni, inv_rho_local, context):
    """Functions::impose_max_total_ni. Returns clipped ni_local."""
    ni_local = jnp.asarray(ni_local)
    mask = (ni_local >= 1e-20) & context
    ni_safe = jnp.where(mask, ni_local, 1.0)
    return jnp.where(mask,
                     ni_local * jnp.minimum(
                         max_total_ni * jnp.asarray(inv_rho_local) / ni_safe, 1.0),
                     ni_local)


def calculate_incloud_mixingratios(qc, qr, qi, qm, nc, nr, ni, bm,
                                   inv_cld_frac_l, inv_cld_frac_i,
                                   inv_cld_frac_r, context):
    """Functions::calculate_incloud_mixingratios. Returns
    (qc_incld, qr_incld, qi_incld, qm_incld,
     nc_incld, nr_incld, ni_incld, bm_incld)."""
    qc, qr, qi, qm = (jnp.asarray(a) for a in (qc, qr, qi, qm))
    nc, nr, ni, bm = (jnp.asarray(a) for a in (nc, nr, ni, bm))
    il, ii, ir = (jnp.asarray(a) for a in
                  (inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r))

    qc_ok = (qc >= c.QSMALL) & context
    qi_ok = (qi >= c.QSMALL) & context
    qim_ok = (qi >= c.QSMALL) & (qm >= c.QSMALL) & context
    qr_ok = (qr >= c.QSMALL) & context

    qc_incld = jnp.where(qc_ok, qc * il, 0.0)
    nc_incld = jnp.where(qc_ok, jnp.maximum(nc * il, 0.0), 0.0)
    qi_incld = jnp.where(qi_ok, qi * ii, 0.0)
    ni_incld = jnp.where(qi_ok, jnp.maximum(ni * ii, 0.0), 0.0)
    qm_incld = jnp.where(qim_ok, qm * ii, 0.0)
    # NB: bm uses the LIQUID inverse cloud fraction in the C++.
    bm_incld = jnp.where(qim_ok, jnp.maximum(bm * il, 0.0), 0.0)
    qr_incld = jnp.where(qr_ok, qr * ir, 0.0)
    nr_incld = jnp.where(qr_ok, jnp.maximum(nr * ir, 0.0), 0.0)

    over = ((qc_incld > c.incloud_limit) | (qi_incld > c.incloud_limit)
            | (qr_incld > c.precip_limit) | (bm_incld > c.incloud_limit)) & context
    qc_incld = jnp.where(over, jnp.minimum(qc_incld, c.incloud_limit), qc_incld)
    qi_incld = jnp.where(over, jnp.minimum(qi_incld, c.incloud_limit), qi_incld)
    bm_incld = jnp.where(over, jnp.minimum(bm_incld, c.incloud_limit), bm_incld)
    qr_incld = jnp.where(over, jnp.minimum(qr_incld, c.precip_limit), qr_incld)

    return (qc_incld, qr_incld, qi_incld, qm_incld,
            nc_incld, nr_incld, ni_incld, bm_incld)
