"""P3 sedimentation, homogeneous freezing and complete ice melting.

Source: micro_p3.F90 cloud/rain/ice_sedimentation,
generalized_sedimentation, calc_first_order_upwind_step,
compute_rain_fall_velocity, homogeneous_freezing, ice_complete_melting
(eam variant). COPIED from scream_jax/p3/sedimentation.py; adaptations:

  * the time-integrated fluxes: cflx (cloud), rflx + precip_liq_flux
    (rain) and sflx (ice) accumulate flux(k)*dt_sub over the band
    [k_qxtop, k_qxbot] each substep and are scaled by inv_dt at the
    end. precip_ice_flux is NOT accumulated anywhere (the Fortran
    leaves it zero). Cloud sedimentation ADDS to precip_liq_surf.
  * compute_rain_fall_velocity uses the EAM get_rain_dsd2
    (p3_max_mean_rain_size).
  * ice sedimentation has no ice_sedimentation_factor (EAMxx-only).
  * homogeneous_freezing tests th/exner (theta updated SEQUENTIALLY:
    the qr test sees the warming from the qc freezing at the same
    level), not the frozen-in T_atm.
  * NEW ice_complete_melting: instantaneous (fractional) melting of
    ice where th/exner > t_snow_melt (f32(273.15)+2), routed to cloud
    droplets (small unrimed crystals) or rain.

Level convention: k=0 model top, kbot=nlev-1 (kdir=-1); the C++-style
while(dt_left>tol) CFL substepping is a lax.while_loop over batched
columns; already-finished columns are masked out of all updates.
"""

import functools

import jax
import jax.numpy as jnp
from jax import lax
from jax.scipy.special import gammaln

from . import constants as c
from .dsd import calc_bulk_rho_rime, get_cloud_dsd2, get_rain_dsd2
from .table_lookups import (
    apply_table3,
    apply_table_ice,
    lookup_ice,
    lookup_table3,
)


def _tgamma(x):
    return jnp.exp(gammaln(x))


def _cbrt(x):
    """bfb_cbrt: x**(1/3) via pow (matches gfortran; NaN for x<0)."""
    return x ** (1.0 / 3.0)


def find_top_bottom(q, small):
    """find top/bottom loops for kdir = -1. Returns
    (k_qxtop, k_qxbot, present)."""
    q = jnp.asarray(q)
    nlev = q.shape[-1]
    has = q >= small
    present = jnp.any(has, axis=-1)
    k_top = jnp.argmax(has, axis=-1)
    k_bot = nlev - 1 - jnp.argmax(has[..., ::-1], axis=-1)
    k_top = jnp.where(present, k_top, 0)
    k_bot = jnp.where(present, k_bot, 0)
    return k_top, k_bot, present


def _upwind_substep(fields, velocities, rho, inv_rho, inv_dz, dt_sub, band):
    """One first-order upwind substep over the band mask
    (calc_first_order_upwind_step, kdir = -1). Returns (updated fields,
    fluxes)."""
    dt_sub = dt_sub[..., None]
    out, fluxes = [], []
    for r, V in zip(fields, velocities):
        flux = jnp.where(band, V * r * rho, 0.0)
        flux_above = jnp.concatenate(
            [jnp.zeros_like(flux[..., :1]), flux[..., :-1]], axis=-1)
        r_new = jnp.where(
            band, r + (flux_above - flux) * inv_dz * dt_sub * inv_rho, r)
        out.append(r_new)
        fluxes.append(flux)
    return tuple(out), tuple(fluxes)


def _generalized_substep(fields, velocities, rho, inv_rho, inv_dz, idx,
                         k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum,
                         active):
    """One iteration of generalized_sedimentation (kdir = -1). Returns
    (fields, fluxes, k_qxbot, dt_left, prt_accum, dt_sub)."""
    nstep = jnp.floor(Co_max + 1.0)
    dt_sub = jnp.where(active, dt_left / nstep, 0.0)
    at_ground = k_qxbot == kbot
    k_temp = jnp.where(at_ground, k_qxbot, k_qxbot + 1)
    band = ((idx >= k_qxtop[..., None]) & (idx <= k_temp[..., None])
            & active[..., None])
    fields, fluxes = _upwind_substep(fields, velocities, rho, inv_rho,
                                     inv_dz, dt_sub, band)
    prt_accum = prt_accum + jnp.where(active & at_ground,
                                      fluxes[0][..., -1] * dt_sub, 0.0)
    k_qxbot = jnp.where(active & ~at_ground, k_qxbot + 1, k_qxbot)
    dt_left = jnp.where(active, dt_left - dt_sub, dt_left)
    return fields, fluxes, k_qxbot, dt_left, prt_accum, dt_sub


def _accum_interface_flux(acc, flux, dt_sub, k_qxtop, k_qxbot, active):
    """cflx/rflx/sflx/precip_liq_flux accumulation:
    acc(k+1) += flux(k)*dt_sub over k in [k_qxtop, k_qxbot] (0-based,
    band AFTER the k_qxbot extension). acc lives on nlev+1 interfaces."""
    nlev1 = acc.shape[-1]
    idx_f = jnp.arange(nlev1)
    fmask = ((idx_f >= (k_qxtop + 1)[..., None])
             & (idx_f <= (k_qxbot + 1)[..., None]) & active[..., None])
    flux_shift = jnp.concatenate([flux[..., :1], flux], axis=-1)
    return jnp.where(fmask, acc + flux_shift * dt_sub[..., None], acc)


@functools.partial(jax.jit, static_argnames=("do_predict_nc",))
def cloud_sedimentation(qc_incld, rho, inv_rho, cld_frac_l, acn, inv_dz,
                        dt, inv_dt, do_predict_nc: bool,
                        qc, nc, nc_incld, mu_c, lamc,
                        precip_liq_surf_in, cflx_in,
                        qc_tend_in, nc_tend_in):
    """cloud_sedimentation. Returns a dict with qc, nc, qc_incld,
    nc_incld, mu_c, lamc, cflx, qc_tend, nc_tend, precip_liq_surf
    (precip_liq_surf ACCUMULATES: in + contribution)."""
    qc, nc = jnp.asarray(qc), jnp.asarray(nc)
    qc_incld, nc_incld = jnp.asarray(qc_incld), jnp.asarray(nc_incld)
    mu_c, lamc = jnp.asarray(mu_c), jnp.asarray(lamc)
    rho, inv_rho = jnp.asarray(rho), jnp.asarray(inv_rho)
    cld_frac_l = jnp.asarray(cld_frac_l)
    acn, inv_dz = jnp.asarray(acn), jnp.asarray(inv_dz)
    cflx = jnp.asarray(cflx_in)

    nlev = qc.shape[-1]
    kbot = nlev - 1
    idx = jnp.arange(nlev)

    k_qxtop, k_qxbot0, present = find_top_bottom(qc, c.qsmall)
    dt_left0 = jnp.where(present, dt, 0.0)
    prt0 = jnp.zeros_like(dt_left0)

    def cond(carry):
        return jnp.any(carry[6] > c.dt_left_tol)

    def body(carry):
        (qc, nc, qc_incld, nc_incld, mu_c, lamc, dt_left, prt_accum,
         k_qxbot, cflx) = carry
        active = dt_left > c.dt_left_tol
        co_band = ((idx >= k_qxtop[..., None]) & (idx <= k_qxbot[..., None])
                   & active[..., None])
        m = co_band & (qc_incld > c.qsmall)

        nc2, mu2, _, lam2, _, _ = get_cloud_dsd2(qc_incld, nc_incld, rho, m)
        nc_incld = jnp.where(m, nc2, nc_incld)
        mu_c = jnp.where(m, mu2, mu_c)
        lamc = jnp.where(m, lam2, lamc)
        nc = jnp.where(m, nc_incld * cld_frac_l, nc)

        lam_safe = jnp.where(m, lamc, 1.0)
        dum = 1.0 / lam_safe ** c.bcn
        V_qc = jnp.where(m, acn * _tgamma(4.0 + c.bcn + mu_c) * dum
                         / (_tgamma(mu_c + 4.0)), 0.0)
        if do_predict_nc:
            V_nc = jnp.where(m, acn * _tgamma(1.0 + c.bcn + mu_c) * dum
                             / (_tgamma(mu_c + 1.0)), 0.0)
            fields, velocities = (qc, nc), (V_qc, V_nc)
        else:
            fields, velocities = (qc,), (V_qc,)

        Co_max = jnp.max(
            jnp.where(co_band, V_qc * dt_left[..., None] * inv_dz, 0.0),
            axis=-1)

        fields, fluxes, k_qxbot, dt_left, prt_accum, dt_sub = \
            _generalized_substep(
                fields, velocities, rho, inv_rho, inv_dz, idx,
                k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, active)
        qc = fields[0]
        if do_predict_nc:
            nc = fields[1]
        cflx = _accum_interface_flux(cflx, fluxes[0], dt_sub,
                                     k_qxtop, k_qxbot, active)

        act = active[..., None]
        qc_incld = jnp.where(act, qc / cld_frac_l, qc_incld)
        nc_incld = jnp.where(act, nc / cld_frac_l, nc_incld)
        return (qc, nc, qc_incld, nc_incld, mu_c, lamc, dt_left, prt_accum,
                k_qxbot, cflx)

    carry = (qc, nc, qc_incld, nc_incld, mu_c, lamc,
             dt_left0, prt0, k_qxbot0, cflx)
    (qc, nc, qc_incld, nc_incld, mu_c, lamc,
     _, prt_accum, _, cflx) = lax.while_loop(cond, body, carry)

    cflx = jnp.where(present[..., None], cflx * inv_dt, cflx)
    precip_liq_surf = (jnp.asarray(precip_liq_surf_in)
                       + jnp.where(present,
                                   prt_accum * c.inv_rho_h2o * inv_dt, 0.0))
    return {
        "qc": qc, "nc": nc, "qc_incld": qc_incld, "nc_incld": nc_incld,
        "mu_c": mu_c, "lamc": lamc, "cflx": cflx,
        "qc_tend": (qc - jnp.asarray(qc_tend_in)) * inv_dt,
        "nc_tend": (nc - jnp.asarray(nc_tend_in)) * inv_dt,
        "precip_liq_surf": precip_liq_surf,
    }


def compute_rain_fall_velocity(vn_table_vals, vm_table_vals, qr_incld,
                               rhofacr, p3_max_mean_rain_size,
                               nr_incld, mu_r, lamr, context):
    """compute_rain_fall_velocity. Returns
    (nr_incld, mu_r, lamr, V_qr, V_nr)."""
    nr2, mu2, lam2, _, _ = get_rain_dsd2(qr_incld, nr_incld,
                                         p3_max_mean_rain_size, context)
    nr_incld = jnp.where(context, nr2, nr_incld)
    mu_r = jnp.where(context, mu2, mu_r)
    lamr = jnp.where(context, lam2, lamr)
    tab = lookup_table3(mu_r, lamr, context)
    V_qr = jnp.where(context, apply_table3(vm_table_vals, tab) * rhofacr, 0.0)
    V_nr = jnp.where(context, apply_table3(vn_table_vals, tab) * rhofacr, 0.0)
    return nr_incld, mu_r, lamr, V_qr, V_nr


@jax.jit
def rain_sedimentation(rho, inv_rho, rhofacr, cld_frac_r, inv_dz, qr_incld,
                       vn_table_vals, vm_table_vals, dt, inv_dt,
                       p3_max_mean_rain_size,
                       qr, nr, nr_incld, mu_r, lamr,
                       precip_liq_surf_in, precip_liq_flux_in, rflx_in,
                       qr_tend_in, nr_tend_in):
    """rain_sedimentation. Returns a dict with qr, nr, qr_incld,
    nr_incld, mu_r, lamr, precip_liq_flux, rflx, qr_tend, nr_tend,
    precip_liq_surf (accumulated)."""
    qr, nr = jnp.asarray(qr), jnp.asarray(nr)
    qr_incld, nr_incld = jnp.asarray(qr_incld), jnp.asarray(nr_incld)
    mu_r, lamr = jnp.asarray(mu_r), jnp.asarray(lamr)
    rho, inv_rho = jnp.asarray(rho), jnp.asarray(inv_rho)
    rhofacr, cld_frac_r = jnp.asarray(rhofacr), jnp.asarray(cld_frac_r)
    inv_dz = jnp.asarray(inv_dz)
    precip_liq_flux = jnp.asarray(precip_liq_flux_in)
    rflx = jnp.asarray(rflx_in)

    nlev = qr.shape[-1]
    kbot = nlev - 1
    idx = jnp.arange(nlev)

    k_qxtop, k_qxbot0, present = find_top_bottom(qr, c.qsmall)
    dt_left0 = jnp.where(present, dt, 0.0)
    prt0 = jnp.zeros_like(dt_left0)

    def cond(carry):
        return jnp.any(carry[7] > c.dt_left_tol)

    def body(carry):
        (qr, nr, qr_incld, nr_incld, mu_r, lamr, precip_liq_flux,
         dt_left, prt_accum, k_qxbot, rflx) = carry
        active = dt_left > c.dt_left_tol
        co_band = ((idx >= k_qxtop[..., None]) & (idx <= k_qxbot[..., None])
                   & active[..., None])
        m = co_band & (qr_incld > c.qsmall)

        nr_incld, mu_r, lamr, V_qr, V_nr = compute_rain_fall_velocity(
            vn_table_vals, vm_table_vals, qr_incld, rhofacr,
            p3_max_mean_rain_size, nr_incld, mu_r, lamr, m)
        nr = jnp.where(m, nr_incld * cld_frac_r, nr)

        Co_max = jnp.max(
            jnp.where(co_band, V_qr * dt_left[..., None] * inv_dz, 0.0),
            axis=-1)

        fields, fluxes, k_qxbot, dt_left, prt_accum, dt_sub = \
            _generalized_substep(
                (qr, nr), (V_qr, V_nr), rho, inv_rho, inv_dz, idx,
                k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, active)
        qr, nr = fields
        precip_liq_flux = _accum_interface_flux(
            precip_liq_flux, fluxes[0], dt_sub, k_qxtop, k_qxbot, active)
        rflx = _accum_interface_flux(rflx, fluxes[0], dt_sub,
                                     k_qxtop, k_qxbot, active)

        act = active[..., None]
        qr_incld = jnp.where(act, qr / cld_frac_r, qr_incld)
        nr_incld = jnp.where(act, nr / cld_frac_r, nr_incld)
        return (qr, nr, qr_incld, nr_incld, mu_r, lamr, precip_liq_flux,
                dt_left, prt_accum, k_qxbot, rflx)

    carry = (qr, nr, qr_incld, nr_incld, mu_r, lamr, precip_liq_flux,
             dt_left0, prt0, k_qxbot0, rflx)
    (qr, nr, qr_incld, nr_incld, mu_r, lamr, precip_liq_flux,
     _, prt_accum, _, rflx) = lax.while_loop(cond, body, carry)

    pm = present[..., None]
    precip_liq_flux = jnp.where(pm, precip_liq_flux * inv_dt,
                                precip_liq_flux)
    rflx = jnp.where(pm, rflx * inv_dt, rflx)
    precip_liq_surf = (jnp.asarray(precip_liq_surf_in)
                       + jnp.where(present,
                                   prt_accum * c.inv_rho_h2o * inv_dt, 0.0))
    return {
        "qr": qr, "nr": nr, "qr_incld": qr_incld, "nr_incld": nr_incld,
        "mu_r": mu_r, "lamr": lamr, "precip_liq_flux": precip_liq_flux,
        "rflx": rflx,
        "qr_tend": (qr - jnp.asarray(qr_tend_in)) * inv_dt,
        "nr_tend": (nr - jnp.asarray(nr_tend_in)) * inv_dt,
        "precip_liq_surf": precip_liq_surf,
    }


@jax.jit
def ice_sedimentation(rho, inv_rho, rhofaci, cld_frac_i, inv_dz,
                      dt, inv_dt,
                      qi, qi_incld, ni, ni_incld, qm, qm_incld,
                      bm, bm_incld, ice_table_vals,
                      precip_ice_surf_in, sflx_in, qi_tend_in, ni_tend_in):
    """ice_sedimentation. Returns a dict with qi, ni, qm, bm, their
    in-cloud values, sflx, qi_tend, ni_tend, precip_ice_surf
    (accumulated). precip_ice_flux is untouched by the Fortran and
    therefore not an output here."""
    qi, ni = jnp.asarray(qi), jnp.asarray(ni)
    qm, bm = jnp.asarray(qm), jnp.asarray(bm)
    qi_incld, ni_incld = jnp.asarray(qi_incld), jnp.asarray(ni_incld)
    qm_incld, bm_incld = jnp.asarray(qm_incld), jnp.asarray(bm_incld)
    rho, inv_rho = jnp.asarray(rho), jnp.asarray(inv_rho)
    rhofaci, cld_frac_i = jnp.asarray(rhofaci), jnp.asarray(cld_frac_i)
    inv_dz = jnp.asarray(inv_dz)
    sflx = jnp.asarray(sflx_in)

    nlev = qi.shape[-1]
    kbot = nlev - 1
    idx = jnp.arange(nlev)

    k_qxtop, k_qxbot0, present = find_top_bottom(qi, c.qsmall)
    dt_left0 = jnp.where(present, dt, 0.0)
    prt0 = jnp.zeros_like(dt_left0)

    def cond(carry):
        return jnp.any(carry[8] > c.dt_left_tol)

    def body(carry):
        (qi, ni, qm, bm, qi_incld, ni_incld, qm_incld, bm_incld,
         dt_left, prt_accum, k_qxbot, sflx) = carry
        active = dt_left > c.dt_left_tol
        co_band = ((idx >= k_qxtop[..., None]) & (idx <= k_qxbot[..., None])
                   & active[..., None])
        m = co_band & (qi_incld > c.qsmall)

        ni_incld = jnp.where(m, jnp.maximum(ni_incld, c.nsmall), ni_incld)
        rhop, qm2, bm2 = calc_bulk_rho_rime(qi_incld, qm_incld, bm_incld, m)
        qm_incld = jnp.where(m, qm2, qm_incld)
        bm_incld = jnp.where(m, bm2, bm_incld)
        qm = jnp.where(m, qm_incld * cld_frac_i, qm)
        bm = jnp.where(m, bm_incld * cld_frac_i, bm)

        ti = lookup_ice(qi_incld, ni_incld, qm_incld, rhop, m)
        t_ni_fallspd = jnp.where(m, apply_table_ice(0, ice_table_vals, ti),
                                 0.0)
        t_qi_fallspd = jnp.where(m, apply_table_ice(1, ice_table_vals, ti),
                                 0.0)
        t_ni_lammax = jnp.where(m, apply_table_ice(6, ice_table_vals, ti),
                                0.0)
        t_ni_lammin = jnp.where(m, apply_table_ice(7, ice_table_vals, ti),
                                0.0)
        ni_incld = jnp.where(m, jnp.minimum(ni_incld,
                                            t_ni_lammax * ni_incld),
                             ni_incld)
        ni_incld = jnp.where(m, jnp.maximum(ni_incld,
                                            t_ni_lammin * ni_incld),
                             ni_incld)
        ni = jnp.where(m, ni_incld * cld_frac_i, ni)

        V_qit = jnp.where(m, t_qi_fallspd * rhofaci, 0.0)
        V_nit = jnp.where(m, t_ni_fallspd * rhofaci, 0.0)

        Co_max = jnp.max(
            jnp.where(co_band, V_qit * dt_left[..., None] * inv_dz, 0.0),
            axis=-1)

        fields, fluxes, k_qxbot, dt_left, prt_accum, dt_sub = \
            _generalized_substep(
                (qi, ni, qm, bm), (V_qit, V_nit, V_qit, V_qit),
                rho, inv_rho, inv_dz, idx,
                k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, active)
        qi, ni, qm, bm = fields
        sflx = _accum_interface_flux(sflx, fluxes[0], dt_sub,
                                     k_qxtop, k_qxbot, active)

        act = active[..., None]
        qi_incld = jnp.where(act, qi / cld_frac_i, qi_incld)
        ni_incld = jnp.where(act, ni / cld_frac_i, ni_incld)
        qm_incld = jnp.where(act, qm / cld_frac_i, qm_incld)
        bm_incld = jnp.where(act, bm / cld_frac_i, bm_incld)
        return (qi, ni, qm, bm, qi_incld, ni_incld, qm_incld, bm_incld,
                dt_left, prt_accum, k_qxbot, sflx)

    carry = (qi, ni, qm, bm, qi_incld, ni_incld, qm_incld, bm_incld,
             dt_left0, prt0, k_qxbot0, sflx)
    (qi, ni, qm, bm, qi_incld, ni_incld, qm_incld, bm_incld,
     _, prt_accum, _, sflx) = lax.while_loop(cond, body, carry)

    sflx = jnp.where(present[..., None], sflx * inv_dt, sflx)
    precip_ice_surf = (jnp.asarray(precip_ice_surf_in)
                       + jnp.where(present,
                                   prt_accum * c.inv_rho_h2o * inv_dt, 0.0))
    return {
        "qi": qi, "ni": ni, "qm": qm, "bm": bm,
        "qi_incld": qi_incld, "ni_incld": ni_incld,
        "qm_incld": qm_incld, "bm_incld": bm_incld, "sflx": sflx,
        "qi_tend": (qi - jnp.asarray(qi_tend_in)) * inv_dt,
        "ni_tend": (ni - jnp.asarray(ni_tend_in)) * inv_dt,
        "precip_ice_surf": precip_ice_surf,
    }


@jax.jit
def homogeneous_freezing(exner, latent_heat_fusion, qc, nc, qr, nr,
                         qi, ni, qm, bm, th):
    """Instantaneous freezing of cloud water and rain below T_homogfrz
    (homogeneous_freezing, EAM form: tests th/exner, SEQUENTIALLY —
    the qr test sees the theta warming from the qc freezing).
    Returns a dict with qc, nc, qr, nr, qi, ni, qm, bm, th."""
    exner = jnp.asarray(exner)
    qc, nc, qr, nr = (jnp.asarray(a) for a in (qc, nc, qr, nr))
    qi, ni, qm, bm = (jnp.asarray(a) for a in (qi, ni, qm, bm))
    th = jnp.asarray(th)

    qc_ge = (qc >= c.qsmall) & ((th / exner) < c.T_homogfrz)
    Qc_nuc = qc
    Nc_nuc = jnp.maximum(nc, c.nsmall)
    qm = jnp.where(qc_ge, qm + Qc_nuc, qm)
    qi = jnp.where(qc_ge, qi + Qc_nuc, qi)
    bm = jnp.where(qc_ge, bm + Qc_nuc * c.inv_rho_rimeMax, bm)
    ni = jnp.where(qc_ge, ni + Nc_nuc, ni)
    th = jnp.where(qc_ge,
                   th + exner * Qc_nuc * latent_heat_fusion * c.inv_cp, th)
    qc = jnp.where(qc_ge, 0.0, qc)
    nc = jnp.where(qc_ge, 0.0, nc)

    # rain test uses the UPDATED theta
    qr_ge = (qr >= c.qsmall) & ((th / exner) < c.T_homogfrz)
    Qr_nuc = qr
    Nr_nuc = jnp.maximum(nr, c.nsmall)
    qm = jnp.where(qr_ge, qm + Qr_nuc, qm)
    qi = jnp.where(qr_ge, qi + Qr_nuc, qi)
    bm = jnp.where(qr_ge, bm + Qr_nuc * c.inv_rho_rimeMax, bm)
    ni = jnp.where(qr_ge, ni + Nr_nuc, ni)
    th = jnp.where(qr_ge,
                   th + exner * Qr_nuc * latent_heat_fusion * c.inv_cp, th)
    qr = jnp.where(qr_ge, 0.0, qr)
    nr = jnp.where(qr_ge, 0.0, nr)

    return {"qc": qc, "nc": nc, "qr": qr, "nr": nr,
            "qi": qi, "ni": ni, "qm": qm, "bm": bm, "th": th}


@jax.jit
def ice_complete_melting(exner, latent_heat_fusion, qi, ni, qm,
                         qr, nr, qc, nc, th):
    """Instantaneous (fractional) melting of ice at th/exner >
    t_snow_melt = f32(273.15)+2 (ice_complete_melting; EAM-only).
    Small unrimed crystals melt to cloud droplets, otherwise to rain.
    Returns a dict with qi, ni, qm, qr, nr, qc, nc, th (bm untouched,
    as in the Fortran)."""
    exner = jnp.asarray(exner)
    qi, ni, qm = jnp.asarray(qi), jnp.asarray(ni), jnp.asarray(qm)
    qr, nr = jnp.asarray(qr), jnp.asarray(nr)
    qc, nc = jnp.asarray(qc), jnp.asarray(nc)
    th = jnp.asarray(th)

    t_here = th / exner
    melt = (qi >= c.qsmall) & (t_here > c.t_snow_melt)

    del_mass = qi
    del_num = jnp.maximum(ni, c.nsmall)

    partial = t_here - del_mass * latent_heat_fusion / c.cp < c.t_snow_melt
    del_safe = jnp.where(melt & (del_mass != 0), del_mass, 1.0)
    equiv_mass = (t_here - c.t_snow_melt) * c.cp / latent_heat_fusion
    frac_mlt = jnp.where(partial,
                         jnp.clip(equiv_mass / del_safe, 0.0, 1.0), 1.0)

    # crystal volume radius [um]; 3*qi/ni/4/pi/900 with RAW ni exactly
    # as the Fortran: ni=0 -> inf -> rain branch; ni<0 -> pow-cbrt NaN
    # -> rain branch (rv<100 false)
    rv_tmp = 3.0 * qi / jnp.where(melt, ni, 1.0) / 4.0 / c.pi / 900.0
    rv = 1.0e6 * rv_tmp ** (1.0 / 3.0)

    qi_safe = jnp.where(melt, qi, 1.0)
    lightly_rimed = (qm / qi_safe) < 0.1
    to_cloud = melt & lightly_rimed & (rv < 100.0)
    to_rain = melt & ~to_cloud

    qi_new = jnp.maximum((1.0 - frac_mlt) * del_mass, 0.0)
    qm_new = jnp.maximum((1.0 - frac_mlt) * del_mass, 0.0)
    ni_new = jnp.maximum((1.0 - frac_mlt) * del_num, 0.0)

    qc = jnp.where(to_cloud, qc + frac_mlt * del_mass, qc)
    nc = jnp.where(to_cloud, nc + frac_mlt * del_num, nc)
    qr = jnp.where(to_rain, qr + frac_mlt * del_mass, qr)
    nr = jnp.where(to_rain, nr + frac_mlt * del_num, nr)
    qi = jnp.where(melt, qi_new, qi)
    qm = jnp.where(melt, qm_new, qm)
    ni = jnp.where(melt, ni_new, ni)
    th = jnp.where(melt,
                   th - exner * frac_mlt * del_mass * latent_heat_fusion
                   / c.cp, th)
    return {"qi": qi, "ni": ni, "qm": qm, "qr": qr, "nr": nr,
            "qc": qc, "nc": nc, "th": th}
