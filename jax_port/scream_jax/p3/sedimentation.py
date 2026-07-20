"""P3 sedimentation: adaptive-substep first-order upwind transport.

Sources (components/eamxx/src/physics/p3/impl/):
  p3_find_impl.hpp (find_top / find_bottom)
  p3_upwind_impl.hpp (calc_first_order_upwind_step, generalized_sedimentation)
  p3_cloud_sed_impl.hpp, p3_rain_sed_impl.hpp, p3_ice_sed_impl.hpp
  (drivers + compute_rain_fall_velocity + homogeneous_freezing)

EAMxx always calls these with kdir = -1 (k=0 at model top, kbot = nlev-1 at
the surface); only that orientation is implemented. The C++ per-column
while(dt_left > tol) CFL substepping has two bit-identical realizations
selected by the static `use_while_loop` flag:

- use_while_loop=True (DEFAULT — the fast primal path for production
  runs without autodiff): the direct lax.while_loop transcription.
  NOT reverse-mode differentiable (jax.grad raises).
- use_while_loop=False (SWITCH TO THIS FOR AUTODIFF, e.g. via
  p3 opts["sed_use_while_loop"] = False): a fixed-length lax.scan over
  batched columns; every step selects per column, with a pure jnp.where
  on each carry component, between the substep result and the incoming
  carry, so columns whose dt_left is already spent pass through bitwise
  unchanged and the result is bit-identical to the while_loop as long
  as the static bound MAX_SEDI_SUBSTEPS_* covers the actual trip count
  (~5.5x slower primal at the default 4x-margin bounds).

Each returned dict carries a per-column "converged" flag (dt_left <= tol
at the end) guarding against silent truncation in scan mode. The
moving active band [k_qxtop, k_qxbot] is realized as boolean masks
against a level-index array; levels outside the band keep their values,
and the zero incoming flux at the band top falls out of fluxes being
zero outside the band.

Columns with no condensate (log_qxpresent false) start with dt_left = 0 so
the loop body never touches them, matching the C++ early exit.
"""

import functools

import jax
import jax.numpy as jnp
from jax import lax

from ..foundation import constants as c
from ..foundation import smoothing
from .dsd import get_cloud_dsd2, get_rain_dsd2, _tgamma
from .table_lookups import lookup_table3, apply_table3, lookup_ice, apply_table_ice
from .main_part3 import calc_bulk_rho_rime

# Static substep bounds for the masked-scan CFL loops. Measured maximum
# GLOBAL trip counts (iterations until every column's dt is spent) on the
# golden states — p3_218x72_dt1800_5steps (218 cols, 5 steps, dt=1800),
# physics_suite_218x72_dt1800_2steps states run at dt=300 and dt=1800,
# and the property-test inputs incl. the thin-layer CFL stress case:
#   cloud:  6  (dt=1800; dt=300: 1, thin-layer dt=1800 stress: 6)
#   rain: 421  (dt=1800; dt=300: 109)
#   ice:  297  (dt=1800; dt=300:  64)
# Bounds are ~4x the observed max (minimum 32). Truncation is guarded by
# the returned "converged" flag and asserted in tests.
MAX_SEDI_SUBSTEPS_CLOUD = 32
MAX_SEDI_SUBSTEPS_RAIN = 1684
MAX_SEDI_SUBSTEPS_ICE = 1188


def _masked_substep_scan(body, carry, length):
    """Fixed-length, reverse-differentiable replacement for the
    while(dt_left > tol) substep loops (dt_left is carry[-3]).

    Runs `body` `length` times under lax.scan. Each step selects, per
    column, between the body's result and the incoming carry with a pure
    jnp.where on every carry component: columns whose dt_left is already
    spent pass through bitwise unchanged, so the result is bit-identical
    to the lax.while_loop whenever `length` covers the actual trip count
    (see MAX_SEDI_SUBSTEPS_*). The step is wrapped in jax.checkpoint so
    reverse-mode memory stays linear in `length` with cheap recompute.
    """
    def step(carry, _):
        active = carry[-3] > c.dt_left_tol
        new = body(carry)
        merged = tuple(
            jnp.where(
                active.reshape(active.shape + (1,) * (n.ndim - active.ndim)),
                n, o)
            for n, o in zip(new, carry))
        return merged, None

    final, _ = lax.scan(jax.checkpoint(step), carry, None, length=length)
    return final


def assert_converged(sed_out, name="sedimentation"):
    """Host-side guard against silent substep truncation (tests/debug):
    raises if any column's CFL loop did not spend its full dt within the
    static substep bound."""
    import numpy as np
    conv = np.asarray(sed_out["converged"])
    if not conv.all():
        raise RuntimeError(
            f"{name}: {int((~conv).sum())} column(s) did not exhaust dt "
            "within MAX_SEDI_SUBSTEPS; raise the bound")


def find_top_bottom(q, small):
    """find_top + find_bottom for kdir = -1 (Functions::find_top/find_bottom).

    Returns (k_qxtop, k_qxbot, present): the smallest and largest level
    index with q >= small, and whether any level qualifies. Where no level
    qualifies the indices are 0 (unused, as in the C++).
    """
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
    (Functions::calc_first_order_upwind_step, kdir = -1).

    fields/velocities are tuples of (..., nlev) arrays. Returns (updated
    fields tuple, fluxes tuple). Incoming flux at the band top is zero
    because fluxes vanish outside the band.
    """
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
    """One iteration of Functions::generalized_sedimentation (kdir = -1).

    Returns (fields, fluxes, k_qxbot, dt_left, prt_accum). The bottom of
    the upwind band extends one level below k_qxbot unless already at the
    ground; when at the ground the surface flux accumulates into
    prt_accum, otherwise the band grows downward.
    """
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
    return fields, fluxes, k_qxbot, dt_left, prt_accum


@functools.partial(jax.jit, static_argnames=("do_predict_nc", "max_substeps",
                                             "use_while_loop"))
def cloud_sedimentation(qc_incld, rho, inv_rho, cld_frac_l, acn, inv_dz,
                        dt, inv_dt, do_predict_nc: bool,
                        qc, nc, nc_incld, mu_c, lamc,
                        qc_tend_in, nc_tend_in, precip_liq_surf_in,
                        max_substeps=None, use_while_loop=True):
    """Functions::cloud_sedimentation. Returns a dict with qc, nc,
    qc_incld, nc_incld, mu_c, lamc, qc_tend, nc_tend, precip_liq_surf,
    converged (per-column: dt spent within the substep bound).

    precip_liq_surf is SET (not accumulated) where cloud water is present,
    as in the C++."""
    qc, nc = jnp.asarray(qc), jnp.asarray(nc)
    qc_incld, nc_incld = jnp.asarray(qc_incld), jnp.asarray(nc_incld)
    mu_c, lamc = jnp.asarray(mu_c), jnp.asarray(lamc)
    rho, inv_rho = jnp.asarray(rho), jnp.asarray(inv_rho)
    cld_frac_l = jnp.asarray(cld_frac_l)
    acn, inv_dz = jnp.asarray(acn), jnp.asarray(inv_dz)

    nlev = qc.shape[-1]
    kbot = nlev - 1
    idx = jnp.arange(nlev)

    k_qxtop, k_qxbot0, present = find_top_bottom(qc, c.QSMALL)
    dt_left0 = jnp.where(present, dt, 0.0)
    prt0 = jnp.zeros_like(dt_left0)

    def cond(carry):
        return jnp.any(carry[-3] > c.dt_left_tol)

    def body(carry):
        qc, nc, qc_incld, nc_incld, mu_c, lamc, dt_left, prt_accum, k_qxbot = carry
        active = dt_left > c.dt_left_tol
        co_band = ((idx >= k_qxtop[..., None]) & (idx <= k_qxbot[..., None])
                   & active[..., None])
        m = co_band & (qc_incld > c.QSMALL)

        nc2, mu2, _, lam2, _, _ = get_cloud_dsd2(qc_incld, nc_incld, rho, m)
        nc_incld = jnp.where(m, nc2, nc_incld)
        mu_c = jnp.where(m, mu2, mu_c)
        lamc = jnp.where(m, lam2, lamc)
        nc = jnp.where(m, nc_incld * cld_frac_l, nc)

        lam_safe = jnp.where(m, lamc, 1.0)
        dum = 1.0 / lam_safe ** c.bcn
        V_qc = jnp.where(m, acn * _tgamma(4.0 + c.bcn + mu_c) * dum
                         / _tgamma(mu_c + 4.0), 0.0)
        if do_predict_nc:
            V_nc = jnp.where(m, acn * _tgamma(1.0 + c.bcn + mu_c) * dum
                             / _tgamma(mu_c + 1.0), 0.0)
            fields, velocities = (qc, nc), (V_qc, V_nc)
        else:
            fields, velocities = (qc,), (V_qc,)

        Co_max = jnp.max(
            jnp.where(m, V_qc * dt_left[..., None] * inv_dz, 0.0), axis=-1)

        fields, _, k_qxbot, dt_left, prt_accum = _generalized_substep(
            fields, velocities, rho, inv_rho, inv_dz, idx,
            k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, active)
        qc = fields[0]
        if do_predict_nc:
            nc = fields[1]

        act = active[..., None]
        qc_incld = jnp.where(act, qc / cld_frac_l, qc_incld)
        nc_incld = jnp.where(act, nc / cld_frac_l, nc_incld)
        return qc, nc, qc_incld, nc_incld, mu_c, lamc, dt_left, prt_accum, k_qxbot

    carry = (qc, nc, qc_incld, nc_incld, mu_c, lamc,
             dt_left0, prt0, k_qxbot0)
    if use_while_loop:
        final = lax.while_loop(cond, body, carry)
    else:
        n = MAX_SEDI_SUBSTEPS_CLOUD if max_substeps is None else max_substeps
        final = _masked_substep_scan(body, carry, n)
    (qc, nc, qc_incld, nc_incld, mu_c, lamc,
     dt_left, prt_accum, _) = final

    precip_liq_surf = jnp.where(
        present, prt_accum * c.INV_RHO_H2O * inv_dt,
        jnp.asarray(precip_liq_surf_in))
    return {
        "qc": qc, "nc": nc, "qc_incld": qc_incld, "nc_incld": nc_incld,
        "mu_c": mu_c, "lamc": lamc,
        "qc_tend": (qc - jnp.asarray(qc_tend_in)) * inv_dt,
        "nc_tend": (nc - jnp.asarray(nc_tend_in)) * inv_dt,
        "precip_liq_surf": precip_liq_surf,
        "converged": dt_left <= c.dt_left_tol,
    }


def compute_rain_fall_velocity(vn_table_vals, vm_table_vals, qr_incld,
                               rhofacr, nr_incld, mu_r, lamr,
                               constant_mu_rain, context):
    """Functions::compute_rain_fall_velocity. Returns
    (nr_incld, mu_r, lamr, V_qr, V_nr) with the persistent arrays updated
    only under context."""
    nr2, mu2, lam2 = get_rain_dsd2(qr_incld, nr_incld, constant_mu_rain,
                                   context)
    nr_incld = jnp.where(context, nr2, nr_incld)
    mu_r = jnp.where(context, mu2, mu_r)
    lamr = jnp.where(context, lam2, lamr)
    tab = lookup_table3(mu_r, lamr, context)
    V_qr = jnp.where(context, apply_table3(vm_table_vals, tab) * rhofacr, 0.0)
    V_nr = jnp.where(context, apply_table3(vn_table_vals, tab) * rhofacr, 0.0)
    return nr_incld, mu_r, lamr, V_qr, V_nr


@functools.partial(jax.jit, static_argnames=("max_substeps",
                                             "use_while_loop"))
def rain_sedimentation(rho, inv_rho, rhofacr, cld_frac_r, inv_dz, qr_incld,
                       vn_table_vals, vm_table_vals, dt, inv_dt,
                       qr, nr, nr_incld, mu_r, lamr,
                       precip_liq_flux_in, qr_tend_in, nr_tend_in,
                       precip_liq_surf_in, opts,
                       max_substeps=None, use_while_loop=True):
    """Functions::rain_sedimentation. Returns a dict with qr, nr,
    qr_incld, nr_incld, mu_r, lamr, precip_liq_flux (nlev+1 interfaces),
    qr_tend, nr_tend, precip_liq_surf (accumulated: in + contribution),
    converged (per-column: dt spent within the substep bound)."""
    qr, nr = jnp.asarray(qr), jnp.asarray(nr)
    qr_incld, nr_incld = jnp.asarray(qr_incld), jnp.asarray(nr_incld)
    mu_r, lamr = jnp.asarray(mu_r), jnp.asarray(lamr)
    rho, inv_rho = jnp.asarray(rho), jnp.asarray(inv_rho)
    rhofacr, cld_frac_r = jnp.asarray(rhofacr), jnp.asarray(cld_frac_r)
    inv_dz = jnp.asarray(inv_dz)
    precip_liq_flux = jnp.asarray(precip_liq_flux_in)

    nlev = qr.shape[-1]
    kbot = nlev - 1
    idx = jnp.arange(nlev)
    idx_f = jnp.arange(nlev + 1)

    k_qxtop, k_qxbot0, present = find_top_bottom(qr, c.QSMALL)
    dt_left0 = jnp.where(present, dt, 0.0)
    prt0 = jnp.zeros_like(dt_left0)

    def cond(carry):
        return jnp.any(carry[-3] > c.dt_left_tol)

    def body(carry):
        (qr, nr, qr_incld, nr_incld, mu_r, lamr, precip_liq_flux,
         dt_left, prt_accum, k_qxbot) = carry
        active = dt_left > c.dt_left_tol
        co_band = ((idx >= k_qxtop[..., None]) & (idx <= k_qxbot[..., None])
                   & active[..., None])
        m = co_band & (qr_incld > c.QSMALL)

        nr_incld, mu_r, lamr, V_qr, V_nr = compute_rain_fall_velocity(
            vn_table_vals, vm_table_vals, qr_incld, rhofacr,
            nr_incld, mu_r, lamr, opts["constant_mu_rain"], m)
        # keep nr consistent with the dsd-limited nr_incld
        nr = jnp.where(m, nr_incld * cld_frac_r, nr)

        Co_max = jnp.max(
            jnp.where(m, V_qr * dt_left[..., None] * inv_dz, 0.0), axis=-1)

        fields, fluxes, k_qxbot, dt_left, prt_accum = _generalized_substep(
            (qr, nr), (V_qr, V_nr), rho, inv_rho, inv_dz, idx,
            k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, active)
        qr, nr = fields

        act = active[..., None]
        qr_incld = jnp.where(act, qr / cld_frac_r, qr_incld)
        nr_incld = jnp.where(act, nr / cld_frac_r, nr_incld)

        # accumulate the qr flux on interfaces [k_qxtop+1, k_qxbot+1]
        # (k_qxbot after the band extension above)
        flux_qx = fluxes[0]
        fmask = ((idx_f >= (k_qxtop + 1)[..., None])
                 & (idx_f <= (k_qxbot + 1)[..., None]) & act)
        flux_shift = jnp.concatenate([flux_qx[..., :1], flux_qx], axis=-1)
        precip_liq_flux = jnp.where(fmask, precip_liq_flux + flux_shift,
                                    precip_liq_flux)
        return (qr, nr, qr_incld, nr_incld, mu_r, lamr, precip_liq_flux,
                dt_left, prt_accum, k_qxbot)

    carry = (qr, nr, qr_incld, nr_incld, mu_r, lamr, precip_liq_flux,
             dt_left0, prt0, k_qxbot0)
    if use_while_loop:
        final = lax.while_loop(cond, body, carry)
    else:
        n = MAX_SEDI_SUBSTEPS_RAIN if max_substeps is None else max_substeps
        final = _masked_substep_scan(body, carry, n)
    (qr, nr, qr_incld, nr_incld, mu_r, lamr, precip_liq_flux,
     dt_left, prt_accum, _) = final

    precip_liq_surf = (jnp.asarray(precip_liq_surf_in)
                       + jnp.where(present,
                                   prt_accum * c.INV_RHO_H2O * inv_dt, 0.0))
    return {
        "qr": qr, "nr": nr, "qr_incld": qr_incld, "nr_incld": nr_incld,
        "mu_r": mu_r, "lamr": lamr, "precip_liq_flux": precip_liq_flux,
        "qr_tend": (qr - jnp.asarray(qr_tend_in)) * inv_dt,
        "nr_tend": (nr - jnp.asarray(nr_tend_in)) * inv_dt,
        "precip_liq_surf": precip_liq_surf,
        "converged": dt_left <= c.dt_left_tol,
    }


@functools.partial(jax.jit, static_argnames=("max_substeps",
                                             "use_while_loop"))
def ice_sedimentation(rho, inv_rho, rhofaci, cld_frac_i, inv_dz,
                      dt, inv_dt,
                      qi, qi_incld, ni, ni_incld, qm, qm_incld,
                      bm, bm_incld, ice_table_vals,
                      qi_tend_in, ni_tend_in, precip_ice_surf_in, opts,
                      max_substeps=None, use_while_loop=True):
    """Functions::ice_sedimentation. Returns a dict with qi, ni, qm, bm,
    their in-cloud values, qi_tend, ni_tend, precip_ice_surf (accumulated),
    converged (per-column: dt spent within the substep bound).
    qm and bm advect with the mass-weighted velocity V_qit."""
    qi, ni = jnp.asarray(qi), jnp.asarray(ni)
    qm, bm = jnp.asarray(qm), jnp.asarray(bm)
    qi_incld, ni_incld = jnp.asarray(qi_incld), jnp.asarray(ni_incld)
    qm_incld, bm_incld = jnp.asarray(qm_incld), jnp.asarray(bm_incld)
    rho, inv_rho = jnp.asarray(rho), jnp.asarray(inv_rho)
    rhofaci, cld_frac_i = jnp.asarray(rhofaci), jnp.asarray(cld_frac_i)
    inv_dz = jnp.asarray(inv_dz)

    nlev = qi.shape[-1]
    kbot = nlev - 1
    idx = jnp.arange(nlev)

    k_qxtop, k_qxbot0, present = find_top_bottom(qi, c.QSMALL)
    dt_left0 = jnp.where(present, dt, 0.0)
    prt0 = jnp.zeros_like(dt_left0)

    def cond(carry):
        return jnp.any(carry[-3] > c.dt_left_tol)

    def body(carry):
        (qi, ni, qm, bm, qi_incld, ni_incld, qm_incld, bm_incld,
         dt_left, prt_accum, k_qxbot) = carry
        active = dt_left > c.dt_left_tol
        co_band = ((idx >= k_qxtop[..., None]) & (idx <= k_qxbot[..., None])
                   & active[..., None])
        m = co_band & (qi_incld > c.QSMALL)

        ni_incld = jnp.where(m, jnp.maximum(ni_incld, c.NSMALL), ni_incld)
        rhop, qm2, bm2 = calc_bulk_rho_rime(qi_incld, qm_incld, bm_incld,
                                            opts, m)
        qm_incld = jnp.where(m, qm2, qm_incld)
        bm_incld = jnp.where(m, bm2, bm_incld)
        qm = jnp.where(m, qm_incld * cld_frac_i, qm)
        bm = jnp.where(m, bm_incld * cld_frac_i, bm)

        ti = lookup_ice(qi_incld, ni_incld, qm_incld, rhop, m)
        t_ni_fallspd = jnp.where(m, apply_table_ice(0, ice_table_vals, ti), 0.0)
        t_qi_fallspd = jnp.where(m, apply_table_ice(1, ice_table_vals, ti), 0.0)
        t_ni_lammax = jnp.where(m, apply_table_ice(6, ice_table_vals, ti), 0.0)
        t_ni_lammin = jnp.where(m, apply_table_ice(7, ice_table_vals, ti), 0.0)
        ni_incld = jnp.where(m, jnp.minimum(ni_incld, t_ni_lammax * ni_incld),
                             ni_incld)
        ni_incld = jnp.where(m, jnp.maximum(ni_incld, t_ni_lammin * ni_incld),
                             ni_incld)
        ni = jnp.where(m, ni_incld * cld_frac_i, ni)

        factor = opts["ice_sedimentation_factor"]
        V_qit = jnp.where(m, factor * t_qi_fallspd * rhofaci, 0.0)
        V_nit = jnp.where(m, factor * t_ni_fallspd * rhofaci, 0.0)

        Co_max = jnp.max(
            jnp.where(m, V_qit * dt_left[..., None] * inv_dz, 0.0), axis=-1)

        fields, _, k_qxbot, dt_left, prt_accum = _generalized_substep(
            (qi, ni, qm, bm), (V_qit, V_nit, V_qit, V_qit),
            rho, inv_rho, inv_dz, idx,
            k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, active)
        qi, ni, qm, bm = fields

        act = active[..., None]
        qi_incld = jnp.where(act, qi / cld_frac_i, qi_incld)
        ni_incld = jnp.where(act, ni / cld_frac_i, ni_incld)
        qm_incld = jnp.where(act, qm / cld_frac_i, qm_incld)
        bm_incld = jnp.where(act, bm / cld_frac_i, bm_incld)
        return (qi, ni, qm, bm, qi_incld, ni_incld, qm_incld, bm_incld,
                dt_left, prt_accum, k_qxbot)

    carry = (qi, ni, qm, bm, qi_incld, ni_incld, qm_incld, bm_incld,
             dt_left0, prt0, k_qxbot0)
    if use_while_loop:
        final = lax.while_loop(cond, body, carry)
    else:
        n = MAX_SEDI_SUBSTEPS_ICE if max_substeps is None else max_substeps
        final = _masked_substep_scan(body, carry, n)
    (qi, ni, qm, bm, qi_incld, ni_incld, qm_incld, bm_incld,
     dt_left, prt_accum, _) = final

    precip_ice_surf = (jnp.asarray(precip_ice_surf_in)
                       + jnp.where(present,
                                   prt_accum * c.INV_RHO_H2O * inv_dt, 0.0))
    return {
        "qi": qi, "ni": ni, "qm": qm, "bm": bm,
        "qi_incld": qi_incld, "ni_incld": ni_incld,
        "qm_incld": qm_incld, "bm_incld": bm_incld,
        "qi_tend": (qi - jnp.asarray(qi_tend_in)) * inv_dt,
        "ni_tend": (ni - jnp.asarray(ni_tend_in)) * inv_dt,
        "precip_ice_surf": precip_ice_surf,
        "converged": dt_left <= c.dt_left_tol,
    }


@functools.partial(jax.jit, static_argnames=("smooth_width",))
def homogeneous_freezing(T_atm, inv_exner, qc, nc, qr, nr, qi, ni, qm, bm,
                         th_atm, smooth_width=0.0):
    """Instantaneous freezing of all cloud water and rain below T_homogfrz
    (Functions::homogeneous_freezing). Returns a dict with the updated
    qc, nc, qr, nr, qi, ni, qm, bm, th_atm."""
    T_atm = jnp.asarray(T_atm)
    inv_exner = jnp.asarray(inv_exner)
    qc, nc, qr, nr = (jnp.asarray(a) for a in (qc, nc, qr, nr))
    qi, ni, qm, bm = (jnp.asarray(a) for a in (qi, ni, qm, bm))
    th_atm = jnp.asarray(th_atm)

    Qc_nuc, Qr_nuc = qc, qr
    Nc_nuc = jnp.maximum(nc, c.NSMALL)
    Nr_nuc = jnp.maximum(nr, c.NSMALL)

    if smooth_width == 0.0:
        t_lt = T_atm < c.T_homogfrz
        qc_ge = t_lt & (qc >= c.QSMALL)
        qr_ge = t_lt & (qr >= c.QSMALL)

        qm = jnp.where(qc_ge, qm + Qc_nuc, qm)
        qi = jnp.where(qc_ge, qi + Qc_nuc, qi)
        bm = jnp.where(qc_ge, bm + Qc_nuc * c.INV_RHO_RIMEMAX, bm)
        ni = jnp.where(qc_ge, ni + Nc_nuc, ni)
        th_atm = jnp.where(qc_ge,
                           th_atm + inv_exner * Qc_nuc * c.LatIce * c.INV_CP,
                           th_atm)

        qm = jnp.where(qr_ge, qm + Qr_nuc, qm)
        qi = jnp.where(qr_ge, qi + Qr_nuc, qi)
        bm = jnp.where(qr_ge, bm + Qr_nuc * c.INV_RHO_RIMEMAX, bm)
        ni = jnp.where(qr_ge, ni + Nr_nuc, ni)
        th_atm = jnp.where(qr_ge,
                           th_atm + inv_exner * Qr_nuc * c.LatIce * c.INV_CP,
                           th_atm)

        qc = jnp.where(qc_ge, 0.0, qc)
        nc = jnp.where(qc_ge, 0.0, nc)
        qr = jnp.where(qr_ge, 0.0, qr)
        nr = jnp.where(qr_ge, 0.0, nr)
    else:
        # PORT_NOTES (smoothing, family "homog"): JUMP — s = T_homogfrz -
        # T_atm, scale = 1 K. At T = T_homogfrz the ENTIRE cloud-water
        # and rain categories (finite Qc_nuc = qc, Qr_nuc = qr) convert
        # to ice with latent heat LatIce*q/cp: value-discontinuous in T.
        # The smoothed form freezes the fraction hT of each category.
        # The qc/qr >= QSMALL sub-gates stay hard (jump size QSMALL,
        # effective kinks).
        hT = smoothing.step(c.T_homogfrz - T_atm, smooth_width, scale=1.0)
        fc = jnp.where(qc >= c.QSMALL, hT, 0.0)
        fr = jnp.where(qr >= c.QSMALL, hT, 0.0)

        qm = qm + fc * Qc_nuc + fr * Qr_nuc
        qi = qi + fc * Qc_nuc + fr * Qr_nuc
        bm = bm + (fc * Qc_nuc + fr * Qr_nuc) * c.INV_RHO_RIMEMAX
        ni = ni + fc * Nc_nuc + fr * Nr_nuc
        th_atm = th_atm + inv_exner * (fc * Qc_nuc + fr * Qr_nuc) \
            * c.LatIce * c.INV_CP

        qc = (1.0 - fc) * qc
        nc = (1.0 - fc) * nc
        qr = (1.0 - fr) * qr
        nr = (1.0 - fr) * nr
    return {"qc": qc, "nc": nc, "qr": qr, "nr": nr,
            "qi": qi, "ni": ni, "qm": qm, "bm": bm, "th_atm": th_atm}
