"""p3_main_part1: pre-microphysics state prep and mass clipping.

Source: micro_p3.F90 p3_main_part1 (eam variant). COPIED from
scream_jax/p3/main_part1.py; adaptations:

  * saturation via the EAM qv_sat (MurphyKoop, p - es denominator).
  * EAM naming: `exner` multiplies T (th = T*exner); the theta
    adjustments use exner and the final T uses inv_exner.
  * prescribed-CCN branch: nc = max(nc, nccn_prescribed) — no cloud
    fraction scaling, no SPA factor/exponent.
  * nccnst is a runtime argument (namelist micro_nccons).
"""

import functools

import jax
import jax.numpy as jnp

from . import constants as c
from .conservation import calculate_incloud_mixingratios
from .saturation import qv_sat


@functools.partial(jax.jit,
                   static_argnames=("do_predict_nc", "do_prescribed_ccn"))
def p3_main_part1(do_predict_nc: bool, do_prescribed_ccn: bool, dt, nccnst,
                  pres, dpres, dz, nc_nuceat_tend, nccn_prescribed,
                  exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i,
                  inv_cld_frac_r, t_atm, qv, th,
                  qc, nc, qr, nr, qi, ni, qm, bm):
    """Returns a dict with the updated state (t_atm, qv, th, qc, nc, qr,
    nr, qi, ni, qm, bm), derived fields (rho, inv_rho, qv_sat_l,
    qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn), the in-cloud
    mixing ratios, and per-column flags (nucleation_possible,
    hydrometeors_present)."""
    t_atm = jnp.asarray(t_atm)
    qv = jnp.asarray(qv)
    th = jnp.asarray(th)
    qc, nc, qr, nr = (jnp.asarray(a) for a in (qc, nc, qr, nr))
    qi, ni, qm, bm = (jnp.asarray(a) for a in (qi, ni, qm, bm))
    pres = jnp.asarray(pres)
    exner = jnp.asarray(exner)
    inv_exner = jnp.asarray(inv_exner)

    rho = jnp.asarray(dpres) / jnp.asarray(dz) / c.g
    inv_rho = 1.0 / rho
    qv_sat_l = qv_sat(t_atm, pres, False)
    qv_sat_i = qv_sat(t_atm, pres, True)
    qv_supersat_i = qv / qv_sat_i - 1.0

    rhofacr = (c.rho_1000mb * inv_rho) ** 0.54
    rhofaci = (c.rho_600mb * inv_rho) ** 0.54
    mu = 1.496e-6 * t_atm ** 1.5 / (t_atm + 120.0)
    acn = c.g * c.rho_h2o / (18.0 * mu)

    nucleation_possible = jnp.any(
        (t_atm < c.T_zerodegc) & (qv_supersat_i >= -0.05), axis=-1)

    # --- mass clipping: cloud water ---
    drymass = qc < c.qsmall
    qv = jnp.where(drymass, qv + qc, qv)
    th = jnp.where(drymass, th - exner * qc * c.latvap * c.inv_cp, th)
    qc = jnp.where(drymass, 0.0, qc)
    nc = jnp.where(drymass, 0.0, nc)
    hydro = jnp.any(~drymass, axis=-1)

    # droplet activation / prescription on non-dry levels
    not_dry = ~drymass
    if do_prescribed_ccn:
        nc = jnp.where(not_dry,
                       jnp.maximum(nc, jnp.asarray(nccn_prescribed)), nc)
    elif do_predict_nc:
        nc = jnp.where(not_dry,
                       jnp.maximum(nc + jnp.asarray(nc_nuceat_tend) * dt,
                                   0.0), nc)
    else:
        nc = jnp.where(not_dry, nccnst * inv_rho, nc)

    # --- rain ---
    drymass = qr < c.qsmall
    qv = jnp.where(drymass, qv + qr, qv)
    th = jnp.where(drymass, th - exner * qr * c.latvap * c.inv_cp, th)
    qr = jnp.where(drymass, 0.0, qr)
    nr = jnp.where(drymass, 0.0, nr)
    hydro = hydro | jnp.any(~drymass, axis=-1)

    # --- ice (also sublimate tiny ice in strongly subsaturated air) ---
    drymass = (qi < c.qsmall) | ((qi < 1.0e-8) & (qv_supersat_i < -0.1))
    qv = jnp.where(drymass, qv + qi, qv)
    th = jnp.where(drymass, th - exner * qi * c.latsub * c.inv_cp, th)
    qi = jnp.where(drymass, 0.0, qi)
    ni = jnp.where(drymass, 0.0, ni)
    qm = jnp.where(drymass, 0.0, qm)
    bm = jnp.where(drymass, 0.0, bm)
    hydro = hydro | jnp.any(~drymass, axis=-1)

    # tiny warm ice melts instantly to rain
    drymass = (qi >= c.qsmall) & (qi < 1.0e-8) & (t_atm >= c.T_zerodegc)
    qr = jnp.where(drymass, qr + qi, qr)
    th = jnp.where(drymass, th - exner * qi * c.latice * c.inv_cp, th)
    qi = jnp.where(drymass, 0.0, qi)
    ni = jnp.where(drymass, 0.0, ni)
    qm = jnp.where(drymass, 0.0, qm)
    bm = jnp.where(drymass, 0.0, bm)

    t_atm = th * inv_exner

    ctx = jnp.ones_like(qc, dtype=bool)
    (qc_incld, qr_incld, qi_incld, qm_incld,
     nc_incld, nr_incld, ni_incld, bm_incld) = calculate_incloud_mixingratios(
        qc, qr, qi, qm, nc, nr, ni, bm,
        inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, ctx)

    return {
        "t_atm": t_atm, "qv": qv, "th": th,
        "qc": qc, "nc": nc, "qr": qr, "nr": nr,
        "qi": qi, "ni": ni, "qm": qm, "bm": bm,
        "rho": rho, "inv_rho": inv_rho,
        "qv_sat_l": qv_sat_l, "qv_sat_i": qv_sat_i,
        "qv_supersat_i": qv_supersat_i,
        "rhofacr": rhofacr, "rhofaci": rhofaci, "acn": acn,
        "qc_incld": qc_incld, "qr_incld": qr_incld, "qi_incld": qi_incld,
        "qm_incld": qm_incld, "nc_incld": nc_incld, "nr_incld": nr_incld,
        "ni_incld": ni_incld, "bm_incld": bm_incld,
        "nucleation_possible": nucleation_possible,
        "hydrometeors_present": hydro,
    }
