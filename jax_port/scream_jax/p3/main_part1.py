"""p3_main_part1: pre-microphysics state prep and mass clipping.

Source: components/eamxx/src/physics/p3/impl/p3_main_impl_part1.hpp

The C++ per-column early-exit booleans (nucleationPossible,
hydrometeorsPresent) become per-column boolean arrays; the driver uses
them to select whether part2+'s results apply (columns that would have
early-returned keep their part1 state, exactly as in the C++).
"""

import functools

import jax
import jax.numpy as jnp

from ..foundation import constants as c
from ..foundation.saturation import SaturationFcn, qv_sat_dry
from .conservation import calculate_incloud_mixingratios


@functools.partial(jax.jit, static_argnames=("predict_nc", "do_prescribed_ccn"))
def p3_main_part1(predict_nc: bool, do_prescribed_ccn: bool, dt,
                  pres, dpres, dz, nc_nuceat_tend, nccn_prescribed,
                  inv_exner, exner, inv_cld_frac_l, inv_cld_frac_i,
                  inv_cld_frac_r, T_atm, qv, th_atm,
                  qc, nc, qr, nr, qi, ni, qm, bm, opts):
    """Returns a dict with the updated state (T_atm, qv, th_atm, qc, nc,
    qr, nr, qi, ni, qm, bm), derived fields (rho, inv_rho, qv_sat_l,
    qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn), the in-cloud mixing
    ratios, and per-column flags (nucleation_possible,
    hydrometeors_present)."""
    T_atm = jnp.asarray(T_atm)
    qv = jnp.asarray(qv)
    th_atm = jnp.asarray(th_atm)
    qc, nc, qr, nr = (jnp.asarray(a) for a in (qc, nc, qr, nr))
    qi, ni, qm, bm = (jnp.asarray(a) for a in (qi, ni, qm, bm))
    pres = jnp.asarray(pres)
    inv_exner = jnp.asarray(inv_exner)

    rho = jnp.asarray(dpres) / jnp.asarray(dz) / c.gravit
    inv_rho = 1.0 / rho
    qv_sat_l = qv_sat_dry(T_atm, pres, False, SaturationFcn.MURPHY_KOOP)
    qv_sat_i = qv_sat_dry(T_atm, pres, True, SaturationFcn.MURPHY_KOOP)
    qv_supersat_i = qv / qv_sat_i - 1.0

    rhofacr = (c.RHO_1000MB * inv_rho) ** 0.54
    rhofaci = (c.RHO_600MB * inv_rho) ** 0.54
    mu = 1.496e-6 * T_atm ** 1.5 / (T_atm + 120.0)
    acn = c.gravit * c.RHO_H2O / (18.0 * mu)

    nucleation_possible = jnp.any(
        (T_atm < c.T_zerodegc) & (qv_supersat_i >= -0.05), axis=-1)

    # --- mass clipping: cloud water ---
    drymass = qc < c.QSMALL
    qv = jnp.where(drymass, qv + qc, qv)
    th_atm = jnp.where(drymass, th_atm - inv_exner * qc * c.LatVap * c.INV_CP,
                       th_atm)
    qc = jnp.where(drymass, 0.0, qc)
    nc = jnp.where(drymass, 0.0, nc)
    hydro = jnp.any(~drymass, axis=-1)

    # droplet activation / prescription on non-dry levels
    not_dry = ~drymass
    if do_prescribed_ccn:
        nccn_scaled = jnp.asarray(nccn_prescribed) / jnp.asarray(inv_cld_frac_l)
        nccn_scaled = nccn_scaled ** opts["spa_ccn_to_nc_exponent"]
        nc = jnp.where(not_dry,
                       jnp.maximum(nc, opts["spa_ccn_to_nc_factor"] * nccn_scaled),
                       nc)
    elif predict_nc:
        nc = jnp.where(not_dry,
                       jnp.maximum(nc + jnp.asarray(nc_nuceat_tend) * dt, 0.0), nc)
    else:
        nc = jnp.where(not_dry, c.NCCNST * inv_rho, nc)

    # --- rain ---
    drymass = qr < c.QSMALL
    qv = jnp.where(drymass, qv + qr, qv)
    th_atm = jnp.where(drymass, th_atm - inv_exner * qr * c.LatVap * c.INV_CP,
                       th_atm)
    qr = jnp.where(drymass, 0.0, qr)
    nr = jnp.where(drymass, 0.0, nr)
    hydro = hydro | jnp.any(~drymass, axis=-1)

    # --- ice (also sublimate tiny ice in strongly subsaturated air) ---
    drymass = (qi < c.QSMALL) | ((qi < 1.0e-8) & (qv_supersat_i < -0.1))
    qv = jnp.where(drymass, qv + qi, qv)
    th_atm = jnp.where(drymass,
                       th_atm - inv_exner * qi * (c.LatVap + c.LatIce) * c.INV_CP,
                       th_atm)
    qi = jnp.where(drymass, 0.0, qi)
    ni = jnp.where(drymass, 0.0, ni)
    qm = jnp.where(drymass, 0.0, qm)
    bm = jnp.where(drymass, 0.0, bm)
    hydro = hydro | jnp.any(~drymass, axis=-1)

    # tiny warm ice melts instantly to rain
    drymass = (qi >= c.QSMALL) & (qi < 1.0e-8) & (T_atm >= c.T_zerodegc)
    qr = jnp.where(drymass, qr + qi, qr)
    th_atm = jnp.where(drymass, th_atm - inv_exner * qi * c.LatIce * c.INV_CP,
                       th_atm)
    qi = jnp.where(drymass, 0.0, qi)
    ni = jnp.where(drymass, 0.0, ni)
    qm = jnp.where(drymass, 0.0, qm)
    bm = jnp.where(drymass, 0.0, bm)

    T_atm = th_atm * jnp.asarray(exner)

    ctx = jnp.ones_like(qc, dtype=bool)
    (qc_incld, qr_incld, qi_incld, qm_incld,
     nc_incld, nr_incld, ni_incld, bm_incld) = calculate_incloud_mixingratios(
        qc, qr, qi, qm, nc, nr, ni, bm,
        inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, ctx)

    return {
        "T_atm": T_atm, "qv": qv, "th_atm": th_atm,
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
