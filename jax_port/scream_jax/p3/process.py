"""EAMxx P3 process step: pre-processing, p3_main, post-processing.

Sources: components/eamxx/src/physics/p3/eamxx_p3_process_interface.hpp
(p3_preamble / p3_postamble functors), eamxx_p3_process_interface.cpp
(field wiring: p3_main's pres/dpres are the DRY pressure/thickness) and
eamxx_p3_run.cpp (run_impl).

This is the unit the EAMxx swap test exchanges and the unit the golden
archives capture: EAMxx fields in (wet mixing ratios), EAMxx fields out.
"""

import functools

import jax
import jax.numpy as jnp

from ..foundation import constants as c
from ..foundation.thermo import (
    calculate_drymmr_from_wetmmr_dp_based,
    calculate_dz,
    calculate_T_from_theta,
    calculate_theta_from_T,
    calculate_wetmmr_from_drymmr_dp_based,
    exner_function,
)
from .main import p3_main

MINCLD = 0.0001

_WATER_KEYS = ("qc", "nc", "qr", "nr", "qi", "ni", "qm", "bm", "qv")


@functools.partial(jax.jit, static_argnames=(
    "predict_nc", "prescribed_ccn", "do_ice_production",
    "use_hetfrz_classnuc", "use_separate_ice_liq_frac",
    "set_cld_frac_l_to_one", "set_cld_frac_i_to_one",
    "set_cld_frac_r_to_one", "sed_use_while_loop"))
def p3_process_step(dt,
                    predict_nc: bool, prescribed_ccn: bool,
                    do_ice_production: bool, use_hetfrz_classnuc: bool,
                    use_separate_ice_liq_frac: bool,
                    set_cld_frac_l_to_one: bool, set_cld_frac_i_to_one: bool,
                    set_cld_frac_r_to_one: bool,
                    # EAMxx input fields
                    T_mid, p_mid, p_dry_mid, pseudo_density,
                    pseudo_density_dry, cldfrac_tot,
                    qv, qc, nc, qr, nr, qi, qm, ni, bm,
                    qv_prev_micro_step, T_prev_micro_step,
                    nc_nuceat_tend, nccn, ni_activated, inv_qc_relvar,
                    precip_liq_surf_mass, precip_ice_surf_mass,
                    tables, opts,
                    cldfrac_liq_in=None, cldfrac_ice_in=None,
                    hetfrz_immersion_nucleation_tend=None,
                    hetfrz_contact_nucleation_tend=None,
                    hetfrz_deposition_nucleation_tend=None,
                    sed_use_while_loop=True):
    """One P3 process step (one AD subcycle). Returns a dict of updated /
    computed EAMxx fields keyed by their EAMxx field names."""
    T_mid = jnp.asarray(T_mid)
    p_mid = jnp.asarray(p_mid)
    pseudo_density = jnp.asarray(pseudo_density)
    pseudo_density_dry = jnp.asarray(pseudo_density_dry)
    cldfrac_tot = jnp.asarray(cldfrac_tot)
    qv = jnp.asarray(qv)

    # ---------------- p3_preamble ----------------
    # dz from FULL pressure/density with the still-wet qv
    dz = calculate_dz(pseudo_density, p_mid, T_mid, qv)

    wet = dict(qc=qc, nc=nc, qr=qr, nr=nr, qi=qi, ni=ni, qm=qm, bm=bm, qv=qv)
    dry = {k: calculate_drymmr_from_wetmmr_dp_based(
        jnp.asarray(v), pseudo_density, pseudo_density_dry)
        for k, v in wet.items()}
    qv_prev_dry = calculate_drymmr_from_wetmmr_dp_based(
        jnp.asarray(qv_prev_micro_step), pseudo_density, pseudo_density_dry)

    inv_exner = 1.0 / exner_function(p_mid)
    th_atm = calculate_theta_from_T(T_mid, p_mid)

    if use_separate_ice_liq_frac:
        cld_frac_l = (jnp.ones_like(cldfrac_tot) if set_cld_frac_l_to_one
                      else jnp.maximum(jnp.asarray(cldfrac_liq_in), MINCLD))
        cld_frac_i = (jnp.ones_like(cldfrac_tot) if set_cld_frac_i_to_one
                      else jnp.maximum(jnp.asarray(cldfrac_ice_in), MINCLD))
    else:
        cld_frac_l = (jnp.ones_like(cldfrac_tot) if set_cld_frac_l_to_one
                      else jnp.maximum(cldfrac_tot, MINCLD))
        cld_frac_i = (jnp.ones_like(cldfrac_tot) if set_cld_frac_i_to_one
                      else jnp.maximum(cldfrac_tot, MINCLD))
    if set_cld_frac_r_to_one:
        cld_frac_r = jnp.ones_like(cldfrac_tot)
    else:
        # max-overlap: rain fraction at k inherits the total cloud
        # fraction of the level above where that is larger
        cld_frac_r = jnp.maximum(cldfrac_tot, MINCLD)
        cld_frac_r = cld_frac_r.at[..., 1:].set(
            jnp.maximum(cld_frac_r[..., 1:], cldfrac_tot[..., :-1]))

    zcol = jnp.zeros_like(T_mid)
    hetfrz = [jnp.asarray(a) if a is not None else zcol
              for a in (hetfrz_immersion_nucleation_tend,
                        hetfrz_contact_nucleation_tend,
                        hetfrz_deposition_nucleation_tend)]
    nccn_arr = jnp.asarray(nccn) if nccn is not None else zcol

    # ---------------- p3_main ----------------
    out = p3_main(
        dt, predict_nc, prescribed_ccn, do_ice_production,
        use_hetfrz_classnuc, use_separate_ice_liq_frac,
        dry["qc"], dry["nc"], dry["qr"], dry["nr"], dry["qi"], dry["qm"],
        dry["ni"], dry["bm"], dry["qv"], th_atm,
        nc_nuceat_tend, nccn_arr, ni_activated, inv_qc_relvar,
        cld_frac_i, cld_frac_l, cld_frac_r,
        p_dry_mid, dz, pseudo_density_dry, inv_exner,
        qv_prev_dry, T_prev_micro_step,
        hetfrz[0], hetfrz[1], hetfrz[2], tables, opts,
        sed_use_while_loop=sed_use_while_loop)

    # ---------------- p3_postamble ----------------
    # rescaled temperature update: T += (T(th_new) - T_before)*dp_dry/dp
    T_new = (T_mid + (calculate_T_from_theta(out["th_atm"], p_mid) - T_mid)
             * pseudo_density_dry / pseudo_density)

    wet_out = {k: calculate_wetmmr_from_drymmr_dp_based(
        out[k], pseudo_density, pseudo_density_dry) for k in _WATER_KEYS}

    precip_liq_mass = (jnp.asarray(precip_liq_surf_mass)
                       + out["precip_liq_surf"] * c.RHO_H2O * dt)
    precip_ice_mass = (jnp.asarray(precip_ice_surf_mass)
                       + out["precip_ice_surf"] * c.RHO_H2O * dt)

    return {
        "T_mid": T_new,
        "T_prev_micro_step": T_new,
        "qv": wet_out["qv"],
        "qc": wet_out["qc"], "nc": wet_out["nc"],
        "qr": wet_out["qr"], "nr": wet_out["nr"],
        "qi": wet_out["qi"], "ni": wet_out["ni"],
        "qm": wet_out["qm"], "bm": wet_out["bm"],
        "qv_prev_micro_step": wet_out["qv"],
        "eff_radius_qc": out["diag_eff_radius_qc"] * 1e6,
        "eff_radius_qi": out["diag_eff_radius_qi"] * 1e6,
        "eff_radius_qr": out["diag_eff_radius_qr"] * 1e6,
        "precip_liq_surf_mass": precip_liq_mass,
        "precip_ice_surf_mass": precip_ice_mass,
        "precip_total_tend": out["precip_total_tend"],
        "nevapr": out["nevapr"],
        "diag_equiv_reflectivity": out["diag_equiv_reflectivity"],
        "micro_liq_ice_exchange": out["liq_ice_exchange"],
        "micro_vap_liq_exchange": out["vap_liq_exchange"],
        "micro_vap_ice_exchange": out["vap_ice_exchange"],
        "rainfrac": cld_frac_r,
        # boundary fluxes for mass/energy conservation checks
        "vapor_flux": jnp.zeros_like(precip_liq_mass),
        "water_flux": out["precip_liq_surf"] + out["precip_ice_surf"],
        "ice_flux": out["precip_ice_surf"],
        "heat_flux": jnp.zeros_like(precip_liq_mass),
    }
