"""p3_main: full P3 microphysics step for a batch of columns.

Source: components/eamxx/src/physics/p3/impl/p3_main_impl.hpp
(p3_main_init + p3_main_internal loop body).

The C++ processes one column per Kokkos team with two early exits:
  - after part1, if neither nucleation is possible nor hydrometeors are
    present, the column keeps its part1 state and the init-value
    diagnostics;
  - after part2, if no hydrometeors remain, the column skips
    sedimentation, homogeneous freezing and part3.
Here every stage runs on all columns and the per-column results are
selected afterwards with the same two masks, which reproduces the C++
control flow exactly (the skipped stages cannot affect the selected
values).

Level convention: k=0 is the model top, kbot = nlev-1 (kdir = -1).
precip_liq_flux / precip_ice_flux live on the nlev+1 interfaces;
precip_ice_flux stays zero, exactly as in the C++ (only the liquid flux
is accumulated by rain sedimentation).
"""

import functools

import jax
import jax.numpy as jnp

from .main_part1 import p3_main_part1
from .main_part2 import p3_main_part2
from .main_part3 import p3_main_part3
from .sedimentation import (cloud_sedimentation, homogeneous_freezing,
                            ice_sedimentation, rain_sedimentation)

# part1 state entries selected column-wise at the first early exit
_STATE_KEYS = ("T_atm", "qv", "th_atm", "qc", "nc", "qr", "nr",
               "qi", "ni", "qm", "bm")


def _sel(col_mask, a, b):
    """Column-wise select: a where col_mask else b (level arrays)."""
    return jnp.where(col_mask[..., None], a, b)


@functools.partial(jax.jit, static_argnames=(
    "predict_nc", "prescribed_ccn", "do_ice_production",
    "use_hetfrz_classnuc", "use_separate_ice_liq_frac"))
def p3_main(dt,
            predict_nc: bool, prescribed_ccn: bool,
            do_ice_production: bool, use_hetfrz_classnuc: bool,
            use_separate_ice_liq_frac: bool,
            # prognostic state
            qc, nc, qr, nr, qi, qm, ni, bm, qv, th_atm,
            # diagnostic inputs
            nc_nuceat_tend, nccn_prescribed, ni_activated, inv_qc_relvar,
            cld_frac_i, cld_frac_l, cld_frac_r, pres, dz, dpres, inv_exner,
            qv_prev, t_prev,
            hetfrz_immersion_nucleation_tend,
            hetfrz_contact_nucleation_tend,
            hetfrz_deposition_nucleation_tend,
            tables, opts):
    """One P3 step. Returns a dict with the updated prognostics and all
    diagnostic outputs (see the return statement)."""
    inv_dt = 1.0 / dt
    inv_exner = jnp.asarray(inv_exner)
    dz = jnp.asarray(dz)
    cld_frac_i = jnp.asarray(cld_frac_i)
    cld_frac_l = jnp.asarray(cld_frac_l)
    cld_frac_r = jnp.asarray(cld_frac_r)

    # ---------------- p3_main_init ----------------
    zcol = jnp.zeros_like(jnp.asarray(qc))
    ze_ice = jnp.full_like(zcol, 1.0e-22)
    ze_rain = jnp.full_like(zcol, 1.0e-22)
    inv_cld_frac_i = 1.0 / cld_frac_i
    inv_cld_frac_l = 1.0 / cld_frac_l
    inv_cld_frac_r = 1.0 / cld_frac_r
    exner = 1.0 / inv_exner
    T_atm = jnp.asarray(th_atm) * exner
    qv = jnp.maximum(jnp.asarray(qv), 0.0)
    inv_dz = 1.0 / dz

    # ---------------- part 1 ----------------
    st = p3_main_part1(
        predict_nc, prescribed_ccn, dt,
        pres, dpres, dz, nc_nuceat_tend, nccn_prescribed,
        inv_exner, exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
        T_atm, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, opts)
    active = st["nucleation_possible"] | st["hydrometeors_present"]

    # ---------------- part 2 ----------------
    g2, diags, hydro2 = p3_main_part2(
        predict_nc, prescribed_ccn, do_ice_production, use_hetfrz_classnuc,
        use_separate_ice_liq_frac, dt, opts["max_total_ni"],
        hetfrz_immersion_nucleation_tend, hetfrz_contact_nucleation_tend,
        hetfrz_deposition_nucleation_tend,
        tables, pres, dpres, dz, nc_nuceat_tend, inv_exner, exner,
        inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
        ni_activated, inv_qc_relvar,
        cld_frac_i, cld_frac_l, cld_frac_r, qv_prev, t_prev, st, opts)

    # first early exit: inactive columns keep their part1 state
    g = {k: _sel(active, g2[k], st[k]) for k in _STATE_KEYS}
    for k in ("qc_incld", "qr_incld", "qi_incld", "qm_incld",
              "nc_incld", "nr_incld", "ni_incld", "bm_incld",
              "mu_c", "lamc", "mu_r", "lamr"):
        base = st[k] if k in st else zcol
        g[k] = _sel(active, g2[k], base)
    diags = {k: _sel(active, v, zcol) for k, v in diags.items()}
    # second exit mask: sed/homog-freezing/part3 run only where
    # hydrometeors remain after part2 (in an active column)
    run3 = active & hydro2

    rho, inv_rho = st["rho"], st["inv_rho"]
    rhofacr, rhofaci, acn = st["rhofacr"], st["rhofaci"], st["acn"]

    # ---------------- sedimentation ----------------
    csed = cloud_sedimentation(
        g["qc_incld"], rho, inv_rho, cld_frac_l, acn, inv_dz,
        dt, inv_dt, predict_nc,
        g["qc"], g["nc"], g["nc_incld"], g["mu_c"], g["lamc"],
        zcol, zcol, jnp.zeros(zcol.shape[:-1]))

    rsed = rain_sedimentation(
        rho, inv_rho, rhofacr, cld_frac_r, inv_dz, g["qr_incld"],
        tables["vn_table_vals"], tables["vm_table_vals"], dt, inv_dt,
        g["qr"], g["nr"], g["nr_incld"], g["mu_r"], g["lamr"],
        jnp.zeros(zcol.shape[:-1] + (zcol.shape[-1] + 1,)), zcol, zcol,
        jnp.zeros(zcol.shape[:-1]), opts)

    ised = ice_sedimentation(
        rho, inv_rho, rhofaci, cld_frac_i, inv_dz, dt, inv_dt,
        g["qi"], g["qi_incld"], g["ni"], g["ni_incld"],
        g["qm"], g["qm_incld"], g["bm"], g["bm_incld"],
        tables["ice_table_vals"], zcol, zcol,
        jnp.zeros(zcol.shape[:-1]), opts)

    sed = dict(g)
    sed["qc"], sed["nc"] = csed["qc"], csed["nc"]
    sed["qr"], sed["nr"] = rsed["qr"], rsed["nr"]
    sed["qi"], sed["ni"] = ised["qi"], ised["ni"]
    sed["qm"], sed["bm"] = ised["qm"], ised["bm"]

    # ---------------- homogeneous freezing ----------------
    if do_ice_production:
        hf = homogeneous_freezing(
            g["T_atm"], inv_exner, sed["qc"], sed["nc"], sed["qr"],
            sed["nr"], sed["qi"], sed["ni"], sed["qm"], sed["bm"],
            g["th_atm"])
        for k in ("qc", "nc", "qr", "nr", "qi", "ni", "qm", "bm"):
            sed[k] = hf[k]
        sed["th_atm"] = hf["th_atm"]

    # ---------------- part 3 ----------------
    p3 = p3_main_part3(
        opts["max_total_ni"], tables["dnu_table_vals"],
        tables["ice_table_vals"], inv_exner,
        cld_frac_l, cld_frac_r, cld_frac_i, rho, inv_rho, rhofaci,
        g["qv"], sed["th_atm"],
        sed["qc"], sed["nc"], sed["qr"], sed["nr"],
        sed["qi"], sed["ni"], sed["qm"], sed["bm"],
        diags["vap_liq_exchange"], ze_rain, ze_ice, opts,
        diag_eff_radius_qc_in=jnp.full_like(zcol, 10.0e-6),
        diag_eff_radius_qr_in=jnp.full_like(zcol, 500.0e-6),
        diag_eff_radius_qi_in=jnp.full_like(zcol, 25.0e-6))

    # second early exit: only run3 columns take the sed/freeze/part3 result
    out_state = {}
    for k in ("qv", "th_atm", "qc", "nc", "qr", "nr", "qi", "ni", "qm", "bm"):
        out_state[k] = _sel(run3, p3[k], g[k])

    vap_liq_exchange = _sel(run3, p3["vap_liq_exchange"],
                            diags["vap_liq_exchange"])

    # diagnostics that keep their init values in skipped columns
    diag_eff_radius_qc = _sel(run3, p3["diag_eff_radius_qc"],
                              jnp.full_like(zcol, 10.0e-6))
    diag_eff_radius_qi = _sel(run3, p3["diag_eff_radius_qi"],
                              jnp.full_like(zcol, 25.0e-6))
    diag_eff_radius_qr = _sel(run3, p3["diag_eff_radius_qr"],
                              jnp.full_like(zcol, 500.0e-6))
    diag_equiv_reflectivity = _sel(run3, p3["diag_equiv_reflectivity"],
                                   jnp.full_like(zcol, -99.0))
    rho_qi = _sel(run3, p3["rho_qi"], zcol)
    diag_vm_qi = _sel(run3, p3["diag_vm_qi"], zcol)
    diag_diam_qi = _sel(run3, p3["diag_diam_qi"], zcol)

    # cloud sed SETS precip_liq_surf, rain sed ADDS to it; both ran with a
    # zero input so the C++ value is the sum of the two contributions
    precip_liq_surf = jnp.where(
        run3, csed["precip_liq_surf"] + rsed["precip_liq_surf"], 0.0)
    precip_ice_surf = jnp.where(run3, ised["precip_ice_surf"], 0.0)
    precip_liq_flux = _sel(run3, rsed["precip_liq_flux"],
                           jnp.zeros(zcol.shape[:-1] + (zcol.shape[-1] + 1,)))
    precip_ice_flux = jnp.zeros_like(precip_liq_flux)

    return {
        **out_state,
        "qv2qi_depos_tend": diags["qv2qi_depos_tend"],
        "precip_total_tend": diags["precip_total_tend"],
        "nevapr": diags["nevapr"],
        "qr_evap_tend": diags["qr_evap_tend"],
        "vap_ice_exchange": diags["vap_ice_exchange"],
        "vap_liq_exchange": vap_liq_exchange,
        "liq_ice_exchange": diags["liq_ice_exchange"],
        "pratot": diags["pratot"],
        "prctot": diags["prctot"],
        "precip_liq_surf": precip_liq_surf,
        "precip_ice_surf": precip_ice_surf,
        "precip_liq_flux": precip_liq_flux,
        "precip_ice_flux": precip_ice_flux,
        "diag_eff_radius_qc": diag_eff_radius_qc,
        "diag_eff_radius_qi": diag_eff_radius_qi,
        "diag_eff_radius_qr": diag_eff_radius_qr,
        "diag_equiv_reflectivity": diag_equiv_reflectivity,
        "rho_qi": rho_qi,
        "diag_vm_qi": diag_vm_qi,
        "diag_diam_qi": diag_diam_qi,
        "qc_sed_tend": jnp.where(run3[..., None], csed["qc_tend"], 0.0),
        "qr_sed_tend": jnp.where(run3[..., None], rsed["qr_tend"], 0.0),
        "qi_sed_tend": jnp.where(run3[..., None], ised["qi_tend"], 0.0),
        # per-column guard against silent CFL-substep truncation in the
        # fixed-length sedimentation scans (see sedimentation.py)
        "sed_converged": (csed["converged"] & rsed["converged"]
                          & ised["converged"]),
    }
