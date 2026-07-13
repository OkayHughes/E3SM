"""p3_main: full EAM P3 microphysics step for a batch of columns.

Source: micro_p3.F90 p3_main (eam variant). COPIED from
scream_jax/p3/main.py; adaptations (see package docstring): EAM
column sequence part1 -> part2 -> sedimentation (with cflx/rflx/sflx)
-> homogeneous_freezing -> ice_complete_melting -> mincdnc floor ->
part3; p3_tend_out reconstruction; EAM naming and outputs.

The Fortran per-column early exits become column masks:
  - after part1, columns with neither nucleation possible nor
    hydrometeors present keep their part1 state (goto 333);
  - after part2, columns without hydrometeors skip sedimentation,
    homogeneous freezing, complete melting and part3.

Level convention: k=0 is the model top, kbot = nlev-1 (kdir = -1).
precip fluxes live on the nlev+1 interfaces; precip_ice_flux stays
zero exactly as in the Fortran (never accumulated).
"""

import functools

import jax
import jax.numpy as jnp

from . import constants as c
from .main_part1 import p3_main_part1
from .main_part2 import p3_main_part2
from .main_part3 import p3_main_part3
from .sedimentation import (
    cloud_sedimentation,
    homogeneous_freezing,
    ice_complete_melting,
    ice_sedimentation,
    rain_sedimentation,
)

_STATE_KEYS = ("t_atm", "qv", "th", "qc", "nc", "qr", "nr",
               "qi", "ni", "qm", "bm")


def _sel(col_mask, a, b):
    """Column-wise select: a where col_mask else b (level arrays)."""
    return jnp.where(col_mask[..., None], a, b)


@functools.partial(jax.jit, static_argnames=(
    "do_predict_nc", "do_prescribed_ccn", "do_precip_off",
    "use_hetfrz_classnuc", "do_cooper"))
def p3_main(dt,
            do_predict_nc: bool, do_prescribed_ccn: bool,
            do_precip_off: bool, use_hetfrz_classnuc: bool,
            do_cooper: bool,
            # prognostic state
            qc, nc, qr, nr, th, qv, qi, qm, ni, bm,
            # diagnostic inputs
            pres, dz, nc_nuceat_tend, nccn_prescribed, ni_activated,
            frzimm, frzcnt, frzdep, inv_qc_relvar, dpres, exner,
            cld_frac_r, cld_frac_l, cld_frac_i, qv_prev, t_prev,
            tables, opts):
    """One P3 step. Returns a dict with the updated prognostics and all
    diagnostic outputs (state fields named as in the golden driver)."""
    inv_dt = 1.0 / dt
    exner = jnp.asarray(exner)
    dz = jnp.asarray(dz)
    cld_frac_i = jnp.asarray(cld_frac_i)
    cld_frac_l = jnp.asarray(cld_frac_l)
    cld_frac_r = jnp.asarray(cld_frac_r)

    # ---------------- init ----------------
    qc = jnp.asarray(qc)
    zcol = jnp.zeros_like(qc)
    zflx = jnp.zeros(qc.shape[:-1] + (qc.shape[-1] + 1,))
    ze_ice = jnp.full_like(zcol, 1.0e-22)
    ze_rain = jnp.full_like(zcol, 1.0e-22)
    inv_cld_frac_i = 1.0 / cld_frac_i
    inv_cld_frac_l = 1.0 / cld_frac_l
    inv_cld_frac_r = 1.0 / cld_frac_r
    inv_exner = 1.0 / exner
    qv = jnp.maximum(jnp.asarray(qv), 0.0)
    inv_dz = 1.0 / dz

    qc_old, nc_old = qc, jnp.asarray(nc)
    qr_old, nr_old = jnp.asarray(qr), jnp.asarray(nr)
    qi_old, ni_old = jnp.asarray(qi), jnp.asarray(ni)
    qv_old, th_old = qv, jnp.asarray(th)

    # ---------------- part 1 ----------------
    st = p3_main_part1(
        do_predict_nc, do_prescribed_ccn, dt, opts["nccnst"],
        pres, dpres, dz, nc_nuceat_tend, nccn_prescribed,
        exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
        th_old * inv_exner, qv, th, qc, nc, qr, nr, qi, ni, qm, bm)
    active = st["nucleation_possible"] | st["hydrometeors_present"]

    # ---------------- part 2 ----------------
    g2, diags, tend2, hydro2 = p3_main_part2(
        do_predict_nc, do_prescribed_ccn, use_hetfrz_classnuc, do_cooper,
        do_precip_off, dt, tables, pres, dpres, dz, nc_nuceat_tend,
        exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
        ni_activated, inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r,
        qv_prev, t_prev, frzimm, frzcnt, frzdep, st, opts)

    # first early exit: inactive columns keep their part1 state
    g = {k: _sel(active, g2[k], st[k]) for k in _STATE_KEYS}
    for k in ("qc_incld", "qr_incld", "qi_incld", "qm_incld",
              "nc_incld", "nr_incld", "ni_incld", "bm_incld"):
        g[k] = _sel(active, g2[k], st[k])
    for k in ("mu_c", "nu", "lamc", "cdist", "cdist1", "cdistr",
              "mu_r", "lamr", "logn0r"):
        g[k] = _sel(active, g2[k], zcol)
    diags = {k: _sel(active, v, zcol) for k, v in diags.items()}
    tend2 = {k: _sel(active, v, zcol) for k, v in tend2.items()}

    # measured process tendencies (slots 42-49), for active columns
    tend2[42] = _sel(active, (g["qc"] - qc_old) * inv_dt, zcol)
    tend2[43] = _sel(active, (g["nc"] - nc_old) * inv_dt, zcol)
    tend2[44] = _sel(active, (g["qr"] - qr_old) * inv_dt, zcol)
    tend2[45] = _sel(active, (g["nr"] - nr_old) * inv_dt, zcol)
    tend2[46] = _sel(active, (g["qi"] - qi_old) * inv_dt, zcol)
    tend2[47] = _sel(active, (g["ni"] - ni_old) * inv_dt, zcol)
    tend2[48] = _sel(active, (g["qv"] - qv_old) * inv_dt, zcol)
    tend2[49] = _sel(active, (g["th"] - th_old) * inv_dt, zcol)

    # second exit mask: the post-part2 sequence runs only where
    # hydrometeors remain in an active column
    run3 = active & hydro2

    rho, inv_rho = st["rho"], st["inv_rho"]
    rhofacr, rhofaci, acn = st["rhofacr"], st["rhofaci"], st["acn"]

    # ---------------- sedimentation ----------------
    csed = cloud_sedimentation(
        g["qc_incld"], rho, inv_rho, cld_frac_l, acn, inv_dz,
        dt, inv_dt, do_predict_nc,
        g["qc"], g["nc"], g["nc_incld"], g["mu_c"], g["lamc"],
        jnp.zeros(zcol.shape[:-1]), zflx, g["qc"], g["nc"])

    rsed = rain_sedimentation(
        rho, inv_rho, rhofacr, cld_frac_r, inv_dz, g["qr_incld"],
        tables["vn_table_vals"], tables["vm_table_vals"], dt, inv_dt,
        opts["p3_max_mean_rain_size"],
        g["qr"], g["nr"], g["nr_incld"], g["mu_r"], g["lamr"],
        csed["precip_liq_surf"], zflx, zflx, g["qr"], g["nr"])

    ised = ice_sedimentation(
        rho, inv_rho, rhofaci, cld_frac_i, inv_dz, dt, inv_dt,
        g["qi"], g["qi_incld"], g["ni"], g["ni_incld"],
        g["qm"], g["qm_incld"], g["bm"], g["bm_incld"],
        tables["ice_table_vals"], jnp.zeros(zcol.shape[:-1]), zflx,
        g["qi"], g["ni"])

    # sedimentation tendency slots 36-41 (init value = pre-sed state)
    tend2[36] = _sel(run3, csed["qc_tend"], zcol)
    tend2[37] = _sel(run3, csed["nc_tend"], zcol)
    tend2[38] = _sel(run3, rsed["qr_tend"], zcol)
    tend2[39] = _sel(run3, rsed["nr_tend"], zcol)
    tend2[40] = _sel(run3, ised["qi_tend"], zcol)
    tend2[41] = _sel(run3, ised["ni_tend"], zcol)

    sed = dict(g)
    sed["qc"], sed["nc"] = csed["qc"], csed["nc"]
    sed["mu_c"], sed["lamc"] = csed["mu_c"], csed["lamc"]
    sed["qr"], sed["nr"] = rsed["qr"], rsed["nr"]
    sed["mu_r"], sed["lamr"] = rsed["mu_r"], rsed["lamr"]
    sed["qi"], sed["ni"] = ised["qi"], ised["ni"]
    sed["qm"], sed["bm"] = ised["qm"], ised["bm"]
    sed["nc_incld"] = csed["nc_incld"]

    # ---------------- homogeneous freezing ----------------
    hf = homogeneous_freezing(
        exner, c.latice, sed["qc"], sed["nc"], sed["qr"], sed["nr"],
        sed["qi"], sed["ni"], sed["qm"], sed["bm"], g["th"])
    for k in ("qc", "nc", "qr", "nr", "qi", "ni", "qm", "bm", "th"):
        sed[k] = hf[k]

    # ---------------- instantaneous complete melting ----------------
    cm = ice_complete_melting(
        exner, c.latice, sed["qi"], sed["ni"], sed["qm"],
        sed["qr"], sed["nr"], sed["qc"], sed["nc"], sed["th"])
    for k in ("qi", "ni", "qm", "qr", "nr", "qc", "nc", "th"):
        sed[k] = cm[k]

    # ---------------- p3_mincdnc floor ----------------
    mincdnc = jnp.asarray(opts["p3_mincdnc"])
    floor = (mincdnc > 0.0) & (sed["qc"] >= c.qsmall)
    sed["nc"] = jnp.where(
        floor, jnp.maximum(sed["nc"], mincdnc * cld_frac_l / rho),
        sed["nc"])
    sed["nc_incld"] = jnp.where(
        floor, jnp.maximum(sed["nc_incld"], mincdnc / rho),
        sed["nc_incld"])

    # ---------------- part 3 ----------------
    p3 = p3_main_part3(
        opts["p3_max_mean_rain_size"], mincdnc,
        exner, cld_frac_l, cld_frac_r, cld_frac_i, rho, inv_rho, rhofaci,
        g["qv"], sed["th"],
        sed["qc"], sed["nc"], sed["qr"], sed["nr"],
        sed["qi"], sed["ni"], sed["qm"], sed["bm"],
        sed["mu_c"], sed["lamc"], sed["mu_r"], sed["lamr"],
        diags["vap_liq_exchange"], ze_rain, ze_ice,
        tables["ice_table_vals"])

    # second early exit: only run3 columns take the post-part2 results
    out_state = {}
    for k in ("qv", "th", "qc", "nc", "qr", "nr", "qi", "ni", "qm", "bm"):
        out_state[k] = _sel(run3, p3[k], g[k])

    vap_liq_exchange = _sel(run3, p3["vap_liq_exchange"],
                            diags["vap_liq_exchange"])
    mu_c_out = _sel(run3, p3["mu_c"], g["mu_c"])
    lamc_out = _sel(run3, p3["lamc"], g["lamc"])

    diag_eff_radius_qc = _sel(run3, p3["diag_eff_radius_qc"],
                              jnp.full_like(zcol, 10.0e-6))
    diag_eff_radius_qi = _sel(run3, p3["diag_eff_radius_qi"],
                              jnp.full_like(zcol, 25.0e-6))
    diag_equiv_reflectivity = _sel(run3, p3["diag_equiv_reflectivity"],
                                   jnp.full_like(zcol, -99.0))
    diag_ze_rain = _sel(run3, p3["diag_ze_rain"],
                        jnp.full_like(zcol, -99.0))
    diag_ze_ice = _sel(run3, p3["diag_ze_ice"], jnp.full_like(zcol, -99.0))
    rho_qi = _sel(run3, p3["rho_qi"], zcol)

    precip_liq_surf = jnp.where(run3, rsed["precip_liq_surf"], 0.0)
    precip_ice_surf = jnp.where(run3, ised["precip_ice_surf"], 0.0)
    precip_liq_flux = _sel(run3, rsed["precip_liq_flux"], zflx)
    rflx = _sel(run3, rsed["rflx"], zflx)
    sflx = _sel(run3, ised["sflx"], zflx)
    cflx = _sel(run3, csed["cflx"], zflx)
    precip_ice_flux = zflx  # never accumulated in the Fortran

    # p3_tend_out: assemble (…, nlev, 49); slots not recorded stay 0
    nlev = zcol.shape[-1]
    tend_out = jnp.zeros(zcol.shape[:-1] + (nlev, 49))
    for slot, arr in tend2.items():
        tend_out = tend_out.at[..., slot - 1].set(arr)

    return {
        **out_state,
        "mu_c": mu_c_out, "lamc": lamc_out,
        "qv2qi_depos_tend": diags["qv2qi_depos_tend"],
        "precip_total_tend": diags["precip_total_tend"],
        "nevapr": diags["nevapr"],
        "qr_evap_tend": diags["qr_evap_tend"],
        "vap_ice_exchange": diags["vap_ice_exchange"],
        "vap_liq_exchange": vap_liq_exchange,
        "liq_ice_exchange": diags["liq_ice_exchange"],
        "precip_liq_surf": precip_liq_surf,
        "precip_ice_surf": precip_ice_surf,
        "precip_liq_flux": precip_liq_flux,
        "precip_ice_flux": precip_ice_flux,
        "rflx": rflx, "sflx": sflx, "cflx": cflx,
        "diag_eff_radius_qc": diag_eff_radius_qc,
        "diag_eff_radius_qi": diag_eff_radius_qi,
        "diag_equiv_reflectivity": diag_equiv_reflectivity,
        "diag_ze_rain": diag_ze_rain, "diag_ze_ice": diag_ze_ice,
        "rho_qi": rho_qi,
        "p3_tend_out": tend_out,
    }
