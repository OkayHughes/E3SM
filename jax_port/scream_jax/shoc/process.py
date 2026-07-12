"""EAMxx SHOC process step: pre-processing, shoc_main, post-processing.

Source: components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.hpp
(SHOCPreprocess/SHOCPostprocess functors) and .cpp (run_impl: z_surf = 0,
hdtime = dt so nadv = 1 per AD subcycle, wtracer_sfc = 0, dx = dy =
cell_length).

This is the unit the EAMxx swap test exchanges and the unit the golden
archives capture: EAMxx fields in, EAMxx fields out.
"""

import functools

import jax
import jax.numpy as jnp

from ..foundation import constants as c
from ..foundation.thermo import (
    calculate_dse,
    calculate_dz,
    calculate_temperature_from_dse,
    calculate_theta_from_T,
    calculate_thetal_from_theta,
    calculate_z_int,
    calculate_z_mid,
    exner_function,
)
from . import constants as sc
from .interp import linear_interp
from .main import shoc_main


@functools.partial(jax.jit, static_argnames=("npbl", "shoc_1p5tke", "extra_diags"))
def shoc_process_step(dt, npbl: int, cell_length,
                      # Runtime parameters (input.yaml shoc block)
                      lambda_low, lambda_high, lambda_slope, lambda_thresh,
                      thl2tune, qw2tune, qwthl2tune, w2tune, length_fac,
                      c_diag_3rd_mom, ckh, ckm,
                      shoc_1p5tke: bool, extra_diags: bool,
                      # EAMxx fields (one AD subcycle's inputs)
                      T_mid, p_mid, p_int, pseudo_density, omega, phis,
                      surf_sens_flux, surf_evap, surf_mom_flux_x,
                      surf_mom_flux_y, qv, qc, tke, u_wind, v_wind,
                      cldfrac_liq, sgs_buoy_flux, eddy_diff_mom):
    """One SHOC process step (= one AD subcycle; nadv = 1).

    Returns a dict of updated/computed EAMxx fields keyed by their EAMxx
    field names.
    """
    T_mid = jnp.asarray(T_mid)
    p_mid = jnp.asarray(p_mid)
    qv = jnp.asarray(qv)
    qc = jnp.asarray(qc)
    tke = jnp.asarray(tke)
    cldfrac_liq = jnp.asarray(cldfrac_liq)
    phis = jnp.asarray(phis)

    # ---------------- SHOCPreprocess ----------------
    cldfrac_liq_prev = cldfrac_liq
    inv_exner = 1.0 / exner_function(p_mid)
    tke = jnp.maximum(sc.mintke, tke)
    qw = qv + qc

    theta = calculate_theta_from_T(T_mid, p_mid)
    thlm = calculate_thetal_from_theta(theta, T_mid, qc)
    thv = theta * (1.0 + c.ZVIR * qv - qc)

    dz = calculate_dz(pseudo_density, p_mid, T_mid, qv)
    rrho = (1.0 / c.gravit) * (jnp.asarray(pseudo_density) / dz)
    wm_zt = -jnp.asarray(omega) / (rrho * c.gravit)

    z_int = calculate_z_int(dz, 0.0)          # z_surf = 0
    z_mid = calculate_z_mid(z_int)
    zt_grid = z_mid - z_int[..., -1:]
    zi_grid = (z_int - z_int[..., -1:]).at[..., -1].set(0.0)

    shoc_s = calculate_dse(T_mid, z_mid, 0.0) + phis[..., None]

    nlev = T_mid.shape[-1]
    rrho_i = linear_interp(zt_grid, zi_grid, rrho, nlev, nlev + 1, 0.0)

    exner_int_sfc = exner_function(jnp.asarray(p_int)[..., -1])
    rrho_sfc = rrho_i[..., -1]
    wpthlp_sfc = (jnp.asarray(surf_sens_flux) / (c.Cpair * rrho_sfc)) / exner_int_sfc
    wprtp_sfc = jnp.asarray(surf_evap) / rrho_sfc
    upwp_sfc = jnp.asarray(surf_mom_flux_x) / rrho_sfc
    vpwp_sfc = jnp.asarray(surf_mom_flux_y) / rrho_sfc

    # Tracer group {qv, qc, tke}: tke enters clipped (the field was updated
    # in place before the copy in the C++); qc/qv as-is. Group-diffused
    # qv/qc/tke are discarded in postprocessing in favor of qw/shoc_ql and
    # the dedicated copies, exactly as in the C++ *_copy dance.
    qtracers = jnp.stack([qv, qc, tke], axis=-2)
    wtracer_sfc = jnp.zeros(qtracers.shape[:-1], dtype=T_mid.dtype)

    # ---------------- shoc_main (nadv = 1) ----------------
    out = shoc_main(
        dt, 1, npbl,
        lambda_low, lambda_high, lambda_slope, lambda_thresh,
        thl2tune, qw2tune, qwthl2tune, w2tune, length_fac,
        c_diag_3rd_mom, ckh, ckm, shoc_1p5tke, extra_diags,
        cell_length, cell_length, zt_grid, zi_grid, p_mid, p_int,
        pseudo_density, thv, wm_zt,
        wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
        wtracer_sfc, inv_exner, phis,
        shoc_s, tke, thlm, qw, jnp.asarray(u_wind), jnp.asarray(v_wind),
        jnp.asarray(sgs_buoy_flux), qtracers, jnp.asarray(eddy_diff_mom),
        cldfrac_liq, qc)

    # ---------------- SHOCPostprocess ----------------
    tke_new = out["tke"]
    qc_new = out["shoc_ql"]
    qv_new = out["qw"] - qc_new
    cldfrac_new = jnp.minimum(out["shoc_cldfrac"], 1.0)

    qc2 = out["shoc_ql2"]
    cond = (qc_new != 0.0) & (qc2 != 0.0)
    qc2_safe = jnp.where(cond, qc2, 1.0)
    inv_qc_relvar = jnp.where(
        cond,
        jnp.minimum(10.0, jnp.maximum(0.001, qc_new ** 2 / qc2_safe)),
        1.0)

    T_new = calculate_temperature_from_dse(
        out["host_dse"] - phis[..., None], z_mid, 0.0)

    return {
        "T_mid": T_new, "qv": qv_new, "qc": qc_new, "tke": tke_new,
        "u_wind": out["u_wind"], "v_wind": out["v_wind"],
        "cldfrac_liq": cldfrac_new, "cldfrac_liq_prev": cldfrac_liq_prev,
        "sgs_buoy_flux": out["wthv_sec"], "eddy_diff_mom": out["tk"],
        "eddy_diff_heat": out["tkh"], "inv_qc_relvar": inv_qc_relvar,
        "pbl_height": out["pblh"], "ustar": out["ustar"],
        "obklen": out["obklen"], "w_variance": out["w_sec"],
        "thl_sec": out["thl_sec"],
    }
