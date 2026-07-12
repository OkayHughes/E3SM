"""SHOC main driver.

Source: components/eamxx/src/physics/shoc/impl/shoc_main_impl.hpp
(shoc_init for npbl, shoc_main_internal for the time loop; the
small-kernels _disp variant is a GPU dispatch detail).

State is threaded functionally: prognostics go in as arrays and come back
updated in a dict, in the exact operation order of the C++.
"""

import functools

import numpy as np

import jax
import jax.numpy as jnp

from . import constants as sc
from .assumed_pdf import shoc_assumed_pdf
from .energy import shoc_energy_fixer, shoc_energy_integrals, update_host_dse
from .grid import shoc_grid
from .length import shoc_length
from .pblintd import pblintd
from .second_moments import diag_second_shoc_moments
from .solver import update_prognostics_implicit
from .surface import shoc_diag_obklen
from .thermo import compute_shoc_temperature, compute_shoc_vapor
from .third_moments import diag_third_shoc_moments
from .tke import check_tke, shoc_tke


def shoc_init(nbot_shoc: int, ntop_shoc: int, pref_mid) -> int:
    """Maximum number of levels in the PBL from the surface
    (Functions::shoc_init): the largest (nbot_shoc - k) over levels
    ntop_shoc <= k < nbot_shoc with reference pressure >= pblmaxp; >= 1.
    Host-side setup (static), matching the C++ one-time computation.
    """
    pref = np.asarray(pref_mid)
    k = np.arange(ntop_shoc, nbot_shoc)
    eligible = pref[ntop_shoc:nbot_shoc] >= sc.pblmaxp
    levels_from_surface = np.where(eligible, nbot_shoc - k, 1)
    return int(max(1, levels_from_surface.max(initial=1)))


@functools.partial(jax.jit, static_argnames=(
    "nadv", "npbl", "shoc_1p5tke", "extra_diags"))
def shoc_main(dtime, nadv: int, npbl: int,
              # Runtime parameters
              lambda_low, lambda_high, lambda_slope, lambda_thresh,
              thl2tune, qw2tune, qwthl2tune, w2tune, length_fac,
              c_diag_3rd_mom, ckh, ckm, shoc_1p5tke: bool, extra_diags: bool,
              # Inputs
              dx, dy, zt_grid, zi_grid, pres, presi, pdel, thv, w_field,
              wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, wtracer_sfc, inv_exner, phis,
              # In/out state
              host_dse, tke, thetal, qw, u_wind, v_wind, wthv_sec, qtracers,
              tk, shoc_cldfrac, shoc_ql):
    """One SHOC host step (nadv internal subcycles of dtime each).

    Array conventions as everywhere in scream_jax: level axis last,
    midpoints nlev, interfaces nlev+1, leading axes batch columns;
    qtracers (..., num_qtracers, nlev); surface scalars (...,).

    Returns a dict with the updated prognostics, the pblh/ustar/obklen
    outputs, and all diagnostic moments.
    """
    thetal = jnp.asarray(thetal)
    qw = jnp.asarray(qw)
    shoc_ql = jnp.asarray(shoc_ql)

    # Total energy integrals before SHOC (for the energy fixer)
    se_b, ke_b, wv_b, wl_b = shoc_energy_integrals(
        host_dse, pdel, qw, shoc_ql, u_wind, v_wind)

    # Loop-carried diagnostics (defined on first iteration; nadv >= 1)
    pblh = ustar = obklen = None
    shoc_mix = brunt = isotropy = tkh = None
    thl_sec = qw_sec = wthl_sec = wqw_sec = qwthl_sec = None
    uw_sec = vw_sec = wtke_sec = w_sec = w3 = None
    wqls_sec = shoc_ql2 = shoc_cond = shoc_evap = None

    for _ in range(nadv):
        tke = check_tke(tke)

        dz_zt, dz_zi, rho_zt = shoc_grid(zt_grid, zi_grid, pdel)

        shoc_qv = compute_shoc_vapor(qw, shoc_ql)
        shoc_tabs = compute_shoc_temperature(thetal, shoc_ql, inv_exner)

        ustar, kbfs, obklen = shoc_diag_obklen(
            uw_sfc, vw_sfc, wthl_sfc, wqw_sfc,
            thetal[..., -1], shoc_ql[..., -1], shoc_qv[..., -1])

        pblh = pblintd(zt_grid, zi_grid, thetal, shoc_ql, shoc_qv,
                       u_wind, v_wind, ustar, obklen, kbfs,
                       shoc_cldfrac, npbl)

        brunt, shoc_mix = shoc_length(length_fac, dx, dy, zt_grid, zi_grid,
                                      dz_zt, tke, thv)

        tke, tk, tkh, isotropy = shoc_tke(
            dtime, lambda_low, lambda_high, lambda_slope, lambda_thresh,
            ckh, ckm, shoc_1p5tke, wthv_sec, shoc_mix, dz_zi, dz_zt,
            pres, shoc_tabs, u_wind, v_wind, brunt, zt_grid, zi_grid,
            pblh, tke, tk)

        thetal, qw, qtracers, tke, u_wind, v_wind = update_prognostics_implicit(
            dtime, dz_zt, dz_zi, rho_zt, zt_grid, zi_grid, tk, tkh,
            uw_sfc, vw_sfc, wthl_sfc, wqw_sfc, wtracer_sfc,
            thetal, qw, qtracers, tke, u_wind, v_wind)

        (thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec,
         wtke_sec, w_sec, _ustar2, _wstar) = diag_second_shoc_moments(
            thl2tune, qw2tune, qwthl2tune, w2tune, shoc_1p5tke,
            thetal, qw, u_wind, v_wind, tke, isotropy, tkh, tk,
            dz_zi, zt_grid, zi_grid, shoc_mix,
            wthl_sfc, wqw_sfc, uw_sfc, vw_sfc)

        w3 = diag_third_shoc_moments(
            c_diag_3rd_mom, shoc_1p5tke, w_sec, thl_sec, wthl_sec,
            isotropy, brunt, thetal, tke, dz_zt, dz_zi, zt_grid, zi_grid)

        (shoc_cldfrac, shoc_ql, wqls_sec, wthv_sec, shoc_ql2,
         shoc_cond, shoc_evap) = shoc_assumed_pdf(
            thetal, qw, w_field, thl_sec, qw_sec, dtime, extra_diags,
            wthl_sec, w_sec, wqw_sec, qwthl_sec, w3, pres,
            zt_grid, zi_grid, shoc_ql)

        tke = check_tke(tke)

    # Update host dry static energy and apply the energy fixer
    host_dse = update_host_dse(thetal, shoc_ql, inv_exner, zt_grid, phis)
    se_a, ke_a, wv_a, wl_a = shoc_energy_integrals(
        host_dse, pdel, qw, shoc_ql, u_wind, v_wind)
    _, _, rho_zt = shoc_grid(zt_grid, zi_grid, pdel)
    host_dse = shoc_energy_fixer(
        dtime, nadv, zt_grid, zi_grid, se_b, ke_b, wv_b, wl_b,
        se_a, ke_a, wv_a, wl_a, wthl_sfc, wqw_sfc, rho_zt, tke, presi,
        host_dse)

    # Final PBL diagnostics (no answer-changing code past this point)
    shoc_qv = compute_shoc_vapor(qw, shoc_ql)
    ustar, kbfs, obklen = shoc_diag_obklen(
        uw_sfc, vw_sfc, wthl_sfc, wqw_sfc,
        thetal[..., -1], shoc_ql[..., -1], shoc_qv[..., -1])
    pblh = pblintd(zt_grid, zi_grid, thetal, shoc_ql, shoc_qv,
                   u_wind, v_wind, ustar, obklen, kbfs, shoc_cldfrac, npbl)

    return {
        "host_dse": host_dse, "tke": tke, "thetal": thetal, "qw": qw,
        "u_wind": u_wind, "v_wind": v_wind, "wthv_sec": wthv_sec,
        "qtracers": qtracers, "tk": tk, "shoc_cldfrac": shoc_cldfrac,
        "shoc_ql": shoc_ql,
        "pblh": pblh, "ustar": ustar, "obklen": obklen,
        "shoc_ql2": shoc_ql2, "tkh": tkh,
        "shoc_cond": shoc_cond, "shoc_evap": shoc_evap,
        "shoc_mix": shoc_mix, "w_sec": w_sec, "thl_sec": thl_sec,
        "qw_sec": qw_sec, "qwthl_sec": qwthl_sec, "wthl_sec": wthl_sec,
        "wqw_sec": wqw_sec, "wtke_sec": wtke_sec, "uw_sec": uw_sec,
        "vw_sec": vw_sec, "w3": w3, "wqls_sec": wqls_sec,
        "brunt": brunt, "isotropy": isotropy,
    }
