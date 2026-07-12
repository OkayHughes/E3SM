"""TKE chain: production, dissipation, isotropy, eddy diffusivities.

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_check_tke_impl.hpp, shoc_integ_column_stability_impl.hpp,
  shoc_compute_shr_prod_impl.hpp, shoc_adv_sgs_tke_impl.hpp,
  shoc_isotropic_ts_impl.hpp, shoc_eddy_diffusivities_impl.hpp,
  shoc_tke_impl.hpp (driver)
"""

import functools

import jax
import jax.numpy as jnp

from ..foundation import constants as c
from . import constants as sc
from .interp import linear_interp


@jax.jit
def check_tke(tke):
    """Clip TKE to its minimum allowed value (Functions::check_tke)."""
    return jnp.maximum(jnp.asarray(tke), sc.mintke)


@jax.jit
def integ_column_stability(dz_zt, pres, brunt):
    """Lower-troposphere (p > 800 hPa) column integral of brunt
    (Functions::integ_column_stability). Returns shape (...,)."""
    troppres = 80000.0
    contrib = jnp.where(jnp.asarray(pres) > troppres,
                        jnp.asarray(dz_zt) * jnp.asarray(brunt), 0.0)
    return jnp.sum(contrib, axis=-1)


@jax.jit
def compute_shr_prod(dz_zi, u_wind, v_wind):
    """TKE shear production on interfaces (Functions::compute_shr_prod),
    Bretherton & Park (2010):

        sterm(k) = Ck_sh * [ (du/dz)^2 + (dv/dz)^2 ]   for k = 1..nlev-1
        sterm(0) = sterm(nlev) = 0

    u/v on midpoints (nlev); dz_zi and the result on interfaces (nlev+1).
    """
    ck_sh = 0.1
    u = jnp.asarray(u_wind)
    v = jnp.asarray(v_wind)
    dz_zi = jnp.asarray(dz_zi)

    band = slice(1, u.shape[-1])
    grid_dz = 1.0 / dz_zi[..., band]
    u_grad = grid_dz * (u[..., :-1] - u[..., 1:])
    v_grad = grid_dz * (v[..., :-1] - v[..., 1:])
    interior = ck_sh * (u_grad ** 2 + v_grad ** 2)

    zeros = jnp.zeros_like(dz_zi[..., :1])
    return jnp.concatenate([zeros, interior, zeros], axis=-1)


@functools.partial(jax.jit, static_argnames=("shoc_1p5tke",))
def adv_sgs_tke(dtime, shoc_1p5tke: bool, shoc_mix, wthv_sec, sterm_zt,
                tk, brunt, tke):
    """Advance the SGS TKE equation one timestep (Functions::adv_sgs_tke).

    Returns (tke_new, a_diss). Operation order matches the C++ exactly:
    the buoyancy production term uses the UNCLIPPED input tke.
    """
    tke = jnp.asarray(tke)
    brunt = jnp.asarray(brunt)

    cs = 0.15
    ck = 0.1
    ce = ck ** 3 / (cs ** 2 * cs ** 2)
    ce1 = ce / 0.7 * 0.19
    ce2 = ce / 0.7 * 0.51
    cee = ce1 + ce2

    if shoc_1p5tke:
        # No SGS variability: buoyant production closed via local moist brunt
        a_prod_bu = -tke * brunt
    else:
        a_prod_bu = (c.gravit / c.basetemp) * jnp.asarray(wthv_sec)

    tke = jnp.maximum(0.0, tke)
    a_prod_sh = jnp.asarray(tk) * jnp.asarray(sterm_zt)
    a_diss = cee / jnp.asarray(shoc_mix) * tke ** 1.5

    tke_new = jnp.maximum(
        sc.mintke,
        tke + dtime * (jnp.maximum(0.0, a_prod_sh + a_prod_bu) - a_diss))
    return jnp.minimum(tke_new, sc.maxtke), a_diss


@jax.jit
def isotropic_ts(lambda_low, lambda_high, lambda_slope, lambda_thresh,
                 brunt_int, tke, a_diss, brunt):
    """Return-to-isotropy timescale (Functions::isotropic_ts).

    brunt_int is per-column (...,); the rest are (..., nlev).
    """
    tke = jnp.asarray(tke)
    brunt = jnp.asarray(brunt)
    brunt_int = jnp.asarray(brunt_int)
    if brunt_int.ndim:
        brunt_int = brunt_int[..., None]

    tscale = 2.0 * tke / jnp.asarray(a_diss)

    lam = lambda_low + (brunt_int / c.gravit - lambda_thresh) * lambda_slope
    lam = jnp.maximum(lambda_low, jnp.minimum(lambda_high, lam))
    lam = jnp.where(brunt <= 0, 0.0, lam)

    return jnp.minimum(sc.maxiso, tscale / (1.0 + lam * brunt * tscale ** 2))


@jax.jit
def eddy_diffusivities(ckh, ckm, pblh, zt_grid, tabs, shoc_mix,
                       sterm_zt, isotropy, tke):
    """Eddy diffusivities of heat (tkh) and momentum (tk)
    (Functions::eddy_diffusivities).

    Near-runaway-cooling columns (surface tabs < 182 K) use undamped
    stable-PBL forms below pblh+200 m; otherwise tkh = Ckh*isotropy*tke etc.
    pblh is per-column (...,). Returns (tkh, tk).
    """
    tabs_crit = 182.0
    pbl_trans = 200.0
    ckh_s = 0.1
    ckm_s = 0.1

    zt_grid = jnp.asarray(zt_grid)
    tabs = jnp.asarray(tabs)
    shoc_mix = jnp.asarray(shoc_mix)
    sterm_zt = jnp.asarray(sterm_zt)
    isotropy = jnp.asarray(isotropy)
    tke = jnp.asarray(tke)
    pblh = jnp.asarray(pblh)
    if pblh.ndim:
        pblh = pblh[..., None]

    condition = (zt_grid < pblh + pbl_trans) & (tabs[..., -1:] < tabs_crit)
    stable_form = shoc_mix ** 2 * jnp.sqrt(sterm_zt)
    tkh = jnp.where(condition, ckh_s * stable_form, ckh * isotropy * tke)
    tk = jnp.where(condition, ckm_s * stable_form, ckm * isotropy * tke)
    return tkh, tk


@functools.partial(jax.jit, static_argnames=("shoc_1p5tke",))
def shoc_tke(dtime, lambda_low, lambda_high, lambda_slope, lambda_thresh,
             ckh, ckm, shoc_1p5tke: bool,
             wthv_sec, shoc_mix, dz_zi, dz_zt, pres, tabs,
             u_wind, v_wind, brunt, zt_grid, zi_grid, pblh, tke, tk):
    """TKE driver (Functions::shoc_tke).

    Returns (tke, tk, tkh, isotropy). tk in the argument list is the
    previous-step momentum diffusivity (used in shear production).
    """
    zt_grid = jnp.asarray(zt_grid)
    nlev = zt_grid.shape[-1]

    brunt_int = integ_column_stability(dz_zt, pres, brunt)
    sterm = compute_shr_prod(dz_zi, u_wind, v_wind)
    sterm_zt = linear_interp(jnp.asarray(zi_grid), zt_grid, sterm,
                             nlev + 1, nlev, 0.0)
    tke, a_diss = adv_sgs_tke(dtime, shoc_1p5tke, shoc_mix, wthv_sec,
                              sterm_zt, tk, brunt, tke)
    isotropy = isotropic_ts(lambda_low, lambda_high, lambda_slope,
                            lambda_thresh, brunt_int, tke, a_diss, brunt)
    tkh, tk = eddy_diffusivities(ckh, ckm, pblh, zt_grid, tabs, shoc_mix,
                                 sterm_zt, isotropy, tke)
    return tke, tk, tkh, isotropy
