"""Second-moment building blocks (fluxes, variances/covariances).

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_calc_shoc_vertflux_impl.hpp, shoc_calc_shoc_varorcovar_impl.hpp
  (the diag_second_* driver chain follows in this module.)

These kernels write only the interior interface levels k = 1..nlev-1;
boundary values are owned by the lbycond/ubycond kernels. They therefore
take the current output array and return a functional update of its
interior band.
"""

import jax
import jax.numpy as jnp


@jax.jit
def calc_shoc_vertflux(tkh_zi, dz_zi, invar, vertflux):
    """Downgradient vertical turbulent flux on interior interfaces
    (Functions::calc_shoc_vertflux):

        vertflux(k) = -tkh_zi(k) * (invar(k-1) - invar(k)) / dz_zi(k)

    invar has nlev entries; tkh_zi, dz_zi, vertflux have nlev+1. Only
    k = 1..nlev-1 of vertflux is replaced.
    """
    tkh_zi = jnp.asarray(tkh_zi)
    dz_zi = jnp.asarray(dz_zi)
    invar = jnp.asarray(invar)
    vertflux = jnp.asarray(vertflux)

    band = slice(1, invar.shape[-1])  # 1..nlev-1 inclusive
    flux = -(tkh_zi[..., band] / dz_zi[..., band]
             * (invar[..., :-1] - invar[..., 1:]))
    return vertflux.at[..., band].set(flux)


@jax.jit
def calc_shoc_varorcovar(tunefac, isotropy_zi, tkh_zi, dz_zi,
                         invar1, invar2, varorcovar):
    """SGS variance/covariance on interior interfaces
    (Functions::calc_shoc_varorcovar):

        varorcovar(k) = tunefac * isotropy_zi(k) * tkh_zi(k)
                        * (d invar1)(d invar2) / dz_zi(k)^2

    invar1/2 have nlev entries; the interface arrays nlev+1. Only
    k = 1..nlev-1 of varorcovar is replaced.
    """
    isotropy_zi = jnp.asarray(isotropy_zi)
    tkh_zi = jnp.asarray(tkh_zi)
    dz_zi = jnp.asarray(dz_zi)
    invar1 = jnp.asarray(invar1)
    invar2 = jnp.asarray(invar2)
    varorcovar = jnp.asarray(varorcovar)

    band = slice(1, invar1.shape[-1])
    d1 = invar1[..., :-1] - invar1[..., 1:]
    d2 = invar2[..., :-1] - invar2[..., 1:]
    val = (tunefac * (isotropy_zi[..., band] * tkh_zi[..., band])
           * d1 * d2 / dz_zi[..., band] ** 2)
    return varorcovar.at[..., band].set(val)


from ..foundation import constants as c  # noqa: E402
from .interp import linear_interp  # noqa: E402
import functools  # noqa: E402


@jax.jit
def shoc_diag_second_moments_srf(wthl_sfc, uw_sfc, vw_sfc):
    """Surface properties for the lower-boundary moments
    (Functions::shoc_diag_second_moments_srf). Returns (ustar2, wstar).

    NB: ustar2 is the surface momentum-flux magnitude sqrt(uw^2+vw^2),
    transcribed as named in the C++/Fortran.
    """
    wthl_sfc = jnp.asarray(wthl_sfc)
    ustar2 = jnp.sqrt(jnp.asarray(uw_sfc) ** 2 + jnp.asarray(vw_sfc) ** 2)
    # where(wthl_sfc >= 0, cbrt(max(0, coef*wthl_sfc)), 0) rewritten with
    # the double-where idiom: cbrt'(0) is inf, so reverse-mode AD NaNs on
    # wthl_sfc < 0 lanes. coef > 0, and cbrt(0) = 0 matches the fallback,
    # so the primal is unchanged.
    warg = jnp.maximum(0.0, (1.0 / c.basetemp) * c.gravit * wthl_sfc)
    warg_pos = warg > 0.0
    wstar = jnp.where(warg_pos,
                      jnp.cbrt(jnp.where(warg_pos, warg, 1.0)), 0.0)
    return ustar2, wstar


@jax.jit
def shoc_diag_second_moments_lbycond(wthl_sfc, wqw_sfc, uw_sfc, vw_sfc,
                                     ustar2, wstar):
    """Lower-boundary (surface) values of the second moments
    (Functions::shoc_diag_second_moments_lbycond), Andre et al. (1978).

    Returns (wthl_sec, wqw_sec, uw_sec, vw_sec, wtke_sec,
             thl_sec, qw_sec, qwthl_sec) surface scalars.
    """
    ufmin = 0.01
    wthl_sfc = jnp.asarray(wthl_sfc)
    wqw_sfc = jnp.asarray(wqw_sfc)
    ustar2 = jnp.asarray(ustar2)
    wstar = jnp.asarray(wstar)

    uf = jnp.maximum(ufmin, jnp.sqrt(ustar2 + 0.3 * wstar * wstar))

    thl_sec = 0.4 * 1.8 * (wthl_sfc / uf) ** 2
    qw_sec = 0.4 * 1.8 * (wqw_sfc / uf) ** 2
    qwthl_sec = 0.2 * 1.8 * (wthl_sfc / uf) * (wqw_sfc / uf)

    wtke_sec = jnp.maximum(jnp.sqrt(ustar2), ufmin) ** 3
    return (wthl_sfc, wqw_sfc, jnp.asarray(uw_sfc), jnp.asarray(vw_sfc),
            wtke_sec, thl_sec, qw_sec, qwthl_sec)


@functools.partial(jax.jit, static_argnames=("shoc_1p5tke",))
def diag_second_moments(thl2tune, qw2tune, qwthl2tune, w2tune,
                        shoc_1p5tke: bool,
                        thetal, qw, u_wind, v_wind, tke, isotropy, tkh, tk,
                        dz_zi, zt_grid, zi_grid, shoc_mix,
                        thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec,
                        uw_sec, vw_sec, wtke_sec):
    """Interior second moments (Functions::diag_second_moments).

    Updates the interior interface band of the *_sec arrays (boundaries are
    owned by the lbycond/ubycond kernels) and computes w_sec on midpoints.
    Returns (isotropy_zi, tkh_zi, tk_zi, thl_sec, qw_sec, wthl_sec, wqw_sec,
    qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec).
    """
    zt_grid = jnp.asarray(zt_grid)
    tke = jnp.asarray(tke)
    nlev = zt_grid.shape[-1]

    isotropy_zi = linear_interp(zt_grid, jnp.asarray(zi_grid),
                                jnp.asarray(isotropy), nlev, nlev + 1, 0.0)
    tkh_zi = linear_interp(zt_grid, jnp.asarray(zi_grid), jnp.asarray(tkh),
                           nlev, nlev + 1, 0.0)
    tk_zi = linear_interp(zt_grid, jnp.asarray(zi_grid), jnp.asarray(tk),
                          nlev, nlev + 1, 0.0)

    # Vertical velocity variance ~ TKE (zero under the 1.5-TKE closure)
    if shoc_1p5tke:
        w_sec = jnp.zeros_like(tke)
        # No SGS variability: zero the FULL variance/covariance arrays
        # (including boundary entries set by lbycond — as in the C++).
        thl_sec = jnp.zeros_like(jnp.asarray(thl_sec))
        qw_sec = jnp.zeros_like(jnp.asarray(qw_sec))
        qwthl_sec = jnp.zeros_like(jnp.asarray(qwthl_sec))
    else:
        w_sec = w2tune * (2.0 / 3.0) * tke
        thl_sec = calc_shoc_varorcovar(thl2tune, isotropy_zi, tkh_zi, dz_zi,
                                       thetal, thetal, thl_sec)
        qw_sec = calc_shoc_varorcovar(qw2tune, isotropy_zi, tkh_zi, dz_zi,
                                      qw, qw, qw_sec)
        qwthl_sec = calc_shoc_varorcovar(qwthl2tune, isotropy_zi, tkh_zi, dz_zi,
                                         thetal, qw, qwthl_sec)

    wthl_sec = calc_shoc_vertflux(tkh_zi, dz_zi, thetal, wthl_sec)
    wqw_sec = calc_shoc_vertflux(tkh_zi, dz_zi, qw, wqw_sec)
    wtke_sec = calc_shoc_vertflux(tkh_zi, dz_zi, tke, wtke_sec)
    uw_sec = calc_shoc_vertflux(tk_zi, dz_zi, u_wind, uw_sec)
    vw_sec = calc_shoc_vertflux(tk_zi, dz_zi, v_wind, vw_sec)

    return (isotropy_zi, tkh_zi, tk_zi, thl_sec, qw_sec, wthl_sec, wqw_sec,
            qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec)


@functools.partial(jax.jit, static_argnames=("shoc_1p5tke",))
def diag_second_shoc_moments(thl2tune, qw2tune, qwthl2tune, w2tune,
                             shoc_1p5tke: bool,
                             thetal, qw, u_wind, v_wind, tke, isotropy,
                             tkh, tk, dz_zi, zt_grid, zi_grid, shoc_mix,
                             wthl_sfc, wqw_sfc, uw_sfc, vw_sfc):
    """Second-moments driver (Functions::diag_second_shoc_moments).

    Sequence: surface properties -> lower bc -> interior band -> upper bc
    (zeros). Returns (thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec,
    vw_sec, wtke_sec, w_sec, ustar2, wstar); *_sec on interfaces, w_sec on
    midpoints.
    """
    zt_grid = jnp.asarray(zt_grid)
    nlevi = zt_grid.shape[-1] + 1
    iface_shape = zt_grid.shape[:-1] + (nlevi,)
    zeros_i = jnp.zeros(iface_shape, dtype=zt_grid.dtype)

    ustar2, wstar = shoc_diag_second_moments_srf(wthl_sfc, uw_sfc, vw_sfc)

    (wthl_b, wqw_b, uw_b, vw_b, wtke_b, thl_b, qw_b, qwthl_b) = \
        shoc_diag_second_moments_lbycond(wthl_sfc, wqw_sfc, uw_sfc, vw_sfc,
                                         ustar2, wstar)

    thl_sec = zeros_i.at[..., -1].set(thl_b)
    qw_sec = zeros_i.at[..., -1].set(qw_b)
    qwthl_sec = zeros_i.at[..., -1].set(qwthl_b)
    wthl_sec = zeros_i.at[..., -1].set(wthl_b)
    wqw_sec = zeros_i.at[..., -1].set(wqw_b)
    uw_sec = zeros_i.at[..., -1].set(uw_b)
    vw_sec = zeros_i.at[..., -1].set(vw_b)
    wtke_sec = zeros_i.at[..., -1].set(wtke_b)

    (_, _, _, thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec,
     uw_sec, vw_sec, wtke_sec, w_sec) = diag_second_moments(
        thl2tune, qw2tune, qwthl2tune, w2tune, shoc_1p5tke,
        thetal, qw, u_wind, v_wind, tke, isotropy, tkh, tk,
        dz_zi, zt_grid, zi_grid, shoc_mix,
        thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec,
        uw_sec, vw_sec, wtke_sec)

    # Upper boundary: all moments zero at the model top interface
    for_zero = (thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec,
                uw_sec, vw_sec, wtke_sec)
    (thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec,
     uw_sec, vw_sec, wtke_sec) = (a.at[..., 0].set(0.0) for a in for_zero)

    return (thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec,
            uw_sec, vw_sec, wtke_sec, w_sec, ustar2, wstar)
