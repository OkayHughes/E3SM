"""PBL height diagnosis (pblintd chain, from CAM's pblintd).

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_pblintd_init_pot_impl.hpp, shoc_pblintd_height_impl.hpp,
  shoc_pblintd_surf_temp_impl.hpp, shoc_pblintd_check_pblh_impl.hpp,
  shoc_pblintd_cldcheck_impl.hpp, shoc_pblintd_impl.hpp (driver)

The C++ per-column `bool& check` in/out flags become boolean arrays
threaded functionally through the chain.
"""

import functools

import jax
import jax.numpy as jnp

from ..foundation import constants as c

_RICR = 0.3
_FAC = 100.0
_TINY = 1e-36
_FAK = 8.5
_BETAM = 15.0
_SFFRAC = 0.1
_BINM = _BETAM * _SFFRAC


@jax.jit
def shoc_pblintd_init_pot(thl, ql, q):
    """Virtual potential temperature from thl, ql, qv
    (Functions::shoc_pblintd_init_pot)."""
    th = jnp.asarray(thl) + (c.LatVap / c.Cpair) * jnp.asarray(ql)
    return th * (1.0 + c.ZVIR * jnp.asarray(q) - jnp.asarray(ql))


@functools.partial(jax.jit, static_argnames=("npbl",))
def pblintd_height(z, u, v, ustar, thv, thv_ref, pblh, rino, check, npbl: int):
    """Richardson-number PBL height sweep (Functions::pblintd_height).

    Scans levels k in [nlev-npbl, nlev-1) for the lowest level with bulk
    Richardson number >= 0.3 and interpolates the crossing height. Columns
    with check=False pass through unchanged; columns where a crossing is
    found get check=False. Returns (pblh, rino, check).
    """
    z = jnp.asarray(z)
    u = jnp.asarray(u)
    v = jnp.asarray(v)
    thv = jnp.asarray(thv)
    thv_ref = jnp.asarray(thv_ref)
    pblh = jnp.asarray(pblh)
    rino = jnp.asarray(rino)
    check = jnp.asarray(check)
    ustar = jnp.asarray(ustar)

    nlev = z.shape[-1]
    band = slice(nlev - npbl, nlev - 1)

    vvk = jnp.maximum(
        _TINY,
        (u[..., band] - u[..., -1:]) ** 2
        + (v[..., band] - v[..., -1:]) ** 2
        + _FAC * (ustar[..., None] ** 2))
    rino_band = (c.gravit * (thv[..., band] - thv_ref[..., None])
                 * (z[..., band] - z[..., -1:]) / (thv[..., -1:] * vvk))
    rino = jnp.where(check[..., None], rino.at[..., band].set(rino_band), rino)

    idx = jnp.arange(nlev)
    in_band = (idx >= nlev - npbl) & (idx < nlev - 1)
    cand = jnp.where(in_band & (rino >= _RICR), idx, -1)
    max_indx = jnp.max(cand, axis=-1)
    found = (max_indx >= 0) & check

    mi = jnp.clip(max_indx, 0, nlev - 2)
    z_mi = jnp.take_along_axis(z, mi[..., None], axis=-1)[..., 0]
    z_mi1 = jnp.take_along_axis(z, mi[..., None] + 1, axis=-1)[..., 0]
    r_mi = jnp.take_along_axis(rino, mi[..., None], axis=-1)[..., 0]
    r_mi1 = jnp.take_along_axis(rino, mi[..., None] + 1, axis=-1)[..., 0]

    denom = r_mi - r_mi1
    denom = jnp.where(denom == 0, 1.0, denom)
    pblh_new = z_mi1 + (_RICR - r_mi1) / denom * (z_mi - z_mi1)

    pblh = jnp.where(found, pblh_new, pblh)
    check = check & ~found
    return pblh, rino, check


@functools.partial(jax.jit, static_argnames=("npbl",))
def pblintd_surf_temp(z, ustar, obklen, kbfs, thv, pblh, check, rino, npbl: int):
    """Effective surface (convective excess) temperature
    (Functions::pblintd_surf_temp). Returns (tlv, pblh, check, rino)."""
    z = jnp.asarray(z)
    thv = jnp.asarray(thv)
    pblh = jnp.asarray(pblh)
    check = jnp.asarray(check)
    rino = jnp.asarray(rino)
    kbfs = jnp.asarray(kbfs)
    ustar = jnp.asarray(ustar)
    obklen = jnp.asarray(obklen)

    nlev = z.shape[-1]
    nlevi = nlev + 1

    # Unresolved columns default to the top of the allowed PBL depth
    pblh = jnp.where(check, z[..., nlevi - npbl - 1], pblh)
    check = kbfs > 0.0

    tlv = thv[..., -1]
    phiminv = jnp.cbrt(1.0 - _BINM * pblh / obklen)
    tlv = jnp.where(check, tlv + kbfs * _FAK / (ustar * phiminv), tlv)
    rino = jnp.where(check[..., None], rino.at[..., -1].set(0.0), rino)
    return tlv, pblh, check, rino


@functools.partial(jax.jit, static_argnames=("npbl",))
def pblintd_check_pblh(z, ustar, check, pblh, npbl: int):
    """Fallback and mechanical-minimum bound (Functions::pblintd_check_pblh)."""
    z = jnp.asarray(z)
    nlev = z.shape[-1]
    nlevi = nlev + 1
    pblh = jnp.where(jnp.asarray(check), z[..., nlevi - npbl - 1],
                     jnp.asarray(pblh))
    return jnp.maximum(pblh, 700.0 * jnp.asarray(ustar))


@jax.jit
def shoc_pblintd_cldcheck(zi_low, cldn_low, pblh):
    """Raise PBL top to just above the lowest interface where cloud is
    present (Functions::shoc_pblintd_cldcheck; the C++ condition is
    cldn >= 0, i.e. effectively always for physical cloud fractions)."""
    return jnp.where(jnp.asarray(cldn_low) >= 0.0,
                     jnp.maximum(jnp.asarray(pblh), jnp.asarray(zi_low) + 50.0),
                     jnp.asarray(pblh))


@functools.partial(jax.jit, static_argnames=("npbl",))
def pblintd(z, zi, thl, ql, q, u, v, ustar, obklen, kbfs, cldn, npbl: int):
    """PBL height driver (Functions::pblintd). Returns pblh (...,).

    z/thl/ql/q/u/v/cldn on midpoints; zi on interfaces; surface scalars
    per column.
    """
    z = jnp.asarray(z)
    thv = shoc_pblintd_init_pot(thl, ql, q)

    rino = jnp.zeros_like(z)
    pblh = z[..., -1]
    check = jnp.ones(z.shape[:-1], dtype=bool)

    pblh, rino, check = pblintd_height(z, u, v, ustar, thv, thv[..., -1],
                                       pblh, rino, check, npbl)
    tlv, pblh, check, rino = pblintd_surf_temp(z, ustar, obklen, kbfs, thv,
                                               pblh, check, rino, npbl)
    pblh, rino, check = pblintd_height(z, u, v, ustar, thv, tlv,
                                       pblh, rino, check, npbl)
    pblh = pblintd_check_pblh(z, ustar, check, pblh, npbl)
    return shoc_pblintd_cldcheck(jnp.asarray(zi)[..., -2],
                                 jnp.asarray(cldn)[..., -1], pblh)
