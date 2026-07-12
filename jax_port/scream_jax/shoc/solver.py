"""Implicit vertical diffusion of the SHOC prognostics.

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_tridiag_solver_impl.hpp (vd_shoc_decomp, vd_shoc_solve),
  shoc_update_prognostics_implicit_impl.hpp

The tridiagonal solve delegates to the foundation Thomas solver
(transcribed from EKAT); the C++'s bfb/cr variants are parallel
decompositions of the same arithmetic.
"""

import jax
import jax.numpy as jnp

from ..foundation import constants as c
from ..foundation.tridiag import thomas
from .grid import compute_tmpi, dp_inverse
from .interp import linear_interp


@jax.jit
def vd_shoc_decomp(kv_term, tmpi, rdp_zt, dtime, flux):
    """Assemble the tridiagonal diffusion matrix (Functions::vd_shoc_decomp).

    kv_term, tmpi are on interfaces (..., nlev+1); rdp_zt on midpoints.
    flux (per-column, (...,) or scalar) enters the surface diagonal
    explicitly. Returns (dl, d, du), each (..., nlev), in the row convention
    dl(k) x(k-1) + d(k) x(k) + du(k) x(k+1) = rhs(k).
    """
    kv_term = jnp.asarray(kv_term)
    tmpi = jnp.asarray(tmpi)
    rdp_zt = jnp.asarray(rdp_zt)
    flux = jnp.asarray(flux)

    du = -kv_term[..., 1:] * tmpi[..., 1:] * rdp_zt
    dl = -kv_term[..., :-1] * tmpi[..., :-1] * rdp_zt
    du = du.at[..., -1].set(0.0)
    dl = dl.at[..., 0].set(0.0)

    d = 1.0 - du - dl
    d = d.at[..., -1].add(flux * dtime * c.gravit * rdp_zt[..., -1])
    return dl, d, du


@jax.jit
def vd_shoc_solve(dl, d, du, rhs):
    """Solve the diffusion system for one or more right-hand sides
    (Functions::vd_shoc_solve). rhs may carry extra leading axes (e.g.
    (..., nrhs, nlev) against (..., nlev) diagonals)."""
    dl, d, du = (jnp.asarray(a) for a in (dl, d, du))
    rhs = jnp.asarray(rhs)
    if rhs.ndim == dl.ndim + 1:
        dl, d, du = (a[..., None, :] for a in (dl, d, du))
    return thomas(dl, d, du, rhs)


@jax.jit
def update_prognostics_implicit(dtime, dz_zt, dz_zi, rho_zt, zt_grid, zi_grid,
                                tk, tkh, uw_sfc, vw_sfc, wthl_sfc, wqw_sfc,
                                wtracer_sfc, thetal, qw, qtracers, tke,
                                u_wind, v_wind):
    """Implicit vertical diffusion of winds, thetal, qw, tke, and tracers
    (Functions::update_prognostics_implicit).

    Args mirror the C++: interface quantities are derived internally by
    linear interpolation; surface fluxes of scalars are applied explicitly
    at the lowest midpoint before the solve; the wind solve embeds the
    surface drag ksrf implicitly in the matrix.

    qtracers has shape (..., num_qtracers, nlev); wtracer_sfc
    (..., num_qtracers). Surface flux args are per-column (...,).

    Returns (thetal, qw, qtracers, tke, u_wind, v_wind), updated.
    """
    zt_grid = jnp.asarray(zt_grid)
    zi_grid = jnp.asarray(zi_grid)
    thetal = jnp.asarray(thetal)
    qw = jnp.asarray(qw)
    qtracers = jnp.asarray(qtracers)
    tke = jnp.asarray(tke)
    u_wind = jnp.asarray(u_wind)
    v_wind = jnp.asarray(v_wind)
    wtracer_sfc = jnp.asarray(wtracer_sfc)
    uw_sfc = jnp.asarray(uw_sfc)
    vw_sfc = jnp.asarray(vw_sfc)

    nlev = zt_grid.shape[-1]

    # Interface diffusivities and density
    tkh_zi = linear_interp(zt_grid, zi_grid, jnp.asarray(tkh), nlev, nlev + 1, 0.0)
    tk_zi = linear_interp(zt_grid, zi_grid, jnp.asarray(tk), nlev, nlev + 1, 0.0)
    rho_zi = linear_interp(zt_grid, zi_grid, jnp.asarray(rho_zt), nlev, nlev + 1, 0.0)

    tmpi = compute_tmpi(dtime, rho_zi, dz_zi)
    rdp_zt = dp_inverse(rho_zt, dz_zt)

    # Implicit surface stress coefficient and surface TKE flux
    wsmin = 1.0
    ksrfmin = 1e-4
    ustarmin = 0.01
    rho_sfc = rho_zi[..., -1]
    taux = rho_sfc * uw_sfc
    tauy = rho_sfc * vw_sfc
    ws = jnp.maximum(jnp.sqrt(u_wind[..., -1] ** 2 + v_wind[..., -1] ** 2), wsmin)
    tau = jnp.sqrt(taux ** 2 + tauy ** 2)
    ksrf = jnp.maximum(tau / ws, ksrfmin)
    ustar = jnp.maximum(jnp.sqrt(jnp.sqrt(uw_sfc ** 2 + vw_sfc ** 2)), ustarmin)
    wtke_sfc = ustar ** 3

    # Explicit surface fluxes for scalars, applied at the lowest midpoint
    cmnfac = dtime * (c.gravit * rho_sfc * rdp_zt[..., -1])
    thetal = thetal.at[..., -1].add(cmnfac * jnp.asarray(wthl_sfc))
    qw = qw.at[..., -1].add(cmnfac * jnp.asarray(wqw_sfc))
    tke = tke.at[..., -1].add(cmnfac * wtke_sfc)
    qtracers = qtracers.at[..., -1].add(cmnfac[..., None] * wtracer_sfc)

    # Momentum solve (kv = tk_zi, implicit surface drag)
    dl, d, du = vd_shoc_decomp(tk_zi, tmpi, rdp_zt, dtime, ksrf)
    wind_rhs = jnp.stack([u_wind, v_wind], axis=-2)
    wind_sol = vd_shoc_solve(dl, d, du, wind_rhs)
    u_wind, v_wind = wind_sol[..., 0, :], wind_sol[..., 1, :]

    # Thermo/tracer solve (kv = tkh_zi, fluxes already applied explicitly)
    dl, d, du = vd_shoc_decomp(tkh_zi, tmpi, rdp_zt, dtime, 0.0)
    thermo_rhs = jnp.concatenate(
        [qtracers, thetal[..., None, :], qw[..., None, :], tke[..., None, :]],
        axis=-2)
    thermo_sol = vd_shoc_solve(dl, d, du, thermo_rhs)
    nq = qtracers.shape[-2]
    qtracers = thermo_sol[..., :nq, :]
    thetal = thermo_sol[..., nq, :]
    qw = thermo_sol[..., nq + 1, :]
    tke = thermo_sol[..., nq + 2, :]

    return thetal, qw, qtracers, tke, u_wind, v_wind
