"""RTE solvers: shortwave two-stream + adding, longwave no-scattering.

Sources: cpp/rte/kernels/mo_rte_solver_kernels.h (Kokkos variants:
sw_two_stream, sw_source_2str, adding, sw_solver_2stream,
lw_source_noscat_stencil, lw_transport_noscat, lw_solver_noscat,
lw_solver_noscat_GaussQuad) and cpp/rte/mo_rte_sw.h / mo_rte_lw.h
(driver logic: albedo/emissivity band->gpt expansion, boundary
conditions, flux reduction). EAMxx always runs with top_at_1 = True
(k=0 at model top) and 1 LW quadrature angle (D=1.66, w=0.5); only that
configuration is implemented.

All fields are (ncol, nlay(+1), ngpt); vertical scans are lax.scan.
"""

import jax.numpy as jnp
from jax import lax

_EPS = 2.220446049250313e-16  # numeric_limits<double>::epsilon()
_PI = 3.14159265358979323846
GAUSS_DS_1 = 1.66
GAUSS_WTS_1 = 0.5


def sw_two_stream(mu0, tau, w0, g):
    """Meador & Weaver two-stream layer coefficients (Zdunkowski PIFM).
    mu0: (ncol,). Returns (Rdif, Tdif, Rdir, Tdir, Tnoscat)."""
    mu0b = jnp.asarray(mu0)[:, None, None]
    mu0_inv = 1.0 / mu0b
    gamma1 = (8.0 - w0 * (5.0 + 3.0 * g)) * 0.25
    gamma2 = 3.0 * (w0 * (1.0 - g)) * 0.25
    gamma3 = (2.0 - 3.0 * mu0b * g) * 0.25
    gamma4 = 1.0 - gamma3

    alpha1 = gamma1 * gamma4 + gamma2 * gamma3
    alpha2 = gamma1 * gamma3 + gamma2 * gamma4
    k = jnp.sqrt(jnp.maximum((gamma1 - gamma2) * (gamma1 + gamma2), 1.0e-12))
    exp_minusktau = jnp.exp(-tau * k)
    exp_minus2ktau = exp_minusktau * exp_minusktau

    RT_term = 1.0 / (k * (1.0 + exp_minus2ktau)
                     + gamma1 * (1.0 - exp_minus2ktau))
    Rdif = RT_term * gamma2 * (1.0 - exp_minus2ktau)
    Tdif = RT_term * 2.0 * k * exp_minusktau
    Tnoscat = jnp.exp(-tau * mu0_inv)

    k_mu = k * mu0b
    k_gamma3 = k * gamma3
    k_gamma4 = k * gamma4
    denom = 1.0 - k_mu * k_mu
    denom = jnp.where(jnp.abs(denom) >= _EPS, denom, _EPS)
    RT_term2 = w0 * RT_term / denom

    Rdir = RT_term2 * ((1.0 - k_mu) * (alpha2 + k_gamma3)
                       - (1.0 + k_mu) * (alpha2 - k_gamma3) * exp_minus2ktau
                       - 2.0 * (k_gamma3 - alpha2 * k_mu) * exp_minusktau
                       * Tnoscat)
    Tdir = -RT_term2 * ((1.0 + k_mu) * (alpha1 + k_gamma4) * Tnoscat
                        - (1.0 - k_mu) * (alpha1 - k_gamma4)
                        * exp_minus2ktau * Tnoscat
                        - 2.0 * (k_gamma4 + alpha1 * k_mu) * exp_minusktau)
    return Rdif, Tdif, Rdir, Tdir, Tnoscat


def sw_source_2str(Rdir, Tdir, Tnoscat, sfc_albedo_dir_gpt, flux_dn_dir_top):
    """sw_source_2str (top_at_1). flux_dn_dir_top: (ncol, ngpt) direct
    beam at the domain top. Returns (source_up, source_dn, source_sfc,
    flux_dn_dir (ncol, nlay+1, ngpt))."""
    def step(fdir, args):
        rdir_l, tdir_l, tnos_l = args
        src_up = rdir_l * fdir
        src_dn = tdir_l * fdir
        fdir_next = tnos_l * fdir
        return fdir_next, (src_up, src_dn, fdir_next)

    args = (jnp.moveaxis(Rdir, 1, 0), jnp.moveaxis(Tdir, 1, 0),
            jnp.moveaxis(Tnoscat, 1, 0))
    _, (src_up, src_dn, fdir_lev) = lax.scan(step, flux_dn_dir_top, args)
    source_up = jnp.moveaxis(src_up, 0, 1)
    source_dn = jnp.moveaxis(src_dn, 0, 1)
    flux_dn_dir = jnp.concatenate([flux_dn_dir_top[:, None, :],
                                   jnp.moveaxis(fdir_lev, 0, 1)], axis=1)
    source_sfc = flux_dn_dir[:, -1, :] * sfc_albedo_dir_gpt
    return source_up, source_dn, source_sfc, flux_dn_dir


def adding(albedo_sfc_gpt, rdif, tdif, src_dn, src_up, src_sfc, flux_dn_top):
    """Shonk & Hogan 2008 adding (top_at_1). flux_dn_top: (ncol, ngpt)
    diffuse down at the top. Returns (flux_up, flux_dn), (ncol,nlay+1,ngpt)."""
    rdif_l = jnp.moveaxis(rdif, 1, 0)   # (nlay, ncol, ngpt)
    tdif_l = jnp.moveaxis(tdif, 1, 0)
    src_dn_l = jnp.moveaxis(src_dn, 1, 0)
    src_up_l = jnp.moveaxis(src_up, 1, 0)

    # bottom-up sweep: albedo/src at level ilev from ilev+1
    def up_step(carry, args):
        alb_below, src_below = carry
        r, t, sdn, sup = args
        denom = 1.0 / (1.0 - r * alb_below)                     # Eq 10
        alb = r + t * t * alb_below * denom                     # Eq 9
        src = sup + t * denom * (src_below + alb_below * sdn)   # Eq 11
        return (alb, src), (alb, src, denom)

    args_rev = (rdif_l[::-1], tdif_l[::-1], src_dn_l[::-1], src_up_l[::-1])
    (_, _), (alb_rev, src_rev, denom_rev) = lax.scan(
        up_step, (albedo_sfc_gpt, src_sfc), args_rev)
    # albedo/src at levels 0..nlay-1 (top..just-above-surface), denom per layer
    alb = alb_rev[::-1]      # level ilev in 0..nlay-1
    src = src_rev[::-1]
    denom = denom_rev[::-1]  # layer ilev

    albedo_sfc_l = albedo_sfc_gpt  # level nlay
    src_sfc_l = src_sfc

    # top-down flux sweep
    flux_up_top = flux_dn_top * alb[0] + src[0]                 # Eq 12

    alb_lev = jnp.concatenate([alb, albedo_sfc_l[None]], axis=0)
    src_lev = jnp.concatenate([src, src_sfc_l[None]], axis=0)

    def dn_step(fdn, args):
        t, r, sdn, dnm, alb_next, src_next = args
        fdn_next = (t * fdn + r * src_next + sdn) * dnm         # Eq 13
        fup_next = fdn_next * alb_next + src_next               # Eq 12
        return fdn_next, (fdn_next, fup_next)

    args_dn = (tdif_l, rdif_l, src_dn_l, denom, alb_lev[1:], src_lev[1:])
    _, (fdn_levs, fup_levs) = lax.scan(dn_step, flux_dn_top, args_dn)

    flux_dn = jnp.concatenate([flux_dn_top[:, None, :],
                               jnp.moveaxis(fdn_levs, 0, 1)], axis=1)
    flux_up = jnp.concatenate([flux_up_top[:, None, :],
                               jnp.moveaxis(fup_levs, 0, 1)], axis=1)
    return flux_up, flux_dn


def sw_solver_2stream(tau, ssa, g, mu0, sfc_alb_dir_gpt, sfc_alb_dif_gpt,
                      inc_flux):
    """sw_solver_2stream (top_at_1): BCs + two-stream + adding.
    inc_flux: (ncol, ngpt) TOA spectral flux. Returns
    (flux_up, flux_dn_total, flux_dn_dir), each (ncol, nlay+1, ngpt)."""
    Rdif, Tdif, Rdir, Tdir, Tnoscat = sw_two_stream(mu0, tau, ssa, g)

    flux_dn_dir_top = inc_flux * jnp.asarray(mu0)[:, None]  # apply_BC w/factor
    source_up, source_dn, source_sfc, flux_dir = sw_source_2str(
        Rdir, Tdir, Tnoscat, sfc_alb_dir_gpt, flux_dn_dir_top)

    flux_dn_top = jnp.zeros_like(inc_flux)  # apply_BC (no incident diffuse)
    flux_up, flux_dn = adding(sfc_alb_dif_gpt, Rdif, Tdif,
                              source_dn, source_up, source_sfc, flux_dn_top)
    flux_dn = flux_dn + flux_dir
    return flux_up, flux_dn, flux_dir


def lw_solver_noscat(tau, lay_source, lev_source_inc, lev_source_dec,
                     sfc_emis_gpt, sfc_src, D=GAUSS_DS_1, weight=GAUSS_WTS_1,
                     flux_dn_top=None):
    """lw_solver_noscat with a single quadrature angle (top_at_1).
    Returns (flux_up, flux_dn) on levels."""
    tau = jnp.asarray(tau)
    ncol, nlay, ngpt = tau.shape
    if flux_dn_top is None:
        flux_dn_top = jnp.zeros((ncol, ngpt), dtype=tau.dtype)

    tau_thresh = jnp.sqrt(jnp.asarray(_EPS, dtype=tau.dtype))

    # top_at_1: lev_source_up = lev_source_dec, lev_source_dn = lev_source_inc
    lev_source_up = lev_source_dec
    lev_source_dn = lev_source_inc

    radn_dn_top = flux_dn_top / (2.0 * _PI * weight)
    sfc_albedo = 1.0 - sfc_emis_gpt
    source_sfc = sfc_emis_gpt * sfc_src

    tau_loc = tau * D
    trans = jnp.exp(-tau_loc)
    # double-where: sanitize the denominator with the SAME predicate as
    # the branch selection so 0 < tau_loc <= tau_thresh never produces a
    # huge 1/tau_loc in the unselected branch (NaN-safe under AD;
    # bit-neutral on the primal).
    use_exact = tau_loc > tau_thresh
    fact = jnp.where(use_exact,
                     (1.0 - trans) / jnp.where(use_exact, tau_loc, 1.0)
                     - trans,
                     tau_loc * (0.5 - 1.0 / 3.0 * tau_loc))
    source_dn = (1.0 - trans) * lev_source_dn \
        + 2.0 * fact * (lay_source - lev_source_dn)
    source_up = (1.0 - trans) * lev_source_up \
        + 2.0 * fact * (lay_source - lev_source_up)

    # downward propagation
    def dn_step(rdn, args):
        tr, sdn = args
        rdn_next = tr * rdn + sdn
        return rdn_next, rdn_next

    trans_l = jnp.moveaxis(trans, 1, 0)
    sdn_l = jnp.moveaxis(source_dn, 1, 0)
    sup_l = jnp.moveaxis(source_up, 1, 0)
    rdn_bot, rdn_levs = lax.scan(dn_step, radn_dn_top, (trans_l, sdn_l))
    radn_dn = jnp.concatenate([radn_dn_top[:, None, :],
                               jnp.moveaxis(rdn_levs, 0, 1)], axis=1)

    # surface reflection + emission, then upward propagation
    radn_up_sfc = rdn_bot * sfc_albedo + source_sfc

    def up_step(rup, args):
        tr, sup = args
        rup_next = tr * rup + sup
        return rup_next, rup_next

    _, rup_levs = lax.scan(up_step, radn_up_sfc,
                           (trans_l[::-1], sup_l[::-1]))
    radn_up = jnp.concatenate([jnp.moveaxis(rup_levs[::-1], 0, 1),
                               radn_up_sfc[:, None, :]], axis=1)

    w = 2.0 * _PI * weight
    return w * radn_up, w * radn_dn


def expand_and_transpose(arr_bnd_T, gpt2band):
    """rte drivers' expand_and_transpose: (nbnd, ncol) -> (ncol, ngpt)."""
    return jnp.asarray(arr_bnd_T).T[:, jnp.asarray(gpt2band)]


def reduce_fluxes(gpt_flux, band2gpt=None):
    """FluxesBroadband/Byband reduce: sum over g-points; optionally also
    per-band sums. gpt_flux (ncol, nlay+1, ngpt)."""
    broadband = gpt_flux.sum(axis=-1)
    if band2gpt is None:
        return broadband
    nband = band2gpt.shape[1]
    parts = [gpt_flux[..., band2gpt[0, b]:band2gpt[1, b] + 1].sum(axis=-1)
             for b in range(nband)]
    return broadband, jnp.stack(parts, axis=-1)


def rte_sw(optics, mu0, inc_flux, sfc_alb_dir_T, sfc_alb_dif_T, gpt2band,
           band2gpt):
    """mo_rte_sw driver (top_at_1): returns dict with broadband
    flux_up/flux_dn/flux_dn_dir and by-band bnd_flux_* arrays."""
    sfc_alb_dir_gpt = expand_and_transpose(sfc_alb_dir_T, gpt2band)
    sfc_alb_dif_gpt = expand_and_transpose(sfc_alb_dif_T, gpt2band)
    fu, fd, fdir = sw_solver_2stream(optics["tau"], optics["ssa"],
                                     optics["g"], mu0, sfc_alb_dir_gpt,
                                     sfc_alb_dif_gpt, inc_flux)
    up, bnd_up = reduce_fluxes(fu, band2gpt)
    dn, bnd_dn = reduce_fluxes(fd, band2gpt)
    dnd, bnd_dnd = reduce_fluxes(fdir, band2gpt)
    return {"flux_up": up, "flux_dn": dn, "flux_dn_dir": dnd,
            "bnd_flux_up": bnd_up, "bnd_flux_dn": bnd_dn,
            "bnd_flux_dn_dir": bnd_dnd}


def rte_lw(optics, sources, sfc_emis_T, gpt2band, band2gpt):
    """mo_rte_lw driver (top_at_1, 1 quadrature angle): returns dict with
    broadband and by-band fluxes."""
    sfc_emis_gpt = expand_and_transpose(sfc_emis_T, gpt2band)
    fu, fd = lw_solver_noscat(optics["tau"], sources["lay_src"],
                              sources["lev_src_inc"], sources["lev_src_dec"],
                              sfc_emis_gpt, sources["sfc_src"])
    up, bnd_up = reduce_fluxes(fu, band2gpt)
    dn, bnd_dn = reduce_fluxes(fd, band2gpt)
    return {"flux_up": up, "flux_dn": dn,
            "bnd_flux_up": bnd_up, "bnd_flux_dn": bnd_dn}
