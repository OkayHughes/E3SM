"""EAMxx RRTMGP interface drivers: rrtmgp_sw, rrtmgp_lw, rrtmgp_main and
the surrounding helpers.

Source: components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_interface.hpp.

Differences in mechanics (not in results): the C++ subsets daytime
columns before the SW solve; here every column is computed (each column
is independent, so day columns get identical results) with mu0
sanitized for night columns and their fluxes zeroed afterward, matching
the C++'s zero-initialized night fluxes exactly.

EAMxx configuration assumed: top_at_1 (k=0 at model top), LW surface
emissivity 0.98, one LW quadrature angle.
"""

import jax.numpy as jnp

from . import gas_optics as go
from . import optical_props as op
from . import cloud_optics as co
from . import mcica
from . import rte

GRAVIT = 9.80665
CPAIR = 1004.64  # scream::physics::Constants Cpair


def compute_band_by_band_surface_albedos(band_lims_wvn, alb_dir_vis,
                                         alb_dir_nir, alb_dif_vis,
                                         alb_dif_nir):
    """-> (sfc_alb_dir, sfc_alb_dif), each (ncol, nswbands). Bands fully
    above 14286 cm^-1 (0.7 um) are visible, fully below are near-IR,
    straddling bands average the two."""
    thresh = 14286.0
    vis1 = band_lims_wvn[0] > thresh
    vis2 = band_lims_wvn[1] > thresh
    w_vis = jnp.where(vis1 & vis2, 1.0, jnp.where(~vis1 & ~vis2, 0.0, 0.5))
    dir_ = (jnp.asarray(alb_dir_vis)[:, None] * w_vis
            + jnp.asarray(alb_dir_nir)[:, None] * (1.0 - w_vis))
    dif_ = (jnp.asarray(alb_dif_vis)[:, None] * w_vis
            + jnp.asarray(alb_dif_nir)[:, None] * (1.0 - w_vis))
    return dir_, dif_


def compute_broadband_surface_fluxes(band_lims_wvn, bnd_flux_dir_sfc,
                                     bnd_flux_dif_sfc):
    """-> (dir_vis, dir_nir, dif_vis, dif_nir), each (ncol,). Inputs are
    the surface-level by-band fluxes (ncol, nswbands)."""
    thresh = 14286.0
    vis1 = band_lims_wvn[0] > thresh
    vis2 = band_lims_wvn[1] > thresh
    w_vis = jnp.where(vis1 & vis2, 1.0, jnp.where(~vis1 & ~vis2, 0.0, 0.5))
    dir_vis = (bnd_flux_dir_sfc * w_vis).sum(axis=-1)
    dir_nir = (bnd_flux_dir_sfc * (1.0 - w_vis)).sum(axis=-1)
    dif_vis = (bnd_flux_dif_sfc * w_vis).sum(axis=-1)
    dif_nir = (bnd_flux_dif_sfc * (1.0 - w_vis)).sum(axis=-1)
    return dir_vis, dir_nir, dif_vis, dif_nir


def mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp):
    """Layer cloud mass per unit area [kg/m2]; in-cloud mixing ratio
    limited to 0.005 kg/kg, cloud fraction floored at 1e-4."""
    q = jnp.asarray(mixing_ratio)
    cf = jnp.asarray(cloud_fraction)
    incld = jnp.minimum(q / jnp.maximum(0.0001, cf), 0.005)
    return jnp.where(cf > 0, incld * jnp.asarray(dp) / GRAVIT, 0.0)


def compute_heating_rate(flux_up, flux_dn, pdel):
    """(ncol, nlay) heating rate from level fluxes."""
    return ((flux_up[:, 1:] - flux_up[:, :-1]
             - flux_dn[:, 1:] + flux_dn[:, :-1])
            * GRAVIT / (CPAIR * jnp.asarray(pdel)))


def compute_cloud_area(pmid, cld_tau_gpt, pmin, pmax):
    """Vertically-projected cloud area for layers with pmin <= p < pmax."""
    cloudy = (jnp.asarray(cld_tau_gpt) > 0) \
        & (jnp.asarray(pmid)[..., None] >= pmin) \
        & (jnp.asarray(pmid)[..., None] < pmax)
    subcol_mask = jnp.any(cloudy, axis=1)  # (ncol, ngpt)
    ngpt = cld_tau_gpt.shape[-1]
    return subcol_mask.astype(jnp.float64).sum(axis=-1) * (1.0 / ngpt)


def compute_aerocom_cloudtop(tmid, pmid, p_del, z_del, qc, qi, rel, rei,
                             cldfrac_tot, nc):
    """AeroCom cloud-top diagnostics; returns a dict of (ncol,) arrays."""
    tmid = jnp.asarray(tmid)
    cf = jnp.asarray(cldfrac_tot)
    qc = jnp.asarray(qc)
    qi = jnp.asarray(qi)
    ncol, nlay = tmid.shape
    thr = 0.001

    outs = {k: jnp.zeros(ncol) for k in
            ("T_mid_at_cldtop", "p_mid_at_cldtop", "cldfrac_ice_at_cldtop",
             "cldfrac_liq_at_cldtop", "cdnc_at_cldtop",
             "eff_radius_qc_at_cldtop", "eff_radius_qi_at_cldtop")}
    clr = jnp.ones(ncol)
    for ilay in range(1, nlay):
        active = ((qc[:, ilay] + qi[:, ilay]) > 0.0) & (cf[:, ilay] > thr)
        tmp = clr * (1.0 - jnp.maximum(cf[:, ilay - 1], cf[:, ilay])) \
            / (1.0 - jnp.minimum(cf[:, ilay - 1], 1.0 - thr))
        wts = jnp.where(active, clr - tmp, 0.0)
        denom = qc[:, ilay] + qi[:, ilay]
        phi = qc[:, ilay] / jnp.where(denom == 0.0, 1.0, denom)
        outs["T_mid_at_cldtop"] += tmid[:, ilay] * wts
        outs["p_mid_at_cldtop"] += jnp.asarray(pmid)[:, ilay] * wts
        outs["cldfrac_ice_at_cldtop"] += (1.0 - phi) * wts
        outs["cldfrac_liq_at_cldtop"] += phi * wts
        cdnc = (jnp.asarray(nc)[:, ilay] * jnp.asarray(p_del)[:, ilay]
                / jnp.asarray(z_del)[:, ilay] / GRAVIT
                / jnp.where(cf[:, ilay] == 0.0, 1.0, cf[:, ilay]))
        outs["cdnc_at_cldtop"] += cdnc * phi * wts
        outs["eff_radius_qc_at_cldtop"] += jnp.asarray(rel)[:, ilay] * phi * wts
        outs["eff_radius_qi_at_cldtop"] += jnp.asarray(rei)[:, ilay] \
            * (1.0 - phi) * wts
        clr = jnp.where(active, tmp, clr)
    outs["cldfrac_tot_at_cldtop"] = 1.0 - clr
    return outs


def rrtmgp_sw(kd_sw, play, plev, tlay, vmr, sfc_alb_dir, sfc_alb_dif, mu0,
              aerosol, clouds_gpt, tsi_scaling,
              extra_clnclrsky_diag=False, extra_clnsky_diag=False):
    """eamxx rrtmgp_sw. aerosol: by-band 2-stream dict; clouds_gpt:
    subsampled by-gpoint 2-stream dict. Returns dict of flux sets."""
    mu0 = jnp.asarray(mu0)
    day = mu0 > 0.0
    mu0_safe = jnp.where(day, mu0, 1.0)
    ncol, nlay = jnp.asarray(play).shape
    nlev = nlay + 1
    zero2 = jnp.zeros((ncol, nlev))

    tlay_lim = jnp.clip(jnp.asarray(tlay), kd_sw["temp_ref_min"],
                        kd_sw["temp_ref_max"])
    optics, toa_flux, _ = go.gas_optics_sw(kd_sw, play, plev, tlay_lim, vmr)
    toa_flux = tsi_scaling * toa_flux

    b2g, g2b = kd_sw["band2gpt"], kd_sw["gpt2band"]
    alb_dir_T = jnp.asarray(sfc_alb_dir).T  # (nbnd, ncol)
    alb_dif_T = jnp.asarray(sfc_alb_dif).T

    def solve(optics_):
        return rte.rte_sw(optics_, mu0_safe, toa_flux, alb_dir_T, alb_dif_T,
                          g2b, b2g)

    def mask(f):
        d = day.reshape((-1,) + (1,) * (f.ndim - 1))
        return jnp.where(d, f, 0.0)

    out = {}
    if extra_clnclrsky_diag:
        r = solve(optics)
        out["clnclrsky"] = {k: mask(r[k]) for k in
                            ("flux_up", "flux_dn", "flux_dn_dir")}
    else:
        out["clnclrsky"] = {k: zero2 for k in
                            ("flux_up", "flux_dn", "flux_dn_dir")}

    aerosol_ds = op.delta_scale_2str(aerosol)
    optics_aer = op.increment_2stream_by_2stream(optics, aerosol_ds,
                                                 gpt2band=g2b)
    r = solve(optics_aer)
    out["clrsky"] = {k: mask(r[k]) for k in
                     ("flux_up", "flux_dn", "flux_dn_dir")}

    clouds_ds = op.delta_scale_2str(clouds_gpt)
    optics_all = op.increment_2stream_by_2stream(optics_aer, clouds_ds)
    r = solve(optics_all)
    out["allsky"] = {k: mask(r[k]) for k in
                     ("flux_up", "flux_dn", "flux_dn_dir")}
    out["allsky"]["bnd_flux_up"] = mask(r["bnd_flux_up"])
    out["allsky"]["bnd_flux_dn"] = mask(r["bnd_flux_dn"])
    out["allsky"]["bnd_flux_dn_dir"] = mask(r["bnd_flux_dn_dir"])

    if extra_clnsky_diag:
        optics_cln = op.increment_2stream_by_2stream(optics, clouds_ds)
        r = solve(optics_cln)
        out["clnsky"] = {k: mask(r[k]) for k in
                         ("flux_up", "flux_dn", "flux_dn_dir")}
    else:
        out["clnsky"] = {k: zero2 for k in
                         ("flux_up", "flux_dn", "flux_dn_dir")}
    return out


def rrtmgp_lw(kd_lw, play, plev, tlay, tlev, vmr, aerosol, clouds_gpt,
              extra_clnclrsky_diag=False, extra_clnsky_diag=False):
    """eamxx rrtmgp_lw (surface emissivity 0.98, 1 quadrature angle)."""
    play = jnp.asarray(play)
    ncol, nlay = play.shape
    nband = kd_lw["nband"]
    zero2 = jnp.zeros((ncol, nlay + 1))

    tlay_lim = jnp.clip(jnp.asarray(tlay), kd_lw["temp_ref_min"],
                        kd_lw["temp_ref_max"])
    tlev_lim = jnp.clip(jnp.asarray(tlev), kd_lw["temp_ref_min"],
                        kd_lw["temp_ref_max"])
    t_sfc = tlev_lim[:, nlay]
    emis_sfc_T = jnp.full((nband, ncol), 0.98)

    optics, sources, _ = go.gas_optics_lw(kd_lw, play, plev, tlay_lim,
                                          t_sfc, vmr, tlev=tlev_lim)
    b2g, g2b = kd_lw["band2gpt"], kd_lw["gpt2band"]

    def solve(optics_):
        return rte.rte_lw(optics_, sources, emis_sfc_T, g2b, b2g)

    out = {}
    if extra_clnclrsky_diag:
        r = solve(optics)
        out["clnclrsky"] = {k: r[k] for k in ("flux_up", "flux_dn")}
    else:
        out["clnclrsky"] = {"flux_up": zero2, "flux_dn": zero2}

    optics_aer = op.increment_1scalar_by_1scalar(optics, aerosol,
                                                 gpt2band=g2b)
    r = solve(optics_aer)
    out["clrsky"] = {k: r[k] for k in ("flux_up", "flux_dn")}

    optics_all = op.increment_1scalar_by_1scalar(optics_aer, clouds_gpt)
    r = solve(optics_all)
    out["allsky"] = {k: r[k] for k in ("flux_up", "flux_dn",
                                       "bnd_flux_up", "bnd_flux_dn")}

    if extra_clnsky_diag:
        optics_cln = op.increment_1scalar_by_1scalar(optics, clouds_gpt)
        r = solve(optics_cln)
        out["clnsky"] = {k: r[k] for k in ("flux_up", "flux_dn")}
    else:
        out["clnsky"] = {"flux_up": zero2, "flux_dn": zero2}
    return out


def rrtmgp_main(kd_sw, kd_lw, co_sw, co_lw,
                play, tlay, plev, tlev, vmr,
                sfc_alb_dir, sfc_alb_dif, mu0,
                lwp, iwp, rel, rei, cldfrac,
                aer_tau_sw, aer_ssa_sw, aer_g_sw, aer_tau_lw,
                tsi_scaling,
                extra_clnclrsky_diag=False, extra_clnsky_diag=False):
    """eamxx rrtmgp_main. lwp/iwp in g/m2; aerosol arrays
    (ncol, nlay, nband). Returns (sw dict, lw dict, cld_tau dicts)."""
    aerosol_sw = {"tau": jnp.asarray(aer_tau_sw),
                  "ssa": jnp.asarray(aer_ssa_sw),
                  "g": jnp.asarray(aer_g_sw)}
    aerosol_lw = {"tau": jnp.asarray(aer_tau_lw)}

    clouds_sw_bnd = co.get_cloud_optics(co_sw, lwp, iwp, rel, rei,
                                        two_stream=True)
    clouds_lw_bnd = co.get_cloud_optics(co_lw, lwp, iwp, rel, rei,
                                        two_stream=False)

    clouds_sw_gpt = mcica.get_subsampled_clouds(
        clouds_sw_bnd, cldfrac, play, kd_sw["gpt2band"], kd_sw["ngpt"],
        two_stream=True)
    clouds_lw_gpt = mcica.get_subsampled_clouds(
        clouds_lw_bnd, cldfrac, play, kd_lw["gpt2band"], kd_lw["ngpt"],
        two_stream=False)

    sw = rrtmgp_sw(kd_sw, play, plev, tlay, vmr, sfc_alb_dir, sfc_alb_dif,
                   mu0, aerosol_sw, clouds_sw_gpt, tsi_scaling,
                   extra_clnclrsky_diag, extra_clnsky_diag)
    lw = rrtmgp_lw(kd_lw, play, plev, tlay, tlev, vmr, aerosol_lw,
                   clouds_lw_gpt, extra_clnclrsky_diag, extra_clnsky_diag)

    cld_tau = {"sw_bnd": clouds_sw_bnd["tau"], "lw_bnd": clouds_lw_bnd["tau"],
               "sw_gpt": clouds_sw_gpt["tau"], "lw_gpt": clouds_lw_gpt["tau"]}
    return sw, lw, cld_tau
