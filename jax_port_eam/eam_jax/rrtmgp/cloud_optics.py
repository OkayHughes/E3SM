"""EAM cloud optics: Conley gamma-distribution liquid + Mitchell ice
(the EAMv3 defaults liqcldoptics='gammadist', icecldoptics='mitchell').

FRESH PORT (EAMxx uses the RRTMGP cloud-optics LUT instead). Sources,
ported verbatim:
  components/eam/src/physics/rrtmgp/cloud_rad_props.F90
    (cloud_rad_props_init table reads incl. the in-code
     ext_sw_liq/abs_lw_liq /= 0.9970449e3 unit fix,
     gammadist_liq_optics_sw/lw [lamc > 0 gate], gam_liquid_sw/lw
     [clwptn < 1e-80 gate], get_mu_lambda_weights,
     mitchell_ice_optics_sw/lw [iciwpth < 1e-80 .or. dei == 0 gate])
  components/eam/src/control/interpolate_data.F90 (lininterp_init with
    extrap_method_bndry — boundary-copy weights, supports the
    DECREASING lambda axis of the table — lininterp1d, lininterp2d1d
    bilinear product weights)
  components/eam/src/physics/rrtmgp/cam_optics.F90
    (get_cloud_optics_sw/lw assembly: ice+liq sums, snow combined via
     combine_properties max-fraction weighting, tau/ssa/asm extraction
     with >0 guards)

Table axes: g_mu increasing; the per-mu lambda row (interpolated to
the target mu first) is DECREASING along the lambda index; g_d_eff
increasing. The lininterp weight construction reproduces the Fortran
interval search exactly, including the <=/>= boundary conventions.

The liquid tables are loaded with the netCDF dims reversed into
Fortran index order: g_lambda (nmu, nlambda), ext_sw_liq
(nmu, nlambda, nswbands), ext_sw_ice (n_g_d, nswbands).
"""

import numpy as np
import jax.numpy as jnp


def load_liq_optics(filename):
    """cloud_rad_props_init, liquid file (incl. the /0.9970449e3 unit
    conversion applied to ext_sw_liq and abs_lw_liq)."""
    import netCDF4
    ds = netCDF4.Dataset(filename)
    v = ds.variables

    def rev(name):
        a = np.asarray(v[name][:], dtype=np.float64)
        return np.ascontiguousarray(
            np.transpose(a, tuple(range(a.ndim))[::-1]))

    liq = {
        "g_mu": np.asarray(v["mu"][:], dtype=np.float64),
        "g_lambda": rev("lambda"),                # (nmu, nlambda)
        "ext_sw_liq": rev("k_ext_sw") / 0.9970449e3,
        "ssa_sw_liq": rev("ssa_sw"),
        "asm_sw_liq": rev("asm_sw"),
        "abs_lw_liq": rev("k_abs_lw") / 0.9970449e3,
    }
    ds.close()
    liq["nmu"] = liq["g_mu"].shape[0]
    liq["nlambda"] = liq["g_lambda"].shape[1]
    return liq


def load_ice_optics(filename):
    """cloud_rad_props_init, ice file."""
    import netCDF4
    ds = netCDF4.Dataset(filename)
    v = ds.variables

    def rev(name):
        a = np.asarray(v[name][:], dtype=np.float64)
        return np.ascontiguousarray(
            np.transpose(a, tuple(range(a.ndim))[::-1]))

    ice = {
        "g_d_eff": np.asarray(v["d_eff"][:], dtype=np.float64),
        "ext_sw_ice": rev("sw_ext"),              # (n_g_d, nswbands)
        "ssa_sw_ice": rev("sw_ssa"),
        "asm_sw_ice": rev("sw_asm"),
        "abs_lw_ice": rev("lw_abs"),
    }
    ds.close()
    return ice


def _lininterp_weights_inc(yin, yout):
    """lininterp_init (extrap_method_bndry) for an INCREASING 1-D grid
    yin; yout any shape. Returns (jjm, jjp, wgts, wgtn), 0-based."""
    yin = jnp.asarray(yin)
    yout = jnp.asarray(yout)
    n = yin.shape[0]
    # interior: yin[jj] < yout <= yin[jj+1]  ->  jj = ss_left - 1
    k = jnp.searchsorted(yin, yout, side="left")
    jj = jnp.clip(k - 1, 0, n - 2)
    dy = yin[jj + 1] - yin[jj]
    wgts = (yin[jj + 1] - yout) / dy
    wgtn = (yout - yin[jj]) / dy
    jjm, jjp = jj, jj + 1
    # boundaries: copy end values
    lo = yout <= yin[0]
    hi = yout > yin[n - 1]
    jjm = jnp.where(lo, 0, jnp.where(hi, n - 1, jjm))
    jjp = jnp.where(lo, 0, jnp.where(hi, n - 1, jjp))
    wgts = jnp.where(lo | hi, 1.0, wgts)
    wgtn = jnp.where(lo | hi, 0.0, wgtn)
    return jjm, jjp, wgts, wgtn


def _lininterp_weights_dec(yin, yout):
    """lininterp_init (extrap_method_bndry) for a DECREASING grid:
    yin (..., n) may vary per output point (the lambda row after mu
    interpolation). Interior: yout <= yin[jj] .and. yout > yin[jj+1]."""
    yin = jnp.asarray(yin)
    yout = jnp.asarray(yout)
    n = yin.shape[-1]
    rev = yin[..., ::-1]  # increasing
    # left insertion point of yout in rev = count of rev values < yout
    # (rev may be a different grid per output point)
    m = jnp.sum(rev < yout[..., None], axis=-1)
    # map back: rev index m corresponds to original index n-1-m
    jj = jnp.clip(n - 1 - m, 0, n - 2)
    if yin.ndim > 1:
        y_j = jnp.take_along_axis(yin, jj[..., None], axis=-1)[..., 0]
        y_j1 = jnp.take_along_axis(yin, (jj + 1)[..., None], axis=-1)[..., 0]
        y_first = yin[..., 0]
        y_last = yin[..., n - 1]
    else:
        y_j = yin[jj]
        y_j1 = yin[jj + 1]
        y_first = yin[0]
        y_last = yin[n - 1]
    dy = y_j1 - y_j
    wgts = (y_j1 - yout) / dy
    wgtn = (yout - y_j) / dy
    jjm, jjp = jj, jj + 1
    lo = yout > y_first          # beyond the (large-value) start
    hi = yout <= y_last          # beyond the (small-value) end
    jjm = jnp.where(lo, 0, jnp.where(hi, n - 1, jjm))
    jjp = jnp.where(lo, 0, jnp.where(hi, n - 1, jjp))
    wgts = jnp.where(lo | hi, 1.0, wgts)
    wgtn = jnp.where(lo | hi, 0.0, wgtn)
    return jjm, jjp, wgts, wgtn


def _gam_liquid(liq, clwptn, lamc, pgam, tables):
    """get_mu_lambda_weights + lininterp2d1d over each table in
    `tables` (list of (nmu, nlambda, nbnd) arrays). clwptn/lamc/pgam
    (ncol, nlev). Returns per-table (ncol, nlev, nbnd) interpolants
    (WITHOUT the clwptn factor)."""
    g_mu = jnp.asarray(liq["g_mu"])
    g_lambda = jnp.asarray(liq["g_lambda"])
    # mu weights (increasing grid)
    im, ip, wm, wn = _lininterp_weights_inc(g_mu, pgam)
    # interpolate the lambda table rows to the target mu
    lam_interp = (g_lambda[im, :] * wm[..., None]
                  + g_lambda[ip, :] * wn[..., None])  # (ncol,nlev,nlambda)
    # lambda weights (decreasing, per-point grid)
    jm, jp, ws, wnn = _lininterp_weights_dec(lam_interp, lamc)
    out = []
    for tab in tables:
        tab = jnp.asarray(tab)
        # lininterp2d1d: arrin(iim,jjm)*wgtw*wgts + arrin(iip,jjm)*wgte*wgts
        #              + arrin(iim,jjp)*wgtw*wgtn + arrin(iip,jjp)*wgte*wgtn
        val = (tab[im, jm, :] * (wm * ws)[..., None]
               + tab[ip, jm, :] * (wn * ws)[..., None]
               + tab[im, jp, :] * (wm * wnn)[..., None]
               + tab[ip, jp, :] * (wn * wnn)[..., None])
        out.append(val)
    return out


def gammadist_liq_optics_sw(liq, iclwp, lamc, pgam):
    """gammadist_liq_optics_sw + gam_liquid_sw. Returns
    (tau, tau_w, tau_w_g, tau_w_f), each (ncol, nlev, nbnd)."""
    ext, ssa, asy = _gam_liquid(liq, iclwp, lamc, pgam,
                                [liq["ext_sw_liq"], liq["ssa_sw_liq"],
                                 liq["asm_sw_liq"]])
    active = ((jnp.asarray(lamc) > 0.0)
              & (jnp.asarray(iclwp) >= 1.0e-80))[..., None]
    tau = jnp.where(active, jnp.asarray(iclwp)[..., None] * ext, 0.0)
    tau_w = tau * jnp.where(active, ssa, 0.0)
    tau_w_g = tau_w * jnp.where(active, asy, 0.0)
    tau_w_f = tau_w_g * jnp.where(active, asy, 0.0)
    return tau, tau_w, tau_w_g, tau_w_f


def gammadist_liq_optics_lw(liq, iclwp, lamc, pgam):
    """gammadist_liq_optics_lw + gam_liquid_lw -> abs_od."""
    (absl,) = _gam_liquid(liq, iclwp, lamc, pgam, [liq["abs_lw_liq"]])
    active = ((jnp.asarray(lamc) > 0.0)
              & (jnp.asarray(iclwp) >= 1.0e-80))[..., None]
    return jnp.where(active, jnp.asarray(iclwp)[..., None] * absl, 0.0)


def mitchell_ice_optics_sw(ice, iciwp, dei):
    """mitchell_ice_optics_sw: (tau, tau_w, tau_w_g, tau_w_f)."""
    im, ip, wm, wn = _lininterp_weights_inc(jnp.asarray(ice["g_d_eff"]),
                                            dei)
    def interp(tab):
        tab = jnp.asarray(tab)
        return tab[im, :] * wm[..., None] + tab[ip, :] * wn[..., None]

    ext = interp(ice["ext_sw_ice"])
    ssa = interp(ice["ssa_sw_ice"])
    asy = interp(ice["asm_sw_ice"])
    active = ((jnp.asarray(iciwp) >= 1.0e-80)
              & (jnp.asarray(dei) != 0.0))[..., None]
    tau = jnp.where(active, jnp.asarray(iciwp)[..., None] * ext, 0.0)
    tau_w = tau * jnp.where(active, ssa, 0.0)
    tau_w_g = tau_w * jnp.where(active, asy, 0.0)
    tau_w_f = tau_w_g * jnp.where(active, asy, 0.0)
    return tau, tau_w, tau_w_g, tau_w_f


def mitchell_ice_optics_lw(ice, iciwp, dei):
    """mitchell_ice_optics_lw -> abs_od."""
    im, ip, wm, wn = _lininterp_weights_inc(jnp.asarray(ice["g_d_eff"]),
                                            dei)
    tab = jnp.asarray(ice["abs_lw_ice"])
    absr = tab[im, :] * wm[..., None] + tab[ip, :] * wn[..., None]
    active = ((jnp.asarray(iciwp) >= 1.0e-80)
              & (jnp.asarray(dei) != 0.0))[..., None]
    return jnp.where(active, jnp.asarray(iciwp)[..., None] * absr, 0.0)


def combine_properties(fraction1, property1, fraction2, property2):
    """cam_optics combine_properties (fractions (ncol,nlev),
    properties (ncol,nlev,nbnd))."""
    f1 = jnp.asarray(fraction1)[..., None]
    f2 = jnp.asarray(fraction2)[..., None]
    combined = jnp.maximum(f1, f2)
    return jnp.where(combined > 0.0,
                     (f1 * property1 + f2 * property2)
                     / jnp.where(combined == 0.0, 1.0, combined),
                     0.0)


def get_cloud_optics_sw(liq, ice, do_snow, cld, cldfsnow, iclwp, iciwp,
                        icswp, lambdac, mu, dei, des):
    """cam_optics get_cloud_optics_sw (gammadist + mitchell).
    Returns dict with tau/ssa/asm (ncol, nlev, nbnd) in the ORIGINAL
    (RRTMG-ordered) bands plus liq_tau/ice_tau/snw_tau diagnostics."""
    ice_tau, ice_tau_ssa, ice_tau_ssa_g, _ = mitchell_ice_optics_sw(
        ice, iciwp, dei)
    liq_tau, liq_tau_ssa, liq_tau_ssa_g, _ = gammadist_liq_optics_sw(
        liq, iclwp, lambdac, mu)

    if do_snow:
        snow_tau, snow_tau_ssa, snow_tau_ssa_g, _ = mitchell_ice_optics_sw(
            ice, icswp, des)
    else:
        snow_tau = jnp.zeros_like(ice_tau)
        snow_tau_ssa = jnp.zeros_like(ice_tau)
        snow_tau_ssa_g = jnp.zeros_like(ice_tau)

    cld_tau = ice_tau + liq_tau
    cld_tau_ssa = ice_tau_ssa + liq_tau_ssa
    cld_tau_ssa_g = ice_tau_ssa_g + liq_tau_ssa_g
    if do_snow:
        combined_tau = combine_properties(cld, cld_tau, cldfsnow, snow_tau)
        combined_tau_ssa = combine_properties(cld, cld_tau_ssa, cldfsnow,
                                              snow_tau_ssa)
        combined_tau_ssa_g = combine_properties(cld, cld_tau_ssa_g,
                                                cldfsnow, snow_tau_ssa_g)
    else:
        combined_tau = cld_tau
        combined_tau_ssa = cld_tau_ssa
        combined_tau_ssa_g = cld_tau_ssa_g

    ssa = jnp.where(combined_tau > 0.0,
                    combined_tau_ssa
                    / jnp.where(combined_tau == 0.0, 1.0, combined_tau),
                    1.0)
    asm_ = jnp.where(combined_tau_ssa > 0.0,
                     combined_tau_ssa_g
                     / jnp.where(combined_tau_ssa == 0.0, 1.0,
                                 combined_tau_ssa),
                     0.0)
    return {"tau": combined_tau, "ssa": ssa, "asm": asm_,
            "liq_tau": liq_tau, "ice_tau": ice_tau, "snw_tau": snow_tau}


def get_cloud_optics_lw(liq, ice, do_snow, cld, cldfsnow, iclwp, iciwp,
                        icswp, lambdac, mu, dei, des):
    """cam_optics get_cloud_optics_lw."""
    ice_tau = mitchell_ice_optics_lw(ice, iciwp, dei)
    liq_tau = gammadist_liq_optics_lw(liq, iclwp, lambdac, mu)
    cld_tau = liq_tau + ice_tau
    if do_snow:
        snow_tau = mitchell_ice_optics_lw(ice, icswp, des)
        combined_tau = combine_properties(cld, cld_tau, cldfsnow, snow_tau)
    else:
        snow_tau = jnp.zeros_like(ice_tau)
        combined_tau = cld_tau
    return {"tau": combined_tau, "liq_tau": liq_tau, "ice_tau": ice_tau,
            "snw_tau": snow_tau}
