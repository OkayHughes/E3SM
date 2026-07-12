"""Cloud optics: LUT loading and liquid/ice cloud optical properties.

Sources: cpp/extensions/cloud_optics/mo_cloud_optics.h (CloudOpticsK:
load [LUT variant], cloud_optics, compute_all_from_table, combine) and
cpp/examples/all-sky/mo_load_cloud_coefficients.cpp (load_cld_lutcoeff
netCDF names). EAMxx uses the LUT path with ice roughness 2
(get_cloud_optics_sw/lw set set_ice_roughness(2)) and limits rel/rei
into the table bounds before the call.
"""

import numpy as np
import jax.numpy as jnp

# conv::epsilon(tau) = std::numeric_limits<double>::epsilon()
_MACH_EPS = 2.220446049250313e-16


def load_cloud_optics(filename):
    """load_cld_lutcoeff: returns a dict of LUT arrays (index order as
    the C++ views: lut_extliq (nsize_liq, nbnd), lut_extice
    (nsize_ice, nbnd, nrghice))."""
    import netCDF4
    ds = netCDF4.Dataset(filename)
    v = ds.variables

    def rev(name):
        a = np.asarray(v[name][:], dtype=np.float64)
        return np.ascontiguousarray(np.transpose(a, tuple(range(a.ndim))[::-1]))

    co = {
        "band_lims_wvn": rev("bnd_limits_wavenumber"),
        "radliq_lwr": float(np.asarray(v["radliq_lwr"][:])),
        "radliq_upr": float(np.asarray(v["radliq_upr"][:])),
        "radice_lwr": float(np.asarray(v["radice_lwr"][:])),
        "radice_upr": float(np.asarray(v["radice_upr"][:])),
        "lut_extliq": rev("lut_extliq"),
        "lut_ssaliq": rev("lut_ssaliq"),
        "lut_asyliq": rev("lut_asyliq"),
        "lut_extice": rev("lut_extice"),
        "lut_ssaice": rev("lut_ssaice"),
        "lut_asyice": rev("lut_asyice"),
    }
    ds.close()
    co["liq_nsteps"] = co["lut_extliq"].shape[0]
    co["ice_nsteps"] = co["lut_extice"].shape[0]
    co["liq_step_size"] = (co["radliq_upr"] - co["radliq_lwr"]) \
        / (co["liq_nsteps"] - 1.0)
    co["ice_step_size"] = (co["radice_upr"] - co["radice_lwr"]) \
        / (co["ice_nsteps"] - 1.0)
    co["nband"] = co["lut_extliq"].shape[1]
    return co


def _compute_all_from_table(mask, wp, re, nsteps, step_size, offset,
                            tau_table, ssa_table, asy_table):
    """compute_all_from_table -> (tau, taussa, taussag), each
    (ncol, nlay, nbnd); zero where mask is false."""
    re = jnp.asarray(re)
    wp = jnp.asarray(wp)
    tau_table = jnp.asarray(tau_table)
    ssa_table = jnp.asarray(ssa_table)
    asy_table = jnp.asarray(asy_table)

    index = jnp.minimum(jnp.floor((re - offset) / step_size) + 1,
                        nsteps - 1.0).astype(jnp.int32) - 1
    index = jnp.maximum(index, 0)  # guard (C++ relies on bounded re)
    fint = (re - offset) / step_size - index.astype(re.dtype)
    t = wp[..., None] * (tau_table[index] + fint[..., None]
                         * (tau_table[index + 1] - tau_table[index]))
    ts = t * (ssa_table[index] + fint[..., None]
              * (ssa_table[index + 1] - ssa_table[index]))
    tsg = ts * (asy_table[index] + fint[..., None]
                * (asy_table[index + 1] - asy_table[index]))
    m = jnp.asarray(mask)[..., None]
    return (jnp.where(m, t, 0.0), jnp.where(m, ts, 0.0),
            jnp.where(m, tsg, 0.0))


def cloud_optics(co, clwp, ciwp, reliq, reice, two_stream: bool,
                 icergh: int = 2):
    """CloudOpticsK::cloud_optics (LUT path). clwp/ciwp in g/m2.

    Returns a 2-stream dict {tau, ssa, g} (SW) or {tau} (LW), by band."""
    clwp = jnp.asarray(clwp)
    ciwp = jnp.asarray(ciwp)
    liqmsk = clwp > 0.0
    icemsk = ciwp > 0.0

    ltau, ltaussa, ltaussag = _compute_all_from_table(
        liqmsk, clwp, reliq, co["liq_nsteps"], co["liq_step_size"],
        co["radliq_lwr"], co["lut_extliq"], co["lut_ssaliq"],
        co["lut_asyliq"])
    itau, itaussa, itaussag = _compute_all_from_table(
        icemsk, ciwp, reice, co["ice_nsteps"], co["ice_step_size"],
        co["radice_lwr"], co["lut_extice"][:, :, icergh - 1],
        co["lut_ssaice"][:, :, icergh - 1],
        co["lut_asyice"][:, :, icergh - 1])

    if two_stream:
        tau = ltau + itau
        taussa = ltaussa + itaussa
        g = (ltaussag + itaussag) / jnp.maximum(_MACH_EPS, taussa)
        ssa = taussa / jnp.maximum(_MACH_EPS, tau)
        return {"tau": tau, "ssa": ssa, "g": g}
    # absorption optical depth = (1-ssa)*tau = tau - taussa
    return {"tau": (ltau - ltaussa) + (itau - itaussa)}


def get_cloud_optics(co, lwp, iwp, rel, rei, two_stream: bool):
    """eamxx get_cloud_optics_sw/lw: limit radii into the LUT bounds,
    ice roughness 2, then cloud_optics."""
    rel_lim = jnp.clip(jnp.asarray(rel), co["radliq_lwr"], co["radliq_upr"])
    rei_lim = jnp.clip(jnp.asarray(rei), co["radice_lwr"], co["radice_upr"])
    return cloud_optics(co, lwp, iwp, rel_lim, rei_lim, two_stream,
                        icergh=2)
