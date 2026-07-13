"""RRTMGP gas optics: interpolation, absorption/Rayleigh optical depths,
Planck sources, solar source.

COPIED UNCHANGED from ../jax_port/scream_jax/rrtmgp/gas_optics.py and
re-validated line-by-line against the EAM FORTRAN sources
(components/eam/src/physics/rrtmgp/external/rrtmgp/kernels/
mo_gas_optics_kernels.F90: interpolation, gas_optical_depths_major/
minor, compute_tau_absorption, compute_tau_rayleigh,
compute_Planck_source, combine_and_reorder_2str, interpolate1D;
external/rrtmgp/mo_gas_optics_rrtmgp.F90: compute_gas_taus, source,
get_col_dry [no-latitude branch: g0 = grav, as EAM calls it];
external/rrtmgp/mo_rrtmgp_constants.F90 — EAM never calls
init_constants, so grav/m_dry/m_h2o/avogad keep their defaults, which
equal the constants below). The Fortran minor-gas kernel restricts each
atmosphere to the contiguous layer_limits range from minloc/maxloc of
play under the tropo mask; the tropo-mask where() here is equivalent
because log(play) crosses press_ref_trop_log monotonically in these
columns.

Array conventions: profile fields are (ncol, nlay); interpolation
outputs keep the C++ leading axes, e.g. fmajor(2,2,2,nflav,ncol,nlay);
tables keep the C++ index order (kmajor(gpt,eta,press,temp)). Optical
properties are returned as (ncol, nlay, ngpt). The C++ accumulates
minor-gas contributions serially per g-point; here they are summed with
a segment-sum, which can reorder FP additions (roundoff-level, non-BFB).
"""

import jax.numpy as jnp
from jax import ops as jops

GRAV = 9.80665
M_DRY = 0.028964
M_H2O = 0.018016
AVOGAD = 6.02214076e23

_TINY = 2.2250738585072014e-308  # numeric_limits<double>::min()


def interpolation(kd, play, tlay, col_gas):
    """mo_gas_optics_kernels interpolation(). Returns a dict with
    jtemp, jpress, tropo (ncol,nlay), jeta/col_mix (2,nflav,ncol,nlay),
    fminor (2,2,nflav,ncol,nlay), fmajor (2,2,2,nflav,ncol,nlay)."""
    play = jnp.asarray(play)
    tlay = jnp.asarray(tlay)
    col_gas = jnp.asarray(col_gas)
    flavor = jnp.asarray(kd["flavor"])
    vmr_ref = jnp.asarray(kd["vmr_ref"])
    temp_ref = jnp.asarray(kd["temp_ref"])
    press_ref_log = jnp.asarray(kd["press_ref_log"])
    ntemp, npres, neta = kd["ntemp"], kd["npres"], kd["neta"]

    jtemp = ((tlay - (kd["temp_ref_min"] - kd["temp_ref_delta"]))
             / kd["temp_ref_delta"]).astype(jnp.int32)
    jtemp = jnp.clip(jtemp, 1, ntemp - 1) - 1
    ftemp = (tlay - temp_ref[jtemp]) / kd["temp_ref_delta"]

    locpress = 1.0 + (jnp.log(play) - press_ref_log[0]) / kd["press_ref_log_delta"]
    jpress = jnp.clip(locpress.astype(jnp.int32), 1, npres - 1)
    fpress = locpress - jpress.astype(play.dtype)

    tropo = jnp.log(play) > kd["press_ref_trop_log"]
    itropo = jnp.where(tropo, 0, 1)  # 0 = lower atmosphere

    # (2, nflav, ncol, nlay) with itemp on the leading axis
    ig1 = flavor[0][None, :, None, None]  # gas indices (may be -1)
    ig2 = flavor[1][None, :, None, None]
    itemp = jnp.arange(2)[:, None, None, None]
    jt = jtemp[None, None] + itemp  # (2, 1->nflav, ncol, nlay) broadcast

    it_b = itropo[None, None]
    ratio_eta_half = (vmr_ref[it_b, ig1 + 1, jt]
                      / vmr_ref[it_b, ig2 + 1, jt])
    # col_gas gathered at the two flavor gases: (1, nflav, ncol, nlay)
    cg1 = jnp.moveaxis(jnp.take(col_gas, flavor[0] + 1, axis=-1), -1, 0)[None]
    cg2 = jnp.moveaxis(jnp.take(col_gas, flavor[1] + 1, axis=-1), -1, 0)[None]

    col_mix = cg1 + ratio_eta_half * cg2
    eta = jnp.where(col_mix > 2.0 * _TINY, cg1 / jnp.where(col_mix == 0.0, 1.0, col_mix), 0.5)
    loceta = eta * (neta - 1.0)
    jeta = jnp.minimum(loceta.astype(jnp.int32) + 1, neta - 1) - 1
    feta = jnp.mod(loceta, 1.0)
    ftemp_term = ((2.0 - (itemp + 1)) + (2.0 * (itemp + 1) - 3.0) * ftemp[None, None])
    fminor = jnp.stack([(1.0 - feta) * ftemp_term, feta * ftemp_term])  # (2,2,nflav,ncol,nlay)
    fp = fpress[None, None]
    # fmajor(i_eta, i_press, itemp, ...):
    #   fmajor(0,0,..)=(1-fpress)*fminor(0,..), fmajor(1,0,..)=(1-fpress)*fminor(1,..)
    #   fmajor(0,1,..)=fpress*fminor(0,..),     fmajor(1,1,..)=fpress*fminor(1,..)
    fmajor = jnp.stack([
        jnp.stack([(1.0 - fp) * fminor[0], fp * fminor[0]]),
        jnp.stack([(1.0 - fp) * fminor[1], fp * fminor[1]]),
    ])  # (eta 2, press 2, itemp 2, nflav, ncol, nlay)

    return {"jtemp": jtemp, "jpress": jpress, "tropo": tropo,
            "jeta": jeta, "col_mix": col_mix, "fminor": fminor,
            "fmajor": fmajor}


def gas_optical_depths_major(kd, interp):
    """gas_optical_depths_major -> tau (ncol, nlay, ngpt)."""
    kmajor = jnp.asarray(kd["kmajor"])
    gpoint_flavor = jnp.asarray(kd["gpoint_flavor"])
    tropo = interp["tropo"]
    itropo = jnp.where(tropo, 0, 1)  # (ncol, nlay)
    ngpt = kd["ngpt"]

    igpt = jnp.arange(ngpt)[None, None, :]              # (1,1,ngpt)
    itropo_b = itropo[..., None]                        # (ncol,nlay,1)
    iflav = gpoint_flavor[itropo_b, igpt]               # (ncol,nlay,ngpt)
    jpress_loc = interp["jpress"][..., None] + itropo_b
    jtemp_loc = interp["jtemp"][..., None]

    ncol, nlay = tropo.shape
    ic = jnp.arange(ncol)[:, None, None]
    il = jnp.arange(nlay)[None, :, None]
    jeta0 = interp["jeta"][0, iflav, ic, il]
    jeta1 = interp["jeta"][1, iflav, ic, il]
    cm0 = interp["col_mix"][0, iflav, ic, il]
    cm1 = interp["col_mix"][1, iflav, ic, il]

    def fmaj(e, p, t):
        return interp["fmajor"][e, p, t, iflav, ic, il]

    tau = (cm0 * (
        fmaj(0, 0, 0) * kmajor[igpt, jeta0, jpress_loc - 1, jtemp_loc] +
        fmaj(1, 0, 0) * kmajor[igpt, jeta0 + 1, jpress_loc - 1, jtemp_loc] +
        fmaj(0, 1, 0) * kmajor[igpt, jeta0, jpress_loc, jtemp_loc] +
        fmaj(1, 1, 0) * kmajor[igpt, jeta0 + 1, jpress_loc, jtemp_loc])
        + cm1 * (
        fmaj(0, 0, 1) * kmajor[igpt, jeta1, jpress_loc - 1, jtemp_loc + 1] +
        fmaj(1, 0, 1) * kmajor[igpt, jeta1 + 1, jpress_loc - 1, jtemp_loc + 1] +
        fmaj(0, 1, 1) * kmajor[igpt, jeta1, jpress_loc, jtemp_loc + 1] +
        fmaj(1, 1, 1) * kmajor[igpt, jeta1 + 1, jpress_loc, jtemp_loc + 1]))
    return tau


def gas_optical_depths_minor(kd, interp, play, tlay, col_gas, atm):
    """gas_optical_depths_minor for atm in {'lower','upper'} ->
    tau contribution (ncol, nlay, ngpt)."""
    kminor = jnp.asarray(kd[f"kminor_{atm}"])
    idx_minor = jnp.asarray(kd[f"idx_minor_{atm}"])
    idx_minor_scaling = jnp.asarray(kd[f"idx_minor_scaling_{atm}"])
    swd = jnp.asarray(kd[f"minor_scales_with_density_{atm}"])
    sbc = jnp.asarray(kd[f"scale_by_complement_{atm}"])
    ent_gpt = jnp.asarray(kd[f"minor_entry_gpt_{atm}"])
    ent_imnr = jnp.asarray(kd[f"minor_entry_imnr_{atm}"])
    ent_k = jnp.asarray(kd[f"minor_entry_k_{atm}"])
    gpoint_flavor = jnp.asarray(kd["gpoint_flavor"])
    idx_tropo = 0 if atm == "lower" else 1
    idx_h2o = kd["idx_h2o"]
    ngpt = kd["ngpt"]

    play = jnp.asarray(play)
    tlay = jnp.asarray(tlay)
    col_gas = jnp.asarray(col_gas)
    nminor = idx_minor.shape[0]
    if nminor == 0:
        return jnp.zeros(play.shape + (ngpt,), dtype=play.dtype)

    # per-(col,lay,minor) scaling factor
    scaling = jnp.take(col_gas, idx_minor + 1, axis=-1)  # (ncol,nlay,nminor)
    dens = (0.01 * play / tlay)[..., None]
    scaled = scaling * dens
    vmr_fact = 1.0 / col_gas[..., 0:1]
    dry_fact = 1.0 / (1.0 + col_gas[..., idx_h2o + 1:idx_h2o + 2] * vmr_fact)
    cg_scl = jnp.take(col_gas, jnp.maximum(idx_minor_scaling, 0) + 1, axis=-1)
    fac = cg_scl * vmr_fact * dry_fact
    scaled_special = jnp.where(sbc[None, None, :], scaled * (1.0 - fac),
                               scaled * fac)
    scaled = jnp.where(idx_minor_scaling[None, None, :] > -1,
                       scaled_special, scaled)
    scaling = jnp.where(swd[None, None, :], scaled, scaling)

    # per-entry interpolation and segment-sum into g-points
    ncol, nlay = play.shape
    ic = jnp.arange(ncol)[:, None, None]
    il = jnp.arange(nlay)[None, :, None]
    iflav = gpoint_flavor[idx_tropo, ent_gpt][None, None, :]  # (1,1,nent)
    jt = interp["jtemp"][..., None]
    jeta0 = interp["jeta"][0, iflav, ic, il]
    jeta1 = interp["jeta"][1, iflav, ic, il]

    def fmin(e, t):
        return interp["fminor"][e, t, iflav, ic, il]

    kloc = (fmin(0, 0) * kminor[ent_k, jeta0, jt] +
            fmin(1, 0) * kminor[ent_k, jeta0 + 1, jt] +
            fmin(0, 1) * kminor[ent_k, jeta1, jt + 1] +
            fmin(1, 1) * kminor[ent_k, jeta1 + 1, jt + 1])
    tau_ent = kloc * scaling[..., ent_imnr]  # (ncol,nlay,nent)

    # active only in the matching atmosphere region (the contiguous
    # tropo/~tropo band, equivalent to the C++ layer_limits ranges)
    active = (interp["tropo"] if idx_tropo == 0 else ~interp["tropo"])
    tau_ent = jnp.where(active[..., None], tau_ent, 0.0)

    tau = jops.segment_sum(jnp.moveaxis(tau_ent, -1, 0), ent_gpt,
                           num_segments=ngpt)  # (ngpt,ncol,nlay)
    return jnp.moveaxis(tau, 0, -1)


def compute_tau_rayleigh(kd, interp, col_gas, col_dry):
    """compute_tau_rayleigh -> tau_rayleigh (ncol, nlay, ngpt)."""
    krayl = jnp.asarray(kd["krayl"])
    gpoint_flavor = jnp.asarray(kd["gpoint_flavor"])
    idx_h2o = kd["idx_h2o"]
    ngpt = kd["ngpt"]
    col_gas = jnp.asarray(col_gas)
    col_dry = jnp.asarray(col_dry)

    tropo = interp["tropo"]
    itropo = jnp.where(tropo, 0, 1)
    ncol, nlay = tropo.shape
    igpt = jnp.arange(ngpt)[None, None, :]
    itropo_b = itropo[..., None]
    iflav = gpoint_flavor[itropo_b, igpt]
    ic = jnp.arange(ncol)[:, None, None]
    il = jnp.arange(nlay)[None, :, None]
    jt = interp["jtemp"][..., None]
    jeta0 = interp["jeta"][0, iflav, ic, il]
    jeta1 = interp["jeta"][1, iflav, ic, il]

    def fmin(e, t):
        return interp["fminor"][e, t, iflav, ic, il]

    k = (fmin(0, 0) * krayl[igpt, jeta0, jt, itropo_b] +
         fmin(1, 0) * krayl[igpt, jeta0 + 1, jt, itropo_b] +
         fmin(0, 1) * krayl[igpt, jeta1, jt + 1, itropo_b] +
         fmin(1, 1) * krayl[igpt, jeta1 + 1, jt + 1, itropo_b])
    return k * (col_gas[..., idx_h2o + 1] + col_dry)[..., None]


def get_col_dry(vmr_h2o, plev):
    """GasOpticsRRTMGPK::get_col_dry (constant gravity path)."""
    vmr_h2o = jnp.asarray(vmr_h2o)
    plev = jnp.asarray(plev)
    delta_plev = jnp.abs(plev[..., :-1] - plev[..., 1:])
    fact = 1.0 / (1.0 + vmr_h2o)
    m_air = (M_DRY + M_H2O * vmr_h2o) * fact
    return 10.0 * delta_plev * AVOGAD * fact / (1000.0 * m_air * 100.0 * GRAV)


def compute_gas_taus(kd, play, plev, tlay, vmr, col_dry=None):
    """GasOpticsRRTMGPK::compute_gas_taus.

    vmr: (ncol, nlay, ngas) volume mixing ratios ordered like
    kd['gas_names']. Returns (optics dict, col_gas, interp dict).
    For SW (krayl present): optics has tau/ssa/g; for LW: tau only.
    """
    play = jnp.asarray(play)
    vmr = jnp.asarray(vmr)
    if col_dry is None:
        col_dry = get_col_dry(vmr[..., kd["idx_h2o"]], plev)
    col_gas = jnp.concatenate([col_dry[..., None], vmr * col_dry[..., None]],
                              axis=-1)

    interp = interpolation(kd, play, tlay, col_gas)
    tau = gas_optical_depths_major(kd, interp)
    tau = tau + gas_optical_depths_minor(kd, interp, play, tlay, col_gas,
                                         "lower")
    tau = tau + gas_optical_depths_minor(kd, interp, play, tlay, col_gas,
                                         "upper")

    if kd["krayl"] is not None:
        tau_rayleigh = compute_tau_rayleigh(kd, interp, col_gas, col_dry)
        t = tau + tau_rayleigh
        ssa = jnp.where(t > 2.0 * _TINY,
                        tau_rayleigh / jnp.where(t == 0.0, 1.0, t), 0.0)
        optics = {"tau": t, "ssa": ssa, "g": jnp.zeros_like(t)}
    else:
        optics = {"tau": tau}
    return optics, col_gas, interp


def _interpolate1d(val, offset, delta, table):
    """interpolate1D: val (...,), table (tab_d1, tab_d2) ->
    (..., tab_d2)."""
    val0 = (val - offset) / delta
    frac = val0 - val0.astype(jnp.int32).astype(val.dtype)
    index = jnp.clip(val0.astype(jnp.int32) + 1, 1,
                     table.shape[0] - 1) - 1
    lo = table[index]           # (..., tab_d2)
    hi = table[index + 1]
    return lo + frac[..., None] * (hi - lo)


def source(kd, play, plev, tlay, tsfc, interp, top_at_1=True, tlev=None):
    """GasOpticsRRTMGPK::source + compute_Planck_source.

    tlev is REQUIRED: EAMxx always provides t_lev (computed in the
    process interface), so the C++'s internal tlev-interpolation branch
    is dead code there — and its ilay==nlay boundary formula reads
    tlay/play out of bounds and divides by zero, so it is not
    transcribable anyway.

    Returns dict with sfc_src (ncol,ngpt), lay_src, lev_src_inc,
    lev_src_dec (ncol,nlay,ngpt)."""
    if tlev is None:
        raise ValueError("source() requires tlev (see docstring)")
    play = jnp.asarray(play)
    plev = jnp.asarray(plev)
    tlay = jnp.asarray(tlay)
    tsfc = jnp.asarray(tsfc)
    pfracin = jnp.asarray(kd["planck_frac"])
    totplnk = jnp.asarray(kd["totplnk"])
    gpoint_bands = jnp.asarray(kd["gpt2band"])
    gpoint_flavor = jnp.asarray(kd["gpoint_flavor"])
    ngpt = kd["ngpt"]
    ncol, nlay = tlay.shape
    tlev = jnp.asarray(tlev)

    sfc_lay = nlay - 1 if top_at_1 else 0

    # pfrac (ncol, nlay, ngpt)
    tropo = interp["tropo"]
    itropo = jnp.where(tropo, 0, 1)
    igpt = jnp.arange(ngpt)[None, None, :]
    itropo_b = itropo[..., None]
    iflav = gpoint_flavor[itropo_b, igpt]
    jpress_loc = interp["jpress"][..., None] + itropo_b
    jt = interp["jtemp"][..., None]
    ic = jnp.arange(ncol)[:, None, None]
    il = jnp.arange(nlay)[None, :, None]
    jeta0 = interp["jeta"][0, iflav, ic, il]
    jeta1 = interp["jeta"][1, iflav, ic, il]

    def fmaj(e, p, t):
        return interp["fmajor"][e, p, t, iflav, ic, il]

    pfrac = (
        fmaj(0, 0, 0) * pfracin[igpt, jeta0, jpress_loc - 1, jt] +
        fmaj(1, 0, 0) * pfracin[igpt, jeta0 + 1, jpress_loc - 1, jt] +
        fmaj(0, 1, 0) * pfracin[igpt, jeta0, jpress_loc, jt] +
        fmaj(1, 1, 0) * pfracin[igpt, jeta0 + 1, jpress_loc, jt] +
        fmaj(0, 0, 1) * pfracin[igpt, jeta1, jpress_loc - 1, jt + 1] +
        fmaj(1, 0, 1) * pfracin[igpt, jeta1 + 1, jpress_loc - 1, jt + 1] +
        fmaj(0, 1, 1) * pfracin[igpt, jeta1, jpress_loc, jt + 1] +
        fmaj(1, 1, 1) * pfracin[igpt, jeta1 + 1, jpress_loc, jt + 1])

    tmin, tdelta = kd["temp_ref_min"], kd["totplnk_delta"]
    planck_sfc = _interpolate1d(tsfc, tmin, tdelta, totplnk)      # (ncol,nbnd)
    planck_lay = _interpolate1d(tlay, tmin, tdelta, totplnk)      # (ncol,nlay,nbnd)
    planck_lev = _interpolate1d(tlev, tmin, tdelta, totplnk)      # (ncol,nlay+1,nbnd)

    band = gpoint_bands  # (ngpt,)
    sfc_src = pfrac[:, sfc_lay, :] * planck_sfc[:, band]
    lay_src = pfrac * planck_lay[:, :, band]
    lev_src_dec = pfrac * planck_lev[:, :-1, :][:, :, band]
    lev_src_inc = pfrac * planck_lev[:, 1:, :][:, :, band]
    return {"sfc_src": sfc_src, "lay_src": lay_src,
            "lev_src_inc": lev_src_inc, "lev_src_dec": lev_src_dec,
            "tlev": tlev}


def gas_optics_sw(kd, play, plev, tlay, vmr, col_dry=None):
    """External-source gas_optics: returns (optics{tau,ssa,g}, toa_flux
    (ncol,ngpt), col_gas)."""
    optics, col_gas, _ = compute_gas_taus(kd, play, plev, tlay, vmr,
                                          col_dry)
    ncol = jnp.asarray(play).shape[0]
    toa_flux = jnp.broadcast_to(jnp.asarray(kd["solar_src"])[None, :],
                                (ncol, kd["ngpt"]))
    return optics, toa_flux, col_gas


def gas_optics_lw(kd, play, plev, tlay, tsfc, vmr, col_dry=None,
                  top_at_1=True, tlev=None):
    """Internal-source gas_optics: returns (optics{tau}, sources dict,
    col_gas)."""
    optics, col_gas, interp = compute_gas_taus(kd, play, plev, tlay, vmr,
                                               col_dry)
    sources = source(kd, play, plev, tlay, tsfc, interp, top_at_1, tlev)
    return optics, sources, col_gas
