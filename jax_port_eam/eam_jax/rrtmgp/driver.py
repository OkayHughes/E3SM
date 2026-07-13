"""EAM radiation driver sequencing (radiation.F90 radiation_tend for a
single diagnostic call, plus the f90 rrtmgp_interface run wrappers).

FRESH PORT. Sources, transcribed verbatim:
  components/eam/src/physics/rrtmgp/radiation_state.F90
    (set_rad_state, set_interface_temperature: RRTMG-style tint from
     lnp weights, surface tint = (lwup/stebol)**0.25, extra layer above
     the model top: tmid/tint copied, pmid = 0.5*pint(top),
     pint = 1.01 Pa)
  components/eam/src/physics/rrtmgp/radiation.F90 (radiation_tend
    sequencing with dosw = dolw = .true.: combined cloud/snow fraction,
    temperature clipping to the k-distribution limits, cloud optics ->
    RRTMG->RRTMGP band reorder -> MCICA sampling -> optics clipping ->
    radiation_driver_sw/lw [day-column compression, empty top level,
    rrtmgp_run_sw/lw, day expansion, heating rates], set_albedo
    [visible/IR split at 14286 cm^-1, straddling bands averaged with a
    DEFAULT-REAL 0.5 which is exact in binary], set_net_fluxes_sw/lw,
    export_surface_fluxes [NIR bands 1-9 + half of band 10, visible
    bands 11-14 + half of band 10])
  components/eam/src/physics/rrtmgp/f90/rrtmgp_interface.F90
    (rrtmgp_run_sw: cloud gpt optics delta-scaled, aerosol band optics
     delta-scaled, gas_optics + tsi scaling, clear sky = gas + aerosol,
     all sky += clouds, both through rte_sw; rrtmgp_run_lw analogous
     with 1-scalar optics [delta scale is a no-op], surface emissivity
     1, t_sfc = tint(:, nlev+1), n_gauss_angles = 1)
  components/eam/src/physics/rrtmgp/radiation_utils.F90
    (calculate_heating_rate: dF * gravit / dp, clip_values = clamp)
  components/eam/src/physics/rrtmgp/radconstants.F90 (wavenumber band
    edges and rrtmg_to_rrtmgp_swbands, copied verbatim)

NOT ported (documented scope): the icall/N_DIAG diagnostic loop (one
call), COSP/history/pbuf plumbing, aerosol optics computation (plain
band-array inputs, as radiation_tend receives them from
set_aerosol_optics_sw in RRTMG band order), cosine solar zenith angle
and orbital factors (coszrs and tsi_scaling are plain inputs),
conserve_energy pdel bookkeeping (outside the radiation call).

Gas VMRs enter as an (ngas, ncol, nlev) array in the radiation.F90
active_gases order H2O CO2 O3 N2O CO CH4 O2 N2 (get_gas_vmr is pbuf
plumbing; its CO/N2 defaults are applied by the caller/test).
"""

import numpy as np
import jax.numpy as jnp

from . import gas_optics, rte, mcica, cloud_optics, optical_props

STEBOL = 5.670374419e-8   # shr_const_stebol via physconst (verbatim)
GRAVIT = 9.80616          # shr_const_g via physconst (verbatim)

ACTIVE_GASES = ["h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"]

# radconstants.F90, verbatim (RRTMG band order)
WAVENUM_SW_LOWER = np.array(
    [2600., 3250., 4000., 4650., 5150., 6150., 7700.,
     8050., 12850., 16000., 22650., 29000., 38000., 820.])
WAVENUM_SW_UPPER = np.array(
    [3250., 4000., 4650., 5150., 6150., 7700., 8050.,
     12850., 16000., 22650., 29000., 38000., 50000., 2600.])
RRTMG_TO_RRTMGP_SWBANDS = np.array(
    [14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]) - 1  # 0-based


def reordered(arr_bnd, new_indexing=RRTMG_TO_RRTMGP_SWBANDS):
    """radiation.F90 reordered() along the last (band) axis."""
    return jnp.asarray(arr_bnd)[..., jnp.asarray(new_indexing)]


def clip_values(x, min_x, max_x):
    """radiation_utils clip_values (warnings aside, a pure clamp)."""
    return jnp.clip(jnp.asarray(x), min_x, max_x)


def set_rad_state(t, pmid, pint, lnpmid, lnpint, lwup):
    """radiation_state set_rad_state + set_interface_temperature.
    Returns (tmid, tint, pmid_rad, pint_rad) on the radiation grid
    (nlev_rad = pver + 1 midpoints, nlev_rad + 1 interfaces)."""
    t = jnp.asarray(t)
    pmid = jnp.asarray(pmid)
    pint = jnp.asarray(pint)
    lnpmid = jnp.asarray(lnpmid)
    lnpint = jnp.asarray(lnpint)

    # interface temperatures on the model grid (pver+1 levels)
    dy = (lnpint[:, 1:-1] - lnpmid[:, 1:]) / (lnpmid[:, :-1] - lnpmid[:, 1:])
    tint_mid = t[:, 1:] - dy * (t[:, 1:] - t[:, :-1])
    tint_sfc = jnp.sqrt(jnp.sqrt(jnp.asarray(lwup) / STEBOL))
    tint = jnp.concatenate([t[:, :1], tint_mid, tint_sfc[:, None]], axis=1)

    tmid_rad = jnp.concatenate([t[:, :1], t], axis=1)
    tint_rad = jnp.concatenate([tint[:, :1], tint], axis=1)
    pmid_rad = jnp.concatenate([0.5 * pint[:, :1], pmid], axis=1)
    pint_rad = jnp.concatenate([jnp.full_like(pint[:, :1], 1.01), pint],
                               axis=1)
    return tmid_rad, tint_rad, pmid_rad, pint_rad


def is_visible(wavenumber):
    return wavenumber > 14286.0


def set_albedo(asdir, asdif, aldir, aldif):
    """radiation.F90 set_albedo. Returns (albedo_dir, albedo_dif),
    each (nswbands, ncol), in RRTMGP band order, clipped to [0,1]."""
    lower = WAVENUM_SW_LOWER[RRTMG_TO_RRTMGP_SWBANDS]
    upper = WAVENUM_SW_UPPER[RRTMG_TO_RRTMGP_SWBANDS]
    asdir = jnp.asarray(asdir)
    asdif = jnp.asarray(asdif)
    aldir = jnp.asarray(aldir)
    aldif = jnp.asarray(aldif)
    dirs, difs = [], []
    for ib in range(lower.shape[0]):
        if is_visible(lower[ib]) and is_visible(upper[ib]):
            dirs.append(asdir)
            difs.append(asdif)
        elif not is_visible(lower[ib]) and not is_visible(upper[ib]):
            dirs.append(aldir)
            difs.append(aldif)
        else:
            dirs.append(0.5 * (aldir + asdir))
            difs.append(0.5 * (aldif + asdif))
    albedo_dir = jnp.clip(jnp.stack(dirs), 0.0, 1.0)
    albedo_dif = jnp.clip(jnp.stack(difs), 0.0, 1.0)
    return albedo_dir, albedo_dif


def calculate_heating_rate(flux_up, flux_dn, pint):
    """radiation_utils calculate_heating_rate: dF * g / dp."""
    return ((flux_up[:, 1:] - flux_up[:, :-1]
             - flux_dn[:, 1:] + flux_dn[:, :-1])
            * GRAVIT / (pint[:, 1:] - pint[:, :-1]))


def _vmr_for_kdist(kd, gas_vmr_rad):
    """(ngas, ncol, nlev) in ACTIVE_GASES order -> (ncol, nlev, ngas)
    ordered like kd['gas_names']."""
    idx = [ACTIVE_GASES.index(g) for g in kd["gas_names"]]
    return jnp.moveaxis(jnp.asarray(gas_vmr_rad)[jnp.asarray(idx)], 0, -1)


def rrtmgp_run_sw(kd_sw, gas_vmr_rad, pmid, tmid, pint, coszrs,
                  albedo_dir, albedo_dif, cld_tau, cld_ssa, cld_asm,
                  aer_tau, aer_ssa, aer_asm, tsi_scaling):
    """f90 rrtmgp_interface rrtmgp_run_sw + mo_rrtmgp_clr_all_sky
    rte_sw. All inputs on the radiation grid, daytime columns only.
    Returns (allsky, clrsky) flux dicts."""
    cld = optical_props.delta_scale_2str(
        {"tau": jnp.asarray(cld_tau), "ssa": jnp.asarray(cld_ssa),
         "g": jnp.asarray(cld_asm)})
    aer = optical_props.delta_scale_2str(
        {"tau": jnp.asarray(aer_tau), "ssa": jnp.asarray(aer_ssa),
         "g": jnp.asarray(aer_asm)})

    vmr = _vmr_for_kdist(kd_sw, gas_vmr_rad)
    gas, toa_flux, _ = gas_optics.gas_optics_sw(kd_sw, pmid, pint, tmid,
                                                vmr)
    toa_flux = toa_flux * tsi_scaling

    clr_optics = optical_props.increment_2stream_by_2stream(
        gas, aer, gpt2band=kd_sw["gpt2band"])
    clrsky = rte.rte_sw(clr_optics, coszrs, toa_flux, albedo_dir,
                        albedo_dif, kd_sw["gpt2band"], kd_sw["band2gpt"])
    all_optics = optical_props.increment_2stream_by_2stream(clr_optics, cld)
    allsky = rte.rte_sw(all_optics, coszrs, toa_flux, albedo_dir,
                        albedo_dif, kd_sw["gpt2band"], kd_sw["band2gpt"])
    return allsky, clrsky


def rrtmgp_run_lw(kd_lw, gas_vmr_rad, pmid, tmid, pint, tint, sfc_emis,
                  cld_tau, aer_tau):
    """f90 rrtmgp_interface rrtmgp_run_lw + mo_rrtmgp_clr_all_sky
    rte_lw (delta scaling of 1-scalar optics is a no-op)."""
    vmr = _vmr_for_kdist(kd_lw, gas_vmr_rad)
    tint = jnp.asarray(tint)
    gas, sources, _ = gas_optics.gas_optics_lw(
        kd_lw, pmid, pint, tmid, tint[:, -1], vmr, tlev=tint)
    clr_optics = optical_props.increment_1scalar_by_1scalar(
        gas, {"tau": jnp.asarray(aer_tau)}, gpt2band=kd_lw["gpt2band"])
    clrsky = rte.rte_lw(clr_optics, sources, sfc_emis, kd_lw["gpt2band"],
                        kd_lw["band2gpt"])
    all_optics = optical_props.increment_1scalar_by_1scalar(
        clr_optics, {"tau": jnp.asarray(cld_tau)})
    allsky = rte.rte_lw(all_optics, sources, sfc_emis, kd_lw["gpt2band"],
                        kd_lw["band2gpt"])
    return allsky, clrsky


def _pad_top(a):
    """Add the empty radiation level above the model top (zeros)."""
    a = jnp.asarray(a)
    return jnp.concatenate([jnp.zeros_like(a[:, :1]), a], axis=1)


def rad_step(kd_sw, kd_lw, liq, ice, do_snow,
             t, pmid, pint, lnpmid, lnpint,
             lwup, asdir, asdif, aldir, aldif, coszrs,
             cld, cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des,
             gas_vmr, aer_tau_sw, aer_ssa_sw, aer_asm_sw, aer_tau_lw,
             tsi_scaling):
    """radiation_tend for one call: returns a dict mirroring the
    Fortran harness drv_rad_step outputs."""
    ncol, nlev = np.shape(t)
    coszrs = np.asarray(coszrs)

    c_cldf = jnp.maximum(jnp.asarray(cld), jnp.asarray(cldfsnow))

    tmid, tint, pmid_rad, pint_rad = set_rad_state(
        t, pmid, pint, lnpmid, lnpint, lwup)
    tmin = min(kd_sw["temp_ref_min"], kd_lw["temp_ref_min"])
    tmax = max(kd_sw["temp_ref_max"], kd_lw["temp_ref_max"])
    tmid = clip_values(tmid, tmin, tmax)
    tint = clip_values(tint, tmin, tmax)

    out = {"diag_tmid": tmid, "diag_tint": tint}

    # ---------------- shortwave ----------------
    albedo_dir, albedo_dif = set_albedo(asdir, asdif, aldir, aldif)
    out["diag_alb"] = jnp.stack([albedo_dir, albedo_dif], axis=-1)

    cldsw = cloud_optics.get_cloud_optics_sw(
        liq, ice, do_snow, cld, cldfsnow, iclwp, iciwp, icswp,
        lambdac, mu, dei, des)
    tau_bnd_sw = reordered(cldsw["tau"])
    ssa_bnd_sw = reordered(cldsw["ssa"])
    asm_bnd_sw = reordered(cldsw["asm"])
    out["diag_cld_tau_bnd_sw"] = tau_bnd_sw

    gpb_sw = kd_sw["gpt2band"]
    tau_gpt_sw, ssa_gpt_sw, asm_gpt_sw = mcica.sample_cloud_optics_sw(
        pmid, cld, cldfsnow, tau_bnd_sw, ssa_bnd_sw, asm_bnd_sw, gpb_sw)

    atau_sw = reordered(aer_tau_sw)
    assa_sw = reordered(aer_ssa_sw)
    aasm_sw = reordered(aer_asm_sw)

    tau_gpt_sw = clip_values(tau_gpt_sw, 0.0, np.finfo(np.float64).max)
    ssa_gpt_sw = clip_values(ssa_gpt_sw, 0.0, 1.0)
    asm_gpt_sw = clip_values(asm_gpt_sw, -1.0, 1.0)
    atau_sw = clip_values(atau_sw, 0.0, np.finfo(np.float64).max)
    assa_sw = clip_values(assa_sw, 0.0, 1.0)
    aasm_sw = clip_values(aasm_sw, -1.0, 1.0)
    out["diag_cld_gpt_sw"] = jnp.stack(
        [tau_gpt_sw, ssa_gpt_sw, asm_gpt_sw], axis=-1)

    day = np.where(coszrs > 0.0)[0]
    nday = day.size
    nlevp2 = nlev + 2
    nswb = tau_bnd_sw.shape[-1]

    sw_names = ["flux_up", "flux_dn", "flux_net", "flux_dn_dir"]
    if nday == 0:
        sw_all = {k: jnp.zeros((ncol, nlevp2)) for k in sw_names}
        sw_clr = {k: jnp.zeros((ncol, nlevp2)) for k in sw_names}
        sw_all_b = {k: jnp.zeros((ncol, nlevp2, nswb)) for k in sw_names}
        sw_clr_b = {k: jnp.zeros((ncol, nlevp2, nswb)) for k in sw_names}
        qrs = jnp.zeros((ncol, nlev))
        qrsc = jnp.zeros((ncol, nlev))
    else:
        gas_vmr = jnp.asarray(gas_vmr)
        gas_vmr_rad = jnp.concatenate([gas_vmr[:, :, :1], gas_vmr], axis=2)
        allsky_d, clrsky_d = rrtmgp_run_sw(
            kd_sw, gas_vmr_rad[:, day, :], pmid_rad[day], tmid[day],
            pint_rad[day], coszrs[day],
            albedo_dir[:, day], albedo_dif[:, day],
            _pad_top(tau_gpt_sw[day]), _pad_top(ssa_gpt_sw[day]),
            _pad_top(asm_gpt_sw[day]),
            _pad_top(atau_sw[day]), _pad_top(assa_sw[day]),
            _pad_top(aasm_sw[day]), tsi_scaling)

        def expand(d):
            bb, bnd = {}, {}
            for k in sw_names:
                full = jnp.zeros((ncol, nlevp2)).at[day].set(d[k])
                bb[k] = full
                bfull = jnp.zeros((ncol, nlevp2, nswb)).at[day].set(
                    d["bnd_" + k])
                bnd[k] = bfull
            return bb, bnd

        sw_all, sw_all_b = expand(allsky_d)
        sw_clr, sw_clr_b = expand(clrsky_d)
        qrs = calculate_heating_rate(sw_all["flux_up"][:, 1:],
                                     sw_all["flux_dn"][:, 1:],
                                     pint_rad[:, 1:])
        qrsc = calculate_heating_rate(sw_clr["flux_up"][:, 1:],
                                      sw_clr["flux_dn"][:, 1:],
                                      pint_rad[:, 1:])

    out["qrs"] = qrs
    out["qrsc"] = qrsc
    out["sw_all"] = jnp.stack([sw_all[k] for k in sw_names], axis=-1)
    out["sw_clr"] = jnp.stack([sw_clr[k] for k in sw_names], axis=-1)

    # surface exports (set_net_fluxes_sw + export_surface_fluxes)
    srf = jnp.zeros((ncol, 11))
    kb = nlevp2 - 1              # kbot+1 in radiation-grid indexing
    kt = 1                       # ktop
    srf = srf.at[:, 0].set(sw_all["flux_dn"][:, kb])                 # fsds
    srf = srf.at[:, 1].set(sw_all["flux_dn"][:, kb]
                           - sw_all["flux_up"][:, kb])               # fsns
    srf = srf.at[:, 2].set(sw_all["flux_dn"][:, kt]
                           - sw_all["flux_up"][:, kt])               # fsnt
    dn_dir_sfc = sw_all_b["flux_dn_dir"][:, kb, :]
    dn_dif_sfc = sw_all_b["flux_dn"][:, kb, :] - dn_dir_sfc
    srf = srf.at[:, 5].set(dn_dir_sfc[:, :9].sum(-1)
                           + 0.5 * dn_dir_sfc[:, 9])                 # soll
    srf = srf.at[:, 6].set(0.5 * dn_dir_sfc[:, 9]
                           + dn_dir_sfc[:, 10:14].sum(-1))           # sols
    srf = srf.at[:, 7].set(dn_dif_sfc[:, :9].sum(-1)
                           + 0.5 * dn_dif_sfc[:, 9])                 # solld
    srf = srf.at[:, 8].set(0.5 * dn_dif_sfc[:, 9]
                           + dn_dif_sfc[:, 10:14].sum(-1))           # solsd
    srf = srf.at[:, 9].set(sw_all["flux_net"][:, kb])                # netsw

    # ---------------- longwave ----------------
    cldlw = cloud_optics.get_cloud_optics_lw(
        liq, ice, do_snow, cld, cldfsnow, iclwp, iciwp, icswp,
        lambdac, mu, dei, des)
    out["diag_cld_tau_bnd_lw"] = cldlw["tau"]
    tau_gpt_lw = mcica.sample_cloud_optics_lw(
        pmid, cld, cldfsnow, cldlw["tau"], kd_lw["gpt2band"])
    atau_lw = clip_values(aer_tau_lw, 0.0, np.finfo(np.float64).max)
    tau_gpt_lw = clip_values(tau_gpt_lw, 0.0, np.finfo(np.float64).max)
    out["diag_cld_gpt_lw"] = tau_gpt_lw

    gas_vmr = jnp.asarray(gas_vmr)
    gas_vmr_rad = jnp.concatenate([gas_vmr[:, :, :1], gas_vmr], axis=2)
    nlwb = np.asarray(aer_tau_lw).shape[-1]
    sfc_emis = jnp.ones((ncol, nlwb)).T   # (nlwbands, ncol)
    allsky, clrsky = rrtmgp_run_lw(
        kd_lw, gas_vmr_rad, pmid_rad, tmid, pint_rad, tint, sfc_emis,
        _pad_top(tau_gpt_lw), _pad_top(atau_lw))

    qrl_rad = calculate_heating_rate(allsky["flux_up"], allsky["flux_dn"],
                                     pint_rad)
    qrlc_rad = calculate_heating_rate(clrsky["flux_up"], clrsky["flux_dn"],
                                      pint_rad)
    out["qrl"] = qrl_rad[:, 1:]
    out["qrlc"] = qrlc_rad[:, 1:]
    lw_names = ["flux_up", "flux_dn", "flux_net"]
    out["lw_all"] = jnp.stack([allsky[k] for k in lw_names], axis=-1)
    out["lw_clr"] = jnp.stack([clrsky[k] for k in lw_names], axis=-1)

    srf = srf.at[:, 3].set(allsky["flux_up"][:, kb]
                           - allsky["flux_dn"][:, kb])               # flns
    srf = srf.at[:, 4].set(allsky["flux_up"][:, kt]
                           - allsky["flux_dn"][:, kt])               # flnt
    srf = srf.at[:, 10].set(allsky["flux_dn"][:, kb])                # flwds
    out["srf"] = srf
    return out
