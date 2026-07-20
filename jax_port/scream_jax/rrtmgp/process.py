"""EAMxx RRTMGP process step: RRTMGPRadiation::run_impl.

Source: components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp
plus foundation thermo (calculate_dz, calculate_vmr_from_mmr) and
eamxx_column_ops (compute_interface_values_linear).

This is the unit the EAMxx swap test exchanges: EAMxx fields in, EAMxx
fields out, for a radiation-update step (the caller decides whether
radiation runs this step via radiation_do; on non-update steps only the
T_mid tendency from the stored rad_heating_pdel applies).

Assumes k=0 = model top (EAMxx standard), do_aerosol_rad with provided
aero fields or zeros, and do_subcol_sampling (default true).
"""

import numpy as np
import jax.numpy as jnp

from ..foundation.thermo import calculate_dz
from . import orbital
from . import interface as ifc
from .orbital import GAS_MOL_WEIGHTS, MWDRY

STEBOL = 5.670374419e-8  # scream::physics::Constants stebol
GAS_NAMES = ["h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"]

DEFAULT_PARAMS = {
    "rad_frequency": 1,
    "orbital_year": -9999,
    "orbital_eccentricity": -9999.0,
    "orbital_obliquity": -9999.0,
    "orbital_mvelp": -9999.0,
    "fixed_total_solar_irradiance": -9999.0,
    "fixed_solar_zenith_angle": -9999.0,
    "co2vmr": 388.717e-6,
    "n2ovmr": 323.141e-9,
    "ch4vmr": 1807.851e-9,
    "f11vmr": 768.7644e-12,
    "f12vmr": 531.2820e-12,
    "n2vmr": 0.7906,
    "covmr": 1.0e-7,
    "do_subcol_sampling": True,
    "extra_clnclrsky_diag": False,
    "extra_clnsky_diag": False,
}


def radiation_do(irad, nstep):
    if irad == 0:
        return False
    return nstep == 0 or (nstep % irad) == 0


def orbit_and_declination(params, year, calday):
    """Orbital-parameter and declination logic from run_impl. Returns
    (delta, eccf)."""
    eccen = params["orbital_eccentricity"]
    obliq = params["orbital_obliquity"]
    mvelp = params["orbital_mvelp"]
    orbital_year = params["orbital_year"]
    if eccen >= 0 and obliq >= 0 and mvelp >= 0:
        orbital_year = orbital.SHR_ORB_UNDEF_INT
    elif orbital_year < 0:
        orbital_year = year
    eccen, obliq, mvelp, obliqr, lambm0, mvelpp = orbital.shr_orb_params(
        orbital_year, eccen if eccen >= 0 else None,
        obliq if obliq >= 0 else None, mvelp if mvelp >= 0 else None)
    delta, eccf = orbital.shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr)
    if params["fixed_total_solar_irradiance"] >= 0:
        eccf = params["fixed_total_solar_irradiance"] / 1360.9
    return delta, eccf


def compute_vmr(qv, o3_vmr, lat_deg, p_mid, params):
    """The gas VMR block of run_impl: (ncol, nlay, 8) ordered like
    GAS_NAMES. Pure jnp on the state-dependent inputs (qv, o3, p_mid) so
    it stays traceable under JAX transforms."""
    qv = jnp.asarray(qv)
    gases = []
    for name in GAS_NAMES:
        if name == "h2o":
            # calculate_vmr_from_mmr(mol_weight, qv, mmr=qv): wet mmr ->
            # vmr (dry mole fraction)
            mw = GAS_MOL_WEIGHTS["h2o"]
            gases.append(MWDRY / mw * qv / (1.0 - qv))
        elif name == "o3":
            gases.append(jnp.broadcast_to(jnp.asarray(o3_vmr), qv.shape))
        elif name == "n2":
            gases.append(jnp.full(qv.shape, params["n2vmr"]))
        elif name == "co":
            gases.append(jnp.full(qv.shape, params["covmr"]))
        else:
            mmr = orbital.trcmix(name, lat_deg, p_mid, params["co2vmr"],
                                 params["n2ovmr"], params["ch4vmr"],
                                 params["f11vmr"], params["f12vmr"])
            gases.append(MWDRY / GAS_MOL_WEIGHTS[name] * mmr)
    return jnp.stack(gases, axis=-1)


def get_wavelength_index(band_lims_wvn, wavelength_m):
    """eamxx get_wavelength_index: band containing the wavelength."""
    wl_bounds = 1.0 / np.asarray(band_lims_wvn)  # wavelength in cm
    target = wavelength_m * 1e2
    idx = -1
    for ibnd in range(wl_bounds.shape[1]):
        lo, hi = wl_bounds[0, ibnd], wl_bounds[1, ibnd]
        if (lo < hi and lo <= target <= hi) or (lo >= hi and hi <= target <= lo):
            idx = max(idx, ibnd)
    return idx


def rrtmgp_process_step(kd_sw, kd_lw, co_sw, co_lw, params,
                        dt, nstep, year, calday, lat_deg, lon_deg,
                        # EAMxx fields
                        T_mid, p_mid, p_int, pseudo_density,
                        sfc_alb_dir_vis, sfc_alb_dir_nir,
                        sfc_alb_dif_vis, sfc_alb_dif_nir,
                        qv, qc, nc, qi, cldfrac_tot,
                        eff_radius_qc, eff_radius_qi, surf_lw_flux_up,
                        o3_volume_mix_ratio, rad_heating_pdel,
                        aero_tau_sw=None, aero_ssa_sw=None, aero_g_sw=None,
                        aero_tau_lw=None, smooth_width=0.0):
    """One radiation process step. Returns a dict of updated/computed
    EAMxx fields keyed by their EAMxx names. Aerosol arrays are
    (ncol, nband, nlay) as in the FM; None means no aerosol radiative
    effects (zeros).

    smooth_width (STATIC Python float, never a traced value): 0.0 ->
    exact binary MCICA subcolumn masks (bitwise original); > 0 ->
    approximation-by-identity fractional masks so that
    d(fluxes)/d(cldfrac_tot) is nonzero and finite under AD (see
    mcica.get_subcolumn_mask)."""
    T_mid = jnp.asarray(T_mid, dtype=jnp.float64)
    p_mid = jnp.asarray(p_mid)
    pdel = jnp.asarray(pseudo_density)
    ncol, nlay = T_mid.shape
    out = {}

    update_rad = radiation_do(params["rad_frequency"], nstep)
    if not update_rad:
        rad_heat = jnp.asarray(rad_heating_pdel) / pdel
        out["T_mid"] = T_mid + rad_heat * dt
        out["rad_heating_pdel"] = jnp.asarray(rad_heating_pdel)
        return out

    delta, eccf = orbit_and_declination(params, year, calday)

    vmr = compute_vmr(qv, o3_volume_mix_ratio, lat_deg, p_mid, params)

    # solar zenith angle. Depends only on (calday, lat, lon, orbital
    # params) -- never on prognostic state -- so it is computed host-side
    # with numpy from concrete lat/lon (which come in concrete from the
    # driver) outside the traced region.
    if params["fixed_solar_zenith_angle"] > 0:
        mu0 = np.full(ncol, params["fixed_solar_zenith_angle"])
    else:
        lat_r = np.asarray(lat_deg, dtype=np.float64) * np.pi / 180.0
        lon_r = np.asarray(lon_deg, dtype=np.float64) * np.pi / 180.0
        mu0 = np.asarray(orbital.shr_orb_cosz(
            calday, lat_r, lon_r, delta,
            params["rad_frequency"] * dt), dtype=np.float64).reshape(ncol)
    out["cosine_solar_zenith_angle"] = mu0

    # dz and T at interfaces (T_int bc: T top, blackbody surface bottom)
    dz = calculate_dz(pdel, p_mid, T_mid, qv)
    bc_bot = jnp.sqrt(jnp.sqrt(jnp.asarray(surf_lw_flux_up) / STEBOL))
    t_int_mid = (T_mid[:, 1:] * dz[:, :-1] + T_mid[:, :-1] * dz[:, 1:]) \
        / (dz[:, :-1] + dz[:, 1:])
    t_int = jnp.concatenate([T_mid[:, :1], t_int_mid, bc_bot[:, None]],
                            axis=1)

    # radiative cloud fraction
    cldfrac_tot = jnp.asarray(cldfrac_tot)
    if params["do_subcol_sampling"]:
        cldfrac_rad = cldfrac_tot
    else:
        cldfrac_rad = jnp.where(cldfrac_tot > 0, 1.0, 0.0)
    out["cldfrac_rad"] = cldfrac_rad

    # cloud mass [g/m2]
    lwp = ifc.mixing_ratio_to_cloud_mass(qc, cldfrac_rad, pdel) * 1e3
    iwp = ifc.mixing_ratio_to_cloud_mass(qi, cldfrac_rad, pdel) * 1e3

    alb_dir, alb_dif = ifc.compute_band_by_band_surface_albedos(
        kd_sw["band_lims_wvn"], sfc_alb_dir_vis, sfc_alb_dir_nir,
        sfc_alb_dif_vis, sfc_alb_dif_nir)

    nswb, nlwb = kd_sw["nband"], kd_lw["nband"]
    if aero_tau_sw is None:
        a_tau_sw = jnp.zeros((ncol, nlay, nswb))
        a_ssa_sw = jnp.zeros((ncol, nlay, nswb))
        a_g_sw = jnp.zeros((ncol, nlay, nswb))
        a_tau_lw = jnp.zeros((ncol, nlay, nlwb))
    else:
        # FM keeps (ncol, band, lay); rrtmgp wants (ncol, lay, band)
        a_tau_sw = jnp.transpose(jnp.asarray(aero_tau_sw), (0, 2, 1))
        a_ssa_sw = jnp.transpose(jnp.asarray(aero_ssa_sw), (0, 2, 1))
        a_g_sw = jnp.transpose(jnp.asarray(aero_g_sw), (0, 2, 1))
        a_tau_lw = jnp.transpose(jnp.asarray(aero_tau_lw), (0, 2, 1))

    sw, lw, cld_tau = ifc.rrtmgp_main(
        kd_sw, kd_lw, co_sw, co_lw,
        p_mid, T_mid, p_int, t_int, vmr,
        alb_dir, alb_dif, mu0, lwp, iwp, eff_radius_qc, eff_radius_qi,
        cldfrac_rad, a_tau_sw, a_ssa_sw, a_g_sw, a_tau_lw, eccf,
        params["extra_clnclrsky_diag"], params["extra_clnsky_diag"],
        smooth_width=smooth_width)

    for key, grp, name in (
            ("SW_flux_up", "allsky", "flux_up"),
            ("SW_flux_dn", "allsky", "flux_dn"),
            ("SW_flux_dn_dir", "allsky", "flux_dn_dir"),
            ("SW_clnclrsky_flux_up", "clnclrsky", "flux_up"),
            ("SW_clnclrsky_flux_dn", "clnclrsky", "flux_dn"),
            ("SW_clnclrsky_flux_dn_dir", "clnclrsky", "flux_dn_dir"),
            ("SW_clrsky_flux_up", "clrsky", "flux_up"),
            ("SW_clrsky_flux_dn", "clrsky", "flux_dn"),
            ("SW_clrsky_flux_dn_dir", "clrsky", "flux_dn_dir"),
            ("SW_clnsky_flux_up", "clnsky", "flux_up"),
            ("SW_clnsky_flux_dn", "clnsky", "flux_dn"),
            ("SW_clnsky_flux_dn_dir", "clnsky", "flux_dn_dir")):
        out[key] = jnp.asarray(sw[grp][name])
    for key, grp, name in (
            ("LW_flux_up", "allsky", "flux_up"),
            ("LW_flux_dn", "allsky", "flux_dn"),
            ("LW_clnclrsky_flux_up", "clnclrsky", "flux_up"),
            ("LW_clnclrsky_flux_dn", "clnclrsky", "flux_dn"),
            ("LW_clrsky_flux_up", "clrsky", "flux_up"),
            ("LW_clrsky_flux_dn", "clrsky", "flux_dn"),
            ("LW_clnsky_flux_up", "clnsky", "flux_up"),
            ("LW_clnsky_flux_dn", "clnsky", "flux_dn")):
        out[key] = jnp.asarray(lw[grp][name])

    # heating rate -> T update and pdel-scaled storage
    sw_heat = ifc.compute_heating_rate(
        out["SW_flux_up"], out["SW_flux_dn"], pdel)
    lw_heat = ifc.compute_heating_rate(
        out["LW_flux_up"], out["LW_flux_dn"], pdel)
    heat = sw_heat + lw_heat
    out["T_mid"] = T_mid + heat * dt
    out["rad_heating_pdel"] = pdel * heat

    # surface fluxes
    bnd_dir = jnp.asarray(sw["allsky"]["bnd_flux_dn_dir"])
    bnd_dn = jnp.asarray(sw["allsky"]["bnd_flux_dn"])
    bnd_dif = bnd_dn - bnd_dir
    dvis, dnir, fvis, fnir = ifc.compute_broadband_surface_fluxes(
        kd_sw["band_lims_wvn"], bnd_dir[:, nlay, :], bnd_dif[:, nlay, :])
    out["sfc_flux_dir_vis"] = dvis
    out["sfc_flux_dir_nir"] = dnir
    out["sfc_flux_dif_vis"] = fvis
    out["sfc_flux_dif_nir"] = fnir
    out["sfc_flux_sw_net"] = out["SW_flux_dn"][:, nlay] \
        - out["SW_flux_up"][:, nlay]
    out["sfc_flux_lw_dn"] = out["LW_flux_dn"][:, nlay]

    # cloud area diagnostics from LW gpt cloud tau
    ct = jnp.asarray(cld_tau["lw_gpt"])
    big = np.finfo(np.float64).max
    out["cldlow"] = ifc.compute_cloud_area(p_mid, ct, 700e2, big)
    out["cldmed"] = ifc.compute_cloud_area(p_mid, ct, 400e2, 700e2)
    out["cldhgh"] = ifc.compute_cloud_area(p_mid, ct, 0.0, 400e2)
    out["cldtot"] = ifc.compute_cloud_area(p_mid, ct, 0.0, big)

    # COSP optical depths and sunlit mask (band indices are static:
    # computed host-side from the coefficient tables only)
    i067 = get_wavelength_index(kd_sw["band_lims_wvn"], 0.67e-6)
    i105 = get_wavelength_index(kd_lw["band_lims_wvn"], 10.5e-6)
    out["dtau067"] = jnp.asarray(cld_tau["sw_bnd"])[:, :, i067]
    out["dtau105"] = jnp.asarray(cld_tau["lw_bnd"])[:, :, i105]
    out["sunlit_mask"] = (out["SW_clrsky_flux_dn"][:, 0] > 0).astype(jnp.int32)

    # AeroCom cloud-top diagnostics
    aero = ifc.compute_aerocom_cloudtop(T_mid, p_mid, pdel, dz, qc, qi,
                                        eff_radius_qc, eff_radius_qi,
                                        cldfrac_rad, nc)
    for k, v in aero.items():
        out[k] = v

    # gas VMR output fields (all but o3/n2/co are Computed)
    for igas, name in enumerate(GAS_NAMES):
        if name not in ("o3",):
            out[f"{name}_volume_mix_ratio"] = vmr[..., igas]
    return out
