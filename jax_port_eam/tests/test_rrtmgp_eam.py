"""Golden replay + Tier-0 property tests for the EAM RRTMGP radiation
port (eam_jax/rrtmgp/).

Golden: harness/gen_rrtmgp_golden.py run in the scream-dev container
against the REAL EAM Fortran RRTMGP stack (external f90 kernels +
driver layer; see harness/build_rrtmgp.py). Coefficient files are the
EAM defaults (the g112/g128 210809 k-distributions — identical files
staged under atm/scream/init — and the atm/cam/physprops cloud-optics
tables); sha256 recorded in the golden metadata.

MEASURED AGREEMENT (asserted tolerances leave ~10x headroom):
  gas optics SW/LW           <= 1.9e-15 rel   (segment-sum reorder of
                                              minor-gas additions)
  cloud optics SW/LW         <= 6.6e-16 rel
  MCICA masks + sampling     bit-exact (KISS RNG reproduced exactly)
  solvers (run_sw / run_lw)  <= 4.1e-13 rel, <= 2e-12 W/m2 abs
  full step fluxes           <= 4.1e-13 rel, <= 6.1e-11 W/m2 abs
  full step heating rates    qrs/qrsc <= 2.7e-12 rel; qrl 1.6e-11,
                             qrlc 1.4e-10 rel (flux-difference
                             cancellation in small clear-sky values;
                             abs <= 2.8e-13)
"""
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from eam_jax.rrtmgp import (cloud_optics as co, coefficients, driver,
                            gas_optics, mcica)

GOLD = Path(__file__).resolve().parents[1] / "golden" / "rrtmgp_golden.npz"
D = np.load(GOLD)
META = json.loads(str(D["__metadata__"]))

INPUTDATA = Path("/Users/ostensiblyowen/development/vibecode_scream/"
                 "e3sm-inputdata")
SW_FILE = INPUTDATA / "atm/scream/init/rrtmgp-data-sw-g112-210809.nc"
LW_FILE = INPUTDATA / "atm/scream/init/rrtmgp-data-lw-g128-210809.nc"
LIQ_FILE = INPUTDATA / ("atm/cam/physprops/"
                        "F_nwvl200_mu20_lam50_res64_t298_c080428.nc")
ICE_FILE = INPUTDATA / "atm/cam/physprops/iceoptics_c080917.nc"

ZEROS = np.zeros_like(D["cld"])


@lru_cache(maxsize=None)
def kdists():
    kd_sw = coefficients.load_kdist(str(SW_FILE), driver.ACTIVE_GASES)
    kd_lw = coefficients.load_kdist(str(LW_FILE), driver.ACTIVE_GASES)
    return kd_sw, kd_lw


@lru_cache(maxsize=None)
def cld_tables():
    return co.load_liq_optics(str(LIQ_FILE)), co.load_ice_optics(str(ICE_FILE))


def gas_vmr_rad():
    g = D["gas_vmr"]
    return np.concatenate([g[:, :, :1], g], axis=2)


@lru_cache(maxsize=None)
def step_out(cfg):
    kd_sw, kd_lw = kdists()
    liq, ice = cld_tables()
    if cfg == "a":
        args = (True, D["coszrs_a"], D["cld"], D["cldfsnow"], D["iclwp"],
                D["iciwp"], D["icswp"], D["lambdac"], D["mu"], D["dei"],
                D["des"])
    elif cfg == "b":
        args = (False, D["coszrs_b"], D["cld"], ZEROS, D["iclwp"],
                D["iciwp"], ZEROS, D["lambdac"], D["mu"], D["dei"],
                D["des"])
    elif cfg == "c":
        args = (True, D["coszrs_a"], ZEROS, ZEROS, ZEROS, ZEROS, ZEROS,
                ZEROS, ZEROS, ZEROS, ZEROS)
    else:  # JAX-only: permanent night
        args = (True, np.full(D["cld"].shape[0], -0.3), D["cld"],
                D["cldfsnow"], D["iclwp"], D["iciwp"], D["icswp"],
                D["lambdac"], D["mu"], D["dei"], D["des"])
    do_snow, coszrs, cld, cldfsnow, iclwp, iciwp, icswp, lambdac, mu, \
        dei, des = args
    out = driver.rad_step(
        kd_sw, kd_lw, liq, ice, do_snow,
        D["t"], D["pmid"], D["pint"], D["lnpmid"], D["lnpint"],
        D["lwup"], D["asdir"], D["asdif"], D["aldir"], D["aldif"], coszrs,
        cld, cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des,
        D["gas_vmr"], D["aer_tau_sw"], D["aer_ssa_sw"], D["aer_asm_sw"],
        D["aer_tau_lw"], float(D["tsi_scaling"]))
    return {k: np.asarray(v) for k, v in out.items()}


# ------------------------------------------------------------------
# structure
# ------------------------------------------------------------------
def test_band_structure():
    kd_sw, kd_lw = kdists()
    assert kd_sw["nband"] == META["nswbands"]
    assert kd_lw["nband"] == META["nlwbands"]
    assert kd_sw["ngpt"] == META["nswgpts"]
    assert kd_lw["ngpt"] == META["nlwgpts"]
    np.testing.assert_array_equal(kd_sw["gpt2band"] + 1, D["gpb_sw"])
    np.testing.assert_array_equal(kd_lw["gpt2band"] + 1, D["gpb_lw"])
    assert kd_sw["temp_ref_min"] == D["temp_limits"][0]
    assert kd_sw["temp_ref_max"] == D["temp_limits"][1]


# ------------------------------------------------------------------
# gas optics
# ------------------------------------------------------------------
def test_gas_optics_sw():
    kd_sw, _ = kdists()
    vmr = driver._vmr_for_kdist(kd_sw, gas_vmr_rad())
    optics, toa, _ = gas_optics.gas_optics_sw(
        kd_sw, D["pmid_rad"], D["pint_rad"], D["a_diag_tmid"], vmr)
    np.testing.assert_allclose(optics["tau"], D["gassw_tau"],
                               rtol=1e-12, atol=1e-280)
    np.testing.assert_allclose(optics["ssa"], D["gassw_ssa"],
                               rtol=1e-12, atol=1e-15)
    np.testing.assert_array_equal(np.asarray(optics["g"]), D["gassw_g"])
    np.testing.assert_array_equal(np.asarray(toa), D["gassw_toa"])


def test_gas_optics_lw():
    _, kd_lw = kdists()
    vmr = driver._vmr_for_kdist(kd_lw, gas_vmr_rad())
    tint = D["a_diag_tint"]
    optics, sources, _ = gas_optics.gas_optics_lw(
        kd_lw, D["pmid_rad"], D["pint_rad"], D["a_diag_tmid"],
        tint[:, -1], vmr, tlev=tint)
    np.testing.assert_allclose(optics["tau"], D["gaslw_tau"],
                               rtol=1e-12, atol=1e-280)
    for k, g in [("lay_src", "gaslw_lay_src"),
                 ("lev_src_inc", "gaslw_lev_src_inc"),
                 ("lev_src_dec", "gaslw_lev_src_dec"),
                 ("sfc_src", "gaslw_sfc_src")]:
        np.testing.assert_allclose(sources[k], D[g], rtol=1e-12,
                                   atol=1e-280)


# ------------------------------------------------------------------
# cloud optics
# ------------------------------------------------------------------
@pytest.mark.parametrize("tag,do_snow", [("snw", True), ("nosnw", False)])
def test_cloud_optics(tag, do_snow):
    liq, ice = cld_tables()
    snowf = D["cldfsnow"] if do_snow else ZEROS
    swpf = D["icswp"] if do_snow else ZEROS
    sw = co.get_cloud_optics_sw(liq, ice, do_snow, D["cld"], snowf,
                                D["iclwp"], D["iciwp"], swpf,
                                D["lambdac"], D["mu"], D["dei"], D["des"])
    for k, g in [("tau", "tau"), ("ssa", "ssa"), ("asm", "asm"),
                 ("liq_tau", "liq_tau"), ("ice_tau", "ice_tau"),
                 ("snw_tau", "snw_tau")]:
        np.testing.assert_allclose(sw[k], D[f"cldsw_{tag}_{g}"],
                                   rtol=1e-12, atol=1e-280)
    lw = co.get_cloud_optics_lw(liq, ice, do_snow, D["cld"], snowf,
                                D["iclwp"], D["iciwp"], swpf,
                                D["lambdac"], D["mu"], D["dei"], D["des"])
    np.testing.assert_allclose(lw["tau"], D[f"cldlw_{tag}_tau"],
                               rtol=1e-12, atol=1e-280)


def test_cloud_optics_tiny_lwp_zero():
    # the iclwp = 1e-82 wisp (col 0, lev 40) hits the < 1e-80 gate
    liq, ice = cld_tables()
    sw = co.get_cloud_optics_sw(liq, ice, True, D["cld"], D["cldfsnow"],
                                D["iclwp"], D["iciwp"], D["icswp"],
                                D["lambdac"], D["mu"], D["dei"], D["des"])
    assert float(np.asarray(sw["liq_tau"])[0, 40].max()) == 0.0


# ------------------------------------------------------------------
# MCICA
# ------------------------------------------------------------------
def test_mcica_masks_bit_exact():
    msw = mcica.mcica_subcol_mask(112, D["pmid"], D["mcica_cldf"],
                                  changeseed=1)
    mlw = mcica.mcica_subcol_mask(128, D["pmid"], D["mcica_cldf"],
                                  changeseed=1)
    np.testing.assert_array_equal(np.asarray(msw).astype(int),
                                  D["mcica_mask_sw"])
    np.testing.assert_array_equal(np.asarray(mlw).astype(int),
                                  D["mcica_mask_lw"])


def test_mcica_sampling_bit_exact():
    kd_sw, kd_lw = kdists()
    RR = driver.RRTMG_TO_RRTMGP_SWBANDS
    tg, sg, ag = mcica.sample_cloud_optics_sw(
        D["pmid"], D["cld"], D["cldfsnow"],
        D["cldsw_snw_tau"][:, :, RR], D["cldsw_snw_ssa"][:, :, RR],
        D["cldsw_snw_asm"][:, :, RR], kd_sw["gpt2band"])
    np.testing.assert_array_equal(np.asarray(tg), D["sample_sw_tau"])
    np.testing.assert_array_equal(np.asarray(sg), D["sample_sw_ssa"])
    np.testing.assert_array_equal(np.asarray(ag), D["sample_sw_asm"])
    tlg = mcica.sample_cloud_optics_lw(D["pmid"], D["cld"], D["cldfsnow"],
                                       D["cldlw_snw_tau"],
                                       kd_lw["gpt2band"])
    np.testing.assert_array_equal(np.asarray(tlg), D["sample_lw_tau"])


def test_mcica_mask_properties():
    # Tier-0: clear layers never cloudy; overcast layers always cloudy
    cldf = np.zeros((4, 12))
    cldf[:, 3] = 1.0
    pmid = np.linspace(5000.0, 99000.0, 12)[None, :] \
        + np.linspace(0.1, 0.7, 4)[:, None]
    mask = np.asarray(mcica.mcica_subcol_mask(16, pmid, cldf))
    assert mask[:, :, 3].all()
    assert not mask[:, :, [0, 1, 2] + list(range(4, 12))].any()


# ------------------------------------------------------------------
# solvers (rrtmgp_run_sw / rrtmgp_run_lw as EAM drives them)
# ------------------------------------------------------------------
def test_run_sw():
    kd_sw, _ = kdists()
    allsky, clrsky = driver.rrtmgp_run_sw(
        kd_sw, D["runsw_gas_vmr"], D["runsw_pmid"], D["runsw_tmid"],
        D["runsw_pint"], D["runsw_coszrs"], D["runsw_alb_dir"],
        D["runsw_alb_dif"], D["runsw_cld_tau"], D["runsw_cld_ssa"],
        D["runsw_cld_asm"], D["runsw_aer_tau"], D["runsw_aer_ssa"],
        D["runsw_aer_asm"], float(D["tsi_scaling"]))
    names = ["flux_up", "flux_dn", "flux_net", "flux_dn_dir"]
    for i, n in enumerate(names):
        np.testing.assert_allclose(allsky[n], D["runsw_flx_all"][:, :, i],
                                   rtol=1e-12, atol=1e-11)
        np.testing.assert_allclose(clrsky[n], D["runsw_flx_clr"][:, :, i],
                                   rtol=1e-12, atol=1e-11)
        np.testing.assert_allclose(allsky["bnd_" + n],
                                   D["runsw_bnd_all"][:, :, :, i],
                                   rtol=1e-12, atol=1e-11)
        np.testing.assert_allclose(clrsky["bnd_" + n],
                                   D["runsw_bnd_clr"][:, :, :, i],
                                   rtol=1e-12, atol=1e-11)


def test_run_lw():
    _, kd_lw = kdists()
    allsky, clrsky = driver.rrtmgp_run_lw(
        kd_lw, D["runlw_gas_vmr"], D["runlw_pmid"], D["runlw_tmid"],
        D["runlw_pint"], D["runlw_tint"], D["runlw_sfc_emis"],
        D["runlw_cld_tau"], D["runlw_aer_tau"])
    for i, n in enumerate(["flux_up", "flux_dn", "flux_net"]):
        np.testing.assert_allclose(allsky[n], D["runlw_flx_all"][:, :, i],
                                   rtol=1e-12, atol=1e-11)
        np.testing.assert_allclose(clrsky[n], D["runlw_flx_clr"][:, :, i],
                                   rtol=1e-12, atol=1e-11)
        np.testing.assert_allclose(allsky["bnd_" + n],
                                   D["runlw_bnd_all"][:, :, :, i],
                                   rtol=1e-12, atol=1e-11)
        np.testing.assert_allclose(clrsky["bnd_" + n],
                                   D["runlw_bnd_clr"][:, :, :, i],
                                   rtol=1e-12, atol=1e-11)


# ------------------------------------------------------------------
# full radiation step
# ------------------------------------------------------------------
STEP_TOLS = {
    "qrs": (1e-11, 1e-11), "qrsc": (1e-11, 1e-11),
    "qrl": (1e-9, 1e-11), "qrlc": (2e-9, 1e-11),
    "sw_all": (1e-11, 1e-9), "sw_clr": (1e-11, 1e-9),
    "lw_all": (1e-11, 1e-9), "lw_clr": (1e-11, 1e-9),
    "srf": (1e-11, 1e-9),
    "diag_tmid": (1e-14, 0.0), "diag_tint": (1e-14, 0.0),
    "diag_alb": (0.0, 0.0),
    "diag_cld_tau_bnd_sw": (1e-12, 1e-280),
    "diag_cld_gpt_sw": (1e-12, 1e-280),
    "diag_cld_tau_bnd_lw": (1e-12, 1e-280),
    "diag_cld_gpt_lw": (1e-12, 1e-280),
}


@pytest.mark.parametrize("cfg", ["a", "b", "c"])
def test_rad_step(cfg):
    out = step_out(cfg)
    for name, (rtol, atol) in STEP_TOLS.items():
        np.testing.assert_allclose(out[name], D[f"{cfg}_{name}"],
                                   rtol=rtol, atol=atol,
                                   err_msg=f"{cfg} {name}")


# ------------------------------------------------------------------
# Tier-0 properties (JAX side alone)
# ------------------------------------------------------------------
def test_heating_consistent_with_flux_divergence():
    # column integral of qrs*dp/g equals the net-flux convergence
    # between the model-top interface and the surface
    out = step_out("a")
    dp = D["pint"][:, 1:] - D["pint"][:, :-1]
    for hr, flx in [(out["qrs"], out["sw_all"]), (out["qrl"], out["lw_all"])]:
        col = (hr * dp).sum(axis=1) / driver.GRAVIT
        net = (flx[:, :, 1] - flx[:, :, 0])          # dn - up
        # qrs_k = -(net_{k+1} - net_k) g / dp_k, so the column integral
        # equals net(model-top interface) - net(surface)
        conv = net[:, 1] - net[:, -1]
        np.testing.assert_allclose(col, conv, rtol=1e-9, atol=1e-9)


def test_clear_sky_olr_bound():
    # clouds can only reduce OLR: all-sky FLUT <= clear-sky FLUT
    out = step_out("a")
    flut_all = out["lw_all"][:, 0, 0]
    flut_clr = out["lw_clr"][:, 0, 0]
    assert np.all(flut_all <= flut_clr + 1e-8)


def test_zero_sun_zero_sw():
    # night columns of config a produce exactly zero SW everywhere
    out = step_out("a")
    night = np.asarray(D["coszrs_a"]) <= 0.0
    assert np.all(out["sw_all"][night] == 0.0)
    assert np.all(out["sw_clr"][night] == 0.0)
    assert np.all(out["qrs"][night] == 0.0)
    # all-night step: the nday == 0 branch zeroes SW but LW still runs
    out_n = step_out("night")
    assert np.all(out_n["sw_all"] == 0.0)
    assert np.all(out_n["qrs"] == 0.0)
    assert np.all(out_n["lw_all"][:, 0, 0] > 0.0)


def test_clear_config_identities():
    # config c: no clouds -> all-sky == clear-sky (not bitwise even in
    # the Fortran: incrementing by a zero-tau/ssa=1 cloud recomputes
    # ssa as (tau*ssa)/tau, a 1-ulp perturbation) and cloud optics zero
    out = step_out("c")
    np.testing.assert_allclose(out["sw_all"], out["sw_clr"],
                               rtol=1e-12, atol=1e-9)
    np.testing.assert_array_equal(out["lw_all"], out["lw_clr"])
    assert np.all(out["diag_cld_tau_bnd_sw"] == 0.0)
    assert np.all(out["diag_cld_gpt_lw"] == 0.0)


def test_surface_partition_sums_to_fsds():
    # soll + sols + solld + solsd == fsds (band partition identity)
    out = step_out("a")
    srf = out["srf"]
    total = srf[:, 5] + srf[:, 6] + srf[:, 7] + srf[:, 8]
    np.testing.assert_allclose(total, srf[:, 0], rtol=1e-12, atol=1e-9)


def test_sw_flux_positivity_and_direct_bound():
    out = step_out("a")
    assert np.all(out["sw_all"][:, :, :2] >= 0.0)         # up, dn
    # direct beam never exceeds total downward
    assert np.all(out["sw_all"][:, :, 3] <= out["sw_all"][:, :, 1] + 1e-9)
