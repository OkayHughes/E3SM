#!/usr/bin/env python3
"""Tier-1 golden generator for the EAM RRTMGP radiation layer (run in
the scream-dev container; see build_rrtmgp.py).

Coefficient files: EAM's namelist defaults are
atm/cam/rad/rrtmgp-data-{sw-g112,lw-g128}-210809.nc; the identical
files (same names/content, staged for SCREAM) live locally in
e3sm-inputdata/atm/scream/init and are used here. Cloud optics files
are the EAM defaults (atm/cam/physprops), fetched from the E3SM
inputdata server.

16 columns x 72 levels of synthetic but physically consistent
profiles: tropical/midlat/polar temperature families, hydrostatic-ish
pressure with per-column surface-pressure noise (the MCICA KISS seeds
are the fractional parts of the bottom-four pmid values), realistic
trace-gas profiles, three cloud decks (low liquid / mixed / cirrus)
plus snow, banded aerosol optics.

Configs (same thermo profiles):
  a: do_snow=1, 11 day + 5 night columns, cloudy   (EAMv3 default path)
  b: do_snow=0, all-day, cloudy (cldfsnow == 0)
  c: do_snow=1, mixed day/night, fully clear sky

Golden groups:
  kernel: gas optics SW/LW on the padded radiation grid, cloud optics
    SW/LW (both snow branches), MCICA masks (112 and 128 gpts),
    MCICA-sampled gpt optics, rrtmgp_run_sw / rrtmgp_run_lw solver
    calls on the day-compressed padded arrays (inputs reconstructed
    in numpy exactly as radiation_driver_sw does: pure copies + one
    0.5*pint multiply).
  step: drv_rad_step (full radiation_tend sequencing for one icall)
    outputs: qrs/qrsc/qrl/qrlc, broadband fluxes, surface exports,
    diagnostic intermediates (tmid/tint post-clip, albedos, band and
    gpt cloud optics).
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_rrtmgp_f  # noqa: E402

d = eam_rrtmgp_f.rrtmgp_driver

INIT = "/work/e3sm-inputdata/atm/scream/init"
PHYS = "/work/e3sm-inputdata/atm/cam/physprops"
SW_FILE = INIT + "/rrtmgp-data-sw-g112-210809.nc"
LW_FILE = INIT + "/rrtmgp-data-lw-g128-210809.nc"
LIQ_FILE = PHYS + "/F_nwvl200_mu20_lam50_res64_t298_c080428.nc"
ICE_FILE = PHYS + "/iceoptics_c080917.nc"

F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731
FI = lambda a: np.asarray(a, dtype=np.int32, order="F")   # noqa: E731

STEBOL = 5.67e-8
ncol, nlev = 16, 72
NGAS = 8

d.w_init(SW_FILE.ljust(256), LW_FILE.ljust(256),
         LIQ_FILE.ljust(256), ICE_FILE.ljust(256))
nswb, nlwb, nswg, nlwg = d.w_dims()
tmin_k, tmax_k = d.w_temp_limits()
gpb_sw, gpb_lw = (np.asarray(a) for a in d.w_gpt_bands(nswg, nlwg))
print(f"dims: nswb={nswb} nlwb={nlwb} nswg={nswg} nlwg={nlwg} "
      f"T=[{tmin_k},{tmax_k}]")

# ------------------------------------------------------------------
# synthetic profiles
# ------------------------------------------------------------------
rng = np.random.default_rng(20260712)

ptop = 10.0                                     # ~0.1 hPa model top
psfc = 1.0e5 * (1 + 0.03 * rng.uniform(-1, 1, ncol)) \
    + rng.uniform(0, 1, ncol)                   # non-degenerate fracs
s = np.linspace(0.0, 1.0, nlev + 1)
pint = ptop + (psfc[:, None] - ptop) * s[None, :] ** 2.2
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
lnpmid = np.log(pmid)
lnpint = np.log(pint)

# temperature: 3 families (tropical / midlat / polar)
t_srf = np.array([300, 296, 288, 282, 275, 265, 255, 302, 297, 290,
                  284, 271, 262, 250, 295, 285], dtype=np.float64)
z_approx = -7600.0 * np.log(pmid / psfc[:, None])
t_trop = np.maximum(t_srf[:, None] - 6.2e-3 * z_approx, 197.0)
t_strat = 197.0 + 12.0 * np.log(np.maximum(20000.0 / pmid, 1.0)) ** 1.1
T = np.where(pmid > 12000.0, t_trop, np.minimum(t_strat, 290.0))
T = T + rng.uniform(-0.4, 0.4, (ncol, nlev))
T = np.clip(T, tmin_k + 2.0, tmax_k - 2.0)

# gas vmrs
def qsat_bolton(T, p):
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    return 0.622 * es / np.maximum(p - es, 1.0)

rh = np.clip(0.75 * (pmid / psfc[:, None]) ** 0.6, 0.02, 0.9)
q = np.minimum(rh * qsat_bolton(T, pmid), 0.025)
h2o_vmr = q / (1.0 - q) * 28.97 / 18.01528
o3 = 3.0e-8 + 8.0e-6 * np.exp(-0.5 * ((np.log(pmid) - np.log(2000.0))
                                      / 0.9) ** 2)
gas_vmr = np.zeros((NGAS, ncol, nlev))
gas_vmr[0] = h2o_vmr
gas_vmr[1] = 397.0e-6                      # CO2
gas_vmr[2] = o3                            # O3
gas_vmr[3] = 320.0e-9 * np.clip(pmid / 5.0e4, 0.3, 1.0)  # N2O
gas_vmr[4] = 1.0e-7                        # CO
gas_vmr[5] = 1.8e-6 * np.clip(pmid / 3.0e4, 0.4, 1.0)    # CH4
gas_vmr[6] = 0.2095                        # O2
gas_vmr[7] = 0.7906                        # N2

# clouds: low liquid deck, mixed-phase deck, cirrus; snow below mixed
cld = np.zeros((ncol, nlev))
cldfsnow = np.zeros((ncol, nlev))
iclwp = np.zeros((ncol, nlev))
iciwp = np.zeros((ncol, nlev))
icswp = np.zeros((ncol, nlev))
lambdac = np.zeros((ncol, nlev))
mu = np.zeros((ncol, nlev))
dei = np.zeros((ncol, nlev))
des = np.zeros((ncol, nlev))
rel = np.full((ncol, nlev), 8.0)     # unused by gammadist/mitchell
rei = np.full((ncol, nlev), 25.0)    # but passed through the interface

for i in range(ncol):
    # low liquid deck 900-800 hPa
    low = (pmid[i] > 8.0e4) & (pmid[i] < 9.3e4)
    cld[i][low] = 0.35 + 0.6 * rng.uniform(0, 1, low.sum())
    iclwp[i][low] = 0.02 + 0.25 * rng.uniform(0, 1, low.sum())
    lambdac[i][low] = 10.0 ** rng.uniform(4.9, 6.8, low.sum())
    mu[i][low] = rng.uniform(1.0, 16.0, low.sum())   # exercises bndry
    # mixed deck 700-500 hPa (liquid + ice + snow)
    mid = (pmid[i] > 5.0e4) & (pmid[i] < 7.0e4) & (i % 3 != 0)
    cld[i][mid] = 0.25 + 0.55 * rng.uniform(0, 1, mid.sum())
    iclwp[i][mid] = 0.01 + 0.08 * rng.uniform(0, 1, mid.sum())
    lambdac[i][mid] = 10.0 ** rng.uniform(5.0, 6.5, mid.sum())
    mu[i][mid] = rng.uniform(2.0, 14.0, mid.sum())
    iciwp[i][mid] = 0.002 + 0.05 * rng.uniform(0, 1, mid.sum())
    dei[i][mid] = rng.uniform(20.0, 180.0, mid.sum())
    cldfsnow[i][mid] = np.maximum(
        cld[i][mid] * rng.uniform(0.2, 1.3, mid.sum()), 0.0)
    icswp[i][mid] = 0.005 + 0.06 * rng.uniform(0, 1, mid.sum())
    des[i][mid] = rng.uniform(60.0, 500.0, mid.sum())
    # cirrus 300-150 hPa
    cir = (pmid[i] > 1.5e4) & (pmid[i] < 3.0e4) & (i % 4 != 1)
    cld[i][cir] = 0.15 + 0.75 * rng.uniform(0, 1, cir.sum())
    iciwp[i][cir] = 0.0005 + 0.02 * rng.uniform(0, 1, cir.sum())
    dei[i][cir] = rng.uniform(15.0, 90.0, cir.sum())
# tiny-value edge: sub-1e-80 path
iclwp[0, 40] = 1.0e-82
lambdac[0, 40] = 2.0e5
mu[0, 40] = 3.0
cld[0, 40] = 0.3

cldfsnow = np.clip(cldfsnow, 0.0, 1.0)

# aerosols (RRTMG band order, as set_aerosol_optics_sw returns them)
zfac = np.exp(-z_approx / 2500.0)
aer_tau_sw = np.zeros((ncol, nlev, nswb))
aer_ssa_sw = np.ones((ncol, nlev, nswb))
aer_asm_sw = np.zeros((ncol, nlev, nswb))
band_fac = 0.3 + 0.7 * rng.uniform(0, 1, nswb)
aer_tau_sw[:] = (0.12 * zfac * (0.3 + rng.uniform(0, 1, (ncol, 1))))[
    :, :, None] * band_fac[None, None, :]
aer_ssa_sw[:] = np.clip(0.92 + 0.07 * rng.uniform(-1, 1,
                                                  (ncol, nlev, nswb)),
                        0.0, 1.0)
aer_asm_sw[:] = np.clip(0.62 + 0.15 * rng.uniform(-1, 1,
                                                  (ncol, nlev, nswb)),
                        -1.0, 1.0)
aer_tau_lw = (0.03 * zfac * (0.3 + rng.uniform(0, 1, (ncol, 1))))[
    :, :, None] * (0.2 + 0.8 * rng.uniform(0, 1, nlwb))[None, None, :]

# surface / solar
lwup = STEBOL * (t_srf + rng.uniform(-4, 8, ncol)) ** 4
asdir = np.clip(rng.uniform(0.04, 0.6, ncol), 0, 1)
asdif = np.clip(asdir + rng.uniform(-0.03, 0.03, ncol), 0, 1)
aldir = np.clip(asdir + rng.uniform(0.0, 0.25, ncol), 0, 1)
aldif = np.clip(aldir + rng.uniform(-0.03, 0.03, ncol), 0, 1)
coszrs_a = rng.uniform(0.05, 0.98, ncol)
coszrs_a[[2, 5, 9, 12, 15]] = [0.0, -0.1, 0.0, -0.4, 0.0]  # night
coszrs_b = rng.uniform(0.08, 0.95, ncol)
tsi_scaling = 1.0123456

zeros = np.zeros((ncol, nlev))

out = {
    "t": T, "pmid": pmid, "pint": pint, "lnpmid": lnpmid,
    "lnpint": lnpint, "gas_vmr": gas_vmr, "lwup": lwup,
    "asdir": asdir, "asdif": asdif, "aldir": aldir, "aldif": aldif,
    "coszrs_a": coszrs_a, "coszrs_b": coszrs_b,
    "cld": cld, "cldfsnow": cldfsnow, "iclwp": iclwp, "iciwp": iciwp,
    "icswp": icswp, "lambdac": lambdac, "mu": mu, "dei": dei,
    "des": des, "rel": rel, "rei": rei,
    "aer_tau_sw": aer_tau_sw, "aer_ssa_sw": aer_ssa_sw,
    "aer_asm_sw": aer_asm_sw, "aer_tau_lw": aer_tau_lw,
    "tsi_scaling": np.float64(tsi_scaling),
    "gpb_sw": gpb_sw, "gpb_lw": gpb_lw,
    "temp_limits": np.array([tmin_k, tmax_k]),
}

# ------------------------------------------------------------------
# full radiation step (configs a, b, c)
# ------------------------------------------------------------------
STEP_OUT = ["qrs", "qrsc", "qrl", "qrlc", "sw_all", "sw_clr", "lw_all",
            "lw_clr", "srf", "diag_tmid", "diag_tint", "diag_alb",
            "diag_cld_tau_bnd_sw", "diag_cld_gpt_sw",
            "diag_cld_tau_bnd_lw", "diag_cld_gpt_lw"]


def run_step(cfg, do_snow, coszrs, cldx, cldfsnowx, iclwpx, iciwpx,
             icswpx, lambdacx, mux, deix, desx):
    res = d.w_rad_step(
        nswg, nlwg, do_snow,
        F(T), F(pmid), F(pint), F(lnpmid), F(lnpint),
        F(lwup), F(asdir), F(asdif), F(aldir), F(aldif), F(coszrs),
        F(cldx), F(cldfsnowx), F(iclwpx), F(iciwpx), F(icswpx),
        F(lambdacx), F(mux), F(deix), F(desx), F(rel), F(rei),
        F(gas_vmr), F(aer_tau_sw), F(aer_ssa_sw), F(aer_asm_sw),
        F(aer_tau_lw), tsi_scaling)
    for name, v in zip(STEP_OUT, res):
        out[f"{cfg}_{name}"] = np.asarray(v)
        assert np.all(np.isfinite(np.asarray(v))), f"{cfg} {name} nonfinite"
    print(f"step {cfg}: FSNT={out[f'{cfg}_srf'][:, 2].max():.1f} "
          f"FLNT={out[f'{cfg}_srf'][:, 4].max():.1f} "
          f"qrl range=[{out[f'{cfg}_qrl'].min():.2e},"
          f"{out[f'{cfg}_qrl'].max():.2e}]")


run_step("a", 1, coszrs_a, cld, cldfsnow, iclwp, iciwp, icswp,
         lambdac, mu, dei, des)
run_step("b", 0, coszrs_b, cld, zeros, iclwp, iciwp, zeros,
         lambdac, mu, dei, des)
run_step("c", 1, coszrs_a, zeros, zeros, zeros, zeros, zeros,
         zeros, zeros, zeros, zeros)

# ------------------------------------------------------------------
# kernel goldens (radiation grid: one empty level above model top)
# ------------------------------------------------------------------
nlev_rad = nlev + 1
pmid_rad = np.concatenate([0.5 * pint[:, :1], pmid], axis=1)
pint_rad = np.concatenate([np.full((ncol, 1), 1.01), pint], axis=1)
tmid_rad = out["a_diag_tmid"]      # post set_rad_state + clip
tint_rad = out["a_diag_tint"]
gas_vmr_rad = np.concatenate([gas_vmr[:, :, :1], gas_vmr], axis=2)
out["pmid_rad"] = pmid_rad
out["pint_rad"] = pint_rad

# gas optics SW + LW
tau, ssa, g, toa = d.w_gas_optics_sw(
    nswg, F(gas_vmr_rad), F(pmid_rad), F(tmid_rad), F(pint_rad))
out["gassw_tau"] = np.asarray(tau)
out["gassw_ssa"] = np.asarray(ssa)
out["gassw_g"] = np.asarray(g)
out["gassw_toa"] = np.asarray(toa)

tau, lay, inc, dec, sfc = d.w_gas_optics_lw(
    nlwg, F(gas_vmr_rad), F(pmid_rad), F(tmid_rad), F(pint_rad),
    F(tint_rad))
out["gaslw_tau"] = np.asarray(tau)
out["gaslw_lay_src"] = np.asarray(lay)
out["gaslw_lev_src_inc"] = np.asarray(inc)
out["gaslw_lev_src_dec"] = np.asarray(dec)
out["gaslw_sfc_src"] = np.asarray(sfc)

# cloud optics (both snow branches)
for tag, dsn, snowf, swpf in (("snw", 1, cldfsnow, icswp),
                              ("nosnw", 0, zeros, zeros)):
    res = d.w_cloud_optics_sw(nswb, dsn, F(cld), F(snowf), F(iclwp),
                              F(iciwp), F(swpf), F(lambdac), F(mu),
                              F(dei), F(des), F(rel), F(rei))
    for name, v in zip(["tau", "ssa", "asm", "liq_tau", "ice_tau",
                        "snw_tau"], res):
        out[f"cldsw_{tag}_{name}"] = np.asarray(v)
    res = d.w_cloud_optics_lw(nlwb, dsn, F(cld), F(snowf), F(iclwp),
                              F(iciwp), F(swpf), F(lambdac), F(mu),
                              F(dei), F(des), F(rei))
    for name, v in zip(["tau", "liq_tau", "ice_tau", "snw_tau"], res):
        out[f"cldlw_{tag}_{name}"] = np.asarray(v)

# MCICA masks + sampling (combined cloud fraction, changeseed=1)
c_cldf = np.maximum(cld, cldfsnow)
out["mcica_cldf"] = c_cldf
out["mcica_mask_sw"] = np.asarray(
    d.w_mcica_mask(nswg, 1, F(pmid), F(c_cldf)))
out["mcica_mask_lw"] = np.asarray(
    d.w_mcica_mask(nlwg, 1, F(pmid), F(c_cldf)))

# band optics in RRTMGP order: exact permutation of the standalone
# cloud-optics golden (radconstants rrtmg_to_rrtmgp_swbands, 0-based)
RRTMG2GP = np.array([14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]) - 1
tau_bnd_sw = out["cldsw_snw_tau"][:, :, RRTMG2GP]
ssa_bnd_sw = out["cldsw_snw_ssa"][:, :, RRTMG2GP]
asm_bnd_sw = out["cldsw_snw_asm"][:, :, RRTMG2GP]
np.testing.assert_array_equal(tau_bnd_sw, out["a_diag_cld_tau_bnd_sw"])

res = d.w_sample_sw(FI(gpb_sw), F(pmid), F(cld),
                    F(cldfsnow), F(tau_bnd_sw), F(ssa_bnd_sw),
                    F(asm_bnd_sw))
for name, v in zip(["tau", "ssa", "asm"], res):
    out[f"sample_sw_{name}"] = np.asarray(v)
tau_bnd_lw = out["cldlw_snw_tau"]     # LW bands are not reordered
out["sample_lw_tau"] = np.asarray(
    d.w_sample_lw(FI(gpb_lw), F(pmid), F(cld), F(cldfsnow),
                  F(tau_bnd_lw)))

# ------------------------------------------------------------------
# solver goldens: rrtmgp_run_sw / rrtmgp_run_lw exactly as
# radiation_driver_sw/lw call them (day compression + padding are
# pure copies, reconstructed here in numpy)
# ------------------------------------------------------------------
day = np.where(coszrs_a > 0.0)[0]
nday = day.size
cld_gpt_sw = out["a_diag_cld_gpt_sw"]        # post-clip (tau,ssa,asm)
cld_gpt_lw = out["a_diag_cld_gpt_lw"]
alb_dir = out["a_diag_alb"][:, :, 0]          # (nswb, ncol)
alb_dif = out["a_diag_alb"][:, :, 1]

pad = lambda a: np.concatenate(  # noqa: E731
    [np.zeros_like(a[:, :1, :]), a], axis=1)
aer_sw_gp = aer_tau_sw[:, :, RRTMG2GP], aer_ssa_sw[:, :, RRTMG2GP], \
    aer_asm_sw[:, :, RRTMG2GP]

runsw_in = dict(
    gas_vmr=gas_vmr_rad[:, day, :], pmid=pmid_rad[day], tmid=tmid_rad[day],
    pint=pint_rad[day], coszrs=coszrs_a[day],
    alb_dir=alb_dir[:, day], alb_dif=alb_dif[:, day],
    cld_tau=pad(cld_gpt_sw[day, :, :, 0]),
    cld_ssa=pad(cld_gpt_sw[day, :, :, 1]),
    cld_asm=pad(cld_gpt_sw[day, :, :, 2]),
    aer_tau=pad(aer_sw_gp[0][day]), aer_ssa=pad(aer_sw_gp[1][day]),
    aer_asm=pad(aer_sw_gp[2][day]))
for k, v in runsw_in.items():
    out[f"runsw_{k}"] = v
res = d.w_run_sw(F(runsw_in["gas_vmr"]),
                 F(runsw_in["pmid"]), F(runsw_in["tmid"]),
                 F(runsw_in["pint"]), F(runsw_in["coszrs"]),
                 F(runsw_in["alb_dir"]), F(runsw_in["alb_dif"]),
                 F(runsw_in["cld_tau"]), F(runsw_in["cld_ssa"]),
                 F(runsw_in["cld_asm"]), F(runsw_in["aer_tau"]),
                 F(runsw_in["aer_ssa"]), F(runsw_in["aer_asm"]),
                 tsi_scaling)
for name, v in zip(["flx_all", "flx_clr", "bnd_all", "bnd_clr"], res):
    out[f"runsw_{name}"] = np.asarray(v)

emis = np.ones((nlwb, ncol))
runlw_in = dict(gas_vmr=gas_vmr_rad, pmid=pmid_rad, tmid=tmid_rad,
                pint=pint_rad, tint=tint_rad, sfc_emis=emis,
                cld_tau=pad(cld_gpt_lw), aer_tau=pad(aer_tau_lw))
for k, v in runlw_in.items():
    out[f"runlw_{k}"] = v
res = d.w_run_lw(F(runlw_in["gas_vmr"]),
                 F(runlw_in["pmid"]), F(runlw_in["tmid"]),
                 F(runlw_in["pint"]), F(runlw_in["tint"]),
                 F(runlw_in["sfc_emis"]), F(runlw_in["cld_tau"]),
                 F(runlw_in["aer_tau"]))
for name, v in zip(["flx_all", "flx_clr", "bnd_all", "bnd_clr"], res):
    out[f"runlw_{name}"] = np.asarray(v)

# consistency: solver goldens must reproduce the step's fluxes
sw_all = out["a_sw_all"]
np.testing.assert_allclose(out["runsw_flx_all"][:, :, 0],
                           sw_all[day, :, 0], rtol=0, atol=0)
np.testing.assert_allclose(out["runlw_flx_all"][:, :, 0],
                           out["a_lw_all"][:, :, 0], rtol=0, atol=0)

# ------------------------------------------------------------------
# metadata + save
# ------------------------------------------------------------------
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {
    "scheme": "EAM RRTMGP radiation layer (f90 path, EAMv3 default)",
    "source_sha": sha, "ncol": ncol, "nlev": nlev, "nlev_rad": nlev_rad,
    "nswbands": int(nswb), "nlwbands": int(nlwb),
    "nswgpts": int(nswg), "nlwgpts": int(nlwg),
    "active_gases": ["H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2", "N2"],
    "coefficients": {
        "sw": {"path": SW_FILE, "sha256": sha256(SW_FILE),
               "note": "identical file to EAM default "
                       "atm/cam/rad/rrtmgp-data-sw-g112-210809.nc"},
        "lw": {"path": LW_FILE, "sha256": sha256(LW_FILE)},
        "liq": {"path": LIQ_FILE, "sha256": sha256(LIQ_FILE),
                "note": "EAM default liqopticsfile (gammadist)"},
        "ice": {"path": ICE_FILE, "sha256": sha256(ICE_FILE),
                "note": "EAM default iceopticsfile (mitchell)"},
    },
    "configs": {
        "a": {"do_snow": True, "night_cols": [2, 5, 9, 12, 15],
              "cloudy": True},
        "b": {"do_snow": False, "night_cols": [], "cloudy": True},
        "c": {"do_snow": True, "night_cols": [2, 5, 9, 12, 15],
              "cloudy": False},
    },
    "step_outputs": STEP_OUT,
    "srf_fields": ["fsds", "fsns", "fsnt", "flns", "flnt", "soll",
                   "sols", "solld", "solsd", "netsw", "flwds"],
    "tsi_scaling": tsi_scaling,
    "icecldoptics": "mitchell", "liqcldoptics": "gammadist",
    "compiler": "gfortran -O2 -ffp-contract=off (fbuild defaults)",
}
gold = Path(__file__).resolve().parents[1] / "golden" / "rrtmgp_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
