#!/usr/bin/env python3
"""Tier-1 golden generator for the full-spectrum (non-orographic)
gravity-wave suite: gw_beres_src (convective), gw_cm_src (frontal),
gw_drag_prof with ngwv = pgwv = 32 (stress profiles, gwd_project_tau,
tendency limiting, gwd_precalc_rhoi constituent/DSE diffusion via the
vdiff_lu_solver), and momentum_energy_conservation, plus a chained
gw_tend-style sequence (prof -> beres src -> drag_prof -> conservation
-> cm src -> drag_prof -> conservation, accumulating ptend like
gw_drag.F90).

EAMv3 defaults (bld/build-namelist + namelist_defaults_eam.xml
phys="default" + gw_drag.F90 parameters): pgwv=32, gw_dc=2.5,
fcrit2=1.0, kwv=6.28e-5, taubgnd=2.5e-3, effgw_cm=1.0,
effgw_beres=0.35, frontgfc=1.25e-15, gw_convect_hcf=10.0,
hdepth_scaling_factor=0.5, gw_convect_hdepth_min=2.5,
gw_convect_storm_speed_min=10.0, gw_convect_plev_src_wind=70000,
use_gw_convect_old=.true., tau_0_ubc=.false., ktop=0,
kbotbg/kfront/k_src_wind from pref_edge at 500/600/700 hPa,
do_latitude_taper=.false. (SE dycore is UNSTRUCTURED). alpha is the
gw_drag.F90 alpha0/palph table interpolated to pref_edge (lininterp
with boundary clamping, as CAM's lininterp extrapolates flat).

The Beres lookup table mfcc is pure input data read from
inputdata atm/waccm/gw/newmfspectra40_dc25.nc in production; that file
is not on this machine, so a smooth synthetic nonnegative table with
the same shape (maxh=20, -maxuh:maxuh=81, -pgwv:pgwv=65) and a
realistic magnitude is used and stored in the archive (the JAX port
takes the table as an input array either way).
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_gw_spectrum_f  # noqa: E402

d = eam_gw_spectrum_f.gw_spectrum_driver

rng = np.random.default_rng(20260712)
ncol, nlev, pgwv, ncnst = 32, 72, 32, 2
nwav = 2 * pgwv + 1
GRAVIT, RAIR, CPAIR = 9.80616, 287.042, 1004.64
DC, FCRIT2, KWV = 2.5, 1.0, 6.28e-5
TAUBGND, FRONTGFC = 2.5e-3, 1.25e-15
EFFGW_BERES, EFFGW_CM = 0.35, 1.0
HCF, HDSF, HDMIN, SSMIN = 10.0, 0.5, 2.5, 10.0
PLEV_SRC_WIND = 70000.0
DT = 1800.0

# ---- state ---------------------------------------------------------------
ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
pref_edge = 225.5 + ai * (1.0e5 - 225.5)
pint = np.broadcast_to(pref_edge, (ncol, nlev + 1)).copy()
pint *= (1 + 0.03 * rng.uniform(-1, 1, (ncol, 1)))
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
dpm = np.diff(pint, axis=1)
rdpm = 1.0 / dpm
piln = np.log(pint)

t = 216.0 + 72.0 * (pmid / pmid[:, -1:]) ** 0.32
t += 8.0 * np.exp(-((np.log(pmid / 1e2)) ** 2))
t += rng.uniform(-3.0, 3.0, t.shape)
zm = (RAIR * t / GRAVIT) * np.log(pint[:, -1:] / pmid)

zj = 1.1e4
u = 28.0 * np.exp(-((zm - zj) / 8e3) ** 2) + rng.uniform(-5, 5, t.shape)
u += rng.uniform(-4, 16, (ncol, 1))
u[::5] *= -0.7                       # easterly columns
v = 8.0 * np.sin(zm / 6e3) + rng.uniform(-4, 4, t.shape)
lat = np.deg2rad(rng.uniform(-80, 80, ncol))

# convective heating [K/s]: cols 0-3 zero; then blobs of varying top/
# depth/magnitude, some with an interior negative patch (melting-level
# case: old vs new gw_heating_depth differ), some too shallow to launch
netdt = np.zeros((ncol, nlev))
for i in range(4, ncol):
    ztop = rng.uniform(4e3, 15e3)
    zbot = rng.uniform(1e3, min(3.5e3, ztop - 1e3))
    amp = rng.uniform(10.0, 90.0) / 86400.0
    prof = np.exp(-((zm[i] - 0.5 * (ztop + zbot)) / (0.35 * (ztop - zbot)))
                  ** 2)
    blob = np.where((zm[i] > zbot) & (zm[i] < ztop), amp * prof, 0.0)
    if i % 3 == 0:                   # interior negative patch
        zmelt = 0.5 * (ztop + zbot)
        blob = np.where(np.abs(zm[i] - zmelt) < 400.0, -0.1 * amp, blob)
    if i % 7 == 0:                   # too shallow to launch
        blob = np.where(zm[i] > zbot + 900.0, 0.0, blob)
    netdt[i] = blob

# frontogenesis function: some columns above frontgfc, some below
frontgf = np.abs(rng.lognormal(-34.5, 1.4, (ncol, nlev)))
frontgf[1::4] *= 1e-3                # keep some columns quiet

# two synthetic tracers + dry static energy
q = np.stack([1e-3 * np.exp(-zm / 2.5e3) * (1 + 0.1 * np.sin(zm / 3e3)),
              1e-6 + 2e-7 * np.cos(zm / 7e3) ** 2], axis=-1)
dse = CPAIR * t + GRAVIT * zm

# newtonian cooling: gw_drag.F90 alpha0/palph table -> pref_edge
ALPHA0 = np.array([
    1.896007, 1.196965, 0.7251356, 0.6397463, 0.5777858, 0.5712274,
    0.6836302, 0.6678557, 0.5683219, 0.4754283, 0.3960519, 0.332022,
    0.2497581, 0.168667, 0.1323903, 0.1257139, 0.1069889, 0.09873954,
    0.09215571, 0.09398635, 0.1061087, 0.1294598, 0.1544743, 0.1648226,
    0.1687332, 0.1691513, 0.1664987, 0.159048, 0.149292, 0.1351563,
    0.1174998, 0.09913579, 0.08300615, 0.0707, 0.0615588, 0.0542623,
    0.0478562, 0.04132157, 0.03454087, 0.02296682, 0.006723819,
    0.02164464, 0.05756261, 0.003844868, 0.02929285, 0.006627098,
    0.04558291, 0.02042176, 0.0, 0.005880283, 0.00689498, 0.01343466,
    0.0, 0.03415992, 0.02855049, 0.01688839, 0.0272628, 0.02772121,
    0.02135626, 0.04863235, 0.04568304, 0.0, 0.009604108, 0.0, 0.0, 0.0])
PALPH = np.array([
    5.11075e-6, 9.8269e-6, 1.620185e-5, 2.671225e-5, 4.4041e-5,
    7.261275e-5, 1.19719e-4, 1.9738e-4, 3.254225e-4, 5.365325e-4,
    8.846025e-4, 0.001458458, 0.002404575, 0.00397825, 0.006556825,
    0.01081382, 0.017898, 0.02955775, 0.04873075, 0.07991075, 0.1282732,
    0.19812, 0.292025, 0.4101675, 0.55347, 0.73048, 0.9559475, 1.244795,
    1.61285, 2.079325, 2.667425, 3.404875, 4.324575, 5.4654, 6.87285,
    8.599725, 10.70705, 13.26475, 16.35175, 20.05675, 24.479, 29.728,
    35.92325, 43.19375, 51.6775, 61.5205, 72.8745, 85.65715, 100.5147,
    118.2503, 139.1154, 163.6621, 192.5399, 226.5132, 266.4812,
    313.5013, 368.818, 433.8952, 510.4553, 600.5242, 696.7963, 787.7021,
    867.1607, 929.6489, 970.5548, 992.5561])
alpha = np.interp(pref_edge, PALPH * 1e2,
                  np.maximum(ALPHA0 / 86400.0, 1.e-6))

# level indices from reference pressures (gw_drag.F90 gw_init logic)
kbotbg = int(np.sum(pref_edge < 5.0e4)) - 1   # spectrum source at 500mb
kfront = int(np.sum(pref_edge < 6.0e4))       # frontogenesis check level
k_src_wind = int(np.sum(pref_edge < PLEV_SRC_WIND))  # steering level

# synthetic Beres source-spectrum table (see module docstring)
maxh, maxuh = 20, 40
hh = np.arange(1, maxh + 1)[:, None, None]
uu = np.arange(-maxuh, maxuh + 1)[None, :, None]
cc = np.arange(-pgwv, pgwv + 1)[None, None, :] * DC
mfcc = 4.0e8 * (hh / maxh) ** 2 * np.exp(-((cc - 0.4 * uu) / 25.0) ** 2) \
    * (0.05 + (cc - 0.4 * uu) ** 2 / 625.0)

d.drv_gw_spec_init(DC, 0, kbotbg, FCRIT2, KWV, GRAVIT, RAIR, alpha, 0,
                   TAUBGND, FRONTGFC, kfront, PLEV_SRC_WIND,
                   np.asfortranarray(mfcc), pref_edge)

out = dict(pmid=pmid, pint=pint, dpm=dpm, rdpm=rdpm, piln=piln, t=t,
           u=u, v=v, zm=zm, lat=lat, netdt=netdt, frontgf=frontgf, q=q,
           dse=dse, alpha=alpha, mfcc=mfcc, pref_edge=pref_edge)

# ---- gw_prof --------------------------------------------------------------
rhoi, ti, nm, ni = d.drv_gw_prof(CPAIR, t, pmid, pint)
out.update(rhoi=rhoi, ti=ti, nm=nm, ni=ni)

# ---- Beres source, both heating-depth variants ----------------------------
for use_old, sfx in ((1, "_old"), (0, "_new")):
    r = d.drv_gw_beres_src(nwav, pgwv, lat, u, v, netdt, zm,
                           HCF, HDSF, HDMIN, SSMIN, use_old)
    keys = ("src_level", "tend_level", "tau", "ubm", "ubi", "xv", "yv",
            "c", "hdepth", "maxq0")
    out.update({f"b{k}{sfx}": v for k, v in zip(keys, r)})

# ---- frontal source -------------------------------------------------------
r = d.drv_gw_cm_src(nwav, pgwv, kbotbg, u, v, frontgf)
for key, val in zip(("src_level", "tend_level", "tau", "ubm", "ubi",
                     "xv", "yv", "c"), r):
    out[f"c{key}"] = val

# ---- gwd_project_tau standalone kernel (on the Beres source tau) ----------
out["ptaucd"] = d.drv_gwd_project_tau(pgwv, out["btend_level_old"],
                                      out["btau_old"], out["bubi_old"],
                                      out["bc_old"], out["bxv_old"],
                                      out["byv_old"])

# ---- gw_drag_prof: Beres (default old), C&M, and C&M with taper -----------
runs = [("bd", "b", "_old", EFFGW_BERES, 0),
        ("cd", "c", "", EFFGW_CM, 0),
        ("cdt", "c", "", EFFGW_CM, 1)]
dragkeys = ("tau", "utgw", "vtgw", "ttgw", "qtgw", "taucd", "egwdffi",
            "gwut", "dttdf", "dttke")
for tag, src, sfx, effgw, taper in runs:
    r = d.drv_gw_drag_prof_spec(
        pgwv, out[f"{src}src_level{sfx}"], out[f"{src}tend_level{sfx}"],
        taper, DT, lat, t, ti, pmid, pint, dpm, rdpm, piln, rhoi, nm,
        ni, out[f"{src}ubm{sfx}"], out[f"{src}ubi{sfx}"],
        out[f"{src}xv{sfx}"], out[f"{src}yv{sfx}"], effgw,
        out[f"{src}c{sfx}"], q, dse, out[f"{src}tau{sfx}"])
    out.update({f"{tag}_{k}": v for k, v in zip(dragkeys, r)})

# ---- momentum_energy_conservation + chained ptend (gw_tend order) ---------
# ptend starts from the Beres tendencies (gw_drag.F90 assigns, then MEC)
ptu, ptv, pts = out["bd_utgw"], out["bd_vtgw"], out["bd_ttgw"]
r = d.drv_momentum_energy_conservation(
    out["btend_level_old"], DT, out["bd_taucd"], pint, dpm, u, v,
    ptu, ptv, pts, out["bd_utgw"], out["bd_vtgw"], out["bd_ttgw"])
ptu, ptv, pts, mutgw, mvtgw, mttgw = r
out.update(mec_b_dudt=ptu, mec_b_dvdt=ptv, mec_b_dsdt=pts,
           mec_b_utgw=mutgw, mec_b_vtgw=mvtgw, mec_b_ttgw=mttgw)
ptq = out["bd_qtgw"].copy()

# add C&M and conserve again
ptu = ptu + out["cd_utgw"]
ptv = ptv + out["cd_vtgw"]
pts = pts + out["cd_ttgw"]
ptq = ptq + out["cd_qtgw"]
r = d.drv_momentum_energy_conservation(
    out["ctend_level"], DT, out["cd_taucd"], pint, dpm, u, v,
    ptu, ptv, pts, out["cd_utgw"], out["cd_vtgw"], out["cd_ttgw"])
out.update(chain_dudt=r[0], chain_dvdt=r[1], chain_dsdt=r[2],
           mec_c_utgw=r[3], mec_c_vtgw=r[4], mec_c_ttgw=r[5],
           chain_dqdt=ptq)

# ---- vdiff_lu_solver standalone kernel ------------------------------------
ksrf = np.abs(rng.normal(0.05, 0.02, ncol))
kv = np.abs(rng.lognormal(0.0, 1.0, (ncol, nlev + 1)))
tmpi = np.abs(rng.lognormal(-2.0, 0.7, (ncol, nlev + 1)))
cc_top = np.abs(rng.normal(0.01, 0.005, ncol))
cd_top = rng.normal(0.0, 1e-4, ncol)
qlu = np.stack([q[:, :, 0], q[:, :, 1], dse], axis=-1)
ntop_lu, nbot_lu = 3, 65
r = d.drv_vd_lu(ksrf, kv, tmpi, rdpm, DT, GRAVIT, cc_top, ntop_lu,
                nbot_lu, qlu, cd_top)
out.update(lu_ksrf=ksrf, lu_kv=kv, lu_tmpi=tmpi, lu_cc_top=cc_top,
           lu_cd_top=cd_top, lu_q=qlu, lu_q_out=r[0], lu_ca=r[1],
           lu_cc=r[2], lu_dnom=r[3], lu_ze=r[4])

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "gw spectrum: gw_beres_src + gw_cm_src + "
                  "gw_drag_prof(ngwv=pgwv) + momentum_energy_conservation"
                  " + vd_lu_decomp/solve",
        "source_sha": sha, "ncol": ncol, "nlev": nlev, "dt": DT,
        "ncnst": ncnst,
        "params": {"pgwv": pgwv, "dc": DC, "fcrit2": FCRIT2, "kwv": KWV,
                   "taubgnd": TAUBGND, "frontgfc": FRONTGFC,
                   "effgw_beres": EFFGW_BERES, "effgw_cm": EFFGW_CM,
                   "gw_convect_hcf": HCF, "hdepth_scaling_factor": HDSF,
                   "gw_convect_hdepth_min": HDMIN,
                   "gw_convect_storm_speed_min": SSMIN,
                   "gw_convect_plev_src_wind": PLEV_SRC_WIND,
                   "use_gw_convect_old": True, "tau_0_ubc": False,
                   "ktop": 0, "kbotbg": kbotbg, "kfront": kfront,
                   "k_src_wind": k_src_wind, "maxh": maxh,
                   "maxuh": maxuh, "gravit": GRAVIT, "rair": RAIR,
                   "cpair": CPAIR, "lu_ntop": ntop_lu,
                   "lu_nbot": nbot_lu,
                   "mfcc": "synthetic (real Beres04 file not on disk)"},
        "flags": "-O2 -fPIC -ffree-line-length-none "
                 "-fallow-argument-mismatch -std=legacy"}
gold = Path(__file__).resolve().parents[1] / "golden" \
    / "gw_spectrum_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
for k in ("bd_utgw", "cd_utgw", "bd_qtgw", "cd_taucd", "bd_egwdffi",
          "chain_dsdt"):
    print(f"  {k}: max|.| = {np.abs(out[k]).max():.3e}")
