#!/usr/bin/env python3
"""Tier-1 golden generator for the ZM precipitation evaporation / snow
production routine (zm_conv.F90 zm_conv_evap, run in the scream-dev
container via the eam_zm_evap_f extension built against the REAL
cloud_fraction.F90 cldfrc_fice).

Realistic configurations (e1/e2) reproduce the exact zm_conv_intr.F90
sequence: the eam_zm_conv_f extension is re-run on the zm_conv golden
profiles (config a: EAMv3 phys="default", is_first_step=T; asserted
bitwise-equal to the archive), the state is advanced exactly as
physics_update does for the 'zm_conv_main' ptend (t1 = t +
dt*heat/cpair; q1 = q + dt*qtnd, then the qneg3 clip at qmin(Q) =
1e-12 -- use_mass_borrower is .false. by default), and zm_conv_evap
runs on (t1, q1) with prdprec = rprd and prec from zm_conv_main, plus
a synthetic pbuf-style cloud fraction:
  e1: old_snow=T (the zm_param_t default; zm_conv_intr only sets
      old_snow=F under zmconv_microp=T, so this is the configuration
      consistent with the harness-wide zm_microp=F scope)
  e2: old_snow=F (the EAMv3 production value, reached via
      zmconv_microp=.true.; with zm_microp=F here prdsnow=0, which
      exercises the ntsnprd=-min(...) branch and the final
      flxsnow<=flxprec protection loop)

Synthetic sweep (e3 old_snow=T dt=1800, e4 old_snow=F dt=900), 40
columns covering: freezing level crossing at varying heights (snow
melt), entirely frozen columns (all snow), warm columns (no snow),
near-saturated columns (evaporation shuts off), cld=1 columns (the
1-cldfrc factor kills evaporation), cld=0 dry columns (max
evaporation), prec_in scaled to 70% of the column production (the
"total evaporation cannot exceed input precipitation" limiter),
negative prdprec patches (downdraft evaporation), and zero-production
no-op columns.

EAMv3 phys="default" evaporation efficiency: zmconv_ke = 2.5E-6
(namelist_defaults_eam.xml phys="default" dyn="se" microphys="p3").
This ZM version has a SINGLE ke: there is no zmconv_ke_lnd /
land-ocean split in zm_conv.F90's zm_conv_evap. omsm = 0.99999 is a
zm_conv.F90 module parameter. PERGRO is off (pergro_active=F).
cldfrc_fice runs with top_lev=1 (see drivers/zm_evap_driver.F90: in
production top_lev sits above the 40 hPa limcnv cap where all ZM
precip fluxes are zero, so it cannot affect any output).
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "fmod"))
import eam_zm_conv_f  # noqa: E402
import eam_zm_evap_f  # noqa: E402

F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731
GRAV = 9.80616
CPAIR = 1.00464e3
QMIN_VAPOR = 1.0e-12   # cnst_add('Q', ..., 1.E-12_r8, ...) in physpkg
KE_V3 = 2.5e-6         # zmconv_ke phys="default" dyn="se" microphys="p3"

de = eam_zm_evap_f.zm_evap_driver
assert "array('d')" in de.drv_zm_conv_evap.__doc__
de.drv_init()

ENAMES = ("tend_s", "tend_q", "tend_s_snwprd", "tend_s_snwevmlt",
          "prec", "snow", "ntprprd", "ntsnprd", "flxprec", "flxsnow")

# --- e1/e2: exact intr sequence on the zm_conv golden profiles ---
g = np.load(ROOT / "golden" / "zm_conv_golden.npz")
zmeta = json.loads(str(g["__metadata__"]))
p = zmeta["params"]["a"]
d = eam_zm_conv_f.zm_conv_driver
d.drv_init()
res = d.drv_zm_conv_main(
    float(g["time_step"]), p["is_first_step"], p["limcnv"],
    p["no_deep_pbl"], p["tau"], p["alfa"], p["ke"], p["dmpdz"],
    p["tpert_fix"], p["tpert_fac"], p["tiedke_add"], p["c0_lnd"],
    p["c0_ocn"], p["num_cin"], p["mx_bot_lyr_adj"], p["trig_dcape"],
    p["trig_ull"], p["clos_dyn_adj"], p["old_snow"],
    F(g["t"]), F(g["q"]), F(g["omega"]), F(g["pmid"]), F(g["pint"]),
    F(g["pdel"]), g["geos"], F(g["zmid"]), F(g["zint"]),
    g["pbl_hgt"], g["tpert"], g["landfrac"], F(g["t_star"]),
    F(g["q_star"]))
NAMES = ("lengath", "gather_index", "msemax_klev_g", "jctop", "jcbot",
         "jt", "prec", "heat", "qtnd", "cape", "dcape", "mcon", "pflx",
         "zdu", "mflx_up", "entr_up", "detr_up", "mflx_dn", "entr_dn",
         "p_del", "dsubcld", "ql", "rliq", "rprd", "dlf")
r = dict(zip(NAMES, (np.asarray(v) for v in res)))
for n in ("heat", "qtnd", "rprd", "prec", "gather_index"):
    assert np.array_equal(r[n], g[f"{n}_a"]), n
ncol, nlev = np.asarray(g["t"]).shape
dt = float(g["time_step"])
lengath = int(r["lengath"])
print(f"zm_conv rerun: lengath={lengath} of ncol={ncol}")

# physics_update for the 'zm_conv_main' ptend (ls, lq(1))
t1 = np.asarray(g["t"]) + dt * r["heat"] / CPAIR
q1_raw = np.asarray(g["q"]) + dt * r["qtnd"]
q1 = np.maximum(q1_raw, QMIN_VAPOR)   # qneg3 clip at qmin(Q)
print(f"qneg3 clip points after zm_conv_main: {(q1_raw < QMIN_VAPOR).sum()}")

# synthetic pbuf-style cloud fraction (deterministic)
rng = np.random.default_rng(20260714)
sig = np.asarray(g["pmid"]) / np.asarray(g["pmid"])[:, -1:]
cld = np.clip(0.25 + 0.55 * np.exp(-0.5 * ((sig - 0.55) / 0.2) ** 2)
              * (1.0 + 0.3 * rng.uniform(-1, 1, (ncol, nlev))), 0.0, 0.9)

zeros2 = np.zeros((ncol, nlev))
out = dict(t1=t1, q1=q1, cld=cld, rprd=r["rprd"], prec_main=r["prec"],
           pmid=np.asarray(g["pmid"]), pdel=np.asarray(g["pdel"]))
for cfg, old_snow in (("e1", 1), ("e2", 0)):
    er = dict(zip(ENAMES, (np.asarray(v) for v in de.drv_zm_conv_evap(
        KE_V3, old_snow, dt, F(g["pmid"]), F(g["pdel"]), F(t1), F(q1),
        F(r["rprd"]), F(cld), F(zeros2), F(zeros2), r["prec"]))))
    out.update({f"{n}_{cfg}": v for n, v in er.items()})
    out[f"old_snow_{cfg}"] = np.array(old_snow)
    out[f"dt_{cfg}"] = np.array(dt)
    assert all(np.isfinite(v).all() for v in er.values())
    # untriggered columns (prec=0, rprd=0) must no-op
    idle = np.asarray(r["prec"]) == 0.0
    assert np.all(er["prec"][idle] == 0.0)
    assert np.all(er["tend_q"][idle] == 0.0)
    # water conservation: flxprec_sfc = int (prdprec - evap) dp/g
    chk = ((r["rprd"] - er["tend_q"]) * np.asarray(g["pdel"])).sum(1) / GRAV
    assert np.abs(chk - er["flxprec"][:, -1]).max() < 1e-18
    print(f"cfg {cfg}: prec_out[mm/day] max="
          f"{er['prec'].max() * 8.64e7:.2f} snow max="
          f"{er['snow'].max() * 8.64e7:.2f} "
          f"evap max={er['tend_q'].max():.2e}")

# all triggered profiles are warm at the surface, so surface snow is
# zero, but the in-column snow flux (and its melt) must be active in
# the old_snow branch; with old_snow=F and zm_microp=F prdsnow=0, so
# flxsnow is identically zero (snow production needs zm_microphysics)
assert out["flxsnow_e1"].max() > 0.0
assert np.all(out["flxsnow_e2"] == 0.0)

# --- e3/e4: synthetic sweep ---
ns, nl = 40, 72
ai = np.linspace(0.0, 1.0, nl + 1) ** 1.7
pint_s = 225.5 + ai * (1.0e5 - 225.5)
pmid_s = np.broadcast_to(0.5 * (pint_s[:-1] + pint_s[1:]), (ns, nl)).copy()
pdel_s = np.broadcast_to(np.diff(pint_s), (ns, nl)).copy()
sig_s = pmid_s / pmid_s[:, -1:]

# surface temperatures sweeping the freezing level from below-ground
# (warm, all rain) to the surface (cold, all snow)
ts = np.linspace(255.0, 310.0, ns)[:, None]
t_s = np.maximum(ts * sig_s ** 0.19, 195.0)
t_s += rng.uniform(-0.4, 0.4, t_s.shape)


def qs_bolton(T, p):
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    return 0.622 * es / (p - es)


rh = np.clip(0.55 * sig_s ** 0.8, 0.02, 0.55)
rh[6::7] = 0.999                       # near-saturated columns
rh[3::11] = 0.03                       # very dry columns
q_s = np.maximum(rh * qs_bolton(t_s, pmid_s), 1.0e-9)

cld_s = np.clip(0.3 + 0.4 * np.sin(7.0 * sig_s)
                * rng.uniform(0.5, 1.0, (ns, nl)), 0.0, 0.95)
cld_s[5::9] = 1.0                      # (1-cldfrc)=0: no evaporation
cld_s[4::9] = 0.0                      # max evaporation

# precip production: mid-level bump, deeper for warmer columns,
# with negative patches (downdraft evaporation) in some columns
prd_s = 4.0e-8 * np.exp(-0.5 * ((sig_s - 0.55) / 0.15) ** 2) \
    * (1.0 + 0.4 * rng.uniform(-1, 1, (ns, nl)))
prd_s[2::5, 60:66] = -1.5e-8           # negative production patches
prd_s[7::13] = 0.0                     # no-op columns
colint = (np.maximum(prd_s, 0.0) * pdel_s).sum(1) / GRAV
prec_in_s = colint / 1000.0
prec_in_s[1::3] *= 0.7                 # engage the prec-evpvint limiter
assert (prd_s < 0).any()

out.update(t_s=t_s, q_s=q_s, cld_s=cld_s, prd_s=prd_s,
           pmid_s=pmid_s, pdel_s=pdel_s, prec_in_s=prec_in_s)
zs = np.zeros((ns, nl))
for cfg, old_snow, dts in (("e3", 1, 1800.0), ("e4", 0, 900.0)):
    er = dict(zip(ENAMES, (np.asarray(v) for v in de.drv_zm_conv_evap(
        KE_V3, old_snow, dts, F(pmid_s), F(pdel_s), F(t_s), F(q_s),
        F(prd_s), F(cld_s), F(zs), F(zs), prec_in_s))))
    out.update({f"{n}_{cfg}": v for n, v in er.items()})
    out[f"old_snow_{cfg}"] = np.array(old_snow)
    out[f"dt_{cfg}"] = np.array(dts)
    assert all(np.isfinite(v).all() for v in er.values())
    if old_snow:
        assert er["snow"].max() > 0.0, "cold columns must produce snow"
    else:
        # zm_microp=F: prdsnow=0, so the new-snow branch never makes snow
        assert np.all(er["flxsnow"] == 0.0)
    assert (er["snow"] <= er["prec"] * (1 + 1e-12) + 1e-30).all()
    print(f"cfg {cfg}: prec_out[mm/day] max="
          f"{er['prec'].max() * 8.64e7:.2f} snow max="
          f"{er['snow'].max() * 8.64e7:.2f} "
          f"evap max={er['tend_q'].max():.2e}")

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "zm_conv zm_conv_evap (ZM precip evaporation / snow "
                  "production) + real cloud_fraction cldfrc_fice",
        "source_sha": sha,
        "ncol_realistic": ncol, "ncol_sweep": ns, "nlev": nlev,
        "zm_microp": False, "pergro_active": False,
        "ke": KE_V3,
        "ke_note": "zmconv_ke=2.5e-6 (namelist_defaults_eam.xml "
                   "phys='default' dyn='se' microphys='p3'); single ke "
                   "-- no zmconv_ke_lnd exists in EAMv3 zm_conv.F90",
        "qmin_vapor": QMIN_VAPOR,
        "state_update": "t1=t+dt*heat/cpair; q1=max(q+dt*qtnd,qmin) "
                        "(physics_update + qneg3; use_mass_borrower=F "
                        "default)",
        "configs": {c: {"old_snow": int(out[f"old_snow_{c}"]),
                        "dt": float(out[f"dt_{c}"])}
                    for c in ("e1", "e2", "e3", "e4")},
        "top_lev": 1,
        "top_lev_note": "cldfrc_fice with module default top_lev=1; "
                        "production top_lev=trop_cloud_top_lev is above "
                        "the 40 hPa limcnv cap where ZM precip fluxes "
                        "are identically zero, so outputs are unaffected",
        "constants": {"omsm": 0.99999, "cpair": CPAIR, "grav": GRAV,
                      "latvap": 2.501e6, "latice": 3.337e5,
                      "tfreez": 273.15,
                      "pergro_perturbation": 8.64e-11}}
gold = ROOT / "golden" / "zm_evap_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
