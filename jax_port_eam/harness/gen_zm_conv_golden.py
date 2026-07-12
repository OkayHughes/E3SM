#!/usr/bin/env python3
"""Tier-1 golden generator for the ZM deep-convection main routine
(zm_conv_main, run in the scream-dev container).

Synthetic hydrostatic 72-level columns, 6 families x 7 columns = 42:
  0-6   strongly unstable tropical  (CAPE O(1000) J/kg, deep plumes)
  7-13  marginally unstable         (small positive CAPE)
  14-20 stable midlatitude          (no trigger -> all outputs no-op)
  21-27 midlatitude convective      (moderate CAPE, land columns)
  28-34 high-CAPE dry boundary layer(hot/dry BL, LCL above the 600 hPa
                                     threshold -> CAPE=0, no trigger)
  35-41 saturated tropical          (LCL at/near launch level)

Three configurations sweep the zm_conv_main branches (zm_microp is
ALWAYS .false. -- see PORTING_PLAN.md scope note; the EAMv3 default
zmconv_microp=.true. path needs zm_microphysics.F90 which cannot be
built standalone):
  a: EAMv3 phys="default" params, is_first_step=T (CAPE>70 trigger,
     no DCAPE second call)
  b: same params, is_first_step=F (DCAPE trigger: second
     compute_dilute_cape call on t_star/q_star with frozen launch
     level; gather needs cape>0 AND dcape>0)
  c: legacy flags, is_first_step=F: trig_dcape/trig_ull/clos_dyn_adj/
     tpert_fix off, num_cin=5, mx_bot_lyr_adj=0, no_deep_pbl=T,
     pre-DCAPE tunings (dmpdz=-1e-3, tiedke_add=0.5, tpert_fac=0,
     c0=0.003, alfa=0.1, ke=3e-6, old_snow=T)

EAMv3 phys="default" values from bld/namelist_files/
namelist_defaults_eam.xml: zmconv_tau=3600 (base default; no
phys-default override), zmconv_alfa=0.14, zmconv_ke=2.5e-6 (dyn=se,
microphys=p3), zmconv_dmpdz=-0.7e-3, zmconv_tpert_fix=.true.,
zmconv_tp_fac=2.0, zmconv_tiedke_add=0.8, zmconv_c0_lnd/ocn=0.0020,
zmconv_cape_cin=1, zmconv_mx_bot_lyr_adj=1, zmconv_trig_dcape=.true.,
zmconv_trig_ull=.true., zmconv_clos_dyn_adj=.true.; old_snow=.false.
(zm_conv_intr sets it when zm_microp; only used by zm_conv_evap which
is not driven). no_deep_pbl=.false. (phys_deepconv_pbl default).
limcnv follows the 40 hPa reference-interface rule (zm_conv_intr).

All level/column indices stored 1-based exactly as Fortran returns
them. Gathered-index outputs (msemax_klev_g, jt, mflx_up, entr_up,
detr_up, mflx_dn, entr_dn, p_del, dsubcld, ql) are stored in gathered
row order; scattered outputs (jctop, jcbot, prec, heat, qtnd, mcon,
pflx, zdu, rprd, dlf, cape, dcape, rliq) in column order.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_zm_conv_f  # noqa: E402

d = eam_zm_conv_f.zm_conv_driver
d.drv_init()

rng = np.random.default_rng(20260712)
nfam, npc = 6, 7
ncol, nlev = nfam * npc, 72
RD, GRAV = 287.042, 9.80616
TIME_STEP = 1800.0

# --- pressure grid [Pa] (zm_cape golden recipe) ---
ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
pint1 = 225.5 + ai * (1.0e5 - 225.5)
pint = np.broadcast_to(pint1, (ncol, nlev + 1)).copy()
pint *= (1 + 0.03 * rng.uniform(-1, 1, (ncol, 1)))
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
pdel = pint[:, 1:] - pint[:, :-1]
ps = pmid[:, -1:]

# limcnv: first 1-based interface k with pref_edge(k) < 40 hPa <=
# pref_edge(k+1) (zm_conv_intr.F90); with 0-based j = argmax(pint1 >=
# 4000 Pa) the crossing is at 1-based k = j.
limcnv = int(np.argmax(pint1 >= 4000.0))


def qs_bolton(T, p):
    """Generator-side approximate qsat (input construction only)."""
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    return 0.622 * es / (p - es)


T = np.zeros((ncol, nlev))
rh = np.zeros((ncol, nlev))
sig = pmid / ps
for f, (ts, ex, rh0) in enumerate([
        (303.0, 0.19, 0.85),   # strongly unstable tropical
        (296.0, 0.16, 0.70),   # marginally unstable
        (282.0, 0.10, 0.50),   # stable midlatitude
        (291.0, 0.18, 0.80),   # midlatitude convective
        (312.0, 0.21, 0.05),   # high-CAPE dry boundary layer
        (300.0, 0.19, 0.99)]):  # saturated tropical
    rows = slice(f * npc, (f + 1) * npc)
    tsv = ts + rng.uniform(-1.0, 2.0, (npc, 1))
    T[rows] = np.maximum(tsv * sig[rows] ** ex, 195.0)
    rh[rows] = np.clip(rh0 * sig[rows] ** 0.8, 0.02, rh0)

T += rng.uniform(-0.5, 0.5, T.shape)
q = np.maximum(rh * qs_bolton(T, pmid), 1.0e-9)
q *= 1 + 0.05 * rng.uniform(-1, 1, q.shape)

# hydrostatic heights ABOVE THE SURFACE (the driver takes z_mid/z_int
# relative to the surface plus geos separately, as zm_conv_tend does)
tv = T * (1 + 0.608 * q)
zint = np.zeros((ncol, nlev + 1))
zmid = np.zeros((ncol, nlev))
for k in range(nlev - 1, -1, -1):
    zmid[:, k] = zint[:, k + 1] + RD * tv[:, k] / GRAV \
        * np.log(pint[:, k + 1] / pmid[:, k])
    zint[:, k] = zint[:, k + 1] + RD * tv[:, k] / GRAV \
        * np.log(pint[:, k + 1] / pint[:, k])

z_srf = rng.uniform(0.0, 500.0, (ncol,))
geos = GRAV * z_srf

pbl_hgt = np.where(np.arange(ncol) % 3 == 0, 2200.0, 1000.0)
tpert = np.where(np.arange(ncol) % 2 == 0, 0.0,
                 rng.uniform(0.2, 1.5, ncol))
landfrac = np.tile([0.0, 1.0, 0.37, 1.0, 0.0, 1.0, 0.5],
                   nfam)[:ncol].astype(np.float64)
omega = rng.uniform(-0.6, 0.6, (ncol, nlev))

# perturbed previous-step state for the DCAPE trigger path
t_star = T + rng.uniform(-0.5, 0.5, T.shape)
q_star = np.maximum(q * (1 + 0.03 * rng.uniform(-1, 1, q.shape)), 1.0e-9)

# parameter sets: EAMv3 phys="default" (a, b) and legacy (c)
P_V3 = dict(tau=3600.0, alfa=0.14, ke=2.5e-6, dmpdz=-0.7e-3,
            tpert_fix=1, tpert_fac=2.0, tiedke_add=0.8,
            c0_lnd=0.0020, c0_ocn=0.0020, num_cin=1,
            mx_bot_lyr_adj=1, trig_dcape=1, trig_ull=1,
            clos_dyn_adj=1, old_snow=0, no_deep_pbl=0)
P_LEG = dict(tau=3600.0, alfa=0.1, ke=3.0e-6, dmpdz=-1.0e-3,
             tpert_fix=0, tpert_fac=0.0, tiedke_add=0.5,
             c0_lnd=0.0030, c0_ocn=0.0030, num_cin=5,
             mx_bot_lyr_adj=0, trig_dcape=0, trig_ull=0,
             clos_dyn_adj=0, old_snow=1, no_deep_pbl=1)
CONFIGS = {"a": (1, P_V3), "b": (0, P_V3), "c": (0, P_LEG)}

F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731
NAMES = ("lengath", "gather_index", "msemax_klev_g", "jctop", "jcbot",
         "jt", "prec", "heat", "qtnd", "cape", "dcape", "mcon", "pflx",
         "zdu", "mflx_up", "entr_up", "detr_up", "mflx_dn", "entr_dn",
         "p_del", "dsubcld", "ql", "rliq", "rprd", "dlf")

out = dict(t=T, q=q, omega=omega, pmid=pmid, pint=pint, pdel=pdel,
           geos=geos, zmid=zmid, zint=zint, pbl_hgt=pbl_hgt,
           tpert=tpert, landfrac=landfrac, t_star=t_star,
           q_star=q_star, limcnv=np.array(limcnv),
           time_step=np.array(TIME_STEP))

for cfg, (first, p) in CONFIGS.items():
    res = d.drv_zm_conv_main(
        TIME_STEP, first, limcnv, p["no_deep_pbl"],
        p["tau"], p["alfa"], p["ke"], p["dmpdz"], p["tpert_fix"],
        p["tpert_fac"], p["tiedke_add"], p["c0_lnd"], p["c0_ocn"],
        p["num_cin"], p["mx_bot_lyr_adj"], p["trig_dcape"],
        p["trig_ull"], p["clos_dyn_adj"], p["old_snow"],
        F(T), F(q), F(omega), F(pmid), F(pint), F(pdel), geos,
        F(zmid), F(zint), pbl_hgt, tpert, landfrac,
        F(t_star), F(q_star))
    r = dict(zip(NAMES, res))
    out.update({f"{n}_{cfg}": np.asarray(v) for n, v in r.items()})
    ng = int(r["lengath"])
    print(f"cfg {cfg}: lengath={ng} "
          f"gathered={np.sort(r['gather_index'][:ng])}")
    print(f"  cape={np.array2string(r['cape'], precision=0)}")
    print(f"  prec[mm/day]={np.array2string(r['prec'] * 8.64e7, precision=1)}")

# sanity: config a triggers most strongly-unstable tropical columns
# with positive precip and heating aloft / drying below
ra = {n: out[f"{n}_a"] for n in NAMES}
gath_a = set(ra["gather_index"][:int(ra["lengath"])] - 1)
trop = [i for i in range(npc) if i in gath_a]
assert len(trop) >= 5, "expected most tropical columns to trigger"
for i in trop:
    assert ra["prec"][i] > 0.0
    assert ra["heat"][i, 30:55].max() > 0.0, "heating aloft expected"
    assert ra["qtnd"][i, 60:].min() < 0.0, "low-level drying expected"
stable = range(2 * npc, 3 * npc)
assert not (set(stable) & gath_a), "stable family must not trigger"
assert np.all(ra["prec"][list(stable)] == 0.0)
dry = range(4 * npc, 5 * npc)
assert np.all(ra["cape"][list(dry)] == 0.0), \
    "dry-BL family should hit the 600 hPa LCL threshold (CAPE=0)"
assert int(out["lengath_b"]) > 0 and int(out["lengath_c"]) > 0

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "zm_conv zm_conv_main (ZM deep convection core, "
                  "zm_microp=.false.)",
        "source_sha": sha, "ncol": ncol, "nlev": nlev, "npc": npc,
        "families": ["tropical_strong", "marginal", "stable_midlat",
                     "midlat_convective", "dry_boundary_layer",
                     "saturated_tropical"],
        "limcnv": limcnv, "msg": limcnv - 1,
        "time_step": TIME_STEP, "index_base": 1,
        "zm_microp": False, "mcsp_enabled": False,
        "params": {"a": {**P_V3, "is_first_step": 1, "limcnv": limcnv},
                   "b": {**P_V3, "is_first_step": 0, "limcnv": limcnv},
                   "c": {**P_LEG, "is_first_step": 0, "limcnv": limcnv}},
        "zm_const": {"zvir": 1.608, "grav": 9.80616, "rdair": RD,
                     "note": "filled by real zm_const_set_to_global"}}
gold = Path(__file__).resolve().parents[1] / "golden" / "zm_conv_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
