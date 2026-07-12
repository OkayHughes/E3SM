#!/usr/bin/env python3
"""Tier-1 golden generator for the ZM dilute CAPE core (run in the
scream-dev container).

Synthetic hydrostatic 72-level columns, 6 families x 7 columns = 42:
  0-6   tropical moist unstable  (CAPE O(1000) J/kg)
  7-13  midlatitude stable       (small/zero CAPE)
  14-20 dry desert               (LCL above 600 hPa threshold -> CAPE 0)
  21-27 saturated tropical       (LCL at/near launch level)
  28-34 inversion-capped marine  (negative buoyancy above capping inv.)
  35-41 cold polar stable        (tmix <= tfreez: freezing branch from
                                  the first ascent level)

Four flag configurations sweep every branch of compute_dilute_cape:
  a: EAMv3 defaults, first call  (calc_msemax_klev=T, trig_ull/dcape=T,
     tpert_fix=T, num_cin=1)     -> ULL launch-level search
  b: dcape second call           (calc_msemax_klev=F, prev launch level
     from (a), perturbed t_star/q_star state)
  c: legacy flags                (trig_ull=F -> pblt-capped search,
     tpert_fix=F, trig_dcape=F, num_cin=5)
  d: use_input_tq_mx=T           (prev level + q_mx/t_mx from (a) applied
     to the perturbed state)

Params: EAMv3 phys="default" from namelist_defaults_eam.xml --
dmpdz=-0.7e-3, tiedke_add=0.8, tpert_fac=2.0 (zmconv_tp_fac),
mx_bot_lyr_adj=1, cape_cin=1. num_msg = limcnv-1 with limcnv from the
40 hPa reference-interface rule in zm_conv_intr.F90.
All level indices stored 1-based exactly as the Fortran returns them.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_zm_cape_f  # noqa: E402

d = eam_zm_cape_f.zm_cape_driver
d.drv_init()

rng = np.random.default_rng(20260712)
nfam, npc = 6, 7
ncol, nlev = nfam * npc, 72
RD, GRAV = 287.042, 9.80616

# --- pressure grid (gen_gw_golden.py recipe), in hPa as zm_conv uses ---
ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
pint1 = (225.5 + ai * (1.0e5 - 225.5)) * 1e-2          # hPa
pint = np.broadcast_to(pint1, (ncol, nlev + 1)).copy()
pint *= (1 + 0.03 * rng.uniform(-1, 1, (ncol, 1)))
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
ps = pmid[:, -1:]

# num_msg = limcnv - 1; limcnv = first interface k (1-based) with
# pref_edge(k) < 40 hPa <= pref_edge(k+1)  (zm_conv_intr.F90). With
# 0-based j = argmax(pint1 >= 40), the crossing is at 1-based k = j.
limcnv = int(np.argmax(pint1 >= 40.0))
num_msg = limcnv - 1


def qs_bolton(T, p):
    """Generator-side approximate qsat (only used to construct inputs)."""
    es = 6.112 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    return 0.622 * es / (p - es)


T = np.zeros((ncol, nlev))
rh = np.zeros((ncol, nlev))
sig = pmid / ps
for f, (ts, ex, rh0) in enumerate([
        (302.0, 0.19, 0.85),   # tropical unstable
        (282.0, 0.11, 0.50),   # midlat stable
        (312.0, 0.21, 0.05),   # dry desert
        (300.0, 0.19, 0.99),   # saturated tropical
        (290.0, 0.13, 0.90),   # inversion-capped marine
        (262.0, 0.10, 0.70)]):  # cold polar stable
    rows = slice(f * npc, (f + 1) * npc)
    tsv = ts + rng.uniform(-1.0, 2.0, (npc, 1))
    T[rows] = np.maximum(tsv * sig[rows] ** ex, 195.0)
    rh[rows] = np.clip(rh0 * sig[rows] ** 0.8, 0.02, rh0)

# capping inversion for the marine family: warm bump 900-850 hPa,
# dry air above the inversion
inv = slice(4 * npc, 5 * npc)
bump = 8.0 * np.exp(-((pmid[inv] - 870.0) / 30.0) ** 2)
T[inv] += bump
rh[inv] = np.where(pmid[inv] < 850.0, 0.30 * sig[inv] ** 0.5, rh[inv])

T += rng.uniform(-0.7, 0.7, T.shape)
q = np.maximum(rh * qs_bolton(T, pmid), 1.0e-9)
q *= 1 + 0.05 * rng.uniform(-1, 1, q.shape)

# hydrostatic heights (surface elevation 0-500 m)
z_srf = rng.uniform(0.0, 500.0, (ncol,))
tv = T * (1 + 0.608 * q)
zint = np.zeros((ncol, nlev + 1))
zmid = np.zeros((ncol, nlev))
zint[:, nlev] = z_srf
for k in range(nlev - 1, -1, -1):
    zmid[:, k] = zint[:, k + 1] + RD * tv[:, k] / GRAV \
        * np.log(pint[:, k + 1] / pmid[:, k])
    zint[:, k] = zint[:, k + 1] + RD * tv[:, k] / GRAV \
        * np.log(pint[:, k + 1] / pint[:, k])

# PBL top index (1-based): ~1 km AGL, some columns deeper (2.2 km) so
# tpert_fix zeroes tpert when the ULL launch sits above the PBL top
pbl_hgt = np.where(np.arange(ncol) % 3 == 0, 2200.0, 1000.0)
pblt = np.array([int(np.argmin(np.abs(zmid[i] - z_srf[i] - pbl_hgt[i]))) + 1
                 for i in range(ncol)], dtype=np.int32)
tpert = np.where(np.arange(ncol) % 2 == 0, 0.0,
                 rng.uniform(0.2, 1.5, ncol))

# perturbed state for the dcape second call (t_star/q_star analogue)
t_star = T + rng.uniform(-0.5, 0.5, T.shape)
q_star = np.maximum(q * (1 + 0.03 * rng.uniform(-1, 1, q.shape)), 1.0e-9)

# EAMv3 phys="default" parameters
DMPDZ, TIEDKE_ADD, TPERT_FAC, MX_BOT_LYR_ADJ = -0.7e-3, 0.8, 2.0, 1

F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731
dummy_prev = np.ones(ncol, dtype=np.int32)
zeros = np.zeros(ncol)


def run(qa, ta, num_cin, tpert_fix, trig_ull, trig_dcape, calc,
        prev, use_input, q_mx_in, t_mx_in):
    return d.drv_compute_dilute_cape(
        num_cin, num_msg, F(qa), F(ta), F(zmid), F(pmid), F(pint),
        pblt, tpert, DMPDZ, TIEDKE_ADD, TPERT_FAC, MX_BOT_LYR_ADJ,
        tpert_fix, trig_ull, trig_dcape, calc,
        np.asarray(prev, dtype=np.int32), use_input, q_mx_in, t_mx_in)


out = dict(q=q, t=T, zmid=zmid, pmid=pmid, pint=pint, pblt=pblt,
           tpert=tpert, t_star=t_star, q_star=q_star,
           num_msg=np.array(num_msg))

names = ("parcel_temp", "parcel_qsat", "msemax_klev", "lcl_temperature",
         "lcl_klev", "eql_klev", "cape", "q_mx", "t_mx")

# a: EAMv3 default first call
res_a = run(q, T, 1, 1, 1, 1, 1, dummy_prev, 0, zeros, zeros)
out.update({f"{n}_a": v for n, v in zip(names, res_a)})
# b: dcape second call, launch level frozen to (a)'s, perturbed state
res_b = run(q_star, t_star, 1, 1, 1, 1, 0, res_a[2], 0, zeros, zeros)
out.update({f"{n}_b": v for n, v in zip(names, res_b)})
# c: legacy flags (pblt-capped search, no tpert_fix, num_cin=5)
res_c = run(q, T, 5, 0, 0, 0, 1, dummy_prev, 0, zeros, zeros)
out.update({f"{n}_c": v for n, v in zip(names, res_c)})
# d: use_input_tq_mx with (a)'s launch level and launch-level T/q
res_d = run(q_star, t_star, 1, 1, 1, 1, 1, res_a[2], 1,
            res_a[7], res_a[8])
out.update({f"{n}_d": v for n, v in zip(names, res_d)})

print("cape_a:", np.array2string(res_a[6], precision=1))
print("cape_b:", np.array2string(res_b[6], precision=1))
print("cape_c:", np.array2string(res_c[6], precision=1))
print("cape_d:", np.array2string(res_d[6], precision=1))
assert res_a[6].max() > 500.0, "expected O(1000) tropical CAPE"
assert (res_a[6][2 * npc:3 * npc] == 0.0).all(), \
    "desert columns should hit the 600 hPa LCL threshold (CAPE=0)"

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "zm_conv_cape compute_dilute_cape (dilute CAPE core)",
        "source_sha": sha, "ncol": ncol, "nlev": nlev,
        "families": ["tropical", "midlat_stable", "desert",
                     "saturated", "inversion", "polar"],
        "npc": npc, "num_msg": num_msg, "limcnv": limcnv,
        "index_base": 1,
        "params": {"dmpdz": DMPDZ, "tiedke_add": TIEDKE_ADD,
                   "tpert_fac": TPERT_FAC,
                   "mx_bot_lyr_adj": MX_BOT_LYR_ADJ,
                   "num_cin": {"a": 1, "b": 1, "c": 5, "d": 1},
                   "tpert_fix": {"a": 1, "b": 1, "c": 0, "d": 1},
                   "trig_ull": {"a": 1, "b": 1, "c": 0, "d": 1},
                   "trig_dcape": {"a": 1, "b": 1, "c": 0, "d": 1},
                   "calc_msemax_klev": {"a": 1, "b": 0, "c": 1, "d": 1},
                   "use_input_tq_mx": {"a": 0, "b": 0, "c": 0, "d": 1},
                   "state": {"a": "t/q", "b": "t_star/q_star",
                             "c": "t/q", "d": "t_star/q_star"},
                   "zvir": 1.608, "grav": 9.80616, "rdair": RD}}
gold = Path(__file__).resolve().parents[1] / "golden" / "zm_cape_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
