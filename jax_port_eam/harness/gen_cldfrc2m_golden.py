#!/usr/bin/env python3
"""Tier-1 golden generator for the cldfrc2m cloud-fraction module (run
in the scream-dev container).

Contents:
 1. Profile set: 40 columns x 72 levels (flattened; every point is
    independent) spanning all three pressure bands (p >= premib,
    premit <= p < premib, p < premit), land/ocean/fractional landfrac,
    snow depths straddling the 1e-6 threshold, RH from 0 through every
    triangular-PDF branch to supersaturated, qi from 0 through the
    minice / qist_min / qist_max limiter regimes, ni over 7 decades.
    Outputs: astG_PDF, astG_RHU (a, G) and aist for every iceopt 1-7
    (re-initializing the module between iceopts) plus an alternate
    (rhmaxi, rhmini) pair for the default iceopt=5.
 2. Dense 1-D single-point sweeps of astG_PDF_single/astG_RHU_single
    (U in [0, 1.2] x 6 pressures x 5 surface types, recording a, G and
    orhmin) and of aist_single at iceopt=5 (qv, qi, and T sweeps).

Params are the EAMv3 phys="default" values from
bld/namelist_files/namelist_defaults_eam.xml:
  cldfrc_rhminl=0.950 (microphys="p3"), cldfrc_rhminl_adj_land=0.100,
  cldfrc_rhminh=0.800, cldfrc_premit=25000 Pa (dyn="se"),
  cldfrc_premib=70000 Pa, cldfrc_iceopt=5, cldfrc_icecrit=0.93,
  cldfrc_minice=1.0e-12, cldfrc2m_rhmini=0.80, cldfrc2m_rhmaxi=1.05
  (clubb_sgs="1" clubb_do_deep="0").
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_cldfrc2m_f  # noqa: E402

d = eam_cldfrc2m_f.cldfrc2m_driver

# EAMv3 phys="default" parameters
RHMINL, RHMINL_ADJ_LAND, RHMINH = 0.950, 0.100, 0.800
PREMIT, PREMIB = 25000.0, 70000.0
ICEOPT, ICECRIT, MINICE = 5, 0.93, 1.0e-12
RHMINI, RHMAXI = 0.80, 1.05
RHMINI_ALT, RHMAXI_ALT = 0.70, 1.00


def init(iceopt=ICEOPT):
    d.drv_init(RHMINL, RHMINL_ADJ_LAND, RHMINH, PREMIT, PREMIB,
               iceopt, ICECRIT, MINICE)


rng = np.random.default_rng(20260712)
ncol, nlev = 40, 72

# ---------------------------------------------------------------- #
# 1. profile set                                                    #
# ---------------------------------------------------------------- #
ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
pint1 = 225.5 + ai * (1.0e5 - 225.5)                     # Pa
pint = np.broadcast_to(pint1, (ncol, nlev + 1)).copy()
pint *= (1 + 0.03 * rng.uniform(-1, 1, (ncol, 1)))
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
sig = pmid / pmid[:, -1:]

T = np.zeros((ncol, nlev))
for i in range(ncol):
    ts = rng.uniform(255.0, 305.0)
    ex = rng.uniform(0.10, 0.21)
    T[i] = np.maximum(ts * sig[i] ** ex, 180.0)
T += rng.uniform(-1.0, 1.0, T.shape)

# relative humidity U: smooth base + structured slabs guaranteeing
# every astG branch in every pressure band
U = np.clip(rng.uniform(0.6, 0.95) * sig ** 0.6
            + rng.uniform(-0.15, 0.15, (ncol, nlev)), 0.0, 1.2)
U[0:4] = np.linspace(0.90, 1.10, nlev)     # dense upper-branch coverage
U[4:8] = np.linspace(0.70, 1.00, nlev)
U[8:10] = np.linspace(0.0, 1.2, nlev)
U[10, :] = 1.0                             # exactly saturated
U[11, ::2] = 0.95                          # exactly rhminl
U[11, 1::2] = 0.80                         # exactly rhminh
U[12, :] = 0.85                            # exactly rhminl - adj_land

landfrac = rng.choice([0.0, 1.0, 0.5, 0.3, 0.7], size=ncol,
                      p=[0.35, 0.35, 0.1, 0.1, 0.1])
snowh = rng.choice([0.0, 1.0e-7, 1.0e-6, 2.0e-6, 0.1], size=ncol,
                   p=[0.4, 0.15, 0.15, 0.15, 0.15])
lf2 = np.broadcast_to(landfrac[:, None], (ncol, nlev)).copy()
sh2 = np.broadcast_to(snowh[:, None], (ncol, nlev)).copy()


def qs_bolton(t, p):
    """Generator-side approximate qsat (input construction only)."""
    es = 611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65))
    return 0.622 * es / np.maximum(p - es, 10.0)


qv = np.maximum(U * qs_bolton(T, pmid), 1.0e-9)

# ice mass / number: cover 0, sub-minice, minice, limiter regimes
qi = 10.0 ** rng.uniform(-14.0, -2.3, (ncol, nlev))
qi.ravel()[::17] = 0.0
qi.ravel()[3::17] = 0.5e-12                # below minice
qi.ravel()[7::17] = MINICE                 # exactly minice
qi.ravel()[11::17] = 10.0 ** rng.uniform(-8.0, -7.0,
                                         qi.ravel()[11::17].shape)
ni = np.maximum(10.0 ** rng.uniform(0.0, 7.0, (ncol, nlev)), 1.0)

flat = lambda x: np.ascontiguousarray(x, dtype=np.float64).ravel()  # noqa: E731
uu, pp, qq, tt = flat(U), flat(pmid), flat(qv), flat(T)
qqi, nni, ll, ss = flat(qi), flat(ni), flat(lf2), flat(sh2)

out = dict(u=uu, p=pp, qv=qq, t=tt, qi=qqi, ni=nni,
           landfrac=ll, snowh=ss)

init()
out["a_pdf"], out["ga_pdf"] = d.drv_astg_pdf(uu, pp, qq, ll, ss)
out["a_rhu"], out["ga_rhu"] = d.drv_astg_rhu(uu, pp, qq, ll, ss)

for opt in (1, 2, 3, 4, 5, 6, 7):
    init(opt)
    out[f"aist_opt{opt}"] = d.drv_aist_vector(
        qq, tt, pp, qqi, nni, ll, ss, RHMAXI, RHMINI)
init()
out["aist_opt5_alt"] = d.drv_aist_vector(
    qq, tt, pp, qqi, nni, ll, ss, RHMAXI_ALT, RHMINI_ALT)

# ---------------------------------------------------------------- #
# 2. dense single-point sweeps                                      #
# ---------------------------------------------------------------- #
u_sweep = np.concatenate([np.linspace(0.0, 1.2, 481),
                          [0.80, 0.85, 0.95, 1.0, 1.05,
                           1.0 - 0.05 / 6.0, 1.0 - 0.2 / 6.0]])
p_cases = [95000.0, PREMIB, 50000.0, PREMIT, 24999.9, 10000.0]
surf_cases = [(0.0, 0.0), (1.0, 0.0), (1.0, 0.5),
              (0.5, 1.0e-6), (0.5, 2.0e-6)]

nu, npr, ns = u_sweep.size, len(p_cases), len(surf_cases)
for name, fn in (("pdf", d.drv_astg_pdf_single),
                 ("rhu", d.drv_astg_rhu_single)):
    a = np.zeros((npr, ns, nu))
    ga = np.zeros((npr, ns, nu))
    rhmin = np.zeros((npr, ns, nu))
    for ip, p in enumerate(p_cases):
        for isf, (lf, sh) in enumerate(surf_cases):
            for iu, u in enumerate(u_sweep):
                a[ip, isf, iu], ga[ip, isf, iu], rhmin[ip, isf, iu] = \
                    fn(u, p, 5.0e-3, lf, sh)
    out[f"sweep_a_{name}"] = a
    out[f"sweep_ga_{name}"] = ga
    out[f"sweep_rhmin_{name}"] = rhmin
out["sweep_u"] = u_sweep
out["sweep_p"] = np.array(p_cases)
out["sweep_landfrac"] = np.array([s[0] for s in surf_cases])
out["sweep_snowh"] = np.array([s[1] for s in surf_cases])

# aist_single sweeps (iceopt=5 default)
qv_sw = np.linspace(0.0, 6.0e-4, 301)
out["aist_qv_sweep_in"] = qv_sw
out["aist_qv_sweep"] = np.array(
    [d.drv_aist_single(q, 230.0, 30000.0, 1.0e-5, 0.0, 0.0,
                       RHMAXI, RHMINI) for q in qv_sw])
qi_sw = np.concatenate([[0.0], 10.0 ** np.linspace(-14.0, -2.0, 241)])
out["aist_qi_sweep_in"] = qi_sw
out["aist_qi_sweep"] = np.array(
    [d.drv_aist_single(2.0e-4, 230.0, 30000.0, q, 0.0, 0.0,
                       RHMAXI, RHMINI) for q in qi_sw])
t_sw = np.linspace(180.0, 273.0, 187)
out["aist_t_sweep_in"] = t_sw
out["aist_t_sweep"] = np.array(
    [d.drv_aist_single(1.0e-4, t, 30000.0, 1.0e-5, 0.0, 0.0,
                       RHMAXI, RHMINI) for t in t_sw])

# ---- sanity: RH near 1 -> fraction near 1; dry -> 0 ----
assert out["a_pdf"][uu >= 1.0].min() == 1.0
assert out["a_pdf"][uu < 0.5].max() == 0.0
assert out["aist_opt5"][qqi < MINICE].max() == 0.0
assert 0.0 <= out["aist_opt5"].max() <= 0.999
print("a_pdf mean/max:", out["a_pdf"].mean(), out["a_pdf"].max())
print("aist_opt5 mean/max:", out["aist_opt5"].mean(),
      out["aist_opt5"].max())
for opt in (1, 2, 3, 4, 5, 6, 7):
    print(f"aist_opt{opt} nonzero frac:",
          (out[f"aist_opt{opt}"] > 0).mean())

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "cldfrc2m (astG_PDF/astG_RHU/aist, single+vector)",
        "source_sha": sha, "ncol": ncol, "nlev": nlev,
        "params": {"rhminl": RHMINL,
                   "rhminl_adj_land": RHMINL_ADJ_LAND,
                   "rhminh": RHMINH, "premit": PREMIT,
                   "premib": PREMIB, "icecrit": ICECRIT,
                   "minice": MINICE, "rhmini": RHMINI,
                   "rhmaxi": RHMAXI, "rhmini_alt": RHMINI_ALT,
                   "rhmaxi_alt": RHMAXI_ALT, "iceopt_default": ICEOPT,
                   "CAMstfrac": False, "freeze_dry": False},
        "notes": ("rhmini/rhmaxi passed explicitly per call (readnl "
                  "path is masterproc-only in the harness); other "
                  "params set via the cloud_fraction stub + real "
                  "cldfrc2m_init. aist single sweeps use iceopt=5; "
                  "vector goldens cover iceopt 1-7.")}
gold = Path(__file__).resolve().parents[1] / "golden" / \
    "cldfrc2m_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
