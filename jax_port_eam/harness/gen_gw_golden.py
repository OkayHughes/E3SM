#!/usr/bin/env python3
"""Tier-1 golden generator for the orographic gravity-wave spine.

Chained realistic states: hydrostatic 72L columns with jet-like wind
profiles, orography sgh from 0 (ocean) to 1500 m, alpha profile
qualitatively like the gw_drag file (decays from top). Captures
gw_prof, gw_oro_src and gw_drag_prof(ngwv=0) inputs/outputs.
EAMv3 defaults: kwv=6.28e-5, fcrit2=1.0, effgw_oro=0.375, ktop=0,
kbotbg = last interface above 500 hPa, orographic_only distinguishes
the two tendency-limiter code paths -- both captured.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_gw_f  # noqa: E402

d = eam_gw_f.gw_driver

rng = np.random.default_rng(2026)
ncol, nlev = 24, 72
GRAVIT, RAIR, CPAIR = 9.80616, 287.042, 1004.64

ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
pint = 225.5 + ai * (1.0e5 - 225.5)
pint = np.broadcast_to(pint, (ncol, nlev + 1)).copy()
pint *= (1 + 0.03 * rng.uniform(-1, 1, (ncol, 1)))
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
dpm = np.diff(pint, axis=1)
rdpm = 1.0 / dpm
piln = np.log(pint)

t = 216.0 + 72.0 * (pmid / pmid[:, -1:]) ** 0.32
t += 8.0 * np.exp(-((np.log(pmid / 1e2)) ** 2))     # warm stratopause bump
t += rng.uniform(-3.0, 3.0, t.shape)

# geopotential-ish midpoint heights (hydrostatic, isothermal-layer)
zm = (RAIR * t / GRAVIT) * np.log(pint[:, -1:] / pmid)

# jet wind profile + shear + some reversed columns (critical levels)
zj = 1.1e4
u = 25.0 * np.exp(-((zm - zj) / 8e3) ** 2) + rng.uniform(-5, 5, t.shape)
u += rng.uniform(-2, 12, (ncol, 1))
v = 6.0 * np.sin(zm / 6e3) + rng.uniform(-3, 3, t.shape)
u[::5] *= -0.6          # some easterly columns
sgh = np.concatenate([np.zeros(4), rng.uniform(5.0, 1500.0, ncol - 4)])
lat = np.deg2rad(rng.uniform(-80, 80, ncol))

# Newtonian cooling profile qualitatively like the gw_drag file
alpha = 1.0e-6 + 3.0e-6 * np.exp(-np.arange(nlev + 1) / 8.0)

kbotbg = int(np.argmax(pint[0] >= 5.0e4)) - 1  # last interface < 500 hPa
dt = 1800.0
effgw = 0.375

out = dict(pmid=pmid, pint=pint, dpm=dpm, rdpm=rdpm, piln=piln, t=t,
           u=u, v=v, zm=zm, sgh=sgh, lat=lat, alpha=alpha,
           kbotbg=np.array(kbotbg))

for oro_only in (0, 1):
    d.drv_gw_init(0, kbotbg, 1.0, 6.28e-5, GRAVIT, RAIR, alpha, oro_only)
    rhoi, ti, nm, ni = d.drv_gw_prof(CPAIR, t, pmid, pint)
    src_level, tend_level, tau0, ubm, ubi, xv, yv = d.drv_gw_oro_src(
        u, v, t, sgh, pmid, pint, dpm, zm, nm)
    tau, utgw, vtgw = d.drv_gw_drag_prof_oro(
        src_level, tend_level, dt, lat, t, ti, pmid, pint, dpm, rdpm,
        piln, rhoi, nm, ni, ubm, ubi, xv, yv, effgw, tau0)
    sfx = f"_oro{oro_only}"
    out.update({f"rhoi{sfx}": rhoi, f"ti{sfx}": ti, f"nm{sfx}": nm,
                f"ni{sfx}": ni, f"src_level{sfx}": src_level,
                f"tend_level{sfx}": tend_level, f"tau0{sfx}": tau0,
                f"ubm{sfx}": ubm, f"ubi{sfx}": ubi, f"xv{sfx}": xv,
                f"yv{sfx}": yv, f"tau{sfx}": tau, f"utgw{sfx}": utgw,
                f"vtgw{sfx}": vtgw})

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "gw_prof + gw_oro_src + gw_drag_prof(ngwv=0)",
        "source_sha": sha, "ncol": ncol, "nlev": nlev, "dt": dt,
        "params": {"kwv": 6.28e-5, "fcrit2": 1.0, "effgw_oro": effgw,
                   "ktop": 0, "kbotbg": kbotbg, "gravit": GRAVIT,
                   "rair": RAIR, "cpair": CPAIR}}
gold = Path(__file__).resolve().parents[1] / "golden" / "gw_oro_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
