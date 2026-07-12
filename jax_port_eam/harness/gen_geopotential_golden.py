#!/usr/bin/env python3
"""Tier-1 golden generator for geopotential_t / geopotential_dse."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_geopotential_f  # noqa: E402

d = eam_geopotential_f.geopotential_driver
rng = np.random.default_rng(7)
ncol, nlev = 16, 72
GRAVIT, RAIR, CPAIR, ZVIR = 9.80616, 287.042, 1004.64, 0.60779262

ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
pint = 225.5 + ai * (1.0e5 - 225.5)
pint = np.broadcast_to(pint, (ncol, nlev + 1)).copy()
pint *= (1 + 0.04 * rng.uniform(-1, 1, (ncol, 1)))
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
pdel = np.diff(pint, axis=1)
rpdel = 1.0 / pdel
piln, pmln = np.log(pint), np.log(pmid)

t = 216.0 + 72.0 * (pmid / pmid[:, -1:]) ** 0.32 \
    + rng.uniform(-4, 4, (ncol, nlev))
q = np.clip(rng.uniform(0.2, 1.0, (ncol, nlev))
            * 2e-2 * (pmid / 1e5) ** 3, 1e-9, None)
phis = rng.uniform(0.0, 3.0e4, ncol)
rair = np.full((ncol, nlev), RAIR)
cpair = np.full((ncol, nlev), CPAIR)
zvir = np.full((ncol, nlev), ZVIR)

zi_t, zm_t = d.drv_geopotential_t(piln, pmln, pint, pmid, pdel, rpdel,
                                  t, q, rair, GRAVIT, zvir)

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "geopotential_t (SE branch, fvdyn=false; "
                  "geopotential_dse is private/'do not use')",
        "source_sha": sha,
        "params": {"gravit": GRAVIT, "rair": RAIR, "cpair": CPAIR,
                   "zvir": ZVIR}}
gold = Path(__file__).resolve().parents[1] / "golden" / \
    "geopotential_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta),
                    piln=piln, pmln=pmln, pint=pint, pmid=pmid,
                    pdel=pdel, rpdel=rpdel, t=t, q=q, phis=phis,
                    zi_t=zi_t, zm_t=zm_t)
print(f"wrote {gold}")
