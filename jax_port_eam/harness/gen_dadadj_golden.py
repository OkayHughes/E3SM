#!/usr/bin/env python3
"""Tier-1 golden generator for dadadj.

States: hydrostatic 72-level columns with a realistic stratospheric
lapse structure, then the top-nlvdry region perturbed into
superadiabatic (unstable) configurations of varying strength so the
sweep exercises: no-op (stable), single-pair adjustment, multi-pair
cascades, and multiple outer iterations. nlvdry is swept over {3, 8}.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_dadadj_f  # noqa: E402

d = eam_dadadj_f.dadadj_driver

rng = np.random.default_rng(123)
ncol, nlev = 16, 72

# hybrid-ish pressure grid, ptop ~ 1 hPa
ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.6
pint = 100.0 + ai * (1.0e5 - 100.0)          # (nlev+1,)
pint = np.broadcast_to(pint, (ncol, nlev + 1)).copy()
pint *= (1 + 0.02 * rng.uniform(-1, 1, (ncol, 1)))
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
pdel = np.diff(pint, axis=1)

out = {"pmid": pmid, "pint": pint, "pdel": pdel}

for nlvdry in (3, 8):
    # smooth reference temperature profile (warm strat top -> cold tropo)
    t = 260.0 + 40.0 * (pmid / pmid[:, -1:]) ** 0.28
    t += rng.uniform(-2.0, 2.0, t.shape)
    # destabilize the top layers with varying strength per column:
    # columns 0..3 stable (no-op), the rest superadiabatic
    for i in range(4, ncol):
        strength = (i - 3) * 2.0
        t[i, :nlvdry + 1] += strength * np.arange(nlvdry + 1, 0, -1)
    q = np.clip(rng.uniform(1e-7, 1e-5, t.shape)
                * (pmid / 1e4) ** 2, 1e-9, 2e-2)

    t_out, q_out = d.drv_dadadj(nlvdry, pmid, pint, pdel, t, q)
    out[f"t_in_{nlvdry}"] = t
    out[f"q_in_{nlvdry}"] = q
    out[f"t_out_{nlvdry}"] = t_out
    out[f"q_out_{nlvdry}"] = q_out

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "dadadj", "source_sha": sha, "ncol": ncol,
        "nlev": nlev, "nlvdry_values": [3, 8],
        "notes": "cam_control_mod default nlvdry=3"}
gold = Path(__file__).resolve().parents[1] / "golden" / "dadadj_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
