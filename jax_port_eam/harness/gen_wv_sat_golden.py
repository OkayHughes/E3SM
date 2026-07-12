#!/usr/bin/env python3
"""Tier-1 golden generator for wv_sat_methods + wv_saturation.

Run in the scream-dev container after build_wv_sat.py:
    docker exec -w /work/E3SM/jax_port_eam/harness scream-dev \
        python3 gen_wv_sat_golden.py

Coverage: full temperature range of the SVP table (127..375 K) for all
four schemes, water/ice/transition; qsat family with derivatives over a
(t, p) grid spanning stratosphere to surface including the es >= p
limiter branch; findsp over sub/super-saturated states in every status
regime (0 converged, 1 unusable, 4 tmin bailout, 8 enthalpy).
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_wv_sat_f  # noqa: E402

d = eam_wv_sat_f.wv_sat_driver
d.drv_init()

out = {}

# --- raw SVP, all schemes, dense temperature sweep ---
t_sweep = np.arange(127.0, 375.01, 0.25)
out["t_sweep"] = t_sweep
for idx, name in [(0, "oldgoffgratch"), (1, "goffgratch"),
                  (2, "murphykoop"), (3, "bolton")]:
    es_w, es_i = d.drv_svp(t_sweep, idx)
    out[f"svp_water_{name}"] = es_w
    out[f"svp_ice_{name}"] = es_i
    out[f"svp_trans_{name}"] = d.drv_svp_trans(t_sweep, idx)

# --- table lookup ---
out["estblf"] = d.drv_estblf(t_sweep)

# --- qsat family with derivatives on a (t, p) product grid ---
t_g = np.arange(150.0, 350.01, 0.5)
p_g = np.array([30.0, 100.0, 1.0e3, 1.0e4, 3.0e4, 5.0e4,
                7.0e4, 8.5e4, 1.0e5, 1.05e5])
tt, pp = np.meshgrid(t_g, p_g, indexing="ij")
tt, pp = tt.ravel(), pp.ravel()
out["qsat_t"], out["qsat_p"] = tt, pp
for fam, fn in [("qsat", d.drv_qsat), ("qsat_water", d.drv_qsat_water),
                ("qsat_ice", d.drv_qsat_ice)]:
    es, qs, gam, dqsdt, enthalpy = fn(tt, pp)
    out[f"{fam}_es"], out[f"{fam}_qs"] = es, qs
    out[f"{fam}_gam"], out[f"{fam}_dqsdt"] = gam, dqsdt
    out[f"{fam}_enthalpy"] = enthalpy

# --- findsp: regimes ---
rng = np.random.default_rng(42)
n = 4000
t_f = rng.uniform(160.0, 340.0, n)
p_f = 10.0 ** rng.uniform(3.0, 5.02, n)
# mix of strongly sub-saturated, near-saturated, and super-saturated q
_, qs0, *_ = d.drv_qsat(t_f, p_f)
fac = rng.choice([0.05, 0.5, 0.95, 1.05, 1.5, 3.0], n)
q_f = np.clip(qs0 * fac, 1e-9, 0.45)
# sprinkle unusable states (status 1: p<=5*es, qs>=0.5) and cold extremes
q_f[:50] = 0.49
t_f[50:100] = 360.0
p_f[100:150] = 20.0
t_f[150:200] = 130.0   # near tmin -> bailout candidates
out["findsp_q"], out["findsp_t"], out["findsp_p"] = q_f, t_f, p_f
for ui, tag in [(1, "ice"), (0, "noice")]:
    tsp, qsp, status = d.drv_findsp(q_f, t_f, p_f, ui)
    out[f"findsp_tsp_{tag}"] = tsp
    out[f"findsp_qsp_{tag}"] = qsp
    out[f"findsp_status_{tag}"] = status

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {
    "scheme": "wv_sat_methods + wv_saturation",
    "source_sha": sha,
    "fortran": "gfortran container, -O2",
    "svp_scheme_indices": {"oldgoffgratch": 0, "goffgratch": 1,
                           "murphykoop": 2, "bolton": 3},
    "notes": "default scheme GoffGratch (idx 1); table 127.16..375.16 K",
}
gold = Path(__file__).resolve().parents[1] / "golden" / "wv_sat_golden.npz"
gold.parent.mkdir(exist_ok=True)
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold} ({len(out)} arrays)")
