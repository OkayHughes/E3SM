#!/usr/bin/env python3
"""Tier-1 golden generator for conv_water_4rad (run in the scream-dev
container). 16 columns x 72 levels of stratified-random per-point
states covering every branch: shallow/deep fractions and in-cloud
water at 0 / below the frac_limit-ic_limit thresholds / significant,
stratiform fraction ast across the frac_limit threshold, condensate
ice fractions 0/interior/1, temperatures across the pergro 243/263 K
nodes, rei across the 13/130 micron clamps, and NaN patches in FICE
(guarded COSP outputs only).

Five configurations sweep the runtime switches:
  a: mode=1, zm_microp=T, P3, pergro=F  -- EAMv3 production default
     (bld/namelist_files/namelist_defaults_eam.xml phys="default":
     conv_water_in_rad=1, zmconv_microp=.true., microp_scheme P3,
     pergro_mods=.false.). mode 0 means cloud_diagnostics never calls
     conv_water_4rad, so only modes 1 and 2 exist inside the routine.
  b: mode=1, zm_microp=F, P3, pergro=F
  c: mode=2, zm_microp=F, P3, pergro=F  (emissivity-weighted average)
  d: mode=1, zm_microp=F, P3, pergro=T  (pergro repartition of wrk1)
  e: mode=2, zm_microp=F, RK, pergro=F  (legacy unclamped kabsi)
The zm_microp branch never reads conv_water_mode, so a covers it for
both modes (asserted below).
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_conv_water_f  # noqa: E402

dcw = eam_conv_water_f.conv_water_driver.drv_conv_water
assert "array('d')" in dcw.__doc__

rng = np.random.default_rng(20260712)
ncol, nlev = 16, 72
npts = (ncol, nlev)

# stratified draws: explicit threshold-straddling values mixed with
# uniform ranges (frac_limit=0.01, ic_limit=1e-12 in conv_water.F90)
frac_choices = np.array([0.0, 0.004, 0.01, 0.011, 0.05, 0.15, 0.35])
icw_choices = np.array([0.0, 5.0e-13, 1.0e-12, 2.0e-12, 5.0e-5,
                        4.0e-4, 2.0e-3])
sh_frac = rng.choice(frac_choices, npts)
dp_frac = rng.choice(frac_choices, npts)
sh_icwmr = rng.choice(icw_choices, npts)
dp_icwmr = rng.choice(icw_choices, npts)
dp_icimr = dp_icwmr * rng.uniform(0.0, 1.0, npts)
ast = rng.choice(np.array([0.0, 0.005, 0.02, 0.2, 0.6, 0.95]), npts)

# stratiform condensate: zeros, pergro-tiny (< 100*ic_limit after the
# /max(0.01, ast) in-cloud conversion), and significant, with pure-ice
# and pure-liquid points for the wrk1 extremes
mag = rng.choice(np.array([0.0, 5.0e-13, 3.0e-8, 2.0e-5, 3.0e-4]),
                 npts)
wfrac = rng.choice(np.array([0.0, 0.13, 0.5, 0.87, 1.0]), npts)
q_cldice = mag * wfrac
q_cldliq = mag * (1.0 - wfrac)

t = rng.uniform(200.0, 305.0, npts)
t.ravel()[:8] = [242.9, 243.0, 243.1, 253.0, 262.9, 263.0, 263.1, 280.0]
pdel = rng.uniform(300.0, 1500.0, npts)
rei = rng.choice(np.array([5.0, 12.9, 13.0, 14.0, 60.0, 129.0, 130.0,
                           131.0, 200.0]), npts)
fice = rng.uniform(0.0, 1.0, npts)
fice[rng.uniform(size=npts) < 0.15] = np.nan

F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731
INP = dict(t=t, pdel=pdel, q_cldliq=q_cldliq, q_cldice=q_cldice,
           sh_icwmr=sh_icwmr, dp_icwmr=dp_icwmr, dp_icimr=dp_icimr,
           fice=fice, sh_frac=sh_frac, dp_frac=dp_frac, ast=ast,
           rei=rei)
ARGS = [F(INP[k]) for k in ("t", "pdel", "q_cldliq", "q_cldice",
                            "sh_icwmr", "dp_icwmr", "dp_icimr",
                            "fice", "sh_frac", "dp_frac", "ast",
                            "rei")]

CONFIGS = {  # mode, zm_microp, is_rk, pergro
    "a": (1, 1, 0, 0),
    "b": (1, 0, 0, 0),
    "c": (2, 0, 0, 0),
    "d": (1, 0, 0, 1),
    "e": (2, 0, 1, 0),
}
NAMES = ("totg_liq", "totg_ice", "sh_cldliq", "sh_cldice")

out = dict(INP)
for cfg, (mode, microp, is_rk, pergro) in CONFIGS.items():
    res = dcw(mode, microp, is_rk, pergro, *ARGS)
    r = dict(zip(NAMES, (np.asarray(v) for v in res)))
    out.update({f"{n}_{cfg}": v for n, v in r.items()})
    for n in NAMES:  # everything finite and nonnegative
        assert np.isfinite(r[n]).all(), (cfg, n)
        assert r[n].min() >= 0.0, (cfg, n)
    print(f"cfg {cfg}: totg_liq mean={r['totg_liq'].mean():.3e} "
          f"nonzero={(r['totg_liq'] > 0).mean():.2f} "
          f"totg_ice mean={r['totg_ice'].mean():.3e}")

# the zm_microp branch ignores conv_water_mode
res2 = dcw(2, 1, 0, 0, *ARGS)
for n, v in zip(NAMES, res2):
    assert np.array_equal(np.asarray(v), out[f"{n}_a"]), n
# pergro must change the partition somewhere (tiny ls_icwmr points)
assert not np.array_equal(out["totg_ice_b"], out["totg_ice_d"])
# RK vs P3 kabsi must matter somewhere in mode 2
assert not np.array_equal(out["totg_liq_c"], out["totg_liq_e"])
# NaN fice zeroes only the COSP outputs
nanm = np.isnan(fice)
assert np.all(out["sh_cldliq_a"][nanm] == 0.0)
assert np.all(out["sh_cldice_a"][nanm] == 0.0)
good = ~nanm
np.testing.assert_allclose(
    out["sh_cldliq_a"][good] + out["sh_cldice_a"][good],
    (sh_icwmr * sh_frac)[good], rtol=1e-12, atol=1e-30)

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "conv_water conv_water_4rad (grid-box condensate "
                  "for radiation)",
        "source_sha": sha, "ncol": ncol, "nlev": nlev,
        "params": {c: {"conv_water_mode": v[0], "zm_microp": v[1],
                       "microp_scheme": "RK" if v[2] else "P3",
                       "pergro_mods": bool(v[3])}
                   for c, v in CONFIGS.items()},
        "defaults": {"conv_water_in_rad": 1, "zmconv_microp": True,
                     "microp_scheme": "P3", "pergro_mods": False,
                     "note": "EAMv3 phys='default'; mode 0 = routine "
                             "not called (cloud_diagnostics.F90)"},
        "constants": {"kabsl": 0.090361, "frac_limit": 0.01,
                      "ic_limit": 1e-12, "gravit": 9.80616}}
gold = Path(__file__).resolve().parents[1] / "golden" \
    / "conv_water_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
