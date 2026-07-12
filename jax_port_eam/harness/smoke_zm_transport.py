#!/usr/bin/env python3
"""Container smoke test for the zm_transport f2py extension: float64
wrapper check, then physically consistent mass fluxes from re-running
eam_zm_conv_f on the zm_conv golden profiles, and sanity properties
(uniform tracer/wind no-op, untriggered columns zero, column-mass
conservation of the tracer tendency)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "fmod"))
import eam_zm_conv_f  # noqa: E402
import eam_zm_transport_f  # noqa: E402

dtt = eam_zm_transport_f.zm_transport_driver.drv_transport_tracer
dtm = eam_zm_transport_f.zm_transport_driver.drv_transport_momentum
assert "array('d')" in dtt.__doc__ and "array('d')" in dtm.__doc__
print("zm_transport wrappers float64 OK")

g = np.load(ROOT / "golden" / "zm_conv_golden.npz")
meta = json.loads(str(g["__metadata__"]))
p = meta["params"]["a"]
F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731

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
(lengath, ideep, mx_g, jctop, jcbot, jt_g, prec, heat, qtnd, cape,
 dcape, mcon, pflx, zdu, mu, eu, du, md, ed, dp, dsubcld, ql, rliq,
 rprd, dlf) = res
lengath = int(lengath)
ncol, nlev = np.asarray(g["t"]).shape
print(f"zm_conv rerun: lengath={lengath}")
assert lengath > 0 and lengath < ncol

# uniform tracer + fracis=1 must be a no-op (mass continuity of the
# plume fluxes); uniform winds must give zero momentum tendency
ncnst = 3
q = np.ones((ncol, nlev, ncnst), order="F")
q[:, :, 2] = 1e-6
fracis = np.ones((ncol, nlev, ncnst), order="F")
dpdry = F(np.zeros((ncol, nlev)))
dqdt = dtt(0, [0, 1, 1], [0, 0, 0], q, F(mu), F(md), F(du), F(eu),
           F(ed), F(dp), jt_g, mx_g, ideep, lengath, fracis, dpdry,
           1800.0)
print("uniform tracer max |dqdt| =", np.abs(dqdt).max())
assert np.abs(dqdt).max() < 1e-18

winds = np.zeros((ncol, nlev, 2), order="F")
winds[:, :, 0] = 7.0
winds[:, :, 1] = -3.0
wt, pgu, pgd, icwu, icwd, seten = dtm(lengath, winds, F(mu), F(md),
                                      F(du), F(eu), F(ed), F(dp),
                                      jt_g, mx_g, ideep, 1800.0)
print("uniform wind max |tend| =", np.abs(wt).max(),
      "max |seten| =", np.abs(seten).max())
assert np.abs(wt).max() < 1e-16 and np.abs(seten).max() < 1e-15

# structured tracer: column-mass conservation + untriggered zeros
rng = np.random.default_rng(7)
sig = np.linspace(0.01, 1.0, nlev)
q[:, :, 1] = 1e-6 * np.exp(3.0 * (sig - 1.0))[None, :] \
    * (1 + 0.2 * rng.uniform(-1, 1, (ncol, nlev)))
q[:, :, 2] = 1e-7 + 1e-6 * np.exp(-0.5 * ((sig - 0.5) / 0.15) ** 2)[None, :]
q = np.asarray(q, order="F")
dqdt = dtt(0, [0, 1, 1], [0, 0, 0], q, F(mu), F(md), F(du), F(eu),
           F(ed), F(dp), jt_g, mx_g, ideep, lengath, fracis, dpdry,
           1800.0)
gset = ideep[:lengath] - 1
idle = sorted(set(range(ncol)) - set(gset.tolist()))
assert np.all(dqdt[idle] == 0.0)
# gathered-column integral of dqdt*dp (dp gathered rows)
for m in (1, 2):
    col = np.abs(sum(dqdt[gset, k, m] * dp[:lengath, k]
                     for k in range(nlev)))
    scale = np.abs(dqdt[gset, :, m] * dp[:lengath]).max()
    print(f"tracer m={m}: max |sum dqdt*dp| = {col.max():.3e} "
          f"(scale {scale:.3e})")
    assert col.max() < 1e-12 * max(scale, 1e-30) or col.max() < 1e-22

# momentum: sheared wind, conservation + untriggered zeros
winds[:, :, 0] = 20.0 * (1 - sig)[None, :] + rng.uniform(-2, 2, (ncol, nlev))
winds[:, :, 1] = 5.0 * np.sin(3 * sig)[None, :]
winds = np.asarray(winds, order="F")
wt, pgu, pgd, icwu, icwd, seten = dtm(lengath, winds, F(mu), F(md),
                                      F(du), F(eu), F(ed), F(dp),
                                      jt_g, mx_g, ideep, 1800.0)
assert np.all(wt[idle] == 0.0) and np.all(seten[idle] == 0.0)
assert np.array_equal(icwu[idle], winds[idle])
for m in (0, 1):
    col = np.abs(sum(wt[gset, k, m] * dp[:lengath, k]
                     for k in range(nlev)))
    print(f"momentum m={m}: max |sum tend*dp| = {col.max():.3e}")
    assert col.max() < 1e-13
print("zm_transport smoke OK")
