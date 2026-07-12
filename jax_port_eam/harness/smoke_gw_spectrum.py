#!/usr/bin/env python3
"""Container smoke test for the gw spectrum extension (float64 wrapper
check + one pass through every driver entry point)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_gw_spectrum_f as f  # noqa: E402

d = f.gw_spectrum_driver
for name in dir(d):
    if name.startswith("drv"):
        doc = getattr(d, name).__doc__
        assert "array('d')" in doc, (name, doc)
print("all drv_* wrappers are float64")

ncol, nlev, pgwv = 4, 72, 32
nwav = 2 * pgwv + 1
ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
pint = 225.5 + ai * (1e5 - 225.5)
pint = np.broadcast_to(pint, (ncol, nlev + 1)).copy()
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
dpm = np.diff(pint, axis=1)
rdpm = 1 / dpm
piln = np.log(pint)
t = 216 + 72 * (pmid / pmid[:, -1:]) ** 0.32
zm = (287.042 * t / 9.80616) * np.log(pint[:, -1:] / pmid)
u = 20 * np.exp(-((zm - 1.1e4) / 8e3) ** 2) + 2
v = 0.3 * u
alpha = np.full(nlev + 1, 1e-6)
pref_edge = pint[0]
kbotbg = int(np.sum(pref_edge < 5e4)) - 1
kfront = int(np.sum(pref_edge < 6e4))
maxh, maxuh = 20, 40
h = np.arange(1, maxh + 1)[:, None, None]
uh = np.arange(-maxuh, maxuh + 1)[None, :, None]
cl = np.arange(-pgwv, pgwv + 1)[None, None, :] * 2.5
mfcc = 4e6 * (h / maxh) ** 2 * np.exp(-((cl - 0.4 * uh) / 25.0) ** 2) \
    * (0.05 + (cl - 0.4 * uh) ** 2 / 625.0)
d.drv_gw_spec_init(2.5, 0, kbotbg, 1.0, 6.28e-5, 9.80616, 287.042,
                   alpha, 0, 2.5e-3, 1.25e-15, kfront, 70000.0,
                   np.asfortranarray(mfcc), pref_edge)
rhoi, ti, nm, ni = d.drv_gw_prof(1004.64, t, pmid, pint)
netdt = np.zeros((ncol, nlev))
netdt[:, 40:60] = 6e-4
lat = np.deg2rad(np.array([0., 15., -30., 45.]))
res = d.drv_gw_beres_src(nwav, pgwv, lat, u, v, netdt, zm,
                         10.0, 0.5, 2.5, 10.0, 1)
src, tend, tau, ubm, ubi, xv, yv, c, hdepth, maxq0 = res
print("beres: src", src, "hdepth", hdepth, "tau max", tau.max(),
      "maxq0", maxq0)
frontgf = np.full((ncol, nlev), 2e-15)
res2 = d.drv_gw_cm_src(nwav, pgwv, kbotbg, u, v, frontgf)
src2, tend2, tau2, ubm2, ubi2, xv2, yv2, c2 = res2
print("cm: src", src2, "tau max", tau2.max())
q = np.stack([1e-3 * np.exp(-zm / 2e3), 1e-6 + 1e-7 * np.sin(zm / 5e3)],
             axis=-1)
dse = 1004.64 * t + 9.80616 * zm
out = d.drv_gw_drag_prof_spec(pgwv, src, tend, 0, 1800.0, lat,
                              t, ti, pmid, pint, dpm, rdpm, piln, rhoi,
                              nm, ni, ubm, ubi, xv, yv, 0.35, c, q, dse,
                              tau)
tau_o, utgw, vtgw, ttgw, qtgw, taucd, egwdffi, gwut, dttdf, dttke = out
print("drag: utgw", np.abs(utgw).max(), "qtgw", np.abs(qtgw).max(),
      "taucd", np.abs(taucd).max(), "egwdffi", egwdffi.max())
mec = d.drv_momentum_energy_conservation(tend, 1800.0, taucd, pint, dpm,
                                         u, v, utgw, vtgw, ttgw,
                                         utgw, vtgw, ttgw)
print("mec dudt", np.abs(mec[0]).max(),
      "finite", all(np.isfinite(a).all() for a in mec))
ksrf = np.zeros(ncol)
kv = np.maximum(egwdffi, 0.0) + 1e-3
tmpi = np.abs(np.random.default_rng(1).normal(1.0, 0.1, (ncol, nlev + 1)))
cc_top = np.zeros(ncol)
qo, ca, cc2, dnom, ze = d.drv_vd_lu(ksrf, kv, tmpi, rdpm, 1800.0,
                                    9.80616, cc_top, 1, kbotbg + 1, q,
                                    np.zeros(ncol))
print("vd_lu dq max", np.abs(qo - q).max(), "dnom", dnom.max())
