#!/usr/bin/env python3
"""Smoke the eam_clubb_f extension (run in the scream-dev container)."""
import numpy as np

from fbuild import import_extension

f = import_extension("eam_clubb_f")
d = f.clubb_driver

nparams = d.drv_get_nparams()
params = d.drv_default_params(nparams)
print("nparams =", nparams)
idx = d.drv_param_indices()
names = ["C1", "C1b", "C1c", "C2rt", "C2thl", "C2rtthl", "C6rt", "C6rtb",
         "C6rtc", "C6thl", "C6thlb", "C6thlc", "C7", "C7b", "C8", "C11",
         "C11b", "C11c", "C14", "beta", "gamma_coef", "gamma_coefb",
         "gamma_coefc", "mu", "nu1", "c_K10", "c_K10h", "wpxp_L_thresh",
         "altitude_threshold", "Skw_denom_coef", "Skw_max_mag"]
for n, i in zip(names, idx):
    print(f"  {n:20s} params[{i}] = {params[i - 1]}")

# EAMv3 overrides (see gen_clubb_golden.py)
ov = {"C1": 2.4, "C1b": 2.8, "C1c": 0.75, "C2rt": 1.75, "C2thl": 1.75,
      "C2rtthl": 1.75 * 1.3, "C6rtb": 7.5, "C6rtc": 0.5, "C6thlb": 7.5,
      "C6thlc": 0.5, "C8": 5.2, "C11": 0.7, "C11b": 0.2, "C11c": 0.85,
      "gamma_coef": 0.12, "gamma_coefb": 0.28, "gamma_coefc": 1.2,
      "mu": 0.0005, "wpxp_L_thresh": 60.0}
for n, v in ov.items():
    params[idx[names.index(n)] - 1] = v

nz = 73
zi = np.zeros(nz)
zi[:] = np.geomspace(25.0, 40000.0, nz)
zi[0] = 0.0
zt = np.empty(nz)
zt[1:] = 0.5 * (zi[1:] + zi[:-1])
zt[0] = -zt[1]

err = d.drv_setup(2, params, zi, zt)
print("setup err =", err)

zm, ztg, dzm, dzt, idzm, idzt, wt2m, wm2t = d.drv_grid_arrays(nz)
assert np.allclose(zm, zi) and np.allclose(ztg, zt)
print("grid ok; dzt[0:3] =", dzt[:3])

azt = 300.0 + 0.01 * zt
azm = np.sin(zm / 5000.0)
o1, o2, o3, o4 = d.drv_grid_ops(azt, azm)
print("zt2zm[-3:] =", o1[-3:])

p = 1.0e5 * np.exp(-zt / 8000.0)
T = np.clip(300.0 - 0.0065 * zt, 180.0, None)
rsl, rsi = d.drv_sat(p, T)
print("rsl[1], rsi[1] =", rsl[1], rsi[1])

# tridiag
rng = np.random.default_rng(0)
n = 12
sub = rng.normal(size=n)
dia = rng.normal(size=n) * 0.5
sup = rng.normal(size=n)
rhs = rng.normal(size=(n, 2))
sol, terr = d.drv_tridag_solve(sup, dia, sub, rhs)
A = np.diag(dia) + np.diag(sub[1:], -1) + np.diag(sup[:-1], 1)
print("tridag err:", terr, "resid:", np.abs(A @ sol - rhs).max())

# pdf_closure on a convective-ish profile
wp2 = np.full(nz, 0.4)
wp3 = np.full(nz, 0.3)
skw = d.drv_skx(wp2, wp3, 2.0e-2)
rtp2 = np.full(nz, 1e-7)
thlp2 = np.full(nz, 0.09)
wprtp = np.full(nz, 5e-5)
wpthlp = np.full(nz, 0.05)
up2 = np.full(nz, 0.3)
vp2 = np.full(nz, 0.3)
upwp = np.full(nz, -0.02)
vpwp = np.full(nz, -0.01)
gamma_skw = 0.28 + (0.12 - 0.28) * np.exp(-0.5 * (skw / 1.2) ** 2)
sig = d.drv_sigma_sqd_w(gamma_skw, wp2, thlp2, rtp2, up2, vp2,
                        wpthlp, wprtp, upwp, vpwp)
exner = (p / 1e5) ** 0.2856
thlm = T / exner
rtm = np.minimum(0.8 * rsl, 0.016)
skthl = d.drv_skx(thlp2, np.full(nz, 0.001), 1e-2)
skrt = d.drv_skx(rtp2, np.full(nz, 1e-11), 1e-8)
mom, pdfp, sig_out, perr = d.drv_pdf_closure(
    p, exner, T / thlm, np.zeros(nz), wp2, wp3, sig, skw, skthl, skrt,
    rtm, rtp2, wprtp, thlm, thlp2, wpthlp,
    np.full(nz, 5.0), up2, upwp, np.full(nz, -3.0), vp2, vpwp,
    np.full(nz, 1e-5))
print("pdf_closure err:", perr)
print("cloud_frac range:", mom[:, 2].min(), mom[:, 2].max())
print("rcm max:", mom[:, 4].max(), "mixt_frac range:",
      pdfp[:, 44].min(), pdfp[:, 44].max())
print("SMOKE OK")
