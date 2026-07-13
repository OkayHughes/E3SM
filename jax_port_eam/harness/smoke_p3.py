#!/usr/bin/env python3
"""Smoke test for the eam_p3_f extension (run in container).

Initializes P3 from the v4.1.1 lookup table (the only version on local
disk; EAMv3's default 4.1.2 lives in atm/cam/physprops which is not in
the local inputdata — the table is pure input data recorded in the
golden, so version choice cancels between Fortran and JAX) and runs
p3_main on a 4-column warm/mixed/cold/clear set.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_p3_f  # noqa: E402

d = eam_p3_f.p3_driver

TABLE_DIR = "/work/e3sm-inputdata/atm/scream/tables"
d.drv_init(TABLE_DIR.ljust(256), "4.1.1".ljust(16), 1, 0)

mu_r, revap, vn, vm = d.drv_get_tables()
print(f"mu_r table: min={mu_r.min()} max={mu_r.max()}")
print(f"vn[0,0]={vn[0,0]:.6e} vm[0,0]={vm[0,0]:.6e} revap[0,0]={revap[0,0]:.6e}")
print(f"vn[250,0]={vn[250,0]:.6e} vm[250,0]={vm[250,0]:.6e}")

qs = d.drv_qv_sat(np.array([250.0, 273.15, 300.0]),
                  np.array([5.0e4, 8.0e4, 1.0e5]), 0)
print("qv_sat(liq):", qs)

ncol, nlev = 4, 72
F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731

zi = np.linspace(1.0, 0.0, nlev + 1) ** 1.4 * 24000.0
zm = 0.5 * (zi[:-1] + zi[1:])
dz = zi[:-1] - zi[1:]
pm = 1.0e5 * np.exp(-zm / 7600.0)
T0 = np.array([300.0, 285.0, 260.0, 290.0])
T = np.maximum(T0[:, None] - 6.5e-3 * zm[None, :], 200.0)
RAIR, CP, GRAV = 287.042311365049, 1004.64, 9.80616
exner = 1.0 / ((pm * 1e-5) ** (RAIR / CP))
th = T * exner[None, :]
qsatw = 0.622 * 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65)) / pm
qv = 0.8 * qsatw
rho = pm[None, :] / (RAIR * T)
dpres = rho * GRAV * dz[None, :]

qc = np.zeros((ncol, nlev)); qr = np.zeros((ncol, nlev))
qi = np.zeros((ncol, nlev)); qm = np.zeros((ncol, nlev))
nc = np.zeros((ncol, nlev)); nr = np.zeros((ncol, nlev))
ni = np.zeros((ncol, nlev)); bm = np.zeros((ncol, nlev))
cldlay = (zm > 1000) & (zm < 4000)
qc[0, :] = np.where(cldlay, 5e-4, 0.0); nc[0, :] = np.where(cldlay, 5e7, 0.0)
qr[0, :] = np.where(zm < 2000, 1e-4, 0.0); nr[0, :] = np.where(zm < 2000, 1e5, 0.0)
mixed = (zm > 3000) & (zm < 7000)
qc[1, :] = np.where(mixed, 3e-4, 0.0); nc[1, :] = np.where(mixed, 8e7, 0.0)
qi[1, :] = np.where(mixed, 2e-4, 0.0); ni[1, :] = np.where(mixed, 5e4, 0.0)
qm[1, :] = 0.3 * qi[1, :]; bm[1, :] = qm[1, :] / 400.0
cold = (zm > 5000) & (zm < 10000)
qi[2, :] = np.where(cold, 4e-4, 0.0); ni[2, :] = np.where(cold, 1e5, 0.0)
qm[2, :] = 0.5 * qi[2, :]; bm[2, :] = qm[2, :] / 300.0
# column 3: clear

ast = np.clip((qc + qi) * 2000.0, 0.0001, 1.0)
cld_frac_l = np.maximum(ast, 0.0001)
cld_frac_i = np.maximum(ast, 0.0001)
cld_frac_r = np.maximum(ast, 0.0001)
for k in range(1, nlev):
    upfall = (qr[:, k - 1] >= 1e-14) | (qi[:, k - 1] >= 1e-14)
    cld_frac_r[:, k] = np.where(upfall,
                                np.maximum(cld_frac_r[:, k - 1],
                                           cld_frac_r[:, k]),
                                cld_frac_r[:, k])

zero = np.zeros((ncol, nlev))
col_loc = np.ones((ncol, 3))
P = dict(autocon=30500.0, accret=117.25, qcauto=3.19, ncauto=-1.10,
         qcaccr=1.15, wbf=1.0, mincdnc=20.0e6, maxrain=0.005,
         embryo=2.5e-5, nccnst=200.0e6)

state, diag, flux, tend, surf = d.drv_p3_main(
    1800.0, 2, 1, 0, 0,
    P["autocon"], P["accret"], P["qcauto"], P["ncauto"], P["qcaccr"],
    P["wbf"], P["mincdnc"], P["maxrain"], P["embryo"], P["nccnst"],
    F(qc), F(nc), F(qr), F(nr), F(th), F(qv), F(qi), F(qm), F(ni), F(bm),
    F(pm[None, :].repeat(ncol, 0)), F(dz[None, :].repeat(ncol, 0)),
    F(zero + 2e4), F(zero), F(zero + 1e4),
    F(zero + 0.01), F(zero + 0.001), F(zero + 0.001),
    F(zero + 1.0), F(dpres), F(exner[None, :].repeat(ncol, 0)),
    F(cld_frac_r), F(cld_frac_l), F(cld_frac_i), F(qv), F(T), F(col_loc))

state = np.asarray(state)
names = ["qc", "nc", "qr", "nr", "th", "qv", "qi", "qm", "ni", "bm"]
for j, n in enumerate(names):
    print(f"{n:3s}: max={state[:, :, j].max():.4e} min={state[:, :, j].min():.4e}")
print("precip_liq_surf:", np.asarray(surf)[:, 0])
print("precip_ice_surf:", np.asarray(surf)[:, 1])
assert np.all(np.isfinite(state))
clear = 3
assert np.asarray(surf)[clear].max() == 0.0, "clear column must not precip"
print("SMOKE OK")
