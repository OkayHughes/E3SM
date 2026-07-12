#!/usr/bin/env python3
"""Smoke test for the tropopause f2py extension (run in container)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_tropopause_f  # noqa: E402

d = eam_tropopause_f.tropopause_driver

# f2py kind check: wrappers must be float64
assert "array('d')" in (d.drv_find.__doc__ or ""), d.drv_find.__doc__
print("drv_find wrapper is float64: OK")

PCOLS, NLEV = 16, 72
RAIR, GRAV = 287.042311365, 9.80616

# ---- grid: hybrid-ish pint from ~225 Pa to 1e5 Pa -------------------
ai = np.linspace(0.0, 1.0, NLEV + 1) ** 1.7
pint1 = 225.5 + ai * (1.0e5 - 225.5)
pint = np.broadcast_to(pint1, (PCOLS, NLEV + 1)).copy()
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])


def build_column(ts, ztrop, gam_strat=0.0):
    """T(z): lapse -6.5 K/km to ztrop, then gam_strat above; z built
    hydrostatically (fixed point, 8 iters)."""
    p, pi = pmid[0], pint[0]
    T = np.full(NLEV, 250.0)
    for _ in range(8):
        zi = np.zeros(NLEV + 1)
        for k in range(NLEV - 1, -1, -1):
            zi[k] = zi[k + 1] + RAIR * T[k] / GRAV * np.log(pi[k + 1] / pi[k])
        zm = zi[1:] + RAIR * T / GRAV * np.log(pi[1:] / p)
        T = np.where(zm <= ztrop, ts - 6.5e-3 * zm,
                     ts - 6.5e-3 * ztrop + gam_strat * (zm - ztrop))
        T = np.maximum(T, 150.0)
    return T, zm, zi


t = np.zeros((PCOLS, NLEV))
zm = np.zeros((PCOLS, NLEV))
zi = np.zeros((PCOLS, NLEV + 1))
ztrops = np.array([16.5e3, 16.0e3, 11.0e3, 10.0e3, 8.5e3, 8.0e3] +
                  [12.0e3] * 10)
tss = np.array([300.0, 299.0, 288.0, 285.0, 258.0, 255.0] + [285.0] * 10)
for i in range(PCOLS):
    t[i], zm[i], zi[i] = build_column(tss[i], ztrops[i], gam_strat=1.0e-3)

# synthetic climatology: constant in lon, per-lat monthly values
nlon, nlat = 4, PCOLS
lon_deg = np.array([0.0, 90.0, 180.0, 270.0])
lat_deg = np.linspace(-75.0, 75.0, nlat)
months = np.arange(12)
clim = (30000.0 - 20000.0 * np.cos(np.deg2rad(lat_deg))[None, :, None] ** 2
        + 3000.0 * np.sin(np.deg2rad(lat_deg))[None, :, None]
        * np.cos(2 * np.pi * months / 12.0)[None, None, :])
tropp3 = np.broadcast_to(clim, (nlon, nlat, 12)).copy()
col_lat = np.deg2rad(lat_deg)
col_lon = np.deg2rad(lon_deg[np.arange(PCOLS) % nlon])

d.drv_init(np.asfortranarray(lon_deg), np.asfortranarray(lat_deg),
           np.asfortranarray(tropp3), col_lat, col_lon)

TWMO, CLIMATE, HYBSTOB, NONE = 5, 3, 7, 1
F = np.asfortranarray
lev, tp, tt, tz = d.drv_find(PCOLS, 100.0, TWMO, NONE, col_lat,
                             F(t), F(pmid), F(pint), F(zm), F(zi))
print("TWMO-only:")
print("  tropLev:", lev)
print("  tropP (hPa):", np.round(tp / 100.0, 1))
print("  tropT (K):", np.round(tt, 1))
print("  tropZ (km):", np.round(tz / 1000.0, 1))

assert (lev != -1).all(), "twmo should find all these"
# tropical ~100 hPa / 16-17 km; polar lower (higher pressure)
assert 80e2 <= tp[0] <= 130e2, tp[0]
assert 15.0e3 <= tz[0] <= 18.0e3, tz[0]
assert tp[4] > tp[0] and tz[4] < tz[0], (tp[4], tp[0])
assert 200e2 <= tp[4] <= 400e2, tp[4]

# climate-only: reproduces node climatology exactly (time-interp'd)
days = np.array([16., 45., 75., 105., 136., 166., 197., 228., 258.,
                 289., 319., 350.])
calday = 100.0
m = np.searchsorted(days, calday, side="right")  # days[m-1] <= calday
dels = (calday - days[m - 1]) / (days[m] - days[m - 1])
expect = clim[0, :, m - 1] + dels * (clim[0, :, m] - clim[0, :, m - 1])
lev2, tp2, tt2, tz2 = d.drv_find(PCOLS, calday, CLIMATE, NONE, col_lat,
                                 F(t), F(pmid), F(pint), F(zm), F(zi))
print("CLIMATE-only tropP (hPa):", np.round(tp2 / 100.0, 1))
err = np.abs(tp2 / expect - 1.0).max()
print("max rel |fortran - expected clim interp|:", err)
# node reproduction is exact (dels=0 caldays replay bitwise); for
# interpolated caldays gfortran -O2/aarch64 contracts a+dels*(b-a)
# into an fma -> <=1 ulp
assert err < 1e-15, "climatology injection not near-exact"

# hybridstobie: always finds
lev3, tp3, tt3, tz3 = d.drv_find(PCOLS, calday, HYBSTOB, NONE, col_lat,
                                 F(t), F(pmid), F(pint), F(zm), F(zi))
print("HYBSTOB tropLev:", lev3)
assert (lev3 != -1).all()

# twmo-failure column: constant steep lapse everywhere -> backup engages
t_fail = t.copy()
for i in range(PCOLS):
    tf = 310.0 - 6.5e-3 * zm[i]
    t_fail[i] = np.maximum(tf, 130.0)
lev4, tp4, _, _ = d.drv_find(PCOLS, calday, TWMO, NONE, col_lat,
                             F(t_fail), F(pmid), F(pint), F(zm), F(zi))
print("TWMO on constant-lapse columns, tropLev:", lev4)
assert (lev4 == -1).all(), "expected twmo failure"
lev5, tp5, _, _ = d.drv_find(PCOLS, calday, TWMO, CLIMATE, col_lat,
                             F(t_fail), F(pmid), F(pint), F(zm), F(zi))
assert (lev5 != -1).all(), "backup should engage"
assert np.abs(tp5 / expect - 1.0).max() < 1e-15
print("backup engages with (near-)exact climatology: OK")
print("SMOKE OK")
