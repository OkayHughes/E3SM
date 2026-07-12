#!/usr/bin/env python3
"""Tier-1 golden generator for the tropopause finder (run in the
scream-dev container).

40 columns x 72 levels, hydrostatically consistent synthetic
profiles:
  0-7   tropical (Ts 297-303 K, tropopause 15.5-17.2 km, stratosphere
        isothermal or inversion up to +2.5 K/km)
  8-15  midlatitude (Ts 282-290 K, tropopause 10-12.5 km)
  16-23 polar (Ts 248-266 K, tropopause 7.5-9.5 km, half with a
        surface inversion)
  24-27 constant -6.5 K/km lapse to the 130 K floor -> twmo FAILS,
        the backup engages
  28-29 fully isothermal columns (twmo finds ~plimu)
  30-33 midlatitude with a mid-troposphere inversion layer
  34-39 randomly perturbed midlatitude profiles

Climatology: 16 lat nodes (-80..80 deg), 4 lons, lon-constant
trop_p(lat, month) with a seasonal cycle + jitter, served through the
pio stub; column lats sit exactly on the lat nodes so the REAL
tropopause_read_file/lininterp path reproduces node values exactly
(verified: dels=0 caldays replay bitwise). All chunks share lchnk=1,
so global column i uses climatology slot i%16.

Runs recorded (TROP_ALG enums: NONE=1, CLIMATE=3, TWMO=5, HYBSTOB=7):
  twmo             TWMO only (backup NONE)
  twmo_climate_cd* TWMO + CLIMATE (production default) per calday
  climate_cd*      CLIMATE only per calday (time-branch coverage:
                   wrap-lo, ==days(1), mid, dels=0 mid, ==days(12),
                   wrap-hi)
  hybstob          HYBSTOB only
  hybstob_climate  HYBSTOB + CLIMATE (production chemistry chain;
                   identical to hybstob since HYBSTOB cannot fail)

tropLev is stored RAW from the Fortran: 1-based, NOTFOUND=-1
(the JAX port returns 0-based levels; tests add 1).
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_tropopause_f  # noqa: E402

d = eam_tropopause_f.tropopause_driver
F = np.asfortranarray

PCOLS, NLEV, NCOL = 16, 72, 40
RAIR, GRAV = 287.042311365, 9.80616
NONE, CLIMATE, TWMO, HYBSTOB = 1, 3, 5, 7
DAYS = np.array([16., 45., 75., 105., 136., 166., 197., 228., 258.,
                 289., 319., 350.])

rng = np.random.default_rng(20260712)

# ---------------------------------------------------------------- #
# grid + profiles                                                   #
# ---------------------------------------------------------------- #
ai = np.linspace(0.0, 1.0, NLEV + 1) ** 1.7
pint1 = 225.5 + ai * (1.0e5 - 225.5)
pint = np.broadcast_to(pint1, (NCOL, NLEV + 1)).copy()
pint *= (1 + 0.03 * rng.uniform(-1, 1, (NCOL, 1)))
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])


def build_column(i, tfun, tfloor=150.0):
    """T = tfun(zm); z built hydrostatically (fixed point, 10 iters)."""
    p, pi = pmid[i], pint[i]
    T = np.full(NLEV, 250.0)
    zi = np.zeros(NLEV + 1)
    zm = np.zeros(NLEV)
    for _ in range(10):
        for k in range(NLEV - 1, -1, -1):
            zi[k] = zi[k + 1] + RAIR * T[k] / GRAV * np.log(pi[k + 1] / pi[k])
        zm = zi[1:] + RAIR * T / GRAV * np.log(pi[1:] / p)
        T = np.maximum(tfun(zm), tfloor)
    return T, zm.copy(), zi.copy()


def std_profile(ts, ztrop, gam_strat, inv_amp=0.0, inv_depth=1500.0,
                blob=None):
    def tfun(z):
        T = np.where(z <= ztrop, ts - 6.5e-3 * z,
                     ts - 6.5e-3 * ztrop + gam_strat * (z - ztrop))
        if inv_amp:
            T = T - inv_amp * np.maximum(0.0, 1.0 - z / inv_depth)
        if blob is not None:
            zc, amp = blob
            T = T + amp * np.exp(-((z - zc) / 800.0) ** 2)
        return T
    return tfun


t = np.zeros((NCOL, NLEV))
zm = np.zeros((NCOL, NLEV))
zi = np.zeros((NCOL, NLEV + 1))
ptype = np.empty(NCOL, dtype="U16")

strat_cycle = [0.0, 1.0e-3, 2.0e-3, 2.5e-3]
for j, i in enumerate(range(0, 8)):
    ptype[i] = "tropical"
    t[i], zm[i], zi[i] = build_column(
        i, std_profile(297.0 + 0.75 * j, 15.5e3 + 212.5 * j,
                       strat_cycle[j % 4]))
for j, i in enumerate(range(8, 16)):
    ptype[i] = "midlat"
    t[i], zm[i], zi[i] = build_column(
        i, std_profile(282.0 + j, 10.0e3 + 312.5 * j,
                       [0.0, 5.0e-4, 1.0e-3][j % 3]))
for j, i in enumerate(range(16, 24)):
    ptype[i] = "polar"
    t[i], zm[i], zi[i] = build_column(
        i, std_profile(248.0 + 2.25 * j, 7.5e3 + 250.0 * j, 1.0e-3,
                       inv_amp=(8.0 if j % 2 else 0.0)))
for j, i in enumerate(range(24, 28)):
    ptype[i] = "twmo_fail"
    ts = 305.0 + 3.0 * j
    t[i], zm[i], zi[i] = build_column(
        i, lambda z, ts=ts: ts - 6.5e-3 * z, tfloor=130.0)
for j, i in enumerate(range(28, 30)):
    ptype[i] = "isothermal"
    tc = [200.0, 230.0][j]
    t[i], zm[i], zi[i] = build_column(i, lambda z, tc=tc: tc + 0.0 * z)
for j, i in enumerate(range(30, 34)):
    ptype[i] = "inversion"
    t[i], zm[i], zi[i] = build_column(
        i, std_profile(286.0, 11.0e3, 5.0e-4,
                       blob=(3.0e3 + 600.0 * j, 4.0 + 1.2 * j)))
for j, i in enumerate(range(34, 40)):
    ptype[i] = "random"
    base = std_profile(283.0 + 1.3 * j, 10.5e3 + 300.0 * j, 8.0e-4)
    pert = np.convolve(rng.normal(0.0, 1.5, NLEV + 8),
                       np.ones(9) / 9.0, mode="valid")

    def tfun(z, base=base, pert=pert):
        return base(z) + pert
    t[i], zm[i], zi[i] = build_column(i, tfun)

# ---------------------------------------------------------------- #
# climatology                                                       #
# ---------------------------------------------------------------- #
nlon, nlat = 4, PCOLS
lon_deg = np.array([0.0, 90.0, 180.0, 270.0])
lat_deg = np.linspace(-80.0, 80.0, nlat)
latr = np.deg2rad(lat_deg)
season = np.cos(2 * np.pi * (DAYS - 16.0) / 365.0)
clim = (28000.0 - 19000.0 * np.cos(latr)[:, None] ** 2
        + 2500.0 * np.sin(latr)[:, None] * season[None, :]
        + rng.uniform(-300.0, 300.0, (nlat, 12)))          # (nlat, 12)
tropp3 = np.broadcast_to(clim[None, :, :], (nlon, nlat, 12)).copy()
col_lat16 = np.deg2rad(lat_deg)                    # on the nodes
col_lon16 = np.deg2rad(lon_deg[np.arange(PCOLS) % nlon])

d.drv_init(F(lon_deg), F(lat_deg), F(tropp3), col_lat16, col_lon16)

slots = np.arange(NCOL) % PCOLS
tropp_clim_cols = clim[slots]                      # (NCOL, 12)
col_lat = col_lat16[slots]

# calday cases: wrap-lo, ==days(1) (dels=0), mid, mid, mid dels=0
# (==days(6)), late mid, ==days(12) (wrap dels=0), wrap-hi
CALDAYS = np.array([1.5, 16.0, 45.75, 100.0, 166.0, 349.0, 350.0,
                    360.25])


def find(primary, backup, calday):
    lev = np.zeros(NCOL, dtype=np.int64)
    tp, tt, tz = (np.zeros(NCOL) for _ in range(3))
    for beg in range(0, NCOL, PCOLS):
        m = min(PCOLS, NCOL - beg)
        sl = slice(beg, beg + m)

        def pad2(a):
            buf = np.repeat(a[beg:beg + 1], PCOLS, axis=0)
            buf[:m] = a[sl]
            return F(buf)
        latb = np.repeat(col_lat[beg:beg + 1], PCOLS)
        latb[:m] = col_lat[sl]
        lb, pb, tb, zb = d.drv_find(m, calday, primary, backup, latb,
                                    pad2(t), pad2(pmid), pad2(pint),
                                    pad2(zm), pad2(zi))
        lev[sl], tp[sl], tt[sl], tz[sl] = lb[:m], pb[:m], tb[:m], zb[:m]
    return lev, tp, tt, tz


out = dict(t=t, pmid=pmid, pint=pint, zm=zm, zi=zi, col_lat=col_lat,
           tropp_clim=tropp_clim_cols, days=DAYS, caldays=CALDAYS,
           ptype=ptype, clim_nodes=clim, clim_lat_deg=lat_deg)

runs = {"twmo": (TWMO, NONE, 100.0),
        "hybstob": (HYBSTOB, NONE, 100.0),
        "hybstob_climate": (HYBSTOB, CLIMATE, 100.0)}
for k, cd in enumerate(CALDAYS):
    runs[f"twmo_climate_cd{k}"] = (TWMO, CLIMATE, cd)
    runs[f"climate_cd{k}"] = (CLIMATE, NONE, cd)

for name, (pa, ba, cd) in runs.items():
    lev, tp, tt, tz = find(pa, ba, cd)
    out[f"{name}_lev"] = lev
    out[f"{name}_p"] = tp
    out[f"{name}_t"] = tt
    out[f"{name}_z"] = tz

# ---------------------------------------------------------------- #
# sanity                                                            #
# ---------------------------------------------------------------- #
lev, tp, tt, tz = (out["twmo_lev"], out["twmo_p"], out["twmo_t"],
                   out["twmo_z"])
fail = ptype == "twmo_fail"
assert (lev[fail] == -1).all(), lev[fail]
assert (lev[~fail] != -1).all(), lev[~fail]
trop, pol = ptype == "tropical", ptype == "polar"
assert 80e2 <= tp[trop].max() <= 130e2, tp[trop]
assert (tz[trop] > 15.0e3).all() and (tz[trop] < 18.0e3).all()
assert (tp[pol] > tp[trop].max()).all()
assert (tz[pol] < 11.0e3).all()
assert (out["twmo_climate_cd3_lev"] != -1).all()
# backup engages exactly on the twmo failures
assert (out["twmo_climate_cd3_p"][~fail] == tp[~fail]).all()
assert (out["twmo_climate_cd3_p"][fail] != tp[fail]).all()
assert (out["hybstob_lev"] != -1).all()
for suff in ("lev", "p", "t", "z"):
    assert (out[f"hybstob_{suff}"]
            == out[f"hybstob_climate_{suff}"]).all()
# climatology injection is exact at dels=0 caldays
np.testing.assert_array_equal(out["climate_cd1_p"],
                              tropp_clim_cols[:, 0])
np.testing.assert_array_equal(out["climate_cd4_p"],
                              tropp_clim_cols[:, 5])
np.testing.assert_array_equal(out["climate_cd6_p"],
                              tropp_clim_cols[:, 11])
print("twmo tropP (hPa):", np.round(tp / 100.0, 1))
print("twmo tropZ (km):", np.round(tz / 1000.0, 2))
print("found fraction (twmo):", (lev != -1).mean())

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": ("tropopause (twmo/climate/hybridstobie + "
                   "tropopause_find dispatch)"),
        "source_sha": sha, "ncol": NCOL, "nlev": NLEV,
        "params": {"gam": -0.002, "plimu": 45000.0, "pliml": 7500.0,
                   "deltaz": 2000.0, "alpha": 0.03,
                   "min_stobie_pressure": 40.0e2,
                   "max_linoz_pressure": 208.0e2,
                   "fillvalue": 1.0e20, "notfound": -1},
        "trop_alg": {"NONE": 1, "CLIMATE": 3, "TWMO": 5, "HYBSTOB": 7},
        "notes": ("tropLev stored 1-based (Fortran), NOTFOUND=-1; the "
                  "JAX port is 0-based (+1 to compare). Climatology "
                  "served through the pio stub with column lats on "
                  "the lat nodes -> tropopause_read_file reproduces "
                  "clim_nodes exactly (dels=0 caldays replay "
                  "bitwise). gfortran -O2/aarch64 contracts "
                  "a+dels*(b-a) to fma: interpolated-calday "
                  "climatology pressures can differ from unfused "
                  "arithmetic by ~1 ulp.")}
gold = Path(__file__).resolve().parents[1] / "golden" / \
    "tropopause_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
