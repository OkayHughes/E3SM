"""Tier-0 smoke tests for the assembled p3_main.

Checks finiteness/positivity, total-water conservation including surface
precipitation, the two early-exit tiers, and diagnostic init values in
skipped columns. Golden replay against EAMxx is the correctness test.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.p3 import DEFAULT_OPTS, tables  # noqa: E402
from scream_jax.p3.main import p3_main  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
TDIR = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "tables"

OPTS = DEFAULT_OPTS


@pytest.fixture(scope="module")
def tbl():
    if not (TDIR / f"p3_lookup_table_1.dat-v{tables.P3_VERSION}").exists():
        pytest.skip("P3 tables not available")
    return tables.p3_init(str(TDIR))


def _setup(ncol=6, nlev=48, seed=0):
    rng = np.random.default_rng(seed)
    shape = (ncol, nlev)
    pres = np.linspace(100e2, 1000e2, nlev)[None] * np.ones((ncol, 1))
    dpres = np.gradient(pres, axis=-1)
    T = 200 + (295 - 200) * (pres / pres.max()) ** 1.2
    dz = dpres * c.Rair * T / (pres * c.gravit)
    exner = (pres / c.P0) ** (c.Rair / c.CP)
    th = T / exner
    cld = np.clip(rng.uniform(0.1, 1.0, shape), 1e-4, 1.0)

    qv = 0.8 * 0.622 * 610.78 * np.exp(
        17.27 * (T - 273.15) / (T - 35.85)) / pres  # ~80% RH wrt liquid
    qc = np.where((rng.uniform(size=shape) < 0.4) & (T > 240), 2e-4, 0.0)
    qr = np.where((rng.uniform(size=shape) < 0.3) & (T > 260), 1e-4, 0.0)
    qi = np.where((rng.uniform(size=shape) < 0.4) & (T < 265), 1e-4, 0.0)
    nc = np.where(qc > 0, 5e7, 0.0)
    nr = np.where(qr > 0, 1e5, 0.0)
    ni = np.where(qi > 0, 5e4, 0.0)
    qm = 0.2 * qi
    bm = qm / 400.0
    return dict(shape=shape, pres=pres, dpres=dpres, T=T, dz=dz, exner=exner,
                th=th, cld=cld, qv=qv, qc=qc, nc=nc, qr=qr, nr=nr,
                qi=qi, ni=ni, qm=qm, bm=bm)


def _run(s, dt=300.0):
    shape = s["shape"]
    z = np.zeros(shape)
    tbl = tables.p3_init(str(TDIR))
    return p3_main(
        dt, True, False, True, False, False,
        s["qc"], s["nc"], s["qr"], s["nr"], s["qi"], s["qm"], s["ni"],
        s["bm"], s["qv"], s["th"],
        z, z, np.full(shape, 5e4), np.ones(shape),
        s["cld"], s["cld"], np.maximum.accumulate(s["cld"], axis=-1),
        s["pres"], s["dz"], s["dpres"], (c.P0 / s["pres"]) ** (c.Rair / c.CP),
        s["qv"], s["T"], z, z, z, tbl, OPTS)


def test_p3_main_smoke_and_conservation(tbl):
    s = _setup()
    dt = 300.0
    out = _run(s, dt)

    for k, v in out.items():
        assert np.all(np.isfinite(np.asarray(v))), k
    for k in ("qv", "qc", "nc", "qr", "nr", "qi", "ni", "qm", "bm"):
        assert np.all(np.asarray(out[k]) >= 0), k

    # total water (vapor+liquid+rain+ice) conserved up to surface precip
    rho_dz = s["dpres"] / c.gravit  # rho*dz = dp/g
    def tot(qv, qc, qr, qi):
        return np.sum((qv + qc + qr + qi) * rho_dz, axis=-1)
    before = tot(np.maximum(s["qv"], 0), s["qc"], s["qr"], s["qi"])
    after = tot(*(np.asarray(out[k]) for k in ("qv", "qc", "qr", "qi")))
    surf = (np.asarray(out["precip_liq_surf"])
            + np.asarray(out["precip_ice_surf"])) * c.RHO_H2O * dt
    # rho used inside sedimentation is dp/(g*dz)*... consistent: closure
    np.testing.assert_allclose(after + surf, before, rtol=5e-7)

    assert np.all(np.asarray(out["precip_liq_flux"]) >= 0)
    assert np.all(np.asarray(out["precip_ice_flux"]) == 0.0)


def test_p3_main_early_exit_column(tbl):
    s = _setup(seed=3)
    # column 0: warm-ish, bone dry, no hydrometeors, subsaturated ->
    # nucleationPossible false and hydrometeorsPresent false
    for k in ("qc", "nc", "qr", "nr", "qi", "ni", "qm", "bm"):
        s[k] = s[k].copy()
        s[k][0] = 0.0
    s["qv"] = s["qv"].copy()
    s["qv"][0] = 1e-6

    out = _run(s)

    # early-exit column: prognostics keep part1 values (qv >= 0, no change
    # from part2/sed), diagnostics keep their init values
    np.testing.assert_array_equal(np.asarray(out["qc"])[0], 0.0)
    np.testing.assert_array_equal(np.asarray(out["diag_eff_radius_qc"])[0],
                                  10.0e-6)
    np.testing.assert_array_equal(np.asarray(out["diag_eff_radius_qi"])[0],
                                  25.0e-6)
    np.testing.assert_array_equal(
        np.asarray(out["diag_equiv_reflectivity"])[0], -99.0)
    np.testing.assert_array_equal(np.asarray(out["precip_liq_flux"])[0], 0.0)
    assert np.asarray(out["precip_liq_surf"])[0] == 0.0
    assert np.asarray(out["precip_ice_surf"])[0] == 0.0
    # theta unchanged in the early-exit column (no processes ran)
    np.testing.assert_allclose(np.asarray(out["th_atm"])[0], s["th"][0],
                               rtol=1e-14)


def test_p3_main_cold_homogeneous_freezing(tbl):
    s = _setup(seed=5)
    # make one column entirely below the homogeneous freezing point with
    # cloud water present: all liquid must exit as ice (no qc/qr survive)
    s["T"] = s["T"].copy()
    s["T"][1] = 220.0
    s["th"] = s["T"] / s["exner"]
    s["qc"] = s["qc"].copy()
    s["qc"][1] = 1e-4
    s["nc"] = s["nc"].copy()
    s["nc"][1] = 5e7

    out = _run(s)
    qc1 = np.asarray(out["qc"])[1]
    qr1 = np.asarray(out["qr"])[1]
    assert np.all(qc1 == 0.0)
    assert np.all(qr1 == 0.0)
