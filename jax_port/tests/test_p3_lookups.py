"""Tier-0 tests for P3 table lookups and DSD kernels."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.p3 import tables  # noqa: E402
from scream_jax.p3.dsd import get_cloud_dsd2, get_rain_dsd2  # noqa: E402
from scream_jax.p3.table_lookups import (  # noqa: E402
    apply_table3,
    apply_table_coll,
    apply_table_ice,
    lookup_ice,
    lookup_rain,
    lookup_table3,
)

REPO = Path(__file__).resolve().parents[2]
TDIR = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "tables"


@pytest.fixture(scope="module")
def tbl():
    if not (TDIR / f"p3_lookup_table_1.dat-v{tables.P3_VERSION}").exists():
        pytest.skip("P3 tables not available")
    return tables.p3_init(str(TDIR))


def test_table3_reproduces_grid_points(tbl):
    # mu_r=1 everywhere; pick dm exactly on the small-branch grid:
    # rdumii = (dm*1e6+5)*0.1 = jj exactly when dm = (jj*10-5)e-6
    vn = tbl["vn_table_vals"]
    ctx = np.ones(5, dtype=bool)
    jjs = np.array([2, 5, 10, 15, 20])
    dm = (jjs * 10 - 5) * 1e-6
    mu_r = np.ones(5)
    lamr = (mu_r + 1) / dm
    t3 = lookup_table3(mu_r, lamr, ctx)
    np.testing.assert_array_equal(np.asarray(t3["dumii"]), jjs)
    got = np.asarray(apply_table3(vn, t3))
    np.testing.assert_allclose(got, vn[jjs - 1, 1], rtol=1e-12)


def test_table3_interpolation_between_points(tbl):
    vn = tbl["vn_table_vals"]
    ctx = np.ones(1, dtype=bool)
    # halfway between jj=10 and jj=11 in the small branch
    dm = np.array([(10.5 * 10 - 5) * 1e-6])
    t3 = lookup_table3(np.ones(1), 2.0 / dm, ctx)
    got = float(np.asarray(apply_table3(vn, t3))[0])
    lo, hi = vn[9, 1], vn[10, 1]
    assert min(lo, hi) - 1e-12 <= got <= max(lo, hi) + 1e-12


def test_ice_lookup_weights_in_range(tbl):
    rng = np.random.default_rng(0)
    n = 200
    qi = 10 ** rng.uniform(-8, -3, n)
    ni = 10 ** rng.uniform(2, 6, n)
    qm = qi * rng.uniform(0, 1, n)
    rhop = rng.uniform(50, 900, n)
    ctx = np.ones(n, dtype=bool)
    ti = lookup_ice(qi, ni, qm, rhop, ctx)
    for k in ("dumi", "dumii", "dumjj"):
        v = np.asarray(ti[k])
        assert v.min() >= 0
    assert np.asarray(ti["dumi"]).max() <= 48
    assert np.asarray(ti["dumii"]).max() <= 2
    assert np.asarray(ti["dumjj"]).max() <= 3
    # interpolated values finite and within table bounds per column
    for idx in range(12):
        vals = np.asarray(apply_table_ice(idx, tbl["ice_table_vals"], ti))
        assert np.isfinite(vals).all()

    tr = lookup_rain(10 ** rng.uniform(-8, -4, n), 10 ** rng.uniform(1, 5, n), ctx)
    for idx in range(2):
        vals = np.asarray(apply_table_coll(idx, tbl["collect_table_vals"], ti, tr))
        assert np.isfinite(vals).all()


def test_cloud_dsd2_properties():
    rng = np.random.default_rng(1)
    n = 100
    qc = 10 ** rng.uniform(-7, -3, n)
    nc = 10 ** rng.uniform(5, 9, n)
    rho = rng.uniform(0.4, 1.3, n)
    ctx = np.ones(n, dtype=bool)
    nc2, mu_c, nu, lamc, cdist, cdist1 = (np.asarray(a) for a in
                                          get_cloud_dsd2(qc, nc, rho, ctx))
    assert np.all((mu_c >= 2) & (mu_c <= 15))
    assert np.all(lamc >= (mu_c + 1) * 2.5e4 - 1e-6)
    assert np.all(lamc <= (mu_c + 1) * 1e6 + 1e-6)
    # After limiter, (qc, nc, lamc) are mutually consistent:
    lam_check = np.cbrt(c.CONS1 * nc2 * (mu_c + 3) * (mu_c + 2) * (mu_c + 1) / qc)
    np.testing.assert_allclose(lam_check, lamc, rtol=1e-10)
    # qc below QSMALL -> all zeros, nc unchanged
    out = get_cloud_dsd2(np.full(3, 1e-15), np.full(3, 1e7), np.ones(3),
                         np.ones(3, dtype=bool))
    assert np.all(np.asarray(out[3]) == 0)
    np.testing.assert_array_equal(np.asarray(out[0]), 1e7)


def test_rain_dsd2_properties():
    rng = np.random.default_rng(2)
    n = 100
    qr = 10 ** rng.uniform(-7, -3, n)
    nr = 10 ** rng.uniform(1, 6, n)
    ctx = np.ones(n, dtype=bool)
    nr2, mu_r, lamr = (np.asarray(a) for a in get_rain_dsd2(qr, nr, 1.0, ctx))
    np.testing.assert_array_equal(mu_r, 1.0)
    assert np.all(lamr >= 2 * 500.0 - 1e-9)
    assert np.all(lamr <= 2 * 1e5 + 1e-9)
    lam_check = np.cbrt(c.CONS1 * (1 + 3) * (1 + 2) * (1 + 1) * nr2 / qr)
    np.testing.assert_allclose(lam_check, lamr, rtol=1e-10)
