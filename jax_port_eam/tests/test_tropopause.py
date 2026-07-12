"""Tier-1 golden replay + Tier-0 properties for the tropopause-finder
port (eam_jax.tropopause).

Golden archive: harness/gen_tropopause_golden.py (run in the
scream-dev container) — 40x72 synthetic profile set (tropical /
midlat / polar / twmo-failure / isothermal / inversion / random)
through TWMO-only, TWMO+CLIMATE (production default), CLIMATE-only at
8 caldays covering every time-interp branch, and HYBSTOB(+CLIMATE)
(production chemistry chain), with a seasonally varying synthetic
climatology injected through the real tropopause_read_file path.

Level indices replay EXACTLY (golden is 1-based Fortran, the port is
0-based: +1). Pressures/temperatures/heights replay at rtol 1e-12
(gfortran -O2/aarch64 fma contraction, <= 1 ulp — see the golden
metadata); fillvalue lanes are compared exactly.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "tropopause_golden.npz"

from eam_jax import tropopause as tr  # noqa: E402

CALDAY_IDX = range(8)


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


def run_port(gold, primary, backup, calday):
    return tr.tropopause_find(
        gold["t"], gold["pmid"], gold["pint"], gold["zm"], gold["zi"],
        tropp_clim=gold["tropp_clim"], days=gold["days"],
        calday=calday, primary=primary, backup=backup)


def check_run(gold, name, state):
    lev, tp, tt, tz = (np.asarray(x) for x in state)
    glev = gold[f"{name}_lev"]
    found = glev != tr.NOTFOUND
    # level indices exact (Fortran 1-based; port 0-based)
    np.testing.assert_array_equal(np.where(lev != tr.NOTFOUND,
                                           lev + 1, tr.NOTFOUND), glev)
    for arr, key in ((tp, "p"), (tt, "t"), (tz, "z")):
        garr = gold[f"{name}_{key}"]
        np.testing.assert_allclose(arr[found], garr[found],
                                   rtol=1e-12, atol=0.0)
        # fillvalue lanes exact
        np.testing.assert_array_equal(arr[~found], garr[~found])


# ---------------------------------------------------------------------------
# Tier-1 golden replay
# ---------------------------------------------------------------------------
def test_twmo_only(gold):
    check_run(gold, "twmo", run_port(gold, "twmo", "none", 100.0))


@pytest.mark.parametrize("k", CALDAY_IDX)
def test_default_chain(gold, k):
    """Production default: TWMO primary + CLIMATE backup."""
    cd = gold["caldays"][k]
    check_run(gold, f"twmo_climate_cd{k}",
              run_port(gold, "twmo", "climate", cd))


@pytest.mark.parametrize("k", CALDAY_IDX)
def test_climate_only(gold, k):
    cd = gold["caldays"][k]
    check_run(gold, f"climate_cd{k}",
              run_port(gold, "climate", "none", cd))


def test_hybridstobie(gold):
    """Production chemistry chain: HYBSTOB primary (+CLIMATE backup,
    which can never engage)."""
    check_run(gold, "hybstob", run_port(gold, "hybridstobie", "none",
                                        100.0))
    check_run(gold, "hybstob_climate",
              run_port(gold, "hybridstobie", "climate", 100.0))


# ---------------------------------------------------------------------------
# Tier-0 properties
# ---------------------------------------------------------------------------
def test_bounds_and_level_consistency(gold):
    """Where found: tropP within [pmid top, pmid bottom], and tropP
    lies inside the pint bounds of the reported layer."""
    for primary, backup in (("twmo", "none"), ("climate", "none"),
                            ("twmo", "climate")):
        lev, tp, _, _ = (np.asarray(x) for x in
                         run_port(gold, primary, backup, 200.5))
        f = lev != tr.NOTFOUND
        assert np.all(tp[f] > gold["pmid"][f, 0])
        assert np.all(tp[f] < gold["pmid"][f, -1])
        # layer k contains tP: pint[k] <= tP < pint[k+1]
        rows = np.nonzero(f)[0]
        assert np.all(gold["pint"][rows, lev[f]] <= tp[f])
        assert np.all(tp[f] < gold["pint"][rows, lev[f] + 1])


def test_tropical_higher_than_polar(gold):
    """Monotone consistency on the synthetic set: tropical tropopause
    higher (lower pressure, greater height) than polar."""
    lev, tp, _, tz = (np.asarray(x) for x in
                      run_port(gold, "twmo", "none", 100.0))
    ptype = gold["ptype"]
    trop, pol = ptype == "tropical", ptype == "polar"
    assert (lev[trop] != tr.NOTFOUND).all()
    assert (lev[pol] != tr.NOTFOUND).all()
    assert tp[trop].max() < tp[pol].min()
    assert tz[trop].min() > tz[pol].max()


def test_backup_engages_exactly_on_notfound(gold):
    """The CLIMATE backup fills exactly the columns where twmo returns
    NOTFOUND, with the climate-only answer, and leaves the twmo answer
    untouched elsewhere."""
    cd = float(gold["caldays"][3])
    prim = run_port(gold, "twmo", "none", cd)
    chain = run_port(gold, "twmo", "climate", cd)
    clim = run_port(gold, "climate", "none", cd)
    nf = np.asarray(prim[0]) == tr.NOTFOUND
    assert nf.any() and not nf.all()
    assert np.array_equal(nf, gold["ptype"] == "twmo_fail")
    for a, b, c in zip(prim, chain, clim):
        a, b, c = np.asarray(a), np.asarray(b), np.asarray(c)
        np.testing.assert_array_equal(b[~nf], a[~nf])
        np.testing.assert_array_equal(b[nf], c[nf])


def test_hybridstobie_cannot_fail(gold):
    lev, tp, tt, tz = (np.asarray(x) for x in
                       run_port(gold, "hybridstobie", "none", 1.5))
    assert (lev != tr.NOTFOUND).all()
    assert (tp != tr.FILLVALUE).all()
    # midpoint outputs by construction
    rows = np.arange(lev.size)
    np.testing.assert_array_equal(tp, gold["pmid"][rows, lev])
    np.testing.assert_array_equal(tt, gold["t"][rows, lev])
    np.testing.assert_array_equal(tz, gold["zm"][rows, lev])


def test_notfound_outputs_are_fillvalue(gold):
    lev, tp, tt, tz = (np.asarray(x) for x in
                       run_port(gold, "twmo", "none", 100.0))
    nf = lev == tr.NOTFOUND
    assert nf.any()
    for arr in (tp, tt, tz):
        assert np.all(arr[nf] == tr.FILLVALUE)
