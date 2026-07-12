"""Tier-0 tests for P3 lookup tables (golden: the C++-written .dat8 caches)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.p3 import tables  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
TDIR = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "tables"


def _tables():
    if not (TDIR / f"p3_lookup_table_1.dat-v{tables.P3_VERSION}").exists():
        pytest.skip(f"P3 tables not found under {TDIR}")
    return tables.p3_init(str(TDIR))


def test_shapes_and_finiteness():
    t = _tables()
    assert t["ice_table_vals"].shape == (5, 4, 50, 12)
    assert t["collect_table_vals"].shape == (5, 4, 50, 30, 2)
    assert t["vn_table_vals"].shape == (300, 10)
    assert t["dnu_table_vals"].shape == (16,)
    for v in t.values():
        assert np.isfinite(v).all()


def test_computed_rain_tables_match_cpp_binaries():
    t = _tables()
    mu_f, vn_f, vm_f, revap_f = tables.read_computed_tables(str(TDIR))
    np.testing.assert_array_equal(t["mu_r_table_vals"], mu_f)
    for mine, ref in ((t["vn_table_vals"], vn_f), (t["vm_table_vals"], vm_f),
                      (t["revap_table_vals"], revap_f)):
        np.testing.assert_allclose(mine, ref, rtol=1e-12)
