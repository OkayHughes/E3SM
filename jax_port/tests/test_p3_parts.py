"""Tier-0 smoke tests for p3_main_part1 and part3."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.p3 import tables  # noqa: E402
from scream_jax.p3.main_part1 import p3_main_part1  # noqa: E402
from scream_jax.p3.main_part3 import p3_main_part3  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
TDIR = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "tables"

OPTS = dict(spa_ccn_to_nc_factor=1.0, spa_ccn_to_nc_exponent=1.0,
            constant_mu_rain=1.0, min_rime_rho=50.0, max_rime_rho=900.0)


def _cols(ncol=4, nlev=24, seed=0):
    rng = np.random.default_rng(seed)
    pres = np.linspace(300e2, 1000e2, nlev)[None] * np.ones((ncol, 1))
    dpres = np.gradient(pres, axis=-1)
    T = 210 + (300 - 210) * (pres / pres.max()) ** 1.5
    dz = dpres * c.Rair * T / (pres * c.gravit)
    exner = (pres / c.P0) ** (c.Rair / c.CP)
    return rng, pres, dpres, T, dz, exner


def test_part1_flags_and_clipping():
    rng, pres, dpres, T, dz, exner = _cols()
    ncol, nlev = T.shape
    shape = (ncol, nlev)
    qc = np.where(rng.uniform(size=shape) < 0.5, 1e-4, 1e-16)  # some dry levels
    out = p3_main_part1(
        True, False, 300.0, pres, dpres, dz, np.zeros(shape), np.zeros(shape),
        1 / exner, exner, np.ones(shape), np.ones(shape), np.ones(shape),
        T, np.full(shape, 5e-3), T / exner,
        qc, np.full(shape, 1e8), np.full(shape, 1e-5), np.full(shape, 1e3),
        np.full(shape, 1e-5), np.full(shape, 1e4), np.full(shape, 5e-6),
        np.full(shape, 1e-8), OPTS)
    # dry levels: qc clipped to zero and vapor increased
    qc2 = np.asarray(out["qc"])
    assert np.all(qc2[qc < c.QSMALL] == 0)
    assert np.all(np.asarray(out["qv"]) >= 5e-3 - 1e-15)
    # cold upper levels are supersaturable -> nucleation possible
    assert np.asarray(out["nucleation_possible"]).all()
    assert np.asarray(out["hydrometeors_present"]).all()
    for k, v in out.items():
        assert np.all(np.isfinite(np.asarray(v))), k


def test_part3_diagnostics_and_cleanup():
    if not (TDIR / f"p3_lookup_table_1.dat-v{tables.P3_VERSION}").exists():
        pytest.skip("P3 tables not available")
    tbl = tables.p3_init(str(TDIR))
    rng, pres, dpres, T, dz, exner = _cols(seed=1)
    ncol, nlev = T.shape
    shape = (ncol, nlev)
    ones = np.ones(shape)
    rho = dpres / dz / c.gravit
    qi = np.where(rng.uniform(size=shape) < 0.5, 2e-4, 0.0)

    out = p3_main_part3(
        740e3, tbl["dnu_table_vals"], tbl["ice_table_vals"],
        1 / exner, ones, ones, ones, rho, 1 / rho,
        (600e2 / (c.Rair * 253.15) / rho) ** 0.54,
        np.full(shape, 5e-3), T / exner,
        np.full(shape, 1e-4), np.full(shape, 1e8),
        np.full(shape, 1e-5), np.full(shape, 1e3),
        qi, np.full(shape, 1e4), qi * 0.3, qi * 0.3 / 400,
        np.zeros(shape), np.zeros(shape), np.zeros(shape), OPTS)

    ice = qi > 0
    assert np.all(np.asarray(out["diag_eff_radius_qi"])[ice] > 0)
    assert np.all(np.asarray(out["rho_qi"])[ice] > 0)
    assert np.all(np.asarray(out["qi"])[~ice] == 0)
    assert np.all(np.asarray(out["diag_eff_radius_qc"]) > 0)   # qc present everywhere
    # reflectivity finite everywhere (floors at 1e-22 m^6/m^3)
    assert np.all(np.isfinite(np.asarray(out["diag_equiv_reflectivity"])))
    for k, v in out.items():
        assert np.all(np.isfinite(np.asarray(v))), k
