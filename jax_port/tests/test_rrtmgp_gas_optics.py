"""Tier-0 tests for RRTMGP k-distribution loading and gas optics."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

REPO = Path(__file__).resolve().parents[2]
DDIR = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "init"
GASES = ["h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"]

from scream_jax.rrtmgp.coefficients import load_kdist  # noqa: E402
from scream_jax.rrtmgp import gas_optics as go  # noqa: E402


@pytest.fixture(scope="module")
def kdists():
    sw = DDIR / "rrtmgp-data-sw-g112-210809.nc"
    lw = DDIR / "rrtmgp-data-lw-g128-210809.nc"
    if not sw.exists():
        pytest.skip("rrtmgp coefficient files not available")
    return load_kdist(str(sw), GASES), load_kdist(str(lw), GASES)


def _atmosphere(ncol=8, nlay=72):
    plev = np.linspace(100.0, 101000.0, nlay + 1)[None, :] * np.ones((ncol, 1))
    play = 0.5 * (plev[:, :-1] + plev[:, 1:])
    tlay = 200 + 100 * (play / play.max()) ** 0.3
    tlev = 200 + 100 * (plev / plev.max()) ** 0.3
    vmr = np.zeros((ncol, nlay, 8))
    vmr[..., 0] = 5e-3 * (play / play.max()) ** 2
    vmr[..., 1] = 400e-6
    vmr[..., 2] = 5e-8 + 1e-5 * np.exp(-((np.log(play) - np.log(2000)) ** 2))
    vmr[..., 3] = 320e-9
    vmr[..., 4] = 1e-7
    vmr[..., 5] = 1.8e-6
    vmr[..., 6] = 0.21
    vmr[..., 7] = 0.78
    return play, plev, tlay, tlev, vmr


def test_kdist_shapes(kdists):
    kd_sw, kd_lw = kdists
    assert (kd_sw["ngpt"], kd_sw["nband"], kd_sw["nflav"]) == (112, 14, 9)
    assert (kd_lw["ngpt"], kd_lw["nband"], kd_lw["nflav"]) == (128, 16, 10)
    assert kd_sw["kmajor"].shape == (112, 9, 60, 14)
    assert kd_lw["planck_frac"].shape == (128, 9, 60, 14)
    assert kd_sw["krayl"].shape == (112, 9, 14, 2)
    # gpt ranges are contiguous and cover all g-points
    for kd in kdists:
        b2g = kd["band2gpt"]
        assert b2g[0, 0] == 0 and b2g[1, -1] == kd["ngpt"] - 1
        assert np.all(b2g[0, 1:] == b2g[1, :-1] + 1)
        # every minor entry indexes a valid kminor row
        for atm in ("lower", "upper"):
            assert kd[f"minor_entry_k_{atm}"].max() \
                == kd[f"kminor_{atm}"].shape[0] - 1


def test_gas_optics_sw(kdists):
    kd_sw, _ = kdists
    play, plev, tlay, tlev, vmr = _atmosphere()
    optics, toa, col_gas = go.gas_optics_sw(kd_sw, play, plev, tlay, vmr)
    tau = np.asarray(optics["tau"])
    ssa = np.asarray(optics["ssa"])
    assert np.isfinite(tau).all() and (tau >= 0).all()
    assert (ssa >= 0).all() and (ssa <= 1 + 1e-15).all()
    # the g112 dataset integrates to the known total solar irradiance
    np.testing.assert_allclose(np.asarray(toa)[0].sum(), 1360.3756,
                               rtol=1e-6)
    # visible band (16000-22650 cm^-1) in the mid-troposphere (away from
    # the synthetic ozone maximum): gas absorption is weak there, so ssa
    # (the Rayleigh fraction) should be near 1
    kd_sw, _ = kdists
    g0, g1 = kd_sw["band2gpt"][:, 10]
    assert ssa[0, 35, g0:g1 + 1].max() > 0.9


def test_gas_optics_lw(kdists):
    _, kd_lw = kdists
    play, plev, tlay, tlev, vmr = _atmosphere()
    tsfc = tlev[:, -1]
    optics, src, col_gas = go.gas_optics_lw(kd_lw, play, plev, tlay, tsfc,
                                            vmr, tlev=tlev)
    tau = np.asarray(optics["tau"])
    assert np.isfinite(tau).all() and (tau >= 0).all()
    for k in ("sfc_src", "lay_src", "lev_src_inc", "lev_src_dec"):
        a = np.asarray(src[k])
        assert np.isfinite(a).all() and (a >= 0).all(), k
    # Planck sources increase toward the (warmer) surface
    lay = np.asarray(src["lay_src"])
    assert lay[0, -1].sum() > lay[0, 0].sum()
    # col_gas: dry air column ~ 2.1e25 molec/cm2 total
    total_col_dry = np.asarray(col_gas)[0, :, 0].sum()
    assert 1e25 < total_col_dry < 3e25


def test_gas_optics_vs_cpp_golden(kdists):
    """BFB-level comparison against the real C++ GasOpticsRRTMGPK run on
    identical inputs (golden/rrtmgp_gasopt_cpp_8x72.npz, generated with
    harness/cpp_dumpers/gasopt_dump.cpp in the dev container)."""
    kd_sw, kd_lw = kdists
    g = REPO / "jax_port" / "golden" / "rrtmgp_gasopt_cpp_8x72.npz"
    if not g.exists():
        pytest.skip("gas optics C++ golden archive not available")
    z = np.load(g)
    tlay_lim = np.clip(z["tlay"], kd_sw["temp_ref_min"], kd_sw["temp_ref_max"])
    tlev_lim = np.clip(z["tlev"], kd_lw["temp_ref_min"], kd_lw["temp_ref_max"])

    optics, toa, _ = go.gas_optics_sw(kd_sw, z["play"], z["plev"], tlay_lim,
                                      z["vmr"])
    optics_lw, src, _ = go.gas_optics_lw(kd_lw, z["play"], z["plev"],
                                         tlay_lim, z["tsfc"], z["vmr"],
                                         tlev=tlev_lim)
    for mine, name in ((optics["tau"], "sw_tau"), (optics["ssa"], "sw_ssa"),
                       (optics["g"], "sw_g"), (toa, "sw_toa"),
                       (optics_lw["tau"], "lw_tau"),
                       (src["lay_src"], "lw_lay_src"),
                       (src["lev_src_inc"], "lw_lev_src_inc"),
                       (src["lev_src_dec"], "lw_lev_src_dec"),
                       (src["sfc_src"], "lw_sfc_src")):
        ref = z[name]
        scale = max(np.abs(ref).max(), 1e-300)
        err = np.abs(np.asarray(mine) - ref).max() / scale
        assert err < 1e-14, f"{name}: {err:.3e}"
