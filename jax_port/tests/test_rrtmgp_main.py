"""Tier-0/BFB tests for the full rrtmgp_main driver chain (cloud optics,
MCICA, RTE solvers) against the real C++ run on identical inputs
(golden/rrtmgp_main_cpp_8x72.npz, generated with
harness/cpp_dumpers/rrtmgpmain_dump.cpp in the dev container)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

REPO = Path(__file__).resolve().parents[2]
DDIR = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "init"
GOLDEN = REPO / "jax_port" / "golden" / "rrtmgp_main_cpp_8x72.npz"
GASES = ["h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"]

from scream_jax.rrtmgp.coefficients import load_kdist  # noqa: E402
from scream_jax.rrtmgp.cloud_optics import load_cloud_optics  # noqa: E402
from scream_jax.rrtmgp.mcica import jsf64_value  # noqa: E402
from scream_jax.rrtmgp import interface as ifc  # noqa: E402


def test_jsf64_matches_cpp():
    import jax.numpy as jnp
    seeds = jnp.array([0, 1, 42, 123456789, 999999999, 987654321987],
                      dtype=jnp.uint64)
    ref = np.array([2.93850194127627862e-01, 6.81447781867085323e-01,
                    6.46265019899146753e-01, 7.14257566607624628e-01,
                    7.36117488556373989e-02, 1.55679042839169581e-02])
    np.testing.assert_array_equal(np.asarray(jsf64_value(seeds)), ref)


def test_rrtmgp_main_vs_cpp_golden():
    if not GOLDEN.exists() or not (DDIR / "rrtmgp-data-sw-g112-210809.nc").exists():
        pytest.skip("rrtmgp golden/coefficient files not available")
    z = np.load(GOLDEN)
    ncol, nlay = 8, 72

    kd_sw = load_kdist(str(DDIR / "rrtmgp-data-sw-g112-210809.nc"), GASES)
    kd_lw = load_kdist(str(DDIR / "rrtmgp-data-lw-g128-210809.nc"), GASES)
    co_sw = load_cloud_optics(str(DDIR / "rrtmgp-cloud-optics-coeffs-sw.nc"))
    co_lw = load_cloud_optics(str(DDIR / "rrtmgp-cloud-optics-coeffs-lw.nc"))

    alb_dir, alb_dif = ifc.compute_band_by_band_surface_albedos(
        kd_sw["band_lims_wvn"], z["adv"], z["anv"], z["afv"], z["afn"])
    np.testing.assert_array_equal(np.asarray(alb_dir), z["albdir"])

    zsw = np.zeros((ncol, nlay, 14))
    zlw = np.zeros((ncol, nlay, 16))
    sw, lw, cld_tau = ifc.rrtmgp_main(
        kd_sw, kd_lw, co_sw, co_lw,
        z["play"], z["tlay"], z["plev"], z["tlev"], z["vmr"],
        alb_dir, alb_dif, z["mu0"],
        z["lwp"], z["iwp"], z["rel"], z["rei"], z["cldfrac"],
        zsw, zsw, zsw, zlw, 1.0,
        extra_clnclrsky_diag=True, extra_clnsky_diag=True)

    checks = [
        (cld_tau["sw_bnd"], "ct_swb", 1e-14), (cld_tau["lw_bnd"], "ct_lwb", 1e-14),
        (cld_tau["sw_gpt"], "ct_swg", 1e-14), (cld_tau["lw_gpt"], "ct_lwg", 1e-14),
        (sw["allsky"]["flux_up"], "swu", 1e-11),
        (sw["allsky"]["flux_dn"], "swd", 1e-11),
        (sw["allsky"]["flux_dn_dir"], "swdd", 1e-11),
        (sw["clrsky"]["flux_up"], "csu", 1e-11),
        (sw["clnclrsky"]["flux_up"], "ccsu", 1e-11),
        (sw["clnsky"]["flux_up"], "cnu", 1e-11),
        (lw["allsky"]["flux_up"], "lwu", 1e-12),
        (lw["allsky"]["flux_dn"], "lwd", 1e-12),
        (lw["clrsky"]["flux_up"], "lcu", 1e-12),
        (lw["clnsky"]["flux_dn"], "lnd", 1e-12),
        (sw["allsky"]["bnd_flux_up"], "swbu", 1e-11),
        (sw["allsky"]["bnd_flux_dn_dir"], "swbdd", 1e-11),
        (lw["allsky"]["bnd_flux_up"], "lwbu", 1e-12),
    ]
    for mine, name, tol in checks:
        ref = z[name]
        scale = max(np.abs(ref).max(), 1e-300)
        err = np.abs(np.asarray(mine) - ref).max() / scale
        assert err < tol, f"{name}: {err:.3e}"
