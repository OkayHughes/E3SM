"""EAMxx adapter for the JAX SHOC implementation.

`py_module_name` target for the shoc atmosphere process (host backend).
The `main` signature must match the py_module_call in the EAMXX_HAS_PYTHON
branch of eamxx_shoc_process_interface.cpp::run_impl on this branch — the
ENTIRE SHOC step (pre-process + shoc_main + post-process) is swapped, which
is exactly the unit validated against golden data by
jax_port/tests/test_shoc_golden.py.

Usage in an EAMxx input.yaml:
    shoc:
      py_module_name: shoc_jax
      py_module_path: <repo>/jax_port/scream_jax/adapters/eamxx
      py_backend: host

No physics lives here — see scream_jax/shoc/process.py.
"""

import sys
from pathlib import Path

import numpy as np

_PKG_ROOT = str(Path(__file__).resolve().parents[3])
if _PKG_ROOT not in sys.path:
    sys.path.insert(0, _PKG_ROOT)

from scream_jax.shoc.process import shoc_process_step  # noqa: E402


def init():
    """Called from SHOCMacrophysics::initialize_impl; nothing to set up."""
    pass


def main(dt, npbl,
         lambda_low, lambda_high, lambda_slope, lambda_thresh,
         thl2tune, qw2tune, qwthl2tune, w2tune, length_fac, c_diag_3rd_mom,
         coeff_kh, coeff_km, shoc_1p5tke,
         cell_length,
         # inputs
         p_mid, p_int, pseudo_density, omega, phis,
         surf_sens_flux, surf_evap, surf_mom_flux,
         # updated (in-place)
         T_mid, qv, qc, tke, horiz_winds, cldfrac_liq,
         sgs_buoy_flux, eddy_diff_mom,
         # computed (in-place)
         pbl_height, inv_qc_relvar, eddy_diff_heat, w_variance,
         cldfrac_liq_prev, ustar, obklen, thl_sec):
    """In-place shim over shoc_process_step (one AD subcycle, nadv = 1)."""
    out = shoc_process_step(
        float(dt), int(npbl), np.asarray(cell_length),
        float(lambda_low), float(lambda_high), float(lambda_slope),
        float(lambda_thresh), float(thl2tune), float(qw2tune),
        float(qwthl2tune), float(w2tune), float(length_fac),
        float(c_diag_3rd_mom), float(coeff_kh), float(coeff_km),
        bool(shoc_1p5tke), False,
        np.asarray(T_mid), np.asarray(p_mid), np.asarray(p_int),
        np.asarray(pseudo_density), np.asarray(omega), np.asarray(phis),
        np.asarray(surf_sens_flux), np.asarray(surf_evap),
        np.asarray(surf_mom_flux)[:, 0], np.asarray(surf_mom_flux)[:, 1],
        np.asarray(qv), np.asarray(qc), np.asarray(tke),
        np.asarray(horiz_winds)[:, 0, :], np.asarray(horiz_winds)[:, 1, :],
        np.asarray(cldfrac_liq), np.asarray(sgs_buoy_flux),
        np.asarray(eddy_diff_mom))

    # Updated fields
    T_mid[...] = np.asarray(out["T_mid"])
    qv[...] = np.asarray(out["qv"])
    qc[...] = np.asarray(out["qc"])
    tke[...] = np.asarray(out["tke"])
    horiz_winds[:, 0, :] = np.asarray(out["u_wind"])
    horiz_winds[:, 1, :] = np.asarray(out["v_wind"])
    cldfrac_liq[...] = np.asarray(out["cldfrac_liq"])
    sgs_buoy_flux[...] = np.asarray(out["sgs_buoy_flux"])
    eddy_diff_mom[...] = np.asarray(out["eddy_diff_mom"])

    # Computed fields
    pbl_height[...] = np.asarray(out["pbl_height"])
    inv_qc_relvar[...] = np.asarray(out["inv_qc_relvar"])
    eddy_diff_heat[...] = np.asarray(out["eddy_diff_heat"])
    w_variance[...] = np.asarray(out["w_variance"])
    cldfrac_liq_prev[...] = np.asarray(out["cldfrac_liq_prev"])
    ustar[...] = np.asarray(out["ustar"])
    obklen[...] = np.asarray(out["obklen"])
    thl_sec[...] = np.asarray(out["thl_sec"])
