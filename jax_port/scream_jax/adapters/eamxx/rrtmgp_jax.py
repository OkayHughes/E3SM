"""EAMxx adapter for the JAX RRTMGP implementation.

`py_module_name` target for the rrtmgp atmosphere process (host
backend). The `main` signature must match the py_module_call in the
EAMXX_HAS_PYTHON branch of
eamxx_rrtmgp_process_interface.cpp::run_impl on this branch — the
ENTIRE radiation step is swapped, which is exactly the unit validated
against golden data by jax_port/tests/test_rrtmgp_golden.py.

No physics lives here — see scream_jax/rrtmgp/process.py.
"""

import sys
from pathlib import Path

import numpy as np

_PKG_ROOT = str(Path(__file__).resolve().parents[3])
if _PKG_ROOT not in sys.path:
    sys.path.insert(0, _PKG_ROOT)

from scream_jax.rrtmgp.coefficients import load_kdist  # noqa: E402
from scream_jax.rrtmgp.cloud_optics import load_cloud_optics  # noqa: E402
from scream_jax.rrtmgp.process import (GAS_NAMES,  # noqa: E402
                                       rrtmgp_process_step)

_KD_SW = _KD_LW = _CO_SW = _CO_LW = None

# field names in the exact order the C++ passes the computed buffers
_COMPUTED = (
    "cldfrac_rad", "cosine_solar_zenith_angle",
    "SW_flux_up", "SW_flux_dn", "SW_flux_dn_dir",
    "LW_flux_up", "LW_flux_dn",
    "SW_clnclrsky_flux_up", "SW_clnclrsky_flux_dn",
    "SW_clnclrsky_flux_dn_dir",
    "SW_clrsky_flux_up", "SW_clrsky_flux_dn", "SW_clrsky_flux_dn_dir",
    "SW_clnsky_flux_up", "SW_clnsky_flux_dn", "SW_clnsky_flux_dn_dir",
    "LW_clnclrsky_flux_up", "LW_clnclrsky_flux_dn",
    "LW_clrsky_flux_up", "LW_clrsky_flux_dn",
    "LW_clnsky_flux_up", "LW_clnsky_flux_dn",
    "sfc_flux_dir_vis", "sfc_flux_dir_nir",
    "sfc_flux_dif_vis", "sfc_flux_dif_nir",
    "sfc_flux_sw_net", "sfc_flux_lw_dn",
    "cldlow", "cldmed", "cldhgh", "cldtot",
    "dtau067", "dtau105", "sunlit_mask",
    "T_mid_at_cldtop", "p_mid_at_cldtop",
    "cldfrac_ice_at_cldtop", "cldfrac_liq_at_cldtop",
    "cldfrac_tot_at_cldtop", "cdnc_at_cldtop",
    "eff_radius_qc_at_cldtop", "eff_radius_qi_at_cldtop",
    "h2o_volume_mix_ratio", "co2_volume_mix_ratio",
    "n2o_volume_mix_ratio", "co_volume_mix_ratio",
    "ch4_volume_mix_ratio", "o2_volume_mix_ratio",
    "n2_volume_mix_ratio",
)


def init(coeff_sw, coeff_lw, cloud_sw, cloud_lw):
    """Called from RRTMGPRadiation::initialize_impl with the four
    coefficient file paths."""
    global _KD_SW, _KD_LW, _CO_SW, _CO_LW
    _KD_SW = load_kdist(str(coeff_sw), GAS_NAMES)
    _KD_LW = load_kdist(str(coeff_lw), GAS_NAMES)
    _CO_SW = load_cloud_optics(str(cloud_sw))
    _CO_LW = load_cloud_optics(str(cloud_lw))


def main(dt, nstep, year, calday,
         rad_frequency, orbital_year, orbital_eccen, orbital_obliq,
         orbital_mvelp, fixed_tsi, fixed_zenith,
         co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, n2vmr, covmr,
         do_subcol_sampling, extra_clnclrsky_diag, extra_clnsky_diag,
         lat, lon,
         # inputs
         p_mid, p_int, pseudo_density,
         sfc_alb_dir_vis, sfc_alb_dir_nir, sfc_alb_dif_vis, sfc_alb_dif_nir,
         qv, qc, nc, qi, cldfrac_tot, eff_radius_qc, eff_radius_qi,
         surf_lw_flux_up, o3_volume_mix_ratio,
         # updated
         T_mid, rad_heating_pdel,
         # computed (in the _COMPUTED order)
         *computed):
    assert _KD_SW is not None, "rrtmgp_jax.init() was not called"
    assert len(computed) == len(_COMPUTED), \
        f"expected {len(_COMPUTED)} computed fields, got {len(computed)}"

    params = {
        "rad_frequency": int(rad_frequency),
        "orbital_year": int(orbital_year),
        "orbital_eccentricity": float(orbital_eccen),
        "orbital_obliquity": float(orbital_obliq),
        "orbital_mvelp": float(orbital_mvelp),
        "fixed_total_solar_irradiance": float(fixed_tsi),
        "fixed_solar_zenith_angle": float(fixed_zenith),
        "co2vmr": float(co2vmr), "n2ovmr": float(n2ovmr),
        "ch4vmr": float(ch4vmr), "f11vmr": float(f11vmr),
        "f12vmr": float(f12vmr), "n2vmr": float(n2vmr),
        "covmr": float(covmr),
        "do_subcol_sampling": bool(do_subcol_sampling),
        "extra_clnclrsky_diag": bool(extra_clnclrsky_diag),
        "extra_clnsky_diag": bool(extra_clnsky_diag),
    }

    out = rrtmgp_process_step(
        _KD_SW, _KD_LW, _CO_SW, _CO_LW, params,
        float(dt), int(nstep), int(year), float(calday),
        np.asarray(lat), np.asarray(lon),
        np.asarray(T_mid), np.asarray(p_mid), np.asarray(p_int),
        np.asarray(pseudo_density),
        np.asarray(sfc_alb_dir_vis), np.asarray(sfc_alb_dir_nir),
        np.asarray(sfc_alb_dif_vis), np.asarray(sfc_alb_dif_nir),
        np.asarray(qv), np.asarray(qc), np.asarray(nc), np.asarray(qi),
        np.asarray(cldfrac_tot),
        np.asarray(eff_radius_qc), np.asarray(eff_radius_qi),
        np.asarray(surf_lw_flux_up), np.asarray(o3_volume_mix_ratio),
        np.asarray(rad_heating_pdel))

    T_mid[...] = np.asarray(out["T_mid"])
    rad_heating_pdel[...] = np.asarray(out["rad_heating_pdel"])
    for buf, name in zip(computed, _COMPUTED):
        if name in out:
            buf[...] = np.asarray(out[name]).astype(np.asarray(buf).dtype)
