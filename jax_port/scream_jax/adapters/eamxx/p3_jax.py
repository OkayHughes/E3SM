"""EAMxx adapter for the JAX P3 implementation.

`py_module_name` target for the p3 atmosphere process (host backend).
The `main` signature must match the py_module_call in the
EAMXX_HAS_PYTHON branch of eamxx_p3_run.cpp::run_impl on this branch —
the ENTIRE P3 step (pre-process + p3_main + post-process) is swapped,
which is exactly the unit validated against golden data by
jax_port/tests/test_p3_golden.py.

Usage in an EAMxx input.yaml:
    p3:
      py_module_name: p3_jax
      py_module_path: <repo>/jax_port/scream_jax/adapters/eamxx
      py_backend: host

No physics lives here — see scream_jax/p3/process.py.
"""

import sys
from pathlib import Path

import numpy as np

_PKG_ROOT = str(Path(__file__).resolve().parents[3])
if _PKG_ROOT not in sys.path:
    sys.path.insert(0, _PKG_ROOT)

from scream_jax.p3 import tables  # noqa: E402
from scream_jax.p3.process import p3_process_step  # noqa: E402

_OPT_NAMES = (
    "max_total_ni", "autoconversion_prefactor", "autoconversion_qc_exponent",
    "autoconversion_nc_exponent", "autoconversion_radius",
    "accretion_prefactor", "accretion_qc_exponent", "accretion_qr_exponent",
    "rain_selfcollection_prefactor", "rain_selfcollection_breakup_diameter",
    "constant_mu_rain", "spa_ccn_to_nc_factor", "spa_ccn_to_nc_exponent",
    "cldliq_to_ice_collection_factor", "rain_to_ice_collection_factor",
    "min_rime_rho", "max_rime_rho", "immersion_freezing_exponent",
    "deposition_nucleation_exponent", "ice_sedimentation_factor",
)

_TBL = None


def init(table_file):
    """Called from P3Microphysics::initialize_impl with the full path of
    the ice lookup table; loads/computes all P3 tables once."""
    global _TBL
    _TBL = tables.p3_init(str(Path(table_file).parent))


def main(dt, predict_nc, do_ice_production,
         set_cld_frac_l_to_one, set_cld_frac_i_to_one, set_cld_frac_r_to_one,
         *args):
    """In-place shim over p3_process_step (one AD subcycle)."""
    assert _TBL is not None, "p3_jax.init() was not called"
    opts = {n: float(v) for n, v in zip(_OPT_NAMES, args[:20])}
    (p_mid, p_dry_mid, pseudo_density, pseudo_density_dry, cldfrac_tot,
     nc_nuceat_tend, ni_activated, inv_qc_relvar,
     T_mid, qv, qc, nc, qr, nr, qi, ni, qm, bm,
     qv_prev_micro_step, T_prev_micro_step,
     precip_liq_surf_mass, precip_ice_surf_mass,
     eff_radius_qc, eff_radius_qi, eff_radius_qr,
     precip_total_tend, nevapr, diag_equiv_reflectivity,
     micro_liq_ice_exchange, micro_vap_liq_exchange,
     micro_vap_ice_exchange, rainfrac) = args[20:]

    out = p3_process_step(
        float(dt),
        bool(predict_nc), False,          # prescribed_ccn rejected in C++
        bool(do_ice_production), False, False,  # hetfrz/sep-frac rejected
        bool(set_cld_frac_l_to_one), bool(set_cld_frac_i_to_one),
        bool(set_cld_frac_r_to_one),
        np.asarray(T_mid), np.asarray(p_mid), np.asarray(p_dry_mid),
        np.asarray(pseudo_density), np.asarray(pseudo_density_dry),
        np.asarray(cldfrac_tot),
        np.asarray(qv), np.asarray(qc), np.asarray(nc), np.asarray(qr),
        np.asarray(nr), np.asarray(qi), np.asarray(qm), np.asarray(ni),
        np.asarray(bm),
        np.asarray(qv_prev_micro_step), np.asarray(T_prev_micro_step),
        np.asarray(nc_nuceat_tend), None, np.asarray(ni_activated),
        np.asarray(inv_qc_relvar),
        np.asarray(precip_liq_surf_mass), np.asarray(precip_ice_surf_mass),
        _TBL, opts)

    # Updated fields
    for buf, key in ((T_mid, "T_mid"), (qv, "qv"), (qc, "qc"), (nc, "nc"),
                     (qr, "qr"), (nr, "nr"), (qi, "qi"), (ni, "ni"),
                     (qm, "qm"), (bm, "bm"),
                     (qv_prev_micro_step, "qv_prev_micro_step"),
                     (T_prev_micro_step, "T_prev_micro_step"),
                     (precip_liq_surf_mass, "precip_liq_surf_mass"),
                     (precip_ice_surf_mass, "precip_ice_surf_mass"),
                     # Computed fields
                     (eff_radius_qc, "eff_radius_qc"),
                     (eff_radius_qi, "eff_radius_qi"),
                     (eff_radius_qr, "eff_radius_qr"),
                     (precip_total_tend, "precip_total_tend"),
                     (nevapr, "nevapr"),
                     (diag_equiv_reflectivity, "diag_equiv_reflectivity"),
                     (micro_liq_ice_exchange, "micro_liq_ice_exchange"),
                     (micro_vap_liq_exchange, "micro_vap_liq_exchange"),
                     (micro_vap_ice_exchange, "micro_vap_ice_exchange"),
                     (rainfrac, "rainfrac")):
        buf[...] = np.asarray(out[key])
