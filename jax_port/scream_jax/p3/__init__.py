"""P3 microphysics (p3 process).

Transcribed kernel-by-kernel from
components/eamxx/src/physics/p3/impl/*_impl.hpp (master header
p3_functions.hpp). Conventions as in scream_jax.shoc: level axis last,
k=0 model top, leading axes batch columns.
"""

from . import tables

# Scalar runtime options with the C++ defaults (P3Runtime in
# p3_functions.hpp). The boolean P3Runtime members are separate static
# arguments of the ported functions, not opts entries.
# NOTE on autodiff: p3_main / p3_process_step take a static kwarg
# sed_use_while_loop (default True = fast lax.while_loop primal, NOT
# reverse-mode differentiable). Pass sed_use_while_loop=False for
# autodiff (bounded masked lax.scan, bit-identical primal, ~5.5x
# slower; see p3/sedimentation.py docstring).
DEFAULT_OPTS = {
    "max_total_ni": 740.0e3,
    "autoconversion_prefactor": 1350.0,
    "autoconversion_qc_exponent": 2.47,
    "autoconversion_nc_exponent": 1.79,
    "autoconversion_radius": 25.0e-6,
    "accretion_prefactor": 67.0,
    "accretion_qc_exponent": 1.15,
    "accretion_qr_exponent": 1.15,
    "rain_selfcollection_prefactor": 5.78,
    "rain_selfcollection_breakup_diameter": 0.00028,
    "constant_mu_rain": 1.0,
    "spa_ccn_to_nc_factor": 1.0,
    "spa_ccn_to_nc_exponent": 1.0,
    "cldliq_to_ice_collection_factor": 0.5,
    "rain_to_ice_collection_factor": 1.0,
    "min_rime_rho": 50.0,
    "max_rime_rho": 900.0,
    "immersion_freezing_exponent": 0.65,
    "deposition_nucleation_exponent": 0.304,
    "ice_sedimentation_factor": 1.0,
}
