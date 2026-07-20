"""P3 microphysics (p3 process).

Transcribed kernel-by-kernel from
components/eamxx/src/physics/p3/impl/*_impl.hpp (master header
p3_functions.hpp). Conventions as in scream_jax.shoc: level axis last,
k=0 model top, leading axes batch columns.
"""

from . import tables

# Approximation-by-identity smoothing (scream_jax.foundation.smoothing):
# p3_process_step / p3_main / p3_main_part2 take static kwargs
# smooth_width (default 0.0 = exact, bitwise-original) and
# smooth_families (None = all families when smooth_width > 0, else a
# tuple drawn from SMOOTH_FAMILIES enabling site groups selectively for
# ablation). Site inventory (switching variable s, per-site scale, and
# jump-vs-kink justification) lives in comments at each adoption site:
#   "homog"  homogeneous_freezing T_homogfrz gate     (sedimentation.py)
#   "tmelt"  Tmelt/T_zerodegc gates in the ice kernels (processes_ice.py)
#   "frz"    T_rainfrz immersion-freezing gates       (processes_warm.py)
#   "evap"   evaporate_rain cloud-presence gate       (processes_ice.py)
#   "rime"   calc_bulk_rho_rime BSMALL snap           (main_part3.py)
#   "nucl"   ice_nucleation activation gate           (processes_warm.py)
SMOOTH_FAMILIES = ("homog", "tmelt", "frz", "evap", "rime", "nucl")


def family_width(smooth_width, smooth_families, family):
    """Static (trace-time) per-family smoothing width: `smooth_width` if
    `family` is enabled, else 0.0 (exact hard path). smooth_families=None
    enables every family."""
    if smooth_width == 0.0:
        return 0.0
    if smooth_families is None or family in smooth_families:
        return smooth_width
    return 0.0


# Scalar runtime options with the C++ defaults (P3Runtime in
# p3_functions.hpp). The boolean P3Runtime members are separate static
# arguments of the ported functions, not opts entries.
# NOTE on autodiff: p3_main / p3_process_step take a static kwarg
# sed_use_while_loop (default True = fast lax.while_loop primal, NOT
# reverse-mode differentiable). Pass sed_use_while_loop=False for
# autodiff (bounded masked lax.scan, bit-identical primal, ~5.5x
# slower; see p3/sedimentation.py docstring). The newer static kwarg
# sed_mode ("while" | "scan" | "uniform") supersedes it when given:
# "uniform" is a SMOOTH fixed-substep sedimentation surrogate (not
# bit-identical) that removes the adaptive CFL controller's
# discreteness. Static kwarg p3_soft_masks=True additionally disables
# the hard active/run3 column gates in p3_main (compute everywhere,
# O(qsmall) contributions from inactive columns).
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
