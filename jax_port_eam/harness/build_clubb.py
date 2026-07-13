#!/usr/bin/env python3
"""Build the EAM CLUBB f2py extension (run in the scream-dev container).

Compiles the ENTIRE unmodified CLUBB stack
(components/eam/src/physics/clubb/*.F90 + mt95.f90) in dependency
order with the exact defines EAM's bld/configure uses for CLUBB:
  -DCLUBB_CAM -DCLUBB_SGS -DCLUBB_REAL_TYPE=dp
(core_rknd = double precision; constants_clubb takes its physical
constants from the REAL shr_const_mod under CLUBB_CAM).

NOT defined (so those paths compile out, exactly as in an EAM build
without them): NETCDF (output_netcdf is an empty shell; EAM compiles
CLUBB the same way -- stats output goes through cam_history, and
l_stats=.false. in this harness anyway), GFDL, MKL, CLUBBND_CAM, SPMD.

Infrastructure stubs (never physics):
  - stubs/infrastructure_stubs.F90: cam_abortutils, spmd_utils
    (masterproc=.false.), mpishorthand, namelist_utils, units
    (only parameters_tunable's clubb_param_readnl / read_parameters
    reference them; with masterproc=.false. no namelist file is read
    and every clubb_* override variable is initialized to init_value,
    which is exactly the state EAM has on non-master ranks before the
    broadcast -- the harness then applies the EAMv3 overrides to the
    params vector on the Python side, see gen_clubb_golden.py).
  - stubs/clubb_stubs.F90: phys_control (use_od_fd=.false., the EAMv3
    default; only advance_windm_edsclrm references it).

Links against the container's reference LAPACK (liblapack) for
CLUBB's lapack_wrap.F90 (dgtsv/dgbsv/dgbsvx and the *svx condition
estimators).

-ffp-contract=off as for the other goldens (FMA contraction is not
reproducible in XLA/numpy).
"""
from pathlib import Path

from fbuild import REPO, build_extension

CLUBB = REPO / "components/eam/src/physics/clubb"
SHARE = REPO / "share/util"
HERE = Path(__file__).resolve().parent

# Dependency-ordered (topologically sorted use-graph; computed from the
# actual `use` statements, see PORTING_PLAN.md row 12).
CLUBB_SOURCES = [
    "clubb_precision.F90",
    "constants_clubb.F90",
    "error_code.F90",
    "file_functions.F90",
    "model_flags.F90",
    "interpolation.F90",
    "grid_class.F90",
    "LY93_pdf.F90",
    "Nc_Ncn_eqns.F90",
    "parameter_indices.F90",
    "parameters_tunable.F90",
    "Skx_module.F90",
    "T_in_K_module.F90",
    "calc_roots.F90",
    "parameters_model.F90",
    "adg1_adg2_3d_luhar_pdf.F90",
    "saturation.F90",
    "stat_file_module.F90",
    "stats_type.F90",
    "stats_type_utilities.F90",
    "stats_variables.F90",
    "advance_helper_module.F90",
    "clip_explicit.F90",
    "diffusion.F90",
    "lapack_wrap.F90",
    "mean_adv.F90",
    "sponge_layer_damping.F90",
    "advance_windm_edsclrm_module.F90",
    "csr_matrix_module.F90",
    "array_index.F90",
    "index_mapping.F90",
    "fill_holes.F90",
    "gmres_cache.F90",
    "gmres_wrap.F90",
    "new_pdf.F90",
    "pdf_parameter_module.F90",
    "new_pdf_main.F90",
    "new_tsdadg_pdf.F90",
    "variables_diagnostic_module.F90",
    "variables_prognostic_module.F90",
    "numerical_check.F90",
    "pdf_utilities.F90",
    "pdf_closure_module.F90",
    "advance_wp2_wp3_module.F90",
    "clip_semi_implicit.F90",
    "mono_flux_limiter.F90",
    "pos_definite_module.F90",
    "turbulent_adv_pdf.F90",
    "advance_xm_wpxp_module.F90",
    "advance_xp2_xpyp_module.F90",
    "advance_xp3_module.F90",
    "mixing_length.F90",
    "sigma_sqd_w_module.F90",
    "calendar.F90",
    "endian.F90",
    "output_grads.F90",
    "output_netcdf.F90",
    "stats_lh_sfc_module.F90",
    "stats_lh_zt_module.F90",
    "stats_rad_zm_module.F90",
    "stats_rad_zt_module.F90",
    "stats_sfc_module.F90",
    "stats_zm_module.F90",
    "stats_zt_module.F90",
    "stats_clubb_utilities.F90",
    "surface_varnce_module.F90",
    "advance_clubb_core_module.F90",
    "input_names.F90",
    "input_reader.F90",
    "matrix_operations.F90",
    "corr_varnce_module.F90",
    "hydromet_pdf_parameter_module.F90",
    "mt95.f90",
    "diagnose_correlations_module.F90",
    "precipitation_fraction.F90",
    "setup_clubb_pdf_params.F90",
    "clubb_api_module.F90",
    "code_timer_module.F90",
]

build_extension(
    name="eam_clubb_f",
    sources=[SHARE / "shr_kind_mod.F90",
             SHARE / "shr_const_mod.F90",
             HERE / "stubs/infrastructure_stubs.F90",
             HERE / "stubs/clubb_stubs.F90"]
            + [CLUBB / s for s in CLUBB_SOURCES],
    driver=HERE / "drivers/clubb_driver.F90",
    includes=[CLUBB],
    fflags=["-ffp-contract=off", "-DCLUBB_CAM", "-DCLUBB_SGS",
            "-DCLUBB_REAL_TYPE=dp"],
    link_flags=["-llapack"],
)
