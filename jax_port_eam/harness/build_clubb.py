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

Verbatim extraction (slice B): advance_clubb_core_module.F90's
PRIVATE subroutine pdf_closure_driver (and the private helpers it
calls: trapezoidal_rule_zt/zm, trapezoid_zt/zm, compute_cloud_cover,
clip_rcm) is re-wrapped as a PUBLIC module procedure so the harness
can golden it directly.  The extraction is generated at build time by
copying the subroutine text VERBATIM from the real source (see
gen_pdf_extract below) into module clubb_pdf_extract in the build
dir, compiled with the same defines -- verbatim by construction, and
additionally validated end-to-end against the real (public)
advance_clubb_core in gen_clubb_golden.py: under EAMv3's
ipdf_call_placement=2 the final pdf_closure_driver call's outputs
pass through advance_clubb_core untouched, and replaying the
extracted routine on the post-advance state reproduces them BITWISE.

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

def _extract_block(lines, start_marker, end_marker):
    """Return the verbatim text lines[start..end] (inclusive) where
    start is the first line beginning with start_marker and end is the
    first subsequent line beginning with end_marker."""
    i0 = None
    for i, line in enumerate(lines):
        if i0 is None and line.startswith(start_marker):
            i0 = i
        elif i0 is not None and line.startswith(end_marker):
            return "".join(lines[i0:i + 1])
    raise RuntimeError(f"markers not found: {start_marker!r} .. "
                       f"{end_marker!r}")


def gen_pdf_extract(dst):
    """Generate clubb_pdf_extract.F90: a module wrapping VERBATIM
    copies of advance_clubb_core_module.F90's private
    pdf_closure_driver + trapezoidal_rule_zt/zm + trapezoid_zt/zm +
    compute_cloud_cover + clip_rcm.  No line of subroutine text is
    edited (the GFDL ifdef blocks stay and compile out exactly as in
    the real build); only the module wrapper is new."""
    src = CLUBB / "advance_clubb_core_module.F90"
    lines = src.read_text().splitlines(keepends=True)
    body = _extract_block(lines, "  subroutine pdf_closure_driver(",
                          "  end subroutine pdf_closure_driver")
    helpers = _extract_block(lines, "    subroutine trapezoidal_rule_zt",
                             "    end subroutine clip_rcm")
    header = (
        "! AUTO-GENERATED by build_clubb.py -- DO NOT EDIT.\n"
        "! Verbatim extraction of the PRIVATE subroutine\n"
        "! pdf_closure_driver (+ its private helpers) from\n"
        f"! {src}\n"
        "! re-wrapped as a public module procedure for the f2py\n"
        "! harness.  See build_clubb.py docstring for the\n"
        "! verbatim-by-construction + end-to-end validation story.\n"
        "module clubb_pdf_extract\n"
        "  implicit none\n"
        "  public :: pdf_closure_driver\n"
        "  private\n"
        "  contains\n\n")
    footer = "\nend module clubb_pdf_extract\n"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(header + body + "\n" + helpers + footer)
    return dst


EXTRACT = gen_pdf_extract(
    Path(__file__).resolve().parent.parent
    / "build/eam_clubb_f/clubb_pdf_extract.F90")

build_extension(
    name="eam_clubb_f",
    sources=[SHARE / "shr_kind_mod.F90",
             SHARE / "shr_const_mod.F90",
             HERE / "stubs/infrastructure_stubs.F90",
             HERE / "stubs/clubb_stubs.F90"]
            + [CLUBB / s for s in CLUBB_SOURCES]
            + [EXTRACT, HERE / "drivers/clubb_driver_helpers.F90"],
    driver=HERE / "drivers/clubb_driver.F90",
    includes=[CLUBB],
    fflags=["-ffp-contract=off", "-DCLUBB_CAM", "-DCLUBB_SGS",
            "-DCLUBB_REAL_TYPE=dp"],
    link_flags=["-llapack"],
)
