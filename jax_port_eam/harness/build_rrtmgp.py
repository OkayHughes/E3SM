#!/usr/bin/env python3
"""Build the EAM RRTMGP radiation f2py extension (run in the
scream-dev container; needs netcdf-fortran, present in the container).

Compiles the REAL EAM RRTMGP stack exactly as bld/configure assembles
it for rad=rrtmgp without -rrtmgpxx (the EAMv3 default):
  external/rte + rte/kernels, external/rrtmgp + rrtmgp/kernels,
  external/extensions, external/examples (mo_simple_netcdf),
  f90/{rrtmgp_interface,mo_load_coefficients},
  the EAM driver layer (radconstants, assertions, radiation_state,
  radiation_utils, cloud_rad_props, ebert_curry, slingo,
  mcica_subcol_gen, cam_optics), control/interpolate_data, and the
  REAL share/RandNum KISS generator (C kernel included).

Infrastructure-only stubs in stubs/rrtmgp_stubs.F90 (own copies of
cam_logfile/cam_abortutils/spmd_utils because this harness needs
masterproc=.true. for the netCDF reads; see the stub file header).

-ffp-contract=off as for the other goldens. -DDSFMT_MEXP=19937 is the
standard dSFMT config (needed to compile shr_RandNum's dSFMT backend,
which the harness links but never calls — MCICA uses KISS).
"""
from pathlib import Path

from fbuild import REPO, build_extension

HERE = Path(__file__).resolve().parent
RAD = REPO / "components/eam/src/physics/rrtmgp"
EXT = RAD / "external"
CTL = REPO / "components/eam/src/control"
RAND = REPO / "share/RandNum/src"

SOURCES = [
    REPO / "share/util/shr_kind_mod.F90",
    REPO / "share/util/shr_const_mod.F90",
    HERE / "stubs/grid_stubs.F90",
    HERE / "stubs/rrtmgp_stubs.F90",
    HERE / "stubs/physconst_stub.F90",
    HERE / "stubs/constituents_stub.F90",
    # share RandNum (real RNG; C kernels)
    RAND / "kissvec/kissvec.c",
    RAND / "kissvec/kissvec_mod.F90",
    RAND / "dsfmt_f03/dSFMT.c",
    RAND / "dsfmt_f03/dSFMT_utils.c",
    RAND / "dsfmt_f03/dSFMT_interface.F90",
    RAND / "mt19937/mersennetwister_mod.F90",
    RAND / "shr_RandNum_mod.F90",
    # RRTMGP external (dependency order)
    EXT / "rte/mo_rte_kind.F90",
    EXT / "rrtmgp/mo_rrtmgp_constants.F90",
    EXT / "rrtmgp/mo_rrtmgp_util_string.F90",
    EXT / "rte/mo_rte_util_array.F90",
    EXT / "rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90",
    EXT / "rrtmgp/mo_rrtmgp_util_reorder.F90",
    EXT / "rte/kernels/mo_optical_props_kernels.F90",
    EXT / "rte/mo_optical_props.F90",
    EXT / "rte/mo_source_functions.F90",
    EXT / "rrtmgp/mo_gas_concentrations.F90",
    EXT / "rrtmgp/mo_gas_optics.F90",
    EXT / "rrtmgp/kernels/mo_gas_optics_kernels.F90",
    EXT / "rrtmgp/mo_gas_optics_rrtmgp.F90",
    EXT / "rte/kernels/mo_fluxes_broadband_kernels.F90",
    EXT / "rte/mo_fluxes.F90",
    EXT / "extensions/mo_fluxes_byband_kernels.F90",
    EXT / "extensions/mo_fluxes_byband.F90",
    EXT / "rte/kernels/mo_rte_solver_kernels.F90",
    EXT / "rte/mo_rte_lw.F90",
    EXT / "rte/mo_rte_sw.F90",
    EXT / "extensions/mo_rrtmgp_clr_all_sky.F90",
    EXT / "examples/mo_simple_netcdf.F90",
    RAD / "f90/mo_load_coefficients.F90",
    # EAM driver layer
    RAD / "assertions.F90",
    RAD / "radconstants.F90",
    CTL / "interpolate_data.F90",
    RAD / "cloud_rad_props.F90",
    RAD / "ebert_curry.F90",
    RAD / "slingo.F90",
    RAD / "mcica_subcol_gen.F90",
    RAD / "radiation_state.F90",
    RAD / "radiation_utils.F90",
    RAD / "cam_optics.F90",
    RAD / "f90/rrtmgp_interface.F90",
    # harness core (derived-type storage + transcribed private helpers)
    HERE / "drivers/rrtmgp_core.F90",
]

build_extension(
    name="eam_rrtmgp_f",
    sources=SOURCES,
    driver=HERE / "drivers/rrtmgp_driver.F90",
    includes=[REPO / "share/RandNum/include", "/usr/include"],
    fflags=["-ffp-contract=off", "-DDSFMT_MEXP=19937"],
    link_flags=["-lnetcdff", "-lnetcdf"],
)
