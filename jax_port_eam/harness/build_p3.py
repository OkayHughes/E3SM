#!/usr/bin/env python3
"""Build the EAM P3 stratiform-microphysics f2py extension (run in the
scream-dev container). Compiles the REAL micro_p3.F90 (eam variant) with
its real dependencies micro_p3_utils.F90, wv_sat_scream.F90 (needs
cam/bfb_math.inc; without SCREAM_CONFIG_IS_CMAKE the bfb_* macros are
the plain gfortran intrinsics, matching the in-file bfb_* wrappers of
micro_p3.F90), physics_utils.F90 and scream_abortutils.F90 (non-SPMD
branch: plain abort()). New infrastructure-only stubs in
stubs/p3_stubs.F90: phys_control (use_hetfrz_classnuc flag + setter)
and debug_info (print-only report_error_info; referenced only by the
never-called polysvp1 error path). See drivers/p3_driver.F90.

-ffp-contract=off as for the other goldens: FMA contraction is not
reproducible in XLA/numpy and P3 is full of knife-edge thresholds.
"""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
P3 = REPO / "components/eam/src/physics/p3/eam"
HERE = Path(__file__).resolve().parent

build_extension(
    name="eam_p3_f",
    sources=SHR_SOURCES + [CAM / "physics_utils.F90",
                           HERE / "stubs/p3_stubs.F90",
                           P3 / "micro_p3_utils.F90",
                           CAM / "scream_abortutils.F90",
                           CAM / "wv_sat_scream.F90",
                           P3 / "micro_p3.F90"],
    driver=HERE / "drivers/p3_driver.F90",
    includes=[CAM],
    fflags=["-ffp-contract=off"],
)
