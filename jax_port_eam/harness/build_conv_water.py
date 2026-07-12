#!/usr/bin/env python3
"""Build the conv_water_4rad f2py extension (run in container).
conv_water.F90 use-associates zm_param from the real zm_conv, so the
whole zm_conv stack from build_zm_conv.py is compiled too; on top of
that come grid_stubs (ppgrid), the constituents stub, and
conv_water_stubs (perf_mod / shr_infnan / phys_control / cam_history /
physics_types / physics_buffer). conv_water.F90 itself is compiled
unmodified. See drivers/conv_water_driver.F90."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
ZM = CAM / "zm"
HERE = Path(__file__).resolve().parent

build_extension(
    name="eam_conv_water_f",
    sources=SHR_SOURCES + [HERE / "stubs/shr_sys_stubs.F90",
                           HERE / "stubs/cloud_fraction_stub.F90",
                           HERE / "stubs/grid_stubs.F90",
                           HERE / "stubs/constituents_stub.F90",
                           HERE / "stubs/conv_water_stubs.F90",
                           CAM / "wv_sat_methods.F90",
                           CAM / "wv_saturation.F90",
                           ZM / "zm_conv_types.F90",
                           ZM / "zm_conv_util.F90",
                           ZM / "zm_conv_cape.F90",
                           ZM / "zm_aero_type.F90",
                           ZM / "zm_microphysics_state.F90",
                           HERE / "stubs/zm_microphysics_stub.F90",
                           ZM / "zm_conv.F90",
                           CAM / "conv_water.F90"],
    driver=HERE / "drivers/conv_water_driver.F90",
    includes=[CAM],
)
