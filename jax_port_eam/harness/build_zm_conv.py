#!/usr/bin/env python3
"""Build the ZM deep-convection main-routine f2py extension (run in
container). Extends the zm_cape recipe with zm_aero_type +
zm_microphysics_state (real, state containers only) and the
zm_microphysics/cldfrc_fice abort-only interface stubs; zm_conv.F90
itself is compiled unmodified. See drivers/zm_conv_driver.F90 for the
zm_microp=.false. scope note."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
ZM = CAM / "zm"
HERE = Path(__file__).resolve().parent

build_extension(
    name="eam_zm_conv_f",
    sources=SHR_SOURCES + [HERE / "stubs/shr_sys_stubs.F90",
                           HERE / "stubs/cloud_fraction_stub.F90",
                           CAM / "wv_sat_methods.F90",
                           CAM / "wv_saturation.F90",
                           ZM / "zm_conv_types.F90",
                           ZM / "zm_conv_util.F90",
                           ZM / "zm_conv_cape.F90",
                           ZM / "zm_aero_type.F90",
                           ZM / "zm_microphysics_state.F90",
                           HERE / "stubs/zm_microphysics_stub.F90",
                           ZM / "zm_conv.F90"],
    driver=HERE / "drivers/zm_conv_driver.F90",
    includes=[CAM],
)
