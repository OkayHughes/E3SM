#!/usr/bin/env python3
"""Build the ZM dilute-CAPE core f2py extension (run in container)."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
ZM = CAM / "zm"
HERE = Path(__file__).resolve().parent

build_extension(
    name="eam_zm_cape_f",
    sources=SHR_SOURCES + [HERE / "stubs/shr_sys_stubs.F90",
                           CAM / "wv_sat_methods.F90",
                           CAM / "wv_saturation.F90",
                           ZM / "zm_conv_types.F90",
                           ZM / "zm_conv_util.F90",
                           ZM / "zm_conv_cape.F90"],
    driver=HERE / "drivers/zm_cape_driver.F90",
    includes=[CAM],
)
