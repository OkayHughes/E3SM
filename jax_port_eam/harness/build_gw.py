#!/usr/bin/env python3
"""Build the gw (orographic spine) f2py extension (container)."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
GW = CAM / "gw"

build_extension(
    name="eam_gw_f",
    sources=SHR_SOURCES + [GW / "gw_utils.F90",
                           CAM / "vdiff_lu_solver.F90",
                           GW / "gw_diffusion.F90",
                           GW / "gw_common.F90", GW / "gw_oro.F90"],
    driver=Path(__file__).parent / "drivers/gw_driver.F90",
    includes=[CAM],
)
