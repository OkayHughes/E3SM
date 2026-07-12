#!/usr/bin/env python3
"""Build the gw spectrum (non-orographic) f2py extension (container)."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
GW = CAM / "gw"
STUBS = Path(__file__).parent / "stubs"

build_extension(
    name="eam_gw_spectrum_f",
    sources=SHR_SOURCES + [STUBS / "ref_pres_stub.F90",
                           GW / "gw_utils.F90",
                           CAM / "vdiff_lu_solver.F90",
                           GW / "gw_diffusion.F90",
                           GW / "gw_common.F90",
                           GW / "gw_front.F90",
                           GW / "gw_convect.F90"],
    driver=Path(__file__).parent / "drivers/gw_spectrum_driver.F90",
    includes=[CAM],
)
