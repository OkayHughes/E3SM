#!/usr/bin/env python3
"""Build the geopotential f2py extension (container)."""
from pathlib import Path

from fbuild import HERE, REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"

build_extension(
    name="eam_geopotential_f",
    sources=SHR_SOURCES + [HERE / "stubs/grid_stubs.F90",
                           CAM / "geopotential.F90"],
    driver=Path(__file__).parent / "drivers/geopotential_driver.F90",
    includes=[CAM],
)
