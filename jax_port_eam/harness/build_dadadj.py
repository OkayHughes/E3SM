#!/usr/bin/env python3
"""Build the dadadj f2py extension (run in the scream-dev container)."""
from pathlib import Path

from fbuild import HERE, REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"

build_extension(
    name="eam_dadadj_f",
    sources=SHR_SOURCES + [HERE / "stubs/grid_stubs.F90",
                           CAM / "dadadj.F90"],
    driver=Path(__file__).parent / "drivers/dadadj_driver.F90",
    includes=[CAM],
)
