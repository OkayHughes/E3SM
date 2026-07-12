#!/usr/bin/env python3
"""Build the tropopause-finder f2py extension (run in container).

Compiles the REAL tropopause.F90 and the REAL interpolate_data.F90
regridder; the climatology-file read is served by the pio stub in
tropopause_stubs.F90 (see driver header)."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
CTL = REPO / "components/eam/src/control"
HERE = Path(__file__).resolve().parent

build_extension(
    name="eam_tropopause_f",
    sources=SHR_SOURCES + [HERE / "stubs/grid_stubs.F90",
                           HERE / "stubs/tropopause_stubs.F90",
                           CTL / "interpolate_data.F90",
                           CAM / "tropopause.F90"],
    driver=HERE / "drivers/tropopause_driver.F90",
    includes=[CAM],
)
