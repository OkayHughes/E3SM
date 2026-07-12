#!/usr/bin/env python3
"""Build the cldfrc2m cloud-fraction f2py extension (run in container)."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
HERE = Path(__file__).resolve().parent

build_extension(
    name="eam_cldfrc2m_f",
    sources=SHR_SOURCES + [HERE / "stubs/grid_stubs.F90",
                           HERE / "stubs/cloud_fraction_stub.F90",
                           CAM / "wv_sat_methods.F90",
                           CAM / "wv_saturation.F90",
                           CAM / "cldfrc2m.F90"],
    driver=HERE / "drivers/cldfrc2m_driver.F90",
    includes=[CAM],
)
