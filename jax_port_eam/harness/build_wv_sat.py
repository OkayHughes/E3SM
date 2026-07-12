#!/usr/bin/env python3
"""Build the wv_sat f2py extension (run with jax_port_eam/.venv python)."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"

build_extension(
    name="eam_wv_sat_f",
    sources=SHR_SOURCES + [CAM / "wv_sat_methods.F90",
                           CAM / "wv_saturation.F90"],
    driver=Path(__file__).parent / "drivers/wv_sat_driver.F90",
    includes=[CAM],
)
