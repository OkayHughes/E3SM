#!/usr/bin/env python3
"""Build the ZM convective transport f2py extension (run in
container). Same source stack as build_zm_conv.py (zm_transport.F90
use-associates zm_param from the real zm_conv, which needs the whole
chain) plus the constituents stub (cnst_get_type_byind for the
dry-mixing-ratio branch) and zm_transport.F90 itself, compiled
unmodified. See drivers/zm_transport_driver.F90."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
ZM = CAM / "zm"
HERE = Path(__file__).resolve().parent

build_extension(
    name="eam_zm_transport_f",
    sources=SHR_SOURCES + [HERE / "stubs/shr_sys_stubs.F90",
                           HERE / "stubs/cloud_fraction_stub.F90",
                           HERE / "stubs/constituents_stub.F90",
                           CAM / "wv_sat_methods.F90",
                           CAM / "wv_saturation.F90",
                           ZM / "zm_conv_types.F90",
                           ZM / "zm_conv_util.F90",
                           ZM / "zm_conv_cape.F90",
                           ZM / "zm_aero_type.F90",
                           ZM / "zm_microphysics_state.F90",
                           HERE / "stubs/zm_microphysics_stub.F90",
                           ZM / "zm_conv.F90",
                           ZM / "zm_transport.F90"],
    driver=HERE / "drivers/zm_transport_driver.F90",
    includes=[CAM],
)
