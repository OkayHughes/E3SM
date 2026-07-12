#!/usr/bin/env python3
"""Build the ZM precipitation-evaporation f2py extension (run in
container). Same zm_conv stack as build_zm_conv.py, except the
cloud_fraction STUB is replaced by the REAL cloud_fraction.F90
(zm_conv_evap calls its cldfrc_fice -- real physics), whose
infrastructure dependencies come from zm_intr_stubs.F90 (ppgrid with
pcols=48, phys_grid, dycore, ref_pres, phys_control, cam_history,
physics_buffer -- compile-time only; nothing but cldfrc_fice ever
executes). See drivers/zm_evap_driver.F90."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
ZM = CAM / "zm"
HERE = Path(__file__).resolve().parent

build_extension(
    name="eam_zm_evap_f",
    sources=SHR_SOURCES + [HERE / "stubs/shr_sys_stubs.F90",
                           HERE / "stubs/zm_intr_stubs.F90",
                           CAM / "wv_sat_methods.F90",
                           CAM / "wv_saturation.F90",
                           CAM / "cloud_fraction.F90",
                           ZM / "zm_conv_types.F90",
                           ZM / "zm_conv_util.F90",
                           ZM / "zm_conv_cape.F90",
                           ZM / "zm_aero_type.F90",
                           ZM / "zm_microphysics_state.F90",
                           HERE / "stubs/zm_microphysics_stub.F90",
                           ZM / "zm_conv.F90"],
    driver=HERE / "drivers/zm_evap_driver.F90",
    includes=[CAM],
)
