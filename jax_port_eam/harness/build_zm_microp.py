#!/usr/bin/env python3
"""Build the ZM convective-microphysics f2py extension (run in
container). Same zm_conv stack as build_zm_evap.py (real
cloud_fraction.F90 + zm_intr_stubs for its infrastructure), except the
abort-only zm_microphysics stub is replaced by the REAL
zm_microphysics.F90 together with its REAL aerosol-activation
dependencies (activate_drop_mam.F90, nucleate_ice_conv.F90) and the
real share/util/shr_spfn_mod.F90. New infrastructure-only stubs in
stubs/zm_microp_stubs.F90: time_manager (get_step_size = driver-set
step size) and an abort-only ndrop_bam (bulk activation; never
executed -- the harness, like EAMv3 production, uses the modal
scheme). See drivers/zm_microp_driver.F90."""
from pathlib import Path

from fbuild import REPO, SHR_SOURCES, build_extension

CAM = REPO / "components/eam/src/physics/cam"
ZM = CAM / "zm"
HERE = Path(__file__).resolve().parent

build_extension(
    name="eam_zm_microp_f",
    sources=SHR_SOURCES + [HERE / "stubs/shr_sys_stubs.F90",
                           HERE / "stubs/zm_intr_stubs.F90",
                           HERE / "stubs/zm_microp_stubs.F90",
                           REPO / "share/util/shr_spfn_mod.F90",
                           CAM / "wv_sat_methods.F90",
                           CAM / "wv_saturation.F90",
                           CAM / "cloud_fraction.F90",
                           CAM / "activate_drop_mam.F90",
                           CAM / "nucleate_ice_conv.F90",
                           ZM / "zm_conv_types.F90",
                           ZM / "zm_conv_util.F90",
                           ZM / "zm_conv_cape.F90",
                           ZM / "zm_aero_type.F90",
                           ZM / "zm_microphysics_state.F90",
                           ZM / "zm_microphysics.F90",
                           ZM / "zm_conv.F90"],
    driver=HERE / "drivers/zm_microp_driver.F90",
    includes=[CAM],
    # gfortran defaults to -ffp-contract=fast; FMA contraction of the
    # hygroscopicity accumulation (hygro + vol*spechygro) shifts the
    # result by 1 ulp across the razor-edge `hygro > 1e-10` activation
    # threshold (MAM pom/bc hygro is exactly 1e-10). XLA/numpy do not
    # fuse, so the golden reference is built without contraction.
    fflags=["-ffp-contract=off"],
)
