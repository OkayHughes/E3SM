"""JAX port of the EAM P3 stratiform microphysics
(components/eam/src/physics/p3/eam/micro_p3.F90, EAMv3 microp_scheme='P3').

Adapted from the validated SCREAM port ../../jax_port/scream_jax/p3
(same scheme family; the C++ descends from this Fortran). Modules were
COPIED and then adapted line-by-line against the EAM Fortran, which is
the truth for this port. PORT_NOTES (differences from the scream_jax
starting point, all traced to the eam-vs-scream Fortran diff):

  * naming: EAM's input `exner` MULTIPLIES T to give theta
    (th = T*exner, exner = (p0/p)^(Rd/cp)); scream calls the same
    quantity `inv_exner`. This port uses the EAM names throughout.
  * saturation: qv_sat = ep_2*es/max(1e-3, p - es) (MurphyKoop es;
    "wet" denominator p - es) everywhere; scream uses qv_sat_dry
    (denominator p).
  * tables: mu_r_constant = 0 (scream: 1), so the rain fallspeed /
    ventilation tables differ numerically; the p3_init_b integration
    weight for the ventilation term uses a SINGLE-precision 1.e-6
    literal (tables.py reproduces it).
  * tuning parameters are namelist inputs (opts): p3_autocon_coeff,
    p3_qc/nc_autocon_expon, p3_accret_coeff, p3_qc_accret_expon,
    p3_wbf_coeff, p3_mincdnc, p3_max_mean_rain_size,
    p3_embryonic_rain_size, nccnst; subgrid_variance_scaling is ACTIVE
    in autoconversion/accretion/immersion freezing (scream hard-set 1).
  * get_rain_dsd2: lammin from p3_max_mean_rain_size, nr recomputed via
    exp(3log lamr + log qr + log G(mu+1) - log G(mu+4))/cons1, logn0r
    from nr (not cdistr); cdistr/logn0r computed inside the routine.
  * ice_deposition_sublimation: EAM (older) form — qidep = epsi/abi*
    (qv - qv_sat_i) with NO min(.., inv_dt) limiter, sublimation at any
    T, Bergeron * p3_wbf_coeff.
  * conservation sequence: ...ni_conservation ->
    prevent_ice_overdepletion (EAM-only, replaces scream's
    prevent_liq_supersaturation) -> ice_supersat_conservation.
  * wet growth: EAM does NOT clamp qccol/qrcol at 0 after shedding
    (scream added max(0,..)).
  * hetfrz classnuc (CNT_couple) qc threshold is qsmall = 1e-14
    (scream ice_classical_nucleation uses 1e-18).
  * ice_nucleation: EAM Cooper branch requires N_nuc >= 1e-20 and has
    the do_Cooper_inP3 add-on branch under predict_nc.
  * update_prognostic_ice additionally returns qi_wetDepos; EAM's
    precip_total_tend = qcacc + qcaut + qcshd + qi_wetDepos (scream
    uses qccol as the last term).
  * post-sedimentation: homogeneous_freezing tests th/exner (updated
    theta, sequentially qc then qr); NEW ice_complete_melting at
    T > 273.15(f32)+2; p3_mincdnc floor on nc before/within part3.
  * sedimentation accumulates time-integrated fluxes cflx/rflx/sflx and
    precip_liq_flux as sum(flux*dt_sub)*inv_dt; precip_ice_flux stays 0
    (never accumulated in the Fortran); cloud sed ADDS to
    precip_liq_surf.
  * part3: diag_ze_rain/diag_ze_ice outputs, eff-radius/vm zeroing in
    the no-condensate branches, no diag_eff_radius_qr, mincdnc floor on
    nc_incld after the cloud dsd call.
  * part1: prescribed-CCN branch is nc = max(nc, nccn_prescribed)
    (no cloud-fraction scaling, no SPA factor/exponent); nccnst is a
    runtime input.
  * p3_tend_out (49 slots) is reconstructed for output parity.
  * constants from micro_p3_utils.F90/physconst verbatim, including
    max_total_ni = 500e3 (scream: 740e3) and the single-precision
    literal promotions incloud_limit = f32(5.1e-3), precip_limit =
    f32(1e-2).

float64 via the package-level jax config in eam_jax/__init__.py.
"""

from . import constants as c  # noqa: F401

# EAMv3 phys="default" namelist values
# (bld/namelist_files/namelist_defaults_eam.xml; recorded in the golden
# metadata). Booleans (do_predict_nc, do_prescribed_CCN, do_precip_off,
# use_hetfrz_classnuc, do_Cooper_inP3) are static arguments.
DEFAULT_OPTS = {
    "p3_autocon_coeff": 30500.0,
    "p3_qc_autocon_expon": 3.19,
    "p3_nc_autocon_expon": -1.10,
    "p3_accret_coeff": 117.25,
    "p3_qc_accret_expon": 1.15,
    "p3_wbf_coeff": 1.0,
    "p3_mincdnc": 20.0e6,
    "p3_max_mean_rain_size": 0.005,
    "p3_embryonic_rain_size": 25.0e-6,
    "nccnst": 200.0e6,
}
