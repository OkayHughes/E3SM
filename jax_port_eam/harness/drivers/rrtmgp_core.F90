! f2py driver for the EAM RRTMGP radiation layer.
!
! Compiled against the REAL, UNMODIFIED EAM sources:
!   components/eam/src/physics/rrtmgp/external/{rte,rrtmgp,extensions}
!   (Fortran RRTMGP as EAM builds it: bld/configure non-rrtmgpxx path),
!   f90/{rrtmgp_interface,mo_load_coefficients}.F90,
!   {radconstants,assertions,radiation_state,radiation_utils,cam_optics,
!    cloud_rad_props,ebert_curry,slingo,mcica_subcol_gen}.F90,
!   control/interpolate_data.F90, share/RandNum (real KISS RNG).
! Infrastructure-only stubs in stubs/rrtmgp_stubs.F90.
!
! Subroutines transcribed VERBATIM from radiation.F90 PRIVATE
! procedures (they cannot be compiled standalone because radiation.F90
! drags in pbuf/history/COSP): reordered(), set_albedo(),
! radiation_driver_sw/lw bodies (day gather + extra-level padding +
! heating rates), set_net_fluxes_sw/lw, export_surface_fluxes,
! set_daynight_indices. Transcriptions are validated end-to-end by the
! full drv_rad_step golden.
module rrtmgp_core
  use mo_rte_kind, only: wp
  use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)

  type(ty_gas_optics_rrtmgp), save :: kdist_sw_drv, kdist_lw_drv

  ! driver-owned k-distributions for kernel-level goldens (the real
  ! rrtmgp_interface holds its own private pair, initialized in
  ! drv_init through the same load_and_init path)
  logical :: kdist_loaded = .false.

contains

  subroutine stop_on_err(msg)
    character(len=*), intent(in) :: msg
    if (trim(msg) /= '') then
      print *, 'RRTMGP DRIVER ERROR: ', trim(msg)
      stop 1
    end if
  end subroutine stop_on_err

  subroutine active_gas_list(gases)
    ! radiation.F90 active_gases, verbatim
    character(len=3), dimension(8), intent(out) :: gases
    gases = (/ 'H2O', 'CO2', 'O3 ', 'N2O', 'CO ', 'CH4', 'O2 ', 'N2 ' /)
  end subroutine active_gas_list

  subroutine drv_init(sw_file, lw_file, liq_file, ice_file)
    use rrtmgp_interface, only: rrtmgp_initialize
    use rad_constituents, only: iceopticsfile, liqopticsfile, &
                                icecldoptics, liqcldoptics
    use cloud_rad_props,  only: cloud_rad_props_init
    use radiation_state,  only: ktop, kbot, nlev_rad
    use ppgrid,           only: pver
    character(len=256), intent(in) :: sw_file, lw_file, liq_file, ice_file
    character(len=3), dimension(8) :: gases

    call active_gas_list(gases)

    ! EAMv3 defaults (namelist_defaults_eam.xml rad="rrtmgp")
    icecldoptics = 'mitchell'
    liqcldoptics = 'gammadist'
    liqopticsfile = trim(liq_file)
    iceopticsfile = trim(ice_file)
    call cloud_rad_props_init()

    call rrtmgp_initialize(size(gases), gases, trim(sw_file), trim(lw_file))

    ! radiation_init: nlev_rad = pver + 1 (no NO_EXTRA_RAD_LEVEL)
    nlev_rad = pver + 1
    ktop = nlev_rad - pver + 1
    kbot = nlev_rad
    kdist_loaded = .true.
  end subroutine drv_init

  subroutine drv_dims(nswb, nlwb, nswg, nlwg)
    use rrtmgp_interface, only: nswbands, nlwbands, nswgpts, nlwgpts
    integer, intent(out) :: nswb, nlwb, nswg, nlwg
    nswb = nswbands; nlwb = nlwbands; nswg = nswgpts; nlwg = nlwgpts
  end subroutine drv_dims

  subroutine drv_gpt_bands(nswg, nlwg, gpb_sw, gpb_lw)
    use rrtmgp_interface, only: get_gpoint_bands_sw, get_gpoint_bands_lw
    integer, intent(in) :: nswg, nlwg
    integer, intent(out) :: gpb_sw(nswg), gpb_lw(nlwg)
    call get_gpoint_bands_sw(gpb_sw)
    call get_gpoint_bands_lw(gpb_lw)
  end subroutine drv_gpt_bands

  subroutine drv_temp_limits(tmin, tmax)
    use rrtmgp_interface, only: get_min_temperature, get_max_temperature
    real(r8), intent(out) :: tmin, tmax
    tmin = get_min_temperature()
    tmax = get_max_temperature()
  end subroutine drv_temp_limits

  ! ------------------------------------------------------------------
  ! gas optics kernels (driver-owned k-dists via the same
  ! load_and_init path the interface uses)
  ! ------------------------------------------------------------------
  subroutine drv_gas_optics_sw(ncol, nlev, ngas, ngpt, gas_vmr, &
       pmid, tmid, pint, tau, ssa, g, toa_flux)
    use mo_gas_concentrations, only: ty_gas_concs
    use mo_optical_props,      only: ty_optical_props_2str
    integer, intent(in) :: ncol, nlev, ngas, ngpt
    real(r8), intent(in) :: gas_vmr(ngas, ncol, nlev)
    real(r8), intent(in) :: pmid(ncol, nlev), tmid(ncol, nlev)
    real(r8), intent(in) :: pint(ncol, nlev + 1)
    real(r8), intent(out) :: tau(ncol, nlev, ngpt), ssa(ncol, nlev, ngpt)
    real(r8), intent(out) :: g(ncol, nlev, ngpt), toa_flux(ncol, ngpt)

    type(ty_gas_concs) :: concs
    type(ty_optical_props_2str) :: optics

    call set_concs(ngas, ncol, nlev, gas_vmr, concs)
    call stop_on_err(optics%alloc_2str(ncol, nlev, kdist_sw_drv))
    call stop_on_err(kdist_sw_drv%gas_optics(pmid, pint, tmid, concs, &
                                        optics, toa_flux))
    tau = optics%tau
    ssa = optics%ssa
    g = optics%g
  end subroutine drv_gas_optics_sw

  subroutine drv_gas_optics_lw(ncol, nlev, ngas, ngpt, gas_vmr, &
       pmid, tmid, pint, tint, tau, lay_src, lev_src_inc, lev_src_dec, &
       sfc_src)
    use mo_gas_concentrations, only: ty_gas_concs
    use mo_optical_props,      only: ty_optical_props_1scl
    use mo_source_functions,   only: ty_source_func_lw
    integer, intent(in) :: ncol, nlev, ngas, ngpt
    real(r8), intent(in) :: gas_vmr(ngas, ncol, nlev)
    real(r8), intent(in) :: pmid(ncol, nlev), tmid(ncol, nlev)
    real(r8), intent(in) :: pint(ncol, nlev + 1), tint(ncol, nlev + 1)
    real(r8), intent(out) :: tau(ncol, nlev, ngpt)
    real(r8), intent(out) :: lay_src(ncol, nlev, ngpt)
    real(r8), intent(out) :: lev_src_inc(ncol, nlev, ngpt)
    real(r8), intent(out) :: lev_src_dec(ncol, nlev, ngpt)
    real(r8), intent(out) :: sfc_src(ncol, ngpt)

    type(ty_gas_concs) :: concs
    type(ty_optical_props_1scl) :: optics
    type(ty_source_func_lw) :: sources

    call set_concs(ngas, ncol, nlev, gas_vmr, concs)
    call stop_on_err(optics%alloc_1scl(ncol, nlev, kdist_lw_drv))
    call stop_on_err(sources%alloc(ncol, nlev, kdist_lw_drv))
    ! mirror rrtmgp_interface rrtmgp_run_lw: t_sfc = tint(:, nlev+1)
    call stop_on_err(kdist_lw_drv%gas_optics(pmid, pint, tmid, &
         tint(:, nlev + 1), concs, optics, sources, tlev=tint))
    tau = optics%tau
    lay_src = sources%lay_source
    lev_src_inc = sources%lev_source_inc
    lev_src_dec = sources%lev_source_dec
    sfc_src = sources%sfc_source
  end subroutine drv_gas_optics_lw

  ! ------------------------------------------------------------------
  ! cloud optics
  ! ------------------------------------------------------------------
  subroutine drv_cloud_optics_sw(ncol, nlev, nbnd, do_snow, cld, &
       cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des, rel, rei, &
       tau_out, ssa_out, asm_out, liq_tau_out, ice_tau_out, snw_tau_out)
    use cam_optics, only: get_cloud_optics_sw
    integer, intent(in) :: ncol, nlev, nbnd, do_snow
    real(r8), intent(in), dimension(ncol, nlev) :: cld, cldfsnow, &
         iclwp, iciwp, icswp, lambdac, mu, dei, des, rel, rei
    real(r8), intent(out), dimension(ncol, nlev, nbnd) :: tau_out, &
         ssa_out, asm_out, liq_tau_out, ice_tau_out, snw_tau_out
    call get_cloud_optics_sw(ncol, nlev, nbnd, do_snow /= 0, cld, &
         cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des, rel, rei, &
         tau_out, ssa_out, asm_out, liq_tau_out, ice_tau_out, snw_tau_out)
  end subroutine drv_cloud_optics_sw

  subroutine drv_cloud_optics_lw(ncol, nlev, nbnd, do_snow, cld, &
       cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des, rei, &
       tau_out, liq_tau_out, ice_tau_out, snw_tau_out)
    use cam_optics, only: get_cloud_optics_lw
    integer, intent(in) :: ncol, nlev, nbnd, do_snow
    real(r8), intent(in), dimension(ncol, nlev) :: cld, cldfsnow, &
         iclwp, iciwp, icswp, lambdac, mu, dei, des, rei
    real(r8), intent(out), dimension(ncol, nlev, nbnd) :: tau_out, &
         liq_tau_out, ice_tau_out, snw_tau_out
    call get_cloud_optics_lw(ncol, nlev, nbnd, do_snow /= 0, cld, &
         cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des, rei, &
         tau_out, liq_tau_out, ice_tau_out, snw_tau_out)
  end subroutine drv_cloud_optics_lw

  ! ------------------------------------------------------------------
  ! MCICA
  ! ------------------------------------------------------------------
  subroutine drv_mcica_mask(ngpt, ncol, nlev, changeseed, pmid, cldfrac, &
                            mask_out)
    use mcica_subcol_gen, only: mcica_subcol_mask
    integer, intent(in) :: ngpt, ncol, nlev, changeseed
    real(r8), intent(in) :: pmid(ncol, nlev), cldfrac(ncol, nlev)
    integer, intent(out) :: mask_out(ngpt, ncol, nlev)
    logical :: iscloudy(ngpt, ncol, nlev)
    call mcica_subcol_mask(ngpt, ncol, nlev, changeseed, pmid, cldfrac, &
                           iscloudy)
    mask_out = merge(1, 0, iscloudy)
  end subroutine drv_mcica_mask

  subroutine drv_sample_sw(ncol, nlev, ngpt, nbnd, gpt2bnd, pmid, cld, &
       cldfsnow, tau_bnd, ssa_bnd, asm_bnd, tau_gpt, ssa_gpt, asm_gpt)
    use cam_optics, only: sample_cloud_optics_sw
    integer, intent(in) :: ncol, nlev, ngpt, nbnd
    integer, intent(in) :: gpt2bnd(ngpt)
    real(r8), intent(in), dimension(ncol, nlev) :: pmid, cld, cldfsnow
    real(r8), intent(in), dimension(ncol, nlev, nbnd) :: tau_bnd, &
         ssa_bnd, asm_bnd
    real(r8), intent(out), dimension(ncol, nlev, ngpt) :: tau_gpt, &
         ssa_gpt, asm_gpt
    call sample_cloud_optics_sw(ncol, nlev, ngpt, gpt2bnd, pmid, cld, &
         cldfsnow, tau_bnd, ssa_bnd, asm_bnd, tau_gpt, ssa_gpt, asm_gpt)
  end subroutine drv_sample_sw

  subroutine drv_sample_lw(ncol, nlev, ngpt, nbnd, gpt2bnd, pmid, cld, &
       cldfsnow, tau_bnd, tau_gpt)
    use cam_optics, only: sample_cloud_optics_lw
    integer, intent(in) :: ncol, nlev, ngpt, nbnd
    integer, intent(in) :: gpt2bnd(ngpt)
    real(r8), intent(in), dimension(ncol, nlev) :: pmid, cld, cldfsnow
    real(r8), intent(in), dimension(ncol, nlev, nbnd) :: tau_bnd
    real(r8), intent(out), dimension(ncol, nlev, ngpt) :: tau_gpt
    call sample_cloud_optics_lw(ncol, nlev, ngpt, gpt2bnd, pmid, cld, &
         cldfsnow, tau_bnd, tau_gpt)
  end subroutine drv_sample_lw

  ! ------------------------------------------------------------------
  ! solvers as EAM drives them (real rrtmgp_run_sw / rrtmgp_run_lw)
  ! ------------------------------------------------------------------
  subroutine drv_run_sw(ngas, nday, nlevrad, ngpt, nbnd, gas_vmr, &
       pmid, tmid, pint, coszrs, alb_dir, alb_dif, &
       cld_tau, cld_ssa, cld_asm, aer_tau, aer_ssa, aer_asm, &
       tsi_scaling, flx_all, flx_clr, bnd_all, bnd_clr)
    use rrtmgp_interface, only: rrtmgp_run_sw
    integer, intent(in) :: ngas, nday, nlevrad, ngpt, nbnd
    real(r8), intent(in) :: gas_vmr(ngas, nday, nlevrad)
    real(r8), intent(in) :: pmid(nday, nlevrad), tmid(nday, nlevrad)
    real(r8), intent(in) :: pint(nday, nlevrad + 1)
    real(r8), intent(in) :: coszrs(nday)
    real(r8), intent(in) :: alb_dir(nbnd, nday), alb_dif(nbnd, nday)
    real(r8), intent(in), dimension(nday, nlevrad, ngpt) :: cld_tau, &
         cld_ssa, cld_asm
    real(r8), intent(in), dimension(nday, nlevrad, nbnd) :: aer_tau, &
         aer_ssa, aer_asm
    real(r8), intent(in) :: tsi_scaling
    real(r8), intent(out) :: flx_all(nday, nlevrad + 1, 4)
    real(r8), intent(out) :: flx_clr(nday, nlevrad + 1, 4)
    real(r8), intent(out) :: bnd_all(nday, nlevrad + 1, nbnd, 4)
    real(r8), intent(out) :: bnd_clr(nday, nlevrad + 1, nbnd, 4)

    call rrtmgp_run_sw(ngas, nday, nlevrad, gas_vmr, pmid, tmid, pint, &
         coszrs, alb_dir, alb_dif, cld_tau, cld_ssa, cld_asm, &
         aer_tau, aer_ssa, aer_asm, &
         flx_all(:, :, 1), flx_all(:, :, 2), flx_all(:, :, 3), flx_all(:, :, 4), &
         bnd_all(:, :, :, 1), bnd_all(:, :, :, 2), bnd_all(:, :, :, 3), bnd_all(:, :, :, 4), &
         flx_clr(:, :, 1), flx_clr(:, :, 2), flx_clr(:, :, 3), flx_clr(:, :, 4), &
         bnd_clr(:, :, :, 1), bnd_clr(:, :, :, 2), bnd_clr(:, :, :, 3), bnd_clr(:, :, :, 4), &
         tsi_scaling)
  end subroutine drv_run_sw

  subroutine drv_run_lw(ngas, ncol, nlevrad, ngpt, nbnd, gas_vmr, &
       pmid, tmid, pint, tint, sfc_emis, cld_tau, aer_tau, &
       flx_all, flx_clr, bnd_all, bnd_clr)
    use rrtmgp_interface, only: rrtmgp_run_lw
    integer, intent(in) :: ngas, ncol, nlevrad, ngpt, nbnd
    real(r8), intent(in) :: gas_vmr(ngas, ncol, nlevrad)
    real(r8), intent(in) :: pmid(ncol, nlevrad), tmid(ncol, nlevrad)
    real(r8), intent(in) :: pint(ncol, nlevrad + 1), tint(ncol, nlevrad + 1)
    real(r8), intent(in) :: sfc_emis(nbnd, ncol)
    real(r8), intent(in) :: cld_tau(ncol, nlevrad, ngpt)
    real(r8), intent(in) :: aer_tau(ncol, nlevrad, nbnd)
    real(r8), intent(out) :: flx_all(ncol, nlevrad + 1, 3)
    real(r8), intent(out) :: flx_clr(ncol, nlevrad + 1, 3)
    real(r8), intent(out) :: bnd_all(ncol, nlevrad + 1, nbnd, 3)
    real(r8), intent(out) :: bnd_clr(ncol, nlevrad + 1, nbnd, 3)

    call rrtmgp_run_lw(ngas, ncol, nlevrad, gas_vmr, pmid, tmid, pint, &
         tint, sfc_emis, cld_tau, aer_tau, &
         flx_all(:, :, 1), flx_all(:, :, 2), flx_all(:, :, 3), &
         bnd_all(:, :, :, 1), bnd_all(:, :, :, 2), bnd_all(:, :, :, 3), &
         flx_clr(:, :, 1), flx_clr(:, :, 2), flx_clr(:, :, 3), &
         bnd_clr(:, :, :, 1), bnd_clr(:, :, :, 2), bnd_clr(:, :, :, 3))
  end subroutine drv_run_lw

  ! ------------------------------------------------------------------
  ! full radiation step, sequenced exactly as radiation_tend +
  ! radiation_driver_sw/lw (aerosol optics and gas vmr as plain inputs,
  ! coszrs precomputed, dosw = dolw = .true., single icall = 0)
  ! ------------------------------------------------------------------
  subroutine drv_rad_step(ncol, nlev, ngas, nswb, nlwb, nswg, nlwg, &
       do_snow, t_in, pmid_in, pint_in, lnpmid_in, lnpint_in, &
       lwup, asdir, asdif, aldir, aldif, coszrs, &
       cld, cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des, rel, rei, &
       gas_vmr, aer_tau_sw, aer_ssa_sw, aer_asm_sw, aer_tau_lw, &
       tsi_scaling, &
       qrs, qrsc, qrl, qrlc, sw_all, sw_clr, lw_all, lw_clr, srf, &
       diag_tmid, diag_tint, diag_alb, diag_cld_tau_bnd_sw, &
       diag_cld_gpt_sw, diag_cld_tau_bnd_lw, diag_cld_gpt_lw)
    use ppgrid,           only: pcols, pver
    use physics_types,    only: physics_state
    use camsrfexch,       only: cam_in_t
    use radiation_state,  only: set_rad_state, ktop, kbot, nlev_rad
    use radiation_utils,  only: clip_values, handle_error, fluxes_t, &
                                initialize_fluxes, reset_fluxes, &
                                free_fluxes, expand_day_fluxes, &
                                calculate_heating_rate
    use rrtmgp_interface, only: rrtmgp_run_sw, rrtmgp_run_lw, &
                                get_min_temperature, get_max_temperature, &
                                get_gpoint_bands_sw, get_gpoint_bands_lw
    use cam_optics,       only: get_cloud_optics_sw, get_cloud_optics_lw, &
                                sample_cloud_optics_sw, sample_cloud_optics_lw
    use radconstants,     only: rrtmg_to_rrtmgp_swbands, &
                                get_sw_spectral_boundaries

    integer, intent(in) :: ncol, nlev, ngas, nswb, nlwb, nswg, nlwg
    integer, intent(in) :: do_snow
    real(r8), intent(in), dimension(ncol, nlev) :: t_in, pmid_in, &
         lnpmid_in, cld, cldfsnow, iclwp, iciwp, icswp, lambdac, mu, &
         dei, des, rel, rei
    real(r8), intent(in), dimension(ncol, nlev + 1) :: pint_in, lnpint_in
    real(r8), intent(in), dimension(ncol) :: lwup, asdir, asdif, aldir, &
         aldif, coszrs
    real(r8), intent(in) :: gas_vmr(ngas, ncol, nlev)
    real(r8), intent(in), dimension(ncol, nlev, nswb) :: aer_tau_sw, &
         aer_ssa_sw, aer_asm_sw
    real(r8), intent(in), dimension(ncol, nlev, nlwb) :: aer_tau_lw
    real(r8), intent(in) :: tsi_scaling

    real(r8), intent(out), dimension(ncol, nlev) :: qrs, qrsc, qrl, qrlc
    real(r8), intent(out) :: sw_all(ncol, nlev + 2, 4)  ! up,dn,net,dndir
    real(r8), intent(out) :: sw_clr(ncol, nlev + 2, 4)
    real(r8), intent(out) :: lw_all(ncol, nlev + 2, 3)  ! up,dn,net
    real(r8), intent(out) :: lw_clr(ncol, nlev + 2, 3)
    real(r8), intent(out) :: srf(ncol, 11)
    real(r8), intent(out) :: diag_tmid(ncol, nlev + 1)
    real(r8), intent(out) :: diag_tint(ncol, nlev + 2)
    real(r8), intent(out) :: diag_alb(nswb, ncol, 2)
    real(r8), intent(out) :: diag_cld_tau_bnd_sw(ncol, nlev, nswb)
    real(r8), intent(out) :: diag_cld_gpt_sw(ncol, nlev, nswg, 3)
    real(r8), intent(out) :: diag_cld_tau_bnd_lw(ncol, nlev, nlwb)
    real(r8), intent(out) :: diag_cld_gpt_lw(ncol, nlev, nlwg)

    ! locals mirroring radiation_tend
    type(physics_state) :: state
    type(cam_in_t) :: cam_in
    type(fluxes_t) :: fluxes_allsky, fluxes_clrsky
    type(fluxes_t) :: fluxes_allsky_day, fluxes_clrsky_day
    real(r8), dimension(ncol, nlev + 1) :: tmid, pmid
    real(r8), dimension(ncol, nlev + 2) :: pint, tint
    real(r8), dimension(nswb, ncol) :: albedo_dir, albedo_dif
    real(r8), dimension(ncol, nlev, nswg) :: cld_tau_gpt_sw, &
         cld_ssa_gpt_sw, cld_asm_gpt_sw
    real(r8), dimension(ncol, nlev, nswb) :: cld_tau_bnd_sw, &
         cld_ssa_bnd_sw, cld_asm_bnd_sw, atau_sw, assa_sw, aasm_sw
    real(r8), dimension(ncol, nlev, nswb) :: liq_tau_bnd_sw, &
         ice_tau_bnd_sw, snw_tau_bnd_sw
    real(r8), dimension(ncol, nlev, nlwb) :: cld_tau_bnd_lw, atau_lw, &
         liq_tau_bnd_lw, ice_tau_bnd_lw, snw_tau_bnd_lw
    real(r8), dimension(ncol, nlev, nlwg) :: cld_tau_gpt_lw
    integer :: gpb_sw(nswg), gpb_lw(nlwg)
    integer :: day_indices(ncol), night_indices(ncol)
    integer :: nday, nnight, icol, ilay, iday
    ! radiation_driver_sw day-compressed / padded arrays
    real(r8), dimension(ncol) :: coszrs_day
    real(r8), dimension(nswb, ncol) :: albedo_dir_day, albedo_dif_day
    real(r8), dimension(ncol, nlev + 1) :: pmid_day, tmid_day
    real(r8), dimension(ncol, nlev + 2) :: pint_day
    real(r8), dimension(ngas, ncol, nlev) :: gas_vmr_day
    real(r8), dimension(ngas, ncol, nlev + 1) :: gas_vmr_rad
    real(r8), dimension(ncol, nlev, nswg) :: cld_tau_gpt_day, &
         cld_ssa_gpt_day, cld_asm_gpt_day
    real(r8), dimension(ncol, nlev, nswb) :: aer_tau_bnd_day, &
         aer_ssa_bnd_day, aer_asm_bnd_day
    real(r8), dimension(ncol, nlev + 1, nswg) :: cld_tau_gpt_rad, &
         cld_ssa_gpt_rad, cld_asm_gpt_rad
    real(r8), dimension(ncol, nlev + 1, nswb) :: aer_tau_bnd_rad, &
         aer_ssa_bnd_rad, aer_asm_bnd_rad
    real(r8), dimension(ncol, nlev + 1, nlwg) :: cld_tau_gpt_rad_lw
    real(r8), dimension(ncol, nlev + 1, nlwb) :: aer_tau_bnd_rad_lw
    real(r8) :: surface_emissivity(nlwb, ncol)
    real(r8), dimension(ncol, nlev + 1) :: qrl_rad, qrlc_rad
    real(r8), dimension(ncol, nlev) :: c_cldf
    real(r8) :: flux_dn_diffuse(nswb)

    if (nlev /= pver) call stop_on_err('drv_rad_step: nlev /= pver')
    if (ncol > pcols) call stop_on_err('drv_rad_step: ncol > pcols')

    ! fill the state / cam_in containers
    state%ncol = ncol
    state%t(1:ncol, :) = t_in
    state%pmid(1:ncol, :) = pmid_in
    state%pint(1:ncol, :) = pint_in
    state%lnpmid(1:ncol, :) = lnpmid_in
    state%lnpint(1:ncol, :) = lnpint_in
    cam_in%lwup(1:ncol) = lwup
    cam_in%asdir(1:ncol) = asdir
    cam_in%asdif(1:ncol) = asdif
    cam_in%aldir(1:ncol) = aldir
    cam_in%aldif(1:ncol) = aldif

    ! combined cloud/snow fraction (radiation_tend)
    do icol = 1, ncol
      do ilay = 1, nlev
        c_cldf(icol, ilay) = max(cld(icol, ilay), cldfsnow(icol, ilay))
      end do
    end do

    qrsc = 0._r8
    qrlc = 0._r8

    call set_rad_state(state, cam_in, tmid(1:ncol, 1:nlev_rad), &
         tint(1:ncol, 1:nlev_rad + 1), pmid(1:ncol, 1:nlev_rad), &
         pint(1:ncol, 1:nlev_rad + 1))

    call handle_error(clip_values(tmid(1:ncol, 1:nlev_rad), &
         get_min_temperature(), get_max_temperature(), 'drv tmid'), &
         fatal=.false., warn=.false.)
    call handle_error(clip_values(tint(1:ncol, 1:nlev_rad + 1), &
         get_min_temperature(), get_max_temperature(), 'drv tint'), &
         fatal=.false., warn=.false.)
    diag_tmid = tmid
    diag_tint = tint

    ! ================= shortwave =================
    call initialize_fluxes(ncol, nlev_rad + 1, nswb, fluxes_allsky, &
                           do_direct=.true.)
    call initialize_fluxes(ncol, nlev_rad + 1, nswb, fluxes_clrsky, &
                           do_direct=.true.)

    call set_albedo_x(cam_in, ncol, nswb, albedo_dir(:, 1:ncol), &
                      albedo_dif(:, 1:ncol))
    diag_alb(:, :, 1) = albedo_dir(:, 1:ncol)
    diag_alb(:, :, 2) = albedo_dif(:, 1:ncol)

    cld_tau_gpt_sw = 0._r8
    cld_ssa_gpt_sw = 0._r8
    cld_asm_gpt_sw = 0._r8
    call get_cloud_optics_sw(ncol, nlev, nswb, do_snow /= 0, cld, &
         cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des, rel, rei, &
         cld_tau_bnd_sw, cld_ssa_bnd_sw, cld_asm_bnd_sw, &
         liq_tau_bnd_sw, ice_tau_bnd_sw, snw_tau_bnd_sw)
    ! reorder bands to RRTMGP order (radiation_tend, via reordered())
    do icol = 1, size(cld_tau_bnd_sw, 1)
      do ilay = 1, size(cld_tau_bnd_sw, 2)
        cld_tau_bnd_sw(icol, ilay, :) = reordered_x( &
             cld_tau_bnd_sw(icol, ilay, :), rrtmg_to_rrtmgp_swbands)
        cld_ssa_bnd_sw(icol, ilay, :) = reordered_x( &
             cld_ssa_bnd_sw(icol, ilay, :), rrtmg_to_rrtmgp_swbands)
        cld_asm_bnd_sw(icol, ilay, :) = reordered_x( &
             cld_asm_bnd_sw(icol, ilay, :), rrtmg_to_rrtmgp_swbands)
      end do
    end do
    diag_cld_tau_bnd_sw = cld_tau_bnd_sw

    call get_gpoint_bands_sw(gpb_sw)
    call sample_cloud_optics_sw(ncol, nlev, nswg, gpb_sw, state%pmid, &
         cld, cldfsnow, cld_tau_bnd_sw, cld_ssa_bnd_sw, cld_asm_bnd_sw, &
         cld_tau_gpt_sw, cld_ssa_gpt_sw, cld_asm_gpt_sw)

    ! aerosol band reorder (radiation_tend does this after
    ! set_aerosol_optics_sw; aerosol optics enter in RRTMG order)
    atau_sw = aer_tau_sw
    assa_sw = aer_ssa_sw
    aasm_sw = aer_asm_sw
    do icol = 1, size(atau_sw, 1)
      do ilay = 1, size(atau_sw, 2)
        atau_sw(icol, ilay, :) = reordered_x(atau_sw(icol, ilay, :), &
                                             rrtmg_to_rrtmgp_swbands)
        assa_sw(icol, ilay, :) = reordered_x(assa_sw(icol, ilay, :), &
                                             rrtmg_to_rrtmgp_swbands)
        aasm_sw(icol, ilay, :) = reordered_x(aasm_sw(icol, ilay, :), &
                                             rrtmg_to_rrtmgp_swbands)
      end do
    end do

    call handle_error(clip_values(cld_tau_gpt_sw(1:ncol, :, :), 0._r8, &
         huge(cld_tau_gpt_sw), 'cld_tau_gpt_sw', tolerance=1e-10_r8))
    call handle_error(clip_values(cld_ssa_gpt_sw(1:ncol, :, :), 0._r8, &
         1._r8, 'cld_ssa_gpt_sw', tolerance=1e-10_r8))
    call handle_error(clip_values(cld_asm_gpt_sw(1:ncol, :, :), -1._r8, &
         1._r8, 'cld_asm_gpt_sw', tolerance=1e-10_r8))
    call handle_error(clip_values(atau_sw(1:ncol, :, :), 0._r8, &
         huge(atau_sw), 'aer_tau_bnd_sw', tolerance=1e-10_r8))
    call handle_error(clip_values(assa_sw(1:ncol, :, :), 0._r8, 1._r8, &
         'aer_ssa_bnd_sw', tolerance=1e-10_r8))
    call handle_error(clip_values(aasm_sw(1:ncol, :, :), -1._r8, 1._r8, &
         'aer_asm_bnd_sw', tolerance=1e-10_r8))
    diag_cld_gpt_sw(:, :, :, 1) = cld_tau_gpt_sw
    diag_cld_gpt_sw(:, :, :, 2) = cld_ssa_gpt_sw
    diag_cld_gpt_sw(:, :, :, 3) = cld_asm_gpt_sw

    ! ---- radiation_driver_sw body ----
    call set_daynight_indices_x(coszrs(1:ncol), day_indices(1:ncol), &
                                night_indices(1:ncol))
    nday = count(day_indices(1:ncol) > 0)
    nnight = count(night_indices(1:ncol) > 0)

    if (nday == 0) then
      call reset_fluxes(fluxes_allsky)
      call reset_fluxes(fluxes_clrsky)
      qrs(1:ncol, 1:nlev) = 0._r8
      qrsc(1:ncol, 1:nlev) = 0._r8
    else
      do iday = 1, nday
        icol = day_indices(iday)
        tmid_day(iday, :) = tmid(icol, :)
        pmid_day(iday, :) = pmid(icol, :)
        pint_day(iday, :) = pint(icol, :)
        albedo_dir_day(:, iday) = albedo_dir(:, icol)
        albedo_dif_day(:, iday) = albedo_dif(:, icol)
        coszrs_day(iday) = coszrs(icol)
        gas_vmr_day(:, iday, :) = gas_vmr(:, icol, :)
        cld_tau_gpt_day(iday, :, :) = cld_tau_gpt_sw(icol, :, :)
        cld_ssa_gpt_day(iday, :, :) = cld_ssa_gpt_sw(icol, :, :)
        cld_asm_gpt_day(iday, :, :) = cld_asm_gpt_sw(icol, :, :)
        aer_tau_bnd_day(iday, :, :) = atau_sw(icol, :, :)
        aer_ssa_bnd_day(iday, :, :) = assa_sw(icol, :, :)
        aer_asm_bnd_day(iday, :, :) = aasm_sw(icol, :, :)
      end do

      call initialize_fluxes(nday, nlev_rad + 1, nswb, fluxes_allsky_day, &
                             do_direct=.true.)
      call initialize_fluxes(nday, nlev_rad + 1, nswb, fluxes_clrsky_day, &
                             do_direct=.true.)

      cld_tau_gpt_rad = 0._r8
      cld_ssa_gpt_rad = 0._r8
      cld_asm_gpt_rad = 0._r8
      cld_tau_gpt_rad(1:nday, ktop:kbot, :) = cld_tau_gpt_day(1:nday, 1:nlev, :)
      cld_ssa_gpt_rad(1:nday, ktop:kbot, :) = cld_ssa_gpt_day(1:nday, 1:nlev, :)
      cld_asm_gpt_rad(1:nday, ktop:kbot, :) = cld_asm_gpt_day(1:nday, 1:nlev, :)
      aer_tau_bnd_rad = 0._r8
      aer_ssa_bnd_rad = 0._r8
      aer_asm_bnd_rad = 0._r8
      aer_tau_bnd_rad(1:nday, ktop:kbot, :) = aer_tau_bnd_day(1:nday, :, :)
      aer_ssa_bnd_rad(1:nday, ktop:kbot, :) = aer_ssa_bnd_day(1:nday, :, :)
      aer_asm_bnd_rad(1:nday, ktop:kbot, :) = aer_asm_bnd_day(1:nday, :, :)
      gas_vmr_rad(:, 1:nday, 1) = gas_vmr_day(:, 1:nday, 1)
      gas_vmr_rad(:, 1:nday, ktop:kbot) = gas_vmr_day(:, 1:nday, 1:nlev)

      call rrtmgp_run_sw(ngas, nday, nlev_rad, &
           gas_vmr_rad(:, 1:nday, 1:nlev_rad), &
           pmid_day(1:nday, 1:nlev_rad), &
           tmid_day(1:nday, 1:nlev_rad), &
           pint_day(1:nday, 1:nlev_rad + 1), &
           coszrs_day(1:nday), &
           albedo_dir_day(1:nswb, 1:nday), &
           albedo_dif_day(1:nswb, 1:nday), &
           cld_tau_gpt_rad(1:nday, 1:nlev_rad, 1:nswg), &
           cld_ssa_gpt_rad(1:nday, 1:nlev_rad, 1:nswg), &
           cld_asm_gpt_rad(1:nday, 1:nlev_rad, 1:nswg), &
           aer_tau_bnd_rad(1:nday, 1:nlev_rad, 1:nswb), &
           aer_ssa_bnd_rad(1:nday, 1:nlev_rad, 1:nswb), &
           aer_asm_bnd_rad(1:nday, 1:nlev_rad, 1:nswb), &
           fluxes_allsky_day%flux_up, fluxes_allsky_day%flux_dn, &
           fluxes_allsky_day%flux_net, fluxes_allsky_day%flux_dn_dir, &
           fluxes_allsky_day%bnd_flux_up, fluxes_allsky_day%bnd_flux_dn, &
           fluxes_allsky_day%bnd_flux_net, fluxes_allsky_day%bnd_flux_dn_dir, &
           fluxes_clrsky_day%flux_up, fluxes_clrsky_day%flux_dn, &
           fluxes_clrsky_day%flux_net, fluxes_clrsky_day%flux_dn_dir, &
           fluxes_clrsky_day%bnd_flux_up, fluxes_clrsky_day%bnd_flux_dn, &
           fluxes_clrsky_day%bnd_flux_net, fluxes_clrsky_day%bnd_flux_dn_dir, &
           tsi_scaling)

      call expand_day_fluxes(fluxes_allsky_day, fluxes_allsky, &
                             day_indices(1:nday))
      call expand_day_fluxes(fluxes_clrsky_day, fluxes_clrsky, &
                             day_indices(1:nday))
      call free_fluxes(fluxes_allsky_day)
      call free_fluxes(fluxes_clrsky_day)

      call calculate_heating_rate( &
           fluxes_allsky%flux_up(1:ncol, ktop:kbot + 1), &
           fluxes_allsky%flux_dn(1:ncol, ktop:kbot + 1), &
           pint(1:ncol, ktop:kbot + 1), qrs(1:ncol, 1:nlev))
      call calculate_heating_rate( &
           fluxes_clrsky%flux_up(1:ncol, ktop:kbot + 1), &
           fluxes_clrsky%flux_dn(1:ncol, ktop:kbot + 1), &
           pint(1:ncol, ktop:kbot + 1), qrsc(1:ncol, 1:nlev))
    end if

    sw_all(:, :, 1) = fluxes_allsky%flux_up(1:ncol, :)
    sw_all(:, :, 2) = fluxes_allsky%flux_dn(1:ncol, :)
    sw_all(:, :, 3) = fluxes_allsky%flux_net(1:ncol, :)
    sw_all(:, :, 4) = fluxes_allsky%flux_dn_dir(1:ncol, :)
    sw_clr(:, :, 1) = fluxes_clrsky%flux_up(1:ncol, :)
    sw_clr(:, :, 2) = fluxes_clrsky%flux_dn(1:ncol, :)
    sw_clr(:, :, 3) = fluxes_clrsky%flux_net(1:ncol, :)
    sw_clr(:, :, 4) = fluxes_clrsky%flux_dn_dir(1:ncol, :)

    ! set_net_fluxes_sw + export_surface_fluxes('shortwave'), verbatim
    srf = 0._r8
    do icol = 1, ncol
      srf(icol, 1) = fluxes_allsky%flux_dn(icol, kbot + 1)           ! fsds
      srf(icol, 2) = fluxes_allsky%flux_dn(icol, kbot + 1) &
                   - fluxes_allsky%flux_up(icol, kbot + 1)           ! fsns
      srf(icol, 3) = fluxes_allsky%flux_dn(icol, ktop) &
                   - fluxes_allsky%flux_up(icol, ktop)               ! fsnt
      ! soll / sols / solld / solsd (bands 1-9 NIR, 10 split, 11-14 vis)
      flux_dn_diffuse = fluxes_allsky%bnd_flux_dn(icol, kbot + 1, :) &
                      - fluxes_allsky%bnd_flux_dn_dir(icol, kbot + 1, :)
      srf(icol, 6) = sum(fluxes_allsky%bnd_flux_dn_dir(icol, kbot + 1, 1:9)) &
                   + 0.5_r8 * fluxes_allsky%bnd_flux_dn_dir(icol, kbot + 1, 10)
      srf(icol, 7) = 0.5_r8 * fluxes_allsky%bnd_flux_dn_dir(icol, kbot + 1, 10) &
                   + sum(fluxes_allsky%bnd_flux_dn_dir(icol, kbot + 1, 11:14))
      srf(icol, 8) = sum(flux_dn_diffuse(1:9)) + 0.5_r8 * flux_dn_diffuse(10)
      srf(icol, 9) = 0.5_r8 * flux_dn_diffuse(10) + sum(flux_dn_diffuse(11:14))
      srf(icol, 10) = fluxes_allsky%flux_net(icol, kbot + 1)         ! netsw
    end do
    call free_fluxes(fluxes_allsky)
    call free_fluxes(fluxes_clrsky)

    ! ================= longwave =================
    call initialize_fluxes(ncol, nlev_rad + 1, nlwb, fluxes_allsky)
    call initialize_fluxes(ncol, nlev_rad + 1, nlwb, fluxes_clrsky)

    cld_tau_gpt_lw = 0._r8
    call get_cloud_optics_lw(ncol, nlev, nlwb, do_snow /= 0, cld, &
         cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des, rei, &
         cld_tau_bnd_lw, liq_tau_bnd_lw, ice_tau_bnd_lw, snw_tau_bnd_lw)
    call get_gpoint_bands_lw(gpb_lw)
    call sample_cloud_optics_lw(ncol, nlev, nlwg, gpb_lw, state%pmid, &
         cld, cldfsnow, cld_tau_bnd_lw, cld_tau_gpt_lw)
    diag_cld_tau_bnd_lw = cld_tau_bnd_lw

    atau_lw = aer_tau_lw
    call handle_error(clip_values(cld_tau_gpt_lw(1:ncol, :, :), 0._r8, &
         huge(cld_tau_gpt_lw), 'cld_tau_gpt_lw', tolerance=1e-10_r8))
    call handle_error(clip_values(atau_lw(1:ncol, :, :), 0._r8, &
         huge(atau_lw), 'aer_tau_bnd_lw', tolerance=1e-10_r8))
    diag_cld_gpt_lw = cld_tau_gpt_lw

    ! ---- radiation_driver_lw body ----
    surface_emissivity(1:nlwb, 1:ncol) = 1.0_r8
    cld_tau_gpt_rad_lw = 0._r8
    cld_tau_gpt_rad_lw(:, ktop:kbot, :) = cld_tau_gpt_lw(:, :, :)
    aer_tau_bnd_rad_lw = 0._r8
    aer_tau_bnd_rad_lw(:, ktop:kbot, :) = atau_lw(:, :, :)
    gas_vmr_rad(:, 1:ncol, 1) = gas_vmr(:, 1:ncol, 1)
    gas_vmr_rad(:, 1:ncol, ktop:kbot) = gas_vmr(:, 1:ncol, :)

    call rrtmgp_run_lw(ngas, ncol, nlev_rad, gas_vmr_rad(:, 1:ncol, :), &
         pmid(1:ncol, 1:nlev_rad), tmid(1:ncol, 1:nlev_rad), &
         pint(1:ncol, 1:nlev_rad + 1), tint(1:ncol, 1:nlev_rad + 1), &
         surface_emissivity(1:nlwb, 1:ncol), &
         cld_tau_gpt_rad_lw(1:ncol, :, :), aer_tau_bnd_rad_lw(1:ncol, :, :), &
         fluxes_allsky%flux_up, fluxes_allsky%flux_dn, fluxes_allsky%flux_net, &
         fluxes_allsky%bnd_flux_up, fluxes_allsky%bnd_flux_dn, &
         fluxes_allsky%bnd_flux_net, &
         fluxes_clrsky%flux_up, fluxes_clrsky%flux_dn, fluxes_clrsky%flux_net, &
         fluxes_clrsky%bnd_flux_up, fluxes_clrsky%bnd_flux_dn, &
         fluxes_clrsky%bnd_flux_net)

    call calculate_heating_rate(fluxes_allsky%flux_up, &
         fluxes_allsky%flux_dn, pint(1:ncol, 1:nlev_rad + 1), &
         qrl_rad(1:ncol, 1:nlev_rad))
    call calculate_heating_rate(fluxes_clrsky%flux_up, &
         fluxes_clrsky%flux_dn, pint(1:ncol, 1:nlev_rad + 1), &
         qrlc_rad(1:ncol, 1:nlev_rad))
    qrl(1:ncol, 1:nlev) = qrl_rad(1:ncol, ktop:kbot)
    qrlc(1:ncol, 1:nlev) = qrlc_rad(1:ncol, ktop:kbot)

    lw_all(:, :, 1) = fluxes_allsky%flux_up(1:ncol, :)
    lw_all(:, :, 2) = fluxes_allsky%flux_dn(1:ncol, :)
    lw_all(:, :, 3) = fluxes_allsky%flux_net(1:ncol, :)
    lw_clr(:, :, 1) = fluxes_clrsky%flux_up(1:ncol, :)
    lw_clr(:, :, 2) = fluxes_clrsky%flux_dn(1:ncol, :)
    lw_clr(:, :, 3) = fluxes_clrsky%flux_net(1:ncol, :)

    ! set_net_fluxes_lw + export_surface_fluxes('longwave'), verbatim
    do icol = 1, ncol
      srf(icol, 4) = fluxes_allsky%flux_up(icol, kbot + 1) &
                   - fluxes_allsky%flux_dn(icol, kbot + 1)           ! flns
      srf(icol, 5) = fluxes_allsky%flux_up(icol, ktop) &
                   - fluxes_allsky%flux_dn(icol, ktop)               ! flnt
      srf(icol, 11) = fluxes_allsky%flux_dn(icol, kbot + 1)          ! flwds
    end do
    call free_fluxes(fluxes_allsky)
    call free_fluxes(fluxes_clrsky)
  end subroutine drv_rad_step

  ! ------------------------------------------------------------------
  ! transcriptions of radiation.F90 private helpers (verbatim logic)
  ! ------------------------------------------------------------------
  function reordered_x(array_in, new_indexing) result(array_out)
    real(r8), intent(in) :: array_in(:)
    integer, intent(in) :: new_indexing(:)
    real(r8), dimension(size(array_in)) :: array_out
    integer :: ii
    do ii = 1, size(new_indexing)
      array_out(ii) = array_in(new_indexing(ii))
    end do
  end function reordered_x

  subroutine set_daynight_indices_x(coszrs, day_indices, night_indices)
    real(r8), intent(in) :: coszrs(:)
    integer, intent(inout) :: day_indices(:), night_indices(:)
    integer :: icol, iday, inight
    day_indices(:) = 0
    night_indices(:) = 0
    iday = 0
    inight = 0
    do icol = 1, size(coszrs)
      if (coszrs(icol) > 0._r8) then
        iday = iday + 1
        day_indices(iday) = icol
      else
        inight = inight + 1
        night_indices(inight) = icol
      end if
    end do
  end subroutine set_daynight_indices_x

  logical function is_visible_x(wavenumber)
    real(r8), intent(in) :: wavenumber
    real(r8), parameter :: visible_wavenumber_threshold = 14286._r8
    is_visible_x = wavenumber > visible_wavenumber_threshold
  end function is_visible_x

  subroutine set_albedo_x(cam_in, ncol, nswbands_in, albedo_dir, albedo_dif)
    ! radiation.F90 set_albedo, verbatim (clip via the real clip_values)
    use camsrfexch,      only: cam_in_t
    use radconstants,    only: get_sw_spectral_boundaries, &
                               rrtmg_to_rrtmgp_swbands, nswbands
    use radiation_utils, only: clip_values, handle_error
    type(cam_in_t), intent(in) :: cam_in
    integer, intent(in) :: ncol, nswbands_in
    real(r8), intent(inout) :: albedo_dir(:, :), albedo_dif(:, :)
    real(r8), dimension(nswbands) :: lower_bounds, upper_bounds
    integer :: iband

    albedo_dir(:, :) = 0._r8
    albedo_dif(:, :) = 0._r8
    call get_sw_spectral_boundaries(lower_bounds, upper_bounds, 'cm^-1')
    lower_bounds = reordered_x(lower_bounds, rrtmg_to_rrtmgp_swbands)
    upper_bounds = reordered_x(upper_bounds, rrtmg_to_rrtmgp_swbands)
    do iband = 1, nswbands
      if (is_visible_x(lower_bounds(iband)) .and. &
          is_visible_x(upper_bounds(iband))) then
        albedo_dir(iband, 1:ncol) = cam_in%asdir(1:ncol)
        albedo_dif(iband, 1:ncol) = cam_in%asdif(1:ncol)
      else if (.not. is_visible_x(lower_bounds(iband)) .and. &
               .not. is_visible_x(upper_bounds(iband))) then
        albedo_dir(iband, 1:ncol) = cam_in%aldir(1:ncol)
        albedo_dif(iband, 1:ncol) = cam_in%aldif(1:ncol)
      else
        ! NOTE 0.5 is a DEFAULT-REAL literal in radiation.F90
        albedo_dir(iband, 1:ncol) = 0.5 * (cam_in%aldir(1:ncol) &
                                           + cam_in%asdir(1:ncol))
        albedo_dif(iband, 1:ncol) = 0.5 * (cam_in%aldif(1:ncol) &
                                           + cam_in%asdif(1:ncol))
      end if
    end do
    call handle_error(clip_values(albedo_dir, 0._r8, 1._r8, &
         'set_albedo: albedo_dir', tolerance=0.01_r8))
    call handle_error(clip_values(albedo_dif, 0._r8, 1._r8, &
         'set_albedo: albedo_dif', tolerance=0.01_r8))
  end subroutine set_albedo_x

  ! ------------------------------------------------------------------
  ! internal helpers
  ! ------------------------------------------------------------------
  subroutine set_concs(ngas, ncol, nlev, gas_vmr, concs)
    ! mirrors rrtmgp_interface's private set_gas_concentrations
    use mo_gas_concentrations, only: ty_gas_concs
    use mo_rrtmgp_util_string, only: lower_case
    integer, intent(in) :: ngas, ncol, nlev
    real(r8), intent(in) :: gas_vmr(ngas, ncol, nlev)
    type(ty_gas_concs), intent(out) :: concs
    character(len=3), dimension(8) :: gases
    character(len=32) :: gas_names_lower(8)
    integer :: igas
    call active_gas_list(gases)
    do igas = 1, ngas
      gas_names_lower(igas) = trim(lower_case(gases(igas)))
    end do
    call stop_on_err(concs%init(gas_names_lower(1:ngas)))
    do igas = 1, ngas
      call stop_on_err(concs%set_vmr(trim(gas_names_lower(igas)), &
                                     gas_vmr(igas, 1:ncol, 1:nlev)))
    end do
  end subroutine set_concs

  subroutine drv_load_kdists(sw_file, lw_file)
    ! load the driver-owned k-dists for the kernel-level goldens
    use mo_load_coefficients,  only: load_and_init
    use mo_gas_concentrations, only: ty_gas_concs
    use mo_rrtmgp_util_string, only: lower_case
    character(len=256), intent(in) :: sw_file, lw_file
    type(ty_gas_concs) :: available
    character(len=3), dimension(8) :: gases
    character(len=32) :: low(8)
    integer :: igas
    call active_gas_list(gases)
    do igas = 1, 8
      low(igas) = trim(lower_case(gases(igas)))
    end do
    call stop_on_err(available%init(low))
    call load_and_init(kdist_sw_drv, trim(sw_file), available)
    call load_and_init(kdist_lw_drv, trim(lw_file), available)
  end subroutine drv_load_kdists

end module rrtmgp_core
