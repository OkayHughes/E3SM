! f2py driver for the ZM precipitation evaporation / snow production
! routine (eam/src/physics/cam/zm/zm_conv.F90 zm_conv_evap):
! below-cloud Sundqvist-type evaporation of convective precip, snow
! melt / production with the latent heat of fusion, and the final
! surface precip / snow fluxes.
!
! zm_conv.F90 is compiled unmodified against the REAL
! cloud_fraction.F90 (zm_conv_evap calls its cldfrc_fice T-based
! ice/snow partition -- real physics, so no stub; the abort-only
! cldfrc_fice placeholder of cloud_fraction_stub.F90 is NOT used
! here). cldfrc_fice runs with the module default top_lev = 1
! (cldfrc_init is not called): in production top_lev =
! trop_cloud_top_lev sits far above the 40 hPa limcnv cap where every
! ZM precip flux is identically zero, so fsnow at those levels only
! ever multiplies zeros and the choice cannot affect any output.
!
! zm_param subset used by zm_conv_evap: ke (zmconv_ke evaporation
! efficiency; EAMv3 has a single ke -- there is no zmconv_ke_lnd in
! this ZM version), old_snow (snow production treatment; the
! zm_param_t default is .true., and zm_conv_intr only sets it .false.
! under zmconv_microp), and zm_microp (FORCED .false. here -- same
! scope as the zm_conv harness; the prdsnow=microp_st%sprd branch
! needs zm_microphysics, PORTING_PLAN.md row 9). zm_const is filled
! by the real zm_const_set_to_global via zm_conv_main_init.
!
! Arrays are padded to ppgrid pcols = 48 (zm_intr_stubs.F90) because
! the real cldfrc_fice declares (pcols, pver) dummies; padding rows
! are zero-initialized and never read (all loops run 1..ncol).
! tend_s/tend_q are intent(inout) accumulators in the Fortran; the
! intr layer hands in a freshly zeroed ptend, which the driver
! reproduces with explicit input arrays.
module zm_evap_driver
  implicit none
  ! local kind param (a use-associated kind is invisible to f2py's
  ! crackfortran; fbuild passes --f2cmap mapping r8 -> double)
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_init()
    ! Build the wv_saturation SVP table / default scheme (GoffGratch);
    ! zm_conv_evap uses wv_saturation's mixed-phase qsat.
    use wv_saturation, only: wv_sat_init
    call wv_sat_init()
  end subroutine drv_init

  subroutine drv_zm_conv_evap(ncol, nlev, ke, old_snow, time_step, &
       p_mid, p_del, t_mid, q_mid, prdprec, cldfrc, &
       tend_s_in, tend_q_in, prec_in, &
       tend_s, tend_q, tend_s_snwprd, tend_s_snwevmlt, &
       prec, snow, ntprprd, ntsnprd, flxprec, flxsnow)
    use ppgrid,                only: pcols, ppver => pver
    use zm_conv,               only: zm_conv_evap, zm_conv_main_init, &
                                     zm_param
    use zm_microphysics_state, only: zm_microp_st
    use cam_abortutils,        only: endrun

    integer, intent(in) :: ncol, nlev
    real(r8), intent(in) :: ke                        ! [1/mm**0.5/s?] zmconv_ke
    integer, intent(in) :: old_snow                   ! 0/1 flag
    real(r8), intent(in) :: time_step                 ! [s]
    real(r8), intent(in) :: p_mid(ncol, nlev)         ! [Pa]
    real(r8), intent(in) :: p_del(ncol, nlev)         ! [Pa]
    real(r8), intent(in) :: t_mid(ncol, nlev)         ! [K]
    real(r8), intent(in) :: q_mid(ncol, nlev)         ! [kg/kg]
    real(r8), intent(in) :: prdprec(ncol, nlev)       ! [kg/kg/s]
    real(r8), intent(in) :: cldfrc(ncol, nlev)        ! cloud fraction
    real(r8), intent(in) :: tend_s_in(ncol, nlev)     ! [J/kg/s]
    real(r8), intent(in) :: tend_q_in(ncol, nlev)     ! [kg/kg/s]
    real(r8), intent(in) :: prec_in(ncol)             ! [m/s]

    real(r8), intent(out) :: tend_s(ncol, nlev)          ! [J/kg/s]
    real(r8), intent(out) :: tend_q(ncol, nlev)          ! [kg/kg/s]
    real(r8), intent(out) :: tend_s_snwprd(ncol, nlev)   ! [J/kg/s]
    real(r8), intent(out) :: tend_s_snwevmlt(ncol, nlev) ! [J/kg/s]
    real(r8), intent(out) :: prec(ncol)                  ! [m/s]
    real(r8), intent(out) :: snow(ncol)                  ! [m/s]
    real(r8), intent(out) :: ntprprd(ncol, nlev)         ! [kg/kg/s]
    real(r8), intent(out) :: ntsnprd(ncol, nlev)         ! [kg/kg/s]
    real(r8), intent(out) :: flxprec(ncol, nlev + 1)     ! [kg/m2/s]
    real(r8), intent(out) :: flxsnow(ncol, nlev + 1)     ! [kg/m2/s]

    ! pcols-padded working copies (cldfrc_fice dummies are
    ! (ppgrid pcols, ppgrid pver))
    real(r8) :: w_pmid(pcols, nlev), w_pdel(pcols, nlev)
    real(r8) :: w_t(pcols, nlev), w_q(pcols, nlev)
    real(r8) :: w_prd(pcols, nlev), w_cld(pcols, nlev)
    real(r8) :: w_ts(pcols, nlev), w_tq(pcols, nlev)
    real(r8) :: w_tssp(pcols, nlev), w_tsse(pcols, nlev)
    real(r8) :: w_prec(pcols), w_snow(pcols)
    real(r8) :: w_ntp(pcols, nlev), w_nts(pcols, nlev)
    real(r8) :: w_fp(pcols, nlev + 1), w_fs(pcols, nlev + 1)
    type(zm_microp_st) :: microp_st  ! untouched: zm_microp=.false.

    if (ncol > pcols) call endrun('drv_zm_conv_evap: ncol > pcols')
    if (nlev /= ppver) call endrun('drv_zm_conv_evap: nlev /= pver')

    zm_param%ke        = ke
    zm_param%old_snow  = old_snow /= 0
    zm_param%zm_microp = .false.
    ! fills zm_const via the real zm_const_set_to_global (limcnv is
    ! not used by zm_conv_evap; masterproc=.false. mutes the print)
    call zm_conv_main_init(1)

    w_pmid = 0._r8; w_pdel = 1._r8; w_t = 200._r8; w_q = 0._r8
    w_prd = 0._r8; w_cld = 0._r8; w_ts = 0._r8; w_tq = 0._r8
    w_prec = 0._r8
    w_pmid(1:ncol, :) = p_mid
    w_pdel(1:ncol, :) = p_del
    w_t(1:ncol, :)    = t_mid
    w_q(1:ncol, :)    = q_mid
    w_prd(1:ncol, :)  = prdprec
    w_cld(1:ncol, :)  = cldfrc
    w_ts(1:ncol, :)   = tend_s_in
    w_tq(1:ncol, :)   = tend_q_in
    w_prec(1:ncol)    = prec_in

    call zm_conv_evap(pcols, ncol, nlev, nlev + 1, time_step, &
         w_pmid, w_pdel, w_t, w_q, w_prd, w_cld, &
         w_ts, w_tq, w_tssp, w_tsse, &
         w_prec, w_snow, w_ntp, w_nts, w_fp, w_fs, microp_st)

    tend_s          = w_ts(1:ncol, :)
    tend_q          = w_tq(1:ncol, :)
    tend_s_snwprd   = w_tssp(1:ncol, :)
    tend_s_snwevmlt = w_tsse(1:ncol, :)
    prec            = w_prec(1:ncol)
    snow            = w_snow(1:ncol)
    ntprprd         = w_ntp(1:ncol, :)
    ntsnprd         = w_nts(1:ncol, :)
    flxprec         = w_fp(1:ncol, :)
    flxsnow         = w_fs(1:ncol, :)
  end subroutine drv_zm_conv_evap

end module zm_evap_driver
