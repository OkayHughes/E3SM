! f2py driver for the ZM deep-convection main routine
! (eam/src/physics/cam/zm/zm_conv.F90 zm_conv_main): dilute CAPE +
! DCAPE trigger, gather, updraft/downdraft plume properties
! (zm_cloud_properties -> zm_calc_fractional_entrainment /
! zm_downdraft_properties), CAPE closure (zm_closure), cloud-base
! mass-flux limiting/scaling, and output tendencies
! (zm_calc_output_tend) + precip/rliq integrals.
!
! Scope: zm_param%zm_microp is FORCED .false. here (the simple
! in-plume condensation path). The EAMv3 default is zmconv_microp =
! .true., but zm_microphysics.F90 cannot be compiled standalone (it
! needs the aerosol activation stack); goldens and the JAX port
! therefore both use the microp=.false. branch of zm_conv_main --
! documented in PORTING_PLAN.md. zm_conv_evap (and hence cldfrc_fice)
! is NOT driven; MCSP and zm_aero paths are inactive by construction
! (aero/microp_st are only touched under zm_microp).
!
! zm_const is filled by the REAL zm_conv_types zm_const_set_to_global
! called from zm_conv_main_init (physconst stub derives every value
! verbatim from share/util/shr_const_mod.F90; zvir = 1.608 hardcoded
! in zm_conv_types "to avoid non-BFB diffs").
!
! All zm_param values used by zm_conv_main are driver inputs and are
! recorded in the golden metadata. limcnv follows the 40 hPa
! reference-interface rule of zm_conv_intr.F90 (computed by the
! generator); msg = limcnv-1.
!
! gather_index is pre-zeroed: zm_conv_main only writes entries
! 1..lengath of that intent(out) array.
module zm_conv_driver
  implicit none
  ! local kind param (a use-associated kind is invisible to f2py's
  ! crackfortran; fbuild passes --f2cmap mapping r8 -> double)
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_init()
    ! Build the wv_saturation SVP table / default scheme (GoffGratch;
    ! the stubbed namelist read leaves the defaults untouched).
    ! Required by qsat_hPa used in zm_conv_cape and zm_cloud_properties.
    use wv_saturation, only: wv_sat_init
    call wv_sat_init()
  end subroutine drv_init

  subroutine drv_zm_conv_main(ncol, nlev, time_step, is_first_step, &
       limcnv, no_deep_pbl, tau, alfa, ke, dmpdz, tpert_fix, &
       tpert_fac, tiedke_add, c0_lnd, c0_ocn, num_cin, &
       mx_bot_lyr_adj, trig_dcape, trig_ull, clos_dyn_adj, old_snow, &
       t_mid, q_mid, omega, p_mid, p_int, p_del, geos, z_mid, z_int, &
       pbl_hgt, tpert, landfrac, t_star_in, q_star_in, &
       lengath, gather_index, msemax_klev_g, jctop, jcbot, jt, &
       prec, heat, qtnd, cape, dcape, mcon, pflx, zdu, &
       mflx_up, entr_up, detr_up, mflx_dn, entr_dn, p_del_out, &
       dsubcld, ql, rliq, rprd, dlf)
    use zm_conv,               only: zm_conv_main, zm_conv_main_init, zm_param
    use zm_aero_type,          only: zm_aero_t
    use zm_microphysics_state, only: zm_microp_st

    integer, intent(in) :: ncol, nlev
    real(r8), intent(in) :: time_step
    integer, intent(in) :: is_first_step               ! 0/1 flag
    integer, intent(in) :: limcnv                      ! 1-based interface
    integer, intent(in) :: no_deep_pbl                 ! 0/1 flag
    real(r8), intent(in) :: tau, alfa, ke, dmpdz
    integer, intent(in) :: tpert_fix                   ! 0/1 flag
    real(r8), intent(in) :: tpert_fac, tiedke_add, c0_lnd, c0_ocn
    integer, intent(in) :: num_cin, mx_bot_lyr_adj
    integer, intent(in) :: trig_dcape, trig_ull        ! 0/1 flags
    integer, intent(in) :: clos_dyn_adj, old_snow      ! 0/1 flags
    real(r8), intent(in) :: t_mid(ncol, nlev)          ! [K]
    real(r8), intent(in) :: q_mid(ncol, nlev)          ! [kg/kg]
    real(r8), intent(in) :: omega(ncol, nlev)          ! [Pa/s]
    real(r8), intent(in) :: p_mid(ncol, nlev)          ! [Pa]
    real(r8), intent(in) :: p_int(ncol, nlev + 1)      ! [Pa]
    real(r8), intent(in) :: p_del(ncol, nlev)          ! [Pa]
    real(r8), intent(in) :: geos(ncol)                 ! [m2/s2]
    real(r8), intent(in) :: z_mid(ncol, nlev)          ! [m above sfc]
    real(r8), intent(in) :: z_int(ncol, nlev + 1)      ! [m above sfc]
    real(r8), intent(in) :: pbl_hgt(ncol)              ! [m]
    real(r8), intent(in) :: tpert(ncol)                ! [K]
    real(r8), intent(in) :: landfrac(ncol)
    real(r8), intent(in) :: t_star_in(ncol, nlev)      ! DCAPE prev T
    real(r8), intent(in) :: q_star_in(ncol, nlev)      ! DCAPE prev q

    integer, intent(out) :: lengath
    integer, intent(out) :: gather_index(ncol)         ! 1-based
    integer, intent(out) :: msemax_klev_g(ncol)        ! 1-based, gathered
    integer, intent(out) :: jctop(ncol), jcbot(ncol)   ! 1-based, scattered
    integer, intent(out) :: jt(ncol)                   ! 1-based, gathered
    real(r8), intent(out) :: prec(ncol)                ! [m/s]
    real(r8), intent(out) :: heat(ncol, nlev)          ! [W/kg] scattered
    real(r8), intent(out) :: qtnd(ncol, nlev)          ! [kg/kg/s] scattered
    real(r8), intent(out) :: cape(ncol), dcape(ncol)   ! [J/kg]
    real(r8), intent(out) :: mcon(ncol, nlev + 1)      ! [mb/s] scattered
    real(r8), intent(out) :: pflx(ncol, nlev + 1)      ! [kg/m2/s] scattered
    real(r8), intent(out) :: zdu(ncol, nlev)           ! [1/s] scattered
    real(r8), intent(out) :: mflx_up(ncol, nlev)       ! [mb/s] gathered
    real(r8), intent(out) :: entr_up(ncol, nlev)       ! [1/mb] gathered
    real(r8), intent(out) :: detr_up(ncol, nlev)       ! [1/mb] gathered
    real(r8), intent(out) :: mflx_dn(ncol, nlev)       ! [mb/s] gathered
    real(r8), intent(out) :: entr_dn(ncol, nlev)       ! [1/mb] gathered
    real(r8), intent(out) :: p_del_out(ncol, nlev)     ! [mb] gathered
    real(r8), intent(out) :: dsubcld(ncol)             ! [mb] gathered
    real(r8), intent(out) :: ql(ncol, nlev)            ! gathered
    real(r8), intent(out) :: rliq(ncol)                ! [m/s]
    real(r8), intent(out) :: rprd(ncol, nlev)          ! [kg/kg/s] scattered
    real(r8), intent(out) :: dlf(ncol, nlev)           ! [kg/kg/s] scattered

    type(zm_aero_t)        :: aero       ! untouched: zm_microp=.false.
    type(zm_microp_st)     :: microp_st  ! untouched: zm_microp=.false.
    real(r8), pointer      :: t_star(:, :), q_star(:, :)

    ! zm_param subset used by zm_conv_main with zm_microp=.false.;
    ! microp tuning + MCSP fields keep their type defaults (inactive)
    zm_param%tau            = tau
    zm_param%alfa           = alfa
    zm_param%ke             = ke
    zm_param%dmpdz          = dmpdz
    zm_param%tpert_fix      = tpert_fix /= 0
    zm_param%tpert_fac      = tpert_fac
    zm_param%tiedke_add     = tiedke_add
    zm_param%c0_lnd         = c0_lnd
    zm_param%c0_ocn         = c0_ocn
    zm_param%num_cin        = num_cin
    zm_param%mx_bot_lyr_adj = mx_bot_lyr_adj
    zm_param%trig_dcape     = trig_dcape /= 0
    zm_param%trig_ull       = trig_ull /= 0
    zm_param%clos_dyn_adj   = clos_dyn_adj /= 0
    zm_param%old_snow       = old_snow /= 0
    zm_param%zm_microp      = .false.
    ! sets zm_param%limcnv/no_deep_pbl and fills zm_const via the real
    ! zm_const_set_to_global (masterproc=.false. mutes the log print)
    call zm_conv_main_init(limcnv, no_deep_pbl_in=(no_deep_pbl /= 0))

    ! zm_conv_main takes t_star/q_star as pointers
    allocate(t_star(ncol, nlev), q_star(ncol, nlev))
    t_star = t_star_in
    q_star = q_star_in

    ! zm_conv_main only writes gather_index(1:lengath)
    gather_index = 0

    call zm_conv_main(ncol, ncol, nlev, nlev + 1, is_first_step /= 0, &
         time_step, t_mid, q_mid, omega, p_mid, p_int, p_del, geos, &
         z_mid, z_int, pbl_hgt, tpert, landfrac, t_star, q_star, &
         lengath, gather_index, msemax_klev_g, jctop, jcbot, jt, &
         prec, heat, qtnd, cape, dcape, mcon, pflx, zdu, &
         mflx_up, entr_up, detr_up, mflx_dn, entr_dn, p_del_out, &
         dsubcld, ql, rliq, rprd, dlf, aero, microp_st)

    deallocate(t_star, q_star)
  end subroutine drv_zm_conv_main

end module zm_conv_driver
