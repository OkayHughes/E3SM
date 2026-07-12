! f2py driver for the ZM deep-convection dilute CAPE core
! (eam/src/physics/cam/zm/zm_conv_cape.F90): compute_dilute_cape =
! find_mse_max + compute_dilute_parcel + compute_cape_from_parcel,
! with entropy/ientropy/qsat_hPa from zm_conv_util.F90.
!
! Constants (zm_const_t) are filled by the REAL zm_conv_types
! zm_const_set_to_global(), compiled against the physconst stub, which
! derives every value verbatim from share/util/shr_const_mod.F90:
!   pi     = 3.14159265358979323846
!   grav   = SHR_CONST_G      = 9.80616
!   rdair  = SHR_CONST_RDAIR  = SHR_CONST_RGAS/28.966
!   rh2o   = SHR_CONST_RWV    = SHR_CONST_RGAS/18.016
!   cpair  = SHR_CONST_CPDAIR = 1.00464e3
!   cpwv   = SHR_CONST_CPWV   = 1.810e3
!   cpliq  = SHR_CONST_CPFW   = 4.188e3
!   tfreez = SHR_CONST_TKFRZ  = 273.15
!   latvap = SHR_CONST_LATVAP = 2.501e6
!   latice = SHR_CONST_LATICE = 3.337e5
!   epsilo = 18.016/28.966
!   zvir   = 1.608 (hardcoded in zm_const_set_to_global "to avoid
!            non-BFB diffs"; NOT rh2o/rdair-1)
!
! Tunable parameters (zm_param_t) are inputs here; EAMv3 defaults from
! bld/namelist_files/namelist_defaults_eam.xml (phys="default"):
!   zmconv_dmpdz          = -0.7e-3   [1/m]
!   zmconv_tiedke_add     =  0.8      [K]
!   zmconv_tp_fac         =  2.0      (tpert_fac)
!   zmconv_mx_bot_lyr_adj =  1
!   zmconv_cape_cin       =  1        (num_cin)
!   zmconv_tpert_fix      = .true.
!   zmconv_trig_ull       = .true.
!   zmconv_trig_dcape     = .true.
! num_msg = limcnv-1, where limcnv is the first interface with
! reference pressure crossing ZM_upper_limit_pref = 40 hPa
! (zm_conv_intr.F90).
!
! compute_dilute_cape takes pcols as a dummy argument (it does not use
! the ppgrid module), so the driver passes pcols = ncol; no padding
! copies are needed. q_mx/t_mx are optional in the Fortran; the driver
! ALWAYS passes them, so the launch-level T/q save branch executes on
! every call (it only writes those two outputs, never feeds back).
module zm_cape_driver
  implicit none
  ! local kind param (a use-associated kind is invisible to f2py's
  ! crackfortran; fbuild passes --f2cmap mapping r8 -> double)
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_init()
    ! Build the wv_saturation SVP table / default scheme (GoffGratch;
    ! the stubbed namelist read leaves the defaults untouched).
    ! Required by qsat_hPa <- qsat_water inside entropy/ientropy.
    use wv_saturation, only: wv_sat_init
    call wv_sat_init()
  end subroutine drv_init

  subroutine drv_compute_dilute_cape(ncol, nlev, num_cin, num_msg, &
       sp_humidity, temperature, zmid, pmid, pint, pblt, tpert, &
       dmpdz, tiedke_add, tpert_fac, mx_bot_lyr_adj, &
       tpert_fix, trig_ull, trig_dcape, &
       calc_msemax_klev, prev_msemax_klev, &
       use_input_tq_mx, q_mx_in, t_mx_in, &
       parcel_temp, parcel_qsat, msemax_klev, lcl_temperature, &
       lcl_klev, eql_klev, cape, q_mx, t_mx)
    use zm_conv_types, only: zm_const_t, zm_param_t, &
                             zm_const_set_to_global
    use zm_conv_cape, only: compute_dilute_cape
    integer, intent(in) :: ncol, nlev
    integer, intent(in) :: num_cin, num_msg
    real(r8), intent(in) :: sp_humidity(ncol, nlev)
    real(r8), intent(in) :: temperature(ncol, nlev)
    real(r8), intent(in) :: zmid(ncol, nlev)   ! [m]
    real(r8), intent(in) :: pmid(ncol, nlev)   ! [hPa] (as in zm_conv)
    real(r8), intent(in) :: pint(ncol, nlev + 1)  ! [hPa]
    integer, intent(in) :: pblt(ncol)          ! 1-based level index
    real(r8), intent(in) :: tpert(ncol)
    real(r8), intent(in) :: dmpdz, tiedke_add, tpert_fac
    integer, intent(in) :: mx_bot_lyr_adj
    integer, intent(in) :: tpert_fix, trig_ull, trig_dcape  ! 0/1 flags
    integer, intent(in) :: calc_msemax_klev                 ! 0/1 flag
    integer, intent(in) :: prev_msemax_klev(ncol)  ! 1-based
    integer, intent(in) :: use_input_tq_mx                  ! 0/1 flag
    real(r8), intent(in) :: q_mx_in(ncol), t_mx_in(ncol)
    real(r8), intent(out) :: parcel_temp(ncol, nlev)
    real(r8), intent(out) :: parcel_qsat(ncol, nlev)
    integer, intent(out) :: msemax_klev(ncol)      ! 1-based
    real(r8), intent(out) :: lcl_temperature(ncol)
    integer, intent(out) :: lcl_klev(ncol), eql_klev(ncol)  ! 1-based
    real(r8), intent(out) :: cape(ncol)
    real(r8), intent(out) :: q_mx(ncol), t_mx(ncol)

    type(zm_const_t) :: zm_const
    type(zm_param_t) :: zm_param
    integer :: prev_k(ncol)

    call zm_const_set_to_global(zm_const)
    ! CAPE-relevant subset of zm_param; the remaining fields are not
    ! referenced by zm_conv_cape and keep their type defaults.
    zm_param%dmpdz          = dmpdz
    zm_param%tiedke_add     = tiedke_add
    zm_param%tpert_fac      = tpert_fac
    zm_param%mx_bot_lyr_adj = mx_bot_lyr_adj
    zm_param%num_cin        = num_cin
    zm_param%tpert_fix      = tpert_fix /= 0
    zm_param%trig_ull       = trig_ull /= 0
    zm_param%trig_dcape     = trig_dcape /= 0

    prev_k = prev_msemax_klev
    q_mx = q_mx_in
    t_mx = t_mx_in
    msemax_klev = 0
    lcl_klev = 0
    eql_klev = 0
    cape = 0._r8
    parcel_temp = 0._r8
    parcel_qsat = 0._r8
    lcl_temperature = 0._r8

    if (use_input_tq_mx /= 0) then
       call compute_dilute_cape(ncol, ncol, nlev, nlev + 1, &
            num_cin, num_msg, sp_humidity, temperature, zmid, pmid, &
            pint, pblt, tpert, parcel_temp, parcel_qsat, msemax_klev, &
            lcl_temperature, lcl_klev, eql_klev, cape, zm_const, &
            zm_param, calc_msemax_klev /= 0, &
            prev_msemax_klev=prev_k, use_input_tq_mx=.true., &
            q_mx=q_mx, t_mx=t_mx)
    else if (calc_msemax_klev /= 0) then
       call compute_dilute_cape(ncol, ncol, nlev, nlev + 1, &
            num_cin, num_msg, sp_humidity, temperature, zmid, pmid, &
            pint, pblt, tpert, parcel_temp, parcel_qsat, msemax_klev, &
            lcl_temperature, lcl_klev, eql_klev, cape, zm_const, &
            zm_param, .true., q_mx=q_mx, t_mx=t_mx)
    else
       call compute_dilute_cape(ncol, ncol, nlev, nlev + 1, &
            num_cin, num_msg, sp_humidity, temperature, zmid, pmid, &
            pint, pblt, tpert, parcel_temp, parcel_qsat, msemax_klev, &
            lcl_temperature, lcl_klev, eql_klev, cape, zm_const, &
            zm_param, .false., prev_msemax_klev=prev_k, &
            q_mx=q_mx, t_mx=t_mx)
    end if
  end subroutine drv_compute_dilute_cape

end module zm_cape_driver
