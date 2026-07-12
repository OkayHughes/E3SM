! f2py driver for conv_water_4rad
! (eam/src/physics/cam/conv_water.F90): grid-box average liquid/ice
! from stratus + cumulus for radiation.
!
! The pbuf/phys_state plumbing becomes plain arrays: the driver fills
! the physics_buffer stub registry (ICWMRSH/ICWMRDP/ICIMRDP/FICE/
! SH_FRAC/DP_FRAC/AST/REI) and a minimal physics_state (t, pdel, q
! with CLDLIQ=2, CLDICE=3 per the constituents stub), then calls the
! real conv_water_register + conv_water_init + conv_water_4rad.
!
! Switches driven per call:
!  - conv_water_mode: 1 (area-weighted arithmetic) or 2 (arithmetic in
!    emissivity). EAMv3 phys="default" conv_water_in_rad = 1; mode 0
!    means conv_water_4rad is never called (cloud_diagnostics.F90), so
!    only 1 and 2 are driven.
!  - zm_param%zm_microp (set directly on the real zm_conv zm_param):
!    selects the convective-microphysics partition branch. EAMv3
!    default zmconv_microp=.true.; BOTH branches are pure per-point
!    arithmetic here (dp_icimr is just a pbuf input), so both are
!    goldened, unlike zm_conv itself.
!  - microp_scheme ('P3' EAMv3 default / 'RK' legacy kabsi formula)
!    and pergro_mods (.false. default) via the phys_control stub.
!
! fice may contain NaN (the real FICE pbuf field can be unset); the
! shr_infnan stub's x/=x test reproduces the guarded COSP outputs
! sh_cldliq/sh_cldice, returned from the pbuf stub after the call.
module conv_water_driver
  implicit none
  ! local kind param (a use-associated kind is invisible to f2py's
  ! crackfortran; fbuild passes --f2cmap mapping r8 -> double)
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_conv_water(ncol, nlev, mode, zm_microp, is_rk, &
       pergro, t, pdel, q_cldliq, q_cldice, sh_icwmr, dp_icwmr, &
       dp_icimr, fice, sh_frac, dp_frac, ast, rei, &
       totg_liq, totg_ice, sh_cldliq, sh_cldice)
    use ppgrid,         only: pcols, pver
    use conv_water,     only: conv_water_register, conv_water_init, &
                              conv_water_4rad
    use zm_conv,        only: zm_param
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc, pbuf_get_index, &
                              pbuf_stub_fill, pbuf_stub_read

    use phys_control,   only: phys_control_stub_set

    integer, intent(in) :: ncol, nlev
    integer, intent(in) :: mode                      ! conv_water_mode
    integer, intent(in) :: zm_microp                 ! 0/1 flag
    integer, intent(in) :: is_rk                     ! 0/1 flag
    integer, intent(in) :: pergro                    ! 0/1 flag
    real(r8), intent(in) :: t(ncol, nlev)            ! [K]
    real(r8), intent(in) :: pdel(ncol, nlev)         ! [Pa]
    real(r8), intent(in) :: q_cldliq(ncol, nlev)     ! [kg/kg] gbx
    real(r8), intent(in) :: q_cldice(ncol, nlev)     ! [kg/kg] gbx
    real(r8), intent(in) :: sh_icwmr(ncol, nlev)     ! [kg/kg] in-cloud
    real(r8), intent(in) :: dp_icwmr(ncol, nlev)     ! [kg/kg] in-cloud
    real(r8), intent(in) :: dp_icimr(ncol, nlev)     ! [kg/kg] in-cloud
    real(r8), intent(in) :: fice(ncol, nlev)         ! may be NaN
    real(r8), intent(in) :: sh_frac(ncol, nlev)
    real(r8), intent(in) :: dp_frac(ncol, nlev)
    real(r8), intent(in) :: ast(ncol, nlev)
    real(r8), intent(in) :: rei(ncol, nlev)          ! [micron]

    real(r8), intent(out) :: totg_liq(ncol, nlev)    ! [kg/kg] gbx
    real(r8), intent(out) :: totg_ice(ncol, nlev)    ! [kg/kg] gbx
    real(r8), intent(out) :: sh_cldliq(ncol, nlev)   ! [kg/kg] gbx COSP
    real(r8), intent(out) :: sh_cldice(ncol, nlev)   ! [kg/kg] gbx COSP

    type(physics_state) :: state
    type(physics_buffer_desc), pointer :: pbuf(:), pbuf2d(:, :)
    real(r8) :: gl(pcols, pver), gi(pcols, pver)
    real(r8) :: buf(pcols, pver)

    if (ncol > pcols .or. nlev /= pver) then
      write(*, *) 'drv_conv_water: need ncol <= ', pcols, &
           ' and nlev == ', pver
      stop 2
    end if

    if (is_rk /= 0) then
      call phys_control_stub_set('RK', pergro /= 0)
    else
      call phys_control_stub_set('P3', pergro /= 0)
    end if
    zm_param%zm_microp = zm_microp /= 0

    ! real register + init (sets the module pbuf/constituent indices;
    ! init also zeroes the stub FICE field, so fill pbuf afterwards)
    allocate(pbuf(1), pbuf2d(1, 1))
    call conv_water_register()
    call conv_water_init(pbuf2d)

    call fill(pbuf_get_index('ICWMRSH'), sh_icwmr)
    call fill(pbuf_get_index('ICWMRDP'), dp_icwmr)
    call fill(pbuf_get_index('ICIMRDP'), dp_icimr)
    call fill(pbuf_get_index('FICE'),    fice)
    call fill(pbuf_get_index('SH_FRAC'), sh_frac)
    call fill(pbuf_get_index('DP_FRAC'), dp_frac)
    call fill(pbuf_get_index('AST'),     ast)
    call fill(pbuf_get_index('REI'),     rei)

    state%lchnk = 1
    state%ncol = ncol
    state%t = 0.0_r8
    state%pdel = 0.0_r8
    state%q = 0.0_r8
    state%t(1:ncol, :) = t
    state%pdel(1:ncol, :) = pdel
    state%q(1:ncol, :, 2) = q_cldliq
    state%q(1:ncol, :, 3) = q_cldice

    gl = 0.0_r8
    gi = 0.0_r8
    call conv_water_4rad(state, pbuf, mode, gl, gi)
    totg_liq = gl(1:ncol, :)
    totg_ice = gi(1:ncol, :)

    call pbuf_stub_read(pbuf_get_index('SH_CLDLIQ1'), buf)
    sh_cldliq = buf(1:ncol, :)
    call pbuf_stub_read(pbuf_get_index('SH_CLDICE1'), buf)
    sh_cldice = buf(1:ncol, :)

    deallocate(pbuf, pbuf2d)

  contains

    subroutine fill(idx, a)
      integer,  intent(in) :: idx
      real(r8), intent(in) :: a(ncol, nlev)
      buf = 0.0_r8
      buf(1:ncol, :) = a
      call pbuf_stub_fill(idx, buf)
    end subroutine fill

  end subroutine drv_conv_water

end module conv_water_driver
