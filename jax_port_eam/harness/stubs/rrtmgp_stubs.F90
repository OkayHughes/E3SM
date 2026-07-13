! Infrastructure-only stubs for the EAM RRTMGP radiation harness
! (METHODOLOGY.md rules: never physics). This build does NOT use
! stubs/infrastructure_stubs.F90 because cloud_rad_props_init and
! mo_load_coefficients read netCDF only on masterproc, so this harness
! needs masterproc = .true. (non-SPMD, single process: every read
! happens in-process; no broadcast required).
!
! Abort-only stubs (documented, never executed by the harness):
!   rad_constituents rad_cnst_get_gas   (gas vmr is a plain input)
!   aer_rad_props aer_rad_props_sw/lw   (aerosol optics are plain inputs)
!   physics_buffer accessors            (no pbuf in the harness)
!   time_manager get_curr_date/calday   (only referenced by the
!                                        never-called MACv2 branch)
!   cam_history addfld/outfld           (history plumbing; outfld is a
!                                        no-op because cam_optics calls
!                                        it inside set_aerosol_optics_sw
!                                        which the harness never calls)

module cam_logfile
  implicit none
  integer, parameter :: iulog = 6
end module cam_logfile

module infnan
  ! infrastructure stub for control/infnan.F90 (which wraps the
  ! genf90-templated shr_infnan_mod): IEEE nan/inf tests used only by
  ! assertions.F90 assert_valid diagnostics.
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
contains
  elemental logical function isnan(x)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    real(r8), intent(in) :: x
    isnan = ieee_is_nan(x)
  end function isnan
  elemental logical function isinf(x)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    real(r8), intent(in) :: x
    isinf = (.not. ieee_is_finite(x)) .and. (.not. ieee_is_nan(x))
  end function isinf
end module infnan

module cam_abortutils
  implicit none
contains
  subroutine endrun(msg)
    character(len=*), intent(in), optional :: msg
    if (present(msg)) print *, 'ENDRUN: ', msg
    stop 1
  end subroutine endrun
end module cam_abortutils

module spmd_utils
  implicit none
  logical, parameter :: masterproc = .true.
end module spmd_utils

module perf_mod
  implicit none
contains
  subroutine t_startf(name)
    character(len=*), intent(in) :: name
  end subroutine t_startf
  subroutine t_stopf(name)
    character(len=*), intent(in) :: name
  end subroutine t_stopf
end module perf_mod

module ioFileMod
  implicit none
contains
  ! getfil: in the model this resolves a file against the inputdata
  ! area; here paths are passed in fully resolved.
  subroutine getfil(fil, locfn, iflag)
    character(len=*), intent(in)  :: fil
    character(len=*), intent(out) :: locfn
    integer, optional, intent(in) :: iflag
    locfn = trim(fil)
  end subroutine getfil
end module ioFileMod

module error_messages
  implicit none
contains
  subroutine handle_ncerr(ret, mes, line)
    ! real behavior: abort when a netCDF call fails (NF90_NOERR = 0)
    integer, intent(in) :: ret
    character(len=*), intent(in) :: mes
    integer, intent(in), optional :: line
    if (ret /= 0) then
      print *, 'NETCDF ERROR: ', mes, ' ret=', ret
      stop 1
    end if
  end subroutine handle_ncerr
end module error_messages

module physics_types
  ! minimal data container consumed by radiation_state set_rad_state
  ! (t, pmid, pint, lnpmid, lnpint, ncol) and referenced (never
  ! executed) by cam_optics' aerosol subroutines (lchnk, phis, zm).
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver, pverp
  implicit none
  type physics_state
    integer :: lchnk = 1
    integer :: ncol = 0
    real(r8) :: t(pcols, pver) = 0._r8
    real(r8) :: pmid(pcols, pver) = 0._r8
    real(r8) :: pint(pcols, pverp) = 0._r8
    real(r8) :: lnpmid(pcols, pver) = 0._r8
    real(r8) :: lnpint(pcols, pverp) = 0._r8
    real(r8) :: phis(pcols) = 0._r8
    real(r8) :: zm(pcols, pver) = 0._r8
  end type physics_state
end module physics_types

module camsrfexch
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols
  implicit none
  type cam_in_t
    real(r8) :: lwup(pcols)  = 0._r8
    real(r8) :: asdir(pcols) = 0._r8
    real(r8) :: asdif(pcols) = 0._r8
    real(r8) :: aldir(pcols) = 0._r8
    real(r8) :: aldif(pcols) = 0._r8
  end type cam_in_t
end module camsrfexch

module physics_buffer
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  type physics_buffer_desc
    integer :: dummy = 0
  end type physics_buffer_desc
contains
  integer function pbuf_get_index(name, errcode) result(idx)
    character(len=*), intent(in) :: name
    integer, optional, intent(out) :: errcode
    print *, 'rrtmgp_stubs pbuf_get_index called: ', name
    stop 1
    idx = -1
  end function pbuf_get_index
  subroutine pbuf_get_field(pbuf, index, field)
    type(physics_buffer_desc), pointer :: pbuf(:)
    integer, intent(in) :: index
    real(r8), pointer :: field(:,:)
    print *, 'rrtmgp_stubs pbuf_get_field called'
    stop 1
  end subroutine pbuf_get_field
  integer function pbuf_old_tim_idx() result(idx)
    print *, 'rrtmgp_stubs pbuf_old_tim_idx called'
    stop 1
    idx = -1
  end function pbuf_old_tim_idx
end module physics_buffer

module rad_constituents
  ! In the model these come from the rad_climate namelist parsing.
  ! Here they are plain public module variables set by the driver
  ! before *_init calls (EAMv3 defaults: gammadist / mitchell and the
  ! atm/cam/physprops optics files).
  implicit none
  character(len=32)  :: icecldoptics = 'mitchell'
  character(len=32)  :: liqcldoptics = 'gammadist'
  character(len=256) :: iceopticsfile = ' '
  character(len=256) :: liqopticsfile = ' '
contains
  subroutine rad_cnst_get_gas(icall, gasname, state, pbuf, mmr)
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    integer, intent(in) :: icall
    character(len=*), intent(in) :: gasname
    type(physics_state), intent(in) :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), pointer :: mmr(:,:)
    print *, 'rrtmgp_stubs rad_cnst_get_gas called: ', gasname
    stop 1
  end subroutine rad_cnst_get_gas
end module rad_constituents

module aer_rad_props
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
contains
  subroutine aer_rad_props_sw(icall, dt, state, pbuf, nnite, idxnite, &
                              is_cmip6_volc, tau, tau_w, tau_w_g, tau_w_f, &
                              clear_rh)
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    integer, intent(in) :: icall
    real(r8), intent(in) :: dt
    type(physics_state), intent(in), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    integer, intent(in) :: nnite, idxnite(:)
    logical, intent(in) :: is_cmip6_volc
    real(r8), intent(out) :: tau(:,0:,:), tau_w(:,0:,:), tau_w_g(:,0:,:), &
                             tau_w_f(:,0:,:)
    real(r8), optional, intent(in) :: clear_rh(:,:)
    print *, 'rrtmgp_stubs aer_rad_props_sw called'
    stop 1
  end subroutine aer_rad_props_sw
  subroutine aer_rad_props_lw(is_cmip6_volc, icall, dt, state, pbuf, odap_aer)
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    logical, intent(in) :: is_cmip6_volc
    integer, intent(in) :: icall
    real(r8), intent(in) :: dt
    type(physics_state), intent(in), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(out) :: odap_aer(:,:,:)
    print *, 'rrtmgp_stubs aer_rad_props_lw called'
    stop 1
  end subroutine aer_rad_props_lw
end module aer_rad_props

module time_manager
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
contains
  subroutine get_curr_date(yr, mon, day, tod, offset)
    integer, intent(out) :: yr, mon, day, tod
    integer, optional, intent(in) :: offset
    print *, 'rrtmgp_stubs get_curr_date called'
    stop 1
    yr = 0; mon = 0; day = 0; tod = 0
  end subroutine get_curr_date
  function get_curr_calday(offset) result(calday)
    use shr_kind_mod, only: r8 => shr_kind_r8
    integer, optional, intent(in) :: offset
    real(r8) :: calday
    print *, 'rrtmgp_stubs get_curr_calday called'
    stop 1
    calday = 0._r8
  end function get_curr_calday
end module time_manager

module cam_history
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  interface outfld
    module procedure outfld_1d, outfld_2d
  end interface outfld
contains
  subroutine addfld()
    ! name-only import target (never called in compiled code paths)
  end subroutine addfld
  subroutine outfld_1d(fname, field, idim, lchnk)
    character(len=*), intent(in) :: fname
    real(r8), intent(in) :: field(:)
    integer, intent(in) :: idim, lchnk
  end subroutine outfld_1d
  subroutine outfld_2d(fname, field, idim, lchnk)
    character(len=*), intent(in) :: fname
    real(r8), intent(in) :: field(:,:)
    integer, intent(in) :: idim, lchnk
  end subroutine outfld_2d
end module cam_history
