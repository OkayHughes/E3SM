! Infrastructure-only stubs for compiling EAM parameterizations
! standalone (see METHODOLOGY.md: stubs never contain physics).
! One file, many tiny modules.

module cam_logfile
  implicit none
  integer, parameter :: iulog = 6
end module cam_logfile

module cam_abortutils
  implicit none
contains
  subroutine endrun(msg)
    character(len=*), intent(in), optional :: msg
    if (present(msg)) write(*, *) 'ENDRUN: ', msg
    stop 1
  end subroutine endrun
end module cam_abortutils

module spmd_utils
  implicit none
  logical, parameter :: masterproc = .false.
end module spmd_utils

module mpishorthand
  implicit none
  integer :: mpicom = 0
  integer :: mpichar = 0
  integer :: mpilog = 0
  integer :: mpiint = 0
  integer :: mpir8 = 0
end module mpishorthand

module namelist_utils
  implicit none
contains
  subroutine find_group_name(unit, group, status)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: group
    integer, intent(out) :: status
    status = 1  ! group never found; callers then skip the read
  end subroutine find_group_name
end module namelist_utils

module units
  implicit none
contains
  integer function getunit(iu)
    integer, intent(in), optional :: iu
    getunit = 99
    if (present(iu)) getunit = iu
  end function getunit
  subroutine freeunit(iu)
    integer, intent(in) :: iu
  end subroutine freeunit
end module units

module shr_assert_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
contains
  subroutine shr_assert_in_domain(var, ge, le, gt, lt, varname, msg)
    use cam_abortutils, only: endrun
    real(r8), intent(in) :: var
    real(r8), intent(in), optional :: ge, le, gt, lt
    character(len=*), intent(in), optional :: varname, msg
    logical :: ok
    ok = .true.
    if (present(ge)) ok = ok .and. (var >= ge)
    if (present(le)) ok = ok .and. (var <= le)
    if (present(gt)) ok = ok .and. (var > gt)
    if (present(lt)) ok = ok .and. (var < lt)
    if (.not. ok) call endrun(msg)
  end subroutine shr_assert_in_domain
end module shr_assert_mod

module error_messages
  implicit none
contains
  subroutine handle_errmsg(errmsg, subname, extra_msg)
    use cam_abortutils, only: endrun
    character(len=*), intent(in) :: errmsg
    character(len=*), intent(in), optional :: subname, extra_msg
    if (trim(errmsg) /= '') then
      if (present(subname)) write(*, *) 'ERROR in ', subname
      call endrun(errmsg)
    end if
  end subroutine handle_errmsg
end module error_messages
