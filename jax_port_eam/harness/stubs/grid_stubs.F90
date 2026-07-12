! Grid/control stubs for column-physics extraction (infrastructure
! only; dims chosen for the harness drivers).

module ppgrid
  implicit none
  public
  integer, parameter :: pcols = 16
  integer, parameter :: psubcols = 1
  integer, parameter :: pver = 72
  integer, parameter :: pverp = 73
end module ppgrid

module phys_grid
  implicit none
contains
  integer function get_lat_p(lchnk, i)
    integer, intent(in) :: lchnk, i
    get_lat_p = 0
  end function get_lat_p
  integer function get_lon_p(lchnk, i)
    integer, intent(in) :: lchnk, i
    get_lon_p = 0
  end function get_lon_p
end module phys_grid

module cam_control_mod
  implicit none
  public
  integer :: nlvdry = 3
end module cam_control_mod

module dycore
  implicit none
contains
  logical function dycore_is(name)
    character(len=*), intent(in) :: name
    ! E3SM uses the SE dycore; 'LR' (finite volume) is false. The
    ! geopotential harness exercises only this (non-FV) branch, matching
    ! all E3SM production configurations.
    dycore_is = .false.
  end function dycore_is
end module dycore
