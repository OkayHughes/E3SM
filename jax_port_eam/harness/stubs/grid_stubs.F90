! Grid/control stubs for column-physics extraction (infrastructure
! only; dims chosen for the harness drivers).

module ppgrid
  implicit none
  public
  integer, parameter :: pcols = 16
  integer, parameter :: psubcols = 1
  integer, parameter :: pver = 72
  integer, parameter :: pverp = 73
  ! single-chunk harness (real ppgrid sets these at runtime)
  integer, parameter :: begchunk = 1
  integer, parameter :: endchunk = 1
end module ppgrid

module phys_grid
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols
  implicit none
  ! settable column coordinates (radians) for stubs that must return
  ! real values (e.g. the tropopause climatology regrid); harness
  ! drivers set them via phys_grid_stub_set_cols.
  integer  :: stub_ncols = pcols
  real(r8) :: stub_rlat(pcols) = 0._r8
  real(r8) :: stub_rlon(pcols) = 0._r8
contains
  subroutine phys_grid_stub_set_cols(ncols, rlat, rlon)
    integer,  intent(in) :: ncols
    real(r8), intent(in) :: rlat(pcols), rlon(pcols)
    stub_ncols = ncols
    stub_rlat  = rlat
    stub_rlon  = rlon
  end subroutine phys_grid_stub_set_cols
  integer function get_lat_p(lchnk, i)
    integer, intent(in) :: lchnk, i
    get_lat_p = 0
  end function get_lat_p
  integer function get_lon_p(lchnk, i)
    integer, intent(in) :: lchnk, i
    get_lon_p = 0
  end function get_lon_p
  integer function get_ncols_p(lchnk)
    integer, intent(in) :: lchnk
    get_ncols_p = stub_ncols
  end function get_ncols_p
  subroutine get_rlat_all_p(lchnk, rlatdim, rlats)
    integer,  intent(in)  :: lchnk, rlatdim
    real(r8), intent(out) :: rlats(rlatdim)
    rlats = stub_rlat(1:rlatdim)
  end subroutine get_rlat_all_p
  subroutine get_rlon_all_p(lchnk, rlondim, rlons)
    integer,  intent(in)  :: lchnk, rlondim
    real(r8), intent(out) :: rlons(rlondim)
    rlons = stub_rlon(1:rlondim)
  end subroutine get_rlon_all_p
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
