! Infrastructure-only stubs for compiling the REAL
! eam/src/physics/cam/cloud_fraction.F90 into the zm_evap harness
! (METHODOLOGY.md: stubs never contain physics). cloud_fraction.F90 is
! compiled UNMODIFIED because zm_conv_evap needs its cldfrc_fice
! (real physics: the T-based ice/snow partition); everything else in
! that module (cldfrc_readnl/register/init and the RK cldfrc scheme)
! merely has to COMPILE and is never invoked by the harness, so the
! modules below satisfy compile-time interfaces only.
!
! This file is self-contained (it is used INSTEAD of grid_stubs.F90 /
! conv_water_stubs.F90, whose module names it reuses) so that ppgrid
! can carry pcols = 48: cldfrc_fice declares its arrays (pcols, pver)
! from ppgrid, and the zm_conv_evap driver passes 48-row padded arrays
! so the golden profiles (42 columns, 72 levels) fit in one call.

module ppgrid
  implicit none
  public
  integer, parameter :: pcols = 48
  integer, parameter :: psubcols = 1
  integer, parameter :: pver = 72
  integer, parameter :: pverp = 73
  integer, parameter :: begchunk = 1
  integer, parameter :: endchunk = 1
end module ppgrid

module phys_grid
  ! only referenced by the never-called cldfrc scheme
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
contains
  subroutine get_rlat_all_p(lchnk, rlatdim, rlats)
    integer,  intent(in)  :: lchnk, rlatdim
    real(r8), intent(out) :: rlats(rlatdim)
    rlats = 0._r8
  end subroutine get_rlat_all_p
  subroutine get_rlon_all_p(lchnk, rlondim, rlons)
    integer,  intent(in)  :: lchnk, rlondim
    real(r8), intent(out) :: rlons(rlondim)
    rlons = 0._r8
  end subroutine get_rlon_all_p
end module phys_grid

module dycore
  implicit none
contains
  logical function dycore_is(name)
    character(len=*), intent(in) :: name
    ! E3SM uses the SE dycore; 'LR' (finite volume) is false
    dycore_is = .false.
  end function dycore_is
  character(len=7) function get_resolution()
    get_resolution = 'ne30np4'
  end function get_resolution
end module dycore

module ref_pres
  ! reference-pressure metadata (grid plumbing, no physics); only
  ! referenced by cldfrc_init, which the harness never calls
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver
  implicit none
  public
  real(r8) :: pref_mid(pver) = 0._r8
  integer  :: trop_cloud_top_lev = 1
end module ref_pres

module phys_control
  ! keyword names must match the real phys_control.F90 phys_getopts;
  ! only referenced by cldfrc_init (never called by the harness)
  implicit none
contains
  subroutine phys_getopts(shallow_scheme_out, eddy_scheme_out, &
                          macrop_scheme_out)
    character(len=*), intent(out), optional :: shallow_scheme_out
    character(len=*), intent(out), optional :: eddy_scheme_out
    character(len=*), intent(out), optional :: macrop_scheme_out
    if (present(shallow_scheme_out)) shallow_scheme_out = 'CLUBB_SGS'
    if (present(eddy_scheme_out))    eddy_scheme_out    = 'CLUBB_SGS'
    if (present(macrop_scheme_out))  macrop_scheme_out  = 'CLUBB_SGS'
  end subroutine phys_getopts
end module phys_control

module cam_history
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private
  public :: addfld, outfld
contains
  subroutine addfld(fname, dimnames, avgflag, units, long_name)
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: dimnames(:)
    character(len=*), intent(in) :: avgflag, units, long_name
  end subroutine addfld
  subroutine outfld(fname, field, idim, lchnk)
    character(len=*), intent(in) :: fname
    integer,          intent(in) :: idim, lchnk
    real(r8),         intent(in) :: field(idim, *)
  end subroutine outfld
end module cam_history

module physics_buffer
  ! minimal pbuf plumbing for cldfrc_register / cldfrc (never called)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver
  implicit none
  private
  save
  public :: physics_buffer_desc, dtype_r8
  public :: pbuf_add_field, pbuf_get_field

  integer, parameter :: dtype_r8 = 1

  type physics_buffer_desc
    integer :: dummy = 0
  end type physics_buffer_desc

  real(r8), target :: pbuf_dummy(pcols, pver) = 0._r8

contains

  subroutine pbuf_add_field(name, scope, dtype, dims, index)
    character(len=*), intent(in)  :: name, scope
    integer,          intent(in)  :: dtype
    integer,          intent(in)  :: dims(:)
    integer,          intent(out) :: index
    index = 1
  end subroutine pbuf_add_field

  subroutine pbuf_get_field(pbuf, idx, field)
    type(physics_buffer_desc), pointer :: pbuf(:)
    integer,  intent(in) :: idx
    real(r8), pointer    :: field(:, :)
    field => pbuf_dummy
  end subroutine pbuf_get_field

end module physics_buffer
