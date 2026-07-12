! Infrastructure-only stubs needed to compile
! eam/src/physics/cam/conv_water.F90 unmodified (METHODOLOGY.md: stubs
! never contain physics). One file, several tiny modules:
!  - perf_mod: timers (empty; conv_water only use-associates it)
!  - shr_infnan_mod: shr_infnan_isnan (IEEE x/=x test; the real module
!    is a genf90 template share/util/shr_infnan_mod.F90.in)
!  - phys_control: phys_getopts keyword plumbing; values injected by
!    the harness driver (the real module just returns namelist reads)
!  - cam_history: addfld/outfld no-ops
!  - physics_types: minimal physics_state container (the subset
!    conv_water_4rad reads: lchnk, ncol, t, pdel, q with vapor=1,
!    CLDLIQ=2, CLDICE=3 as returned by the constituents stub)
!  - physics_buffer: fixed-index registry backed by one module array;
!    the driver fills fields with pbuf_stub_fill and reads them back
!    with pbuf_stub_read

module perf_mod
  implicit none
contains
  subroutine t_startf(event)
    character(len=*), intent(in) :: event
  end subroutine t_startf
  subroutine t_stopf(event)
    character(len=*), intent(in) :: event
  end subroutine t_stopf
end module perf_mod

module shr_infnan_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private
  public :: shr_infnan_isnan
  interface shr_infnan_isnan
    module procedure shr_infnan_isnan_r8
  end interface
contains
  elemental logical function shr_infnan_isnan_r8(x)
    real(r8), intent(in) :: x
    shr_infnan_isnan_r8 = (x /= x)
  end function shr_infnan_isnan_r8
end module shr_infnan_mod

module phys_control
  implicit none
  private
  save
  public :: phys_getopts, phys_control_stub_set

  character(len=16) :: stub_microp_scheme = 'P3'
  logical :: stub_pergro_mods = .false.
  logical :: stub_use_mmf = .false.

contains

  subroutine phys_control_stub_set(microp_scheme, pergro_mods)
    character(len=*), intent(in) :: microp_scheme
    logical,          intent(in) :: pergro_mods
    stub_microp_scheme = microp_scheme
    stub_pergro_mods   = pergro_mods
  end subroutine phys_control_stub_set

  ! keyword names must match the real phys_control.F90 phys_getopts
  subroutine phys_getopts(microp_scheme_out, pergro_mods_out, &
                          use_MMF_out)
    character(len=*), intent(out), optional :: microp_scheme_out
    logical,          intent(out), optional :: pergro_mods_out
    logical,          intent(out), optional :: use_MMF_out
    if (present(microp_scheme_out)) microp_scheme_out = stub_microp_scheme
    if (present(pergro_mods_out))   pergro_mods_out   = stub_pergro_mods
    if (present(use_MMF_out))       use_MMF_out       = stub_use_mmf
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

module physics_types
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver
  implicit none
  private
  public :: physics_state

  type physics_state
    integer :: lchnk = 1
    integer :: ncol = pcols
    real(r8) :: t(pcols, pver) = 0.0_r8
    real(r8) :: pdel(pcols, pver) = 0.0_r8
    real(r8) :: q(pcols, pver, 3) = 0.0_r8
  end type physics_state
end module physics_types

module physics_buffer
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver
  implicit none
  private
  save

  public :: physics_buffer_desc, dtype_r8
  public :: pbuf_add_field, pbuf_get_index, pbuf_get_field
  public :: pbuf_set_field, pbuf_old_tim_idx
  public :: pbuf_stub_fill, pbuf_stub_read

  integer, parameter :: dtype_r8 = 1

  type physics_buffer_desc
    integer :: dummy = 0
  end type physics_buffer_desc

  ! fixed name -> index registry for the fields conv_water touches
  integer, parameter :: nflds = 10
  real(r8), target :: pbuf_data(pcols, pver, nflds) = 0.0_r8

  interface pbuf_get_field
    module procedure pbuf_get_field_2d
  end interface
  interface pbuf_set_field
    module procedure pbuf_set_field_scalar
  end interface

contains

  integer function name_to_index(name) result(idx)
    use cam_abortutils, only: endrun
    character(len=*), intent(in) :: name
    select case (trim(name))
    case ('ICWMRSH');    idx = 1
    case ('ICWMRDP');    idx = 2
    case ('ICIMRDP');    idx = 3
    case ('FICE');       idx = 4
    case ('SH_FRAC');    idx = 5
    case ('DP_FRAC');    idx = 6
    case ('AST');        idx = 7
    case ('REI');        idx = 8
    case ('SH_CLDLIQ1'); idx = 9
    case ('SH_CLDICE1'); idx = 10
    case default
      idx = -1
      call endrun('pbuf stub: unknown field '//trim(name))
    end select
  end function name_to_index

  integer function pbuf_get_index(name) result(idx)
    character(len=*), intent(in) :: name
    idx = name_to_index(name)
  end function pbuf_get_index

  subroutine pbuf_add_field(name, scope, dtype, dims, index)
    character(len=*), intent(in)  :: name, scope
    integer,          intent(in)  :: dtype
    integer,          intent(in)  :: dims(:)
    integer,          intent(out) :: index
    index = name_to_index(name)
  end subroutine pbuf_add_field

  integer function pbuf_old_tim_idx()
    pbuf_old_tim_idx = 1
  end function pbuf_old_tim_idx

  subroutine pbuf_set_field_scalar(pbuf2d, idx, value)
    type(physics_buffer_desc), pointer :: pbuf2d(:, :)
    integer,  intent(in) :: idx
    real(r8), intent(in) :: value
    pbuf_data(:, :, idx) = value
  end subroutine pbuf_set_field_scalar

  ! start/kount accepted (conv_water passes them for the time-indexed
  ! AST field) but ignored: the stub stores a single time level
  subroutine pbuf_get_field_2d(pbuf, idx, field, start, kount)
    type(physics_buffer_desc), pointer :: pbuf(:)
    integer,  intent(in) :: idx
    real(r8), pointer    :: field(:, :)
    integer,  intent(in), optional :: start(:), kount(:)
    field => pbuf_data(:, :, idx)
  end subroutine pbuf_get_field_2d

  subroutine pbuf_stub_fill(idx, field)
    integer,  intent(in) :: idx
    real(r8), intent(in) :: field(pcols, pver)
    pbuf_data(:, :, idx) = field
  end subroutine pbuf_stub_fill

  subroutine pbuf_stub_read(idx, field)
    integer,  intent(in)  :: idx
    real(r8), intent(out) :: field(pcols, pver)
    field = pbuf_data(:, :, idx)
  end subroutine pbuf_stub_read

end module physics_buffer
