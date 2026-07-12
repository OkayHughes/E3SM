! Infrastructure-only stubs needed to compile the UNMODIFIED
! eam/src/physics/cam/tropopause.F90 standalone (METHODOLOGY.md rules:
! stubs never contain physics).
!
! Runtime-relevant stubs:
!  - physics_types: plain data container mirroring the physics_state
!    fields tropopause.F90 reads (lchnk/ncol/lat/t/pmid/pint/zm/zi/q).
!  - pio + cam_pio_utils + ioFileMod: a fake climatology "file". The
!    REAL tropopause_read_file runs against it (with the REAL
!    interpolate_data.F90 regridder), so the module-private
!    tropp_p_loc/days get filled through the production code path.
!    The driver injects lon/lat/tropp via pio_stub_set_climo.
!  - time_manager: get_calday reproduces the noleap day-of-year that
!    ESMF returns for the mmdd codes tropopause_read_file passes
!    (116 -> 16.0, ..., 1216 -> 350.0); get_curr_calday returns a
!    driver-settable value (tm_stub_set_calday).
!  - cam_history: no-op addfld/add_default/outfld (tropopause_hybridstobie
!    calls outfld at runtime; output is not part of the golden).
!  - cam_history_support: fillvalue = 1.e+20_r8 VERBATIM from
!    eam/src/control/cam_history_support.F90.
!  - mo_chem_utls: get_spc_ndx -> -1 (no E90 tracer; disables the E90
!    paths exactly as in an EAMv3 run without E90).
!
! Compile-only stubs, NEVER executed (guarded by get_spc_ndx() > 0 or
! masterproc, both false here); abort-only placeholders per
! METHODOLOGY.md:
!  - constituents (cnst_get_ind), chem_mods (adv_mass), dyn_grid
!    (get_dyn_grid_parm is called but its plon/plat results are unused
!    by tropopause_read_file).

module cam_history_support
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  ! VERBATIM from eam/src/control/cam_history_support.F90
  real(r8), parameter, public :: fillvalue = 1.e+20_r8
end module cam_history_support

module cam_history
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  character(len=16), parameter :: horiz_only = 'horiz_only'
  interface addfld
     module procedure addfld_hor
     module procedure addfld_nd
  end interface addfld
  interface outfld
     module procedure outfld_1d
     module procedure outfld_2d
  end interface outfld
contains
  subroutine addfld_hor(fname, dims, avgflag, units, long_name, flag_xyfill)
    character(len=*), intent(in) :: fname, dims, avgflag, units, long_name
    logical, intent(in), optional :: flag_xyfill
  end subroutine addfld_hor
  subroutine addfld_nd(fname, dimnames, avgflag, units, long_name, flag_xyfill)
    character(len=*), intent(in) :: fname, avgflag, units, long_name
    character(len=*), intent(in) :: dimnames(:)
    logical, intent(in), optional :: flag_xyfill
  end subroutine addfld_nd
  subroutine add_default(name, tindex, flag)
    character(len=*), intent(in) :: name, flag
    integer, intent(in) :: tindex
  end subroutine add_default
  subroutine outfld_1d(fname, field, idim, c)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: idim, c
    real(r8), intent(in) :: field(idim)
  end subroutine outfld_1d
  subroutine outfld_2d(fname, field, idim, c)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: idim, c
    real(r8), intent(in) :: field(idim, *)
  end subroutine outfld_2d
end module cam_history

module physics_types
  ! plain-data mirror of the physics_state fields tropopause.F90 uses
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver, pverp
  implicit none
  type physics_state
     integer  :: lchnk = 1
     integer  :: ncol = pcols
     real(r8) :: lat(pcols) = 0._r8
     real(r8) :: t(pcols, pver) = 0._r8
     real(r8) :: pmid(pcols, pver) = 0._r8
     real(r8) :: pint(pcols, pverp) = 0._r8
     real(r8) :: zm(pcols, pver) = 0._r8
     real(r8) :: zi(pcols, pverp) = 0._r8
     real(r8) :: q(pcols, pver, 1) = 0._r8
  end type physics_state
end module physics_types

module time_manager
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  real(r8) :: stub_calday = 1._r8
contains
  subroutine tm_stub_set_calday(calday)
    real(r8), intent(in) :: calday
    stub_calday = calday
  end subroutine tm_stub_set_calday
  real(r8) function get_curr_calday()
    get_curr_calday = stub_calday
  end function get_curr_calday
  real(r8) function get_calday(ymd, tod)
    ! noleap day-of-year for mmdd codes (Jan 1 00Z -> 1.0), matching
    ! the ESMF NO_LEAP dayOfYear_r8 the real time_manager returns
    integer, intent(in) :: ymd, tod
    integer, parameter :: doy0(12) = &
         (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
    integer :: month, day
    month = ymd / 100
    day   = mod(ymd, 100)
    get_calday = real(doy0(month) + day, r8) + real(tod, r8) / 86400._r8
  end function get_calday
end module time_manager

module mo_chem_utls
  implicit none
contains
  integer function get_spc_ndx(spc_name)
    character(len=*), intent(in) :: spc_name
    get_spc_ndx = -1   ! no chemistry tracers in the harness
  end function get_spc_ndx
end module mo_chem_utls

module constituents
  implicit none
contains
  subroutine cnst_get_ind(name, ind, abrtf)
    ! compile-only: unreachable (guarded by get_spc_ndx('E90') > 0)
    use cam_abortutils, only: endrun
    character(len=*), intent(in) :: name
    integer, intent(out) :: ind
    logical, intent(in), optional :: abrtf
    ind = -1
    call endrun('cnst_get_ind stub called')
  end subroutine cnst_get_ind
end module constituents

module chem_mods
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  ! compile-only: adv_mass is read solely under get_spc_ndx('E90') > 0
  real(r8) :: adv_mass(1) = 0._r8
end module chem_mods

module dyn_grid
  implicit none
contains
  integer function get_dyn_grid_parm(name)
    ! called by tropopause_read_file, but its plon/plat results are
    ! never used there
    character(len=*), intent(in) :: name
    get_dyn_grid_parm = 0
  end function get_dyn_grid_parm
end module dyn_grid

module ioFileMod
  implicit none
contains
  subroutine getfil(fulpath, locfn, iflag)
    character(len=*), intent(in)  :: fulpath
    character(len=*), intent(out) :: locfn
    integer, intent(in) :: iflag
    locfn = fulpath
  end subroutine getfil
end module ioFileMod

module pio
  ! fake netcdf climatology "file": serves lon/lat (degrees) and
  ! trop_p(nlon,nlat,12) injected by the driver via pio_stub_set_climo
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  integer, parameter :: pio_nowrite = 0
  type file_desc_t
     integer :: fh = 0
  end type file_desc_t
  type var_desc_t
     integer :: id = 0
  end type var_desc_t
  integer :: stub_nlon = 0
  integer :: stub_nlat = 0
  real(r8), allocatable :: stub_lon(:)        ! degrees
  real(r8), allocatable :: stub_lat(:)        ! degrees
  real(r8), allocatable :: stub_tropp(:,:,:)  ! (nlon,nlat,12) Pa
  interface pio_get_var
     module procedure pio_get_var_1d
     module procedure pio_get_var_3d
  end interface pio_get_var
contains
  subroutine pio_stub_set_climo(nlon, nlat, lon_deg, lat_deg, tropp)
    integer, intent(in) :: nlon, nlat
    real(r8), intent(in) :: lon_deg(nlon), lat_deg(nlat)
    real(r8), intent(in) :: tropp(nlon, nlat, 12)
    if (allocated(stub_lon)) deallocate(stub_lon, stub_lat, stub_tropp)
    stub_nlon = nlon
    stub_nlat = nlat
    allocate(stub_lon(nlon), stub_lat(nlat), stub_tropp(nlon, nlat, 12))
    stub_lon = lon_deg
    stub_lat = lat_deg
    stub_tropp = tropp
  end subroutine pio_stub_set_climo
  integer function pio_inq_dimid(file, name, dimid)
    type(file_desc_t), intent(in) :: file
    character(len=*), intent(in) :: name
    integer, intent(out) :: dimid
    select case (trim(name))
    case ('time')
       dimid = 1
    case ('lat')
       dimid = 2
    case ('lon')
       dimid = 3
    case default
       dimid = -1
    end select
    pio_inq_dimid = 0
  end function pio_inq_dimid
  integer function pio_inq_dimlen(file, dimid, dimlen)
    type(file_desc_t), intent(in) :: file
    integer, intent(in) :: dimid
    integer, intent(out) :: dimlen
    select case (dimid)
    case (1)
       dimlen = 12
    case (2)
       dimlen = stub_nlat
    case (3)
       dimlen = stub_nlon
    case default
       dimlen = -1
    end select
    pio_inq_dimlen = 0
  end function pio_inq_dimlen
  integer function pio_inq_varid(file, name, vid)
    type(file_desc_t), intent(in) :: file
    character(len=*), intent(in) :: name
    type(var_desc_t), intent(out) :: vid
    select case (trim(name))
    case ('lat')
       vid%id = 2
    case ('lon')
       vid%id = 3
    case ('trop_p')
       vid%id = 4
    case default
       vid%id = -1
    end select
    pio_inq_varid = 0
  end function pio_inq_varid
  integer function pio_get_var_1d(file, vid, arr)
    type(file_desc_t), intent(in) :: file
    type(var_desc_t), intent(in) :: vid
    real(r8), intent(out) :: arr(:)
    use_id: select case (vid%id)
    case (2)
       arr = stub_lat
    case (3)
       arr = stub_lon
    end select use_id
    pio_get_var_1d = 0
  end function pio_get_var_1d
  integer function pio_get_var_3d(file, vid, start, count, arr)
    type(file_desc_t), intent(in) :: file
    type(var_desc_t), intent(in) :: vid
    integer, intent(in) :: start(:), count(:)
    real(r8), intent(out) :: arr(:,:,:)
    arr = stub_tropp(start(1):start(1)+count(1)-1, &
                     start(2):start(2)+count(2)-1, &
                     start(3):start(3)+count(3)-1)
    pio_get_var_3d = 0
  end function pio_get_var_3d
  subroutine pio_closefile(file)
    type(file_desc_t), intent(in) :: file
  end subroutine pio_closefile
end module pio

module cam_pio_utils
  use pio, only: file_desc_t
  implicit none
contains
  subroutine cam_pio_openfile(file, fname, mode)
    type(file_desc_t), intent(inout) :: file
    character(len=*), intent(in) :: fname
    integer, intent(in) :: mode
    file%fh = 1
  end subroutine cam_pio_openfile
end module cam_pio_utils
