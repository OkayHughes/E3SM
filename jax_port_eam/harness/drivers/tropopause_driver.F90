! f2py driver for the EAM tropopause finder
! (eam/src/physics/cam/tropopause.F90). EAMv3 production invokes only
!   tropopause_find(state, tropLev)                              [TWMO+CLIMATE
!     default: aer_rad_props, tropopause_output]
!   tropopause_find(..., primary=TROP_ALG_TWMO,    backup=TROP_ALG_CLIMATE)
!     [prescribed_volcaero]
!   tropopause_find(..., primary=TROP_ALG_HYBSTOB, backup=TROP_ALG_CLIMATE)
!     [mozart chemistry, modal_aero_wateruptake]
! so this driver exposes tropopause_find with caller-selected
! primary/backup enums (NONE=1 disables the backup).
!
! The climatology enters through drv_init: the pio stub
! (tropopause_stubs.F90) serves a synthetic lon/lat/trop_p "file" and
! the REAL tropopause_init -> tropopause_read_file path fills the
! module-private tropp_p_loc/days using the REAL interpolate_data
! regridder. Column lats are placed exactly on the climatology lat
! nodes (and lons on lon nodes) so the regrid reproduces node values
! exactly (lininterp weights become exact {0,1}); the golden generator
! relies on that to know the per-column climatology bit-for-bit.
! days(:) = noleap day-of-year of mmdd codes 116..1216
!         = 16,45,75,105,136,166,197,228,258,289,319,350.
!
! The physics_state coupling becomes plain arrays: the driver fills a
! stub physics_state (pcols=16, pver=72; ppgrid in grid_stubs.F90) and
! chunks arbitrary column counts through it. lchnk is always 1, so all
! chunks share the same 16 climatology column slots (global column i
! uses slot mod(i-1,16)+1; the generator mirrors this).
module tropopause_driver
  implicit none
  ! local kind param (a use-associated kind is invisible to f2py's
  ! crackfortran; fbuild passes --f2cmap mapping r8 -> double)
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_init(nlon, nlat, lon_deg, lat_deg, tropp, &
                      col_lat_rad, col_lon_rad)
    use pio,        only: pio_stub_set_climo
    use phys_grid,  only: phys_grid_stub_set_cols
    use tropopause, only: tropopause_init
    integer,  intent(in) :: nlon, nlat
    real(r8), intent(in) :: lon_deg(nlon), lat_deg(nlat)
    real(r8), intent(in) :: tropp(nlon, nlat, 12)
    real(r8), intent(in) :: col_lat_rad(16), col_lon_rad(16)
    call pio_stub_set_climo(nlon, nlat, lon_deg, lat_deg, tropp)
    call phys_grid_stub_set_cols(16, col_lat_rad, col_lon_rad)
    call tropopause_init()
  end subroutine drv_init

  subroutine drv_find(ncol, calday, primary, backup, lat, t, pmid, &
                      pint, zm, zi, troplev, tropp, tropt, tropz)
    use physics_types, only: physics_state
    use time_manager,  only: tm_stub_set_calday
    use tropopause,    only: tropopause_find
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: calday
    integer,  intent(in)  :: primary, backup   ! TROP_ALG_* enums
    real(r8), intent(in)  :: lat(16)           ! radians
    real(r8), intent(in)  :: t(16, 72), pmid(16, 72), zm(16, 72)
    real(r8), intent(in)  :: pint(16, 73), zi(16, 73)
    integer,  intent(out) :: troplev(16)
    real(r8), intent(out) :: tropp(16), tropt(16), tropz(16)
    type(physics_state) :: st
    st%lchnk = 1
    st%ncol  = ncol
    st%lat   = lat
    st%t     = t
    st%pmid  = pmid
    st%pint  = pint
    st%zm    = zm
    st%zi    = zi
    st%q     = 0._r8
    call tm_stub_set_calday(calday)
    call tropopause_find(st, troplev, tropP=tropp, tropT=tropt, &
                         tropZ=tropz, primary=primary, backup=backup)
  end subroutine drv_find

end module tropopause_driver
