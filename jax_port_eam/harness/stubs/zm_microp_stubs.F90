! Infrastructure-only stubs for compiling the REAL
! eam/src/physics/cam/zm/zm_microphysics.F90 into the zm_microp
! harness (METHODOLOGY.md: stubs never contain physics).
!
! time_manager: zm_mphy calls get_step_size() for the model time step
! (pure infrastructure -- the value is a driver input, set through
! set_step_size below and recorded in the golden metadata).
!
! ndrop_bam: abort-only placeholder for the BULK-aerosol droplet
! activation (eam/src/physics/cam/ndrop_bam.F90 ndrop_bam_run).
! zm_microphysics.F90 use-associates it at module scope so the symbol
! is needed to COMPILE, but the harness (like EAMv3 production, which
! runs MAM modal aerosols) only exercises aero%scheme == 'modal', so
! it is never CALLED: the modal path uses the REAL activate_drop_mam
! and nucleate_ice_conv modules, which are compiled unmodified.
! Documented in PORTING_PLAN.md row 9; per METHODOLOGY.md this stub
! only refuses to run (endrun).

module time_manager
  implicit none
  private
  public :: get_step_size, set_step_size
  integer :: harness_step_size = 1800
contains
  subroutine set_step_size(step_size_in)
    integer, intent(in) :: step_size_in
    harness_step_size = step_size_in
  end subroutine set_step_size

  integer function get_step_size()
    get_step_size = harness_step_size
  end function get_step_size
end module time_manager

module ndrop_bam
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private
  public :: ndrop_bam_run
contains
  subroutine ndrop_bam_run(wbar, tair, rhoair, na, pmode, nmode, ma, nact)
    use cam_abortutils, only: endrun
    integer,  intent(in) :: pmode, nmode
    real(r8), intent(in) :: wbar, tair, rhoair
    real(r8), intent(in) :: na(pmode), ma(pmode)
    real(r8), intent(out) :: nact
    nact = 0._r8
    call endrun('ndrop_bam_run stub called: zm_microp harness only '// &
                'supports the modal aerosol scheme (EAMv3 default)')
  end subroutine ndrop_bam_run
end module ndrop_bam
