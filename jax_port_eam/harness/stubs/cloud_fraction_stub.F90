! Infrastructure-only stub for eam/src/physics/cam/cloud_fraction.F90.
! cldfrc2m_init() fetches its tuning parameters through
! cldfrc_getparams(), which in the real module merely returns values
! read from the cldfrc_nl namelist -- parameter plumbing, not physics.
! The stub stores values injected by the harness driver
! (cldfrc_stub_set) and returns them through the real, unmodified
! cldfrc_getparams keyword interface. Nothing else of cloud_fraction
! (the RK cldfrc scheme itself) is provided or needed by cldfrc2m.

module cloud_fraction
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private
  save

  public :: cldfrc_getparams, cldfrc_stub_set

  real(r8) :: rhminl          = 0.0_r8
  real(r8) :: rhminl_adj_land = 0.0_r8
  real(r8) :: rhminh          = 0.0_r8
  real(r8) :: premit          = 0.0_r8
  real(r8) :: premib          = 0.0_r8
  integer  :: iceopt          = 0
  real(r8) :: icecrit         = 0.0_r8
  real(r8) :: minice          = 0.0_r8

contains

  subroutine cldfrc_stub_set(rhminl_in, rhminl_adj_land_in, rhminh_in, &
       premit_in, premib_in, iceopt_in, icecrit_in, minice_in)
    real(r8), intent(in) :: rhminl_in, rhminl_adj_land_in, rhminh_in
    real(r8), intent(in) :: premit_in, premib_in
    integer,  intent(in) :: iceopt_in
    real(r8), intent(in) :: icecrit_in, minice_in
    rhminl          = rhminl_in
    rhminl_adj_land = rhminl_adj_land_in
    rhminh          = rhminh_in
    premit          = premit_in
    premib          = premib_in
    iceopt          = iceopt_in
    icecrit         = icecrit_in
    minice          = minice_in
  end subroutine cldfrc_stub_set

  ! Signature copied from the real cloud_fraction.F90 (keyword names
  ! must match cldfrc2m_init's call).
  subroutine cldfrc_getparams(rhminl_out, rhminl_adj_land_out, &
       rhminh_out, premit_out, premib_out, iceopt_out, icecrit_out, &
       minice_out)
    real(r8), intent(out), optional :: rhminl_out
    real(r8), intent(out), optional :: rhminl_adj_land_out
    real(r8), intent(out), optional :: rhminh_out
    real(r8), intent(out), optional :: premit_out
    real(r8), intent(out), optional :: premib_out
    integer,  intent(out), optional :: iceopt_out
    real(r8), intent(out), optional :: icecrit_out
    real(r8), intent(out), optional :: minice_out

    if (present(rhminl_out))          rhminl_out          = rhminl
    if (present(rhminl_adj_land_out)) rhminl_adj_land_out = rhminl_adj_land
    if (present(rhminh_out))          rhminh_out          = rhminh
    if (present(premit_out))          premit_out          = premit
    if (present(premib_out))          premib_out          = premib
    if (present(iceopt_out))          iceopt_out          = iceopt
    if (present(icecrit_out))         icecrit_out         = icecrit
    if (present(minice_out))          minice_out          = minice
  end subroutine cldfrc_getparams

end module cloud_fraction
