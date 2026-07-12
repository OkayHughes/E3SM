! Infrastructure-only stub for eam/src/physics/cam/constituents.F90
! (the constituent registry: pure name/index/type bookkeeping, no
! physics). Only the pieces the compiled schemes touch:
!  - conv_water.F90 use-associates cnst_add/pcnst (never calls
!    cnst_add) and calls cnst_get_ind('CLDICE'/'CLDLIQ') in
!    conv_water_init; the harness physics_state stub carries q with
!    vapor=1, CLDLIQ=2, CLDICE=3, matching the indices returned here.
!  - zm_transport.F90 calls cnst_get_type_byind(m) to select the
!    dry-mixing-ratio branch; the per-constituent type is injected by
!    the harness driver via cnst_stub_set_type (the real module reads
!    it from the constituent registry filled at model init).

module constituents
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private
  save

  public :: pcnst, cnst_add, cnst_get_ind, cnst_get_type_byind
  public :: cnst_stub_set_type

  integer, parameter :: pcnst = 8      ! upper bound for stub arrays
  integer, parameter :: max_cnst = 64
  character(len=3) :: stub_cnst_type(max_cnst) = 'wet'

contains

  subroutine cnst_stub_set_type(ind, ttype)
    integer,          intent(in) :: ind
    character(len=*), intent(in) :: ttype
    stub_cnst_type(ind) = ttype
  end subroutine cnst_stub_set_type

  ! never called by the harness (conv_water_register only imports it);
  ! abort if it ever is
  subroutine cnst_add(name, mwc, cpc, qminc, ind)
    use cam_abortutils, only: endrun
    character(len=*), intent(in)  :: name
    real(r8),         intent(in)  :: mwc, cpc, qminc
    integer,          intent(out) :: ind
    ind = -1
    call endrun('cnst_add stub called: registry not available')
  end subroutine cnst_add

  subroutine cnst_get_ind(name, ind, abrtf)
    use cam_abortutils, only: endrun
    character(len=*),  intent(in)  :: name
    integer,           intent(out) :: ind
    logical, optional, intent(in)  :: abrtf
    select case (trim(name))
    case ('Q');      ind = 1
    case ('CLDLIQ'); ind = 2
    case ('CLDICE'); ind = 3
    case default
      ind = -1
      call endrun('cnst_get_ind stub: unknown constituent '//trim(name))
    end select
  end subroutine cnst_get_ind

  character*3 function cnst_get_type_byind(ind)
    integer, intent(in) :: ind
    cnst_get_type_byind = stub_cnst_type(ind)
  end function cnst_get_type_byind

end module constituents
