! Infrastructure-only stub of ref_pres (reference-pressure grid
! metadata; no physics). gw_convect_init scans pref_edge to locate the
! steering-flow level; the harness driver fills it via the setter.

module ref_pres
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  public
  real(r8), allocatable :: pref_edge(:)
contains
  subroutine ref_pres_stub_set_edge(n, edges)
    integer,  intent(in) :: n
    real(r8), intent(in) :: edges(n)
    if (allocated(pref_edge)) deallocate(pref_edge)
    allocate(pref_edge(n))
    pref_edge = edges
  end subroutine ref_pres_stub_set_edge
end module ref_pres
