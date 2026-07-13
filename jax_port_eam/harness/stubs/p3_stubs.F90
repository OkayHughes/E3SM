! Infrastructure-only stubs for the EAM P3 harness (METHODOLOGY.md:
! stubs never contain physics).
!
! - phys_control: the real eam/src/physics/cam/phys_control.F90 is a
!   large namelist-plumbing module; micro_p3.F90 only imports the
!   logical use_hetfrz_classnuc from it (a pure configuration flag,
!   EAMv3 phys="default" sets it .true.). The stub exposes the flag
!   plus a setter so the golden generator can drive BOTH branches.
!
! - debug_info: only referenced by wv_sat_scream.F90's polysvp1 error
!   path (report_error_info); polysvp1 itself is never called (P3 uses
!   MurphyKoop_svp). Print-only placeholder, documented as never
!   executed in normal operation.

module phys_control
  implicit none
  public
  logical :: use_hetfrz_classnuc = .true.  ! EAMv3 phys="default"
contains
  subroutine phys_control_set_hetfrz(flag)
    logical, intent(in) :: flag
    use_hetfrz_classnuc = flag
  end subroutine phys_control_set_hetfrz
end module phys_control

module debug_info
  implicit none
  public
contains
  subroutine report_error_info(msg, subname)
    character(len=*), intent(in) :: msg, subname
    write(*, *) 'DEBUG_INFO(', trim(subname), '): ', trim(msg)
  end subroutine report_error_info
end module debug_info
