! Infrastructure-only stubs for the CLUBB f2py harness (eam_clubb_f).
! (METHODOLOGY.md: stubs never contain physics.)
!
! phys_control: CLUBB's advance_windm_edsclrm_module.F90 does
!   `use phys_control, only: use_od_fd` (turbulent orographic form
!   drag switch).  In EAM the flag comes from the phys_ctl_nl
!   namelist; its default (physics/cam/phys_control.F90) is
!   use_od_fd = .false., which is also the EAMv3 phys="default"
!   value (namelist_defaults_eam.xml does not override it).  A
!   setter is provided in case a later slice needs the .true. branch.

module phys_control
  implicit none
  logical, public :: use_od_fd = .false.
contains
  subroutine set_use_od_fd(flag)
    logical, intent(in) :: flag
    use_od_fd = flag
  end subroutine set_use_od_fd
end module phys_control
