! Infrastructure-only stub for share/util/shr_sys_mod.F90 (the real
! module drags in shr_log_mod/shr_abort_mod/MPI). Only the pieces the
! ported schemes touch: zm_conv_types uses shr_sys_flush in its
! log-print helpers (never called by the harness drivers).

module shr_sys_mod
  implicit none
contains
  subroutine shr_sys_flush(unit)
    integer, intent(in) :: unit
    flush(unit)
  end subroutine shr_sys_flush

  subroutine shr_sys_abort(string)
    character(len=*), intent(in), optional :: string
    if (present(string)) write(*, *) 'shr_sys_abort: ', string
    stop 1
  end subroutine shr_sys_abort
end module shr_sys_mod

! zm_conv_types calls mpibcast (a bare external in EAM's mpishorthand)
! from zm_param_mpi_broadcast. Single-process harness: no-op. Untyped
! legacy interface on purpose (callers pass r8/integer/logical buffers;
! compiled with -fallow-argument-mismatch like the rest of the build).
subroutine mpibcast(buf, count, datatype, root, comm)
  implicit none
  integer :: buf(*)
  integer :: count, datatype, root, comm
end subroutine mpibcast
