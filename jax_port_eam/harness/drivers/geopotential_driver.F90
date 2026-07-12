! f2py driver for geopotential_t / geopotential_dse.
module geopotential_driver
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_geopotential_t(ncol, nlev, piln, pmln, pint, pmid, &
                                pdel, rpdel, t, q, rair, gravit, zvir, &
                                zi, zm)
    use geopotential, only: geopotential_t
    integer, intent(in) :: ncol, nlev
    real(r8), intent(in) :: piln(ncol, nlev + 1), pmln(ncol, nlev)
    real(r8), intent(in) :: pint(ncol, nlev + 1), pmid(ncol, nlev)
    real(r8), intent(in) :: pdel(ncol, nlev), rpdel(ncol, nlev)
    real(r8), intent(in) :: t(ncol, nlev), q(ncol, nlev)
    real(r8), intent(in) :: rair(ncol, nlev), zvir(ncol, nlev)
    real(r8), intent(in) :: gravit
    real(r8), intent(out) :: zi(ncol, nlev + 1), zm(ncol, nlev)
    call geopotential_t(piln, pmln, pint, pmid, pdel, rpdel, t, q, &
         rair, gravit, zvir, zi, zm, ncol)
  end subroutine drv_geopotential_t

end module geopotential_driver
