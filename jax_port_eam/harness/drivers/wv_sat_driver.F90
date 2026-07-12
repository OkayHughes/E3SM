! f2py driver for wv_sat_methods + wv_saturation (plain array
! interfaces around the elemental/optional-argument originals).
module wv_sat_driver
  implicit none
  ! local kind param (f2py can evaluate this; a use-associated kind from
  ! shr_kind_mod is invisible to f2py's crackfortran and silently
  ! becomes real(4) in the wrappers)
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_init()
    use wv_saturation, only: wv_sat_init
    call wv_sat_init()
  end subroutine drv_init

  ! raw scheme SVP formulae by index (0=OldGoffGratch, 1=GoffGratch,
  ! 2=MurphyKoop, 3=Bolton)
  subroutine drv_svp(n, t, idx, es_water, es_ice)
    use wv_sat_methods, only: wv_sat_svp_water, wv_sat_svp_ice
    integer, intent(in) :: n, idx
    real(r8), intent(in) :: t(n)
    real(r8), intent(out) :: es_water(n), es_ice(n)
    integer :: i
    do i = 1, n
      es_water(i) = wv_sat_svp_water(t(i), idx)
      es_ice(i) = wv_sat_svp_ice(t(i), idx)
    end do
  end subroutine drv_svp

  subroutine drv_svp_trans(n, t, idx, es)
    use wv_sat_methods, only: wv_sat_svp_trans
    integer, intent(in) :: n, idx
    real(r8), intent(in) :: t(n)
    real(r8), intent(out) :: es(n)
    integer :: i
    do i = 1, n
      es(i) = wv_sat_svp_trans(t(i), idx)
    end do
  end subroutine drv_svp_trans

  ! mixed-phase table qsat with all derivatives
  subroutine drv_qsat(n, t, p, es, qs, gam, dqsdt, enthalpy)
    use wv_saturation, only: qsat
    integer, intent(in) :: n
    real(r8), intent(in) :: t(n), p(n)
    real(r8), intent(out) :: es(n), qs(n), gam(n), dqsdt(n), enthalpy(n)
    call qsat(t, p, es, qs, gam=gam, dqsdt=dqsdt, enthalpy=enthalpy)
  end subroutine drv_qsat

  subroutine drv_qsat_water(n, t, p, es, qs, gam, dqsdt, enthalpy)
    use wv_saturation, only: qsat_water
    integer, intent(in) :: n
    real(r8), intent(in) :: t(n), p(n)
    real(r8), intent(out) :: es(n), qs(n), gam(n), dqsdt(n), enthalpy(n)
    call qsat_water(t, p, es, qs, gam=gam, dqsdt=dqsdt, enthalpy=enthalpy)
  end subroutine drv_qsat_water

  subroutine drv_qsat_ice(n, t, p, es, qs, gam, dqsdt, enthalpy)
    use wv_saturation, only: qsat_ice
    integer, intent(in) :: n
    real(r8), intent(in) :: t(n), p(n)
    real(r8), intent(out) :: es(n), qs(n), gam(n), dqsdt(n), enthalpy(n)
    call qsat_ice(t, p, es, qs, gam=gam, dqsdt=dqsdt, enthalpy=enthalpy)
  end subroutine drv_qsat_ice

  subroutine drv_estblf(n, t, es)
    use wv_saturation, only: estblf
    integer, intent(in) :: n
    real(r8), intent(in) :: t(n)
    real(r8), intent(out) :: es(n)
    es = estblf(t)
  end subroutine drv_estblf

  subroutine drv_findsp(n, q, t, p, use_ice, tsp, qsp, status)
    use wv_saturation, only: findsp
    integer, intent(in) :: n
    real(r8), intent(in) :: q(n), t(n), p(n)
    integer, intent(in) :: use_ice
    real(r8), intent(out) :: tsp(n), qsp(n)
    integer, intent(out) :: status(n)
    call findsp(q, t, p, use_ice /= 0, tsp, qsp, status)
  end subroutine drv_findsp

end module wv_sat_driver
