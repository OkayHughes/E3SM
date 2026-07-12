! f2py driver for the orographic gravity-wave spine:
! gw_common_init + gw_prof + gw_oro_src + gw_drag_prof (ngwv = 0).
module gw_driver
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_gw_init(nlev, ktop_in, kbotbg_in, fcrit2_in, kwv_in, &
                         gravit_in, rair_in, alpha, orographic_only_in)
    use gw_common, only: gw_common_init
    use gw_oro, only: gw_oro_init
    integer, intent(in) :: nlev, ktop_in, kbotbg_in
    real(r8), intent(in) :: fcrit2_in, kwv_in, gravit_in, rair_in
    real(r8), intent(in) :: alpha(nlev + 1)
    integer, intent(in) :: orographic_only_in
    real(r8) :: cref(1)
    character(len=256) :: errstring
    cref = 0._r8
    call gw_common_init(nlev, 0, 0._r8, cref, &
         orographic_only_in /= 0, .false., .false., huge(1), &
         ktop_in, kbotbg_in, fcrit2_in, kwv_in, gravit_in, rair_in, &
         alpha, errstring)
    if (trim(errstring) /= '') stop 'gw_common_init failed'
    call gw_oro_init(errstring)
    if (trim(errstring) /= '') stop 'gw_oro_init failed'
  end subroutine drv_gw_init

  subroutine drv_gw_prof(ncol, nlev, cpair, t, pmid, pint, &
                         rhoi, ti, nm, ni)
    use gw_common, only: gw_prof
    integer, intent(in) :: ncol, nlev
    real(r8), intent(in) :: cpair
    real(r8), intent(in) :: t(ncol, nlev), pmid(ncol, nlev)
    real(r8), intent(in) :: pint(ncol, nlev + 1)
    real(r8), intent(out) :: rhoi(ncol, nlev + 1), ti(ncol, nlev + 1)
    real(r8), intent(out) :: nm(ncol, nlev), ni(ncol, nlev + 1)
    call gw_prof(ncol, cpair, t, pmid, pint, rhoi, ti, nm, ni)
  end subroutine drv_gw_prof

  subroutine drv_gw_oro_src(ncol, nlev, u, v, t, sgh, pmid, pint, dpm, &
                            zm, nm, src_level, tend_level, tau, ubm, &
                            ubi, xv, yv)
    use gw_oro, only: gw_oro_src
    integer, intent(in) :: ncol, nlev
    real(r8), intent(in) :: u(ncol, nlev), v(ncol, nlev), t(ncol, nlev)
    real(r8), intent(in) :: sgh(ncol)
    real(r8), intent(in) :: pmid(ncol, nlev), pint(ncol, nlev + 1)
    real(r8), intent(in) :: dpm(ncol, nlev), zm(ncol, nlev)
    real(r8), intent(in) :: nm(ncol, nlev)
    integer, intent(out) :: src_level(ncol), tend_level(ncol)
    real(r8), intent(out) :: tau(ncol, nlev + 1)
    real(r8), intent(out) :: ubm(ncol, nlev), ubi(ncol, nlev + 1)
    real(r8), intent(out) :: xv(ncol), yv(ncol)
    real(r8) :: tau3(ncol, 1, 0:nlev), c(ncol, 0:0)
    tau3 = 0._r8
    call gw_oro_src(ncol, u, v, t, sgh, pmid, pint, dpm, zm, nm, &
         src_level, tend_level, tau3, ubm, ubi, xv, yv, c)
    tau = tau3(:, 1, :)
  end subroutine drv_gw_oro_src

  subroutine drv_gw_drag_prof_oro(ncol, nlev, src_level, tend_level, &
       dt, lat, t, ti, pmid, pint, dpm, rdpm, piln, rhoi, nm, ni, &
       ubm, ubi, xv, yv, effgw, tau_in, tau_out, utgw, vtgw)
    use gw_common, only: gw_drag_prof
    integer, intent(in) :: ncol, nlev
    integer, intent(in) :: src_level(ncol), tend_level(ncol)
    real(r8), intent(in) :: dt, effgw
    real(r8), intent(in) :: lat(ncol)
    real(r8), intent(in) :: t(ncol, nlev), ti(ncol, nlev + 1)
    real(r8), intent(in) :: pmid(ncol, nlev), pint(ncol, nlev + 1)
    real(r8), intent(in) :: dpm(ncol, nlev), rdpm(ncol, nlev)
    real(r8), intent(in) :: piln(ncol, nlev + 1)
    real(r8), intent(in) :: rhoi(ncol, nlev + 1)
    real(r8), intent(in) :: nm(ncol, nlev), ni(ncol, nlev + 1)
    real(r8), intent(in) :: ubm(ncol, nlev), ubi(ncol, nlev + 1)
    real(r8), intent(in) :: xv(ncol), yv(ncol)
    real(r8), intent(in) :: tau_in(ncol, nlev + 1)
    real(r8), intent(out) :: tau_out(ncol, nlev + 1)
    real(r8), intent(out) :: utgw(ncol, nlev), vtgw(ncol, nlev)

    real(r8) :: tau3(ncol, 1, 0:nlev), c(ncol, 0:0)
    real(r8) :: kvtt(ncol, 0:nlev), q(ncol, nlev, 1), dse(ncol, nlev)
    real(r8) :: ttgw(ncol, nlev), qtgw(ncol, nlev, 1)
    real(r8) :: taucd(ncol, 0:nlev, 4), egwdffi(ncol, 0:nlev)
    real(r8) :: gwut(ncol, nlev, 1), dttdf(ncol, nlev)
    real(r8) :: dttke(ncol, nlev)

    tau3(:, 1, :) = tau_in
    c = 0._r8
    kvtt = 0._r8
    q = 0._r8
    dse = 0._r8
    call gw_drag_prof(ncol, 0, src_level, tend_level, .false., dt, &
         lat, t, ti, pmid, pint, dpm, rdpm, piln, rhoi, nm, ni, &
         ubm, ubi, xv, yv, effgw, c, kvtt, q, dse, tau3, utgw, vtgw, &
         ttgw, qtgw, taucd, egwdffi, gwut, dttdf, dttke)
    tau_out = tau3(:, 1, :)
  end subroutine drv_gw_drag_prof_oro

end module gw_driver
