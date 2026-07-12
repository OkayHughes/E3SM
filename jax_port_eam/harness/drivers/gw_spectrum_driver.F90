! f2py driver for the full-spectrum (non-orographic) gravity-wave
! suite: gw_common_init with pgwv > 0 + gw_front_init + gw_convect_init,
! gw_beres_src (convective source), gw_cm_src (frontal source),
! gw_drag_prof with ngwv = pgwv (stress profiles, gwd_project_tau,
! tendencies, gwd_precalc_rhoi constituent/DSE diffusion),
! momentum_energy_conservation, and the vdiff_lu_solver kernels.
!
! ubm/ubi driver outputs are zero-initialized before the src calls:
! the Fortran sources leave slots above/below the projection range
! unwritten (intent(out) => undefined); zeroing makes goldens
! deterministic. Those slots are never read by the downstream kernels.
module gw_spectrum_driver
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_gw_spec_init(nlev, npgwv, dc_in, ktop_in, kbotbg_in, &
       fcrit2_in, kwv_in, gravit_in, rair_in, alpha, tau_0_ubc_in, &
       taubgnd_in, frontgfc_in, kfront_in, plev_src_wind, &
       nh, nuh, mfcc_in, pref_edge_in)
    use gw_common,  only: gw_common_init
    use gw_front,   only: gw_front_init
    use gw_convect, only: gw_convect_init
    use ref_pres,   only: ref_pres_stub_set_edge
    integer, intent(in) :: nlev, npgwv, ktop_in, kbotbg_in, kfront_in
    integer, intent(in) :: tau_0_ubc_in, nh, nuh
    real(r8), intent(in) :: dc_in, fcrit2_in, kwv_in, gravit_in, rair_in
    real(r8), intent(in) :: alpha(nlev + 1)
    real(r8), intent(in) :: taubgnd_in, frontgfc_in, plev_src_wind
    real(r8), intent(in) :: mfcc_in(nh, 2 * nuh + 1, 2 * npgwv + 1)
    real(r8), intent(in) :: pref_edge_in(nlev + 1)
    real(r8) :: cref(-npgwv:npgwv)
    character(len=256) :: errstring
    integer :: l
    do l = -npgwv, npgwv
       cref(l) = dc_in * l
    end do
    call gw_common_init(nlev, npgwv, dc_in, cref, .false., .false., &
         tau_0_ubc_in /= 0, huge(1), ktop_in, kbotbg_in, fcrit2_in, &
         kwv_in, gravit_in, rair_in, alpha, errstring)
    if (trim(errstring) /= '') stop 'gw_common_init failed'
    call ref_pres_stub_set_edge(nlev + 1, pref_edge_in)
    call gw_front_init(taubgnd_in, frontgfc_in, kfront_in, errstring)
    if (trim(errstring) /= '') stop 'gw_front_init failed'
    call gw_convect_init(plev_src_wind, mfcc_in, errstring)
    if (trim(errstring) /= '') stop 'gw_convect_init failed'
  end subroutine drv_gw_spec_init

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

  subroutine drv_gw_beres_src(ncol, nlev, nwav, npgwv, lat, u, v, &
       netdt, zm, hcf, hdsf, hdmin, ssmin, use_old, &
       src_level, tend_level, tau, ubm, ubi, xv, yv, c, hdepth, maxq0)
    use gw_convect, only: gw_beres_src
    integer, intent(in) :: ncol, nlev, nwav, npgwv, use_old
    real(r8), intent(in) :: lat(ncol), u(ncol, nlev), v(ncol, nlev)
    real(r8), intent(in) :: netdt(ncol, nlev), zm(ncol, nlev)
    real(r8), intent(in) :: hcf, hdsf, hdmin, ssmin
    integer, intent(out) :: src_level(ncol), tend_level(ncol)
    real(r8), intent(out) :: tau(ncol, nwav, nlev + 1)
    real(r8), intent(out) :: ubm(ncol, nlev), ubi(ncol, nlev + 1)
    real(r8), intent(out) :: xv(ncol), yv(ncol), c(ncol, nwav)
    real(r8), intent(out) :: hdepth(ncol), maxq0(ncol)
    if (nwav /= 2 * npgwv + 1) stop 'nwav /= 2*pgwv+1'
    ubm = 0._r8
    ubi = 0._r8
    call gw_beres_src(ncol, npgwv, lat, u, v, netdt, zm, src_level, &
         tend_level, tau, ubm, ubi, xv, yv, c, hdepth, maxq0, &
         hcf, hdsf, hdmin, ssmin, use_old /= 0)
  end subroutine drv_gw_beres_src

  subroutine drv_gw_cm_src(ncol, nlev, nwav, npgwv, kbot, u, v, &
       frontgf, src_level, tend_level, tau, ubm, ubi, xv, yv, c)
    use gw_front, only: gw_cm_src
    integer, intent(in) :: ncol, nlev, nwav, npgwv, kbot
    real(r8), intent(in) :: u(ncol, nlev), v(ncol, nlev)
    real(r8), intent(in) :: frontgf(ncol, nlev)
    integer, intent(out) :: src_level(ncol), tend_level(ncol)
    real(r8), intent(out) :: tau(ncol, nwav, nlev + 1)
    real(r8), intent(out) :: ubm(ncol, nlev), ubi(ncol, nlev + 1)
    real(r8), intent(out) :: xv(ncol), yv(ncol), c(ncol, nwav)
    if (nwav /= 2 * npgwv + 1) stop 'nwav /= 2*pgwv+1'
    ubm = 0._r8
    ubi = 0._r8
    call gw_cm_src(ncol, npgwv, kbot, u, v, frontgf, &
         src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  end subroutine drv_gw_cm_src

  subroutine drv_gw_drag_prof_spec(ncol, nlev, nwav, npgwv, ncnst, &
       src_level, tend_level, do_taper, dt, lat, t, ti, pmid, pint, &
       dpm, rdpm, piln, rhoi, nm, ni, ubm, ubi, xv, yv, effgw, c, q, &
       dse, tau_in, tau_out, utgw, vtgw, ttgw, qtgw, taucd, egwdffi, &
       gwut, dttdf, dttke)
    use gw_common, only: gw_drag_prof
    integer, intent(in) :: ncol, nlev, nwav, npgwv, ncnst, do_taper
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
    real(r8), intent(in) :: c(ncol, nwav)
    real(r8), intent(in) :: q(ncol, nlev, ncnst), dse(ncol, nlev)
    real(r8), intent(in) :: tau_in(ncol, nwav, nlev + 1)
    real(r8), intent(out) :: tau_out(ncol, nwav, nlev + 1)
    real(r8), intent(out) :: utgw(ncol, nlev), vtgw(ncol, nlev)
    real(r8), intent(out) :: ttgw(ncol, nlev)
    real(r8), intent(out) :: qtgw(ncol, nlev, ncnst)
    real(r8), intent(out) :: taucd(ncol, nlev + 1, 4)
    real(r8), intent(out) :: egwdffi(ncol, nlev + 1)
    real(r8), intent(out) :: gwut(ncol, nlev, nwav)
    real(r8), intent(out) :: dttdf(ncol, nlev), dttke(ncol, nlev)
    real(r8) :: kvtt(ncol, 0:nlev)
    if (nwav /= 2 * npgwv + 1) stop 'nwav /= 2*pgwv+1'
    kvtt = 0._r8
    tau_out = tau_in
    call gw_drag_prof(ncol, npgwv, src_level, tend_level, &
         do_taper /= 0, dt, lat, t, ti, pmid, pint, dpm, rdpm, piln, &
         rhoi, nm, ni, ubm, ubi, xv, yv, effgw, c, kvtt, q, dse, &
         tau_out, utgw, vtgw, ttgw, qtgw, taucd, egwdffi, gwut, &
         dttdf, dttke)
  end subroutine drv_gw_drag_prof_spec

  subroutine drv_gwd_project_tau(ncol, nlev, nwav, npgwv, tend_level, &
       tau, ubi, c, xv, yv, taucd)
    use gw_common, only: gwd_project_tau
    integer, intent(in) :: ncol, nlev, nwav, npgwv
    integer, intent(in) :: tend_level(ncol)
    real(r8), intent(in) :: tau(ncol, nwav, nlev + 1)
    real(r8), intent(in) :: ubi(ncol, nlev + 1), c(ncol, nwav)
    real(r8), intent(in) :: xv(ncol), yv(ncol)
    real(r8), intent(out) :: taucd(ncol, nlev + 1, 4)
    if (nwav /= 2 * npgwv + 1) stop 'nwav /= 2*pgwv+1'
    taucd = 0._r8
    call gwd_project_tau(ncol, npgwv, tend_level, tau, ubi, c, xv, &
         yv, taucd)
  end subroutine drv_gwd_project_tau

  subroutine drv_momentum_energy_conservation(ncol, nlev, tend_level, &
       dt, taucd, pint, pdel, u, v, dudt_in, dvdt_in, dsdt_in, &
       utgw_in, vtgw_in, ttgw_in, dudt, dvdt, dsdt, utgw, vtgw, ttgw)
    use gw_common, only: momentum_energy_conservation
    integer, intent(in) :: ncol, nlev
    integer, intent(in) :: tend_level(ncol)
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: taucd(ncol, nlev + 1, 4)
    real(r8), intent(in) :: pint(ncol, nlev + 1), pdel(ncol, nlev)
    real(r8), intent(in) :: u(ncol, nlev), v(ncol, nlev)
    real(r8), intent(in) :: dudt_in(ncol, nlev), dvdt_in(ncol, nlev)
    real(r8), intent(in) :: dsdt_in(ncol, nlev)
    real(r8), intent(in) :: utgw_in(ncol, nlev), vtgw_in(ncol, nlev)
    real(r8), intent(in) :: ttgw_in(ncol, nlev)
    real(r8), intent(out) :: dudt(ncol, nlev), dvdt(ncol, nlev)
    real(r8), intent(out) :: dsdt(ncol, nlev)
    real(r8), intent(out) :: utgw(ncol, nlev), vtgw(ncol, nlev)
    real(r8), intent(out) :: ttgw(ncol, nlev)
    dudt = dudt_in
    dvdt = dvdt_in
    dsdt = dsdt_in
    utgw = utgw_in
    vtgw = vtgw_in
    ttgw = ttgw_in
    call momentum_energy_conservation(ncol, tend_level, dt, taucd, &
         pint, pdel, u, v, dudt, dvdt, dsdt, utgw, vtgw, ttgw)
  end subroutine drv_momentum_energy_conservation

  subroutine drv_vd_lu(ncol, nlev, ncnst, ksrf, kv, tmpi, rpdel, &
       ztodt, gravit_in, cc_top, ntop, nbot, q_in, cd_top, &
       q_out, ca, cc, dnom, ze)
    use vdiff_lu_solver, only: vd_lu_decomp, vd_lu_solve, lu_decomp
    integer, intent(in) :: ncol, nlev, ncnst, ntop, nbot
    real(r8), intent(in) :: ksrf(ncol), kv(ncol, nlev + 1)
    real(r8), intent(in) :: tmpi(ncol, nlev + 1), rpdel(ncol, nlev)
    real(r8), intent(in) :: ztodt, gravit_in
    real(r8), intent(in) :: cc_top(ncol), cd_top(ncol)
    real(r8), intent(in) :: q_in(ncol, nlev, ncnst)
    real(r8), intent(out) :: q_out(ncol, nlev, ncnst)
    real(r8), intent(out) :: ca(ncol, nlev), cc(ncol, nlev)
    real(r8), intent(out) :: dnom(ncol, nlev), ze(ncol, nlev)
    type(lu_decomp) :: decomp
    real(r8) :: qbuf(ncol, nlev)
    integer :: m
    call vd_lu_decomp(ncol, nlev, ncol, ksrf, kv, tmpi, rpdel, ztodt, &
         gravit_in, cc_top, ntop, nbot, decomp)
    ca = decomp%ca
    cc = decomp%cc
    dnom = decomp%dnom
    ze = decomp%ze
    do m = 1, ncnst
       qbuf = q_in(:, :, m)
       call vd_lu_solve(ncol, nlev, ncol, qbuf, decomp, ntop, nbot, &
            cd_top)
       q_out(:, :, m) = qbuf
    end do
    call decomp%finalize()
  end subroutine drv_vd_lu

end module gw_spectrum_driver
