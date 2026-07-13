! f2py driver for the EAM P3 stratiform microphysics
! (components/eam/src/physics/p3/eam/micro_p3.F90, compiled UNMODIFIED,
! plus its real dependencies micro_p3_utils.F90, wv_sat_scream.F90,
! physics_utils.F90, scream_abortutils.F90).
!
! drv_init mirrors micro_p3_interface.F90 micro_p3_init exactly:
! micro_p3_utils_init with the physconst values (the physconst stub is
! derived verbatim from the real shr_const_mod), then p3_init reading
! p3_lookup_table_1.dat-v<version> from lookup_dir and generating the
! rain tables (p3_init_b with mu_r_constant = 0).
!
! use_hetfrz_classnuc (phys_control stub) and do_Cooper_inP3 (real
! micro_p3_utils module variable) are runtime-settable so the golden
! covers both branches of each.
!
! drv_p3_main packs outputs:
!   state_out(ncol,nlev,10): qc, nc, qr, nr, th, qv, qi, qm, ni, bm
!   diag_out(ncol,nlev,14): diag_eff_radius_qc, diag_eff_radius_qi,
!     rho_qi, qv2qi_depos_tend, precip_total_tend, nevapr,
!     qr_evap_tend, mu_c, lamc, liq_ice_exchange, vap_liq_exchange,
!     vap_ice_exchange, diag_equiv_reflectivity, diag_ze_rain
!     (+ diag_ze_ice as 15th... see DIAG_FIELDS in gen_p3_golden.py)
!   flux_out(ncol,nlev+1,5): precip_liq_flux, precip_ice_flux, rflx,
!     sflx, cflx
!   tend_out(ncol,nlev,49): the p3_tend_out catch-all
!   surf_out(ncol,2): precip_liq_surf, precip_ice_surf
module p3_driver
  implicit none
  ! local kind param (use-associated kinds are invisible to f2py's
  ! crackfortran; fbuild passes --f2cmap mapping r8 -> double)
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_init(lookup_dir, version, use_hetfrz, do_cooper)
    use micro_p3_utils, only: micro_p3_utils_init, do_Cooper_inP3
    use micro_p3,       only: p3_init
    use phys_control,   only: phys_control_set_hetfrz
    use physconst,      only: cpair, rair, rh2o, rhoh2o, mwh2o, mwdry, &
                              gravit, latvap, latice, cpliq, tmelt, pi
    character(len=256), intent(in) :: lookup_dir
    character(len=16), intent(in)  :: version
    integer, intent(in) :: use_hetfrz   ! 0/1
    integer, intent(in) :: do_cooper    ! 0/1

    call micro_p3_utils_init(cpair, rair, rh2o, rhoh2o, mwh2o, mwdry, &
         gravit, latvap, latice, cpliq, tmelt, pi, 6, .false.)
    do_Cooper_inP3 = (do_cooper /= 0)
    call phys_control_set_hetfrz(use_hetfrz /= 0)
    call p3_init(trim(lookup_dir), version)
  end subroutine drv_init

  subroutine drv_set_flags(use_hetfrz, do_cooper)
    ! toggle the configuration flags between golden configs
    use micro_p3_utils, only: do_Cooper_inP3
    use phys_control,   only: phys_control_set_hetfrz
    integer, intent(in) :: use_hetfrz, do_cooper
    call phys_control_set_hetfrz(use_hetfrz /= 0)
    do_Cooper_inP3 = (do_cooper /= 0)
  end subroutine drv_set_flags

  subroutine drv_get_tables(mu_r_user, revap_user, vn_user, vm_user)
    ! dump the p3_init_b-generated rain tables for the JAX side
    use micro_p3, only: p3_get_tables
    real(r8), intent(out) :: mu_r_user(150)
    real(r8), intent(out) :: revap_user(300, 10)
    real(r8), intent(out) :: vn_user(300, 10), vm_user(300, 10)
    call p3_get_tables(mu_r_user, revap_user, vn_user, vm_user)
  end subroutine drv_get_tables

  subroutine drv_qv_sat(n, t_atm, p_atm, i_wrt, qsat)
    ! direct qv_sat calls for unit bisection of the saturation port
    use wv_sat_scream, only: qv_sat
    integer, intent(in) :: n
    real(r8), intent(in) :: t_atm(n), p_atm(n)
    integer, intent(in) :: i_wrt
    real(r8), intent(out) :: qsat(n)
    integer :: i
    do i = 1, n
      qsat(i) = qv_sat(t_atm(i), p_atm(i), i_wrt)
    end do
  end subroutine drv_qv_sat

  subroutine drv_cloud_sed(ncol, nlev, dt, do_predict_nc, &
       qc_incld, rho, inv_rho, cld_frac_l, acn, inv_dz, &
       qc, nc, nc_incld, mu_c, lamc, &
       qc_out, nc_out, nc_incld_out, mu_c_out, lamc_out, &
       precip_liq_surf, cflx, qc_tend, nc_tend)
    ! direct cloud_sedimentation call for unit bisection (per column)
    use micro_p3,       only: cloud_sedimentation
    use micro_p3_utils, only: dnu
    integer, intent(in) :: ncol, nlev
    real(r8), intent(in) :: dt
    integer, intent(in) :: do_predict_nc
    real(r8), intent(in), dimension(ncol, nlev) :: qc_incld, rho, &
         inv_rho, cld_frac_l, acn, inv_dz, qc, nc, nc_incld, mu_c, lamc
    real(r8), intent(out), dimension(ncol, nlev) :: qc_out, nc_out, &
         nc_incld_out, mu_c_out, lamc_out, qc_tend, nc_tend
    real(r8), intent(out) :: precip_liq_surf(ncol)
    real(r8), intent(out), dimension(ncol, nlev + 1) :: cflx
    real(r8) :: qc_incld_l(nlev)
    integer :: i

    qc_out = qc; nc_out = nc; nc_incld_out = nc_incld
    mu_c_out = mu_c; lamc_out = lamc
    precip_liq_surf = 0._r8
    cflx = 0._r8
    qc_tend = qc
    nc_tend = nc
    do i = 1, ncol
      qc_incld_l = qc_incld(i, :)
      call cloud_sedimentation(1, nlev, 1, nlev, -1, &
           qc_incld_l, rho(i, :), inv_rho(i, :), cld_frac_l(i, :), &
           acn(i, :), inv_dz(i, :), dt, 1._r8 / dt, dnu, &
           do_predict_nc /= 0, &
           qc_out(i, :), nc_out(i, :), nc_incld_out(i, :), &
           mu_c_out(i, :), lamc_out(i, :), precip_liq_surf(i), &
           cflx(i, :), qc_tend(i, :), nc_tend(i, :))
    end do
  end subroutine drv_cloud_sed

  subroutine drv_rain_sed(ncol, nlev, dt, p3_max_mean_rain_size, &
       qr_incld, rho, inv_rho, rhofacr, cld_frac_r, inv_dz, &
       qr, nr, nr_incld, mu_r, lamr, &
       qr_out, nr_out, nr_incld_out, mu_r_out, lamr_out, &
       precip_liq_surf, precip_liq_flux, rflx, qr_tend, nr_tend)
    ! direct rain_sedimentation call for unit bisection (per column)
    use micro_p3, only: rain_sedimentation
    integer, intent(in) :: ncol, nlev
    real(r8), intent(in) :: dt, p3_max_mean_rain_size
    real(r8), intent(in), dimension(ncol, nlev) :: qr_incld, rho, &
         inv_rho, rhofacr, cld_frac_r, inv_dz, qr, nr, nr_incld, mu_r, lamr
    real(r8), intent(out), dimension(ncol, nlev) :: qr_out, nr_out, &
         nr_incld_out, mu_r_out, lamr_out, qr_tend, nr_tend
    real(r8), intent(out) :: precip_liq_surf(ncol)
    real(r8), intent(out), dimension(ncol, nlev + 1) :: precip_liq_flux, rflx
    real(r8) :: qr_incld_l(nlev)
    integer :: i

    qr_out = qr; nr_out = nr; nr_incld_out = nr_incld
    mu_r_out = mu_r; lamr_out = lamr
    precip_liq_surf = 0._r8
    precip_liq_flux = 0._r8
    rflx = 0._r8
    qr_tend = qr
    nr_tend = nr
    do i = 1, ncol
      qr_incld_l = qr_incld(i, :)
      call rain_sedimentation(1, nlev, 1, nlev, -1, &
           qr_incld_l, rho(i, :), inv_rho(i, :), rhofacr(i, :), &
           cld_frac_r(i, :), inv_dz(i, :), dt, 1._r8 / dt, &
           p3_max_mean_rain_size, &
           qr_out(i, :), nr_out(i, :), nr_incld_out(i, :), &
           mu_r_out(i, :), lamr_out(i, :), precip_liq_surf(i), &
           precip_liq_flux(i, :), rflx(i, :), qr_tend(i, :), nr_tend(i, :))
    end do
  end subroutine drv_rain_sed

  subroutine drv_ice_sed(ncol, nlev, dt, &
       rho, inv_rho, rhofaci, cld_frac_i, inv_dz, &
       qi, qi_incld, ni, ni_incld, qm, qm_incld, bm, bm_incld, &
       qi_out, ni_out, qm_out, bm_out, &
       qi_incld_out, ni_incld_out, qm_incld_out, bm_incld_out, &
       precip_ice_surf, sflx, qi_tend, ni_tend)
    ! direct ice_sedimentation call for unit bisection (per column)
    use micro_p3, only: ice_sedimentation
    integer, intent(in) :: ncol, nlev
    real(r8), intent(in) :: dt
    real(r8), intent(in), dimension(ncol, nlev) :: rho, inv_rho, &
         rhofaci, cld_frac_i, inv_dz, qi, qi_incld, ni, ni_incld, qm, &
         qm_incld, bm, bm_incld
    real(r8), intent(out), dimension(ncol, nlev) :: qi_out, ni_out, &
         qm_out, bm_out, qi_incld_out, ni_incld_out, qm_incld_out, &
         bm_incld_out, qi_tend, ni_tend
    real(r8), intent(out) :: precip_ice_surf(ncol)
    real(r8), intent(out), dimension(ncol, nlev + 1) :: sflx
    real(r8) :: precip_ice_flux(nlev + 1)
    integer :: i

    qi_out = qi; ni_out = ni; qm_out = qm; bm_out = bm
    qi_incld_out = qi_incld; ni_incld_out = ni_incld
    qm_incld_out = qm_incld; bm_incld_out = bm_incld
    precip_ice_surf = 0._r8
    sflx = 0._r8
    qi_tend = qi
    ni_tend = ni
    do i = 1, ncol
      precip_ice_flux = 0._r8
      call ice_sedimentation(1, nlev, 1, nlev, -1, &
           rho(i, :), inv_rho(i, :), rhofaci(i, :), cld_frac_i(i, :), &
           inv_dz(i, :), dt, 1._r8 / dt, &
           qi_out(i, :), qi_incld_out(i, :), ni_out(i, :), &
           qm_out(i, :), qm_incld_out(i, :), bm_out(i, :), &
           bm_incld_out(i, :), ni_incld_out(i, :), &
           precip_ice_surf(i), precip_ice_flux, sflx(i, :), &
           qi_tend(i, :), ni_tend(i, :))
    end do
  end subroutine drv_ice_sed

  subroutine drv_p3_main(ncol, nlev, dt, it, &
       do_predict_nc, do_prescribed_ccn, do_precip_off, &
       p3_autocon_coeff, p3_accret_coeff, p3_qc_autocon_expon, &
       p3_nc_autocon_expon, p3_qc_accret_expon, p3_wbf_coeff, &
       p3_mincdnc, p3_max_mean_rain_size, p3_embryonic_rain_size, &
       nccnst, &
       qc, nc, qr, nr, th, qv, qi, qm, ni, bm, &
       pres, dz, nc_nuceat_tend, nccn_prescribed, ni_activated, &
       frzimm, frzcnt, frzdep, inv_qc_relvar, dpres, exner, &
       cld_frac_r, cld_frac_l, cld_frac_i, qv_prev, t_prev, &
       col_location, &
       state_out, diag_out, flux_out, tend_out, surf_out)
    use physics_utils, only: rtype, btype
    use micro_p3,      only: p3_main

    integer, intent(in) :: ncol, nlev
    real(r8), intent(in) :: dt
    integer, intent(in) :: it
    integer, intent(in) :: do_predict_nc, do_prescribed_ccn, do_precip_off
    real(r8), intent(in) :: p3_autocon_coeff, p3_accret_coeff
    real(r8), intent(in) :: p3_qc_autocon_expon, p3_nc_autocon_expon
    real(r8), intent(in) :: p3_qc_accret_expon, p3_wbf_coeff
    real(r8), intent(in) :: p3_mincdnc, p3_max_mean_rain_size
    real(r8), intent(in) :: p3_embryonic_rain_size, nccnst
    real(r8), intent(in), dimension(ncol, nlev) :: qc, nc, qr, nr, th, &
         qv, qi, qm, ni, bm
    real(r8), intent(in), dimension(ncol, nlev) :: pres, dz, &
         nc_nuceat_tend, nccn_prescribed, ni_activated, frzimm, frzcnt, &
         frzdep, inv_qc_relvar, dpres, exner, cld_frac_r, cld_frac_l, &
         cld_frac_i, qv_prev, t_prev
    real(r8), intent(in), dimension(ncol, 3) :: col_location

    real(r8), intent(out) :: state_out(ncol, nlev, 10)
    real(r8), intent(out) :: diag_out(ncol, nlev, 15)
    real(r8), intent(out) :: flux_out(ncol, nlev + 1, 5)
    real(r8), intent(out) :: tend_out(ncol, nlev, 49)
    real(r8), intent(out) :: surf_out(ncol, 2)

    ! local copies (p3_main updates its prognostics in place)
    real(rtype), dimension(ncol, nlev) :: lqc, lnc, lqr, lnr, lth, lqv, &
         lqi, lqm, lni, lbm
    real(rtype), dimension(ncol) :: precip_liq_surf, precip_ice_surf
    real(rtype), dimension(ncol, nlev) :: diag_eff_radius_qc, &
         diag_eff_radius_qi, rho_qi, qv2qi_depos_tend, &
         precip_total_tend, nevapr, qr_evap_tend, mu_c, lamc, &
         liq_ice_exchange, vap_liq_exchange, vap_ice_exchange, &
         diag_equiv_reflectivity, diag_ze_rain, diag_ze_ice
    real(rtype), dimension(ncol, nlev + 1) :: precip_liq_flux, &
         precip_ice_flux, rflx, sflx, cflx
    real(rtype), dimension(ncol, nlev, 49) :: p3_tend_out
    logical(btype) :: l_predict_nc, l_prescribed_ccn, l_precip_off

    lqc = qc; lnc = nc; lqr = qr; lnr = nr; lth = th
    lqv = qv; lqi = qi; lqm = qm; lni = ni; lbm = bm
    l_predict_nc     = (do_predict_nc /= 0)
    l_prescribed_ccn = (do_prescribed_ccn /= 0)
    l_precip_off     = (do_precip_off /= 0)

    call p3_main(lqc, lnc, lqr, lnr, lth, lqv, dt, lqi, lqm, lni, lbm, &
         pres, dz, nc_nuceat_tend, nccn_prescribed, ni_activated, &
         frzimm, frzcnt, frzdep, inv_qc_relvar, it, &
         precip_liq_surf, precip_ice_surf, 1, ncol, 1, nlev, &
         diag_eff_radius_qc, diag_eff_radius_qi, rho_qi, &
         l_predict_nc, l_prescribed_ccn, &
         p3_autocon_coeff, p3_accret_coeff, p3_qc_autocon_expon, &
         p3_nc_autocon_expon, p3_qc_accret_expon, p3_wbf_coeff, &
         p3_mincdnc, p3_max_mean_rain_size, p3_embryonic_rain_size, &
         dpres, exner, qv2qi_depos_tend, precip_total_tend, nevapr, &
         qr_evap_tend, precip_liq_flux, precip_ice_flux, rflx, sflx, &
         cflx, cld_frac_r, cld_frac_l, cld_frac_i, p3_tend_out, mu_c, &
         lamc, liq_ice_exchange, vap_liq_exchange, vap_ice_exchange, &
         qv_prev, t_prev, col_location, l_precip_off, nccnst, &
         diag_equiv_reflectivity, diag_ze_rain, diag_ze_ice)

    state_out(:, :, 1)  = lqc
    state_out(:, :, 2)  = lnc
    state_out(:, :, 3)  = lqr
    state_out(:, :, 4)  = lnr
    state_out(:, :, 5)  = lth
    state_out(:, :, 6)  = lqv
    state_out(:, :, 7)  = lqi
    state_out(:, :, 8)  = lqm
    state_out(:, :, 9)  = lni
    state_out(:, :, 10) = lbm

    diag_out(:, :, 1)  = diag_eff_radius_qc
    diag_out(:, :, 2)  = diag_eff_radius_qi
    diag_out(:, :, 3)  = rho_qi
    diag_out(:, :, 4)  = qv2qi_depos_tend
    diag_out(:, :, 5)  = precip_total_tend
    diag_out(:, :, 6)  = nevapr
    diag_out(:, :, 7)  = qr_evap_tend
    diag_out(:, :, 8)  = mu_c
    diag_out(:, :, 9)  = lamc
    diag_out(:, :, 10) = liq_ice_exchange
    diag_out(:, :, 11) = vap_liq_exchange
    diag_out(:, :, 12) = vap_ice_exchange
    diag_out(:, :, 13) = diag_equiv_reflectivity
    diag_out(:, :, 14) = diag_ze_rain
    diag_out(:, :, 15) = diag_ze_ice

    flux_out(:, :, 1) = precip_liq_flux
    flux_out(:, :, 2) = precip_ice_flux
    flux_out(:, :, 3) = rflx
    flux_out(:, :, 4) = sflx
    flux_out(:, :, 5) = cflx

    tend_out = p3_tend_out
    surf_out(:, 1) = precip_liq_surf
    surf_out(:, 2) = precip_ice_surf
  end subroutine drv_p3_main

end module p3_driver
