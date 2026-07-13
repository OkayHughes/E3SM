! f2py driver for the EAM CLUBB harness (slice 1: grid + saturation +
! tridiag solver + Skx/sigma_sqd_w helpers + pdf_closure, ADG1 path).
!
! Setup mirrors EAM's clubb_intr.F90 exactly:
!   clubb_ini_cam:  set_clubb_debug_level_api(0);
!                   read_parameters_api(-99,"",params) [defaults + any
!                   clubb_* overrides -- here the overrides are applied
!                   to the params vector on the Python side, equivalent
!                   by construction, see drv_default_params];
!                   setup_clubb_core_api(pverp, theta0=300,
!                   ts_nudge=86400, hydromet_dim=0, sclr_dim=0,
!                   sclr_tol=[], edsclr_dim, params,
!                   l_host_applies_sfc_fluxes=F, l_uv_nudge=F, "flatau",
!                   l_input_fields=F, l_implemented=T, grid_type=3, ...)
!   clubb_tend_cam (per column, heights change):
!                   setup_grid_heights_api + setup_parameters_api.
module clubb_driver
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)
  logical, save :: l_core_setup_done = .false.

contains

  !---------------------------------------------------------------------
  ! Default tunable parameter vector, exactly as EAM obtains it on a
  ! non-masterproc rank: clubb_param_readnl initializes every clubb_*
  ! override to init_value (no file is read; stub masterproc=.false.),
  ! then read_parameters(-99, "", params) packs the compiled-in CLUBB
  ! defaults.  EAMv3 namelist overrides are applied to this vector by
  ! the golden generator (indices from drv_param_indices).
  subroutine drv_get_nparams(nparams_out)
    use parameter_indices, only: nparams
    integer, intent(out) :: nparams_out
    nparams_out = nparams
  end subroutine drv_get_nparams

  subroutine drv_default_params(np_in, params)
    use parameter_indices, only: nparams
    use parameters_tunable, only: clubb_param_readnl, read_parameters
    integer, intent(in) :: np_in
    real(r8), intent(out) :: params(np_in)
    if (np_in /= nparams) stop 'drv_default_params: np_in /= nparams'
    call clubb_param_readnl('unused_no_masterproc')
    call read_parameters(-99, "", params)
  end subroutine drv_default_params

  ! 1-based indices into params for the clubb_param_nl overrides (order
  ! matches the EAMv3 namelist names documented in gen_clubb_golden.py)
  ! plus Skw_denom_coef, Skw_max_mag, beta, gamma coefs for the port.
  subroutine drv_param_indices(idx)
    use parameter_indices, only: &
        iC1, iC1b, iC1c, iC2rt, iC2thl, iC2rtthl, iC6rt, iC6rtb, &
        iC6rtc, iC6thl, iC6thlb, iC6thlc, iC7, iC7b, iC8, iC11, &
        iC11b, iC11c, iC14, ibeta, igamma_coef, igamma_coefb, &
        igamma_coefc, imu, inu1, ic_K10, ic_K10h, iwpxp_L_thresh, &
        ialtitude_thresh, iSkw_denom_coef, iSkw_max_mag
    integer, intent(out) :: idx(31)
    idx = (/ iC1, iC1b, iC1c, iC2rt, iC2thl, iC2rtthl, iC6rt, iC6rtb, &
             iC6rtc, iC6thl, iC6thlb, iC6thlc, iC7, iC7b, iC8, iC11, &
             iC11b, iC11c, iC14, ibeta, igamma_coef, igamma_coefb, &
             igamma_coefc, imu, inu1, ic_K10, ic_K10h, iwpxp_L_thresh, &
             ialtitude_thresh, iSkw_denom_coef, iSkw_max_mag /)
  end subroutine drv_param_indices

  !---------------------------------------------------------------------
  ! Full setup for one column's grid (call before any kernel drv_*).
  ! zi_g/zt_g are CLUBB-oriented heights (index 1 = surface ghost),
  ! exactly the arrays clubb_tend_cam builds (zi_g(1)=0 at the surface,
  ! zt_g(1) = -zt_g(2)).  edsclr_dim as in EAM (pcnst-aerosols+2; value
  ! only matters for later slices).
  subroutine drv_setup(nz, edsclr_dim_in, np_in, params_in, zi_g, zt_g, &
                       err_out)
    use clubb_api_module, only: setup_clubb_core_api, &
        setup_grid_heights_api, setup_parameters_api, &
        set_clubb_debug_level_api
    use parameter_indices, only: nparams
    integer, intent(in) :: nz, edsclr_dim_in, np_in
    real(r8), intent(in) :: params_in(np_in)
    real(r8), intent(in) :: zi_g(nz), zt_g(nz)
    integer, intent(out) :: err_out
    real(r8) :: sclr_tol_loc(0)
    real(r8) :: params_loc(nparams)
    integer :: err_code_loc

    if (np_in /= nparams) stop 'drv_setup: np_in /= nparams'
    params_loc = params_in(1:nparams)
    err_out = 0

    call set_clubb_debug_level_api(0)

    if (.not. l_core_setup_done) then
      call setup_clubb_core_api( &
           nz, 300._r8, 86400._r8, &
           0, 0, &
           sclr_tol_loc, edsclr_dim_in, params_loc, &
           .false., &
           .false., "flatau", .false., &
           .true., 3, zi_g(2), zi_g(1), zi_g(nz), &
           zi_g, zt_g, zi_g(1), &
           err_code_loc)
      if (err_code_loc /= 0) then
        err_out = err_code_loc
        return
      end if
      l_core_setup_done = .true.
    else
      ! heights (and, harmlessly, parameters) refresh: the EAM
      ! clubb_tend_cam per-column sequence
      call setup_grid_heights_api(.true., 3, zi_g(2), zi_g(1), &
                                  zi_g, zt_g)
      call setup_parameters_api(zi_g(2), params_loc, nz, 3, &
                                zi_g, zt_g, err_code_loc)
      if (err_code_loc /= 0) then
        err_out = err_code_loc
        return
      end if
    end if
  end subroutine drv_setup

  !---------------------------------------------------------------------
  subroutine drv_grid_arrays(nz, zm, zt, dzm, dzt, invrs_dzm, invrs_dzt, &
                             w_zt2zm, w_zm2zt)
    use grid_class, only: gr
    integer, intent(in) :: nz
    real(r8), intent(out) :: zm(nz), zt(nz), dzm(nz), dzt(nz)
    real(r8), intent(out) :: invrs_dzm(nz), invrs_dzt(nz)
    real(r8), intent(out) :: w_zt2zm(2, nz), w_zm2zt(2, nz)
    if (nz /= gr%nz) stop 'drv_grid_arrays: nz mismatch'
    zm = gr%zm
    zt = gr%zt
    dzm = gr%dzm
    dzt = gr%dzt
    invrs_dzm = gr%invrs_dzm
    invrs_dzt = gr%invrs_dzt
    w_zt2zm = gr%weights_zt2zm
    w_zm2zt = gr%weights_zm2zt
  end subroutine drv_grid_arrays

  subroutine drv_grid_ops(nz, azt, azm, zt2zm_out, zm2zt_out, &
                          ddzt_out, ddzm_out)
    use grid_class, only: zt2zm, zm2zt, ddzt, ddzm
    integer, intent(in) :: nz
    real(r8), intent(in) :: azt(nz), azm(nz)
    real(r8), intent(out) :: zt2zm_out(nz), zm2zt_out(nz)
    real(r8), intent(out) :: ddzt_out(nz), ddzm_out(nz)
    zt2zm_out = zt2zm(azt)
    zm2zt_out = zm2zt(azm)
    ddzt_out = ddzt(azt)   ! derivative of zt-field, on zm levels
    ddzm_out = ddzm(azm)   ! derivative of zm-field, on zt levels
  end subroutine drv_grid_ops

  !---------------------------------------------------------------------
  subroutine drv_sat(nz, p_in_pa, t_in_k, rsl, rsi)
    use saturation, only: sat_mixrat_liq, sat_mixrat_ice
    integer, intent(in) :: nz
    real(r8), intent(in) :: p_in_pa(nz), t_in_k(nz)
    real(r8), intent(out) :: rsl(nz), rsi(nz)
    rsl = sat_mixrat_liq(p_in_pa, t_in_k)
    rsi = sat_mixrat_ice(p_in_pa, t_in_k)
  end subroutine drv_sat

  !---------------------------------------------------------------------
  subroutine drv_skx(nz, xp2, xp3, x_tol, skx)
    use Skx_module, only: Skx_func
    integer, intent(in) :: nz
    real(r8), intent(in) :: xp2(nz), xp3(nz), x_tol
    real(r8), intent(out) :: skx(nz)
    skx = Skx_func(xp2, xp3, x_tol)
  end subroutine drv_skx

  subroutine drv_sigma_sqd_w(nz, gamma_skw_fnc, wp2, thlp2, rtp2, &
                             up2, vp2, wpthlp, wprtp, upwp, vpwp, &
                             sigma_sqd_w)
    use sigma_sqd_w_module, only: compute_sigma_sqd_w
    integer, intent(in) :: nz
    real(r8), intent(in), dimension(nz) :: gamma_skw_fnc, wp2, thlp2, &
        rtp2, up2, vp2, wpthlp, wprtp, upwp, vpwp
    real(r8), intent(out) :: sigma_sqd_w(nz)
    sigma_sqd_w = compute_sigma_sqd_w(gamma_skw_fnc, wp2, thlp2, rtp2, &
                                      up2, vp2, wpthlp, wprtp, upwp, vpwp)
  end subroutine drv_sigma_sqd_w

  !---------------------------------------------------------------------
  ! CLUBB's LAPACK tridiagonal wrapper (dgtsv).  supd/diag/subd/rhs are
  ! copied because tridag_solve/dgtsv overwrite them.
  subroutine drv_tridag_solve(ndim, nrhs, supd, diag, subd, rhs, &
                              solution, err_out)
    use lapack_wrap, only: tridag_solve
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: ndim, nrhs
    real(r8), intent(in) :: supd(ndim), diag(ndim), subd(ndim)
    real(r8), intent(in) :: rhs(ndim, nrhs)
    real(r8), intent(out) :: solution(ndim, nrhs)
    integer, intent(out) :: err_out
    real(r8) :: supd_c(ndim), diag_c(ndim), subd_c(ndim)
    real(r8) :: rhs_c(ndim, nrhs)
    supd_c = supd
    diag_c = diag
    subd_c = subd
    rhs_c = rhs
    err_code = clubb_no_error
    call tridag_solve("drv_tridag", ndim, nrhs, supd_c, diag_c, subd_c, &
                      rhs_c, solution)
    err_out = err_code
    err_code = clubb_no_error
  end subroutine drv_tridag_solve

  !---------------------------------------------------------------------
  ! pdf_closure (iiPDF_ADG1 compile-time path), hydromet_dim = 0,
  ! sclr_dim = 0, exactly as EAM calls it.
  !
  ! moments(:,17) slots:
  !   1 wp2rtp   2 wp2thlp  3 cloud_frac  4 ice_supersat_frac  5 rcm
  !   6 wpthvp   7 wp2thvp  8 rtpthvp     9 thlpthvp          10 wprcp
  !  11 wp2rcp  12 rtprcp  13 thlprcp    14 rcp2              15 uprcp
  !  16 vprcp   17 rc_coef
  ! (wp4, wprtp2, wpthlp2, wprtpthlp are NOT computed under the EAM
  ! configuration: l_explicit_turbulent_adv_wp3/xpyp = .false. and all
  ! stats indices are 0 with l_stats=.false. -- intent(out) garbage.)
  !
  ! pdfp(:,47) slots: the 47 pdf_params fields in the order of
  ! pdf_parameter_module's init_pdf_params (w_1, w_2, varnce_w_1,
  ! varnce_w_2, rt_1, rt_2, varnce_rt_1, varnce_rt_2, thl_1, thl_2,
  ! varnce_thl_1, varnce_thl_2, corr_w_rt_1, corr_w_rt_2, corr_w_thl_1,
  ! corr_w_thl_2, corr_rt_thl_1, corr_rt_thl_2, alpha_thl, alpha_rt,
  ! crt_1, crt_2, cthl_1, cthl_2, chi_1, chi_2, stdev_chi_1,
  ! stdev_chi_2, stdev_eta_1, stdev_eta_2, covar_chi_eta_1,
  ! covar_chi_eta_2, corr_w_chi_1, corr_w_chi_2, corr_w_eta_1,
  ! corr_w_eta_2, corr_chi_eta_1, corr_chi_eta_2, rsatl_1, rsatl_2,
  ! rc_1, rc_2, cloud_frac_1, cloud_frac_2, mixt_frac,
  ! ice_supersat_frac_1, ice_supersat_frac_2).
  subroutine drv_pdf_closure(nz, p_in_pa, exner, thv_ds, wm, wp2, wp3, &
                             sigma_sqd_w, skw, skthl, skrt, &
                             rtm, rtp2, wprtp, thlm, thlp2, wpthlp, &
                             um, up2, upwp, vm, vp2, vpwp, rtpthlp, &
                             moments, pdfp, sigma_sqd_w_out, err_out)
    use pdf_closure_module, only: pdf_closure
    use pdf_parameter_module, only: pdf_parameter, implicit_coefs_terms, &
        init_pdf_params
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: nz
    real(r8), intent(in), dimension(nz) :: p_in_pa, exner, thv_ds, wm, &
        wp2, wp3, sigma_sqd_w, skw, skthl, skrt, rtm, rtp2, wprtp, &
        thlm, thlp2, wpthlp, um, up2, upwp, vm, vp2, vpwp, rtpthlp
    real(r8), intent(out) :: moments(nz, 17)
    real(r8), intent(out) :: pdfp(nz, 47)
    real(r8), intent(out) :: sigma_sqd_w_out(nz)
    integer, intent(out) :: err_out

    real(r8), dimension(nz) :: sigw_loc
    real(r8), dimension(nz, 0) :: sclrm, wpsclrp, sclrp2, sclrprtp, &
        sclrpthlp, wpsclrprtp, wpsclrp2, sclrpthvp, wpsclrpthlp, &
        sclrprcp, wp2sclrp
    real(r8), dimension(nz, 0) :: wphydrometp, wp2hmp, rtphmp, thlphmp
    real(r8), dimension(nz) :: wp4, wprtp2, wp2rtp, wpthlp2, wp2thlp, &
        wprtpthlp, cloud_frac, ice_supersat_frac, rcm, wpthvp, wp2thvp, &
        rtpthvp, thlpthvp, wprcp, wp2rcp, rtprcp, thlprcp, rcp2, &
        uprcp, vprcp, rc_coef
    real(r8), dimension(nz) :: f_w, f_rt, f_thl, min_f_w, max_f_w, &
        min_f_rt, max_f_rt, min_f_thl, max_f_thl
    type(pdf_parameter) :: pdf_params
    type(implicit_coefs_terms), dimension(nz) :: pdf_ict

    call init_pdf_params(nz, pdf_params)

    sigw_loc = sigma_sqd_w
    err_code = clubb_no_error

    ! Pre-fill the outputs that the EAM configuration leaves
    ! uncomputed, so their slots are deterministic (not recorded).
    wp4 = -9999._r8
    wprtp2 = -9999._r8
    wpthlp2 = -9999._r8
    wprtpthlp = -9999._r8

    call pdf_closure(0, p_in_pa, exner, thv_ds, &
                     wm, wp2, wp3, sigw_loc, &
                     skw, skthl, skrt, &
                     rtm, rtp2, wprtp, &
                     thlm, thlp2, wpthlp, &
                     um, up2, upwp, &
                     vm, vp2, vpwp, &
                     rtpthlp, &
                     sclrm, wpsclrp, sclrp2, &
                     sclrprtp, sclrpthlp, &
                     wphydrometp, wp2hmp, &
                     rtphmp, thlphmp, &
                     wp4, wprtp2, wp2rtp, &
                     wpthlp2, wp2thlp, wprtpthlp, &
                     cloud_frac, ice_supersat_frac, &
                     rcm, wpthvp, wp2thvp, rtpthvp, &
                     thlpthvp, wprcp, wp2rcp, rtprcp, &
                     thlprcp, rcp2, &
                     uprcp, vprcp, &
                     pdf_params, pdf_ict, &
                     f_w, f_rt, f_thl, &
                     min_f_w, max_f_w, &
                     min_f_rt, max_f_rt, &
                     min_f_thl, max_f_thl, &
                     wpsclrprtp, wpsclrp2, sclrpthvp, &
                     wpsclrpthlp, sclrprcp, wp2sclrp, &
                     rc_coef)

    err_out = err_code
    err_code = clubb_no_error

    moments(:, 1) = wp2rtp
    moments(:, 2) = wp2thlp
    moments(:, 3) = cloud_frac
    moments(:, 4) = ice_supersat_frac
    moments(:, 5) = rcm
    moments(:, 6) = wpthvp
    moments(:, 7) = wp2thvp
    moments(:, 8) = rtpthvp
    moments(:, 9) = thlpthvp
    moments(:, 10) = wprcp
    moments(:, 11) = wp2rcp
    moments(:, 12) = rtprcp
    moments(:, 13) = thlprcp
    moments(:, 14) = rcp2
    moments(:, 15) = uprcp
    moments(:, 16) = vprcp
    moments(:, 17) = rc_coef

    sigma_sqd_w_out = sigw_loc

    pdfp(:, 1) = pdf_params%w_1
    pdfp(:, 2) = pdf_params%w_2
    pdfp(:, 3) = pdf_params%varnce_w_1
    pdfp(:, 4) = pdf_params%varnce_w_2
    pdfp(:, 5) = pdf_params%rt_1
    pdfp(:, 6) = pdf_params%rt_2
    pdfp(:, 7) = pdf_params%varnce_rt_1
    pdfp(:, 8) = pdf_params%varnce_rt_2
    pdfp(:, 9) = pdf_params%thl_1
    pdfp(:, 10) = pdf_params%thl_2
    pdfp(:, 11) = pdf_params%varnce_thl_1
    pdfp(:, 12) = pdf_params%varnce_thl_2
    pdfp(:, 13) = pdf_params%corr_w_rt_1
    pdfp(:, 14) = pdf_params%corr_w_rt_2
    pdfp(:, 15) = pdf_params%corr_w_thl_1
    pdfp(:, 16) = pdf_params%corr_w_thl_2
    pdfp(:, 17) = pdf_params%corr_rt_thl_1
    pdfp(:, 18) = pdf_params%corr_rt_thl_2
    pdfp(:, 19) = pdf_params%alpha_thl
    pdfp(:, 20) = pdf_params%alpha_rt
    pdfp(:, 21) = pdf_params%crt_1
    pdfp(:, 22) = pdf_params%crt_2
    pdfp(:, 23) = pdf_params%cthl_1
    pdfp(:, 24) = pdf_params%cthl_2
    pdfp(:, 25) = pdf_params%chi_1
    pdfp(:, 26) = pdf_params%chi_2
    pdfp(:, 27) = pdf_params%stdev_chi_1
    pdfp(:, 28) = pdf_params%stdev_chi_2
    pdfp(:, 29) = pdf_params%stdev_eta_1
    pdfp(:, 30) = pdf_params%stdev_eta_2
    pdfp(:, 31) = pdf_params%covar_chi_eta_1
    pdfp(:, 32) = pdf_params%covar_chi_eta_2
    pdfp(:, 33) = pdf_params%corr_w_chi_1
    pdfp(:, 34) = pdf_params%corr_w_chi_2
    pdfp(:, 35) = pdf_params%corr_w_eta_1
    pdfp(:, 36) = pdf_params%corr_w_eta_2
    pdfp(:, 37) = pdf_params%corr_chi_eta_1
    pdfp(:, 38) = pdf_params%corr_chi_eta_2
    pdfp(:, 39) = pdf_params%rsatl_1
    pdfp(:, 40) = pdf_params%rsatl_2
    pdfp(:, 41) = pdf_params%rc_1
    pdfp(:, 42) = pdf_params%rc_2
    pdfp(:, 43) = pdf_params%cloud_frac_1
    pdfp(:, 44) = pdf_params%cloud_frac_2
    pdfp(:, 45) = pdf_params%mixt_frac
    pdfp(:, 46) = pdf_params%ice_supersat_frac_1
    pdfp(:, 47) = pdf_params%ice_supersat_frac_2
  end subroutine drv_pdf_closure

end module clubb_driver
