! f2py driver for the EAM CLUBB harness.
! Slice A: grid + saturation + tridiag solver + Skx/sigma_sqd_w
!          helpers + pdf_closure (ADG1 path).
! Slice B: pdf_closure_driver (the zt+zm double pdf_closure call,
!          trapezoidal-rule averaging, compute_cloud_cover, clip_rcm)
!          via the verbatim extraction module clubb_pdf_extract (see
!          build_clubb.py), plus the full advance_clubb_core for
!          end-to-end validation of that extraction.
! Slice C: advance_xp2_xpyp (PUBLIC -- driven directly, no extraction
!          needed) + clip_covars_denom, plus getters for the extra
!          tunables (C4/C5/C14/c_K2/c_K9/nu2/nu9), the
!          setup_parameters-derived nu2/nu9_vert_res_dep profiles, and
!          the model_flags/sponge configuration the routine reads.
! Slice D: Lscale/tau infrastructure.  compute_mixing_length,
!          calc_brunt_vaisala_freq_sqd, calc_stability_correction,
!          term_wp2_splat/term_wp3_splat and calc_surface_varnce are
!          all PUBLIC and driven directly; drv_lscale_tau_segment is a
!          documented VERBATIM TRANSCRIPTION of the inline
!          advance_clubb_core segment (em -> thvm -> Lscale -> tau ->
!          Kh -> splat -> surface variances -> stability-corrected
!          tau), validated end-to-end against the REAL
!          advance_clubb_core through its khzm/khzt diagnostics
!          (bitwise, asserted in gen_clubb_golden.py).
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

  !---------------------------------------------------------------------
  ! Set the two model_flags module variables EAM's clubb_intr assigns
  ! outside of setup_clubb_core_api: ipdf_call_placement (EAMv3
  ! phys="default" namelist: clubb_ipdf_call_placement = 2 =
  ! ipdf_post_advance_fields) and l_do_expldiff_rtm_thlm
  ! (= do_expldiff, clubb_expldiff = .true. in EAMv3).  Only
  ! advance_clubb_core reads them; pdf_closure_driver itself is
  ! placement-independent.
  subroutine drv_set_eam_flags(ipdf, l_expldiff)
    use model_flags, only: ipdf_call_placement, l_do_expldiff_rtm_thlm
    integer, intent(in) :: ipdf
    logical, intent(in) :: l_expldiff
    ipdf_call_placement = ipdf
    l_do_expldiff_rtm_thlm = l_expldiff
  end subroutine drv_set_eam_flags

  !---------------------------------------------------------------------
  ! pdf_closure_driver exactly as advance_clubb_core invokes it under
  ! the EAMv3 configuration (hydromet_dim=0, sclr_dim=0, l_stats=F,
  ! flags: l_gamma_Skw=T, l_call_pdf_closure_twice=T,
  ! l_trapezoidal_rule_zt/zm=T, l_use_cloud_cover=T,
  ! l_use_ice_latent=F, l_rcm_supersat_adj=F, l_rtm_nudge=F), via the
  ! verbatim extraction module clubb_pdf_extract.
  !
  ! Level placement of the inputs (CLUBB orientation, index 1 =
  ! surface): momentum-level: wprtp, wpthlp, rtp2, thlp2, rtpthlp,
  ! wp2, up2, upwp, vp2, vpwp, wm_zm, thv_ds_zm; thermo-level: thlm,
  ! rtm, rtp3, thlp3, wp3, wm_zt, um, vm, p_in_pa, exner, thv_ds_zt,
  ! rfrzm.
  !
  ! outs(:,29) slots (only outputs DEFINED under the EAM config; the
  ! intent(out) garbage slots wp4/wprtp2/wpthlp2/wprtpthlp/
  ! Skw_velocity/rtm_frz/thlm_frz are not returned):
  !   1 rcm            2 cloud_frac    3 ice_supersat_frac
  !   4 wprcp          5 sigma_sqd_w   6 wpthvp     7 wp2thvp
  !   8 rtpthvp        9 thlpthvp     10 rc_coef   11 rcm_in_layer
  !  12 cloud_cover   13 rcp2_zt      14 thlprcp   15 rc_coef_zm
  !  16 wp2rtp        17 wp2thlp      18 wp2rcp    19 rtprcp
  !  20 rcp2          21 uprcp        22 vprcp     23 cloud_frac_zm
  !  24 ice_supersat_frac_zm  25 rtm_zm  26 thlm_zm  27 rcm_zm
  !  28 rcm_supersat_adj      29 sigma_sqd_w_zt (module diagnostic
  !                              set by the routine)
  subroutine drv_pdf_closure_driver(nz, dt, wprtp, thlm, wpthlp, &
      rtp2, rtp3, thlp2, thlp3, rtpthlp, wp2, wp3, wm_zm, wm_zt, &
      um, up2, upwp, vm, vp2, vpwp, p_in_pa, exner, thv_ds_zm, &
      thv_ds_zt, rfrzm, rtm, rtm_out, outs, pdfp_zt, pdfp_zm, err_out)
    use clubb_pdf_extract, only: pdf_closure_driver
    use clubb_driver_helpers, only: pack_pdf_params
    use pdf_parameter_module, only: pdf_parameter, implicit_coefs_terms, &
        init_pdf_params
    use variables_diagnostic_module, only: sigma_sqd_w_zt
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: nz
    real(r8), intent(in) :: dt
    real(r8), intent(in), dimension(nz) :: wprtp, thlm, wpthlp, rtp2, &
        rtp3, thlp2, thlp3, rtpthlp, wp2, wp3, wm_zm, wm_zt, um, up2, &
        upwp, vm, vp2, vpwp, p_in_pa, exner, thv_ds_zm, thv_ds_zt, &
        rfrzm, rtm
    real(r8), intent(out) :: rtm_out(nz)
    real(r8), intent(out) :: outs(nz, 29)
    real(r8), intent(out) :: pdfp_zt(nz, 47), pdfp_zm(nz, 47)
    integer, intent(out) :: err_out

    real(r8), dimension(nz) :: rtm_loc
    real(r8), dimension(nz, 0) :: hydromet, wphydrometp, wp2hmp, &
        rtphmp_zt, thlphmp_zt, sclrm, wpsclrp, sclrp2, sclrprtp, &
        sclrpthlp, sclrpthvp, wp2sclrp, wpsclrp2, sclrprcp, &
        wpsclrprtp, wpsclrpthlp
    real(r8), dimension(nz) :: rcm, cloud_frac, ice_supersat_frac, &
        wprcp, sigma_sqd_w, wpthvp, wp2thvp, rtpthvp, thlpthvp, &
        rc_coef, rcm_in_layer, cloud_cover, rcp2_zt, thlprcp, &
        rc_coef_zm, rtm_frz, thlm_frz, wp4, wp2rtp, wprtp2, wp2thlp, &
        wpthlp2, wprtpthlp, wp2rcp, rtprcp, rcp2, uprcp, vprcp, &
        skw_velocity, cloud_frac_zm, ice_supersat_frac_zm, rtm_zm, &
        thlm_zm, rcm_zm, rcm_supersat_adj
    type(pdf_parameter) :: pdf_params, pdf_params_frz, pdf_params_zm
    type(implicit_coefs_terms), dimension(nz) :: pdf_ict

    call init_pdf_params(nz, pdf_params)
    call init_pdf_params(nz, pdf_params_frz)
    call init_pdf_params(nz, pdf_params_zm)

    rtm_loc = rtm
    err_code = clubb_no_error

    ! Sentinel-fill the outputs the EAM configuration never computes,
    ! so they are deterministic locals (not recorded).
    wp4 = -9999._r8
    wprtp2 = -9999._r8
    wpthlp2 = -9999._r8
    wprtpthlp = -9999._r8
    skw_velocity = -9999._r8
    rtm_frz = -9999._r8
    thlm_frz = -9999._r8

    call pdf_closure_driver( dt, 0, wprtp,                 & ! Intent(in)
                             thlm, wpthlp, rtp2, rtp3,     & ! Intent(in)
                             thlp2, thlp3, rtpthlp, wp2,   & ! Intent(in)
                             wp3, wm_zm, wm_zt,            & ! Intent(in)
                             um, up2, upwp,                & ! Intent(in)
                             vm, vp2, vpwp,                & ! Intent(in)
                             p_in_pa, exner,               & ! Intent(in)
                             thv_ds_zm, thv_ds_zt,         & ! Intent(in)
                             rfrzm, hydromet, wphydrometp, & ! Intent(in)
                             wp2hmp, rtphmp_zt, thlphmp_zt,& ! Intent(in)
                             sclrm, wpsclrp, sclrp2,       & ! Intent(in)
                             sclrprtp, sclrpthlp,          & ! Intent(in)
                             .false.,                      & ! Intent(in)
                             rtm_loc,                      & ! Intent(i/o)
                             rcm, cloud_frac,              & ! Intent(out)
                             ice_supersat_frac, wprcp,     & ! Intent(out)
                             sigma_sqd_w, wpthvp, wp2thvp, & ! Intent(out)
                             rtpthvp, thlpthvp, rc_coef,   & ! Intent(out)
                             rcm_in_layer, cloud_cover,    & ! Intent(out)
                             rcp2_zt, thlprcp, rc_coef_zm, & ! Intent(out)
                             rtm_frz, thlm_frz, sclrpthvp, & ! Intent(out)
                             wp4, wp2rtp, wprtp2, wp2thlp, & ! Intent(out)
                             wpthlp2, wprtpthlp, wp2rcp,   & ! Intent(out)
                             rtprcp, rcp2,                 & ! Intent(out)
                             uprcp, vprcp,                 & ! Intent(out)
                             skw_velocity,                 & ! Intent(out)
                             cloud_frac_zm,                & ! Intent(out)
                             ice_supersat_frac_zm,         & ! Intent(out)
                             rtm_zm, thlm_zm, rcm_zm,      & ! Intent(out)
                             rcm_supersat_adj,             & ! Intent(out)
                             wp2sclrp, wpsclrp2, sclrprcp, & ! Intent(out)
                             wpsclrprtp, wpsclrpthlp,      & ! Intent(out)
                             pdf_params, pdf_params_frz,   & ! Intent(out)
                             pdf_params_zm,                & ! Intent(out)
                             pdf_ict )                       ! Intent(out)

    err_out = err_code
    err_code = clubb_no_error

    rtm_out = rtm_loc
    outs(:, 1) = rcm
    outs(:, 2) = cloud_frac
    outs(:, 3) = ice_supersat_frac
    outs(:, 4) = wprcp
    outs(:, 5) = sigma_sqd_w
    outs(:, 6) = wpthvp
    outs(:, 7) = wp2thvp
    outs(:, 8) = rtpthvp
    outs(:, 9) = thlpthvp
    outs(:, 10) = rc_coef
    outs(:, 11) = rcm_in_layer
    outs(:, 12) = cloud_cover
    outs(:, 13) = rcp2_zt
    outs(:, 14) = thlprcp
    outs(:, 15) = rc_coef_zm
    outs(:, 16) = wp2rtp
    outs(:, 17) = wp2thlp
    outs(:, 18) = wp2rcp
    outs(:, 19) = rtprcp
    outs(:, 20) = rcp2
    outs(:, 21) = uprcp
    outs(:, 22) = vprcp
    outs(:, 23) = cloud_frac_zm
    outs(:, 24) = ice_supersat_frac_zm
    outs(:, 25) = rtm_zm
    outs(:, 26) = thlm_zm
    outs(:, 27) = rcm_zm
    outs(:, 28) = rcm_supersat_adj
    outs(:, 29) = sigma_sqd_w_zt

    call pack_pdf_params(nz, pdf_params, pdfp_zt)
    call pack_pdf_params(nz, pdf_params_zm, pdfp_zm)
  end subroutine drv_pdf_closure_driver

  !---------------------------------------------------------------------
  ! The REAL (public) advance_clubb_core, driven exactly as EAM's
  ! clubb_tend_cam does (see clubb_intr.F90; edsclrm_forcing = 0 and
  ! wpedsclrp_sfc = 0 as in EAM, pert pointers unassociated).  Used to
  ! validate the clubb_pdf_extract transcription END-TO-END: under
  ! ipdf_call_placement = 2 (EAMv3) the final pdf_closure_driver call
  ! inside advance_clubb_core produces outputs that pass through to
  ! the caller untouched, and its inputs are exactly the advanced
  ! prognostic state returned here -- so replaying
  ! drv_pdf_closure_driver on prog_out must reproduce them bitwise.
  !
  ! prog(:,23) slots (in and out):
  !   1 um   2 vm   3 upwp  4 vpwp  5 up2  6 vp2  7 thlm  8 rtm
  !   9 wprtp  10 wpthlp  11 wp2  12 wp3  13 rtp2  14 rtp3
  !  15 thlp2  16 thlp3  17 rtpthlp  18 rcm  19 cloud_frac
  !  20 wpthvp  21 wp2thvp  22 rtpthvp  23 thlpthvp
  ! diag(:,8) slots:
  !   1 khzm  2 khzt  3 qclvar(=rcp2_zt)  4 thlprcp_out  5 wprcp
  !   6 ice_supersat_frac  7 rcm_in_layer  8 cloud_cover
  subroutine drv_advance_clubb_core(nz, ned, dt, fcor, sfc_elevation, &
      wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, host_dx, host_dy, &
      thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
      wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
      rtpthlp_forcing, wm_zm, wm_zt, p_in_pa, rho_zm, rho, exner, &
      rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, &
      thv_ds_zm, thv_ds_zt, rfrzm, radf, prog_in, edsclrm_in, &
      prog_out, edsclrm_out, diag_out, pdfp_zt, pdfp_zm, err_out)
    use clubb_api_module, only: advance_clubb_core_api
    use clubb_driver_helpers, only: pack_pdf_params
    use pdf_parameter_module, only: pdf_parameter, init_pdf_params
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: nz, ned
    real(r8), intent(in) :: dt, fcor, sfc_elevation, wpthlp_sfc, &
        wprtp_sfc, upwp_sfc, vpwp_sfc, host_dx, host_dy
    real(r8), intent(in), dimension(nz) :: thlm_forcing, rtm_forcing, &
        um_forcing, vm_forcing, wprtp_forcing, wpthlp_forcing, &
        rtp2_forcing, thlp2_forcing, rtpthlp_forcing, wm_zm, wm_zt, &
        p_in_pa, rho_zm, rho, exner, rho_ds_zm, rho_ds_zt, &
        invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
        rfrzm, radf
    real(r8), intent(in) :: prog_in(nz, 23)
    real(r8), intent(in) :: edsclrm_in(nz, ned)
    real(r8), intent(out) :: prog_out(nz, 23)
    real(r8), intent(out) :: edsclrm_out(nz, ned)
    real(r8), intent(out) :: diag_out(nz, 8)
    real(r8), intent(out) :: pdfp_zt(nz, 47), pdfp_zm(nz, 47)
    integer, intent(out) :: err_out

    real(r8), dimension(nz) :: um, vm, upwp, vpwp, up2, vp2, thlm, &
        rtm, wprtp, wpthlp, wp2, wp3, rtp2, rtp3, thlp2, thlp3, &
        rtpthlp, rcm, cloud_frac, wpthvp, wp2thvp, rtpthvp, thlpthvp
    real(r8), dimension(nz) :: khzm, khzt, qclvar, thlprcp_out, &
        wprcp, ice_supersat_frac, rcm_in_layer, cloud_cover
    real(r8), dimension(nz, ned) :: edsclrm, edsclrm_forcing
    real(r8), dimension(ned) :: wpedsclrp_sfc
    real(r8), dimension(nz, 0) :: hydromet, wphydrometp, wp2hmp, &
        rtphmp_zt, thlphmp_zt, sclrm, wpsclrp, sclrp2, sclrprtp, &
        sclrpthlp, sclrm_forcing, sclrpthvp
    real(r8), dimension(0) :: wpsclrp_sfc
    type(pdf_parameter) :: pdf_params, pdf_params_zm
    integer :: err_code_api
    real(r8), pointer :: upwp_sfc_pert => null(), vpwp_sfc_pert => null()
    real(r8), pointer, dimension(:) :: um_pert => null(), &
        vm_pert => null(), upwp_pert => null(), vpwp_pert => null()

    call init_pdf_params(nz, pdf_params)
    call init_pdf_params(nz, pdf_params_zm)

    um = prog_in(:, 1)
    vm = prog_in(:, 2)
    upwp = prog_in(:, 3)
    vpwp = prog_in(:, 4)
    up2 = prog_in(:, 5)
    vp2 = prog_in(:, 6)
    thlm = prog_in(:, 7)
    rtm = prog_in(:, 8)
    wprtp = prog_in(:, 9)
    wpthlp = prog_in(:, 10)
    wp2 = prog_in(:, 11)
    wp3 = prog_in(:, 12)
    rtp2 = prog_in(:, 13)
    rtp3 = prog_in(:, 14)
    thlp2 = prog_in(:, 15)
    thlp3 = prog_in(:, 16)
    rtpthlp = prog_in(:, 17)
    rcm = prog_in(:, 18)
    cloud_frac = prog_in(:, 19)
    wpthvp = prog_in(:, 20)
    wp2thvp = prog_in(:, 21)
    rtpthvp = prog_in(:, 22)
    thlpthvp = prog_in(:, 23)
    edsclrm = edsclrm_in
    edsclrm_forcing = 0._r8
    wpedsclrp_sfc = 0._r8

    err_code = clubb_no_error
    err_code_api = clubb_no_error

    call advance_clubb_core_api( &
        .true., dt, fcor, sfc_elevation, 0,                & ! intent(in)
        thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! intent(in)
        sclrm_forcing, edsclrm_forcing, wprtp_forcing,     & ! intent(in)
        wpthlp_forcing, rtp2_forcing, thlp2_forcing,       & ! intent(in)
        rtpthlp_forcing, wm_zm, wm_zt,                     & ! intent(in)
        wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,         & ! intent(in)
        wpsclrp_sfc, wpedsclrp_sfc,                        & ! intent(in)
        p_in_pa, rho_zm, rho, exner,                       & ! intent(in)
        rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,             & ! intent(in)
        invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet,   & ! intent(in)
        rfrzm, radf,                                       & ! intent(in)
        wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt,        & ! intent(in)
        host_dx, host_dy,                                  & ! intent(in)
        um, vm, upwp, vpwp, up2, vp2,                      & ! intent(inout)
        thlm, rtm, wprtp, wpthlp,                          & ! intent(inout)
        wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp,       & ! intent(inout)
        sclrm,                                             & ! intent(inout)
        sclrp2, sclrprtp, sclrpthlp,                       & ! intent(inout)
        wpsclrp, edsclrm, err_code_api,                    & ! intent(inout)
        rcm, cloud_frac,                                   & ! intent(inout)
        wpthvp, wp2thvp, rtpthvp, thlpthvp,                & ! intent(inout)
        sclrpthvp,                                         & ! intent(inout)
        pdf_params, pdf_params_zm,                         & ! intent(inout)
        khzm, khzt, qclvar, thlprcp_out,                   & ! intent(out)
        wprcp, ice_supersat_frac,                          & ! intent(out)
        rcm_in_layer, cloud_cover,                         & ! intent(out)
        upwp_sfc_pert, vpwp_sfc_pert,                      & ! intent(in)
        um_pert, vm_pert, upwp_pert, vpwp_pert )             ! intent(inout)

    err_out = err_code_api
    err_code = clubb_no_error

    prog_out(:, 1) = um
    prog_out(:, 2) = vm
    prog_out(:, 3) = upwp
    prog_out(:, 4) = vpwp
    prog_out(:, 5) = up2
    prog_out(:, 6) = vp2
    prog_out(:, 7) = thlm
    prog_out(:, 8) = rtm
    prog_out(:, 9) = wprtp
    prog_out(:, 10) = wpthlp
    prog_out(:, 11) = wp2
    prog_out(:, 12) = wp3
    prog_out(:, 13) = rtp2
    prog_out(:, 14) = rtp3
    prog_out(:, 15) = thlp2
    prog_out(:, 16) = thlp3
    prog_out(:, 17) = rtpthlp
    prog_out(:, 18) = rcm
    prog_out(:, 19) = cloud_frac
    prog_out(:, 20) = wpthvp
    prog_out(:, 21) = wp2thvp
    prog_out(:, 22) = rtpthvp
    prog_out(:, 23) = thlpthvp
    edsclrm_out = edsclrm

    diag_out(:, 1) = khzm
    diag_out(:, 2) = khzt
    diag_out(:, 3) = qclvar
    diag_out(:, 4) = thlprcp_out
    diag_out(:, 5) = wprcp
    diag_out(:, 6) = ice_supersat_frac
    diag_out(:, 7) = rcm_in_layer
    diag_out(:, 8) = cloud_cover

    call pack_pdf_params(nz, pdf_params, pdfp_zt)
    call pack_pdf_params(nz, pdf_params_zm, pdfp_zm)
  end subroutine drv_advance_clubb_core

  !---------------------------------------------------------------------
  ! Slice C getters ----------------------------------------------------
  ! 1-based indices into params for the additional tunables
  ! advance_xp2_xpyp reads (beta is already in drv_param_indices).
  subroutine drv_param_indices_xp2(idx)
    use parameter_indices, only: iC4, iC5, iC14, ic_K2, ic_K9, inu2, inu9
    integer, intent(out) :: idx(7)
    idx = (/ iC4, iC5, iC14, ic_K2, ic_K9, inu2, inu9 /)
  end subroutine drv_param_indices_xp2

  ! The setup_parameters-derived background-diffusivity profiles
  ! (l_adj_low_res_nu=.true. compile-time: nu * scalar mult_factor from
  ! the average grid spacing) and the model_flags / sponge-settings
  ! configuration advance_xp2_xpyp + clip_covars_denom read.  iflags:
  !   1 l_single_C2_Skw   2 l_explicit_turbulent_adv_xpyp
  !   3 l_upwind_xpyp_ta  4 l_min_xp2_from_corr_wx  5 l_C2_cloud_frac
  !   6 l_hole_fill       7 l_tke_aniso
  !   8 up2_vp2_sponge_damp_settings%l_sponge_damping (never assigned
  !     in an EAM build -- static storage, gfortran zero-init = F)
  subroutine drv_xp2_config(nz, nu2_out, nu9_out, iflags)
    use grid_class, only: gr
    use parameters_tunable, only: nu2_vert_res_dep, nu9_vert_res_dep
    use model_flags, only: l_single_C2_Skw, &
        l_explicit_turbulent_adv_xpyp, l_upwind_xpyp_ta, &
        l_min_xp2_from_corr_wx, l_C2_cloud_frac, l_hole_fill, l_tke_aniso
    use sponge_layer_damping, only: up2_vp2_sponge_damp_settings
    integer, intent(in) :: nz
    real(r8), intent(out) :: nu2_out(nz), nu9_out(nz)
    integer, intent(out) :: iflags(8)
    if (nz /= gr%nz) stop 'drv_xp2_config: nz mismatch'
    nu2_out = nu2_vert_res_dep
    nu9_out = nu9_vert_res_dep
    iflags = 0
    if (l_single_C2_Skw) iflags(1) = 1
    if (l_explicit_turbulent_adv_xpyp) iflags(2) = 1
    if (l_upwind_xpyp_ta) iflags(3) = 1
    if (l_min_xp2_from_corr_wx) iflags(4) = 1
    if (l_C2_cloud_frac) iflags(5) = 1
    if (l_hole_fill) iflags(6) = 1
    if (l_tke_aniso) iflags(7) = 1
    if (up2_vp2_sponge_damp_settings%l_sponge_damping) iflags(8) = 1
  end subroutine drv_xp2_config

  !---------------------------------------------------------------------
  ! The REAL (public) advance_xp2_xpyp, exactly as advance_clubb_core
  ! calls it under the EAMv3 configuration: sclr_dim = 0 (all scalar
  ! arguments are (nz,0)), l_iter = l_iter_xp2_xpyp = .true. (a
  ! compile-time parameter in advance_clubb_core_module),
  ! pdf_implicit_coefs_terms untouched on the compile-time iiPDF_ADG1
  ! path (zero-filled here to keep it a deterministic local).
  ! wprtp2/wpthlp2/wprtpthlp are unused under
  ! l_explicit_turbulent_adv_xpyp=.false. but are real arguments.
  subroutine drv_advance_xp2_xpyp(nz, dt, l_iter_in, tau_zm, wm_zm, &
      rtm, wprtp, thlm, wpthlp, wpthvp, um, vm, wp2, wp2_zt, wp3, &
      upwp, vpwp, sigma_sqd_w, skw_zm, wprtp2, wpthlp2, wprtpthlp, &
      kh_zt, rtp2_forcing, thlp2_forcing, rtpthlp_forcing, rho_ds_zm, &
      rho_ds_zt, invrs_rho_ds_zm, thv_ds_zm, cloud_frac, lscale, &
      wp3_on_wp2, wp3_on_wp2_zt, wp2_splat, rtp2_in, thlp2_in, &
      rtpthlp_in, up2_in, vp2_in, rtp2_out, thlp2_out, rtpthlp_out, &
      up2_out, vp2_out, err_out)
    use advance_xp2_xpyp_module, only: advance_xp2_xpyp
    use pdf_parameter_module, only: implicit_coefs_terms
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: nz
    real(r8), intent(in) :: dt
    logical, intent(in) :: l_iter_in
    real(r8), intent(in), dimension(nz) :: tau_zm, wm_zm, rtm, wprtp, &
        thlm, wpthlp, wpthvp, um, vm, wp2, wp2_zt, wp3, upwp, vpwp, &
        sigma_sqd_w, skw_zm, wprtp2, wpthlp2, wprtpthlp, kh_zt, &
        rtp2_forcing, thlp2_forcing, rtpthlp_forcing, rho_ds_zm, &
        rho_ds_zt, invrs_rho_ds_zm, thv_ds_zm, cloud_frac, lscale, &
        wp3_on_wp2, wp3_on_wp2_zt, wp2_splat
    real(r8), intent(in), dimension(nz) :: rtp2_in, thlp2_in, &
        rtpthlp_in, up2_in, vp2_in
    real(r8), intent(out), dimension(nz) :: rtp2_out, thlp2_out, &
        rtpthlp_out, up2_out, vp2_out
    integer, intent(out) :: err_out

    real(r8), dimension(nz) :: rtp2, thlp2, rtpthlp, up2, vp2
    real(r8), dimension(nz, 0) :: sclrm, wpsclrp, wpsclrp2, &
        wpsclrprtp, wpsclrpthlp, sclrp2, sclrprtp, sclrpthlp
    type(implicit_coefs_terms), dimension(nz) :: pdf_ict
    real(r8), parameter :: z = 0.0_r8

    pdf_ict = implicit_coefs_terms(z, z, z, z, z, z, z, z, z)

    rtp2 = rtp2_in
    thlp2 = thlp2_in
    rtpthlp = rtpthlp_in
    up2 = up2_in
    vp2 = vp2_in
    err_code = clubb_no_error

    call advance_xp2_xpyp( tau_zm, wm_zm, rtm, wprtp, thlm,        & ! In
                           wpthlp, wpthvp, um, vm, wp2, wp2_zt,    & ! In
                           wp3, upwp, vpwp, sigma_sqd_w, skw_zm,   & ! In
                           wprtp2, wpthlp2, wprtpthlp,             & ! In
                           kh_zt, rtp2_forcing, thlp2_forcing,     & ! In
                           rtpthlp_forcing, rho_ds_zm, rho_ds_zt,  & ! In
                           invrs_rho_ds_zm, thv_ds_zm, cloud_frac, & ! In
                           lscale, wp3_on_wp2, wp3_on_wp2_zt,      & ! In
                           pdf_ict,                                & ! In
                           l_iter_in, dt,                          & ! In
                           sclrm, wpsclrp,                         & ! In
                           wpsclrp2, wpsclrprtp, wpsclrpthlp,      & ! In
                           wp2_splat,                              & ! In
                           rtp2, thlp2, rtpthlp, up2, vp2,         & ! Inout
                           sclrp2, sclrprtp, sclrpthlp)              ! Inout

    err_out = err_code
    err_code = clubb_no_error

    rtp2_out = rtp2
    thlp2_out = thlp2
    rtpthlp_out = rtpthlp
    up2_out = up2
    vp2_out = vp2
  end subroutine drv_advance_xp2_xpyp

  !---------------------------------------------------------------------
  ! The REAL clip_covars_denom, exactly as advance_clubb_core calls it
  ! after advance_xp2_xpyp (cl_num = 2/2/2/1/1 under EAMv3
  ! l_predict_upwp_vpwp=F) and after advance_wp2_wp3 (3/3/3/2/2).
  ! sclr_dim = 0; the pert pointers stay unassociated (as in EAM).
  subroutine drv_clip_covars_denom(nz, dt, rtp2, thlp2, up2, vp2, &
      wp2, wprtp_cl_num, wpthlp_cl_num, wpsclrp_cl_num, upwp_cl_num, &
      vpwp_cl_num, wprtp_in, wpthlp_in, upwp_in, vpwp_in, &
      wprtp_out, wpthlp_out, upwp_out, vpwp_out, err_out)
    use clip_explicit, only: clip_covars_denom
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: nz
    real(r8), intent(in) :: dt
    real(r8), intent(in), dimension(nz) :: rtp2, thlp2, up2, vp2, wp2
    integer, intent(in) :: wprtp_cl_num, wpthlp_cl_num, &
        wpsclrp_cl_num, upwp_cl_num, vpwp_cl_num
    real(r8), intent(in), dimension(nz) :: wprtp_in, wpthlp_in, &
        upwp_in, vpwp_in
    real(r8), intent(out), dimension(nz) :: wprtp_out, wpthlp_out, &
        upwp_out, vpwp_out
    integer, intent(out) :: err_out

    real(r8), dimension(nz) :: wprtp, wpthlp, upwp, vpwp
    real(r8), dimension(nz, 0) :: sclrp2, wpsclrp
    real(r8), pointer, dimension(:) :: upwp_pert => null(), &
        vpwp_pert => null()

    wprtp = wprtp_in
    wpthlp = wpthlp_in
    upwp = upwp_in
    vpwp = vpwp_in
    err_code = clubb_no_error

    call clip_covars_denom( dt, rtp2, thlp2, up2, vp2, wp2,           & ! in
                            sclrp2, wprtp_cl_num, wpthlp_cl_num,      & ! in
                            wpsclrp_cl_num, upwp_cl_num, vpwp_cl_num, & ! in
                            wprtp, wpthlp, upwp, vpwp, wpsclrp,       & ! i/o
                            upwp_pert, vpwp_pert )                      ! i/o

    err_out = err_code
    err_code = clubb_no_error

    wprtp_out = wprtp
    wpthlp_out = wpthlp
    upwp_out = upwp
    vpwp_out = vpwp
  end subroutine drv_clip_covars_denom

  !---------------------------------------------------------------------
  ! Slice D getters ----------------------------------------------------
  ! 1-based indices into params for the additional tunables the
  ! Lscale/tau infrastructure reads.
  subroutine drv_param_indices_lscale(idx)
    use parameter_indices, only: ic_K, itaumin, itaumax, &
        iLscale_mu_coef, iLscale_pert_coef, ilmin_coef, iC_wp2_splat, &
        ilambda0_stability_coef, iup2_vp2_factor
    integer, intent(out) :: idx(9)
    idx = (/ ic_K, itaumin, itaumax, iLscale_mu_coef, &
             iLscale_pert_coef, ilmin_coef, iC_wp2_splat, &
             ilambda0_stability_coef, iup2_vp2_factor /)
  end subroutine drv_param_indices_lscale

  ! Derived/module state the slice-D routines read:
  !   lmin  = lmin_coef * 40 m (setup_parameters; NOT in params vector)
  !   T0    = theta0 from setup_clubb_core (parameters_model)
  ! iflags:
  !   1 l_stability_correct_tau_zm   2 l_diag_Lscale_from_tau
  !   3 l_use_C7_Richardson          4 l_use_C11_Richardson
  !   5 l_use_wp3_pr3                6 l_Lscale_plume_centered
  !   7 l_use_ice_latent             8 l_brunt_vaisala_freq_moist
  !   9 l_use_thvm_in_bv_freq       10 l_sat_mixrat_lookup
  !  11 l_tke_aniso
  subroutine drv_lscale_config(lmin_out, t0_out, iflags)
    use parameters_tunable, only: lmin
    use parameters_model, only: T0
    use model_flags, only: l_stability_correct_tau_zm, &
        l_diag_Lscale_from_tau, l_use_C7_Richardson, &
        l_use_C11_Richardson, l_use_wp3_pr3, l_Lscale_plume_centered, &
        l_use_ice_latent, l_brunt_vaisala_freq_moist, &
        l_use_thvm_in_bv_freq, l_sat_mixrat_lookup, l_tke_aniso
    real(r8), intent(out) :: lmin_out, t0_out
    integer, intent(out) :: iflags(11)
    lmin_out = lmin
    t0_out = T0
    iflags = 0
    if (l_stability_correct_tau_zm) iflags(1) = 1
    if (l_diag_Lscale_from_tau) iflags(2) = 1
    if (l_use_C7_Richardson) iflags(3) = 1
    if (l_use_C11_Richardson) iflags(4) = 1
    if (l_use_wp3_pr3) iflags(5) = 1
    if (l_Lscale_plume_centered) iflags(6) = 1
    if (l_use_ice_latent) iflags(7) = 1
    if (l_brunt_vaisala_freq_moist) iflags(8) = 1
    if (l_use_thvm_in_bv_freq) iflags(9) = 1
    if (l_sat_mixrat_lookup) iflags(10) = 1
    if (l_tke_aniso) iflags(11) = 1
  end subroutine drv_lscale_config

  ! Toggle the Brunt-Vaisala variant flags (module variables) so all
  ! three formula branches can be goldened; the generator restores the
  ! EAMv3 defaults (F, F) afterwards.
  subroutine drv_set_bv_flags(l_moist, l_use_thvm)
    use model_flags, only: l_brunt_vaisala_freq_moist, &
        l_use_thvm_in_bv_freq
    logical, intent(in) :: l_moist, l_use_thvm
    l_brunt_vaisala_freq_moist = l_moist
    l_use_thvm_in_bv_freq = l_use_thvm
  end subroutine drv_set_bv_flags

  !---------------------------------------------------------------------
  ! The REAL (public) compute_mixing_length, exactly as
  ! advance_clubb_core calls it (l_implemented=.true.; mu = the
  ! module tunable newmu -- passed explicitly here so entrainment-rate
  ! sweeps can be goldened).
  subroutine drv_compute_mixing_length(nz, thvm, thlm, rtm, em, &
      lscale_max_in, mu_in, p_in_pa, exner, thv_ds, &
      lscale, lscale_up, lscale_down, err_out)
    use mixing_length, only: compute_mixing_length
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: nz
    real(r8), intent(in), dimension(nz) :: thvm, thlm, rtm, em, &
        p_in_pa, exner, thv_ds
    real(r8), intent(in) :: lscale_max_in, mu_in
    real(r8), intent(out), dimension(nz) :: lscale, lscale_up, &
        lscale_down
    integer, intent(out) :: err_out
    err_code = clubb_no_error
    call compute_mixing_length( thvm, thlm,                        & ! in
                         rtm, em, lscale_max_in, p_in_pa,          & ! in
                         exner, thv_ds, mu_in, .true.,             & ! in
                         lscale, lscale_up, lscale_down )            ! out
    err_out = err_code
    err_code = clubb_no_error
  end subroutine drv_compute_mixing_length

  !---------------------------------------------------------------------
  ! The REAL (public) calc_brunt_vaisala_freq_sqd under the CURRENT
  ! model_flags (see drv_set_bv_flags for the variants).
  subroutine drv_brunt_vaisala(nz, thlm, exner, rtm, rcm, p_in_pa, &
                               thvm, bv_out)
    use advance_helper_module, only: calc_brunt_vaisala_freq_sqd
    integer, intent(in) :: nz
    real(r8), intent(in), dimension(nz) :: thlm, exner, rtm, rcm, &
        p_in_pa, thvm
    real(r8), intent(out) :: bv_out(nz)
    call calc_brunt_vaisala_freq_sqd( thlm, exner, rtm, rcm, &
                                      p_in_pa, thvm, bv_out )
  end subroutine drv_brunt_vaisala

  ! The REAL (public) calc_stability_correction (uses the EAMv3 dry-T0
  ! Brunt-Vaisala branch and lambda0_stability_coef).
  subroutine drv_stability_correction(nz, thlm, lscale, em, exner, &
      rtm, rcm, p_in_pa, thvm, sc_out)
    use advance_helper_module, only: calc_stability_correction
    integer, intent(in) :: nz
    real(r8), intent(in), dimension(nz) :: thlm, lscale, em, exner, &
        rtm, rcm, p_in_pa, thvm
    real(r8), intent(out) :: sc_out(nz)
    sc_out = calc_stability_correction( thlm, lscale, em, exner, &
                                        rtm, rcm, p_in_pa, thvm )
  end subroutine drv_stability_correction

  ! The REAL (public) term_wp2_splat + term_wp3_splat.  C_wp2_splat is
  ! an explicit argument of both routines (EAMv3 passes the tunable,
  ! default 0.0 -- a nonzero sweep goldens the clip branch).
  subroutine drv_term_splat(nz, dt, c_wp2_splat_in, wp2, wp2_zt, wp3, &
                            tau_zm, tau_zt, wp2_splat, wp3_splat)
    use advance_helper_module, only: term_wp2_splat, term_wp3_splat
    integer, intent(in) :: nz
    real(r8), intent(in) :: dt, c_wp2_splat_in
    real(r8), intent(in), dimension(nz) :: wp2, wp2_zt, wp3, tau_zm, &
        tau_zt
    real(r8), intent(out), dimension(nz) :: wp2_splat, wp3_splat
    call term_wp2_splat( c_wp2_splat_in, nz, dt, wp2, wp2_zt, tau_zm, &
                         wp2_splat )
    call term_wp3_splat( c_wp2_splat_in, nz, dt, wp2, wp3, tau_zt, &
                         wp3_splat )
  end subroutine drv_term_splat

  !---------------------------------------------------------------------
  ! The REAL (public) calc_surface_varnce, exactly as
  ! advance_clubb_core calls it (sclr_dim = 0; l_andre_1978 = .false.
  ! is compile-time in surface_varnce_module).
  ! outs slots: 1 wp2  2 up2  3 vp2  4 thlp2  5 rtp2  6 rtpthlp
  subroutine drv_calc_surface_varnce(upwp_sfc, vpwp_sfc, wpthlp_sfc, &
      wprtp_sfc, um_sfc, vm_sfc, lscale_up_sfc, wp2_splat_sfc, &
      tau_zm_sfc, outs, err_out)
    use surface_varnce_module, only: calc_surface_varnce
    use error_code, only: err_code, clubb_no_error
    real(r8), intent(in) :: upwp_sfc, vpwp_sfc, wpthlp_sfc, &
        wprtp_sfc, um_sfc, vm_sfc, lscale_up_sfc, wp2_splat_sfc, &
        tau_zm_sfc
    real(r8), intent(out) :: outs(6)
    integer, intent(out) :: err_out
    real(r8), dimension(0) :: wpsclrp_sfc, sclrp2_sfc, sclrprtp_sfc, &
        sclrpthlp_sfc
    err_code = clubb_no_error
    call calc_surface_varnce( upwp_sfc, vpwp_sfc, wpthlp_sfc,       & ! in
                         wprtp_sfc, um_sfc, vm_sfc, lscale_up_sfc,  & ! in
                         wpsclrp_sfc, wp2_splat_sfc, tau_zm_sfc,    & ! in
                         outs(1), outs(2), outs(3),                 & ! out
                         outs(4), outs(5), outs(6),                 & ! out
                         sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc )    ! out
    err_out = err_code
    err_code = clubb_no_error
  end subroutine drv_calc_surface_varnce

  !---------------------------------------------------------------------
  ! VERBATIM TRANSCRIPTION of the inline Lscale/tau segment of
  ! advance_clubb_core (advance_clubb_core_module.F90, EAMv3 path):
  !   wp2_zt (line ~1006), thvm (~1083), em l_tke_aniso=T (~1094),
  !   sqrt_em_zt (~1097), the single compute_mixing_length call
  !   (~1244; l_avg_Lscale=.false. is a COMPILE-TIME parameter, so the
  !   perturbed-Lscale calls and the averaging are dead), tau_zt/tau_zm
  !   (~1271), Kh_zt/Kh_zm (~1336), term_wp2_splat/term_wp3_splat
  !   (~1346), the surface-variance branch (~1360; sfc_elevation
  !   comparison verbatim), calc_stability_correction (~1483),
  !   tau_N2_zm = tau_zm/stability_correction under
  !   l_stability_correct_tau_zm=.true. (~1497), and
  !   Cx_fnc_Richardson = 0 (~1515; l_use_C7_Richardson =
  !   l_use_C11_Richardson = l_use_wp3_pr3 = .false.).
  !   set_Lscale_max is private; its l_implemented=T body is one line,
  !   transcribed at the top.
  ! Validated end-to-end: Kh_zm/Kh_zt here must match the khzm/khzt
  ! diagnostics of the REAL drv_advance_clubb_core on the same entry
  ! state BITWISE (asserted in gen_clubb_golden.py).
  !
  ! seg(:,14) slots:
  !   1 em        2 thvm      3 sqrt_em_zt  4 Lscale  5 Lscale_up
  !   6 Lscale_down  7 tau_zt  8 tau_zm  9 Kh_zt  10 Kh_zm
  !  11 wp2_splat  12 wp3_splat  13 stability_correction  14 tau_N2_zm
  ! sfc(6): wp2(1), up2(1), vp2(1), thlp2(1), rtp2(1), rtpthlp(1)
  ! after the surface-variance update.
  subroutine drv_lscale_tau_segment(nz, dt, sfc_elevation, host_dx, &
      host_dy, wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, thlm, rtm, &
      rcm, wp2, wp3, up2, vp2, um, vm, p_in_pa, exner, thv_ds_zt, &
      seg, sfc, err_out)
    use grid_class, only: gr, zm2zt, zt2zm
    use mixing_length, only: compute_mixing_length
    use advance_helper_module, only: calc_stability_correction, &
        term_wp2_splat, term_wp3_splat
    use surface_varnce_module, only: calc_surface_varnce
    use parameters_tunable, only: mu, c_K, taumax, C_wp2_splat
    use constants_clubb, only: w_tol_sqd, em_min, ep1, ep2, Lv, Cp, &
        zero_threshold, eps, thl_tol, rt_tol
    use model_flags, only: l_stability_correct_tau_zm
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: nz
    real(r8), intent(in) :: dt, sfc_elevation, host_dx, host_dy, &
        wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc
    real(r8), intent(in), dimension(nz) :: thlm, rtm, rcm, wp2, wp3, &
        up2, vp2, um, vm, p_in_pa, exner, thv_ds_zt
    real(r8), intent(out) :: seg(nz, 14)
    real(r8), intent(out) :: sfc(6)
    integer, intent(out) :: err_out

    real(r8), dimension(nz) :: wp2_zt, thvm, em, sqrt_em_zt, &
        Lscale, Lscale_up, Lscale_down, tau_zt, tau_zm, Kh_zt, Kh_zm, &
        wp2_splat, wp3_splat, stability_correction, tau_N2_zm
    real(r8), dimension(0) :: wpsclrp_sfc, sclrp2_sfc, sclrprtp_sfc, &
        sclrpthlp_sfc
    real(r8) :: Lscale_max, wp2_1, up2_1, vp2_1, thlp2_1, rtp2_1, &
        rtpthlp_1

    if (nz /= gr%nz) stop 'drv_lscale_tau_segment: nz mismatch'
    err_code = clubb_no_error

    ! set_Lscale_max, l_implemented=.true. branch
    Lscale_max = 0.25_r8 * min( host_dx, host_dy )

    wp2_zt = max( zm2zt( wp2 ), w_tol_sqd )

    thvm = thlm + ep1 * thv_ds_zt * rtm &
                + ( Lv/(Cp*exner) - ep2 * thv_ds_zt ) * rcm

    ! l_tke_aniso = .true.
    em = 0.5_r8 * ( wp2 + vp2 + up2 )

    sqrt_em_zt = SQRT( MAX( em_min, zm2zt( em ) ) )

    call compute_mixing_length( thvm, thlm,                        & ! in
                         rtm, em, Lscale_max, p_in_pa,             & ! in
                         exner, thv_ds_zt, mu, .true.,             & ! in
                         Lscale, Lscale_up, Lscale_down )            ! out

    tau_zt = MIN( Lscale / sqrt_em_zt, taumax )
    tau_zm = MIN( ( MAX( zt2zm( Lscale ), zero_threshold )  &
                   / SQRT( MAX( em_min, em ) ) ), taumax )

    Kh_zt = c_K * Lscale * sqrt_em_zt
    Kh_zm = c_K * max( zt2zm( Lscale ), zero_threshold )  &
                * sqrt( max( em, em_min ) )

    call term_wp2_splat( C_wp2_splat, gr%nz, dt, wp2, wp2_zt, tau_zm, &
                         wp2_splat )
    call term_wp3_splat( C_wp2_splat, gr%nz, dt, wp2, wp3, tau_zt, &
                         wp3_splat )

    if ( abs(gr%zm(1)-sfc_elevation) &
         <= abs(gr%zm(1)+sfc_elevation)*eps/2 ) then
      call calc_surface_varnce( upwp_sfc, vpwp_sfc, wpthlp_sfc,     & ! in
                           wprtp_sfc, um(2), vm(2), Lscale_up(2),   & ! in
                           wpsclrp_sfc, wp2_splat(1), tau_zm(1),    & ! in
                           wp2_1, up2_1, vp2_1,                     & ! out
                           thlp2_1, rtp2_1, rtpthlp_1,              & ! out
                           sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc )  ! out
    else
      wp2_1     = w_tol_sqd
      up2_1     = w_tol_sqd
      vp2_1     = w_tol_sqd
      thlp2_1   = thl_tol**2
      rtp2_1    = rt_tol**2
      rtpthlp_1 = 0.0_r8
    end if

    stability_correction = calc_stability_correction( thlm, Lscale, &
        em, exner, rtm, rcm, p_in_pa, thvm )

    if ( l_stability_correct_tau_zm ) then
      tau_N2_zm = tau_zm / stability_correction
    else
      tau_N2_zm = -999._r8   ! unused_var
    end if

    err_out = err_code
    err_code = clubb_no_error

    seg(:, 1) = em
    seg(:, 2) = thvm
    seg(:, 3) = sqrt_em_zt
    seg(:, 4) = Lscale
    seg(:, 5) = Lscale_up
    seg(:, 6) = Lscale_down
    seg(:, 7) = tau_zt
    seg(:, 8) = tau_zm
    seg(:, 9) = Kh_zt
    seg(:, 10) = Kh_zm
    seg(:, 11) = wp2_splat
    seg(:, 12) = wp3_splat
    seg(:, 13) = stability_correction
    seg(:, 14) = tau_N2_zm
    sfc(1) = wp2_1
    sfc(2) = up2_1
    sfc(3) = vp2_1
    sfc(4) = thlp2_1
    sfc(5) = rtp2_1
    sfc(6) = rtpthlp_1
  end subroutine drv_lscale_tau_segment

  !---------------------------------------------------------------------
  ! Slice E: advance_xm_wpxp ---------------------------------------------
  ! 1-based indices into params for the additional tunables
  ! advance_xm_wpxp reads (C6rt/C6rtb/C6rtc/C6thl*/C7/C7b/wpxp_L_thresh/
  ! altitude_threshold are already in drv_param_indices).
  subroutine drv_param_indices_xm_wpxp(idx)
    use parameter_indices, only: iC7c, ic_K6, inu6, &
        iC6rt_Lscale0, iC6thl_Lscale0, iC7_Lscale0
    integer, intent(out) :: idx(6)
    idx = (/ iC7c, ic_K6, inu6, iC6rt_Lscale0, iC6thl_Lscale0, &
             iC7_Lscale0 /)
  end subroutine drv_param_indices_xm_wpxp

  ! Model-flag / module-state configuration advance_xm_wpxp reads.
  ! iflags:
  !   1 l_clip_semi_implicit          2 l_explicit_turbulent_adv_wpxp
  !   3 l_upwind_wpxp_ta              4 l_predict_upwp_vpwp
  !   5 l_uv_nudge                    6 l_pos_def
  !   7 l_hole_fill                   8 l_clip_turb_adv
  !   9 l_diffuse_rtm_and_thlm       10 l_stability_correct_Kh_N2_zm
  !  11 rtm_sponge_damp_settings%l_sponge_damping   (zero-init in EAM)
  !  12 thlm_sponge_damp_settings%l_sponge_damping  (zero-init in EAM)
  !  13 l_tke_aniso
  subroutine drv_xm_wpxp_config(nz, nu6_out, iflags)
    use grid_class, only: gr
    use parameters_tunable, only: nu6_vert_res_dep
    use model_flags, only: l_clip_semi_implicit, &
        l_explicit_turbulent_adv_wpxp, l_upwind_wpxp_ta, &
        l_predict_upwp_vpwp, l_uv_nudge, l_pos_def, l_hole_fill, &
        l_clip_turb_adv, l_diffuse_rtm_and_thlm, &
        l_stability_correct_Kh_N2_zm, l_tke_aniso
    use sponge_layer_damping, only: rtm_sponge_damp_settings, &
        thlm_sponge_damp_settings
    integer, intent(in) :: nz
    real(r8), intent(out) :: nu6_out(nz)
    integer, intent(out) :: iflags(13)
    if (nz /= gr%nz) stop 'drv_xm_wpxp_config: nz mismatch'
    nu6_out = nu6_vert_res_dep
    iflags = 0
    if (l_clip_semi_implicit) iflags(1) = 1
    if (l_explicit_turbulent_adv_wpxp) iflags(2) = 1
    if (l_upwind_wpxp_ta) iflags(3) = 1
    if (l_predict_upwp_vpwp) iflags(4) = 1
    if (l_uv_nudge) iflags(5) = 1
    if (l_pos_def) iflags(6) = 1
    if (l_hole_fill) iflags(7) = 1
    if (l_clip_turb_adv) iflags(8) = 1
    if (l_diffuse_rtm_and_thlm) iflags(9) = 1
    if (l_stability_correct_Kh_N2_zm) iflags(10) = 1
    if (rtm_sponge_damp_settings%l_sponge_damping) iflags(11) = 1
    if (thlm_sponge_damp_settings%l_sponge_damping) iflags(12) = 1
    if (l_tke_aniso) iflags(13) = 1
  end subroutine drv_xm_wpxp_config

  !---------------------------------------------------------------------
  ! CLUBB's LAPACK band-diagonal wrapper (band_solve -> dgbsv), exactly
  ! as xm_wpxp_solve calls it (nsup=nsub=2).  lhs is in CLUBB's
  ! row-oriented band storage (nsup+nsub+1, ndim): lhs(3,i) is the main
  ! diagonal of equation i, lhs(1:2,i) couple variables i+2/i+1,
  ! lhs(4:5,i) couple variables i-1/i-2.  rhs copied (dgbsv overwrites).
  subroutine drv_band_solve(ndim, nrhs, lhs, rhs, solution, err_out)
    use lapack_wrap, only: band_solve
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: ndim, nrhs
    real(r8), intent(in) :: lhs(5, ndim)
    real(r8), intent(in) :: rhs(ndim, nrhs)
    real(r8), intent(out) :: solution(ndim, nrhs)
    integer, intent(out) :: err_out
    real(r8) :: rhs_c(ndim, nrhs)
    rhs_c = rhs
    err_code = clubb_no_error
    call band_solve("drv_band", 2, 2, ndim, nrhs, lhs, rhs_c, solution)
    err_out = err_code
    err_code = clubb_no_error
  end subroutine drv_band_solve

  !---------------------------------------------------------------------
  ! The REAL (public) calc_turb_adv_range (mono_flux_limiter.F90):
  ! turbulent-advection level-range for the monotonic flux limiter
  ! (l_constant_thickness=.false. compile-time path, driven by the
  ! zm-level PDF component means/variances).
  subroutine drv_calc_turb_adv_range(nz, dt, w_1_zm, w_2_zm, &
      varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm, low_lev, high_lev)
    use grid_class, only: gr
    use mono_flux_limiter, only: calc_turb_adv_range
    integer, intent(in) :: nz
    real(r8), intent(in) :: dt
    real(r8), intent(in), dimension(nz) :: w_1_zm, w_2_zm, &
        varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm
    integer, intent(out), dimension(nz) :: low_lev, high_lev
    if (nz /= gr%nz) stop 'drv_calc_turb_adv_range: nz mismatch'
    call calc_turb_adv_range( dt, w_1_zm, w_2_zm, varnce_w_1_zm, &
                              varnce_w_2_zm, mixt_frac_zm, &
                              low_lev, high_lev )
  end subroutine drv_calc_turb_adv_range

  !---------------------------------------------------------------------
  ! The REAL (public) monotonic_turbulent_flux_limit, exactly as
  ! xm_wpxp_clipping_and_stats calls it (l_implemented=.true.,
  ! l_mfl_xm_imp_adj=.true. compile-time).  solve_type: 1=thlm, 2=rtm.
  subroutine drv_mono_flux_limiter(nz, solve_type, dt, xm_old, xp2, &
      wm_zt, xm_forcing, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
      invrs_rho_ds_zt, xp2_threshold, xm_tol, low_lev, high_lev, &
      xm_in, wpxp_in, xm_out, wpxp_out, err_out)
    use grid_class, only: gr
    use mono_flux_limiter, only: monotonic_turbulent_flux_limit
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: nz, solve_type
    real(r8), intent(in) :: dt, xp2_threshold, xm_tol
    real(r8), intent(in), dimension(nz) :: xm_old, xp2, wm_zt, &
        xm_forcing, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
        invrs_rho_ds_zt
    integer, intent(in), dimension(nz) :: low_lev, high_lev
    real(r8), intent(in), dimension(nz) :: xm_in, wpxp_in
    real(r8), intent(out), dimension(nz) :: xm_out, wpxp_out
    integer, intent(out) :: err_out
    real(r8), dimension(nz) :: xm, wpxp
    if (nz /= gr%nz) stop 'drv_mono_flux_limiter: nz mismatch'
    xm = xm_in
    wpxp = wpxp_in
    err_code = clubb_no_error
    call monotonic_turbulent_flux_limit( solve_type, dt, xm_old, &
                                         xp2, wm_zt, xm_forcing, &
                                         rho_ds_zm, rho_ds_zt, &
                                         invrs_rho_ds_zm, &
                                         invrs_rho_ds_zt, &
                                         xp2_threshold, .true., &
                                         low_lev, high_lev, &
                                         xm, xm_tol, wpxp )
    err_out = err_code
    err_code = clubb_no_error
    xm_out = xm
    wpxp_out = wpxp
  end subroutine drv_mono_flux_limiter

  !---------------------------------------------------------------------
  ! The REAL (public) fill_holes_vertical on the "zt" grid, exactly as
  ! xm_wpxp_clipping_and_stats calls it (num_draw_pts=2).
  subroutine drv_fill_holes_zt(nz, threshold, rho_ds_zt, rho_ds_zm, &
                               field_in, field_out)
    use grid_class, only: gr
    use fill_holes, only: fill_holes_vertical
    integer, intent(in) :: nz
    real(r8), intent(in) :: threshold
    real(r8), intent(in), dimension(nz) :: rho_ds_zt, rho_ds_zm, &
        field_in
    real(r8), intent(out), dimension(nz) :: field_out
    real(r8), dimension(nz) :: field
    if (nz /= gr%nz) stop 'drv_fill_holes_zt: nz mismatch'
    field = field_in
    call fill_holes_vertical( 2, threshold, "zt", rho_ds_zt, &
                              rho_ds_zm, field )
    field_out = field
  end subroutine drv_fill_holes_zt

  !---------------------------------------------------------------------
  ! The REAL (public) advance_xm_wpxp, exactly as advance_clubb_core
  ! calls it under the EAMv3 configuration: sclr_dim = 0,
  ! l_implemented = .true., Cx_fnc_Richardson = 0 (the l_use_*
  ! Richardson flags are all false), pdf_implicit_coefs_terms
  ! zero-filled (only read on the iiPDF_new path; EAM is compile-time
  ! ADG1), pert pointers unassociated, and the um/vm prediction inputs
  ! zero/dummy (l_predict_upwp_vpwp = .false.: um/vm/upwp/vpwp pass
  ! through untouched -- asserted by the golden generator).
  subroutine drv_advance_xm_wpxp(nz, dt, sigma_sqd_w, wm_zm, wm_zt, &
      wp2, Lscale, em, wp3_on_wp2, wp3_on_wp2_zt, Kh_zt, Kh_zm, &
      tau_C6_zm, Skw_zm, wp2rtp, rtpthvp, rtm_forcing, wprtp_forcing, &
      wp2thlp, thlpthvp, thlm_forcing, wpthlp_forcing, rho_ds_zm, &
      rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, rtp2, &
      thlp2, w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, &
      mixt_frac_zm, exner, rcm, p_in_pa, thvm, &
      rtm_in, wprtp_in, thlm_in, wpthlp_in, &
      rtm_out, wprtp_out, thlm_out, wpthlp_out, err_out)
    use grid_class, only: gr
    use advance_xm_wpxp_module, only: advance_xm_wpxp
    use pdf_parameter_module, only: implicit_coefs_terms
    use error_code, only: err_code, clubb_no_error
    integer, intent(in) :: nz
    real(r8), intent(in) :: dt
    real(r8), intent(in), dimension(nz) :: sigma_sqd_w, wm_zm, wm_zt, &
        wp2, Lscale, em, wp3_on_wp2, wp3_on_wp2_zt, Kh_zt, Kh_zm, &
        tau_C6_zm, Skw_zm, wp2rtp, rtpthvp, rtm_forcing, &
        wprtp_forcing, wp2thlp, thlpthvp, thlm_forcing, &
        wpthlp_forcing, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
        invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm, &
        varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm, exner, rcm, &
        p_in_pa, thvm
    real(r8), intent(in), dimension(nz) :: rtm_in, wprtp_in, thlm_in, &
        wpthlp_in
    real(r8), intent(out), dimension(nz) :: rtm_out, wprtp_out, &
        thlm_out, wpthlp_out
    integer, intent(out) :: err_out

    real(r8), dimension(nz) :: rtm, wprtp, thlm, wpthlp
    real(r8), dimension(nz) :: um, vm, upwp, vpwp
    real(r8), dimension(nz) :: zeros, cx_fnc_richardson
    real(r8), dimension(nz, 0) :: wp2sclrp, sclrpthvp, sclrm_forcing, &
        sclrp2, sclrm, wpsclrp
    type(implicit_coefs_terms), dimension(nz) :: pdf_ict
    real(r8), pointer, dimension(:) :: um_pert => null(), &
        vm_pert => null(), upwp_pert => null(), vpwp_pert => null()
    real(r8), parameter :: z = 0.0_r8

    if (nz /= gr%nz) stop 'drv_advance_xm_wpxp: nz mismatch'
    pdf_ict = implicit_coefs_terms(z, z, z, z, z, z, z, z, z)
    zeros = 0.0_r8
    cx_fnc_richardson = 0.0_r8
    um = 0.0_r8
    vm = 0.0_r8
    upwp = 0.0_r8
    vpwp = 0.0_r8

    rtm = rtm_in
    wprtp = wprtp_in
    thlm = thlm_in
    wpthlp = wpthlp_in
    err_code = clubb_no_error

    call advance_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                          Lscale, wp3_on_wp2, wp3_on_wp2_zt, Kh_zt, &
                          Kh_zm, tau_C6_zm, Skw_zm, wp2rtp, rtpthvp, &
                          rtm_forcing, wprtp_forcing, zeros, wp2thlp, &
                          thlpthvp, thlm_forcing, wpthlp_forcing, &
                          zeros, rho_ds_zm, rho_ds_zt, &
                          invrs_rho_ds_zm, invrs_rho_ds_zt, &
                          thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm, &
                          varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm, &
                          .true., em, wp2sclrp, sclrpthvp, &
                          sclrm_forcing, sclrp2, exner, rcm, p_in_pa, &
                          thvm, cx_fnc_richardson, pdf_ict, &
                          zeros, zeros, zeros, zeros, zeros, &
                          1.0e-4_r8, zeros, zeros, zeros, zeros, &
                          zeros, zeros, zeros, &
                          rtm, wprtp, thlm, wpthlp, &
                          sclrm, wpsclrp, um, upwp, vm, vpwp, &
                          um_pert, vm_pert, upwp_pert, vpwp_pert )

    err_out = err_code
    err_code = clubb_no_error

    if (any(um /= 0.0_r8) .or. any(vm /= 0.0_r8) .or. &
        any(upwp /= 0.0_r8) .or. any(vpwp /= 0.0_r8)) then
      stop 'drv_advance_xm_wpxp: um/vm touched (l_predict_upwp_vpwp?)'
    end if

    rtm_out = rtm
    wprtp_out = wprtp
    thlm_out = thlm
    wpthlp_out = wpthlp
  end subroutine drv_advance_xm_wpxp

end module clubb_driver
