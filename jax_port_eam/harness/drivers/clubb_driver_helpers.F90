! Helpers for the CLUBB f2py driver that take derived-type arguments
! (f2py cannot wrap those, so they live in a separate plain module
! compiled as an ordinary source; clubb_driver use-associates them).
module clubb_driver_helpers
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)
  public :: pack_pdf_params
  private

contains

  ! Pack the 47 pdf_parameter fields into pdfp(nz, 47), in the order of
  ! pdf_parameter_module's init_pdf_params (the same slot order as
  ! drv_pdf_closure in clubb_driver.F90 and PDFP_SLOTS in
  ! tests/test_clubb.py).
  subroutine pack_pdf_params(nz, pp, pdfp)
    use pdf_parameter_module, only: pdf_parameter
    integer, intent(in) :: nz
    type(pdf_parameter), intent(in) :: pp
    real(r8), intent(out) :: pdfp(nz, 47)
    pdfp(:, 1) = pp%w_1
    pdfp(:, 2) = pp%w_2
    pdfp(:, 3) = pp%varnce_w_1
    pdfp(:, 4) = pp%varnce_w_2
    pdfp(:, 5) = pp%rt_1
    pdfp(:, 6) = pp%rt_2
    pdfp(:, 7) = pp%varnce_rt_1
    pdfp(:, 8) = pp%varnce_rt_2
    pdfp(:, 9) = pp%thl_1
    pdfp(:, 10) = pp%thl_2
    pdfp(:, 11) = pp%varnce_thl_1
    pdfp(:, 12) = pp%varnce_thl_2
    pdfp(:, 13) = pp%corr_w_rt_1
    pdfp(:, 14) = pp%corr_w_rt_2
    pdfp(:, 15) = pp%corr_w_thl_1
    pdfp(:, 16) = pp%corr_w_thl_2
    pdfp(:, 17) = pp%corr_rt_thl_1
    pdfp(:, 18) = pp%corr_rt_thl_2
    pdfp(:, 19) = pp%alpha_thl
    pdfp(:, 20) = pp%alpha_rt
    pdfp(:, 21) = pp%crt_1
    pdfp(:, 22) = pp%crt_2
    pdfp(:, 23) = pp%cthl_1
    pdfp(:, 24) = pp%cthl_2
    pdfp(:, 25) = pp%chi_1
    pdfp(:, 26) = pp%chi_2
    pdfp(:, 27) = pp%stdev_chi_1
    pdfp(:, 28) = pp%stdev_chi_2
    pdfp(:, 29) = pp%stdev_eta_1
    pdfp(:, 30) = pp%stdev_eta_2
    pdfp(:, 31) = pp%covar_chi_eta_1
    pdfp(:, 32) = pp%covar_chi_eta_2
    pdfp(:, 33) = pp%corr_w_chi_1
    pdfp(:, 34) = pp%corr_w_chi_2
    pdfp(:, 35) = pp%corr_w_eta_1
    pdfp(:, 36) = pp%corr_w_eta_2
    pdfp(:, 37) = pp%corr_chi_eta_1
    pdfp(:, 38) = pp%corr_chi_eta_2
    pdfp(:, 39) = pp%rsatl_1
    pdfp(:, 40) = pp%rsatl_2
    pdfp(:, 41) = pp%rc_1
    pdfp(:, 42) = pp%rc_2
    pdfp(:, 43) = pp%cloud_frac_1
    pdfp(:, 44) = pp%cloud_frac_2
    pdfp(:, 45) = pp%mixt_frac
    pdfp(:, 46) = pp%ice_supersat_frac_1
    pdfp(:, 47) = pp%ice_supersat_frac_2
  end subroutine pack_pdf_params

end module clubb_driver_helpers
