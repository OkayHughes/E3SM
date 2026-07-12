! f2py driver for dadadj (dry adiabatic adjustment).
module dadadj_driver
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_dadadj(ncol, nlev, nlvdry_in, pmid, pint, pdel, &
                        t_in, q_in, t_out, q_out)
    use ppgrid, only: pcols, pver, pverp
    use cam_control_mod, only: nlvdry
    integer, intent(in) :: ncol, nlev, nlvdry_in
    real(r8), intent(in) :: pmid(ncol, nlev), pint(ncol, nlev + 1)
    real(r8), intent(in) :: pdel(ncol, nlev)
    real(r8), intent(in) :: t_in(ncol, nlev), q_in(ncol, nlev)
    real(r8), intent(out) :: t_out(ncol, nlev), q_out(ncol, nlev)

    real(r8) :: pmid_f(pcols, pver), pint_f(pcols, pverp)
    real(r8) :: pdel_f(pcols, pver), t_f(pcols, pver), q_f(pcols, pver)

    if (nlev /= pver .or. ncol > pcols) stop 'drv_dadadj: bad dims'
    nlvdry = nlvdry_in

    pmid_f = 1.0e5_r8
    pint_f = 1.0e5_r8
    pdel_f = 1.0_r8
    t_f = 300.0_r8
    q_f = 0.0_r8
    pmid_f(1:ncol, :) = pmid
    pint_f(1:ncol, :) = pint
    pdel_f(1:ncol, :) = pdel
    t_f(1:ncol, :) = t_in
    q_f(1:ncol, :) = q_in

    call dadadj(0, ncol, pmid_f, pint_f, pdel_f, t_f, q_f)

    t_out = t_f(1:ncol, :)
    q_out = q_f(1:ncol, :)
  end subroutine drv_dadadj

end module dadadj_driver
