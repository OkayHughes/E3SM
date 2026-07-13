! Thin f2py wrapper around drivers/rrtmgp_core.F90 (which holds the
! derived-type k-dist storage f2py's crackfortran cannot parse). All
! real work (REAL EAM sources plus documented transcriptions of
! radiation.F90 private helpers) lives in rrtmgp_core; see that file.
module rrtmgp_driver
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine w_init(sw_file, lw_file, liq_file, ice_file)
    use rrtmgp_core, only: drv_init, drv_load_kdists
    character(len=256), intent(in) :: sw_file, lw_file, liq_file, ice_file
    call drv_init(sw_file, lw_file, liq_file, ice_file)
    call drv_load_kdists(sw_file, lw_file)
  end subroutine w_init

  subroutine w_dims(nswb, nlwb, nswg, nlwg)
    use rrtmgp_core, only: drv_dims
    integer, intent(out) :: nswb, nlwb, nswg, nlwg
    call drv_dims(nswb, nlwb, nswg, nlwg)
  end subroutine w_dims

  subroutine w_gpt_bands(nswg, nlwg, gpb_sw, gpb_lw)
    use rrtmgp_core, only: drv_gpt_bands
    integer, intent(in) :: nswg, nlwg
    integer, intent(out) :: gpb_sw(nswg), gpb_lw(nlwg)
    call drv_gpt_bands(nswg, nlwg, gpb_sw, gpb_lw)
  end subroutine w_gpt_bands

  subroutine w_temp_limits(tmin, tmax)
    use rrtmgp_core, only: drv_temp_limits
    real(r8), intent(out) :: tmin, tmax
    call drv_temp_limits(tmin, tmax)
  end subroutine w_temp_limits

  subroutine w_gas_optics_sw(ncol, nlev, ngas, ngpt, gas_vmr, pmid, &
       tmid, pint, tau, ssa, g, toa_flux)
    use rrtmgp_core, only: drv_gas_optics_sw
    integer, intent(in) :: ncol, nlev, ngas, ngpt
    real(r8), intent(in) :: gas_vmr(ngas, ncol, nlev)
    real(r8), intent(in) :: pmid(ncol, nlev), tmid(ncol, nlev)
    real(r8), intent(in) :: pint(ncol, nlev + 1)
    real(r8), intent(out) :: tau(ncol, nlev, ngpt), ssa(ncol, nlev, ngpt)
    real(r8), intent(out) :: g(ncol, nlev, ngpt), toa_flux(ncol, ngpt)
    call drv_gas_optics_sw(ncol, nlev, ngas, ngpt, gas_vmr, pmid, tmid, &
                           pint, tau, ssa, g, toa_flux)
  end subroutine w_gas_optics_sw

  subroutine w_gas_optics_lw(ncol, nlev, ngas, ngpt, gas_vmr, pmid, &
       tmid, pint, tint, tau, lay_src, lev_src_inc, lev_src_dec, sfc_src)
    use rrtmgp_core, only: drv_gas_optics_lw
    integer, intent(in) :: ncol, nlev, ngas, ngpt
    real(r8), intent(in) :: gas_vmr(ngas, ncol, nlev)
    real(r8), intent(in) :: pmid(ncol, nlev), tmid(ncol, nlev)
    real(r8), intent(in) :: pint(ncol, nlev + 1), tint(ncol, nlev + 1)
    real(r8), intent(out) :: tau(ncol, nlev, ngpt)
    real(r8), intent(out) :: lay_src(ncol, nlev, ngpt)
    real(r8), intent(out) :: lev_src_inc(ncol, nlev, ngpt)
    real(r8), intent(out) :: lev_src_dec(ncol, nlev, ngpt)
    real(r8), intent(out) :: sfc_src(ncol, ngpt)
    call drv_gas_optics_lw(ncol, nlev, ngas, ngpt, gas_vmr, pmid, tmid, &
         pint, tint, tau, lay_src, lev_src_inc, lev_src_dec, sfc_src)
  end subroutine w_gas_optics_lw

  subroutine w_cloud_optics_sw(ncol, nlev, nbnd, do_snow, cld, cldfsnow, &
       iclwp, iciwp, icswp, lambdac, mu, dei, des, rel, rei, &
       tau_out, ssa_out, asm_out, liq_tau_out, ice_tau_out, snw_tau_out)
    use rrtmgp_core, only: drv_cloud_optics_sw
    integer, intent(in) :: ncol, nlev, nbnd, do_snow
    real(r8), intent(in), dimension(ncol, nlev) :: cld, cldfsnow, &
         iclwp, iciwp, icswp, lambdac, mu, dei, des, rel, rei
    real(r8), intent(out), dimension(ncol, nlev, nbnd) :: tau_out, &
         ssa_out, asm_out, liq_tau_out, ice_tau_out, snw_tau_out
    call drv_cloud_optics_sw(ncol, nlev, nbnd, do_snow, cld, cldfsnow, &
         iclwp, iciwp, icswp, lambdac, mu, dei, des, rel, rei, &
         tau_out, ssa_out, asm_out, liq_tau_out, ice_tau_out, snw_tau_out)
  end subroutine w_cloud_optics_sw

  subroutine w_cloud_optics_lw(ncol, nlev, nbnd, do_snow, cld, cldfsnow, &
       iclwp, iciwp, icswp, lambdac, mu, dei, des, rei, &
       tau_out, liq_tau_out, ice_tau_out, snw_tau_out)
    use rrtmgp_core, only: drv_cloud_optics_lw
    integer, intent(in) :: ncol, nlev, nbnd, do_snow
    real(r8), intent(in), dimension(ncol, nlev) :: cld, cldfsnow, &
         iclwp, iciwp, icswp, lambdac, mu, dei, des, rei
    real(r8), intent(out), dimension(ncol, nlev, nbnd) :: tau_out, &
         liq_tau_out, ice_tau_out, snw_tau_out
    call drv_cloud_optics_lw(ncol, nlev, nbnd, do_snow, cld, cldfsnow, &
         iclwp, iciwp, icswp, lambdac, mu, dei, des, rei, &
         tau_out, liq_tau_out, ice_tau_out, snw_tau_out)
  end subroutine w_cloud_optics_lw

  subroutine w_mcica_mask(ngpt, ncol, nlev, changeseed, pmid, cldfrac, &
                          mask_out)
    use rrtmgp_core, only: drv_mcica_mask
    integer, intent(in) :: ngpt, ncol, nlev, changeseed
    real(r8), intent(in) :: pmid(ncol, nlev), cldfrac(ncol, nlev)
    integer, intent(out) :: mask_out(ngpt, ncol, nlev)
    call drv_mcica_mask(ngpt, ncol, nlev, changeseed, pmid, cldfrac, &
                        mask_out)
  end subroutine w_mcica_mask

  subroutine w_sample_sw(ncol, nlev, ngpt, nbnd, gpt2bnd, pmid, cld, &
       cldfsnow, tau_bnd, ssa_bnd, asm_bnd, tau_gpt, ssa_gpt, asm_gpt)
    use rrtmgp_core, only: drv_sample_sw
    integer, intent(in) :: ncol, nlev, ngpt, nbnd
    integer, intent(in) :: gpt2bnd(ngpt)
    real(r8), intent(in), dimension(ncol, nlev) :: pmid, cld, cldfsnow
    real(r8), intent(in), dimension(ncol, nlev, nbnd) :: tau_bnd, &
         ssa_bnd, asm_bnd
    real(r8), intent(out), dimension(ncol, nlev, ngpt) :: tau_gpt, &
         ssa_gpt, asm_gpt
    call drv_sample_sw(ncol, nlev, ngpt, nbnd, gpt2bnd, pmid, cld, &
         cldfsnow, tau_bnd, ssa_bnd, asm_bnd, tau_gpt, ssa_gpt, asm_gpt)
  end subroutine w_sample_sw

  subroutine w_sample_lw(ncol, nlev, ngpt, nbnd, gpt2bnd, pmid, cld, &
       cldfsnow, tau_bnd, tau_gpt)
    use rrtmgp_core, only: drv_sample_lw
    integer, intent(in) :: ncol, nlev, ngpt, nbnd
    integer, intent(in) :: gpt2bnd(ngpt)
    real(r8), intent(in), dimension(ncol, nlev) :: pmid, cld, cldfsnow
    real(r8), intent(in), dimension(ncol, nlev, nbnd) :: tau_bnd
    real(r8), intent(out), dimension(ncol, nlev, ngpt) :: tau_gpt
    call drv_sample_lw(ncol, nlev, ngpt, nbnd, gpt2bnd, pmid, cld, &
                       cldfsnow, tau_bnd, tau_gpt)
  end subroutine w_sample_lw

  subroutine w_run_sw(ngas, nday, nlevrad, ngpt, nbnd, gas_vmr, pmid, &
       tmid, pint, coszrs, alb_dir, alb_dif, cld_tau, cld_ssa, cld_asm, &
       aer_tau, aer_ssa, aer_asm, tsi_scaling, flx_all, flx_clr, &
       bnd_all, bnd_clr)
    use rrtmgp_core, only: drv_run_sw
    integer, intent(in) :: ngas, nday, nlevrad, ngpt, nbnd
    real(r8), intent(in) :: gas_vmr(ngas, nday, nlevrad)
    real(r8), intent(in) :: pmid(nday, nlevrad), tmid(nday, nlevrad)
    real(r8), intent(in) :: pint(nday, nlevrad + 1)
    real(r8), intent(in) :: coszrs(nday)
    real(r8), intent(in) :: alb_dir(nbnd, nday), alb_dif(nbnd, nday)
    real(r8), intent(in), dimension(nday, nlevrad, ngpt) :: cld_tau, &
         cld_ssa, cld_asm
    real(r8), intent(in), dimension(nday, nlevrad, nbnd) :: aer_tau, &
         aer_ssa, aer_asm
    real(r8), intent(in) :: tsi_scaling
    real(r8), intent(out) :: flx_all(nday, nlevrad + 1, 4)
    real(r8), intent(out) :: flx_clr(nday, nlevrad + 1, 4)
    real(r8), intent(out) :: bnd_all(nday, nlevrad + 1, nbnd, 4)
    real(r8), intent(out) :: bnd_clr(nday, nlevrad + 1, nbnd, 4)
    call drv_run_sw(ngas, nday, nlevrad, ngpt, nbnd, gas_vmr, pmid, &
         tmid, pint, coszrs, alb_dir, alb_dif, cld_tau, cld_ssa, cld_asm, &
         aer_tau, aer_ssa, aer_asm, tsi_scaling, flx_all, flx_clr, &
         bnd_all, bnd_clr)
  end subroutine w_run_sw

  subroutine w_run_lw(ngas, ncol, nlevrad, ngpt, nbnd, gas_vmr, pmid, &
       tmid, pint, tint, sfc_emis, cld_tau, aer_tau, flx_all, flx_clr, &
       bnd_all, bnd_clr)
    use rrtmgp_core, only: drv_run_lw
    integer, intent(in) :: ngas, ncol, nlevrad, ngpt, nbnd
    real(r8), intent(in) :: gas_vmr(ngas, ncol, nlevrad)
    real(r8), intent(in) :: pmid(ncol, nlevrad), tmid(ncol, nlevrad)
    real(r8), intent(in) :: pint(ncol, nlevrad + 1), tint(ncol, nlevrad + 1)
    real(r8), intent(in) :: sfc_emis(nbnd, ncol)
    real(r8), intent(in) :: cld_tau(ncol, nlevrad, ngpt)
    real(r8), intent(in) :: aer_tau(ncol, nlevrad, nbnd)
    real(r8), intent(out) :: flx_all(ncol, nlevrad + 1, 3)
    real(r8), intent(out) :: flx_clr(ncol, nlevrad + 1, 3)
    real(r8), intent(out) :: bnd_all(ncol, nlevrad + 1, nbnd, 3)
    real(r8), intent(out) :: bnd_clr(ncol, nlevrad + 1, nbnd, 3)
    call drv_run_lw(ngas, ncol, nlevrad, ngpt, nbnd, gas_vmr, pmid, &
         tmid, pint, tint, sfc_emis, cld_tau, aer_tau, flx_all, flx_clr, &
         bnd_all, bnd_clr)
  end subroutine w_run_lw

  subroutine w_rad_step(ncol, nlev, ngas, nswb, nlwb, nswg, nlwg, &
       do_snow, t_in, pmid_in, pint_in, lnpmid_in, lnpint_in, &
       lwup, asdir, asdif, aldir, aldif, coszrs, &
       cld, cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des, rel, rei, &
       gas_vmr, aer_tau_sw, aer_ssa_sw, aer_asm_sw, aer_tau_lw, &
       tsi_scaling, &
       qrs, qrsc, qrl, qrlc, sw_all, sw_clr, lw_all, lw_clr, srf, &
       diag_tmid, diag_tint, diag_alb, diag_cld_tau_bnd_sw, &
       diag_cld_gpt_sw, diag_cld_tau_bnd_lw, diag_cld_gpt_lw)
    use rrtmgp_core, only: drv_rad_step
    integer, intent(in) :: ncol, nlev, ngas, nswb, nlwb, nswg, nlwg
    integer, intent(in) :: do_snow
    real(r8), intent(in), dimension(ncol, nlev) :: t_in, pmid_in, &
         lnpmid_in, cld, cldfsnow, iclwp, iciwp, icswp, lambdac, mu, &
         dei, des, rel, rei
    real(r8), intent(in), dimension(ncol, nlev + 1) :: pint_in, lnpint_in
    real(r8), intent(in), dimension(ncol) :: lwup, asdir, asdif, aldir, &
         aldif, coszrs
    real(r8), intent(in) :: gas_vmr(ngas, ncol, nlev)
    real(r8), intent(in), dimension(ncol, nlev, nswb) :: aer_tau_sw, &
         aer_ssa_sw, aer_asm_sw
    real(r8), intent(in), dimension(ncol, nlev, nlwb) :: aer_tau_lw
    real(r8), intent(in) :: tsi_scaling
    real(r8), intent(out), dimension(ncol, nlev) :: qrs, qrsc, qrl, qrlc
    real(r8), intent(out) :: sw_all(ncol, nlev + 2, 4)
    real(r8), intent(out) :: sw_clr(ncol, nlev + 2, 4)
    real(r8), intent(out) :: lw_all(ncol, nlev + 2, 3)
    real(r8), intent(out) :: lw_clr(ncol, nlev + 2, 3)
    real(r8), intent(out) :: srf(ncol, 11)
    real(r8), intent(out) :: diag_tmid(ncol, nlev + 1)
    real(r8), intent(out) :: diag_tint(ncol, nlev + 2)
    real(r8), intent(out) :: diag_alb(nswb, ncol, 2)
    real(r8), intent(out) :: diag_cld_tau_bnd_sw(ncol, nlev, nswb)
    real(r8), intent(out) :: diag_cld_gpt_sw(ncol, nlev, nswg, 3)
    real(r8), intent(out) :: diag_cld_tau_bnd_lw(ncol, nlev, nlwb)
    real(r8), intent(out) :: diag_cld_gpt_lw(ncol, nlev, nlwg)
    call drv_rad_step(ncol, nlev, ngas, nswb, nlwb, nswg, nlwg, &
         do_snow, t_in, pmid_in, pint_in, lnpmid_in, lnpint_in, &
         lwup, asdir, asdif, aldir, aldif, coszrs, &
         cld, cldfsnow, iclwp, iciwp, icswp, lambdac, mu, dei, des, rel, rei, &
         gas_vmr, aer_tau_sw, aer_ssa_sw, aer_asm_sw, aer_tau_lw, &
         tsi_scaling, &
         qrs, qrsc, qrl, qrlc, sw_all, sw_clr, lw_all, lw_clr, srf, &
         diag_tmid, diag_tint, diag_alb, diag_cld_tau_bnd_sw, &
         diag_cld_gpt_sw, diag_cld_tau_bnd_lw, diag_cld_gpt_lw)
  end subroutine w_rad_step

end module rrtmgp_driver
