! f2py driver for the cldfrc2m RH-based stratiform cloud-fraction
! module (eam/src/physics/cam/cldfrc2m.F90): astG_PDF[_single]
! (triangular-PDF liquid stratus fraction + inverse gradient G=dU/da),
! astG_RHU[_single] (CAM3.5 quadratic formula) and aist_single /
! aist_vector (ice stratus fraction, iceopt 1-7).
!
! Runtime parameters are injected through drv_init:
!  - rhminl/rhminl_adj_land/rhminh/premit/premib/iceopt/icecrit/minice
!    go into the cloud_fraction STUB (cloud_fraction_stub.F90), and the
!    REAL cldfrc2m_init() pulls them into its private module variables
!    via the real cldfrc_getparams keyword interface.
!  - rhmini/rhmaxi CANNOT be set this way: cldfrc2m_readnl assigns the
!    protected rhmini_const/rhmaxi_const only under masterproc (stubbed
!    .false.) and the SPMD broadcast is compiled out. The aist drivers
!    therefore ALWAYS pass rhmaxi_in/rhmini_in explicitly (every use of
!    rhmini/rhmaxi in aist_* is overridden when present).
!
! EAMv3 phys="default" values (namelist_defaults_eam.xml, recorded in
! the golden metadata):
!   cldfrc_rhminl          = 0.950  (phys="default" microphys="p3")
!   cldfrc_rhminl_adj_land = 0.100
!   cldfrc_rhminh          = 0.800
!   cldfrc_premit          = 25000.0 Pa (dyn="se")
!   cldfrc_premib          = 70000.0 Pa (phys="default")
!   cldfrc_iceopt          = 5      (phys="default")
!   cldfrc_icecrit         = 0.93   (phys="default")
!   cldfrc_minice          = 1.0e-12
!   cldfrc2m_rhmini        = 0.80
!   cldfrc2m_rhmaxi        = 1.05   (phys="default" clubb_sgs="1"
!                                    clubb_do_deep="0")
!
! The vector routines take pcols-sized (=16, grid_stubs.F90) assumed
! arrays with an ncol count, so the drivers chunk arbitrary-length
! inputs through pcols-wide buffers (columns are independent; chunking
! is exact). wv_sat_init builds the SVP table needed by
! qsat_water/svp_water/svp_ice inside aist_*.
module cldfrc2m_driver
  implicit none
  ! local kind param (a use-associated kind is invisible to f2py's
  ! crackfortran; fbuild passes --f2cmap mapping r8 -> double)
  integer, parameter :: r8 = selected_real_kind(12)
  logical :: wv_sat_ready = .false.

contains

  subroutine drv_init(rhminl, rhminl_adj_land, rhminh, premit, premib, &
                      iceopt, icecrit, minice)
    use cloud_fraction, only: cldfrc_stub_set
    use cldfrc2m,       only: cldfrc2m_init
    use wv_saturation,  only: wv_sat_init
    real(r8), intent(in) :: rhminl, rhminl_adj_land, rhminh
    real(r8), intent(in) :: premit, premib
    integer,  intent(in) :: iceopt
    real(r8), intent(in) :: icecrit, minice
    call cldfrc_stub_set(rhminl, rhminl_adj_land, rhminh, premit, &
                         premib, iceopt, icecrit, minice)
    call cldfrc2m_init()
    ! GoffGratch default scheme; the stubbed namelist read leaves the
    ! wv_saturation defaults untouched. Once only: wv_sat_init
    ! allocates its SVP table and aborts if called twice, while
    ! drv_init is re-invoked per iceopt by the golden generator.
    if (.not. wv_sat_ready) then
       call wv_sat_init()
       wv_sat_ready = .true.
    end if
  end subroutine drv_init

  subroutine drv_astg_pdf_single(u, p, qv, landfrac, snowh, a, ga, rhmin)
    use cldfrc2m, only: astG_PDF_single
    real(r8), intent(in)  :: u, p, qv, landfrac, snowh
    real(r8), intent(out) :: a, ga, rhmin
    call astG_PDF_single(u, p, qv, landfrac, snowh, a, ga, orhmin=rhmin)
  end subroutine drv_astg_pdf_single

  subroutine drv_astg_rhu_single(u, p, qv, landfrac, snowh, a, ga, rhmin)
    use cldfrc2m, only: astG_RHU_single
    real(r8), intent(in)  :: u, p, qv, landfrac, snowh
    real(r8), intent(out) :: a, ga, rhmin
    call astG_RHU_single(u, p, qv, landfrac, snowh, a, ga, orhmin=rhmin)
  end subroutine drv_astg_rhu_single

  subroutine drv_astg_pdf(n, u, p, qv, landfrac, snowh, a, ga)
    use ppgrid,   only: pcols
    use cldfrc2m, only: astG_PDF
    integer,  intent(in)  :: n
    real(r8), intent(in)  :: u(n), p(n), qv(n), landfrac(n), snowh(n)
    real(r8), intent(out) :: a(n), ga(n)
    real(r8), dimension(pcols) :: ub, pb, qvb, lfb, shb, ab, gab
    integer :: beg, m
    do beg = 1, n, pcols
       m = min(pcols, n - beg + 1)
       ub = 0._r8; pb = 1.e5_r8; qvb = 0._r8; lfb = 0._r8; shb = 0._r8
       ub(1:m)  = u(beg:beg+m-1)
       pb(1:m)  = p(beg:beg+m-1)
       qvb(1:m) = qv(beg:beg+m-1)
       lfb(1:m) = landfrac(beg:beg+m-1)
       shb(1:m) = snowh(beg:beg+m-1)
       call astG_PDF(ub, pb, qvb, lfb, shb, ab, gab, m)
       a(beg:beg+m-1)  = ab(1:m)
       ga(beg:beg+m-1) = gab(1:m)
    end do
  end subroutine drv_astg_pdf

  subroutine drv_astg_rhu(n, u, p, qv, landfrac, snowh, a, ga)
    use ppgrid,   only: pcols
    use cldfrc2m, only: astG_RHU
    integer,  intent(in)  :: n
    real(r8), intent(in)  :: u(n), p(n), qv(n), landfrac(n), snowh(n)
    real(r8), intent(out) :: a(n), ga(n)
    real(r8), dimension(pcols) :: ub, pb, qvb, lfb, shb, ab, gab
    integer :: beg, m
    do beg = 1, n, pcols
       m = min(pcols, n - beg + 1)
       ub = 0._r8; pb = 1.e5_r8; qvb = 0._r8; lfb = 0._r8; shb = 0._r8
       ub(1:m)  = u(beg:beg+m-1)
       pb(1:m)  = p(beg:beg+m-1)
       qvb(1:m) = qv(beg:beg+m-1)
       lfb(1:m) = landfrac(beg:beg+m-1)
       shb(1:m) = snowh(beg:beg+m-1)
       call astG_RHU(ub, pb, qvb, lfb, shb, ab, gab, m)
       a(beg:beg+m-1)  = ab(1:m)
       ga(beg:beg+m-1) = gab(1:m)
    end do
  end subroutine drv_astg_rhu

  subroutine drv_aist_single(qv, t, p, qi, landfrac, snowh, &
                             rhmaxi, rhmini, aist)
    use cldfrc2m, only: aist_single
    real(r8), intent(in)  :: qv, t, p, qi, landfrac, snowh
    real(r8), intent(in)  :: rhmaxi, rhmini
    real(r8), intent(out) :: aist
    call aist_single(qv, t, p, qi, landfrac, snowh, aist, &
                     rhmaxi_in=rhmaxi, rhmini_in=rhmini)
  end subroutine drv_aist_single

  subroutine drv_aist_vector(n, qv, t, p, qi, ni, landfrac, snowh, &
                             rhmaxi, rhmini, aist)
    use ppgrid,   only: pcols
    use cldfrc2m, only: aist_vector
    integer,  intent(in)  :: n
    real(r8), intent(in)  :: qv(n), t(n), p(n), qi(n), ni(n)
    real(r8), intent(in)  :: landfrac(n), snowh(n)
    real(r8), intent(in)  :: rhmaxi, rhmini
    real(r8), intent(out) :: aist(n)
    real(r8), dimension(pcols) :: qvb, tb, pb, qib, nib, lfb, shb, ab
    real(r8), dimension(pcols) :: rmib
    integer :: beg, m
    rmib = rhmini
    do beg = 1, n, pcols
       m = min(pcols, n - beg + 1)
       ! benign pads (aist_vector only touches 1:ncol)
       qvb = 1.e-6_r8; tb = 250._r8; pb = 5.e4_r8
       qib = 1.e-9_r8; nib = 1.e3_r8; lfb = 0._r8; shb = 0._r8
       qvb(1:m) = qv(beg:beg+m-1)
       tb(1:m)  = t(beg:beg+m-1)
       pb(1:m)  = p(beg:beg+m-1)
       qib(1:m) = qi(beg:beg+m-1)
       nib(1:m) = ni(beg:beg+m-1)
       lfb(1:m) = landfrac(beg:beg+m-1)
       shb(1:m) = snowh(beg:beg+m-1)
       call aist_vector(qvb, tb, pb, qib, nib, lfb, shb, ab, m, &
                        rhmaxi_in=rhmaxi, rhmini_in=rmib)
       aist(beg:beg+m-1) = ab(1:m)
    end do
  end subroutine drv_aist_vector

end module cldfrc2m_driver
