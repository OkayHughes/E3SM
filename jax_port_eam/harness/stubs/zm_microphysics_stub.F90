! Abort-only interface stub for
! eam/src/physics/cam/zm/zm_microphysics.F90 (zm_mphy +
! zm_microphysics_adjust). zm_conv.F90 use-associates both at module
! scope, so the symbols are needed to COMPILE, but with the harness
! configuration zm_param%zm_microp = .false. neither is ever CALLED
! (the ZM golden scope is documented in PORTING_PLAN.md). The real
! module cannot be compiled standalone: it pulls in the aerosol
! activation stack (activate_drop_mam, ndrop_bam, nucleate_ice_conv)
! and CAM time_manager/pbuf infrastructure. Per METHODOLOGY.md, stubs
! never *supply* physics -- these only refuse to run (endrun).
!
! Signatures are copied argument-for-argument from the real
! zm_microphysics.F90 so the explicit module interfaces seen by
! zm_conv.F90 match its call sites; 2D arrays are assumed-size to
! avoid dragging in ppgrid.

module zm_microphysics
  use shr_kind_mod,          only: r8 => shr_kind_r8
  use zm_aero_type,          only: zm_aero_t
  use zm_conv_types,         only: zm_const_t
  use zm_microphysics_state, only: zm_microp_st
  implicit none
  private
  public :: zm_mphy, zm_microphysics_adjust

contains

  subroutine zm_mphy(pcols, il2g, msg, grav, cp, rd, auto_fac, accr_fac, dcs, jb, jt, jlcl, &
                     su, qu, mu, du, eu, zf, pm, te, qe, gamhat, eps0, cmel, cmei, aero, &
                     qc,    qi,   nc,   ni,   qcde,  qide,  ncde,  nide, rprd, sprd, frz,       &
                     wu,    qr,   qni,  nr,   ns,    qg,    ng,    qnide,nsde,                  &
                     autolm, accrlm, bergnm, fhtimm, fhtctm,    &
                     fhmlm,  hmpim,  accslm, dlfm,   autoln, accrln, bergnn, fhtimn, fhtctn,    &
                     fhmln,  accsln, activn, dlfn,   autoim, accsim, difm,   nuclin, autoin,    &
                     accsin, hmpin,  difn,   trspcm, trspcn, trspim, trspin, lamc,   pgam,      &
                     accgrm, accglm, accgslm,accgsrm,accgirm,accgrim,accgrsm,accgsln,accgsrn,   &
                     accgirn,accsrim,acciglm,accigrm,accsirm,accigln,accigrn,accsirn,accgln,    &
                     accgrn ,accilm, acciln ,fallrm ,fallsm ,fallgm ,fallrn ,fallsn ,fallgn,    &
                     fhmrm  ,dsfm, dsfn)
    use cam_abortutils, only: endrun
    integer,  intent(in) :: pcols, il2g, msg
    real(r8), intent(in) :: grav, cp, rd, auto_fac, accr_fac, dcs
    integer,  intent(in) :: jb(pcols), jt(pcols), jlcl(pcols)
    real(r8), intent(in) :: su(pcols,*), qu(pcols,*), mu(pcols,*)
    real(r8), intent(in) :: du(pcols,*), eu(pcols,*), zf(pcols,*)
    real(r8), intent(in) :: pm(pcols,*), te(pcols,*), qe(pcols,*)
    real(r8), intent(in) :: gamhat(pcols,*), eps0(pcols)
    real(r8), intent(in) :: cmel(pcols,*), cmei(pcols,*)
    type(zm_aero_t), intent(in) :: aero
    real(r8) :: qc(pcols,*), qi(pcols,*), nc(pcols,*), ni(pcols,*)
    real(r8) :: qcde(pcols,*), qide(pcols,*), ncde(pcols,*), nide(pcols,*)
    real(r8) :: rprd(pcols,*), sprd(pcols,*), frz(pcols,*)
    real(r8) :: wu(pcols,*), qr(pcols,*), qni(pcols,*), nr(pcols,*)
    real(r8) :: ns(pcols,*), qg(pcols,*), ng(pcols,*), qnide(pcols,*)
    real(r8) :: nsde(pcols,*)
    real(r8) :: autolm(pcols,*), accrlm(pcols,*), bergnm(pcols,*)
    real(r8) :: fhtimm(pcols,*), fhtctm(pcols,*), fhmlm(pcols,*)
    real(r8) :: hmpim(pcols,*), accslm(pcols,*), dlfm(pcols,*)
    real(r8) :: autoln(pcols,*), accrln(pcols,*), bergnn(pcols,*)
    real(r8) :: fhtimn(pcols,*), fhtctn(pcols,*), fhmln(pcols,*)
    real(r8) :: accsln(pcols,*), activn(pcols,*), dlfn(pcols,*)
    real(r8) :: autoim(pcols,*), accsim(pcols,*), difm(pcols,*)
    real(r8) :: nuclin(pcols,*), autoin(pcols,*), accsin(pcols,*)
    real(r8) :: hmpin(pcols,*), difn(pcols,*)
    real(r8) :: trspcm(pcols,*), trspcn(pcols,*), trspim(pcols,*)
    real(r8) :: trspin(pcols,*), lamc(pcols,*), pgam(pcols,*)
    real(r8) :: accgrm(pcols,*), accglm(pcols,*), accgslm(pcols,*)
    real(r8) :: accgsrm(pcols,*), accgirm(pcols,*), accgrim(pcols,*)
    real(r8) :: accgrsm(pcols,*), accgsln(pcols,*), accgsrn(pcols,*)
    real(r8) :: accgirn(pcols,*), accsrim(pcols,*), acciglm(pcols,*)
    real(r8) :: accigrm(pcols,*), accsirm(pcols,*), accigln(pcols,*)
    real(r8) :: accigrn(pcols,*), accsirn(pcols,*), accgln(pcols,*)
    real(r8) :: accgrn(pcols,*), accilm(pcols,*), acciln(pcols,*)
    real(r8) :: fallrm(pcols,*), fallsm(pcols,*), fallgm(pcols,*)
    real(r8) :: fallrn(pcols,*), fallsn(pcols,*), fallgn(pcols,*)
    real(r8) :: fhmrm(pcols,*), dsfm(pcols,*), dsfn(pcols,*)
    call endrun('zm_mphy stub called: harness only supports zm_microp=.false.')
  end subroutine zm_mphy

  subroutine zm_microphysics_adjust(pcols, ncol, pver, jt, msg, delt, zm_const, &
                                    dp, qh, dl, dsdt, dqdt, prdprec, microp_st)
    use cam_abortutils, only: endrun
    integer,  intent(in) :: pcols, ncol, pver, msg
    integer,  intent(in) :: jt(pcols)
    real(r8), intent(in) :: delt
    type(zm_const_t), intent(in) :: zm_const
    real(r8) :: dp(pcols,*), qh(pcols,*), dl(pcols,*)
    real(r8) :: dsdt(pcols,*), dqdt(pcols,*), prdprec(pcols,*)
    type(zm_microp_st) :: microp_st
    call endrun('zm_microphysics_adjust stub called: harness only supports zm_microp=.false.')
  end subroutine zm_microphysics_adjust

end module zm_microphysics
