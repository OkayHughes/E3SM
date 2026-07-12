! f2py driver for the ZM convective transport routines
! (eam/src/physics/cam/zm/zm_transport.F90): zm_transport_tracer
! (tracer transport on the ZM mu/md/du/eu/ed mass fluxes) and
! zm_transport_momentum (momentum transport with the Gregory et al.
! pressure-gradient term and the Boville & Bretherton KE-dissipation
! heating).
!
! Gathered-array convention exactly as zm_conv_intr.F90 passes them:
! mu/md [mb/s], du/eu/ed [1/s], dp/dpdry [mb] are GATHERED rows
! i = 1..lengath; q/fracis/wind_in and the scattered outputs are full
! column arrays indexed through ideep. jt/mx/ideep are 1-based.
!
! zm_transport_tracer's only zm_param use is the zm_microp flag
! gating the negative-tracer conservation fixer; the driver sets it
! directly (both branches are driven -- the fixer is plain arithmetic
! and does not need zm_microphysics.F90). Constituent types ('wet' /
! 'dry') are injected through the constituents stub; dqdt slices the
! routine leaves untouched (m=1 water vapor and doconvtran=F) are
! zeroed by the driver, matching the zeroed ptend%q the intr layer
! hands in.
!
! zm_transport_momentum has NO runtime switch: the pressure-gradient
! constants momcu = momcd = 0.4 are compile-time parameters in
! zm_transport.F90 (no zmconv_mom* namelist), and nwind = 2 always
! (the KE fixer hard-codes components 1 and 2).
module zm_transport_driver
  implicit none
  ! local kind param (a use-associated kind is invisible to f2py's
  ! crackfortran; fbuild passes --f2cmap mapping r8 -> double)
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_transport_tracer(ncol, nlev, ncnst, zm_microp, &
       doconvtran, is_dry, q, mu, md, du, eu, ed, dp, jt, mx, ideep, &
       lengath, fracis, dpdry, dt, dqdt)
    use zm_transport, only: zm_transport_tracer, btype
    use zm_conv,      only: zm_param
    use constituents, only: cnst_stub_set_type

    integer, intent(in) :: ncol, nlev, ncnst
    integer, intent(in) :: zm_microp                 ! 0/1 flag
    integer, intent(in) :: doconvtran(ncnst)         ! 0/1 flags
    integer, intent(in) :: is_dry(ncnst)             ! 0/1 flags
    real(r8), intent(in) :: q(ncol, nlev, ncnst)     ! ungathered
    real(r8), intent(in) :: mu(ncol, nlev)           ! [mb/s] gathered
    real(r8), intent(in) :: md(ncol, nlev)           ! [mb/s] gathered
    real(r8), intent(in) :: du(ncol, nlev)           ! [1/s] gathered
    real(r8), intent(in) :: eu(ncol, nlev)           ! [1/s] gathered
    real(r8), intent(in) :: ed(ncol, nlev)           ! [1/s] gathered
    real(r8), intent(in) :: dp(ncol, nlev)           ! [mb] gathered
    integer, intent(in) :: jt(ncol)                  ! 1-based gathered
    integer, intent(in) :: mx(ncol)                  ! 1-based gathered
    integer, intent(in) :: ideep(ncol)               ! 1-based
    integer, intent(in) :: lengath
    real(r8), intent(in) :: fracis(ncol, nlev, ncnst)! ungathered
    real(r8), intent(in) :: dpdry(ncol, nlev)        ! [mb] gathered
    real(r8), intent(in) :: dt                       ! [s]
    real(r8), intent(out) :: dqdt(ncol, nlev, ncnst) ! ungathered

    logical(btype) :: doconv(ncnst)
    integer :: m

    zm_param%zm_microp = zm_microp /= 0
    do m = 1, ncnst
      if (is_dry(m) /= 0) then
        call cnst_stub_set_type(m, 'dry')
      else
        call cnst_stub_set_type(m, 'wet')
      end if
    end do
    doconv = doconvtran /= 0

    call zm_transport_tracer(ncol, nlev, doconv, q, ncnst, &
         mu, md, du, eu, ed, dp, jt, mx, ideep, 1, lengath, &
         fracis, dqdt, dpdry, dt)

    ! slices the routine never writes (m=1 and doconvtran=F) are
    ! undefined intent(out) memory: zero them deterministically, as
    ! the zeroed ptend%q of zm_conv_intr effectively does
    do m = 1, ncnst
      if (m == 1 .or. .not. doconv(m)) dqdt(:, :, m) = 0.0_r8
    end do
  end subroutine drv_transport_tracer

  subroutine drv_transport_momentum(ncol, nlev, lengath, wind_in, &
       mu, md, du, eu, ed, dp, jt, mx, ideep, dt, &
       wind_tend, pguall, pgdall, icwu, icwd, seten)
    use zm_transport, only: zm_transport_momentum

    integer, intent(in) :: ncol, nlev, lengath
    real(r8), intent(in) :: wind_in(ncol, nlev, 2)   ! [m/s] ungathered
    real(r8), intent(in) :: mu(ncol, nlev)           ! [mb/s] gathered
    real(r8), intent(in) :: md(ncol, nlev)           ! [mb/s] gathered
    real(r8), intent(in) :: du(ncol, nlev)           ! [1/s] gathered
    real(r8), intent(in) :: eu(ncol, nlev)           ! [1/s] gathered
    real(r8), intent(in) :: ed(ncol, nlev)           ! [1/s] gathered
    real(r8), intent(in) :: dp(ncol, nlev)           ! [mb] gathered
    integer, intent(in) :: jt(ncol)                  ! 1-based gathered
    integer, intent(in) :: mx(ncol)                  ! 1-based gathered
    integer, intent(in) :: ideep(ncol)               ! 1-based
    real(r8), intent(in) :: dt                       ! [s]
    real(r8), intent(out) :: wind_tend(ncol, nlev, 2)! [m/s2] ungathered
    real(r8), intent(out) :: pguall(ncol, nlev, 2)   ! ungathered
    real(r8), intent(out) :: pgdall(ncol, nlev, 2)   ! ungathered
    real(r8), intent(out) :: icwu(ncol, nlev, 2)     ! ungathered
    real(r8), intent(out) :: icwd(ncol, nlev, 2)     ! ungathered
    real(r8), intent(out) :: seten(ncol, nlev)       ! [J/kg/s] ungathered

    call zm_transport_momentum(ncol, ncol, nlev, nlev + 1, wind_in, 2, &
         mu, md, du, eu, ed, dp, jt, mx, ideep, 1, lengath, &
         wind_tend, pguall, pgdall, icwu, icwd, dt, seten)
  end subroutine drv_transport_momentum

end module zm_transport_driver
