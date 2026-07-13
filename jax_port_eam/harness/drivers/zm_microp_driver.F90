! f2py driver for the ZM convective microphysics
! (eam/src/physics/cam/zm/zm_microphysics.F90, compiled UNMODIFIED)
! and for zm_conv.F90 zm_conv_main with zm_param%zm_microp = .true.
! (the EAMv3 default path, zmconv_microp=.true.).
!
! AEROSOL ACTIVATION STRATEGY (PORTING_PLAN.md row 9): the REAL
! activation modules are compiled unmodified --
!   cam/activate_drop_mam.F90  (actdrop_mam_calc, AR-G 2000 droplet
!                               activation; init via the real
!                               actdrop_mam_init exactly as
!                               ndrop.F90 calls it)
!   cam/nucleate_ice_conv.F90  (nucleati_conv, Liu & Penner 2005)
! Only the BULK activation (ndrop_bam) is an abort-only stub
! (stubs/zm_microp_stubs.F90): EAMv3 runs MAM modal aerosols, so
! aero%scheme is always 'modal' here and ndrop_bam_run is never
! executed. shr_spfn_mod is the real share/util module; under
! gfortran its erf/gamma resolve to the compiler intrinsics.
!
! The zm_aero_t modal object is filled from plain driver arrays the
! same way zm_aero.F90 zm_aero_init + zm_conv_intr fill it from
! rad_constituents + pbuf: mode indices, per-mode species counts,
! specdens/spechygro, voltonumblo/hi, sigmag_aitken, and per-column
! num/mmr/dgnum profiles. All values are generator inputs recorded in
! the golden metadata (rad_constituents only supplies DATA, not
! physics, so prescribing it keeps the physics untouched).
!
! Compiled WITHOUT MODAL_AERO_4MODE_MOM / MODAL_AERO_5MODE /
! RAIN_EVAP_TO_COARSE_AERO, so the coarse-mode dust weight in
! zm_mphy is wght = dmc/(ssmc+dmc+so4mc) (the 3-species variant); the
! JAX port implements the same variant.
!
! deltat inside zm_mphy comes from the time_manager stub
! (get_step_size), set via drv_init.
!
! Output packing: zm_mphy's 22 state/inout arrays and 60 diagnostic
! tendency arrays are returned stacked along the last dimension in
! the documented order below (STATE_FIELDS / DIAG_FIELDS in
! gen_zm_microp_golden.py must match):
!
! state (22): qc, qi, nc, ni, qcde, qide, qnide, ncde, nide, nsde,
!             qni, qr, ns, nr, qg, ng, rprd, sprd, frz, wu, lamc, pgam
! diag (60):  autolm, accrlm, bergnm, fhtimm, fhtctm, fhmlm, hmpim,
!             accslm, dlfm, autoln, accrln, bergnn, fhtimn, fhtctn,
!             fhmln, accsln, activn, dlfn, autoim, accsim, difm,
!             nuclin, autoin, accsin, hmpin, difn, trspcm, trspcn,
!             trspim, trspin, accgrm, accglm, accgslm, accgsrm,
!             accgirm, accgrim, accgrsm, accgsln, accgsrn, accgirn,
!             accsrim, acciglm, accigrm, accsirm, accigln, accigrn,
!             accsirn, accgln, accgrn, accilm, acciln, fallrm,
!             fallsm, fallgm, fallrn, fallsn, fallgn, fhmrm, dsfm,
!             dsfn
!
! zm_conv_main microp_st packing (27 2D fields, column-scattered
! space, + rice separately):
!             wu, qliq, qice, qrain, qsnow, qgraupel, qnl, qni, qnr,
!             qns, qng, sprd, mudpcu, lambdadpcu, qcde, qide, qsde,
!             ncde, nide, nsde, dif, dsf, dnlf, dnif, dnsf, frz, cmel,
!             cmei  -- actually 28; see MICROP_FIELDS below.
module zm_microp_driver
  implicit none
  ! local kind param (a use-associated kind is invisible to f2py's
  ! crackfortran; fbuild passes --f2cmap mapping r8 -> double)
  integer, parameter :: r8 = selected_real_kind(12)

contains

  subroutine drv_init(step_size, nmodes, sigmag_amode)
    ! wv_saturation SVP table (GoffGratch defaults), the real
    ! actdrop_mam_init exactly as ndrop.F90 aerosol_activate_init
    ! calls it, zm_mphyi (derived microphysics constants), and the
    ! time_manager stub step size.
    use wv_saturation,     only: wv_sat_init
    use activate_drop_mam, only: actdrop_mam_init
    use zm_microphysics,   only: zm_mphyi
    use time_manager,      only: set_step_size
    use physconst,         only: pi, rair, r_universal, rh2o, rhoh2o, &
                                 mwh2o, epsilo, latvap, cpair, gravit
    integer, intent(in) :: step_size          ! [s] model time step
    integer, intent(in) :: nmodes             ! number of aerosol modes
    real(r8), intent(in) :: sigmag_amode(nmodes) ! geometric stddev per mode

    call wv_sat_init()
    call set_step_size(step_size)
    call actdrop_mam_init(6, pi, rair, r_universal, rh2o, &
                          rhoh2o, mwh2o, epsilo, latvap, cpair, &
                          gravit, nmodes, sigmag_amode)
    call zm_mphyi()
  end subroutine drv_init

  subroutine drv_actdrop(np, nmodes, wbar, sigw, wdiab, wminf, wmaxf, &
       tair, rhoair, na, volume, hygro, in_cloud, smax_f, &
       fn, fm, fluxn, fluxm, flux_fullact)
    ! direct actdrop_mam_calc call for unit bisection (drv_init must
    ! have been called first)
    use activate_drop_mam, only: actdrop_mam_calc
    integer, intent(in) :: np, nmodes
    real(r8), intent(in) :: wbar(np), sigw(np), wdiab(np)
    real(r8), intent(in) :: wminf(np), wmaxf(np)
    real(r8), intent(in) :: tair(np), rhoair(np)
    real(r8), intent(in) :: na(np, nmodes), volume(np, nmodes)
    real(r8), intent(in) :: hygro(np, nmodes)
    integer, intent(in) :: in_cloud(np)         ! 0/1
    real(r8), intent(in) :: smax_f(np)
    real(r8), intent(out) :: fn(np, nmodes), fm(np, nmodes)
    real(r8), intent(out) :: fluxn(np, nmodes), fluxm(np, nmodes)
    real(r8), intent(out) :: flux_fullact(np)
    integer :: i
    do i = 1, np
      call actdrop_mam_calc(wbar(i), sigw(i), wdiab(i), wminf(i), &
           wmaxf(i), tair(i), rhoair(i), na(i, :), nmodes, &
           volume(i, :), hygro(i, :), in_cloud(i) /= 0, smax_f(i), &
           fn(i, :), fm(i, :), fluxn(i, :), fluxm(i, :), &
           flux_fullact(i))
    end do
  end subroutine drv_actdrop

  subroutine drv_nucleati(np, wbar, tair, relhum, cldn, qc, nfice, &
       rhoair, so4_num, dst_num, soot_num, nuci, onihf, oniimm, &
       onidep, onimey)
    ! direct nucleati_conv call (zm_microp=.true.) for unit bisection
    use nucleate_ice_conv, only: nucleati_conv
    integer, intent(in) :: np
    real(r8), intent(in) :: wbar(np), tair(np), relhum(np), cldn(np)
    real(r8), intent(in) :: qc(np), nfice(np), rhoair(np)
    real(r8), intent(in) :: so4_num(np), dst_num(np), soot_num(np)
    real(r8), intent(out) :: nuci(np), onihf(np), oniimm(np)
    real(r8), intent(out) :: onidep(np), onimey(np)
    integer :: i
    do i = 1, np
      call nucleati_conv(wbar(i), tair(i), relhum(i), cldn(i), qc(i), &
           nfice(i), rhoair(i), so4_num(i), dst_num(i), soot_num(i), &
           .true., nuci(i), onihf(i), oniimm(i), onidep(i), onimey(i))
    end do
  end subroutine drv_nucleati

  subroutine drv_zm_mphy(ncol, nlev, nmodes, nspecmx, msg, &
       auto_fac, accr_fac, dcs, lambdadpcu0, mudpcu0, &
       jb, jt, jlcl, su, qu, mu, du, eu, zf, pm, te, qe, gamhat, &
       eps0, cmel, cmei, &
       nspec, m_accum, m_aitken, m_coarse, l_dust, l_nacl, l_so4, &
       sigmag_aitken, specdens, spechygro, voltonumblo, voltonumbhi, &
       numg_a, mmrg_a, dgnumg, &
       state_out, diag_out)
    ! Direct zm_mphy kernel call on gathered plume profiles.
    ! grav/cp/rd are the zm_const values (zm_const_set_to_global).
    use ppgrid,          only: ppver => pver
    use zm_conv,         only: zm_conv_main_init, zm_const
    use zm_microphysics, only: zm_mphy
    use zm_aero_type,    only: zm_aero_t
    use cam_abortutils,  only: endrun

    integer, intent(in) :: ncol, nlev, nmodes, nspecmx
    integer, intent(in) :: msg                     ! levels skipped at top
    real(r8), intent(in) :: auto_fac, accr_fac     ! zmconv_auto_fac/accr_fac
    real(r8), intent(in) :: dcs                    ! zmconv_micro_dcs [m]
    real(r8), intent(in) :: lambdadpcu0, mudpcu0   ! lamc/pgam seed values
    integer, intent(in) :: jb(ncol), jt(ncol), jlcl(ncol)  ! 1-based
    real(r8), intent(in) :: su(ncol, nlev)         ! updraft DSE/cp [K]
    real(r8), intent(in) :: qu(ncol, nlev)         ! updraft q [kg/kg]
    real(r8), intent(in) :: mu(ncol, nlev)         ! updraft mass flux [normalized]
    real(r8), intent(in) :: du(ncol, nlev)         ! detrainment [1/m]
    real(r8), intent(in) :: eu(ncol, nlev)         ! entrainment [1/m]
    real(r8), intent(in) :: zf(ncol, nlev + 1)     ! interface height [m]
    real(r8), intent(in) :: pm(ncol, nlev)         ! env pressure [hPa]
    real(r8), intent(in) :: te(ncol, nlev)         ! env T [K]
    real(r8), intent(in) :: qe(ncol, nlev)         ! env q [kg/kg]
    real(r8), intent(in) :: gamhat(ncol, nlev)     ! L/cp dq*/dT at interface
    real(r8), intent(in) :: eps0(ncol)             ! lambda_max
    real(r8), intent(in) :: cmel(ncol, nlev)       ! liq condensation rate
    real(r8), intent(in) :: cmei(ncol, nlev)       ! ice deposition rate
    integer, intent(in) :: nspec(nmodes)
    integer, intent(in) :: m_accum, m_aitken, m_coarse  ! 1-based mode idx
    integer, intent(in) :: l_dust, l_nacl, l_so4        ! coarse species idx
    real(r8), intent(in) :: sigmag_aitken
    real(r8), intent(in) :: specdens(nspecmx, nmodes)   ! [kg/m3]
    real(r8), intent(in) :: spechygro(nspecmx, nmodes)
    real(r8), intent(in) :: voltonumblo(nmodes), voltonumbhi(nmodes)
    real(r8), intent(in) :: numg_a(ncol, nlev, nmodes)  ! [#/kg]
    real(r8), intent(in) :: mmrg_a(ncol, nlev, nspecmx, nmodes) ! [kg/kg]
    real(r8), intent(in) :: dgnumg(ncol, nlev, nmodes)  ! [m]

    real(r8), intent(out) :: state_out(ncol, nlev, 22)
    real(r8), intent(out) :: diag_out(ncol, nlev, 60)

    type(zm_aero_t) :: aero

    if (nlev /= ppver) call endrun('drv_zm_mphy: nlev /= ppgrid pver')

    ! fills zm_const via the real zm_const_set_to_global (limcnv value
    ! irrelevant for zm_mphy; masterproc=.false. mutes the print)
    call zm_conv_main_init(1)

    aero%scheme = 'modal'
    aero%nmodes = nmodes
    allocate(aero%nspec(nmodes))
    aero%nspec = nspec
    aero%mode_accum_idx  = m_accum
    aero%mode_aitken_idx = m_aitken
    aero%mode_coarse_idx = m_coarse
    aero%coarse_dust_idx = l_dust
    aero%coarse_nacl_idx = l_nacl
    aero%coarse_so4_idx  = l_so4
    aero%sigmag_aitken   = sigmag_aitken
    allocate(aero%specdens(nspecmx, nmodes), aero%spechygro(nspecmx, nmodes))
    aero%specdens  = specdens
    aero%spechygro = spechygro
    allocate(aero%voltonumblo(nmodes), aero%voltonumbhi(nmodes))
    aero%voltonumblo = voltonumblo
    aero%voltonumbhi = voltonumbhi
    allocate(aero%numg_a(ncol, nlev, nmodes))
    allocate(aero%mmrg_a(ncol, nlev, nspecmx, nmodes))
    allocate(aero%dgnumg(ncol, nlev, nmodes))
    aero%numg_a = numg_a
    aero%mmrg_a = mmrg_a
    aero%dgnumg = dgnumg

    state_out = 0._r8
    diag_out  = 0._r8
    ! lamc/pgam are intent(inout): seeded exactly as zm_conv_main
    ! seeds loc_microp_st%lambdadpcu/mudpcu
    state_out(:, :, 21) = lambdadpcu0
    state_out(:, :, 22) = mudpcu0

    call zm_mphy(ncol, ncol, msg, &
         zm_const%grav, zm_const%cpair, zm_const%rdair, &
         auto_fac, accr_fac, dcs, jb, jt, jlcl, &
         su, qu, mu, du, eu, zf, pm, te, qe, gamhat, eps0, &
         cmel, cmei, aero, &
         state_out(:, :, 1),  state_out(:, :, 2),  state_out(:, :, 3),  state_out(:, :, 4),  &
         state_out(:, :, 5),  state_out(:, :, 6),  state_out(:, :, 8),  state_out(:, :, 9),  &
         state_out(:, :, 17), state_out(:, :, 18), state_out(:, :, 19), state_out(:, :, 20), &
         state_out(:, :, 12), state_out(:, :, 11), state_out(:, :, 14), state_out(:, :, 13), &
         state_out(:, :, 15), state_out(:, :, 16), state_out(:, :, 7),  state_out(:, :, 10), &
         diag_out(:, :, 1),  diag_out(:, :, 2),  diag_out(:, :, 3),  diag_out(:, :, 4),  &
         diag_out(:, :, 5),  diag_out(:, :, 6),  diag_out(:, :, 7),  diag_out(:, :, 8),  &
         diag_out(:, :, 9),  diag_out(:, :, 10), diag_out(:, :, 11), diag_out(:, :, 12), &
         diag_out(:, :, 13), diag_out(:, :, 14), diag_out(:, :, 15), diag_out(:, :, 16), &
         diag_out(:, :, 17), diag_out(:, :, 18), diag_out(:, :, 19), diag_out(:, :, 20), &
         diag_out(:, :, 21), diag_out(:, :, 22), diag_out(:, :, 23), diag_out(:, :, 24), &
         diag_out(:, :, 25), diag_out(:, :, 26), diag_out(:, :, 27), diag_out(:, :, 28), &
         diag_out(:, :, 29), diag_out(:, :, 30), &
         state_out(:, :, 21), state_out(:, :, 22), &
         diag_out(:, :, 31), diag_out(:, :, 32), diag_out(:, :, 33), diag_out(:, :, 34), &
         diag_out(:, :, 35), diag_out(:, :, 36), diag_out(:, :, 37), diag_out(:, :, 38), &
         diag_out(:, :, 39), diag_out(:, :, 40), diag_out(:, :, 41), diag_out(:, :, 42), &
         diag_out(:, :, 43), diag_out(:, :, 44), diag_out(:, :, 45), diag_out(:, :, 46), &
         diag_out(:, :, 47), diag_out(:, :, 48), diag_out(:, :, 49), diag_out(:, :, 50), &
         diag_out(:, :, 51), diag_out(:, :, 52), diag_out(:, :, 53), diag_out(:, :, 54), &
         diag_out(:, :, 55), diag_out(:, :, 56), diag_out(:, :, 57), diag_out(:, :, 58), &
         diag_out(:, :, 59), diag_out(:, :, 60))
  end subroutine drv_zm_mphy

  subroutine drv_zm_conv_main_microp(ncol, nlev, nmodes, nspecmx, &
       time_step, is_first_step, &
       limcnv, no_deep_pbl, tau, alfa, ke, dmpdz, tpert_fix, &
       tpert_fac, tiedke_add, c0_lnd, c0_ocn, num_cin, &
       mx_bot_lyr_adj, trig_dcape, trig_ull, clos_dyn_adj, old_snow, &
       auto_fac, accr_fac, micro_dcs, &
       t_mid, q_mid, omega, p_mid, p_int, p_del, geos, z_mid, z_int, &
       pbl_hgt, tpert, landfrac, t_star_in, q_star_in, &
       nspec, m_accum, m_aitken, m_coarse, l_dust, l_nacl, l_so4, &
       sigmag_aitken, specdens, spechygro, voltonumblo, voltonumbhi, &
       num_a_in, mmr_a_in, dgnum_in, &
       lengath, gather_index, msemax_klev_g, jctop, jcbot, jt, &
       prec, heat, qtnd, cape, dcape, mcon, pflx, zdu, &
       mflx_up, entr_up, detr_up, mflx_dn, entr_dn, p_del_out, &
       dsubcld, ql, rliq, rprd, dlf, microp_out, rice)
    ! Full zm_conv_main with zm_param%zm_microp = .true. -- identical
    ! plumbing to drv_zm_conv_main (zm_conv_driver.F90) plus the
    ! microphysics tuning params, the modal aerosol object (num/mmr/
    ! dgnum in ungathered column space, gathered internally by
    ! zm_conv_main), and the scattered microp_st outputs.
    !
    ! microp_out field order (MICROP_FIELDS in the generator):
    !  1 wu       2 qliq    3 qice    4 qrain   5 qsnow   6 qgraupel
    !  7 qnl      8 qni     9 qnr    10 qns    11 qng    12 sprd
    ! 13 mudpcu  14 lambdadpcu       15 qcde   16 qide   17 qsde
    ! 18 ncde    19 nide   20 nsde   21 dif    22 dsf    23 dnlf
    ! 24 dnif    25 dnsf   26 frz    27 cmel   28 cmei   29 autolm
    ! 30 accrlm  31 bergnm 32 fhtimm 33 fhtctm 34 fhmlm  35 hmpim
    ! 36 accslm  37 dlfm   38 dsfm   39 autoln 40 accrln 41 bergnn
    ! 42 fhtimn  43 fhtctn 44 fhmln  45 accsln 46 activn 47 dlfn
    ! 48 dsfn    49 autoim 50 accsim 51 difm   52 nuclin 53 autoin
    ! 54 accsin  55 hmpin  56 difn   57 trspcm 58 trspcn 59 trspim
    ! 60 trspin  61 accgrm 62 accglm 63 accgslm 64 accgsrm 65 accgirm
    ! 66 accgrim 67 accgrsm 68 accgsln 69 accgsrn 70 accgirn
    ! 71 accsrim 72 acciglm 73 accigrm 74 accsirm 75 accigln
    ! 76 accigrn 77 accsirn 78 accgln 79 accgrn 80 accilm 81 acciln
    ! 82 fallrm  83 fallsm 84 fallgm 85 fallrn 86 fallsn 87 fallgn
    ! 88 fhmrm
    use ppgrid,                only: ppver => pver
    use zm_conv,               only: zm_conv_main, zm_conv_main_init, zm_param
    use zm_aero_type,          only: zm_aero_t
    use zm_microphysics_state, only: zm_microp_st, zm_microp_st_alloc, &
                                     zm_microp_st_dealloc
    use cam_abortutils,        only: endrun

    integer, intent(in) :: ncol, nlev, nmodes, nspecmx
    real(r8), intent(in) :: time_step
    integer, intent(in) :: is_first_step               ! 0/1 flag
    integer, intent(in) :: limcnv                      ! 1-based interface
    integer, intent(in) :: no_deep_pbl                 ! 0/1 flag
    real(r8), intent(in) :: tau, alfa, ke, dmpdz
    integer, intent(in) :: tpert_fix                   ! 0/1 flag
    real(r8), intent(in) :: tpert_fac, tiedke_add, c0_lnd, c0_ocn
    integer, intent(in) :: num_cin, mx_bot_lyr_adj
    integer, intent(in) :: trig_dcape, trig_ull        ! 0/1 flags
    integer, intent(in) :: clos_dyn_adj, old_snow      ! 0/1 flags
    real(r8), intent(in) :: auto_fac, accr_fac         ! microp tuning
    real(r8), intent(in) :: micro_dcs                  ! [m]
    real(r8), intent(in) :: t_mid(ncol, nlev)          ! [K]
    real(r8), intent(in) :: q_mid(ncol, nlev)          ! [kg/kg]
    real(r8), intent(in) :: omega(ncol, nlev)          ! [Pa/s]
    real(r8), intent(in) :: p_mid(ncol, nlev)          ! [Pa]
    real(r8), intent(in) :: p_int(ncol, nlev + 1)      ! [Pa]
    real(r8), intent(in) :: p_del(ncol, nlev)          ! [Pa]
    real(r8), intent(in) :: geos(ncol)                 ! [m2/s2]
    real(r8), intent(in) :: z_mid(ncol, nlev)          ! [m above sfc]
    real(r8), intent(in) :: z_int(ncol, nlev + 1)      ! [m above sfc]
    real(r8), intent(in) :: pbl_hgt(ncol)              ! [m]
    real(r8), intent(in) :: tpert(ncol)                ! [K]
    real(r8), intent(in) :: landfrac(ncol)
    real(r8), intent(in) :: t_star_in(ncol, nlev)      ! DCAPE prev T
    real(r8), intent(in) :: q_star_in(ncol, nlev)      ! DCAPE prev q
    integer, intent(in) :: nspec(nmodes)
    integer, intent(in) :: m_accum, m_aitken, m_coarse
    integer, intent(in) :: l_dust, l_nacl, l_so4
    real(r8), intent(in) :: sigmag_aitken
    real(r8), intent(in) :: specdens(nspecmx, nmodes)
    real(r8), intent(in) :: spechygro(nspecmx, nmodes)
    real(r8), intent(in) :: voltonumblo(nmodes), voltonumbhi(nmodes)
    real(r8), intent(in) :: num_a_in(ncol, nlev, nmodes)          ! [#/kg]
    real(r8), intent(in) :: mmr_a_in(ncol, nlev, nspecmx, nmodes) ! [kg/kg]
    real(r8), intent(in) :: dgnum_in(ncol, nlev, nmodes)          ! [m]

    integer, intent(out) :: lengath
    integer, intent(out) :: gather_index(ncol)         ! 1-based
    integer, intent(out) :: msemax_klev_g(ncol)        ! 1-based, gathered
    integer, intent(out) :: jctop(ncol), jcbot(ncol)   ! 1-based, scattered
    integer, intent(out) :: jt(ncol)                   ! 1-based, gathered
    real(r8), intent(out) :: prec(ncol)                ! [m/s]
    real(r8), intent(out) :: heat(ncol, nlev)          ! [W/kg] scattered
    real(r8), intent(out) :: qtnd(ncol, nlev)          ! [kg/kg/s] scattered
    real(r8), intent(out) :: cape(ncol), dcape(ncol)   ! [J/kg]
    real(r8), intent(out) :: mcon(ncol, nlev + 1)      ! [mb/s] scattered
    real(r8), intent(out) :: pflx(ncol, nlev + 1)      ! [kg/m2/s] scattered
    real(r8), intent(out) :: zdu(ncol, nlev)           ! [1/s] scattered
    real(r8), intent(out) :: mflx_up(ncol, nlev)       ! [mb/s] gathered
    real(r8), intent(out) :: entr_up(ncol, nlev)       ! [1/mb] gathered
    real(r8), intent(out) :: detr_up(ncol, nlev)       ! [1/mb] gathered
    real(r8), intent(out) :: mflx_dn(ncol, nlev)       ! [mb/s] gathered
    real(r8), intent(out) :: entr_dn(ncol, nlev)       ! [1/mb] gathered
    real(r8), intent(out) :: p_del_out(ncol, nlev)     ! [mb] gathered
    real(r8), intent(out) :: dsubcld(ncol)             ! [mb] gathered
    real(r8), intent(out) :: ql(ncol, nlev)            ! gathered
    real(r8), intent(out) :: rliq(ncol)                ! [m/s]
    real(r8), intent(out) :: rprd(ncol, nlev)          ! [kg/kg/s] scattered
    real(r8), intent(out) :: dlf(ncol, nlev)           ! [kg/kg/s] scattered
    real(r8), intent(out) :: microp_out(ncol, nlev, 88) ! scattered microp_st
    real(r8), intent(out) :: rice(ncol)                ! [m/s?] microp_st%rice

    type(zm_aero_t)        :: aero
    type(zm_microp_st)     :: microp_st
    real(r8), pointer      :: t_star(:, :), q_star(:, :)
    real(r8), allocatable, target :: tgt_num(:, :, :)
    real(r8), allocatable, target :: tgt_mmr(:, :, :, :)
    real(r8), allocatable, target :: tgt_dgn(:, :, :)
    integer :: l, m

    if (nlev /= ppver) &
        call endrun('drv_zm_conv_main_microp: nlev /= ppgrid pver')

    ! zm_param subset used by zm_conv_main; zm_conv_intr sets
    ! old_snow = .false. under zmconv_microp (driver input anyway)
    zm_param%tau            = tau
    zm_param%alfa           = alfa
    zm_param%ke             = ke
    zm_param%dmpdz          = dmpdz
    zm_param%tpert_fix      = tpert_fix /= 0
    zm_param%tpert_fac      = tpert_fac
    zm_param%tiedke_add     = tiedke_add
    zm_param%c0_lnd         = c0_lnd
    zm_param%c0_ocn         = c0_ocn
    zm_param%num_cin        = num_cin
    zm_param%mx_bot_lyr_adj = mx_bot_lyr_adj
    zm_param%trig_dcape     = trig_dcape /= 0
    zm_param%trig_ull       = trig_ull /= 0
    zm_param%clos_dyn_adj   = clos_dyn_adj /= 0
    zm_param%old_snow       = old_snow /= 0
    zm_param%zm_microp      = .true.
    zm_param%auto_fac       = auto_fac
    zm_param%accr_fac       = accr_fac
    zm_param%micro_dcs      = micro_dcs
    call zm_conv_main_init(limcnv, no_deep_pbl_in=(no_deep_pbl /= 0))

    ! modal aerosol object: ungathered pointers (num_a/mmr_a/dgnum)
    ! plus the gathered work arrays zm_conv_main fills itself
    aero%scheme = 'modal'
    aero%nmodes = nmodes
    allocate(aero%nspec(nmodes))
    aero%nspec = nspec
    aero%mode_accum_idx  = m_accum
    aero%mode_aitken_idx = m_aitken
    aero%mode_coarse_idx = m_coarse
    aero%coarse_dust_idx = l_dust
    aero%coarse_nacl_idx = l_nacl
    aero%coarse_so4_idx  = l_so4
    aero%sigmag_aitken   = sigmag_aitken
    allocate(aero%specdens(nspecmx, nmodes), aero%spechygro(nspecmx, nmodes))
    aero%specdens  = specdens
    aero%spechygro = spechygro
    allocate(aero%voltonumblo(nmodes), aero%voltonumbhi(nmodes))
    aero%voltonumblo = voltonumblo
    aero%voltonumbhi = voltonumbhi
    allocate(tgt_num(ncol, nlev, nmodes))
    allocate(tgt_mmr(ncol, nlev, nspecmx, nmodes))
    allocate(tgt_dgn(ncol, nlev, nmodes))
    tgt_num = num_a_in
    tgt_mmr = mmr_a_in
    tgt_dgn = dgnum_in
    allocate(aero%num_a(nmodes), aero%mmr_a(nspecmx, nmodes), aero%dgnum(nmodes))
    do m = 1, nmodes
      aero%num_a(m)%val => tgt_num(:, :, m)
      aero%dgnum(m)%val => tgt_dgn(:, :, m)
      do l = 1, nspecmx
        aero%mmr_a(l, m)%val => tgt_mmr(:, :, l, m)
      end do
    end do
    allocate(aero%numg_a(ncol, nlev, nmodes))
    allocate(aero%mmrg_a(ncol, nlev, nspecmx, nmodes))
    allocate(aero%dgnumg(ncol, nlev, nmodes))
    aero%numg_a = 0._r8
    aero%mmrg_a = 0._r8
    aero%dgnumg = 0._r8

    call zm_microp_st_alloc(microp_st, ncol, nlev)

    allocate(t_star(ncol, nlev), q_star(ncol, nlev))
    t_star = t_star_in
    q_star = q_star_in

    gather_index = 0

    call zm_conv_main(ncol, ncol, nlev, nlev + 1, is_first_step /= 0, &
         time_step, t_mid, q_mid, omega, p_mid, p_int, p_del, geos, &
         z_mid, z_int, pbl_hgt, tpert, landfrac, t_star, q_star, &
         lengath, gather_index, msemax_klev_g, jctop, jcbot, jt, &
         prec, heat, qtnd, cape, dcape, mcon, pflx, zdu, &
         mflx_up, entr_up, detr_up, mflx_dn, entr_dn, p_del_out, &
         dsubcld, ql, rliq, rprd, dlf, aero, microp_st)

    microp_out(:, :, 1)  = microp_st%wu
    microp_out(:, :, 2)  = microp_st%qliq
    microp_out(:, :, 3)  = microp_st%qice
    microp_out(:, :, 4)  = microp_st%qrain
    microp_out(:, :, 5)  = microp_st%qsnow
    microp_out(:, :, 6)  = microp_st%qgraupel
    microp_out(:, :, 7)  = microp_st%qnl
    microp_out(:, :, 8)  = microp_st%qni
    microp_out(:, :, 9)  = microp_st%qnr
    microp_out(:, :, 10) = microp_st%qns
    microp_out(:, :, 11) = microp_st%qng
    microp_out(:, :, 12) = microp_st%sprd
    microp_out(:, :, 13) = microp_st%mudpcu
    microp_out(:, :, 14) = microp_st%lambdadpcu
    microp_out(:, :, 15) = microp_st%qcde
    microp_out(:, :, 16) = microp_st%qide
    microp_out(:, :, 17) = microp_st%qsde
    microp_out(:, :, 18) = microp_st%ncde
    microp_out(:, :, 19) = microp_st%nide
    microp_out(:, :, 20) = microp_st%nsde
    microp_out(:, :, 21) = microp_st%dif
    microp_out(:, :, 22) = microp_st%dsf
    microp_out(:, :, 23) = microp_st%dnlf
    microp_out(:, :, 24) = microp_st%dnif
    microp_out(:, :, 25) = microp_st%dnsf
    microp_out(:, :, 26) = microp_st%frz
    microp_out(:, :, 27) = microp_st%cmel
    microp_out(:, :, 28) = microp_st%cmei
    microp_out(:, :, 29) = microp_st%autolm
    microp_out(:, :, 30) = microp_st%accrlm
    microp_out(:, :, 31) = microp_st%bergnm
    microp_out(:, :, 32) = microp_st%fhtimm
    microp_out(:, :, 33) = microp_st%fhtctm
    microp_out(:, :, 34) = microp_st%fhmlm
    microp_out(:, :, 35) = microp_st%hmpim
    microp_out(:, :, 36) = microp_st%accslm
    microp_out(:, :, 37) = microp_st%dlfm
    microp_out(:, :, 38) = microp_st%dsfm
    microp_out(:, :, 39) = microp_st%autoln
    microp_out(:, :, 40) = microp_st%accrln
    microp_out(:, :, 41) = microp_st%bergnn
    microp_out(:, :, 42) = microp_st%fhtimn
    microp_out(:, :, 43) = microp_st%fhtctn
    microp_out(:, :, 44) = microp_st%fhmln
    microp_out(:, :, 45) = microp_st%accsln
    microp_out(:, :, 46) = microp_st%activn
    microp_out(:, :, 47) = microp_st%dlfn
    microp_out(:, :, 48) = microp_st%dsfn
    microp_out(:, :, 49) = microp_st%autoim
    microp_out(:, :, 50) = microp_st%accsim
    microp_out(:, :, 51) = microp_st%difm
    microp_out(:, :, 52) = microp_st%nuclin
    microp_out(:, :, 53) = microp_st%autoin
    microp_out(:, :, 54) = microp_st%accsin
    microp_out(:, :, 55) = microp_st%hmpin
    microp_out(:, :, 56) = microp_st%difn
    microp_out(:, :, 57) = microp_st%trspcm
    microp_out(:, :, 58) = microp_st%trspcn
    microp_out(:, :, 59) = microp_st%trspim
    microp_out(:, :, 60) = microp_st%trspin
    microp_out(:, :, 61) = microp_st%accgrm
    microp_out(:, :, 62) = microp_st%accglm
    microp_out(:, :, 63) = microp_st%accgslm
    microp_out(:, :, 64) = microp_st%accgsrm
    microp_out(:, :, 65) = microp_st%accgirm
    microp_out(:, :, 66) = microp_st%accgrim
    microp_out(:, :, 67) = microp_st%accgrsm
    microp_out(:, :, 68) = microp_st%accgsln
    microp_out(:, :, 69) = microp_st%accgsrn
    microp_out(:, :, 70) = microp_st%accgirn
    microp_out(:, :, 71) = microp_st%accsrim
    microp_out(:, :, 72) = microp_st%acciglm
    microp_out(:, :, 73) = microp_st%accigrm
    microp_out(:, :, 74) = microp_st%accsirm
    microp_out(:, :, 75) = microp_st%accigln
    microp_out(:, :, 76) = microp_st%accigrn
    microp_out(:, :, 77) = microp_st%accsirn
    microp_out(:, :, 78) = microp_st%accgln
    microp_out(:, :, 79) = microp_st%accgrn
    microp_out(:, :, 80) = microp_st%accilm
    microp_out(:, :, 81) = microp_st%acciln
    microp_out(:, :, 82) = microp_st%fallrm
    microp_out(:, :, 83) = microp_st%fallsm
    microp_out(:, :, 84) = microp_st%fallgm
    microp_out(:, :, 85) = microp_st%fallrn
    microp_out(:, :, 86) = microp_st%fallsn
    microp_out(:, :, 87) = microp_st%fallgn
    microp_out(:, :, 88) = microp_st%fhmrm
    rice = microp_st%rice

    call zm_microp_st_dealloc(microp_st)
    deallocate(t_star, q_star)
    deallocate(tgt_num, tgt_mmr, tgt_dgn)
  end subroutine drv_zm_conv_main_microp

end module zm_microp_driver
