# EAMv3 Parameterization Inventory & Porting Order

EAMv3 default physics (from `bld/namelist_files/namelist_defaults_eam.xml`
`phys="default"` and `physpkg.F90`): CLUBB (macroturbulence/shallow),
P3 (stratiform microphysics), ZM with convective microphysics
(`zmconv_microp=.true.`) + MCSP (deep convection), RRTMGP (radiation),
MAM4/5 (aerosols, via chemistry), orographic + convective + frontal
gravity-wave drag, dry adiabatic adjustment, tropopause, `cldfrc2m`
(cloud fraction), `conv_water`, vertical diffusion support layers.

Order below is feasibility x impact. "Lines" = Fortran source lines.

| # | Scheme | Source (components/eam/src/physics/) | Lines | Status |
|---|--------|--------------------------------------|-------|--------|
| 1 | Saturation vapor pressure foundation | cam/wv_sat_methods.F90, cam/wv_saturation.F90 | 1300 | **ported + kernel-golden** (eam_jax/wv_sat.py) |
| 2 | Dry adiabatic adjustment | cam/dadadj.F90 | 135 | **ported + kernel-golden** (eam_jax/dadadj.py) |
| 3 | ZM deep convection core | cam/zm/zm_conv.F90, zm_conv_cape.F90, zm_conv_util.F90, zm_conv_types.F90 | 3400 | **ported + kernel-golden**: dilute CAPE core (eam_jax/zm_cape.py) AND the full zm_conv_main driver (eam_jax/zm_conv.py: DCAPE trigger + gather, zm_cloud_properties updraft/downdraft plumes with zm_calc_fractional_entrainment + zm_downdraft_properties, zm_closure, mass-flux limiting/scaling, zm_calc_output_tend, prec/rliq). SCOPE at the time of this row: zm_microp=.false. (goldens + port used the in-plume condensation branch; Fortran harness ran the same flag, abort-only interface stub documented it). SCOPE SINCE LIFTED by row 9: eam_jax/zm_conv.py now also validates zm_microp=.true. (the EAMv3 default) against the eam_zm_microp_f harness, which compiles the REAL zm_microphysics + activation stack. MCSP + zm_aero inactive by construction. Replay <= 1e-10 rel everywhere except heat at 2e-9 (measured 1-ulp container-vs-host libm exp() through the (exp(x)-1)/x mass-flux cancellation; see tests/test_zm_conv.py). **zm_conv_evap also ported + kernel-golden** (eam_jax/zm_conv.py zm_conv_evap + the REAL cloud_fraction.F90 cldfrc_fice, compiled unmodified in the eam_zm_evap_f harness with new infrastructure-only stubs, stubs/zm_intr_stubs.F90): Sundqvist below-cloud evaporation, snow melt/production + fusion heating, surface prec/snow fluxes; BOTH old_snow branches goldened (old_snow=T is the zm_param_t default consistent with zm_microp=F; old_snow=F has prdsnow=0 without zm_microphysics, so zero snow flux and an unreachable flxsnow<=flxprec fixer, ported verbatim + validated as no-op). zmconv_ke=2.5e-6 (phys="default" dyn=se microphys=p3; single ke, no zmconv_ke_lnd exists). Goldens: exact intr sequence (zm_conv rerun config a -> physics_update with qneg3 at qmin(Q)=1e-12 -> evap) on the 42 golden columns + 40-column synthetic sweep (freezing-level crossings, saturated/dry, cld=0/1, prec-limited, negative prdprec). Replay 1e-12 rel (atol floors at ulp-cancellation residue, measured <=1.4e-20); Tier-0: rain-flux-in == precip+evap out water closure, evaporation ceases at saturation, snow flux in [0, flxprec], no-rain no-op (tests/test_zm_evap.py) |
| 4 | ZM transport + intr layer | cam/zm/zm_transport.F90, zm_conv_intr.F90 | 1600 | **zm_transport.F90 ported + kernel-golden** (eam_jax/zm_transport.py): zm_transport_tracer (chat geometric/arithmetic interface averages, in-cloud updraft/downdraft recursions, min(chat,const)-limited flux divergence, dry-mixing-ratio dpdry branch, zm_microp negative-tracer conservation fixer — BOTH fixer branches goldened, it needs no zm_microphysics code) and zm_transport_momentum (Gregory PGF momcu=momcd=0.4 compile-time constants, no namelist switch; B&B03 KE-dissipation heating; verbatim odd parenthesization of the k=2 downdraft seed). Goldens drive the real routines with physically consistent plume fluxes from re-running eam_zm_conv_f on the zm_conv golden profiles (config a, 18/42 columns triggered) + 6 synthetic tracers (positive-definite contrasting structures, dry-type, negative-valued, uniform no-op, fixer-bait) + sheared/uniform winds. Replay: tracer rtol 1e-12 (atol 1e-24 denormal netflux-clip residue), momentum icw/pg 2e-12, wind_tend/seten 1e-11 (measured 1.9e-12; cancelling mu*(wiu-wi) differences; see tests/test_zm_transport.py). Tier-0: tracer column-mass conservation, untriggered columns exact zero, fixer positivity, uniform tracer/wind no-ops, KE energy closure. **zm_conv_intr.F90 zm_conv_tend sequence ported** as the pure-JAX driver eam_jax/zm_intr.py zm_tend on plain arrays: zm_conv_main (incl. the two compute_dilute_cape DCAPE-trigger calls) -> physics_update (t+=dt*s/cpair, vapor qneg3=max(q,1e-12); use_mass_borrower=F default) -> zm_conv_evap on state1 -> physics_update -> zm_transport_momentum on the untouched winds -> zm_transport_tracer_1 (lq=cnst_is_convtran1, fake_dpdry=0) -> summed ptend_all, plus the mcon mb/s->kg/m2/s conversion (midlevels only, verbatim) and the pbuf outputs as plain arrays. SKIPPED (documented in zm_intr.py): MCSP (phys="default" would set zmconv_MCSP_heat_coeff=0.3; not ported, mcsp_enabled=F scope), aerosol/convproc paths + zm_conv_tend_2, pbuf/history plumbing, zm_microp=F throughout. **Tier-1.5 chained replay green** (harness/gen_zm_chain_golden.py + tests/test_zm_chain.py): 42 columns x 3 steps dt=1800 of the full Fortran chain (eam_zm_conv_f + eam_zm_evap_f + eam_zm_transport_f + eam_geopotential_f) with tendency feedback, end-of-step T_STAR/Q_STAR (physpkg bookkeeping), destabilizing inter-step forcing, geopotential_t zmid/zint refresh; JAX replays the identical chain through zm_tend. Trigger sets exact all steps (18/21/24 columns); final state t 1.6e-12, zmid 7.2e-14, tr2 4.6e-11 rel (<=1e-10 target); qv 1.1e-10, tr1 1.9e-10, u 2.2e-10, v 1.6e-10 (asserted 1e-9): accumulated feedback of the documented 1-ulp libm exp() heat noise (per-step s_tend 7.6e-10/3.4e-9/2.2e-9), not an interface error -- prec <=1.9e-12, snow bitwise |
| 5 | Gravity-wave drag | cam/gw/{gw_common,gw_oro,gw_convect,gw_front,gw_diffusion,gw_utils}.F90 + cam/vdiff_lu_solver.F90 | 2200 | **ported + kernel-golden**, both branches. Orographic spine (gw_prof + gw_oro_src + gw_drag_prof ngwv=0) in eam_jax/gw.py; full spectrum (ngwv=pgwv=32) in eam_jax/gw_spectrum.py: gw_beres_src (convective, mfcc lookup table as plain input; synthetic table in the golden since the Beres04 file isn't on disk), gw_cm_src (frontal, fav init), gwd_compute_stress/tendencies over -pgwv:pgwv, gwd_project_tau, gwd_precalc_rhoi with gw_ediff + gw_diff_tend on the exact vd_lu_decomp/solve port, momentum_energy_conservation. EAMv3 defaults recorded in golden metadata (dc=2.5, taubgnd=2.5e-3, effgw_beres=0.35, effgw_cm=1.0, frontgfc=1.25e-15, hcf=10, hdsf=0.5, use_gw_convect_old=.true.; both heating-depth variants goldened). do_molec_diff=.false. only (no WACCM in v3 defaults). Replay <= 1e-12 rel per kernel + chained beres->cm gw_tend sequence; sole atol exception: chained qtgw at abs 1e-24 ((qnew-q)/dt cancellation noise floor ulp(q)/dt ~ 3e-25, measured 1.2e-25). Tier-0: spectral symmetry, MEC momentum/energy closure, tracer column-mass conservation, zero-source => zero tendency (tests/test_gw_spectrum.py) |
| 6 | Cloud fraction | cam/cldfrc2m.F90 | 1100 | **ported + kernel-golden** (eam_jax/cldfrc2m.py): astG_PDF/astG_RHU (single+vector bodies are identical, ported once) and aist for all iceopt 1-7 (EAMv3 default 5); rhmini/rhmaxi injected per call (readnl is masterproc-only), other params via cloud_fraction stub + real cldfrc2m_init |
| 7 | Tropopause finder | cam/tropopause.F90 | 1700 | **ported + kernel-golden** (eam_jax/tropopause.py): twmo core (Reichler WMO lapse-rate), climatology fallback (climo field becomes a plain (ncol,12)+days+calday input; the harness fills the real module-private tropp_p_loc through tropopause_read_file via a pio stub + the real interpolate_data regridder, column lats on climo nodes for exact injection), hybridstobie, and the tropopause_find primary+backup dispatch. Exactly the EAMv3 production combinations: default TWMO+CLIMATE (aer_rad_props, prescribed_volcaero, tropopause_output) and HYBSTOB+CLIMATE (mozart chemistry, modal_aero_wateruptake). NOT ported (no EAMv3 callers): analytic, stobie, wmo, e90/e90_3d, findChemTrop. Replay: level indices exact, P/T/Z <= 1.8e-15 rel (fma/libm ulp; tol 1e-12) |
| 8 | Convective cloud water | cam/conv_water.F90 | ~500 | **ported + kernel-golden** (eam_jax/conv_water.py): conv_water_4rad, all supported modes — conv_water_in_rad=1 (EAMv3 phys="default"; arithmetic average) and 2 (emissivity log-average); mode 0 = routine never called (cloud_diagnostics.F90). Both zm_param%zm_microp branches goldened (default zmconv_microp=.true.; the branch is plain arithmetic on dp_icimr), plus pergro_mods repartition and RK-vs-P3 kabsi, and the NaN-FICE-guarded COSP sh_cldliq/sh_cldice outputs. pbuf/phys_state become plain arrays (physics_buffer/physics_types/phys_control/cam_history/constituents stubs; harness compiles the full zm_conv stack for zm_param). Replay: mode 1 + microp at 1e-13 rel; mode 2 rtol 1e-12 + atol 1e-17 (measured 2.25e-8 rel on 2/1152 points of ~7e-14 kg/kg: 1-ulp cross-libm exp/log amplified by the log(1+alpha*w)/alpha cancellation near ic_limit; abs error bounded by ulp(1)/|alpha| < 4e-19). Tier-0: nonnegativity, microp mode-independence, no-convection = stratiform limit, COSP sum identity + NaN guard, pergro scope (tests/test_conv_water.py) |
| 9 | ZM convective microphysics | cam/zm/zm_microphysics.F90 | 3300 | **ported + kernel-golden** (eam_jax/zm_microphysics.py), AND zm_conv_main now validates the full zmconv_microp=.true. path (EAMv3 default), lifting the row-3 microp=.false. scope. ACTIVATION STRATEGY: the REAL aerosol-activation modules compile cleanly against the existing infrastructure stubs, so the harness (harness/build_zm_microp.py, eam_zm_microp_f) compiles cam/activate_drop_mam.F90 (AR-G 2000 droplet activation; real actdrop_mam_init exactly as ndrop.F90 calls it) and cam/nucleate_ice_conv.F90 (Liu-Penner 2005) UNMODIFIED and both are ported outright — no prescribed-activation shim anywhere. Only bulk ndrop_bam is an abort-only stub (never executed: EAMv3 runs MAM modal aerosols; the port raises on scheme != 'modal'). New infrastructure stubs (stubs/zm_microp_stubs.F90): time_manager get_step_size (driver-set dt) + the ndrop_bam placeholder. MAM4 mode/species data (sigmag, densities, hygroscopicities, dgnum lo/hi -> voltonumb) are physprop stand-ins recorded in golden metadata (pure data, replayed identically on both sides). CPP variant: no MODAL_AERO_4MODE_MOM/RAIN_EVAP_TO_COARSE_AERO (3-species coarse dust weight). Ported: zm_mphyi, zm_mphy (2-internal-iteration in-plume 2-moment scheme: KK2000 autoconversion/accretion with auto_fac=7/accr_fac=1.5, Ferrier ice autoconv at micro_dcs=150e-6, riming/graupel conversion, Hallett-Mossop, Bigg immersion + dust contact freezing, Rotstayn Bergeron, homogeneous freezing of cloud/rain at -40C, size-distribution clamps, plume vertical velocity (ECMWF KE), per-species conservation scalings, vertical integration with falling-precip re-averaging, niadj/ncadj redistribution), actdrop_mam_calc (single-updraft branch only — zm_mphy calls with sigw=0; maxsat scan reproduced exactly; only fn feeds back), nucleati_conv (hetero/hf), zm_microphysics_adjust, and the zm_conv.F90 microp branches (itnum=2, freezing-coupled h_upd/cu, tot_frz top feedback, evp<=rprd, sprd/frz evp removal, pflxs<=pflx fixer, jt>=jlcl disable, latice*frz heating + dif/dnlf/dnif/dsf/dnsf detrainment, mx-jt<2 closure cut, prec/rliq/rice). BFB LESSONS (measured): gfortran default -ffp-contract=fast FMA-fuses hygro accumulation across the razor-edge hygro>1e-10 threshold (MAM pom/bc hygro is EXACTLY 1e-10) -> harness builds with -ffp-contract=off; XLA rewrites division by scalar constants AND broadcast divisors into reciprocal-multiplies (incorrectly rounded) -> _divs/_bdiv helpers force true IEEE division at every inexact-divisor site; mode/species reductions made sequential (Fortran loop order, no pairwise); single-precision literal promotions reproduced (mg0=1.6E-10, qric/qgic 1.e-8 thresholds, the two post-integration (rhosu/rho)**0.54 exponents). Goldens (gen_zm_microp_golden.py): config k = direct zm_mphy on 40 synthetic plumes (warm/mixed/ice families — ice family goldens the faithful glaciated-from-base quirk where kqi never fires and qi stays 0 — aerosol sweep 0.05..20x, dynamics sweep incl. eps0=0 and cond=0 no-ops and a cmei-gap libase-branch column); configs m/n = full zm_conv_main zm_microp=.true. on the 42-column zm_conv sounding set (CAPE + DCAPE triggers, 18/11 columns, column-varying aerosol incl. dust-heavy). Replay: kernel default 1e-12 rel (qc/qi/nc/ni/wu/lamc/pgam/frz + most diags); precip chain (qr/nr/qni/ns/qg/ng + partners) 5e-10 (measured 2.1e-10: sub-QSMALL cancellation residues through fallspeed pow chains and conservation-ratio branches); fhmrm atol 5e-10 (denormal-qr sign flips at the -40C rain-freeze branch, measured zero impact on frz/sprd/qg); bergnn/trspcm 5e-12 (kernel) / 1e-10+5e-9 (main; flux-difference cancellation); main outputs at the zm_conv tolerances (heat 2e-9 libm exp). Tier-0: warm plumes liquid-only, positivity, no-op columns exact, aerosol->droplet monotonicity, no freezing above 0C, number-mass consistency, water closure prec+rliq==-int(qtnd), column snow<=precip production, rice==detrained ice integral (tests/test_zm_microp.py, 309 tests). NOT ported: bulk-aerosol branches (ndrop_bam + bulk immersion/contact freezing; unreachable in EAMv3), the spectral-updraft branch of actdrop_mam_calc (sigw>0; zm_mphy never uses it), fm/fluxn/fluxm outputs (unused), zm_microphysics_register/history (pbuf/history plumbing). NEXT STEP (documented): zm_intr.py zm_tend still runs zm_microp=F — wiring microp through the intr chain (old_snow=F evap with prdsnow=sprd, dnlf/dnif/dnsf pbuf outputs, chained replay) is the remaining integration |
| 10 | P3 stratiform microphysics | p3/ (eam variant) | large | **ported + kernel-golden** (eam_jax/p3/: main, main_part1/2/3, sedimentation, processes_warm/ice, conservation, update, dsd, cell_average, table_lookups, tables, saturation, constants), COPIED from ../jax_port/scream_jax/p3 and adapted line-by-line against p3/eam/micro_p3.F90 (the truth; every adaptation in the package docstring): EAM exner naming (multiplies T), qv_sat = MurphyKoop es/(p-es) (scream: /p), mu_r_constant=0 rain tables (scream: 1; p3_init_b's ventilation weight has a SINGLE-precision 1.e-6 literal), namelist tuning params with ACTIVE subgrid_variance_scaling (autoconv 3.19 / accretion 1.15 / immersion 2.0 exponents), EAM get_rain_dsd2 (lammin from p3_max_mean_rain_size=5mm, log/gamma nr recompute, logn0r from nr), EAM ice_deposition_sublimation (no inv_dt limiter; qiberg*p3_wbf_coeff), conservation sequence ending prevent_ice_overdepletion (EAM-only) -> ice_supersat_conservation, wet growth WITHOUT the max(0,..) qccol/qrcol clamp, CNT_couple at qsmall=1e-14 (use_hetfrz_classnuc=T is the EAMv3 default; frzimm/frzcnt/frzdep pure inputs), Cooper N_nuc>=1e-20 + do_Cooper_inP3 add-on, update_prognostic_ice with qi_wetDepos (feeds precip_total_tend), post-sed homogeneous_freezing on th/exner (sequential qc-then-qr), NEW ice_complete_melting at f32(273.15)+2, p3_mincdnc=20e6 floor, cflx/rflx/sflx dt_sub-integrated fluxes (precip_ice_flux stays 0), max_total_ni=500e3 (scream: 740e3), f32-promoted incloud_limit/precip_limit, bfb_cbrt = pow(x,1/3) NOT libm cbrt, 49-slot p3_tend_out reconstruction. HARNESS (build_p3.py, eam_p3_f): REAL micro_p3.F90 + micro_p3_utils + wv_sat_scream + physics_utils + scream_abortutils compiled unmodified; new infrastructure stubs (stubs/p3_stubs.F90): phys_control (use_hetfrz_classnuc flag+setter) and print-only debug_info (only referenced by the never-called polysvp1 error path). Table: p3_lookup_table_1.dat-v4.1.1 (EAMv3's default 4.1.2 expects atm/cam/physprops, not on local disk; pure input data recorded in the golden — reader identical). GOLDEN (gen_p3_golden.py): 40 cols x 72 lev, 5 regime families (warm rain / mixed-phase / ice-only / heavy precip / evap+edge incl. clear, nucleation-only, sub-qsmall wisps, T>2C melt, T<-40C homog-freeze columns), 5 configs x 8 chained steps: a=EAMv3 default (hetfrz=T, predict_nc=T, 3 steps), b=hetfrz F (Bigg), c=predict_nc F (Cooper+nccnst), d=do_Cooper_inP3, e=prescribed CCN; EAMv3 phys="default" params in metadata; plus qv_sat probe grids, Fortran rain tables, and drv_cloud/rain/ice_sed sub-kernel goldens. REPLAY (tests/test_p3_eam.py, 281 tests): process rates (tend slots 2-35) <=8.6e-8 rel (worst slot 26 ni_sublim, ratio of leading-edge tinies; others <=4.1e-10); sed sub-kernels 5e-12 (measured 4e-13); liquid-path fields strict 1e-10 rel (measured <=1.7e-11: qc/nc/qv/mu_c/lamc/eff_radius_qc/exchanges/cflx/rflx/precip_liq_flux/precip_liq_surf). Ice-path full-chain fields carry a documented KNIFE-EDGE ENVELOPE (>=70-85% of points at 1e-9 + max-abs <= per-field fraction of field scale, measured maxima half the bounds): sedimentation's nstep=int(Co_max+1) and qsmall front gates amplify sub-1e-11 part2 noise (libm tgamma/pow, container glibc vs XLA) to O(1e-3) locally; PROVEN not a port bug — the Fortran sed fed the JAX part2 state reproduces the golden mismatch magnitudes exactly (9.05e-4/1.34e-3 qi/ni_sed a/0) while JAX-vs-Fortran sed on identical inputs agrees to 4e-13. Tier-0: column water closure vs surface precip <1e-12, cp*T-Lv*(qc+qr)-Ls*qi column energy closure <1e-10, positivity, clear-column exact no-op, nucleation-only column gains ice. NOT ported: do_precip_off=T goldens (branch ported, iop-only option), p3_main's debug check_values (print-only), SCREAM_CONFIG_IS_CMAKE paths, prescribed-CCN file reader (get_prescribed_CCN is interface-level; nccn_prescribed is a plain input). NEXT: micro_p3_interface (micro_p3_tend pbuf/ptend layer + get_cloud_fraction) and chained CLUBB->P3 replay |
| 11 | RRTMGP driver | rrtmgp/ | large | **ported + kernel-golden** (eam_jax/rrtmgp/): the full EAM radiation layer on the FORTRAN path EAMv3 builds (bld/configure rad=rrtmgp without -rrtmgpxx: external/{rte,rrtmgp} f90 kernels + f90/rrtmgp_interface). COPIED from ../jax_port/scream_jax/rrtmgp and re-validated against the Fortran: coefficients.py (loader identical; netCDF names/index conventions match mo_load_coefficients), gas_optics.py (kernels identical; get_col_dry no-latitude branch g0=grav; mo_rrtmgp_constants defaults never overridden), rte.py (+ flux_net = dn - up outputs, mo_fluxes), optical_props.py (delta_scale ADAPTED to the Fortran kernel: no tau>eps gate, max(eps,.) division guards; LW 1-scalar delta-scale is a no-op). FRESH: mcica.py (EAM's KISS RNG — kissvec.c + ShrKissRandGen + mcica_subcol_gen.F90 max-random overlap, seeds = frac(pmid) of the bottom FOUR layers, changeseed=1 for both SW and LW — completely different from EAMxx's per-cell JSF64), cloud_optics.py (EAMv3 defaults gammadist liquid [mu/lambda bilinear, in-code /0.9970449e3 unit fix, lamc>0 and 1e-80 gates] + mitchell ice [d_eff grid, dei==0 gate], CAM lininterp with extrap_method_bndry incl. the DECREASING per-mu lambda axis, cam_optics snow combine_properties), driver.py (radiation_tend sequencing: set_rad_state RRTMG tint + (lwup/stebol)^0.25 surface + extra level [pmid=0.5*pint_top, pint=1.01 Pa], temperature clip to k-dist limits, RRTMG->RRTMGP band reorder, MCICA sampling, optics clipping, day compression, set_albedo 14286 cm^-1 split, clr=gas+aer / all+=cloud through rte, dF*g/dp heating, surface exports with the hard-coded band-10 NIR/vis split). HARNESS (build_rrtmgp.py, eam_rrtmgp_f, container has netcdf-fortran): REAL external kernels + interface + driver layer (radconstants, assertions, radiation_state/utils, cloud_rad_props, ebert_curry, slingo, mcica_subcol_gen, cam_optics, interpolate_data) + REAL share/RandNum (C KISS) compiled unmodified; new infrastructure stubs stubs/rrtmgp_stubs.F90 (own masterproc=.true. copies for the netCDF reads; abort-only pbuf/aer_rad_props/rad_cnst_get_gas/time_manager; physics_state/cam_in_t data containers; infnan via ieee_arithmetic); drivers/rrtmgp_core.F90 holds documented VERBATIM TRANSCRIPTIONS of radiation.F90 PRIVATE helpers (reordered, set_albedo, set_daynight_indices, radiation_driver_sw/lw bodies, set_net_fluxes, export_surface_fluxes) validated end-to-end by the full-step golden. COEFFICIENTS: EAM defaults rrtmgp-data-{sw-g112,lw-g128}-210809.nc (identical files staged in e3sm-inputdata/atm/scream/init) + atm/cam/physprops {F_nwvl200_mu20_lam50_res64_t298_c080428, iceoptics_c080917}.nc fetched from the LCRC inputdata server; sha256 in golden metadata. GOLDEN (gen_rrtmgp_golden.py): 16 cols x 72 lev synthetic profiles (3 T families, 3 cloud decks + snow, banded aerosols, 11 day/5 night), configs a (snow, EAMv3 path) / b (no snow, all-day) / c (clear-sky); kernel goldens (gas optics SW/LW on the padded rad grid, cloud optics both snow branches, MCICA masks 112+128 gpt, sampled gpt optics, rrtmgp_run_sw/lw solver calls) + full drv_rad_step. REPLAY (tests/test_rrtmgp_eam.py, 20 tests): gas optics <=1.9e-15 rel, cloud optics <=6.6e-16 rel, MCICA masks + sampling BIT-EXACT, solvers <=4.1e-13 rel, step fluxes <=4.1e-13 rel / 6.1e-11 W/m2 abs, qrs/qrsc <=2.7e-12 rel, qrl 1.6e-11 / qrlc 1.4e-10 rel (clear-sky flux-difference cancellation; abs <=2.8e-13); asserted with ~10x headroom. Tier-0: column-integral heating == net-flux convergence, all-sky OLR <= clear-sky, zero-sun => zero SW (incl. the nday==0 branch), clear-sky==all-sky in the clear config (1-ulp increment note), surface soll+sols+solld+solsd == fsds, direct <= total down, MCICA clear/overcast limits. NOT ported (documented in driver.py): aerosol optics computation (MAM/aer_rad_props — plain band-array inputs as radiation_tend receives them), gas vmr pbuf plumbing (get_gas_vmr; plain (8,ncol,nlev) input in active_gases order), coszrs/orbital factors (plain inputs; shr_orb port exists in ../jax_port/scream_jax/rrtmgp/orbital.py), icall diagnostic loop, spectralflux pbuf, COSP/HIRS/history, radheat_tend. NEXT: chained CLUBB/P3->radiation replay once cloud_diagnostics (iclwp/lambdac pbuf layer) is ported |
| 12 | CLUBB | clubb/ | very large | **slices A+B+C delivered** (A: grid + saturation + dgtsv tridiag + pdf_closure ADG1; B: pdf_closure_driver — zt+zm double call, trapezoidal rule, clip_rcm, compute_cloud_cover — kernel-golden vs a verbatim extraction validated BITWISE end-to-end against advance_clubb_core; C: advance_xp2_xpyp — ADG1 upwind semi-implicit tridiagonal advance of rtp2/thlp2/rtpthlp/up2/vp2 with pr1/pr2, hole filling and all clips — plus clip_covars_denom, goldened DIRECTLY against the real public routine; clip_covars_denom bitwise, advance <= 8.5e-14 rel); next slice: **D** (Lscale/tau infrastructure); see the CLUBB slicing plan below |
| 13 | Vertical diffusion stack | cam/{vertical_diffusion,eddy_diff,hb_diff,diffusion_solver,trb_mtn_stress}.F90 | ~4000 | pending |
| 14 | Aerosol activation/het. freezing | cam/{activate_drop_mam,hetfrz_classnuc}.F90 | ~1500 | pending |

## CLUBB slicing plan (row 12)

SURVEY (2026-07): EAM builds CLUBB with `-DCLUBB_CAM -DCLUBB_SGS
-DCLUBB_REAL_TYPE=dp` (bld/configure); clubb_intr.F90 config:
grid_type=3 + l_implemented=T (host heights, zi_g index 1 = surface,
zt ghost below), sclr_dim=0, hydromet_dim=0, theta0=300,
ts_nudge=86400, saturation="flatau", l_uv_nudge=F, l_input_fields=F,
l_host_applies_sfc_fluxes=F, edsclr_dim = pcnst - MAM modes/species -
NUMLIQ + 2 (clubb_expldiff=T default), debug_level 0, l_stats=F.
clubb_readnl defaults: clubb_vert_avg_closure=T -> l_vert_avg_closure=
l_trapezoidal_rule_zt/zm = l_call_pdf_closure_twice = T;
clubb_ipdf_call_placement: CORRECTED in slice B — the EAMv3
phys="default" namelist value is **2** (ipdf_post_advance_fields;
namelist_defaults_eam.xml line ~1944 overrides the generic default 1),
so advance_clubb_core runs the sub-advances FIRST and calls
pdf_closure_driver LAST (its outputs pass through to clubb_intr
untouched — which is what made the slice-B end-to-end validation
possible); pdf_closure_driver itself is placement-independent.
clubb_timestep=300 (dtime; hdtime = host dt / cld_macmic_num_steps=6).
iiPDF_type=iiPDF_ADG1 is COMPILE-TIME in pdf_closure_module.
advance_clubb_core sequence (advance_clubb_core_module.F90) under
EAMv3 ipdf_call_placement=2 (ipdf_post_advance_fields; the pre-call
block at line ~961 is SKIPPED): inline sigma_sqd_w (Skw/gamma +
zt2zm(zm2zt) smoothing, line ~1012) -> compute_mixing_length/Lscale
(+ perturbed Lscale avg) + brunt_vaisala + tau/Cx_fnc_Richardson +
calc_surface_varnce -> advance_xm_wpxp (band solver, mono_flux_limiter,
uses the PREVIOUS step's pdf_params_zm/wp2rtp/rtpthvp/... inouts)
-> clip_rcm -> advance_xp2_xpyp (tridiag) -> clip_covars_denom ->
advance_wp2_wp3 (band solver) -> [advance_xp3: l_advance_xp3=F ->
xp3_LG_2005_ansatz] -> advance_windm_edsclrm (tridiag + fill_holes)
-> pdf_closure_driver LAST (Skx/sigma_sqd_w -> pdf_closure at zt
[+ at zm: l_call_pdf_closure_twice] -> trapezoidal rule -> clip_rcm
-> compute_cloud_cover -> l_use_cloud_cover substitution; its outputs
pass through to clubb_intr untouched).
Tunables: 83-element params vector, CLUBB defaults + EAMv3
phys="default" clubb_param_nl overrides (recorded in goldens;
Skw_denom_coef=0 is the CLUBB_CAM compile-time default).

HARNESS (slice A, reused by ALL later slices): eam_clubb_f compiles
the ENTIRE unmodified clubb/ directory (75 files, topo-sorted in
harness/build_clubb.py) + real shr_kind/shr_const + infrastructure
stubs (masterproc=F skips all namelist reads; stubs/clubb_stubs.F90
adds phys_control use_od_fd=F) + `-llapack`; NETCDF/GFDL/MKL/SPMD
undefined so those paths compile out.  drivers/clubb_driver.F90
exposes the EAM-exact setup (read_parameters(-99) via
clubb_param_readnl init -> setup_clubb_core_api -> per-column
setup_grid_heights_api + setup_parameters_api) and per-kernel entry
points; any internal routine can be goldened by adding a drv_.

Slices, in order (A delivered by this row's update):
- **A (DONE)**: grid_class (setup_grid_heights l_implemented path +
  zt2zm/zm2zt/ddzt/ddzm) -> eam_jax/clubb/grid.py; flatau saturation
  -> saturation.py; tridag_solve = exact LAPACK dgtsv port ->
  tridiag.py; Skx_func + gamma_Skw_fnc + compute_sigma_sqd_w +
  pdf_closure (ADG1: ADG1_w_closure, responder params, binormal comp
  corrs, chi/eta transform, cloud/ice-supersat fractions,
  x'rc'/th_v moments, rcp2) -> pdf_closure.py.  Goldens
  golden/clubb_{grid,tridag,sat,pdf_closure}.npz (40 regime-sweeping
  columns x 73 levels: stable/convective-cloudy/extreme-skew/
  degenerate-wp2/degenerate-xp2/saturated/cirrus/randomized).
  MEASURED: pdf_closure BITWISE on everything not passing through
  erf/exp (all 47 pdf_params but rc/cloud-fraction, wp2rtp/wp2thlp,
  rc_coef, mixt_frac); erf/exp-tail fields cloud_frac abs<=2.3e-16,
  rc_i abs<=1.7e-19 (rel 4.5e-9 at ~2e-11 leading edges), downstream
  x'rc'/thv moments at the same abs level (tests at rtol 1e-12 +
  ~10x-measured atol floors).  Grid/sat replay <=1e-14 rel; tridag:
  nrhs=1 mostly bitwise, multi-RHS <=4.7e-14 elementwise (Ubuntu
  liblapack built with FMA contraction; ill-conditioned probe family
  amplifies to 5.2e-11 solution-scaled, asserted at 1e-9).  Tier-0:
  realizability (variances>=0, |corr|<=0.99, fractions in [0,1],
  rcm/rcp2>=0, mixt_frac clip), clear-column exact zeros, saturated
  limit cf~1, convective moist-updraft wprcp>0 + buoyancy-flux
  enhancement, Skw-negation mirror symmetry, degenerate-wp2 collapse
  (tests/test_clubb.py, 19 tests).
- **B (DONE)**: pdf_closure_driver (Skx/gamma/sigma_sqd_w assembly
  with the zt2zm(zm2zt(.)) smoothing, the zt+zm double pdf_closure
  call l_call_pdf_closure_twice=T, trapezoidal_rule_zt/zm, clip_rcm,
  compute_cloud_cover + l_use_cloud_cover substitution;
  rcm_supersat_adj=0 and the rel-humidity diag are dead under EAMv3
  l_rcm_supersat_adj=F) -> eam_jax/clubb/pdf_closure_driver.py.
  HARNESS: build_clubb.py generates clubb_pdf_extract.F90 at build
  time — a VERBATIM-BY-CONSTRUCTION extraction (mechanical text copy)
  of the private pdf_closure_driver + trapezoid/cloud-cover/clip
  helpers rewrapped as a public module; drivers/clubb_driver.F90 adds
  drv_pdf_closure_driver, drv_advance_clubb_core (the REAL public
  advance_clubb_core driven as clubb_tend_cam does) and
  drv_set_eam_flags(ipdf=2, expldiff=T); drivers/
  clubb_driver_helpers.F90 packs pdf_parameter (f2py can't take
  derived-type args).  VALIDATION: with EAMv3 placement=2 the final
  internal pdf call's outputs pass through advance_clubb_core
  untouched and its inputs are the returned advanced prognostics —
  replaying the extraction on that state reproduces rcm, cloud_frac,
  wpthvp, wp2thvp, rtpthvp, thlpthvp, rcp2_zt, thlprcp, wprcp,
  ice_supersat_frac, rcm_in_layer, cloud_cover and both 47-field pdf
  param sets BITWISE (asserted for 8 columns inside
  gen_clubb_golden.py).  GOLDEN golden/clubb_pdf_driver.npz: 40
  regime columns (moments on zm, means on zt; fam 4 = one-level dry
  notch + large rtp2 that TRIGGERS clip_rcm and the
  compute_cloud_cover cloud-top/base branches) + 8 advanced-state
  cases; 29 defined outputs (wp4/wprtp2/wpthlp2/wprtpthlp/
  Skw_velocity/rtm_frz/thlm_frz are intent(out) garbage under EAM
  config and not goldened; rtm is bitwise-passthrough, asserted).
  MEASURED: rc_coef/rcm_supersat_adj bitwise; everything else 1-2 ulp
  of its own scale — grid-op FMA-contraction ulps at thl scale
  (thlm_zm/thl_1/thl_2 abs<=1.2e-13 = 1.4e-16 rel; Fortran objects
  -ffp-contract=off, XLA contracts) and slice-A erf/exp libm tails
  amplified by knife-edge zeta sensitivity in the zm pdf call +
  the cloud-cover division (cloud fraction leading edges
  abs<=1.7e-11 at values <=1e-6; rcm-family abs<=1.3e-17,
  thvp-family abs<=1.1e-14); tests at rtol 1e-12 + ~10x-measured
  atol floors.  Tier-0: trapezoid convex-combination bounds +
  constant invariance, clip_rcm exactness (rtm - 2^-52 floor at 0),
  cloud-cover boost/realizability/interior-no-op, driver-level
  realizability incl. rcm==rcm_in_layer + cf==min(1,cloud_cover)
  substitution identities, bone-dry-column exact zeros, clip_rcm
  branch-coverage assertion (tests/test_clubb_driver.py, 13 tests).
- **C (DONE)**: advance_xp2_xpyp + clip_covars_denom ->
  eam_jax/clubb/advance_xp2_xpyp.py.  EAMv3-active path (flags
  asserted against the Fortran module state via drv_xp2_config):
  l_iter_xp2_xpyp=T, l_upwind_xpyp_ta=T, ADG1 semi-implicit ta
  (l_explicit_turbulent_adv_xpyp=F), l_single_C2_Skw=F,
  l_C2_cloud_frac=F, l_min_xp2_from_corr_wx=F, l_hole_fill=T,
  l_clip_large_rtp2=T (rtp2 <= 0.5*rtm^2), l_tke_aniso=T,
  gamma_over_implicit_ts=1.5, up2/vp2 sponge damping OFF (settings
  never assigned in an EAM build; static zero-init, read back via the
  driver).  Ported: xp2_xpyp_lhs/rhs/uv_rhs (diffusion_zm +
  term_ma_zm + upwind xpyp_term_ta_pdf lhs/rhs + dp1 + the
  over-implicit RHS balance terms, term_tp/term_pr1/term_pr2 with
  Fortran-exact fp order incl. vp2's swapped shear-term order),
  slice-A dgtsv tridiag solves (up2/vp2 share one matrix, 2 RHS),
  pos_definite_variances -> fill_holes_vertical "zm"
  (sequential 5-point windows + global pass, strict left-to-right
  vertical_avg accumulation), clip_variance / large-rtp2 cap / 1000
  cap / clip_covar, clip_covars_denom.  DEAD Fortran args dropped and
  documented (incl. the discovery that under l_upwind_xpyp_ta=T the
  ENTIRE zt-side ADG1 coefficient set — a1_zt, wp2_zt, wp3_on_wp2_zt,
  all coef/term zt arrays — is computed but never read).
  HARNESS: drv_advance_xp2_xpyp drives the REAL PUBLIC routine
  directly (no transcription; pdf_implicit_coefs_terms zero-filled,
  unread under ADG1) + drv_clip_covars_denom + drv_xp2_config /
  drv_param_indices_xp2 getters (nu2/nu9_vert_res_dep captured from
  the module).  GOLDEN golden/clubb_xp2.npz: 16 cases whose moments/Kh
  are REAL one-step drv_advance_clubb_core-advanced states + 24
  synthetic stress columns (degenerate variances, Cauchy-Schwarz
  violations, sign-alternating wp3_on_wp2 for both upwind arms,
  0.5*rtm^2 and 1000 caps, hole bait via strongly negative forcings,
  advection/dissipation heavy); coverage asserted at generation
  (rtp2 floor 384 pts / cap 188 / up2 cap 3 / rtpthlp clipped 1669 /
  65 clip_covars columns).  MEASURED: clip_covars_denom BITWISE;
  advance outputs 59-73% bitwise, rest few-ulp (max rel: rtpthlp
  8.5e-14 at ~4e-8 leading edges, others <= 3.5e-15; XLA FMA
  contraction through the term assemblies + dgtsv); tests at
  rtol=1e-12, atol=0 — NO loosened fields.  Tier-0: post-clip
  realizability (floors/caps/corr bound), dissipation-to-threshold
  fixed point (incl. the up2/vp2 dp1+pr1 algebra x*=w_tol^2),
  production==dissipation steady state (thlp2 = thl_tol^2 +
  P*tau/C2thl invariant), exact no-op on dead paths (hole-fill gate,
  in-bound clips), fill_holes mass conservation + untouched
  boundaries, clip lands exactly on the bound
  (tests/test_clubb_xp2.py, 15 tests).
  NAMELIST CORRECTIONS recorded in this slice (EAMV3_OVERRIDES in
  gen_clubb_golden.py; all clubb goldens regenerated, A/B outputs
  unaffected since none of these feed pdf_closure): clubb_C14 = 2.5
  (phys="default"; CLUBB default 1.0 — REQUIRED by slice C's up2/vp2
  terms), c_K10 = c_K10h = 0.35 (defaults 0.6/1.0; first used slice
  D), and wpxp_L_thresh = 100.0 NOT 60.0 (build-namelist's defaults
  lookup is case-insensitive, so <clubb_wpxp_l_thresh phys="default">
  100.0 beats the attribute-less 60.0 entry; first used slice E).
- **D (NEXT)**: Lscale/tau infrastructure: compute_mixing_length (+
  the l_avg_Lscale perturbed calls), calc_brunt_vaisala_freq_sqd,
  compute_Cx_Fnc_Richardson, calc_surface_varnce, tau_zm/tau_N2
  (l_stability_correct_tau_zm=T).  Note c_K10/c_K10h = 0.35 (EAMv3
  overrides recorded in slice C).
- **E**: advance_xm_wpxp: band solver (dgbsv/dgbsvx port — the *svx
  condition-estimate path needs dgbtrf/dgbcon/dgbtrs, check which
  l_use_* flag EAM hits), C6/C7 skewness functions, semi-implicit
  matrix, monotonic flux limiter (mono_flux_limiter.F90, tridiag).
- **F**: advance_wp2_wp3 (5-band solver, l_damp_wp2_using_em=F path).
- **G**: advance_windm_edsclrm (tridiag, l_do_expldiff_rtm_thlm=T
  extra 2 eddy scalars, fill_holes_vertical) + wind nudging off.
- **H**: advance_clubb_core shell (sponge/splat, sequencing,
  clipping order) + multi-substep Tier-1.5 chained replay; then the
  clubb_tend_cam intr layer (state conversion, macmic subcycling,
  pdf pbuf persistence) as eam_jax clubb_intr driver.

Out of first scope: MAM4/5 chemistry proper, COSP, WACCM/x, CRM
(MMF), iondrag/qbo, subcolumn machinery, deprecated schemes EAMv3 does
not use in the default configuration (MG1/MG2, UW shallow, HK, RRTMG,
SHOC-in-EAM, stratiform.F90 path).

Shared-with-SCREAM note: EAMxx's SHOC/P3/RRTMGP are C++ descendants of
EAM-family Fortran; where the algorithms coincide, the JAX ports under
`../jax_port/scream_jax/` are the starting point and the f2py golden
distinguishes genuine v3 differences from shared code.
