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
| 3 | ZM deep convection core | cam/zm/zm_conv.F90, zm_conv_cape.F90, zm_conv_util.F90, zm_conv_types.F90 | 3400 | pending |
| 4 | ZM transport + intr layer | cam/zm/zm_transport.F90, zm_conv_intr.F90 | 1600 | pending |
| 5 | Gravity-wave drag | cam/gw/{gw_common,gw_oro,gw_convect,gw_front,gw_diffusion,gw_utils}.F90 | 2200 | orographic spine (gw_prof + gw_oro_src + gw_drag_prof ngwv=0) **ported + kernel-golden** (eam_jax/gw.py); spectrum branch (project_tau, LU diffusion, convect/front sources) pending |
| 6 | Cloud fraction | cam/cldfrc2m.F90 | 1100 | pending |
| 7 | Tropopause finder | cam/tropopause.F90 | 1700 | pending |
| 8 | Convective cloud water | cam/conv_water.F90 | ~500 | pending |
| 9 | ZM convective microphysics | cam/zm/zm_microphysics.F90 | 3300 | pending (zmconv_microp=.true. in v3 default) |
| 10 | P3 stratiform microphysics | p3/ (eam variant) | large | pending — start from ../jax_port/scream_jax/p3 (same scheme family; EAM Fortran P3 has aerosol coupling differences) |
| 11 | RRTMGP driver | rrtmgp/ | large | pending — gas optics/solvers reusable from ../jax_port/scream_jax/rrtmgp (same k-distribution algorithm); EAM driver layer differs |
| 12 | CLUBB | clubb/ | very large | pending (single-moment 1.5-order closure; own generator) |
| 13 | Vertical diffusion stack | cam/{vertical_diffusion,eddy_diff,hb_diff,diffusion_solver,trb_mtn_stress}.F90 | ~4000 | pending |
| 14 | Aerosol activation/het. freezing | cam/{activate_drop_mam,hetfrz_classnuc}.F90 | ~1500 | pending |

Out of first scope: MAM4/5 chemistry proper, COSP, WACCM/x, CRM
(MMF), iondrag/qbo, subcolumn machinery, deprecated schemes EAMv3 does
not use in the default configuration (MG1/MG2, UW shallow, HK, RRTMG,
SHOC-in-EAM, stratiform.F90 path).

Shared-with-SCREAM note: EAMxx's SHOC/P3/RRTMGP are C++ descendants of
EAM-family Fortran; where the algorithms coincide, the JAX ports under
`../jax_port/scream_jax/` are the starting point and the f2py golden
distinguishes genuine v3 differences from shared code.
