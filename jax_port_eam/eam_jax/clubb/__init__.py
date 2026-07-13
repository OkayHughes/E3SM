"""JAX port of EAMv3 CLUBB — slice 1 (see PORTING_PLAN.md row 12).

Delivered in this slice, all validated against the eam_clubb_f f2py
harness (the ENTIRE unmodified components/eam/src/physics/clubb stack
compiled with EAM's defines -DCLUBB_CAM -DCLUBB_SGS
-DCLUBB_REAL_TYPE=dp; see harness/build_clubb.py):

- grid.py         CLUBB vertical grid (setup_grid_heights for the EAM
                  host-model path l_implemented=T/grid_type=3) and the
                  four grid operators zt2zm/zm2zt/ddzt/ddzm
                  (l_cubic_interp=.false. linear versions).
- saturation.py   Flatau et al. (1992) SVP polynomials and
                  sat_mixrat_liq/ice (saturation_formula="flatau",
                  the EAM clubb_intr setting).
- tridiag.py      tridag_solve: an exact port of LAPACK dgtsv
                  (reference LAPACK, non-blocked partial-pivoting LU),
                  which is what CLUBB's lapack_wrap.F90 calls for the
                  xp2_xpyp / windm_edsclrm tridiagonal solves.
- pdf_closure.py  The heart of CLUBB: Skx_func, gamma_Skw_fnc,
                  compute_sigma_sqd_w, and pdf_closure on the
                  compile-time iiPDF_ADG1 path (ADG1_w_closure,
                  ADG1/ADG2 responder params, binormal component
                  correlations, chi/eta transform, cloud fraction /
                  liquid water / ice supersaturation fraction, the
                  x'rc' contributions and the th_v moments).

EAM configuration baked into the port scope (all verbatim from
clubb_intr.F90 + model_flags.F90 defaults; goldens use the same):
sclr_dim=0, hydromet_dim=0, iiPDF_type=iiPDF_ADG1 (compile-time),
l_explicit_turbulent_adv_wp3/xpyp=.false., l_stats=.false. (so wp4,
wprtp2, wpthlp2, wprtpthlp are never computed by pdf_closure in EAM
and are not returned here), debug_level=0, l_predict_upwp_vpwp=.false.,
l_cubic_interp=.false., saturation_formula="flatau".

Tunable parameters come from CLUBB's compiled-in defaults with the
EAMv3 phys="default" namelist overrides applied (see
harness/gen_clubb_golden.py for the list and the C2rt->C2thl/C2rtthl
coupling read_parameters applies); the goldens record the full packed
params vector.
"""

from . import grid, pdf_closure, saturation, tridiag  # noqa: F401
