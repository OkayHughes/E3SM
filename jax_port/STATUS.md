# JAX Port — Status Ledger

> **P3 TIER-2 SWAP-TESTED (2026-07-12):** `p3_standalone_cpp_vs_jax_strict`
> and `_ice_feedback` pass in-container: the whole P3 step swapped for
> scream_jax via the embedded-Python bridge over a 5-step coupled run.
> T_mid/qv/rainfrac stay within 1e-6 field-scale error at >99.9% of
> points (max 9e-7 for T_mid); hydrometeor fields carry the documented
> ice-sed knife-edge feedback, bounded (qc max 1.4e-3 at 0.5% of points).
> All 14 p3 standalone tests pass (incl. the original np1-np4 BFB family).
>
> **P3 TIER-1 VALIDATED (2026-07-12):** The complete P3 port (tables,
> ~30 process kernels, parts 1/2/3, sedimentation, homogeneous freezing,
> p3_main, process pre/post) replays the EAMxx golden archive (5 steps,
> dt=1800, 218x72). Warm/thermo fields match at 1e-12..1e-14 field-scale
> relative error at EVERY point (T_mid 7e-14, qc/qr/nc/nr ~1e-13,
> exchanges ~1e-12). Ice-path fields carry bounded knife-edge outliers
> (qi max 6e-4 at 2.4% of points, qm/bm max 1.8e-2): verified NOT a port
> bug — the real C++ ice_sedimentation fed two part2 states differing
> only by 1e-13 FP noise reproduces the same divergence (discontinuous
> rime-density/qsmall branches + exponential cloud-top drainage).
> Direct C++ cross-validation via p3_test_data host wrappers in-container:
> my part1/part2/cloud-rain-ice sed/homog freezing/part3 are each
> **bit-level identical** to the C++ on identical inputs (<= 3e-15,
> ice sed 5e-21). Two bugs found & fixed this way: inverted
> ice_nucleation branch (C++ `any_if_not_log` naming trap) and part3
> eff-radius init values. Run: `pytest jax_port/tests/test_p3_golden.py -s`.
>
> **SHOC TIER-2 SWAP-TESTED (2026-07-12):** `shoc_standalone_cpp_vs_jax`
> passes in-container: the whole SHOC step swapped for scream_jax via the
> embedded-Python bridge, matching the C++ standalone run at <= 2e-5 max
> relative error with <0.004% of points beyond 1e-6 (fraction-aware
> comparator, `jax_port/harness/compare_nc_tolerant.py`). All 9 shoc
> standalone tests and all 5 cld_fraction tests pass.
>
> **SHOC TIER-1 VALIDATED (2026-07-12):** The complete SHOC port (all ~57
> kernels, shoc_main, process pre/post) replays the EAMxx golden archive —
> 5 host steps x 6 subcycles — with max relative error <= 7e-8 on every
> field (most 1e-9..1e-14; single-subcycle agreement is 1e-13..1e-16). The
> only caveat is a documented knife-edge `shoc_ql2 != 0` branch that flips
> on <0.1% of points. Run: `pytest jax_port/tests/test_shoc_golden.py -s`.
>
> **M0 COMPLETE (2026-07-11):** The Docker environment is verified end-to-end
> (configure → full build → ctest), EAMXX_ENABLE_PYTHON works in-container,
> and `cldfrac_standalone_cpp_vs_jax` **passes with bitwise-identical (cprnc)
> output** — the JAX swap toolchain is proven. Recipes: jax_port/dev/README.md.
> Tier-0 pytest suite: 35 passing (`uv venv && uv pip install jax pytest
> numpy; pytest jax_port/tests/`). Next: golden-data generators (M1), then
> shoc/p3 kernel porting.

One row per ported file. **A row here is a claim about provenance and
validation state, nothing more** — `draft` code has been reviewed against the
C++ by its translator but has passed no numerical comparison.

Validation states (cumulative, see [TEST_HARNESS_DESIGN.md](TEST_HARNESS_DESIGN.md)):

| State | Meaning |
|---|---|
| `draft` | translated, self-reviewed, no numerical validation |
| `kernel-golden` | passes Tier-0 golden-data / property tests |
| `swap-tested` | passes Tier-2 in-situ single-process swap comparison |

`Source @` is the E3SM commit the C++/Fortran was read at; re-check rows whose
source files changed upstream before trusting them.

## foundation/

| File | Source file(s) | Source @ | Translator | State |
|---|---|---|---|---|
| `scream_jax/foundation/constants.py` | `components/eamxx/src/share/physics/physics_constants.hpp` | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/foundation/saturation.py` | `components/eamxx/src/share/physics/physics_saturation_impl.hpp` (+ decls in `physics_functions.hpp`) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/foundation/thermo.py` | `components/eamxx/src/share/physics/eamxx_common_physics_functions_impl.hpp` | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/foundation/column_ops.py` | `components/eamxx/src/share/util/eamxx_column_ops.hpp` (pack_size==1 paths) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/foundation/tridiag.py` | `externals/ekat/src/algorithm/ekat_tridiag.hpp` (`impl::thomas_a1x1`) | d957a16d34 | Claude (Fable 5) | draft |

## cld_fraction/ (pilot)

| File | Source file(s) | Source @ | Translator | State |
|---|---|---|---|---|
| `scream_jax/cld_fraction/main.py` | `components/eamxx/src/physics/cld_fraction/cld_fraction_main_impl.hpp` | d957a16d34 | Claude (Fable 5) | **swap-tested** (`cldfrac_standalone_cpp_vs_jax`, cprnc-BFB) |
| `scream_jax/adapters/eamxx/cld_fraction_jax.py` | calling convention of `cld_fraction_numpy.py` / `eamxx_cld_fraction_process_interface.cpp` | d957a16d34 | Claude (Fable 5) | **swap-tested** |

## tms/

| File | Source file(s) | Source @ | Translator | State |
|---|---|---|---|---|
| `scream_jax/tms/main.py` | `components/eamxx/src/physics/tms/impl/compute_tms_impl.hpp` | d957a16d34 | Claude (Fable 5) | draft |

## shoc/ (in progress — kernels ported bottom-up, see PORTING_PLAN.md §4.3)

| File | Source file(s) (impl/ prefix = `components/eamxx/src/physics/shoc/impl/`) | Source @ | Translator | State |
|---|---|---|---|---|
| `scream_jax/shoc/constants.py` | `shoc_constants.hpp` | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/grid.py` | impl/ `shoc_grid`, `shoc_dp_inverse`, `shoc_compute_tmpi` | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/thermo.py` | impl/ `shoc_compute_shoc_vapor`, `shoc_compute_shoc_temperature` | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/interp.py` | impl/ `shoc_linear_interp` | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/tke.py` | impl/ `shoc_check_tke`, `shoc_integ_column_stability`, `shoc_compute_shr_prod`, `shoc_adv_sgs_tke`, `shoc_isotropic_ts`, `shoc_eddy_diffusivities`, `shoc_tke` (driver) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/length.py` | impl/ `shoc_compute_brunt_shoc_length`, `shoc_compute_l_inf_shoc_length`, `shoc_compute_shoc_mix_shoc_length`, `shoc_check_length_scale_shoc_length`, `shoc_length` (driver) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/energy.py` | impl/ `shoc_energy_integrals`, `shoc_energy_fixer`, `shoc_update_host_dse` | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/surface.py` | impl/ `shoc_diag_obklen` | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/second_moments.py` | impl/ `shoc_calc_shoc_vertflux`, `shoc_calc_shoc_varorcovar`, `shoc_diag_second_moments{,_srf,_lbycond,_ubycond}`, `shoc_diag_second_shoc_moments` (driver) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/solver.py` | impl/ `shoc_tridiag_solver` (`vd_shoc_decomp`/`vd_shoc_solve`), `shoc_update_prognostics_implicit` | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/third_moments.py` | impl/ `shoc_compute_diag_third_shoc_moment`, `shoc_clipping_diag_third_shoc_moments`, `shoc_diag_third_shoc_moments` (driver) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/assumed_pdf.py` | impl/ `shoc_assumed_pdf` (driver) + 12 `shoc_assumed_pdf_*` kernels | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/pblintd.py` | impl/ `shoc_pblintd_init_pot`, `pblintd_height`, `pblintd_surf_temp`, `pblintd_check_pblh`, `shoc_pblintd_cldcheck`, `pblintd` (driver) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/main.py` | impl/ `shoc_main` (`shoc_init` + `shoc_main_internal` loop) | d957a16d34 | Claude (Fable 5) | draft (whole-scheme water/energy budget invariants pass) |

| `scream_jax/shoc/process.py` | `eamxx_shoc_process_interface.hpp` (SHOCPreprocess/SHOCPostprocess) + run_impl setup | d957a16d34 | Claude (Fable 5) | **swap-tested** |
| `scream_jax/adapters/eamxx/shoc_jax.py` | py_module_call marshalling in `eamxx_shoc_process_interface.cpp` (this branch) | d957a16d34 | Claude (Fable 5) | **swap-tested** |

The whole SHOC package is **swap-tested** (Tier-2, see banner) on top of
the Tier-1 golden replay. The C++ enablement branch swaps the ENTIRE step
(pre+main+post), matching the validated unit; it rejects the unsupported
extra_shoc_diags/apply_tms/check_flux_state_consistency configs.
Whole-scheme validation target: `jax_port/golden/shoc_218x72_dt1800_5steps.npz`.

## p3/ (in progress)

| File | Source file(s) | Source @ | Translator | State |
|---|---|---|---|---|
| `scream_jax/p3/tables.py` | `impl/p3_init_impl.hpp` (ice table text parser, rain-table numerical integration, dnu) | d957a16d34 | Claude (Fable 5) | **kernel-golden** (computed rain tables match the C++ .dat8 binaries to 1e-14) |
| `scream_jax/p3/table_lookups.py` | impl/ `p3_table3`, `p3_table_ice` (lookup + tri/quadrilinear apply) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/p3/dsd.py` | impl/ `p3_dsd2` (`get_cloud_dsd2`, `get_rain_dsd2`) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/p3/processes_warm.py` | impl/ autoconversion, cloud_rain_acc, droplet/rain self-collection, subgrid_variance_scaling, ice_nucleation, ice_classical_nucleation, cldliq/rain_imm_freezing, calc_rime_density | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/p3/processes_ice.py` | impl/ ice_collection (3), ice_melting, ice_cldliq_wet_growth, ice_deposition_sublimation, ice/liq relaxation timescales, evaporate_rain (+3 helpers) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/p3/conservation.py` | impl/ q/n conservation (6), ice_supersat_conservation, prevent_liq_supersaturation, impose_max_total_ni, incloud_mixingratios | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/p3/cell_average.py` | impl/ back_to_cell_average, get_time_space_phys_variables | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/p3/update.py` | impl/ update_prognostic_ice, update_prognostic_liquid | d957a16d34 | Claude (Fable 5) | **kernel-golden** (via part2 C++ cross-validation) |
| `scream_jax/p3/main_part1.py` | impl/ `p3_main_impl_part1.hpp` | d957a16d34 | Claude (Fable 5) | **kernel-golden** (BFB vs `p3_main_part1_host`) |
| `scream_jax/p3/main_part2.py` | impl/ `p3_main_impl_part2.hpp` (593-line k-loop) | d957a16d34 | Claude (Fable 5) | **kernel-golden** (<=3e-15 vs `p3_main_part2_host` on identical inputs) |
| `scream_jax/p3/main_part3.py` | impl/ `p3_main_impl_part3.hpp` + calc_bulk_rho_rime | d957a16d34 | Claude (Fable 5) | **kernel-golden** (BFB vs `p3_main_part3_host`) |
| `scream_jax/p3/sedimentation.py` | impl/ `p3_find`, `p3_upwind`, `p3_cloud_sed`, `p3_rain_sed`, `p3_ice_sed` (incl. homogeneous_freezing) | d957a16d34 | Claude (Fable 5) | **kernel-golden** (ice sed 5e-21 vs `ice_sedimentation_host`; cloud/rain sed BFB in stage pipeline) |
| `scream_jax/p3/main.py` | impl/ `p3_main_impl.hpp` (init + orchestration + early exits) | d957a16d34 | Claude (Fable 5) | **kernel-golden** (whole-main vs `p3_main_host`: warm fields <=1e-13, ice knife-edge bounded) |
| `scream_jax/p3/process.py` | `eamxx_p3_process_interface.hpp/.cpp`, `eamxx_p3_run.cpp` (preamble/postamble, wet<->dry, cld-frac max-overlap) | d957a16d34 | Claude (Fable 5) | **swap-tested** (`p3_standalone_cpp_vs_jax_*`) |
| `scream_jax/adapters/eamxx/p3_jax.py` | py_module_call in `eamxx_p3_run.cpp` (this branch) | d957a16d34 | Claude (Fable 5) | **swap-tested** |

Tier-0/Tier-1 status: `pytest jax_port/tests/` — 119 passing, including
`test_p3_golden.py` (Tier-1 replay) and `test_p3_main.py` /
`test_p3_sedimentation.py` (conservation & property tests).

Notes:
- `p3/process.py` mirrors the C++: p3_main's `pres`/`dpres` are the DRY
  pressure/thickness; state converts wet->dry before and dry->wet after;
  dz uses full pseudo_density with still-wet qv; rainfrac gets the
  max-overlap of cldfrac_tot from the level above.
- `precip_ice_flux` stays zero — the C++ zero-initializes it and never
  accumulates it (only the liquid flux is filled, by rain sed).
- C++ cross-validation dumpers used for the above (`sed_dump`,
  `part2_dump`, `main_dump`, `stage_dump`) live in the container at
  `/work/*.cpp`, compiled against `p3_test_infra` — rebuild recipe in the
  respective compile commands (link line from
  `CMakeFiles/p3_tests.dir/link.txt`).

P3 is complete through Tier-2.

> **PHYSICS SUITE ASSEMBLED + pySEs BRIDGE (2026-07-12):**
> `scream_jax.driver.ScreamPhysics` assembles the validated chain in
> EAMxx AD order — [mac_mic(shoc -> cld_fraction -> spa -> p3, 6
> subcycles) + rrtmgp] with prescribed CCN (SPA -> P3) and aerosol
> optics (SPA -> RRTMGP) — and replays a multi-process pyeamxx golden
> run within the suite's inherent divergence envelope (control: the C++
> suite vs itself under a 1e-6 T perturbation shows the same
> binary-cldfrac flip growth; a wiring error shows a systematic
> ~100%-of-points signature, which is how this test caught a missing
> aerosol-optics connection). `scream_jax.pyses_bridge` converts pySEs
> state (dry mass/dry mixing ratios, GLL layout) <-> EAMxx fields,
> holds physics-internal persistent state, and returns FT/FU/FV/FQ
> tendencies for `advance_coupling_step`; end-to-end coupler smoke
> passes on real IC thermodynamics. pyeamxx gained group-process
> support (one-line factory registration). Local pytest: 133 passing.
>
> **SPA TIER-2 SWAP-TESTED (2026-07-12):** `spa_standalone_cpp_vs_jax`
> passes in-container at 1e-10 tolerance (all 9 spa standalone tests
> pass). Tier-1 golden replay matches all five prescribed-aerosol
> fields to <=6e-15. The port covers the modern DataInterpolation path:
> yearly-periodic time interpolation, Dynamic3DRef vertical remap
> (p = PS*hybm + P0*hyam, ekat LinInterp + P0 extrapolation) and the
> repairable within-interval postcondition clamps (the data files
> contain g > 1 points that the C++ silently repairs).
>
> **RRTMGP TIER-2 SWAP-TESTED (2026-07-12):** `rrtmgp_standalone_cpp_vs_jax`
> passes in-container: the whole radiation step swapped for scream_jax via
> the embedded-Python bridge over the standalone multi-step run
> (rad_frequency 3, so both update and no-update steps are exercised),
> within 1e-8 field-scale error at >99.9% of points. All 11 rrtmgp
> standalone tests pass (incl. the chunked/not_chunked BFB family), and
> all four swap tests (cld_fraction, shoc, p3, rrtmgp) pass together.
>
> **RRTMGP TIER-1 VALIDATED (2026-07-12):** 3-step golden replay matches
> every computed field to <=2.6e-11 field-scale error (most bit-exact,
> MCICA masks identical, zero knife-edge outliers). Two constant
> mismatches found via replay: SCREAM gravit (9.80616) vs RRTMGP-internal
> grav (9.80665), and stebol 5.670374419e-8.

## rrtmgp/ (complete through Tier-2)

All C++ comparisons below are against the REAL C++/Kokkos RTE+RRTMGP
compiled in the dev container, on identical inputs, via the dumpers in
`harness/cpp_dumpers/` (gasopt_dump.cpp, rrtmgpmain_dump.cpp).

| File | Source file(s) | Source @ | Translator | State |
|---|---|---|---|---|
| `scream_jax/rrtmgp/coefficients.py` | cpp/examples/mo_load_coefficients.h + GasOpticsRRTMGPK::load/init_abs_coeffs (conv::SimpleNetCDF conventions: reversed dims, ints shifted 0-based) | d957a16d34 | Claude (Fable 5) | **kernel-golden** |
| `scream_jax/rrtmgp/gas_optics.py` | cpp/rrtmgp/kernels/mo_gas_optics_kernels.h (Kokkos), mo_gas_optics_rrtmgp.h (compute_gas_taus/source/get_col_dry) | d957a16d34 | Claude (Fable 5) | **kernel-golden** (<=5.6e-16 vs C++, golden/rrtmgp_gasopt_cpp_8x72.npz) |
| `scream_jax/rrtmgp/optical_props.py` | cpp/rte/kernels/mo_optical_props_kernels.h (delta scale, increments) | d957a16d34 | Claude (Fable 5) | **kernel-golden** (via rrtmgp_main comparison) |
| `scream_jax/rrtmgp/cloud_optics.py` | cpp/extensions/cloud_optics/mo_cloud_optics.h (CloudOpticsK LUT path) | d957a16d34 | Claude (Fable 5) | **kernel-golden** (1e-16 vs C++) |
| `scream_jax/rrtmgp/mcica.py` | rrtmgp_conversion.h conv::Random (JSF64) + eamxx get_subcolumn_mask/get_subsampled_clouds | d957a16d34 | Claude (Fable 5) | **kernel-golden** (RNG bit-exact; subsampled tau identical to C++) |
| `scream_jax/rrtmgp/rte.py` | cpp/rte/kernels/mo_rte_solver_kernels.h (SW 2-stream+adding, LW noscat), mo_rte_sw/lw.h drivers | d957a16d34 | Claude (Fable 5) | **kernel-golden** (SW <=4e-13, LW <=2e-15 vs C++) |
| `scream_jax/rrtmgp/interface.py` | eamxx_rrtmgp_interface.hpp (rrtmgp_sw/lw/main + helpers; day-mask instead of day-subset) | d957a16d34 | Claude (Fable 5) | **kernel-golden** (golden/rrtmgp_main_cpp_8x72.npz) |
| `scream_jax/rrtmgp/orbital.py` | share/util/shr_orb_mod.F90 (Berger series, decl, cosz/avg_cosz) + eamxx_trcmix.cpp | d957a16d34 | Claude (Fable 5) | **kernel-golden** (~1e-15 vs Fortran in-container) |
| `scream_jax/rrtmgp/process.py` | eamxx_rrtmgp_process_interface.cpp run_impl | d957a16d34 | Claude (Fable 5) | **swap-tested** (`rrtmgp_standalone_cpp_vs_jax`) |
| `scream_jax/adapters/eamxx/rrtmgp_jax.py` | py_module_call in eamxx_rrtmgp_process_interface.cpp (this branch) | d957a16d34 | Claude (Fable 5) | **swap-tested** |

RRTMGP is complete through Tier-2.

## spa/ (complete through Tier-2)

| File | Source file(s) | Source @ | Translator | State |
|---|---|---|---|---|
| `scream_jax/spa/process.py` | `eamxx_spa_process_interface.cpp`, `share/algorithm/eamxx_data_interpolation.cpp`, `share/remap/vertical_remapper.cpp`, ekat LinInterp | d957a16d34 | Claude (Fable 5) | **swap-tested** (`spa_standalone_cpp_vs_jax`, 1e-10 tol; Tier-1 replay 6e-15) |
| `scream_jax/adapters/eamxx/spa_jax.py` | py_module_call in `eamxx_spa_process_interface.cpp` (this branch) | d957a16d34 | Claude (Fable 5) | **swap-tested** |

## driver / pySEs bridge

| File | Source / role | Translator | State |
|---|---|---|---|
| `scream_jax/driver.py` | EAMxx AD group semantics for [mac_mic + rrtmgp] (tests/multi-process/physics_only/shoc_cld_spa_p3_rrtmgp) | Claude (Fable 5) | **suite-golden** (`tests/test_suite_golden.py` vs multi-process pyeamxx run) |
| `scream_jax/pyses_bridge.py` | pySEs <-> EAMxx state/forcing conversion + persistent coupler (new design, not a transcription) | Claude (Fable 5) | smoke-tested end-to-end (`tests/test_pyses_bridge.py`) |
| `pyses_ext/finite_volume_grid_operational.py` | operational pg-N layer from `components/homme/src/share/gllfvremap_mod.F90` (dp-weighted remaps, CAAS limiter, theta-form T, hydrostatic dp_fv, D-tensor vector remap, tendency/state-asymmetric drivers), written to be appended to pySEs `dynamical_cores/finite_volume_grid.py` | Claude (Fable 5) | **property-tested** (`tests/test_operational_fv.py`, 18 tests, incl. proof that HOMME's constrained-projection FV->GLL operator equals pySEs' reference operator for nf>=2) |

Winds need NO coordinate conversion: pySEs' `horizontal_wind` is
physical lon-lat (u, v) m/s (verified in pySEs `initialization.py`
`wind = jnp.stack((u, v), axis=-1)`, the Coriolis term in
`explicit_terms_theta.py`, and the operator signatures — divergence/
vorticity call `physical_to_contravariant` internally), so components
map 1:1 onto EAMxx `horiz_winds`.

Remaining before production pySEs runs: surface fluxes/albedos are
prescribed (no surface model); omega from pySEs vertical motion; a real
coupled pySEs+scream_jax integration run (requires the pySEs
environment). Physics-grid placement is now an option:
`PysesScreamCoupler` runs physics on GLL columns (np4);
`PysesScreamCouplerPg2` runs it on the pg-N FV grid the way SCREAM
does operationally, via `pyses_ext/finite_volume_grid_operational.py`
(smoke-tested end-to-end in `tests/test_pyses_bridge_pg2.py`; forcing
is returned pre-DSS — supply the `dss` hook or project on the pySEs
side). The pg2 path remaps the DRY state (documented divergence from
SCREAM's wet remap, consistent with pySEs' dry-mass prognostics), and
`ScreamPhysics` gained `spa_col_indices` to source SPA data for a
physics grid that differs from the data file's grid.

## Pending (next in port order — see PORTING_PLAN.md §5)

## Infrastructure

| Item | State |
|---|---|
| Docker build/test environment (`jax_port/dev/`) | **verified**: configure + full build + ctest, C++ and Python paths |
| EAMxx Python swap path (`EAMXX_ENABLE_PYTHON` + adapter) | **verified** via `cldfrac_standalone_cpp_vs_jax` (BFB) |
| cld_fraction JAX swap test (`eamxx/tests/single-process/cld_fraction`) | added on this branch, passing |
| pyeamxx driver (`EAMXX_ENABLE_PYSCREAM`) | **repaired & verified** on this branch (see below) and used for golden data |
| Golden data (`jax_port/golden/*.npz`, via `harness/gen_golden.py`) | **p3 and shoc captured**: 218 cols x 72 lev, dt=1800, 5 steps, all fields, per-step snapshots; validated NaN-free with plausible evolution |
| Tier-0 pytest scaffolding (`jax_port/tests/`) | 35 property/cross-check tests passing |

**pyeamxx repairs on this branch (candidates for upstreaming):** stale include
path; `FieldRequest`/`OutputManager`/`FieldHeader` API drift; **nanobind
stride-units bug** (strides passed in bytes where DLPack requires elements —
corrupted memory on any in-place write through `Field.get()`); field/group
creation now routed through a real `FieldManager` so processes with
(monolithic) group requests work (shoc); single-process tracer-group
defaulting mirroring `pre_process_tracer_requests`; **subfield data-pointer
bug** (`get()` used the parent block's base pointer for monolithic-group
subfields, so qv/qc/tke snapshots all aliased the block start — data
pointer now taken from the typed subview).

## Automatic differentiation hardening (2026-07 session)

Measured baseline and progress live in `harness/ad_probe.py` (grad/jvp/
FD per process on golden states). Current state — every process step
now has running, NaN-free reverse-mode gradients, with primals
bit-identical to the validated code (all changes primal-bit-neutral by
construction; suite at 185 tests green):

| Process | reverse grad | notes |
|---|---|---|
| cld_fraction | runs, ≡0 by physics (binary) | `smooth_width>0` opt-in surrogate gives finite gradients (foundation/smoothing.py, approximation-by-identity: width=0 is bitwise the hard op) |
| shoc | NaN-free, dot-product identity vs FD-validated jvp at 1e-12 | 4 unsafe-where sqrt/cbrt sites fixed |
| p3 | NaN-free through the full step | sedimentation while_loops -> bounded masked scans (bitwise-equal primal; `use_while_loop=True` escape hatch, scan default ~5.5x slower primal; `sed_converged` guards truncation) |
| rrtmgp | NaN-free; FD-vs-jvp 2e-5 | traced path pure jnp (52/52 outputs bit-identical incl. an XLA reciprocal-multiply pitfall); MCICA mask gradients structurally zero w.r.t. cldfrac (expected) |

Remaining AD roadmap: smoothing adoption at further jump sites (P3
qsmall gates, SHOC branches — physics-judgment pass); lax.custom_linear_solve
for tridiag solves; custom_root wrappers for iterative solvers as they
enter (CAAS projection via its scalar dual root per Blondel et al.
2022); an expected-overlap differentiable radiation mode if gradients
w.r.t. cloud fraction through radiation are needed; suite-level
checkpointing strategy for reverse mode through the 6-substep mac_mic
loop.

### B1 uniform-substep sedimentation: implemented; hypothesis refuted
`sed_mode` in {"while","scan","uniform"} (+ `p3_soft_masks`): the
uniform mode removes ALL adaptive-controller discreteness (fixed
dt/M substeps, full-column band, per-substep surface accumulation).
Verified: default bitwise (suite 204), conservation/envelope/
M-refinement tests in tests/test_p3_sed_uniform.py. Stage C
measurement (p3_jump_gap.py --stage c): uniform == scan to 4 digits
for ALL objectives and soft masks change nothing -- the precip
boundary term does NOT live in the sedimentation controller or the
column masks. Remaining candidates: part2's unsmoothed qv/T/ni-driven
jumps (ssat gates, wet-growth qm:=qi snap, ni-gated epsi) and
joint-flip interactions; note also precip's E2 truth is not
FD-converged (strong eps-dependence), so its gap magnitude is itself
uncertain. Recommendation: score-function/ES estimation for
precip-like objectives; smoothing suffices for heat sign and warm.
