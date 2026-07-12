# JAX Port — Status Ledger

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
| `scream_jax/p3/process.py` | `eamxx_p3_process_interface.hpp/.cpp`, `eamxx_p3_run.cpp` (preamble/postamble, wet<->dry, cld-frac max-overlap) | d957a16d34 | Claude (Fable 5) | **kernel-golden** (Tier-1 golden replay passes) |

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

Next for P3: Tier-2 swap test (C++ EAMXX_HAS_PYTHON branch in
`eamxx_p3_process_interface.cpp` + `p3_jax.py` adapter + single-process
CMake test, mirroring SHOC's).

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
