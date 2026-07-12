# JAX Port — Status Ledger

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
| `scream_jax/shoc/second_moments.py` | impl/ `shoc_calc_shoc_vertflux`, `shoc_calc_shoc_varorcovar` (driver chain pending) | d957a16d34 | Claude (Fable 5) | draft |
| `scream_jax/shoc/solver.py` | impl/ `shoc_tridiag_solver` (`vd_shoc_decomp`/`vd_shoc_solve`), `shoc_update_prognostics_implicit` | d957a16d34 | Claude (Fable 5) | draft |

Remaining SHOC kernels (in port order): second/third
moment drivers + boundary conditions, assumed-PDF chain (13 kernels),
pblintd chain (6), `shoc_main` driver, process-interface pre/post.
Whole-scheme validation target: `jax_port/golden/shoc_218x72_dt1800_5steps.npz`.

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
defaulting mirroring `pre_process_tracer_requests`.
