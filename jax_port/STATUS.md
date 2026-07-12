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

## Pending (next in port order — see PORTING_PLAN.md §5)

## Infrastructure

| Item | State |
|---|---|
| Docker build/test environment (`jax_port/dev/`) | **verified**: configure + full build + ctest, C++ and Python paths |
| EAMxx Python swap path (`EAMXX_ENABLE_PYTHON` + adapter) | **verified** via `cldfrac_standalone_cpp_vs_jax` (BFB) |
| cld_fraction JAX swap test (`eamxx/tests/single-process/cld_fraction`) | added on this branch, passing |
| Tier-0 pytest scaffolding (`jax_port/tests/`) | 35 property/cross-check tests passing; golden-data generators not yet written |
