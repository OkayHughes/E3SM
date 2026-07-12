# SCREAMv1 → JAX Porting Plan

**Status:** Phase-1 planning document (inventory + dependency ordering).
**Branch:** `jax_scream` (fork of E3SM).
**Companion document:** [TEST_HARNESS_DESIGN.md](TEST_HARNESS_DESIGN.md) — how each ported
component is validated inside SCREAM/EAMxx.

## 1. Goal and scope

Port the SCREAMv1 atmospheric physics suite (the C++/Kokkos EAMxx implementation) to
JAX/Python so the ported parameterizations can be assembled, on the pySEs side, into a
model that is *nearly identical* to SCREAM (bit-for-bit equivalence is explicitly **not**
a goal). The port is incremental: each parameterization is validated against its EAMxx
counterpart by swapping the Python implementation into EAMxx at runtime (see the test
harness document) before it is composed with the pySEs dynamical core.

This document is written to be self-contained for collaborators (and their LLM
assistants) who have not read the EAMxx source. All paths are relative to the E3SM repo
root unless noted.

### 1.1 The SCREAMv1 process suite (what "the physics" is)

The atmosphere process list is defined in
`components/eamxx/cime_config/namelist_defaults_eamxx.xml` (note: there is no
`namelist_defaults_scream.xml`; that name is obsolete). The full atmosphere is:

```
sc_import → homme → physics → sc_export
```

where the `physics` group defaults to:

```
physics       = mac_aero_mic, rrtmgp
mac_aero_mic  = shoc, cld_fraction, spa, p3            # np4 (GLL) grids
mac_aero_mic  = tms, shoc, cld_fraction, spa, p3       # pg2 grids (production SCREAM)
```

`mac_aero_mic` is **subcycled** relative to the driver timestep (`number_of_subcycles`
per grid, e.g. 24 at ne4, 12 at ne30, 1 at ne1024 — same XML file, lines ~589–607).
Radiation runs at reduced cadence (`rad_frequency`).

**Therefore the port scope, in order of the run loop, is:**

| # | Process | EAMxx name | In default suite? | Port? |
|---|---------|-----------|-------------------|-------|
| 0 | Shared foundation (constants, saturation, thermo, column ops) | — | (dependency of all) | **Yes, first** |
| 1 | Turbulent mountain stress | `tms` | pg2 grids only | **Yes** (small) |
| 2 | SHOC turbulence/macrophysics | `shoc` | Yes | **Yes** |
| 3 | Cloud fraction diagnosis | `cld_fraction` | Yes | **Yes** (pilot) |
| 4 | Prescribed aerosol optics | `spa` | Yes (non-noAero) | **Yes** (data plumbing, no kernels) |
| 5 | P3 microphysics | `p3` | Yes | **Yes** |
| 6 | RRTMGP radiation | `rrtmgp` | Yes | **Yes** (largest) |
| — | HOMME dynamics | `homme` | Yes | **No — replaced by pySEs** |
| — | Surface coupling import/export | `sc_import`/`sc_export` | Yes | **No** — pySEs-side responsibility; but the *field contract* (§7) must be honored |
| — | ZM deep convection, gravity-wave drag, COSP, MAM4xx aerosols, nudging, prescribed chemistry, IOP forcing | `zm`, `gw`, `Cosp`, `mam4_*`, `Nudging`, `SPC`, `iop_forcing` | No (opt-in compsets) | **Deferred** (§9) |

### 1.2 What is explicitly out of scope

- **HOMME and the physics↔dynamics remap** (`components/eamxx/src/dynamics/homme/**`):
  pySEs provides the dycore. Note the pg2 path (`eamxx_homme_fv_phys.cpp`, Homme
  `GllFvRemap`) is where GLL↔FV(pg2) remapping happens — pySEs must supply an
  equivalent if targeting pg2; on np4 grids physics runs directly on GLL points.
- **MCT coupler glue** (`components/eamxx/src/mct_coupling/**`).
- **Scorpio/PIO I/O** (`components/eamxx/src/share/io/**`) — used by the test harness,
  not ported.
- **EKAT** (`externals/ekat`; source layout is `externals/ekat/src/{algorithm,core,
  kokkos,pack,...}` — e.g. packs at `src/pack/ekat_pack.hpp`, workspace at
  `src/kokkos/ekat_workspace.hpp`, units at `src/core/ekat_units.hpp`): packs, workspace
  manager, Kokkos types all disappear in JAX. The **one algorithmic exception** is the
  tridiagonal solver `externals/ekat/src/algorithm/ekat_tridiag.hpp` (namespace
  `ekat::tridiag`; three interchangeable solvers: `thomas` — team-parallel and serial
  variants, `cr` — cyclic reduction, `bfb`) used by SHOC's implicit diffusion.
  Reimplement as a Thomas-algorithm `jax.lax.scan` (or `jax.scipy.linalg` banded solve);
  the serial `thomas` overload in that header is the cleanest reference.
- **The AtmosphereDriver / FieldManager infrastructure**: on the pySEs side, replicate
  only the *semantics* — named fields shared between sequentially-run processes,
  process groups, subcycling — not the C++ machinery.

### 1.3 Repo prerequisites

Two submodules hold reference code the port needs to read: `externals/ekat` and
`components/eam/src/physics/rrtmgp/external` (rte-rrtmgp). They are checked out on the
`jax_scream` working tree; on a fresh clone run:

```bash
git submodule update --init --recursive --depth=1
```

## 2. Architecture pattern common to P3/SHOC/TMS/GW/ZM

Each of these packages is a header-only Kokkos template library with an identical layout;
knowing it makes every package legible:

- `<pkg>_functions.hpp` — a single `Functions<Scalar,Device>` struct declaring every
  kernel plus the runtime-options and state structs; `#include`s every impl at the bottom.
- `impl/*_impl.hpp` — **all the physics lives here**. One kernel per file, pure functions
  of Kokkos views. These are the files to translate.
- `eti/*.cpp` — explicit template instantiations. **No physics; ignore.**
- `disp/*.cpp` — GPU "small kernels" dispatch splitting (register pressure). **No unique
  physics; ignore** — JAX handles fusion itself.
- `eamxx_<pkg>_process_interface.{hpp,cpp}` — the `AtmosphereProcess` wrapper: field
  registration (names, layouts, units), unit/variable conversions before/after the kernel
  call, runtime-option parsing from YAML. **Must be read carefully**: the pre/post
  conversions (e.g. wet↔dry mixing ratios, T↔θl) are part of the physics contract and
  must be reproduced in the Python port (either inside the ported process or in the
  pySEs assembly layer).
- `tests/` — Catch2 per-kernel tests + a whole-scheme `*_run_and_cmp` baseline driver;
  `tests/infra/*_test_data.{hpp,cpp}` defines per-kernel I/O structs (these enumerate the
  exact argument list of every kernel — useful as a porting spec).

Kokkos/EKAT concepts and their JAX replacements:

| C++ concept | JAX replacement |
|---|---|
| `ekat::Pack<Real,N>`, `Mask`, `scalarize` | plain `jnp` arrays (packs are a SIMD detail) |
| Team parallelism over columns | `jax.vmap` over the column axis |
| `ekat::WorkspaceManager` scratch | ordinary temporaries |
| `ekat::tridiag::{thomas,cr,bfb}` | Thomas via `lax.scan` (see §1.2) |
| `Kokkos::parallel_scan` (column integrals) | `jnp.cumsum` / `lax.associative_scan` |
| `ekat::units` dimensional constants | plain floats |
| Explicit adaptive sub-stepping loops (sedimentation) | `lax.while_loop` / fixed-trip `lax.fori_loop` |

## 3. Layer 0 — shared foundation (port first; everything depends on it)

All in `components/eamxx/src/share/`:

| File | Provides | Port notes |
|---|---|---|
| `physics/physics_constants.hpp` | All physical constants (Cp=1004.64, Rd=287.042, Rv=461.505, g=9.80616, LatVap=2.501e6, LatIce=3.337e5, Tmelt=273.15, P0=1e5, Karman=0.4, plus P3-specific rime/QSMALL constants) | Transcribe **verbatim** into a constants module. This is the single most-shared file. |
| `physics/physics_functions.hpp` + `physics/physics_saturation_impl.hpp` | Saturation vapor pressure: `polysvp1` (Flatau 1992) and `MurphyKoop_svp` (Murphy & Koop 2005); `qv_sat_dry`, `qv_sat_wet` | Coefficient arrays are load-bearing — copy verbatim. Both P3 and SHOC's assumed PDF call these. |
| `physics/eamxx_common_physics_functions.hpp` (+`_impl.hpp`) | ~30 thermo conversions: exner, θ↔T, θl, virtual T, DSE↔T, dry↔wet mmr, dz, z_mid/z_int, density, vmr↔mmr, psl, rayleigh friction | Used by every process interface's pre/post loops. |
| `util/eamxx_column_ops.hpp` | `ColumnOps`: midpoint↔interface interpolation, midpoint deltas, column scans | Small; maps to slicing + `cumsum`. |
| `physics/eamxx_trcmix.hpp/.cpp` | Background trace-gas profiles (CO2/CH4/N2O/CFC) for radiation | Needed by RRTMGP only. |
| `physics/physics_share.hpp/.cpp` | Thin math wrappers (pow/gamma/erf/…) | Map directly to `jnp`/`jax.scipy.special`; no port needed beyond choosing equivalents. |

*Not ported from share:* `physics_share_f2c.F90` (Fortran-parity test shim),
`physics_test_data.*` (C++ test harness), `util/eamxx_repro_sum_mod.F90` (BFB global sums).

## 4. Per-package inventory

### 4.1 cld_fraction — pilot component (trivial physics, harness already exists)

Dir: `components/eamxx/src/physics/cld_fraction/`

- Physics: `cld_fraction_main_impl.hpp` + `cld_fraction_functions.hpp` — ice cloud
  fraction by thresholding `qi`, total = max(ice, liquid), plus `_for_analysis` variants
  with a looser threshold. A few dozen lines of math.
- Interface: `eamxx_cld_fraction_process_interface.{cpp,hpp}`.
- **Already has NumPy and CuPy ports** (`cld_fraction_numpy.py`, `cld_fraction_cupy.py`)
  swapped in at runtime via YAML (`py_module_name`/`py_module_path`/`py_backend`) — the
  exact mechanism the JAX port uses (see harness doc). Also an ML-emulator sibling
  (`cld_frac_net/`) demonstrating the pytorch path.
- Fortran counterpart: loosely `components/eam/src/physics/cam/cldfrc2m.F90` (the EAMxx
  scheme is one branch of it); no BFB test exists or is needed.
- Tests: `components/eamxx/tests/single-process/cld_fraction/` — already contains
  cpp-vs-numpy and cpp-vs-cupy cprnc comparisons. **Adding `cld_fraction_jax.py` here is
  the first milestone of the whole project.**

### 4.2 tms — turbulent mountain stress (smallest real kernel)

Dir: `components/eamxx/src/physics/tms/`

- Physics: `impl/compute_tms_impl.hpp` (single kernel: surface stress from subgrid
  orography); `tms_functions.hpp`; interface `eamxx_tms_process_interface.{cpp,hpp}`.
- Fortran counterpart: `components/eam/src/physics/cam/trb_mtn_stress.F90`
  (`init_tms`, `compute_tms`) — the C++ is a direct port of it.
- Tests: `components/eamxx/src/physics/tms/tests/compute_tms_tests.cpp` (+`infra/`).
  No single-process AD test exists yet; one must be added for harness Tier-2 (or TMS can
  be validated purely at kernel level).
- Only active on pg2 grids (`surf_drag_coeff_tms` is consumed by SHOC).

### 4.3 shoc — turbulence / macrophysics (~57 kernels)

Dir: `components/eamxx/src/physics/shoc/`

- Master header: `shoc_functions.hpp` (all signatures, `SHOCRuntime` options,
  I/O structs); constants in `shoc_constants.hpp` (mintke, maxlen, w3clip, ustar_min, …).
- Entry point `shoc_main` → `impl/shoc_main_impl.hpp`. Per-column time loop
  (`nadv` subcycles):
  1. `shoc_energy_integrals` (pre)
  2. per subcycle: `check_tke → shoc_grid → compute_shoc_vapor →
     compute_shoc_temperature → shoc_diag_obklen → pblintd → shoc_length → shoc_tke →
     update_prognostics_implicit → diag_second_shoc_moments → diag_third_shoc_moments →
     shoc_assumed_pdf → check_tke`
  3. `update_host_dse → shoc_energy_integrals (post) → shoc_energy_fixer →
     compute_shoc_vapor → shoc_diag_obklen → pblintd`
- Kernel groups in `impl/` (full file list is in the directory; one file per kernel):
  - grid/thermo: `shoc_grid`, `compute_shoc_vapor`, `compute_shoc_temperature`,
    `shoc_linear_interp`
  - energy: `shoc_energy_integrals`, `shoc_energy_fixer`, `update_host_dse`
  - length scale: `shoc_length`, `compute_brunt_shoc_length`, `compute_l_inf_shoc_length`,
    `compute_shoc_mix_shoc_length`, `check_length_scale_shoc_length`
  - TKE: `shoc_tke`, `adv_sgs_tke`, `compute_shr_prod`, `integ_column_stability`,
    `isotropic_ts`, `eddy_diffusivities`, `check_tke`
  - implicit diffusion: `update_prognostics_implicit`, `shoc_tridiag_solver`
    (`vd_shoc_decomp`/`vd_shoc_solve` → **needs the tridiag primitive**),
    `compute_tmpi`, `dp_inverse`
  - second/third moments: `diag_second_shoc_moments` (+`_srf`, `_lbycond`, `_ubycond`,
    core `diag_second_moments`, `calc_shoc_vertflux`, `calc_shoc_varorcovar`),
    `diag_third_shoc_moments`, `compute_diag_third_shoc_moment`,
    `clipping_diag_third_shoc_moments`
  - assumed double-Gaussian PDF (cloud fraction, ql, buoyancy flux):
    `shoc_assumed_pdf` + 12 `shoc_assumed_pdf_*` helper kernels (thl/qw/vv parameters,
    in-plume correlations, tilde-to-real, temperature, qs, s, sgs liquid, liquid variance,
    liquid flux, buoyancy flux)
  - PBL height: `pblintd` + `pblintd_{init_pot,height,surf_temp,check_pblh,cldcheck}`
- Fortran counterpart: `components/eam/src/physics/cam/shoc.F90` (+`shoc_intr.F90`).
  No Fortran remains in EAMxx SHOC; BFB unit tests compare against **stored baseline data
  files** (see harness doc).
- Tests: `components/eamxx/src/physics/shoc/tests/` — ~50 per-kernel tests each with
  `run_property()` (invariants) and `run_bfb()` (baseline compare); whole-scheme
  `shoc_run_and_cmp.cpp`; standalone AD test `components/eamxx/tests/single-process/shoc/`.

### 4.4 p3 — microphysics (~40 kernels + lookup tables)

Dir: `components/eamxx/src/physics/p3/`

- Master header: `p3_functions.hpp` (all signatures; `P3Runtime` — 20 tunables + 7
  flags; `P3PrognosticState/DiagnosticInputs/DiagnosticOutputs/Infrastructure/
  LookupTables/Temporaries` structs).
- Entry point `p3_main` (`impl/p3_main_impl.hpp`) staged as:
  `p3_main_init` → `p3_main_part1` (pre-microphysics: saturation, in-cloud mixing
  ratios, flags) → `p3_main_part2` (main k-loop of all process rates) →
  `cloud_sedimentation` / `rain_sedimentation` / `ice_sedimentation` →
  `homogeneous_freezing` → `p3_main_part3` (conservation, effective radii,
  reflectivity).
- Kernel inventory in `impl/` (one file per kernel): autoconversion,
  back_to_cell_average, calc_liq_relaxation_timescale, calc_rime_density,
  cldliq_imm_freezing, cloud_rain_acc, cloud_sed, conservation (mass/number clipping),
  droplet_self_coll, dsd2 (cloud & rain DSD parameters), evaporate_rain, find (table
  index search), get_time_space_phys_variables, ice_classical_nucleation,
  ice_cldliq_wet_growth, ice_collection, ice_deposition_sublimation, ice_melting,
  ice_nucleation, ice_relaxation_timescale, ice_sed, ice_supersat_conservation,
  impose_max_total_ni, incloud_mixingratios, nc/ni/nr_conservation,
  prevent_liq_supersaturation, rain_imm_freezing, rain_self_collection,
  subgrid_variance_scaling, table3 (rain tables), table_ice (ice table interpolation),
  update_prognostics, upwind (generalized sedimentation flux solver).
- **Data dependencies** (staged from `scream/tables/` by CMake `GetInputFile`):
  - `p3_lookup_table_1.dat-v4.1.1` — the 4-D ice lookup table (read at init)
  - `mu_r_table_vals_v2.dat{4,8}`, `vn_table_vals_v2.dat{4,8}`,
    `vm_table_vals_v2.dat{4,8}`, `revap_table_vals_v2.dat{4,8}` — rain tables. These are
    *computable* (`compute_tables` in `impl/p3_init_impl.hpp`) — the Python port can
    either read the `.dat` files or regenerate them; regeneration is preferred (pure
    function) with a one-time check against the files.
  - hard-coded 16-element `dnu` array (`compute_dnu`).
- Fortran counterpart: `components/eam/src/physics/p3/scream/micro_p3.F90`
  (+`micro_p3_utils.F90`, `micro_p3_interface.F90`). The C++ was ported from this. (An
  EAM-native variant sits in `components/eam/src/physics/p3/eam/`.) No live Fortran
  bridge remains; BFB tests use stored baselines.
- Tests: `components/eamxx/src/physics/p3/tests/` — ~34 per-kernel unit tests,
  `p3_run_and_cmp.cpp` whole-scheme baseline driver, `infra/p3_ic_cases.cpp` canned
  column ICs; standalone AD test `components/eamxx/tests/single-process/p3/` (includes a
  per-step BFB hash extractor); **Python-driven standalone already exists:**
  `components/eamxx/tests/python/pyp3/p3_standalone.py`.

### 4.5 spa — prescribed aerosol (data plumbing, no physics kernels)

Dir: `components/eamxx/src/physics/spa/`

- `eamxx_spa_process_interface.{cpp,hpp}` only. Reads prescribed aerosol optics from
  NetCDF and time-interpolates (`yearly_periodic`) + vertically remaps onto the model
  grid via `components/eamxx/src/share/algorithm/eamxx_data_interpolation.*`. Outputs
  `aero_g_sw`, `aero_ssa_sw`, `aero_tau_sw`, `aero_tau_lw` consumed by RRTMGP.
- No Fortran counterpart (replaces EAM's prescribed-aerosol pathway).
- Port = a Python data reader + time/vertical interpolation. On the pySEs side this is
  I/O infrastructure more than physics; keep it isolated from the numerical suite.
- Tests: `components/eamxx/tests/single-process/spa/`.

### 4.6 rrtmgp — radiation (largest single work item)

Dir: `components/eamxx/src/physics/rrtmgp/`

- **The actual radiative-transfer kernels are NOT in EAMxx.** They live in the external
  submodule `components/eam/src/physics/rrtmgp/external` (repo
  `E3SM-Project/rte-rrtmgp`), consumed via its C++ port under `external/cpp/`.
  With the submodule checked out, the numeric core turns out to be compact
  (~5.4 kLOC total, and much of that is class plumbing):
  - `cpp/rrtmgp/kernels/mo_gas_optics_kernels.cpp` (580 lines) — gas-optics
    interpolation/absorption kernels
  - `cpp/rte/kernels/mo_rte_solver_kernels.cpp` (639) — the SW/LW two-stream and
    no-scattering solvers
  - `cpp/rte/kernels/mo_optical_props_kernels.cpp` (597) — optical-property
    increment/delta-scale ops
  - `cpp/rte/kernels/mo_fluxes_broadband_kernels.cpp` (85) — spectral flux reduction
  - `cpp/extensions/cloud_optics/mo_cloud_optics.h` (1080) — cloud optics LUT/Padé
  - `cpp/rrtmgp/mo_gas_optics_rrtmgp.h` (2447) — gas-optics class: coefficient-file
    ingestion, T/p interpolation setup, source functions
  - plus small classes: `cpp/rrtmgp/mo_gas_concentrations.h`, `cpp/rte/mo_rte_sw.h`,
    `cpp/rte/mo_rte_lw.h`, `cpp/rte/mo_optical_props.h`, `cpp/rte/mo_fluxes.h`,
    `cpp/rte/mo_source_functions.h`, coefficient loaders under `cpp/examples/`.
  The reference Fortran RTE+RRTMGP (Pincus et al.) sits at the submodule top level and
  mirrors this structure 1:1 — `rrtmgp/mo_gas_optics_rrtmgp.F90` +
  `rrtmgp/kernels/mo_gas_optics_kernels.F90`, `rte/mo_rte_{sw,lw}.F90` +
  `rte/kernels/mo_rte_solver_kernels.F90` — useful as a second reference when the C++
  is opaque. EAM's own Fortran driver is `components/eam/src/physics/rrtmgp/*.F90`
  (`radiation.F90`, `cam_optics.F90`, `mcica_subcol_gen.F90`, …).
- EAMxx-owned code is orchestration:
  - `eamxx_rrtmgp_interface.hpp` (~67 KB, header-only): init/load coefficients,
    `rrtmgp_main`, band-by-band surface albedos, broadband fluxes,
    `get_cloud_optics_sw/lw`, `get_subsampled_clouds` (**MCICA random subcolumn
    sampling — a determinism hazard for comparison; see harness doc**),
    `compute_cloud_area`.
  - `eamxx_rrtmgp_process_interface.cpp` (~65 KB): field registration, column chunking,
    gas-concentration population (uses `trcmix`), SPA aerosol optics ingestion,
    heating-rate computation, cloud diagnostics.
  - `rrtmgp_utils.hpp`: heating rate, radiation cadence (`radiation_do`), range checks.
  - `shr_orb_mod_c2f.F90/.hpp`: bridge to CIME's `share/util/shr_orb_mod.F90` — orbital
    parameters, solar declination, cos(zenith). **The Python port must reimplement
    `shr_orb_params` / `shr_orb_decl` / `shr_orb_cosz` from that Fortran** (small,
    self-contained).
- **Data dependencies** (NetCDF, staged from `scream/init/`): k-distributions
  `rrtmgp-data-sw-g112-210809.nc` / `rrtmgp-data-lw-g128-210809.nc` (production) and
  g224/g256 (reference tests); cloud optics `rrtmgp-cloud-optics-coeffs-{sw,lw}.nc`;
  all-sky reference `rrtmgp-allsky.nc`. Filenames are YAML params
  (`rrtmgp_coefficients_file_sw`, etc.).
- **Port strategy decision (flagged, not decided here):** (a) port the RTE+RRTMGP
  kernels themselves to JAX — given the inventory above (~2 kLOC of true numeric
  kernels plus coefficient-file ingestion) this is comparable in size to P3 and is
  mechanical, with a clean 1:1 Fortran cross-reference — or (b) adopt/adapt an existing
  third-party JAX/Python RRTMGP port if one proves compatible with these coefficient
  files. Either way the EAMxx-side orchestration (cloud optics assembly, MCICA, heating
  rates, orbital geometry) must be ported from the files above. Recommend deferring
  RRTMGP to last; the g112/g128 k-distribution files define the spectral discretization
  and must be read as-is.
- Tests: `components/eamxx/src/physics/rrtmgp/tests/` (unit tests + `rrtmgp_run_and_cmp`
  against a generated `rrtmgp-allsky-baseline.nc`); standalone AD test
  `components/eamxx/tests/single-process/rrtmgp/`.

## 5. Dependency graph and recommended port order

```
Layer 0: constants → saturation → common_physics_functions → column_ops → tridiag(JAX)
            │
            ├─ 1. cld_fraction (pilot: exercises the harness end-to-end, ~days)
            ├─ 2. tms          (one kernel; unblocks pg2-suite fidelity)
            ├─ 3. shoc         (needs tridiag, saturation, column_ops)
            ├─ 4. p3           (needs saturation, tables; independent of shoc)
            ├─ 5. spa          (data interpolation only; needed before rrtmgp-with-aerosol)
            └─ 6. rrtmgp       (needs spa outputs, trcmix, orbital code, external kernels)
```

Notes on ordering:
- **cld_fraction first** — not for science value but because the Python-swap harness for
  it already exists and a JAX version validates the entire toolchain (build flags,
  pybind arg marshalling, cprnc comparison) with trivial physics.
- **shoc and p3 are independent of each other** and can be ported in parallel by
  different contributors once Layer 0 merges. Within each, port kernels bottom-up in the
  call-tree order listed in §4.3/§4.4, validating each kernel against golden data before
  assembling `shoc_main`/`p3_main`.
- **Inter-process couplings to respect** (field contract, §7): shoc → (tke, eddy_diff,
  cldfrac_liq, inv_qc_relvar) → p3; cld_fraction → cldfrac_tot → p3/rrtmgp;
  spa → aero_* → rrtmgp; tms → surf_drag_coeff_tms → shoc; p3 → effective radii → rrtmgp.
- The **process-interface pre/post conversions** (wet↔dry mmr, θl↔T, DSE updates) are
  part of each port item, not an afterthought — they are where several historical
  SCREAM bugs lived.

## 6. C++ ↔ Fortran counterpart summary

| EAMxx package | Language of truth | Fortran counterpart (reference/lineage) | Live bridge? |
|---|---|---|---|
| p3 | C++ (`impl/*.hpp`) | `components/eam/src/physics/p3/scream/micro_p3*.F90` | No (stored baselines) |
| shoc | C++ | `components/eam/src/physics/cam/shoc.F90` | No (stored baselines) |
| cld_fraction | C++ (+NumPy/CuPy) | `components/eam/src/physics/cam/cldfrc2m.F90` (loose) | No |
| tms | C++ | `components/eam/src/physics/cam/trb_mtn_stress.F90` | No |
| rrtmgp | external C++/Kokkos (`rte-rrtmgp` submodule) | reference Fortran in same submodule; EAM driver `components/eam/src/physics/rrtmgp/*.F90`; orbital `share/util/shr_orb_mod.F90` (**live**, via `shr_orb_mod_c2f.F90`) | Orbital only |
| spa | C++ (data prescriber) | none | No |
| zm (deferred) | C++ port + **live Fortran bridge** | `components/eam/src/physics/cam/zm/*.F90` via `components/eamxx/src/physics/zm/fortran_bridge/` | **Yes** |
| gw (deferred) | C++ kernels (process is a placeholder, not registered) | `components/eam/src/physics/cam/gw/*.F90` | No |
| cosp (deferred) | Fortran (external COSP2) | `components/eam/src/physics/cosp2/**` via `cosp_c2f.F90` (**live**) | **Yes** |
| mam4xx (deferred) | C++ (`externals/mam4xx`) | lineage: `components/eam/src/chemistry/modal_aero` | No |

Fortran remaining inside `components/eamxx/src` is only: homme interface glue, MCT
coupler glue, scorpio interface, repro-sum, and the ZM/COSP/orbital bridges listed above
— i.e. **nothing in the default suite's numerics is Fortran except the orbital functions**.

## 7. The field contract (physics ↔ pySEs)

Fields the physics consumes from / returns to the dycore+coupler (from `set_grids` in the
process interfaces; layouts: `mid` = nlev midpoints, `int` = nlev+1 interfaces):

- **From dynamics (pySEs must provide):** `T_mid`, `horiz_winds` (u,v; mid),
  `omega` (Pa/s), `p_mid`, `p_int`, `p_dry_mid`, `pseudo_density` (Δp, Pa),
  `pseudo_density_dry`, `phis`, and all advected tracers: `qv, qc, qr, qi, qm, nc, nr,
  ni, bm, tke` (+ the `turbulence_advected_tracers` group).
- **From surface coupler:** `surf_sens_flux`, `surf_evap`, `surf_mom_flux`,
  `sfc_alb_dir/dif_vis/nir` (rad), `surf_lw_flux_up` (rad), `snow_depth_land`/fractions
  as applicable.
- **Back to dynamics:** updated `T_mid`, `horiz_winds`, tracers (HOMME applies them as
  forcings `FT/FM/FQ`; pySEs decides its own application scheme).
- **Physics-internal couplings** are listed per package in §4/§5.

Porting rule: **keep EAMxx field names, units, and (ncol, nlev[+1]) array conventions in
the Python state dict** so that (a) the EAMxx harness adapters are trivial and (b) the
pySEs assembly reads like the EAMxx `atm_procs_list`.

## 8. Design constraints so the Python assembles cleanly on the pySEs side

1. **Functional core / imperative shell.** Every ported process is a pure JAX function
   `new_state, diags = <pkg>_main(state, params, dt)` operating on `(ncol, nlev)` arrays
   with no in-place mutation and no I/O. A thin adapter (harness doc §3) wraps it for the
   in-place numpy calling convention EAMxx's pybind bridge expects. pySEs composes the
   pure functions directly.
2. **float64 by default** (`jax.config.update("jax_enable_x64", True)`). SCREAM's
   production build is double precision; comparisons assume it.
3. **Runtime options as frozen dataclasses/pytrees** mirroring the YAML parameter names
   in `namelist_defaults_eamxx.xml` (e.g. `SHOCRuntime.lambda_low`, `P3Runtime.
   max_total_ni`) so a SCREAM input.yaml can be translated mechanically.
4. **Reproduce the sequential-splitting semantics**: processes run in `atm_procs_list`
   order, each seeing the previous one's updated state; `mac_aero_mic` subcycles with
   dt/number_of_subcycles; radiation on `rad_frequency` cadence. Provide a tiny
   `ProcessGroup` composition helper in Python that mirrors
   `AtmosphereProcessGroup::run_sequential`.
5. **No pack/workspace artifacts** in ported code; column dimension is a `vmap`/leading
   axis, vertical loops become scans where order matters (sedimentation, tridiag,
   integrals) and vectorized ops elsewhere.
6. **Table/data loading isolated** in a `data/` module (P3 tables, RRTMGP NetCDF, SPA
   files) with content-hash checks against the canonical files listed in §4.
7. **One repo-visible Python package** (suggested: `jax_port/scream_jax/` with
   subpackages `foundation/`, `shoc/`, `p3/`, `cld_fraction/`, `tms/`, `spa/`,
   `rrtmgp/`, `adapters/eamxx/`) importable both by the EAMxx test harness and,
   unchanged, by pySEs. Keep adapters (EAMxx in-place shims) strictly out of the
   numerics subpackages.
8. **Provider-agnostic docs**: each subpackage gets a `PORTING_NOTES.md` recording the
   source files translated, deviations taken (e.g. tridiag algorithm choice), and the
   validation evidence — so contributors using any LLM/tooling can pick up mid-stream.

## 9. Deferred packages (inventoried, not scheduled)

- **zm** (Zhang-McFarlane deep convection; `SCREAM%ZM` compsets): full C++ ETI port
  exists (`zm_functions.hpp` + ~22 impl kernels) *plus* a live Fortran bridge used for
  BFB validation. If ported later, the same kernel-level approach as p3/shoc applies.
- **gw** (gravity-wave drag): kernels fully ported to C++ with per-kernel tests, but the
  process interface is an unregistered placeholder — not runnable in the suite yet.
- **cosp** (satellite simulator diagnostics): thin C++ wrapper over the external Fortran
  COSP2 library; diagnostic-only, opt-in.
- **mam4xx** (interactive aerosols; `SCREAM%MAM4xx` compsets): 7 processes over the
  `externals/mam4xx` C++ library — a separate project on its own.
- **nudging, spc, iop_forcing**: data-driven/config processes, port on demand.
