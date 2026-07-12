# EAMv3 → JAX Porting Methodology

Port of the E3SMv3 atmosphere (EAM, Fortran) parameterizations to
Python/JAX, adapted from the SCREAMv1 port methodology in `../jax_port/`
(see `../jax_port/TEST_HARNESS_DESIGN.md`). This directory is
independent of the SCREAM port; nothing here is imported by
`jax_port/` and vice versa (except where a scheme is shared, noted in
PORTING_PLAN.md).

## How this differs from the SCREAM port

The SCREAM port validated against C++/Kokkos through three tiers:
property tests, golden replay via pyeamxx, and in-situ swap via
EAMXX_HAS_PYTHON ctests. EAM has none of that infrastructure — no
Python driver, no per-process swap hooks, no single-process ctests.
What it has instead is **Fortran + f2py**: the decisive technique from
the SCREAM port (compiling reference-kernel dumpers against the real
build) becomes the *primary* harness here, and improves on it — the
reference kernel is directly callable in-process from Python, so golden
generation, regression sweeps, and debugging all happen in one REPL.

A second difference cuts the other way: SCREAM's C++ and our JAX both
descend from Fortran, so cross-language divergence there was knife-edge
noise through discontinuous branches. Here the port is
**like-for-like** (float64 Fortran → float64 JAX, same formulas), so
Tier-1 tolerances are near-BFB (relative ~1e-12, with documented
exceptions where associativity of reductions differs), and any larger
mismatch is a bug, not an envelope.

## Tiers

**Tier-0 — property tests.** Physical invariants of the JAX port alone:
conservation (energy/water/momentum where the scheme claims it), bounds
(positivity, saturation limits), limiting behavior (zero forcing → zero
tendency, single-regime columns), symmetry/monotonicity where
applicable. Written before or alongside the port; they catch wiring
errors goldens can miss (a golden replay passes if you replay the same
mistake the harness made).

**Tier-1 — kernel-golden via f2py.** The *actual, unmodified* EAM
Fortran source files are compiled with `f2py` into Python extensions,
with minimal **stub modules** replacing infrastructure-only
dependencies (logging, abort, MPI, pbuf, grid queries — never physics).
A generator script drives the Fortran kernel with

  1. *synthetic sweeps* — states constructed to cover every regime
     branch (stable/unstable, saturated/sub-saturated, all
     hydrometeor combinations, top/bottom boundary cases), and
  2. *realistic profiles* — columns from the SCREAM 72-level IC file
     (`screami_unit_tests_ne2np4L72_20220822.nc`; same atmosphere,
     valid thermodynamic states) interpolated to EAM's 72/80-level grid
     as needed,

and records every input and output into a `golden/*.npz` archive with
full metadata (params used, compiler flags, source git SHA). The JAX
port replays the archive under `tests/`.

**Tier-1.5 — chained replay.** Once several schemes from the same
`tphysbc`/`tphysac` sequence are ported, drive the Fortran chain and
the JAX chain from the same initial state for several steps. This
catches interface misunderstandings between schemes (field aliasing,
pbuf persistence, wet/dry mixing-ratio conventions, tendency
application order) that per-kernel goldens cannot.

**Tier-2 — in-situ (deferred).** True in-model validation requires
running EAM on a supported machine (see AGENTS.md; F1850 compset,
ne4pg2_oQU480). The plan: instrument a real run to write each
parameterization's in/out state (nstep ≤ 5), then replay offline — the
same discipline as SCREAM Tier-2 without needing swap hooks. Not
executable on this workstation; every scheme's STATUS notes it as
pending.

## Fortran-specific rules

- **Never modify parameterization source.** Extraction is done by
  compiling the real files plus stubs. If a file cannot compile without
  touching physics code, that scheme's harness notes why and what was
  stubbed instead.
- **Stubs are infrastructure-only**: `shr_kind_mod`, `ppgrid`,
  `cam_logfile`, `cam_abortutils` (abort → `stop`), `spmd_utils`,
  `phys_grid`, pbuf accessors when unavoidable. Physical constants
  (`physconst`, `shr_const_mod`) are copied **verbatim from the E3SM
  source**, never retyped from memory.
- **Runtime parameters** come from
  `bld/namelist_files/namelist_defaults_eam.xml` with `phys="default"`
  attributes (EAMv3 default physics) and are recorded in each golden's
  metadata. Where a scheme's `*_init` routine derives internal
  constants from namelist values, the harness calls the real init.
- **Compilation**: gfortran, `-O2 -ffree-line-length-none
  -fallow-argument-mismatch -std=legacy` (flags recorded in metadata).
  No `-ffast-math` ever. f2py via numpy's meson backend (numpy ≥ 2.0
  needs `meson`+`ninja` in the venv).
- **Array conventions**: EAM arrays are `(pcols, pver)` with level 1 =
  model top; f2py exposes them Fortran-ordered — the harness stores
  goldens C-ordered `(ncol, nlev)` (top→surface) and the JAX ports use
  that layout, matching `jax_port/` conventions.
- **float64 everywhere**: `jax.config.update("jax_enable_x64", True)`
  at import in the package `__init__`.
- **Randomness/iteration**: schemes with internal iteration (e.g.
  saturation adjustment Newton loops) port the *same* iteration count
  and convergence test, not a mathematically equivalent one — that is
  what keeps tolerances near-BFB.

## Workflow per scheme (the loop that worked for SCREAM)

1. Read the Fortran top-to-bottom; write PORT_NOTES in the module
   docstring: inputs/outputs, regime branches, anything surprising.
2. Build the f2py extension (`harness/build_<scheme>.py`), smoke it
   from Python against a hand-checkable case.
3. Generate goldens (`harness/gen_<scheme>_golden.py`): synthetic
   sweeps + IC-file profiles.
4. Port to JAX (`eam_jax/<scheme>/`), replay goldens
   (`tests/test_<scheme>_golden.py`) at 1e-12, plus Tier-0 property
   tests.
5. When a mismatch appears: bisect *inside* the Fortran by calling the
   f2py kernel on sub-inputs (the in-process advantage), not by staring
   at the port.
6. Update STATUS.md, commit. One scheme (or coherent sub-kernel) per
   commit.

## Environment

`uv venv` in `jax_port_eam/`; `uv pip install numpy jax pytest netCDF4
meson ninja`. gfortran from Homebrew (15.x verified). Everything runs
locally — no container needed (the scream-dev container remains
available if a build ever needs the E3SM CMake infrastructure).
