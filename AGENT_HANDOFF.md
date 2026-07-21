# Handoff: JAX ports of SCREAM & EAMv3 physics (branch `jax_scream`)

Orientation for a coding agent continuing this work (written for a
transfer to a new machine/agent). Everything below is elaborated in the
per-directory docs — read those before editing anything; this file is
the map.

## What this branch contains

Two independent ports of E3SM atmospheric physics to Python/JAX, plus
their validation harnesses, on a fork of E3SM (upstream sources under
`components/` are READ-ONLY except for the already-committed
EAMXX_HAS_PYTHON swap hooks in `components/eamxx`):

1. **`jax_port/` — SCREAMv1 (EAMxx C++) suite**: SHOC, cld_fraction,
   SPA, P3, RRTMGP, assembled driver (`scream_jax/driver.py`), and a
   coupling bridge to the pySEs dynamical core
   (`scream_jax/pyses_bridge.py`, np4 and pg2 physics-grid options;
   the operational GLL↔FV remap lives in `pyses_ext/`).
   Validated in tiers: property tests → golden replay vs pyeamxx →
   in-situ swap ctests (JAX physics embedded in the EAMxx executable).
   Cross-language bit-for-bit was NEVER the goal; divergences are
   characterized against a measured knife-edge/perturbation envelope
   (binary branches amplify roundoff — proven inherent by control
   experiments). Docs: `jax_port/STATUS.md` (per-file provenance +
   AD state), `TEST_HARNESS_DESIGN.md`, `FORMULATION.md` (math),
   `PORTING_PLAN.md`, `dev/README.md` (container recipes).
2. **`jax_port_eam/` — EAMv3 (Fortran) parameterizations**: wv_sat,
   dadadj, geopotential, full gravity-wave suite, complete ZM deep
   convection (incl. convective microphysics), cldfrc2m, tropopause,
   conv_water, EAM P3, EAM RRTMGP, CLUBB slices A–D (E is WIP —
   Fortran golden trusted, JAX files unvalidated; see PORTING_PLAN
   row 12). Validation: f2py-compiled UNMODIFIED Fortran as the golden
   generator (near-BFB, rtol 1e-12 defaults with measured documented
   exceptions). Docs: `jax_port_eam/STATUS.md` (how to resume),
   `METHODOLOGY.md` (rules + environment traps), `PORTING_PLAN.md`.

## Differentiability layer (jax_port; recent work)

- Every process step has NaN-free reverse-mode gradients; defaults are
  **bitwise identical** to the validated primal (non-negotiable
  invariant — the full pytest suite is the guard, currently 214
  passing non-slow: `cd jax_port && .venv/bin/python -m pytest tests/
  -q -m "not slow"`).
- **Smoothing** (`scream_jax/foundation/smoothing.py`): approximation-
  by-identity operators; width=0 returns the exact hard op via a
  Python-level branch. Adopted in cld_fraction, 10 P3 jump families,
  and MCICA (smooth masks → d(flux)/d(cldfrac)). Theory + measurements
  (jump vs kink; the boundary-term/estimator theory) in
  `harness/smoothing_estimator_demo.py` and `harness/p3_jump_gap.py`
  (Stages A/B/C; key results recorded in STATUS.md and commit
  messages `d83e746e83`, `46c628684d`, `e2a445b306`).
- **Code-path flags are STATIC kwargs, never opts/traced values**
  (jit rule): `sed_use_while_loop` / `sed_mode` / `p3_soft_masks` /
  `smooth_width` / `smooth_families`.
- **`ScreamPhysics.step(ad_mode=True)`**: end-to-end traceable suite
  step, bit-exact primal, jax.checkpoint at process/subcycle
  granularity. Measured (harness/rollout_ad_profile.py): grad ≈ 3.5×
  primal, memory flat ~2–2.6 GB to 24 steps (the 6–12 h ML window),
  gradient-norm growth ~×1.09/step. Unsafe-`where` fixes use the
  double-where idiom (sanitize operands BEFORE risky ops; primal-bit-
  neutral).
- Deferred by explicit decision: ES/score-function estimators (kept
  for future adjoint-sensitivity reference tests; design notes in
  STATUS.md), B2 implicit sedimentation (speed-only), expected-overlap
  radiation solver (smooth-MCICA suffices for training gradients).

## Working agreements (do not relax these)

1. Default-path numerics are sacred: any refactor must leave committed
   golden/suite tests bit-stable. One scheme/slice per commit, with
   measured validation evidence in the commit message.
2. `components/` physics sources are read-only. Harness stubs are
   infrastructure-only (abort-only placeholders allowed if provably
   never executed, documented).
3. Physical constants verbatim from source files, never from memory
   (a from-memory Stefan–Boltzmann constant has already caused a bug).
4. Jump-vs-kink triage before smoothing anything: gates on rates that
   vanish continuously at the threshold are kinks — leave hard.
5. Tolerances: state the measured number and its root cause (FMA,
   libm ulp, cancellation floor) whenever looser than 1e-12.

## Environment

- Host: venvs via `uv` in `jax_port/.venv` and `jax_port_eam/.venv`
  (numpy, jax x64, pytest, netCDF4; + scipy, meson, ninja for eam).
  Data at `../e3sm-inputdata/atm/scream` (+ `atm/cam/physprops`).
- Container `scream-dev` (recipes in `jax_port/dev/`): EAMxx builds
  with `EAMXX_ENABLE_PYTHON`/`EAMXX_ENABLE_PYSCREAM` (swap ctests,
  pyeamxx golden generation) and ALL `jax_port_eam` f2py builds
  (`docker exec -w /work/E3SM/jax_port_eam/harness scream-dev python3
  build_<x>.py`). Never build f2py on macOS (dlopen/kind traps —
  documented in jax_port_eam/METHODOLOGY.md).

## Next steps on a supported E3SM machine (agreed plan)

Phase 0: build EAMxx with the Python swap on the machine; re-run
`test-all-eamxx` and the `*_standalone_cpp_vs_jax` ctests there.
Phase 1: EAM Tier-2 golden capture — instrument a real F1850/F2010 run
to dump per-parameterization in/out states, replay through `eam_jax`
(the deferred plan in jax_port_eam/METHODOLOGY.md).
Phase 2: **PGN** (`cime/CIME/SystemTests/pgn.py`) with the JAX-swapped
executable at ne4pg2 — institutional form of our perturbation-envelope
claim.
Phase 3: **MVK** (`cime/CIME/SystemTests/mvk.py`;
`MVK_PS.ne4pg2_oQU480.F2010` in `cime_config/tests.py` ~line 285, via
`cime/scripts/create_test`, report by evv4esm) — the statistical-
climate-equivalence verdict, stock vs JAX-swapped SCREAM. Run TSC
alongside. These validate the DEFAULT hard path only (smoothed/AD
modes are training-time surrogates, not the climate model under test).

## Reading order for a new agent

1. This file → 2. `jax_port/STATUS.md` → 3. `jax_port/FORMULATION.md`
→ 4. `jax_port/TEST_HARNESS_DESIGN.md` + `dev/README.md` →
5. `jax_port_eam/STATUS.md` + `METHODOLOGY.md` → 6. the harness
scripts named above, which are executable documentation of every
measured claim.
