# foundation/ — Porting Notes

Source commit: `d957a16d34` (jax_scream branch). See [../../STATUS.md](../../STATUS.md)
for validation state.

## Deviations from the C++ (all deliberate)

| Where | Deviation | Rationale |
|---|---|---|
| everywhere | EKAT Pack/Mask/WorkspaceManager dropped; `range_mask` args dropped | pack padding doesn't exist in JAX |
| everywhere | `sp(x)` literals treated as double | port targets double precision only |
| `constants.py` | EKAT unit annotations dropped; only SI magnitudes kept | units are compile-time-only in the C++ |
| `saturation.py` | `check_temperature` is a host-side helper, not an in-kernel abort | no aborts inside jit; use in adapters/tests |
| `saturation.py` | branch selection via `jnp.where` instead of masked `Pack::set` | identical arithmetic, both branches evaluated |
| `thermo.py` | team overloads not ported | they are per-level loops over the scalar math |
| `thermo.py::apply_rayleigh_friction` | returns updated (u, v, T) instead of mutating | pure-functional convention |
| `thermo.py::calculate_psl` | `phi_safe` guard inside the `where` | avoids inf/nan in the unselected branch (values match C++; gradients near phi=0 are still not meaningful) |
| `column_ops.py` | CombineMode (y = beta*y + alpha*f(x)) not ported | output-blending is an in-place-output detail; compose instead |
| `tridiag.py` | only the serial `thomas` ordering ported | `cr`/team variants are parallel decompositions of the same math |

## Traps for future kernel ports (learned here)

- **Match binary FP arithmetic, not decimal values**: `Tmelt - 40` is
  233.14999999999998; write expressions, not rounded literals.
- **k=0 is the model top** everywhere; `calculate_z_int` scans from the
  *bottom* (`FromTop=false` with `s0=z_surf`).
- **`column_scan(from_top=False)`** yields interface values whose forward
  delta is **−dx**, not dx.
- **polysvp1 clamps** dt at −80 (constant es below 193.15 K) — intentional.
- **Transcribe quirks verbatim** (e.g. Rayleigh friction's direct
  application of tendency coefficients); golden-data agreement with EAMxx is
  the only correctness criterion, physical plausibility arguments are not.
