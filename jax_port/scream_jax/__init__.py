"""scream_jax — JAX port of the SCREAMv1 (EAMxx) atmospheric physics suite.

Layout (see jax_port/PORTING_PLAN.md):
    foundation/     shared constants, saturation, thermodynamics, column ops
    <pkg>/          one subpackage per parameterization (shoc, p3, ...)
    adapters/       EAMxx-facing in-place shims (py_module_name targets)

All numerics are pure functions of (ncol, nlev[+1]) float64 arrays; nothing in
this package performs I/O or in-place mutation outside adapters/.
"""

import jax

# SCREAM production builds are double precision; every comparison tolerance in
# the test harness assumes it (TEST_HARNESS_DESIGN.md section 3).
jax.config.update("jax_enable_x64", True)
