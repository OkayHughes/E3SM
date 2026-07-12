"""P3 microphysics (p3 process).

Transcribed kernel-by-kernel from
components/eamxx/src/physics/p3/impl/*_impl.hpp (master header
p3_functions.hpp). Conventions as in scream_jax.shoc: level axis last,
k=0 model top, leading axes batch columns.
"""

from . import tables
