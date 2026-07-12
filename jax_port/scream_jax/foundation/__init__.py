"""Layer-0 foundation: shared physics constants and functions.

Everything here is transcribed from components/eamxx/src/share/ (see each
module's docstring for the exact source file). Port order and rationale:
jax_port/PORTING_PLAN.md section 3.
"""

from . import column_ops, constants, thermo, tridiag
from .saturation import (
    SaturationFcn,
    murphy_koop_svp,
    polysvp1,
    qv_sat_dry,
    qv_sat_wet,
)
