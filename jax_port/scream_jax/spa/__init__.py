"""SPA (Simple Prescribed Aerosols).

The modern EAMxx SPA process is a thin wrapper over the generic
DataInterpolation machinery: yearly-periodic linear time interpolation
between monthly file slices, followed by a Dynamic3DRef vertical remap
(source pressure p = PS*hybm + P0*hyam) onto the model p_mid with
P0 (constant) extrapolation outside the source pressure range.
"""
