"""Port of eam/src/physics/cam/geopotential.F90 (geopotential_t).

PORT_NOTES
----------
- Computes zi (interfaces, pver+1) and zm (midpoints) above the
  surface from t/q hydrostatically, bottom-up. zi[surface] = 0.
- Two hydrostatic-weight branches: fvdyn (Lin-Rood, hkl from log
  pressures) and the default SE/Eulerian branch (hkl = pdel/pmid,
  hkk = hkl/2). E3SM runs the SE dycore, so fvdyn=False is the
  production path; both are ported, only the SE branch is
  golden-validated (the Fortran's fvdyn comes from dycore_is('LR')).
- geopotential_dse exists in the Fortran but is private and marked
  "Do not use this subroutine" — not ported.
- The bottom-up recurrence zi[k] = zi[k+1] + rog*tv*hkl is a pure
  reversed cumulative sum — vectorized here (identical FP result:
  same additions in the same order per column).
- rair/zvir are per-point arrays in the Fortran signature (variable
  composition support); scalars broadcast fine.
"""

import jax.numpy as jnp


def geopotential_t(piln, pmln, pint, pmid, pdel, rpdel, t, q, rair,
                   gravit, zvir, fvdyn=False):
    """Heights above the surface from temperature (hydrostatic).
    Returns (zi, zm): (ncol, pver+1) and (ncol, pver)."""
    t = jnp.asarray(t, dtype=jnp.float64)
    rair = jnp.broadcast_to(jnp.asarray(rair, dtype=jnp.float64),
                            t.shape)
    zvir = jnp.broadcast_to(jnp.asarray(zvir, dtype=jnp.float64),
                            t.shape)
    rog = rair / gravit

    if fvdyn:
        hkl = piln[:, 1:] - piln[:, :-1]
        hkk = 1.0 - pint[:, :-1] * hkl * rpdel
    else:
        hkl = pdel / pmid
        hkk = 0.5 * hkl

    tv = t * (1.0 + zvir * q)
    dzi = rog * tv * hkl                       # interface-to-interface
    # zi[k] = sum of dzi below (reverse cumsum); zi[pver] = 0
    zi_above = jnp.cumsum(dzi[:, ::-1], axis=1)[:, ::-1]
    zi = jnp.concatenate([zi_above, jnp.zeros_like(zi_above[:, :1])],
                         axis=1)
    zm = zi[:, 1:] + rog * tv * hkk
    return zi, zm
