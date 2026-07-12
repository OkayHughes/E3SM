"""Port of the orographic gravity-wave spine of
eam/src/physics/cam/gw/: gw_prof + gw_oro_src (gw_oro.F90) +
gw_drag_prof with ngwv=0 (gw_common.F90).

PORT_NOTES
----------
- Index convention: Fortran interfaces run 0:pver and midpoints 1:pver.
  Here interface arrays have pver+1 entries (same interface index k);
  Fortran midpoint k maps to python k-1.
- ngwv=0 restricts to the single c=0 orographic wave; the spectrum
  branch (gwd_project_tau, gwd_precalc_rhoi with the LU diffusion
  solver) is NOT ported yet — see PORTING_PLAN.md.
- `orographic_only` changes three things, all ported: tndmax (500 vs
  400 m/s/day), the extra WKB saturation limiter on the tendency, and
  whether efficiency/taper is applied inside the wave loop (spectrum
  style) or on the projected tendency (vanilla-CAM answer
  preservation).
- The stress build (bottom-up) and tendency sweep (top-down, with the
  stress-rewrite feedback tau[k] = tau[k-1] + ubtl*dpm/g) are
  order-critical sequential loops, ported as unrolled masked updates.
- gw_oro_src's source-region average accumulates layers while
  2*sgh > sqrt(zm[k]*zm[k+1]) scanning upward; ported with per-level
  masks (identical result, including non-contiguity if it ever arose).
- EAMv3 defaults (namelist/gw_drag.F90): kwv=6.28e-5, fcrit2=1.0
  (namelist, must be set), effgw_oro=0.375, ktop=0, and alpha is the
  Newtonian cooling profile from the gw_drag file (input here).

Parameters live in a dict from make_gw_params(); alpha is (pver+1,).
"""

import jax.numpy as jnp

DBACK = 0.05        # background diffusivity
TAUMIN = 1.0e-10
UMCFAC = 0.5
UBMC2MN = 0.01
N2MIN = 1.0e-8
OROHMIN = 10.0      # min surface displacement for oro waves
OROVMIN = 2.0       # min wind speed for oro waves


def make_gw_params(alpha, kbotbg, pver, fcrit2=1.0, kwv=6.28e-5,
                   gravit=9.80616, rair=287.042, ktop=0,
                   orographic_only=True, tau_0_ubc=False):
    tndmax = 500.0 / 86400.0 if orographic_only else 400.0 / 86400.0
    return dict(alpha=jnp.asarray(alpha, dtype=jnp.float64),
                kbotbg=int(kbotbg), pver=int(pver), fcrit2=fcrit2,
                kwv=kwv, effkwv=kwv * fcrit2, gravit=gravit, rair=rair,
                rog=rair / gravit, ktop=int(ktop),
                orographic_only=bool(orographic_only),
                tau_0_ubc=bool(tau_0_ubc), tndmax=tndmax,
                oroko2=0.5 * kwv)


def midpoint_interp(a):
    return 0.5 * (a[:, :-1] + a[:, 1:])


def get_unit_vector(u, v):
    mag = jnp.sqrt(u * u + v * v)
    ok = mag > 0.0
    return (jnp.where(ok, u / jnp.where(ok, mag, 1.0), 0.0),
            jnp.where(ok, v / jnp.where(ok, mag, 1.0), 0.0), mag)


def _sign(a, b):
    """Fortran sign(a, b): |a| with the sign of b (+ for b == 0)."""
    return jnp.where(b >= 0.0, jnp.abs(a), -jnp.abs(a))


def gw_prof(cpair, t, pmid, pint, params):
    """Interface densities/temperatures and Brunt-Vaisala frequencies.
    t/pmid: (ncol, pver); pint: (ncol, pver+1).
    Returns (rhoi, ti, nm, ni) — interface arrays (ncol, pver+1)."""
    rair, gravit = params["rair"], params["gravit"]
    ti_top = t[:, :1]
    ti_int = midpoint_interp(t)
    ti_bot = t[:, -1:]
    ti = jnp.concatenate([ti_top, ti_int, ti_bot], axis=1)
    rhoi = pint / (rair * ti)

    ni_top = jnp.sqrt(gravit * gravit / (cpair * ti_top))
    dtdp = (t[:, 1:] - t[:, :-1]) / (pmid[:, 1:] - pmid[:, :-1])
    n2 = gravit * gravit / ti_int * (1.0 / cpair - rhoi[:, 1:-1] * dtdp)
    ni_int = jnp.sqrt(jnp.maximum(N2MIN, n2))
    ni = jnp.concatenate([ni_top, ni_int, ni_int[:, -1:]], axis=1)
    nm = midpoint_interp(ni)
    return rhoi, ti, nm, ni


def gw_oro_src(u, v, t, sgh, pmid, pint, dpm, zm, nm, params):
    """Orographic (c=0) wave source, McFarlane 1987.
    Returns (src_level, tend_level, tau, ubm, ubi, xv, yv); tau is
    (ncol, pver+1) for the single wave; *_level are interface indices."""
    pver = params["pver"]
    rair = params["rair"]
    ncol = u.shape[0]

    hdsp = 2.0 * sgh

    # source-region dp-weighted averages, scanning up from the bottom
    src_level = jnp.full(ncol, pver - 1, dtype=jnp.int64)
    rsrc = pmid[:, pver - 1] / (rair * t[:, pver - 1]) * dpm[:, pver - 1]
    usrc = u[:, pver - 1] * dpm[:, pver - 1]
    vsrc = v[:, pver - 1] * dpm[:, pver - 1]
    nsrc = nm[:, pver - 1] * dpm[:, pver - 1]
    for kf in range(pver - 1, pver // 2 - 1, -1):   # Fortran k
        kk = kf - 1                                  # python midpoint
        m = hdsp > jnp.sqrt(zm[:, kk] * zm[:, kk + 1])
        src_level = jnp.where(m, kf - 1, src_level)
        rsrc = rsrc + jnp.where(
            m, pmid[:, kk] / (rair * t[:, kk]) * dpm[:, kk], 0.0)
        usrc = usrc + jnp.where(m, u[:, kk] * dpm[:, kk], 0.0)
        vsrc = vsrc + jnp.where(m, v[:, kk] * dpm[:, kk], 0.0)
        nsrc = nsrc + jnp.where(m, nm[:, kk] * dpm[:, kk], 0.0)

    dpsrc = pint[:, pver] - jnp.take_along_axis(
        pint, src_level[:, None], axis=1)[:, 0]
    rsrc, usrc = rsrc / dpsrc, usrc / dpsrc
    vsrc, nsrc = vsrc / dpsrc, nsrc / dpsrc

    xv, yv, ubi_sfc = get_unit_vector(usrc, vsrc)
    ubm = u * xv[:, None] + v * yv[:, None]
    ubi = jnp.concatenate(
        [ubm[:, :1], midpoint_interp(ubm), ubi_sfc[:, None]], axis=1)

    on = (ubi_sfc > OROVMIN) & (hdsp > OROHMIN)
    sghmax = params["fcrit2"] * (ubi_sfc / nsrc) ** 2
    tauoro = jnp.where(
        on, params["oroko2"] * jnp.minimum(hdsp ** 2, sghmax)
        * rsrc * nsrc * ubi_sfc, 0.0)
    src_level = jnp.where(on, src_level, pver)

    ki = jnp.arange(pver + 1)
    tau = jnp.where((ki[None, :] >= src_level[:, None])
                    & (ki[None, :] <= pver), tauoro[:, None], 0.0)

    tend_level = jnp.full(ncol, pver, dtype=jnp.int64)
    return src_level, tend_level, tau, ubm, ubi, xv, yv


def gw_drag_prof_oro(src_level, tend_level, dt, lat, t, ti, pmid, pint,
                     dpm, rdpm, piln, rhoi, nm, ni, ubm, ubi, xv, yv,
                     effgw, tau, params, do_taper=False):
    """gw_drag_prof for the single c=0 orographic wave (ngwv=0).
    Returns (tau, utgw, vtgw). All sequential level sweeps are ported
    with their exact update order and masks."""
    pver, ktop = params["pver"], params["ktop"]
    alpha = params["alpha"]
    effkwv, kwv = params["effkwv"], params["kwv"]
    rog = params["rog"]
    gravit = params["gravit"]
    oro_only = params["orographic_only"]

    tau = jnp.asarray(tau, dtype=jnp.float64)
    if params["tau_0_ubc"]:
        tau = tau.at[:, 0].set(0.0)

    # ---- stress profile, bottom-up (c = 0 so ubmc = ubi) ----
    for k in range(pver - 1, ktop - 1, -1):
        m = src_level > k
        ubmc = ubi[:, k]
        samesign = ubmc * ubi[:, k + 1] > 0.0
        tausat = jnp.where(
            samesign,
            jnp.abs(effkwv * rhoi[:, k] * ubmc ** 3 / (2.0 * ni[:, k])),
            0.0)
        tausat = jnp.where(tausat <= TAUMIN, 0.0, tausat)

        dsat = (ubmc / ni[:, k]) ** 2 * (
            effkwv * ubmc ** 2 / (rog * ti[:, k] * ni[:, k]) - alpha[k])
        dscal = jnp.minimum(1.0, tau[:, k + 1] / (tausat + TAUMIN))
        dd = jnp.maximum(DBACK, dscal * dsat)

        ubmc2 = jnp.maximum(ubmc ** 2, UBMC2MN)
        mi = ni[:, k] / (2.0 * kwv * ubmc2) * (
            alpha[k] + ni[:, k] ** 2 / ubmc2 * dd)
        wrk = -2.0 * mi * rog * t[:, k] * (piln[:, k + 1] - piln[:, k])
        taudmp = tau[:, k + 1] * jnp.exp(wrk)
        taudmp = jnp.where(taudmp <= TAUMIN, 0.0, taudmp)
        tau = tau.at[:, k].set(
            jnp.where(m, jnp.minimum(taudmp, tausat), tau[:, k]))

    # ---- tendencies from stress divergence, top-down ----
    ptaper = jnp.cos(lat) if do_taper else jnp.ones_like(lat)
    utgw = jnp.zeros_like(ubm)
    vtgw = jnp.zeros_like(ubm)

    for k in range(ktop + 1, pver + 1):     # Fortran interface/mid k
        km = k - 1                           # python midpoint
        m = k <= tend_level
        ubtl = gravit * (tau[:, k] - tau[:, k - 1]) * rdpm[:, km]
        if oro_only:
            ubtlsat = effkwv * jnp.abs((-ubm[:, km]) ** 3) / (
                2.0 * rog * t[:, km] * nm[:, km])
            ubtl = jnp.minimum(ubtl, ubtlsat)
        ubtl = jnp.minimum(ubtl, UMCFAC * jnp.abs(-ubm[:, km]) / dt)
        ubtl = jnp.minimum(ubtl, params["tndmax"])

        gwut_k = _sign(ubtl, -ubm[:, km]) * effgw * ptaper
        if oro_only:
            ubt_k = _sign(ubtl, -ubm[:, km])
        else:
            ubt_k = gwut_k
        ubt_k = jnp.where(m, ubt_k, 0.0)

        tau = tau.at[:, k].set(
            jnp.where(m, tau[:, k - 1] + ubtl * dpm[:, km] / gravit,
                      tau[:, k]))

        scale = effgw * ptaper if oro_only else 1.0
        utgw = utgw.at[:, km].set(
            jnp.where(m, ubt_k * xv * scale, utgw[:, km]))
        vtgw = vtgw.at[:, km].set(
            jnp.where(m, ubt_k * yv * scale, vtgw[:, km]))

    return tau, utgw, vtgw
