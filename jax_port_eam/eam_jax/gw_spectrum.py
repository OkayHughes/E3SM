"""Port of the full-spectrum (non-orographic) branch of the EAMv3
gravity-wave suite, eam/src/physics/cam/gw/: gw_convect.F90
(gw_beres_src, Beres 2004 convective source), gw_front.F90 (gw_cm_src,
Charron & Manzini frontal source), the ngwv > 0 paths of
gw_common.F90 gw_drag_prof (gwd_compute_stress_profiles_and_
diffusivities, gwd_project_tau, gwd_compute_tendencies_from_stress_
divergence, gwd_precalc_rhoi) with gw_diffusion.F90 (gw_ediff +
gw_diff_tend) and cam/vdiff_lu_solver.F90 (vd_lu_decomp/vd_lu_solve),
and momentum_energy_conservation. Extends eam_jax.gw (orographic
spine); gw_prof and the shared helpers are imported from there.

PORT_NOTES
----------
- Index conventions follow eam_jax.gw: interface arrays have pver+1
  entries with the same interface index k as the Fortran 0:pver;
  Fortran midpoint k maps to python k-1. src_level/tend_level keep
  their Fortran interface values. The wave axis is 0-based: Fortran
  wavenumber l in -pgwv:pgwv maps to python l+pgwv; tau is
  (ncol, nwav, pver+1), c is (ncol, nwav), cref[l+pgwv] = dc*l.
  taucd is (ncol, pver+1, 4) with the last axis in Fortran order
  (west, east, south, north) = indices (0, 1, 2, 3).
- do_molec_diff is .false. in EAMv3 default physics (no WACCM);
  the kvtt/nbot_molec branches of the stress loop collapse to the
  non-molecular path and are ported that way (kvtt is not an input).
- Sequential level sweeps (stress bottom-up, tendency top-down with
  the tau rewrite, LU forward/backward recursions, the conservation
  dz/dE integrals, heating-depth scans) are python loops in the exact
  Fortran order. Reductions whose order affects the last bits (the
  ubt wave sum, taub/tauf in project_tau, egwdffm/dttke wave sums,
  the dE level sum) accumulate term-by-term in Fortran loop order;
  order-independent reductions (the max over waves in the diffusivity
  d) are vectorized.
- gw_convect_init's k_src_wind and gw_front_init's fav are computed
  in make_beres_params/make_front_params from the same inputs
  (pref_edge scan; trapezoid average of the gaussian over each bin,
  fav(0) forced to 0). The tau assignment in gw_cm_src uses fav(|l|)
  for both signs, so the port indexes fav at |l| rather than assuming
  numerical symmetry of the init loop.
- The Beres lookup table mfcc is a plain input array (production
  reads it from the Beres04 netCDF file); index mapping
  mfcc[h-1, uh+maxuh, l+pgwv]. Fortran NINT (round half away from
  zero) and int() (truncate) are ported explicitly; cshift's circular
  left shift becomes a modular gather.
- gw_beres_src/gw_cm_src leave ubm/ubi (and ubm above kbot for the
  frontal source) unset outside the projection range; the harness
  driver zero-initializes them, and the port writes zeros in the same
  slots. ubi[:, k_src_wind] is first set to the source-wind magnitude
  and then overwritten by the midpoint interpolation, exactly as the
  Fortran assignment order does.
- gw_heating_depth: both the use_gw_convect_old=.true. (EAMv3
  default; CAM4-6/EAMv1-3 behavior, melting-level bug included) and
  the corrected variant are ported. Columns are assumed to have
  zm >= 20 km somewhere (true for any EAM grid), so mini/maxi are
  always found.
- vd_lu_decomp/vd_lu_solve are ported exactly, with ntop/nbot as the
  Fortran 1-based midpoint indices; the decomp is a dict of
  (ca, cc, dnom, ze), each (ncol, pver) with Fortran column k at
  python k-1.
- momentum_energy_conservation divides by the sub-source-level mass;
  columns with tend_level == pver would divide by zero exactly as the
  Fortran does (not reachable from these sources on real grids).
- EAMv3 default spectrum settings (bld/build-namelist +
  namelist_defaults_eam.xml phys="default"): pgwv=32, dc=2.5,
  fcrit2=1.0, kwv=6.28e-5, taubgnd=2.5e-3, frontgfc=1.25e-15,
  effgw_cm=1.0, effgw_beres=0.35, gw_convect_hcf=10.0,
  hdepth_scaling_factor=0.5, hdepth_min=2.5, storm_speed_min=10.0,
  plev_src_wind=70000 Pa, use_gw_convect_old=.true.,
  tau_0_ubc=.false., ktop=0, kbotbg/kfront from pref_edge at 500/600
  hPa, do_taper=.false. (SE dycore).
"""

import numpy as np
import jax.numpy as jnp

from .gw import (DBACK, TAUMIN, UMCFAC, UBMC2MN, get_unit_vector,
                 make_gw_params, midpoint_interp, _sign)
from .gw import gw_prof  # noqa: F401  (re-exported for chained use)

# gw_convect.F90 fixed parameters
HEATING_ALTITUDE_MAX = np.float32(20e3)   # single-precision literal
TAU_AVG_LENGTH = np.float32(100e3)        # single-precision literal
# gw_diffusion.F90
PRNDL = 0.25
EGWDFFI_MAX = 150.0
# gw_front.F90 fav integration
FAV_DCA = 0.1
FAV_C0 = 30.0


def _nint(x):
    """Fortran NINT: round to nearest, halves away from zero."""
    return jnp.where(x >= 0.0, jnp.floor(x + 0.5),
                     jnp.ceil(x - 0.5)).astype(jnp.int64)


def make_gw_spectrum_params(alpha, pgwv, dc, kbotbg, pver, fcrit2=1.0,
                            kwv=6.28e-5, gravit=9.80616, rair=287.042,
                            ktop=0, tau_0_ubc=False):
    """gw_common_init for the spectrum configuration (pgwv > 0,
    orographic_only=.false., do_molec_diff=.false.)."""
    p = make_gw_params(alpha, kbotbg, pver, fcrit2=fcrit2, kwv=kwv,
                       gravit=gravit, rair=rair, ktop=ktop,
                       orographic_only=False, tau_0_ubc=tau_0_ubc)
    p["pgwv"] = int(pgwv)
    p["dc"] = float(dc)
    p["nwav"] = 2 * int(pgwv) + 1
    p["cref"] = jnp.asarray(
        dc * np.arange(-pgwv, pgwv + 1, dtype=np.float64))
    return p


def make_front_params(params, taubgnd=2.5e-3, frontgfc=1.25e-15,
                      kfront=1):
    """gw_front_init: frontogenesis threshold + fav (bin-averaged
    gaussian source spectrum scaled by taubgnd, wavenumber 0
    prohibited). kfront is the Fortran 1-based midpoint index of the
    600 hPa check level."""
    pgwv, dc = params["pgwv"], params["dc"]
    cref = np.asarray(params["cref"])
    fav = np.zeros(2 * pgwv + 1)
    if pgwv > 0:
        n_sub = int(np.rint(dc / FAV_DCA))
        for li, l in enumerate(range(-pgwv, pgwv + 1)):
            cmn = cref[li] - 0.5 * dc
            cmx = cref[li] + 0.5 * dc
            f = 0.5 * FAV_DCA * (np.exp(-(cmn / FAV_C0) ** 2)
                                 + np.exp(-(cmx / FAV_C0) ** 2))
            for n in range(1, n_sub):
                f += FAV_DCA * np.exp(-((cmn + n * FAV_DCA) / FAV_C0)
                                      ** 2)
            fav[li] = f
        fav = fav / dc
    else:
        fav[0] = 1.0
    fav = taubgnd * fav
    fav[pgwv] = 0.0
    return dict(taubgnd=float(taubgnd), frontgfc=float(frontgfc),
                kfront=int(kfront), fav=jnp.asarray(fav))


def make_beres_params(params, mfcc, pref_edge, plev_src_wind=70000.0,
                      maxq0_conversion_factor=10.0,
                      hdepth_scaling_factor=0.5, hdepth_min=2.5,
                      storm_speed_min=10.0, use_gw_convect_old=True):
    """gw_convect_init: steering-flow level from pref_edge plus the
    Beres source-spectrum table (mfcc[h-1, uh+maxuh, l+pgwv])."""
    mfcc = jnp.asarray(mfcc, dtype=jnp.float64)
    maxh = mfcc.shape[0]
    maxuh = (mfcc.shape[1] - 1) // 2
    pver = params["pver"]
    pref_edge = np.asarray(pref_edge)
    k_src_wind = 0
    for k in range(0, pver + 1):
        if pref_edge[k] < plev_src_wind:   # Fortran pref_edge(k+1)
            k_src_wind = k + 1
    return dict(mfcc=mfcc, maxh=maxh, maxuh=maxuh,
                k_src_wind=int(k_src_wind),
                maxq0_conversion_factor=float(maxq0_conversion_factor),
                hdepth_scaling_factor=float(hdepth_scaling_factor),
                hdepth_min=float(hdepth_min),
                storm_speed_min=float(storm_speed_min),
                use_gw_convect_old=bool(use_gw_convect_old))


# ---------------------------------------------------------------------------
# gw_convect.F90: Beres deep-convection source
# ---------------------------------------------------------------------------
def _convect_project_winds(u, v, params, bp):
    pver = params["pver"]
    ksw = bp["k_src_wind"]                 # Fortran 1-based midpoint
    xv, yv, mag = get_unit_vector(u[:, ksw - 1], v[:, ksw - 1])
    ubm = u * xv[:, None] + v * yv[:, None]
    ubi = jnp.zeros((u.shape[0], pver + 1), dtype=jnp.float64)
    ubi = ubi.at[:, ksw].set(mag)          # overwritten below, as in F
    ubi = ubi.at[:, 0].set(ubm[:, 0])
    ubi = ubi.at[:, 1:pver].set(midpoint_interp(ubm))
    return xv, yv, ubm, ubi


def _heating_depth(zm, netdt, params, bp):
    pver = params["pver"]
    ncol = zm.shape[0]
    mini = jnp.zeros(ncol, dtype=jnp.int64)
    maxi = jnp.zeros(ncol, dtype=jnp.int64)
    hmax = jnp.float64(HEATING_ALTITUDE_MAX)
    if bp["use_gw_convect_old"]:
        # original CAM4-6/EAMv1-3 scan (melting-level bug preserved)
        for k in range(pver, 0, -1):
            above = zm[:, k - 1] >= hmax
            pos = netdt[:, k - 1] > 0.0
            m0 = mini == 0
            x0 = (~m0) & (maxi == 0)
            maxi = jnp.where(m0 & above, k,
                             jnp.where(x0 & (above | ~pos), k, maxi))
            mini = jnp.where(m0 & (above | pos), k, mini)
    else:
        for k in range(pver, 0, -1):
            below = zm[:, k - 1] < hmax
            pos = netdt[:, k - 1] > 0.0
            mini = jnp.where(below & pos & (mini == 0), k,
                             jnp.where(~below & (mini == 0), k, mini))
            maxi = jnp.where(below & pos, k,
                             jnp.where(~below & (maxi == 0), k, maxi))

    zmaxi = jnp.take_along_axis(zm, maxi[:, None] - 1, axis=1)[:, 0]
    zmini = jnp.take_along_axis(zm, mini[:, None] - 1, axis=1)[:, 0]
    hdepth = (zmaxi - zmini) / 1000.0
    hdepth = jnp.minimum(hdepth, float(bp["maxh"]))
    hdepth = hdepth * bp["hdepth_scaling_factor"]

    maxq0 = jnp.zeros(ncol, dtype=jnp.float64)
    for k in range(int(maxi.min()), int(mini.max()) + 1):
        m = (k >= maxi) & (k <= mini)
        maxq0 = jnp.where(m, jnp.maximum(maxq0, netdt[:, k - 1]), maxq0)
    maxq0_out = maxq0 * 24.0 * 3600.0
    maxq0 = maxq0 * bp["maxq0_conversion_factor"]
    return mini, maxi, hdepth, maxq0_out, maxq0


def _storm_speed(ubm, mini, maxi, params, bp):
    pgwv, dc = params["pgwv"], params["dc"]
    maxuh = bp["maxuh"]
    ub_src = ubm[:, bp["k_src_wind"] - 1]
    storm_speed = jnp.trunc(_sign(
        jnp.maximum(jnp.abs(ub_src) - bp["storm_speed_min"], 0.0),
        ub_src)).astype(jnp.int64)

    uh = jnp.zeros(ubm.shape[0], dtype=jnp.float64)
    for k in range(int(maxi.min()), int(mini.max()) + 1):
        m = (k >= maxi) & (k <= mini)
        uh = uh + jnp.where(m, ubm[:, k - 1] / (mini - maxi + 1), 0.0)
    uh = uh - storm_speed
    uh = jnp.minimum(uh, float(maxuh))
    uh = jnp.maximum(uh, -float(maxuh))

    umin = jnp.full(ubm.shape[0], pgwv * dc, dtype=jnp.float64)
    umax = jnp.full(ubm.shape[0], -pgwv * dc, dtype=jnp.float64)
    for k in range(int(maxi.min()), int(mini.max()) + 1):
        m = (k >= maxi) & (k <= mini)
        umin = jnp.where(m, jnp.minimum(umin, ubm[:, k - 1]), umin)
        umax = jnp.where(m, jnp.maximum(umax, ubm[:, k - 1]), umax)
    return storm_speed, uh, umin, umax


def gw_beres_src(lat, u, v, netdt, zm, params, bp):
    """Beres 2004 convective gravity-wave source (gw_beres_src).
    Returns (src_level, tend_level, tau, ubm, ubi, xv, yv, c, hdepth,
    maxq0_out); tau is (ncol, nwav, pver+1)."""
    pver, pgwv, dc = params["pver"], params["pgwv"], params["dc"]
    nwav = params["nwav"]
    maxuh = bp["maxuh"]
    ncol = u.shape[0]

    xv, yv, ubm, ubi = _convect_project_winds(u, v, params, bp)
    mini, maxi, hdepth, maxq0_out, maxq0 = _heating_depth(
        zm, netdt, params, bp)
    storm_speed, uh, umin, umax = _storm_speed(ubm, mini, maxi, params,
                                               bp)

    launch = (hdepth >= bp["hdepth_min"]) & (jnp.abs(lat) < np.pi / 2)
    hg = jnp.clip(_nint(hdepth), 1, bp["maxh"]) - 1
    ug = jnp.clip(_nint(uh), -maxuh, maxuh) + maxuh
    tau0 = bp["mfcc"][hg, ug, :]                     # (ncol, nwav)
    # cshift(tau0, shift): circular left shift, relative to the ground
    shift = -_nint(storm_speed / dc)
    idx = (jnp.arange(nwav)[None, :] + shift[:, None]) % nwav
    tau0 = jnp.take_along_axis(tau0, idx, axis=1)
    tau0 = tau0 * (maxq0 * maxq0 / jnp.float64(TAU_AVG_LENGTH))[:, None]
    # critical-level filtering
    umini = jnp.maximum(_nint(umin / dc), -pgwv)
    umaxi = jnp.minimum(_nint(umax / dc), pgwv)
    li = jnp.arange(-pgwv, pgwv + 1)
    filt = ((umaxi > umini)[:, None] & (li[None, :] >= umini[:, None])
            & (li[None, :] <= umaxi[:, None]))
    tau0 = jnp.where(filt, 0.0, tau0)

    ki = jnp.arange(pver + 1)
    tau = jnp.where(
        launch[:, None, None]
        & (ki[None, None, :] == maxi[:, None, None]),
        tau0[:, :, None], 0.0)

    c = jnp.broadcast_to(params["cref"][None, :], (ncol, nwav)).copy()
    return (maxi, maxi, tau, ubm, ubi, xv, yv, c, hdepth, maxq0_out)


# ---------------------------------------------------------------------------
# gw_front.F90: frontogenesis source
# ---------------------------------------------------------------------------
def gw_cm_src(u, v, frontgf, params, fp):
    """Charron & Manzini frontal gravity-wave source (gw_cm_src) with
    kbot = kbotbg. Returns (src_level, tend_level, tau, ubm, ubi, xv,
    yv, c)."""
    pver, pgwv, nwav = params["pver"], params["pgwv"], params["nwav"]
    kbot = params["kbotbg"]
    ncol = u.shape[0]

    usrc = 0.5 * (u[:, kbot] + u[:, kbot - 1])   # F: u(:,kbot+1)+u(:,kbot)
    vsrc = 0.5 * (v[:, kbot] + v[:, kbot - 1])
    xv, yv, mag = get_unit_vector(usrc, vsrc)
    ubm = jnp.zeros((ncol, pver), dtype=jnp.float64)
    ubm = ubm.at[:, :kbot].set(u[:, :kbot] * xv[:, None]
                               + v[:, :kbot] * yv[:, None])
    ubi = jnp.zeros((ncol, pver + 1), dtype=jnp.float64)
    ubi = ubi.at[:, kbot].set(mag)
    ubi = ubi.at[:, 0].set(ubm[:, 0])
    ubi = ubi.at[:, 1:kbot].set(midpoint_interp(ubm[:, :kbot]))

    launch = frontgf[:, fp["kfront"] - 1] > fp["frontgfc"]
    lsym = jnp.abs(jnp.arange(-pgwv, pgwv + 1))
    fav_used = fp["fav"][pgwv + lsym]            # fav(|l|), fav(0)=0
    tau = jnp.zeros((ncol, nwav, pver + 1), dtype=jnp.float64)
    tau = tau.at[:, :, kbot].set(
        jnp.where(launch[:, None], fav_used[None, :], 0.0))

    src_level = jnp.full(ncol, kbot, dtype=jnp.int64)
    c = params["cref"][None, :] + jnp.abs(ubi[:, kbot])[:, None]
    return src_level, src_level, tau, ubm, ubi, xv, yv, c


# ---------------------------------------------------------------------------
# gw_common.F90: full-spectrum gw_drag_prof
# ---------------------------------------------------------------------------
def _stress_profiles(src_level, ubi, c, rhoi, ni, t, ti, piln, tau,
                     params):
    """gwd_compute_stress_profiles_and_diffusivities
    (do_molec_diff=.false. path)."""
    ktop = params["ktop"]
    alpha = params["alpha"]
    effkwv, kwv, rog = params["effkwv"], params["kwv"], params["rog"]

    for k in range(int(src_level.max()) - 1, ktop - 1, -1):
        m = src_level > k
        ubmc = ubi[:, k:k + 1] - c                       # (ncol, nwav)
        samesign = ubmc * (ubi[:, k + 1:k + 2] - c) > 0.0
        tausat = jnp.where(
            samesign,
            jnp.abs(effkwv * rhoi[:, k:k + 1] * ubmc ** 3
                    / (2.0 * ni[:, k:k + 1])), 0.0)
        tausat = jnp.where(tausat <= TAUMIN, 0.0, tausat)

        dsat = (ubmc / ni[:, k:k + 1]) ** 2 * (
            effkwv * ubmc ** 2 / (rog * ti[:, k:k + 1] * ni[:, k:k + 1])
            - alpha[k])
        dscal = jnp.minimum(1.0, tau[:, :, k + 1] / (tausat + TAUMIN))
        dd = jnp.maximum(DBACK, jnp.max(dscal * dsat, axis=1))

        ubmc2 = jnp.maximum(ubmc ** 2, UBMC2MN)
        mi = ni[:, k:k + 1] / (2.0 * kwv * ubmc2) * (
            alpha[k] + ni[:, k:k + 1] ** 2 / ubmc2 * dd[:, None])
        wrk = -2.0 * mi * rog * t[:, k:k + 1] * (
            piln[:, k + 1:k + 2] - piln[:, k:k + 1])
        taudmp = tau[:, :, k + 1] * jnp.exp(wrk)
        taudmp = jnp.where(taudmp <= TAUMIN, 0.0, taudmp)
        tau = tau.at[:, :, k].set(
            jnp.where(m[:, None], jnp.minimum(taudmp, tausat),
                      tau[:, :, k]))
    return tau


def _project_tau(tend_level, tau, ubi, c, xv, yv, params):
    """gwd_project_tau: tau in the four cardinal directions."""
    pver, ktop, nwav = params["pver"], params["ktop"], params["nwav"]
    ncol = tau.shape[0]
    ki = jnp.arange(pver + 1)
    mask_k = ((ki[None, :] <= tend_level[:, None])
              & (ki[None, :] >= ktop))
    ubi_tend = jnp.take_along_axis(ubi, tend_level[:, None],
                                   axis=1)[:, 0]

    taub = jnp.zeros((ncol, pver + 1), dtype=jnp.float64)
    tauf = jnp.zeros((ncol, pver + 1), dtype=jnp.float64)
    for l in range(nwav):                    # Fortran accumulation order
        tausg = _sign(tau[:, l, :], (c[:, l:l + 1] - ubi))
        behind = (c[:, l] < ubi_tend)[:, None]
        forward = (c[:, l] > ubi_tend)[:, None]
        taub = taub + jnp.where(mask_k & behind, tausg, 0.0)
        tauf = tauf + jnp.where(mask_k & forward, tausg, 0.0)

    xvp = mask_k & (xv[:, None] > 0.0)
    xvn = mask_k & (xv[:, None] < 0.0)
    yvp = mask_k & (yv[:, None] > 0.0)
    yvn = mask_k & (yv[:, None] < 0.0)
    east = jnp.where(xvp, tauf * xv[:, None],
                     jnp.where(xvn, taub * xv[:, None], 0.0))
    west = jnp.where(xvp, taub * xv[:, None],
                     jnp.where(xvn, tauf * xv[:, None], 0.0))
    north = jnp.where(yvp, tauf * yv[:, None],
                      jnp.where(yvn, taub * yv[:, None], 0.0))
    south = jnp.where(yvp, taub * yv[:, None],
                      jnp.where(yvn, tauf * yv[:, None], 0.0))
    return jnp.stack([west, east, south, north], axis=-1)


def _tendencies_from_stress_divergence(tend_level, do_taper, dt, effgw,
                                       lat, dpm, rdpm, c, ubm, t, nm,
                                       xv, yv, tau, params):
    """gwd_compute_tendencies_from_stress_divergence
    (orographic_only=.false. path)."""
    pver, ktop, nwav = params["pver"], params["ktop"], params["nwav"]
    gravit, tndmax = params["gravit"], params["tndmax"]
    ncol = ubm.shape[0]

    ptaper = jnp.cos(lat) if do_taper else jnp.ones_like(lat)
    if params["tau_0_ubc"]:
        tau = tau.at[:, :, 0].set(0.0)

    gwut = jnp.zeros((ncol, pver, nwav), dtype=jnp.float64)
    utgw = jnp.zeros((ncol, pver), dtype=jnp.float64)
    vtgw = jnp.zeros((ncol, pver), dtype=jnp.float64)

    for k in range(ktop + 1, int(tend_level.max()) + 1):
        km = k - 1                            # python midpoint
        mk = k <= tend_level
        cmu = c - ubm[:, km:km + 1]           # (ncol, nwav)
        ubtl = gravit * (tau[:, :, k] - tau[:, :, k - 1]) \
            * rdpm[:, km:km + 1]
        ubtl = jnp.minimum(ubtl, UMCFAC * jnp.abs(cmu) / dt)
        ubtl = jnp.minimum(ubtl, tndmax)

        gwut_k = _sign(ubtl, cmu) * effgw * ptaper[:, None]
        gwut = gwut.at[:, km, :].set(
            jnp.where(mk[:, None], gwut_k, gwut[:, km, :]))

        ubt = jnp.zeros(ncol, dtype=jnp.float64)
        for l in range(nwav):                 # Fortran accumulation order
            ubt = ubt + jnp.where(mk, gwut_k[:, l], 0.0)

        tau = tau.at[:, :, k].set(
            jnp.where(mk[:, None],
                      tau[:, :, k - 1] + ubtl * dpm[:, km:km + 1]
                      / gravit, tau[:, :, k]))
        utgw = utgw.at[:, km].set(jnp.where(mk, ubt * xv, utgw[:, km]))
        vtgw = vtgw.at[:, km].set(jnp.where(mk, ubt * yv, vtgw[:, km]))
    return tau, gwut, utgw, vtgw


def vd_lu_decomp(ksrf, kv, tmpi, rpdel, ztodt, gravit, cc_top, ntop,
                 nbot):
    """vdiff_lu_solver.F90 vd_lu_decomp. ntop/nbot are the Fortran
    1-based midpoint bounds; kv/tmpi are (ncol, pver+1), Fortran
    column j at python j-1. Returns dict(ca, cc, dnom, ze)."""
    ncol, pver = rpdel.shape
    ca = jnp.zeros((ncol, pver), dtype=jnp.float64)
    cc = jnp.zeros((ncol, pver), dtype=jnp.float64)
    dnom = jnp.zeros((ncol, pver), dtype=jnp.float64)
    ze = jnp.zeros((ncol, pver), dtype=jnp.float64)

    kvt = kv[:, ntop:nbot] * tmpi[:, ntop:nbot]  # F kv(:,k+1)*tmpi(:,k+1)
    ca = ca.at[:, ntop - 1:nbot - 1].set(kvt * rpdel[:, ntop - 1:nbot - 1])
    cc = cc.at[:, ntop:nbot].set(kvt * rpdel[:, ntop:nbot])
    # ca(:,nbot) = 0 -- already zero.

    dnom = dnom.at[:, nbot - 1].set(1.0 / (
        1.0 + cc[:, nbot - 1]
        + ksrf * ztodt * gravit * rpdel[:, nbot - 1]))
    ze = ze.at[:, nbot - 1].set(cc[:, nbot - 1] * dnom[:, nbot - 1])
    for k in range(nbot - 1, ntop, -1):          # F k = nbot-1 .. ntop+1
        j = k - 1
        dnom = dnom.at[:, j].set(1.0 / (
            1.0 + ca[:, j] + cc[:, j] - ca[:, j] * ze[:, j + 1]))
        ze = ze.at[:, j].set(cc[:, j] * dnom[:, j])
    dnom = dnom.at[:, ntop - 1].set(1.0 / (
        1.0 + ca[:, ntop - 1] + cc_top
        - ca[:, ntop - 1] * ze[:, ntop]))
    return dict(ca=ca, cc=cc, dnom=dnom, ze=ze)


def vd_lu_solve(q, decomp, ntop, nbot, cd_top):
    """vdiff_lu_solver.F90 vd_lu_solve (Richtmyer & Morton back
    substitution); q is (ncol, pver), returns the solved q."""
    ca, ze, dnom = decomp["ca"], decomp["ze"], decomp["dnom"]
    q = jnp.asarray(q, dtype=jnp.float64)
    zf = jnp.zeros_like(q)
    zf = zf.at[:, nbot - 1].set(q[:, nbot - 1] * dnom[:, nbot - 1])
    for k in range(nbot - 1, ntop, -1):          # F k = nbot-1 .. ntop+1
        j = k - 1
        zf = zf.at[:, j].set((q[:, j] + ca[:, j] * zf[:, j + 1])
                             * dnom[:, j])
    zf = zf.at[:, ntop - 1].set(
        (q[:, ntop - 1] + cd_top + ca[:, ntop - 1] * zf[:, ntop])
        * dnom[:, ntop - 1])
    q = q.at[:, ntop - 1].set(zf[:, ntop - 1])
    for k in range(ntop + 1, nbot + 1):          # F k = ntop+1 .. nbot
        j = k - 1
        q = q.at[:, j].set(zf[:, j] + ze[:, j] * q[:, j - 1])
    return q


def _gw_ediff(tend_level, gwut, ubm, nm, rho, dt, pmid, rdpm, c,
              params):
    """gw_diffusion.F90 gw_ediff with kbot=kbotbg, ktop; rho is the
    1-based rhoi_kludge (ncol, pver+1). Returns (egwdffi, decomp)."""
    pver, ktop, kbot = params["pver"], params["ktop"], params["kbotbg"]
    nwav, gravit = params["nwav"], params["gravit"]
    ncol = ubm.shape[0]

    egwdffm = jnp.zeros((ncol, pver), dtype=jnp.float64)
    sl = slice(ktop, kbot)                    # F midpoints ktop+1..kbot
    for l in range(nwav):                     # Fortran accumulation order
        egwdffm = egwdffm.at[:, sl].add(
            PRNDL * 0.5 * gwut[:, sl, l] * (c[:, l:l + 1] - ubm[:, sl])
            / nm[:, sl] ** 2)

    egwdffi = jnp.zeros((ncol, pver + 1), dtype=jnp.float64)
    egwdffi = egwdffi.at[:, ktop + 1:kbot].set(
        midpoint_interp(egwdffm[:, sl]))
    egwdffi = jnp.minimum(EGWDFFI_MAX, egwdffi)
    ki = jnp.arange(pver + 1)
    zero_m = ((ki[None, :] >= tend_level[:, None])
              & (ki[None, :] >= ktop) & (ki[None, :] <= kbot))
    egwdffi = jnp.where(zero_m, 0.0, egwdffi)

    tmpi2 = jnp.zeros((ncol, pver + 1), dtype=jnp.float64)
    tmpi2 = tmpi2.at[:, ktop + 1:kbot + 1].set(
        dt * (gravit * rho[:, ktop + 1:kbot + 1]) ** 2
        / (pmid[:, ktop + 1:kbot + 1] - pmid[:, ktop:kbot]))

    zero = jnp.zeros(ncol, dtype=jnp.float64)
    decomp = vd_lu_decomp(zero, egwdffi, tmpi2, rdpm, dt, gravit, zero,
                          ktop + 1, kbot + 1)
    return egwdffi, decomp


def _gw_diff_tend(q, dt, decomp, ntop, nbot):
    """gw_diffusion.F90 gw_diff_tend: dq = (solve(q) - q)/dt."""
    ncol = q.shape[0]
    zero = jnp.zeros(ncol, dtype=jnp.float64)
    qnew = vd_lu_solve(q, decomp, ntop, nbot, zero)
    return (qnew - q) / dt


def _precalc_rhoi(dt, tend_level, pmid, pint, t, gwut, ubm, nm, rdpm,
                  c, q, dse, params):
    """gwd_precalc_rhoi: recomputed interface density feeding gw_ediff
    plus the constituent/DSE diffusion and the KE->heat term."""
    pver, ktop, kbot = params["pver"], params["ktop"], params["kbotbg"]
    rair, nwav = params["rair"], params["nwav"]

    rhoi_kludge = jnp.concatenate([
        pint[:, :1] / (rair * t[:, :1]),
        pint[:, 1:pver] * 2.0 / (rair * (t[:, 1:pver]
                                         + t[:, :pver - 1])),
        pint[:, pver:pver + 1] / (rair * t[:, pver - 1:pver])], axis=1)

    egwdffi, decomp = _gw_ediff(tend_level, gwut, ubm, nm, rhoi_kludge,
                                dt, pmid, rdpm, c, params)

    qtgw = jnp.stack(
        [_gw_diff_tend(q[:, :, m], dt, decomp, ktop + 1, kbot + 1)
         for m in range(q.shape[2])], axis=-1)
    dttdf = _gw_diff_tend(dse, dt, decomp, ktop + 1, kbot + 1)

    dttke = jnp.zeros_like(dttdf)
    sl = slice(ktop, kbot)                    # F midpoints ktop+1..kbotbg
    for l in range(nwav):                     # Fortran accumulation order
        dttke = dttke.at[:, sl].add(c[:, l:l + 1] * gwut[:, sl, l])

    ttgw = dttke + dttdf
    return egwdffi, qtgw, dttdf, dttke, ttgw


def gw_drag_prof(src_level, tend_level, dt, lat, t, ti, pmid, pint,
                 dpm, rdpm, piln, rhoi, nm, ni, ubm, ubi, xv, yv,
                 effgw, c, q, dse, tau, params, do_taper=False):
    """gw_common.F90 gw_drag_prof with ngwv = pgwv > 0. q is
    (ncol, pver, ncnst). Returns (tau, utgw, vtgw, ttgw, qtgw, taucd,
    egwdffi, gwut, dttdf, dttke)."""
    tau = jnp.asarray(tau, dtype=jnp.float64)
    tau = _stress_profiles(src_level, ubi, c, rhoi, ni, t, ti, piln,
                           tau, params)
    taucd = _project_tau(tend_level, tau, ubi, c, xv, yv, params)
    tau, gwut, utgw, vtgw = _tendencies_from_stress_divergence(
        tend_level, do_taper, dt, effgw, lat, dpm, rdpm, c, ubm, t, nm,
        xv, yv, tau, params)
    egwdffi, qtgw, dttdf, dttke, ttgw = _precalc_rhoi(
        dt, tend_level, pmid, pint, t, gwut, ubm, nm, rdpm, c,
        jnp.asarray(q, dtype=jnp.float64),
        jnp.asarray(dse, dtype=jnp.float64), params)
    return (tau, utgw, vtgw, ttgw, qtgw, taucd, egwdffi, gwut, dttdf,
            dttke)


# ---------------------------------------------------------------------------
# gw_common.F90: momentum & energy conservation (C.-C. Chen)
# ---------------------------------------------------------------------------
def momentum_energy_conservation(tend_level, dt, taucd, pint, pdel, u,
                                 v, dudt, dvdt, dsdt, utgw, vtgw, ttgw,
                                 params):
    """Adds the below-source compensation tendencies; returns updated
    (dudt, dvdt, dsdt, utgw, vtgw, ttgw)."""
    pver, gravit = params["pver"], params["gravit"]
    ncol = u.shape[0]
    dudt, dvdt, dsdt, utgw, vtgw, ttgw = (
        jnp.asarray(a, dtype=jnp.float64)
        for a in (dudt, dvdt, dsdt, utgw, vtgw, ttgw))
    ktl = int(tend_level.min())

    dz = jnp.zeros(ncol, dtype=jnp.float64)
    for k in range(pver, ktl, -1):            # F k = pver .. min+1
        dz = dz + jnp.where(k > tend_level, pdel[:, k - 1] / gravit,
                            0.0)

    tl = tend_level[:, None, None]
    taucd_tl = jnp.take_along_axis(
        taucd, jnp.broadcast_to(tl, (ncol, 1, 4)), axis=1)[:, 0, :]
    ut_dz = -(taucd_tl[:, 1] + taucd_tl[:, 0]) / dz    # east + west
    vt_dz = -(taucd_tl[:, 3] + taucd_tl[:, 2]) / dz    # north + south

    for k in range(ktl + 1, pver + 1):
        mk = k > tend_level
        dudt = dudt.at[:, k - 1].add(jnp.where(mk, ut_dz, 0.0))
        dvdt = dvdt.at[:, k - 1].add(jnp.where(mk, vt_dz, 0.0))
        utgw = utgw.at[:, k - 1].add(jnp.where(mk, ut_dz, 0.0))
        vtgw = vtgw.at[:, k - 1].add(jnp.where(mk, vt_dz, 0.0))

    dE = jnp.zeros(ncol, dtype=jnp.float64)
    for k in range(1, pver + 1):              # Fortran summation order
        j = k - 1
        dE = dE + pdel[:, j] * (
            dsdt[:, j]
            + dudt[:, j] * (u[:, j] + dudt[:, j] * 0.5 * dt)
            + dvdt[:, j] * (v[:, j] + dvdt[:, j] * 0.5 * dt))
    pint_tl = jnp.take_along_axis(pint, tend_level[:, None],
                                  axis=1)[:, 0]
    dE = dE / (pint[:, pver] - pint_tl)

    for k in range(ktl + 1, pver + 1):
        mk = k > tend_level
        dsdt = dsdt.at[:, k - 1].add(jnp.where(mk, -dE, 0.0))
        ttgw = ttgw.at[:, k - 1].add(jnp.where(mk, -dE, 0.0))
    return dudt, dvdt, dsdt, utgw, vtgw, ttgw
