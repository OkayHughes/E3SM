"""Port of the EAMv3 tropopause finder, eam/src/physics/cam/
tropopause.F90: the twmo core (Reichler et al. [2003] WMO lapse-rate
tropopause), the climatology fallback, the hybrid Stobie-Linoz
algorithm, and the tropopause_find primary+backup dispatch.

PORT_NOTES
----------
- EAMv3 production invokes exactly three method combinations (grep of
  components/eam/src): the default TWMO primary + CLIMATE backup
  (aer_rad_props, tropopause_output), explicit TWMO+CLIMATE
  (prescribed_volcaero) and HYBSTOB+CLIMATE (mozart chemistry,
  modal_aero_wateruptake). Those three methods and the dispatch are
  ported; ANALYTIC/STOBIE/WMO/E90 and tropopause_findChemTrop have no
  EAMv3 callers and are not ported.
- Index convention: arrays are C-ordered (ncol, nlev), level 0 = model
  top; pint/zi carry nlev+1 interfaces. The returned tropLev is
  0-BASED (Fortran tropLev - 1); NOTFOUND = -1 as in the Fortran.
  Found tropLev is the layer k with pint[k] <= tropP < pint[k+1]
  (Fortran `do k = pver,2,-1; if (tP >= pint(i,k))` scan), so level 0
  can never be returned by twmo/climate; hybridstobie indexes pmid
  directly and defaults to Fortran level 1 (here 0).
- physics_state/pbuf coupling becomes plain arrays. The climatology
  (module tropp_p_loc, read from the trop_climo file in
  tropopause_init) enters as tropp_clim (ncol, 12) [Pa] — the
  already-column-interpolated monthly values — together with days(12)
  (noleap day-of-year of the file's mid-month dates, 16, 45, ..., 350)
  and the scalar calday from get_curr_calday().
- Constants, verbatim: gam = -0.002 K/m, plimu = 45000 Pa,
  pliml = 7500 Pa, deltaz = 2000 m (twmo); ALPHA = 0.03,
  min_Stobie_Pressure = 40.E2 Pa, max_Linoz_Pressure = 208.E2 Pa
  (hybridstobie); fillvalue = 1.e+20 (cam_history_support);
  cnst_kap = cappa, cnst_faktor = -gravit/rair (physconst).
- twmo: exact sequential port. Pair index q in [0, nlev-2] holds the
  half-level quantities of Fortran pair (j-1, j), j = q+2 (1-based);
  the main loop starts at the pair ABOVE the bottom pair and, because
  the pair quantities are recomputed unconditionally each iteration,
  dtdz0/pmk0 are always the (q+1) pair values (no carry). A candidate
  needs dtdz > gam and pm <= plimu; ptph interpolates the gam crossing
  in p**kap space only when dtdz0 < gam, else ptph = pm; then
  pliml <= ptph <= plimu; then the 2-km test scans pairs qq = q..0
  (skip pm2 > ptph, accept on pm2 < p2km, else accumulate dtdz and
  discard if any running mean <= gam; running out of pairs also
  accepts). First success scanning up wins; trp stays -99.0 on
  failure. Masked-lane guards (log argument, 0-count mean) never
  change a lane the Fortran evaluates: dtdz0 < gam < dtdz implies the
  crossing pressure is strictly inside the pair interval.
- Fortran quirks kept: tropopause_twmo writes tropP (and interpolated
  tropT/tropZ) whenever trp > 0 even if the pint level scan failed
  (tropLev stays NOTFOUND; unreachable for pliml = 7500 on any real
  grid) — the port clips the gather index for that unreachable case.
  tropopause_climate writes the outputs for every previously-NOTFOUND
  column unconditionally. interpolateT falls back to the midpoint
  value at the (unreachable) column top/bottom, where the Fortran
  would return an undefined local; interpolateZ has no such guards
  (all its indices are valid).
- Climatology time interpolation: wrap branches for calday < days(1)
  and calday >= days(12) use the 365-day (noleap) year; otherwise
  last = max{m in 1..11: calday >= days(m)}; dels clipped to [0, 1].
- hybridstobie: levels with pmid < 40 hPa are skipped entirely (the
  cycle also skips the Linoz check); ltrop_trop is the bottom-most
  level attaining the minimum of ALPHA*t - log10(pmid) (strict <
  keeps the first minimum found scanning up); ltrop_linoz is the
  bottom-most non-skipped level with pmid < 208 hPa; both default to
  Fortran level 1 (top). tropLev = min of the two; tropP/T/Z are the
  midpoint values (no interpolation). It cannot fail, so a backup
  after it never engages.
- tropopause_find applies the backup unconditionally (the Fortran's
  `any(tropLev == NOTFOUND)` guard is an optimization only — methods
  touch nothing but NOTFOUND columns).
- Like-for-like float64; the only measured deviation source is fused
  multiply-add contraction by gfortran -O2 on aarch64 (e.g. the
  climatology a + dels*(b-a)), <= 1 ulp, inside the 1e-12 replay
  tolerance.
"""

import jax
import jax.numpy as jnp

from .constants import CAPPA, GRAVIT, RAIR

NOTFOUND = -1
FILLVALUE = 1.0e20          # cam_history_support fillvalue
ALPHA = 0.03

# twmo parameters (tropopause_twmo locals, verbatim)
GAM = -0.002                # K/m
PLIMU = 45000.0             # Pa
PLIML = 7500.0              # Pa
DELTAZ = 2000.0             # m

# hybridstobie parameters (verbatim)
MIN_STOBIE_PRESSURE = 40.0e2   # Pa
MAX_LINOZ_PRESSURE = 208.0e2   # Pa

# physical constants as tropopause_init derives them
CNST_KAP = CAPPA
CNST_FAKTOR = -GRAVIT / RAIR
CNST_KA1 = CNST_KAP - 1.0


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _gather(a, lev):
    """a[:, lev] per column (lev (ncol,) int)."""
    return jnp.take_along_axis(a, lev[:, None], axis=1)[:, 0]


def _find_level(pint, tP):
    """Fortran `do k = pver, 2, -1; if (tP >= pint(i,k)) exit` — the
    bottom-most 1-based k in [2, pver] with pint(k) <= tP. Returns the
    0-BASED level, NOTFOUND if none."""
    nlevp = pint.shape[1]
    ks = jnp.arange(2, nlevp)                       # 1-based k = 2..pver
    cond = tP[:, None] >= pint[:, 1:nlevp - 1]      # pint(i, k), k=2..pver
    lev_f = jnp.max(jnp.where(cond, ks[None, :], NOTFOUND), axis=1)
    return jnp.where(lev_f > 0, lev_f - 1, NOTFOUND)


def _interpolate_t(t, pmid, lev0, tP):
    """tropopause_interpolateT: T linear in log(p) between midpoints.
    Falls back to t[lev] at the (unreachable) top/bottom where the
    Fortran result is undefined."""
    nlev = t.shape[1]
    t_l = _gather(t, lev0)
    pm_l = _gather(pmid, lev0)
    t_up = _gather(t, jnp.maximum(lev0 - 1, 0))
    pm_up = _gather(pmid, jnp.maximum(lev0 - 1, 0))
    t_dn = _gather(t, jnp.minimum(lev0 + 1, nlev - 1))
    pm_dn = _gather(pmid, jnp.minimum(lev0 + 1, nlev - 1))

    dtdlogp_a = (t_l - t_up) / (jnp.log(pm_l) - jnp.log(pm_up))
    val_a = t_l + (jnp.log(tP) - jnp.log(pm_l)) * dtdlogp_a
    dtdlogp_b = (t_dn - t_l) / (jnp.log(pm_dn) - jnp.log(pm_l))
    val_b = t_l + (jnp.log(tP) - jnp.log(pm_l)) * dtdlogp_b

    above = tP < pm_l
    return jnp.where(
        tP == pm_l, t_l,
        jnp.where(above,
                  jnp.where(lev0 > 0, val_a, t_l),
                  jnp.where(lev0 < nlev - 1, val_b, t_l)))


def _interpolate_z(zm, zi, pmid, pint, lev0, tP):
    """tropopause_interpolateZ: Z linear in log(p) between the midpoint
    and the interface on the tP side (no level guards in the
    Fortran)."""
    zm_l = _gather(zm, lev0)
    pm_l = _gather(pmid, lev0)
    zi_a = _gather(zi, lev0)                # zi(tropLev): top interface
    pi_a = _gather(pint, lev0)
    zi_b = _gather(zi, lev0 + 1)            # zi(tropLev+1): bottom
    pi_b = _gather(pint, lev0 + 1)

    dzdlogp_a = (zm_l - zi_a) / (jnp.log(pm_l) - jnp.log(pi_a))
    dzdlogp_b = (zm_l - zi_b) / (jnp.log(pm_l) - jnp.log(pi_b))
    dzdlogp = jnp.where(tP < pm_l, dzdlogp_a, dzdlogp_b)
    return jnp.where(tP == pm_l, zm_l,
                     zm_l + (jnp.log(tP) - jnp.log(pm_l)) * dzdlogp)


# ---------------------------------------------------------------------------
# twmo core (subroutine twmo, exact sequential port)
# ---------------------------------------------------------------------------
def twmo(t, p, plimu=PLIMU, pliml=PLIML, gam=GAM):
    """Reichler et al. [2003] tropopause pressure per column.

    t, p: (ncol, nlev) temperature/pressure midpoints, level 0 = top.
    Returns trp (ncol,): tropopause pressure [Pa], -99.0 if not found.
    """
    t = jnp.asarray(t, dtype=jnp.float64)
    p = jnp.asarray(p, dtype=jnp.float64)
    ncol, nlev = t.shape

    # half-level (pair) quantities; pair q holds Fortran pair (j-1, j)
    # with j = q+2 (1-based), i.e. 0-based levels (q, q+1)
    pk = p ** CNST_KAP
    pmk = 0.5 * (pk[:, :-1] + pk[:, 1:])
    pm = pmk ** (1.0 / CNST_KAP)
    a = (t[:, :-1] - t[:, 1:]) / (pk[:, :-1] - pk[:, 1:])
    b = t[:, 1:] - a * pk[:, 1:]
    tm = a * pmk + b
    dtdp = a * CNST_KAP * pm ** CNST_KA1
    dtdz = CNST_FAKTOR * dtdp * pm / tm

    # main_loop pairs: q = nlev-3 .. 0 (Fortran j = level-1 .. 2);
    # in_loop pairs:   qq = q .. 0    (Fortran jj = j .. 2)
    qs = jnp.arange(nlev - 3, -1, -1)
    pm_in = pm.T[qs]                 # (nq, ncol): inner xs (all pairs
    dz_in = dtdz.T[qs]               # a candidate can ever visit)

    def outer(carry, x):
        found, trp = carry
        q, pmk_c, pmk0, pm_c, tm_c, dz_c, dz0 = x

        # dt/dz valid? (cycle conditions inverted)
        ok1 = (dz_c > gam) & (pm_c <= plimu)
        # tropopause pressure: interpolate the gam crossing in p**kap
        # space when the pair below is still sub-critical
        ag = (dz_c - dz0) / (pmk_c - pmk0)
        bg = dz0 - ag * pmk0
        arg = (gam - bg) / ag
        # guard: on lanes with dz0 < gam < dz_c (the only ones that
        # select this branch) arg is strictly inside the pair interval,
        # hence > 0; the where only sanitizes masked lanes
        ptph = jnp.where(dz0 < gam,
                         jnp.exp(jnp.log(jnp.where(arg > 0.0, arg, 1.0))
                                 / CNST_KAP),
                         pm_c)
        ok2 = ok1 & (ptph >= pliml) & (ptph <= plimu)

        # 2nd test: mean dt/dz over the 2 km above must stay > gam
        p2km = ptph + DELTAZ * (pm_c / tm_c) * CNST_FAKTOR

        def inner(icarry, ix):
            asum, icount, discard, exited = icarry
            pm2, dz2, qq = ix
            active = (qq <= q) & ~discard & ~exited
            skip = pm2 > ptph                    # doesn't happen
            exi = pm2 < p2km                     # ptropo is valid
            acc = active & ~skip & ~exi
            asum = asum + jnp.where(acc, dz2, 0.0)
            icount = icount + jnp.where(acc, 1, 0)
            aquer = asum / jnp.maximum(icount, 1)  # guarded: only read
            discard = discard | (acc & (aquer <= gam))  # where acc
            exited = exited | (active & ~skip & exi)
            return (asum, icount, discard, exited), None

        init = (jnp.zeros(ncol), jnp.zeros(ncol, dtype=jnp.int32),
                jnp.zeros(ncol, dtype=bool), jnp.zeros(ncol, dtype=bool))
        (_, _, discard, _), _ = jax.lax.scan(inner, init,
                                             (pm_in, dz_in, qs))

        passed = ok2 & ~discard
        newly = passed & ~found
        trp = jnp.where(newly, ptph, trp)
        return (found | passed, trp), None

    xs = (qs, pmk.T[qs], pmk.T[qs + 1], pm.T[qs], tm.T[qs],
          dtdz.T[qs], dtdz.T[qs + 1])
    init = (jnp.zeros(ncol, dtype=bool), jnp.full(ncol, -99.0))
    (_, trp), _ = jax.lax.scan(outer, init, xs)
    return trp


# ---------------------------------------------------------------------------
# per-method wrappers (each fills only previously-NOTFOUND columns)
# ---------------------------------------------------------------------------
def tropopause_twmo(t, pmid, pint, zm, zi, tropLev, tropP, tropT, tropZ):
    trp = twmo(t, pmid)
    ok = trp > 0.0
    upd = (tropLev == NOTFOUND) & ok
    lev0 = _find_level(pint, trp)
    safe = jnp.maximum(lev0, 0)     # unreachable NOTFOUND-lev gather
    tropLev = jnp.where(upd & (lev0 != NOTFOUND), lev0, tropLev)
    tropP = jnp.where(upd, trp, tropP)
    tropT = jnp.where(upd, _interpolate_t(t, pmid, safe, trp), tropT)
    tropZ = jnp.where(upd, _interpolate_z(zm, zi, pmid, pint, safe, trp),
                      tropZ)
    return tropLev, tropP, tropT, tropZ


def tropopause_climate(t, pmid, pint, zm, zi, tropp_clim, days, calday,
                       tropLev, tropP, tropT, tropZ):
    """tropp_clim: (ncol, 12) monthly climatological tropopause
    pressures at the columns [Pa]; days: (12,) noleap day-of-year of
    the climatology months; calday: scalar current day of year."""
    tropp_clim = jnp.asarray(tropp_clim, dtype=jnp.float64)
    days = jnp.asarray(days, dtype=jnp.float64)
    calday = jnp.asarray(calday, dtype=jnp.float64)

    wrap_lo = calday < days[0]
    wrap_hi = calday >= days[11]
    wrap = wrap_lo | wrap_hi
    # mid branch: last = max{m in 1..11 : calday >= days(m)} (1-based)
    mm = jnp.sum(calday >= days[:11]).astype(jnp.int32)
    mm_s = jnp.clip(mm, 1, 11)      # gather-safe on wrap lanes
    dels_wrap = jnp.where(
        wrap_lo,
        (365.0 + calday - days[11]) / (365.0 + days[0] - days[11]),
        (calday - days[11]) / (365.0 + days[0] - days[11]))
    dels_mid = (calday - days[mm_s - 1]) / (days[mm_s] - days[mm_s - 1])
    dels = jnp.where(wrap, dels_wrap, dels_mid)
    dels = jnp.maximum(jnp.minimum(1.0, dels), 0.0)
    last0 = jnp.where(wrap, 11, mm_s - 1)   # 0-based month indices
    next0 = jnp.where(wrap, 0, mm_s)

    tP = (tropp_clim[:, last0]
          + dels * (tropp_clim[:, next0] - tropp_clim[:, last0]))

    upd = tropLev == NOTFOUND       # outputs written unconditionally
    lev0 = _find_level(pint, tP)    # for previously-unfound columns
    safe = jnp.maximum(lev0, 0)
    tropLev = jnp.where(upd & (lev0 != NOTFOUND), lev0, tropLev)
    tropP = jnp.where(upd, tP, tropP)
    tropT = jnp.where(upd, _interpolate_t(t, pmid, safe, tP), tropT)
    tropZ = jnp.where(upd, _interpolate_z(zm, zi, pmid, pint, safe, tP),
                      tropZ)
    return tropLev, tropP, tropT, tropZ


def tropopause_hybridstobie(t, pmid, zm, tropLev, tropP, tropT, tropZ):
    """Stobie-Linoz hybrid; cannot fail (defaults to level 0)."""
    t = jnp.asarray(t, dtype=jnp.float64)
    pmid = jnp.asarray(pmid, dtype=jnp.float64)
    nlev = t.shape[1]
    k0 = jnp.arange(nlev)

    considered = pmid >= MIN_STOBIE_PRESSURE   # `< min -> cycle`
    stobie = ALPHA * t - jnp.log10(pmid)
    vals = jnp.where(considered, stobie, jnp.inf)
    vmin = jnp.min(vals, axis=1)
    # strict < while scanning k = pver..1: the bottom-most minimum wins
    ltrop_trop = jnp.max(jnp.where(vals == vmin[:, None], k0[None, :],
                                   NOTFOUND), axis=1)
    ltrop_trop = jnp.where(jnp.any(considered, axis=1), ltrop_trop, 0)

    linoz_ok = considered & (pmid < MAX_LINOZ_PRESSURE)
    ltrop_linoz = jnp.max(jnp.where(linoz_ok, k0[None, :], NOTFOUND),
                          axis=1)
    ltrop_linoz = jnp.where(jnp.any(linoz_ok, axis=1), ltrop_linoz, 0)

    lev0 = jnp.minimum(ltrop_trop, ltrop_linoz)
    upd = tropLev == NOTFOUND
    tropLev = jnp.where(upd, lev0, tropLev)
    tropP = jnp.where(upd, _gather(pmid, lev0), tropP)
    tropT = jnp.where(upd, _gather(t, lev0), tropT)
    tropZ = jnp.where(upd, _gather(jnp.asarray(zm, dtype=jnp.float64),
                                   lev0), tropZ)
    return tropLev, tropP, tropT, tropZ


# ---------------------------------------------------------------------------
# dispatch (tropopause_find; default primary/backup as in the Fortran)
# ---------------------------------------------------------------------------
def tropopause_find(t, pmid, pint, zm, zi, tropp_clim=None, days=None,
                    calday=None, primary="twmo", backup="climate"):
    """Primary+backup tropopause search. Returns (tropLev, tropP,
    tropT, tropZ); tropLev is 0-based with NOTFOUND=-1, the others are
    FILLVALUE where nothing was found. Algorithms: "twmo", "climate",
    "hybridstobie", "none" (only the combinations EAMv3 uses are
    ported; see PORT_NOTES)."""
    t = jnp.asarray(t, dtype=jnp.float64)
    pmid = jnp.asarray(pmid, dtype=jnp.float64)
    pint = jnp.asarray(pint, dtype=jnp.float64)
    zm = jnp.asarray(zm, dtype=jnp.float64)
    zi = jnp.asarray(zi, dtype=jnp.float64)
    ncol = t.shape[0]

    state = (jnp.full(ncol, NOTFOUND, dtype=jnp.int32),
             jnp.full(ncol, FILLVALUE),
             jnp.full(ncol, FILLVALUE),
             jnp.full(ncol, FILLVALUE))

    for alg in (primary, backup):
        if alg == "none":
            continue
        elif alg == "twmo":
            state = tropopause_twmo(t, pmid, pint, zm, zi, *state)
        elif alg == "climate":
            if tropp_clim is None or days is None or calday is None:
                raise ValueError(
                    "climate algorithm needs tropp_clim/days/calday")
            state = tropopause_climate(t, pmid, pint, zm, zi,
                                       tropp_clim, days, calday, *state)
        elif alg == "hybridstobie":
            state = tropopause_hybridstobie(t, pmid, zm, *state)
        else:
            raise ValueError(f"unsupported algorithm {alg!r}")
    return state
