"""Port of the ZM convective transport routines
(eam/src/physics/cam/zm/zm_transport.F90): zm_transport_tracer
(convective tracer transport on the ZM plume mass fluxes) and
zm_transport_momentum (convective momentum transport with the
Gregory et al. pressure-gradient term and the Boville & Bretherton
2003 kinetic-energy-dissipation heating).

PORT_NOTES
----------
- Gathered-array convention exactly as zm_conv_intr.F90 passes it:
  mu/md [mb/s], du/eu/ed [1/s], dp/dpdry [mb] are GATHERED rows
  0..lengath-1 (the direct zm_conv_main outputs mflx_up/mflx_dn/
  detr_up/entr_up/entr_dn/p_del); q/fracis/wind_in and all scattered
  outputs are full-column arrays indexed through ideep. Indices here
  are 0-based (goldens store Fortran 1-based; tests shift).
- ideep/jt/mx are concrete numpy indices (eager-mode gather, same
  style as eam_jax.zm_conv; not jit-compatible as written).
- Tracer loop starts at m=1 (Fortran m=2): water vapor is never
  transported here; slices with doconvtran False stay exactly zero
  (matching the zeroed ptend%q of the intr layer / harness driver).
- 'dry' constituents (cnst_get_type_byind == 'dry') swap dp for the
  gathered dpdry and rescale du/eu/ed by dp/dpdry, exactly as the
  Fortran does per constituent.
- chat interface interpolation: geometric (log) average when the
  normalized layer difference exceeds 1e-6 and both layers are
  nonnegative, else arithmetic mean; the geometric branch floors both
  values at maxc*1e-12. bfb_log == log (SCREAM_CONFIG_IS_CMAKE
  undefined in EAM builds).
- The updraft (bottom->top) and downdraft (top->bottom) in-cloud
  recursions are sequential python level loops with per-column masks
  (mu + du*dp > mbsth resp. md < -mbsth), data flow as in the
  Fortran; everything else is vectorized over (column, level).
- The flux limiter block min(chat,const) is antisymmetric between
  adjacent layers, so column mass sum(dqdt*dptmp) telescopes to the
  (zero) fluxes at the cloud-top interface: conservation to roundoff.
  The |netflux| < max(fluxin,fluxout)*1e-12 clip is kept verbatim.
- zm_param%zm_microp gates ONLY the negative-tracer conservation
  fixer (a per-column sequential borrow up/down of the negative
  mass), which deliberately ADDS tracer when the column cannot absorb
  the deficit (`dcondt(k) - negadt`); it is ported as concrete numpy
  loops in the exact Fortran order. Unlike zm_conv itself, both
  branches are goldened (the fixer needs no zm_microphysics code).
- zm_transport_momentum has NO runtime switch: momcu = momcd = 0.4
  are compile-time parameters (no zmconv_mom* namelist) and nwind = 2
  always (the KE fixer hard-codes components 1 and 2). The
  boundary-level downdraft seed keeps the Fortran's odd
  parenthesization verbatim: only the pgd term is divided by md
  (`(-ed*wm*dp) - pgd*dp/md`), unlike the in-loop recursion.
- windf (used only by the KE fixer) is zero-initialized and only
  filled for k >= ktm, exactly as in the Fortran; wind0 is the full
  gathered wind.
- Constants verbatim: mbsth=1e-15, small=1e-36, cdifr_min=1e-6,
  maxc_factor=1e-12, flux_factor=1e-12, momcu=momcd=0.4.
- Everything float64.
"""

import numpy as np
import jax.numpy as jnp

__all__ = ["zm_transport_tracer", "zm_transport_momentum"]

MBSTH = 1.0e-15        # mass-flux threshold [mb/s]
SMALL = 1.0e-36        # divide-by-zero guard
CDIFR_MIN = 1.0e-6     # geometric-average threshold
MAXC_FACTOR = 1.0e-12  # floor factor for the geometric average
FLUX_FACTOR = 1.0e-12  # netflux round-off clip
MOMCU = 0.4            # updraft pressure-gradient constant
MOMCD = 0.4            # downdraft pressure-gradient constant


def _chat(const):
    """Environment interface values (geometric/arithmetic average).
    Interface k sits between midpoints k-1 and k; k=0 uses km1=0."""
    n, pver = const.shape
    km1 = np.maximum(np.arange(pver) - 1, 0)
    a = const[:, km1]           # const(km1)
    b = const                   # const(k)
    minc = jnp.minimum(a, b)
    maxc = jnp.maximum(a, b)
    cdifr = jnp.where(minc < 0.0, 0.0,
                      jnp.abs(b - a) / jnp.maximum(maxc, SMALL))
    geo = cdifr > CDIFR_MIN
    cabv = jnp.maximum(a, maxc * MAXC_FACTOR)
    cbel = jnp.maximum(b, maxc * MAXC_FACTOR)
    cabv_s = jnp.where(geo, cabv, 2.0)
    cbel_s = jnp.where(geo, cbel, 1.0)
    diff_s = jnp.where(geo, cabv - cbel, 1.0)
    geo_val = jnp.log(cabv_s / cbel_s) / diff_s * cabv_s * cbel_s
    return jnp.where(geo, geo_val, 0.5 * (b + a))


def zm_transport_tracer(q, fracis, doconvtran, is_dry, mu, md, du,
                        eu, ed, dp, dpdry, jt, mx, ideep, dt,
                        zm_microp=False):
    """Convective tracer transport of ncnst constituents (m=0 is
    water vapor, always skipped). q, fracis: (ncol, pver, ncnst)
    ungathered; mu/md [mb/s], du/eu/ed [1/s], dp/dpdry [mb] gathered;
    jt/mx 0-based gathered plume top/bottom; ideep 0-based gather
    indices (length lengath); dt [s]. Returns dqdt (ncol, pver,
    ncnst) [1/s], zero outside gathered columns and for untransported
    constituents."""
    q = jnp.asarray(q, dtype=jnp.float64)
    fracis = jnp.asarray(fracis, dtype=jnp.float64)
    ncol, pver, ncnst = q.shape
    ideep = np.asarray(ideep)
    n = ideep.size
    dqdt = jnp.zeros((ncol, pver, ncnst))
    if n == 0:
        return dqdt

    jt = np.asarray(jt)[:n]
    mx = np.asarray(mx)[:n]
    mu = jnp.asarray(mu, dtype=jnp.float64)[:n]
    md = jnp.asarray(md, dtype=jnp.float64)[:n]
    du = jnp.asarray(du, dtype=jnp.float64)[:n]
    eu = jnp.asarray(eu, dtype=jnp.float64)[:n]
    ed = jnp.asarray(ed, dtype=jnp.float64)[:n]
    dp = jnp.asarray(dp, dtype=jnp.float64)[:n]
    dpdry = jnp.asarray(dpdry, dtype=jnp.float64)[:n]

    ktm = int(np.min(jt))
    kbm = int(np.min(mx))
    rows = np.arange(n)

    for m in range(1, ncnst):
        if not doconvtran[m]:
            continue
        if is_dry[m]:
            dptmp = dpdry
            dutmp = du * dp / dpdry
            eutmp = eu * dp / dpdry
            edtmp = ed * dp / dpdry
        else:
            dptmp = dp
            dutmp = du
            eutmp = eu
            edtmp = ed

        const = q[ideep, :, m]
        fisg = fracis[ideep, :, m]

        chat = _chat(const)
        conu = jnp.array(chat)
        cond = jnp.array(chat)
        dcondt = jnp.zeros((n, pver))

        mupdudp = mu + dutmp * dptmp

        # levels adjacent to top and bottom (Fortran k=2, kk=pver)
        mk = mupdudp[:, pver - 1] > MBSTH
        conu = conu.at[:, pver - 1].set(jnp.where(
            mk,
            (eutmp[:, pver - 1] * fisg[:, pver - 1]
             * const[:, pver - 1] * dptmp[:, pver - 1])
            / jnp.where(mk, mupdudp[:, pver - 1], 1.0),
            conu[:, pver - 1]))
        mk = md[:, 1] < -MBSTH
        cond = cond.at[:, 1].set(jnp.where(
            mk,
            (-edtmp[:, 0] * fisg[:, 0] * const[:, 0] * dptmp[:, 0])
            / jnp.where(mk, md[:, 1], 1.0),
            cond[:, 1]))

        # updraft from bottom to top
        for kk in range(pver - 2, -1, -1):
            mk = mupdudp[:, kk] > MBSTH
            num = (mu[:, kk + 1] * conu[:, kk + 1]
                   + eutmp[:, kk] * fisg[:, kk] * const[:, kk]
                   * dptmp[:, kk])
            conu = conu.at[:, kk].set(jnp.where(
                mk, num / jnp.where(mk, mupdudp[:, kk], 1.0),
                conu[:, kk]))

        # downdraft from top to bottom
        for k in range(2, pver):
            mk = md[:, k] < -MBSTH
            num = (md[:, k - 1] * cond[:, k - 1]
                   - edtmp[:, k - 1] * fisg[:, k - 1]
                   * const[:, k - 1] * dptmp[:, k - 1])
            cond = cond.at[:, k].set(jnp.where(
                mk, num / jnp.where(mk, md[:, k], 1.0), cond[:, k]))

        # limited flux divergence, k = ktm..pver-1 (vectorized)
        ks = np.arange(ktm, pver)
        km1 = np.maximum(ks - 1, 0)
        kp1 = np.minimum(ks + 1, pver - 1)
        fluxin = (mu[:, kp1] * conu[:, kp1]
                  + mu[:, ks] * jnp.minimum(chat[:, ks], const[:, km1])
                  - (md[:, ks] * cond[:, ks]
                     + md[:, kp1] * jnp.minimum(chat[:, kp1],
                                                const[:, kp1])))
        fluxout = (mu[:, ks] * conu[:, ks]
                   + mu[:, kp1] * jnp.minimum(chat[:, kp1],
                                              const[:, ks])
                   - (md[:, kp1] * cond[:, kp1]
                      + md[:, ks] * jnp.minimum(chat[:, ks],
                                                const[:, ks])))
        netflux = fluxin - fluxout
        netflux = jnp.where(
            jnp.abs(netflux) < jnp.maximum(fluxin, fluxout)
            * FLUX_FACTOR, 0.0, netflux)
        dcondt = dcondt.at[:, ktm:].set(netflux / dptmp[:, ks])

        # cloud-base layer and below, k = kbm..pver-1
        kb = np.arange(kbm, pver)
        km1b = np.maximum(kb - 1, 0)
        at_mx = kb[None, :] == mx[:, None]
        below = kb[None, :] > mx[:, None]
        fluxin_b = (mu[:, kb] * jnp.minimum(chat[:, kb],
                                            const[:, km1b])
                    - md[:, kb] * cond[:, kb])
        fluxout_b = (mu[:, kb] * conu[:, kb]
                     - md[:, kb] * jnp.minimum(chat[:, kb],
                                               const[:, kb]))
        net_b = fluxin_b - fluxout_b
        net_b = jnp.where(
            jnp.abs(net_b) < jnp.maximum(fluxin_b, fluxout_b)
            * FLUX_FACTOR, 0.0, net_b)
        dcondt = dcondt.at[:, kbm:].set(
            jnp.where(at_mx, net_b / dptmp[:, kb],
                      jnp.where(below, 0.0, dcondt[:, kbm:])))

        # conservation fixer for ZM microphysics (sequential borrow)
        if zm_microp:
            dc = np.asarray(dcondt).copy()
            cn = np.asarray(const)
            dpt = np.asarray(dptmp)
            for i in range(n):
                for k in range(jt[i], mx[i] + 1):
                    if dc[i, k] * dt + cn[i, k] < 0.0:
                        negadt = dc[i, k] + cn[i, k] / dt
                        dc[i, k] = -cn[i, k] / dt
                        for kk in range(k + 1, mx[i] + 1):
                            if negadt < 0.0 and \
                                    dc[i, kk] * dt + cn[i, kk] > 0.0:
                                qtmp = dc[i, kk] + negadt \
                                    * dpt[i, k] / dpt[i, kk]
                                if qtmp * dt + cn[i, kk] > 0.0:
                                    dc[i, kk] = qtmp
                                    negadt = 0.0
                                else:
                                    negadt = negadt \
                                        + (cn[i, kk] / dt + dc[i, kk]) \
                                        * dpt[i, kk] / dpt[i, k]
                                    dc[i, kk] = -cn[i, kk] / dt
                        for kk in range(k - 1, jt[i] - 1, -1):
                            if negadt < 0.0 and \
                                    dc[i, kk] * dt + cn[i, kk] > 0.0:
                                qtmp = dc[i, kk] + negadt \
                                    * dpt[i, k] / dpt[i, kk]
                                if qtmp * dt + cn[i, kk] > 0.0:
                                    dc[i, kk] = qtmp
                                    negadt = 0.0
                                else:
                                    negadt = negadt \
                                        + (cn[i, kk] / dt + dc[i, kk]) \
                                        * dpt[i, kk] / dpt[i, k]
                                    dc[i, kk] = -cn[i, kk] / dt
                        if negadt < 0.0:
                            dc[i, k] = dc[i, k] - negadt
            dcondt = jnp.asarray(dc)

        dqdt = dqdt.at[ideep, :, m].set(dcondt)

    return dqdt


def zm_transport_momentum(wind_in, mu, md, du, eu, ed, dp, jt, mx,
                          ideep, dt):
    """Convective momentum transport of nwind=2 components (u, v).
    wind_in: (ncol, pver, 2) ungathered; fluxes gathered as in
    zm_transport_tracer; dt [s] (ztodt). Returns a dict with
    wind_tend [m/s2], pguall/pgdall (apparent PG force), icwu/icwd
    (in-cloud winds) all (ncol, pver, 2), and seten (ncol, pver)
    [J/kg/s] dry static energy tendency."""
    wind_in = jnp.asarray(wind_in, dtype=jnp.float64)
    ncol, pver, nwind = wind_in.shape
    assert nwind == 2
    ideep = np.asarray(ideep)
    n = ideep.size

    wind_tend = jnp.zeros((ncol, pver, 2))
    pguall = jnp.zeros((ncol, pver, 2))
    pgdall = jnp.zeros((ncol, pver, 2))
    icwu = jnp.array(wind_in)
    icwd = jnp.array(wind_in)
    seten = jnp.zeros((ncol, pver))
    if n == 0:
        return dict(wind_tend=wind_tend, pguall=pguall, pgdall=pgdall,
                    icwu=icwu, icwd=icwd, seten=seten)

    jt = np.asarray(jt)[:n]
    mx = np.asarray(mx)[:n]
    mu = jnp.asarray(mu, dtype=jnp.float64)[:n]
    md = jnp.asarray(md, dtype=jnp.float64)[:n]
    du = jnp.asarray(du, dtype=jnp.float64)[:n]
    dp = jnp.asarray(dp, dtype=jnp.float64)[:n]
    eu = jnp.asarray(eu, dtype=jnp.float64)[:n]
    ed = jnp.asarray(ed, dtype=jnp.float64)[:n]

    ktm = int(np.min(jt))
    kbm = int(np.min(mx))

    wind0 = jnp.zeros((n, pver, 2))
    windf = jnp.zeros((n, pver, 2))
    mflux = jnp.zeros((n, pver + 1, 2))
    gseten = jnp.zeros((n, pver))
    mupdudp = mu + du * dp

    for m in range(2):
        wm = wind_in[ideep, :, m]
        wind0 = wind0.at[:, :, m].set(wm)

        # interfaces: arithmetic mean (k=0 uses km1=0)
        km1 = np.maximum(np.arange(pver) - 1, 0)
        wi = 0.5 * (wm + wm[:, km1])
        wiu = jnp.array(wi)
        wid = jnp.array(wi)

        # pressure perturbation terms
        pgu = jnp.zeros((n, pver))
        pgd = jnp.zeros((n, pver))
        ks = np.arange(1, pver - 1)
        mududp = (mu[:, ks] * (wm[:, ks] - wm[:, ks - 1])
                  / dp[:, ks - 1]
                  + mu[:, ks + 1] * (wm[:, ks + 1] - wm[:, ks])
                  / dp[:, ks])
        mddudp = (md[:, ks] * (wm[:, ks] - wm[:, ks - 1])
                  / dp[:, ks - 1]
                  + md[:, ks + 1] * (wm[:, ks + 1] - wm[:, ks])
                  / dp[:, ks])
        pgu = pgu.at[:, 1:pver - 1].set(-MOMCU * 0.5 * mududp)
        pgd = pgd.at[:, 1:pver - 1].set(-MOMCD * 0.5 * mddudp)
        kb = pver - 1
        mududp_b = mu[:, kb] * (wm[:, kb] - wm[:, kb - 1]) \
            / dp[:, kb - 1]
        mddudp_b = md[:, kb] * (wm[:, kb] - wm[:, kb - 1]) \
            / dp[:, kb - 1]
        pgu = pgu.at[:, kb].set(-MOMCU * mududp_b)
        pgd = pgd.at[:, kb].set(-MOMCD * mddudp_b)

        # in-cloud velocities: boundary levels (Fortran k=2, kk=pver)
        mk = mupdudp[:, pver - 1] > MBSTH
        wiu = wiu.at[:, pver - 1].set(jnp.where(
            mk,
            (eu[:, pver - 1] * wm[:, pver - 1] * dp[:, pver - 1]
             + pgu[:, pver - 1] * dp[:, pver - 1])
            / jnp.where(mk, mupdudp[:, pver - 1], 1.0),
            wiu[:, pver - 1]))
        mk = md[:, 1] < -MBSTH
        # verbatim Fortran parenthesization: only the pgd term is
        # divided by md
        wid = wid.at[:, 1].set(jnp.where(
            mk,
            (-ed[:, 0] * wm[:, 0] * dp[:, 0])
            - pgd[:, 0] * dp[:, 0] / jnp.where(mk, md[:, 1], 1.0),
            wid[:, 1]))

        # updraft from bottom to top
        for kk in range(pver - 2, -1, -1):
            mk = mupdudp[:, kk] > MBSTH
            num = (mu[:, kk + 1] * wiu[:, kk + 1]
                   + eu[:, kk] * wm[:, kk] * dp[:, kk]
                   + pgu[:, kk] * dp[:, kk])
            wiu = wiu.at[:, kk].set(jnp.where(
                mk, num / jnp.where(mk, mupdudp[:, kk], 1.0),
                wiu[:, kk]))

        # downdraft from top to bottom
        for k in range(2, pver):
            mk = md[:, k] < -MBSTH
            num = (md[:, k - 1] * wid[:, k - 1]
                   - ed[:, k - 1] * wm[:, k - 1] * dp[:, k - 1]
                   - pgd[:, k - 1] * dp[:, k - 1])
            wid = wid.at[:, k].set(jnp.where(
                mk, num / jnp.where(mk, md[:, k], 1.0), wid[:, k]))

        # momentum tendency, k = ktm..pver-1
        ks = np.arange(ktm, pver)
        kp1 = np.minimum(ks + 1, pver - 1)
        wtt = jnp.zeros((n, pver))
        wtt = wtt.at[:, ktm:].set(
            (mu[:, kp1] * (wiu[:, kp1] - wi[:, kp1])
             - mu[:, ks] * (wiu[:, ks] - wi[:, ks])
             + md[:, kp1] * (wid[:, kp1] - wi[:, kp1])
             - md[:, ks] * (wid[:, ks] - wi[:, ks])) / dp[:, ks])

        # cloud-base layer
        kb2 = np.arange(kbm, pver)
        at_mx = kb2[None, :] == mx[:, None]
        wtt = wtt.at[:, kbm:].set(jnp.where(
            at_mx,
            (-mu[:, kb2] * (wiu[:, kb2] - wi[:, kb2])
             - md[:, kb2] * (wid[:, kb2] - wi[:, kb2]))
            * (1.0 / dp[:, kb2]),
            wtt[:, kbm:]))

        # scatter
        wind_tend = wind_tend.at[ideep, :, m].set(wtt)
        pguall = pguall.at[ideep, :, m].set(-pgu)
        pgdall = pgdall.at[ideep, :, m].set(-pgd)
        icwu = icwu.at[ideep, :, m].set(wiu)
        icwd = icwd.at[ideep, :, m].set(wid)

        # momentum flux and end-of-step winds (k = ktm..pver-1)
        mflux = mflux.at[:, ktm:pver, m].set(
            -mu[:, ks] * (wiu[:, ks] - wi[:, ks])
            - md[:, ks] * (wid[:, ks] - wi[:, ks]))
        windf = windf.at[:, ktm:, m].set(
            wm[:, ks] - (mflux[:, ks + 1, m] - mflux[:, ks, m])
            * dt / dp[:, ks])

    # KE dissipation heating (Boville & Bretherton 2003)
    ks = np.arange(ktm, pver)
    km1 = np.maximum(ks - 1, 0)
    kp1 = np.minimum(ks + 1, pver - 1)
    utop = (wind0[:, ks, 0] + wind0[:, km1, 0]) / 2.0
    vtop = (wind0[:, ks, 1] + wind0[:, km1, 1]) / 2.0
    ubot = (wind0[:, kp1, 0] + wind0[:, ks, 0]) / 2.0
    vbot = (wind0[:, kp1, 1] + wind0[:, ks, 1]) / 2.0
    fket = utop * mflux[:, ks, 0] + vtop * mflux[:, ks, 1]
    fkeb = ubot * mflux[:, ks + 1, 0] + vbot * mflux[:, ks + 1, 1]
    ketend_cons = (fket - fkeb) / dp[:, ks]
    ketend = ((windf[:, ks, 0]**2 + windf[:, ks, 1]**2)
              - (wind0[:, ks, 0]**2 + wind0[:, ks, 1]**2)) * 0.5 / dt
    gseten = gseten.at[:, ktm:].set(ketend_cons - ketend)
    seten = seten.at[ideep].set(gseten)

    return dict(wind_tend=wind_tend, pguall=pguall, pgdall=pgdall,
                icwu=icwu, icwd=icwd, seten=seten)
