"""Port of the zm_conv_intr.F90 zm_conv_tend per-timestep sequence as
a pure-JAX driver over the already-ported ZM kernels: the primary
Zhang-McFarlane interface that CAM's convect_deep calls each step.

PORT_NOTES
----------
Sequence, exactly as zm_conv_tend performs it (all zm_param/zm_const
handling and the DCAPE trigger's two compute_dilute_cape calls live
inside eam_jax.zm_conv.zm_conv_main):

  1. zm_conv_main on the incoming state (t, q(:,:,0)=vapor, omega,
     pmid/pint/pdel, phis, zm/zi, pblh, tpert, landfrac, t_star,
     q_star) -> ptend_loc(s=heat, q1=qtnd) + prec, rprd, dlf, plume
     fluxes, cape/dcape, jctop/jcbot.
  2. mcon is converted from mb/s to kg/m2/s over levels 0..pver-1
     ONLY (the intr loop is 1:pver; the last interface entry stays
     unconverted -- verbatim).
  3. physics_update of the local state copy (ptend 'zm_conv_main',
     ls + lq(1)): t1 = t + dt*heat/cpair; q1 = q + dt*qtnd followed
     by the qneg3 clip at qmin(Q) = 1e-12 (physpkg cnst_add('Q', ...,
     1.E-12); use_mass_borrower is .false. by default in
     phys_control, so qneg3 -- an elementwise max against qmin -- is
     the production path). physics_update also refreshes zm/zi/s via
     geopotential_t, but nothing downstream of it in THIS routine
     reads them (evap uses t/q/pmid/pdel; the transports use u/v/q),
     so the refresh is omitted here; chained stepping recomputes
     geopotential between steps (see tests/test_zm_chain.py).
  4. zm_conv_evap on the updated state1 with prdprec = rprd, the
     pbuf cloud fraction, and prec from step 1 -> evap heating /
     moistening, final prec + snow, flxprec/flxsnow (the DP_FLXPRC /
     DP_FLXSNW pbuf outputs), evapcdp = evap tend_q (NEVAPR_DPCU).
  5. physics_update again (ptend 'zm_conv_evap'): t2, q2.
  6. zm_transport_momentum on state1 winds (u/v are untouched so far,
     so winds = the input u/v) with the gathered plume fluxes ->
     wind tendencies + KE-dissipation heating seten; physics_update
     (u/v/t updated; nothing downstream reads them).
  7. zm_transport_tracer ('zm_transport_tracer_1') on state1%q --
     vapor slot updated through q2 but never transported (lq(1) =
     .false.); lq(2:) = cnst_is_convtran1(2:) arrives as the
     doconvtran argument; fake_dpdry = 0 since convtran1 species are
     wet mixing ratios (verbatim).
  8. ptend_all = sum of the four ptend_locs: s = heat + evap_s +
     seten; q(:,:,0) = qtnd + evap_q; q(:,:,m>0) = tracer transport;
     u/v = momentum transport. The CALLER applies ptend_all (with
     name 'zm_conv_tend') to the model state.

SKIPPED (documented scope, matching the kernel ports and goldens):
  - zm_microp=.false. only (EAMv3 default is zmconv_microp=.true.;
    zm_microphysics.F90 is PORTING_PLAN.md row 9). Hence dlftot =
    dlf, rice = 0, old_snow keeps its zm_param_t default .true.
    (zm_conv_intr sets .false. only under zm_microp), and the
    microp_st pbuf plumbing (qi/dif/dsf/dnlf/...) is absent.
  - MCSP: EAMv3 phys="default" sets zmconv_MCSP_heat_coeff = 0.3 in
    namelist_defaults_eam.xml, but the MCSP tendency path is NOT
    ported (zm_conv_mcsp.F90); this driver corresponds to the
    all-zero-coefficient configuration (mcsp_enabled = .false.), the
    same scope as the zm_conv golden.
  - aerosol convective transport (zm_conv_tend_2 with convproc /
    species_class logic) and the aero(lchnk) object: EAMv3 runs them
    via chemistry; out of scope with zm_microp=F.
  - pbuf/cam_history plumbing: pbuf fields become plain array inputs
    (cld, t_star, q_star, fracis) and outputs (flxprec, flxsnow,
    evapcdp, dlf, ql=ICWMRDP, rprd=RPRDDP, prec=PREC_DP, snow=
    SNOW_DP); outfld diagnostics (FREQZM, PCONVT/B, MAXI, mu_out/
    md_out scatter) are omitted.
  - is_first_step handling of T_STAR/Q_STAR pbuf initialization is
    the caller's job (physpkg records them at the end of tphysac;
    zm_conv_tend seeds them with the current state on nstep 0 --
    tests/test_zm_chain.py replicates exactly that bookkeeping).

Namelist-dependent constants applied by the intr layer, EAMv3
phys="default" (namelist_defaults_eam.xml): zmconv_ke = 2.5e-6
(dyn="se" microphys="p3"; single evaporation efficiency -- EAMv3
zm_conv.F90 has no zmconv_ke_lnd) enters via zp["ke"]; ztodt is the
full physics step (1800 s for ne30pg2-class grids); qmin(Q) = 1e-12.
The 40 hPa ZM_upper_limit_pref rule that fixes limcnv is grid
metadata and enters through zp["limcnv"] exactly as in the zm_conv
port. Everything float64, eager-mode (concrete gathers), plain-array
in/out.
"""

import jax.numpy as jnp
import numpy as np

from .zm_conv import zm_conv_evap, zm_conv_main
from .zm_transport import zm_transport_momentum, zm_transport_tracer

__all__ = ["zm_tend", "QMIN_VAPOR"]

QMIN_VAPOR = 1.0e-12   # physpkg cnst_add('Q', ..., 1.E-12_r8, ...)
GRAVIT_MB = 100.0      # mb -> Pa in the mcon conversion (100/gravit)


def zm_tend(t, q, u, v, omega, p_mid, p_int, p_del, geos, z_mid,
            z_int, pbl_hgt, tpert, landfrac, cld, fracis, doconvtran,
            t_star, q_star, ztodt, is_first_step, zc, zp):
    """zm_conv_tend sequence on plain arrays (level 0 = model top).

    t (ncol,pver) [K]; q (ncol,pver,ncnst) [kg/kg] with slot 0 =
    water vapor; u/v (ncol,pver) [m/s]; omega [Pa/s]; p_mid/p_int/
    p_del [Pa]; geos (ncol) [m2/s2]; z_mid/z_int [m above surface];
    pbl_hgt [m]; tpert [K]; landfrac; cld (ncol,pver) cloud fraction
    (pbuf 'CLD'); fracis (ncol,pver,ncnst); doconvtran (ncnst,) bool
    = cnst_is_convtran1 (slot 0 ignored); t_star/q_star DCAPE
    previous-state fields; ztodt [s]; zc/zp from
    zm_conv.make_zm_const / make_zm_param.

    Returns a dict: total tendencies (s_tend [J/kg/s], q_tend
    (ncol,pver,ncnst) [1/s], u_tend, v_tend [m/s2]), the per-process
    pieces (heat, qtnd, evap_s, evap_q, seten, wind_tends), and the
    zm_conv_tend outputs (prec, snow, rliq, cape, dcape, jctop,
    jcbot, mcon [kg/m2/s], dlftot, pflx, zdu, ql, rprd, evapcdp,
    flxprec, flxsnow, ntprprd, ntsnprd, tend_s_snwprd,
    tend_s_snwevmlt, plus the gathered plume arrays mu/eu/du/md/ed/
    dp/jt/maxg/ideep/lengath/dsubcld for zm_conv_tend_2-style reuse).
    """
    t = jnp.asarray(t, dtype=jnp.float64)
    q = jnp.asarray(q, dtype=jnp.float64)
    ncol, pver, ncnst = q.shape
    cpair = zc["cpair"]

    # ---- 1. primary ZM parameterization on the incoming state ----
    main = zm_conv_main(t, q[:, :, 0], omega, p_mid, p_int, p_del,
                        geos, z_mid, z_int, pbl_hgt, tpert, landfrac,
                        t_star, q_star, ztodt, is_first_step, zc, zp)
    heat = main["heat"]
    qtnd = main["qtnd"]
    lengath = main["lengath"]
    ideep = np.asarray(main["gather_index"][:lengath])
    jt_g = np.asarray(main["jt"][:lengath])
    mx_g = np.asarray(main["msemax_klev_g"][:lengath])
    dlftot = main["dlf"]                      # zm_microp=F: dlf only

    # ---- 2. mcon mb/s -> kg/m2/s over midlevel rows only (verbatim:
    # the intr loop is 1:pver, so the bottom interface entry keeps
    # its mb/s value) ----
    mcon = main["mcon"]
    mcon = mcon.at[:, :pver].set(mcon[:, :pver] * GRAVIT_MB
                                 / zc["grav"])

    # ---- 3. physics_update for ptend 'zm_conv_main' ----
    t1 = t + ztodt * heat / cpair
    q1 = jnp.maximum(q[:, :, 0] + ztodt * qtnd, QMIN_VAPOR)  # qneg3

    # ---- 4. convective rain evaporation on state1 ----
    ev = zm_conv_evap(p_mid, p_del, t1, q1, main["rprd"], cld,
                      main["prec"], ztodt, zc, zp)
    evap_s = ev["tend_s"]
    evap_q = ev["tend_q"]

    # ---- 5. physics_update for ptend 'zm_conv_evap' ----
    # (t2 has no consumer inside this routine; q2 feeds the tracer
    # transport's state1%q vapor slot)
    q2 = jnp.maximum(q1 + ztodt * evap_q, QMIN_VAPOR)

    # ---- 6. convective momentum transport (u/v untouched so far) ---
    winds = jnp.stack([jnp.asarray(u, dtype=jnp.float64),
                       jnp.asarray(v, dtype=jnp.float64)], axis=2)
    mom = zm_transport_momentum(winds, main["mflx_up"][:lengath],
                                main["mflx_dn"][:lengath],
                                main["detr_up"][:lengath],
                                main["entr_up"][:lengath],
                                main["entr_dn"][:lengath],
                                main["p_del"][:lengath],
                                jt_g, mx_g, ideep, ztodt)
    seten = mom["seten"]

    # ---- 7. convective tracer transport (zm_transport_tracer_1) ----
    q_state1 = q.at[:, :, 0].set(q2)
    fake_dpdry = jnp.zeros((max(lengath, 1), pver))
    doconv = np.asarray(doconvtran, dtype=bool).copy()
    doconv[0] = False                          # vapor never here
    dqdt = zm_transport_tracer(
        q_state1, fracis, doconv, np.zeros(ncnst, dtype=bool),
        main["mflx_up"][:lengath], main["mflx_dn"][:lengath],
        main["detr_up"][:lengath], main["entr_up"][:lengath],
        main["entr_dn"][:lengath], main["p_del"][:lengath],
        fake_dpdry, jt_g, mx_g, ideep, ztodt,
        zm_microp=zp["zm_microp"])

    # ---- 8. ptend_all ----
    s_tend = heat + evap_s + seten
    q_tend = dqdt.at[:, :, 0].set(qtnd + evap_q)
    u_tend = mom["wind_tend"][:, :, 0]
    v_tend = mom["wind_tend"][:, :, 1]

    return dict(
        s_tend=s_tend, q_tend=q_tend, u_tend=u_tend, v_tend=v_tend,
        # per-process pieces
        heat=heat, qtnd=qtnd, evap_s=evap_s, evap_q=evap_q,
        seten=seten, wind_tends=mom["wind_tend"],
        pguall=mom["pguall"], pgdall=mom["pgdall"],
        icwu=mom["icwu"], icwd=mom["icwd"],
        # zm_conv_tend outputs / pbuf fields
        prec=ev["prec"], snow=ev["snow"], rliq=main["rliq"],
        rice=jnp.zeros(ncol), cape=main["cape"], dcape=main["dcape"],
        jctop=main["jctop"], jcbot=main["jcbot"], mcon=mcon,
        dlftot=dlftot, pflx=main["pflx"], zdu=main["zdu"],
        ql=main["ql"], rprd=main["rprd"], evapcdp=evap_q,
        flxprec=ev["flxprec"], flxsnow=ev["flxsnow"],
        ntprprd=ev["ntprprd"], ntsnprd=ev["ntsnprd"],
        tend_s_snwprd=ev["tend_s_snwprd"],
        tend_s_snwevmlt=ev["tend_s_snwevmlt"],
        # gathered plume arrays (zm_conv_tend also returns these for
        # the later zm_conv_tend_2 call)
        mu=main["mflx_up"], eu=main["entr_up"], du=main["detr_up"],
        md=main["mflx_dn"], ed=main["entr_dn"], dp=main["p_del"],
        dsubcld=main["dsubcld"], jt=main["jt"],
        maxg=main["msemax_klev_g"], ideep=main["gather_index"],
        lengath=lengath,
        # intermediate states (chain/debug)
        t1=t1, q1=q1, q2=q2)
