# SCREAM Physics in `scream_jax`: Mathematical Formulation and Numerical Machinery

This document explains, for each parameterization in the ported SCREAMv1
physics suite, (1) the general continuous process being modeled, (2) the
closure assumptions adopted, (3) the resulting discrete equations, and
(4) how they are solved — followed by how the schemes are coupled. A
final annotation per scheme records its differentiability structure,
since this port treats that as a first-class property (see STATUS.md,
"Automatic differentiation hardening").

The four-step pattern fits the process models (SHOC, P3, RRTMGP)
naturally. Two components are *not* process models — cld_fraction is a
diagnostic map and SPA is a data interpolator — and are described in
their own terms. Sources of truth are the C++ implementations under
`components/eamxx/src/physics/` at the commit recorded in STATUS.md;
this document describes what the code *does*, in the notation of the
underlying literature where one exists.

Conventions used throughout: columns are independent (all horizontal
coupling happens in the dynamical core); vertical index k = 0 is the
model top; `dp` (pseudo_density) is the wet hydrostatic layer thickness
in Pa; mixing ratios are wet unless suffixed `_dry`; Exner function
Π = (p/p₀)^κ links temperature T and potential temperature θ = T/Π.

---

## 1. The host frame: operator-split column physics

The suite advances the column state
S = (T, qv, qc, nc, qr, nr, qi, qm, ni, bm, u, v, tke, …) by
first-order sequential operator splitting, mirroring EAMxx's
AtmosphereDriver group:

```
for n in 1..6:                       # mac_mic subcycle, Δt/6
    S ← SHOC(S, Δt/6)                # turbulence & mixing
    S ← CldFraction(S)               # diagnostic
    S ← SPA(t_end_of_substep)        # prescribed aerosol
    S ← P3(S, Δt/6)                  # stratiform microphysics
S ← RRTMGP(S, Δt, t_start_of_step)   # radiation, once per step
```

Each process sees the state left by its predecessor ("time-split");
no attempt is made at higher-order splitting. Two pieces of process
memory persist across steps: P3 keeps (qv, T) from the previous
micro-step for its supersaturation relaxation, and RRTMGP keeps
`rad_heating_pdel` so its heating can be re-applied on non-radiation
steps (the port calls radiation every step). SPA is evaluated at the
end-of-substep timestamp, RRTMGP at the start-of-step timestamp —
these conventions are load-bearing for reproducing EAMxx.

---

## 2. SHOC — Simplified Higher-Order Closure (turbulence, shallow clouds)

### 2.1 Governing process

Reynolds-averaged anelastic boundary-layer turbulence. Decomposing
fields into mean and fluctuation (φ = φ̄ + φ′), the exact moment
hierarchy for the conserved thermodynamic pair — liquid water potential
temperature θl and total water qw — couples each moment to the next:

∂ₜ tke = P + B − ε − ∂z(w′e)  with
P = −(u′w′ ∂z ū + v′w′ ∂z v̄)  (shear production),
B = (g/θv) w′θv′  (buoyancy production),
ε (dissipation),

plus prognostic/diagnostic equations for the second moments
(w′θl′, w′qw′, θl′², qw′², θl′qw′, w′²) and the third moment w′³.
The closure problem: B and the cloud quantities depend on the *joint
distribution* of (w, θl, qw) within the grid cell, not just its
moments.

### 2.2 Assumptions

1. **1.5-order closure**: only tke is prognostic; all second moments
   are diagnosed via downgradient forms with eddy diffusivities
   K_m = c_k ℓ √tke (momentum/tke) and K_h (heat/moisture), where the
   **mixing length** ℓ blends an integral (Blackadar-type) scale
   L_inf computed from ∫ z √tke dz / ∫ √tke dz, a stability
   (Brunt–Väisälä) limitation ∝ √tke/N, and a near-surface vK-scaling,
   with tunable length_fac and hard clip bounds.
2. **Assumed joint PDF**: within each cell, (w, θl, qw) follows a
   **double-Gaussian mixture**: two plumes with weights (a, 1−a) set by
   the skewness of w; per-plume scalar means/variances follow from
   linear regression of scalars on w (correlations w-θl, w-qw from the
   diagnosed fluxes). This converts the unclosed buoyancy and cloud
   terms into *analytic Gaussian integrals*.
3. **Subgrid condensation** within each plume: define the saturation
   excess s = (qw − q_sat(T_l, p))/(1 + β) with the psychrometric
   factor β = (L/cp)(∂q_sat/∂T); s is Gaussian within a plume, so
   cloud fraction C = ½(1 + erf(s̄/√2 σ_s)) and liquid
   ql = s̄·C + σ_s·φ(s̄/σ_s) — closed forms per plume, mixture-summed.
   The buoyancy flux w′θv′ follows analytically from the same PDF
   (including the latent contribution through s).
4. **Dissipation** ε = tke^{3/2}/(c_ε ℓ) with stability-dependent
   coefficients; isotropization time scale τ = 2ℓ/√tke feeding the
   return-to-isotropy terms.
5. Implicit (backward-Euler) treatment of vertical diffusion for
   unconditional stability; surface fluxes enter as lower boundary
   conditions of the implicit solve.

### 2.3 Resulting equations

Per substep of length Δt_s ≤ 300 s:

- tke update: explicit source S_e = P + B − ε (with P from diagnosed
  shear, B from the PDF), then implicit diffusion:
  (I − Δt_s ∂z ρK_m ∂z / ρ) tke^{n+1} = tke^n + Δt_s S_e, clipped to
  [mintke, maxtke].
- Prognostic thermodynamics/momentum: the same implicit diffusion
  operator applied to (u, v, θl, qw, tracers) with K_m/K_h and
  surface-flux BCs; afterwards (θl, qw) are inverted back to (T, qv,
  qc) through the PDF's condensation diagnosis.
- Diagnostic chain each substep: grid/thermo setup → Obukhov length
  and surface scales → tke advance → mixing length → assumed-PDF
  (cloud fraction, ql, w′θv′, higher moments) → second/third moments →
  PBL height (first Richardson-number crossing, interpolated).
- A column **energy fixer** enforces conservation of the column moist
  static energy integral against the implicit solve's discretization
  by a uniform additive temperature correction below the PBL.

### 2.4 Solution machinery

The implicit diffusion is a per-column **tridiagonal system**
(Thomas algorithm, `foundation/tridiag.py`, ported from ekat's
`thomas_a1x1`), one solve per prognostic bundle. Everything else is
explicit pointwise algebra plus vertical cumulative sums (length-scale
integrals, PBL search). The scheme subcycles internally so that its
effective time step respects its own stability/accuracy limit
(6 subcycles at Δt = 1800 s).

### 2.5 Differentiability

Smooth *by physical design*: the assumed-PDF closure already is the
"approximation by identity" of the binary saturation decision — cloud
quantities are erf/Gaussian integrals, not thresholds. Remaining
non-smoothness is kinks (clips, mintke floor, lambda ramps) plus the
PBL-height first-crossing selection (can jump layers for non-monotone
Ri). Reverse-mode gradients validate against forward/FD at ~2.5e-4.

---

## 3. Cloud fraction (diagnostic; the pattern degenerates)

There is no governing process: EAMxx's cld_fraction is the map

C_ice = 𝟙[qi > τ],  C_tot = max(C_ice, C_liq),  τ = 1e-12 (internal)
and a looser τ = 1e-5 for output diagnostics; C_liq comes from SHOC's
PDF. The implicit assumption is binary in-cloud homogeneity for ice.
The port adds an opt-in smooth surrogate (sigmoid in (qi−τ)/τ,
smooth-max), which is exactly the statement "qi has a subgrid
distribution of relative width w about its mean" — the width→0 limit
recovers the binary scheme bitwise.

---

## 4. SPA — Simple Prescribed Aerosol (data interpolator; pattern degenerates)

No dynamics: SPA supplies CCN number (nccn → P3) and aerosol optical
properties (τ, ω̃, g per SW band; τ per LW band → RRTMGP) from a
monthly climatology file. Machinery:

1. **Time interpolation**: yearly-periodic linear interpolation in
   fractional day-of-year between bracketing monthly slices
   (end-of-substep timestamp).
2. **Vertical remap**: the file's fields live on hybrid levels
   p_src = PS_src·hybm + P₀·hyam; they are remapped to the model's
   p_mid by piecewise-linear interpolation (ekat LinInterp semantics:
   `searchsorted` bracket + lerp) with constant extrapolation beyond
   the source range (boundary clamps).
3. **Repair clamps** mirroring EAMxx's "repairable postcondition"
   checks (e.g. g ≤ 1 enforced; the data file genuinely contains
   g > 1 points).

Differentiability: value-continuous piecewise-linear in p_mid — kinks
only.

---

## 5. P3 — Predicted Particle Properties (stratiform microphysics)

### 5.1 Governing process

The kinetic (Smoluchowski-type) equation for each hydrometeor
population's size distribution N(D, x, t):

∂ₜN + ∂z(V(D) N) = (∂ₜN)_growth + (∂ₜN)_collection + (∂ₜN)_phase + (∂ₜN)_nucleation

— gravitational sedimentation with size-dependent fall speed V(D),
diffusional growth/evaporation (supersaturation-driven), pairwise
collection (accretion, self-collection, riming), phase changes
(freezing, melting, deposition/sublimation), and nucleation/activation
sources, coupled to the thermodynamic state through latent heating and
vapor exchange.

### 5.2 Assumptions

1. **Bulk two-moment closure**: prognose only mass and number mixing
   ratios per species — cloud liquid (qc, nc), rain (qr, nr), and a
   **single ice category** (qi, ni) carrying two extra prognostics:
   rime mass qm and rime volume bm. This is P3's defining idea:
   instead of discrete ice classes (snow/graupel/hail), the *predicted
   properties* rime fraction F_rim = qm/qi and rime density
   ρ_rim = qm/bm vary continuously, and all ice process rates and fall
   speeds depend on (qi, ni, F_rim, ρ_rim).
2. **Gamma size distributions**: N(D) = N₀ D^μ e^{−λD} per species;
   (N₀, λ) from the two moments with diagnostic shape μ: for rain a
   fixed/diagnosed μ_r, for cloud μ_c(λ) from the observed
   dispersion–λ relation, for ice μ_i(λ) fits. Limiters keep λ within
   physical mean-size bounds (re-deriving number when clipped).
3. **Process rates as distribution integrals**: warm-phase rates use
   closed-form power laws (Khairoutdinov–Kogan-style autoconversion
   ∝ qc^2.47 nc^−1.79, accretion ∝ (qc qr)^1.15, self-collection,
   breakup switching on mean raindrop size); ice rates and fall speeds
   use **precomputed lookup tables**: integrals of the gamma
   distribution against mass–size, area–size and fall-speed laws that
   depend on (normalized qi, F_rim, ρ_rim) — 3-parameter tables for
   ice properties and a 2D table pair for rain self/ice-rain
   collection, interpolated multilinearly at run time.
4. **Saturation**: Flatau polynomial fits for e_sat over liquid/ice
   (Murphy–Koop available); condensation/evaporation as relaxation of
   supersaturation with psychrometric correction, using the (qv, T)
   from the previous micro-step to build the relaxation target (the
   `qv_prev/T_prev` memory).
5. **Subgrid partition**: processes act on **in-cloud** values
   q_incld = q/C with the appropriate cloud fraction (liquid, ice,
   rain — rain fraction from max overlap of the layers above), then
   tendencies are mapped back to cell means; mincld floors avoid
   division blowups.
6. **Conservation limiting**: within a substep, if the sum of sinks of
   a species exceeds its available mass, all sinks are scaled by the
   ratio (a positional projection onto the non-negativity constraint —
   continuous in the inputs).
7. Homogeneous freezing of all liquid at T < 233.15 K; Cooper curve
   ice nucleation (with activation cap); Hallett–Mossop rime
   splintering; wet-growth switching of riming collection; melting and
   shedding above 0°C.

### 5.3 Resulting equations

Per micro substep (Δt/6), the sequence is:

- **Part 1**: thermodynamic setup (Exner, ρ, saturation), in-cloud
  values, initial DSD parameters.
- **Part 2**: the process-rate network evaluated explicitly — a
  directed graph of ~30 rates (autoconversion, accretion,
  self-collection, rain evap, deposition/sublimation, riming of
  cloud/rain, melting, shedding, nucleation, Bergeron, …) — assembled
  into tendencies for the 9 prognostics + θ, with the conservation
  scalings applied per species, iterated over the substep.
- **Sedimentation** (per species): the flux-form conservative upwind
  discretization of ∂ₜq = −∂z(ρ V_q q)/ρ with distinct
  number/mass-weighted fall speeds (V_n, V_m) from the DSD tables.
- **Homogeneous freezing**, then **Part 3**: final DSD parameters,
  effective radii (from μ, λ), precip diagnostics.

### 5.4 Solution machinery

Everything is explicit; the stiffness is handled by (a) the mac_mic
subcycle and (b) **adaptive CFL substepping inside sedimentation**:
per column, n_step = ⌊Co_max⌋+1 substeps of Δt_left/n_step over a
moving "active band" [k_top, k_bot] that grows downward as the
precipitation front advances, accumulating surface flux when the band
reaches the ground. Table lookups are multilinear interpolation with
integer index snapping. The port provides three bit-compatible or
envelope-compatible realizations of the sedimentation loop (adaptive
`while`, masked fixed-length `scan`, and a uniform-substep smooth mode
— see §8 and sedimentation.py).

### 5.5 Differentiability

The jump structure was measured (harness/p3_jump_gap.py): point-gate
jumps (freezing/melting-layer switches — smoothed opt-in via
`smooth_width`/`smooth_families`) and *controller* discreteness
(adaptive substep counts, band indices, column early-exits) which
dominates the precipitation sensitivity and motivates the uniform
integrator mode. Details and measurements in STATUS.md and the commit
history.

---

## 6. RRTMGP — radiation

### 6.1 Governing process

The plane-parallel, azimuthally-averaged radiative transfer equation
per monochromatic frequency ν:

μ dI_ν/dτ_ν = I_ν − S_ν,   S_ν = (1−ω̃)B_ν(T) + (ω̃/2)∫ P(μ,μ′) I dμ′

with optical depth dτ_ν = (k_gas + k_cloud + k_aerosol + k_Rayleigh) ρ dz,
solar beam boundary condition at TOA (SW) and thermal emission with
surface BCs (LW). The frequency integral over ~10⁵ lines is the
central difficulty.

### 6.2 Assumptions

1. **Correlated-k distribution**: within each of 14 SW / 16 LW bands,
   reorder absorption coefficients by magnitude and integrate over the
   cumulative distribution with a fixed Gaussian quadrature — 112 SW /
   128 LW g-points. k(g; T, p, η) is tabulated; η is the binary
   mixing fraction of each band's two major species, interpolated with
   2×2×2 (temperature × pressure × η) stencils, plus "minor" species
   contributions per atmospheric layer regime and Rayleigh scattering
   (SW). Planck source fractions per g-point are tabulated likewise.
2. **Cloud optics**: bulk lookup tables in effective radius (liquid:
   gamma distribution assumption; ice: roughened aggregate habit),
   giving (τ, ω̃, g) per band from water paths (LWP/IWP computed from
   qc, qi and cloud fraction); **delta-scaling** removes the forward
   peak: τ′ = (1−ω̃g²)τ, etc.
3. **Subgrid cloud overlap — the stochastic element**: cloud
   configuration within a column is a random variable with
   maximum-random overlap statistics. **MCICA** samples one binary
   subcolumn per g-point (JSF64 counter-based PRNG, seeds from the
   state) and the g-point quadrature doubles as the Monte-Carlo
   average — an unbiased estimator of the overlap expectation of the
   flux with O(1) cost.
4. **Two-stream closure (SW)**: Meador–Weaver/PIFM coefficients give
   per-layer reflectance/transmittance (direct and diffuse); layers
   are combined by the **adding method** (interaction principle)
   recursions.
5. **No-scatter emission integral (LW)**: ω̃ ≈ 0 in the thermal
   bands; a single secant μ̄ = 1/1.66 quadrature of the Schwarzschild
   integral with linear-in-τ Planck source within each layer.
6. Gas concentrations: qv from state, prescribed O₃, well-mixed
   CO₂/CH₄/N₂O/CFCs with tropopause-shaped `trcmix` profiles; solar
   geometry from Berger orbital series (host-side, calendar-only).

### 6.3 Resulting equations and solution

Per column: build τ_gas(g, k) by table interpolation → add aerosol
(clear-sky) and sampled cloud (all-sky) optics → SW: two-stream
closed forms per layer + adding sweeps (top-down/bottom-up scans) for
direct/diffuse fluxes with band-partitioned surface albedos (VIS/NIR
split at 0.7 µm) → LW: downward then upward Schwarzschild recursions
with surface emission given by the prescribed upward flux → band/g
sums → heating rate Q = (g/cp) ∂F_net/∂p, stored as Q·dp for re-use.
All recursions are linear scans in k; there are no iterative solves.

### 6.4 Differentiability

Traceable end-to-end (pure jnp on the state path); gradients through
the Planck/gas path are validated (FD 2e-5). The MCICA masks are
binary samples: gradients w.r.t. cloud fraction *through the sampling*
are structurally zero — the roadmap item is an expected-overlap mode
that evaluates the overlap expectation deterministically.

---

## 7. Coupling

### 7.1 Within the suite

The driver (scream_jax/driver.py) is a faithful transcription of the
EAMxx group semantics of §1, exchanging fields through a flat dict
keyed by EAMxx field names. Conversions at process boundaries follow
EAMxx's pre/post steps: wet ↔ dry mixing ratios via
q_wet = q_dry·dp_dry/dp_wet with dp_wet = dp_dry(1+Σq_dry); pressure
levels from cumulative dp with hybrid interfaces. cld_fraction's
output feeds both SPA-era P3 (cldfrac_tot for the in-cloud partition)
and RRTMGP (cloud optics weighting); SPA's nccn enters P3 as the
prescribed-CCN source; P3's effective radii and cloud masses feed
RRTMGP's cloud optics; RRTMGP's heating persists via rad_heating_pdel.

### 7.2 To the dynamical core (pySEs bridge)

pySEs advances (horizontal_wind [physical (u,v)], thermodynamic
variable, dry mass d_mass, dry tracers) on GLL points. The bridge
(scream_jax/pyses_bridge.py):

- builds the EAMxx state (wet dp, wet q, p levels from ptop + Σdp) by
  the identities above,
- runs the suite,
- returns physics_forcing tendencies FT = ΔT/Δt, FU/FV, and dry
  tracer tendencies FQ computed with the *updated* wet mass implied by
  the new water (dp_wet' = dp_dry/(1−Σq_wet')), which pySEs applies in
  its coupling step,
- persists the physics-internal prognostics (tke, second moments,
  hydrometeor numbers, precip accumulators, the qv/T memory) that the
  dycore does not carry.

**Physics-grid placement** is an option: np4 (physics on GLL columns
directly) or pg2 via the operational HOMME `gllfvremap` machinery
(pyses_ext/finite_volume_grid_operational.py): GLL→FV by exact
subcell-integral averaging; FV→GLL by the mass-constrained L2
projection (identical to the Hannah et al. reference operator for
nf ≥ 2); temperature remapped in θ-form (T·(p₀/p)^κ, dp-weighted)
with FV pressure reconstructed hydrostatically from remapped ps;
winds through reference-element map tensors; tracers as *states* with
the CAAS shape-preserving limiter — itself a projection onto
{box bounds} ∩ {mass hyperplane}, i.e. a singly-constrained bounded
QP. T/wind updates return as tendencies, tracers as limited states
(HOMME's asymmetric convention), pre-DSS.

### 7.3 Surface

Surface fluxes (sensible, evaporation, momentum), upward LW flux and
albedos are prescribed inputs in this port (no surface model);
they enter SHOC's implicit-solve BCs and RRTMGP's boundary
conditions respectively.

---

## 8. Cross-cutting numerical machinery (index)

| Machinery | Where | Notes |
|---|---|---|
| Tridiagonal (Thomas) solve | SHOC implicit diffusion | `foundation/tridiag.py`; roadmap: `lax.custom_linear_solve` |
| Conservative upwind + adaptive CFL substepping | P3 sedimentation | three realizations: `while` (fast), masked `scan` (AD), uniform (smooth) |
| Multilinear table interpolation | P3 ice/rain tables; RRTMGP k-tables, cloud optics | integer index snap ⇒ kinks |
| Saturation vapor pressure fits | SHOC, P3 | Flatau polynomials (MurphyKoop alt.) |
| erf/Gaussian integrals | SHOC assumed PDF | the model's native smooth closure |
| Ratio conservation limiting | P3 part2 | continuous projection onto positivity |
| Adding-method / Schwarzschild scans | RRTMGP | linear recursions, no iteration |
| Counter-based PRNG subcolumn sampling | RRTMGP MCICA | JSF64, bit-exact vs C++ |
| Periodic/linear interpolation, hybrid-level remap | SPA | searchsorted + lerp + clamps |
| GLL↔FV remap, CAAS projection | pg2 bridge | constrained L2 + clip-and-sum |
| Smoothing operators (approximation by identity) | foundation/smoothing.py | width=0 bitwise-exact; adopted in cld_fraction, P3 jump families |

For validation methodology (golden tiers), see TEST_HARNESS_DESIGN.md;
for per-file provenance and the AD state, STATUS.md.
