"""MCICA subcolumn sampling: JSF64 PRNG, subcolumn cloud masks and
band->gpoint cloud subsampling.

Sources: cpp/rrtmgp_conversion.h (conv::Random — Bob Jenkins' small fast
64-bit PRNG: a=0xf1ea5eed, b=c=d=seed, 20 warm-up rounds, then one
gen() per value) and eamxx_rrtmgp_interface.hpp (get_subcolumn_mask,
get_subsampled_clouds). The C++ seeds a FRESH generator per
(col, lay, gpt) with seed = seeds(icol) + ilay*ngpt + igpt and draws a
single value — fully deterministic and reproduced here exactly with
vectorized uint64 arithmetic (requires jax_enable_x64).

Seeds come from the decimal part of the layer pressure nearest the
surface: SW uses p_lay(:, nlay-1), LW uses p_lay(:, nlay-2), via
seeds = int(1e9 * (p - int(p))).

Differentiability (approximation-by-identity): the subcolumn mask is a
composition of hard ``>`` comparisons on the random deviates, so
d(mask)/d(cldfrac) is identically zero under AD. The static kwarg
``smooth_width`` (a Python float, never traced) relaxes every one of
those comparisons via scream_jax.foundation.smoothing with scale=1.0
(the deviates and cloud fractions live in [0, 1]): width 0 keeps the
original bitwise int mask through an explicit Python-level branch;
width > 0 yields a FRACTIONAL float mask in [0, 1] that multiplies the
cloud optics, converging to the binary scheme as width -> 0.
"""

import jax.numpy as jnp

from ..foundation import smoothing

_U64_MAX_FLOAT = 1.8446744073709552e19  # (double)std::numeric_limits<u64>::max()


def _rot(x, k):
    return (x << jnp.uint64(k)) | (x >> jnp.uint64(64 - k))


def jsf64_value(seed):
    """conv::Random(seed).genFP<double>(): 20 warm-up rounds + 1 draw.

    seed: uint64 array. Returns float64 in [0, 1)."""
    a = jnp.full_like(seed, 0xF1EA5EED, dtype=jnp.uint64)
    b = seed.astype(jnp.uint64)
    c = seed.astype(jnp.uint64)
    d = seed.astype(jnp.uint64)
    for _ in range(21):  # 20 warm-ups + the drawn value
        e = a - _rot(b, 7)
        a = b ^ _rot(c, 13)
        b = c + _rot(d, 37)
        c = d + e
        d = e + a
    return d.astype(jnp.float64) / _U64_MAX_FLOAT


def get_subcolumn_mask(cldf, seeds, ngpt: int, overlap_option: int = 1,
                       smooth_width: float = 0.0):
    """eamxx get_subcolumn_mask (maximum-random overlap, eq. 14 of
    Raisanen et al. 2004).

    cldf: (ncol, nlay) radiative cloud fraction; seeds: (ncol,) int.
    smooth_width (STATIC Python float): 0.0 -> exact bitwise int mask
    (original comparisons); > 0 -> fractional float mask in [0, 1] with
    every binary comparison replaced by smoothing.step/blend at
    scale=1.0 (the natural [0, 1] range of the deviates).

    Returns an int mask (ncol, nlay, ngpt) at width 0, a float
    fractional mask at width > 0."""
    cldf = jnp.asarray(cldf)
    ncol, nlay = cldf.shape
    if overlap_option == 0:
        # no switching variable here: every subcolumn is cloudy
        return jnp.ones((ncol, nlay, ngpt), dtype=jnp.int32)

    lay = jnp.arange(nlay, dtype=jnp.uint64)[None, :, None]
    gpt = jnp.arange(ngpt, dtype=jnp.uint64)[None, None, :]
    seed = (jnp.asarray(seeds, dtype=jnp.uint64)[:, None, None]
            + lay * jnp.uint64(ngpt) + gpt)
    cldx = jsf64_value(seed)  # (ncol, nlay, ngpt)

    # top-down max-random overlap rewrite:
    #   if cldx(l-1) > 1-cldf(l-1): cldx(l) = cldx(l-1)   (max overlap)
    #   else:                        cldx(l) = cldx(l)*(1-cldf(l-1))
    # sequential in l — small nlay, unrolled python loop over layers.
    # Smoothed site 1 (recursion select): switching variable is the
    # "cloudy above" exceedance cldx(l-1) - (1 - cldf(l-1)); at width>0
    # the hard select becomes smoothing.blend of the two branch values.
    cldf_e = cldf[..., None]
    cols = [cldx[:, 0, :]]
    for l in range(1, nlay):
        above = cols[l - 1]
        cf_above = cldf_e[:, l - 1, 0][:, None]
        rand_below = cldx[:, l, :] * (1.0 - cf_above)
        if smooth_width == 0.0:
            cloudy_above = above > 1.0 - cf_above
            cols.append(jnp.where(cloudy_above, above, rand_below))
        else:
            cols.append(smoothing.blend(above - (1.0 - cf_above),
                                        above, rand_below,
                                        width=smooth_width, scale=1.0))
    cldx = jnp.stack(cols, axis=1)

    # Smoothed site 2 (mask threshold): switching variable is the
    # exceedance cldx - (1 - cldf); at width>0 the Heaviside becomes
    # smoothing.step, giving a fractional mask.
    if smooth_width == 0.0:
        return (cldx > 1.0 - cldf_e).astype(jnp.int32)
    return smoothing.step(cldx - (1.0 - cldf_e), smooth_width, scale=1.0)


def compute_seeds(p_lay, nlay_from_bottom: int):
    """seeds(icol) = 1e9 * frac(p_lay(icol, nlay - nlay_from_bottom)).
    SW: nlay_from_bottom=1, LW: 2."""
    p = jnp.asarray(p_lay)[:, -nlay_from_bottom]
    return (1.0e9 * (p - jnp.floor(p))).astype(jnp.int64)
    # NOTE: the C++ writes `1e9 * (p - int(p))` into an int — truncation


def get_subsampled_clouds(cloud_optics_bnd, cldfrac, p_lay, gpt2band,
                          ngpt: int, two_stream: bool,
                          smooth_width: float = 0.0):
    """eamxx get_subsampled_clouds (SW: 2-stream, LW: 1-scalar).

    cloud_optics_bnd: by-band dict from cloud_optics.get_cloud_optics.
    smooth_width (STATIC): 0.0 -> binary mask zeroing (bitwise
    original); > 0 -> the fractional subcolumn mask multiplies the
    by-band optics (the direct relaxation of where(mask==1, x, 0) ==
    mask * x for a binary mask).
    Returns a by-gpoint dict of the same kind."""
    tau_bnd = jnp.asarray(cloud_optics_bnd["tau"])
    cldfrac = jnp.asarray(cldfrac)
    gpt2band = jnp.asarray(gpt2band)

    # "radiative cloud fraction": zero wherever the cloud has no optical
    # properties in any band. Kept hard at all widths: the switching
    # variable is tau presence (a structural gate on the optics tables),
    # not the cloud fraction; jnp.where passes d/d(cldfrac) through the
    # selected branch, so no smoothing is needed for the gradient goal.
    cldfrac_rad = jnp.where(jnp.any(tau_bnd > 0.0, axis=-1), cldfrac, 0.0)

    seeds = compute_seeds(p_lay, 1 if two_stream else 2)
    cldmask = get_subcolumn_mask(cldfrac_rad, seeds, ngpt,
                                 smooth_width=smooth_width)

    if smooth_width == 0.0:
        keep = cldmask == 1
        out = {"tau": jnp.where(keep, tau_bnd[..., gpt2band], 0.0)}
        if two_stream:
            out["ssa"] = jnp.where(keep,
                                   jnp.asarray(cloud_optics_bnd["ssa"])[..., gpt2band],
                                   0.0)
            out["g"] = jnp.where(keep,
                                 jnp.asarray(cloud_optics_bnd["g"])[..., gpt2band],
                                 0.0)
    else:
        out = {"tau": cldmask * tau_bnd[..., gpt2band]}
        if two_stream:
            out["ssa"] = cldmask * jnp.asarray(cloud_optics_bnd["ssa"])[..., gpt2band]
            out["g"] = cldmask * jnp.asarray(cloud_optics_bnd["g"])[..., gpt2band]
    return out
