"""SPA process step: time interpolation + Dynamic3DRef vertical remap.

Sources:
  components/eamxx/src/physics/spa/eamxx_spa_process_interface.cpp
  components/eamxx/src/share/algorithm/eamxx_data_interpolation.cpp
    (run: fields = alpha*end + (1-alpha)*beg with
     alpha = (ts - t_beg)/interval_len on a YearlyPeriodic timeline;
     p_data = PS*hybm + P0*hyam)
  components/eamxx/src/share/remap/vertical_remapper.cpp + ekat LinInterp
    (per-column piecewise-linear in pressure with edge-slope handling,
     then P0/constant extrapolation strictly outside the source range,
     applied to the bottom half of target levels for p>p_max and the
     top half for p<p_min)

The standalone SPA configuration has no horizontal remap (data on the
model grid) and yearly-periodic time interpolation; only that
configuration is implemented. All arrays are numpy (this process is
data movement, not differentiable physics).
"""

import numpy as np

P0 = 100000.0  # scream::physics::Constants P0
DAYS_PER_YEAR = 365.0

FIELD_MAP = {
    "nccn": "CCN3",
    "aero_g_sw": "AER_G_SW",
    "aero_ssa_sw": "AER_SSA_SW",
    "aero_tau_sw": "AER_TAU_SW",
    "aero_tau_lw": "AER_TAU_LW",
}


def load_spa_data(filename):
    """Load all monthly slices of the SPA fields (+PS, hyam, hybm) into
    memory. Field arrays keep the file layout (time, ncol[, band], lev)."""
    import netCDF4
    ds = netCDF4.Dataset(filename)
    data = {
        "time_doy": np.mod(np.asarray(ds.variables["time"][:],
                                      dtype=np.float64), DAYS_PER_YEAR),
        "PS": np.asarray(ds.variables["PS"][:], dtype=np.float64),
        "hyam": np.asarray(ds.variables["hyam"][:], dtype=np.float64),
        "hybm": np.asarray(ds.variables["hybm"][:], dtype=np.float64),
    }
    for out_name, file_name in FIELD_MAP.items():
        data[out_name] = np.asarray(ds.variables[file_name][:],
                                    dtype=np.float64)
    ds.close()
    return data


def time_interp_coeffs(doy, slice_doys):
    """Yearly-periodic interval [beg, end] containing day-of-year `doy`
    (0-based fractional). Returns (i_beg, i_end, alpha)."""
    n = len(slice_doys)
    # find the slice with the largest doy <= target; wrap if before first
    i_beg = int(np.searchsorted(slice_doys, doy, side="right")) - 1
    if i_beg < 0:
        i_beg = n - 1
    i_end = (i_beg + 1) % n
    t_beg = slice_doys[i_beg]
    t_end = slice_doys[i_end]
    length = (t_end - t_beg) % DAYS_PER_YEAR
    elapsed = (doy - t_beg) % DAYS_PER_YEAR
    return i_beg, i_end, elapsed / length


def _lin_interp_column(x_src, x_tgt, y_src):
    """ekat::LinInterp along the last axis: k1 = last src index with
    x_src[k1] <= x_tgt (clamped), slope from (k1, k1+1) except at the
    last point where it is (k1-1, k1). y_src may have extra leading
    band axes; x arrays are (nlev_src,) / (nlev_tgt,)."""
    km1 = x_src.shape[0]
    k1 = np.searchsorted(x_src, x_tgt, side="right") - 1
    k1 = np.clip(k1, 0, km1 - 1)
    k1ph = np.where(k1 == km1 - 1, k1 - 1, k1 + 1)
    x1 = x_src[k1]
    x1ph = x_src[k1ph]
    y1 = y_src[..., k1]
    y1ph = y_src[..., k1ph]
    return y1 + (y1ph - y1) * (x_tgt - x1) / (x1ph - x1)


def vertical_remap(field_src, p_src, p_tgt):
    """VerticalRemapper (P0 extrapolation both ends). field_src is
    (ncol[, nband], nlev_src); p_src/p_tgt are (ncol, nlev)."""
    ncol = p_src.shape[0]
    nlev_tgt = p_tgt.shape[1]
    mid = nlev_tgt // 2
    out = np.empty(field_src.shape[:-1] + (nlev_tgt,))
    for icol in range(ncol):
        xs = p_src[icol]
        xt = p_tgt[icol]
        y = _lin_interp_column(xs, xt, field_src[icol])
        # P0 extrapolation strictly outside the source range
        bot = (np.arange(nlev_tgt) >= mid) & (xt > xs[-1])
        top = (np.arange(nlev_tgt) < mid) & (xt < xs[0])
        y[..., bot] = field_src[icol][..., -1:]
        y[..., top] = field_src[icol][..., 0:1]
        out[icol] = y
    return out


def spa_process_step(data, doy_end_of_step, p_mid):
    """One SPA step: returns dict with nccn, aero_g_sw, aero_ssa_sw,
    aero_tau_sw, aero_tau_lw on the model grid/levels.

    doy_end_of_step: 0-based fractional day-of-year of the END-of-step
    timestamp (DataInterpolation::run receives end_of_step_ts)."""
    i_beg, i_end, alpha = time_interp_coeffs(doy_end_of_step,
                                             data["time_doy"])
    ps = (1.0 - alpha) * data["PS"][i_beg] + alpha * data["PS"][i_end]
    p_src = ps[:, None] * data["hybm"][None, :] + P0 * data["hyam"][None, :]

    out = {}
    for name in FIELD_MAP:
        f = (1.0 - alpha) * data[name][i_beg] + alpha * data[name][i_end]
        out[name] = vertical_remap(f, p_src, np.asarray(p_mid))

    # the SPA process adds REPAIRABLE FieldWithinIntervalCheck
    # postconditions, which clamp out-of-bounds values (the input files
    # do contain e.g. g > 1 at some points)
    out["nccn"] = np.clip(out["nccn"], 0.0, 1e11)
    for name in ("aero_g_sw", "aero_ssa_sw", "aero_tau_sw", "aero_tau_lw"):
        out[name] = np.clip(out[name], 0.0, 1.0)
    return out
