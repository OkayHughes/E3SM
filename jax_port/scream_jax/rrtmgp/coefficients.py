"""k-distribution coefficient loading and reduction (host-side, numpy).

Sources:
  cpp/examples/mo_load_coefficients.h (Kokkos load_and_init)
  cpp/rrtmgp_conversion.h (conv::SimpleNetCDF::read: netCDF dims are
    REVERSED into the view extents and every int array is shifted to
    0-based)
  cpp/rrtmgp/mo_gas_optics_rrtmgp.h (GasOpticsRRTMGPK::load /
    init_abs_coeffs / reduce_minor_arrays / create_idx_minor /
    create_flavor / create_gpoint_flavor)

The result is a plain dict of numpy arrays ("kdist") with the same
index order the C++ kernels use, e.g. kmajor(gpt, eta, press, temp).
All content indices are 0-based; gpt ranges are inclusive.
"""

import numpy as np


def _read_strings(ds, name):
    """Read a netCDF char array as a list of stripped strings."""
    raw = ds.variables[name][:]
    out = []
    for row in np.asarray(raw):
        s = b"".join(bytes(c) for c in row.filled(b" ") if c is not np.ma.masked) \
            if np.ma.isMaskedArray(row) else b"".join(row)
        out.append(s.decode().strip().lower())
    return out


def _rev(a):
    """Transpose a C-order netCDF array so index order matches the C++
    view (reversed dims), as conv::SimpleNetCDF::read produces."""
    return np.ascontiguousarray(np.transpose(a, tuple(range(a.ndim))[::-1]))


def _string_loc(s, arr):
    """string_loc_in_array: 0-based index or -1."""
    s = s.strip().lower()
    for i, item in enumerate(arr):
        if item.strip().lower() == s:
            return i
    return -1


def _reduce_minor_arrays(available_gases, gas_names, gas_minor,
                         identifier_minor, kminor_atm, minor_gases_atm,
                         minor_limits_gpt_atm, minor_scales_with_density_atm,
                         scaling_gas_atm, scale_by_complement_atm,
                         kminor_start_atm):
    """GasOpticsRRTMGPK::reduce_minor_arrays (host, numpy)."""
    nm = len(minor_gases_atm)
    gas_is_present = np.zeros(nm, dtype=bool)
    tot_g = 0
    for i in range(nm):
        idx_mnr = _string_loc(minor_gases_atm[i], identifier_minor)
        gas_is_present[i] = _string_loc(gas_minor[idx_mnr],
                                        available_gases) >= 0
        if gas_is_present[i]:
            tot_g += minor_limits_gpt_atm[1, i] - minor_limits_gpt_atm[0, i] + 1

    red_nm = int(gas_is_present.sum())
    red = {
        "minor_gases": [], "scaling_gas": [],
        "scales_with_density": np.zeros(red_nm, dtype=bool),
        "scale_by_complement": np.zeros(red_nm, dtype=bool),
        "kminor_start": np.zeros(red_nm, dtype=np.int64),
        "minor_limits_gpt": np.zeros((2, red_nm), dtype=np.int64),
        "kminor": np.zeros((tot_g, kminor_atm.shape[1], kminor_atm.shape[2])),
    }
    if red_nm == nm:
        red["minor_gases"] = list(minor_gases_atm)
        red["scaling_gas"] = list(scaling_gas_atm)
        red["kminor"] = kminor_atm.copy()
        red["minor_limits_gpt"] = minor_limits_gpt_atm.copy()
        red["scales_with_density"] = minor_scales_with_density_atm.copy()
        red["scale_by_complement"] = scale_by_complement_atm.copy()
        red["kminor_start"] = kminor_start_atm.copy()
        return red

    slot = 0
    for i in range(nm):
        if gas_is_present[i]:
            red["minor_gases"].append(minor_gases_atm[i])
            red["scaling_gas"].append(scaling_gas_atm[i])
            red["scales_with_density"][slot] = minor_scales_with_density_atm[i]
            red["scale_by_complement"][slot] = scale_by_complement_atm[i]
            red["kminor_start"][slot] = kminor_start_atm[i]
            slot += 1

    slot = -1
    n_elim = 0
    for i in range(nm):
        ng = minor_limits_gpt_atm[1, i] - minor_limits_gpt_atm[0, i] + 1
        if gas_is_present[i]:
            slot += 1
            red["minor_limits_gpt"][:, slot] = minor_limits_gpt_atm[:, i]
            red["kminor_start"][slot] = kminor_start_atm[i] - n_elim
            s = red["kminor_start"][slot]
            red["kminor"][s:s + ng] = kminor_atm[kminor_start_atm[i]:
                                                 kminor_start_atm[i] + ng]
        else:
            n_elim += ng
    return red


def _create_flavor(key_species_red):
    """GasOpticsRRTMGPK::create_flavor. key_species_red is (2,2,nbnd)."""
    nbnd = key_species_red.shape[2]
    pairs = []
    for ibnd in range(nbnd):
        for iatm in range(2):
            p = [key_species_red[0, iatm, ibnd], key_species_red[1, iatm, ibnd]]
            if p[0] == -1 and p[1] == -1:
                p = [1, 1]
            pairs.append(tuple(p))
    flavors = []
    for p in pairs:
        if p not in flavors:
            flavors.append(p)
    return np.array(flavors, dtype=np.int64).T  # (2, nflav)


def _create_gpoint_flavor(key_species_red, gpt2band, flavor):
    """GasOpticsRRTMGPK::create_gpoint_flavor. Returns (2, ngpt)."""
    ngpt = gpt2band.shape[0]
    gpoint_flavor = np.full((2, ngpt), -1, dtype=np.int64)
    nflav = flavor.shape[1]
    for igpt in range(ngpt):
        for iatm in range(2):
            k1 = key_species_red[0, iatm, gpt2band[igpt]]
            k2 = key_species_red[1, iatm, gpt2band[igpt]]
            if k1 == -1 and k2 == -1:
                k1, k2 = 1, 1
            for iflav in range(nflav):
                if flavor[0, iflav] == k1 and flavor[1, iflav] == k2:
                    gpoint_flavor[iatm, igpt] = iflav
    return gpoint_flavor


def _minor_entries(minor_limits_gpt, kminor_start):
    """Flatten the per-minor-gas gpt ranges into static entry arrays for
    the vectorized minor-species kernel: for each stored absorption
    coefficient row k, which minor gas and which gpt it belongs to."""
    ent_gpt, ent_imnr, ent_k = [], [], []
    for imnr in range(minor_limits_gpt.shape[1]):
        gptS, gptE = minor_limits_gpt[0, imnr], minor_limits_gpt[1, imnr]
        for igpt in range(gptS, gptE + 1):
            ent_gpt.append(igpt)
            ent_imnr.append(imnr)
            ent_k.append(kminor_start[imnr] + (igpt - gptS))
    return (np.array(ent_gpt, dtype=np.int64),
            np.array(ent_imnr, dtype=np.int64),
            np.array(ent_k, dtype=np.int64))


def load_kdist(filename, available_gases):
    """load_and_init + GasOpticsRRTMGPK::load for one k-distribution file.

    available_gases: list of gas names provided by the host model.
    Returns the kdist dict. Longwave files (with totplnk) yield
    'internal' source data (totplnk/planck_frac); shortwave files yield
    'external' (solar_src).
    """
    import netCDF4
    ds = netCDF4.Dataset(filename)
    v = ds.variables

    gas_names_file = _read_strings(ds, "gas_names")
    gas_minor = _read_strings(ds, "gas_minor")
    identifier_minor = _read_strings(ds, "identifier_minor")
    minor_gases_lower = _read_strings(ds, "minor_gases_lower")
    minor_gases_upper = _read_strings(ds, "minor_gases_upper")
    scaling_gas_lower = _read_strings(ds, "scaling_gas_lower")
    scaling_gas_upper = _read_strings(ds, "scaling_gas_upper")

    def iarr(name):  # int arrays: reversed dims, shifted to 0-based
        return _rev(np.asarray(v[name][:], dtype=np.int64)) - 1

    def barr(name):  # bool arrays (stored as int 0/1)
        return _rev(np.asarray(v[name][:], dtype=np.int64)) == 1

    def rarr(name):
        return _rev(np.asarray(v[name][:], dtype=np.float64))

    key_species = iarr("key_species")           # (2, 2, nbnd)
    band2gpt = iarr("bnd_limits_gpt")           # (2, nbnd), 0-based incl.
    band_lims_wvn = rarr("bnd_limits_wavenumber")  # (2, nbnd)
    press_ref = np.asarray(v["press_ref"][:], dtype=np.float64)
    temp_ref = np.asarray(v["temp_ref"][:], dtype=np.float64)
    press_ref_trop = float(np.asarray(v["press_ref_trop"][:]))
    kminor_lower_file = rarr("kminor_lower")    # (contrib, eta, temp)
    kminor_upper_file = rarr("kminor_upper")
    minor_limits_gpt_lower = iarr("minor_limits_gpt_lower")  # (2, nminor)
    minor_limits_gpt_upper = iarr("minor_limits_gpt_upper")
    minor_scales_with_density_lower = barr("minor_scales_with_density_lower")
    minor_scales_with_density_upper = barr("minor_scales_with_density_upper")
    scale_by_complement_lower = barr("scale_by_complement_lower")
    scale_by_complement_upper = barr("scale_by_complement_upper")
    kminor_start_lower = iarr("kminor_start_lower")  # 0-based
    kminor_start_upper = iarr("kminor_start_upper")
    vmr_ref = rarr("vmr_ref")                   # (atmos_layer=2, absorber_ext, temp)
    kmajor = rarr("kmajor")                     # (gpt, eta, press+1, temp)

    has_rayl = "rayl_lower" in v
    rayl_lower = rarr("rayl_lower") if has_rayl else None  # (gpt, eta, temp)
    rayl_upper = rarr("rayl_upper") if has_rayl else None

    kd = {}

    # ---- OpticalProps base init ----
    nband = band2gpt.shape[1]
    ngpt = int(band2gpt.max()) + 1
    gpt2band = np.zeros(ngpt, dtype=np.int64)
    for ibnd in range(nband):
        gpt2band[band2gpt[0, ibnd]:band2gpt[1, ibnd] + 1] = ibnd
    kd["band2gpt"] = band2gpt
    kd["gpt2band"] = gpt2band
    kd["band_lims_wvn"] = band_lims_wvn

    # ---- init_abs_coeffs ----
    avail = [g.strip().lower() for g in available_gases]
    gas_names = [g for g in gas_names_file if g in avail]
    kd["gas_names"] = gas_names
    ngas = len(gas_names)

    # vmr_ref reduction: keep gas 0 (col_dry) plus the present gases
    vmr_ref_red = np.zeros((vmr_ref.shape[0], ngas + 1, vmr_ref.shape[2]))
    vmr_ref_red[:, 0, :] = vmr_ref[:, 0, :]
    for i, g in enumerate(gas_names):
        idx = _string_loc(g, gas_names_file)
        vmr_ref_red[:, i + 1, :] = vmr_ref[:, idx + 1, :]
    kd["vmr_ref"] = vmr_ref_red

    for atm, (kmin, mgas, mlim, mswd, sgas, sbc, kstart) in {
        "lower": (kminor_lower_file, minor_gases_lower,
                  minor_limits_gpt_lower, minor_scales_with_density_lower,
                  scaling_gas_lower, scale_by_complement_lower,
                  kminor_start_lower),
        "upper": (kminor_upper_file, minor_gases_upper,
                  minor_limits_gpt_upper, minor_scales_with_density_upper,
                  scaling_gas_upper, scale_by_complement_upper,
                  kminor_start_upper),
    }.items():
        # presence is checked against the HOST-provided gas list (the C++
        # passes available_gases.gas_name), not the reduced intersection
        red = _reduce_minor_arrays(avail, gas_names_file, gas_minor,
                                   identifier_minor, kmin, mgas, mlim,
                                   mswd, sgas, sbc, kstart)
        kd[f"kminor_{atm}"] = red["kminor"]
        kd[f"minor_limits_gpt_{atm}"] = red["minor_limits_gpt"]
        kd[f"minor_scales_with_density_{atm}"] = red["scales_with_density"]
        kd[f"scale_by_complement_{atm}"] = red["scale_by_complement"]
        kd[f"kminor_start_{atm}"] = red["kminor_start"]
        idx_minor = np.array(
            [_string_loc(gas_minor[_string_loc(m, identifier_minor)],
                         gas_names) for m in red["minor_gases"]],
            dtype=np.int64)
        idx_minor_scaling = np.array(
            [_string_loc(s, gas_names) for s in red["scaling_gas"]],
            dtype=np.int64)
        kd[f"idx_minor_{atm}"] = idx_minor
        kd[f"idx_minor_scaling_{atm}"] = idx_minor_scaling
        (kd[f"minor_entry_gpt_{atm}"], kd[f"minor_entry_imnr_{atm}"],
         kd[f"minor_entry_k_{atm}"]) = _minor_entries(
             red["minor_limits_gpt"], red["kminor_start"])

    kd["press_ref"] = press_ref
    kd["temp_ref"] = temp_ref
    kd["kmajor"] = kmajor
    if has_rayl:
        kd["krayl"] = np.stack([rayl_lower, rayl_upper], axis=-1)
    else:
        kd["krayl"] = None

    kd["press_ref_log"] = np.log(press_ref)
    kd["press_ref_trop_log"] = float(np.log(press_ref_trop))

    # key species reduce + flavors
    key_species_red = np.empty_like(key_species)
    for ip in range(key_species.shape[0]):
        for ia in range(key_species.shape[1]):
            for it in range(key_species.shape[2]):
                ks = key_species[ip, ia, it]
                key_species_red[ip, ia, it] = (
                    -1 if ks == -1 else _string_loc(gas_names_file[ks],
                                                    gas_names))
    flavor = _create_flavor(key_species_red)
    kd["flavor"] = flavor
    kd["gpoint_flavor"] = _create_gpoint_flavor(key_species_red, gpt2band,
                                                flavor)

    kd["temp_ref_min"] = float(temp_ref[0])
    kd["temp_ref_max"] = float(temp_ref[-1])
    kd["press_ref_min"] = float(press_ref[-1])
    kd["press_ref_max"] = float(press_ref[0])
    kd["press_ref_log_delta"] = float(
        (np.log(kd["press_ref_min"]) - np.log(kd["press_ref_max"]))
        / (press_ref.shape[0] - 1))
    kd["temp_ref_delta"] = float(
        (kd["temp_ref_max"] - kd["temp_ref_min"]) / (temp_ref.shape[0] - 1))

    is_key = np.zeros(ngas, dtype=bool)
    for j in range(flavor.shape[1]):
        for i in range(2):
            if flavor[i, j] != -1:
                is_key[flavor[i, j]] = True
    kd["is_key"] = is_key

    kd["idx_h2o"] = _string_loc("h2o", gas_names)

    # ---- source data ----
    if "totplnk" in v:
        kd["totplnk"] = rarr("totplnk")           # (nPlanckTemp, nbnd)
        kd["planck_frac"] = rarr("plank_fraction")  # (gpt, eta, press+1, temp)
        kd["totplnk_delta"] = float(
            (kd["temp_ref_max"] - kd["temp_ref_min"])
            / (kd["totplnk"].shape[0] - 1))
        kd["solar_src"] = None
    else:
        if "solar_source" in v:
            kd["solar_src"] = np.asarray(v["solar_source"][:],
                                         dtype=np.float64)
        else:
            kd["solar_src"] = np.asarray(v["solar_source_quiet"][:],
                                         dtype=np.float64)
        kd["totplnk"] = None
        kd["planck_frac"] = None
        kd["totplnk_delta"] = 0.0

    # sizes
    kd["ngas"] = ngas
    kd["nflav"] = flavor.shape[1]
    kd["neta"] = kmajor.shape[1]
    kd["npres"] = kmajor.shape[2] - 1
    kd["ntemp"] = kmajor.shape[3]
    kd["ngpt"] = ngpt
    kd["nband"] = nband
    ds.close()
    return kd
