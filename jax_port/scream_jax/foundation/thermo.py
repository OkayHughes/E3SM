"""Common thermodynamic conversions (PhysicsFunctions).

Transcribed from
components/eamxx/src/share/physics/eamxx_common_physics_functions_impl.hpp
(struct scream::PhysicsFunctions; declarations in
eamxx_common_physics_functions.hpp). Only the scalar overloads are ported —
the C++ "team" overloads are per-level Kokkos loops over the same scalar
math, which jnp broadcasting subsumes.

All functions are elementwise over arrays of any shape, except
calculate_z_int / calculate_z_mid which operate along the last (level) axis
with the EAMxx convention k=0 = model top.
"""

import jax.numpy as jnp

from . import constants as c
from .column_ops import column_scan, compute_midpoint_values


def calculate_dx_from_area(area, lat):
    """Grid length [m] from cell area [sr] and latitude [deg] (WGS84 expansion)."""
    lat_in_rad = jnp.asarray(lat) * (c.Pi / 180.0)
    m_per_degree_lat = (c.earth_ellipsoid1
                        - c.earth_ellipsoid2 * jnp.cos(2.0 * lat_in_rad)
                        + c.earth_ellipsoid3 * jnp.cos(4.0 * lat_in_rad))
    return m_per_degree_lat * jnp.sqrt(jnp.asarray(area)) * (180.0 / c.Pi)


def calculate_density(pseudo_density, dz):
    """Air density [kg/m3] = dp / (dz * g)."""
    return jnp.asarray(pseudo_density) / jnp.asarray(dz) / c.gravit


def calculate_vertical_velocity(omega, density):
    """Vertical velocity w [m/s] = -omega / (rho * g)."""
    return -jnp.asarray(omega) / (jnp.asarray(density) * c.gravit)


def exner_function(pressure):
    """Exner function (p/p0)^(Rd/cp)."""
    return jnp.power(jnp.asarray(pressure) / c.P0, c.RD * c.INV_CP)


def calculate_theta_from_T(temperature, pressure):
    """Potential temperature theta = T / exner(p)."""
    return jnp.asarray(temperature) / exner_function(pressure)


def calculate_T_from_theta(theta, pressure):
    """Temperature T = theta * exner(p)."""
    return jnp.asarray(theta) * exner_function(pressure)


def calculate_thetal_from_theta(theta, temperature, qc):
    """Liquid-water potential temperature:
    thetal = theta - (theta/T) * (Lv/cp) * qc."""
    theta = jnp.asarray(theta)
    return theta - (theta / jnp.asarray(temperature)) * (c.LatVap / c.Cpair) * jnp.asarray(qc)


# c1 = -1 + 1/ep_2, shared by the virtual-temperature pair below.
_C1_VIRTUAL = -c.ONE + c.ONE / c.ep_2


def calculate_virtual_temperature(temperature, qv):
    """Virtual temperature Tv = T * (1 + c1*qv), c1 = -1 + 1/ep_2."""
    return jnp.asarray(temperature) * (c.ONE + _C1_VIRTUAL * jnp.asarray(qv))


def calculate_temperature_from_virtual_temperature(T_virtual, qv):
    """Inverse of calculate_virtual_temperature."""
    return jnp.asarray(T_virtual) / (c.ONE + _C1_VIRTUAL * jnp.asarray(qv))


def calculate_dse(temperature, z, surf_geopotential):
    """Dry static energy: dse = cp*T + g*z + phis."""
    return c.CP * jnp.asarray(temperature) + c.gravit * jnp.asarray(z) + surf_geopotential


def calculate_temperature_from_dse(dse, z, surf_geopotential):
    """Inverse of calculate_dse."""
    return (jnp.asarray(dse) - c.gravit * jnp.asarray(z) - surf_geopotential) / c.CP


def calculate_wetmmr_from_drymmr(drymmr, qv_dry):
    """Wet mmr from dry mmr: wet = dry / (1 + qv_dry)."""
    return jnp.asarray(drymmr) / (1 + jnp.asarray(qv_dry))


def calculate_drymmr_from_wetmmr(wetmmr, qv_wet):
    """Dry mmr from wet mmr: dry = wet / (1 - qv_wet)."""
    return jnp.asarray(wetmmr) / (1 - jnp.asarray(qv_wet))


def calculate_wetmmr_from_drymmr_dp_based(drymmr, pseudo_density, pseudo_density_dry):
    """Wet mmr from dry mmr via pseudo-density ratio."""
    return jnp.asarray(drymmr) * jnp.asarray(pseudo_density_dry) / jnp.asarray(pseudo_density)


def calculate_drymmr_from_wetmmr_dp_based(wetmmr, pseudo_density, pseudo_density_dry):
    """Dry mmr from wet mmr via pseudo-density ratio."""
    return jnp.asarray(wetmmr) * jnp.asarray(pseudo_density) / jnp.asarray(pseudo_density_dry)


def calculate_dz(pseudo_density, p_mid, T_mid, qv):
    """Layer thickness [m] from hydrostatic balance: dz = (Rd/g)*dp*Tv/p."""
    T_virtual = calculate_virtual_temperature(T_mid, qv)
    return (c.RD / c.gravit) * jnp.asarray(pseudo_density) * T_virtual / jnp.asarray(p_mid)


def calculate_z_int(dz, z_surf):
    """Interface geometric heights from layer thicknesses (level axis last).

    Bottom interface (k=nlev) is prescribed at z_surf; heights accumulate
    upward (C++ calls column_scan<FromTop=false>).
    """
    return column_scan(dz, z_surf, from_top=False)


def calculate_z_mid(z_int):
    """Midpoint heights = adjacent-interface averages."""
    return compute_midpoint_values(z_int)


def calculate_vmr_from_mmr(gas_mol_weight, qv, mmr):
    """Volume mixing ratio from mass mixing ratio (moist-air correction)."""
    return jnp.asarray(mmr) / (1.0 - jnp.asarray(qv)) * c.MWdry / gas_mol_weight


def calculate_mmr_from_vmr(gas_mol_weight, qv, vmr):
    """Mass mixing ratio from volume mixing ratio (moist-air correction)."""
    return (gas_mol_weight / c.MWdry) * jnp.asarray(vmr) * (1.0 - jnp.asarray(qv))


def calculate_surface_air_T(T_mid_bot, z_mid_bot):
    """Air temperature at the ground from the lowest midpoint, assuming a
    6.5 K/km lapse rate. Only intended for calculate_psl (see C++ comment)."""
    return jnp.asarray(T_mid_bot) + 0.0065 * jnp.asarray(z_mid_bot)


def lapse_T_for_psl(T_ground, phi_ground):
    """Lapse rate and effective ground T for the sea-level-pressure reduction.

    Returns (lapse, T_ground_tmp). Branch structure mirrors the C++ if/elif
    chain exactly (including its acknowledged crudeness for cold cases).
    Note: evaluated with jnp.where, so phi_ground ~ 0 produces inf/nan in the
    *unselected* warm branch; calculate_psl guards this with a safe phi.
    """
    T_ground = jnp.asarray(T_ground)
    phi_ground = jnp.asarray(phi_ground)
    T_sl = T_ground + 0.0065 * phi_ground / c.gravit

    lapse_make_290 = c.gravit / phi_ground * (290.5 - T_ground)

    cond1 = (T_ground <= 290.5) & (T_sl > 290.5)   # cap T_sl at 290.5
    cond2 = (T_ground > 290.5) & (T_sl > 290.5)    # hot: no lapse, smooth T
    cond3 = T_ground < 255.0                       # cold: EAM's crude fix

    lapse = jnp.where(cond1, lapse_make_290,
            jnp.where(cond2, 0.0, 0.0065))
    T_ground_tmp = jnp.where(cond1, T_ground,
                   jnp.where(cond2, 0.5 * (290.5 + T_ground),
                   jnp.where(cond3, 0.5 * (255.0 + T_ground), T_ground)))
    return lapse, T_ground_tmp


def calculate_psl(T_ground, p_ground, phi_ground):
    """Sea-level pressure, EAM-style reduction (see C++ docstring and
    components/eamxx/docs for the derivation)."""
    T_ground = jnp.asarray(T_ground)
    p_ground = jnp.asarray(p_ground)
    phi_ground = jnp.asarray(phi_ground)

    near_sea_level = jnp.abs(phi_ground / c.gravit) < 1e-4
    # Avoid div-by-zero inside the unselected branch of the where below.
    phi_safe = jnp.where(near_sea_level, 1.0, phi_ground)

    lapse, T_ground_tmp = lapse_T_for_psl(T_ground, phi_safe)
    alpha = lapse * c.Rair / c.gravit
    beta = phi_safe / (c.Rair * T_ground_tmp)
    psl_lapse = p_ground * jnp.exp(beta * (1 - alpha * beta / 2
                                           + (alpha * beta) ** 2 / 3))
    return jnp.where(near_sea_level, p_ground, psl_lapse)


def apply_rayleigh_friction(dt, otau, u_wind, v_wind, T_mid):
    """Rayleigh friction with frictional heating; returns (u, v, T) updated.

    Pure-functional version of the C++ in-place update. NOTE: the C++
    applies the EAM tendency coefficients (c1, c3) directly to the state
    without a dt factor — transcribed here verbatim, quirks included. Do
    not "fix" the scheme; validate against golden data.
    """
    otau = jnp.asarray(otau)
    u_wind = jnp.asarray(u_wind)
    v_wind = jnp.asarray(v_wind)
    T_mid = jnp.asarray(T_mid)

    dt_inv = 1.0 / dt
    c2 = 1.0 / (1.0 + otau * dt)
    c1 = -1.0 * otau * c2
    c3 = 0.5 * (1.0 - c2 * c2) * dt_inv

    u2 = u_wind * u_wind
    v2 = v_wind * v_wind

    u_new = u_wind + c1 * u_wind
    v_new = v_wind + c1 * v_wind
    T_new = T_mid + c3 * (u2 + v2) / c.CP
    return u_new, v_new, T_new
