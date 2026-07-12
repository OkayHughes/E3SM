"""Orbital mechanics (shr_orb_mod) and trace-gas profiles (trcmix).

Sources: share/util/shr_orb_mod.F90 (shr_orb_params Berger-1978 series,
shr_orb_decl, shr_orb_cosz / shr_orb_avg_cosz [Zhou et al. 2015]) and
components/eamxx/src/share/physics/eamxx_trcmix.cpp.

shr_orb_params/decl are scalar host-side numpy (called once per step);
cosz and trcmix are vectorized over columns.
"""

import numpy as np

_PI = np.pi
_PSECDEG = 1.0 / 3600.0
_DEGRAD = _PI / 180.0

# Berger 1978 series coefficients (shr_orb_params)
_OBAMP = np.array([
    -2462.2214466, -857.3232075, -629.3231835, -414.2804924, -311.7632587,
    308.9408604, -162.5533601, -116.1077911, 101.1189923, -67.6856209,
    24.9079067, 22.5811241, -21.1648355, -15.6549876, 15.3936813,
    14.6660938, -11.7273029, 10.2742696, 6.4914588, 5.8539148,
    -5.4872205, -5.4290191, 5.1609570, 5.0786314, -4.0735782, 3.7227167,
    3.3971932, -2.8347004, -2.6550721, -2.5717867, -2.4712188, 2.4625410,
    2.2464112, -2.0755511, -1.9713669, -1.8813061, -1.8468785, 1.8186742,
    1.7601888, -1.5428851, 1.4738838, -1.4593669, 1.4192259, -1.1818980,
    1.1756474, -1.1316126, 1.0896928])
_OBRATE = np.array([
    31.609974, 32.620504, 24.172203, 31.983787, 44.828336, 30.973257,
    43.668246, 32.246691, 30.599444, 42.681324, 43.836462, 47.439436,
    63.219948, 64.230478, 1.010530, 7.437771, 55.782177, 0.373813,
    13.218362, 62.583231, 63.593761, 76.438310, 45.815258, 8.448301,
    56.792707, 49.747842, 12.058272, 75.278220, 65.241008, 64.604291,
    1.647247, 7.811584, 12.207832, 63.856665, 56.155990, 77.448840,
    6.801054, 62.209418, 20.656133, 48.344406, 55.145460, 69.000539,
    11.071350, 74.291298, 11.047742, 0.636717, 12.844549])
_OBPHAS = np.array([
    251.9025, 280.8325, 128.3057, 292.7252, 15.3747, 263.7951, 308.4258,
    240.0099, 222.9725, 268.7809, 316.7998, 319.6024, 143.8050, 172.7351,
    28.9300, 123.5968, 20.2082, 40.8226, 123.4722, 155.6977, 184.6277,
    267.2772, 55.0196, 152.5268, 49.1382, 204.6609, 56.5233, 200.3284,
    201.6651, 213.5577, 17.0374, 164.4194, 94.5422, 131.9124, 61.0309,
    296.2073, 135.4894, 114.8750, 247.0691, 256.6114, 32.1008, 143.6804,
    16.8784, 160.6835, 27.5932, 348.1074, 82.6496])
_ECAMP = np.array([
    0.01860798, 0.01627522, -0.01300660, 0.00988829, -0.00336700,
    0.00333077, -0.00235400, 0.00140015, 0.00100700, 0.00085700,
    0.00064990, 0.00059900, 0.00037800, -0.00033700, 0.00027600,
    0.00018200, -0.00017400, -0.00012400, 0.00001250])
_ECRATE = np.array([
    4.2072050, 7.3460910, 17.8572630, 17.2205460, 16.8467330, 5.1990790,
    18.2310760, 26.2167580, 6.3591690, 16.2100160, 3.0651810, 16.5838290,
    18.4939800, 6.1909530, 18.8677930, 17.4255670, 6.1860010, 18.4174410,
    0.6678630])
_ECPHAS = np.array([
    28.620089, 193.788772, 308.307024, 320.199637, 279.376984, 87.195000,
    349.129677, 128.443387, 154.143880, 291.269597, 114.860583,
    332.092251, 296.414411, 145.769910, 337.237063, 152.092288,
    126.839891, 210.667199, 72.108838])
_MVAMP = np.array([
    7391.0225890, 2555.1526947, 2022.7629188, -1973.6517951, 1240.2321818,
    953.8679112, -931.7537108, 872.3795383, 606.3544732, -496.0274038,
    456.9608039, 346.9462320, -305.8412902, 249.6173246, -199.1027200,
    191.0560889, -175.2936572, 165.9068833, 161.1285917, 139.7878093,
    -133.5228399, 117.0673811, 104.6907281, 95.3227476, 86.7824524,
    86.0857729, 70.5893698, -69.9719343, -62.5817473, 61.5450059,
    -57.9364011, 57.1899832, -57.0236109, -54.2119253, 53.2834147,
    52.1223575, -49.0059908, -48.3118757, -45.4191685, -42.2357920,
    -34.7971099, 34.4623613, -33.8356643, 33.6689362, -31.2521586,
    -30.8798701, 28.4640769, -27.1960802, 27.0860736, -26.3437456,
    24.7253740, 24.6732126, 24.4272733, 24.0127327, 21.7150294,
    -21.5375347, 18.1148363, -16.9603104, -16.1765215, 15.5567653,
    15.4846529, 15.2150632, 14.5047426, -14.3873316, 13.1351419,
    12.8776311, 11.9867234, 11.9385578, 11.7030822, 11.6018181,
    -11.2617293, -10.4664199, 10.4333970, -10.2377466, 10.1934446,
    -10.1280191, 10.0289441, -10.0034259])
_MVRATE = np.array([
    31.609974, 32.620504, 24.172203, 0.636717, 31.983787, 3.138886,
    30.973257, 44.828336, 0.991874, 0.373813, 43.668246, 32.246691,
    30.599444, 2.147012, 10.511172, 42.681324, 13.650058, 0.986922,
    9.874455, 13.013341, 0.262904, 0.004952, 1.142024, 63.219948,
    0.205021, 2.151964, 64.230478, 43.836462, 47.439436, 1.384343,
    7.437771, 18.829299, 9.500642, 0.431696, 1.160090, 55.782177,
    12.639528, 1.155138, 0.168216, 1.647247, 10.884985, 5.610937,
    12.658184, 1.010530, 1.983748, 14.023871, 0.560178, 1.273434,
    12.021467, 62.583231, 63.593761, 76.438310, 4.280910, 13.218362,
    17.818769, 8.359495, 56.792707, 8.448301, 1.978796, 8.863925,
    0.186365, 8.996212, 6.771027, 45.815258, 12.002811, 75.278220,
    65.241008, 18.870667, 22.009553, 64.604291, 11.498094, 0.578834,
    9.237738, 49.747842, 2.147012, 1.196895, 2.133898, 0.173168])
_MVPHAS = np.array([
    251.9025, 280.8325, 128.3057, 348.1074, 292.7252, 165.1686, 263.7951,
    15.3747, 58.5749, 40.8226, 308.4258, 240.0099, 222.9725, 106.5937,
    114.5182, 268.7809, 279.6869, 39.6448, 126.4108, 291.5795, 307.2848,
    18.9300, 273.7596, 143.8050, 191.8927, 125.5237, 172.7351, 316.7998,
    319.6024, 69.7526, 123.5968, 217.6432, 85.5882, 156.2147, 66.9489,
    20.2082, 250.7568, 48.0188, 8.3739, 17.0374, 155.3409, 94.1709,
    221.1120, 28.9300, 117.1498, 320.5095, 262.3602, 336.2148, 233.0046,
    155.6977, 184.6277, 267.2772, 78.9281, 123.4722, 188.7132, 180.1364,
    49.1382, 152.5268, 98.2198, 97.4808, 221.5376, 168.2438, 161.1199,
    55.0196, 262.6495, 200.3284, 201.6651, 294.6547, 99.8233, 213.5577,
    154.1631, 232.7153, 138.3034, 204.6609, 106.5938, 250.4676, 332.3345,
    27.3039])

SHR_ORB_UNDEF_INT = 2000000000


def shr_orb_params(iyear_AD, eccen=None, obliq=None, mvelp=None):
    """shr_orb_params. If iyear_AD == SHR_ORB_UNDEF_INT, eccen/obliq/mvelp
    must be given (fixed-parameter mode); otherwise they are computed
    from the Berger 1978 series. Returns (eccen, obliq, mvelp, obliqr,
    lambm0, mvelpp)."""
    if iyear_AD == SHR_ORB_UNDEF_INT:
        eccen2 = eccen * eccen
        eccen3 = eccen2 * eccen
    else:
        years = -(1950.0 - float(iyear_AD))
        obsum = np.sum(_OBAMP * _PSECDEG * np.cos(
            (_OBRATE * _PSECDEG * years + _OBPHAS) * _DEGRAD))
        obliq = 23.320556 + obsum

        cossum = np.sum(_ECAMP * np.cos(
            (_ECRATE * _PSECDEG * years + _ECPHAS) * _DEGRAD))
        sinsum = np.sum(_ECAMP * np.sin(
            (_ECRATE * _PSECDEG * years + _ECPHAS) * _DEGRAD))
        eccen2 = cossum * cossum + sinsum * sinsum
        eccen = np.sqrt(eccen2)
        eccen3 = eccen2 * eccen

        if abs(cossum) <= 1.0e-8:
            if sinsum == 0.0:
                fvelp = 0.0
            elif sinsum < 0.0:
                fvelp = 1.5 * _PI
            else:
                fvelp = 0.5 * _PI
        elif cossum < 0.0:
            fvelp = np.arctan(sinsum / cossum) + _PI
        else:
            if sinsum < 0.0:
                fvelp = np.arctan(sinsum / cossum) + 2.0 * _PI
            else:
                fvelp = np.arctan(sinsum / cossum)

        mvsum = np.sum(_MVAMP * _PSECDEG * np.sin(
            (_MVRATE * _PSECDEG * years + _MVPHAS) * _DEGRAD))
        mvelp = fvelp / _DEGRAD + 50.439273 * _PSECDEG * years \
            + 3.392506 + mvsum
        while mvelp < 0.0:
            mvelp += 360.0
        while mvelp >= 360.0:
            mvelp -= 360.0

    obliqr = obliq * _DEGRAD
    mvelpp = (mvelp + 180.0) * _DEGRAD
    beta = np.sqrt(1.0 - eccen2)
    lambm0 = 2.0 * ((0.5 * eccen + 0.125 * eccen3) * (1.0 + beta)
                    * np.sin(mvelpp)
                    - 0.250 * eccen2 * (0.5 + beta) * np.sin(2.0 * mvelpp)
                    + 0.125 * eccen3 * (1.0 / 3.0 + beta)
                    * np.sin(3.0 * mvelpp))
    return eccen, obliq, mvelp, obliqr, lambm0, mvelpp


def shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr):
    """shr_orb_decl -> (delta, eccf)."""
    dayspy = 365.0
    ve = 80.5
    lambm = lambm0 + (calday - ve) * 2.0 * _PI / dayspy
    lmm = lambm - mvelpp
    sinl = np.sin(lmm)
    lamb = lambm + eccen * (2.0 * sinl + eccen * (1.25 * np.sin(2.0 * lmm)
                            + eccen * ((13.0 / 12.0) * np.sin(3.0 * lmm)
                                       - 0.25 * sinl)))
    invrho = (1.0 + eccen * np.cos(lamb - mvelpp)) / (1.0 - eccen * eccen)
    delta = np.arcsin(np.sin(obliqr) * np.sin(lamb))
    eccf = invrho * invrho
    return delta, eccf


def shr_orb_cosz(jday, lat, lon, declin, dt_avg=0.0):
    """shr_orb_cosz, vectorized over lat/lon (radians). dt_avg != 0
    selects the Zhou et al. 2015 time-averaged form."""
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    if dt_avg == 0.0:
        return (np.sin(lat) * np.sin(declin)
                - np.cos(lat) * np.cos(declin)
                * np.cos((jday - np.floor(jday)) * 2.0 * _PI + lon))
    return _shr_orb_avg_cosz(jday, lat, lon, declin, dt_avg)


def _shr_orb_avg_cosz(jday, lat, lon, declin, dt_avg):
    piover2 = _PI / 2.0
    twopi = 2.0 * _PI

    delta_lat = np.where(lat == piover2, lat - 1.0e-5,
                         np.where(lat == -piover2, lat + 1.0e-5, lat))
    phi = declin - 1.0e-5 if declin == piover2 else (
        declin + 1.0e-5 if declin == -piover2 else declin)

    cos_h = -np.tan(delta_lat) * np.tan(phi)
    h = np.where(cos_h <= -1.0, _PI,
                 np.where(cos_h >= 1.0, 0.0,
                          np.arccos(np.clip(cos_h, -1.0, 1.0))))

    t1 = (jday - int(jday)) * twopi + lon - _PI
    t1 = np.where(t1 >= _PI, t1 - twopi, np.where(t1 < -_PI, t1 + twopi, t1))
    dt = dt_avg / 86400.0 * twopi
    t2 = t1 + dt

    aa = np.sin(lat) * np.sin(declin)
    bb = np.cos(lat) * np.cos(declin)

    case1 = (t2 >= _PI) & (t1 <= _PI) & (_PI - h <= dt)
    case2 = (t2 >= -_PI) & (t1 <= -_PI) & (_PI - h <= dt)

    # default case
    tt2_d = np.where(t2 > _PI, np.clip(t2 - twopi, -h, h),
                     np.where(t2 < -_PI, np.clip(t2 + twopi, -h, h),
                              np.clip(t2, -h, h)))
    tt1_d = np.where(t1 > _PI, np.clip(t1 - twopi, -h, h),
                     np.where(t1 < -_PI, np.clip(t1 + twopi, -h, h),
                              np.clip(t1, -h, h)))
    tt3_d = np.zeros_like(h)
    tt4_d = np.zeros_like(h)

    tt2 = np.where(case1, h, np.where(case2, -twopi + h, tt2_d))
    tt1 = np.where(case1, np.clip(t1, -h, h),
                   np.where(case2, np.clip(t1, -twopi - h, -twopi + h),
                            tt1_d))
    tt4 = np.where(case1, np.clip(t2, twopi - h, twopi + h),
                   np.where(case2, np.clip(t2, -h, h), tt4_d))
    tt3 = np.where(case1, twopi - h, np.where(case2, -h, tt3_d))

    integ = (aa * (tt2 - tt1) + bb * (np.sin(tt2) - np.sin(tt1))) / dt \
        + (aa * (tt4 - tt3) + bb * (np.sin(tt4) - np.sin(tt3))) / dt
    return np.where((tt2 > tt1) | (tt4 > tt3), integ, 0.0)


# ---- trcmix ----
# gas molecular weights: Constants<T>::get_gas_mol_weight
GAS_MOL_WEIGHTS = {"h2o": 18.016, "co2": 44.0095, "o3": 47.9982,
                   "n2o": 44.0128, "co": 28.0101, "ch4": 16.04246,
                   "o2": 31.998, "n2": 28.0134, "cfc11": 136.0,
                   "cfc12": 120.0}
_MW = GAS_MOL_WEIGHTS
MWDRY = 28.966
O2MMR = 0.23143


def trcmix(name, clat_deg, pmid, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr):
    """eamxx trcmix: (ncol, nlay) mass mixing ratio for the named gas.
    clat_deg in degrees."""
    clat_r = np.asarray(clat_deg) * _PI / 180.0
    pmid = np.asarray(pmid)

    if name == "o2":
        return np.full_like(pmid, O2MMR)
    if name == "co2":
        return np.full_like(pmid, _MW["co2"] / MWDRY * co2vmr)

    params = {
        "ch4": (_MW["ch4"] / MWDRY * ch4vmr, 0.2353, 0.0, 0.2353, 0.0225489),
        "n2o": (_MW["n2o"] / MWDRY * n2ovmr, 0.3478, 0.00116, 0.4000,
                0.013333),
        "cfc11": (_MW["cfc11"] / MWDRY * f11vmr, 0.7273, 0.00606, 1.0,
                  0.013333),
        "cfc12": (_MW["cfc12"] / MWDRY * f12vmr, 0.4000, 0.00222, 0.5,
                  0.024444),
    }
    trop_mmr, s1b, s1f, s2b, s2f = params[name]
    dlat = np.abs(57.2958 * clat_r)[:, None]
    scale = np.where(dlat <= 45.0, s1b + s1f * dlat, s2b + s2f * (dlat - 45.0))
    ptrop = 250.0e2 - 150.0e2 * np.cos(clat_r)[:, None] ** 2
    return np.where(pmid >= ptrop, trop_mmr,
                    trop_mmr * (pmid / ptrop) ** scale)
