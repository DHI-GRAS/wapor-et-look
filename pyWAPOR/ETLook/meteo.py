"""
    The meteo module contains all functions related to meteorological
    variables. All meteorological functions can be calculated on a daily
    or instantaneous basis. Base functions are available also. The daily
    functions have a 'daily' extension, instantaneous functions have a 'inst'
    extension

"""
from pyWAPOR.ETLook import constants as c
import numpy as np

def air_temperature_celcius(t_air_k):
    r"""
    Converts air temperature from Kelvin to Celcius, where 0 degrees Celcius
    is equal to 273.15 degrees Kelvin

    Parameters
    ----------
    t_air_k : float
        air temperature
        :math:`T_a`
        [K]

    Returns
    -------
    t_air_c : float
        air temperature
        :math:`T_a`
        [C]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.air_temperature_celcius(12.5)
    285.65
    """
    return t_air_k - c.zero_celcius


def air_temperature_kelvin(t_air):
    r"""
    Converts air temperature from Celcius to Kelvin, where 0 degrees Celcius
    is equal to 273.15 degrees Kelvin

    Parameters
    ----------
    t_air : float
        air temperature
        :math:`T_a`
        [C]

    Returns
    -------
    t_air_k : float
        air temperature
        :math:`T_a`
        [K]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.air_temperature_kelvin(12.5)
    285.65
    """
    return t_air + c.zero_celcius


def air_temperature_kelvin_daily(t_air_24):
    """Like :func:`air_temperature_kelvin` but as a daily average

    Parameters
    ----------
    t_air_24 : float
        daily air temperature
        :math:`T_{a,24}`
        [C]

    Returns
    -------
    t_air_k_24 : float
        daily air temperature
        :math:`T_{a,24}`
        [K]
    """
    return air_temperature_kelvin(t_air_24)


def air_temperature_kelvin_inst(t_air_i):
    """Like :func:`air_temperature_kelvin` but as an instantaneous value

    Parameters
    ----------
    t_air_i : float
        instantaneous air temperature
        :math:`T_{a,i}`
        [C]

    Returns
    -------
    t_air_k_i : float
        instantaneous air temperature
        :math:`T_{a,i}`
        [K]
    """
    return air_temperature_kelvin(t_air_i)


def wet_bulb_temperature_kelvin_inst(t_wet_i):
    """
    Converts wet bulb temperature from Celcius to Kelvin, where 0
    degrees Celcius is equal to 273.15 degrees Kelvin

    Parameters
    ----------
    t_wet_i : float
        instantaneous wet bulb temperature
        :math:`T_{w,i}`
        [C]

    Returns
    -------
    t_wet_k_i : float
        instantaneous wet bulb temperature
        :math:`T_{w,i}`
        [K]
    """
    return air_temperature_kelvin(t_wet_i)


def disaggregate_air_temperature(t_air_coarse, z, z_coarse, lapse=-0.006):
    r"""
    Disaggregates GEOS or MERRA or another coarse scale air temperature using
    two digital elevation models. One DEM for the target resolution, another
    DEM smoothed from the original air temperature resolution to the target
    resolution.

    .. math ::
        T_{a}=T_{a,c}+(z-z_{c})L_{T}-T_{K,0}

    where the following constant is used

    * :math:`T_{K,0}` = 273.15 K is equal to 0 degrees Celsius

    Parameters
    ----------
    t_air_coarse : float
        air temperature at coarse resolution
        :math:`T_{a,c}`
        [K]
    z : float
        elevation
        :math:`z`
        [m]
    z_coarse : float
        elevation at coarse resolution
        :math:`z_{c}`
        [m]
    lapse : float
        lapse rate
        :math:`L_{T}`
        [K m-1]

    Returns
    -------
    t_air : float
        air temperature
        :math:`T_{a}`
        [C]

    Notes
    -----
    The input air temperature is specified in Kelvin. The output air
    temperature is specified in C.

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.disaggregate_air_temperature(24.5+273.15, 10, 5)
    24.47
    """
    return t_air_coarse + ((z - z_coarse) * lapse) - c.zero_celcius


def disaggregate_air_temperature_daily(t_air_24_coarse, z, z_coarse, lapse=-0.006):
    r"""Like :func:`disaggregate_air_temperature` but as a daily average

    Parameters
    ----------
    t_air_24_coarse : float
        daily air temperature at coarse resolution
        :math:`T_{a,24,c}`
        [K]
    z : float
        elevation
        :math:`z`
        [m]
    z_coarse : float
        elevation at coarse resolution
        :math:`z_{c}`
        [m]
    lapse : float
        lapse rate
        :math:`L`
        [K m-1]

    Returns
    -------
    t_air_24 : float
        daily air temperature
        :math:`T_{a,24}`
        [C]

    Notes
    -----
    The input air temperature is specified in Kelvin. The output air
    temperature is specified in C.
    """
    return disaggregate_air_temperature(t_air_24_coarse, z, z_coarse, lapse)


def disaggregate_air_temperature_inst(t_air_i_coarse, z, z_coarse, lapse=-0.006):
    r"""Like :func:`disaggregate_air_temperature` but as a instantaneous value

    Parameters
    ----------
    t_air_i_coarse : float
        instantaneous air temperature at coarse resolution
        :math:`T_{a,i,c}`
        [K]
    z : float
        elevation
        :math:`z`
        [m]
    z_coarse : float
        elevation at coarse resolution
        :math:`z_{c}`
        [m]
    lapse : float
        lapse rate
        :math:`L`
        [K m-1]

    Returns
    -------
    t_air_i : float
        instantaneous air temperature
        :math:`T_{a,i}`
        [C]

    Notes
    -----
    The input air temperature is specified in Kelvin. The output air
    temperature is specified in C.
    """
    return disaggregate_air_temperature(t_air_i_coarse, z, z_coarse, lapse)


def disaggregate_dew_point_temperature_inst(t_dew_coarse_i, z, z_coarse, lapse_dew=-0.002):
    r"""
    Disaggregates geos dew point temperature using lapse rate and difference between
    smoothed coarse scale DEM and fine scale DEM

    Parameters
    ----------
    t_dew_coarse_i : float
        coarse instantaneous dew point temperature
        :math:`T_{dew,coarse}`
        [C]
    z : float
        elevation
        :math:`z`
        [m]
    z_coarse : float
        smoothed elevation at coarse resolution
        :math:`z`
        [m]
    lapse_dew : float
        lapse rate
        :math:`L`
        [K m-1]

    Returns
    -------
    t_dew_i : float
        instantaneous dew point temperature
        :math:`T_{dew,i}`
        [C]
    """
    return t_dew_coarse_i + ((z - z_coarse) * lapse_dew)


def vapour_pressure_from_specific_humidity(qv, p_air):
    r"""
    Computes the vapour pressure :math:`e_a` in [mbar] using specific humidity
    and surface pressure

    .. math ::
        e_{a}=\frac{q_{v}P}{\varepsilon}

    where the following constant is used

    * :math:`\varepsilon` = ratio of molecular weight of water to
      dry air = 0.622 [-]


    Parameters
    ----------
    qv : float
        specific humidity
        :math:`q_{v}`
        [kg/kg]
    p_air : float
        air pressure
        :math:`P`
        [mbar]

    Returns
    -------
    vp : float
        vapour pressure
        :math:`e_{a}`
        [mbar]
    """
    return (qv * p_air) / c.r_mw


def vapour_pressure_from_specific_humidity_daily(qv_24, p_air_24):
    """Like :func:`vapour_pressure_from_specific_humidity` but as a daily average

    Parameters
    ----------
    qv_24 : float
        daily specific humidity
        :math:`q_{v,24}`
        [kg/kg]
    p_air_24 : float
        daily air pressure
        :math:`P_{24}`
        [mbar]

    Returns
    -------
    vp_24 : float
        daily vapour pressure
        :math:`e_{a,24}`
        [mbar]
    """
    return vapour_pressure_from_specific_humidity(qv_24, p_air_24)


def saturated_vapour_pressure_minimum(t_air_min_coarse):
    """Like :func:`saturated_vapour_pressure` but based on daily minimum air temperature. This
    is only relevant for reference ET calculations

    Parameters
    ----------
    t_air_min_coarse : float
        daily minimum air temperature
        :math:`T_{a,min}`
        [C]

    Returns
    -------
    svp_24_min : float
        daily saturated vapour pressure based on minimum air temperature
        :math:`e_{s,min}`
        [mbar]

    """
    return saturated_vapour_pressure(t_air_min_coarse)


def saturated_vapour_pressure_maximum(t_air_max_coarse):
    """Like :func:`saturated_vapour_pressure` but based on daily maximum air temperature. This
    is only relevant for reference ET calculations

    Parameters
    ----------
    t_air_max_coarse : float
        daily maximum air temperature
        :math:`T_{a,max}`
        [C]

    Returns
    -------
    svp_24_max : float
        daily saturated vapour pressure based on maximum air temperature
        :math:`e_{s,max}`
        [mbar]

    """
    return saturated_vapour_pressure(t_air_max_coarse)


def saturated_vapour_pressure_average(svp_24_max, svp_24_min):
    """
    Average saturated vapour pressure based on two saturated vapour pressure values
    calculated using minimum and maximum air temperature respectively. This is preferable
    to calculating saturated vapour pressure using the average air temperature, because
    of the strong non-linear relationship between saturated vapour pressure and air
    temperature

    .. math ::
        e_{s}=\frac{e^{0}\left(T_{max}\right)+e^{0}\left(T_{mon}\right)}{2}

    Parameters
    ----------
    svp_24_max : float
        daily saturated vapour pressure based on maximum air temperature
        :math:`e_{s,max}`
        [mbar]
    svp_24_min : float
        daily saturated vapour pressure based on minimum air temperature
        :math:`e_{s,min}`
        [mbar]

    Returns
    -------
    svp_24 : float
        daily saturated vapour pressure
        :math:`e_{s,24}`
        [mbar]

    """
    return (svp_24_max + svp_24_min)/2


def vapour_pressure_from_specific_humidity_inst(qv_i, p_air_i):
    """Like :func:`vapour_pressure_from_specific_humidity` but as an instantaneous value

    Parameters
    ----------
    qv_i : float
        instantaneous specific humidity
        :math:`q_{v,i}`
        [kg/kg]
    p_air_i : float
        instantaneous air pressure
        :math:`P_{i}`
        [mbar]

    Returns
    -------
    vp_i : float
        instantaneous vapour pressure
        :math:`e_{a,i}`
        [mbar]
    """
    return vapour_pressure_from_specific_humidity(qv_i, p_air_i)


def saturated_vapour_pressure(t_air):
    r"""
    Computes saturated vapour pressure :math:`e_s` [mbar], it provides
    the vapour pressure when the air is fully saturated with water. It is
    related to air temperature :math:`T_a` [C]:

    .. math ::
        e_{s}=6.108\exp\left[\frac{17.27T_{a}}{T_{a}+237.3}\right]

    Parameters
    ----------
    t_air : float
        air temperature
        :math:`T_a`
        [C]

    Returns
    -------
    svp : float
        saturated vapour pressure
        :math:`e_s`
        [mbar]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.saturated_vapour_pressure(20)
    23.382812709274457

    .. plot:: pyplots/meteo/plot_saturated_vapour_pressure.py
    """
    return 6.108 * np.exp(((17.27 * t_air) / (237.3 + t_air)))


def saturated_vapour_pressure_daily(t_air_24):
    """Like :func:`saturated_vapour_pressure` but as a daily average

    Parameters
    ----------
    t_air_24 : float
        daily air temperature
        :math:`T_{a,24}`
        [C]

    Returns
    -------
    svp_24 : float
        daily saturated vapour pressure
        :math:`e_{s,24}`
        [mbar]

    """
    return saturated_vapour_pressure(t_air_24)


def saturated_vapour_pressure_inst(t_air_i):
    """Like :func:`saturated_vapour_pressure` but as an instantaneous value

    Parameters
    ----------
    t_air_i : float
        instantaneous air temperature
        :math:`T_{a,i}`
        [C]

    Returns
    -------
    svp_i : float
        instantaneous saturated vapour pressure
        :math:`e_{s,i}`
        [mbar]

    """
    return saturated_vapour_pressure(t_air_i)


def vapour_pressure(svp, rh):
    r"""
    Computes the vapour pressure :math:`e_a` in [mbar]

    .. math ::
        e_{a}=\frac{\phi}{100}e_{s}

    Parameters
    ----------
    svp : float
        saturated vapour pressure
        :math:`e_s`
        [mbar]
    rh : float
        relative humidity
        :math:`\phi`
        [%]

    Returns
    -------
    vp : float
        vapour pressure
        :math:`e_{a}`
        [mbar]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.vapour_pressure(rh=75, svp=meteo.saturated_vapour_pressure(20))
    17.537109531955842

    """
    return 0.01 * rh * svp


def slope_saturated_vapour_pressure(t_air):
    r"""
    Computes the rate of change of vapour pressure :math:`\Delta` in [mbar K-1]
    for a given air temperature :math:`T_a`. It is a function of the air
    temperature :math:`T_a` and the saturated vapour pressure :math:`e_s`
    [mbar] which in itself is a function of :math:`T_a`.

    .. math ::
        \Delta=\frac{4098e_{s}}{\left(237.3+T_{a}\right)^{2}}

    for :math:`e_s` see :func:`saturated_vapour_pressure`

    Parameters
    ----------
    t_air : float
       air temperature
       :math:`T_a`
       [C]

    Returns
    -------
    ssvp : float
       slope of saturated vapour pressure curve
       :math:`\Delta`
       [mbar K-1]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.slope_saturated_vapour_pressure(20)
    1.447401881124136

    .. plot:: pyplots/meteo/plot_slope_saturated_vapour_pressure.py
    """
    svp = saturated_vapour_pressure(t_air)
    return (4098 * svp) / (237.3 + t_air) ** 2


def slope_saturated_vapour_pressure_daily(t_air_24):
    """Like :func:`slope_saturated_vapour_pressure` but as a daily average

    Parameters
    ----------
    t_air_24 : float
       daily air temperature
       :math:`T_{a,24}`
       [C]

    Returns
    -------
    ssvp_24 : float
       daily slope of saturated vapour pressure curve
       :math:`\Delta_{24}`
       [mbar K-1]


    """
    return slope_saturated_vapour_pressure(t_air_24)


def slope_saturated_vapour_pressure_inst(t_air_i):
    """Like :func:`slope_saturated_vapour_pressure` but as an instantaneous
    value

    Parameters
    ----------
    t_air_i : float
        instantaneous air temperature
        :math:`T_{a,i}`
        [C]

    Returns
    -------
    ssvp_i : float
        instantaneous slope of saturated vapour pressure curve
        :math:`e_{s,i}`
        [mbar]
    """
    return slope_saturated_vapour_pressure(t_air_i)


def vapour_pressure_deficit(svp, vp):
    r"""
    Computes the vapour pressure deficit :math:`\Delta_e` in [mbar]

    .. math ::
        \Delta_e=e_s-e_a

    Parameters
    ----------
    svp : float
        saturated vapour pressure
        :math:`e_s`
        [mbar]
    vp : float
        actual vapour pressure
        :math:`e_a`
        [mbar]

    Returns
    -------
    vpd : float
       vapour pressure deficit
       :math:`\Delta_e`
       [mbar]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.vapour_pressure_deficit(12.5, 5.4)
    7.1
    >>> meteo.vapour_pressure_deficit(vp=5.4, svp=12.3)
    6.9

    """
    vpd = svp - vp
    vpd[vpd < 0] = 0

    return vpd


def vapour_pressure_deficit_daily(svp_24, vp_24):
    """Like :func:`vapour_pressure_deficit` but as a daily average

    Parameters
    ----------
    svp_24 : float
        daily saturated vapour pressure
        :math:`e_{s,24}`
        [mbar]
    vp_24 : float
        daily actual vapour pressure
        :math:`e_{a,24}`
        [mbar]

    Returns
    -------
    vpd_24 : float
       daily vapour pressure deficit
       :math:`\Delta_{e,24}`
       [mbar]
    """
    return vapour_pressure_deficit(svp_24, vp_24)


def air_pressure(z, p_air_0=1013.25):
    r"""
    Computes air pressure :math:`P` at a certain elevation derived from the
    air pressure at sea level :math:`P_0`. Air pressure decreases with
    increasing elevation.

    .. math ::
        P=P_{0}\left(\frac{T_{ref,0,K}-\alpha_{1}\left(z-z_{0}\right)}
        {T_{ref,0,K}}\right)^{\frac{g}{\alpha{}_{1^{R}}}}

    where the following constants are used

    * :math:`P_0` = air pressure [mbar] at sea level :math:`z_0` = 1013.25 mbar
    * :math:`T_{ref,0,K}` = reference temperature [K] at sea level
      :math:`z_0` = 293.15 K
    * :math:`g` = gravitational acceleration = 9.807 [m/s2]
    * :math:`R` = specific gas constant = 287.0 [J kg-1 K-1]
    * :math:`\alpha_1` = constant lapse rate for moist air = 0.0065 [K m-1]

    Parameters
    ----------
    z : float
        elevation
        :math:`z`
        [m]
    p_air_0 : float
        air pressure at sea level
        :math:`P_0`
        [mbar]

    Returns
    -------
    p_air : float
        air pressure
        :math:`P`
        [mbar]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.air_pressure(z=1000)
    900.5832172948869

    .. plot:: pyplots/meteo/plot_air_pressure.py
    """
    return p_air_0 * ((c.t_ref + c.lapse * (z - c.z_ref)) / c.t_ref) ** c.power


def air_pressure_daily(z, p_air_0_24=1013.25):
    r"""Like :func:`air_pressure` but as a daily average

    Parameters
    ----------
    z : float
        elevation
        :math:`z`
        [m]
    p_air_0_24 : float
        daily air pressure at sea level
        :math:`P_{0,24}`
        [mbar]

    Returns
    -------
    p_air_24 : float
        daily air pressure
        :math:`P_{24}`
        [mbar]

    """
    return air_pressure(z, p_air_0_24)


def dry_air_density(p_air, vp, t_air_k):
    r"""
    Computes dry air density :math:`\rho_{d}` in [kg m-3]

    .. math ::
        \rho_{d}=\frac{P-e_{a}}{\Re T_{a,K}}

    where the following constants are used

    * :math:`\Re` = gas constant for dry air = 2.87 mbar K-1 m3 kg-1

    Parameters
    ----------
    p_air : float
        air pressure
        :math:`P`
        [mbar]
    vp : float
        vapour pressure
        :math:`e_{a}`
        [mbar]
    t_air_k : float
        daily air temperature
        :math:`T_{a}`
        [K]

    Returns
    -------
    ad_dry : float
        dry air density
        :math:`\rho_{d}`
        [kg m-3]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.dry_air_density(p_air=900, vp=17.5, t_air_k=293.15)
    1.0489213344656534

    .. plot:: pyplots/meteo/plot_dry_air_density.py

    """
    return (p_air - vp) / (t_air_k * c.gc_dry)


def dry_air_density_daily(p_air_24, vp_24, t_air_k_24):
    r"""
    Like :func:`dry_air_density` but as a daily average

    Parameters
    ----------
    p_air_24 : float
        daily air pressure
        :math:`P_{24}`
        [mbar]
    vp_24 : float
        daily vapour pressure
        :math:`e_{a,24}`
        [mbar]
    t_air_k_24 : float
        daily air temperature
        :math:`T_{a,24}`
        [K]

    Returns
    -------
    ad_dry_24 : float
        daily dry air density
        :math:`\rho_{d,24}`
        [kg m-3]
    """
    return dry_air_density(p_air_24, vp_24, t_air_k_24)


def dry_air_density_inst(p_air_i, vp_i, t_air_k_i):
    r"""
    Like :func:`dry_air_density` but as an instantaneous value

    Parameters
    ----------
    p_air_i : float
        instantaneous air pressure
        :math:`P_{i}`
        [mbar]
    vp_i : float
        instantaneous vapour pressure
        :math:`e_{a,i}`
        [mbar]
    t_air_k_i : float
        instantaneous air temperature
        :math:`T_{a,i}`
        [K]

    Returns
    -------
    ad_dry_i : float
        instantaneous dry air density
        :math:`\rho_{d,i}`
        [kg m-3]
    """
    return dry_air_density(p_air_i, vp_i, t_air_k_i)


def moist_air_density(vp, t_air_k):
    r"""
    Computes moist air density :math:`\rho_{s}` in [kg m-3]

    .. math ::
        \rho_{s}=\frac{e_{a}}{R_{v}T_{a,K}}

    where the following constants are used

    * :math:`R_v` = gas constant for moist air = 4.61 mbar K-1 m3 kg-1

    Parameters
    ----------
    vp : float
        vapour pressure
        :math:`e_{a}`
        [mbar]
    t_air_k : float
        air temperature
        :math:`T_{a,K}`
        [K]

    Returns
    -------
    ad_moist : float
        moist air density
        :math:`\rho_{s}`
        [kg m-3]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.moist_air_density(vp=17.5, t_air_k = 293.15)
    0.012949327800393881

    .. plot:: pyplots/meteo/plot_moist_air_density.py

    """
    return vp / (t_air_k * c.gc_moist)


def moist_air_density_daily(vp_24, t_air_k_24):
    r"""
    Like :func:`moist_air_density` but as a daily average

    Parameters
    ----------
    vp_24 : float
        daily vapour pressure
        :math:`e_{a,24}`
        [mbar]
    t_air_k_24 : float
        daily air temperature
        :math:`T_{a,K,24}`
        [K]

    Returns
    -------
    ad_moist_24 : float
        daily moist air density
        :math:`\rho_{s,24}`
        [kg m-3]

    """
    return moist_air_density(vp_24, t_air_k_24)


def moist_air_density_inst(vp_i, t_air_k_i):
    r"""
    Like :func:`moist_air_density` but as an instantaneous value

    Parameters
    ----------
    vp_i : float
        instantaneous vapour pressure
        :math:`e_{a,i}`
        [mbar]
    t_air_k_i : float
        instantaneous air temperature
        :math:`T_{a,K,i}`
        [K]

    Returns
    -------
    ad_moist_i : float
        instantaneous moist air density
        :math:`\rho_{s,i}`
        [kg m-3]

    """
    return moist_air_density(vp_i, t_air_k_i)


def air_density(ad_dry, ad_moist):
    r"""
    Computes air density :math:`\rho` in [kg m-3]

    .. math ::
        \rho=\rho_{s}+\rho_{d}

    Parameters
    ----------
    ad_dry : float
        dry air density
        :math:`\rho_{d}`
        [kg m-3]
    ad_moist : float
        moist air density
        :math:`\rho_{s}`
        [kg m-3]

    Returns
    -------
    ad : float
        air density
        :math:`\rho`
        [kg m-3]

    Examples
    --------
    >>> from ETLook import meteo
    >>> ad_moist = meteo.moist_air_density(vp=17.5, t_air_k = 293.15)
    >>> ad_dry = meteo.dry_air_density(p_air=900, vp=17.5, t_air_k=293.15)
    >>> meteo.air_density(ad_dry=ad_dry, ad_moist=ad_moist)
    1.0618706622660472

    .. plot:: pyplots/meteo/plot_air_density.py

    """
    return ad_dry + ad_moist


def air_density_daily(ad_dry_24, ad_moist_24):
    r"""

    Like :func:`air_density` but as a daily average

    Parameters
    ----------
    ad_dry_24 : float
        daily dry air density
        :math:`\rho_{d,24}`
        [kg m-3]
    ad_moist_24 : float
        daily moist air density
        :math:`\rho_{s,24}`
        [kg m-3]

    Returns
    -------
    ad_24 : float
        daily air density
        :math:`\rho_{24}`
        [kg m-3]

    """
    return air_density(ad_dry_24, ad_moist_24)


def air_density_inst(ad_dry_i, ad_moist_i):
    r"""
    Like :func:`air_density` but as a instantaneous value

    Parameters
    ----------
    ad_dry_i : float
        instantaneous dry air density
        :math:`\rho_{d,i}`
        [kg m-3]
    ad_moist_i : float
        instantaneous moist air density
        :math:`\rho_{s,i}`
        [kg m-3]

    Returns
    -------
    ad_i : float
        instantaneous air density
        :math:`\rho_{i}`
        [kg m-3]

    """
    return air_density(ad_dry_i, ad_moist_i)


def latent_heat(t_air):
    r"""
    Computes latent heat of evaporation :math:`\lambda` [J kg-1], describing
    the amount of energy needed to evaporate one kg of water at constant
    pressure and temperature. At higher temperatures less energy will be
    required than at lower temperatures.

    .. math ::

        \lambda=(\lambda_0 + \Delta_\lambda T_{a})

    where the following constants are used

    * :math:`\lambda_0` = latent heat of evaporation at 0 C = 2501000 [J kg-1]
    * :math:`\Delta_\lambda` = rate of change of latent heat with respect to temperature = -2361 [J Kg-1 C-1]

    Parameters
    ----------
    t_air : float
        air temperature
        :math:`T_a`
        [C]

    Returns
    -------
    lh : float
        latent heat of evaporation
        :math:`\lambda`
        [J/kg]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.latent_heat(20)
    2453780.0

    .. plot:: pyplots/meteo/plot_latent_heat.py
    """
    return c.lh_0 + c.lh_rate * t_air


def latent_heat_daily(t_air_24):
    """Like :func:`latent_heat` but as a daily average

    Parameters
    ----------
    t_air_24 : float
        daily air temperature
        :math:`T_{a,24}`
        [C]

    Returns
    -------
    lh_24 : float
        daily latent heat of evaporation
        :math:`\lambda_{24}`
        [J/kg]

    """
    return latent_heat(t_air_24)


def psychrometric_constant(p_air, lh):
    r"""
    Computes the psychrometric constant :math:`\gamma` [mbar K-1] which
    relates the partial pressure of water in air to the air temperature

    .. math ::

        \gamma=\frac{Pc_{p}}{\varepsilon\lambda}

    where the following constants are used

    * :math:`c_{p}` = specific heat for dry air = 1004 [J Kg-1 K-1]
    * :math:`\varepsilon` = ratio of molecular weight of water to
      dry air = 0.622 [-]

    Parameters
    ----------
    p_air : float
        air pressure
        :math:`P`
        [mbar]
    lh : float
        latent heat of evaporation
        :math:`\lambda`
        [J/kg]

    Returns
    -------
    psy : float
        psychrometric constant
        :math:`\gamma`
        [mbar K-1]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.psychrometric_constant(p_air = 1003.0, lh = 2500000.0)
    0.6475961414790997
    >>> meteo.psychrometric_constant(1003.0, 2500000.0)
    0.6475961414790997

    .. plot:: pyplots/meteo/plot_psychrometric_constant.py

    """
    return (c.sh * p_air) / (c.r_mw * lh)


def psychrometric_constant_daily(p_air_24, lh_24):
    """Like :func:`psychrometric_constant` but as a daily average

    Parameters
    ----------
    p_air_24 : float
        daily air pressure
        :math:`P_{24}`
        [mbar]
    lh_24 : float
        daily latent heat of evaporation
        :math:`\lambda_{24}`
        [J/kg]

    Returns
    -------
    psy_24 : float
        daily psychrometric constant
        :math:`\gamma_{24}`
        [mbar K-1]

    """
    return psychrometric_constant(p_air_24, lh_24)


def wind_speed_blending_height(u, z_obs=2, z_b=100):
    r"""
    Computes the wind speed at blending height :math:`u_{b}` [m/s] using the
    logarithmic wind profile

    .. math ::
        u_{b}=\frac{u_{obs}\ln\left(\frac{z_{b}}{z_{0,m}}\right)}
        {\ln\left(\frac{z_{obs}}{z_{0,m}}\right)}

    Parameters
    ----------
    u : float
        wind speed at observation height
        :math:`u_{obs}`
        [m/s]
    z_obs : float
        observation height of wind speed
        :math:`z_{obs}`
        [m]
    z_b : float
        blending height
        :math:`z_{b}`
        [m]

    Returns
    -------
    u_b : float
        wind speed at blending height
        :math:`u_{b}`
        [m/s]

    Examples
    --------
    >>> from ETLook import meteo
    >>> meteo.wind_speed_blending_height(u=3.0, z_obs=2, z_b=100)
    5.4646162953650572

    .. plot:: pyplots/meteo/plot_wind_speed.py

    """

    z0m = 0.0171

    ws = (c.k * u) / np.log(z_obs / z0m) * np.log(z_b / z0m) / c.k
    ws = np.clip(ws, 1, 150)

    return ws


def wind_speed_blending_height_daily(u_24, z_obs=2, z_b=100):
    """Like :func:`wind_speed_blending_height` but as a daily average

    Parameters
    ----------
    u_24 : float
        daily wind speed at observation height
        :math:`u_{obs,24}`
        [m/s]
    z_obs : float
        observation height of wind speed
        :math:`z_{obs}`
        [m]
    z_b : float
        blending height
        :math:`z_{b}`
        [m]

    Returns
    -------
    u_b_24 : float
        daily wind speed at blending height
        :math:`u_{b, 24}`
        [m/s]

    """
    return wind_speed_blending_height(u_24, z_obs, z_b)


def air_pressure_kpa2mbar(p_air_kpa):
    """Like :func:`p_air`

    Parameters
    ----------
    p_air_kpa : float
        air pressure
        :math:`Pair_{a}`
        [kpa]

    Returns
    -------
    p_air_mbar : float
        air pressure
        :math:`Pair_{a}`
        [mbar]
    """
    return p_air_kpa * 10