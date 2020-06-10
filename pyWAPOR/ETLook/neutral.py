from ETLook import constants as c
from numpy import log


def initial_soil_aerodynamic_resistance(u_24, z_obs=2):
    r"""
    Computes the aerodynamic resistance for soil without stability corrections
    :math:`r_{a,soil}^{0}`.

    .. math ::

        r_{a,soil}^{0}=\frac{\ln\left(\frac{z_{obs}}{z_{0,soil}}\right)
                       \ln\left(\frac{z_{obs}}{0.1z_{0,soil}}\right)}
                       {k^{2}u_{obs}}


    where the following constants are used

    * :math:`z_{0,soil}` soil roughness = 0.001 [m]
    * :math:`k` = karman constant = 0.41 [-]

    The factor 0.1 is the ratio between the surface roughness for momentum and
    heat.

    Parameters
    ----------
    u_24 : float
        daily wind speed at observation height
        :math:`u_obs`
        [m/s]

    z_obs : float
        observation height
        :math:`z_{obs}`
        [m]

    Returns
    -------
    ra_soil_init : float
        aerodynamic resistance without stability corrections
        :math:`r_{a,soil}^{0}`
        [s/m]

    """
    return (log(z_obs / c.z0_soil) * log(z_obs / (0.1 * c.z0_soil))) / (c.k**2 * u_24)


def initial_daily_evaporation(rn_24_soil, g0_24, ssvp_24, ad_24, vpd_24,
                              psy_24, r_soil, ra_soil_init):
    r"""
    Computes the soil evaporation based on the Penman Monteith equation
    adapted for soil.

    .. math ::

        E^{0}=\frac{\Delta\left(Q_{soil}^{*}-G\right)+\rho c_{p}
        \frac{\Delta_{e}}{r_{a,soil}}}{\Delta+
        \gamma\left(1+\frac{r_{soil}}{r_{a,soil}}\right)}

    where the following constants are used

    * :math:`c_{p}` specific heat for dry air = 1004 [J kg-1 K-1]
    * :math:`k` = karman constant = 0.41 [-]

    Parameters
    ----------
    rn_24_soil : float
        daily net radiation for soil
        :math:`Q_{soil}^{*}`
        [W m-2]
    g0_24 : float
        daily soil heat flux
        :math:`G`
        [W m-2]
    ssvp_24 : float
       daily slope of saturated vapour pressure curve
       :math:`\Delta`
       [mbar K-1]
    ad_24 : float
        daily air density
        :math:`\rho`
        [kg m-3]
    vpd_24 : float
        daily vapour pressure deficit
        :math:`\Delta_{e}`
        [mbar]
    psy_24 : float
        daily psychrometric constant
        :math:`\gamma`
        [mbar K-1]
    r_soil : float
        soil resistance
        :math:`r_{soil}`
        [m s-1]
    ra_soil_init : float
        initial soil aerodynamic resistance
        :math:`r_{a,soil}`
        [m s-1]

    Returns
    -------
    e_24_init : float
        initial estimate radiation equivalent daily evaporation
        :math:`E^{0}`
        [W m-2]

    """
    numerator = (ssvp_24 * (rn_24_soil - g0_24) +
                 ad_24 * c.sh * (vpd_24 / ra_soil_init))
    denominator = (ssvp_24 + psy_24 * (1 + r_soil / ra_soil_init))
    return numerator / denominator


def initial_daily_evaporation_mm(e_24_init, lh_24):
    r"""
    Computes the soil evaporation based on the Penman Monteith equation
    adapted for soil.

    .. math ::

        E^{0}=E^{0}d_{sec}\lambda_{24}

    where the following constants are used

    * :math:`d_{sec}` seconds in the day = 86400 [s]

    Parameters
    ----------
    e_24_init : float
        initial estimate daily evaporation
        :math:`E^{0}`
        [W m-2]
    lh_24 : float
        daily latent heat of evaporation
        :math:`\lambda_{24}`
        [J/kg]

    Returns
    -------
    e_24_init_mm : float
        initial estimate daily evaporation in mm
        :math:`E^{0}`
        [mm d-1]
    """
    return e_24_init * c.day_sec / lh_24


def initial_canopy_aerodynamic_resistance(u_24, z0m, z_obs=2):
    r"""
    Computes the aerodynamic resistance for a canopy soil without stability
    corrections :math:`r_{a,}^{0}`.

    .. math ::

        r_{a,canopy}^{0}=\frac{\ln\left(\frac{z_{obs}}{z_{0,m}}\right)\ln
        \left(\frac{z_{obs}}{0.1z_{0,m}}\right)}{k^{2}u_{obs}}

    where the following constants are used

    * :math:`k` = karman constant = 0.41 [-]

    The factor 0.1 is the ratio between the surface roughness for momentum and
    heat.

    Parameters
    ----------
    u_24 : float
        daily wind speed at observation height
        :math:`u_obs`
        [m/s]
    z0m : float
        roughness length
        :math:`z_{0,m}`
        [m]
    z_obs : float
        observation height
        :math:`z_{obs}`
        [m]

    Returns
    -------
    ra_canopy_init : float
        canopy resistance without stability corrections
        :math:`r_{a,canopy}^{0}`
        [s/m]
    """
    return (log(z_obs / z0m) * log(z_obs / (0.1 * z0m))) / (c.k**2 * u_24)


def initial_daily_transpiration(rn_24_canopy, ssvp_24, ad_24, vpd_24,
                                psy_24, r_canopy, ra_canopy_init):
    r"""
    Computes the soil evaporation based on the Penman Monteith equation adapted
    for soil.

    .. math ::

        T_{0}=\frac{\Delta\left(Q_{canopy}^{*}\right)
        +\rho c_{p}\frac{\Delta_{e}}
        {r_{a,canopy}}}{\Delta+
        \gamma\left(1+\frac{r_{canopy}}{r_{a,canopy}}\right)}

    where the following constants are used

    * :math:`c_{p}` specific heat for dry air = 1004 [J kg-1 K-1]
    * :math:`k` = karman constant = 0.41 [-]

    Parameters
    ----------
    rn_24_canopy : float
        daily net radiation for the canopy
        :math:`Q_{soil}^{*}`
        [W m-2]
    ssvp_24 : float
       daily slope of saturated vapour pressure curve
       :math:`\Delta`
       [mbar K-1]
    ad_24 : float
        daily air density
        :math:`\rho`
        [kg m-3]
    vpd_24 : float
        daily vapour pressure deficit
        :math:`\Delta_{e}`
        [mbar]
    psy_24 : float
        daily psychrometric constant
        :math:`\gamma`
        [mbar K-1]
    r_canopy : float
        canopy resistance
        :math:`r_{canopy}`
        [m s-1]
    ra_canopy_init : float
        initial canopy aerodynamic resistance
        :math:`r_{a,canopy}`
        [m s-1]

    Returns
    -------
    t_24_init : float
        initial estimate radiation equivalent daily transpiration
        :math:`T^{0}`
        [W m-2]
    """
    numerator = (ssvp_24 * rn_24_canopy + ad_24 *
                 c.sh * (vpd_24 / ra_canopy_init))
    denominator = (ssvp_24 + psy_24 * (1 + r_canopy / ra_canopy_init))
    return numerator / denominator


def initial_daily_transpiration_mm(t_24_init, lh_24):
    r"""
    Computes the canopy transpiration based on the Penman Monteith equation
    adapted for canopy.

    .. math ::

        T^{0}=T^{0}d_{sec}\lambda_{24}

    where the following constants are used

    * :math:`d_{sec}` seconds in the day = 86400 [s]

    Parameters
    ----------
    t_24_init : float
        initial estimate daily transpiration
        :math:`E^{0}`
        [W m-2]
    lh_24 : float
        daily latent heat of evaporation
        :math:`\lambda_{24}`
        [J/kg]

    Returns
    -------
    t_24_init_mm : float
        initial estimate daily transpiration in mm
        :math:`T^{0}`
        [mm d-1]
    """
    return t_24_init * c.day_sec / lh_24
