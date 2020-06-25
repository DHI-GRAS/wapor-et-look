import numpy as np
from pyWAPOR.ETLook import constants as c

def interception_wm2(int_mm, lh_24):
    r"""
    Computes the energy equivalent for the interception in Wm-2 if it
    is provide in mm/day

    .. math ::
        I = \frac{\lambda I^*}{86400}

    Parameters
    ----------
    int_mm : float
        interception
        :math:`I^*`
        [mm day-1]

    lh_24 : float
        daily latent heat for evaporation
        :math:`\lambda`
        [J kg-1]

    Returns
    -------
    int_wm2 : float
        interception
        :math:`I`
        [W m-2]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> import ETLook.meteo as meteo
    >>> lh = meteo.latent_heat_daily(20.0)
    >>> rad.interception_wm2(1.0, lh)
    28.40023148148148

    """
    return int_mm * (lh_24/c.day_sec)


def soil_fraction(lai):
    r"""
    Computes the effect of the vegetation has in separating the net radiation
    into a soil and canopy component. If the canopy has a full cover almost
    no radiation reaches the soil.

    .. math ::
        s_f = \exp\left(-0.6*I_{lai}\right)

    Parameters
    ----------
    lai : float
        leaf area index
        :math:`I_{lai}`
        [-]

    Returns
    -------
    sf_soil : float
        soil fraction
        :math:`s_f`
        [-]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.soil_fraction(3.0)
    0.16529888822158656
    """
    return np.exp(-0.6*lai)


def longwave_radiation_fao_etref(t_air_k_24, vp_24, trans_24):
    r"""
    Computes the net longwave radiation according to the FAO 56 manual. For the
    reference ET calculation the values for vp_slope, vp_offset, lw_slope and
    lw_offset are being provided as defaults

    .. math ::
        L^{*}=\sigma\left(T_{a,K}\right)^{4}
        \left(vp_off-vp_slp\sqrt{0.1e_{a}}
        \right)\left(lw_slp\frac{\tau}{0.75}+lw_off\right)

    where the following constant is used

    * :math:`\sigma` = Stefan Boltzmann constant = 5.67 e-8 J s-1 m-2 K-4

    Parameters
    ----------
    t_air_k_24 : float
        daily air temperature in Kelvin
        :math:`T_{a,K}`
        [-]
    vp_24 : float
        daily vapour pressure
        :math:`e_{a}`
        [mbar]
    trans_24 : float
        daily atmospheric transmissivity
        :math:`\tau`
        [-]

    Returns
    -------
    l_net : float
        daily net longwave radiation
        :math:`L^{*}`
        [Wm-2]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.longwave_radiation_fao_etref(t_air_k=302.5, vp=10.3, trans_24=0.6)
    68.594182173686306
    """

    vp_slope=0.14
    vp_offset=0.34
    lw_slope=1.35
    lw_offset=-0.35
    return longwave_radiation_fao(t_air_k_24, vp_24, trans_24, vp_slope, vp_offset, lw_slope, lw_offset)


def longwave_radiation_fao(t_air_k_24, vp_24, trans_24, vp_slope=0.14, vp_offset=0.34,
                           lw_slope=1.35, lw_offset=-0.35):
    r"""
    Computes the net longwave radiation according to the FAO 56 manual.

    .. math ::
        L^{*}=\sigma\left(T_{a,K}\right)^{4}
        \left(vp_{off}-vp_{slope}\sqrt{0.1e_{a}}
        \right)\left(lw_{slope}\frac{\tau}{0.75}+lw_{off}\right)

    where the following constant is used

    * :math:`\sigma` = Stefan Boltzmann constant = 5.67 e-8 J s-1 m-2 K-4

    Parameters
    ----------
    t_air_k_24 : float
        daily air temperature in Kelvin
        :math:`T_{a,K}`
        [-]
    vp_24 : float
        daily vapour pressure
        :math:`e_{a}`
        [mbar]
    trans_24 : float
        daily atmospheric transmissivity
        :math:`\tau`
        [-]
    vp_slope : float
        slope of the vp-term in the FAO-56 longwave radiation relationship
        :math:`vp_{slope}`
        [-]
    vp_offset : float
        offset of the vp-term in the FAO-56 longwave radiation relationship
        :math:`vp_{off}`
        [-]
    lw_slope : float
        slope of the tau-term in the FAO-56 longwave radiation relationship
        :math:`lw_{slope}`
        [-]
    lw_offset : float
        offset of the tau-term in the FAO-56 longwave radiation relationship
        :math:`lw_{off}`
        [-]

    Returns
    -------
    l_net : float
        daily net longwave radiation
        :math:`L^{*}`
        [Wm-2]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.longwave_radiation_fao(t_air_k=302.5, vp=10.3, trans_24=0.6)
    68.594182173686306
    """

    return c.sb*t_air_k_24**4*(vp_offset-vp_slope*np.sqrt(0.1*vp_24))*(lw_offset + lw_slope*(trans_24/0.75))


def net_radiation(r0, ra_24, l_net, int_wm2):
    r"""
    Computes the net radiation

    .. math ::
        Q^{*} = \left[\left(1-\alpha_{0}\right)S^{\downarrow}-L^{*}-I\right]

    Parameters
    ----------
    r0 : float
        albedo
        :math:`\alpha_{0}`
        [-]
    ra_24 : float
        daily solar radiation
        :math:`S^{\downarrow}`
        [Wm-2]
    l_net : float
        daily net longwave radiation
        :math:`L^{*}`
        [wm-2]
    int_wm2 : float
        interception
        :math:`I`
        [Wm-2]

    Returns
    -------
    rn_24 : float
        daily net radiation
        :math:`Q^{*}`
        [Wm-2]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.net_radiation(r0=0.10, ra_24=123., l_net=24., int_wm2=0)
    86.7
    """
    return (1-r0)*ra_24-l_net-int_wm2


def net_radiation_canopy(rn_24, sf_soil):
    r"""
    Computes the net radiation for the canopy

    .. math ::
        Q^{*}_{canopy} = \left(1-s_f\right) Q^{*}

    Parameters
    ----------
    rn_24 : float
        net radiation
        :math:`Q^{*}`
        [Wm-2]
    sf_soil : float
        soil fraction
        :math:`s_f`
        [-]

    Returns
    -------
    rn_24_canopy : float
        net radiation for the canopy
        :math:`Q^{*}_{canopy}`
        [Wm-2]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.net_radiation_canopy(rn_24=200, sf_soil=0.4)
    120.0

    """
    return rn_24 * (1-sf_soil)


def net_radiation_soil(rn_24, sf_soil):
    """
    Computes the net radiation for the soil

    .. math ::
        Q^{*}_{soil} = s_f Q^{*}

    Parameters
    ----------
    rn_24 : float
        net radiation
        :math:`Q^{*}`
        [Wm-2]
    sf_soil : float
        soil fraction
        :math:`s_f`
        [-]

    Returns
    -------
    rn_24_soil : float
        net radiation for the soil
        :math:`Q^{*}_{soil}`
        [Wm-2]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.net_radiation_soil(rn_24=200, sf_soil=0.4)
    80.0
    """
    return rn_24 * sf_soil


def net_radiation_grass(ra_24, l_net, r0_grass=0.23):
    r"""
    Computes the net radiation for reference grass

    .. math ::
        Q^{*} = \left[\left(1-\alpha_{0, grass}\right)S^{\downarrow}-L^{*}-I\right]

    Parameters
    ----------
    ra_24 : float
        daily solar radiation
        :math:`S^{\downarrow}`
        [Wm-2]
    l_net : float
        daily net longwave radiation
        :math:`L^{*}`
        [wm-2]
    r0_grass : float
        albedo for reference grass
        :math:`\alpha_{0, grass}`
        [-]

    Returns
    -------
    rn_24_grass : float
        daily net radiation for reference grass
        :math:`Q^{*}`
        [Wm-2]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.net_radiation_grass(ra_24=123., l_net=24.)
    70.7
    """
    return (1-r0_grass)*ra_24-l_net


def volumetric_heat_capacity(se_top=1.0, porosity=0.4):
    r"""
    Computes the volumetric heat capacity of the soil

    .. math ::
        \rho c_{p}=10e^{6}\left[\left(1-\phi\right)^{2}+
        2.5\phi+4.2\phi S_{e,top}\right]

    Parameters
    ----------
    se_top : float
        effective saturation of the topsoil
        :math:`S_{e,top}`
        [-]
    porosity : float
        porosity of the soil
        :math:`\phi`
        [-]

    Returns
    -------
    vhc : float
        volumetric heat capacity
        :math:`\rho c_{p}`
        [J m-3 K-1]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.volumetric_heat_capacity(se_top=0.4, porosity = 0.5)
    23400000.0
    """
    return ((1-porosity)**2+2.5*porosity+4.2*porosity*se_top)*10**6


def soil_thermal_conductivity(se_top):
    r"""
    Computes the soil thermal conductivity

    .. math ::
        k=0.15+18.5S_{e,top}

    Parameters
    ----------
    se_top : float
        effective saturation of the topsoil
        :math:`S_{e,top}`
        [-]

    Returns
    -------
    stc : float
        soil thermal conductivity
        :math:`k`
        [W m-1 K-1]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.soil_thermal_conductivity(se_top=0.4)
    0.8900000000000001
    """
    return 0.15 + 1.85 * se_top


def damping_depth(stc, vhc):
    r"""
    Computes the damping depth

    .. math ::
        z_{d}=\sqrt{\frac{2kP}{2\pi\rho c_{p}}}

    with the following constant

    * :math:`P` period (seconds within a year)

    Parameters
    ----------
    stc : float
        soil thermal conductivity
        :math:`k`
        [W m-1 K-1]
    vhc : float
        volumetric heat capacity
        :math:`\rho c_{p}`
        [J m-3 K-1]

    Returns
    -------
    dd : float
        damping depth
        :math:`z_{d}`
        [m]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.damping_depth(stc=0.9, vhc=volumetric_heat_capacity())
    0.54514600029013294
    """
    return np.sqrt((2*stc*c.year_sec)/(vhc*2*np.pi))


#TODO north-south transition with regard to latitude
def bare_soil_heat_flux(doy, dd, stc, t_amp_year, lat):
    r"""
    Computes the bare soil heat flux

    .. math ::
        G_{0}=\frac{\sqrt{2}A_{t,year}k\sin\left(\frac{2\pi J}{P}-
              \frac{\pi}{4}\right)}{z_{d}}

    where the following constant is used

    * :math:`P` period (seconds within a year)

    The term :math:`-\frac{\pi}{4}` is a phase shift for northern latitudes.
    For southern latitudes the phase shift will be :math:`-\frac{\pi}{4}+\pi`

    Parameters
    ----------
    stc : float
        soil thermal conductivity
        :math:`k`
        [W m-1 K-1]
    dd : float
        damping depth
        :math:`z_{d}`
        [m]
    t_amp_year : float
        yearly air temperature amplitude
        :math:`A_{t,year}`
        [m]
    doy : float
        julian day of the year
        :math:`J`
        [-]
    lat : float
        latitude
        :math:`\lambda`
        [rad]

    Returns
    -------
    g0_bs : float
        bare soil heat flux
        :math:`G_{0}`
        [m]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> stc = rad.soil_thermal_conductivity(se_top=1.0)
    >>> vhc = rad.volumetric_heat_capacity(se_top=1.0)
    >>> dd = damping_depth(stc,vhc)
    >>> rad.bare_soil_heat_flux(126, dd, stc, t_amp_year=13.4, lat=40*(math.pi/180.0))
    array([ 45.82350561])
    """

    phase = np.zeros_like(lat)
    phase = np.where(lat > 0, -np.pi/4.0, phase)
    phase = np.where(lat <= 0, -np.pi/4.0+np.pi, phase)

    out = (np.sqrt(2.0)*t_amp_year*stc*np.sin(2*np.pi/c.year_sec*doy*c.day_sec+phase))/dd

    return out


def soil_heat_flux(g0_bs, sf_soil, land_mask, rn_24_soil, trans_24, ra_24, l_net, rn_slope=0.92, rn_offset=-61.0):
    r"""
    Computes the soil heat flux

    .. math ::
        G=s_f G_{0}

    Parameters
    ----------
    g0_bs : float
        bare soil heat flux
        :math:`G_{0}`
        [W m-2]
    sf_soil : float
        soil fraction
        :math:`s_f`
        [-]
    land_mask : int
        land use classification
        :math:`l`
        [-]
    rn_24_soil : float
        net radiation for the soil
        :math:`Q^{*}_{soil}`
        [Wm-2]
    trans_24 : float
        daily atmospheric transmissivity
        :math:`\tau`
        [-]
    rn_slope : float
        slope rn/g0 relation water
        :math:`lws`
        [-]
    rn_offset : float
        offset rn/g0 relation water
        :math:`lwo`
        [-]
    ra_24 : float
        daily solar radiation
        :math:`S^{\downarrow}`
        [Wm-2]
    l_net : float
        daily net longwave radiation
        :math:`L^{*}`
        [wm-2]

    Returns
    -------
    g0_24 : float
        daily soil heat flux
        :math:`G`
        [W m-2]

    Examples
    --------
    >>> import ETLook.radiation as rad
    >>> rad.soil_heat_flux(g0_bs=12.4, sf_soil=0.4)
    4.960000000000001
    """
    def land_city_func(g0_bs, sf_soil):
        return g0_bs * sf_soil

    def water_func(ra_24, trans_24, l_net, rn_slope, rn_offset, rn_24_soil):
        rn_24_clear = 0.95 * ra_24 / trans_24 - l_net
        g0_24_clear = rn_24_clear * rn_slope + rn_offset
        g0_24_clear = np.minimum(g0_24_clear, 0.5 * rn_24_clear)

        # adjust water heat storage to current net radiation conditions
        g0_24 = g0_24_clear * rn_24_soil / rn_24_clear

        return g0_24

    g0 = np.zeros_like(land_mask)
    g0 = np.where(land_mask == 1, land_city_func(g0_bs, sf_soil), g0)
    g0 = np.where(land_mask == 2, water_func(ra_24, trans_24, l_net, rn_slope, rn_offset, rn_24_soil), g0)
    g0 = np.where(land_mask == 3, land_city_func(g0_bs, sf_soil), g0)

    return g0
