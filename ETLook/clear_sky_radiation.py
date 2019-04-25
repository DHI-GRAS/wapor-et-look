"""
    The clear_sky_radiation module contains all functions related to the
    calculation of (instantaneous) clear sky radation. Most of these
    functions are based upon Suri et Hofierka, 2004. [SuHi04]_

"""
from math import pi, sin, cos, asin, exp, log
import numpy as np

def extraterrestrial_irradiance_normal(I0, ied):
    r"""
    Computes the extraterrestrial irradiance normal to the solar beam

    .. math ::
        G_{0}=I_{0}\varepsilon

    Parameters
    ----------
    I0 : float
        solar constant
        :math:`I_{0}`
        [W m\ :sup:`-2`\]
    ied : float
        inverse earth sun distance
        :math:`\varepsilon`
        [AU\ :sup:`-1`\]

    Returns
    -------
    G0 : float
        ext rad normal to solar bea,
        :math:`G_{0}`
        [W m\ :sup:`-2`\]

    """
    return I0*ied


def inverse_earth_sun_distance(day_angle):
    r"""
    Computes the inverse earth sun distance

    .. math ::
        \varepsilon=1+0.03344\cos\left(j^{\prime}-0.048869\right)

    Parameters
    ----------
    day_angle : float
        day angle
        :math:`j^{\prime}`
        [-]

    Returns
    -------
    ied : float
        inverse earth sun distance
        :math:`\varepsilon`
        [AU]

    """
    return 1 + 0.03344 * cos(day_angle - 0.048869)


def day_angle(doy):
    r"""

    Computes the day angle. 0 is january 1\ :sup:`st`\,
    2\ :math:`\pi` is december 31\ :sup:`st`\.

    .. math ::
        j^{\prime}=\frac{2\pi j}{365.25}

    Parameters
    ----------
    doy : float
        day of year
        :math:`j`
        [-]

    Returns
    -------
    day_angle : float
        day angle
        :math:`j^{\prime}`
        [rad]

    """
    return 2 * pi * doy / 365.25


def solar_constant():
    r"""
    Returns the solar constant. The solar constant is defined as the
    flux density of solar radiation at the mean distance from Sun
    to Earth. The solar constant is estimated to be 1367 W m\ :sup:`-2`\

    Parameters
    ----------
       
    Returns
    -------
    I0 : float
        solar constant
        :math:`I_{0}`
        [W m\ :sup:`-2`\]
    """
    return 1367.


def declination(day_angle):
    r"""
    Computes the solar declination. The solar declination is computed according to
    Gruter (1984) [Gr84]_

    .. math ::
       \delta=\arcsin\left(0.3978\sin\left(j^{\prime}-1.4+0.0355\sin
              \left(j^{\prime}-0.0489\right)\right)\right)

    Parameters
    ----------
    day_angle : float
        day angle
        :math:`j^{\prime}`
        [rad]

    Returns
    -------
    decl : float
        declination
        :math:`\delta`
        [rad]

    References
    ----------
    .. [Gr84] Gruter, J. W. (ed), 1984, Radiation Nomenclature, Brussels, CEC, Second
        Solar Energy Programme, Project F, Solar Radiation Data
    """
    return asin(0.3978 * sin(day_angle - 1.4 + 0.0355 * sin(day_angle - 0.0489)))


def relative_optical_airmass(p_air_i, p_air_0_i, h0ref):
    r"""    
    Computes the relative optical air mass. It is calculated according to Kasten and Young 1989 [KaYo89]_

    .. math ::
       m=\left(\frac{p}{p_{0}}\right)/\left(\sin h_{0}^{ref}+0.50572
         \left(h_{0}^{ref}+6.07995\right)^{-1.6364}\right)

    Parameters
    ----------
    p_air_i : float
        actual instantaneous air pressure
        :math:`p`
        [hPa]
    p_air_0_i : float
        air pressure at sea level
        :math:`p_{0}`
        [-]
    h0ref : float
        solar elevation angle corrected for refraction
        :math:`h_{0}^{ref}`
        [degrees]
        
    Returns
    -------
    m : float
        relative optical airmass
        :math:`m`
        [-]

    References
    ----------
    .. [KaYo89] Kasten F., and T.A. Young. 1989, Revised optical air mass tables
       and approximation formula, Applied Optics 28:4735-8
    """
    
    h0ref_rad = h0ref * np.pi/180.
    m = (p_air_i/p_air_0_i)/(np.sin(h0ref_rad) + 0.50572 * (h0ref + 6.07995)**-1.6364)
    m[h0ref_rad <= 0] = 64

    return m


def solar_elevation_angle_refracted(h0):
    r"""    
    Computes the solar elevation angle corrected for refraction

    .. math ::
       h_{0}^{ref} = h_{0} + \Delta h_{0}^{ref}

    where

    .. math ::
        \Delta h_{0}^{ref}=0.61359\left(0.1594+1.123h_{0}+
        0.065656h_{0}^{2}\right)/\left(
        1+28.9344h_{0}+277.3971h_{0}^{2}\right)

    Parameters
    ----------
    h0 : float
        solar elevation angle
        :math:`h_{0}`
        [degrees]
        
    Returns
    -------
    h0ref : float
        solar elevation angle corrected for refrection
        :math:`h_{0}^{ref}`
        [degrees]
    """
    delta_h0ref = 0.061359 * (0.1594 + 1.123*h0 + 0.065656*h0**2)/(1 + 28.9344*h0 + 277.3971*h0**2)
    h0ref = h0 + delta_h0ref
    
    return h0ref


def hour_angle(solar_time):
    r"""    
    Computes the solar hour angle

    .. math ::
       T=\frac{\pi}{12}\left(t-12\right)

    Parameters
    ----------
    solar_time : float
        solar_time
        :math:`t`
        [hours]         
        
    Returns
    -------
    ha : float
        solar hour angle
        :math:`T`
        [rad]

    """
    ha = pi / 12 * (solar_time - 12)
    
    return ha


def solar_elevation_angle(lat, decl, ha):
    r"""    
    Computes the solar elevation angle

    .. math ::
        h_{0}=\arcsin\left(C_{31}\cos T+C_{33}\right)

    where

    .. math ::
        C_{31}=\cos\varphi\cos\delta\;;\;C_{33}=\sin\varphi\sin\delta

    Parameters
    ----------
    lat : float
        latitude
        :math:`\varphi`
        [rad]
    decl : float
        declination
        :math:`\delta`
        [rad]
    ha : float
        solar hour angle
        :math:`T`
        [rad]       
        
    Returns
    -------
    h0 : float
        solar elevation angle
        :math:`h_0`
        [degrees]

    """
    C31 = np.cos(lat) * np.cos(decl)
    C33 = np.sin(lat) * np.sin(decl)

    sin_h0 = C31 * np.cos(ha) + C33
    
    h0_rad = np.arcsin(sin_h0)
    
    h0 = h0_rad * 180 / np.pi
    
    return h0


def rayleigh_optical_thickness(m):
    r"""    
    Computes the Rayleigh optical thickness at airmass :math:`m`. It is calculated according
    to the improved formula by Kasten (1996) [Ka96]_

    if :math:`m` > 20:

    .. math ::

        \delta_{R}(m)=1/\left(6.6296+1.7513m-0.1202m^{2}+0.0065m^{3}-0.00013m^{4}\right)

    if :math:`m` < 20:

    .. math ::

        \delta_{R}(m)=1/\left(10.4+0.718m\right)

    Parameters
    ----------
    m : float
        relative optical airmass
        :math:`m`
        [-]
        
    Returns
    -------
    rotm : float
        Rayleigh optical thickness at airmass m
        :math:`\delta_{R}`
        [-]

    References
    ----------
    .. [Ka96] Kasten F. 1996, The Linke turbidity factor based on improved values of the
       integral Rayleigh optical thickness. Solar Energy 56: 239-44
    """
    rotm = 1 / (10.4 + 0.718 * m)
    rotm[m <= 20] = 1/(6.6296 + 1.7513*m[m <= 20] - 0.1202*m[m <= 20]**2 + 0.0065*m[m <= 20]**3 - 0.00013*m[m <= 20]**4)
    
    return rotm


def beam_irradiance_normal_clear(G0, Tl2, m, rotm, h0):
    r"""    
    Computes the clear sky beam irradiance normal to the solar beam

    .. math ::

        B_{0c}=G_{0}\exp\left(-0.8662T_{LK}m\delta_{R}(m)\right)

    Parameters
    ----------
    G0 : float
        ext rad normal to solar beam
        :math:`G_0`
        [W/m2]
    Tl2 : float
        airmass  2 Linke atmospheric turbidity factor
        :math:`T_{LK}`
        [-]
    m : float
        relative optical airmass
        :math:`m`
        [-]
    rotm : float
        Rayleigh optical thickness at airmass m
        :math:`\delta_{R}`
        [-]
    h0 : float
        solar elevation angle
        :math:`h_0`
        [degrees]
        
    Returns
    -------
    B0c : float
        beam irradiance normal to the solar beam
        :math:`B_{0c}`
        [W/m2]

    """
    B0c = G0 * np.exp(-0.8662 * Tl2 * m * rotm)
    B0c[h0 < 0] = 0

    return B0c


def beam_irradiance_horizontal_clear(B0c, h0):
    r"""    
    Computes the clear sky beam irradiance on a horizontal surface

    .. math ::
        B_{hc}=B_{0c}\sin h_{0}

    Parameters
    ----------
    B0c : float
        beam irradiance normal to the solar beam
        :math:`B_{0c}`
        [W/m2]
    h0 : float
        solar elevation angle
        :math:`h_0`
        [degrees]
        
    Returns
    -------
    Bhc : float
        beam irradiance at a horizontal surface
        :math:`B_{hc}`
        [W/m2]

    """
    Bhc = B0c * np.sin(h0 * np.pi / 180)
    Bhc[h0 < 0] = 0

    return Bhc


def linke_turbidity(wv_i, aod550_i, p_air_i, p_air_0_i):
    r"""    
    Computes the air mass 2 Linke atmospheric turbidity factor

    .. math ::

        p_{rel}=\frac{p}{p_{0}}

        T_{LK}=3.91\tau_{550}\exp\left(0.689p_{rel}\right)+0.376\ln\left(TCWV\right)+\left(2+0.54p_{rel}-0.34p_{rel}^{2}\right)

    Parameters
    ----------
    wv_i : float
        total column atmospheric water vapor
        :math:`TCWV`
        [kg m-2]
    aod550_i : float
        Aerosol optical depth at 550nm
        :math:`aod550`
        [-]
    p_air_i : float
        actual instantaneous air pressure
        :math:`p`
        [hPa]
    p_air_0_i : float
        air pressure at sea level
        :math:`p_0`
        [-]

    Returns
    -------
    Tl2 : float
        Airmass 2 Linke atmospheric turbidity factor
        :math:`T_{LK}`
        [-]

    Examples
    --------
    """
    # prel = p0 / p # Papers mixes p/p0 and p0/p????
    prel = p_air_i / p_air_0_i

    term1 = 3.91 * np.exp(0.689 * prel) * aod550_i
    term2 = 0.376 * np.log(wv_i)
        
    Tl2 = term1 + term2 + (2 + 0.54 * prel - 0.5 * prel**2 + 0.16 * prel**2)    
        
    return Tl2


def diffuse_irradiance_horizontal_clear(G0, Tl2, h0):
    r"""    
    Computes the clear sky beam irradiance on a horizontal surface

    .. math ::

        D_{hc}=G_{0}Tn\left(T_{LK}\right)F_{d}\left(h_{0}\right)

    For the estimation of the transmission function :math:`Tn\left(T_{LK}\right)`
    the following function is used:

    .. math ::

        Tn\left(T_{LK}\right)=-0.015843+0.030543T_{LK}+0.0003797T_{LK}^{2}

    The solar altitude function :math:`F_{d}\left(h_{0}\right)` is evaluated using the expression:

    .. math ::

        F_{d}\left(h_{0}\right)=A_{1}+A_{2}\sin h_{0}+A_{3}sin^{2}h_{0}

    with:

    .. math ::

        A_{1}^{\prime}=0.26463-0.061581T_{LK}+0.0031408T_{LK}^{2}

        A_{1}=0.0022/Tn\left(T_{LK}\right) \: \text{if} \: A_{1}^{\prime}Tn\left(T_{LK}\right)<0.0022

        A_{1}=A_{1}^{\prime} \: \text{if} \: A_{1}^{\prime}Tn\left(T_{LK}\right)\geq0.0022

        A_{2}=2.04020+0.018945T_{LK}-0.011161T_{LK}^{2}

        A_{3}=-1.3025+0.039231T_{LK}+0.0085079T_{LK}^{2}



    Parameters
    ----------
    G0 : float
        ext rad normal to solar beam
        :math:`G_0`
        [W/m2]
    Tl2 : float
        Airmass 2 Linke atmospheric turbidity factor
        :math:`T_{LK}`
        [-]
    h0 : float
        solar elevation angle
        :math:`h_0`
        [degrees]
        
    Returns
    -------
    Dhc : float
        Diffuse irradiance at a horizontal surface
        :math:`D_{hc}`
        [W/m2]

    Examples
    --------
    """
    h0_rad = h0 * np.pi / 180.
    TnTl2 = -0.015843 + 0.030543 * Tl2 + 0.0003797 * Tl2**2

    A1 = 0.26463 - 0.061581 * Tl2 + 0.0031408 * Tl2**2
    A1 = np.where((A1 * TnTl2) < 0.0022, 0.0022 / TnTl2, A1)
    A2 = 2.04020 + 0.018945 * Tl2 - 0.011161 * Tl2**2
    A3 = -1.3025 + 0.039231 * Tl2 + 0.0085079 * Tl2**2

    FdH0 = A1 + A2 * np.sin(h0_rad) + A3 * np.sin(h0_rad)**2

    Dhc = G0 * TnTl2 * FdH0
    
    Dhc[Dhc < 0] = 0
    
    return Dhc


def ra_clear_horizontal(Bhc, Dhc):
    r"""    
    Computes the clear sky beam irradiance on a horizontal surface

    .. math ::
        G_{hc}=B_{hc}+D_{hc}

    Parameters
    ----------
    Bhc : float
        beam irradiance at a horizontal surface
        :math:`B_{hc}`
        [W/m2]
    Dhc : float
        Diffuse irradiance at a horizontal surface
        :math:`D_{hc}`
        [W/m2]
        
    Returns
    -------
    ra_hor_clear_i : float
        Total clear-sky irradiance on a horizontal surface
        :math:`G_{hc}`
        [W/m2]

    Examples
    --------
    """
    ra_clear_hor = Bhc + Dhc
    
    return ra_clear_hor
