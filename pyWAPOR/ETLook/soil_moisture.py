"""
    The soil_moisture module contains all functions related to soil moisture data components.

"""
import math
from pyWAPOR.ETLook import constants as c
from pyWAPOR.ETLook import unstable
import numpy as np

def wet_bulb_temperature_inst(t_air_i, t_dew_i):
    r"""
    Computes the instantaneous wet bulb temperature.

    Parameters
    ----------
    t_air_i : float
        instantaneous air temperature
        :math:`T_{a}`
        [C]
    t_dew_i : float
        instantaneous dew point temperature
        :math:`Td_{a}`
        [C]

    Returns
    -------
    t_wet_i : float
        instantaneous wet bulb temperature
        :math:`Tw_{a}`
        [C]
    """
    tw = wetbulb_temperature_iter(t_air_i, t_dew_i)

    return tw


def dew_point_temperature_inst(vp_i):
    r"""
    Computes the instantaneous dew point temperature.

    Parameters
    ----------
    vp_i : float
        instantaneous vapour pressure
        :math:`e_{a}`
        [mbar]

    Returns
    -------
    t_dew_i : float
        instantaneous dew point temperature
        :math:`Td_{a}`
        [K]
    """
    t_dew_i = (237.3 * np.log(vp_i / 6.108)) / (17.27 - np.log(vp_i / 6.108))

    return t_dew_i


def dew_point_temperature_coarse_inst(vp_i):
    r"""
    Computes the instantaneous dew point temperature.

    Parameters
    ----------
    vp_i : float
        instantaneous vapour pressure
        :math:`e_{a}`
        [mbar]

    Returns
    -------
    t_dew_coarse_i : float
        instantaneous dew point temperature
        :math:`Td_{a}`
        [K]
    """
    t_dew_i = (237.3 * np.log(vp_i / 6.108)) / (17.27 - np.log(vp_i / 6.108))

    return t_dew_i


# only used internally
def latent_heat_iter(t):
    # temperature in celcius
    lv = 1000 * (2501 - 2.361 * t)

    return lv


# only used internally
def psychometric_constant_iter(lv, p=1013.25, cp=1004, rm=0.622):
    # lv: latent heat of vaporization (J/s?)
    # p:  pressure in hPa
    # cp: specific heat
    # rm: ratio of molecular weight

    psy = (cp * p) / (lv * rm)

    return psy


# only used internally
def vapor_pressure_iter(t):
    # t: temperature in Celcius

    vp = 6.108 * np.exp((17.27 * t) / (237.3 + t))

    return vp


# only used internally
def wetbulb_temperature_iter(ta, td):
    
    maxiter = 1000
    tol = 1e-6
    pressure = 1013.25

    lv = latent_heat_iter(ta)
    psy = psychometric_constant_iter(lv, p=pressure)

    tw = td + ((ta - td) / 3)
    ea_ta = vapor_pressure_iter(td)

    n = 0

    prev_dir = np.zeros(ta.shape)
    step = (ta - td) / 5.

    while abs(np.nanmax(step)) > tol:

        ea_tw = vapor_pressure_iter(tw) - psy * (ta - tw)

        direction = (-1) ** ((ea_tw - ea_ta) > 0)

        step = np.where(prev_dir != direction, step * 0.5, step)

        tw += step * direction
        prev_dir = direction

        n += 1

        if n >= maxiter:
            return tw

    return tw


# some internal functions for the stability correction based on Brutsaert
def psi_m(y):
    r"""
    Computes the stability correction for momentum based on
    Brutsaert (1999) [2]_.

    .. math ::
        \Psi_{M}(y)=\ln(a+y)-3by^{\frac{1}{3}}+ \\
        \frac{ba^{\frac{1}{3}}}{2}\ln[\frac{(1+x)^{2}}{(1-x+x^{2})}]+\\
        \sqrt{3}ba^{\frac{1}{3}}\arctan[\frac{(2x-1)}{\sqrt{3}}]+\Psi_{0}

    where the following constants are used

    * :math:`a` = 0.33
    * :math:`b` = 0.41

    in which

    .. math ::
        x = (\frac{y}{a})^{\frac{1}{3}}

    and

    .. math ::
        y = \frac{-(z-d)}{L}

    where :math:`L` is the monin obukhov length defined by
    :func:`ETLook.unstable.monin_obukhov_length`,
    :math:`z` and :math:`d` are the
    measurement height and displacement height respectively. All aforementioned
    parameters are different for the bare soil and full canopy solutions.

    The symbol :math:`\Psi_{0}` denotes a constant of integration, given by

    .. math ::
       \Psi_{0}=-\ln{a}+\sqrt{3}ba^{\frac{1}{3}}\frac{\pi}{6}

    .. plot:: pyplots/soil_moisture/plot_psi_m.py

    Notes
    -----
    This function should not be used as an input function for a ETLook tool.
    This function is used internally by :func:`aerodynamical_resistance_bare`
    and :func:`aerodynamical_resistance_full` and :func:`wind_speed_soil`.

    References
    ----------
    .. [2] Brutsaert, W., Aspect of bulk atmospheric boundary layer similarity
        under free-convective conditions,
        Reviews of Geophysics, 1999, 37(4), 439-451.
    """
    a = 0.33
    b = 0.41
    x = (y / a) ** (1. / 3.)
    phi_0 = -np.log(a) + np.sqrt(3) * b * a ** (1. / 3.) * np.pi / 6.
    res = (
        np.log(a + y)
        - 3 * b * y ** (1. / 3.)
        + (b * a ** (1. / 3.)) / 2. * np.log((1 + x) ** 2 / (1 - x + x ** 2))
        + np.sqrt(3) * b * a ** (1. / 3.) * np.arctan((2 * x - 1) / np.sqrt(3))
        + phi_0
    )
    return res


def psi_h(y):
    r"""
    Computes the stability correction for momentum based on
    Brutsaert (1999) [2]_.

    .. math ::
        \Psi_{H}(y)=[\frac{(1-d)}{n}]\ln{\frac{(c+y^n)}{c}}

    where the following constants are used

    * :math:`c` = 1.00
    * :math:`d` = 0.057
    * :math:`n` = 0.78

    in which

    .. math ::
        y = \frac{-(z-d)}{L}

    where :math:`L` is the monin obukhov length defined by
    :func:`ETLook.unstable.monin_obukhov_length`,
    :math:`z` and :math:`d` are the
    measurement height and displacement height respectively. All aforementioned
    parameters are different for the bare soil and full canopy solutions.

    .. plot:: pyplots/soil_moisture/plot_psi_h.py

    Notes
    -----
    This function should not be used as an input function for a tool.
    This function is used internally by :func:`aerodynamical_resistance_bare`
    and :func:`aerodynamical_resistance_full` and :func:`wind_speed_soil`.

    References
    ----------
    .. [2] Brutsaert, W., Aspect of bulk atmospheric boundary layer similarity
        under free-convective conditions,
        Reviews of Geophysics, 1999, 37(4), 439-451.
    """
    c = 0.33
    d = 0.057
    n = 0.78
    return ((1 - d) / n) * np.log((c + y ** n) / c)


def initial_friction_velocity_inst(u_b_i, z0m, disp, z_b=100):
    r"""
    Computes the initial instantaneous friction velocity without stability
    corrections.

    .. math ::
        u_{*}=\frac{ku_{b}}{\ln\left(\frac{z_{b}-d}{z_{0,m}}\right)}

    Parameters
    ----------
    u_b_i : float
        instantaneous wind speed at blending height
        :math:`u_{b}`
        [m s-1]
    z0m : float
        surface roughness
        :math:`z_{0,m}`
        [m]
    disp : float
        displacement height
        :math:`d`
        [m]
    z_b : float
        blending height
        :math:`z_{b}`
        [m]


    Returns
    -------
    u_star_i_init : float
        initial estimate of the instantaneous friction velocity
        :math:`u_{*,i}`
        [m s-1]
    """
    return (c.k * u_b_i) / (np.log((z_b - disp) / z0m))


def atmospheric_emissivity_inst(vp_i, t_air_k_i):
    r"""
    Computes the atmospheric emissivity according to Brutsaert [1]_.

    .. math ::
        \varepsilon_{a}=a\left(\frac{e_{a}}{T_{a}}\right)^{b}

    where the following constants are used

    * :math:`a` = 1.24
    * :math:`b` = 1/7

    Parameters
    ----------
    vp_i : float
        instantaneous vapour pressure
        :math:`e_{a}`
        [mbar]
    t_air_k_i : float
        instantaneous air temperature
        :math:`T_{a}`
        [K]

    Returns
    -------
    emiss_atm_i : float
        instantaneous atmospheric emissivity
        :math:`\varepsilon_{a}`
        [-]

    References
    ----------
    .. [1] Brutsaert, W., On a derivable formula for long-wave radiation
        from clear skies, Water Resour. Res, 1975, 11, 742-744.
    """
    return 1.24 * (vp_i / t_air_k_i) ** (1. / 7.)


def net_radiation_bare(ra_hor_clear_i, emiss_atm_i, t_air_k_i, lst, r0_bare=0.38):
    r"""
    Computes the net radiation for the bare soil with zero evaporation

    .. math ::

        Q_{bare}^{*}=\left(1-\alpha_{0,bare}\right)S_{d}+\varepsilon_{s}\varepsilon_{a}\sigma T_{a}^{4}-\varepsilon_{s}\sigma T_{s}^{4}

    Parameters
    ----------
    ra_hor_clear_i : float
        Total clear-sky irradiance on a horizontal surface
        :math:`S_{d}`
        [W/m2]
    emiss_atm_i : float
        instantaneous atmospheric emissivity
        :math:`\varepsilon_{a}`
        [-]
    t_air_k_i : float
        instantaneous air temperature
        :math:`T_{a}`
        [K]
    lst : float
        surface temperature
        :math:`T_{0}`
        [K]
    r0_bare : float
        dry bare soil surface albedo
        :math:`\alpha_{0, bare}`
        [-]


    Returns
    -------
    rn_bare : float
        net radiation bare soil
        :math:`Q^*_{bare}`
        [Wm-2]
    """
    emiss_bare = 0.95
    rn_bare = (
        (1 - r0_bare) * ra_hor_clear_i
        + emiss_atm_i * emiss_bare * c.sb * (t_air_k_i) ** 4
        - emiss_bare * c.sb * (lst) ** 4
    )

    return rn_bare


def net_radiation_full(ra_hor_clear_i, emiss_atm_i, t_air_k_i, lst, r0_full=0.18):
    r"""
    Computes the net radiation at full canopy with zero evaporation

    .. math ::

        Q_{full}^{*}=\left(1-\alpha_{0,full}\right)S_{d}+\varepsilon_{c}\varepsilon_{a}\sigma T_{a}^{4}-\varepsilon_{c}\sigma T_{s}^{4}

    Parameters
    ----------
    ra_hor_clear_i : float
        Total clear-sky irradiance on a horizontal surface
        :math:`ra_hor_clear_i`
        [W/m2]
    emiss_atm_i : float
        instantaneous atmospheric emissivity
        :math:`P`
        [-]
    t_air_k_i : float
        instantaneous air temperature
        :math:`T_{a}`
        [K]
    lst : float
        surface temperature
        :math:`T_{0}`
        [K]
    r0_full : float
        surface albedo full vegetation
        :math:`\alpha_{0, full}`
        [-]

    Returns
    -------
    rn_full : float
        net radiation full vegetation
        :math:`Q^*_{full}`
        [Wm-2]
    """
    emiss_full = 0.99
    rn_full = (
        (1 - r0_full) * ra_hor_clear_i
        + emiss_atm_i * emiss_full * c.sb * (t_air_k_i) ** 4
        - emiss_full * c.sb * (lst) ** 4
    )

    return rn_full


def sensible_heat_flux_bare(rn_bare, fraction_h_bare=0.65):
    r"""
    Computes the bare soil sensible heat flux

    .. math ::

        H_{bare} = H_{f, bare}Q^*_{bare}

    Parameters
    ----------
    rn_bare : float
        net radiation bare soil
        :math:`Q^*_{bare}`
        [Wm-2]
    fraction_h_bare : float
        fraction of H of net radiation bare soil
        :math:`H_{f, bare}`
        [-]

    Returns
    -------
    h_bare : float
        sensible heat flux bare soil
        :math:`H_{bare}`
        [Wm-2]
     """
    return rn_bare * fraction_h_bare


def sensible_heat_flux_full(rn_full, fraction_h_full=0.95):
    r"""
    Computes the full canopy sensible heat flux

    .. math ::

        H_{full} = H_{f, full}Q^*_{full}

    Parameters
    ----------
    rn_full : float
        net radiation full vegetation
        :math:`Q^*_{full}`
        [Wm-2]
    fraction_h_full : float
        fraction of H of net radiation full vegetation
        :math:`H_{f, full}`
        [-]

    Returns
    -------
    h_full : float
        sensible heat flux full vegetation
        :math:`H_{full}`
        [Wm-2]
     """
    return rn_full * fraction_h_full


def wind_speed_blending_height_bare(u_i, z0m_bare=0.001, z_obs=10, z_b=100):
    r"""
    Computes the wind speed at blending height :math:`u_{b}` [m/s] using the
    logarithmic wind profile

    .. math ::
        u_{b}=\frac{u_{obs}\ln\left(\frac{z_{b}}{z_{0,m}}\right)}
        {\ln\left(\frac{z_{obs}}{z_{0,m}}\right)}

    Parameters
    ----------
    u_i : float
        instantaneous wind speed at observation height
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
    z0m_bare : float
        surface roughness bare soil
        :math:`z_{0,m}`
        m

    Returns
    -------
    u_b_i_bare : float
        instantaneous wind speed at blending height for bare soil
        :math:`u_{b,i,bare}`
        [m/s]
    """
    ws = (c.k * u_i) / np.log(z_obs / z0m_bare) * np.log(z_b / z0m_bare) / c.k
    ws = ws.clip(1,150)
    return ws


def wind_speed_blending_height_full_inst(u_i, z0m_full=0.1, z_obs=10, z_b=100):
    r"""
    Computes the wind speed at blending height :math:`u_{b}` [m/s] using the
    logarithmic wind profile

    .. math ::
        u_{b}=\frac{u_{obs}\ln\left(\frac{z_{b}}{z_{0,m}}\right)}
        {\ln\left(\frac{z_{obs}}{z_{0,m}}\right)}

    Parameters
    ----------
    u_i : float
        instantaneous wind speed at observation height
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
    z0m_full : float
        surface roughness vegetation
        :math:`z_{0,m}`
        [m]

    Returns
    -------
    u_b_i_full : float
        instantaneous wind speed at blending height for full vegetation
        :math:`u_{b,i,full}`
        [m s-1]
    """
    ws = (c.k * u_i) / np.log(z_obs / z0m_full) * np.log(z_b / z0m_full) / c.k
    ws = ws.clip(1, 150)

    return ws


def friction_velocity_full_inst(u_b_i_full, z0m_full=0.1, disp_full=0.667, z_b=100):
    r"""
    Like :func:`initial_friction_velocity_inst` but with full vegetation parameters

    Parameters
    ----------
    u_b_i_full : float
        instantaneous wind speed blending height for full vegetation
        :math:`u_{b,d}`
        [m s-1]
    z0m_full : float
        surface roughness vegetation
        :math:`z_{0,m,b}`
        [m]
    disp_full : float
        displacement height vegetation
        :math:`d^{b}`
        [m]
    z_b : float
        blending height
        :math:`z_b`
        [m]

    Returns
    -------
    u_star_i_full : float
        instantaneous friction velocity vegetation
        :math:`u_{f}^{*}`
        [m s-1]

    """
    return initial_friction_velocity_inst(u_b_i_full, z0m_full, disp_full, z_b=100)


def friction_velocity_bare_inst(u_b_i_bare, z0m_bare=0.001, disp_bare=0.0, z_b=100):
    r"""
    Like :func:`initial_friction_velocity_inst` but with bare soil parameters

    Parameters
    ----------
    u_b_i_bare : float
        instantaneous wind speed blending height bare soil
        :math:`u_{b,d}`
        [W m-2]
    z0m_bare : float
        surface roughness bare soil
        :math:`z_{0,m,b}`
        [m]
    disp_bare : float
        displacement height bare soil
        :math:`d^{b}`
        [m]
    z_b : float
        blending height
        :math:`z_b`
        [m]

    Returns
    -------
    u_star_i_bare : float
        instantaneous friction velocity bare soil
        :math:`u_{b}^{*}`
        [m s-1]

    """
    return initial_friction_velocity_inst(u_b_i_bare, z0m_bare, disp_bare, z_b=100)


def monin_obukhov_length_bare(h_bare, ad_i, u_star_i_bare, t_air_k_i):
    r"""
    Like :func:`unstable.monin_obukhov_length` but with bare soil parameters

    Parameters
    ----------
    h_bare : float
        sensible heat flux for dry bare soil
        :math:`H_{b,d}`
        [W m-2]
    ad_i : float
        instantaneous air density
        :math:`\rho`
        [k g m-3]
    u_star_i_bare : float
        instantaneous friction velocity bare soil
        :math:`u^{*}_{b}`
        [m s-1]
    t_air_k_i : float
        instantaneous air temperature
        :math:`T_{a}`
        [K]

    Returns
    -------
    L_bare : float
        monin obukhov length dry vegetation
        :math:`L_{b,d}`
        [m]

    """
    return unstable.monin_obukhov_length(h_bare, ad_i, u_star_i_bare, t_air_k_i)


def monin_obukhov_length_full(h_full, ad_i, u_star_i_full, t_air_k_i):
    r"""
    Like :func:`unstable.monin_obukhov_length` but with full canopy parameters

    Parameters
    ----------
    h_full : float
        sensible heat flux for dry full vegetation
        :math:`H_{f,d}`
        [W m-2]
    ad_i : float
        instantaneous air density
        :math:`\rho`
        [k g m-3]
    u_star_i_full : float
        instantaneous friction velocity vegetation
        :math:`u^{*}_{b}`
        [m s-1]
    t_air_k_i : float
        instantaneous air temperature
        :math:`T_{a}`
        [K]

    Returns
    -------
    L_full : float
        monin obukhov length dry vegetation
        :math:`L_{f,d}`
        [m]

    """
    return unstable.monin_obukhov_length(h_full, ad_i, u_star_i_full, t_air_k_i)


def aerodynamical_resistance_full(u_i, L_full, z0m_full=0.1, disp_full=0.667, z_obs=10):
    r"""
    Computes the aerodynamical resistance for a full canopy.

    .. math ::
        z_{1} = \frac{z_{obs}-d}{z_{0,m}}

        z_{2} = \frac{z_{obs}-d}{L}

        z_{3} = \frac{z_{0,m}}{L}

        z_{4} = \frac{z_{obs}-d}{\frac{z_{0,m}}{7}}

        z_{5} = \frac{\frac{z_{0,m}}{7}}{L}

        r_{a,c}=\frac{(\ln(z_{1})-\phi_{m}(-z_{2})+\phi_{m}(-z_{3}))(\ln(z_{4})-\phi_{h}(-z_{2})+\phi_{h}(-z_{5}))}{k^{2}u}

    Parameters
    ----------
    u_i : float
        instantaneous wind speed at observation height
        :math:`u_{obs}`
        [m/s]
    z_obs : float
        observation height of wind speed
        :math:`z_{obs}`
        [m]
    disp_full : float
        displacement height
        :math:`d`
        [m]
    z0m_full : float
        surface roughness
        :math:`z_{0,m}`
        [m]
    L_full : float
        monin obukhov length
        :math:`L`
        [m]

    Returns
    -------
    rac : float
        aerodynamical resistance canopy
        :math:`r_{a,c}`
        [sm-1]

    """
    z1 = (z_obs - disp_full) / z0m_full
    z2 = (z_obs - disp_full) / L_full
    z3 = z0m_full / L_full
    z4 = (z_obs - disp_full) / (z0m_full / 7)
    z5 = (z0m_full / 7) / L_full
    res = (
        (np.log(z1) - psi_m(-z2) + psi_m(-z3)) * (np.log(z4) - psi_h(-z2) + psi_h(-z5))
    ) / (c.k ** 2 * u_i)

    return res


def aerodynamical_resistance_bare(u_i, L_bare, z0m_bare=0.001, disp_bare=0.0, z_obs=10):
    r"""
    Computes the aerodynamical resistance for a dry bare soil.

    .. math ::
        z_{1} = \frac{z_{obs}-d}{z_{0,b,m}}

        z_{2} = \frac{z_{obs}-d}{L_{b}}

        r_{a,a}=\frac{(\ln(z_{1})-\phi_{m}(-z_{2}))(\ln(z_{1})-\phi_{h}(-z_{2}))}{k^{2}u}

    Parameters
    ----------
    u_i : float
        instantaneous wind speed at observation height
        :math:`u_{obs}`
        [m/s]
    z_obs : float
        observation height of wind speed
        :math:`z_{obs}`
        [m]
    disp_bare : float
        displacement height
        :math:`d`
        [m]
    z0m_bare : float
        surface roughness
        :math:`z_{0,b,m}`
        [m]
    L_bare : float
        monin obukhov length
        :math:`L_{b}`
        [m]

    Returns
    -------
    raa : float
        aerodynamical resistance dry surface
        :math:`r_{a,a}`
        [sm-1]

    """

    z1 = (z_obs - disp_bare) / z0m_bare
    z2 = (z_obs - disp_bare) / L_bare
    res = ((np.log(z1) - psi_m(-z2)) * (np.log(z1) - psi_h(-z2))) / (c.k ** 2 * u_i)

    return res


def wind_speed_soil_inst(u_i, L_bare, z_obs=10):
    r"""
    Computes the instantaneous wind speed at soil surface

    .. math ::

        u_{i,s}=u_{obs}\frac{\ln\left(\frac{z_{obs}}{z_{0}}\right)}
              {\ln\left(\frac{z_{obs}}{z_{0,s}}\right)-\psi_{m}\left(\frac{-z_{0}}{L}\right)}

    Parameters
    ----------
    u_i : float
        wind speed at observation height
        :math:`u_{obs}`
        [m/s]
    z_obs : float
        observation height of wind speed
        :math:`z_{obs}`
        [m]
    L_bare : float
        monin obukhov length
        :math:`L`
        [m]

    Returns
    -------
    u_i_soil : float
        instantaneous wind speed just above soil surface
        :math:`u_{i,s}`
        [ms-1]

    """
    z0_soil = 0.01
    z0_free = 0.1
    return u_i * (
        (np.log(z0_free / z0_soil))
        / (np.log(z_obs / z0_soil) - psi_m(-z0_free / L_bare))
    )


def aerodynamical_resistance_soil(u_i_soil):
    r"""
    Computes the aerodynamical resistance of the soil

    .. math ::
        r_{a,s}=\frac{1}{\left(0.0025T_{dif}^{\frac{1}{3}}+0.012u_{i,s}\right)}

    Parameters
    ----------
    u_i_soil : float
        instantaneous wind speed just above soil surface
        :math:`u_{i,s}`
        [m s-1]

    Returns
    -------
    ras : float
        aerodynamical resistance
        :math:`r_{a,s}`
        [sm-1]

    """
    Tdif = 10.0
    return 1. / (0.0025 * (Tdif) ** (1. / 3.) + 0.012 * u_i_soil)


def maximum_temperature_full(
    ra_hor_clear_i, emiss_atm_i, t_air_k_i, ad_i, rac, r0_full=0.18
):
    r"""
    Computes the maximum temperature under fully vegetated conditions

    .. math ::

        T_{c,max}=\frac{\left(1-\alpha_{c}\right)S_{d}+\varepsilon_{c}\varepsilon_{a}\sigma
                  T_{a}^{4}-\varepsilon_{c}\sigma T_{a}^{4}}{4\varepsilon_{s}\sigma
                  T_{a}^{3}+\rho C_{p}/r_{a,c}}+T_{a}

    Parameters
    ----------
    ra_hor_clear_i : float
        Total clear-sky irradiance on a horizontal surface
        :math:`ra_hor_clear_i`
        [W/m2]
    emiss_atm_i : float
        instantaneous atmospheric emissivity
        :math:`P`
        [-]
    t_air_k_i : float
        instantaneous air temperature
        :math:`T_{a}`
        [K]
    rac : float
        aerodynamic resistance canopy
        :math:`r_{a,c}`
        [sm-1]
    ad_i : float
        instantaneous air density
        :math:`\rho`
        [kg m-3]
    r0_full : float
        surface albedo full vegetation cover
        :math:`\alpha_{0, full}`
        [-]

    Returns
    -------
    t_max_full : float
        maximum temperature at full vegetation cover
        :math:`T_{c,max}`
        [K]

    """
    emiss_full = 0.99

    tc_max_num = (
        (1 - r0_full) * ra_hor_clear_i
        + emiss_full * emiss_atm_i * c.sb * (t_air_k_i) ** 4
        - emiss_full * c.sb * (t_air_k_i) ** 4
    )
    tc_max_denom = 4 * emiss_full * c.sb * (t_air_k_i) ** 3 + (ad_i * c.sh) / rac
    tc_max = tc_max_num / tc_max_denom + t_air_k_i

    return tc_max


def maximum_temperature_bare(
    ra_hor_clear_i, emiss_atm_i, t_air_k_i, ad_i, raa, ras, r0_bare=0.38
):
    r"""
    Computes the maximum temperature under dry bare soil conditions

    .. math ::

        T_{s,max}=\frac{\left(1-\alpha_{s}\right)S_{d}+\varepsilon_{s}\varepsilon_{a}\sigma
        T_{a}^{4}-\varepsilon_{s}\sigma T_{a}^{4}}{4\varepsilon_{s}\sigma T_{a}^{3}+
        \rho C_{p}/\left[\left(r_{a,a}+r_{a,s}\right)\left(1-G/R_{n,s}\right)\right]}+T_{a}

    Parameters
    ----------
    ra_hor_clear_i : float
        Total clear-sky irradiance on a horizontal surface
        :math:`ra_hor_clear_i`
        [W/m2]
    emiss_atm_i : float
        instantaneous atmospheric emissivity
        :math:`P`
        [-]
    t_air_k_i : float
        instantaneous air temperature
        :math:`T_{a}`
        [K]
    ad_i : float
        instantaneous air density
        :math:`\rho`
        [kg m-3]
    raa : float
        aerodynamical resistance
        :math:`r_{a,a}`
        [sm-1]
    ras : float
        aerodynamical resistance
        :math:`r_{a,s}`
        [sm-1]
    r0_bare : float
        dry bare soil surface albedo
        :math:`\alpha_{0, bare}`
        [-]

    Returns
    -------
    t_max_bare : float
        maximum temperature at bare soil
        :math:`T_{c,max}`
        [K]

    """
    emiss_bare = 0.95
    ts_max_num = (
        (1 - r0_bare) * ra_hor_clear_i
        + emiss_bare * emiss_atm_i * c.sb * (t_air_k_i) ** 4
        - emiss_bare * c.sb * (t_air_k_i) ** 4
    )
    ts_max_denom = 4 * emiss_bare * c.sb * (t_air_k_i) ** 3 + (ad_i * c.sh) / (
        (raa + ras) * (1 - 0.35)
    )
    return ts_max_num / ts_max_denom + t_air_k_i


def maximum_temperature(t_max_bare, t_max_full, vc):
    r"""
    Computes the maximum temperature at dry conditions

    .. math ::

        T_{0,max} = c_{veg}(T_{c,max}-T_{s,max})+T_{s,max}


    Parameters
    ----------
    t_max_bare : float
        maximum temperature at bare soil
        :math:`T_{s,max}`
        [K]
    t_max_full : float
        maximum temperature at full dry vegetation
        :math:`T_{c,max}`
        [K]
    vc : float
        vegetation cover
        :math:`c_{veg}`
        [-]


    Returns
    -------
    lst_max : float
        maximum temperature at dry conditions
        :math:`T_{0,max}`
        [K]

    """
    return vc * (t_max_full - t_max_bare) + t_max_bare


def soil_moisture_from_maximum_temperature(lst_max, lst, t_wet_k_i):
    r"""
    Computes the relative root zone soil moisture based on estimates of
    maximum temperature and wet bulb temperature and measured land
    surface temperature

    .. math ::

        \Theta = \frac{T_{0}-T_{w}}{T_{0,max}-T_{w}}

    Parameters
    ----------
    lst : float
        land surface temperature
        :math:`T_{0}`
        [K]
    lst_max : float
        maximum temperature at dry conditions
        :math:`T_{0,max}`
        [K]
    t_wet_k_i : float
        instantaneous wet bulb temperature
        :math:`T_{w}`
        [K]


    Returns
    -------
    se_root : float
        relative root zone soil moisture (based on LST)
        :math:`\Theta`
        [%]

    """
    ratio = (lst - t_wet_k_i) / (lst_max - t_wet_k_i)
    #ratio[ratio < 0] = 0
    #ratio[ratio > 1] = 1
    ratio = ratio.clip(0, 1)
    
    return 1 - ratio
