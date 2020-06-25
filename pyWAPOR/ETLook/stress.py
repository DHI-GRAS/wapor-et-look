import numpy as np

def stress_radiation(ra_24):
    r"""
    Computes the stress for plants when there is not sufficient radiation

    .. math ::
        S_{r}=\frac{S^{\downarrow}}{\left(S^{\downarrow}+60.\right)}
        \left(1+\frac{60}{500}\right)

    Parameters
    ----------
    ra_24 : float
        daily solar radiation
        :math:`S^{\downarrow}`
        [Wm-2]

    Returns
    -------
    stress_rad : float
        stress factor for radiation
        :math:`S_{r}`
        [-]


    Examples
    --------
    >>> import ETLook.stress as stress
    >>> stress.stress_radiation()
    0.0
    >>> stress.stress_radiation(500)
    1.0
    >>> stress.stress_radiation(700)
    1.0
    >>> stress.stress_radiation(250)
    0.90322580645161288
    """
    stress = ra_24/(ra_24 + 60.)*(1 + 60./500.)
    stress = np.clip(stress, 0, 1)

    return stress


def stress_moisture(se_root, tenacity=1.5):
    r"""
    Computes the stress for plants when there is not sufficient soil
    moisture in the root zone

    .. math ::
        S_{m}=K_{sf}S_{e,root}-\frac{\sin\left(2\pi S_{e,root}\right)}{2\pi}

    The tenacity factor :math:`K_{sf}` ranges from 1 for sensitive plants to
    1.5 for moderately sensitive plants to 3 for insensitive
    (tenacious plants).

    Parameters
    ----------
    se_root : float
        effective saturation root zone moisture
        `S_{e,root}`
        [-]
    tenacity : float
        tenacity factor
        `K_{sf}`
        [-]

    Returns
    -------
    stress_moist : float
        stress factor for root zone moisture
        :math:`S_{m}`
        [-]

    Examples
    --------
    >>> import ETLook.stress as stress
    >>> stress.stress_moisture(0.5)
    0.75
    >>> stress.stress_moisture(0.5, tenacity = 1)
    0.5
    >>> stress.stress_moisture(0.5, tenacity = 3)
    1.0
    """
    stress = tenacity*se_root - (np.sin(2*np.pi*se_root))/(2*np.pi)
    stress = np.clip(stress, 0, 1)

    return stress


def stress_temperature(t_air_24, t_opt=25.0, t_min=0.0, t_max=50.0):
    r"""
    Computes the stress for plants when it is too cold or hot

    .. math ::
        f=\frac{T_{max}-T_{opt}}{T_{opt}-T_{min}}

    .. math ::
        s_{T}=\frac{\left(T_{a}-T_{min}\right)\left(T_{max}-T_{a}\right)^{f}}
            {\left(T_{opt}-T_{min}\right)\left(T_{max}-T_{opt}\right)^{f}}

    Parameters
    ----------
    t_air_24 : float
        daily air temperature
        :math:`T_{a}`
        [C]
    t_opt : float
        optimum air temperature for plant growth
        :math:`T_{opt}`
        [C]
    t_min : float
        minimum air temperature for plant growth
        :math:`T_{min}`
        [C]
    t_max : float
        maximum air temperature for plant growth
        :math:`T_{max}`
        [C]

    Returns
    -------
    stress_temp : float
        stress factor for air temperature
        :math:`S_{T}`
        [-]

    Examples
    --------
    >>> import ETLook.stress as stress
    >>> stress.stress_temperature(15)
    0.83999999999999997
    >>> stress.stress_temperature(15, t_opt =20)
    0.9451080185178129
    >>> stress.stress_temperature(15, t_opt =20, t_min=10)
    0.79398148148148151
    >>> stress.stress_temperature(15, t_opt =20, t_min=10, t_max=30)
    0.75

    """
    # f = float((t_max - t_opt))/float((t_opt - t_min))
    f = (t_max - t_opt)/(t_opt - t_min)
    x = (t_air_24 - t_min) * (t_max - t_air_24)**f
    y = (t_opt - t_min) * (t_max - t_opt)**f

    stress = x/y
    stress = np.clip(stress, 0, 1)

    return stress


def stress_vpd(vpd_24, vpd_slope=-0.3):
    r"""
    Computes the stress for plants if the vpd increases too much. With
    lower slopes the stress increases faster. The slope of the curve
    is between -0.3 and -0.7

    .. math ::
        S_{v}=m\ln(0.1\Delta_{e}+\frac{1}{2})+1

    Parameters
    ----------
    vpd_24 : float
        daily vapour pressure deficit
        :math:`\Delta_{e}`
        [mbar]
    vpd_slope : float
        vapour pressure stress curve slope
        :math:`m`
        [mbar-1]

    Returns
    -------
    stress_vpd : float
        stress factor for vapour pressure deficit
        :math:`S_{v}`
        [-]

    Examples
    --------
    >>> import ETLook.stress as stress
    >>> stress.stress_vpd(15)
    0.79205584583201638
    >>> stress.stress_vpd(15, vpd_slope=-0.7)
    0.51479697360803833
    >>> stress.stress_vpd(15, vpd_slope=-0.3)
    0.79205584583201638

    """
    stress = vpd_slope * np.log(vpd_24/10. + 0.5) + 1
    stress = np.clip(stress, 0, 1)

    return stress
