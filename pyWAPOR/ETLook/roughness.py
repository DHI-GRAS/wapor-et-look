"""
    The roughness module contains all functions related to surface roughness

"""
import numpy as np

def orographic_roughness(slope, dem_resolution=1000.0):
    r"""
    Computes the orographic roughness, which is a function of the slope of the
    terrain. Steeper slopes will increase the orographic roughness

    .. math ::
        z_{oro}=0.002\left(\frac{\left(x_{res}\frac{\Delta}{100}\right)}
        {x_{res}}^{2}\right)

    Parameters
    ----------
    slope : float
        slope
        :math:`\Delta`
        [rad]
    dem_resolution : float
        resolution of dem
        :math:`x_{res}`
        [m]

    Returns
    -------
    z_oro : float
        orographic roughness
        :math:`z_{oro}`
        [m]

    """

    slope_per = np.tan(slope)
    return 0.002*(((dem_resolution*slope_per)**2)/dem_resolution)

def roughness_length(lai, z_oro, z_obst, z_obst_max, land_mask=1):
    r"""
    Computes the surface roughness length. The roughness length is related to
    the roughness characteristics. For the logarithmic wind-profile the
    surface roughness length is the height at which the wind speed is zero.
    The roughness length is calculated differently for different types of land use

    Land use is classified as follows:

    0. no data
    1. land
    2. water
    3. urban

    .. math ::

        z_{0,m}=\begin{cases}
        \begin{array}{cc}
        0 & l=0\\
        z_{0,m} & l=1\\
        0.0001 & l=2\\
        \frac{1}{7}z_{obst,max}+z_{oro} & l=3
        \end{array}\end{cases}


    Parameters
    ----------
    lai : float
        leaf area index
        :math:`I_{lai}`
        [-]
    z_oro : float
        orographic roughness
        :math:`z_{oro}`
        [m]
    z_obst : float
        obstacle height
        :math:`z_{obst}`
        [m]
    z_obst_max : float
        maximum obstacle height
        :math:`z_{obst,max}`
        [m]
    land_mask : int
        land use classification
        :math:`l`
        [-]

    Returns
    -------
    z0m : float
        roughness length
        :math:`z_{0,m}`
        [m]

    Examples
    --------
    >>> import ETLook.roughness as roughness
    >>> roughness.roughness_length(0.4)
    0.34179999999999999

    """

    def veg_roughness(zo, d, zo_max, l, oro):
        veg = 0.193
        z_dif = zo-d

        term1 = np.minimum(0.41**2/((np.log(z_dif/(0.002*zo_max))+veg)**2), 1.0) + 0.35*l/2
        term2 = np.exp(0.41/np.minimum(np.sqrt(term1), 0.3)-veg)

        return z_dif/term2+oro

    # roughness length specific displacement height
    disp = z_obst * (1 - (1 - np.exp(-np.sqrt(12 * lai))) / (np.sqrt(12 * lai)))
    disp[lai == 0] = 0

    z0m = (1. / 7.) * z_obst_max + z_oro
    z0m[land_mask == 0] = 0
    z0m = np.where(land_mask == 1, veg_roughness(z_obst, disp, z_obst_max, lai, z_oro), z0m)
    z0m[land_mask == 2] = 0.0001

    return z0m


def obstacle_height(ndvi, z_obst_max, ndvi_obs_min=0.25,
                    ndvi_obs_max=0.75, obs_fr=0.25):
    r"""
    Computes the obstacle height. The ndvi is used to limit the obstacle
    height.

    .. math ::
        z_{obst}	=	\begin{cases}
        \begin{array}{cc}
        f_{obs}z_{obst,max} & I_{ndvi}\leq I_{ndvi,obs,min}\\
        z_{obst,max}\left(f_{obs}+\left(1-f_{obs}\right)\left
        (\frac{I_{ndvi}-I_{ndvi,obs,min}}
        {I_{ndvi,obs,max}-I_{ndvi,obs,min}}\right)\right) &
        I_{ndvi}>I_{ndvi,obs,min}\&I_{ndvi}<I_{ndvi,obs,max}\\
        z_{obst,max} & I_{ndvi}\geq I_{ndvi,obs,max}
        \end{array}\end{cases}

    Parameters
    ----------
    ndvi : float
        normalized difference vegetation index
        :math:`I_{ndvi}`
        [-]
    ndvi_obs_min : float
        normalized difference vegetation index @ min obstacle height
        :math:`I_{ndvi,obs,min}`
        [-]
    ndvi_obs_max : float
        normalized difference vegetation index @ max obstacle height
        :math:`I_{ndvi,obs,max}`
        [-]
    obs_fr : float
        ratio of minimum and maximum obstacle height
        :math:`f_{obs}`
        [-]
    z_obst_max : float
        maximum obstacle height
        :math`z_{obst,max}`
        [m]

    Returns
    -------
    z_obst : float
        obstacle height
        :math:`z_{obst}`
        [m]

    Examples
    --------
    >>> import ETLook.roughness as roughness
    >>> roughness.obstacle_height(0.4, 2.0)
    0.95

    """

    cond = [ndvi <= ndvi_obs_min, (ndvi > ndvi_obs_min) &
                 (ndvi < ndvi_obs_max), ndvi >= ndvi_obs_max]

    def frac_func(n):
        return (obs_fr + (1-obs_fr)*(n-ndvi_obs_min)/
                (ndvi_obs_max-ndvi_obs_min))

    obs_height = z_obst_max.copy()
    obs_height = np.where(ndvi <= ndvi_obs_min, obs_fr*z_obst_max,  obs_height)
    obs_height = np.where(np.logical_and(ndvi > ndvi_obs_min, ndvi < ndvi_obs_max), frac_func(ndvi)*z_obst_max, obs_height)
    obs_height = np.where(ndvi >= ndvi_obs_max, z_obst_max, obs_height)

    return obs_height

def displacement_height(lai, z_obst, land_mask=1, c1=1):
    r"""
    Computes the displacement height. The lai is used to limit the displacement
    height. It is defined differently for different types of landuse.

    Land use is classified as follows:

    0. no data
    1. land
    2. water
    3. urban

    .. math ::

        z_{disp}=\begin{cases}
        \begin{array}{cc}
        0 & l=0\\
        z_{obst}\left(1-\frac{1-\exp\left(-\sqrt{c_{1}I_{lai}}\right)}
        {\sqrt{c_{1}I_{lai}}}\right) & l=1\\
        0 & l=2\\
        \frac{2}{3}z_{obst} & l=3
        \end{array}\end{cases}

    Parameters
    ----------
    lai : float
        leaf area index
        :math:`I_{lai}`
        [-]
    z_obst : float
        obstacle height
        :math:`z_{obst}`
        [m]
    land_mask : int
        land use classification
        :math:`l`
        [-]
    c1 : float
        exponential growth rate displacement height function
        :math:`c_1`
        [-]


    Returns
    -------
    disp : float
        displacement height
        :math:`{disp}`
        [m]

    Examples
    --------
    >>> import ETLook.roughness as roughness
    >>> roughness.displacement_height(0.4, 2.0)
    0.51779495

    """

    def disp_func(l):
        disp = z_obst * (1-(1-np.exp(-np.sqrt(c1*l)))/np.sqrt(c1*l))
        disp[l <= 0] = (2./3.)*z_obst[l <= 0]
        return disp

    disp = np.zeros_like(land_mask)
    disp = np.where(land_mask == 1, disp_func(lai), disp)
    disp = np.where(land_mask == 2, 0, disp)
    disp = np.where(land_mask == 3, (2./3.)*z_obst, disp)

    return disp
