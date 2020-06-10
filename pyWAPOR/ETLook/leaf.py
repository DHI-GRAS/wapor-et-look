"""
    The leaf module contains all functions related to estimating vegetation
    cover. These functions only work on an instantaneous basis.

"""
import math


def vegetation_cover(ndvi, nd_min=0.125, nd_max=0.8, vc_pow=0.7):
    r"""
    Computes vegetation cover based on NDVI

    .. math ::
        c_{veg} =
        \begin{cases}
        \begin{array}{cc}
        0 & I_{NDVI}\leq I_{NDVI,min}\\
        1-\left(\frac{I_{NDVI,max}-I_{NDVI}}{I_{NDVI,max}-I_{NDVI,min}}\right)^
        {a} & I_{NDVI,min}<I_{NDVI}<I_{NDVI,max}\\
        1 & I_{NDVI}\geq I_{NDVI,max}
        \end{array}\end{cases}

    Parameters
    ----------
    ndvi : float
        Normalized Difference Vegetation Index
        :math:`I_{NDVI}`
        [-]
    nd_min : float
        NDVI value where vegetation cover is 0
        :math:`I_{NDVI,min}`
        [-]
    nd_max : float
        NDVI value where vegetation cover is 1
        :math:`I_{NDVI,max}`
        [-]
    vc_pow : float
        Exponential power used in vegetation cover function
        :math:`a`
        [-]

    Returns
    -------
    vc : float
        vegetation cover
        :math:`c_{veg}`
        [-]

    Examples
    --------
    >>> from ETLook import leaf
    >>> leaf.vegetation_cover(0.1, nd_min=0.2)
    0
    >>> leaf.vegetation_cover(0.5)
    0.4331446663885373
    >>> leaf.vegetation_cover(0.85)
    1

    .. plot:: pyplots/leaf/plot_vegetation_cover.py

    """
    if ndvi <= nd_min:
        res = 0
    if (ndvi > nd_min) & (ndvi < nd_max):
        res = 1 - ((nd_max - ndvi) / (nd_max - nd_min)) ** vc_pow
    if ndvi >= nd_max:
        res = 1
    return res


def leaf_area_index(vc, vc_min=0.0, vc_max=vegetation_cover(0.795), lai_pow=-0.45):
    r"""
    Computes leaf area index based on vegetation cover. It is based on the
    Kustas formulation of LAI vs NDVI.

    .. math ::
        I_{lai} =
        \begin{cases}
        \begin{array}{cc}
        0 & c_{veg}\leq c_{veg,min}\\
        \frac{\ln\left(-\left(c_{veg}-1\right)\right)}{b} &
        c_{veg,min}<c_{veg}\leq c_{veg,max}\\
        \frac{\ln\left(-\left(c_{veg,max}-1\right)\right)}{b} &
        c_{veg}>c_{veg,max}
        \end{array}\end{cases}

    Parameters
    ----------
    vc : float
        vegetation cover
        :math:`c_{veg}`
        [-]
    vc_min : float
        vegetation cover where LAI is 0
        :math:`c_{veg,min}`
        [-]
    vc_max : float
        vegetation cover at maximum LAI
        :math:`c_{veg,max}`
        [-]
    lai_pow : float
        exponential factor used in LAI function
        :math:`b`
        [-]

    Returns
    -------
    lai : float
        leaf area index
        :math:`I_{lai}`
        [-]

    Examples
    --------
    >>> from ETLook import leaf
    >>> leaf.leaf_area_index(0.0)
    0
    >>> leaf.leaf_area_index(0.5)
    1.5403270679109895
    >>> leaf.leaf_area_index(1.0)
    7.6304274331264414

    .. plot:: pyplots/leaf/plot_leaf_area_index.py

    """
    if vc <= vc_min:
        res = 0
    if (vc > vc_min) & (vc < vc_max):
        res = math.log(-(vc - 1)) / lai_pow
    if vc >= vc_max:
        res = math.log(-(vc_max - 1)) / lai_pow
    return res


def effective_leaf_area_index(lai):
    r"""
    Computes effective leaf area index, this describes the leaf area which
    actively participates in transpiration. It is based on the actual leaf
    area index and an extinction function. So with a higher leaf area index the
    effective leaf area index is a smaller percentage of the total leaf area
    index.

    .. math ::
        I_{lai,eff}=\frac{I_{lai}}{0.3I_{lai}+1.2}

    Parameters
    ----------
    lai : float
        Leaf area index
        :math:`I_{lai}`
        [-]

    Returns
    -------
    lai_eff : float
        effective leaf area index
        :math:`I_{lai,eff}`
        [-]

    Examples
    --------
    >>> from ETLook import leaf
    >>> leaf.effective_leaf_area_index(3.0)
    1.4285714285714288
    >>> leaf.effective_leaf_area_index(5.0)
    1.8518518518518516

    .. plot:: pyplots/leaf/plot_effective_leaf_area_index.py

    """
    lai_eff = lai / ((0.3 * lai) + 1.2)
    return lai_eff
