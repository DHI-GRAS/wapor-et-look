# -*- coding: utf-8 -*-
"""
Swets NDVI filtering
author: Laust Færch @ DHI GRAS
Created on 2020/08/29

Based on the article:
Swets, D.L, Reed, B.C., Rowland, J.D., Marko, S.E., 1999. A weighted least-squares approach to temporal
NDVI smoothing. In: Proceedings of the 1999 ASPRS Annual Conference, Portland, Oregon, pp. 526-536.

"""

import numpy as np
from numba import jit
from scipy import interpolate
from scipy.ndimage.filters import generic_filter


def _interpolate1d(data):

    # we need at least 2 non-nan elements
    if np.sum(~np.isnan(data)*1) < 2:
        return data

    good = ~np.isnan(data)

    # scipy interpolation
    finterp = interpolate.interp1d(np.flatnonzero(good),
                                   data[good], kind='linear',
                                   fill_value=np.nan,
                                   bounds_error=False)

    yinterp = finterp(np.arange(data.shape[0]))

    return yinterp


# calculate the weight of each sample based on neighbours
def _calc_weights(y):
    # for class local_peak, sloping_points, local_valley
    class_weights = [1.5, 0.5, 0.005] # weights defined in article

    left_shift = (y - np.roll(y, -1)) >= 0
    right_shift = (y - np.roll(y, 1)) >= 0

    peaks = left_shift & right_shift
    valleys = (~left_shift) & (~right_shift)
    slopes = (~peaks) & (~valleys)

    weights = np.zeros_like(y)

    weights[peaks] = class_weights[0]
    weights[slopes] = class_weights[1]
    weights[valleys] = class_weights[2]

    return weights


# calculate the weighted linear regression
@jit(nopython=True)
def _calc_linreg(x, y, w):
    eps = 1e-8

    sw = np.sum(w)
    sy = np.sum(w * y)
    sx = np.sum(w * x)
    sxy = np.sum(w * x * y)
    sx2 = np.sum(w * x ** 2)

    num = (sw * sxy - sx * sy)
    denom = (sw * sx2 - sx ** 2)

    if denom == 0:
        b = 0
    else:
        b = num / denom

    a = (sy - b * sx) / (sw + eps)

    return a, b


@jit(nopython=True)
def _calc_linreg_wrapper_a(xyw):
    n = int(np.round(xyw.shape[0] / 3))
    xyw = xyw.reshape(3, n)

    a, b = _calc_linreg(xyw[0, :], xyw[1, :], xyw[2, :])

    return a


@jit(nopython=True)
def _calc_linreg_wrapper_b(xyw):
    n = int(np.round(xyw.shape[0] / 3))
    xyw = xyw.reshape(3, n)

    a, b = _calc_linreg(xyw[0, :], xyw[1, :], xyw[2, :])

    return b


def _piecewise_linreg(xyw, window_width=3):
    n = int(np.round((window_width - 1) / 2))

    piece_a = generic_filter(xyw.T, _calc_linreg_wrapper_a, size=(3, window_width), mode='nearest')
    piece_b = generic_filter(xyw.T, _calc_linreg_wrapper_b, size=(3, window_width), mode='nearest')

    # pad array
    piece_a = np.pad(piece_a[1, :], n, 'edge')
    piece_b = np.pad(piece_b[1, :], n, 'edge')

    smooth_a = np.convolve(piece_a, np.ones(window_width) / window_width, mode='valid')
    smooth_b = np.convolve(piece_b, np.ones(window_width) / window_width, mode='valid')

    y_est = smooth_b * xyw[:, 0] + smooth_a

    return y_est


# Apply the swets filter on 1d array
def _apply_swets1d(y):
    window_width = 3  # window width

    # dont smooth if all nan
    if np.all(np.isnan(y)):
        return y

    y_smoothed = np.zeros_like(y)
    y_smoothed[:] = np.nan

    x = np.flatnonzero(~np.isnan(y))
    w = _calc_weights(y[x])
    xyw = np.stack((x, y[x], w), axis=1)
    y_smoothed[x] = _piecewise_linreg(xyw, window_width)

    return y_smoothed

# slow: 3024x3024x7 array takes approx 1 hour on my machine
def swets_filter(data, do_interpolate=True, axis=2, invert=False):
    """
    :param data: np.array((y,x,t))
    NDVI raster timeseries. Each image in timeseries is shape y,x. Images are stacked along t (time dimension)
    nan values are initially ignored by swets. At the last step, nan values are replaced by linear interpolation
    :param do_interpolate: bool
    True if we want to apply interpolation of nan values
    :param axis: int
    Axis for the time-dimension in the array. Filtering/interpolation will be apply along this axis
    :param invert: bool
    Inversion of the data and output. This is useful for albedo where bad clous maksing will force the values up,
    instead of down in NDVI.

    :return: y_smoothed: np.array(y,x,t)
    """

    if invert:
        data = data * -1

    print('running filter...')
    y_smoothed = np.apply_along_axis(_apply_swets1d, axis, data)

    if do_interpolate:
        print('running interpolation...')
        y_smoothed = np.apply_along_axis(_interpolate1d, axis, y_smoothed)

    if invert:
        y_smoothed = y_smoothed * -1

    return y_smoothed.astype(data.dtype)
