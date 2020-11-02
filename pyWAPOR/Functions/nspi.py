# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 11:09:19 2020

@author: rmgu
"""

import numpy as np
from numba import njit, jit


# Based on https://www.sciencedirect.com/science/article/pii/S0034425710003482
# A simple and effective method for filling gaps in Landsat ETM+ SLC-off images

def target_pixel_value(input_image_window, target_image_window, pixel, similarity_threshold,
                       geographic_distance, required_pixels):
    target_pixel_value = np.empty(input_image_window.shape[2]) + np.NaN

    if (np.sum(np.isfinite(input_image_window)) < required_pixels or
        np.sum(np.isfinite(target_image_window)) < required_pixels):
        return target_pixel_value

    target_pixel = input_image_window[pixel[0], pixel[1], :].copy()
    if np.isnan(target_pixel[0]):
        return target_pixel_value
    input_image_window = input_image_window.copy()
    input_image_window[pixel[0], pixel[1], :] = np.NaN

    # Eq 1
    rmsd = np.sqrt(np.nanmean((input_image_window - target_pixel)**2, axis=2))

    similar_pixels = np.where(np.logical_and.reduce((rmsd < similarity_threshold,
                                                     rmsd > 0,
                                                     np.isfinite(target_image_window[:, :, 0]))))
    target_similar_pixels = target_image_window[similar_pixels]
    input_similar_pixels = input_image_window[similar_pixels]

    pixel_num = input_similar_pixels.shape[0]
    if pixel_num < required_pixels:
        return target_pixel_value

    distance_weights = geographic_distance[similar_pixels]
    spectral_weights = rmsd[similar_pixels]

    return _target_pixel_value(target_pixel, distance_weights, spectral_weights,
                               target_similar_pixels, input_similar_pixels, pixel_num)


@jit
def _target_pixel_value(target_pixel, distance_weights, spectral_weights, target_similar_pixels,
                        input_similar_pixels, pixel_num):

    # Eq 5 and 6
    cd = distance_weights * spectral_weights
    weights = (1 / cd) / np.sum(1 / cd)
    weights = np.expand_dims(weights, axis=1)

    # Eq 7
    spectral_similarity_value = nb_nansum_axis_0(weights * target_similar_pixels)
    # Eq 8
    temporal_similarity_value = target_pixel + nb_nansum_axis_0(weights *
                                                                (target_similar_pixels -
                                                                 input_similar_pixels))

    # Eq 9
    landscape_homogeneity = np.sum(spectral_weights) / pixel_num
    # Eq 10
    temporal_homogeneity = np.sum(np.sqrt(nb_mean_axis_0((target_similar_pixels -
                                                          input_similar_pixels)**2))) / pixel_num

    # Eq 11
    t_1 = (1/landscape_homogeneity) / (1/landscape_homogeneity + 1/temporal_homogeneity)
    t_2 = (1/temporal_homogeneity) / (1/landscape_homogeneity + 1/temporal_homogeneity)

    # Eq 12
    target_pixel_value = t_1 * spectral_similarity_value + t_2 * temporal_similarity_value

    return target_pixel_value


def histogram_matching(input_image_window, target_image_window, pixel):

    target_pixel = input_image_window[pixel[0], pixel[1], :]
    if np.isnan(target_pixel[0]):
        return np.empty(input_image_window.shape[2]) + np.NaN

    common_pixels = np.where(np.logical_and(np.isfinite(input_image_window[:, :, 0]),
                                            np.isfinite(target_image_window[:, :, 0])))
    target_common_pixels = target_image_window[common_pixels]
    input_common_pixels = input_image_window[common_pixels]
    if target_common_pixels.size == 0:
        return np.empty(input_image_window.shape[2]) + np.NaN

    # Eq 13
    gain = nb_std_axis_0(target_common_pixels) / nb_std_axis_0(input_common_pixels)
    bias = nb_mean_axis_0(target_common_pixels) - nb_mean_axis_0(input_common_pixels) * gain

    # Eq 14
    target_pixel_value = gain * target_pixel + bias

    return target_pixel_value


# Eq 4
def calc_geographic_distance(window_size):
    ones = np.ones((window_size, window_size))
    all_indices = np.argwhere(ones == 1)
    distance_weights = []
    ci = int(window_size/2)
    for i in all_indices:
        distance_weights.append(np.sqrt((i[0] - ci)**2 + (i[1] - ci)**2))
    geographic_distance = np.array(distance_weights).reshape(ones.shape)
    return geographic_distance


def fill_missing_pixels(input_image, target_image, missing_pixels_mask, similarity_threshold, required_pixels,
                        window_size, geographic_distance, is_max_window):

    filled_image = target_image.copy()
    f = int(window_size/2)
    width, height = input_image.shape[0:2]
    # missing_pixels = np.argwhere(np.isnan(target_image[:, :, 0]))
    missing_pixels = np.argwhere(missing_pixels_mask)

    for i in range(len(missing_pixels)):
        missing_pixel = missing_pixels[i]
        window = [[max(missing_pixel[0]-f, 0), min(missing_pixel[0]+f+1, width)],
                  [max(missing_pixel[1]-f, 0), min(missing_pixel[1]+f+1, height)]]
        window_pixel = [missing_pixel[0] - window[0][0], missing_pixel[1] - window[1][0]]
        filled_image[missing_pixel[0], missing_pixel[1], :] = target_pixel_value(
            input_image[window[0][0]:window[0][1], window[1][0]:window[1][1], :],
            target_image[window[0][0]:window[0][1], window[1][0]:window[1][1], :],
            window_pixel,
            similarity_threshold,
            geographic_distance,
            required_pixels)

        if np.isnan(filled_image[missing_pixel[0], missing_pixel[1], 0]) and is_max_window:
            filled_image[missing_pixel[0], missing_pixel[1], :] = histogram_matching(
                input_image[window[0][0]:window[0][1], window[1][0]:window[1][1], :],
                target_image[window[0][0]:window[0][1], window[1][0]:window[1][1], :],
                window_pixel)
    return filled_image


# Section 2.1
def nspi(input_image, target_image, missing_pixels_mask, num_classes, required_pixels, max_window_size=31):

    # Eq 2
    std_dev_image = np.nanstd(input_image, axis=(0, 1))
    similarity_threshold = (np.sum(std_dev_image) * 2 / num_classes) / np.size(std_dev_image)

    # Eq 3
    initial_window_size = int((required_pixels**0.5 + 1) / 2) * 2 + 1

    is_max_window = False
    for window_size in range(initial_window_size, max_window_size+2, 2):
        print("Window size: %i" % window_size)

        if np.all(np.isfinite(target_image)):
            break

        # Eq 4
        geographic_distance = calc_geographic_distance(window_size)

        if window_size == max_window_size:
            is_max_window = True

        target_image = fill_missing_pixels(input_image, target_image, missing_pixels_mask, similarity_threshold,
                                           required_pixels, window_size, geographic_distance,
                                           is_max_window)

    return target_image


###########################
# Numba helper functions, taken from
# https://github.com/numba/numba/issues/1269#issuecomment-702665837

@jit
def apply_along_axis_0(func1d, arr):
    """Like calling func1d(arr, axis=0)"""
    if arr.size == 0:
        raise RuntimeError("Must have arr.size > 0")
    ndim = arr.ndim
    if ndim == 0:
        raise RuntimeError("Must have ndim > 0")
    elif 1 == ndim:
        return func1d(arr)
    else:
        result_shape = arr.shape[1:]
        out = np.empty(result_shape, arr.dtype)
        _apply_along_axis_0(func1d, arr, out)
        return out


@jit
def _apply_along_axis_0(func1d, arr, out):
    """Like calling func1d(arr, axis=0, out=out). Require arr to be 2d or bigger."""
    ndim = arr.ndim
    if ndim < 2:
        raise RuntimeError("_apply_along_axis_0 requires 2d array or bigger")
    elif ndim == 2:  # 2-dimensional case
        for i in range(len(out)):
            out[i] = func1d(arr[:, i])
    else:  # higher dimensional case
        for i, out_slice in enumerate(out):
            _apply_along_axis_0(func1d, arr[:, i], out_slice)


@jit
def nb_nanmean_axis_0(arr):
    return apply_along_axis_0(np.nanmean, arr)


@jit
def nb_nansum_axis_0(arr):
    return apply_along_axis_0(np.nansum, arr)


@jit
def nb_mean_axis_0(arr):
    return apply_along_axis_0(np.mean, arr)


@jit
def nb_std_axis_0(arr):
    return apply_along_axis_0(np.std, arr)
