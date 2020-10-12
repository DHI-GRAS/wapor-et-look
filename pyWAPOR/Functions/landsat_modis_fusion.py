# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 11:38:16 2020

@author: rmgu
"""

import numpy as np

# Based on https://www.spiedigitallibrary.org/journals/journal-of-applied-remote-sensing/volume-9/issue-1/096095/Spatiotemporal-image-fusion-model-for-enhancing-the-temporal-resolution-of/10.1117/1.JRS.9.096095.full?SSO=1
# Spatiotemporal image-fusion model for enhancing the temporal resolution of Landsat-8 surface
# reflectance images using MODIS images

# To be run on 10-day MODIS NDVI composites nad Landsat images
# TO DO
# 1. Read images
# 2. Resample MODIS to Landsat projection, extent and resolution
# 3. Run predict_image
# 4. Save the output


# Section 3
def predict_image(hr_image_t1, lr_image_t1, lr_image_t2):

    rate_of_change = lr_image_t2 / lr_image_t1
    change = lr_image_t2 - lr_image_t1

    negligible_change = np.abs(change) <= 0.15 * rate_of_change
    negative_change = np.logical_and(change < 0, np.abs(change) > 0.15 * rate_of_change)
    positive_change = np.logical_and(change > 0, np.abs(change) > 0.15 * rate_of_change)

    negligible_change_func = np.poly1d(np.polyfit(lr_image_t1[negligible_change],
                                                  lr_image_t2[negligible_change],
                                                  1))
    negative_change_func = np.poly1d(np.polyfit(lr_image_t1[negative_change],
                                                lr_image_t2[negative_change],
                                                1))
    positive_change_func = np.poly1d(np.polyfit(lr_image_t1[positive_change],
                                                lr_image_t2[positive_change],
                                                1))

    predicted_image = np.empty(hr_image_t1.shape)
    predicted_image[negligible_change] = negligible_change_func(hr_image_t1[negligible_change])
    predicted_image[negative_change] = negative_change_func(hr_image_t1[negative_change])
    predicted_image[positive_change] = positive_change_func(hr_image_t1[positive_change])

    return predicted_image
