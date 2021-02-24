
"""
Author: Laust Færch
Module: Collect/PROBAV
"""

import re
import os
import warnings
import rasterio
import rasterio.mask
import rasterio.merge

import numpy as np
import xarray as xr
import vito_download as vito

from tqdm import tqdm
from pathlib import Path
from geojson import Polygon
from requests.exceptions import HTTPError
from datetime import datetime, timedelta
from aiohttp import ClientResponseError, ServerDisconnectedError
from asyncio import TimeoutError
# Required for Python 3.6 and 3.7
import nest_asyncio
nest_asyncio.apply()


def download_data(download_dir, start_date, end_date, latitude_extent, longitude_extent, username,
                  password):

    # setup
    max_retries = 5
    delete_hdf5 = False
    delete_tif = True
    dataset = 'Proba-V-S5-TOC'

    if not os.path.isdir(download_dir):
        os.mkdir(download_dir)

    start_date = datetime.strptime(start_date, "%Y-%m-%d")
    end_date = datetime.strptime(end_date, "%Y-%m-%d")

    # Buffer the download to make sure that a decadal composite can be produced regarding of start
    # and end dates
    start_date = start_date - timedelta(days=14)
    end_date = end_date + timedelta(days=8)
    delta = end_date - start_date

    # Loop over all dates
    for i in tqdm(range(delta.days + 1)):
        date = start_date + timedelta(days=i)

        # retrieve vito URL
        url = vito.build_url(product=dataset, year=date.year, month=date.month, day=date.day,
                             extent={'xmin': longitude_extent[0], 'xmax': longitude_extent[1],
                                     'ymin': latitude_extent[0], 'ymax': latitude_extent[1]})

        no_of_attempts = 0
        download_success = False
        downloaded_files = []

        # sometimes vito.download fails, often it works to retry.
        while no_of_attempts <= max_retries:
            try:
                # download all matching files
                local_files = vito.download_data(url, username=username, password=password,
                                                 download_dir=download_dir, include='*.HDF5',
                                                 download_jobs=4)
                downloaded_files = list(local_files)

                for file in downloaded_files:
                    # If file is corrupted, delete it and retry
                    try:
                        with xr.open_dataset(file, engine='netcdf4') as src:
                            crs = src.crs # all files should have crs
                    except OSError:
                        os.remove(file)
                        raise RuntimeError

                download_success = True

            except (RuntimeError, HTTPError,
                    ClientResponseError, ServerDisconnectedError, TimeoutError):
                no_of_attempts += 1
                continue
            break

        if not download_success:
            warnings.warn(f'vito download failed for date: {date.strftime("%Y-%m-%d")}')
            break

        # convert downloaded HDF5 files to tif files
        for hdf_file in downloaded_files:
            da = _hdf5_to_dataarray(hdf_file, 'LEVEL3/NDVI', 'NDVI')
            _dataarray_to_tif(da, str(Path(hdf_file).parent / Path(hdf_file).stem) + '_NDVI.tif')

            band_list = ['BLUE', 'NIR', 'RED', 'SWIR']
            # read all bands and save as individual tifs
            for band in band_list:
                da = _hdf5_to_dataarray(hdf_file, f'LEVEL3/RADIOMETRY/{band}', 'TOC')
                _dataarray_to_tif(da, str(Path(hdf_file).parent /
                                          Path(hdf_file).stem) + f'_{band}.tif')

            angle_band_list = ['VNIR', 'SWIR']
            # read all angle bands and save as individual tifs
            for band in angle_band_list:
                da = _hdf5_to_dataarray(hdf_file, f'LEVEL3/GEOMETRY/{band}', 'VZA')
                _dataarray_to_tif(da, str(Path(hdf_file).parent /
                                          Path(hdf_file).stem) + f'_{band}-VZA.tif')

            # read and save quality mask
            da = _hdf5_to_dataarray(hdf_file, 'LEVEL3/QUALITY', 'SM')
            _dataarray_to_tif(da, str(Path(hdf_file).parent / Path(hdf_file).stem) + '_SM.tif')

        # loop over all tif files in download folder
        input_files = []
        for tif_file in download_dir.glob('*.tif'):
            date_str = max(re.findall('[0-9]+', str(tif_file)), key=len)
            # collect all files that matches current date
            if date_str == date.strftime('%Y%m%d'):
                input_files.append(str(tif_file))

        for product in ['NDVI', 'ALBEDO']:
            output_file = str(download_dir /
                              ('%s/%s_%s.tif' % (product, product, date.strftime("%Y-%m-%d"))))
            data_files = _preprocess_inputs(input_files, product)

            # merge files and clip to extent, save as tif
            if data_files:
                _merge_and_save_tifs(data_files, output_file, latitude_extent, longitude_extent,
                                     delete_input=delete_tif)

        if delete_tif:
            for file in input_files:
                os.remove(file)

    if delete_hdf5:
        for file in list(download_dir.glob('*.HDF5')):
            os.remove(file)

    return()


# save xarray DataArrays as tif
def _dataarray_to_tif(da, filename):

    meta = {'driver': 'GTiff',
            'height': da.shape[0],
            'width': da.shape[1],
            'dtype': str(da.dtype),
            'count': 1,
            'transform': rasterio.Affine.from_gdal(*da.attrs['affine']),
            'crs': da.crs}

    with rasterio.open(filename, 'w', **meta) as dst:
        dst.write(da.values, 1)


# open hdf5 file as xarray DataArray
# if this function doesn't work: try to update your xarray and netcdf libraries.
def _hdf5_to_dataarray(filename, group, dataset_name):
    # read the group data
    with xr.open_dataset(filename, group=group, engine='netcdf4') as src:
        da = src[dataset_name]

    # fetch metadata
    with xr.open_dataset(filename) as src:
        meta = src.copy()

    da.attrs['crs'] = meta.crs.spatial_ref
    da.attrs['affine'] = [float(i) for i in meta.crs.GeoTransform.split(' ')[0:6]]

    return da


# merge tif files and save as one
def _merge_and_save_tifs(input_files, output_file, latitude_extent, longitude_extent,
                         delete_input=True):
    extent_poly = Polygon([(longitude_extent[0], latitude_extent[0]),
                           (longitude_extent[1], latitude_extent[0]),
                           (longitude_extent[1], latitude_extent[1]),
                           (longitude_extent[0], latitude_extent[1])])

    bbox = rasterio.features.bounds(extent_poly)
    src_files_to_mosaic = []
    crs = None

    for file in list(input_files):
        # use crs from first file as crs for mosaic
        if crs is None:
            with rasterio.open(file) as src:
                meta = src.profile
                crs = meta['crs']

        src_files_to_mosaic.append(rasterio.open(file))

    mosaic, mosaic_trans = rasterio.merge.merge(src_files_to_mosaic, bounds=bbox)
    mosaic = np.squeeze(mosaic)

    # Do not check if we have data in aoi


    meta = {
        'driver': 'GTiff',
        'width': mosaic.shape[1],
        'height': mosaic.shape[0],
        'count': 1,
        'dtype': str(mosaic.dtype),
        'crs': crs,
        'transform': mosaic_trans
    }

    if not os.path.exists(os.path.dirname(output_file)):
        os.mkdir(os.path.dirname(output_file))
    with rasterio.open(output_file, 'w', **meta) as dst:
        dst.write(mosaic, 1)

    if delete_input:
        for src in src_files_to_mosaic:
            src.close()

        for file in input_files:
            os.remove(file)


# calculate albedo from band files
# Albedo should be calculated for each downloaded tile and saved as tif
def _preprocess_inputs(input_files, product):
    tiles = set([os.path.basename(file).split('_')[3] for file in input_files])
    updated_input_files = []

    for tile in list(tiles):
        matching_files = [elem for elem in input_files if tile in elem]

        all_data = []
        band_names = []
        for file in matching_files:
            with rasterio.open(str(file)) as src:
                all_data.append(src.read())
                meta = src.profile
            band_names.append(str(Path(file).stem).split('_')[-1])

        all_bands = np.asarray(all_data).squeeze()

        if product == 'NDVI':
            data = all_bands[band_names.index('NDVI'), ...]
        elif product == 'ALBEDO':
            data = (0.429 * all_bands[band_names.index('BLUE'), ...]
                    + 0.333 * all_bands[band_names.index('RED'), ...]
                    + 0.133 * all_bands[band_names.index('NIR'), ...]
                    + 0.105 * all_bands[band_names.index('SWIR'), ...])

        # Mask clouds
        quality_band = all_bands[band_names.index('SM'), ...]
        flag_mask = [1, 2, 3, 4]  # {1: shadow, 2: cloud, 3: undefined, 4: ice+snow}
        bit_mask_array = np.bitwise_or.reduce(flag_mask)*np.ones_like(quality_band)
        cloudmask = np.bitwise_and(quality_band.astype(np.int), bit_mask_array.astype(np.int)) > 0
        data[cloudmask] = None

        # Mask Viewing Angle
        max_angle = 35

        swir_vza = all_bands[band_names.index('SWIR-VZA'), ...]
        vnir_vza = all_bands[band_names.index('VNIR-VZA'), ...]

        swir_mask = swir_vza >= max_angle
        vnir_mask = vnir_vza >= max_angle

        if product == 'NDVI':
            data[vnir_mask] = None
        elif product == 'ALBEDO':
            data[vnir_mask & swir_mask] = None

        data_filename = str(Path(matching_files[0]).parent /
                            Path(str('_').join(Path(matching_files[0]).stem.split('_')[0:-1]))) + '_temp.tif'

        with rasterio.open(data_filename, 'w', **meta) as dst:
            dst.write(data, 1)

        updated_input_files.append(data_filename)

    return updated_input_files
