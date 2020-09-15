
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

# Required for Python 3.6 and 3.7
import nest_asyncio
nest_asyncio.apply()


def download_data(download_dir, start_date, end_date, latitude_extent, longitude_extent, username,
                  password, product):

    # setup
    max_retries = 5
    delete_hdf5 = False
    delete_tif = True

    if not os.path.isdir(download_dir):
        os.mkdir(download_dir)

    start_date = datetime.strptime(start_date, "%Y-%m-%d")
    end_date = datetime.strptime(end_date, "%Y-%m-%d")

    # download for the date interval +/- the time_buffer
    time_buffer = timedelta(days=7)
    delta = (end_date + time_buffer) - (start_date - time_buffer)

    # Loop over all dates
    for i in tqdm(range(delta.days + 1)):
        date = (start_date - time_buffer) + timedelta(days=i)

        # retrieve vito URL
        url = vito.build_url(product=product, year=date.year, month=date.month, day=date.day,
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

            except (RuntimeError, HTTPError):
                no_of_attempts += 1
                continue
            break

        if not download_success:
            warnings.warn(f'vito download failed for date: {date.strftime("%Y-%m-%d")}')
            break

        # convert downloaded HDF5 files to tif files
        for hdf_file in downloaded_files:
            if product == 'Proba-V-S5-TOC-NDVI':
                da = _hdf5_to_dataarray(hdf_file, 'LEVEL3/NDVI', product)
                _dataarray_to_tif(da, str(Path(hdf_file).with_suffix('.tif')))

            elif product == 'Proba-V-S5-TOC':
                band_list = ['BLUE', 'NIR', 'RED', 'SWIR']
                # read all bands and save as individual tifs
                for band in band_list:
                    da = _hdf5_to_dataarray(hdf_file, f'LEVEL3/RADIOMETRY/{band}', product)
                    _dataarray_to_tif(da, str(Path(hdf_file).parent / Path(hdf_file).stem) + f'_{band}.tif')
                da = _hdf5_to_dataarray(hdf_file, 'LEVEL3/QUALITY', 'quality')
                _dataarray_to_tif(da, str(Path(hdf_file).parent / Path(hdf_file).stem) + '_SM.tif')

        # loop over all tif files in download folder
        input_files = []
        for tif_file in download_dir.glob('*.tif'):
            date_str = max(re.findall('[0-9]+', str(tif_file)), key=len)
            # collect all files that matches current date
            if date_str == date.strftime('%Y%m%d'):
                input_files.append(str(tif_file))

        if product == 'Proba-V-S5-TOC-NDVI':
            output_file = str(download_dir / f'NDVI_{date.strftime("%Y-%m-%d")}.tif')
        elif product == 'Proba-V-S5-TOC':
            output_file = str(download_dir / f'ALBEDO_{date.strftime("%Y-%m-%d")}.tif')
            input_files = _process_and_save_albedo(input_files, delete_tif)

        # merge files and clip to extent, save as tif
        if input_files:
            _merge_and_save_tifs(input_files, output_file, latitude_extent, longitude_extent,
                                 delete_input=delete_tif)

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
def _hdf5_to_dataarray(filename, group, product):
    # read the group data
    with xr.open_dataset(filename, group=group, engine='netcdf4') as src:

        if product == 'Proba-V-S5-TOC-NDVI':
            da = src.NDVI
        elif product == 'Proba-V-S5-TOC':
            da = src.TOC
        elif product == 'quality':
            da = src.SM

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

    # check if we have data in aoi
    if np.any(~((mosaic == 0) | np.isnan(mosaic))):

        meta = {
            'driver': 'GTiff',
            'width': mosaic.shape[1],
            'height': mosaic.shape[0],
            'count': 1,
            'dtype': str(mosaic.dtype),
            'crs': crs,
            'transform': mosaic_trans
        }

        with rasterio.open(output_file, 'w', **meta) as dst:
            dst.write(mosaic, 1)

    if delete_input:
        for src in src_files_to_mosaic:
            src.close()

        for file in input_files:
            os.remove(file)


# calculate albedo from band files
# Albedo should be calculated for each downloaded tile and saved as tif
def _process_and_save_albedo(input_files, delete_input=True):
    tiles = set([file.split('_')[4] for file in input_files])
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

        surface_albedo = (0.429 * all_bands[band_names.index('BLUE'), ...]
                          + 0.333 * all_bands[band_names.index('RED'), ...]
                          + 0.133 * all_bands[band_names.index('NIR'), ...]
                          + 0.105 * all_bands[band_names.index('SWIR'), ...])

        # Mask clouds
        quality_band = all_bands[band_names.index('SM'), ...]
        flag_mask = [1, 2, 3, 4]  # {1: shadow, 2: cloud, 3: undefined, 4: ice+snow}
        bit_mask_array = np.bitwise_or.reduce(flag_mask)*np.ones_like(quality_band)
        mask = np.bitwise_and(quality_band.astype(np.int), bit_mask_array.astype(np.int)) > 0
        surface_albedo[mask] = None

        albedo_filename = str(Path(matching_files[0]).parent /
                              Path(str('_').join(Path(matching_files[0]).stem.split('_')[0:-1]))) + '_ALBEDO.tif'

        with rasterio.open(albedo_filename, 'w', **meta) as dst:
            dst.write(surface_albedo, 1)

        updated_input_files.append(albedo_filename)

    if delete_input:
        for file in input_files:
            os.remove(file)

    return updated_input_files
