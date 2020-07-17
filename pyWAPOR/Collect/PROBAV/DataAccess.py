
"""
Author: Laust Færch
Module: Collect/PROBAV
"""

import glob
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
from datetime import datetime, timedelta

import nest_asyncio
nest_asyncio.apply()


def download_data(download_dir, start_date, end_date, latitude_extent, longitude_extent, username,
                  password):

    download_dir = Path(os.path.join(download_dir, "ProbaV"))
    if not os.path.isdir(download_dir):
        os.mkdir(download_dir)

    start_date = datetime.strptime(start_date, "%Y-%m-%d")
    end_date = datetime.strptime(end_date, "%Y-%m-%d")

    # download for the date interval +/- the time_buffer
    time_buffer = timedelta(days=14)
    delta = (end_date + time_buffer) - (start_date - time_buffer)

    # Loop over all dates
    for i in tqdm(range(delta.days + 1)):
        date = (start_date - time_buffer) + timedelta(days=i)

        if (glob.glob(os.path.join(download_dir, "*%s.tif" % date.strftime("%Y%m%d"))) or
            glob.glob(os.path.join(download_dir, "*%s*.HDF5" % date.strftime("%Y%m%d")))):
            continue

        # retrieve vito URL
        url = vito.build_url(product='Proba-V-NDVI', year=date.year, month=date.month,
                             day=date.day,
                             extent={'xmin': longitude_extent[0], 'xmax': longitude_extent[1],
                                     'ymin': latitude_extent[0], 'ymax': latitude_extent[1]})

        max_retries = 4
        no_of_attempts = 0
        download_success = False

        # TODO: maybe put all this in a function called _download_hdf5 (line 49-70)
        # sometimes vito fails, sometimes it works to retry.
        while no_of_attempts <= max_retries:
            try:
                # TODO: make check to see if *.HDF5 files exists already
                # TODO: we need to figure out why vito sometimes fails
                # download all matching files
                local_files = vito.download_data(url, username=username, password=password,
                                                 download_dir=download_dir, include='*.HDF5',
                                                 download_jobs=4)
                downloaded_files = list(local_files)

                download_success = True

            except:
                no_of_attempts += 1
                continue
            break

        if not download_success:
            warnings.warn(f'vito download failed for date: {date.strftime("%Y-%m-%d")}')
            break

        # convert downloaded HDF5 files to tif files
        for file in downloaded_files:
            # TODO: make it possible to delete hdf5 files
            da = _hdf5_to_dataarray(file, 'LEVEL3/NDVI')
            _dataarray_to_tif(da, str(Path(file).with_suffix('.tif')))

        # loop over all tif files in download folder
        input_files = []
        for file in download_dir.glob('*.tif'):
            date_str = max(re.findall('[0-9]+', str(file)), key=len)
            current_date_str = date.strftime('%Y%m%d')
            # collect all files that matches current date
            if date_str == current_date_str:
                input_files.append(str(file))

        output_file = str(download_dir / f'NDVI_{date.strftime("%Y%m%d")}.tif')

        # merge files and clip to extent, save as tif
        if input_files:
            _merge_and_save_tifs(input_files, output_file, latitude_extent, longitude_extent,
                                 delete_input=True)

    return()


# save xarray DataArrays as tif
def _dataarray_to_tif(da, filename):

    # TODO: we assume crs is EPSG:4326 - change to dynamic crs later
    meta = {'driver': 'GTiff',
            'height': da.shape[0],
            'width': da.shape[1],
            'dtype': str(da.dtype),
            'count': 1,
            'transform': rasterio.Affine.from_gdal(*da.attrs['affine']),
            'crs': rasterio.crs.CRS.from_string('EPSG:4326')}

    with rasterio.open(filename, 'w', **meta) as dst:
        dst.write(da.values, 1)


# open hdf5 file as xarray DataArray
# if this function doesn't work: try to update your xarray and netcdf libraries.
def _hdf5_to_dataarray(filename, group):
    # read the group data
    with xr.open_dataset(filename, group=group) as src:
        da = src.NDVI
    # fetch metadata
    with xr.open_dataset(filename) as src:
        meta = src.copy()

    da.attrs['crs'] = meta.crs.spatial_ref
    da.attrs['affine'] = [float(i) for i in meta.crs.GeoTransform.split(' ')[0:6]]

    return da


# merge tif files and save as one
def _merge_and_save_tifs(input_files, output_file, latitude_extent, longitude_extent,
                         delete_input=True):
    extent_poly = Polygon([[(longitude_extent[0], latitude_extent[0]),
                            (longitude_extent[1], latitude_extent[0]),
                            (longitude_extent[1], latitude_extent[1]),
                            (longitude_extent[0], latitude_extent[1])]])

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
