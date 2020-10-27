

import os
import shutil
import tarfile

import numpy as np
import rasterio as rio

from tqdm import tqdm
from pathlib import Path
from datetime import datetime
from rasterio.vrt import WarpedVRT
from pyWAPOR.Functions.nspi import nspi


def PreprocessLandsat(landsat_dir):

    if isinstance(landsat_dir, str):
        landsat_dir = Path(landsat_dir)

    # unpack the *.tar.gz Landsat files
    L7_files = list((landsat_dir / Path('L7')).glob('*.tar.gz'))
    L8_files = list((landsat_dir / Path('L8')).glob('*.tar.gz'))

    print('Unpacking *.tar.gz files...')
    for file in tqdm(L7_files + L8_files):
        _unpack_and_save(file, delete_input=False)

    # merge the individual landsat bands into multiband files
    L7_dirs = [directory for directory in list((landsat_dir / Path('L7')).glob('*')) if os.path.isdir(directory)]
    L8_dirs = [directory for directory in list((landsat_dir / Path('L8')).glob('*')) if os.path.isdir(directory)]

    filename_list = []
    bandnames_list = []

    print('Merging Landsat bands...')
    for directory in tqdm(L7_dirs+L8_dirs):
        filename, band_names = _merge_and_save_landsat(directory, delete_input=False)
        filename_list.append(filename)
        bandnames_list.append(band_names)

    # apply nspi gap-filling on the landsat-7 data (slow!)
    _apply_nspi(landsat_dir, filename_list, band_names)

    # calculate NDVI/Albedo

    # For NDVI:
        # Fill gaps using Weiss et. al. 2014
        # Smoothe using the Savitzky-Golay filter
        # Calculate dekadal NDVI as mean of the dekade

    # For Albedo
        # Smoothe using the Savitzky-Golay filter
        # Calculate dekadal NDVI as mean of the dekade


def _apply_nspi(landsat_dir, filename_list, bandnames_list):
    L7_idx = [i for i, file in enumerate(filename_list) if str(file.split('_')[0]) == 'LE07']
    L8_idx = [i for i, file in enumerate(filename_list) if str(file.split('_')[0]) == 'LC08']

    L7_tifs = [filename_list[idx] + '.tif' for idx in L7_idx]
    L8_tifs = [filename_list[idx] + '.tif' for idx in L8_idx]

    L7_bandnames = [bandnames_list[idx] for idx in L7_idx]
    L8_bandnames = [bandnames_list[idx] for idx in L8_idx]

    L7_dates = [datetime.strptime(str(f).split('_')[3], '%Y%m%d') for f in L7_tifs]
    L8_dates = [datetime.strptime(str(f).split('_')[3], '%Y%m%d') for f in L8_tifs]

    print('Fill gaps (NSPI)...')
    for L7_idx, L7_file in enumerate(tqdm(L7_tifs)):

        L7_date = L7_dates[L7_idx]
        date_diff = [np.abs(L7_date - L8_date) for L8_date in L8_dates]
        L8_idx = np.argmin(date_diff)
        L8_file = L8_tifs[L8_idx]

        # open target image
        with rio.open(str(landsat_dir / Path('L7') / Path(L7_file))) as master_src:
            meta = master_src.profile
            L7_image = master_src.read()

        master_dict = {'transform': master_src.transform, 'height': master_src.height, 'width': master_src.width,
                       'crs': master_src.crs}

        # open input image
        with rio.open(str(landsat_dir / Path('L8') / Path(L8_file))) as slave_src:
            with WarpedVRT(slave_src, **master_dict) as vrt:
                L8_image = vrt.read()

        bands = ['sr_band1', 'sr_band2', 'sr_band3', 'sr_band4', 'sr_band5', 'sr_band7', ]

        L7_band_idx = [L7_bandnames[L7_idx].index(band) for band in bands]
        L8_band_idx = [L8_bandnames[L8_idx].index(band) for band in bands]

        # reshape images for correct NSPI format
        L7_image_reshaped = np.transpose(L7_image[L7_band_idx, ...], axes=[1, 2, 0])
        L8_image_reshaped = np.transpose(L8_image[L8_band_idx, ...], axes=[1, 2, 0])

        # add nans instead of nodata
        L7_image_reshaped = np.where(L7_image_reshaped == -9999, np.nan, L7_image_reshaped)
        L8_image_reshaped = np.where(L8_image_reshaped == -9999, np.nan, L8_image_reshaped)

        # calculate the missing pixels mask
        missing_pixels_mask = np.isnan(L7_image_reshaped[..., 0]) & ~np.isnan(L8_image_reshaped[..., 0])

        # get the bitmask
        L7_pixel_qa = L7_image[L7_bandnames[L7_idx].index('pixel_qa'), ...]
        L7_pixel_qa = np.where(L7_pixel_qa == -9999, 0, L7_pixel_qa)

        # calculate cloud mask
        L7_pixel_cloudmask = _landsat_cloudmask(L7_pixel_qa)

        # mask clouds
        for n in range(0, len(bands)):
            L7_image_reshaped[..., n] = np.where(L7_pixel_cloudmask, np.nan, L7_image_reshaped[..., n])

        out_image_reshaped = nspi(L8_image_reshaped, L7_image_reshaped, missing_pixels_mask,
                                  num_classes=5, required_pixels=20, max_window_size=15)

        out_image = L7_image.copy()

        out_image[L7_band_idx, ...] = np.transpose(out_image_reshaped, axes=[2, 0, 1])

        out_image = np.where(np.isnan(out_image), -9999, out_image).astype(np.int16)

        output_filename = str(Path(L7_file).parent / Path(L7_file).stem) + '_gap-filled.tif'

        with rio.open(output_filename, 'w', **meta) as dst:
            dst.write(out_image)


# merges individual landsat bands and saves as single tif
def _merge_and_save_landsat(directory, delete_input=False):
    master_file = list(directory.glob('*pixel_qa.tif'))[0]
    slave_files = [f for f in list(directory.glob('*.tif')) if 'pixel_qa' not in str(f)]

    band_names = ['_'.join(master_file.stem.split('_')[-2:])]

    # open master and add to array
    with rio.open(str(master_file)) as master_src:
        meta = master_src.profile
        pixel_qa = master_src.read().squeeze()

    master_dict = {'transform': master_src.transform, 'height': master_src.height, 'width': master_src.width,
                   'crs': master_src.crs}

    nodata_mask = np.where(pixel_qa == 1, True, False)

    data = np.zeros((len(slave_files) + 1, pixel_qa.shape[0], pixel_qa.shape[1])).astype(np.int16)

    data[0, ...] = pixel_qa

    # for each slave-file, open and append to array
    for i, file in enumerate(slave_files):
        with rio.open(file) as slave_src:
            with WarpedVRT(slave_src, **master_dict) as vrt:
                data[i + 1, ...] = vrt.read().squeeze()

        band_names.append('_'.join(file.stem.split('_')[-2:]))

    data[:, nodata_mask] = -9999

    # save file as tif in root
    meta.update({'count': data.shape[0],
                 'dtype': str(data.dtype),
                 'blockxsize': 256,
                 'blockysize': 256,
                 'tiled': True,
                 'compress': 'lzw',
                 'interleave': 'pixel',
                 'nodata': -9999})

    output_filename = master_file.parents[1] / Path('_'.join(master_file.stem.split('_')[0:-2]) + '.tif')

    with rio.open(output_filename, 'w', **meta) as dst:
        dst.write(data)

    if delete_input:
        shutil.rmtree(directory)

    return output_filename.stem, band_names


# unpack and saves *.tar.gz files
def _unpack_and_save(file, delete_input=False):
    path = Path(file).parent / Path(Path(file).stem).stem

    if ~os.path.isdir(path):
        os.mkdir(path)

    tar = tarfile.open(file, "r:gz")
    tar.extractall(path=path)
    tar.close()

    if delete_input:
        os.remove(file)


# returns cloud mask given quality band as input
def _landsat_cloudmask(quality_band):
    # if clouds (bit 5) and low/medium/high probability (bit 6 and 7) then clouds
    clouds = ((quality_band & (1 << 5)) > 1) & ((quality_band & ((1 << 6) | (1 << 7))) > 1)
    # if shadows (pixel 3)
    shadows = (quality_band & (1 << 3)) > 1

    return clouds | shadows
