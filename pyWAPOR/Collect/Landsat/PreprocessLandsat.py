

import os
import shutil
import tarfile

import numpy as np
import rasterio as rio

from tqdm import tqdm
from pathlib import Path
from datetime import datetime
from datetime import timedelta
from pyWAPOR.Functions.nspi import nspi
from pyWAPOR.Functions.SavGol_Filter import savgol_reconstruct
from pyWAPOR.Pre_ETLook import _get_dekadal_date

def PreprocessLandsat(landsat_dir, output_dir):

    if isinstance(landsat_dir, str):
        landsat_dir = Path(landsat_dir)

    # unpack the *.tar.gz Landsat files
    L7_files = list((landsat_dir / Path('L7')).glob('*.tar'))
    L8_files = list((landsat_dir / Path('L8')).glob('*.tar'))

    print('Unpacking *.tar files...')
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


    # apply nspi gap-filling on the landsat-7 data (slow!)[
    _apply_nspi(landsat_dir, filename_list, bandnames_list)

    # calculate NDVI/ABEDO and save
    _process_and_save(landsat_dir, filename_list, bandnames_list, output_dir, delete_input=False)


def _process_and_save(landsat_dir, filename_list, bandnames_list, output_folder, delete_input=False):

    NDVI = []
    ALBEDO = []

    master_src = None
    dates = [datetime.strptime(str(f).split('_')[3], '%Y%m%d') for f in filename_list]

    # we need to sort according to date before filtering the timeseries
    sorted_dates = [dates[idx] for idx in np.argsort(dates)]
    sorted_filenames = [filename_list[idx] for idx in np.argsort(dates)]
    sorted_bandnames = [bandnames_list[idx] for idx in np.argsort(dates)]

    print('Calculating NDVI/ALBEDO...')
    # TODO: # enable delete inputs
    
    # Apply NDVI and albedo calculations and smoothing in tiles otherwise all the landsat data
    # might not fit in memory
    filename = _image_path(landsat_dir, sorted_filenames[0])[0]
    with rio.open(filename) as src:
        block_windows = src.block_window(1)[1]

    for window in block_windows:
        for i, file in enumerate(tqdm(sorted_filenames)):

            filename, sensor = _image_path(landsat_dir, file)
            bandnames = sorted_bandnames[i]

            # use the first file as master
            if master_src is None:
                # open with rasterio
                with rio.open(filename) as master_src:
                    data = master_src.read(window=window)
                    meta = master_src.profile

                master_dict = {'transform': master_src.transform,
                               'height': master_src.height,
                               'width': master_src.width,
                               'crs': master_src.crs}

            # all files after master are reprojected to match master
            else:
                with rio.open(filename) as src:
                    with rio.vrt.WarpedVRT(src, **master_dict) as vrt:
                        data = vrt.read(window=window)

            # calculate NDVI and Albedo
            NDVI.append(_calc_ndvi(data, bandnames, sensor))
            ALBEDO.append(_calc_albedo(data, bandnames, sensor))

            del data

        NDVI = np.asarray(NDVI)
        ALBEDO = np.asarray(ALBEDO)

        print('Applying SavGol filter...')
        NDVI_smooth, _ = savgol_reconstruct(NDVI)
        del NDVI
        ALBEDO_smooth, _ = savgol_reconstruct(ALBEDO, invert=True)
        del ALBEDO

        dekadal_dates = [_get_dekadal_date(date) for date in sorted_dates]
        unique_dekadal_dates = np.unique(dekadal_dates)

        meta.update({'dtype': 'float64',
                     'nodata': np.nan})

        # merge dekadal images
        for dekadal_date in unique_dekadal_dates:

            datestring = dekadal_date.strftime('%Y%m%d')
            idx = np.argwhere(np.isin(dekadal_dates, dekadal_date))[:, 0]

            if idx.shape[0] == 0:
                break

            ndvi_composite_array = NDVI_smooth[idx, ...]
            albedo_composite_array = ALBEDO_smooth[idx, ...]

            ndvi_dekadal_composite = np.nanmean(ndvi_composite_array, axis=0)
            albedo_dekadal_composite = np.nanmean(albedo_composite_array, axis=0)

            # save dekadal images
            ndvi_filename = output_folder / Path(datestring) / Path('NDVI_' + datestring + '.tif')
            albedo_filename = output_folder / Path(datestring) / Path('ALBEDO_' + datestring + '.tif')

            if not os.path.exists(output_folder / Path(datestring)):
                os.makedirs(output_folder / Path(datestring))

            if window == block_windows[1]:
                with rio.open(str(ndvi_filename), 'w', **meta) as dst:
                    dst.write(ndvi_dekadal_composite, 1, window=window)
                with rio.open(str(albedo_filename), 'w', **meta) as dst:
                    dst.write(albedo_dekadal_composite, 1, window=window)
            else:
                with rio.open(str(ndvi_filename), 'r+') as dst:
                    dst.write(ndvi_dekadal_composite, 1, window=window)
                with rio.open(str(albedo_filename), 'r+') as dst:
                    dst.write(albedo_dekadal_composite, 1, window=window)


def _image_path(directory, filename):
    sensor = str(filename.split('_')[0])
    if sensor == 'LE07':
        filename = filename+'_gap-filled'
        filename = str(directory/Path('L7')/Path(filename)) + '.tif'
    elif sensor == 'LC08':
        filename = str(directory/Path('L8')/Path(filename)) + '.tif'
    return filename, sensor


def _calc_ndvi(data, bandnames, sensor):
    if sensor == 'LE07':
        bands = ['SR_B3', 'SR_B4']

    elif sensor == 'LC08':
        bands = ['SR_B4', 'SR_B5']

    red = data[bandnames.index(bands[0]), ...].astype(np.float)
    nir = data[bandnames.index(bands[1]), ...].astype(np.float)

    ndvi = (nir - red) / (nir + red)

    # remove too large or too small values
    ndvi[ndvi > 1] = np.nan
    ndvi[ndvi < -1] = np.nan

    return ndvi


def _calc_albedo(data, bandnames, sensor):
    albedo_Mp = 2.75e-5  # multiplicative scaling factor for Collection 2
    albedo_Ap = -0.2  # additive scaling factor for Collection 2

    # ESUN values: [Blue, Green, Red, NIR, SWIR-1, SWIR-2]
    if sensor == 'LE07':
        ESUN_values = np.array([1970, 1842, 1547, 1044, 225.7, 82.06])
        bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']

    elif sensor == 'LC08':
        ESUN_values = np.array([1991, 1812, 1549, 972.6, 214.7, 80.7])
        bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']

    band_idx = [bandnames.index(band) for band in bands]

    BGRNS = albedo_Mp * data[band_idx, ...] + albedo_Ap

    albedo = np.sum(BGRNS * np.expand_dims(ESUN_values, (1, 2)), axis=0) / np.sum(ESUN_values)

    # remove too large or too small values
    albedo[albedo > 1] = np.nan
    albedo[albedo < 0] = np.nan

    return albedo


def _apply_nspi(landsat_dir, filename_list, bandnames_list, overwrite=False):
    L7_idx = [i for i, file in enumerate(filename_list) if str(file.split('_')[0]) == 'LE07']
    L8_idx = [i for i, file in enumerate(filename_list) if str(file.split('_')[0]) == 'LC08']

    L7_tifs = [filename_list[idx] + '.tif' for idx in L7_idx]
    L8_tifs = [filename_list[idx] + '.tif' for idx in L8_idx]

    L7_bandnames = [bandnames_list[idx] for idx in L7_idx]
    L8_bandnames = [bandnames_list[idx] for idx in L8_idx]

    L7_dates = [datetime.strptime(str(f).split('_')[3], '%Y%m%d') for f in L7_tifs]
    L8_dates = [datetime.strptime(str(f).split('_')[3], '%Y%m%d') for f in L8_tifs]

    L7_bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']
    L8_bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']

    target_folder = Path('L7')

    print('Fill gaps (NSPI)...')
    for target_idx, target_file in enumerate(tqdm(L7_tifs)):

        target_date = L7_dates[target_idx]
        date_diff = [np.abs(target_date - input_date) for input_date in L7_dates + L8_dates]

        output_filename = str(landsat_dir / target_folder / Path(target_file).stem) + '_gap-filled.tif'

        if os.path.isfile(output_filename) and not overwrite:
            continue

        input_idx = np.argpartition(date_diff, 1)[1]

        input_file = (L7_tifs + L8_tifs)[input_idx]

        input_sensor = str(input_file.split('_')[0])

        if input_sensor == 'LE07':
            input_folder = Path('L7')
            input_bands = L7_bands
        elif input_sensor == 'LC08':
            input_folder = Path('L8')
            input_bands = L8_bands

        # open target image
        with rio.open(str(landsat_dir / target_folder / Path(target_file))) as master_src:
            meta = master_src.profile
            target_image = master_src.read()

        master_dict = {'transform': master_src.transform, 'height': master_src.height, 'width': master_src.width,
                       'crs': master_src.crs}

        if date_diff[input_idx] > timedelta(32):
            # save results as tif
            with rio.open(output_filename, 'w', **meta) as dst:
                dst.write(target_image)

            continue

        # open input image
        with rio.open(str(landsat_dir / input_folder / Path(input_file))) as slave_src:
            with rio.vrt.WarpedVRT(slave_src, **master_dict) as vrt:
                input_image = vrt.read()

        target_band_idx = [L7_bandnames[target_idx].index(band) for band in L7_bands]
        input_band_idx = [(L7_bandnames + L8_bandnames)[input_idx].index(band) for band in input_bands]

        # reshape images for correct NSPI format
        target_image_reshaped = np.transpose(target_image[target_band_idx, ...], axes=[1, 2, 0])
        input_image_reshaped = np.transpose(input_image[input_band_idx, ...], axes=[1, 2, 0])

        # add nans instead of nodata
        target_image_reshaped = np.where(target_image_reshaped == -9999, np.nan, target_image_reshaped)
        input_image_reshaped = np.where(input_image_reshaped == -9999, np.nan, input_image_reshaped)

        # calculate the missing pixels mask
        missing_pixels_mask = np.isnan(target_image_reshaped[..., 0]) & ~np.isnan(input_image_reshaped[..., 0])

        # get the bitmask
        target_pixel_qa = target_image[L7_bandnames[target_idx].index('QA_PIXEL'), ...]
        target_pixel_qa = np.where(target_pixel_qa == -9999, 0, target_pixel_qa)

        input_pixel_qa = input_image[(L7_bandnames + L8_bandnames)[input_idx].index('QA_PIXEL'), ...]
        input_pixel_qa = np.where(input_pixel_qa == -9999, 0, input_pixel_qa)

        # calculate cloud mask
        target_pixel_cloudmask = _landsat_cloudmask(target_pixel_qa)
        input_pixel_cloudmask = _landsat_cloudmask(input_pixel_qa)

        # mask clouds
        for n in range(0, len(L7_bands)):
            target_image_reshaped[..., n] = np.where(target_pixel_cloudmask, np.nan, target_image_reshaped[..., n])
            input_image_reshaped[..., n] = np.where(input_pixel_cloudmask, np.nan, input_image_reshaped[..., n])

        out_image_reshaped = nspi(input_image_reshaped, target_image_reshaped, missing_pixels_mask,
                                  num_classes=5, required_pixels=20, max_window_size=15)

        # replace the gap-filled bands in the original image
        out_image = target_image.astype(np.float64)
        out_image[target_band_idx, ...] = np.transpose(out_image_reshaped, axes=[2, 0, 1])

        # convert nan-values to nodata
        out_image = np.where(np.isnan(out_image), -9999, out_image).astype(np.int16)

        # save results as tif
        with rio.open(output_filename, 'w', **meta) as dst:
            dst.write(out_image)


# merges individual landsat bands and saves as single tif
def _merge_and_save_landsat(directory, delete_input=False, overwrite=False):

    master_file = list(directory.glob('*QA_PIXEL.TIF'))[0]
    slave_files = [f for f in list(directory.glob('*.TIF')) if 'QA_PIXEL' not in str(f)]

    output_filename = master_file.parents[1] / Path(
        '_'.join(master_file.stem.split('_')[0:-2]) + '.tif')

    band_names = ['_'.join(master_file.stem.split('_')[-2:])]

    if os.path.isfile(output_filename) and not overwrite:
        # for each slave-file, open and append to array
        for i, file in enumerate(slave_files):
            band_names.append('_'.join(file.stem.split('_')[-2:]))
        return output_filename.stem, band_names

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
            with rio.vrt.WarpedVRT(slave_src, **master_dict) as vrt:
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

    with rio.open(output_filename, 'w', **meta) as dst:
        dst.write(data)

    if delete_input:
        shutil.rmtree(directory)

    return output_filename.stem, band_names


# unpack and saves *.tar files
def _unpack_and_save(file, delete_input=False, overwrite=False):
    path = Path(file).parent / Path(Path(file).stem).stem

    # If folder exists assume that tar has already been uncompressed
    if os.path.isdir(path):
        if overwrite:
            path.rmdir()
        else:
            return

    os.mkdir(path)

    tar = tarfile.open(file, "r:*")
    tar.extractall(path=path)
    tar.close()

    if delete_input:
        os.remove(file)


# returns cloud mask given quality band as input
# Updated to match Collection 2 QA_PIXEL
# https://www.usgs.gov/media/files/landsat-8-9-olitirs-collection-2-level-2-data-format-control-book
# https://www.usgs.gov/media/files/landsat-7-etm-collection-2-level-2-data-format-control-book
def _landsat_cloudmask(quality_band):
    # if clouds (bit 3) and low/medium/high probability (bit 8 and 9) then clouds
    clouds = ((quality_band & (1 << 3)) > 1) & ((quality_band & ((1 << 8) | (1 << 9))) > 1)
    # if shadows (bit 4) and low/medium/high probability shadows (bit 10 and 11) then shadows
    shadows = ((quality_band & (1 << 4)) > 1) & ((quality_band & ((1 << 10) | (1 << 11))) > 1)
    # if cirrus (bit 2) and low/medium/high probability shadows (bit 14 and 15) then shadows
    cirrus = ((quality_band & (1 << 2)) > 1) & ((quality_band & ((1 << 14) | (1 << 15))) > 1)

    return np.logical_or.reduce((clouds, shadows, cirrus))


