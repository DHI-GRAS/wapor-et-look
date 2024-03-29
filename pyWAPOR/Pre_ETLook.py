# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Thu Feb 21 18:57:09 2019
"""
import os
import shutil
import datetime
from osgeo import gdal
import pandas as pd
import numpy as np
import rasterio as rio

import pyWAPOR
import pyWAPOR.Functions.Processing_Functions as PF

from pyWAPOR.Functions.Swets_Filter import swets_filter

from pathlib import Path

def prepare_level1(output_folder, startdate, enddate, latlim, lonlim, username, password,
                   landcover="GlobCover"):

    # Define the input folders
    folders_input_RAW = os.path.join(output_folder, "RAW")
    folder_input_ETLook = os.path.join(output_folder, "ETLook_input", "level_1")

    # Create folders if not exists
    if not os.path.exists(folders_input_RAW):
        os.makedirs(folders_input_RAW)
    if not os.path.exists(folder_input_ETLook):
        os.makedirs(folder_input_ETLook)

    # Define the dates
    dates = pd.date_range(startdate, enddate, freq = "D")

    # Only NDVI and albedo are different between Level 1 and Level 2 so download and process them
    # first

    # Extend the days for NDVI data with +8 from both sides
    startdate_NDVI = datetime.datetime.strptime(startdate, "%Y-%m-%d") - datetime.timedelta(days = 8)
    enddate_NDVI = datetime.datetime.strptime(enddate, "%Y-%m-%d") + datetime.timedelta(days = 8)

    startdate_NDVI_str = datetime.datetime.strftime(startdate_NDVI, "%Y-%m-%d")
    enddate_NDVI_str = datetime.datetime.strftime(enddate_NDVI, "%Y-%m-%d")

    # Download NDVI and ALBEDO data
    pyWAPOR.Collect.MOD13.NDVI(folders_input_RAW, startdate_NDVI_str, enddate_NDVI_str, latlim, lonlim, username, password)
    pyWAPOR.Collect.MYD13.NDVI(folders_input_RAW, startdate_NDVI_str, enddate_NDVI_str, latlim, lonlim, username, password)
    pyWAPOR.Collect.MCD43.ALBEDO(folders_input_RAW, startdate, enddate, latlim, lonlim, username, password)

    # Create the inputs of MODIS NDVI and albedo for all the Dates
    template_file = None
    for date in dates:

        try:
            # Define output folder
            folder_input_ETLook_Date = os.path.join(folder_input_ETLook, "%d%02d%02d" %(date.year, date.month, date.day))
            if not os.path.exists(folder_input_ETLook_Date):
                os.makedirs(folder_input_ETLook_Date)

            # Find nearest date for NDVI
            Startdate_year = "%d-01-01" %date.year
            Enddate_year = "%d-12-31" %date.year

            # Create MODIS NDVI dataset
            Dates_eight_daily_year = pd.date_range(Startdate_year, Enddate_year, freq = "8D")

            # find nearest NDVI date
            Date_nearest = min(Dates_eight_daily_year, key=lambda Dates_eight_daily_year: abs(Dates_eight_daily_year - date))

            # Create NDVI files for ETLook

            # try MOD13 and MYD13
            NDVI_file = os.path.join(folder_input_ETLook_Date, "NDVI_%d%02d%02d.tif" %(date.year, date.month, date.day))
            if not os.path.exists(NDVI_file):
                folder_RAW_file_NDVI = os.path.join(folders_input_RAW, "MODIS", "{v}13")
                filename_NDVI = "NDVI_{v}13Q1_-_16-daily_%d.%02d.%02d.tif" %(Date_nearest.year, Date_nearest.month, Date_nearest.day)

                if os.path.exists(os.path.join(folder_RAW_file_NDVI, filename_NDVI).format(v="MOD")):
                    shutil.copy(os.path.join(folder_RAW_file_NDVI, filename_NDVI).format(v="MOD"),
                                folder_input_ETLook_Date)
                    os.rename(os.path.join(folder_input_ETLook_Date, filename_NDVI).format(v="MOD"), NDVI_file)

                elif os.path.exists(os.path.join(folder_RAW_file_NDVI, filename_NDVI).format(v="MYD")):
                    shutil.copy(os.path.join(folder_RAW_file_NDVI, filename_NDVI).format(v="MYD"),
                                folder_input_ETLook_Date)
                    os.rename(os.path.join(folder_input_ETLook_Date, filename_NDVI).format(v="MYD"), NDVI_file)

                else:
                    print("NDVI is not available for date: %d%02d%02d" %(date.year, date.month, date.day))

            # Get example files
            if not template_file and os.path.exists(NDVI_file):
                template_file = NDVI_file
                dest_ex = gdal.Open(NDVI_file)
                geo_ex = dest_ex.GetGeoTransform()
                proj_ex = dest_ex.GetProjection()
                dest_ex = None

            # Create ALBEDO files for ETLook

            # try MCD43
            ALBEDO_file = os.path.join(folder_input_ETLook_Date, "ALBEDO_%d%02d%02d.tif" %(date.year, date.month, date.day))
            if not os.path.exists(ALBEDO_file):
                folder_RAW_file_ALBEDO = os.path.join(folders_input_RAW, "MODIS", "MCD43")
                filename_ALBEDO = "Albedo_MCD43A3_-_daily_%d.%02d.%02d.tif" %(date.year, date.month, date.day)

                if os.path.exists(os.path.join(folder_RAW_file_ALBEDO, filename_ALBEDO)):
                    destalbedo = PF.reproject_dataset_example(os.path.join(folder_RAW_file_ALBEDO, filename_ALBEDO), template_file, method=1)
                    albedo = destalbedo.GetRasterBand(1).ReadAsArray()
                    albedo[albedo<=-0.4] = -9999
                    PF.Save_as_tiff(ALBEDO_file, albedo, geo_ex, proj_ex)
                else:
                    print("ALBEDO is not available for date: %d%02d%02d" %(date.year, date.month, date.day))

        except:
            print("No ETLook input dataset for %s" %date)

    # Download and create all other Level 1 inputs
    prepare_level1_level2(output_folder, startdate, enddate, latlim, lonlim, username, password,
                          template_file, "level_1", LandCover=landcover)


def prepare_level2(output_folder, startdate, enddate, latlim, lonlim, username_vito, password_vito,
                   username_earthdata, password_earthdata, landcover="GlobCover"):

    # Define the input folders
    folder_input_RAW = os.path.join(output_folder, "RAW")
    folder_input_ETLook = os.path.join(output_folder, "ETLook_input", "level_2")

    # Create folders if not exists
    if not os.path.exists(folder_input_RAW):
        os.makedirs(folder_input_RAW)
    if not os.path.exists(folder_input_ETLook):
        os.makedirs(folder_input_ETLook)

    # Define the dates
    dates = pd.date_range(startdate, enddate, freq="D")

    # Only NDVI and albedo are different between Level 1 and Level 2 so download and process them
    # first

    # Download 5-day PROBA-V composites
    pyWAPOR.Collect.PROBAV.PROBAV_S5(folder_input_RAW, startdate, enddate, latlim, lonlim,
                                     username_vito, password_vito)

    # Create PROBA-V 10-day composites

    # get the dates at which we have data
    folder_input_albedo = (Path(folder_input_RAW) / Path('PROBAV/ALBEDO'))
    folder_input_ndvi = (Path(folder_input_RAW) / Path('PROBAV/NDVI'))

    s5_startdates = [datetime.datetime.strptime(file.stem[-10:], '%Y-%m-%d') for file in
                     list(folder_input_albedo.glob('*.tif'))]
    s5_centerdates = [startdate + datetime.timedelta(days=3) for startdate in s5_startdates]

    # allocate large timeseries array
    albedo_dekadal_timeseries_list = []
    ndvi_dekadal_timeseries_list = []
    datetime_timeseries = []
    first_iteration = True

    # loop over all dates in daterange
    for year in np.unique(dates.year):
        for month in np.unique(dates.month):
            for day in [1, 11, 21]:
                # find all files in relevant daterange
                current_datetime = datetime.datetime(year, month, day)
                composite_daterange = _get_dekadal_daterange(current_datetime)
                matching_dates = list(set(composite_daterange) & set(s5_centerdates))  # we match on center dates in S5

                albedo_composite_list = []
                ndvi_composite_list = []

                # open files
                for s5_centerdate in matching_dates:

                    s5_startdate = s5_centerdate - datetime.timedelta(days=3)  # go back to startdate in S5
                    s5_datestring = s5_startdate.strftime('%Y-%m-%d')

                    albedo_filename = list(folder_input_albedo.glob(f'*{s5_datestring}*.tif'))[0]
                    ndvi_filename = list(folder_input_ndvi.glob(f'*{s5_datestring}*.tif'))[0]

                    # open file and append to mosaic-ready array
                    with rio.open(str(albedo_filename)) as src:
                        albedo_composite_list.append(src.read().squeeze())
                        # save a template of the tif metadata
                        if first_iteration:
                            meta = src.profile
                            first_iteration = False
                    with rio.open(str(ndvi_filename)) as src:
                        ndvi_composite_list.append(src.read().squeeze())

                if albedo_composite_list:
                    ndvi_composite_array = np.asarray(ndvi_composite_list)
                    albedo_composite_array = np.asarray(albedo_composite_list)

                    # remove pixels with nan in all bands (because np.nanargmax can't handle these)
                    nan_mask = np.repeat(np.expand_dims(np.mean(np.isnan(ndvi_composite_array), axis=0) == 1,
                                                        axis=0), ndvi_composite_array.shape[0], axis=0)

                    ndvi_composite_array[nan_mask] = 0
                    albedo_composite_array[nan_mask] = 0

                    # contrained max-composite
                    composite_idx = np.nanargmax(ndvi_composite_array, axis=0)

                    ndvi_dekadal_composite = _numeric_nd_indexing(ndvi_composite_array, composite_idx)
                    albedo_dekadal_composite = _numeric_nd_indexing(albedo_composite_array, composite_idx)

                    # add nans again
                    ndvi_dekadal_composite[nan_mask[0, ...]] = np.nan
                    albedo_dekadal_composite[nan_mask[0, ...]] = np.nan

                    # append to timeseries
                    albedo_dekadal_timeseries_list.append(albedo_dekadal_composite.astype(albedo_composite_array.dtype))
                    ndvi_dekadal_timeseries_list.append(ndvi_dekadal_composite.astype(ndvi_composite_array.dtype))
                    datetime_timeseries.append(current_datetime)

    # transpose for correct format for swets (ROWxCOLxTIME)
    albedo_dekadal_timeseries = np.transpose(np.asarray(albedo_dekadal_timeseries_list), axes=[1, 2, 0])
    ndvi_dekadal_timeseries = np.transpose(np.asarray(ndvi_dekadal_timeseries_list), axes=[1, 2, 0])

    albedo_dekadal_timeseries_smoothed = swets_filter(albedo_dekadal_timeseries, invert=True)
    ndvi_dekadal_timeseries_smoothed = swets_filter(ndvi_dekadal_timeseries)

    # save one composite for each individual date
    template_file = None
    for date in dates:
        # get the dekadal date corresponding with current date
        s10_startdate = _get_dekadal_date(date)

        # get the dekadal image corresponding with current date
        current_dekadal_date_idx = datetime_timeseries.index(s10_startdate)

        current_dekadal_albedo_array = albedo_dekadal_timeseries_smoothed[..., current_dekadal_date_idx]
        current_dekadal_ndvi_array = ndvi_dekadal_timeseries_smoothed[..., current_dekadal_date_idx]

        datestring = date.strftime('%Y%m%d')

        folder_input_ETLook_Date = Path(folder_input_ETLook) / Path(datestring)

        if not os.path.exists(folder_input_ETLook_Date):
            os.makedirs(folder_input_ETLook_Date)

        albedo_filename = folder_input_ETLook_Date / Path(f'ALBEDO_{datestring}.tif')
        ndvi_filename = folder_input_ETLook_Date / Path(f'NDVI_{datestring}.tif')

        # save tif
        with rio.open(str(albedo_filename), 'w', **meta) as dst:
            dst.write(current_dekadal_albedo_array, 1)

        with rio.open(str(ndvi_filename), 'w', **meta) as dst:
            dst.write(current_dekadal_ndvi_array, 1)

        if not template_file and os.path.exists(ndvi_filename):
            template_file = str(ndvi_filename)
            dest_ex = gdal.Open(str(ndvi_filename))
            geo_ex = dest_ex.GetGeoTransform()
            proj_ex = dest_ex.GetProjection()
            dest_ex = None

    # Download and create all other Level 2 inputs
    prepare_level1_level2(output_folder, startdate, enddate, latlim, lonlim, username_earthdata,
                          password_earthdata, template_file, "level_2", LandCover=landcover)


# Preparation of input data which is common for Level 1 and Level 2
def prepare_level1_level2(output_folder, Startdate, Enddate, latlim, lonlim, username, password,
                          template_file, level, LandCover="GlobCover"):

    # Define the input folders
    folders_input_RAW = os.path.join(output_folder, "RAW")
    folder_input_ETLook = os.path.join(output_folder, "ETLook_input", level)

    # Create folders if not exists
    if not os.path.exists(folders_input_RAW):
        os.makedirs(folders_input_RAW)
    if not os.path.exists(folder_input_ETLook):
        os.makedirs(folder_input_ETLook)

    # Define the dates
    Dates = pd.date_range(Startdate, Enddate, freq = "D")

    ######################### Download LST MODIS data #############################

    # Download LST data
    pyWAPOR.Collect.MOD11.LST(folders_input_RAW, Startdate, Enddate, latlim, lonlim, username, password)
    pyWAPOR.Collect.MYD11.LST(folders_input_RAW, Startdate, Enddate, latlim, lonlim, username, password)
    Combine_LST(folders_input_RAW, Startdate, Enddate)

    ######################## Download Rainfall Data ###############################

    # Download CHIRPS data
    pyWAPOR.Collect.CHIRPS.daily(folders_input_RAW, Startdate, Enddate, latlim, lonlim)

    ########################### Download DEM data #################################

    # Download DEM data
    pyWAPOR.Collect.SRTM.DEM(folders_input_RAW, latlim, lonlim)

    ############################ Download Landuse #################################
    if LandCover == "GlobCover":
        # Download Globcover data
        pyWAPOR.Collect.Globcover.Landuse(folders_input_RAW, latlim, lonlim)
    if LandCover == "WAPOR":
        # Download Globcover data
        pyWAPOR.Collect.WAPOR.LandCover(folders_input_RAW, "%s-01-01" % (Startdate.split("-")[0]), "%s-12-31" % (Enddate.split("-")[0]), latlim, lonlim)
    ############### Loop over days for the dynamic data ###############################

    # Create the inputs of MODIS for all the Dates
    dest_ex = gdal.Open(template_file)
    geo_ex = dest_ex.GetGeoTransform()
    proj_ex = dest_ex.GetProjection()
    size_x_ex = dest_ex.RasterXSize
    size_y_ex = dest_ex.RasterYSize
    for Date in Dates:

        try:
            # Define output folder
            folder_input_ETLook_Date = os.path.join(folder_input_ETLook, "%d%02d%02d" %(Date.year, Date.month, Date.day))
            if not os.path.exists(folder_input_ETLook_Date):
                os.makedirs(folder_input_ETLook_Date)

            # Create LST files for ETLook
            LST_file = os.path.join(folder_input_ETLook_Date, "LST_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            Time_file = os.path.join(folder_input_ETLook_Date, "Time_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            VZA_file = os.path.join(folder_input_ETLook_Date, "VZA_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            if not os.path.exists(LST_file):
                folder_RAW_file_LST = os.path.join(folders_input_RAW, "MODIS", "LST")
                filename_LST = "LST_MCD11A1_K_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)
                filename_Time = "Time_MCD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)
                filename_VZA = "VZA_MCD11A1_degrees_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)
                if os.path.exists(os.path.join(folder_RAW_file_LST, filename_LST)):
                    destLST = PF.reproject_dataset_example(os.path.join(folder_RAW_file_LST, filename_LST), template_file, method=2)
                    LST = destLST.GetRasterBand(1).ReadAsArray()
                    LST[LST==0.0] = -9999
                    PF.Save_as_tiff(LST_file, LST, geo_ex, proj_ex)

                    destTime = PF.reproject_dataset_example(os.path.join(folder_RAW_file_LST, filename_Time), template_file, method=1)
                    Time = destTime.GetRasterBand(1).ReadAsArray()
                    Time[Time==0.0] = -9999
                    PF.Save_as_tiff(Time_file, Time, geo_ex, proj_ex)
                    
                    destVZA = PF.reproject_dataset_example(os.path.join(folder_RAW_file_LST, filename_VZA), template_file, method=1)
                    VZA = destVZA.GetRasterBand(1).ReadAsArray()
                    VZA[LST==0.0] = -9999
                    PF.Save_as_tiff(VZA_file, VZA, geo_ex, proj_ex)
                else:
                    print("LST is not available for date: %d%02d%02d" %(Date.year, Date.month, Date.day))

            ####################### Create lat and lon rasters ############################

            Lon_file = os.path.join(folder_input_ETLook_Date, "Lon_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            Lat_file = os.path.join(folder_input_ETLook_Date, "Lat_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            if not (os.path.exists(Lon_file) or os.path.exists(Lat_file)):
                lon_deg = np.array([geo_ex[0] + np.arange(0,size_x_ex) * geo_ex[1]]*size_y_ex)
                lat_deg = np.array([geo_ex[3] + np.arange(0,size_y_ex) * geo_ex[5]]*size_x_ex).transpose()

                # save as tiff
                PF.Save_as_tiff(Lon_file, lon_deg, geo_ex, proj_ex)
                PF.Save_as_tiff(Lat_file, lat_deg, geo_ex, proj_ex)

            ########################## Create Time rasters ################################

            # calculate overall time
            dest_time = gdal.Open(Time_file)
            Time_array = dest_time.GetRasterBand(1).ReadAsArray()
            Time_array[Time_array==-9999] = np.nan
            dtime = np.nanmean(Time_array)
            if np.isnan(dtime):
                dtime = 12
            NowTime = datetime.datetime(Date.year, Date.month, Date.day, int(np.floor(dtime)), int((dtime - np.floor(dtime))*60))

            # Get DOY
            doy = int(Date.strftime("%j"))

            ####################### Create DEM rasters ############################

            # Create DEM files for ETLook
            DEM_file = os.path.join(folder_input_ETLook_Date, "DEM.tif")
            if not os.path.exists(DEM_file):
                folder_RAW_file_DEM = os.path.join(folders_input_RAW, "SRTM", "DEM")
                filename_DEM = "DEM_SRTM_m_3s.tif"
                if os.path.exists(os.path.join(folder_RAW_file_DEM, filename_DEM)):
                    destDEM = PF.reproject_dataset_example(os.path.join(folder_RAW_file_DEM, filename_DEM), template_file, method=4)
                    DEM = destDEM.GetRasterBand(1).ReadAsArray()
                    PF.Save_as_tiff(DEM_file, DEM, geo_ex, proj_ex)

                else:
                    print("DEM is not available")

            ##################### Calculate SLope and Aspect ##############################
            Slope_file = os.path.join(folder_input_ETLook_Date, "Slope.tif")
            Aspect_file = os.path.join(folder_input_ETLook_Date, "Aspect.tif")
            if not (os.path.exists(Slope_file) and os.path.exists(Aspect_file)):

                # open DEM
                destDEM = gdal.Open(DEM_file)
                DEM = destDEM.GetRasterBand(1).ReadAsArray()

                # constants
                pixel_spacing = 1000
                deg2rad = np.pi / 180.0  # Factor to transform from degree to rad
                rad2deg = 180.0 / np.pi  # Factor to transform from rad to degree

                # Calculate slope
                x, y = np.gradient(DEM, pixel_spacing, pixel_spacing)
                hypotenuse_array = np.hypot(x,y)
                slope = np.arctan(hypotenuse_array) * rad2deg

                # calculate aspect
                aspect = np.arctan2(y/pixel_spacing, -x/pixel_spacing) * rad2deg
                aspect = 180 + aspect

                # Save as tiff files
                PF.Save_as_tiff(Slope_file, slope, geo_ex, proj_ex)
                PF.Save_as_tiff(Aspect_file, aspect, geo_ex, proj_ex)
            ######################### Create Rainfall file ################################

            P_file = os.path.join(folder_input_ETLook_Date, "Precipitation_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            if not os.path.exists(P_file):
                folder_RAW_file_P = os.path.join(folders_input_RAW, "Precipitation", "CHIRPS")
                filename_P = "P_CHIRPS.v2.0_mm-day-1_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)
                destP = PF.reproject_dataset_example(os.path.join(folder_RAW_file_P, filename_P), template_file, method=2)
                P = destP.GetRasterBand(1).ReadAsArray()
                PF.Save_as_tiff(P_file, P, geo_ex, proj_ex)


            ############################# Download METEO ##################################

            # Define the startdates for the METEO
            StartTime = datetime.datetime(Date.year, Date.month, Date.day, 0, 0)
            EndTime = datetime.datetime(Date.year, Date.month, Date.day, 23, 59)

            if (Date >= datetime.datetime(2016,1,1) and Date < datetime.datetime(2017,12,1)):
                # find nearest Meteo time
                DateTime = pd.date_range(StartTime, EndTime, freq="H") + pd.offsets.Minute(30)
                Time_nearest = min(DateTime, key=lambda DateTime: abs(DateTime - NowTime))
                Period = np.argwhere(DateTime ==Time_nearest)[0][0] + 1

            else:
                # find nearest Meteo time
                DateTime = pd.date_range(StartTime, EndTime, freq="3H") + pd.offsets.Minute(90)
                Time_nearest = min(DateTime, key=lambda DateTime: abs(DateTime - NowTime))
                Period = np.argwhere(DateTime ==Time_nearest)[0][0] + 1


            # Download METEO data
            if Date < datetime.datetime(2016,1,1):
                pyWAPOR.Collect.MERRA.daily(folders_input_RAW, ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp'], StartTime, EndTime, latlim, lonlim)
                pyWAPOR.Collect.MERRA.three_hourly(folders_input_RAW, ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp'], StartTime, EndTime, latlim, lonlim, [int(Period)])
                str_METEO = "MERRA"
                inst_name = "three_hourly"
                day_name = "daily"
                hour_steps = 3
                file_time_inst = "3-hourly"

            elif (Date >= datetime.datetime(2016,1,1) and Date < datetime.datetime(2017,12,1)):
                pyWAPOR.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp'],StartTime, EndTime, latlim, lonlim, username, password)
                pyWAPOR.Collect.MERRA.hourly_MERRA2(folders_input_RAW, ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp'], StartTime, EndTime, latlim, lonlim, [int(Period)], username, password)
                str_METEO = "MERRA"
                inst_name = "hourly_MERRA2"
                day_name = "daily_MERRA2"
                hour_steps = 1
                file_time_inst = "hourly"

            else:
                pyWAPOR.Collect.GEOS.daily(folders_input_RAW, ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp'],StartTime, EndTime, latlim, lonlim)
                pyWAPOR.Collect.GEOS.three_hourly(folders_input_RAW, ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp'], StartTime, EndTime, latlim, lonlim, [int(Period)])
                str_METEO = "GEOS"
                inst_name = "three_hourly"
                day_name = "daily"
                hour_steps = 3
                file_time_inst = "3-hourly"

            # Air pressure
            pair_inst_file = os.path.join(folder_input_ETLook_Date, "Pair_inst_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            pair_inst_0_file = os.path.join(folder_input_ETLook_Date, "Pair_inst_0_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            pair_24_0_file = os.path.join(folder_input_ETLook_Date, "Pair_24_0_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))

            if not (os.path.exists(pair_inst_file) and os.path.exists(pair_inst_0_file) and os.path.exists(pair_24_0_file)):
                folder_RAW_file_pair_inst = os.path.join(folders_input_RAW, str_METEO, "Surface_Pressure", inst_name)
                folder_RAW_file_pair_inst_0 = os.path.join(folders_input_RAW, str_METEO, "Sea_Level_Pressure", inst_name)
                folder_RAW_file_pair_24_0 = os.path.join(folders_input_RAW, str_METEO, "Sea_Level_Pressure", day_name)
                HourPeriod = hour_steps * (Period - 1)

                filename_pair_inst = "ps_%s_kpa_%s_%d.%02d.%02d_H%02d.M00.tif"  %(str_METEO, file_time_inst, Date.year, Date.month, Date.day, HourPeriod)
                if os.path.exists(os.path.join(folder_RAW_file_pair_inst, filename_pair_inst)):
                    destPairInst = PF.reproject_dataset_example(os.path.join(folder_RAW_file_pair_inst, filename_pair_inst), template_file, method=2)
                    Pair_inst = destPairInst.GetRasterBand(1).ReadAsArray()
                    PF.Save_as_tiff(pair_inst_file, Pair_inst, geo_ex, proj_ex)

                else:
                    print("Pair instantenious is not available")

                filename_pair_inst_sea = "slp_%s_kpa_%s_%d.%02d.%02d_H%02d.M00.tif"  %(str_METEO, file_time_inst, Date.year, Date.month, Date.day, HourPeriod)
                if os.path.exists(os.path.join(folder_RAW_file_pair_inst_0, filename_pair_inst_sea)):
                    destPairInstSea = PF.reproject_dataset_example(os.path.join(folder_RAW_file_pair_inst_0, filename_pair_inst_sea), template_file, method=2)
                    Pair_inst_sea = destPairInstSea.GetRasterBand(1).ReadAsArray()
                    PF.Save_as_tiff(pair_inst_0_file, Pair_inst_sea, geo_ex, proj_ex)

                else:
                    print("Pair sea level instantenious is not available")

                filename_pair_24_sea = "slp_%s_kpa_daily_%d.%02d.%02d.tif"  %(str_METEO, Date.year, Date.month, Date.day)
                if os.path.exists(os.path.join(folder_RAW_file_pair_24_0, filename_pair_24_sea)):
                    destPair24Sea = PF.reproject_dataset_example(os.path.join(folder_RAW_file_pair_24_0, filename_pair_24_sea), template_file, method=2)
                    Pair_24_sea = destPair24Sea.GetRasterBand(1).ReadAsArray()
                    PF.Save_as_tiff(pair_24_0_file, Pair_24_sea, geo_ex, proj_ex)

                else:
                    print("Pair sea level daily is not available")

            # Specific Humidity
            qv_inst_file = os.path.join(folder_input_ETLook_Date, "qv_inst_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            qv_24_file = os.path.join(folder_input_ETLook_Date, "qv_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))

            if not (os.path.exists(qv_inst_file) and os.path.exists(qv_24_file)):
                folder_RAW_file_qv_inst = os.path.join(folders_input_RAW, str_METEO, "Specific_Humidity", inst_name)
                folder_RAW_file_qv_24 = os.path.join(folders_input_RAW, str_METEO, "Specific_Humidity", day_name)
                HourPeriod = hour_steps * (Period - 1)
                if str_METEO == "MERRA":
                    para = "q2m"
                else:
                    para = "qv2m"

                filename_qv_inst = "%s_%s_kg-kg-1_%s_%d.%02d.%02d_H%02d.M00.tif" %(para, str_METEO, file_time_inst, Date.year, Date.month, Date.day, HourPeriod)
                if os.path.exists(os.path.join(folder_RAW_file_qv_inst, filename_qv_inst)):
                    destqvInst = PF.reproject_dataset_example(os.path.join(folder_RAW_file_qv_inst, filename_qv_inst), template_file, method=2)
                    qv_inst = destqvInst.GetRasterBand(1).ReadAsArray()
                    PF.Save_as_tiff(qv_inst_file, qv_inst, geo_ex, proj_ex)

                else:
                    print("qv instantenious is not available")


                filename_qv_24 = "%s_%s_kg-kg-1_daily_%d.%02d.%02d.tif"  %(para, str_METEO, Date.year, Date.month, Date.day)
                if os.path.exists(os.path.join(folder_RAW_file_qv_24, filename_qv_24)):
                    destqv24 = PF.reproject_dataset_example(os.path.join(folder_RAW_file_qv_24, filename_qv_24), template_file, method=2)
                    qv_24 = destqv24.GetRasterBand(1).ReadAsArray()
                    PF.Save_as_tiff(qv_24_file, qv_24, geo_ex, proj_ex)

                else:
                    print("daily qv is not available")

            # Air temperature
            Tair_inst_file = os.path.join(folder_input_ETLook_Date, "tair_inst_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            Tair_24_file = os.path.join(folder_input_ETLook_Date, "tair_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))

            if not (os.path.exists(Tair_inst_file) and os.path.exists(Tair_24_file)):
                folder_RAW_file_tair_inst = os.path.join(folders_input_RAW, str_METEO, "Air_Temperature", inst_name)
                folder_RAW_file_tair_24 = os.path.join(folders_input_RAW, str_METEO, "Air_Temperature", day_name)
                HourPeriod = hour_steps * (Period - 1)

                filename_tair_inst = "t2m_%s_K_%s_%d.%02d.%02d_H%02d.M00.tif"  %(str_METEO, file_time_inst, Date.year, Date.month, Date.day, HourPeriod)
                if os.path.exists(os.path.join(folder_RAW_file_tair_inst, filename_tair_inst)):
                    tair_inst = lapse_rate_temp(os.path.join(folder_RAW_file_tair_inst, filename_tair_inst), DEM_file)
                    PF.Save_as_tiff(Tair_inst_file, tair_inst, geo_ex, proj_ex)

                else:
                    print("Tair instantenious is not available")


                filename_tair_24 = "t2m_%s_K_daily_%d.%02d.%02d.tif"  %(str_METEO, Date.year, Date.month, Date.day)
                if os.path.exists(os.path.join(folder_RAW_file_tair_24, filename_tair_24)):
                    tair_24 = lapse_rate_temp(os.path.join(folder_RAW_file_tair_24, filename_tair_24), DEM_file)
                    PF.Save_as_tiff(Tair_24_file, tair_24, geo_ex, proj_ex)
                # Also save the maximum and minimum daily temperatures
                Tair_24_file = os.path.join(folder_input_ETLook_Date, "tair_max_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
                filename_tair_max_24 = "t2m_%s_K_daily_max_%d.%02d.%02d.tif"  %(str_METEO, Date.year, Date.month, Date.day)
                if os.path.exists(os.path.join(folder_RAW_file_tair_24, filename_tair_max_24)):
                    tair_24 = lapse_rate_temp(os.path.join(folder_RAW_file_tair_24, filename_tair_max_24), DEM_file)
                    PF.Save_as_tiff(Tair_24_file, tair_24, geo_ex, proj_ex)
                Tair_24_file = os.path.join(folder_input_ETLook_Date, "tair_min_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
                filename_tair_min_24 = "t2m_%s_K_daily_min_%d.%02d.%02d.tif"  %(str_METEO, Date.year, Date.month, Date.day)
                if os.path.exists(os.path.join(folder_RAW_file_tair_24, filename_tair_min_24)):
                    tair_24 = lapse_rate_temp(os.path.join(folder_RAW_file_tair_24, filename_tair_min_24), DEM_file)
                    PF.Save_as_tiff(Tair_24_file, tair_24, geo_ex, proj_ex)

                else:
                    print("daily Tair is not available")

            # Wind Speed
            wind_inst_file = os.path.join(folder_input_ETLook_Date, "wind_inst_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
            wind_24_file = os.path.join(folder_input_ETLook_Date, "wind_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))

            if not (os.path.exists(wind_inst_file) and os.path.exists(wind_24_file)):
                folder_RAW_file_u2_inst = os.path.join(folders_input_RAW, str_METEO, "Eastward_Wind", inst_name)
                folder_RAW_file_u2_24 = os.path.join(folders_input_RAW, str_METEO, "Eastward_Wind", day_name)
                folder_RAW_file_v2_inst = os.path.join(folders_input_RAW, str_METEO, "Northward_Wind", inst_name)
                folder_RAW_file_v2_24 = os.path.join(folders_input_RAW, str_METEO, "Northward_Wind", day_name)

                HourPeriod = hour_steps * (Period - 1)

                filename_u2_inst = "u2m_%s_m-s-1_%s_%d.%02d.%02d_H%02d.M00.tif"  %(str_METEO, file_time_inst, Date.year, Date.month, Date.day, HourPeriod)
                filename_v2_inst = "v2m_%s_m-s-1_%s_%d.%02d.%02d_H%02d.M00.tif"  %(str_METEO, file_time_inst, Date.year, Date.month, Date.day, HourPeriod)
                if (os.path.exists(os.path.join(folder_RAW_file_u2_inst, filename_u2_inst)) and os.path.exists(os.path.join(folder_RAW_file_v2_inst, filename_v2_inst))):
                    destu2inst = PF.reproject_dataset_example(os.path.join(folder_RAW_file_u2_inst, filename_u2_inst), template_file, method=2)
                    destv2inst = PF.reproject_dataset_example(os.path.join(folder_RAW_file_v2_inst, filename_v2_inst), template_file, method=2)
                    u2_inst = destu2inst.GetRasterBand(1).ReadAsArray()
                    v2_inst = destv2inst.GetRasterBand(1).ReadAsArray()
                    wind_inst = np.sqrt(u2_inst**2 + v2_inst **2)
                    PF.Save_as_tiff(wind_inst_file, wind_inst, geo_ex, proj_ex)

                else:
                    print("Wind instantenious is not available")

                filename_u2_24 = "u2m_%s_m-s-1_daily_%d.%02d.%02d.tif"  %(str_METEO, Date.year, Date.month, Date.day)
                filename_v2_24 = "v2m_%s_m-s-1_daily_%d.%02d.%02d.tif"  %(str_METEO, Date.year, Date.month, Date.day)
                if (os.path.exists(os.path.join(folder_RAW_file_u2_24, filename_u2_24)) and os.path.exists(os.path.join(folder_RAW_file_v2_24, filename_v2_24))):
                    destu224 = PF.reproject_dataset_example(os.path.join(folder_RAW_file_u2_24, filename_u2_24), template_file, method=2)
                    destv224 = PF.reproject_dataset_example(os.path.join(folder_RAW_file_v2_24, filename_v2_24), template_file, method=2)
                    u2_24 = destu224.GetRasterBand(1).ReadAsArray()
                    v2_24 = destv224.GetRasterBand(1).ReadAsArray()
                    wind_24 = np.sqrt(u2_24**2 + v2_24 **2)
                    PF.Save_as_tiff(wind_24_file, wind_24, geo_ex, proj_ex)

                else:
                    print("daily Wind is not available")

            # Precipitable Water Vapor
            wv_inst_file = os.path.join(folder_input_ETLook_Date, "wv_inst_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))

            if not os.path.exists(wv_inst_file):
                folder_RAW_file_wv_inst = os.path.join(folders_input_RAW, str_METEO, "Total_Precipitable_Water_Vapor", inst_name)
                HourPeriod = hour_steps * (Period - 1)
                if str_METEO == "MERRA":
                    para = "tpw"
                else:
                    para = "tqv"

                filename_wv_inst = "%s_%s_mm_%s_%d.%02d.%02d_H%02d.M00.tif"  %(para, str_METEO, file_time_inst, Date.year, Date.month, Date.day, HourPeriod)
                if os.path.exists(os.path.join(folder_RAW_file_wv_inst, filename_wv_inst)):
                    destwvinst = PF.reproject_dataset_example(os.path.join(folder_RAW_file_wv_inst, filename_wv_inst), template_file, method=2)
                    wv_inst = destwvinst.GetRasterBand(1).ReadAsArray()
                    PF.Save_as_tiff(wv_inst_file, wv_inst, geo_ex, proj_ex)

                else:
                    print("Total Precipitable Water Vapour instantenious is not available")

            ##################### Calculate Landmask ##############################

            LM_file = os.path.join(folder_input_ETLook_Date, "LandMask.tif")
            Bulk_file = os.path.join(folder_input_ETLook_Date, "Bulk_Stomatal_resistance.tif")
            MaxObs_file = os.path.join(folder_input_ETLook_Date, "Maximum_Obstacle_Height.tif")
            if not (os.path.exists(LM_file) and os.path.exists(Bulk_file) and os.path.exists(MaxObs_file)):

                if LandCover == "GlobCover":
                    folder_RAW_file_LC = os.path.join(folders_input_RAW, "GlobCover", "Landuse")
                    filename_LC = "LC_GLOBCOVER_V2.3.tif"
                if LandCover == "WAPOR":
                    folder_RAW_file_LC = os.path.join(folders_input_RAW, "WAPOR", "LandCover")
                    filename_LC = "LC_WAPOR_%s.01.01.tif" %(Date.year)

                if os.path.exists(os.path.join(folder_RAW_file_LC, filename_LC)):
                    destLC = PF.reproject_dataset_example(os.path.join(folder_RAW_file_LC, filename_LC), template_file, method=1)
                    LC = destLC.GetRasterBand(1).ReadAsArray()
                    LC[np.isnan(LC)] = -9999

                    # import list with numbers to convert globcover into other maps
                    import pyWAPOR.Functions.LandCover_Converter as LCC

                    if LandCover == "GlobCover":
                        # Get conversion between globcover and landmask
                        LU_LM_Classes = LCC.Globcover_LM()
                        LU_Bulk_Classes = LCC.Globcover_Bulk()
                        LU_MaxObs_Classes = LCC.Globcover_MaxObs()

                    if LandCover == "WAPOR":
                        # Get conversion between globcover and landmask
                        LU_LM_Classes = LCC.WAPOR_LM_LM()
                        LU_Bulk_Classes = LCC.WAPOR_Bulk_Bulk()
                        LU_MaxObs_Classes = LCC.WAPOR_MaxObs_MaxObs()

                    # Create Array for LandMask
                    LM = np.ones([size_y_ex, size_x_ex]) * np.nan
                    Bulk = np.ones([size_y_ex, size_x_ex]) * np.nan
                    MaxObs = np.ones([size_y_ex, size_x_ex]) * np.nan

                    # Create LandMask
                    for LU_LM_Class in LU_LM_Classes.keys():
                        Value_LM = LU_LM_Classes[LU_LM_Class]
                        Value_Bulk = LU_Bulk_Classes[LU_LM_Class]
                        Value_MaxObs = LU_MaxObs_Classes[LU_LM_Class]
                        LM[LC == LU_LM_Class] = Value_LM
                        Bulk[LC == LU_LM_Class] = Value_Bulk
                        MaxObs[LC  == LU_LM_Class] = Value_MaxObs

                    # Save as tiff files
                    PF.Save_as_tiff(LM_file, LM, geo_ex, proj_ex)
                    PF.Save_as_tiff(Bulk_file, Bulk, geo_ex, proj_ex)
                    PF.Save_as_tiff(MaxObs_file, MaxObs, geo_ex, proj_ex)

                else:
                    print("LandCover is not available")

            ########################### Download amplitude ################################

            pyWAPOR.Collect.MERRA.yearly_T_Amplitude(folders_input_RAW, [Date.year],latlim, lonlim)

            # yearly amplitude temperature air
            Tair_amp_file = os.path.join(folder_input_ETLook_Date, "Tair_amp_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))

            if not os.path.exists(Tair_amp_file):
                folder_RAW_file_Tair_amp = os.path.join(folders_input_RAW, "MERRA", "Temperature_Amplitude", "yearly")

                filename_Tair_amp = "Tamp_MERRA_K_yearly_%d.01.01.tif"  %(Date.year)
                if os.path.exists(os.path.join(folder_RAW_file_Tair_amp, filename_Tair_amp)):
                    desttairamp = PF.reproject_dataset_example(os.path.join(folder_RAW_file_Tair_amp, filename_Tair_amp), template_file, method=2)
                    tair_amp = desttairamp.GetRasterBand(1).ReadAsArray()
                    PF.Save_as_tiff(Tair_amp_file, tair_amp, geo_ex, proj_ex)

                else:
                    print("Yearly Tair amplitude is not available")

            ######################## Download Transmissivity ##############################

            # Download MSGCPP data
            if Date < datetime.datetime(2016,1,1):
                pyWAPOR.Collect.MERRA.daily(folders_input_RAW, ['swgnet'],StartTime, EndTime, latlim, lonlim)
                str_TRANS = "MERRA"
                day_name = "daily"

            # MSG CPP archive seems to go back only 3 calendar years
            elif (Date >= datetime.datetime(2016,1,1) and Date < datetime.datetime(2019,1,1)):
                pyWAPOR.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['swgnet'],StartTime, EndTime, latlim, lonlim, username, password)
                str_TRANS = "MERRA"
                day_name = "daily_MERRA2"

            else:
                pyWAPOR.Collect.MSGCPP.SDS(folders_input_RAW, StartTime, EndTime, latlim, lonlim)
                str_TRANS = "MSGCPP"
                day_name = "daily"

            # yearly amplitude temperature air
            Trans_file = os.path.join(folder_input_ETLook_Date, "Trans_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))

            if not os.path.exists(Trans_file):

                # Calculate the extraterrestrial daily radiation
                destLat = gdal.Open(Lat_file)
                lat = destLat.GetRasterBand(1).ReadAsArray()

                Gsc = 1367        # Solar constant (W / m2)
                deg2rad = np.pi / 180.0
                # Computation of Hour Angle (HRA = w)
                B = 360./365 * (doy-81)           # (degrees)
                # Computation of cos(theta), where theta is the solar incidence angle
                # relative to the normal to the land surface
                delta=np.arcsin(np.sin(23.45*deg2rad)*np.sin(np.deg2rad(B))) # Declination angle (radians)

                phi = lat * deg2rad                                     # latitude of the pixel (radians)

                dr = 1 + 0.033 * np.cos(doy*2*np.pi/365)

                # Daily 24 hr radiation - For flat terrain only !
                ws_angle = np.arccos(-np.tan(phi)*np.tan(delta))   # Sunset hour angle ws

                # Extraterrestrial daily radiation, Ra (W/m2):
                Ra24_flat = (Gsc/np.pi * dr * (ws_angle * np.sin(phi) * np.sin(delta) +
                                np.cos(phi) * np.cos(delta) * np.sin(ws_angle)))

                if str_TRANS == "MERRA":
                    folder_RAW_file_trans = os.path.join(folders_input_RAW, str_TRANS, "Surface_Net_Downward_Shortwave_Flux", day_name)
                    filename_trans = "swgnet_MERRA_W-m-2_daily_%d.%02d.%02d.tif"  %(Date.year, Date.month, Date.day)
                if str_TRANS == "MSGCPP":
                    folder_RAW_file_trans = os.path.join(folders_input_RAW, str_TRANS, "SDS", "15min")
                    filename_trans = "SDS_MSGCPP_W-m-2_15min_%d.%02d.%02d_H{hour}.M{minutes}.tif"  %(Date.year, Date.month, Date.day)

                if str_TRANS == "MSGCPP":
                    import glob
                    os.chdir(folder_RAW_file_trans)
                    files = glob.glob(filename_trans.format(hour = "*", minutes = "*"))
                    i = 0

                    # Open all the 15 minutes files
                    for file in files:
                        file_in = os.path.join(folder_RAW_file_trans, file)
                        destswgone = gdal.Open(file_in)
                        try:
                            swgnet_one = destswgone.GetRasterBand(1).ReadAsArray()
                            swgnet_one[swgnet_one<0] = 0
                            if not "geo_trans" in locals():
                                swgnet = np.ones([destswgone.RasterYSize, destswgone.RasterXSize, len(files)]) * np.nan
                                geo_trans = destswgone.GetGeoTransform()
                                proj_trans = destswgone.GetProjection()
                            swgnet[:,:,i] = swgnet_one
                        except:
                            pass
                        i+=1

                    # Calculate the daily mean
                    swgnet_mean = np.nanmean(swgnet, 2)
                    dest_swgnet_mean = PF.Save_as_MEM(swgnet_mean,geo_trans, proj_trans)
                    destswgnet = PF.reproject_dataset_example(dest_swgnet_mean, template_file, method=2)

                else:
                    destswgnet = PF.reproject_dataset_example(os.path.join(folder_RAW_file_trans, filename_trans), template_file, method=2)

                swgnet = destswgnet.GetRasterBand(1).ReadAsArray()
                trans = swgnet / Ra24_flat
                PF.Save_as_tiff(Trans_file, trans, geo_ex, proj_ex)
                del geo_trans

        except:
            print("No ETLook input dataset for %s" %Date)

    return()

def lapse_rate_temp(tair_file, dem_file):

    import pyWAPOR
    import pyWAPOR.Functions.Processing_Functions as PF

    destT_down = PF.reproject_dataset_example(tair_file, dem_file, 2)
    destDEM_up = PF.reproject_dataset_example(dem_file, tair_file, 4)
    destDEM_down = gdal.Open(dem_file)
    destDEM_up_down = PF.reproject_dataset_example(destDEM_up, dem_file, 2)

    # Open Arrays
    T = destT_down.GetRasterBand(1).ReadAsArray()
    DEM_down = destDEM_down.GetRasterBand(1).ReadAsArray()
    DEM_up_ave = destDEM_up_down.GetRasterBand(1).ReadAsArray()

    # correct wrong values
    DEM_down[DEM_down<=0]=0
    DEM_up_ave[DEM_up_ave<=0]=0

    #
    Tdown = pyWAPOR.ETLook.meteo.disaggregate_air_temperature(T, DEM_down, DEM_up_ave)

    return(Tdown)

def Combine_LST(folders_input_RAW, Startdate, Enddate):

    import pyWAPOR.Functions.Processing_Functions as PF

    Dates = pd.date_range(Startdate, Enddate, freq = "D")

    output_folder_end = os.path.join(folders_input_RAW, "MODIS", "LST")

    if not os.path.exists(output_folder_end):
        os.makedirs(output_folder_end)

    for Date in Dates:

        LST_file = os.path.join(output_folder_end, "LST_MCD11A1_K_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
        Time_file = os.path.join(output_folder_end, "Time_MCD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
        VZA_file = os.path.join(output_folder_end, "VZA_MCD11A1_degrees_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))

        if not (os.path.exists(Time_file) and os.path.exists(LST_file) and os.path.exists(VZA_file)):
            filename_angle_mod = os.path.join(folders_input_RAW, "MODIS", "MOD11", "Angle_MOD11A1_degrees_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_angle_myd = os.path.join(folders_input_RAW, "MODIS", "MYD11", "Angle_MYD11A1_degrees_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_time_mod = os.path.join(folders_input_RAW, "MODIS", "MOD11", "Time_MOD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_time_myd = os.path.join(folders_input_RAW, "MODIS", "MYD11", "Time_MYD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_lst_mod = os.path.join(folders_input_RAW, "MODIS", "MOD11", "LST_MOD11A1_K_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_lst_myd = os.path.join(folders_input_RAW, "MODIS", "MYD11", "LST_MYD11A1_K_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))

            dest_angle_mod = gdal.Open(filename_angle_mod)
            dest_angle_myd = gdal.Open(filename_angle_myd)
            dest_time_mod = gdal.Open(filename_time_mod)
            dest_time_myd = gdal.Open(filename_time_myd)
            dest_lst_mod = gdal.Open(filename_lst_mod)
            dest_lst_myd = gdal.Open(filename_lst_myd)

            Array_angle_mod = dest_angle_mod.GetRasterBand(1).ReadAsArray()
            Array_angle_myd = dest_angle_myd.GetRasterBand(1).ReadAsArray()
            Array_time_mod = dest_time_mod.GetRasterBand(1).ReadAsArray()
            Array_time_myd = dest_time_myd.GetRasterBand(1).ReadAsArray()
            Array_lst_mod = dest_lst_mod.GetRasterBand(1).ReadAsArray()
            Array_lst_myd = dest_lst_myd.GetRasterBand(1).ReadAsArray()

            LST = Array_lst_mod
            Time = Array_time_mod
            VZA = Array_angle_mod
            LST = np.where(np.abs(Array_angle_myd)<np.abs(Array_angle_mod), Array_lst_myd, LST)
            Time = np.where(np.abs(Array_angle_myd)<np.abs(Array_angle_mod), Array_time_myd, Time)
            VZA = np.where(np.abs(Array_angle_myd)<np.abs(Array_angle_mod), Array_angle_myd, VZA)

            proj_ex = dest_angle_mod.GetProjection()
            geo_ex = dest_angle_mod.GetGeoTransform()

            PF.Save_as_tiff(LST_file, LST, geo_ex, proj_ex)
            PF.Save_as_tiff(Time_file, Time, geo_ex, proj_ex)
            PF.Save_as_tiff(VZA_file, VZA, geo_ex, proj_ex)

    return()


# get the last day in the month in date
def _get_last_day_in_month(currentdate):
    next_month = currentdate.replace(day=28) + datetime.timedelta(days=4)
    return next_month - datetime.timedelta(days=next_month.day)


# get dekadal daterange
def _get_dekadal_daterange(currentdate):
    if currentdate.day < 21:
        return [currentdate + datetime.timedelta(days=x) for x in range(10)]
    else:
        delta = _get_last_day_in_month(currentdate) - currentdate
        return [currentdate + datetime.timedelta(days=x) for x in range(delta.days+1)]


# get the dekadal composite date matching the current date
def _get_dekadal_date(currentdate):
    dekadal_start_dates = [datetime.datetime(currentdate.year, currentdate.month, day) for day in [1, 11, 21]]
    dekadal_center_dates = [datetime.datetime(currentdate.year, currentdate.month, day) for day in [5, 15, 25]]
    nearest_arg = np.argmin(np.abs([currentdate - dekadal_center_date for dekadal_center_date in dekadal_center_dates]))
    return dekadal_start_dates[nearest_arg]


# numeric indexing of nd array
def _numeric_nd_indexing(array, idx):
    composite = np.zeros(array.shape[1:])
    for n in range(0, array.shape[0]):
        composite[idx == n] = array[n, ...][idx == n]
    return composite
