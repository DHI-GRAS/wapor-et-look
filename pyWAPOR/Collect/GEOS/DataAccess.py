# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Tue Feb 19 08:47:49 2019
"""

import os
import numpy as np
import pandas as pd
import urllib

def DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, TimeStep, Period, Waitbar):

	# WAPOR modules
    import pyWAPOR.Functions.Processing_Functions as PF

    # Check the latitude and longitude and otherwise set lat or lon on greatest extent
    if latlim[0] < -90 or latlim[1] > 90:
        print('Latitude above 90N or below 90S is not possible. Value set to maximum')
        latlim[0] = np.max(latlim[0], -90)
        latlim[1] = np.min(latlim[1], 90)
    if lonlim[0] < -180 or lonlim[1] > 180:
        print('Longitude must be between 180E and 180W. Now value is set to maximum')
        lonlim[0] = np.max(lonlim[0], -180)
        lonlim[1] = np.min(lonlim[1], 180)

    # Get information of the parameter
    VarInfo = VariablesInfo(TimeStep)
    Parameter = VarInfo.names[Var]
    unit  = VarInfo.units[Var]
    types  = VarInfo.types[Var]

    # Create output folder
    output_folder = os.path.join(Dir, "GEOS", Parameter, TimeStep)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Define IDs
    IDx = [np.floor((lonlim[0] + 180)/0.3125), np.ceil((lonlim[1] + 180)/0.3125)]
    IDy = [np.floor((latlim[0] + 90)/0.25), np.ceil((latlim[1] + 90)/0.25)]

    # Create output geo transform
    Xstart = -180 + 0.3125 * IDx[0]
    Ystart = -90 + 0.25 * IDy[1]
    geo_out = tuple([Xstart, 0.3125, 0, Ystart, 0, -0.25])
    proj = "WGS84"

    Dates = pd.date_range(Startdate, Enddate, freq = "D")

    # Create Waitbar
    if Waitbar == 1:
        import pyWAPOR.Functions.WaitbarConsole as WaitbarConsole
        total_amount = len(Dates)
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    for Date in Dates:

        # Define the IDz
        if TimeStep == "three_hourly":
            IDz_start = IDz_end = int(((Date - pd.Timestamp("2017-12-01")).days) * 8) + (Period - 1)
            Hour = int((Period - 1) * 3)
            output_name = os.path.join(output_folder, "%s_GEOS_%s_3-hourly_%d.%02d.%02d_H%02d.M00.tif"%(Var, unit, Date.year, Date.month, Date.day, Hour))

        if TimeStep == "daily":
            IDz_start = int(((Date - pd.Timestamp("2017-12-01")).days) * 8)
            IDz_end = IDz_start + 7
            output_name = os.path.join(output_folder, "%s_GEOS_%s_daily_%d.%02d.%02d.tif"%(Var, unit, Date.year, Date.month, Date.day))


        if not os.path.exists(output_name):

            # define total url
            url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/inst3_2d_asm_Nx."
            url_GEOS = url_start + 'ascii?%s[%s:1:%s][%s:1:%s][%s:1:%s]' %(Var, IDz_start,IDz_end, int(IDy[0]),int(IDy[1]),int(IDx[0]),int(IDx[1]))

            # Reset the begin parameters for downloading
            downloaded = 0
            N = 0

            # if not downloaded try to download file
            while downloaded == 0:
                try:

                    # download data (first save as text file)
                    pathtext = os.path.join(output_folder,'temp%s.txt' %str(IDz_start))

                    # Download the data
                    urllib.request.urlretrieve(url_GEOS, filename=pathtext)

                    # Reshape data
                    datashape = [int(IDy[1] - IDy[0] + 1), int(IDx[1] - IDx[0] + 1)]
                    data_start = np.genfromtxt(pathtext,dtype = float,skip_header = 1,skip_footer = 6,delimiter = ',')
                    data_list = np.asarray(data_start[:,1:])
                    if TimeStep == "daily":
                        data_end = np.resize(data_list,(8, datashape[0], datashape[1]))
                    if TimeStep == "three_hourly":
                        data_end = np.resize(data_list,(datashape[0], datashape[1]))
                    os.remove(pathtext)

                    # Set no data value
                    data_end[data_end>1000000] = -9999

                    if TimeStep == "daily":
                        if Var == "t2m":
                            data_end_max = np.nanmax(data_end, 0)
                            data_end_min = np.nanmin(data_end, 0)
                        if types == "state":
                            data_end = np.nanmean(data_end, 0)
                        else:
                            data_end = np.nansum(data_end, 0)

                    # Add the VarFactor
                    if VarInfo.factors[Var] < 0:
                        data_end[data_end != -9999] = data_end[data_end != -9999] + VarInfo.factors[Var]
                    else:
                        data_end[data_end != -9999] = data_end[data_end != -9999] * VarInfo.factors[Var]
                    data_end[data_end < -9999] = -9999

                    # Download was succesfull
                    downloaded = 1

                    # twist the data
                    data_end = np.flipud(data_end)

                    # Save as tiff file
                    PF.Save_as_tiff(output_name, data_end, geo_out, proj)
                    
                    if TimeStep == "daily" and Var == "t2m":
                        # Add the VarFactor
                        if VarInfo.factors[Var] < 0:
                            data_end_max[data_end != -9999] = data_end_max[data_end != -9999] + VarInfo.factors[Var]
                            data_end_min[data_end != -9999] = data_end_min[data_end != -9999] + VarInfo.factors[Var]
                        else:
                            data_end_max[data_end != -9999] = data_end_max[data_end != -9999] * VarInfo.factors[Var]
                            data_end_min[data_end != -9999] = data_end_min[data_end != -9999] * VarInfo.factors[Var]
                        data_end_max[data_end < -9999] = -9999
                        data_end_min[data_end < -9999] = -9999
                        
                        # twist the data
                        data_end_min = np.flipud(data_end_min)
                        data_end_max = np.flipud(data_end_max)

                        # Save as tiff file
                        output_name = os.path.join(output_folder, "%s_GEOS_%s_daily_max_%d.%02d.%02d.tif"%(Var, unit, Date.year, Date.month, Date.day))
                        PF.Save_as_tiff(output_name, data_end_max, geo_out, proj)
                        output_name = os.path.join(output_folder, "%s_GEOS_%s_daily_min_%d.%02d.%02d.tif"%(Var, unit, Date.year, Date.month, Date.day))
                        PF.Save_as_tiff(output_name, data_end_min, geo_out, proj)

                # If download was not succesfull
                except:

                    # Try another time
                    N = N + 1

                    # Stop trying after 10 times
                    if N == 10:
                        print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                        downloaded = 1

        if Waitbar == 1:
            amount += 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    return()

class VariablesInfo:
    """
    This class contains the information about the GEOS variables
    """
    names = {'t2m': 'Air_Temperature',
             'u2m': 'Eastward_Wind',
             'v2m': 'Northward_Wind',
             'qv2m': 'Specific_Humidity',
             'tqv': 'Total_Precipitable_Water_Vapor',
             'ps': 'Surface_Pressure',
             'slp': 'Sea_Level_Pressure'}

    descriptions = {'t2m': '2m Air Temperature',
             'u2m': '2m Eastward wind',
             'v2m': '2m Northward wind',
             'qv2m': '2m Specific Humidity',
             'tqv': 'Total Precipitable Water Vapor',
             'ps': 'Surface Pressure',
             'slp': 'Sea Level Pressure'}

    factors = {'t2m': 1,
             'u2m': 1,
             'v2m': 1,
             'qv2m': 1,
             'tqv': 1,
             'ps': 0.001,
             'slp':  0.001}

    types = {'t2m': 'state',
             'u2m': 'state',
             'v2m': 'state',
             'qv2m': 'state',
             'tqv': 'state',
             'ps': 'state',
             'slp': 'state'}

    def __init__(self, step):
        if step == 'three_hourly':
            self.units = {'t2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'qv2m': 'kg-kg-1',
             'tqv': 'mm',
             'ps': 'kpa',
             'slp': 'kpa'}

        elif step == 'daily':
            self.units = {'t2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'qv2m': 'kg-kg-1',
             'tqv': 'mm',
             'ps': 'kpa',
             'slp': 'kpa'}

        else:
            raise KeyError("The input time step is not supported")

