# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Mon Feb 18 20:50:36 2019
"""

import os
import urllib
import gdal
import shutil
import numpy as np

def DownloadData(output_folder, latlim, lonlim):
    
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
        
    # Check output path and create if not exists
    output_folder_end = os.path.join(output_folder, "GlobCover", "Landuse")
    if not os.path.exists(output_folder_end):
        os.makedirs(output_folder_end)

    # Define end output
    filename_out_tiff = os.path.join(output_folder_end, "LC_GLOBCOVER_V2.3.tif")   

    if not os.path.exists(filename_out_tiff):

        # Define url where to download the globcover data
        url = r"http://due.esrin.esa.int/files/Globcover2009_V2.3_Global_.zip"
        
        # Create temp folder
        output_folder_temp = os.path.join(output_folder, "GlobCover", "Landuse", "Temp")
        if not os.path.exists(output_folder_temp):
            os.makedirs(output_folder_temp)
        
        # define output layer
        filename_out = os.path.join(output_folder_temp, "Globcover2009_V2.3_Global_.zip")
        
        # Download the data
        urllib.request.urlretrieve(url, filename=filename_out)
        
        # Extract data
        PF.Extract_Data(filename_out, output_folder_temp)
        
        # Define extracted tiff file
        globcover_filename = os.path.join(output_folder_temp, "GLOBCOVER_L4_200901_200912_V2.3.tif")
        
        # Open extract file
        dest = gdal.Open(globcover_filename)
        Array = dest.GetRasterBand(1).ReadAsArray()
        
        # Get information of geotransform and projection
        geo = dest.GetGeoTransform()
        proj = "WGS84"
        
        # define the spatial ids
        Xid = [np.floor((-geo[0] + lonlim[0])/geo[1]), np.ceil((-geo[0] + lonlim[1])/geo[1])]
        Yid = [np.floor((geo[3] - latlim[1])/-geo[5]), np.ceil((geo[3] - latlim[0])/-geo[5])]
        
        # Define the geotransform
        Xstart = geo[0] + Xid[0] * geo[1]
        Ystart = geo[3] + Yid[0] * geo[5]
        geo_out = tuple([Xstart, geo[1], 0, Ystart, 0, geo[5]])
        
        # Clip data out
        Array_end = Array[int(Yid[0]):int(Yid[1]), int(Xid[0]):int(Xid[1])]

        # Save data as tiff
        PF.Save_as_tiff(filename_out_tiff, Array_end, geo_out, proj)
        dest = None
        
        # remove temporary folder
        shutil.rmtree(output_folder_temp)
        
        