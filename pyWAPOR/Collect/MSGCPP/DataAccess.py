# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Tue Feb 19 07:10:44 2019
"""
import os
import urllib
import datetime
import pandas as pd
import numpy as np

def DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Time = '', Waitbar = 1):
    
    # Check the latitude and longitude and otherwise set lat or lon on greatest extent
    if latlim[0] < -82 or latlim[1] > 82:
        print('Latitude above 82N or below 82S is not possible. Value set to maximum')
        latlim[0] = np.max(latlim[0], -82)
        latlim[1] = np.min(latlim[1], 82)
    if lonlim[0] < -80 or lonlim[1] > 80:
        print('Longitude must be between 80E and 80W. Now value is set to maximum')
        lonlim[0] = np.max(lonlim[0], -80)
        lonlim[1] = np.min(lonlim[1], 80)
    
    output_folder = os.path.join(Dir, "MSGCPP", "SDS", "15min")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    if Time == '':
        Dates = pd.date_range(Startdate, Enddate, freq = "15min")
    else:
        Dates = pd.date_range(Startdate, Enddate, freq = "D")
    
    # Create Waitbar
    if Waitbar == 1:
        import pyWAPOR.Functions.WaitbarConsole as WaitbarConsole
        total_amount = len(Dates)
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
    # Loop over dates
    for Date in Dates:
        
        if Time is not '':
            Hour = int(Time.split(':')[0])
            Minute = int(Time.split(':')[1])        
            Date = datetime.datetime(Date.year, Date.month, Date.day, Hour, Minute)
        
        filename_out = os.path.join(output_folder, "SDS_MSGCPP_W-m-2_15min_%d.%02d.%02d_H%02d.M%02d.tif" %(Date.year, Date.month, Date.day, Date.hour, Date.minute))

        if not os.path.exists(filename_out):
        
            # define url
            url = r"http://msgcpp-ogc-archive.knmi.nl/msgar.cgi?&service=wcs&version=1.0.0&request=getcoverage&coverage=surface_downwelling_shortwave_flux_in_air&FORMAT=GeoTIFF&CRS=EPSG%%3A4326&BBOX=%s,%s,%s,%s&RESX=0.04310344827586207&RESY=0.04418103448275862&time=%d-%02d-%02dT%02d%%3A%02d%%3A00Z" %(lonlim[0],latlim[0], lonlim[1], latlim[1], Date.year, Date.month, Date.day, Date.hour, Date.minute)
            urllib.request.urlretrieve(url, filename=filename_out)
        
        if Waitbar == 1:
            amount += 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
        