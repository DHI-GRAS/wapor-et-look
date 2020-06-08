# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Tue Feb 19 08:47:49 2019
"""

import os
import numpy as np
import pandas as pd
import datetime
import urllib
import requests
from netCDF4 import Dataset

def DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, TimeStep, Period, username, password, Waitbar):

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
    
    if TimeStep == "yearly":
        Parameter = "Temperature_Amplitude"
    
    # Create output folder
    output_folder = os.path.join(Dir, "MERRA", Parameter, TimeStep) 
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if TimeStep.split("_")[-1] == "MERRA2":
        corrx = 0.625 * 0.5
        corry = 0.5 * 0.5
    else:
        corrx = 0
        corry = 0
    
    # Define IDs
    IDx = [np.floor((lonlim[0] + corrx + 180)/0.625), np.ceil((lonlim[1] + corrx + 180)/0.625)]
    IDy = [np.floor((latlim[0] + corry + 90)/0.5), np.ceil((latlim[1] + corry + 90)/0.5)]
    
    # Create output geo transform
    Xstart = -180 + 0.625 * IDx[0] - corrx
    Ystart = -90 + 0.5 * IDy[1] - corry

        
    geo_out = tuple([Xstart, 0.625, 0, Ystart, 0, -0.5])
    proj = "WGS84"
    
    if TimeStep == "yearly":
        Dates = pd.date_range(Startdate, Enddate, freq = "AS")
    else:
        Dates = pd.date_range(Startdate, Enddate, freq = "D")
        
    # Create Waitbar
    if Waitbar == 1:
        import pyWAPOR.Functions.WaitbarConsole as WaitbarConsole
        total_amount = len(Dates)
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
    for Date in Dates:
        
        # Define the IDz
        if TimeStep == "hourly_MERRA2":
            Hour = int((Period - 1) * 1)
            output_name = os.path.join(output_folder, "%s_MERRA_%s_hourly_%d.%02d.%02d_H%02d.M00.tif"%(Var, unit, Date.year, Date.month, Date.day, Hour))
            output_folder_temp = os.path.join(Dir, "MERRA", "Temp")
            if not os.path.exists(output_folder_temp):
                os.makedirs(output_folder_temp)
            year = Date.year
            month = Date.month
            day = Date.day
            
        if TimeStep == "daily_MERRA2":
            output_name = os.path.join(output_folder, "%s_MERRA_%s_daily_%d.%02d.%02d.tif" %(Var, unit, Date.year, Date.month, Date.day))
            output_folder_temp = os.path.join(Dir, "MERRA", "Temp")
            if not os.path.exists(output_folder_temp):
                os.makedirs(output_folder_temp)            
            year = Date.year
            month = Date.month
            day = Date.day
    
        if TimeStep == "three_hourly":
            IDz_start = IDz_end = int(((Date - pd.Timestamp("2002-07-01")).days) * 8) + (Period - 1)
            Hour = int((Period - 1) * 3)
            output_name = os.path.join(output_folder, "%s_MERRA_%s_3-hourly_%d.%02d.%02d_H%02d.M00.tif"%(Var, unit, Date.year, Date.month, Date.day, Hour))
    
        if TimeStep == "daily":
            IDz_start = int(((Date - pd.Timestamp("2002-07-01")).days) * 8) 
            IDz_end = IDz_start + 7
            output_name = os.path.join(output_folder, "%s_MERRA_%s_daily_%d.%02d.%02d.tif"%(Var, unit, Date.year, Date.month, Date.day))

        if TimeStep == "yearly":
            
             # Change date for now, there is no OpenDAP system for those years.
            if Date >= datetime.datetime(2015,1,1):
                Data_data = datetime.datetime(2014,1,1)    
            else:
                Data_data = Date
                
            IDz_start = (Data_data.year - pd.Timestamp("2002-07-01").year) * 12 + Data_data.month - pd.Timestamp("2002-07-01").month 
            IDz_end = IDz_start + 11
            output_name = os.path.join(output_folder, "Tamp_MERRA_%s_yearly_%d.%02d.%02d.tif"%(unit, Date.year, Date.month, Date.day))
     
        if not os.path.exists(output_name):
            
            if (TimeStep == "hourly_MERRA2" or TimeStep == "daily_MERRA2"):
                
                if Date < datetime.datetime(1992,1,1):
                    number = 1
                elif (Date >= datetime.datetime(1992,1,1) and Date < datetime.datetime(2001,1,1)):
                    number = 2
                elif (Date >= datetime.datetime(2001,1,1) and Date < datetime.datetime(2011,1,1)):
                    number = 3
                else:
                    number = 4
                               
                if Var == "swgnet":
                    url_MERRA = r"https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXRAD.5.12.4/%d/%02d/MERRA2_%s00.tavg1_2d_rad_Nx.%d%02d%02d.nc4" %(year, month, number, year, month, day)
                    
                else:    
                    url_MERRA = r"https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/%d/%02d/MERRA2_%s00.inst1_2d_asm_Nx.%d%02d%02d.nc4" %(year, month, number, year, month, day)
            
            if (TimeStep == "three_hourly" or TimeStep == "daily"):                
                
                # define total url
                if (Var == "ps" or Var == "slp"):
                    url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/MERRAero/hourly/inst3hr_3d_asm_Nv."
                else:
                    url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/MERRAero/hourly/tavg3hr_2d_asm_Nx."
                
                url_MERRA = url_start + 'ascii?%s[%s:1:%s][%s:1:%s][%s:1:%s]' %(Var, IDz_start,IDz_end, int(IDy[0]),int(IDy[1]),int(IDx[0]),int(IDx[1]))
            
            if TimeStep == "yearly":
                url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/MERRAero/hourly/tavg3hr_2d_asm_Nx."
                url_MERRA = url_start + 'ascii?%s[%s:1:%s][%s:1:%s][%s:1:%s]' %(Var, IDz_start,IDz_end, int(IDy[0]),int(IDy[1]),int(IDx[0]),int(IDx[1]))
                url_MERRA = url_MERRA.replace("hourly", "monthly")
                

                    
            # Reset the begin parameters for downloading
            downloaded = 0
            N = 0

            # if not downloaded try to download file
            while downloaded == 0:
                try:
                    if (TimeStep == "hourly_MERRA2" or TimeStep == "daily_MERRA2"):
                                          
                        # Define the output name that is downloaded
                        file_name = os.path.join(output_folder_temp, url_MERRA.split("/")[-1])
                        if not os.path.exists(file_name):
                        
                            # make contact with server
                            x = requests.get(url_MERRA, allow_redirects = False)
                            try:
                                y = requests.get(x.headers['location'], auth = (username, password))
                            except:
                                from requests.packages.urllib3.exceptions import InsecureRequestWarning
                                requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                                y = requests.get(x.headers['location'], auth = (username, password), verify = False)
    
                            
                            # Write the download in the output directory
                            z = open(file_name, 'wb')
                            z.write(y.content)
                            z.close()
                            statinfo = os.stat(file_name)
                            # Say that download was succesfull
                            if int(statinfo.st_size) > 1000:
                                 downloaded = 1  
                        else:
                            downloaded = 1                             
                        
                        data_end = Get_NC_data_end(file_name,Var, TimeStep, Period, IDy, IDx, VarInfo)
                        #os.remove(file_name)
                                  
                    else:    
                            
                        # download data (first save as text file)
                        pathtext = os.path.join(output_folder,'temp%s.txt' %str(IDz_start))
                        
                        # Download the data
                        urllib.request.urlretrieve(url_MERRA, filename=pathtext)
        
                        # Reshape data
                        datashape = [int(IDy[1] - IDy[0] + 1), int(IDx[1] - IDx[0] + 1)]
                        data_start = np.genfromtxt(pathtext,dtype = float,skip_header = 1,skip_footer = 6,delimiter = ',')
                        data_list = np.asarray(data_start[:,1:])
                        if TimeStep == "yearly":
                            data_end = np.resize(data_list,(12, datashape[0], datashape[1]))
                        if TimeStep == "daily":
                            data_end = np.resize(data_list,(8, datashape[0], datashape[1]))
                        if TimeStep == "three_hourly":
                            data_end = np.resize(data_list,(datashape[0], datashape[1]))
                        os.remove(pathtext)
        
                        # Set no data value
                        data_end[data_end>1000000] = -9999
                        
                        if TimeStep == "daily":
                            if types == "state":
                                data_end = np.nanmean(data_end, 0)
                            else:
                                data_end = np.nansum(data_end, 0)   
                                
                        if TimeStep == "yearly":
                            data_min = np.nanmin(data_end, 0)
                            data_max = np.nanmax(data_end, 0) 
                            data_end = data_max - data_min
                                
             
                        # Download was succesfull
                        downloaded = 1

                    # Add the VarFactor
                    if VarInfo.factors[Var] < 0:
                        data_end[data_end != -9999] = data_end[data_end != -9999] + VarInfo.factors[Var]
                    else:
                        data_end[data_end != -9999] = data_end[data_end != -9999] * VarInfo.factors[Var]
                    data_end[data_end < -9999] = -9999
        
                    # twist the data                
                    data_end = np.flipud(data_end)
                
                    # Save as tiff file
                    PF.Save_as_tiff(output_name, data_end, geo_out, proj)
            
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


def Get_NC_data_end(file_name,Var, TimeStep, Period, IDy, IDx, VarInfo):
                        
    dict_para = {'t2m': 'T2M',
         'u2m': 'U2M',
         'v2m': 'V2M',
         'q2m': 'QV2M',
         'tpw': 'TQV',
         'ps': 'PS',
         'slp': 'SLP',
         'swgnet': 'SWGDN'}  
      
    types  = VarInfo.types[Var]
    if TimeStep == "hourly_MERRA2":
        data_end = Dataset(file_name)["%s" %dict_para[Var]][int(Period-1),int(IDy[0]):int(IDy[1]),int(IDx[0]):int(IDx[1])]

    else:
        data = Dataset(file_name)["%s" %dict_para[Var]][:,int(IDy[0]):int(IDy[1]),int(IDx[0]):int(IDx[1])]
        if types == "state":
            data_end = np.nanmean(data, 0)
        else:
            data_end = np.nansum(data, 0)   
            
    return(data_end)

class VariablesInfo:
    """
    This class contains the information about the GLDAS variables
    """
    names = {'t2m': 'Air_Temperature',
             'u2m': 'Eastward_Wind',
             'v2m': 'Northward_Wind',
             'q2m': 'Specific_Humidity',
             'tpw': 'Total_Precipitable_Water_Vapor',
             'ps': 'Surface_Pressure',
             'slp': 'Sea_Level_Pressure',
             'swgnet': 'Surface_Net_Downward_Shortwave_Flux'
             }
    
    descriptions = {'t2m': '2m Air Temperature',
             'u2m': '2m Eastward wind',
             'v2m': '2m Northward wind',
             'q2m': '2m Specific Humidity',
             'tpw': 'Total Precipitable Water Vapor',
             'ps': 'Surface Pressure',
             'slp': 'Sea Level Pressure',
             'swgnet': 'Surface Net Downward Shortwave Flux'
             }
    
    factors = {'t2m': 1,
             'u2m': 1,
             'v2m': 1,
             'q2m': 1,
             'tpw': 1,
             'ps': 0.001,
             'slp':  0.001,
             'swgnet':  1
             }
    
    types = {'t2m': 'state',
             'u2m': 'state',
             'v2m': 'state',
             'q2m': 'state',
             'tpw': 'state',
             'ps': 'state',
             'slp': 'state',
             'swgnet': 'state'
             }

    def __init__(self, step):
        if step == 'three_hourly':
            self.units = {'t2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'q2m': 'kg-kg-1',
             'tpw': 'mm',
             'ps': 'kpa',
             'slp': 'kpa',
             'swgnet': 'W-m-2'
             }
            
        elif step == 'hourly_MERRA2':
            self.units = {'t2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'q2m': 'kg-kg-1',
             'tpw': 'mm',
             'ps': 'kpa',
             'slp': 'kpa',
             'swgnet': 'W-m-2'
             }    
            
        elif (step == 'daily' or step == 'daily_MERRA2'):
            self.units = {'t2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'q2m': 'kg-kg-1',
             'tpw': 'mm',
             'ps': 'kpa',
             'slp': 'kpa',
             'swgnet': 'W-m-2'
             }
            
        elif step == 'yearly':
            self.units = {'t2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'q2m': 'kg-kg-1',
             'tpw': 'mm',
             'ps': 'kpa',
             'slp': 'kpa',
             'swgnet': 'W-m-2'
             }
        else:
            raise KeyError("The input time step is not supported")

