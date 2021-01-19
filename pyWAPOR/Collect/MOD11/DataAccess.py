# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MOD11
"""

# import general python modules
import os
import numpy as np
import pandas as pd
from osgeo import gdal
import urllib
from bs4 import BeautifulSoup
import re
import glob
import requests
import sys
if sys.version_info[0] == 3:
    import urllib.parse
if sys.version_info[0] == 2:
    import urlparse
    import urllib2

def DownloadData(Dir, Startdate, Enddate, latlim, lonlim, username, password, Waitbar, hdf_library, remove_hdf):
    """
    This function downloads MOD11 daily LST data

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax] (values must be between -90 and 90)
    lonlim -- [xmin, xmax] (values must be between -180 and 180)
    Waitbar -- 1 (Default) will print a waitbar
    hdf_library -- string, if all the hdf files are already stored on computer
                    define directory to the data here
    remove_hdf -- 1 (Default), if 1 remove all the downloaded hdf files in the end    
    """
	
    # Check start and end date and otherwise set the date to max
    if not Startdate:
        Startdate = pd.Timestamp('2000-02-18')
    if not Enddate:
        Enddate = pd.Timestamp('Now')

    # Make an array of the days of which the LST is taken
    Dates = pd.date_range(Startdate, Enddate, freq = 'D')

    # Create Waitbar
    if Waitbar == 1:
        import pyWAPOR.Functions.WaitbarConsole as WaitbarConsole
        total_amount = len(Dates)
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    # Check the latitude and longitude and otherwise set lat or lon on greatest extent
    if latlim[0] < -90 or latlim[1] > 90:
        print('Latitude above 90N or below 90S is not possible. Value set to maximum')
        latlim[0] = np.max(latlim[0], -90)
        latlim[1] = np.min(latlim[1], 90)
    if lonlim[0] < -180 or lonlim[1] > 180:
        print('Longitude must be between 180E and 180W. Now value is set to maximum')
        lonlim[0] = np.max(lonlim[0], -180)
        lonlim[1] = np.min(lonlim[1], 180)

    # Make directory for the MODIS NDVI data
    Dir = Dir.replace("/", os.sep)
    output_folder = os.path.join(Dir, 'MODIS', 'MOD11')
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Define which MODIS tiles are required
    TilesVertical, TilesHorizontal = Get_tiles_from_txt(latlim, lonlim) 
    
    # Pass variables to parallel function and run
    args = [output_folder, TilesVertical, TilesHorizontal,lonlim, latlim, username, password, hdf_library]
    for Date in Dates:
        RetrieveData(Date, args)
        if Waitbar == 1:
            amount += 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    # Remove all .hdf files
    if remove_hdf == 1:
        try:
            os.chdir(output_folder)
            files = glob.glob("*.hdf")
            for f in files:
                os.remove(os.path.join(output_folder, f))
    
        except:
            pass

    return()

def RetrieveData(Date, args):
    """
    This function retrieves MOD11 LST data for a given date from the
    https://e4ftl01.cr.usgs.gov/ server.

    Keyword arguments:
    Date -- 'yyyy-mm-dd'
    args -- A list of parameters defined in the DownloadData function.
    """
    
	# WAPOR modules
    import pyWAPOR.Functions.Processing_Functions as PF
    
    # Argument
    [output_folder, TilesVertical, TilesHorizontal,lonlim, latlim, username, password, hdf_library] = args
    
    # Define output names
    LSTfileName = os.path.join(output_folder, 'LST_MOD11A1_K_daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')
    TimefileName = os.path.join(output_folder, 'Time_MOD11A1_hour_daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')    
    OnsangfileName = os.path.join(output_folder, 'Angle_MOD11A1_degrees_daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')    

    if not (os.path.exists(LSTfileName) and os.path.exists(TimefileName) and os.path.exists(OnsangfileName)):

        # Collect the data from the MODIS webpage and returns the data and lat and long in meters of those tiles
        try:
            Collect_data(TilesHorizontal, TilesVertical, Date, username, password, output_folder, hdf_library)
        
            # Define the output name of the collect data function
            name_collect = os.path.join(output_folder, 'Merged.tif')
            name_collect_time = os.path.join(output_folder, 'Merged_Time.tif')
            name_collect_obsang = os.path.join(output_folder, 'Merged_Obsang.tif')
            
            # Reproject the MODIS product to epsg_to
            epsg_to ='4326'
            name_reprojected = PF.reproject_MODIS(name_collect, epsg_to)
            name_reprojected_time = PF.reproject_MODIS(name_collect_time, epsg_to)
            name_reprojected_obsang = PF.reproject_MODIS(name_collect_obsang, epsg_to)
            
            # Clip the data to the users extend
            data, geo = PF.clip_data(name_reprojected, latlim, lonlim)
            data_time, geo = PF.clip_data(name_reprojected_time, latlim, lonlim)
            data_obsang, geo = PF.clip_data(name_reprojected_obsang, latlim, lonlim)
            
            # remove wrong values
            data[data==0.] = -9999
            data_time[data_time==25.5] = -9999
           
            # Save results as Gtiff 
            PF.Save_as_tiff(name=LSTfileName, data=data, geo=geo, projection='WGS84')
            PF.Save_as_tiff(name=TimefileName, data=data_time, geo=geo, projection='WGS84')
            PF.Save_as_tiff(name=OnsangfileName, data=data_obsang, geo=geo, projection='WGS84')
             
            # remove the side products
            os.remove(name_collect)
            os.remove(name_reprojected)
            os.remove(name_collect_time)
            os.remove(name_reprojected_time) 
            os.remove(name_collect_obsang)
            os.remove(name_reprojected_obsang)             
        except:
            print("Was not able to download the file")
    
    return True

def Collect_data(TilesHorizontal, TilesVertical, Date, username, password, output_folder, hdf_library):
    '''
    This function downloads all the needed MODIS tiles from https://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.006/ as a hdf file.

    Keywords arguments:
    TilesHorizontal -- [TileMin,TileMax] max and min horizontal tile number
    TilesVertical -- [TileMin,TileMax] max and min vertical tile number
    Date -- 'yyyy-mm-dd'
    output_folder -- 'C:/file/to/path/'
    '''

    # Make a new tile for the data
    sizeX = int((TilesHorizontal[1] - TilesHorizontal[0] + 1) * 1200)
    sizeY = int((TilesVertical[1] - TilesVertical[0] + 1) * 1200)
    DataTot = np.zeros((sizeY, sizeX))
    DataTot_Time = np.zeros((sizeY, sizeX))
    DataTot_ObsAng = np.zeros((sizeY, sizeX))
    
    # Create the Lat and Long of the MODIS tile in meters
    for Vertical in range(int(TilesVertical[0]), int(TilesVertical[1])+1):
        
        Distance = 4*231.65635826395834 # resolution of a MODIS pixel in meter
        countY=(TilesVertical[1] - TilesVertical[0] + 1) - (Vertical - TilesVertical[0])

        for Horizontal in range(int(TilesHorizontal[0]), int(TilesHorizontal[1]) + 1):
            countX=Horizontal - TilesHorizontal[0] + 1

            # Create the URL to the LST MODIS data
            url = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.006/' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '/'

		    # Reset the begin parameters for downloading
            downloaded = 0
            N=0

			# Check the library given by user if the file is already there
            if hdf_library is not None:
                os.chdir(hdf_library)
                hdf_name = glob.glob("MOD11A1.A%s%03s.h%02dv%02d.*" %(Date.strftime('%Y'), Date.strftime('%j'), Horizontal, Vertical))

                if len(hdf_name) == 1:
                    hdf_file = os.path.join(hdf_library, hdf_name[0])

                    if os.path.exists(hdf_file):
                        downloaded = 1
                        file_name = hdf_file
                        
            # If the file is not existing, download the data
            if not downloaded == 1:

                # Get files on FTP server
                if sys.version_info[0] == 3:
                    f = urllib.request.urlopen(url)

                if sys.version_info[0] == 2:
                    f = urllib2.urlopen(url)

                # Sum all the files on the server
                soup = BeautifulSoup(f, "lxml")
                for i in soup.findAll('a', attrs = {'href': re.compile('(?i)(hdf)$')}):

                    # Find the file with the wanted tile number
                    Vfile=str(i)[30:32]
                    Hfile=str(i)[27:29]
                    
                    if int(Vfile) is int(Vertical) and int(Hfile) is int(Horizontal):

                        # Define the whole url name
                        if sys.version_info[0] == 3:
                            full_url = urllib.parse.urljoin(url, i['href'])

                        if sys.version_info[0] == 2:
                            full_url = urlparse.urljoin(url, i['href'])

                        # if not downloaded try to download file
                        while downloaded == 0:

                            try: # open http and download whole .hdf
                                nameDownload = full_url
                                file_name = os.path.join(output_folder,nameDownload.split('/')[-1])
                                if os.path.isfile(file_name):
                                    print("file ", file_name, " already exists")
                                    downloaded = 1
                                else:
                                    x = requests.get(nameDownload, allow_redirects = False)
                                    try:
                                        y = requests.get(x.headers['location'], auth = (username, password))
                                    except:
                                        from requests.packages.urllib3.exceptions import InsecureRequestWarning
                                        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

                                        y = requests.get(x.headers['location'], auth = (username, password), verify = False)
                                        
                                    z = open(file_name, 'wb')
                                    z.write(y.content)
                                    z.close()
                                    statinfo = os.stat(file_name)
                                    # Say that download was succesfull
                                    if int(statinfo.st_size) > 10000:
                                         downloaded = 1

                            # If download was not succesfull
                            except:

                                # Try another time
                                N = N + 1

    				    # Stop trying after 10 times
                        if N == 10:
                            print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                            downloaded = 1
                            
            try:
                # Open .hdf only band with LST and collect all tiles to one array
                dataset = gdal.Open(file_name)
                sdsdict = dataset.GetMetadata('SUBDATASETS')
                sdslist = [sdsdict[k] for k in sdsdict.keys() if (('SUBDATASET_1_NAME') in k or ('SUBDATASET_3_NAME') in k or ('SUBDATASET_4_NAME') in k)]
                sds = []
                sds_time = []
                sds_obsang = []
                                
                full_layer = [i for i in sdslist if 'LST_Day_1km' in i]
                full_layer_time = [i for i in sdslist if 'Day_view_time' in i]
                full_layer_obsang = [i for i in sdslist if 'Day_view_angl' in i]
                idx = sdslist.index(full_layer[0])
                idx_time = sdslist.index(full_layer_time[0])
                idx_obsang = sdslist.index(full_layer_obsang[0])
                sds.append(gdal.Open(sdslist[idx]))
                sds_time.append(gdal.Open(sdslist[idx_time]))     
                sds_obsang.append(gdal.Open(sdslist[idx_obsang]))                  
                if Horizontal == TilesHorizontal[0] and Vertical == TilesVertical[0]:
                    geo_t = sds[idx].GetGeoTransform()

                    # get the projection value
                    proj = sds[idx].GetProjection()

                data = sds[0].ReadAsArray()
                data_time = sds_time[0].ReadAsArray()
                data_obsang = sds_obsang[0].ReadAsArray()

                
                countYdata = (TilesVertical[1] - TilesVertical[0] + 2) - countY
                DataTot[int((countYdata - 1) * 1200):int(countYdata * 1200), int((countX - 1) * 1200):int(countX * 1200)]=data * 0.02
                DataTot_Time[int((countYdata - 1) * 1200):int(countYdata * 1200), int((countX - 1) * 1200):int(countX * 1200)]=data_time * 0.1
                DataTot_ObsAng[int((countYdata - 1) * 1200):int(countYdata * 1200), int((countX - 1) * 1200):int(countX * 1200)]= data_obsang
                del data, data_time, data_obsang

            # if the tile not exists or cannot be opened, create a nan array with the right projection
            except:
                if Horizontal==TilesHorizontal[0] and Vertical==TilesVertical[0]:
                     x1 = (TilesHorizontal[0] - 19) * 1200 * Distance
                     x4 = (TilesVertical[0] - 9) * 1200 * -1 * Distance
                     geo = 	[x1, Distance, 0.0, x4, 0.0, -Distance]
                     geo_t=tuple(geo)

                proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
                data=np.ones((1200,1200)) * (-9999/0.02)
                data_time=np.ones((1200,1200)) * (-9999/0.1)
                data_obsang=np.ones((1200,1200)) * 255
                countYdata=(TilesVertical[1] - TilesVertical[0] + 2) - countY
                DataTot[(countYdata - 1) * 1200:countYdata * 1200,(countX - 1) * 1200:countX * 4800] = data * 0.02
                DataTot_Time[(countYdata - 1) * 1200:countYdata * 1200,(countX - 1) * 1200:countX * 4800] = data_time * 0.1
                DataTot_ObsAng[(countYdata - 1) * 1200:countYdata * 1200,(countX - 1) * 1200:countX * 4800] = data_obsang
                
    # Set limits    
    DataTot_Time[np.logical_and(DataTot_Time < 0., DataTot_Time>24.)] = -9999
    DataTot[DataTot < 1] = -9999   
    DataTot_ObsAng[DataTot_ObsAng<255] = DataTot_ObsAng[DataTot_ObsAng<255] - 65
    DataTot_ObsAng[DataTot_ObsAng==255] = -9999
    DataTot_ObsAng[DataTot_ObsAng>255] = DataTot_ObsAng[DataTot_ObsAng>255] - 256
    
    # Make geotiff file of the LST data
    name_out_LST = os.path.join(output_folder, 'Merged.tif')
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(name_out_LST, DataTot.shape[1], DataTot.shape[0], 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
    try:
         dst_ds.SetProjection(proj)
    except:
        proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
        x1 = (TilesHorizontal[0] - 18) * 1200 * Distance
        x4 = (TilesVertical[0] - 9) * 1200 * -1 * Distance
        geo = [x1, Distance, 0.0, x4, 0.0, -Distance]
        geo_t = tuple(geo)
        dst_ds.SetProjection(proj)

    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo_t)
    dst_ds.GetRasterBand(1).WriteArray(DataTot)
    dst_ds = None
    sds = None

    # Make geotiff file of the time data
    name_out_time = os.path.join(output_folder, 'Merged_Time.tif')
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(name_out_time, DataTot_Time.shape[1], DataTot_Time.shape[0], 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
    try:
         dst_ds.SetProjection(proj)
    except:
        proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
        x1 = (TilesHorizontal[0] - 18) * 1200 * Distance
        x4 = (TilesVertical[0] - 9) * 1200 * -1 * Distance
        geo = [x1, Distance, 0.0, x4, 0.0, -Distance]
        geo_t = tuple(geo)
        dst_ds.SetProjection(proj)

    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo_t)
    dst_ds.GetRasterBand(1).WriteArray(DataTot_Time)
    dst_ds = None
    sds = None    

    # Make geotiff file of the time data
    name_out_obsang = os.path.join(output_folder, 'Merged_Obsang.tif')
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(name_out_obsang, DataTot_ObsAng.shape[1], DataTot_ObsAng.shape[0], 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
    try:
         dst_ds.SetProjection(proj)
    except:
        proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
        x1 = (TilesHorizontal[0] - 18) * 1200 * Distance
        x4 = (TilesVertical[0] - 9) * 1200 * -1 * Distance
        geo = [x1, Distance, 0.0, x4, 0.0, -Distance]
        geo_t = tuple(geo)
        dst_ds.SetProjection(proj)

    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo_t)
    dst_ds.GetRasterBand(1).WriteArray(DataTot_ObsAng)
    dst_ds = None
    sds = None    
    
    return()
    
def Get_tiles_from_txt(latlim, lonlim):
    '''
    Opens a text file storing the extents of every MODIS tile

    Keywords arguments:
    tiletext2 -- 'C:/file/to/path/' to path of the txt file with all the MODIS tiles extents
    lonlim1 -- [ymin, ymax] (longitude limits of the chunk or whole image)
    latlim1 -- [ymin, ymax] (latitude limits of the chunk or whole image)
    '''

    # Get the list (txt file from the internet) which includes the lat and lon information of the MODIS tiles
    path = os.path.abspath(os.path.dirname(__file__))
    file_nametext = os.path.join(path, "sn_gring_10deg.txt")

    # Open text file with tiles which is downloaded before
    tiletext=np.genfromtxt(file_nametext,skip_header=7,skip_footer=1,usecols=(0,1,2,3,4,5,6,7,8,9))
    tiletext2=tiletext[tiletext[:,2]>=-900,:]

    # This function converts the values in the text file into horizontal and vertical number of the tiles which must be downloaded to cover the extent defined by the user
    TilesVertical, TilesHorizontal = Tiles_to_download(tiletext2=tiletext2,lonlim1=lonlim,latlim1=latlim)

    return(TilesVertical, TilesHorizontal)

def Tiles_to_download(tiletext2,lonlim1,latlim1):
    '''
    Defines the MODIS tiles that must be downloaded in order to cover the latitude and longitude limits

    Keywords arguments:
    tiletext2 -- 'C:/file/to/path/' to path of the txt file with all the MODIS tiles extents
    lonlim1 -- [ymin, ymax] (longitude limits of the chunk or whole image)
    latlim1 -- [ymin, ymax] (latitude limits of the chunk or whole image)
    '''
    # calculate min and max longitude and latitude
    # lat down    lat up      lon left     lon right
    tiletextExtremes = np.empty([len(tiletext2),6])
    tiletextExtremes[:,0] = tiletext2[:,0]
    tiletextExtremes[:,1] = tiletext2[:,1]
    tiletextExtremes[:,2] = np.minimum(tiletext2[:,3], tiletext2[:,9])
    tiletextExtremes[:,3] = np.maximum(tiletext2[:,5], tiletext2[:,7])
    tiletextExtremes[:,4] = np.minimum(tiletext2[:,2], tiletext2[:,4])
    tiletextExtremes[:,5] = np.maximum(tiletext2[:,6], tiletext2[:,8])

    # Define the upper left tile
    latlimtiles1UL = tiletextExtremes[np.logical_and(tiletextExtremes[:,2] <= latlim1[1], tiletextExtremes[:,3] >= latlim1[1])]
    latlimtilesUL = latlimtiles1UL[np.logical_and(latlimtiles1UL[:,4] <= lonlim1[0], latlimtiles1UL[:,5] >= lonlim1[0])]

    # Define the lower right tile
    latlimtiles1LR = tiletextExtremes[np.logical_and(tiletextExtremes[:,2] <= latlim1[0], tiletextExtremes[:,3] >= latlim1[0])]
    latlimtilesLR = latlimtiles1LR[np.logical_and(latlimtiles1LR[:,4]<=lonlim1[1],latlimtiles1LR[:,5]>=lonlim1[1])]

    # Define the total tile
    TotalTiles = np.vstack([latlimtilesUL, latlimtilesLR])

    # Find the minimum horizontal and vertical tile value and the maximum horizontal and vertical tile value
    TilesVertical = [TotalTiles[:,0].min(), TotalTiles[:,0].max()]
    TilesHorizontal = [TotalTiles[:,1].min(), TotalTiles[:,1].max()]
    
    return(TilesVertical, TilesHorizontal)    
