# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/SRTM
"""

# General modules
import numpy as np
import os
import urllib
import shutil
import gdal
import sys

def DownloadData(output_folder, latlim, lonlim):
    """
    This function downloads DEM data from SRTM

    Keyword arguments:
    output_folder -- directory of the result
	latlim -- [ymin, ymax] (values must be between -60 and 60)
    lonlim -- [xmin, xmax] (values must be between -180 and 180)
    """
    
    # WA+ modules
    import pyWAPOR.Functions.Processing_Functions as PF
    
   # Check the latitude and longitude and otherwise set lat or lon on greatest extent
    if latlim[0] < -60 or latlim[1] > 60:
        print('Latitude above 60N or below 60S is not possible. Value set to maximum')
        latlim[0] = np.max(latlim[0], -60)
        latlim[1] = np.min(latlim[1], 60)
    if lonlim[0] < -180 or lonlim[1] > 180:
        print('Longitude must be between 180E and 180W. Now value is set to maximum')
        lonlim[0] = np.max(lonlim[0], -180)
        lonlim[1] = np.min(lonlim[1], 180)

    # converts the latlim and lonlim into names of the tiles which must be
    # downloaded
    name, rangeLon, rangeLat = Find_Document_Names(latlim, lonlim)

    # Memory for the map x and y shape (starts with zero)
    size_X_tot = 0
    size_Y_tot = 0

    nameResults = []
    
    # Create a temporary folder for processing
    output_folder_trash = os.path.join(output_folder, "Temp")
    if not os.path.exists(output_folder_trash):
        os.makedirs(output_folder_trash)

    # Download, extract, and converts all the files to tiff files
    for nameFile in name:

        try:
            # Download the data from
            # http://earlywarning.usgs.gov/hydrodata/
            output_file, file_name = Download_Data(nameFile, output_folder_trash)

            # extract zip data
            PF.Extract_Data(output_file, output_folder_trash)

            # The input is the file name and in which directory the data must be stored
            file_name_tiff = file_name.replace(".zip", ".tif")
            output_tiff = os.path.join(output_folder_trash, file_name_tiff)

            # convert data from adf to a tiff file
            dest_SRTM = gdal.Open(output_tiff)
            geo_out = dest_SRTM.GetGeoTransform()
            size_X = dest_SRTM.RasterXSize
            size_Y = dest_SRTM.RasterYSize            
   
            
            if (int(size_X) != int(6001) or int(size_Y) != int(6001)):
                data = np.ones((6001, 6001)) * -9999

                # Create the latitude bound
                Vfile = nameFile.split("_")[2][0:2]
                Bound2 = 60 - 5*(int(Vfile)-1)

                # Create the longitude bound
                Hfile = nameFile.split("_")[1]
                Bound1 = -180 + 5 * (int(Hfile)-1)

                Expected_X_min = Bound1
                Expected_Y_max = Bound2

                Xid_start = int(np.round((geo_out[0] - Expected_X_min)/geo_out[1]))
                Xid_end = int(np.round(((geo_out[0] + size_X * geo_out[1]) - Expected_X_min)/geo_out[1]))
                Yid_start = int(np.round((Expected_Y_max - geo_out[3])/(-geo_out[5])))
                Yid_end = int(np.round((Expected_Y_max - (geo_out[3] + (size_Y * geo_out[5])))/(-geo_out[5])))
                
                data_SRTM = dest_SRTM.GetRasterBand(1).ReadAsArray()

                data[Yid_start:Yid_end,Xid_start:Xid_end] = data_SRTM
                if np.max(data)==255:
                    data[data==255] = -9999
                data[data<-9999] = -9999

                geo_in = [Bound1 - 0.5 * 0.00083333333333333, 0.00083333333333333, 0.0, Bound2 + 0.5 * 0.00083333333333333,
                          0.0, -0.0008333333333333333333]

                # save chunk as tiff file
                PF.Save_as_tiff(name=output_tiff, data=data, geo=geo_in,
                             projection="WGS84")
                
                dest_SRTM = None 

        except:

            # If tile not exist create a replacing zero tile (sea tiles)
            file_name_tiff = file_name.replace(".zip", ".tif")
            output_tiff = os.path.join(output_folder_trash, file_name_tiff)
            file_name = nameFile
            data = np.ones((6001, 6001)) * -9999
            data = data.astype(np.float32)

            # Create the latitude bound
            Vfile = nameFile.split("_")[2][0:2]
            Bound2 = 60 - 5*(int(Vfile)-1)

            # Create the longitude bound
            Hfile = nameFile.split("_")[1]
            Bound1 = -180 + 5 * (int(Hfile)-1)

            # Geospatial data for the tile
            geo_in = [Bound1 - 0.5 * 0.00083333333333333, 0.00083333333333333, 0.0, Bound2 + 0.5 * 0.00083333333333333,
                          0.0, -0.0008333333333333333333]

            # save chunk as tiff file
            PF.Save_as_tiff(name=output_tiff, data=data, geo=geo_in,
                         projection="WGS84")


        # clip data
        Data, Geo_data = PF.clip_data(output_tiff, latlim, lonlim)
        size_Y_out = int(np.shape(Data)[0])
        size_X_out = int(np.shape(Data)[1])

        # Total size of the product so far
        size_Y_tot = int(size_Y_tot + size_Y_out)
        size_X_tot = int(size_X_tot + size_X_out)

        if nameFile is name[0]:
            Geo_x_end = Geo_data[0]
            Geo_y_end = Geo_data[3]
        else:
            Geo_x_end = np.min([Geo_x_end,Geo_data[0]])
            Geo_y_end = np.max([Geo_y_end,Geo_data[3]])

        # create name for chunk
        FileNameEnd = "%s_temporary.tif" % (file_name)
        nameForEnd = os.path.join(output_folder_trash, FileNameEnd)
        nameResults.append(str(nameForEnd))

        # save chunk as tiff file
        PF.Save_as_tiff(name=nameForEnd, data=Data, geo=Geo_data,
                      projection="WGS84")


    size_X_end = int(size_X_tot/len(rangeLat)) + 1 #!
    size_Y_end = int(size_Y_tot/len(rangeLon)) + 1 #!

    # Define the georeference of the end matrix
    geo_out = [Geo_x_end, Geo_data[1], 0, Geo_y_end, 0, Geo_data[5]]

    latlim_out = [geo_out[3] + geo_out[5] * size_Y_end, geo_out[3]]
    lonlim_out = [geo_out[0], geo_out[0] + geo_out[1] * size_X_end]


    # merge chunk together resulting in 1 tiff map
    datasetTot = Merge_DEM(latlim_out, lonlim_out, nameResults, size_Y_end,
                                size_X_end)

    datasetTot[datasetTot<-9999] = -9999


    # name of the end result
    output_DEM_name = "DEM_SRTM_m_3s.tif" 
    
    Save_name = os.path.join(output_folder, output_DEM_name)

    # Make geotiff file
    PF.Save_as_tiff(name=Save_name, data=datasetTot, geo=geo_out, projection="WGS84")
    os.chdir(output_folder)

    # Delete the temporary folder
    try:
        shutil.rmtree(output_folder_trash)
    except:
        pass
    
    return()


def Merge_DEM(latlim, lonlim, nameResults, size_Y_tot, size_X_tot):
    """
    This function will merge the tiles

    Keyword arguments:
    latlim -- [ymin, ymax], (values must be between -50 and 50)
    lonlim -- [xmin, xmax], (values must be between -180 and 180)
    nameResults -- ['string'], The directories of the tiles which must be
                   merged
    size_Y_tot -- integer, the width of the merged array
    size_X_tot -- integer, the length of the merged array
    """
    # Define total size of end dataset and create zero array
    datasetTot = np.ones([size_Y_tot, size_X_tot])*-9999.

    # Put all the files in the datasetTot (1 by 1)
    for nameTot in nameResults:
        f = gdal.Open(nameTot)
        dataset = np.array(f.GetRasterBand(1).ReadAsArray())
        dataset = np.flipud(dataset)
        geo_out = f.GetGeoTransform()
        BoundChunk1 = int(round((geo_out[0]-lonlim[0])/geo_out[1]))
        BoundChunk2 = BoundChunk1 + int(dataset.shape[1])
        BoundChunk4 = size_Y_tot + int(round((geo_out[3] -
                                       latlim[1]) / geo_out[1]))
        BoundChunk3 = BoundChunk4 - int(dataset.shape[0])
        datasetTot[BoundChunk3:BoundChunk4, BoundChunk1:BoundChunk2] = dataset
        f = None
    datasetTot = np.flipud(datasetTot)
    return(datasetTot)

def Find_Document_Names(latlim, lonlim):
    """
    This function will translate the latitude and longitude limits into
    the filenames that must be downloaded from the hydroshed webpage

    Keyword Arguments:
    latlim -- [ymin, ymax] (values must be between -60 and 60)
    lonlim -- [xmin, xmax] (values must be between -180 and 180)
    """
    # find tiles that must be downloaded
    startLat = np.floor((60 - latlim[1])/ 5) + 1
    startLon = np.floor((180 + lonlim[0])/ 5) + 1
    endLat = np.ceil((60 - latlim[0])/ 5.0) + 1
    endLon = np.ceil((180 + lonlim[1])/ 5.0) + 1
    rangeLon = np.arange(startLon, endLon, 1)
    rangeLat = np.arange(startLat, endLat, 1)

    name = []

    # make the names of the files that must be downloaded
    for lonname in rangeLon:
        for latname in rangeLat:
            name.append(str("srtm_%02d_%02d.zip" %(lonname, latname)))
            
    return(name, rangeLon, rangeLat)


def Download_Data(nameFile, output_folder_trash):
    """
    This function downloads the DEM data from the HydroShed website

    Keyword Arguments:
    nameFile -- name, name of the file that must be downloaded
    output_folder_trash -- Dir, directory where the downloaded data must be
                           stored
    """

    # download data from the internet
    url= "http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/%s" %(nameFile) 
    
    try:    
        file_name = url.split('/')[-1]
        output_file = os.path.join(output_folder_trash, file_name)
        if sys.version_info[0] == 3:
            urllib.request.urlretrieve(url, output_file)
        if sys.version_info[0] == 2:
            urllib.urlretrieve(url, output_file)
    except:
        pass

    return(output_file, file_name)


def Find_Document_names_15s_30s(latlim, lonlim, parameter, resolution):

    continents = ['na','ca','sa','eu','af','as','au']
    continents_download = []

    for continent in continents:
        extent = DEM_15s_extents.Continent[continent]
        if (extent[0] < lonlim[0] and extent[1] > lonlim[0] and extent[2] < latlim[0] and extent[3] > latlim[0]) and (extent[0] < lonlim[1] and extent[1] > lonlim[1] and extent[2] < latlim[1] and extent[3] > latlim[1]) == True:
            if resolution == "15s":
                name = '%s_%s_%s_grid.zip' %(continent, parameter, resolution)
            if resolution =="30s":
                name = '%s_%s_%s_bil.zip' %(continent, parameter, resolution)                
            continents_download = np.append(continents_download,name)

    return(continents_download)

class DEM_15s_extents:
        Continent = { 'na' : [-138, -52, 24, 60],
					'ca' : [-119, -60, 5, 39],
                      'sa' : [-93, -32, -56, 15],
                      'eu' : [-14, 70, 12, 62],
                      'af' : [-19, 55, -35, 38],
                      'as' : [57, 180, -12, 61],
                      'au' : [112, 180, -56, -10]}