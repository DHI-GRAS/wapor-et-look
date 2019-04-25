# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Fri Feb 15 09:14:33 2019
"""

import os
import osr
import gdal
import gzip
import zipfile
import numpy as np

def Extract_Data_gz(zip_filename, outfilename):
    """
    This function extract the zip files

    Keyword Arguments:
    zip_filename -- name, name of the file that must be unzipped
    outfilename -- Dir, directory where the unzipped data must be
                           stored
    """

    with gzip.GzipFile(zip_filename, 'rb') as zf:
        file_content = zf.read()
        save_file_content = open(outfilename, 'wb')
        save_file_content.write(file_content)
    save_file_content.close()
    zf.close()
    os.remove(zip_filename)
    
def Save_as_MEM(data='', geo='', projection=''):
    """
    This function save the array as a memory file

    Keyword arguments:
    data -- [array], dataset of the geotiff
    geo -- [minimum lon, pixelsize, rotation, maximum lat, rotation,
            pixelsize], (geospatial dataset)
    projection -- interger, the EPSG code
    """
    # save as a geotiff
    driver = gdal.GetDriverByName("MEM")
    dst_ds = driver.Create('', int(data.shape[1]), int(data.shape[0]), 1,
                           gdal.GDT_Float32)
    srse = osr.SpatialReference()
    if projection == '':
        srse.SetWellKnownGeogCS("WGS84")

    else:
        try:
            if not srse.SetWellKnownGeogCS(projection) == 6:
                srse.SetWellKnownGeogCS(projection)
            else:
                try:
                    srse.ImportFromEPSG(int(projection))
                except:
                    srse.ImportFromWkt(projection)
        except:
            try:
                srse.ImportFromEPSG(int(projection))
            except:
                srse.ImportFromWkt(projection)
    dst_ds.SetProjection(srse.ExportToWkt())
    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo)
    dst_ds.GetRasterBand(1).WriteArray(data)
    return(dst_ds)   
    
def Save_as_tiff(name='', data='', geo='', projection=''):
    """
    This function save the array as a geotiff

    Keyword arguments:
    name -- string, directory name
    data -- [array], dataset of the geotiff
    geo --  [minimum lon, pixelsize, rotation, maximum lat, rotation,
            pixelsize], (geospatial dataset)
    projection -- integer, the EPSG code
    """
    # Change no data values
    data[np.isnan(data)] = -9999
    
    # save as a geotiff
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(name, int(data.shape[1]), int(data.shape[0]), 1,
                           gdal.GDT_Float32, ['COMPRESS=LZW'])
    srse = osr.SpatialReference()
    if projection == '':
        srse.SetWellKnownGeogCS("WGS84")

    # Set the projection, which can be an EPSG code or a well known GeogCS
    else:
        try:
            if not srse.SetWellKnownGeogCS(projection) == 6:
                srse.SetWellKnownGeogCS(projection)
            else:
                try:
                    srse.ImportFromEPSG(int(projection))
                except:
                    srse.ImportFromWkt(projection)
        except:
            try:
                srse.ImportFromEPSG(int(projection))
            except:
                srse.ImportFromWkt(projection)
                
    # Save the tiff file
    dst_ds.SetProjection(srse.ExportToWkt())
    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo)
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds = None
    
    return()
    
def reproject_MODIS(input_name, epsg_to):
    '''
    Reproject the merged data file by using gdalwarp. The input projection must be the MODIS projection.
    The output projection can be defined by the user.

    Keywords arguments:
    input_name -- 'C:/file/to/path/file.tif'
        string that defines the input tiff file
    epsg_to -- integer
        The EPSG code of the output dataset
    '''
    
    # Define the output name
    name_out = ''.join(input_name.split(".")[:-1]) + '_reprojected.tif'
   
    src_ds = gdal.Open(input_name)
    
    # Define target SRS
    dst_srs = osr.SpatialReference()
    dst_srs.ImportFromEPSG(int(epsg_to))
    dst_wkt = dst_srs.ExportToWkt()
    
    error_threshold = 0.125  # error threshold --> use same value as in gdalwarp
    resampling = gdal.GRA_NearestNeighbour
    
    # Call AutoCreateWarpedVRT() to fetch default values for target raster dimensions and geotransform
    tmp_ds = gdal.AutoCreateWarpedVRT( src_ds,
                                   None, # src_wkt : left to default value --> will use the one from source
                                   dst_wkt,
                                   resampling,
                                   error_threshold )
    dst_ds = gdal.GetDriverByName('GTiff').CreateCopy(name_out, tmp_ds)
    dst_ds = None 

    return(name_out)
    
def clip_data(input_file, latlim, lonlim):
    """
    Clip the data to the defined extend of the user (latlim, lonlim)

    Keyword Arguments:
    input_file -- output data, output of the clipped dataset
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    """
    try:
        if input_file.split('.')[-1] == 'tif':
            dest_in = gdal.Open(input_file)
        else:
            dest_in = input_file
    except:
        dest_in = input_file

    # Open Array
    data_in = dest_in.GetRasterBand(1).ReadAsArray()

    # Define the array that must remain
    Geo_in = dest_in.GetGeoTransform()
    Geo_in = list(Geo_in)
    Start_x = np.max([int(np.floor(((lonlim[0]) - Geo_in[0])/ Geo_in[1])),0])
    End_x = np.min([int(np.ceil(((lonlim[1]) - Geo_in[0])/ Geo_in[1])),int(dest_in.RasterXSize)])

    Start_y = np.max([int(np.floor((Geo_in[3] - latlim[1])/ -Geo_in[5])),0])
    End_y = np.min([int(np.ceil(((latlim[0]) - Geo_in[3])/Geo_in[5])), int(dest_in.RasterYSize)])

    #Create new GeoTransform
    Geo_in[0] = Geo_in[0] + Start_x * Geo_in[1]
    Geo_in[3] = Geo_in[3] + Start_y * Geo_in[5]
    Geo_out = tuple(Geo_in)

    data = np.zeros([End_y - Start_y, End_x - Start_x])

    data = data_in[Start_y:End_y,Start_x:End_x]
    dest_in = None

    return(data, Geo_out)    
    
    
def Extract_Data(input_file, output_folder):
    """
    This function extract the zip files

    Keyword Arguments:
    output_file -- name, name of the file that must be unzipped
    output_folder -- Dir, directory where the unzipped data must be
                           stored
    """
    # extract the data
    z = zipfile.ZipFile(input_file, 'r')
    z.extractall(output_folder)
    z.close()    
    
def reproject_dataset_example(dataset, dataset_example, method=1):
    """
    A sample function to reproject and resample a GDAL dataset from within
    Python. The user can define the wanted projection and shape by defining an example dataset.

    Keywords arguments:
    dataset -- 'C:/file/to/path/file.tif' or a gdal file (gdal.Open(filename))
        string that defines the input tiff file or gdal file
    dataset_example -- 'C:/file/to/path/file.tif' or a gdal file (gdal.Open(filename))
        string that defines the input tiff file or gdal file
    method -- 1,2,3,4 default = 1
        1 = Nearest Neighbour, 2 = Bilinear, 3 = lanzcos, 4 = average
    """
    # open dataset that must be transformed
    try:
        if os.path.splitext(dataset)[-1] == '.tif':
            g = gdal.Open(dataset)
        else:
            g = dataset
    except:
            g = dataset
    epsg_from = Get_epsg(g)

    #exceptions
    if epsg_from == 9001:
        epsg_from = 5070

    # open dataset that is used for transforming the dataset
    try:
        gland = gdal.Open(dataset_example)
        epsg_to = Get_epsg(gland)
    except:
        gland = dataset_example
        epsg_to = Get_epsg(gland)

    # Set the EPSG codes
    osng = osr.SpatialReference()
    osng.ImportFromEPSG(epsg_to)
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(epsg_from)

    # Get shape and geo transform from example
    geo_land = gland.GetGeoTransform()
    col=gland.RasterXSize
    rows=gland.RasterYSize

    # Create new raster
    mem_drv = gdal.GetDriverByName('MEM')
    dest = mem_drv.Create('', col, rows, 1, gdal.GDT_Float32)
    dest.SetGeoTransform(geo_land)
    dest.SetProjection(osng.ExportToWkt())

    # Perform the projection/resampling
    if method is 1:
        gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_NearestNeighbour)
    if method is 2:
        gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Bilinear)
    if method is 3:
        gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Lanczos)
    if method is 4:
        gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Average)
    return(dest)

def Get_epsg(g):
    """
    This function reads the projection of a GEOGCS file or tiff file

    Keyword arguments:
    g -- string
        Filename to the file that must be read
    extension -- tiff or GEOGCS
        Define the extension of the dataset (default is tiff)
    """
    try:
        # Get info of the dataset that is used for transforming
        g_proj = g.GetProjection()
        Projection=g_proj.split('EPSG","')
        epsg_to=int((str(Projection[-1]).split(']')[0])[0:-1])            
    except:
        epsg_to=4326

    return(epsg_to)

def Show_tif(image_file, Limits = None, Color = None):
    """
    This function plot a tiff array in the console

    Keyword arguments:
    image_file -- string
        Filename to the file that must be shown
    Limits -- [min, max] (Default = min and max image)
        User can define the limits of the colorbar
    Color -- string (Default = "viridis")
        User can define the wanted colormap, all options are listed here:
        https://matplotlib.org/examples/color/colormaps_reference.html
    """    
    dest = gdal.Open(image_file)
    Array = dest.GetRasterBand(1).ReadAsArray()
    Array[Array==-9999] = np.nan
    if Limits == None:
        Limits = [np.nanmin(Array), np.nanmax(Array)]
    
    if Color == None:
        Color = "viridis"
    
    import matplotlib.pyplot as plt
    
    plt.imshow(Array, cmap = Color, vmin=Limits[0], vmax=Limits[1])
    plt.colorbar()
    plt.show()
    
    return()

    
    
    
    
