import sys
from pyWAPOR.Collect.MYD13.DataAccess import DownloadData


def main(Dir, Startdate, Enddate, latlim, lonlim, username, password, Waitbar = 1, hdf_library = None, remove_hdf = 1):
    """
    This function downloads MYD13 16-daily data for the specified time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    username -- "" string giving the username of your NASA account (https://urs.earthdata.nasa.gov/)
    password -- "" string giving the password of your NASA account    
    Waitbar -- 1 (Default) will print a waitbar
    hdf_library -- string, if all the hdf files are already stored on computer
                    define directory to the data here
    remove_hdf -- 1 (Default), if 1 remove all the downloaded hdf files in the end    
    """
    print('\nDownload 16-daily MODIS NDVI data for period %s till %s' %(Startdate, Enddate))
    DownloadData(Dir, Startdate, Enddate, latlim, lonlim, username, password, Waitbar, hdf_library, remove_hdf)

if __name__ == '__main__':
    main(sys.argv)