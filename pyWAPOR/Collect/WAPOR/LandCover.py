# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/SRTM
"""
from pyWAPOR.Collect.WAPOR.DataAccess import DownloadData
import sys

def main(Dir, Startdate, Enddate, latlim, lonlim, Waitbar = 1):
    """
    Downloads LandCover data from https://io.apps.fao.org/gismgr/api/v1/query/.

    The following keyword arguments are needed:
    Dir -- 'C:/file/to/path/'
	Startdate -- 'yyyy-mm-dd'
	Enddate -- 'yyyy-mm-dd'    
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    Waitbar -- '1' if you want a waitbar (Default = 1)
    """

    print('\nDownloading Yearly LandCover maps from WAPOR for the period %s till %s' %(Startdate, Enddate))

    # Download data
    DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar)

if __name__ == '__main__':
    main(sys.argv)
