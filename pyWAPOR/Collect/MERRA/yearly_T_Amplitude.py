# -*- coding: utf-8 -*-
import sys
from pyWAPOR.Collect.MERRA.DataAccess import DownloadData


def main(Dir, Years, latlim, lonlim, Waitbar = 1):
    """
    This function downloads MERRA temperature for the defined years

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Vars -- ['t2m', 'v2m'] 
    Years -- [2009, 2010]
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    Waitbar -- 1 (Default) Will print a waitbar
    """
    for Year in Years:

        Startdate = "%s-01-01" %(Year)
        Enddate = "%s-12-31" %(Year)
        
        if Waitbar == 1:
            print('\nDownloading yearly MERRA t2m amplitude data for the period %s till %s' %(Startdate, Enddate))

        # Download data
        DownloadData(Dir, 't2m', Startdate, Enddate, latlim, lonlim, "yearly", Period = '', username = '', password = '', Waitbar = 1)

if __name__ == '__main__':
    main(sys.argv)
