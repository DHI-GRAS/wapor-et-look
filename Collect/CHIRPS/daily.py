# -*- coding: utf-8 -*-
import sys
from pyWAPOR.Collect.CHIRPS.DataAccess import DownloadData


def main(Dir, Startdate='', Enddate='',
			latlim=[-50, 50], lonlim=[-180, 180], Waitbar = 1):
	"""
	This function downloads daily CHIRPS data

	Keyword arguments:
	Dir -- 'C:/file/to/path/'
	Startdate -- 'yyyy-mm-dd'
	Enddate -- 'yyyy-mm-dd'
	latlim -- [ymin, ymax] (values must be between -50 and 50)
	lonlim -- [xmin, xmax] (values must be between -180 and 180)
	Waitbar -- 1 (Default) will print a waitbar
	"""
	print('\nDownloading daily CHIRPS rainfall data for the period %s till %s' %(Startdate, Enddate))

	# Download data
	DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar)

if __name__ == '__main__':
    main(sys.argv)
