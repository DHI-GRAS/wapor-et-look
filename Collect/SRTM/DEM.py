# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/SRTM
"""
import os
from pyWAPOR.Collect.SRTM.DataAccess import DownloadData
import sys


def main(Dir, latlim, lonlim, Waitbar = 1):
    """
    Downloads HydroSHED data from http://srtm.csi.cgiar.org/download

    this data includes a Digital Elevation Model (DEM)
    The spatial resolution is 90m (3s)

    The following keyword arguments are needed:
    Dir -- 'C:/file/to/path/'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    Waitbar -- '1' if you want a waitbar (Default = 1)
    """

    # Create directory if not exists for the output
    output_folder = os.path.join(Dir, 'SRTM', 'DEM')
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Define the output map and create this if not exists
    nameEnd = os.path.join(output_folder, 'DEM_SRTM_m_3s.tif')

    if not os.path.exists(nameEnd):

        # Create Waitbar
        if Waitbar == 1:
            print('\nDownload SRTM altitude map with a resolution of 3s')
            import pyWAPOR.Functions.WaitbarConsole as WaitbarConsole
            total_amount = 1
            amount = 0
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

        # Download and process the data
        DownloadData(output_folder, latlim, lonlim)

        if Waitbar == 1:
            amount = 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    else:
        if Waitbar == 1:
            print("\nSRTM altitude map (3s) already exists in output folder")

if __name__ == '__main__':
    main(sys.argv)
