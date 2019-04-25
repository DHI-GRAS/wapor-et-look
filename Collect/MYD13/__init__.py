# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MYD13

Description:
This module downloads MYD13 NDVI data from
http://e4ftl01.cr.usgs.gov/. Use the MYD13.NDVI function to
download and create 16 daily NDVI images in Gtiff format.
The data is available between 2000-02-18 till present.

Examples:
from pyWAPOR.Collect import MYD13
MYD13.NDVI(Dir='C:/Temp3/', Startdate='2003-12-01', Enddate='2003-12-20',
           latlim=[41, 45], lonlim=[-8, -5])
"""

from .NDVI import main as NDVI

__all__ = ['NDVI']

__version__ = '0.1'
