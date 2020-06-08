# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MCD43

Description:
This module downloads MCD43 ALBEDO data from
https://e4ftl01.cr.usgs.gov/. Use the MCD43.ALBEDO function to
download and create daily ALBEDO images in Gtiff format.
The data is available between 2000-02-18 till present.

Examples:
from pyWAPOR.Collect import MCD43
MCD43.ALBEDO(Dir='C:/Temp3/', Startdate='2003-12-01', Enddate='2003-12-20',
           latlim=[41, 45], lonlim=[-8, -5])
"""

from .ALBEDO import main as ALBEDO

__all__ = ['ALBEDO']

__version__ = '0.1'
