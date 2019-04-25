# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MOD04

Description:
This module downloads MOD04 AOD550 data from
https://ladsweb.modaps.eosdis.nasa.gov. Use the MOD04.AOD550 function to
download and create daily AOD550 images in Gtiff format.
The data is available between 2000-02-18 till present.

Examples:
from pyWAPOR.Collect import MOD04
MOD04.AOD550(Dir='C:/Temp3/', Startdate='2003-12-01', Enddate='2003-12-20',
           latlim=[41, 45], lonlim=[-8, -5])
"""

from .AOD550 import main as AOD550

__all__ = ['AOD550']

__version__ = '0.1'
