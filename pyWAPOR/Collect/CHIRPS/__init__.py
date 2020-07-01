# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Contact: timhessels@hotmail.com
Repository: WAPOR
Module: WAPOR/Collect/CHIRPS

Description:
This module downloads daily CHIRPS 2.0 data from
ftp://chg-ftpout.geog.ucsb.edu server. The CHIRPS data is available since 1981-01-01 till the present.
The datasets will be stored in the user defined outputfolder in GEOTIFF format.


Examples:
from pyWAPOR.Collect import CHIRPS
CHIRPS.daily(Dir='C:/Temp/', Startdate='1999-02-01', Enddate='1999-02-03',
             latlim=[-10, 30], lonlim=[-20, 120])
"""

from .daily import main as daily

__all__ = ['daily']

__version__ = '0.1'
