# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Contact: timhessels@hotmail.com
Repository: WAPOR
Module: WAPOR/Collect/GEOS

Description:
This module downloads daily or instantanious GEOS data from
https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim server. The GEOS data is available since 2007-12-01 till the present.
The datasets will be stored in the user defined outputfolder in GEOTIFF format.


Examples:
from pyWAPOR.Collect import GEOS
GEOS.daily(Dir='C:/Temp/', Vars =["t2m", "v2m"], Startdate='1999-02-01', Enddate='1999-02-03',
             latlim=[-10, 30], lonlim=[-20, 120])
GEOS.three_hourly(Dir='C:/Temp/', Vars =["t2m", "v2m"], Startdate='1999-02-01', Enddate='1999-02-03',
             latlim=[-10, 30], lonlim=[-20, 120], Period = [1,2,3,4,5,6,7,8])			 
"""

from .daily import main as daily
from .three_hourly import main as three_hourly

__all__ = ['daily', 'three_hourly']

__version__ = '0.1'
