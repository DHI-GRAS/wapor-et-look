# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Contact: timhessels@hotmail.com
Repository: WAPOR
Module: WAPOR/Collect/MERRA

Description:
This module downloads daily or instantanious MERRA data from
https://opendap.nccs.nasa.gov/dods/GEOS-5/MERRAero server. The MERRA data is available since 1979-01-01 till the present.
The datasets will be stored in the user defined outputfolder in GEOTIFF format.


Examples:
from pyWAPOR.Collect import MERRA
MERRA.daily(Dir='C:/Temp/', Vars = ["t2m", "v2m"], Startdate='2002-02-01', Enddate='2002-02-03',
             latlim=[-10, 30], lonlim=[-20, 120])
MERRA.three_hourly(Dir='C:/Temp/', Vars =["t2m", "v2m"], Startdate='2002-02-01', Enddate='2002-02-03',
             latlim=[-10, 30], lonlim=[-20, 120], Periods = [1,2,3,4,5])			 
"""

from .hourly_MERRA2 import main as hourly_MERRA2
from .daily_MERRA2 import main as daily_MERRA2
from .daily import main as daily
from .three_hourly import main as three_hourly
from. yearly_T_Amplitude import main as yearly_T_Amplitude

__all__ = ['daily', 'three_hourly', 'yearly_T_Amplitude']

__version__ = '0.1'
