# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Contact: timhessels@hotmail.com
Repository: WAPOR
Module: WAPOR/Collect/MSGCPP

Description:
This module downloads daily or instantanious MSGCPP data from
http://msgcpp.knmi.nl/portal/ server. The MSGCPP data is available since 2017-01-01 till the present.
The datasets will be stored in the user defined outputfolder in GEOTIFF format.


Examples:
from pyWAPOR.Collect import MSGCPP
MSGCPP.SDS(Dir='C:/Temp/', Startdate='1999-02-01', Enddate='1999-02-03',
             latlim=[-10, 30], lonlim=[-20, 120], Time = "13:45")		 
"""

from .SDS import main as SDS

__all__ = ['SDS']

__version__ = '0.1'
