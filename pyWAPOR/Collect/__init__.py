# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect

Description:
This module contains scripts used to download Level 1 data (data directly from web).

Products                      Dates                             Password
CHIRPS (daily)                1981/01/01-now                    -
SRTM                           -                                 -
GLDAS                         2000/01/01-now                    NASA
MCD43 (daily)                 2000/02/24-now                    NASA
MOD11 (daily)                 2000/02/24-now                    NASA
MYD11 (daily)                 2000/02/24-now                    NASA
MYD13 (16-daily)              2000/02/18-now                    NASA
MOD13 (16-daily)              2000/02/18-now                    NASA
MERRA
MSGCPP
GEOS

Examples:
from pyWAPOR import Collect
help(Collect)
dir(Collect)
"""

from pyWAPOR.Collect import CHIRPS, MOD11, MYD11, MCD43, MOD13, MYD13, MSGCPP, MERRA, Globcover, GEOS, SRTM, WAPOR

__all__ = ['CHIRPS', 'MOD11', 'MYD11', 'MCD43', 'MOD13', 'MYD13', 'MSGCPP' 'MERRA', 'SRTM', 'Globcover', 'GEOS', 'WAPOR']

__version__ = '0.1'
