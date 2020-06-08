# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/SRTM


Description:
This module downloads DEM data from http://earlywarning.usgs.gov/hydrodata/.
Use the DEM functions to download and create DEM images in Gtiff format.

Examples:
from pyWAPOR.Collect import SRTM
SRTM.DEM(Dir='C:/TempDEM4/', latlim=[29, 32], lonlim=[-113, -109])
"""
from .DEM import main as DEM

__all__ = ['DEM']

__version__ = '0.1'
