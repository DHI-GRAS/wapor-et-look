# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/WAPOR


Description:
This module downloads WAPOR data from https://io.apps.fao.org/gismgr/api/v1/query/.
Use the LandCover functions to download and create LC images in Gtiff format.

Examples:
from pyWAPOR.Collect import WAPOR
WAPOR.LandCover(Dir='C:/TempDEM4/', "2008-01-01", "2017-01-01",  latlim=[29, 32], lonlim=[-113, -109])
"""
from .LandCover import main as LandCover

__all__ = ['LandCover']

__version__ = '0.1'
