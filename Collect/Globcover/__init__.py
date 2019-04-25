# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/Globcover


Description:
This module downloads Landuse data from http://due.esrin.esa.int/page_globcover.php.
Use the Globcover functions to download and create landuse images in Gtiff format.

Examples:
from pyWAPOR.Collect import Globcover
Globcover.Landuse(Dir='C:/TempDEM4/', latlim=[29, 32], lonlim=[-113, -109])
"""
from .Landuse import main as Landuse

__all__ = ['Landuse']

__version__ = '0.1'
