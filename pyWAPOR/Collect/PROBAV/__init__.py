# -*- coding: utf-8 -*-
"""
Authors: Laust Færch
Module: Collect/PROBAV

Description:
This module downloads Proba-V NDVI data from
https://www.vito-eodata.be/PDF/datapool/.
You need to register as a user (which is free) to use this module

"""

from .NDVI import main as NDVI
from .ALBEDO import main as ALBEDO

__all__ = ['NDVI', 'ALBEDO']

__version__ = '0.1'
