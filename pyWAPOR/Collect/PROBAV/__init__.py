# -*- coding: utf-8 -*-
"""
Authors: Laust Færch
Module: Collect/PROBAV

Description:
This module downloads Proba-V NDVI data from
https://www.vito-eodata.be/PDF/datapool/.
You need to register as a user (which is free) to use this module

"""
from .PROBAV_S5 import main as PROBAV_S5

__all__ = ["PROBAV_S5"]

__version__ = '0.1'
