# -*- coding: utf-8 -*-
"""
Authors: Laust Færch
Module: Collect/PROBAV

Description:
This module preprocesses Landsat(7+8) data downloaded externaly

"""
from .PrepareLandsat import main as PrepareLandsat

__all__ = ["PrepareLandsat"]

__version__ = '0.1'
