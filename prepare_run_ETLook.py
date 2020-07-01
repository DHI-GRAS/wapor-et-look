# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 09:31:32 2020

@author: rmgu
"""

import sys

from pyWAPOR import ConfigFileInterface

if __name__ == '__main__':
    args = sys.argv
    if len(args) > 1:
        config_file = args[1]
    else:
        config_file = "example_config.toml"

    ConfigFileInterface.prepare_and_run_ETLook(config_file)
