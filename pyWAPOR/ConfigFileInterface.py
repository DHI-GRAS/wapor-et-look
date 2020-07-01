# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 09:50:54 2020

@author: rmgu
"""

import os
import pandas as pd
import toml

from pyWAPOR import Pre_ETLook, ETLook


def prepare_and_run_ETLook(config_file):

    # Parse config file
    config = toml.load(config_file)
    try:
        data_folder = config["data_folder"]
        start_date = config["start_date"]
        end_date = config["end_date"]
        latitude_extent = config["latitude_extent"]
        longitude_extent = config["longitude_extent"]
        username_earthdata = config["username_earthdata"]
        password_earthdata = config["password_earthdata"]
        landcover = config["landcover"]
        prepare_data = config["prepare_data"]
        run_ETLook = config["run_ETLook"]
    except KeyError as e:
        print("Missing configuration parameter %s" % str(e))
        return

    # Prepare input data
    if prepare_data:
        Pre_ETLook.main(data_folder,
                        start_date,
                        end_date,
                        latitude_extent,
                        longitude_extent,
                        username_earthdata,
                        password_earthdata,
                        landcover)

    # Run the model
    if run_ETLook:
        print(" ")
        for date in pd.date_range(start_date, end_date, freq="D"):
            print("Running ETLook for %s " % date.strftime("%Y-%m-%d"))
            ETLook.ETLook_code.main(os.path.join(data_folder, "ETLook_input"),
                                    os.path.join(data_folder, "ETLook_output"),
                                    date)
