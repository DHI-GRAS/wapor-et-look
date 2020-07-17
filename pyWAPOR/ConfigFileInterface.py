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
        username_vito = config["username_vito"]
        password_vito = config["password_vito"]
        landcover = config["landcover"]
        prepare_level1_data = config["prepare_level1_data"]
        prepare_level2_data = config["prepare_level2_data"]
        run_level1_ETLook = config["run_level1_ETLook"]
        run_level2_ETLook = config["run_level2_ETLook"]
    except KeyError as e:
        print("Missing configuration parameter %s" % str(e))
        return

    # Prepare input data
    if prepare_level1_data:
        Pre_ETLook.prepare_level1(data_folder,
                                  start_date,
                                  end_date,
                                  latitude_extent,
                                  longitude_extent,
                                  username_earthdata,
                                  password_earthdata,
                                  landcover)
    if prepare_level2_data:
        Pre_ETLook.prepare_level2(data_folder,
                                  start_date,
                                  end_date,
                                  latitude_extent,
                                  longitude_extent,
                                  username_vito,
                                  password_vito,
                                  username_earthdata,
                                  password_earthdata,
                                  landcover)

    # Run the model
    if run_level1_ETLook:
        print(" ")
        for date in pd.date_range(start_date, end_date, freq="D"):
            print("Running ETLook Level 1 for %s " % date.strftime("%Y-%m-%d"))
            ETLook.ETLook_code.main(os.path.join(data_folder, "ETLook_input", "level_1"),
                                    os.path.join(data_folder, "ETLook_output", "level_1"),
                                    date)
    if run_level2_ETLook:
        print(" ")
        for date in pd.date_range(start_date, end_date, freq="D"):
            print("Running ETLook Level 2 for %s " % date.strftime("%Y-%m-%d"))
            ETLook.ETLook_code.main(os.path.join(data_folder, "ETLook_input", "level_2"),
                                    os.path.join(data_folder, "ETLook_output", "level_2"),
                                    date)
