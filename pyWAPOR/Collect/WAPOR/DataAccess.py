# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/WAPOR
"""

import requests, json, time
import pandas as pd
import os
import urllib
import gdal
import numpy as np

def DownloadData(output_folder, Startdate, Enddate, latlim, lonlim, Waitbar = 1):

    import pyWAPOR.Functions.Processing_Functions as PF
    
    output_folder_Tot = os.path.join(output_folder, "WAPOR", "LandCover")
    if not os.path.exists(output_folder_Tot):
        os.makedirs(output_folder_Tot)

    # Define dates
    Dates = pd.date_range(Startdate, Enddate, freq = "AS")
 
    # Create Waitbar
    if Waitbar == 1:
        import pyWAPOR.Functions.WaitbarConsole as WaitbarConsole
        total_amount = len(Dates)
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)    
       
    # Loop over dates
    for Date in Dates:
        
        file_name = os.path.join(output_folder_Tot, "LC_WAPOR_%s.01.01.tif" %Date.year)
     
        if not os.path.exists(file_name):
        
            # Get the data
            try:
                Array, geo, proj, file_name_temp = Download_data(output_folder_Tot, Date, latlim, lonlim)
            except:
                print("Was not able to create %s" %file_name)
            # Save the data
            PF.Save_as_tiff(file_name, Array, geo, proj)
            
            # remove temp folder
            os.remove(file_name_temp)
                    
        if Waitbar == 1:
            amount += 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)            
                        
    return()

def Download_data(output_folder_Tot, Date, latlim, lonlim):
    
    polygon = [[[lonlim[0],latlim[0]],[lonlim[0],latlim[1]],[lonlim[1],latlim[1]],[lonlim[1],latlim[0]],[lonlim[0],latlim[0]]]]
    url = 'https://io.apps.fao.org/gismgr/api/v1/query/'
      
    header = {
          "Content-type": "application/json;charset=UTF-8",
          "Accept": "application/json"
      } 
    
    # Create request
    payload ={
      "type": "CropRaster",
      "params": {
        "properties": {
          "outputFileName": "LC_WAPOR_%s.01.01.tif" %Date.year,
          "cutline": True,
          "tiled": True,
          "compressed": True,
          "overviews": True
        },
        "cube": {
          "code": "L1_LCC_A",
          "workspaceCode": "WAPOR",
          "language": "en"
        },
        "dimensions": [
          {
            "code": "YEAR",
            "values": [
              "[%s-01-01,%s-01-01)" %(int(Date.year), int(Date.year)+1)
            ]
          }
        ],
        "measures": [
          "LCC"
        ],
        "shape": {
          "type": "Polygon",
          "coordinates": polygon
        }
      }
    }
    
    # Download the data
    response = requests.post(url, data=json.dumps(payload), headers=header)
    response.raise_for_status()      
      
    response_json = response.json()
    result = response_json['response']
    
    job_url = result['links'][0]['href']
          
    for tries in range(20):

        file_name_temp = os.path.join(output_folder_Tot, "LC_WAPOR_%s.01.01_temp.tif" %Date.year)
        
        if not os.path.exists(file_name_temp):
            time.sleep(5)
            job_response = requests.get(job_url, headers=header)
            if job_response.status_code == 200:
                try:
                    job_result = job_response.json()['response']['output']['downloadUrl']
                    urllib.request.urlretrieve(job_result, file_name_temp)  
                    dest = gdal.Open(file_name_temp)
                    geo = dest.GetGeoTransform()
                    proj = dest.GetProjection()
                    Array = dest.GetRasterBand(1).ReadAsArray()
                    Array = np.float_(Array)        
                except:
                    pass
        
    return(Array, geo, proj, file_name_temp)