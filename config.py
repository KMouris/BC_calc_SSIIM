"""
Module contains the python3 and external libraries needed to run the "main_bc.py" file, as well as the input data
needed.
"""

try:
    import os
    import sys
    import glob
    import calendar
    import datetime

except ModuleNotFoundError as b:
    print('ModuleNotFoundError: Missing basic libraries (required: os, sys')
    print(b)

try:
    import numpy as np
    import pandas as pd

except ModuleNotFoundError as e:
    print('ModuleNotFoundError: Missing fundamental packages (required: numpy, pandas')
    print(e)

"""
INPUT DATA: 
q_path: str, path where the .b16 WaSim data (detail input format)
q_storage: str, path where water level is linked to water storage volume
sy_folder: str, path of folder where the .txt files with the total sediment yield data for each sub-catchment 
    It must have a .txt file for the following sub-catchments: Devoll, Holta, Zalli and Skebices
    Each .txt file name must contain the name of the sub-catchment. 
turbine_capacity: float, maximum turbine capacity

time_interval: int, with the following possible values: 
    * 0 to keep original time interval
    * 1 to resample original flow data to daily data
    * 2 to resample original flow data to monthly data
    
winter_threshold: list of int, with the months where the winter water level is applied
summer_threshold: list of int, with the months where the summer water level is applied

upstream_wl_winter: float, with the constant upstream water level to be used in winter
downstream_wl_winter: float, with the constant downstream water level to be used in winter
upstream_wl_summer: float, with the constant upstream water level to be used in summer
downstream_wl_summer: float, with the constant downstream water level to be used in summer
    
CONSTANT INPUT DATA: that must only be changed if more sub-catchments are to be added, or the order in the final files
 should be changed: 

catchment_order: list of strings, with the names of the sub-catchments to consider, in the order they should be 
considered

sediment_density: float, sediment density (kg/m3) to be used to calculate monthly volume concentration (m3/m3)
"""
q_path = r'/home/yendras/hiwi/01_BC_Calculation/Example_Files/WaSim/y_opt.b16'
q_storage = r'/home/yendras/hiwi/01_BC_Calculation/Python/Input/storage_curve_2019.txt'
sy_folder = os.path.abspath(r'../Input/Soil_Input')

results_folder = os.path.abspath(r'../results')

log_outflow_data = True

turbine_capacity = 108.02

time_interval = 0

winter_threshold = [11, 12, 1, 2, 3, 4]
summer_threshold = [5, 6, 7, 8, 9, 10]

upstream_wl_winter = 168
downstream_wl_winter = 168
upstream_wl_summer = 173
downstream_wl_summer = 173


#timei_date_start = datetime.datetime(2020, 1, 1)
#timei_date_end = datetime.datetime(2025, 12, 31)

catchment_order = ['Devoll', 'Holta', 'Zalli', 'Skebices']
sediment_density = 2650