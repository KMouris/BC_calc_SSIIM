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
    from dateutil import parser
    import matplotlib.pyplot as mp
    from matplotlib.dates import DateFormatter

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
results_folder: str, path where the results are safed to.

plot_outflow_data: boolean, whether a plot for outflow is generated
plot_water_level: boolean, whether a plot for water level is generated
plot_fig_size: array, of two integers [width, height] for the size of the plots. It is recommended to have a larger 
width when handling large time frames to improve the quality of the plots.

turbine_capacity: float, maximum turbine capacity

time_interval: int, with the following possible values: 
    * 0 to keep original time interval
    * 1 to resample original flow data to daily data
    * 2 to resample original flow data to monthly data - Less accurate in seasonal water level calculations
    
wl_threshold: list of int, with the target water levels for each month. This array reads only the first 12 values
target_wl_upper_boundary: int, with the relative upper boundary to the target water level where the water level is kept 
steady
target_wl_lower_boundary: int, with the relative lower boundary to the target water level where the inflow is stored. 
target_wl_maximum: int, with the maximum water level of the reservoir. If the water level rises above the excess 
inflow is released through the spillway.
If the water level is below target_wl_lower_boundary, turbine is turned off. Turbine is run at half of the inflow level,
if the water level is between target_wl_lower_boundary and target water level from wl_threshold. When the level is in 
between the target water level and target_wl_upper_boundary the turbine is kept at inflow level or max capacity if 
exceeded. For water levels above the turbine is ran at full throttle. For water levels above the target_wl_maximum all 
excess inflow is dumped to keep the water level at a maximum level of target_wl_maximum.

restrict_timei_date: boolean, if true the timei is trimed to the timeframe of timei_date_start and timei_date_end
timei_date_start: str, with the start date for the timei file
timei_date_end: str, with the end date for the timei file

CONSTANT INPUT DATA: that must only be changed if more sub-catchments are to be added, or the order in the final files
 should be changed: 

catchment_order: list of strings, with the names of the sub-catchments to consider, in the order they should be 
considered

grsz_const: boolean, if true a constant concentration for all sediment fraction is implemented. If False, there will 
be only one fraction sediment_density: float, sediment density (kg/m3) to be used to calculate monthly volume 
concentration (m3/m3) """
q_path = r'C:\Users\Mouris\Desktop\BC_SSIIM\Test\Input\qgkobanja.b16'
q_storage = r'C:\Users\Mouris\Desktop\BC_SSIIM\Test\Input\storage_curve_2019.txt'
sy_folder = r'C:\Users\Mouris\Desktop\BC_SSIIM\Test\Input\Sediment_input'

results_folder = r'C:\Users\Mouris\Desktop\BC_SSIIM\Test\results'

plot_outflow_data = True
plot_water_level = True
plot_fig_size = [50, 25]

turbine_capacity = 108.02

time_interval = 1

wl_threshold = [168, 168, 168, 168, 173, 173, 173, 173, 173, 173, 168, 168]
# Boundaries are within target_volume ( + 0,1 wl)
# Emergency boundaries are within target_volume ( - 2 wl & 175 m asl)
target_wl_upper_boundary = 0.1
target_wl_lower_boundary = -2
target_wl_maximum = 175

restrict_timei_date = True
timei_date_start = "01-01-2018"
timei_date_end = "31-12-2020"

catchment_order = ['Devoll', 'Holta', 'Zalli', 'Skebices']

grsz_const = True
sediment_density = 2650
