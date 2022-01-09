"""
Module contains the python3 and external libraries needed to run the "main_bc.py" file, as well as the input data
needed.
"""

try:
    import os
    import sys
    import glob
except ModuleNotFoundError as b:
    print('ModuleNotFoundError: Missing basic libraries (required: os, sys')
    print(b)

try:
    import numpy as np
    import pandas as pd
    import math

except ModuleNotFoundError as e:
    print('ModuleNotFoundError: Missing fundamental packages (required: numpy, pandas')
    print(e)

"""
INPUT DATA: 
q_path: str, path where the .b16 WaSim data (detail input format)

sy_folder: str, path of folder where the .txt files with the total soil yield data for each sub-catchment 
    It must have a .txt file for the following sub-catchments: Devoll, Holta, Zalli and Skebices
    Each .txt file name must contain the name of the sub-catchment. 
turbine_capacity: float, maximum turbine capacity

time_interval: int, with the following possible values: 
    * 0 to keep original time interval
    * 1 to resample original flow data to daily data
    * 2 to resample original flow data to monthly data
    
upstream_wl: float, with the constant upstream water level to be used (for all time steps)
downstream_wl: float, with the constant upstream water level to be used (for all time steps)
    * NOTE: If user wants to use time-dependent upstream or downstream water levels, the input file must be a .txt file
    with the water level value and the corresponding time step assigned to it. A function must be added (or modify 
    'upstream_downstream_data' function) to read the data, re-discretize it to the simulation time step (if needed) and
    and assign it to the timei file. 
    
CONSTANT INPUT DATA: that must only be changed if more sub-catchments are to be added, or the order in the final files
 should be changed: 

catchment_order: list of strings, with the names of the sub-catchments to consider, in the order they should be 
considered

sediment_density: float, sediment density to be used to calculate monthly volume concentration (m3/m3)

FOR SEASONAL OUTFLOW: to consider a seasonal outflow, based on water level
seasonal_wl: boolean, is False, only the 'turbine_capacity' variable would be needed, and it considers a constant water
level, if True, it considers a seasonal outflow, in which case the following variables are needed: 
    storage_curve_path: str, path (.txt) where the reservoir storage curve is saved, in a .txt, tab-delimited file, 
    with 2 columns: 1: water level, 2: storage volume
    initial_wl: float, initial water level (at the beginning of the simulation)
    h_max = float, maximum allowable water level (above this, the spillway begins to work)

    NOTE: this configuration is a suggestion, and any method, with any input data could be used, to consider a seasonal
    outflow. 

"""
q_path = r'Y:\Abt1\hiwi\Oreamuno\Tasks\01_BC_Calculation\Example_Files\WaSim\y_opt.b16'
sy_folder = os.path.abspath(r'../Input/Soil_Input')

resuls_folder = os.path.abspath(r'../Results')

turbine_capacity = 108.02

time_interval = 0

upstream_wl = 0
downstream_wl = 175

catchment_order = ['Devoll', 'Holta', 'Zalli', 'Skebices']
sediment_density = 2650

# For seasonal water level:
seasonal_wl = True
storage_curve_path = os.path.abspath(r'../Input/storage_curve_2019.txt')
initial_wl = 171.3
h_max = 175
