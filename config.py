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

except ModuleNotFoundError as e:
    print('ModuleNotFoundError: Missing fundamental packages (required: numpy, pandas')
    print(e)

"""
q_path: str, path where the .b16 WaSim data (detail input format)
sy_folder: str, path of folder where the .txt files with the total soil yield data for each sub-catchment 
    It must have a .txt file for the following sub-catchments: Devoll, Holta, Zalli and Skebices
    Each .txt file name must contain the name of the sub-catchment. 
turbine_capacity: float, maximum turbine capacity

time_interval: int, with the following possible values: 
    * 0 to keep original time interval
    * 1 to resample original flow data to daily data
    * 2 to resample original flow data to monthly data
    
*Constant input data, that must only be changed if more subcatchments are to be added, or the order in the final files
 should be changed: 

catchment_order: list of strings, with the names of the sub-catchments to consider, in the order they should be 
considered
soil_density: float, soil density to be used to calculate monthly volume concentration (m3/m3)
"""
q_path = r'Y:\Abt1\hiwi\Oreamuno\Tasks\01_BC_Calculation\Example_Files\WaSim\y_opt.b16'
sy_folder = r'Y:\Abt1\hiwi\Oreamuno\Tasks\01_BC_Calculation\Python\Soil_Input'

turbine_capacity = 108.02

time_interval = 0

catchment_order = ['Devoll', 'Holta', 'Zalli', 'Skebices']
soil_density = 2650