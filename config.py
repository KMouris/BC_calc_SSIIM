"""
Module contains the python3 and external libraries needed to run the "main_bc.py" file, as well as the input data
needed.
"""

try:
    import os
    import sys

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

turbine_capacity: float, maximum turbine capacity. 
"""
q_path = r'Y:\Abt1\hiwi\Oreamuno\Tasks\01_BC_Calculation\Example_Files\WaSim\y_opt.b16'

turbine_capacity = 108.02
