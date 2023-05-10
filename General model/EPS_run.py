from EPS_ import Solar_Array
import numpy as np
import matplotlib.pyplot as plt
import toml
import pandas as pd
import os
import sys
import time

config_data_path = os.path.join(os.path.dirname(__file__), './Inputs/Config_eps.toml')

# Initialise Solar Array
SA = Solar_Array(config_data_path)

AS = SA.Array_size()
print(AS)
SA.plot_SA_efficiency_vs_Array_size()
print(SA.Power_SA())
