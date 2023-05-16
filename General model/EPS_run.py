from EPS_ import Solar_Array
import numpy as np
import matplotlib.pyplot as plt
import toml
import pandas as pd
import os
import sys
import time
from colored import fg
green = fg('green')
white = fg('white')

config_data_path = os.path.join(os.path.dirname(__file__), './Inputs/Config_eps.toml')
sim_data_path = os.path.join(os.path.dirname(__file__), './Outputs/30_Equatorial_Orbits.csv')
# Initialise Solar Array
SA = Solar_Array(config_data_path, sim_data_path)

AS = SA.Array_size()
# print(green + 'Array size = %.5f m^2' % AS + white)
Xtr = SA.Xtra_AR_from_SA()
# print(Xtr)
BAT = SA.E_bat()
# print(BAT/3600)
# print('------------------')
Size_new = SA.Array_size_improved(ALL_FIX=False)
print(Size_new)
# SA.plot_SA_efficiency_vs_Array_size()
# print(SA.Power_SA())

# dt = pd.read_csv(sim_data_path)
# print(np.shape(np.where(dt["Eclipse"] == 1))[1])

# print(r'Scenario 1 = %lf m^2' % (SA.Array_size_worst()))