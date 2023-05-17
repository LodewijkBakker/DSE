####################################################################################################
# Description : Running module for the Solar Array Sizing                                          #
####################################################################################################

# ----------------------------------- Imports -------------------------------------------------------#
import os 
import numpy as np
from colored import fg
green = fg('green')
white = fg('white')
red = fg('red')
from Solar_Arrays import Solar_Array

# ---------------------------------- Controls -------------------------------------------------------#
PLOT = False
# ----------------------------------- Main ----------------------------------------------------------#

# Initialise Solar Array
param_data_path = os.path.join(os.path.dirname(__file__), './Inputs/SA_Config.toml')
SA = Solar_Array(param_data_path)


# EPS power consumption
EPS_check = SA.EPS_P()
if abs(EPS_check - SA.P_EPS[0][0]) > 0.0001:
    print(f'{red}ERROR: EPS power consumption is {SA.P_EPS[0][0]} W, but should be {EPS_check} W{white}')
    exit()

# ---- Scenario Description ---- #
SA.disp_sc()

# ---- Scenario 1: See Drawing ---- #
Area_1 = SA.SA_Area_I1(PLOT)
P_avg_1 = SA.P_avg(scenario=1)
print(f'Scenario 1 : Average Power to generate per orbit = {green}', P_avg_1, f'W{white}')
# efficiency of 0.3 for comparison
idx_1 = np.where((Area_1[0]['cell efficiency'] <= 0.3+0.0001) & (Area_1[0]['cell efficiency'] >= 0.3-0.0001))[0][0]
print(f'             An efficiency of{green} {Area_1[0]["cell efficiency"][idx_1]} % :{white} Solar panel area{green} {Area_1[1]["Area"][idx_1]} m\u00b2 {white}')

# ---- Scenario 2: See Drawing ---- #
Area_2 = SA.SA_Area_I2(PLOT)
P_avg_2 = SA.P_avg(scenario=2)
print(f'\nScenario 2 : Average Power to generate per orbit = {green}', P_avg_2, f'W{white}')
# efficiency of 0.3 for comparison
idx_1 = np.where((Area_2[0]['cell efficiency'] <= 0.3+0.0001) & (Area_2[0]['cell efficiency'] >= 0.3-0.0001))[0][0]
print(f'             An efficiency of{green} {Area_2[0]["cell efficiency"][idx_1]} % :{white} Solar panel area{green} {Area_2[1]["Area"][idx_1]} m\u00b2 {white}')


# ---- Scenario 3: See Drawing ---- #
Area_3 = SA.SA_Area_I3(PLOT)
P_avg_3 = SA.P_avg(scenario=3)
print(f'\nScenario 3 : Average Power to generate per orbit = {green}', P_avg_3, f'W{white}')
# efficiency of 0.3 for comparison
idx_3 = np.where((Area_3[0]['cell efficiency'] <= 0.3+0.0001) & (Area_3[0]['cell efficiency'] >= 0.3-0.0001))[0][0]
print(f'             An efficiency of{green} {Area_3[0]["cell efficiency"][idx_3]} % :{white} Solar panel area{green} {Area_3[1]["Area"][idx_3]} m\u00b2 {white}')


# ---- Scenario 4: See Drawing ---- #
Area_4 = SA.SA_Area_I4(PLOT)
P_avg_4 = SA.P_avg(scenario=4)
print(f'\nScenario 4 : Average Power to generate per orbit = {green}', P_avg_4, f'W{white}')
# efficiency of 0.3 for comparison
idx_4 = np.where((Area_4[0]['cell efficiency'] <= 0.3+0.0001) & (Area_4[0]['cell efficiency'] >= 0.3-0.0001))[0][0]
print(f'             An efficiency of{green} {Area_4[0]["cell efficiency"][idx_4]} % :{white} Solar panel area{green} {Area_4[1]["Area"][idx_4]} m\u00b2 {white}')

# ---- Scenario 5: See Drawing ---- #
Area_5 = SA.SA_Area_I5(PLOT)
P_avg_5 = SA.P_avg(scenario=5)
print(f'\nScenario 5 : Average Power to generate per orbit = {green}', P_avg_5, f'W{white}')
# efficiency of 0.3 for comparison
idx_5 = np.where((Area_5[0]['cell efficiency'] <= 0.3+0.0001) & (Area_5[0]['cell efficiency'] >= 0.3-0.0001))[0][0]
print(f'             An efficiency of{green} {Area_5[0]["cell efficiency"][idx_5]} % :{white} Solar panel area{green} {Area_5[1]["Area"][idx_5]} m\u00b2 {white}')


# ---- Battery sizing ---- #
Bat = SA.Bat_size_visual()



