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
idx_2 = np.where((Area_2[0]['cell efficiency'] <= 0.3+0.0001) & (Area_2[0]['cell efficiency'] >= 0.3-0.0001))[0][0]
print(f'             An efficiency of{green} {Area_2[0]["cell efficiency"][idx_2]} % :{white} Solar panel area{green} {Area_2[1]["Area"][idx_2]} m\u00b2 {white}')


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


# ---- Configuration Computation ---- #
A_top = 305 * 360 / 1e6
A_back = 220 * 350 / 1e6
A_side = 200 * 360 / 1e6
# Cong. I, Config. II
A_back_I_II = Area_5[1]['Area'][idx_5]
A_side_I_II = 0.5 * (Area_1[1]['Area'][idx_1] - A_back_I_II - A_top)

# Cong. III, Config. IV
A_back_III_IV = Area_5[1]['Area'][idx_5]
A_side_III_IV = 0.5* (Area_3[1]['Area'][idx_3] - (A_top + A_back_III_IV)*2/np.pi * 0.81)

print(f'Back area required for Cong. I, Config. II = {green} ', A_back_I_II, f'm\u00b2 {white}')
print(f'Side area required for Cong. I, Config. II = {green} ', A_side_I_II, f'm\u00b2 {white}')
print(f'Back area required for Cong. III, Config. IV = {green} ', A_back_III_IV, f'm\u00b2 {white}')
print(f'Side area required for Cong. III, Config. IV = {green} ', A_side_III_IV, f'm\u00b2 {white}')

# ---- Equivalent Folds ---- #
# Cong. I
Folds_back_I = A_back_I_II / A_back
Folds_side_I = A_side_I_II / A_side
if Folds_back_I > 2:
    txt_I = f'Required folds for configuration I = {red} {Folds_back_I} {white} > 2'
else:
    txt_I = f'Required folds for configuration I = {green} {Folds_back_I} {white} < 2'
if Folds_side_I > 3:
    txt_I += f'\nRequired folds for configuration I = {red} {Folds_side_I} {white} > 3'
else:
    txt_I += f'\nRequired folds for configuration I = {green} {Folds_side_I} {white} < 3'
print(txt_I)

# Config. II
Folds_back_II = A_back_I_II / A_back
Folds_side_II = A_side_I_II / A_side
if Folds_back_II > 3:
    txt_II = f'Required folds for configuration II = {red} {Folds_back_II} {white} > 3'
else:
    txt_II = f'Required folds for configuration II = {green} {Folds_back_II} {white} < 3'
if Folds_side_II > 3:
    txt_II += f'\nRequired folds for configuration II = {red} {Folds_side_II} {white} > 3'
else:
    txt_II += f'\nRequired folds for configuration II = {green} {Folds_side_II} {white} < 3'
print(txt_II)

# Cong. III
Folds_back_III = A_back_III_IV / A_back
Folds_side_III = A_side_III_IV / A_side
if Folds_back_III > 2:
    txt_III = f'Required folds for configuration III = {red} {Folds_back_III} {white} > 2'
else:
    txt_III = f'Required folds for configuration III = {green} {Folds_back_III} {white} < 2'
if Folds_side_III > 3:
    txt_III += f'\nRequired folds for configuration III = {red} {Folds_side_III} {white} > 3'
else:
    txt_III += f'\nRequired folds for configuration III = {green} {Folds_side_III} {white} < 3'
print(txt_III)

# Config. IV
Folds_back_IV = A_back_III_IV / A_back
Folds_side_IV = A_side_III_IV / A_side
if Folds_back_IV > 3:
    txt_IV = f'Required folds for configuration IV = {red} {Folds_back_IV} {white} > 3'
else:  
    txt_IV = f'Required folds for configuration IV = {green} {Folds_back_IV} {white} < 3'
if Folds_side_IV > 3:
    txt_IV += f'\nRequired folds for configuration IV = {red} {Folds_side_IV} {white} > 3'
else:
    txt_IV += f'\nRequired folds for configuration IV = {green} {Folds_side_IV} {white} < 3'
print(txt_IV)

# ---- Battery sizing ---- #
Bat = SA.Bat_size_visual(CONT=True)



