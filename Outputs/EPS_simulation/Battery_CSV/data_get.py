####################################################################################################
# ---------------------------------- -Data Organisation -------------------------------------------#
# Function which compiles the most promissing configurations from battery_sizings.csv and outputs  #
# the mass of solar panels, batteries and harness estimation for each of the promissing configs.   #
####################################################################################################
from colored import fg
red = fg('red')
green = fg('green')
white = fg('white')

import pandas as pd
import numpy as np
import csv


''' Important visual inspection of the data has revealed that : 
        -> At least 1.5 side folds are required to cover 6h-12h orbits
        -> Largest battery masses occur at 12h and 9h
    >>> Thus the researched data from column 1 are the keys '1.5' and '12', '9' 
'''

# Importing the data from the csv file (remove header i.e row 0) and only keep the rows of interest (col 0 contains '1.5' and '12' or '9')

# Read the CSV file into a pandas DataFrame
df = pd.read_csv('./battery_sizings.csv')

# Filter the DataFrame based on the condition
condition = df['Config'].str.contains('1.5') & (df['Config'].str.contains('12') | df['Config'].str.contains('9') | df['Config'].str.contains('6'))
filtered_df = df[condition]

# Write the filtered DataFrame back to the CSV file
filtered_df.to_csv('./updated_promisssing_values.csv', index=False)

# Read the CSV file into a pandas DataFrame
df = pd.read_csv('./updated_promisssing_values.csv')

# Find the row corresponding to the most constraining propulsion battery capacity & remaining capacity
id_max_C_BOL_prop = df['C_BOL_prop'].idxmax()
id_max_C_BOL_remain = df['C_BOL_remain'].idxmax()

# --- Display the final results of the Simulation --- #

row_max_C_BOL_prop = df.iloc[id_max_C_BOL_prop]
row_max_C_BOL_remain = df.iloc[id_max_C_BOL_remain]

m_prop_bat = row_max_C_BOL_prop['M_prop_bat']
m_remain_bat = row_max_C_BOL_remain['M_remain_bat']

v_prop_bat = row_max_C_BOL_prop['V_prop_bat']
v_remain_bat = row_max_C_BOL_remain['V_remain_bat']

# Total mass and volume of the batteries
m_tot_bat = m_prop_bat + m_remain_bat
v_tot_bat = v_prop_bat + v_remain_bat

# Mass of the solar panels 
mass_p_cell = 53/2/1000     # [kg] Mass of a single solar panel cell   
cells_top = 17
cells_side = 19
cell_bacl = 21

num_side = 1.5
num_top = 1
num_back = 3

m_top_SP = cells_top*num_top*mass_p_cell
m_side_SP = 2*cells_side*num_side*mass_p_cell
m_back_SP = cell_bacl*num_back*mass_p_cell
m_tot_SP = m_top_SP + m_side_SP + m_back_SP

# Mass of the harness = 25% of the total mass of the EPS
mass_pcdu = 0.3

m_harness = 0.25*(m_tot_bat + m_tot_SP + mass_pcdu)

# --- Display the final results of the Simulation --- #
print('The total mass of the batteries is: ', m_tot_bat, 'kg')
print('                                   >>> Battery mass for propulsion system maintenance : ', m_prop_bat, 'kg')
print('                                   >>> Battery mass for remaining system maintenance  : ', m_remain_bat, 'kg')    
print('The total volume of the batteries is: ', v_tot_bat, 'U')
print('                                   >>> Battery volume for propulsion system maintenance : ', v_prop_bat, 'm^3')
print('                                   >>> Battery volume for remaining system maintenance  : ', v_remain_bat, 'm^3')
print('The total mass of the solar panels is: ', m_tot_SP, 'kg')
print('                                   >>> Mass of the top solar panels (1 fold)       : ', m_top_SP, 'kg')
print('                                   >>> Mass of the side solar panels (2x1.5 folds) : ', m_side_SP, 'kg')
print('                                   >>> Mass of the front solar panels (1x3 folds)  : ', m_back_SP, 'kg')
print('The total mass of the harness is: ', m_harness, 'kg')

print(green+'--------------------------------------------------------------------------------')
print('The total mass of the system is: ', m_tot_bat + m_tot_SP + m_harness + mass_pcdu, 'kg')
print('--------------------------------------------------------------------------------'+white)
