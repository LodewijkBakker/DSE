############################################################################
#        All Constants Required for EPS sizing are defined here            #       
############################################################################

#--------------------------- Libraries ------------------------------------#

import numpy as np
from EPS_com import *

#------------------------- Orbit Related ----------------------------------#

Mission_life = 5                         # Mission Lifetime in years
T_orbit = 1*3600+30*60+31                # Orbital Period in seconds
# print('Orbital Period = ', T_orbit, 'seconds') # 5431 seconds
perc_eclipse = 0.404228                  # Percentage of Orbital Period in Eclipse
T_eclipse = perc_eclipse*T_orbit         # Eclipse Period in seconds
T_sunlight = T_orbit - T_eclipse         # Sunlight Period in seconds
S_flux = 1322                            # Solar Flux in W/m^2 - COLD from THERMO SIMULATION

#------------------------- Power Related ----------------------------------#

Peak_Power_range = np.arange(100, 1010, 10)  # Range of Peak Power in Watts

if PSPEPP:
    P_sun = Peak_Power_range
    P_eclipse = Peak_Power_range

if not PSPEPP:
    n_eclipse_PSPPEPP = 0.8                      # Chosen at random for now
    P_sun = Peak_Power_range
    P_eclipse = n_eclipse_PSPPEPP*Peak_Power_range     

#------------------------- Efficiencies, Errors ----------------------------#

n_cell = 0.338                                  # Solar Cell Efficiency NEW SMAT book tab.21-13
n_packing = 0.9                                 # Packing Efficiency according to eq.10-10 Space System Engineering third edition
dtheta_arr = np.radians(1)                      # Array Pointing Error in radians eq.10-10 Space System Engineering third edition
A_degradation_silicon = 3.75/100                # Array Degradation Factor per year according to NEW SMAT book p647 step 5 (Higher estimate for LEO silicon arrays)
A_degradation_galium_arsenite = 2.75/100        # Array Degradation Factor per year according to NEW SMAT book p647 step 5 (Higher estimate for LEO galium arsenite arrays)
D_arr = (1-A_degradation_silicon)**Mission_life # Array Degradation Factor over Mission Lifetime

X_e = 0.6                                       # Eclipse Factor accoring to SMAD book p. 643 step 2           
X_d = 0.8                                       # Day Factor accoring to SMAD book p. 643 step 2 
I_d = 0.72                                      # Total Inherent Degradation using SMAD tab.21-14




