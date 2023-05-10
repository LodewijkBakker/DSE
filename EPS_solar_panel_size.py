############################################################################
#                     Different methods solar panel sizing                 #       
############################################################################

#--------------------------- Libraries ------------------------------------#

from EPS_com import *
from EPS_cst import *
import numpy as np
import matplotlib.pyplot as plt

#-------------------------- First Method ----------------------------------#
### Assuming same power demamd in eclipse and sunlight  & (n_BDR, n_BCR, n_AR) -> (1, 1, 1) ###

P_array = P_sun * (T_orbit/T_sunlight) 
A_array = P_array/(S_flux*n_cell*n_packing*np.cos(dtheta_arr)*(1-D_arr))
print(D_arr)

#-------------------------- Second Method ---------------------------------#
### SMAT MANUAL METHOD ###

P_array_2 = ((P_eclipse*T_eclipse/X_e)+(P_sun*T_sunlight/X_d))/T_orbit
P_o = S_flux * n_cell                                                   # Power output with sun perpendicular to array
P_BOL_fixed = P_o * I_d * np.cos(45)                                    # Worst case scenario for the angle of incidence, assumming the average of a sinusoidal curve, with fixed arrays.
P_BOL_rotate = P_o * I_d * 1                                            # Best case scenario in case of rotating arrays
L_d = (1-A_degradation_galium_arsenite)**Mission_life                   # Degradation factor over mission lifetime
print(L_d)
P_EOL_fixed = P_BOL_fixed * L_d                                         # Power output at EOL with fixed arrays
P_EOL_rotate = P_BOL_rotate * L_d                                       # Power output at EOL with rotating arrays 
print('P_EOL_fixed = ', P_EOL_fixed, r'W/m^2')
print('P_EOL_rotate = ', P_EOL_rotate, r'W/m^2')
A_array_fixed = P_array_2/P_EOL_fixed                                     # Fixed arrays               
A_array_rotate = P_array_2/P_EOL_rotate                                   # Rotating arrays

#--------------------------- Graphing -----------------------------------#
if PLOTTING_SA:
### Plotting the Area of the Solar Panel vs Peak Power P_array split in two graphs ###
    # First half of P_array
    plt.figure()
    plt.plot(P_array[0:50], A_array[0:50],marker='x',markersize=5, label = 'Method 1')
    plt.plot(P_array_2[0:50], A_array_fixed[0:50],marker='x',markersize=5,  label = 'Method 2 - Fixed')
    plt.plot(P_array_2[0:50], A_array_rotate[0:50],marker='x',markersize=5, label = 'Method 2 - Rotating')
    plt.xlabel(r'Power Array $[W]$')
    plt.ylabel(r'Area of Solar Panel $[m^2]$')
    plt.title('Area of Solar Panel vs Array Power')
    plt.legend()
    plt.grid()
    plt.show()

    # Second half of P_array
    plt.figure()
    plt.plot(P_array[50:100], A_array[50:100],marker='x',markersize=5, label = 'Method 1')
    plt.plot(P_array_2[50:100], A_array_fixed[50:100],marker='x',markersize=5,  label = 'Method 2 - Fixed')
    plt.plot(P_array_2[50:100], A_array_rotate[50:100],marker='x',markersize=5, label = 'Method 2 - Rotating')
    plt.xlabel(r'Power Array $[W]$')
    plt.ylabel(r'Area of Solar Panel $[m^2]$')
    plt.title('Area of Solar Panel vs Array Power')
    plt.legend()
    plt.grid()
    plt.show()
