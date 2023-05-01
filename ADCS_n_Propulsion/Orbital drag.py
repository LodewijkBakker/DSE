import numpy as np

# Environmental constants
rho_median= 2.3E-11
rho_max = 4.39E-11
r = (6371.8+300)*10**3
mu = 3.986004418E14

# Satellite parameters
A_extended = 0.330*(0.350+0.368)/2+0.220*0.350
A_cube = 0.350*0.340
Cd_standard = 2.2
Cd_max = 3.04
m = 54

v = np.sqrt(mu/r)
D_median = Cd_standard*0.5*rho_median*v**2*A_extended
D_max = Cd_standard*0.5*rho_max*v**2*A_extended

dv_median = D_median/m*5*365.25*24*3600
dv_max = D_max/m*5*365.25*24*3600
print(f"Drag median {D_median}N, Delta v median {dv_median}")
print(f"Drag max {D_max}N, Delta v max {dv_max}")
