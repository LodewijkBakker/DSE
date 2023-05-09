import numpy as np
from scipy.special import erf
def Cd_calculator(theta, Vinc, A, A_ref):
    R = 8.314
    Ti = 1000
    sigma_a = 0.95
    Tw = 400
    m=28.97E-3
    S = 7# between 1 and 13 also Vinf/cm cm=sqrt(2R/mTi)
    S = Vinc / np.sqrt(2 * R * Ti / m)

    P = 1 / S * np.exp(-S ** 2 * np.cos(theta) ** 2)
    G = 1 / 2 / S ** 2
    Q = 1 + G
    Z = 1 + erf(S * np.cos(theta))



    Vre_Vinc=np.sqrt(0.5*(1+sigma_a*4*R*Tw/Vinc**2-1))
    Cd=(P/np.sqrt(np.pi)+Q*Z*np.cos(theta)+np.cos(theta)/2*Vre_Vinc*(np.sqrt(np.pi)*Z*np.cos(theta)+P))*A/A_ref
    return Cd


#Environmental constants
rho_median=2.30E-11
rho_max=4.39E-11
r=(6378.14+300)*10**3
mu=	3.986004418E14

#Satellite parameters
A_extended=0.330*(0.350+0.368)/2+0.220*0.350
A_cube=0.350*0.340
Cd_standard=2.2
Cd_max=3.04
m=54
print(A_extended)

v_inf=np.sqrt(mu/r)

Cd_standard=Cd_calculator(0,v_inf,A_extended,A_extended)
D_median=Cd_standard*0.5*rho_median*v_inf**2*A_extended*1.2
D_max=Cd_standard*0.5*rho_max*v_inf**2*A_extended*1.2

dv_median=D_median/m*5*365.25*24*3600
dv_max=D_max/m*5*365.25*24*3600
print(Cd_standard)
print(f"Drag median {D_median}N, Delta v_inf median {dv_median}")
print(f"Drag max {D_max}N, Delta v_inf max {dv_max}")
