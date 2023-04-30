import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

Q_gen = 30  # [watt] Amount of heat generated inside the spacecraft
Q_ir = 228  # beta angle dependent 228 for < 30 218 above
a_albedo = 0.14  # beta angle dependent 0.14 for < 30 0.19 above
Q_sun_hot = 1414  # [W * m−2]
Q_sun_cold = 1322  # [W * m−2]
Q_albedo_hot = Q_sun_hot*a_albedo
Q_albedo_cold = Q_sun_cold*a_albedo

A_rad = 1e-2  # [m^2] Surface of radiations
sigma_boltzmann = 5.670374419*10**(-8)  # https://www.britannica.com/science/Stefan-Boltzmann-law check e notation
alpha_aborptivity = 0.96
epsilon_emissivity = 0.9
eps_ = 0.9
m_sat = 5  # [kg] mass of satellite
specific_heat = 800  # not a constant value specifically
alt_orbit = (6371+300)*1e3  # [m] nominal orbit altitude from center of earth
mu_earth = 3.986004418e14  # [km^3/s^-2] https://ssd.jpl.nasa.gov/astro_par.html
v_orbit = ((mu_earth/alt_orbit)**0.5)
print(v_orbit)
t_orbit = 2*np.pi*alt_orbit/v_orbit
print(t_orbit)
f_orbit = 0.2  # fraction in eclipse

Q_in_hot = Q_gen + Q_ir + Q_sun_hot + Q_albedo_hot
Q_in_cold = Q_gen + Q_ir + Q_sun_cold + Q_albedo_cold
Q_eclipse = Q_gen + Q_ir

T_t0 = 293.15
t_final = 15*60*60
# Q =  c_p * m * T_dot
# T_dot = Q/(c_p * m)
#
def rhs_thermal(t, y):
    # starting just out of eclipse:
    # to do it faster convert it to gate form probably
    Q_in = Q_in_hot
    if (t%t_orbit)/t_orbit > (1-f_orbit):
        Q_in = Q_eclipse

    Q = Q_in - A_rad*sigma_boltzmann*epsilon_emissivity*y[0]**4
    T_dot = Q/(specific_heat*m_sat)
    return T_dot

sol = solve_ivp(rhs_thermal, [0, t_final], [T_t0], max_step=1)
# max step is here to not force out vibrations artificially
plt.plot(sol.t, sol.y.T)
plt.show()