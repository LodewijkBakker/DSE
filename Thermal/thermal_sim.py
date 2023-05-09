import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

Q_gen = 238  # [watt] Amount of heat generated inside the spacecraft
Q_ir_hot = 228  # beta angle dependent 228 for < 30 218 above
Q_ir_cold = 218
a_albedo_hot = 0.14  # beta angle dependent 0.14 for < 30 0.19 above
a_albedo_cold = 0.19
Q_sun_hot = 1414  # [W * m−2]
Q_sun_cold = 1322  # [W * m−2]
Q_albedo_hot = Q_sun_hot * a_albedo_hot
Q_albedo_cold = Q_sun_cold*a_albedo_cold

A_solar = 0
A_rad_cube = 2*0.36*0.34 + 2*0.36*0.35 + 2*0.35*0.34 + A_solar  # [m^2] Surface of radiations
A_rad_combined = 2*0.36*(0.34-0.12) + 2*0.35*(0.34-0.12) + 0.35*0.36 + (0.368+0.35)/2*np.pi*0.33 + A_solar

# max in to get max heat, change this to actual max!!!
A_in_cube = 0.2122 + A_solar    # 0.36*0.34 + A_solar
A_in_combined = 0.259 + A_solar     # 0.36*(0.34-0.12) + (0.368+0.35)/2 * 0.33 + A_solar
A_in_cube_min = 0.35*(0.34-0.12) + (0.368+0.35)/2 * 0.33
A_in_combined_min = 0.35*0.34


sigma_boltzmann = 5.670374419*10**(-8)  # https://www.britannica.com/science/Stefan-Boltzmann-law check e notation
alpha_absorptivity = 0.25
epsilon_emissivity = 0.88
m_sat = 54  # [kg] mass of satellite
specific_heat = 890  # (for aluminium) not a constant value specifically
alt_orbit = (6371+300)*1e3  # [m] nominal orbit altitude from center of earth
mu_earth = 3.986004418e14  # [km^3/s^-2] https://ssd.jpl.nasa.gov/astro_par.html
v_orbit = ((mu_earth/alt_orbit)**0.5)
t_orbit = 2*np.pi*alt_orbit/v_orbit
f_orbit = 0.5  # fraction in eclipse, update this to actual value should at least be less than 0.5

Q_in_hot_cube = Q_gen + Q_ir_hot * A_in_cube + Q_sun_hot * A_in_cube * alpha_absorptivity + Q_albedo_hot * A_in_cube * alpha_absorptivity
Q_in_cold_cube = Q_gen + Q_ir_cold*A_in_cube + Q_sun_cold*A_in_cube*alpha_absorptivity + Q_albedo_cold*A_in_cube*alpha_absorptivity
Q_eclipse_cube = Q_gen + Q_ir_hot * A_in_cube_min  # absolute minimum in for combined

Q_in_hot_combined = Q_gen + Q_ir_hot * A_in_combined + Q_sun_hot * A_in_combined * alpha_absorptivity + Q_albedo_hot * A_in_combined * alpha_absorptivity
Q_in_cold_combined = Q_gen + Q_ir_cold*A_in_combined + Q_sun_cold*A_in_combined*alpha_absorptivity + Q_albedo_cold*A_in_combined*alpha_absorptivity
Q_eclipse_combined = Q_gen + Q_ir_cold*A_in_combined_min  # absolute minimum in for combined

# Constant approximation
T_hot_cube = (Q_in_hot_cube/(A_rad_cube*sigma_boltzmann*epsilon_emissivity))**(1/4)
T_cold_cube = (Q_eclipse_cube/(A_rad_cube*sigma_boltzmann*epsilon_emissivity))**(1/4)
print(T_hot_cube, T_cold_cube, 'Temperature range for cube payload')

T_hot_combined = (Q_in_hot_combined/(A_rad_combined*sigma_boltzmann*epsilon_emissivity))**(1/4)
T_cold_combined = (Q_eclipse_combined/(A_rad_combined*sigma_boltzmann*epsilon_emissivity))**(1/4)

print(T_hot_combined, T_cold_combined, 'Temperature range for telescope payload')

T_t0 = 288  # [K] initial temperature
t_final = 50*3600  # [s] final time
# Q =  c_p * m * T_dot
# T_dot = Q/(c_p * m)
#
def rhs_thermal(t, y, Q_in_hot, Q_eclipse, A_rad):
    # starting just out of eclipse:
    # to do it faster convert it to gate form probably
    Q_in = Q_in_hot
    if (t%t_orbit)/t_orbit > (1-f_orbit):
        Q_in = Q_eclipse

    Q = Q_in - A_rad*sigma_boltzmann*epsilon_emissivity*y[0]**4
    T_dot = Q/(specific_heat*m_sat)
    return T_dot

sol_cube_hot = solve_ivp(rhs_thermal, [0, t_final], [T_t0], max_step=10, min_step=3, args=(Q_in_hot_cube, Q_eclipse_cube, A_rad_cube))
sol_combined_hot = solve_ivp(rhs_thermal, [0, t_final], [T_t0], max_step=10, min_step=3, args=(Q_in_hot_combined, Q_eclipse_combined, A_rad_combined))

sol_cube_cold = solve_ivp(rhs_thermal, [0, t_final], [T_t0], max_step=10, min_step=3, args=(Q_in_cold_cube, Q_eclipse_cube, A_rad_cube))
sol_combined_cold = solve_ivp(rhs_thermal, [0, t_final], [T_t0], max_step=10, min_step=3, args=(Q_in_cold_combined, Q_eclipse_combined, A_rad_combined))

# max step is here to not force out vibrations artificially
plt.plot(sol_cube_hot.t, sol_cube_hot.y.T)
plt.plot(sol_combined_hot.t, sol_combined_hot.y.T)
plt.plot(sol_cube_cold.t, sol_cube_cold.y.T)
plt.plot(sol_combined_cold.t, sol_combined_cold.y.T)
plt.hlines([280, 310], 0, t_final, linestyles='dashed')
plt.legend(['Cube hot', 'Combined hot', 'Cube cold', 'Combined cold'])
plt.show()

