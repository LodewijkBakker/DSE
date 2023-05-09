import numpy as np

# Environmental constants
rho = 2.3E-11
r = (6371.8+300)*10**3
mu = 3.986004418E14
g_o = 9.81

# Satellite parameters
A_extended = 0.330*(0.350+0.368)/2+0.220*0.350
A_cube = 0.350*0.340
Cd = 2.2
m_dry = 54
l = 2.7 # length of S/C
T_life = 5*365*86400 # 5 years lifetime
rho_p = 867 # Density of Propellant

v = np.sqrt(mu/r)
F_d = Cd*0.5*rho*v**2*A_extended # force due to drag
J = F_d*T_life # Total Impulse

# Chosen parameters
F_t = 0.03
Isp = 1000

print(f"The aerodynamic drag is {np.round(F_d, 6)} N")
print(f"The total impulse required is {np.round(J, 3)} N s")

D = F_d / F_t # Duty cycle

# Mass of Propulsion System
m_p = J / (g_o * Isp) # Propellant mass
# m_s = 10 + 0.3*m_p # Propulsion structure mass for HET/RF
m_s = 2 + 0.15*m_p
m_wet = m_s + m_p # wet mass of just propulsion system
m_ratio = m_wet / m_dry # Mass ratio

#Volume of Propulsion System
V_p = m_p / rho_p # Volume of propellant
V_s = 0.016 # Volume of structure
V_wet = V_s + V_p
V_ratio = V_wet / l**3 # Volume ratio

print(m_p, m_wet)
