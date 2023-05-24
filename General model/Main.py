
import numpy as np

# Main satellite file to update all other values
vol_sat = 2e-3  # [m^3]
mass_sat = vol_sat*1000  # [kg], water method
vol_sat = 0.35*0.36*0.34  # [m^3]
mass_sat = 54  # [kg]

# First order delta V calculation
a_frontal = 1e-3  # [m^2] frontal area
alt_orbit = (6371+300)*1e3  # [m] nominal orbit altitude from center of earth
mu_earth = 398600.435507*1e3  # [km^3/s^-2] https://ssd.jpl.nasa.gov/astro_par.html
v_orbit = (alt_orbit/mu_earth)**0.5
a_frontal_cube = 0.360*0.34
a_frontal_cyl = 0.36*0.22+0.33*0.35+0.018*0.33/2  # [m^2] frontal area
a_frontal = max([a_frontal_cyl, a_frontal_cube])
alt_orbit = (6371.8+300)*1e3  # [m] nominal orbit altitude from center of earth
mu_earth = 3.98600435507*1e14  # [m^3/s^-2] https://ssd.jpl.nasa.gov/astro_par.html
v_orbit = (mu_earth/alt_orbit)**0.5
t_life = 5*365*24*60*60  # [s] nominal life time, 5 years
c_d_satellite = 3  # flat surfaces
rho = 1e-5  # [kg/m^3]
c_d_satellite = 2.2  # flat surfaces
rho = 4.39e-11  # [kg/m^3]

f_aero = rho*a_frontal*c_d_satellite*v_orbit**2  # [n] aerodynamic force at nominal orbit
f_aero = 0.5*rho*a_frontal*c_d_satellite*v_orbit**2  # [n] aerodynamic force at nominal orbit
print(f_aero)
# https://drive.google.com/drive/u/0/folders/1XzHVbNNy5GXSy6uyikV365bol00h7116
# assuming mass is constant, frontal area does not change and no change in altitude or other density factors
safety_factor = 1.5
safety_factor = 1  # 1.5
delta_v = f_aero / mass_sat * t_life * safety_factor # [m/s] delta v

print(delta_v)
# ADCS disturbance torques:
mag_field_max = 0
mag_dipole = 0
w_used = 20  # power used
mag_dist_torque = 0
mag_field_max = 57587e-9  # max of modelling of very low earth orbit source
mag_dipole = 3.5e-3  # 19690020961 pdf, it's based on amount of mass
mag_dist_torque = mag_dipole*mass_sat*mag_field_max  #

largest_arm = 0.01
a_largest_frontal_area = 1
largest_arm = 0.36/2
a_largest_frontal_area = 0.24

q_sun_max = 1414
#a_albedo = 0.14  # beta angle dependent 0.14 for < 30 0.19 above
c_light = 299,792,458  # [m/s] https://physics.nist.gov/cgi-bin/cuu/Value?c
a_reflectivity = 0.9
f_solar_radiation = q_sun_max*c_light*a_reflectivity*a_frontal
c_light = 299792458  # [m/s] https://physics.nist.gov/cgi-bin/cuu/Value?c
a_reflectivity = 0.6  # SOURCE, use as representative value smad but also
# https://www.sciencedirect.com/science/article/pii/S0273117718306197
f_solar_radiation = q_sun_max/c_light*(1+a_reflectivity)*a_frontal
sun_dist_torque = f_solar_radiation*largest_arm

f_aero = rho*a_largest_frontal_area*c_d_satellite*v_orbit**2
f_aero = 0.5*rho*a_largest_frontal_area*c_d_satellite*v_orbit**2
#print(f_aero)
# [n] aerodynamic force at nominal orbit, maximum around any axis
aero_dist_torque = f_aero*largest_arm

gravity_dist_torque = 0
# body mounted panels dominate moment of inertia so this has to be taken with a grain of salt
m_telescope = 30
m_total_cube = 54
m_bottom_cube = m_total_cube - m_telescope

h_bot_cube = 0.22
h_tot_cube = 0.34
w_tot_cube = 0.36
l_tot_cube = 0.35
r_dst = (0.368+0.35)/4  # check
h_dst = 0.33  #???

Izz_combined = (1/12)*m_bottom_cube*(w_tot_cube**2+l_tot_cube**2) + (1/2)*m_telescope*r_dst**2  # yaw
Iyy_combined = (1/12)*m_bottom_cube*(l_tot_cube**2+h_bot_cube**2) + (1/12)*m_telescope*(3*r_dst**2 + h_dst)  # pitch
Ixx_combined = (1/12)*m_bottom_cube*(w_tot_cube**2+h_bot_cube**2) + (1/12)*m_telescope*(3*r_dst**2 + h_dst) # roll


Izz_cube = (1/12)*m_total_cube*(w_tot_cube**2+l_tot_cube**2)
Iyy_cube = (1/12)*m_total_cube*(l_tot_cube**2+h_tot_cube**2)
Ixx_cube = (1/12)*m_total_cube*(w_tot_cube**2+h_tot_cube**2)

print(max([Izz_combined, Iyy_combined, Ixx_combined, Izz_cube, Iyy_cube, Ixx_cube]), 'max inertia')
#max_I_diff = max([abs(Izz-Iyy), abs(Izz-Ixx), abs(Iyy-Ixx)])  # largest moment of inertia
# difference
print([abs(Izz_combined-Iyy_combined), abs(Izz_combined-Ixx_combined)], 'combined')
print([abs(Izz_cube-Iyy_cube), abs(Izz_cube-Ixx_cube)], 'cube')
max_I_diff_combined = max([abs(Izz_combined-Iyy_combined), abs(Izz_combined-Ixx_combined)])  # largest moment of inertia
# difference
max_I_diff_cube = max([abs(Izz_cube-Iyy_cube), abs(Izz_cube-Ixx_cube)])   # largest moment of inertia



theta = 30/360 *(2*np.pi)  # 30 degree offset
gravity_dist_torque_comb = 3*mu_earth/(alt_orbit**3) * max_I_diff_combined*np.sin(theta)
gravity_dist_torque_cube = 3*mu_earth/(alt_orbit**3) * max_I_diff_cube*np.sin(theta)
# only for when it is in the incorrect orientation
print(gravity_dist_torque_comb, gravity_dist_torque_cube)
gravity_dist_torque = max([gravity_dist_torque_cube, gravity_dist_torque_comb])
total_torque = gravity_dist_torque + aero_dist_torque + sun_dist_torque + mag_dist_torque
print("Tot_torque: ", total_torque, "gravity_dist_torque: ", gravity_dist_torque,
      "aero_dist_torque: ", aero_dist_torque, "sun_dist_torque: ", sun_dist_torque,
      "mag_dist_torque: ", mag_dist_torque)
print("Tot_torque mNm: ", total_torque*1000)
f_aero = 0.5*rho*a_frontal*c_d_satellite*v_orbit**2
aero_dist_torque = f_aero*largest_arm

mag_field_max = 57587e-9  # max of modelling of very low earth orbit source
mag_dipole = 3.5e-3  # 19690020961 pdf, it's based on amount of mass
mag_dist_torque = mag_dipole*mass_sat*mag_field_max  #

theta = 30/(2*np.pi)
mag_field_max = 57587e-9  # max of modelling of very low earth orbit source
mag_dipole = 3.5e-3  # 19690020961 pdf, it's based on amount of mass
mag_dist_torque = mag_dipole*mass_sat*mag_field_max  #
print(mag_dipole*mass_sat)
print(mag_dist_torque, 'sss')

#gravity_dist_torque = 3*mu_earth/(alt_orbit**2) * max([Izz-Iyy, Ixx-Izz, Ixx-Izz])*theta
nominal_torque = aero_dist_torque + mag_dist_torque + sun_dist_torque + gravity_dist_torque
print(nominal_torque, 'Nominal Torque [Nm]', nominal_torque*1000, '[mNm]')
impulse_moment = nominal_torque*t_life
print(impulse_moment, '[Nms]')