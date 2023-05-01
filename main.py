

# Main satellite file to update all other values
vol_sat = 2e-3  # [m^3]
mass_sat = vol_sat*1000  # [kg], water method

# First order delta V calculation
a_frontal = 1e-3  # [m^2] frontal area
alt_orbit = (6371+300)*1e3  # [m] nominal orbit altitude from center of earth
mu_earth = 398600.435507*1e3  # [km^3/s^-2] https://ssd.jpl.nasa.gov/astro_par.html
v_orbit = (alt_orbit/mu_earth)**0.5
t_life = 5*365*24*60*60  # [s] nominal life time, 5 years
c_d_satellite = 3  # flat surfaces
rho = 1e-5  # [kg/m^3]

f_aero = rho*a_frontal*c_d_satellite*v_orbit**2  # [n] aerodynamic force at nominal orbit
# https://drive.google.com/drive/u/0/folders/1XzHVbNNy5GXSy6uyikV365bol00h7116
# assuming mass is constant, frontal area does not change and no change in altitude or other density factors
safety_factor = 1.5
delta_v = f_aero / mass_sat * t_life * safety_factor # [m/s] delta v

