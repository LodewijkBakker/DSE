import multiprocessing as mp

import ADCS_n_Propulsion.project_area_backend as pab
from stl import mesh
import numpy as np

def torque_calculator(filename, rot_angles):
    #busy_text = Label(window, text="Don't worry the program is busy", width="50")
    #busy_text.grid(column=0, row=8, columnspan="2")

    stl_path = filename
    stl_model_main = mesh.Mesh.from_file(stl_path)

    list_grav_centers = [[0, 0, 0], [1, 0, 0]]

    # general parameters:
    mu_earth = 3.986004418e14
    rad_earth = 6371800

    # satellite values, take out to toml after
    alt_sat = 300000
    total_mass = 54
    c_d = 2.027
    thrust = 2.5e-3
    v_orbit = 5000  # !!!!! change
    rho = 0.00001  # !!!!!! change
    #incl = 30

    alt_orbit = rad_earth + alt_sat

    m_telescope = 30
    m_total_cube = 54
    m_bottom_cube = m_total_cube - m_telescope

    h_bot_cube = 0.22
    h_tot_cube = 0.34
    w_tot_cube = 0.36
    l_tot_cube = 0.35
    r_dst = (0.368 + 0.35) / 4  # check
    h_dst = 0.33  # ???

    Izz_combined = (1 / 12) * m_bottom_cube * (w_tot_cube ** 2 + l_tot_cube ** 2) + (
                1 / 2) * m_telescope * r_dst ** 2  # yaw
    Iyy_combined = (1 / 12) * m_bottom_cube * (l_tot_cube ** 2 + h_bot_cube ** 2) + (
                1 / 12) * m_telescope * (3 * r_dst ** 2 + h_dst)  # pitch
    Ixx_combined = (1 / 12) * m_bottom_cube * (w_tot_cube ** 2 + h_bot_cube ** 2) + (
                1 / 12) * m_telescope * (3 * r_dst ** 2 + h_dst)  # roll

    Izz_cube = (1 / 12) * m_total_cube * (w_tot_cube ** 2 + l_tot_cube ** 2)
    Iyy_cube = (1 / 12) * m_total_cube * (l_tot_cube ** 2 + h_tot_cube ** 2)
    Ixx_cube = (1 / 12) * m_total_cube * (w_tot_cube ** 2 + h_tot_cube ** 2)

    ###########

    # gravity torque
    # theta = 30 / 360 * (2 * np.pi)  # 30 degree offset
    # gravity_dist_torque_comb = 3 * mu_earth / (alt_orbit ** 3) * max_I_diff_combined * np.sin(
    #     theta)
    # gravity_dist_torque_cube = 3 * mu_earth / (alt_orbit ** 3) * max_I_diff_cube * np.sin(
    #     theta)

    # Pressure torques
    for rot_angle in rot_angles:
        res = pab._calculate_projected_area(rot_angle, stl_model_main, list_grav_centers)
        for single_grav_arm_area in res: # single grav arm is a area * (arm vector cross trajectory)
            # Aerodynamic Torque
            t_aero = 0.5 * rho * c_d * v_orbit**2 * single_grav_arm_area
            # Solar torque

    #for rot_angle in rot_angles:
        # Sum torques per rot angles
        #t_total =

        # Check instability for all the torques summed, except for the roll angle!!

    # busy_text.config(text="Thanks for your patience it's done!")

    return res

# Function for opening the
# file explorer window