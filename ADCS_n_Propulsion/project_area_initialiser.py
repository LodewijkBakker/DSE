import multiprocessing as mp

import ADCS_n_Propulsion.project_area_backend as pab
from stl import mesh
import os
import numpy as np
from rot_angles_generators import sphere_fib_user_angle, sphere_fib_octant


def torque_calculator(filename):
    #busy_text = Label(window, text="Don't worry the program is busy", width="50")
    #busy_text.grid(column=0, row=8, columnspan="2")

    # TODO MAKE CLASS OR TOML FILE OUT OF THIS!!!!
    # general parameters:
    q_sun_max = 1414
    a_albedo = 0.19
    q_total_max = q_sun_max + a_albedo*q_sun_max
    c_light = 299792458

    mu_earth = 3.986004418e14
    rad_earth = 6371800
    B_orbit_min_orbit_avg = [50000, 0, 0]  # sizing magnetorquers on this represents the minimum size they have to be
    B_orbit_abs_min = [0, 0, 0]  # this is relative to a earth with north and south as axis
    B_orbit_max = [500000, 0, 0]

    # satellite values, take out to toml after
    alt_sat = 300000
    total_mass = 54
    c_d = 2.027
    thrust = 2.5e-3  # no thrust when tumbling
    v_orbit = 5000  # !!!!! change
    rho = 0.00001  # !!!!!! change
    a_reflectivity = 0.94  # jose thermal
    #incl = 30

    alt_orbit = rad_earth + alt_sat

    m_telescope = 30
    m_total_cube = 54
    m_bottom_cube = m_total_cube - m_telescope

    A_kg_dipole = 3.5e-3
    A_dipole = A_kg_dipole * m_total_cube

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

    Izz = Izz_combined
    Iyy = Iyy_combined
    Ixx = Ixx_combined
    ###########

    stl_path = filename
    stl_model_main = mesh.Mesh.from_file(stl_path)

    list_grav_centers = [[0, 0, 0], [1, 0, 0]]

    # gravity gradient
    gravity_t_constant = 3 / 2 * mu_earth / (alt_orbit ** 3)

    ng = 5000
    # store aero, gravity, sun, magnetic torque,
    t_store_extreme = np.zeros((ng, len(list_grav_centers), 4, 3))

    # Extreme case:
    # creating fibonacci angles
    rot_angles = sphere_fib_octant(ng, [True, True, True, True, True, True, True, True])
    for i, rot_angle in enumerate(rot_angles):
        res = pab._calculate_projected_area(rot_angle, stl_model_main, list_grav_centers)
        for j, single_grav_arm_area in enumerate(res[1]): # single grav arm is a area * (arm vector cross trajectory)
            t_aero = 0.5 * rho * c_d * v_orbit ** 2 * single_grav_arm_area
            t_store_extreme[i, j, 0] = t_aero

            # should be checked again this bottom one, again moment in xyz
            gravity_dist_torque = [gravity_t_constant * (Izz - Iyy) * np.sin(-2 * rot_angles[2]),
                                   gravity_t_constant * (Izz - Ixx) * np.sin(-2 * rot_angles[1]), 0]

            t_store_extreme[i, j, 1] = gravity_dist_torque

            t_sun = q_sun_max/c_light*(1+a_reflectivity)*single_grav_arm_area
            t_store_extreme[i, j, 2] = t_sun

            # Magnetic torque given a certain inclination
            t_store_extreme[i, j, 3] = [np.nan, np.nan, np.nan]

    # Nominal case:
    # creating rot angles for nominal case:
    rot_angles = sphere_fib_user_angle(1)
    ng_user_angle = len(rot_angles)  # set by user
    # store aero, gravity, sun, magnetic torque,
    t_store_nominal = np.zeros((ng_user_angle, len(list_grav_centers), 4, 3))

    for i, rot_angle in enumerate(rot_angles):
        res = pab._calculate_projected_area(rot_angle, stl_model_main, list_grav_centers)
        for j, single_grav_arm_area in enumerate(res[1]):  # single grav arm is a area * (arm vector cross trajectory)
            # Aerodynamic Torque
            t_aero = 0.5 * rho * c_d * v_orbit**2 * single_grav_arm_area
            t_store_nominal[i, j, 0] = t_aero

            # should be checked again this bottom one, again moment in xyz
            gravity_dist_torque = [gravity_t_constant * (Izz - Iyy) * np.sin(-2*rot_angles[2]),
                                   gravity_t_constant * (Izz - Ixx) * np.sin(-2*rot_angles[1]), 0]
            t_store_nominal[i, j, 1] = gravity_dist_torque

            # Sun torque
            # initialise as nan, in nominal most extreme case is assumed and that is calculated above
            t_store_nominal[i, j, 2] = [np.nan, np.nan, np.nan]

            # Magnetic torque given a certain inclination
            t_store_nominal[i, j, 3] = [np.nan, np.nan, np.nan]

    return t_store_extreme, t_store_nominal, rot_angles


def adjust_torque_direction(torque_vector, rotation_vector):
    adjusted_torque = np.copy(torque_vector)
    rotation_signs = np.sign(rotation_vector)

    for i in range(1, len(rotation_signs)):  # starting at 1 such to not remove z in any case
        print(2-i)
        if rotation_signs[i] == 1:
            adjusted_torque[2-i] = min(0, adjusted_torque[2-i])  # Set torque to 0 if rotation is positive
        elif rotation_signs[i] == -1:
            adjusted_torque[2-i] = max(0, adjusted_torque[2-i])  # Set torque to 0 if rotation is negative

    return adjusted_torque


def res_parser(t_store_extreme, t_store_nominal, user_rot_angles):
    # result max torque per disturbance per gravitie center per axis, nominal, also give absolute max
    # result max torque per disturbance per gravitie center per axis, extreme, also give absolute max
    # result max torque per gravity center per axis w/o magnetic, nominal (sun torque is used from max torque extreme case)
    # result max torque per gravity center per axis w/o magnetic, extreme
    # result avg torque per gravity center per axis w/o magnetic, extreme

    ## result max torque grav + aero per axis, nominal
    ## result max torque grav + aero per axis, extreme
    n_grav_centers = np.shape(t_store_extreme)[1]

    # --------------------------------------------
    # ---------Extreme----------------

    # result max torque per disturbance per centre of gravity per axis, extreme
    max_torque_per_disturbance_extreme = []
    for i in range(n_grav_centers):
        max_torque_per_disturbance_f_grav_center_extreme = [[np.max(t_store_extreme[:, :, 0, 0]), np.max(t_store_extreme[:, :, 0, 1]), np.max(t_store_extreme[:, :, 0, 2])],
                                                            [np.max(t_store_extreme[:, :, 1, 0]), np.max(t_store_extreme[:, :, 1, 1]), np.max(t_store_extreme[:, :, 1, 2])],
                                                            [np.max(t_store_extreme[:, :, 2, 0]), np.max(t_store_extreme[:, :, 2, 1]), np.max(t_store_extreme[:, :, 2, 2])],
                                                            [np.max(t_store_extreme[:, :, 3, 0]), np.max(t_store_extreme[:, :, 3, 1]), np.max(t_store_extreme[:, :, 3, 2])]]
        max_torque_per_disturbance_extreme.append(max_torque_per_disturbance_f_grav_center_extreme)

    # result max torque per gravity center per axis w/o magnetic, extreme
    max_torque_total_extreme_wo_mag = []  # without magnetic
    # result avg torque per gravity center per axis w/o magnetic, extreme
    avg_torque_extreme_wo_mag = []  # without magnetic
    for i in range(n_grav_centers):
        # sun always in worst position
        max_torque_total_extreme_wo_mag.append([np.max(np.sum(t_store_extreme[:, i, :2, 0], axis=2)) + max_torque_per_disturbance_extreme[i][2][0],
                                                np.max(np.sum(t_store_extreme[:, i, :2, 1], axis=2)) + max_torque_per_disturbance_extreme[i][2][1],
                                                np.max(np.sum(t_store_extreme[:, i, :2, 2], axis=2)) + max_torque_per_disturbance_extreme[i][2][2]])
        # sun always in worst position
        avg_torque_extreme_wo_mag.append([np.average(np.sum(t_store_extreme[:, i, :2, 0], axis=2)) + max_torque_per_disturbance_extreme[i][2][0],
                                          np.average(np.sum(t_store_extreme[:, i, :2, 1], axis=2)) + max_torque_per_disturbance_extreme[i][2][1],
                                          np.average(np.sum(t_store_extreme[:, i, :2, 2], axis=2)) + max_torque_per_disturbance_extreme[i][2][2]])

    # --------------------------------------------
    # ---------Nominal----------------

    # result max torque per disturbance for all centre of gravities per axis, nominal
    max_torque_per_disturbance_nominal = []

    # result max torque per gravity center per axis w/o magnetic, nominal (sun torque is used from max torque extreme case)
    max_torque_total_nominal_wo_mag = []

    for i in range(n_grav_centers):
        max_torque_per_disturbance_f_grav_center_nominal = [[np.max(t_store_nominal[:, :, 0, 0]), np.max(t_store_nominal[:, :, 0, 1]), np.max(t_store_nominal[:, :, 0, 2])],
                                                            [np.max(t_store_nominal[:, :, 1, 0]), np.max(t_store_nominal[:, :, 1, 1]), np.max(t_store_nominal[:, :, 1, 2])],
                                                            [np.max(t_store_nominal[:, :, 2, 0]), np.max(t_store_nominal[:, :, 2, 1]), np.max(t_store_nominal[:, :, 2, 2])],
                                                            [np.max(t_store_nominal[:, :, 3, 0]), np.max(t_store_nominal[:, :, 3, 1]), np.max(t_store_nominal[:, :, 3, 2])]]
        max_torque_per_disturbance_nominal.append(max_torque_per_disturbance_f_grav_center_nominal)

        # sun always in worst position, nominal doesn't really bring sun with it
        max_torque_total_nominal_wo_mag.append([np.max(np.sum(t_store_nominal[:, i, :2, 0], axis=2)) + max_torque_per_disturbance_extreme[i][2][0],
                                                np.max(np.sum(t_store_nominal[:, i, :2, 1], axis=2)) + max_torque_per_disturbance_extreme[i][2][1],
                                                np.max(np.sum(t_store_nominal[:, i, :2, 2], axis=2)) + max_torque_per_disturbance_extreme[i][2][2]])

    # ---------------------------
    # --------Unstable Nominal -----------
    # unstable torque per gravity center per axis w/o magnetic nominal (sun torque is used from max torque extreme case)
    unstable_rot_angles = []
    unstable_torques_total_nominal = []

    for i, rot_angle in enumerate(user_rot_angles):
        for j in range(n_grav_centers):
            total_t_vec = [np.sum(t_store_nominal[i, j, :2, 0]) + max_torque_per_disturbance_extreme[i][2][0],
                           np.sum(t_store_nominal[i, j, :2, 1]) + max_torque_per_disturbance_extreme[i][2][1],
                           np.sum(t_store_nominal[i, j, :2, 2]) + max_torque_per_disturbance_extreme[i][2][2]]
            adjust_total_t_vec = adjust_torque_direction(total_t_vec, rot_angle) # adjust for pitch and yaw stability
            if np.linalg.norm(adjust_total_t_vec) > 0:
                unstable_rot_angles.append(rot_angle)
                unstable_torques_total_nominal.append(adjust_total_t_vec)


    return max_torque_per_disturbance_extreme, max_torque_total_extreme_wo_mag, avg_torque_extreme_wo_mag, \
        max_torque_per_disturbance_nominal, max_torque_total_nominal_wo_mag


def unit_test_parser():
    z = [[[[0., 0., 0.],
           [0., 7., 0.],
           [0., 0., 0.],
           [4., 4., 4.]],

          [[0., 10., 0.],
           [0., 0., 0.],
           [5., 0., 0.],
           [0., 0., 0.]]],

         [[[10., 0., 0.],
           [7., 0., 0.],
           [0., 5., 0.],
           [0., 0., 0.]],

          [[0., 0., 10.],
           [0., 6., 7.],
           [0., 0., 5.],
           [0., 0., 0.]]]]

    z_max_p_disturbance = [[10, 10, 10],
                           [7, 7, 7],
                           [4, 4, 4],
                           [5, 5, 5]]

    print(res_parser(z, False))

def unit_tests():
    stl_path = os.getcwd() + '\cube_ascii_centered.stl'  # 2 by 2 by 2 cube
    stl_model_main_centered = mesh.Mesh.from_file(stl_path)

    stl_path = os.getcwd() + '\cube_ascii_offset.stl'  # 2 by 2 by 2 cube
    stl_model_main_offset = mesh.Mesh.from_file(stl_path)

    list_grav_centers = [[0, 0, 0], [1, 0, 0]]
    rot_angles = [[0, 0, 0],
                  [45, 0, 0],
                  [0, 45, 0],
                  [0, 0, 45],
                  [90, 0, 0],
                  [0, 90, 0],
                  [21, 32, 2],
                  [5, -30, 15]]

    for rot_angle in rot_angles:
        res = pab._calculate_projected_area(rot_angle, stl_model_main_centered, list_grav_centers)
        print(res, 'centered')
        res = pab._calculate_projected_area(rot_angle, stl_model_main_offset, list_grav_centers)
        print(res, 'offset')

if __name__ == "__main__":
    unit_tests()

# Function for opening the
# file explorer window