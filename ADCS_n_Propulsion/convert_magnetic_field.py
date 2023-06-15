import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.signal import find_peaks
from tqdm import tqdm
from time import perf_counter


def mag_field_creator():
    """
    :return mag_field: magnetic field in body axis [tesla]
    :return time_step: time step of magnetic field (resolution)
    """

    """
    Time step - 30s
    N_orbits - 30
    """
    time_step = 30

    F = np.array(pd.read_csv('mag_res/intensity.csv').columns.values).astype(float)
    I = np.array(pd.read_csv('mag_res/inclination.csv').columns.values).astype(float) * np.pi / 180
    D = np.array(pd.read_csv('mag_res/declination.csv').columns.values).astype(float) * np.pi / 180
    FID = np.array([F, I, D]).T
    LLA = np.array(pd.read_csv('mag_res/lla.csv').values)   # latitide, longitude, altitude

    BF = []
    FID_cart = []

    def cartesian_convert(arr):   # intensity, inclination, declination
        return [arr[0] * np.sin(arr[2]) * np.cos(arr[1]),
                arr[0] * np.sin(arr[2]) * np.sin(arr[1]),
                arr[0] * np.cos(arr[2])]

    for arr in FID:
        FID_cart.append(cartesian_convert(arr))

    FID_cart = np.array(FID_cart)

    def matrix(long, lat):
        m = np.array([[-np.cos(long) * np.sin(lat), -np.sin(long) * np.sin(lat), np.cos(lat)],
                      [-np.sin(long), np.cos(long), 0],
                      [-np.cos(long) * np.cos(lat), -np.sin(long) * np.cos(lat), -np.sin(lat)]])
        return m

    for i, arr in enumerate(FID_cart[1:, :]):
        BF.append(list(matrix(LLA[i, 1] * np.pi / 180, LLA[i, 0] * np.pi/180) @ np.array([arr[0], arr[1], arr[2]])))

    BF = np.array(BF)
    x = BF[:, 0]
    y = BF[:, 1]
    z = BF[:, 2]

    yaw = np.interp(np.linspace(0, len(x), num=1000), np.arange(0, len(x)), x)
    roll = np.interp(np.linspace(0, len(y), num=1000), np.arange(0, len(y)), y)
    pitch = np.interp(np.linspace(0, len(z), num=1000), np.arange(0, len(z)), z)

    #Turn arrays to absolute values
    mag_field = 1e-9*np.column_stack((np.absolute(roll), np.absolute(pitch), np.absolute(yaw)))

    # plt.plot(np.absolute(x[:200]))
    # plt.show()
    # plt.plot(np.absolute(y[:200]))
    # plt.show()
    # plt.plot(np.absolute(z[:200]))
    # plt.show()
    # print("avg_mag_field", np.average(mag_field, axis=0))  # in [tesla])
    # plt.plot(mag_field[:, 0])
    # plt.show()
    # plt.plot(mag_field[:, 1])
    # plt.show()
    # plt.plot(mag_field[:, 2])
    # plt.show()
    return mag_field, time_step


def avg_torque_calc(nom_torque: np.ndarray, prop_torque: np.ndarray, t_prop_on: float, t_repeat_prop: float):

    # new average
    avg_torque_w_p = (nom_torque*t_repeat_prop + (prop_torque-nom_torque)*t_prop_on)/t_repeat_prop
    avg_torque = np.maximum(nom_torque, avg_torque_w_p)  # should still be high even if t propulsion brings it down
    return avg_torque

def sizing_dipole(mag_field: np.ndarray, avg_torque: np.ndarray, I_sat: np.ndarray, t_orbits_detumbling: float):
    """
    sizes the dipole of the magnetic torquer
    :param mag_field: roll, pitch yaw magnetic field in [tesla]
    :param avg_torque: roll, pitch, yaw average torque over orbit [Nm]
    :param I_sat: moment of inertia of the satellite along body axis (x, y, z)
    :return dip_moment: dipole moment in [A/m]
    """

    avg_mag_field = np.average(mag_field, axis=0)  # in [tesla]
    # Nominal operation dip_moment
    dip_moment_nom = np.divide(avg_torque, avg_mag_field)
    # print(dip_moment_nom, 'nominal operation dip moment')

    # detumbling dip_moment required
    detumbling_speed = 300/180*np.pi  # 200 deg/s

    detumbling_torque = detumbling_speed/t_orbits_detumbling * I_sat
    dip_moment_tumb = np.divide(detumbling_torque, avg_mag_field)

    dip_moment = np.maximum(dip_moment_nom, dip_moment_tumb)
    # print(dip_moment, 'Dip used for calculation')

    return dip_moment

def angular_momentum_calc(mag_field: np.ndarray, avg_torque: np.ndarray, dip_moment: np.ndarray, time_step: float):
    #@jit(nopython=True)
    def res_torques_calc(mag_field: np.ndarray, avg_torque: np.ndarray, dip_moment: np.ndarray):
        """
        :param mag_field: roll, pitch yaw magnetic field in [nanotesla]
        :param avg_torque: roll, pitch, yaw average torque over orbit [Nm]
        :return res_torques: roll, pitch, yaw resultant torques over orbit with desaturation [Nm]
        """

        mag_torque = np.multiply(dip_moment, mag_field)
        res_torques = avg_torque - mag_torque
        # # Add window of propulsion at one point
        # result_torque_left_pitch[400:415] += (entry_torque_straight_prop[1]-entry_torque[1])
        return res_torques

    #@jit(nopython=True)
    def integrate_torques(res_torques: np.ndarray, time_step: float):
        """
        :param res_torques: roll, pitch, yaw resultant torques over orbit with desaturation [Nm]
        :param time_step: time step resolution of the torques
        :return: angular_momentum stored velocity of spacecraft in reaction wheels [Nms]
        """
        # Integration
        angular_momentum = integrate.cumtrapz(res_torques, axis=0)*time_step  # min to sec

        return angular_momentum  # this should be 0 on average I think


    def angular_momentum_realism_creator(angular_momentum: np.ndarray):
        """
        :param angular_momentum: stored velocity of spacecraft in reaction wheels [Nms]
        :return angular_momentum: stored velocity of spacecraft in reaction wheels but now realistic [Nms]:
        """
        # Find The peaks at minimum points under 0
        roll_peaks, _ = find_peaks(-angular_momentum[:, 0], height=0)
        pitch_peaks, _ = find_peaks(-angular_momentum[:, 1], height=0)
        yaw_peaks, _ = find_peaks(-angular_momentum[:, 2], height=0)

        # Update peaks
        for r_peak in roll_peaks:
            if angular_momentum[r_peak, 0] < 0:
                angular_momentum[r_peak:, 0] -= angular_momentum[r_peak, 0]

        for p_peak in pitch_peaks:
            if angular_momentum[p_peak, 1] < 0:
                angular_momentum[p_peak:, 1] -= angular_momentum[p_peak, 1]

        for y_peak in yaw_peaks:
            if angular_momentum[y_peak, 2] < 0:
                angular_momentum[y_peak:, 2] -= angular_momentum[y_peak, 2]

        # all items that are below zero then reset to 0
        angular_momentum = np.clip(angular_momentum, 0, None)
        return angular_momentum

    #@jit(nopython=True)
    def get_sizing_from_angular_momentum(angular_momentum: np.ndarray):
        """
        :param angular_momentum: stored velocity of spacecraft in reaction wheels [Nms]
        :return design_max_angular_momentum_mNms: maximum momentum in mNms to size cmgs
        """

        # Find Maximum
        max_angular_momentum = np.max(angular_momentum, axis=0)
        safety_factor = 2
        design_max_angular_momentum_Nms = max_angular_momentum*safety_factor

        return design_max_angular_momentum_Nms

    i_t_t1 = perf_counter()
    res_torques = res_torques_calc(mag_field, avg_torque, dip_moment)
    i_t_t2 = perf_counter()
    angular_momentum = integrate_torques(res_torques, time_step)
    i_t_t3 = perf_counter()
    angular_momentum = angular_momentum_realism_creator(angular_momentum)
    i_t_t4 = perf_counter()
    design_max_angular_momentum_Nms = get_sizing_from_angular_momentum(angular_momentum)
    i_t_t5 = perf_counter()
    #print(i_t_t5 - i_t_t4, i_t_t4 - i_t_t3, i_t_t3 - i_t_t2, i_t_t2 - i_t_t1)
    return design_max_angular_momentum_Nms


def sizing_angular_momentum_calc(max_angular_momentum):
    # https://arc.aiaa.org/doi/10.2514/1.G006461
    beta_angle = 2*np.arctan(max_angular_momentum[2]/(2*max(max_angular_momentum[0], max_angular_momentum[1])))  # look out for arctan2 values
    sizing_angular_momentum = max(max_angular_momentum[0], max_angular_momentum[1])/(2*(1+np.cos(beta_angle)))
    return sizing_angular_momentum

def sizing_cmg(max_angular_momentum, r_wheel=0.02):
    """
    :param max_angular_momentum: angular momentum along body axis (roll, pitch, yaw)
    :param r_wheel: radius of the wheel
    :return:
    """
    sizing_angular_momentum = sizing_angular_momentum_calc(max_angular_momentum)
    #print(sizing_angular_momentum*1000)  # varies from okay to really really shit so look out
    # TODO investigate 4 placement
    rho_cmg = 0.8989617244  # [u / kg]
    delta_angular_velocity = 2500*2*np.pi/60  # from statistics
    #r_range = 0.02  # actually 0.015 - 0.025
    inertia_needed = sizing_angular_momentum/delta_angular_velocity
    n_cmg = 4
    m_disk = n_cmg * inertia_needed/r_wheel**2
    m_full = n_cmg * (1+2/3) * m_disk
    v_cmg = n_cmg * m_full*rho_cmg  # [u]

    return m_full, v_cmg

def sizing_magnetorquer(dip_moment):
    m_magnetorquer = 0
    v_magnetorquer = 0
    for dip_moment_axis in dip_moment:
        assert(dip_moment_axis > 0)
        m_magnetorquer += dip_moment_axis*0.15
        v_magnetorquer += m_magnetorquer*2.1  # U/kg

    return m_magnetorquer, v_magnetorquer


def optimum_sizer(dip_moment_orig, avg_torque, mag_field, time_step):
    r_range = np.linspace(0.015, 0.025, 5)
    dip_moment_roll = np.linspace(dip_moment_orig[0], dip_moment_orig[0]+15, 10)
    dip_moment_pitch = np.linspace(dip_moment_orig[1], dip_moment_orig[0]+15, 10)
    dip_moment_yaw = np.linspace(dip_moment_orig[2], dip_moment_orig[0]+15, 10)

    dip_min = dip_moment_orig
    m_total_min = [np.inf, np.inf, np.inf]  # total mass, mass cmg, mass magnetorquer
    v_total_min = [np.inf, np.inf, np.inf]  # total volume, volume cmg, volume magnetorquer
    for r in r_range:
        m_total_min = [np.inf, np.inf, np.inf]  # total mass, mass cmg, mass magnetorquer
        v_total_min = [np.inf, np.inf, np.inf]  # total volume, volume cmg, volume magnetorquer
        for dip_moment_roll_i in tqdm(dip_moment_roll, disable=False):
            for dip_moment_pitch_j in dip_moment_pitch:
                for dip_moment_yaw_k in dip_moment_yaw:
                    dip_moment = np.array([dip_moment_roll_i, dip_moment_pitch_j, dip_moment_yaw_k])
                    design_max_angular_momentum_Nms = angular_momentum_calc(mag_field, avg_torque, dip_moment, time_step)
                    #print(design_max_angular_momentum_Nms*1000)  # Fuck pitch is high!!
                    m1, v1 = sizing_cmg(design_max_angular_momentum_Nms, r)
                    m2, v2 = sizing_magnetorquer([dip_moment_roll_i, dip_moment_pitch_j, dip_moment_yaw_k])
                    m_total = m1 + m2
                    v_total = v1 + v2
                    if m_total_min[0] > m_total:
                        m_total_min = [m_total, m1, m2]
                        v_total_min = [v_total, v1, v2]
                        dip_min = dip_moment

        print(dip_min, m_total_min, v_total_min, r)

if __name__ == "__main__":
    # Avg torque calc tester
    assert (np.all(np.isclose(avg_torque_calc(np.array([0, 0, 0]), np.array([0, 0, 0]), 0, 1000), np.array([0, 0, 0]))))
    assert (np.all(np.isclose(avg_torque_calc(np.array([1, 1, 1]), np.array([0.5, 0, 0.2]), 200, 1000), np.array([1, 1, 1]))))
    assert (np.all(np.isclose(avg_torque_calc(np.array([1, 1, 1]), np.array([2, 0, 2]), 1000, 1000), np.array([2, 1, 2]))))

    # Dipole sizing tester
    t_orbits_detumbling = 40*5423  # time available for detumbling until deorbit
    assert np.all(np.isclose(sizing_dipole(np.array([[0, 0, 0], [2, 2, 2]]), np.array([1, 1, 1]), np.array([1, 1, 1]),
                                           t_orbits_detumbling), np.array([1, 1, 1])))
    m = 0.5  # tesla
    v_r = 300/180*np.pi  # rad /s
    torque_required = v_r/(t_orbits_detumbling*m)
    assert np.all(np.isclose(sizing_dipole(np.array([[0, 0, 0], [1, 1, 1]]), np.array([0, 0, 0]), np.array([1, 1, 1]),
                                           t_orbits_detumbling), np.array([torque_required, torque_required,
                                                                           torque_required])))

    # Cmg h tester
    H_roll = 83.8
    H_pitch = 68.1
    H_yaw = 38.4
    H = [H_roll, H_pitch, H_yaw]
    np.testing.assert_almost_equal(sizing_angular_momentum_calc(H), 22, decimal=1)

    # angular momentum is not correct now it seems. (way to high (cause maybe average magnetic))

    mag_field, time_step = mag_field_creator()
    t_orbit = 5432
    t_repeat_prop = t_orbit*3
    t_prop_on = 0
    nom_torque_1 = np.array([1.71E-05, 7.91E-05, 2.58E-05])
    prop_torque_1 = np.array([1.75E-05,	8.22E-05, 2.30E-05])
    avg_torque = avg_torque_calc(nom_torque_1, prop_torque_1, t_prop_on, t_repeat_prop)
    I_sat = np.array([0.33, 0.42, 0.72])  # moment of inertia of the satellite
    dip_moment_1 = sizing_dipole(mag_field, avg_torque, I_sat, t_orbits_detumbling)
    print(dip_moment_1, 'dip start')
    optimum_sizer(dip_moment_1, avg_torque, mag_field, time_step)

    # r_t_1 = res_torques_calc(mag_field, avg_torque, dip_moment_1)
    # int_torques = integrate_torques(r_t_1, time_step)
    # angular_momentum_1 = angular_momentum_realism_creator(int_torques)
    # print(get_sizing_from_angular_momentum(angular_momentum_1))

    #integrate_torques