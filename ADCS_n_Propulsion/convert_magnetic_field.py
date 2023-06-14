import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.signal import find_peaks

"""
Time step - 30s
N_orbits - 30
"""

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

x = np.interp(np.linspace(0, len(x), num=1000), np.arange(0, len(x)), x)
y = np.interp(np.linspace(0, len(y), num=1000), np.arange(0, len(y)), y)
z = np.interp(np.linspace(0, len(z), num=1000), np.arange(0, len(z)), z)

#Turn arrays to absolute values
mag_field = np.array([np.absolute(x), np.absolute(y), np.absolute(z)])


fig, ax = plt.subplots(3)
ax[1].plot(roll)
ax[2].plot(pitch)
ax[0].plot(yaw)
plt.show()

def avg_torque_calc():
    # torque for roll , pitch, yaw HAS TO BE POSITIVE
    entry_torque = np.array([1.75E-05, 8.22E-05, 2.58E-05])
    # average over 3 orbits so including propulsion
    entry_torque_straight_prop = np.array([1.75E-05, 3.51E-04, 2.58E-05])
    # average over 3 orbits so including propulsion

    t_orbit = 5431  # time of an orbit
    # new average
    # entry_torque[1] = (entry_torque[1]*t_orbit*3 + (entry_torque_straight_prop[1]-entry_torque[1])*900)/(t_orbit*3) #

    avg_torque = np.array([1.75E-05, 8.22E-05, 2.58E-05])
    return avg_torque

def sizing_dipole(mag_field, avg_torque, I_sat):
    """
    sizes the dipole of the magnetic torquer
    :param mag_field: roll, pitch yaw magnetic field in [nanotesla]
    :param avg_torque: roll, pitch, yaw average torque over orbit [Nm]
    :param I_sat: moment of inertia of the satellite along body axis (x, y, z)
    :return dip_moment: dipole moment in [A/m]
    """

    avg_mag_field = 1e-9 * np.average(mag_field, axis=1)  # in [tesla]
    # Nominal operation dip_moment
    dip_moment_nom = np.divide(avg_torque, avg_mag_field)
    # print(dip_moment_nom, 'nominal operation dip moment')

    # detumbling dip_moment required
    detumbling_speed = 300/180*np.pi  # 200 deg/s
    t_orbits_detumbling = 40*5423  # time available for detumbling until deorbit

    detumbling_torque = detumbling_speed/t_orbits_detumbling * I_sat
    dip_moment_tumb = np.divide(detumbling_torque, avg_mag_field)
    # print(dip_moment_tumb, 'if sized on detumbling')

    dip_moment = np.maximum(dip_moment_nom, dip_moment_tumb)
    # print(dip_moment, 'Dip used for calculation')

    return dip_moment

avg_torque = avg_torque_calc()
I_sat = np.array([0.33, 0.42, 0.72])  # moment of inertia of the satellite
dip_moment = sizing_dipole(mag_field, avg_torque, I_sat)


def res_torques_calc(mag_field, avg_torque):
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


def integrate_torques(res_torques):
    """
    :param res_torques: roll, pitch, yaw resultant torques over orbit with desaturation [Nm]
    :return:
    """
    # Integration
    min_to_sec = 30  # half a minute to seconds
    angular_momentum = integrate.cumtrapz(res_torques, axis=1)*min_to_sec  # min to sec

    return angular_momentum

def angular_momentum_realism_creator(angular_momentum):
    """
    :param angular_momentum: stored velocity of spacecraft in reaction wheels [Nms]
    :return:
    """
    # Find The peaks at minimum points under 0
    roll_peaks, _ = find_peaks(-angular_momentum[0], height=0)
    pitch_peaks, _ = find_peaks(-angular_momentum[1], height=0)
    yaw_peaks, _ = find_peaks(-angular_momentum[2], height=0)

    # Add the positive slopes

    for i in range(0, len(roll_peaks)):
        # print(roll_peaks[i])
        if angular_momentum[0][roll_peaks[i]] < 0:
            # print(angular_momentum_roll[roll_peaks[i]])
            angular_momentum[0][roll_peaks[i]:] -= angular_momentum[0][roll_peaks[i]]
        # ideally use representative length,
        # not pitch yaw or something specifically

    for i in range(0, len(pitch_peaks)):
        # print(pitch_peaks[i])
        if angular_momentum[1][pitch_peaks[i]] < 0:
            # print(angular_momentum_pitch[pitch_peaks[i]])
            angular_momentum[1][pitch_peaks[i]:] -= angular_momentum[1][pitch_peaks[i]]
        # ideally use representative length,
        # not pitch yaw or something specifically

    for i in range(0, len(yaw_peaks)):
        # print(yaw_peaks[i])
        if angular_momentum[2][yaw_peaks[i]] < 0:
            # print(angular_momentum_yaw[yaw_peaks[i]])
            angular_momentum[2][yaw_peaks[i]:] -= angular_momentum[2][yaw_peaks[i]]
        # ideally use representative length,
        # not pitch yaw or something specifically

# Find Maximum 
max_roll = np.max(angular_momentum_roll)
max_pitch = np.max(angular_momentum_pitch)
max_yaw = np.max(angular_momentum_yaw)

print(max_roll, max_pitch, max_yaw, "Roll, Pitch, Yaw")

    #To ensure no negative values

    modified_roll = []
    for i in range(len(angular_momentum[0])):
        if angular_momentum[0][i] >= 0:
            modified_roll.append(angular_momentum[0][i])
        else:
            modified_roll.append(0)

    modified_pitch = []
    for i in range(len(angular_momentum[1])):
        if angular_momentum[i][1] >= 0:
            modified_pitch.append(angular_momentum[i][1])
        else:
            modified_pitch.append(0)

    modified_yaw = []
    for i in range(len(angular_momentum[2])):
        if angular_momentum[i][2] >= 0:
            modified_yaw.append(angular_momentum[i][2])
        else:
            modified_yaw.append(0)

    plt.plot(angular_momentum[:][0], label="roll")
    plt.plot(modified_roll, label="mod_roll")
    leg = plt.legend(loc='upper right')
    plt.show()

    plt.plot(angular_momentum[:][1], label="pitch")
    plt.plot(modified_pitch, label="mod _pitch")
    leg = plt.legend(loc='upper right')
    plt.show()

    plt.plot(angular_momentum[:][2], label="yaw")
    plt.plot(modified_yaw, label="mod_yaw")
    leg = plt.legend(loc='upper right')
    plt.show()

    average_test = np.average(angular_momentum, axis=1)

    safety_factor = 2
    print(np.max(angular_momentum_roll)*1000*safety_factor, np.max(angular_momentum_pitch)*1000*safety_factor, np.max(angular_momentum_yaw)*1000*safety_factor, 'Angular momentum roll pitch yaw mNms, with safety factor 2')

# print('average', np.average(x), np.average(y), np.average(z)) #nanotesla
# print('max', np.max(x), np.max(y), np.max(z))
# print('min', np.min(x), np.min(y), np.min(z))
# print(dip_moment)
# print(result_a)
# print(result_b)
# print(result_c)
# print(magn_torque_yaw[0])
# print(arrayA[0])


