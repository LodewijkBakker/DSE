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
yaw = np.absolute(x)
roll = np.absolute(y)
pitch = np.absolute(z)

# fig, ax = plt.subplots(3)
# ax[1].plot(roll)
# ax[2].plot(pitch)
# ax[0].plot(yaw)
# plt.show()

#torque for roll , pitch, yaw HAS TO BE POSITIVE
entry_torque = np.array([1.06e-5, 8.73e-5, 5.4e-6])
  
#Average Dipole Moment
dip_moment = np.divide(entry_torque, 1e-9 * np.array([np.average(roll), np.average(pitch), np.average(yaw)]))

#avg torque per timestep
magn_torque_roll = list(map(lambda entry: (1e-9 * entry) * dip_moment[0], roll))
magn_torque_pitch = list(map(lambda entry: (1e-9 * entry) * dip_moment[1], pitch))
magn_torque_yaw = list(map(lambda entry: (1e-9 * entry) * dip_moment[2], yaw))


# print(np.average(magn_torque_yaw), 'mty')
# print(np.average(magn_torque_roll))
# print(np.average(magn_torque_pitch))
# print(magn_torque_yaw)
# print(magn_torque_roll)
# print(magn_torque_pitch)

#get given array a , b , c
# new_torque = a - yaw , b-roll, c-pitch  
entry_roll_torque_cont = np.ones(len(x))*entry_torque[0]
entry_pitch_torque_cont = np.ones(len(y))*entry_torque[1]
entry_yaw_torque_cont = np.ones(len(z))*entry_torque[2]

# plt.plot(magn_torque_roll)
# plt.plot(magn_torque_pitch)
# plt.plot(magn_torque_yaw)
# plt.plot(entry_roll_torque_cont)
# plt.plot(entry_pitch_torque_cont)
# plt.plot(entry_yaw_torque_cont)
# plt.show()

result_torque_left_roll = list(map(lambda b, roll: b - roll, entry_roll_torque_cont, magn_torque_roll))
result_torque_left_pitch = list(map(lambda c, pitch: c - pitch, entry_pitch_torque_cont, magn_torque_pitch))
result_torque_left_yaw = list(map(lambda a, yaw: a - yaw, entry_yaw_torque_cont, magn_torque_yaw))

# plt.plot(entry_roll_torque_cont)
# plt.plot(entry_pitch_torque_cont)
# plt.plot(entry_yaw_torque_cont)
# plt.plot(result_torque_left_roll)
# plt.plot(result_torque_left_pitch)
# plt.plot(result_torque_left_yaw)
# plt.show()

#Integration

angular_momentum_roll = integrate.cumtrapz(result_torque_left_roll)
angular_momentum_pitch = integrate.cumtrapz(result_torque_left_pitch)
angular_momentum_yaw = integrate.cumtrapz(result_torque_left_yaw)

# Find The peaks at minimum points under 0

roll_peaks, _ = find_peaks(-angular_momentum_roll, height=0)
pitch_peaks, _ = find_peaks(-angular_momentum_pitch, height=0)
yaw_peaks, _ = find_peaks(-angular_momentum_yaw, height=0)

# plt.plot(angular_momentum_roll)
# plt.plot(roll_peaks, angular_momentum_roll[roll_peaks], "x")
# plt.plot(np.zeros_like(angular_momentum_roll), "--", color="gray")
# plt.show()

# plt.plot(angular_momentum_pitch)
# plt.plot(pitch_peaks, angular_momentum_pitch[pitch_peaks], "x")
# plt.plot(np.zeros_like(angular_momentum_pitch), "--", color="gray")
# plt.show()

# plt.plot(angular_momentum_yaw)
# plt.plot(yaw_peaks, angular_momentum_yaw[yaw_peaks], "x")
# plt.plot(np.zeros_like(angular_momentum_yaw), "--", color="gray")
# plt.show()

# Add the positive slopes

for i in range(0, len(roll_peaks)):  
    # print(roll_peaks[i])
    if angular_momentum_roll[roll_peaks[i]] < 0:
        # print(angular_momentum_roll[roll_peaks[i]])
        angular_momentum_roll[roll_peaks[i]:] -= angular_momentum_roll[roll_peaks[i]]
    # ideally use representative length, 
    # not pitch yaw or something specifically

for i in range(0, len(pitch_peaks)):  
    # print(pitch_peaks[i])
    if angular_momentum_pitch[pitch_peaks[i]] < 0:
        # print(angular_momentum_pitch[pitch_peaks[i]])
        angular_momentum_pitch[pitch_peaks[i]:] -= angular_momentum_pitch[pitch_peaks[i]]
    # ideally use representative length, 
    # not pitch yaw or something specifically

for i in range(0, len(yaw_peaks)):  
    # print(yaw_peaks[i])
    if angular_momentum_yaw[yaw_peaks[i]] < 0:
        # print(angular_momentum_yaw[yaw_peaks[i]])
        angular_momentum_yaw[yaw_peaks[i]:] -= angular_momentum_yaw[yaw_peaks[i]]
    # ideally use representative length, 
    # not pitch yaw or something specifically


#To ensure no negative values 

modified_roll = []
for i in range(len(angular_momentum_roll)):
    if angular_momentum_roll[i] >= 0:
        modified_roll.append(angular_momentum_roll[i])
    else:
        modified_roll.append(0)

modified_pitch = []
for i in range(len(angular_momentum_pitch)):
    if angular_momentum_pitch[i] >= 0:
        modified_pitch.append(angular_momentum_pitch[i])
    else:
        modified_pitch.append(0)

modified_yaw = []
for i in range(len(angular_momentum_yaw)):
    if angular_momentum_yaw[i] >= 0:
        modified_yaw.append(angular_momentum_yaw[i])
    else:
        modified_yaw.append(0)

plt.plot(angular_momentum_roll, label="roll")
plt.plot(modified_roll, label="mod_roll")
leg = plt.legend(loc='upper right')
plt.show()

plt.plot(angular_momentum_pitch, label="pitch")
plt.plot(modified_pitch, label="mod _pitch")
leg = plt.legend(loc='upper right')
plt.show()

plt.plot(angular_momentum_yaw, label="yaw")
plt.plot(modified_yaw, label="mod_yaw")
leg = plt.legend(loc='upper right')
plt.show()

average_test = [np.average(angular_momentum_roll), np.average(angular_momentum_pitch), np.average(angular_momentum_yaw)]


#print(angular_momentum_yaw*1000, angular_momentum_roll*1000, angular_momentum_pitch*1000, 'Angular momentum yaw roll pitch mNms')

# print('average', np.average(x), np.average(y), np.average(z)) #nanotesla
# print('max', np.max(x), np.max(y), np.max(z))
# print('min', np.min(x), np.min(y), np.min(z))
# print(dip_moment)
# print(result_a)
# print(result_b)
# print(result_c)
# print(magn_torque_yaw[0])
# print(arrayA[0])


