import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


F = np.array(pd.read_csv('intensity.csv').columns.values).astype(float)
I = np.array(pd.read_csv('inclination.csv').columns.values).astype(float) * np.pi / 180
D = np.array(pd.read_csv('declination.csv').columns.values).astype(float) * np.pi / 180
FID = np.array([F, I, D]).T
LLA = np.array(pd.read_csv('lla.csv').values)   # latitide, longitude, altitude

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

#Turn arrays to absolute values
yaw = np.absolute(x)
roll = np.absolute(y)
pitch = np.absolute(z)

fig, ax = plt.subplots(3)
ax[0].plot(yaw)
ax[1].plot(roll)
ax[2].plot(pitch)
plt.show()

#torque for roll , pitch, yaw
torque = np.array([1.06e-5, 8.73e-5, 5.4e-6])


#Average Dipole Moment
dip_moment = np.divide(torque, 1e-9 * np.array([np.average(roll), np.average(pitch), np.average(yaw)]))

#avg torque per timestep
avg_torque_yaw = list(map(lambda entry: entry * dip_moment[2], x))
avg_torque_roll = list(map(lambda entry: entry * dip_moment[0], y))
avg_torque_pitch = list(map(lambda entry: entry * dip_moment[1], z))

# print(avg_torque_yaw)
# print(avg_torque_roll)
# print(avg_torque_pitch)

#get given array a , b , c
# new_torque = a - yaw , b-roll, c-pitch

result_a = list(map(lambda a, avg_torque_yaw: a - avg_torque_yaw, arrayA, avg_torque_yaw))
result_b = list(map(lambda b, avg_torque_roll: b - avg_torque_roll, arrayB, avg_torque_roll))
result_c = list(map(lambda c, avg_torque_pitch: c - avg_torque_pitch, arrayC, avg_torque_pitch))


# print('average', np.average(x), np.average(y), np.average(z)) #nanotesla
# print('max', np.max(x), np.max(y), np.max(z))
# print('min', np.min(x), np.min(y), np.min(z))
print(dip_moment)
# print(new_torque)


#nano tesla absolute avg $
#avg torque per axis + nano tesla $
#Dipole moment = t / n $

# dipole moment * absolute nanotesla for every axis at every time step $

# new different torque given single value for each axis $
# per timestep ->  torque given - torque possible magnetorquer $

# lijst van times ?
# lijst van torque given - torque possible ?
# simpson integrtor ?


