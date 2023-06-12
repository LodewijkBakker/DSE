import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

fig, ax = plt.subplots(3)
ax[0].plot(x)
ax[1].plot(y)
ax[2].plot(z)
plt.show()

print('average', np.average(x), np.average(y), np.average(z))
print('max', np.max(x), np.max(y), np.max(z))
print('min', np.min(x), np.min(y), np.min(z))
