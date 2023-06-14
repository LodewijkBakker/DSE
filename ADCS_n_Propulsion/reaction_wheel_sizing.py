import numpy as np
import scipy.optimize as sp

angular_momentum = 5/1000
RPM = 8000
w = RPM * 2 * np.pi / 60
rho = 2700

R = np.linspace(0.01, 0.02, 1000)
T = np.linspace(0.005, 0.04, 1000)

def Vol(R, t):
    return np.pi * R**2 * t

v = []
V = []
for r in R:
    for t in T:
        Ixx = 0.5 * Vol(r,t) * rho * r**2
        if Ixx * w >= angular_momentum and Ixx * w <= 5.5/1000:
                v.append([Vol(r, t), r, t])
                V.append(Vol(r, t))


for x in v:
    if x[0] == min(V):
        print(x)
