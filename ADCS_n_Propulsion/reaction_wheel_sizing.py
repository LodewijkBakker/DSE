import numpy as np


angular_momentum = 5/1000
RPM = 21
w = RPM * 2 * np.pi / 60
rho = 7800

R = np.linspace(0.04, 0.05, 100)
T = np.linspace(0.02, 0.1, 1000)

def Vol(r, t):
    return np.pi * r**2 * t

def mass(r, t):
    return rho * Vol(r, t)
c = 0
v = []
V = []
for r in R:
    for t in T:
        Ixx = 0.5 * Vol(r,t) * rho * r**2
        if Ixx * w >= angular_momentum and Ixx * w <= angular_momentum*1.1:
            c += 1
            v.append([mass(r, t), r, t])
            V.append(mass(r, t))


for x in v:
    if x[0] == min(V):
        print(c, x)
