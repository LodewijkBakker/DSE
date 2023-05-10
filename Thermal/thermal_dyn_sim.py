import numpy as np
import numba as nb
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class ThermalNode:
    def __init__(self, node_id, params: dict):
        self.node_id = node_id
        self.area = params['area']
        self.mass = params['mass']
        self.radiation_area = params['radiation_area']
        self.absorptivity = params['absorptivity']
        self.emissivity = params['emissivity']
        self.heat_generated = params['heat_generated']
        self.thermal_conductivity = params['thermal_conductivity']
        self.thermal_capacitance = params['thermal_capacitance']
        self.interior = params['interior']


class ThermalModel:
    def __init__(self, nodes, connections, env, t_sim):
        self.nodes = nodes
        self.connections: list[tuple[int, int]] = connections
        self.env = env
        self.t_sim = t_sim

    def Qs(self, t):
        if (t % self.env['t_orbit']) / self.env['t_orbit'] > (1 - self.env['f_orbit']):
            return 0
        else:
            return 1414

    def Qa(self, t):
        if (t % self.env['t_orbit']) / self.env['t_orbit'] > (1 - self.env['f_orbit']):
            return 0
        else:
            return 1414 * 0.14

    def Qir(self, t):
        return 228


    def Tdot(self, t, T):
        Q = np.empty((len(self.nodes), 1))
        Tdot = np.empty((len(self.nodes), 1))
        for i, node in enumerate(self.nodes):

            cond = sum(sum([node.thermal_conductivity * node.area * (T[connect, 0] - T[i, 0]) for connect in self.connections if connect[0] == i]))

            Q[i,0] = (node.absorptivity * node.area * (self.Qs(t) + self.Qa(t)) +
                      node.emissivity * node.area * self.Qir(t) +
                      node.heat_generated + cond) -\
                      node.radiation_area * node.emissivity * self.env['sigma_boltzmann'] * (T[i, 0] ** 4)

            Tdot[i,0] = Q[i,0] / (node.thermal_capacitance * node.mass)

        return Tdot


sigma_boltzmann = 5.670374419*10**(-8)
m_sat = 54
alt_orbit = (6371+300)*1e3
mu_earth = 3.986004418e14
v_orbit = ((mu_earth/alt_orbit)**0.5)
t_orbit = 2*np.pi*alt_orbit/v_orbit
f_orbit = 0.404

nodes = [ThermalNode(0, {'area': 0.35 * 0.34, 'mass': 4, 'absorptivity': 0.25, 'emissivity': 0.88, 'heat_generated': 200, 'thermal_conductivity': 167, 'thermal_capacitance': 896, 'radiation_area': 0.6, 'interior': False}),
         ThermalNode(1, {'area': 0.35 * 0.2, 'mass': 4, 'absorptivity': 0.25, 'emissivity': 0.88, 'heat_generated': 180, 'thermal_conductivity': 167, 'thermal_capacitance': 896, 'radiation_area': 0.6, 'interior': False}),
         ThermalNode(2, {'area': 0.35 * 0.2, 'mass': 2, 'absorptivity': 0.25, 'emissivity': 0.88, 'heat_generated': 300, 'thermal_conductivity': 50, 'thermal_capacitance': 896, 'radiation_area': 0.6, 'interior': False})]

TM = ThermalModel(nodes, [(0,1), (1,0), (1, 2), (2,1)], {'Qs': 1414, 'Qa': 1414*0.14, 'Qir': 228, 't_orbit': t_orbit, 'f_orbit': f_orbit, 'sigma_boltzmann': sigma_boltzmann}, 50*3600)
T0 = [288.15, 288.15, 288.15]

sol = solve_ivp(fun=TM.Tdot, t_span=(0, 10*t_orbit), y0=T0, method='RK45', max_step=10, min_step=3, vectorized=True)
plt.plot(sol.t, sol.y[0], sol.t, sol.y[1], sol.t, sol.y[2])
plt.show()


