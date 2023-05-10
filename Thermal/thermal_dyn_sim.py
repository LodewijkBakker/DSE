import numpy as np
import numba as nb
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from env import Env
from materials import Material


class ThermalNode:
    """
    This class initializes a thermal node, taking in a node id and a dictionary of parameters.
    """

    def __init__(self, node_id, name, params: dict, mat: Material):
        self.node_id = node_id
        self.name = name
        self.area = params['area']
        self.radiation_area = params['radiation_area']
        self.contact_area = 0
        self.mass = params['mass']
        self.absorptivity = mat.absorptivity
        self.emissivity = mat.emissivity
        self.thermal_conductivity = mat.thermal_conductivity
        self.thermal_capacitance = mat.thermal_capacitance
        self.heat_generated = params['heat_generated']
        self.interior = params['interior']
        self.solution = None


class ThermalModel:
    """
    This class initializes the thermal model, taking in a list of nodes and a list of connections:
    :param nodes: list of ThermalNode objects
    :param connections: list of tuples of node ids that link two nodes together thermally (for conduction)
    :param env: dictionary of environmental parameters -> may turn into ENV class
    :param t_sim: simulation time in seconds
    :param init_temp: list of initial temperature of each node in Kelvin
    """

    def __init__(self, nodes, connections, env, t_sim, init_temp):
        self.nodes = nodes
        self.connections: list[tuple[int, int]] = connections
        self.env = env
        self.t_sim = t_sim
        self.init_temp = init_temp

    def Qs(self, t):
        """
        Solar flux as function of time
        """
        if (t % self.env.t_orbit) / self.env.t_orbit > (1 - self.env.f_orbit):
            return 0
        else:
            return 1414

    def Qa(self, t):
        """
        Albedo flux as function of time
        """
        if (t % self.env.t_orbit) / self.env.t_orbit > (1 - self.env.f_orbit):
            return 0
        else:
            return 1414 * 0.14

    def Qir(self, t):
        """
        Infrared flux as function of time
        """
        return 228

    def Tdot(self, t, T):
        """
        Calculates the temperature derivative of each node at a given time, for solve_ivp to solve temperature
        """
        Q = np.empty((len(self.nodes), 1))
        Tdot = np.empty((len(self.nodes), 1))

        for i, node in enumerate(self.nodes):

            cond = sum([node.thermal_conductivity * node.area * (T[connect, 0] - T[i, 0]) for connect in self.connections[i]])

            if not node.interior:
                Q[i, 0] = (node.absorptivity * node.area * (self.Qs(t) + self.Qa(t)) +
                           node.emissivity * node.area * self.Qir(t) +
                           node.heat_generated + cond) - \
                          node.radiation_area * node.emissivity * self.env.sigma_boltzmann * (T[i, 0] ** 4)
            else:
                Q[ i, 0] = node.heat_generated + cond - node.radiation_area * node.emissivity * self.env.sigma_boltzmann * (
                            T[i, 0] ** 4)

            Tdot[i, 0] = Q[i, 0] / (node.thermal_capacitance * node.mass)

        return Tdot

    def solve(self):
        """
        Solves the temperature of each node over the simulation time
        """
        sol = solve_ivp(self.Tdot, (0, self.t_sim), self.init_temp, method='RK45', max_step=10, vectorized=True)
        self.solution = sol

    def plot(self, node_id=None):
        """
        Plots the temperature of each node over the simulation time,
        can also specify which nodes to plot
        """
        if not node_id:
            for i, node in enumerate(self.nodes):
                plt.plot(self.solution.t, self.solution.y[i], label=f'{node.name}')
        else:
            for i in node_id:
                plt.plot(self.solution.t, self.solution.y[i], label=f'{self.nodes[i].name}')
        plt.legend()
        plt.show()

def make_sat_1():
    nodes = [ThermalNode(0, 'struc1', {'area': 0.35 * 0.34, 'mass': 0.83, 'heat_generated': 2, 'radiation_area': 0.72, 'interior': False}, Material().anodized_aluminium()),
             ThermalNode(1, 'struc2', {'area': 0.35 * 0.34, 'mass': 0.83, 'heat_generated': 2, 'radiation_area': 0.72, 'interior': False}, Material().anodized_aluminium()),
             ThermalNode(2, 'struc3', {'area': 0.35 * 0.34, 'mass': 0.83, 'heat_generated': 2, 'radiation_area': 0.72, 'interior': False}, Material().anodized_aluminium()),
             ThermalNode(3, 'struc4', {'area': 0.35 * 0.34, 'mass': 0.83, 'heat_generated': 2, 'radiation_area': 0.72, 'interior': False}, Material().anodized_aluminium()),
             ThermalNode(4, 'struc5', {'area': 0.35 * 0.34, 'mass': 0.83, 'heat_generated': 2, 'radiation_area': 0.72, 'interior': False}, Material().anodized_aluminium()),
             ThermalNode(5, 'struc6', {'area': 0.35 * 0.36, 'mass': 0.83, 'heat_generated': 2, 'radiation_area': 0.72, 'interior': False}, Material().anodized_aluminium()),
             ThermalNode(6, 'battery', {'area': 0.1 * 0.1, 'mass': 2, 'heat_generated': 5, 'radiation_area': 0.1, 'interior': True}, Material().battery()),
             ThermalNode(7, 'DST', {'area': 0.3, 'mass': 25, 'heat_generated': 0, 'radiation_area': 0.49, 'interior': False}, Material().carbon_fibre()),
             ThermalNode(8, 'inst_box', {'area': 0.022, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.08, 'interior': True}, Material().aluminium())]

    connections = np.array([[1,3,4,5],
                            [0,2,4,5],
                            [1,3,4,5],
                            [0,2,4,5],
                            [0,1,2,3],
                            [0,1,2,3,7,8],
                            [5],
                            [5,8],
                            [5,7]])

    ENV = Env()
    TM = ThermalModel(nodes, connections, ENV, 30 * ENV.t_orbit, [280] * 9)
    TM.solve()
    TM.plot()



if __name__ == '__main__':
    pass
    make_sat_1()
