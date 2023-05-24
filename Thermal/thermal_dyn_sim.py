import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from DSE.Thermal.env import Env
from DSE.Thermal.materials import Material


class ThermalNode:
    """
    This class initializes a thermal node, taking in a node id and a dictionary of parameters.
    """

    def __init__(self, node_id, name, params: dict, mat: Material):
        self.node_id = node_id
        self.name = name
        self.area = params['area']
        self.radiation_area = params['radiation_area']
        self.contact_area = params['contact_area']
        self.mass = params['mass']
        self.absorptivity = mat.absorptivity
        self.emissivity = mat.emissivity
        self.thermal_conductivity = mat.thermal_conductivity
        self.thermal_capacitance = mat.thermal_capacitance
        self.heat_generated = params['heat_generated']


class ThermalModel:
    """
    This class initializes the thermal model, taking in a list of nodes and a list of connections:
    :param nodes: list of ThermalNode objects
    :param connections: list of tuples of node ids that link two nodes together thermally (for conduction)
    :param env: dictionary of environmental parameters -> may turn into ENV class
    :param t_sim: simulation time in seconds
    :param init_temp: list of initial temperature of each node in Kelvin
    """

    def __init__(self, nodes, connections, env, t_sim, init_temp, with_TCS=True):
        self.nodes: list[ThermalNode] = nodes
        self.connections: list[tuple[int, int]] = connections
        self.env = env
        self.t_sim = t_sim
        self.init_temp = init_temp
        self.solution = None
        self.TCS = with_TCS

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
            return 1414 * self.env.albedo

    def Qir(self, t):
        """
        Infrared flux as function of time
        """
        return 228

    def Qin(self, t, i, node):
        if i==5:    # nadir facing
            Fnad = 1
            q = (Fnad + self.env.albedo) * node.area * self.Qs(t) * node.absorptivity + node.emissivity * node.area * self.Qir(t)
        elif i==4:    # zenith facing
            Fzen = 1
            q = Fzen * node.area * self.Qs(t) * node.absorptivity
            q = (1 + self.env.albedo) * node.area * self.Qs(
                t) * node.absorptivity + node.emissivity * node.area * self.Qir(t)
        elif i==0:    # North facing
            FN = 1
            q = FN * node.area * self.Qs(t) * node.absorptivity
            q = (1 + self.env.albedo) * node.area * self.Qs(
                t) * node.absorptivity + node.emissivity * node.area * self.Qir(t)
        elif i==2:    # South facing
            q = 0
        elif i==1:    # +v facing
            Fpv = 1
            q = Fpv * node.area * self.Qs(t) * node.absorptivity
            q = (1 + self.env.albedo) * node.area * self.Qs(
                t) * node.absorptivity + node.emissivity * node.area * self.Qir(t)
        elif i==3:    # -v facing
            Fmv = 1
            q = Fmv * node.area * self.Qs(t) * node.absorptivity
            q = (1 + self.env.albedo) * node.area * self.Qs(
                t) * node.absorptivity + node.emissivity * node.area * self.Qir(t)
        elif i==7:    # DST
            q = node.area * self.Qs(t) * node.absorptivity
            q = (1 + self.env.albedo) * node.area * self.Qs(t) * node.absorptivity + node.emissivity * node.area * self.Qir(t)
        elif i==8:      # power needed to keep instrument box at 200K
            q = -140 if self.TCS else 0
        elif i==12:      # radiator panel
            q = (1 + self.env.albedo) * node.area * self.Qs(t) * node.absorptivity + node.emissivity * node.area * self.Qir(t)
        elif i==9:     # solar panel
            q = (1 + self.env.albedo) * node.area * self.Qs(t) * node.absorptivity + node.emissivity * node.area * self.Qir(t)
        else:
            q = 0
        return q


    def Tdot(self, t, T):
        """
        Calculates the temperature derivative of each node at a given time, for solve_ivp to solve temperature
        """
        Q = np.empty(len(self.nodes))
        Tdot = np.empty(len(self.nodes))

        for i, node in enumerate(self.nodes):
            cond = sum([node.thermal_conductivity * node.contact_area * (T[connect] - T[i]) for connect in self.connections[i]])

            Q[i] = (self.Qin(t, i, node) + node.heat_generated + cond) - \
                      node.radiation_area * node.emissivity * self.env.sigma_boltzmann * (T[i] ** 4)

            Tdot[i] = Q[i] / (node.thermal_capacitance * node.mass)

        return Tdot

    def solve(self):
        """
        Solves the temperature of each node over the simulation time
        """
        sol = solve_ivp(self.Tdot, (0, self.t_sim), self.init_temp, method='RK45', max_step=10, vectorized=True)
        self.solution = sol

    def plot(self, node_id=None, with_legend=False):
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
        if with_legend:
            plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
            plt.subplots_adjust(right=0.74)
        plt.xlim(7*self.env.t_orbit, 9*self.env.t_orbit)
        plt.savefig('temp_TCS.png')
        plt.show()

    def save_csv(self, filename):
        np.savetxt(filename, np.transpose(self.solution.y), delimiter=',')


def make_sat_1():
    nodes = [ThermalNode(0, 'struc1', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
             ThermalNode(1, 'struc2', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
             ThermalNode(2, 'struc3', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
             ThermalNode(3, 'struc4', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
             ThermalNode(4, 'struc5', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
             ThermalNode(5, 'struc6', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
             ThermalNode(6, 'battery', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 2, 'heat_generated': 5+8, 'radiation_area': 0.1}, Material().battery()),
             ThermalNode(7, 'DST', {'area': 0.3, 'contact_area': 0.1, 'mass': 25, 'heat_generated': 0, 'radiation_area': 0.49}, Material().carbon_fibre()),
             ThermalNode(8, 'inst_box', {'area': 0.022, 'contact_area': 0.02, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.2}, Material().aluminium()),
             ThermalNode(9, 'radiator', {'area': 0.09, 'contact_area': 0, 'mass': 1, 'heat_generated': 0, 'radiation_area': 0.09}, Material().white_coating()),
             ThermalNode(10, 'solar_panel', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 2, 'heat_generated': 0, 'radiation_area': 0.35 * 0.36 * 2}, Material().solar_panel()),
             ThermalNode(11, 'OBC', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 0.3, 'heat_generated': 1+6, 'radiation_area': 0.02}, Material().battery()),
             ThermalNode(12, 'XeTank', {'area': 0.1 * 0.1, 'contact_area': 0.05*0.05, 'mass': 2, 'heat_generated': 0+7, 'radiation_area': 0.13}, Material().aluminium())
    ]

    connections = np.array([[1,3,4,5],
                            [0,2,4,5],
                            [1,3,4,5],
                            [0,2,4,5],
                            [0,1,2,3,10],
                            [0,1,2,3,7,8],
                            [5],
                            [5],
                            [5],
                            [],
                            [4],
                            [5],
                            [3]])

    ENV = Env()
    TM = ThermalModel(nodes, connections, ENV, 10 * ENV.t_orbit, [260] * len(nodes))
    TM.solve()
    TM.plot(with_legend=True)


def make_sat_2():
    nodes = [ThermalNode(0, 'struc1', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(1, 'struc2', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(2, 'struc3', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(3, 'struc4', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(4, 'struc5', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(5, 'structure #6', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(6, 'battery', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 2, 'heat_generated': 5+6, 'radiation_area': 0.05}, Material().battery()),
             ThermalNode(7, 'DST', {'area': 0.3, 'contact_area': 0.1, 'mass': 25, 'heat_generated': 0, 'radiation_area': 0.49}, Material().carbon_fibre()),
             ThermalNode(8, 'inst_box', {'area': 0.022, 'contact_area': 0.02, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.2}, Material().aluminium()),
             ThermalNode(9, 'solar panel', {'area': 3*0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.35 * 0.36 * 6}, Material().solar_panel()),
             ThermalNode(10, 'OBC', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 0.3, 'heat_generated': 1+3, 'radiation_area': 0.02}, Material().PCB()),
             ThermalNode(11, 'Xe Tank', {'area': 0.1 * 0.1, 'contact_area': 0.05*0.05, 'mass': 2, 'heat_generated': 0, 'radiation_area': 0.13}, Material().aluminium()),
             ThermalNode(12, 'radiator', {'area': 0.36*0.22*2, 'contact_area': 0.36*0.22, 'mass': 1.1, 'heat_generated': 0, 'radiation_area': 0.36*0.22*2}, Material().white_coating())]

    connections = np.array([[1,3,4,5,8],
                            [0,2,4,5],
                            [1,3,4,5],
                            [0,2,4,5],
                            [0,1,2,3,10],
                            [0,1,2,3,7],
                            [5],
                            [5],
                            [0],
                            [4],
                            [5],
                            [3],
                            []])

    ENV = Env()
    TM = ThermalModel(nodes, connections, ENV, 10 * ENV.t_orbit, [273] * len(nodes))
    TM.solve()
    TM.plot([5,6,7,8,9,10,11,12], with_legend=True)


def make_nice_plot():
    ENV = Env()

    nodes_none = [ThermalNode(0, 'struc1', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(1, 'struc2', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(2, 'struc3', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(3, 'struc4', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(4, 'struc5', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(5, 'structure', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(6, 'battery', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 2, 'heat_generated': 5, 'radiation_area': 0.05}, Material().battery()),
             ThermalNode(7, 'DST', {'area': 0.3, 'contact_area': 0.1, 'mass': 25, 'heat_generated': 0, 'radiation_area': 0.49}, Material().carbon_fibre()),
             ThermalNode(8, 'inst_box', {'area': 0.022, 'contact_area': 0.02, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.2}, Material().aluminium()),
             ThermalNode(9, 'solar panel', {'area': 3*0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.35 * 0.36 * 6}, Material().solar_panel()),
             ThermalNode(10, 'OBC', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 0.3, 'heat_generated': 1, 'radiation_area': 0.02}, Material().PCB()),
             ThermalNode(11, 'Xe Tank', {'area': 0.1 * 0.1, 'contact_area': 0.05*0.05, 'mass': 2, 'heat_generated': 0, 'radiation_area': 0.13}, Material().aluminium())]
    connections_none = np.array([[1,3,4,5,8],
                            [0,2,4,5],
                            [1,3,4,5],
                            [0,2,4,5],
                            [0,1,2,3,10],
                            [0,1,2,3,7],
                            [5],
                            [5],
                            [0],
                            [4],
                            [5],
                            [3]])
    TM1 = ThermalModel(nodes_none, connections_none, ENV, 10 * ENV.t_orbit, [273] * len(nodes_none), with_TCS=False)
    TM1.solve()

    nodes_TCS = [ThermalNode(0, 'struc1', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(1, 'struc2', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(2, 'struc3', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(3, 'struc4', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(4, 'struc5', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(5, 'structure', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(6, 'battery', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 2, 'heat_generated': 5+6, 'radiation_area': 0.05}, Material().battery()),
             ThermalNode(7, 'DST', {'area': 0.3, 'contact_area': 0.1, 'mass': 25, 'heat_generated': 0, 'radiation_area': 0.49}, Material().carbon_fibre()),
             ThermalNode(8, 'inst_box', {'area': 0.022, 'contact_area': 0.02, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.2}, Material().aluminium()),
             ThermalNode(9, 'solar panel', {'area': 3*0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.35 * 0.36 * 6}, Material().solar_panel()),
             ThermalNode(10, 'OBC', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 0.3, 'heat_generated': 1+3, 'radiation_area': 0.02}, Material().PCB()),
             ThermalNode(11, 'Xe Tank', {'area': 0.1 * 0.1, 'contact_area': 0.05*0.05, 'mass': 2, 'heat_generated': 0, 'radiation_area': 0.13}, Material().aluminium()),
             ThermalNode(12, 'radiator', {'area': 0.36*0.22*2, 'contact_area': 0.36*0.22, 'mass': 1.1, 'heat_generated': 0, 'radiation_area': 0.36*0.22*2}, Material().white_coating())]
    connections_TCS = np.array([[1,3,4,5,8],
                            [0,2,4,5],
                            [1,3,4,5],
                            [0,2,4,5],
                            [0,1,2,3,10],
                            [0,1,2,3,7],
                            [5],
                            [5],
                            [0],
                            [4],
                            [5],
                            [3],
                            []])
    TM2 = ThermalModel(nodes_TCS, connections_TCS, ENV, 10 * ENV.t_orbit, [273] * len(nodes_TCS))
    TM2.solve()

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))

    for i, node in enumerate(TM1.nodes[5:], start=5):
        a1.plot(TM1.solution.t, TM1.solution.y[i], label=f'{node.name}')
        a1.set_xlim(7 * TM1.env.t_orbit, 9 * TM1.env.t_orbit)
        a1.set_xlabel('Time [s]')
        a1.set_ylabel('Temperature [K]')
        # a1.legend(bbox_to_anchor=(0.97, 0.5), loc="center left")

    for i, node in enumerate(TM2.nodes[5:], start=5):
        a2.plot(TM2.solution.t, TM2.solution.y[i], label=f'{node.name}')
        a2.set_xlim(7 * TM2.env.t_orbit, 9 * TM2.env.t_orbit)
        a2.set_xlabel('Time [s]')
        a2.set_ylabel('Temperature [K]')
        a2.legend(bbox_to_anchor=(0.97, 0.5), loc="center left")

    plt.savefig('TCS_comparison.png')




if __name__ == '__main__':
    make_nice_plot()
