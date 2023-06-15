"""
------------------------------------------------------------------------------------------------------------------------
                                        SIMULATION OF THERMAL DYNAMICS
------------------------------------------------------------------------------------------------------------------------
This file contains the code that simulates the thermal dynamics of the satellite.
To create a satellite configuration, refer to the configurations.py file in the same folder.

The model is defined by a list of nodes (ThermalNode) with their respective parameters, and a list of connections.
For model solves the following equation:
    m * C * dT/dt = (Q_generated + Q_cond + Q_solar + Q_albedo + Q_EIR) - A * e * sigma * T^4

    Q_generated - heat generated by the components
    Q_cond - sum of heat conducted from adjacent nodes
    Q_solar - heat from solar radiation => Q_solar = A * alpha * q_solar
    Q_albedo - heat from albedo radiation => Q_albedo = A * alpha * q_solar * albedo
    Q_EIR - heat from Earth IR => Q_EIR = A * epsilon * q_EIR

The function Q_in(t, i, node) define the environmental fluxes (solar, albedo, Earth IR) for each node, which are
different for the exterior nodes depending on orientation to Sun, and Q_in = 0 for interior nodes.

The function Tdot(t, T) defines the differential equation to be solved by solve().
"""


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from DSE.Thermal.materials import Material, Coating
from DSE.Thermal.esatan_reader import heat_flows, interpolate_points, prepare_heat_flows


class ThermalNode:
    """
    This class initializes a thermal node, taking in a node id, a dictionary of geometric parameters, a material and
    a coating.
    params -> dictionary e.g. {'area': 0.5, 'radiation_area': 0.8, 'contact_area': [0.2, 0.0, 0.3], 'mass': 3,
    'heat_generated': 3}
    material -> Material() object [see materials.py file]
    coating -> Coating() object [see materials.py file]
    """

    def __init__(self, node_id: int, name: str, params: dict, mat: Material, coat: Coating):
        self.node_id = node_id
        self.name = name
        self.area = params['area']      # incident radiation area [m^2]
        self.radiation_area = params['radiation_area']      # radiation emission area [m^2]
        self.contact_area = params['contact_area']      # list of contact areas to connecting nodes [m^2]
        self.mass = params['mass']      # mass of node [kg]
        self.absorptivity = coat.absorptivity
        self.emissivity = coat.emissivity
        self.thermal_conductivity = mat.thermal_conductivity
        self.thermal_capacitance = mat.thermal_capacitance
        self.heat_generated = params['heat_generated']  # heat generated by node [W]


class ThermalModel:
    """
    This class initializes the thermal model, taking in a list of nodes and a list of connections:
    @param nodes: list of ThermalNode objects
    @param connections: list of tuples of node ids that link two nodes together thermally (for conduction)
    @param env: dictionary of environmental parameters -> may turn into ENV class
    @param t_sim: simulation time in seconds
    @param init_temp: list of initial temperature of each node in Kelvin
    @param n_orbits: number of orbits to be generated (defuault value 1)
    @param: Q_ESATAN: list of heat values for each node, taken from ESATAN Radiative simulation, instead of the heat
    values assumed in Qin()
    """

    def __init__(self, nodes, connections, env, init_temp, n_orbits=1, sun_shield=True, unit='K', ESATAN=True, Q_ESATAN=None):
        self.nodes: list[ThermalNode] = nodes
        self.connections: list[tuple[int, int]] = connections
        self.env = env
        self.n_orbits = n_orbits
        self.t_sim = self.n_orbits * self.env.t_orbit - 1
        self.init_temp = init_temp
        self.solution = None
        self.sun_shield = sun_shield
        self.unit = unit
        self.ESATAN = ESATAN
        self.Q_ESATAN = Q_ESATAN
        self.Tdot_arr = [[] for i in range(len(self.nodes))]

    def Qs(self, t, i):
        """
        Solar flux as function of time (takes into consideration eclipse fraction)
        """
        if (t % self.env.t_orbit) / self.env.t_orbit > (1 - self.env.f_orbit):
            return 0
        else:
            return 1414

    def Qa(self, t):
        """
        Albedo flux as function of time (takes into consideration eclipse fraction)
        """
        if (t % self.env.t_orbit) / self.env.t_orbit > (1 - self.env.f_orbit):
            return 0
        else:
            return 1414 * self.env.albedo

    def Qir(self, t, i):
        """
        Infrared flux as function of time
        """
        if i==4:
            return 0
        else:
            return 228

    def Qin(self, t, i, node):
        """
        Absorbed environmental heat values for each node
        t - time (used in solve_ivp)
        i - node id
        node - node object
        """
        if i==5:    # nadir facing
            q = self.env.albedo * node.area * self.Qs(t, i) * node.absorptivity + node.emissivity * node.area * self.Qir(t,i)

        elif i==4:    # zenith facing
            q = node.area * self.Qs(t, i) * node.absorptivity

        elif i==0:    # North facing
            q = node.area * self.env.albedo * self.Qs(t, i) * node.absorptivity + node.emissivity * node.area * self.Qir(t,i)

        elif i==2:    # South facing
            q = node.area * self.env.albedo * self.Qs(t, i) * node.absorptivity + node.emissivity * node.area * self.Qir(t,i)

        elif i==1:    # +v facing
            q = (1 + self.env.albedo) * node.area * self.Qs(t, i) * node.absorptivity + node.emissivity * node.area * self.Qir(t,i)

        elif i==3:    # -v facing
            q = (1 + self.env.albedo) * node.area * self.Qs(t, i) * node.absorptivity + node.emissivity * node.area * self.Qir(t,i)

        elif i in [7, 13, 14]:   # solar panels
            q = (1 + self.env.albedo) * node.area * self.Qs(t, i) * node.absorptivity + node.emissivity * node.area * self.Qir(t,i)

        elif i ==10:   # radiator
            if self.sun_shield:
                q = node.area * self.Qs(t, i) * node.absorptivity
            else:
                q = (1 + self.env.albedo) * node.area * self.Qs(t, i) * node.absorptivity + node.emissivity * node.area * self.Qir(t,i)

        elif i ==12:   # DST baffle
            q = (1 + self.env.albedo) * node.area * self.Qs(t, i) * node.absorptivity + node.emissivity * node.area * self.Qir(t,i)

        else:      # battery, OBC, DST, propellant
            q = 0

        return q


    def plotting_Q(self):
        Q = np.empty((len(self.nodes), self.env.t_orbit))
        for t in range(0, self.env.t_orbit):
            for i in range(len(self.nodes)):
                Q[i, t] = self.Qin(t, i, self.nodes[i])
        for i in range(len(self.nodes)):
            plt.plot(Q[i, :], label=self.nodes[i].name)
        plt.legend()
        plt.show()
        plt.close()


    def Tdot(self, t, T):
        """
        Calculates the temperature derivative of each node at a given time, for solve_ivp to solve temperature
        """
        Q = np.empty(len(self.nodes))
        Tdot = np.empty(len(self.nodes))

        for i, node in enumerate(self.nodes):
            if i == 11:     # this is used to account for the high thermal conductivity required between the instrument box and radiator
                cond = sum([self.nodes[connect].thermal_conductivity * node.contact_area[connect] * (T[connect] - T[i]) for connect in
                            self.connections[i]])
            else:        # normal conductive heat transfer between nodes
                cond = sum([node.thermal_conductivity * node.contact_area[connect] * (T[connect] - T[i]) for connect in self.connections[i]])

            if not self.ESATAN:
                Q[i] = (self.Qin(t, i, node) + node.heat_generated + cond) - \
                          node.radiation_area * node.emissivity * self.env.sigma_boltzmann * (T[i] ** 4)
            else:
                Q[i] = (self.Q_ESATAN[i][int(t)] + node.heat_generated + cond) - \
                          node.radiation_area * node.emissivity * self.env.sigma_boltzmann * (T[i] ** 4)

            Tdot[i] = Q[i] / (node.thermal_capacitance * node.mass)
            self.Tdot_arr[i].append(Tdot[i])

        return Tdot

    def solve(self):
        """
        Solves the temperature of each node over the simulation time
        """
        # self.plotting_Q()
        sol = solve_ivp(self.Tdot, (0, self.t_sim), self.init_temp, method='RK45', t_eval=np.arange(0, self.t_sim+1, 1),
                        vectorized=True)
        self.solution = sol

    def t_solver(self, Qin, Qgen, A_rad, emissivity):
        """
        Calculates the steady state temperature, given the incident and generated heat for a node
        """
        return ((Qin + Qgen) / (A_rad * emissivity * 5.67 * 10 ** -8)) ** 0.25

    def plot(self, node_id=None, with_legend=False, save=None):
        """
        Plots the temperature of each node over the simulation time,
        can also specify which nodes to plot
        """
        if self.unit != 'K':
            # plt.plot(self.solution.t, np.mean(self.solution.y[:6]-273, axis=0), label='structure')
            if not node_id:
                for i, node in enumerate(self.nodes):
                    plt.plot(self.solution.t, self.solution.y[i]-273, label=f'{node.name}')
            else:
                for i in node_id:
                    plt.plot(self.solution.t, self.solution.y[i]-273, label=f'{self.nodes[i].name}')
        else:
            # plt.plot(self.solution.t, np.mean(self.solution.y[:6], axis=0), label='structure')
            if not node_id:
                for i, node in enumerate(self.nodes):
                    plt.plot(self.solution.t, self.solution.y[i], label=f'{node.name}')
            else:
                for i in node_id:
                    plt.plot(self.solution.t, self.solution.y[i], label=f'{self.nodes[i].name}')
            plt.plot(self.solution.t, np.ones(len(self.solution.t)) * self.t_solver(0, 0.8, 0.016, 0.035),
                     label='pipes')
            plt.plot(self.solution.t, np.ones(len(self.solution.t)) * self.t_solver(0, 1.35, 0.035, 0.035),
                     label='propellant')

        if with_legend:
            plt.legend(bbox_to_anchor=(1, 0.5), loc="center left")
            plt.subplots_adjust(right=0.74)
        plt.xlabel('Time [s]')
        plt.ylabel(f'Temperature [{self.unit}]')
        plt.xlim(9*self.env.t_orbit, 10*self.env.t_orbit)
        # plt.xticks(np.arange(8*self.env.t_orbit, 9*self.env.t_orbit, 1000), np.arange(0, int(self.env.t_orbit), 1000))
        plt.title('Nodal Temperatues')
        if save:
            plt.savefig(save)
        # plt.show()
        plt.close()
        print('CMG', self.t_solver(0, 0.1, 0.06, 0.035))

    def save_csv(self, filename):
        if self.solution is not None:
            np.savetxt(filename, np.transpose(self.solution.y), delimiter=',')
