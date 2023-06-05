"""
------------------------------------------------------------------------------------------------------------------------
                                CONFIGURATION FILE TO SET UP THE THERMAL MODEL
------------------------------------------------------------------------------------------------------------------------
To run the thermal model, you need to first define some parameter:
    - Nodes: a list of ThermalNode objects
    - Connections: a array of node ids that link two nodes together thermally (for conduction)
    - Environment: a dictionary of environmental parameters
    - Simulation time: the time in seconds that the simulation will run for
    - Initial temperature: a list of initial temperature of each node in Kelvin

Then, you can instantiate a ThermalModel object with the above parameters and solve the simulation (and plot it), e.g.
    model = ThermalModel(nodes, connections, env, t_sim, init_temp)
    model.solve()
    model.plot(with_legend=True)
"""

import numpy as np
import toml
from DSE.Thermal.env import Env
from DSE.Thermal.materials import Material
from DSE.Thermal.thermal_dyn_sim import ThermalNode, ThermalModel




# def make_sat_1():
#     nodes = [ThermalNode(0, 'struc1', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
#              ThermalNode(1, 'struc2', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
#              ThermalNode(2, 'struc3', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
#              ThermalNode(3, 'struc4', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
#              ThermalNode(4, 'struc5', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
#              ThermalNode(5, 'struc6', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.56}, Material().anodized_aluminium()),
#              ThermalNode(6, 'battery', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 2, 'heat_generated': 5+8, 'radiation_area': 0.1}, Material().battery()),
#              ThermalNode(7, 'DST', {'area': 0.3, 'contact_area': 0.1, 'mass': 25, 'heat_generated': 0, 'radiation_area': 0.49}, Material().carbon_fibre()),
#              ThermalNode(8, 'inst_box', {'area': 0.022, 'contact_area': 0.02, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.2}, Material().aluminium()),
#              ThermalNode(9, 'radiator', {'area': 0.09, 'contact_area': 0, 'mass': 1, 'heat_generated': 0, 'radiation_area': 0.09}, Material().white_coating()),
#              ThermalNode(10, 'solar_panel', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 2, 'heat_generated': 0, 'radiation_area': 0.35 * 0.36 * 2}, Material().solar_panel()),
#              ThermalNode(11, 'OBC', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 0.3, 'heat_generated': 1+6, 'radiation_area': 0.02}, Material().battery()),
#              ThermalNode(12, 'XeTank', {'area': 0.1 * 0.1, 'contact_area': 0.05*0.05, 'mass': 2, 'heat_generated': 0+7, 'radiation_area': 0.13}, Material().aluminium())
#     ]
#
#     connections = np.array([[1,3,4,5],
#                             [0,2,4,5],
#                             [1,3,4,5],
#                             [0,2,4,5],
#                             [0,1,2,3,10],
#                             [0,1,2,3,7,8],
#                             [5],
#                             [5],
#                             [5],
#                             [],
#                             [4],
#                             [5],
#                             [3]])
#
#     ENV = Env()
#     TM = ThermalModel(nodes, connections, ENV, 10 * ENV.t_orbit, [260] * len(nodes))
#     TM.solve()
#     TM.plot(with_legend=True)
# def make_sat_2():
#     nodes = [ThermalNode(0, 'struc1', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(1, 'struc2', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(2, 'struc3', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(3, 'struc4', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(4, 'struc5', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(5, 'structure #6', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(6, 'battery', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 2, 'heat_generated': 5+6, 'radiation_area': 0.05}, Material().battery()),
#              ThermalNode(7, 'DST', {'area': 0.3, 'contact_area': 0.1, 'mass': 25, 'heat_generated': 0, 'radiation_area': 0.49}, Material().carbon_fibre()),
#              ThermalNode(8, 'inst_box', {'area': 0.022, 'contact_area': 0.02, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.2}, Material().aluminium()),
#              ThermalNode(9, 'solar panel', {'area': 3*0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.35 * 0.36 * 6}, Material().solar_panel()),
#              ThermalNode(10, 'OBC', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 0.3, 'heat_generated': 1+3, 'radiation_area': 0.02}, Material().PCB()),
#              ThermalNode(11, 'Xe Tank', {'area': 0.1 * 0.1, 'contact_area': 0.05*0.05, 'mass': 2, 'heat_generated': 0, 'radiation_area': 0.13}, Material().aluminium()),
#              ThermalNode(12, 'radiator', {'area': 0.36*0.22*2, 'contact_area': 0.36*0.22, 'mass': 1.1, 'heat_generated': 0, 'radiation_area': 0.36*0.22*2}, Material().white_coating())]
#
#     connections = np.array([[1,3,4,5,8],
#                             [0,2,4,5],
#                             [1,3,4,5],
#                             [0,2,4,5],
#                             [0,1,2,3,10],
#                             [0,1,2,3,7],
#                             [5],
#                             [5],
#                             [0],
#                             [4],
#                             [5],
#                             [3],
#                             []])
#
#     ENV = Env()
#     TM = ThermalModel(nodes, connections, ENV, 10 * ENV.t_orbit, [273] * len(nodes))
#     TM.solve()
#     TM.plot([5,6,7,8,9,10,11,12], with_legend=True)
# def make_nice_plot():
#     ENV = Env()
#
#     nodes_none = [ThermalNode(0, 'struc1', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(1, 'struc2', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(2, 'struc3', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(3, 'struc4', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(4, 'struc5', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(5, 'structure', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(6, 'battery', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 2, 'heat_generated': 5, 'radiation_area': 0.05}, Material().battery()),
#              ThermalNode(7, 'DST', {'area': 0.3, 'contact_area': 0.1, 'mass': 25, 'heat_generated': 0, 'radiation_area': 0.49}, Material().carbon_fibre()),
#              ThermalNode(8, 'inst_box', {'area': 0.022, 'contact_area': 0.02, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.2}, Material().aluminium()),
#              ThermalNode(9, 'solar panel', {'area': 3*0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.35 * 0.36 * 6}, Material().solar_panel()),
#              ThermalNode(10, 'OBC', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 0.3, 'heat_generated': 1, 'radiation_area': 0.02}, Material().PCB()),
#              ThermalNode(11, 'Xe Tank', {'area': 0.1 * 0.1, 'contact_area': 0.05*0.05, 'mass': 2, 'heat_generated': 0, 'radiation_area': 0.13}, Material().aluminium())]
#     connections_none = np.array([[1,3,4,5,8],
#                             [0,2,4,5],
#                             [1,3,4,5],
#                             [0,2,4,5],
#                             [0,1,2,3,10],
#                             [0,1,2,3,7],
#                             [5],
#                             [5],
#                             [0],
#                             [4],
#                             [5],
#                             [3]])
#     TM1 = ThermalModel(nodes_none, connections_none, ENV, 10 * ENV.t_orbit, [273] * len(nodes_none), with_TCS=False)
#     TM1.solve()
#
#     nodes_TCS = [ThermalNode(0, 'struc1', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(1, 'struc2', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(2, 'struc3', {'area': 0.36 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(3, 'struc4', {'area': 0.35 * 0.22, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(4, 'struc5', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(5, 'structure', {'area': 0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 0.83, 'heat_generated': 8, 'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(6, 'battery', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 2, 'heat_generated': 5+6, 'radiation_area': 0.05}, Material().battery()),
#              ThermalNode(7, 'DST', {'area': 0.3, 'contact_area': 0.1, 'mass': 25, 'heat_generated': 0, 'radiation_area': 0.49}, Material().carbon_fibre()),
#              ThermalNode(8, 'inst_box', {'area': 0.022, 'contact_area': 0.02, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.2}, Material().aluminium()),
#              ThermalNode(9, 'solar panel', {'area': 3*0.35 * 0.36, 'contact_area': 0.35*0.05, 'mass': 5, 'heat_generated': 5, 'radiation_area': 0.35 * 0.36 * 6}, Material().solar_panel()),
#              ThermalNode(10, 'OBC', {'area': 0.1 * 0.1, 'contact_area': 0.1*0.1, 'mass': 0.3, 'heat_generated': 1+3, 'radiation_area': 0.02}, Material().PCB()),
#              ThermalNode(11, 'Xe Tank', {'area': 0.1 * 0.1, 'contact_area': 0.05*0.05, 'mass': 2, 'heat_generated': 0, 'radiation_area': 0.13}, Material().aluminium()),
#              ThermalNode(12, 'radiator', {'area': 0.36*0.22*2, 'contact_area': 0.36*0.22, 'mass': 1.1, 'heat_generated': 0, 'radiation_area': 0.36*0.22*2}, Material().white_coating())]
#     connections_TCS = np.array([[1,3,4,5,8],
#                             [0,2,4,5],
#                             [1,3,4,5],
#                             [0,2,4,5],
#                             [0,1,2,3,10],
#                             [0,1,2,3,7],
#                             [5],
#                             [5],
#                             [0],
#                             [4],
#                             [5],
#                             [3],
#                             []])
#     TM2 = ThermalModel(nodes_TCS, connections_TCS, ENV, 10 * ENV.t_orbit, [273] * len(nodes_TCS))
#     TM2.solve()
#
#     # fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
#     #
#     # for i, node in enumerate(TM1.nodes[5:], start=5):
#     #     a1.plot(TM1.solution.t, TM1.solution.y[i], label=f'{node.name}')
#     #     a1.set_xlim(7 * TM1.env.t_orbit, 9 * TM1.env.t_orbit)
#     #     a1.set_xlabel('Time [s]')
#     #     a1.set_ylabel('Temperature [K]')
#     #     a1.legend(bbox_to_anchor=(0.97, 0.5), loc="center left")
#     #
#     # for i, node in enumerate(TM2.nodes[5:], start=5):
#     #     a2.plot(TM2.solution.t, TM2.solution.y[i], label=f'{node.name}')
#     #     a2.set_xlim(7 * TM2.env.t_orbit, 9 * TM2.env.t_orbit)
#     #     a2.set_xlabel('Time [s]')
#     #     a2.set_ylabel('Temperature [K]')
#     #     a2.legend(bbox_to_anchor=(0.97, 0.5), loc="center left")
#
#     for i, node in enumerate(TM1.nodes[5:], start=5):
#         plt.plot(TM1.solution.t, TM1.solution.y[i], label=f'{node.name}')
#         plt.xlim(7 * TM1.env.t_orbit, 9 * TM1.env.t_orbit)
#         plt.xlabel('Time [s]')
#         plt.ylabel('Temperature [K]')
#         plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
#         plt.subplots_adjust(right=0.74)
#
#     plt.savefig('TCS_none.png')

def make_TCS(output=None):
    ENV = Env()

    nodes = [ThermalNode(0, 'North',
                              {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.22}, Material().white_coating()),
             ThermalNode(1, '+V',
                              {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(2, 'South',
                              {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(3, '-V',
                              {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003, 0], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(4, 'Zenith',
                              {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0, 0], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(5, 'Nadir',
                              {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0, 0, 0,0,0,0,0,0, 0.003], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, Material().white_coating()),
             ThermalNode(6, 'battery',
                              {'area': 0.35 * 0.36, 'contact_area': [0,0,0,0.2*0.2*0.5,0,0,0], 'mass': 2.5, 'heat_generated': 5,
                               'radiation_area': 0.2*0.2*6}, Material().battery()),
             ThermalNode(7, 'solar',
                              {'area': 0.63, 'contact_area': [], 'mass': 0.5, 'heat_generated': 0,
                               'radiation_area': 0.63*2}, Material().coated_solar_panel()),
             ThermalNode(8, 'OBC',
                              {'area': 0.1, 'contact_area': [0,0,0,0.1*0.1*0.5,0,0,0], 'mass': 0.25, 'heat_generated': 5,
                               'radiation_area': 0.1*0.1*3}, Material().PCB()),
             ThermalNode(9, 'propellant',
                              {'area': 0.1, 'contact_area': [0,0,0,0,0,0.3*0.1,0,0,0], 'mass': 0.25, 'heat_generated': 0,
                               'radiation_area': 0.1}, Material().MLI()),
             ThermalNode(10, 'radiator',
                              {'area': 0.35*0.22, 'contact_area': [0,0,0,0,0,0,0,0,0,0,0,0.01], 'mass': 0.6, 'heat_generated': 0,
                               'radiation_area': 0.35*0.22}, Material().OSR()),
             ThermalNode(11, 'DST box',
                              {'area': 0.11, 'contact_area': [0,0,0,0,0,0.19*0.12,0,0,0,0,0.01,0], 'mass': 5, 'heat_generated': 2,
                               'radiation_area': 0.11}, Material().MLI()),
             ThermalNode(12, 'DST baffle',
                             {'area': 0.36/2, 'contact_area': [0,0,0,0,0,0.55], 'mass': 25, 'heat_generated': 0,
                              'radiation_area': 0.36*2}, Material().white_coating())]

    connections = np.array([[1,3,4,5],
                            [0,2,4,5],
                            [1,3,4,5],
                            [0,2,4,5],
                            [0,1,2,3],
                            [0,1,2,3],
                            [3],
                            [],
                            [3],
                            [5],
                            [],
                            [5,10],
                            [5]])

    TM = ThermalModel(nodes, connections, ENV, 10 * ENV.t_orbit, [270] * len(nodes))
    TM.solve()
    TM.plot([6,7,8,9,10,11,12], with_legend=True)

    if output:
        res = {node.name: {'max': np.max(TM.solution.y[i][1000:]), 'min': np.min(TM.solution.y[i][1000:])}
               for i, node in enumerate(TM.nodes)}
        with open(output, "w") as toml_file:
            toml.dump(res, toml_file)


def make_no_TCS():
    ENV = Env()

    nodes = [ThermalNode(0, 'North',
                              {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.22}, Material().aluminium()),
             ThermalNode(1, '+V',
                              {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(2, 'South',
                              {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(3, '-V',
                              {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003, 0], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(4, 'Zenith',
                              {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0, 0], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(5, 'Nadir',
                              {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0, 0, 0,0,0,0,0, 0.003], 'mass': 0.62, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, Material().aluminium()),
             ThermalNode(6, 'battery',
                              {'area': 0.35 * 0.36, 'contact_area': [0,0,0,0.2*0.2*0.5,0,0,0], 'mass': 2.5, 'heat_generated': 5,
                               'radiation_area': 0.2*0.2*6}, Material().battery()),
             ThermalNode(7, 'solar',
                              {'area': 0.63, 'contact_area': [], 'mass': 0.5, 'heat_generated': 0,
                               'radiation_area': 0.63*2}, Material().solar_panel()),
             ThermalNode(8, 'OBC',
                              {'area': 0.1, 'contact_area': [0,0,0,0.1*0.1*0.5,0,0,0], 'mass': 0.25, 'heat_generated': 5,
                               'radiation_area': 0.1*0.1*3}, Material().PCB()),
             ThermalNode(9, 'propellant',
                              {'area': 0.1, 'contact_area': [0,0,0,0,0,0.3*0.1,0,0,0], 'mass': 0.25, 'heat_generated': 0,
                               'radiation_area': 0.1}, Material().aluminium()),
             ThermalNode(10, 'DST box',
                              {'area': 0.11, 'contact_area': [0,0,0,0,0,0.19*0.12,0,0,0,0,0], 'mass': 5, 'heat_generated': 2,
                               'radiation_area': 0.11}, Material().aluminium()),
             ThermalNode(11, 'DST baffle',
                             {'area': 0.36/2, 'contact_area': [0,0,0,0,0,0.55], 'mass': 25, 'heat_generated': 0,
                              'radiation_area': 0.36*2}, Material().aluminium())]

    connections = np.array([[1,3,4,5],
                            [0,2,4,5],
                            [1,3,4,5],
                            [0,2,4,5],
                            [0,1,2,3],
                            [0,1,2,3],
                            [3],
                            [],
                            [3],
                            [5],
                            [5],
                            [5]])

    TM = ThermalModel(nodes, connections, ENV, 10 * ENV.t_orbit, [350] * len(nodes))
    TM.solve()
    TM.plot([6,7,8,9,10,11], with_legend=True)
    for i, node in enumerate(TM.nodes):
        print(node.name, np.max(TM.solution.y[i][1000:]), np.min(TM.solution.y[i][1000:]))




if __name__ == '__main__':
    #make_no_TCS()
    make_TCS('result.toml')
