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
import matplotlib.pyplot as plt
from DSE.Thermal.env import Env
from DSE.Thermal.materials import *
from DSE.Thermal.thermal_dyn_sim import ThermalNode, ThermalModel
from DSE.Thermal.esatan_reader import heat_flows, interpolate_points, prepare_heat_flows


def make_TCS(output=None, plot=False, unit='K', Q=None, file=None):
    ENV = Env()

    nodes = [ThermalNode(0, 'North',
                              {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.3, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.22}, aluminium(), white_paint()),
             ThermalNode(1, 'East',
                              {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003], 'mass': 0.3, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.22}, aluminium(), white_paint()),
             ThermalNode(2, 'South',    # radiator mounted on this face
                              {'area': 0.36 * 0.11, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.3, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.11}, aluminium(), white_paint()),
             ThermalNode(3, 'West',
                              {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003, 0], 'mass': 0.3, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.22}, aluminium(), white_paint()),
             ThermalNode(4, 'Zenith',    # this is actually Zenith, leave as structure for plotting purposes
                              {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0, 0], 'mass': 0.3, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, aluminium(), coated_solar_panel()),
             ThermalNode(5, 'Nadir',
                              {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003,0,0,0,0,0,0,0,0,0,0,0,0.1*0.1], 'mass': 0.3, 'heat_generated': 0,
                               'radiation_area': 0.35 * 0.36}, aluminium(), white_paint()),
             ThermalNode(6, 'battery',
                              {'area': 0, 'contact_area': [0,0,0,0,0,0,0,0,0.2*0.2], 'mass': 2.5, 'heat_generated': 4,
                               'radiation_area': 0.2*0.2*5}, battery(), MLI_coat()),
             ThermalNode(7, 'solar_panel_1',
                              {'area': 0.35*0.66, 'contact_area': [], 'mass': 1.5, 'heat_generated': 0,
                               'radiation_area': 0.35*0.66*2}, solar_panel(), coated_solar_panel()),
             ThermalNode(8, 'OBC',
                              {'area': 0, 'contact_area': [0,0,0,0,0.1*0.1*0.5,0,0], 'mass': 0.25, 'heat_generated': 5,
                               'radiation_area': 0.1*0.1*3}, PCB(), MLI_coat()),
             ThermalNode(9, 'propellant',
                              {'area': 0, 'contact_area': [], 'mass': 5, 'heat_generated': 0,
                               'radiation_area': 0.1}, hastealloy(), MLI_coat()),
             ThermalNode(10, 'radiator',
                              {'area': 0.35*0.11, 'contact_area': [], 'mass': 0.6, 'heat_generated': 0,
                               'radiation_area': 0.35*0.11}, OSR_mat(), OSR_coat()),
             ThermalNode(11, 'DST box',
                              {'area': 0, 'contact_area': [0,0,0,0,0,0.19*0.12,0,0,0,0,0.01,0.01, 0.01], 'mass': 5, 'heat_generated': 2,
                               'radiation_area': 0.11}, aluminium(), MLI_coat()),
             ThermalNode(12, 'DST baffle',
                             {'area': 0.37, 'contact_area': [0,0,0,0,0,0.1], 'mass': 25, 'heat_generated': 0,
                              'radiation_area': 0.36*2}, aluminium(), white_paint()),
             ThermalNode(13, 'solar_panel_2',
                         {'area': 0.5*0.35, 'contact_area': [], 'mass': 1.2, 'heat_generated': 0,
                          'radiation_area': 0.5*0.35*2}, solar_panel(), coated_solar_panel()),
             ThermalNode(14, 'solar_panel_3',
                         {'area': 0.5*0.35, 'contact_area': [], 'mass': 1.2, 'heat_generated': 0,
                          'radiation_area': 0.5*0.35*2}, solar_panel(), coated_solar_panel()),
             ThermalNode(15, 'propulsion',
                         {'area': 0, 'contact_area': [0,0,0,0.01,0,0.01], 'mass': 1.7, 'heat_generated': 2,
                          'radiation_area': 0.08*0.08*5}, low_conductance_tape(), bare_Al()),
             ThermalNode(16, 'pipes',
                         {'area': 0, 'contact_area': [], 'mass': 0.45, 'heat_generated': 0.8,
                          'radiation_area': 0.016}, aluminium(), MLI_coat())]

    connections = np.array([[1,3,4,5],  # 0
                            [0,2,4,5],  # 1
                            [1,3,4,5],  # 2
                            [0,2,4,5],  # 3
                            [0,1,2,3],  # 4
                            [0,1,2,3],  # 5
                            [8],        # 6
                            [],         # 7
                            [4],        # 8
                            [],         # 9
                            [],         # 10
                            [10,12],    # 11
                            [5],        # 12
                            [],         # 13
                            [],         # 14
                            [3,5],      # 15
                            []])        # 16

    TM = ThermalModel(nodes, connections, ENV, [270] * len(nodes), n_orbits=10, unit=unit, ESATAN=True, Q_ESATAN=Q)
    TM.solve()
    if plot:
        TM.plot([0,1,2,3,4,5,6,7,8,10,11,12,13,14,15], with_legend=True, save=file)
        TM.plotting_Q()

    if output:
        res = {node.name: {'max': np.max(TM.solution.y[i][3000:]), 'min': np.min(TM.solution.y[i][3000:])}
               for i, node in enumerate(TM.nodes)}
        with open(output, "w") as toml_file:
            toml.dump(res, toml_file)

    return TM

# def make_no_TCS(output=None):
#     ENV = Env()
#
#     nodes = [ThermalNode(0, 'North',
#                               {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.22}, Material().aluminium()),
#              ThermalNode(1, '+V',
#                               {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(2, 'South',
#                               {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(3, '-V',
#                               {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003, 0], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(4, 'Zenith',
#                               {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0, 0], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(5, 'Nadir',
#                               {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0, 0, 0,0,0,0,0, 0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, Material().aluminium()),
#              ThermalNode(6, 'battery',
#                               {'area': 0.35 * 0.36, 'contact_area': [0,0,0,0.2*0.2*0.5,0,0,0], 'mass': 2.5, 'heat_generated': 5,
#                                'radiation_area': 0.2*0.2*6}, Material().battery()),
#              ThermalNode(7, 'solar',
#                               {'area': 0.63, 'contact_area': [], 'mass': 0.5, 'heat_generated': 0,
#                                'radiation_area': 0.63*2}, Material().solar_panel()),
#              ThermalNode(8, 'OBC',
#                               {'area': 0.1, 'contact_area': [0,0,0,0.1*0.1*0.5,0,0,0], 'mass': 0.25, 'heat_generated': 5,
#                                'radiation_area': 0.1*0.1*3}, Material().PCB()),
#              ThermalNode(9, 'propellant',
#                               {'area': 0.1, 'contact_area': [0,0,0,0,0,0.3*0.1,0,0,0], 'mass': 0.25, 'heat_generated': 0,
#                                'radiation_area': 0.1}, Material().aluminium()),
#              ThermalNode(10, 'DST box',
#                               {'area': 0.11, 'contact_area': [0,0,0,0,0,0.19*0.12,0,0,0,0,0], 'mass': 5, 'heat_generated': 2,
#                                'radiation_area': 0.11}, Material().aluminium()),
#              ThermalNode(11, 'DST baffle',
#                              {'area': 0.36/2, 'contact_area': [0,0,0,0,0,0.55], 'mass': 25, 'heat_generated': 0,
#                               'radiation_area': 0.36*1.5}, Material().aluminium())]
#
#     connections = np.array([[1,3,4,5],
#                             [0,2,4,5],
#                             [1,3,4,5],
#                             [0,2,4,5],
#                             [0,1,2,3],
#                             [0,1,2,3],
#                             [3],
#                             [],
#                             [3],
#                             [5],
#                             [5],
#                             [5]])
#
#     TM = ThermalModel(nodes, connections, ENV, 10 * ENV.t_orbit, [350] * len(nodes))
#     TM.solve()
#     # TM.plot([6,7,8,9,10,11], with_legend=True, save=True)
#     if output:
#         res = {node.name: {'max': np.max(TM.solution.y[i][1000:]), 'min': np.min(TM.solution.y[i][1000:])}
#                for i, node in enumerate(TM.nodes)}
#         with open(output, "w") as toml_file:
#             toml.dump(res, toml_file)
#
#     return TM
#
# def make_no_rad_TCS(output=None):
#     ENV = Env()
#
#     nodes = [ThermalNode(0, 'North',
#                               {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.22}, Material().white_coating()),
#              ThermalNode(1, '+V',
#                               {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(2, 'South',
#                               {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(3, '-V',
#                               {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003, 0], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(4, 'Zenith',
#                               {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0, 0], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(5, 'Nadir',
#                               {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0,0,0,0,0,0,0,0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, Material().white_coating()),
#              ThermalNode(6, 'battery',
#                               {'area': 0.35 * 0.36, 'contact_area': [0,0,0,0.2*0.2*0.5], 'mass': 2.5, 'heat_generated': 5,
#                                'radiation_area': 0.2*0.2*6}, Material().battery()),
#              ThermalNode(7, 'solar',
#                               {'area': 0.63, 'contact_area': [], 'mass': 0.5, 'heat_generated': 0,
#                                'radiation_area': 0.63*2}, Material().coated_solar_panel()),
#              ThermalNode(8, 'OBC',
#                               {'area': 0.1, 'contact_area': [0,0,0,0.1*0.1*0.5], 'mass': 0.25, 'heat_generated': 5,
#                                'radiation_area': 0.1*0.1*3}, Material().PCB()),
#              ThermalNode(9, 'propellant',
#                               {'area': 0.1, 'contact_area': [0,0,0,0,0,0.3*0.1], 'mass': 0.25, 'heat_generated': 0,
#                                'radiation_area': 0.1}, Material().MLI()),
#              ThermalNode(10, 'DST box',
#                               {'area': 0.11, 'contact_area': [0,0,0,0,0,0.19*0.12,0,0,0,0,0], 'mass': 5, 'heat_generated': 2,
#                                'radiation_area': 0.11}, Material().MLI()),
#              ThermalNode(11, 'DST baffle',
#                              {'area': 0.36/2, 'contact_area': [0,0,0,0,0,0.55], 'mass': 25, 'heat_generated': 0,
#                               'radiation_area': 0.36*1.5}, Material().white_coating())]
#
#     connections = np.array([[1,3,4,5],
#                             [0,2,4,5],
#                             [1,3,4,5],
#                             [0,2,4,5],
#                             [0,1,2,3],
#                             [0,1,2,3],
#                             [3],
#                             [],
#                             [3],
#                             [5],
#                             [5],
#                             [5]])
#
#     TM = ThermalModel(nodes, connections, ENV, 10 * ENV.t_orbit, [270] * len(nodes))
#     TM.solve()
#     # TM.plot([6,7,8,9,10,11], with_legend=True, save=True)
#
#     if output:
#         res = {node.name: {'max': np.max(TM.solution.y[i][1000:]), 'min': np.min(TM.solution.y[i][1000:])}
#                for i, node in enumerate(TM.nodes)}
#         with open(output, "w") as toml_file:
#             toml.dump(res, toml_file)
#
#     return TM
#
# def make_no_sun_sheild(output=None):
#     ENV = Env()
#
#     nodes = [ThermalNode(0, 'North',
#                               {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.22}, aluminium(), white_paint()),
#              ThermalNode(1, '+V',
#                               {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, aluminium(), white_paint()),
#              ThermalNode(2, 'South',
#                               {'area': 0.36 * 0.22, 'contact_area': [0, 0.22*0.003, 0, 0.22*0.003, 0.36*0.003, 0.36*0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, aluminium(), white_paint()),
#              ThermalNode(3, '-V',
#                               {'area': 0.35 * 0.22, 'contact_area': [0.22*0.003, 0, 0.22*0.003, 0, 0.36*0.003, 0.36*0.003, 0], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, aluminium(), white_paint()),
#              ThermalNode(4, 'Zenith',
#                               {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0, 0], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, aluminium(), white_paint()),
#              ThermalNode(5, 'Nadir',
#                               {'area': 0.35 * 0.36, 'contact_area': [0.22*0.003, 0.22*0.003, 0.22*0.003, 0.22*0.003, 0, 0, 0,0,0,0,0,0, 0.003], 'mass': 0.62, 'heat_generated': 0,
#                                'radiation_area': 0.35 * 0.36}, aluminium(), white_paint()),
#              ThermalNode(6, 'battery',
#                               {'area': 0.35 * 0.36, 'contact_area': [0,0,0,0.2*0.2*0.5,0,0,0], 'mass': 2.5, 'heat_generated': 5,
#                                'radiation_area': 0.2*0.2*6},battery(), MLI_coat()),
#              ThermalNode(7, 'solar',
#                               {'area': 0.63, 'contact_area': [], 'mass': 0.5, 'heat_generated': 0,
#                                'radiation_area': 0.63*2}, solar_panel(), coated_solar_panel()),
#              ThermalNode(8, 'OBC',
#                               {'area': 0.1, 'contact_area': [0,0,0,0.1*0.1*0.5,0,0,0], 'mass': 0.25, 'heat_generated': 5,
#                                'radiation_area': 0.1*0.1*3}, PCB(), bare_Al()),
#              ThermalNode(9, 'propellant',
#                               {'area': 0.1, 'contact_area': [0,0,0,0,0,0.3*0.1,0,0,0], 'mass': 0.25, 'heat_generated': 0,
#                                'radiation_area': 0.1}, aluminium(), MLI_coat()),
#              ThermalNode(10, 'radiator',
#                               {'area': 0.35*0.22, 'contact_area': [0,0,0,0,0,0,0,0,0,0,0,0.01], 'mass': 0.6, 'heat_generated': 0,
#                                'radiation_area': 0.35*0.22}, OSR_mat(), OSR_coat()),
#              ThermalNode(11, 'DST box',
#                               {'area': 0.11, 'contact_area': [0,0,0,0,0,0.19*0.12,0,0,0,0,0.01,0], 'mass': 5, 'heat_generated': 2,
#                                'radiation_area': 0.11}, aluminium(), MLI_coat()),
#              ThermalNode(12, 'DST baffle',
#                              {'area': 0.36/2, 'contact_area': [0,0,0,0,0,0.55], 'mass': 25, 'heat_generated': 0,
#                               'radiation_area': 0.36*1.5}, aluminium(), white_paint())]
#
#     connections = np.array([[1,3,4,5],
#                             [0,2,4,5],
#                             [1,3,4,5],
#                             [0,2,4,5],
#                             [0,1,2,3],
#                             [0,1,2,3],
#                             [3],
#                             [],
#                             [3],
#                             [5],
#                             [],
#                             [5,10],
#                             [5]])
#
#     TM = ThermalModel(nodes, connections, ENV, [270] * len(nodes), n_orbits=10 ,sun_shield=True, ESATAN=True)
#     TM.solve()
#     TM.plot(with_legend=True, save=True)
#
#     if output:
#         res = {node.name: {'max': np.max(TM.solution.y[i][1000:]), 'min': np.min(TM.solution.y[i][1000:])}
#                for i, node in enumerate(TM.nodes)}
#         with open(output, "w") as toml_file:
#             toml.dump(res, toml_file)
#
#     return TM
#
# def nice_plots(save=False):
#     TM1 = make_no_TCS("no_TCS.toml")
#     TM2 = make_no_rad_TCS("no_rad_TCS.toml")
#     TM3 = make_no_sun_sheild("no_sun_sheild.toml")
#     TM4 = make_TCS("TCS.toml")
#
#     I = [4,6,7,8,9,10,11,12]
#
#     fig, ax = plt.subplots(2, 2, figsize=(11, 10))
#     for i, node in enumerate(TM1.nodes):
#         if i in I: ax[0, 0].plot(TM1.solution.t, TM1.solution.y[i], label=node.name)
#     for i, node in enumerate(TM2.nodes):
#         if i in I:  ax[0, 1].plot(TM2.solution.t, TM2.solution.y[i], label=node.name)
#     for i, node in enumerate(TM3.nodes):
#         if i in I:  ax[1, 0].plot(TM3.solution.t, TM3.solution.y[i], label=node.name)
#     for i, node in enumerate(TM4.nodes):
#         if i in I:  ax[1, 1].plot(TM4.solution.t, TM4.solution.y[i], label=node.name)
#
#     ax[0, 0].set_title('a) No TCS')
#     ax[0, 1].set_title('b) No radiator')
#     ax[1, 0].set_title('c) No sun sheild')
#     ax[1, 1].set_title('d) Full TCS')
#
#     for i, j in [(0,0), (0,1), (1,0), (1,1)]:
#         ax[i, j].set_xlabel('Time [s]')
#         ax[i, j].set_ylabel('Temperature [K]')
#         # ax[i, j].set_xlim(6 * TM1.env.t_orbit, 9 * TM1.env.t_orbit)
#         # ax[i, j].set_xticks(np.arange(8 * TM1.env.t_orbit, 9 * TM1.env.t_orbit, 1000), np.arange(0, int(TM1.env.t_orbit), 1000))
#
#     ax[1, 1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#
#     fig.tight_layout()
#     if save:
#         plt.savefig('nice_plots_1.svg')
#     plt.show()

if __name__ == '__main__':
    AN15 = prepare_heat_flows(interpolate_points(heat_flows(), Env().t_orbit, 10))

    for i, arr in enumerate(AN15):
        plt.plot(arr, label=f'{i}')
        plt.xlim(8*Env().t_orbit, 9*Env().t_orbit)
        plt.legend(loc='center right')
    plt.show()
    plt.close()

    make_TCS('results/TCS1.toml', plot=True, Q=AN15, file='results/TCS15.png')
