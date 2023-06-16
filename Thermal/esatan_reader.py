import numpy as np
import os
import matplotlib.pyplot as plt
from DSE.Thermal.env import Env


def heat_flows():
    """
    Returns the heat flows for the thermal simulation for each outer face of the satellite.
    Q is a dictionary with the names of the faces as keys and a 2d array of time and heat flows as values.
    """

    files = os.listdir('ESATAN')
    Q = {}
    for file in files:
        arr = np.genfromtxt('ESATAN/'+file, delimiter=';')
        arr = np.clip(arr, 0, np.max(arr))
        Q[file[:-4]] = arr

    return Q


def interpolate_points(Q: dict, orbit_time: int, n: int) -> dict:
    """
    Interpolates the data points to 1 second intervals
    """

    Q_new = {}
    for key in Q.keys():
        Q_new[key] = np.interp(np.arange(0, orbit_time, 1), Q[key][:, 0], Q[key][:, 1])
        if n > 1: Q_new[key] = np.tile(Q_new[key], n)

    return Q_new


def prop_heat():
    z = np.zeros(15442)
    h = np.ones(845) * 60
    res = np.hstack([np.tile(np.hstack([z,h]), 3), np.zeros(5429)])
    return res


def prepare_heat_flows(Q: dict):
    """
    return list of heat flows for each face in order
    """
    Q_res = {}

    Q_res['North'] = Q['north']
    Q_res['East'] = Q['east']
    Q_res['South'] = Q['south']
    Q_res['West'] = Q['west']
    Q_res['Zenith'] = Q['zenith']
    Q_res['Nadir'] = Q['nadir']
    Q_res['battery'] = np.zeros(len(Q['north']))
    Q_res['solar_panel_1'] = Q['back_solar_top'] + Q['back_solar_bottom']
    Q_res['OBC'] = np.zeros(len(Q['north']))
    Q_res['propellant'] = np.zeros(len(Q['north']))
    Q_res['radiator'] = Q['radiator'] * 2
    Q_res['DST box'] = np.zeros(len(Q['north']))
    Q_res['DST baffle'] = Q['baffle_inner'] + Q['baffle_outer']
    Q_res['solar_panel_2'] = Q['side_solar_1_top'] + Q['side_solar_1_bottom']
    Q_res['solar_panel_3'] = Q['side_solar_2_top'] + Q['side_solar_2_bottom']
    Q_res['propulsion'] = Q['west'] * 0.2 + prop_heat()
    Q_res['pipes'] = np.zeros(len(Q['north']))

    return list(Q_res.values())

def plot_ESATAN_heat_values():
    AN15 = prepare_heat_flows(interpolate_points(heat_flows(), Env().t_orbit, 10))

    l = ['North', 'East', 'South', 'West', 'Zenith', 'Nadir', 'Batteries', 'Solar panel 1', 'OBC', 'Propellant tanks',
         'Radiator', 'DST box', 'DST baffle', 'Solar panel 2', 'Solar panel 3', 'Thruster', 'Propellant lines']
    idx = [0,1,2,3,4,5,7,10,12,13,14]
    for i in idx:
        plt.plot(AN15[i], label=f'{l[i]}')
        plt.xlim(8*Env().t_orbit, 9*Env().t_orbit)
        plt.xticks(np.arange(8*Env().t_orbit, 9*Env().t_orbit, 1000), np.arange(0, int(Env().t_orbit), 1000))
        plt.legend(bbox_to_anchor=(1, 0.5), loc="center left")
        plt.subplots_adjust(right=0.74)
    plt.xlabel('Time [s]')
    plt.ylabel('Heat [W]')
    plt.title('ESATAN Heat Values')
    plt.savefig('plots/ESATAN_Q_vals')
    plt.close()

plot_ESATAN_heat_values()
