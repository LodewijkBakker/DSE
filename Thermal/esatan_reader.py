import numpy as np
import os


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
    Q_res['radiator'] = Q['radiator']
    Q_res['DST box'] = np.zeros(len(Q['north']))
    Q_res['DST baffle'] = Q['baffle_inner'] + Q['baffle_outer']
    Q_res['solar_panel_2'] = Q['side_solar_1_top'] + Q['side_solar_1_bottom']
    Q_res['solar_panel_3'] = Q['side_solar_2_top'] + Q['side_solar_2_bottom']
    Q_res['propulsion'] = Q['west'] * 0.2 + prop_heat()
    Q_res['pipes'] = np.zeros(len(Q['north']))

    return list(Q_res.values())
