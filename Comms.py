import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from tqdm import tqdm
import re
from datetime import datetime
from scipy import integrate
import pickle

# NOTE: Standard values:
delta           = 0.1       # [-]   Step size of elevation. 0.1 is good, anything above 1 is bad and lower than 0.1 is not needed
link_margin     = 3         # [dB]  As advised
boltzmann       = -228.6    # [dBW/K/Hz] Logarithmic representation
earth_radius    = 6378      # [km]
orbit_alt       = 300       # [km]  Orbital altitude    IF CHANGED, THE ORBITAL PASS DATA FILES MUST BE CHANGED AS WELL
inclination     = 96.7      # [deg] Orbital inclination IF CHANGED, THE ORBITAL PASS DATA FILES MUST BE CHANGED AS WELL

# NOTE: Design choices that stay constant:
frequency       = 14000  # [MHz] 8000 for X-Band
gs_gain         = 43.2  # [dB]  43.2 for 2.4m X-band GS
l_atm_zenith    = -0.3  # [dB]  -0.3 for X-Band
l_line          = -1    # [dB]  Line losses for antenna
l_line_gs       = -1    # [dB]  Based on delfi ground station
l_pol           = -3    # [dB]  Due to mismatch of polarization of channel and antenna. For circular -> linear, a 3 dB loss can be taken.
l_filter        = 0     # [dB]  If transmission power needs to be distributed to several antennae: ~1 dB loss per split in two
t_noise         = 180   # [K]   Based on EWI GS
coding_gain     = 0     # [dB]  Coding gain, when using the required Es/N0 from a datum datasheet, a coding gain of 0 should be taken
pl_data_rate    = 4.1402 * 0.2  # [Mb/s] Generated data that should be sent to earth

# NOTE: Standard calculations that only have to be done once
wavelength      = 299.792458 / frequency
orbit_radius    = earth_radius + orbit_alt
orbit_time      = 2 * np.pi * np.sqrt(((orbit_radius * 1000) ** 3) / (6.67430e-11 * 5.9736e+24))

#   Open the visibility times file and load the list (Visibility_times is a list of lists, with each sublist corresponding to a maximum pass elevation)
with open(f'Visibility_times\\Visibility_times_{delta}.txt', 'rb') as f:
    visibility_times = pickle.load(f)


#   Create antenna class, which contains all information needed on the antenna
class Antenna:
    def __init__(self, name: str, gains: list, power: float, gimballed: bool, rotating: bool):
        antenna_angles = [-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90]

        self.name = name
        self.gains = gains
        self.gains += reversed(gains[0:-1])
        self.gains_inter = interpolate.interp1d(antenna_angles, gains, kind='cubic')
        self.power = power
        self.gimballed = gimballed
        self.rotating = rotating


#   Create Modulation class, which contains all information needed on the modulation scheme
class Modulation:
    def __init__(self, name: str, code_rate: float, snr: float, bits_hz: float):
        self.name = name
        self.code_rate = code_rate
        self.snr = snr
        self.bits_hz = bits_hz


#   Create an antenna class for all the different antennas
omnidirectional = Antenna("Omnidirectional", [0, 0, 0, 0, 0, 0, 0], power=2, gimballed=False, rotating=False)  # TODO fix this bullshit
low_gain_patch = Antenna("Low gain patch", [-4.3, -2, 0.6, 3, 4.9, 6.1, 6.7], power=4, gimballed=True, rotating=True)
high_gain_patch = Antenna("High gain patch", [-15.7, -10.8, -6, -0.5, 6.3, 10.5, 11.5], power=3, gimballed=True, rotating=True)
comb_patch = Antenna("Combined patch", [-4.3, -2, 0.6, 3, 6.3, 10.5, 11.5], power=3, gimballed=True, rotating=True)
reflect_array = Antenna("Reflectarray", [-10, -10, -10, -10, -8, 4, 24], power=5, gimballed=True, rotating=True)
isoflux = Antenna("Isoflux", [-8, -3, 2, 5, 5.8, 3.7, 2.9], power=5, gimballed=False, rotating=True)
phase_array = Antenna("Phased array", [-20, -15, 15, 27, 30, 30, 30], power=1, gimballed=False, rotating=False)
parabolic_10 = Antenna("10cm-Parabolic", [-30, -25, -20, -8.6, 2.528, 12.24, 16.3], power=4, gimballed=True, rotating=True)
parabolic_5 = Antenna("5cm-Parabolic", [-14.6, -9, -3.5, 2.22, 6.22, 9.22, 10.22], power=4, gimballed=True, rotating=True)
# https://www.satsig.net/pointing/antenna-beamwidth-calculator.htm
# Use these if you want to find the most power efficient antenna (this means the actual data values are only to be used for comparing, since they are not representative.)
# omnidirectional = Antenna("Omnidirectional",    [0,     0,      0,      0,      0,      0,      0],     power=3,    gimballed=False,  rotating=False)
# low_gain_patch  = Antenna("Low gain patch",     [-4.3,  -2,     0.6,    3,      4.9,    6.1,    6.7],   power=3,    gimballed=True,   rotating=True)
# high_gain_patch = Antenna("High gain patch",    [-15.7, -10.8,  -6,     -0.5,   6.3,    10.5,   11.5],  power=3,    gimballed=True,   rotating=True)
# comb_patch      = Antenna("Combined patch",     [-4.3,  -2,     0.6,    3,      6.3,    10.5,   11.5],  power=3,    gimballed=True,   rotating=True)
# reflect_array   = Antenna("Reflectarray",       [-10,   -10,    -10,    -10,    -8,     4,      24],    power=3,    gimballed=False,  rotating=True)
# isoflux         = Antenna("Isoflux",            [-8,    -3, 	2,  	5,  	5.8,	3.7,	2.9],   power=3,    gimballed=False,  rotating=True)
# phase_array     = Antenna("Phased array",       [-20,   -20,	30,	    30, 	30, 	30, 	30],    power=3,    gimballed=False,  rotating=False)

# list of antennas to be used when plotting "all antennas"
antenna_list = [omnidirectional, low_gain_patch, high_gain_patch, reflect_array, isoflux, parabolic_10, parabolic_5, phase_array]


def slant_dist_calc(elev_angle: float) -> float:
    """
    Calculates the distance between the satellite and the ground station. Assumes a 300km orbit

    :param elev_angle: [deg] Elevation angle of satellite, as seen from the ground station
    :return: [km] distance between satellite and ground station
    """
    slant_dist = earth_radius * (
            (((orbit_radius ** 2) / (earth_radius ** 2) - (np.cos(elev_angle * np.pi / 180) ** 2)) ** 0.5) - np.sin(elev_angle * np.pi / 180))

    return slant_dist


def EIRP_calc(power: float, gain: float) -> float:
    """
    Calculates the equivalent isotropically radiated power from the antenna. Takes into account 1dB of line losses.

    :param power: [W] output power of antenna (NOT THE SAME AS INPUT POWER)
    :param gain: [-] The gain of the antenna for this specific angle
    :return: EIRP [dB]
    """
    power_dbw = 10 * np.log10(power)
    return power_dbw + l_line + l_filter + gain


def link_budget(elev_angle: float, power: float, gain: float, modulation: Modulation = Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74),
                print_values: bool = False, bandwidth: float = 100) -> float:
    """
    Calculates the maximum data rate possible for the inserted values.
    Assumes an X-Band communication (8000MHz), 4dB of losses due to line and polarization losses at ground station, 2.4m GS, and a 180K system noise temperature.

    :param bandwidth: [MHz] Maximum available bandwidth
    :param elev_angle: [deg] Elevation angle of the satellite, as seen from the ground station
    :param power: [W] Output power of the antenna
    :param gain: [dB] Gain of the antenna
    :param modulation: Mudulation class
    :param print_values: Print interesting values from link budget calculation
    :return: Maximum data-rate [bit/s]
    """

    #   Modulation dependent: # NOTE: DVB-S2X - ESA ACM (https://datumsystems.com/product/dvbs2x/)
    req_snr = modulation.snr  # [dB]
    code_rate = modulation.code_rate  # [-]

    max_bandwidth = bandwidth  # [MHz]
    max_data_rate = max_bandwidth * modulation.bits_hz * 1000000
    # print(modulation.name, "max data rate = ", max_data_rate)

    slant_dist = slant_dist_calc(elev_angle)
    eirp = EIRP_calc(power, gain)
    l_path = -20 * np.log10(4 * np.pi * slant_dist * (1000 / wavelength))
    l_atm = l_atm_zenith / np.sin(elev_angle * np.pi / 180)
    l_total_channel = l_atm + l_path

    signal_level_gs = eirp + l_total_channel
    t_noise_db = 10 * np.log10(t_noise)
    snr = link_margin + req_snr - coding_gain
    snr_power_density = signal_level_gs + gs_gain + l_pol + l_line_gs - t_noise_db - boltzmann
    data_rate = 10 ** ((snr_power_density - snr) / 10)
    data_rate *= code_rate
    if print_values:
        print("slant_dist       =", slant_dist)
        print("eirp             =", eirp)
        print("wavelength       =", wavelength)
        print("l_path           =", l_path)
        print("l_atm            =", l_atm)
        print("l_total_channel  =", l_total_channel)
        print("snr              =", snr)
        print("snr_power_density=", snr_power_density)
        print("data_rate        =", data_rate)

    if data_rate > max_data_rate:
        data_rate = max_data_rate

    return data_rate


def data_vs_elev(max_elev: int, antenna: Antenna, rotation_axes: int, modulation: Modulation = Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), bandwidth: float = 100) -> list:
    """
    Calculates the maximum data rate in Mbit per second for all elevation angles (steps of 1 deg), starting at 5 deg.
    For example: if max_elev=7, then it will generate: [DR@5deg, DR@6deg, DR@7deg, DR@6deg, DR@5deg]
    :param bandwidth: [MHz] Maximum available bandwidth
    :param max_elev: [deg] maximum elevation of the pass
    :param antenna: Antenna class
    :param rotation_axes: 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)
    :param modulation: Modulation class
    :return: [Mb/s] list of data rates for each elevation angle. (back and forth)
    """
    power = antenna.power
    gain_inter_func = antenna.gains_inter

    elev_angles = np.hstack((np.arange(5, max_elev + 1), np.flip(np.arange(5, max_elev))))  # same as elev_angles = [*range(5, max_elev + 1, 1)] + [*range(max_elev - 1, 4, -1)]

    data_rates = np.zeros(len(elev_angles))
    for i in range(0, len(elev_angles)):
        if rotation_axes == 2:
            gain = gain_inter_func(0)
        elif rotation_axes == 1:
            raise NotImplementedError("single axis rotation has not been implemented yet")
        elif rotation_axes == 0:
            slant_dist = slant_dist_calc(elev_angles[i])
            nadir_angle = np.arccos(((slant_dist ** 2) + (orbit_radius ** 2) - (earth_radius ** 2)) / (2 * slant_dist * orbit_radius)) * 180 / np.pi
            gain = gain_inter_func(nadir_angle)
        else:
            raise ValueError("Invalid rotation axes value. Please insert 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)")

        data_rate = link_budget(elev_angles[i], power, gain, modulation=modulation, bandwidth=bandwidth) / 1000000
        data_rates[i] = data_rate

    return data_rates


def vis_time_calc(max_elev: int, min_elev: int) -> float:  # https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=805436
    """
    Calculates the time between min_elev and min_elev after passing the assigned maximum elevation.
    :param max_elev: [deg] maximum elevation of the pass
    :param min_elev: [deg] elevation to start and stop counting at
    :return: [sec] time between minimum elevation and minimum elevation again whilst passing over the maximum elevation
    """
    a = 2 / (np.deg2rad(360 / orbit_time) - np.deg2rad(360 / 86400) * np.cos(np.deg2rad(inclination)))
    b = np.arccos(np.cos(np.arccos((earth_radius / orbit_radius) * np.cos(np.deg2rad(min_elev))) - np.deg2rad(min_elev)) /
                  np.cos(np.arccos((earth_radius / orbit_radius) * np.cos(np.deg2rad(max_elev))) - np.deg2rad(max_elev)))
    return a * b


def data_vs_time(max_elev: int, antenna: Antenna, rotation_axes: int, modulation: Modulation = Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74),
                 plot: bool = True, bandwidth: float = 100) -> list:
    """
    Creates a list of time stamps and their respective data rates from a pass

    :param bandwidth: [MHz] Maximum available bandwidth
    :param max_elev: [deg] Maximum elevation of the pass
    :param antenna: Antenna class
    :param rotation_axes: 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)
    :param modulation: Modulation class
    :param plot: If the maximum data rate vs the elevation and the maximum data rate vs time should be plotted
    :return: times [secs], data_rate_vs_time [Mb/s]: list of times for each step of elevation, list of data rates according to the times
    """
    data_elev = data_vs_elev(max_elev, antenna, rotation_axes, modulation=modulation, bandwidth=bandwidth)
    elevation = np.arange(-max_elev + 5, max_elev - 4, 1)

    elev_new = np.arange(-max_elev + 5, max_elev - 5 + delta, delta)
    if plot:
        elev_data_inter = interpolate.interp1d(elevation, data_elev, kind='linear', fill_value=(data_elev[0], data_elev[-1]), bounds_error=False)
        plt.plot(elevation, data_elev, 'o')
        plt.plot(elev_new, elev_data_inter(elev_new), '--', label=antenna.name)
        plt.xlabel("Elevation [deg]")
        plt.ylabel("Max Data-rate [Mbit/s]")
        plt.show()

    times = visibility_times[int(max_elev) - 7]
    time_elev_inter = interpolate.interp1d(times, elev_new, kind='linear', fill_value=(times[0], times[-1]), bounds_error=False)

    if plot:
        plt.xlabel("Time passed [s]")
        plt.ylabel("Elevation [deg]")
        plt.plot(times, elev_new, 'o')
        plt.plot(times, time_elev_inter(times), '--')
        plt.show()

    # Calculate data rates vs time
    data_rate_vs_time = np.interp(elev_new, elevation, data_elev)

    return times, data_rate_vs_time


def total_data_over_pass(max_elev: int, antenna: Antenna, rotation_axes: int, modulation: Modulation = Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74),
                         bandwidth: float = 100) -> float:
    """
    Calculates the total data that could be send over a specific pass

    :param bandwidth: [MHz] Maximum available bandwidth
    :param max_elev: [deg] Maximum elevation of the pass
    :param antenna: Antenna class
    :param rotation_axes: 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)
    :param modulation: Modulation class
    :return: [Mbit] Total data sent over the pass
    """
    times, data_rates = data_vs_time(max_elev, antenna, rotation_axes, modulation, plot=False, bandwidth=bandwidth)
    return integrate.simpson(data_rates, times)


def pass_info(lat: int, data_type: str) -> list:
    """
    Gives the pass info by reading the files stored in the "Orbital pass data" folder
    :param lat: [deg] latitude of the ground station, should be 0, 30, 60, 78 or 90 degrees
    :param data_type: "elevation" or "time"
    :return: list with either maximum elevation or time of each pass
    """
    if lat not in [0, 30, 60, 78, 90]:
        raise ValueError("Latitude should be 0, 30, 60, 78 or 90 degrees")
    filename = f"90d_lat{lat}"
    with open("Orbital pass data/" + filename + ".txt") as file:
        count = 0
        passes = []
        for line in file.readlines():
            count += 1
            if count >= 13 and line != "\n":
                line = re.sub(' +', ' ', line)
                passes.append(line.split(" "))

    data = []
    for i in range(0, len(passes), 3):
        # print(passes[i + 1][4])
        if data_type == "elevation":
            data.append(float(passes[i + 1][4]))
        elif data_type == "time":
            start = passes[i][0] + " " + passes[i][1]
            end = passes[i + 2][0] + " " + passes[i + 2][1]
            date_time_start = datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
            date_time_end = datetime.strptime(end, "%Y-%m-%d %H:%M:%S")

            data.append([date_time_start, date_time_end])
        else:
            raise ValueError("data type should either be elevation or time")

    return data


def total_data(lat: int, antenna: Antenna, rotation_axes: int, modulation: Modulation = Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), bandwidth: float = 100) -> float:
    """
    Calculates the maximum amount data sent over 90 days
    :param bandwidth: [MHz] Maximum available bandwidth
    :param modulation: Modulation class
    :param lat: [deg] latitude of the ground station
    :param antenna: Antenna class
    :param rotation_axes: 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)
    :return: [Gb] maximum amount of data sent by antenna over 90 days
    """
    elevations = pass_info(lat, "elevation")
    total = 0
    for max_elev in tqdm(elevations, f"{antenna.name}"):
        if max_elev > 7:
            total += total_data_over_pass(int(max_elev), antenna, rotation_axes, modulation=modulation, bandwidth=bandwidth)
    print(f"{antenna.name} total data = {total / 1000}Gb")
    time.sleep(0.1)
    return total


def read_modulation_file() -> list:
    """
    Reads the modulations file
    :return: list of all modulation classes
    """
    with open("Modulations.csv") as file:
        lines = file.readlines()
        modulation_list = []

        # Only iteration on modulations
        for i, line in enumerate(lines):
            if line == 'ModCod,Code Rate,Es/N0,Bits/Hz\n':
                pass
            else:
                line = line.split(",")
                name = line[0]
                code_rate = float(line[1])
                snr = float(line[2])
                bits_hz = float(line[3][:-1])

                # Class creation, modulation
                mod = Modulation(name, code_rate, snr, bits_hz)
                modulation_list.append(mod)
    return modulation_list


def plot_gains():
    """
    plots the gains of all antennas in antenna_list
    :return: Shows plot
    """
    angles = np.arange(-90, 91, 1)
    antenna_angles = [-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90]
    for antenna in antenna_list:
        plt.plot(angles, antenna.gains_inter(angles), label=antenna.name)
        plt.plot(antenna_angles, antenna.gains, 'o')
    plt.title(f"Antenna gain")
    plt.xlabel(f'Antenna angle offset [deg]')
    plt.ylabel(f'Gain [dB]')
    plt.legend()
    plt.show()


def plot_data_vs_elev(max_elev: int, antenna: Antenna, rotation_axes: int, plot_label: str = '', modulation: Modulation = Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74),
                      bandwidth: float = 100):
    """
    Plots data vs elevation graph for a specific pass, antenna and modulation scheme
    :param bandwidth: [MHz] Maximum available bandwidth
    :param max_elev: [deg] Maximum elevation of the pass
    :param antenna: Antenna class
    :param rotation_axes: 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)
    :param plot_label: string that will be added to each plot in the legend. For example: f"{antenna.name} with {modulation.name} as modulation"
    :param modulation: Modulation class
    :return: plot (not shown)
    """
    x = [*range(-max_elev + 5, max_elev - 4, 1)]
    data_rates = data_vs_elev(max_elev=max_elev, antenna=antenna, rotation_axes=rotation_axes, modulation=modulation, bandwidth=bandwidth)
    plt.plot(x, data_rates, label=plot_label)
    plt.xlabel(f"Angle from maximum elevation [deg]")
    plt.ylabel(f"Data rate [Mb/s]")
    plt.title(f'Maximum data rate vs elevation angle')


def plot_data_vs_elev_all_antennas(max_elev: int, modulation: Modulation = Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), bandwidth: float = 100):
    """
    Plots data vs elevation graph for all antennas in antenna_list
    :param bandwidth: [MHz] Maximum available bandwidth
    :param max_elev: [deg] Maximum elevation of the pass
    :param modulation: Modulation class
    :return: Shows plots
    """
    for antenna in antenna_list:
        if antenna.gimballed:
            plot_data_vs_elev(max_elev=max_elev, antenna=antenna, rotation_axes=2, plot_label=antenna.name, modulation=modulation, bandwidth=bandwidth)
        else:
            plot_data_vs_elev(max_elev=max_elev, antenna=antenna, rotation_axes=0, plot_label=antenna.name, modulation=modulation, bandwidth=bandwidth)

    plt.xlabel(f"Angle from maximum elevation [deg]")
    plt.ylabel(f"Data rate [Mb/s]")
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.show()


def plot_data_vs_elev_all_modulations(max_elev: int, rotation_axes: int, antenna: Antenna, bandwidth: float = 100):
    """
    Plots data vs elevation graph for all modulation schemes
    :param bandwidth: [MHz] Maximum available bandwidth
    :param max_elev: [deg] Maximum elevation of the pass
    :param rotation_axes: 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)
    :param antenna: Antenna class
    :return: Shows plots
    """
    modulation_list = read_modulation_file()
    for modulation in modulation_list:
        plot_data_vs_elev(max_elev=max_elev, antenna=antenna, rotation_axes=rotation_axes, plot_label=modulation.name, modulation=modulation, bandwidth=bandwidth)
    plt.xlabel(f"Angle from maximum elevation [deg]")
    plt.ylabel(f"Data rate [Mb/s]")
    plt.legend()
    plt.show()


def plot_data_vs_time_all_antennas(max_elev: int, modulation: Modulation = Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), bandwidth: float = 100):
    """
    Generates data-rate vs time plots for all antennas in antenna_list
    :param bandwidth: [MHz] Maximum available bandwidth
    :param max_elev: [deg] Maximum elevation of the pass
    :param rotation_axes: 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)
    :param modulation: Modulation class
    :return: Shows plots
    """
    for antenna in antenna_list:
        if antenna.gimballed:
            times, data_rate = data_vs_time(max_elev=max_elev, antenna=antenna, rotation_axes=2, modulation=modulation, plot=False, bandwidth=bandwidth)
        else:
            times, data_rate = data_vs_time(max_elev=max_elev, antenna=antenna, rotation_axes=0, modulation=modulation, plot=False, bandwidth=bandwidth)
        plt.plot(times, data_rate, label=f"{antenna.name}, \n Total = {int(integrate.simpson(data_rate, times))}Mb")
    plt.xlabel("Time passed [s]")
    plt.ylabel("Max data-rate [Mbit/s]")
    plt.title("Maximum data-rate vs time")
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.show()


def plot_data_vs_time_all_elevations(antenna: Antenna, rotation_axes: int, modulation: Modulation = Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), bandwidth: float = 100):
    """
    Generates data-rate vs time plots for different maximum elevations
    :param bandwidth: [MHz] Maximum available bandwidth
    :param antenna: Antenna class
    :param rotation_axes: 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)
    :param modulation: Modulation class
    :return: Shows plots
    """
    for elevation in [90, 80, 70, 60, 50, 40, 30, 20]:
        times, data_rate = data_vs_time(max_elev=elevation, antenna=antenna, rotation_axes=rotation_axes, modulation=modulation, plot=False, bandwidth=bandwidth)
        plt.plot(times, data_rate, label=f"{elevation}°, Total = {int(integrate.simpson(data_rate, times))}Mb")
    plt.xlabel("Time passed [s]")
    plt.ylabel("Max Data-rate [Mbit/s]")
    plt.title(f"{antenna.name} Antenna")
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.show()


def plot_pass_elev_dist(lat: int):
    """
    Shows the distribution of passes over 90 days for a specific latitude of the ground station
    :param lat: [deg] latitude of the ground station
    :return: Shows plot
    """
    data = pass_info(lat, "elevation")
    filename = f"90d_lat{lat}"

    weights = np.ones_like(data) / len(data)
    plt.hist(data, bins=15, range=[5, 90], weights=weights)
    x_limits = plt.xlim()  # Get the x-limits of the histogram
    x_center = (x_limits[1] - x_limits[0]) / 2
    plt.title(f"Distribution of tracking passes over 90 days [Lat {lat}°]")
    plt.text(x_center, -6, f'Total number of passes: {int(len(data))}', size=12, ha='center')
    plt.savefig("Orbital_pass_graphs/hist_" + filename + "_normalized.png", dpi=400)
    plt.show()

    plt.hist(data, bins=15, range=[5, 90])
    plt.title(f"Distribution of tracking passes over 90 days [Lat {lat}°]")
    x_limits = plt.xlim()  # Get the x-limits of the histogram
    x_center = (x_limits[1] - x_limits[0]) / 2
    plt.text(x_center, -160, f'Total number of passes: {int(len(data))}', size=12, ha='center')
    # Adjust the plot layout
    plt.tight_layout()
    plt.savefig("Orbital_pass_graphs/hist_" + filename + ".png", dpi=400)
    plt.show()


def plot_accumulated_data(lat: int, antenna: Antenna, rotation_axes: int, modulation: Modulation = Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), disable_tqdm: bool = False,
                          bandwidth: float = 100, single_antenna: bool = True):
    """
    Generates a graph for accumalated data in Gb vs time passed in days

    :param lat: [deg] Latitude of the ground station
    :param antenna: Antenna class
    :param rotation_axes: 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)
    :param modulation: Modulation class
    :param disable_tqdm: Disable the progress bar
    :param bandwidth: [MHz] Maximum bandwidth that can be used
    :param single_antenna: If the final plot contains 1 antenna or multiple. If True: Shows antenna name in title, if False: Shows antenna name in legend
    :return: plot (not shown!)
    """

    prev_time = 0
    total_accumalated_data = 0
    backlog = 0
    backlog_history = [0]
    elevations = pass_info(lat, data_type="elevation")
    times = pass_info(lat, data_type="time")
    capacity = 0

    data_over_pass_list = np.zeros(91 - 6)
    for max_elev in tqdm(range(7, 91, 1), f"{antenna.name}", disable=disable_tqdm):
        data_over_pass_list[max_elev - 8] = (total_data_over_pass(max_elev, antenna, rotation_axes, modulation=modulation, bandwidth=bandwidth))

    time_data = [0]
    # Summation of data over pass list
    for i, pass_time_list in enumerate(times):  # go through each pass, and take index and the time
        if prev_time != 0:  # max elevation needs to be higher than 7, otherwise the calculation cannot be made
            # print(pass_time_list)
            time_diff = pass_time_list[0] - prev_time[1]
            time_data.append(time_diff.total_seconds() + time_data[-1])
            # print("time_diff=", time_diff)
            # print("time_data=", time_data)

            pass_duration = pass_time_list[1] - pass_time_list[0]
            time_data.append(pass_duration.total_seconds() + time_data[-1])
            # print("pass_duration=", pass_duration)
            # print("time_data=", time_data)

            gathered_data = pl_data_rate * time_diff.total_seconds()
            total_accumalated_data += gathered_data
            backlog += gathered_data
            backlog_history.append(backlog)
            if int(elevations[i]) >= 7:
                sent_data = data_over_pass_list[int(elevations[i]) - 7]
            else:
                sent_data = 0
            capacity += sent_data
            gathered_data = pl_data_rate * pass_duration.total_seconds()
            total_accumalated_data += gathered_data
            backlog += gathered_data
            backlog -= sent_data  # calculates the max data sent over a certain pass

            if backlog < 0:
                backlog = 0
            backlog_history.append(backlog)

        prev_time = pass_time_list

    time_data = [i / 86400 for i in time_data]  # turn x axis into days
    backlog_history = [i / 1000 for i in backlog_history]  # turn y axis into Gbit

    rotation_name = "ERROR IN ROTATION AXIS"
    if rotation_axes == 0:
        rotation_name = "Fixed"
    elif rotation_axes == 1:
        rotation_name = "Single-axis-rotating"
    elif rotation_axes == 2:
        rotation_name = "Gimballed"

    if single_antenna:
        plt.plot(time_data, backlog_history, label=f"{modulation.name}: \nTotal={int(capacity) / 1000}Gbit")
        plt.title(f"{rotation_name}-{antenna.name}: Total acc data = {int(total_accumalated_data / 1000)}[Gb], latitude = {lat}, bandwidth = {bandwidth}MHz")
    else:
        plt.plot(time_data, backlog_history, label=f"{rotation_name}-{antenna.name}-{modulation.name}: \nTotal={int(capacity) / 1000}Gbit")
        plt.title(f"Total acc data = {int(total_accumalated_data / 1000)}[Gb], latitude = {lat}, bandwidth = {bandwidth}MHz")

    plt.ylabel("Data backlog [Gb]")
    plt.xlabel("Time passed [days]")
    return capacity


def plot_accumulated_data_all_mods(lat: int, antenna: Antenna, rotation_axes: int, bandwidth: float = 100):
    """
    Generates a graph for accumalated data in Gb vs time passed in days for all different modulation schemes in Modulations.csv
    :param lat: [deg] Latitude of the ground station
    :param antenna: Antenna class
    :param rotation_axes: 0 for fixed, 1 for single axis rotation or 2 for dual axis rotation (gimballed)
    :param bandwidth: bandwidth: [MHz] Maximum bandwidth that can be used
    :return: shows plot
    """
    modulation_list = read_modulation_file()
    # Call calculation for total modulation and Plot the data
    highest_capacity = 0
    best_modulation = Modulation("ERROR", 0, 0, 0)

    for modulation in tqdm(modulation_list, f'{antenna.name}, lat = {lat}deg'):
        capacity = plot_accumulated_data(lat, antenna, rotation_axes, modulation, disable_tqdm=True, single_antenna=True, bandwidth=bandwidth)
        if capacity > highest_capacity:
            highest_capacity = capacity
            best_modulation = modulation
    print(f'Best modulation = {best_modulation.name}, Capacity = {int(highest_capacity / 1000)}Gb')
    time.sleep(0.1)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # NOTE: Plot the gain graphs of all antennas
    # plot_gains()

    # NOTE: Plot data vs elevation graph
    # plot_data_vs_elev(max_elev=70, antenna=isoflux, rotation_axes=0, modulation=Modulation("8PSK - 3/5", 0.6, 5.5, 1.74), plot_label=f'Fixed-{isoflux.name}', bandwidth=100)
    # plt.legend()
    # plt.show()

    # NOTE: Plot data vs elevation graphs for different passes
    # for elevation in range(10, 91, 10):
    #     plot_data_vs_elev(elevation, isoflux, 0, f'{elevation}°', modulation=Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), bandwidth=100)
    # plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
    # plt.tight_layout(rect=[0, 0, 1, 1])
    # plt.show()

    # NOTE: Plot data vs elevation graphs for all antennas
    # plot_data_vs_elev_all_antennas(max_elev=60, modulation=Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), bandwidth=100)

    # NOTE: Plot data vs elevation graphs for all modulation schemes
    # plot_data_vs_elev_all_modulations(max_elev=90, rotation_axes=0, antenna=isoflux, bandwidth=100)

    # NOTE: plot data vs time graphs for all antennas
    plot_data_vs_time_all_antennas(max_elev=60, modulation=Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), bandwidth=100)

    # NOTE: Plot data vs time graphs for different passes
    plot_data_vs_time_all_elevations(antenna=isoflux, rotation_axes=0, modulation=Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), bandwidth=100)

    # NOTE: Calculate the total data that can be sent for all antennas in the antenna list for a specific latitude and modulation scheme
    # for antenna_class in antenna_list:
    #     total_data(lat=60, antenna=antenna_class, rotation_axes=0, modulation=Modulation("QPSK 3/4", 0.75, 4.03, 1.45), bandwidth=np.inf)

    # NOTE: Plot the accumulated data for all antennas with their best modulation scheme for 78deg latitude and 50MHz bandwidth
    # plot_accumulated_data(78, isoflux,            0, bandwidth=50, single_antenna=False, modulation=Modulation("QPSK 3/4", 0.75, 4.03, 1.45))
    # plot_accumulated_data(78, low_gain_patch,     0, bandwidth=50, single_antenna=False, modulation=Modulation("QPSK 11/20", 0.55, 1.45, 1.06))
    # plot_accumulated_data(78, high_gain_patch,    2, bandwidth=50, single_antenna=False, modulation=Modulation("QPSK 4/5", 0.80, 4.68, 1.55))
    # plot_accumulated_data(78, parabolic_5,        2, bandwidth=50, single_antenna=False, modulation=Modulation("QPSK 4/5", 0.80, 4.68, 1.55))
    # plot_accumulated_data(78, parabolic_10,       2, bandwidth=50, single_antenna=False, modulation=Modulation("16APSK 2/3-L", 0.67, 8.53, 2.57))
    # plt.legend()
    # plt.show()

    # NOTE: plot the accumulated data for all antennas with their best modulation scheme for 60deg latitude and 50MHz bandwidth
    # plot_accumulated_data(60, isoflux,          0, bandwidth=50, single_antenna=False, modulation=Modulation("QPSK 2/3", 0.67, 3.10, 1.29))
    # plot_accumulated_data(60, low_gain_patch,   0, bandwidth=50, single_antenna=False, modulation=Modulation("QPSK 11/20", 0.55, 1.45, 1.06))
    # plot_accumulated_data(60, high_gain_patch,  0, bandwidth=50, single_antenna=False, modulation=Modulation("QPSK 11/20", 0.55, 1.45, 1.06))
    # plot_accumulated_data(60, high_gain_patch,  2, bandwidth=50, single_antenna=False, modulation=Modulation("QPSK 4/5", 0.80, 4.68, 1.55))
    # plot_accumulated_data(60, parabolic_5,      2, bandwidth=50, single_antenna=False, modulation=Modulation("QPSK 4/5", 0.80, 4.68, 1.55))
    # plot_accumulated_data(60, parabolic_10,     2, bandwidth=50, single_antenna=False, modulation=Modulation("16APSK 2/3-L", 0.67, 8.53, 2.57))
    # plt.legend()
    # plt.show()

    # NOTE: plot all modulation schemes for a specific antenna and do this for all latitudes in the list
    # for latitude in [78]:
    #     plot_accumulated_data_all_mods(latitude, phase_array, 0, bandwidth=500)

    pass
