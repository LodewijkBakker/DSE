import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from tqdm import tqdm
import re
from datetime import datetime
from scipy import integrate
import pickle


with open(f'Visibility_times\\Visibility_times_0.1.txt', 'rb') as f:
    visibility_times = pickle.load(f)

class Antenna:
    def __init__(self, name, gains, power, gimballed: bool, rotating: bool):
        antenna_angles = [-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90]

        self.name = name
        self.gains = gains
        self.gains += reversed(gains[0:-1])
        self.gains_inter = interpolate.interp1d(antenna_angles, gains, kind='cubic')
        self.power = power
        self.gimballed = gimballed
        self.rotating = rotating


class Modulation:
    def __init__(self, name, code_rate, snr, bits_hz):
        self.name = name
        self.code_rate = code_rate
        self.snr = snr
        self.bits_hz = bits_hz


omnidirectional = Antenna("Omnidirectional", [0, 0, 0, 0, 0, 0, 0], power=2, gimballed=False, rotating=False)
low_gain_patch = Antenna("Low gain patch", [-4.3, -2, 0.6, 3, 4.9, 6.1, 6.7], power=4, gimballed=True, rotating=True)
high_gain_patch = Antenna("High gain patch", [-15.7, -10.8, -6, -0.5, 6.3, 10.5, 11.5], power=3, gimballed=True, rotating=True)
comb_patch = Antenna("Combined patch", [-4.3, -2, 0.6, 3, 6.3, 10.5, 11.5], power=3, gimballed=True, rotating=True)
reflect_array = Antenna("Reflectarray", [-10, -10, -10, -10, -8, 4, 24], power=5, gimballed=False, rotating=True)
isoflux = Antenna("Isoflux", [-8, -3, 2, 5, 5.8, 3.7, 2.9], power=5, gimballed=False, rotating=True)
phase_array = Antenna("Phased array", [-20, -20, 30, 30, 30, 30, 30], power=2, gimballed=False, rotating=False)

# omnidirectional = Antenna("Omnidirectional",    [0,     0,      0,      0,      0,      0,      0],     power=3,    gimballed=False,  rotating=False)
# low_gain_patch  = Antenna("Low gain patch",     [-4.3,  -2,     0.6,    3,      4.9,    6.1,    6.7],   power=3,    gimballed=True,   rotating=True)
# high_gain_patch = Antenna("High gain patch",    [-15.7, -10.8,  -6,     -0.5,   6.3,    10.5,   11.5],  power=3,    gimballed=True,   rotating=True)
# comb_patch      = Antenna("Combined patch",     [-4.3,  -2,     0.6,    3,      6.3,    10.5,   11.5],  power=3,    gimballed=True,   rotating=True)
# reflect_array   = Antenna("Reflectarray",       [-10,   -10,    -10,    -10,    -8,     4,      24],    power=3,    gimballed=False,  rotating=True)
# isoflux         = Antenna("Isoflux",            [-8,    -3, 	2,  	5,  	5.8,	3.7,	2.9],   power=3,    gimballed=False,  rotating=True)
# phase_array     = Antenna("Phased array",       [-20,   -20,	30,	    30, 	30, 	30, 	30],    power=3,    gimballed=False,  rotating=False)

antenna_list = [omnidirectional, low_gain_patch, high_gain_patch, reflect_array, isoflux]


def slant_dist_calc(elev_angle):
    earth_radius = 6378  # [km]
    orbit_alt = 300  # [km]
    orbit_radius = earth_radius + orbit_alt
    slant_dist = earth_radius * (
            (((orbit_radius ** 2) / (earth_radius ** 2) - (np.cos(elev_angle * np.pi / 180) ** 2)) ** 0.5) - np.sin(elev_angle * np.pi / 180))

    return slant_dist


def EIRP_calc(power, gain):
    power_dbw = 10 * np.log10(power)
    l_line = -1  # [dB]  TODO-CHECK VALUE
    l_filter = 0  # [dB]  TODO-CHECK VALUE
    return power_dbw + l_line + l_filter + gain


def link_budget(elev_angle, power, gain, modulation=Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), print_values=False):
    frequency = 8000  # [MHz] X-Band
    gs_gain = 43.2  # [dB]  Based on 2.4m X-band
    l_pol = -3  # [dB]  TODO-CHECK VALUE
    l_line = -1  # [dB]  Based on delfi ground station
    t_noise = 180  # [K]   Based on EWI GS
    link_margin = 3  # [dB]  As advised
    l_atm_zenith = -0.3  # [dB]  TODO-CHECK VALUE, zenith total atmospheric attenuation
    boltzmann = -228.6  # [dBW/K/Hz] Logarithmic representation

    #   Modulation dependent: # TODO DVB-S2X - ESA ACM (https://datumsystems.com/product/dvbs2x/)
    req_snr = modulation.snr  # [dB]
    code_rate = modulation.code_rate  # [-]
    coding_gain = 0  # [dB]

    max_bandwidth = 3000  # [MHz] TODO GET VALUE
    max_data_rate = max_bandwidth * modulation.bits_hz * 1000000

    slant_dist = slant_dist_calc(elev_angle)
    eirp = EIRP_calc(power, gain)
    wavelength = 299.792458 / frequency
    l_path = -20 * np.log10(4 * np.pi * slant_dist * (1000 / wavelength))
    l_atm = l_atm_zenith / np.sin(elev_angle * np.pi / 180)
    l_total_channel = l_atm + l_path

    signal_level_gs = eirp + l_total_channel
    t_noise_db = 10 * np.log10(t_noise)
    snr = link_margin + req_snr - coding_gain
    snr_power_density = signal_level_gs + gs_gain + l_pol + l_line - t_noise_db - boltzmann
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


def data_vs_elev(max_elev, antenna, rotation_axes, modulation=Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74)):
    power = antenna.power
    gain_inter_func = antenna.gains_inter

    elev_angles = np.hstack((np.arange(5, max_elev + 1), np.flip(np.arange(5, max_elev))))  # same as elev_angles = [*range(5, max_elev + 1, 1)] + [*range(max_elev - 1, 4, -1)]
    earth_radius = 6378  # [km]
    orbit_alt = 300  # [km]
    orbit_radius = earth_radius + orbit_alt

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

        data_rate = link_budget(elev_angles[i], power, gain, modulation) / 1000000
        data_rates[i] = data_rate

    return data_rates


def vis_time_calc(max_elev, min_elev, orbit_alt=300, inclination=96.7):  # https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=805436
    earth_radius = 6378  # [km]
    orbit_radius = earth_radius + orbit_alt
    orbit_time = 2 * np.pi * np.sqrt(((orbit_radius * 1000) ** 3) / (6.67430e-11 * 5.9736e+24))

    a = 2 / (np.deg2rad(360 / orbit_time) - np.deg2rad(360 / 86400) * np.cos(np.deg2rad(inclination)))
    b = np.arccos(np.cos(np.arccos((earth_radius / orbit_radius) * np.cos(np.deg2rad(min_elev))) - np.deg2rad(min_elev)) /
                  np.cos(np.arccos((earth_radius / orbit_radius) * np.cos(np.deg2rad(max_elev))) - np.deg2rad(max_elev)))
    return a * b


def data_vs_time(max_elev, antenna, rotation_axes, modulation=Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), delta=0.1, plot=True, disable_tqdm=True):
    start = time.perf_counter()

    data_elev = data_vs_elev(max_elev, antenna, rotation_axes, modulation)
    elevation = np.arange(-max_elev + 5, max_elev - 4, 1)

    elev_data_inter = interpolate.interp1d(elevation, data_elev, kind='linear', fill_value=(data_elev[0], data_elev[-1]), bounds_error=False)
    elev_new = np.arange(-max_elev + 5, max_elev - 5 + delta, delta)
    plot = False
    if plot:
        plt.plot(elevation, data_elev, 'o')
        plt.plot(elev_new, elev_data_inter(elev_new), '--', label=antenna.name)
        plt.xlabel("Elevation [deg]")
        plt.ylabel("Max Data-rate [Mbit/s]")
        plt.show()

    times = visibility_times[int(max_elev) - 7]

    time_elev_inter = interpolate.interp1d(times, elev_new, kind='linear', fill_value=(times[0], times[-1]), bounds_error=False)

    inter = time.perf_counter() - start
    print("interim time:", inter)
    if plot:
        plt.xlabel("Time passed [s]")
        plt.ylabel("Elevation [deg]")
        plt.plot(times, elev_new, 'o')
        plt.plot(times, time_elev_inter(times), '--')
        plt.show()

    # Calculate data rates vs time
    new_data_rate = np.zeros(len(times))
    for i, elev in enumerate(tqdm(elev_new, f"{antenna.name} Data-rate", disable=disable_tqdm)):
        new_data_rate[i] = elev_data_inter(elev)

    # new_data_rate = np.inter(times, )

    print("final time:", time.perf_counter() - inter - start)

    return times, new_data_rate


def data_over_pass_antenna(max_elev, antenna, rotation_axes, delta,
                           modulation=Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74)):  # same as data_over_pass, but gatheres the times and data_rates lists first

    times, data_rates = data_vs_time(max_elev, antenna, rotation_axes, modulation, delta, plot=False)
    return integrate.simpson(data_rates, times)


def pass_info(lat, data_type):
    if lat not in [0, 30, 60, 78, 90]:
        raise ValueError("Latitude should be 0, 30, 60, 78 or 90 degrees")
    filename = f"90d_lat{lat}"
    with open("Orbital pass data/" + filename + ".txt") as f:
        count = 0
        passes = []
        for line in f.readlines():
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


def total_data(lat, antenna, rotation_axes, delta):
    elevations = pass_info(lat, "elevation")
    total = 0
    for max_elev in tqdm(elevations, f"{antenna.name}"):
        if max_elev > 7:
            total += data_over_pass_antenna(int(max_elev), antenna, rotation_axes, delta)
    print(f"{antenna.name} total data = {total / 1000}Gb")
    print(" ")
    return total


def plot_data_vs_elev(max_elev, antenna, rotation_axes, plot_label=''):
    x = [*range(-max_elev + 5, max_elev - 4, 1)]
    data_rates = data_vs_elev(max_elev, antenna, rotation_axes)
    plt.plot(x, data_rates, label=plot_label)


def plot_gains():
    angles = np.arange(-90, 91, 1)
    antenna_angles = [-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90]
    for antenna in antenna_list:
        plt.plot(angles, antenna.gains_inter(angles), label=antenna.name)
        plt.plot(antenna_angles, antenna.gains, 'o')
    plt.legend()
    plt.show()


def plot_all_elev(max_elev, rotation_axes):
    for antenna in antenna_list:
        plot_data_vs_elev(max_elev, antenna, rotation_axes, antenna.name)
    plt.suptitle(f'fixed, max elevation = {max_elev}')
    plt.xlabel("elevation angle [deg]")
    plt.ylabel("Max data rate [Mbit/s]")
    plt.legend()
    plt.show()


def plot_all_data_vs_time(max_elev, rotation_axes, delta):
    for antenna in antenna_list:
        times, data_rate = data_vs_time(max_elev, antenna, rotation_axes, delta, plot=False, disable_tqdm=False)
        plt.plot(times, data_rate, '--', label=f"{antenna.name},    Total = {integrate.simpson(times, data_rate)}")
    plt.xlabel("Time passed [s]")
    plt.ylabel("Max Data-rate [Mbit/s]")
    plt.legend()
    plt.show()


def plot_single_data_vs_time(antenna, rotation_axes, delta, modulation=Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74)):
    for elevation in [90, 80, 70, 60, 50, 40, 30, 20]:
        times, data_rate = data_vs_time(elevation, antenna, rotation_axes, modulation, delta, plot=False, disable_tqdm=False)
        plt.plot(times, data_rate, '--', label=f"Max elevation = {elevation}deg,    Total = {integrate.simpson(times, data_rate)}")
    plt.xlabel("Time passed [s]")
    plt.ylabel("Max Data-rate [Mbit/s]")
    plt.title(f"{antenna.name}")
    plt.legend()
    plt.show()


def plot_pass_elev_dist(lat):
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


def plot_accumalated_data(lat, antenna, rotation_axes, delta, modulation=Modulation("8PSK - Generic 3/5", 0.6, 5.5, 1.74), disable_tqdm=False):
    payload_data_rate = 0.5175368314 * 0.2 * 8  # [Mb/s], !! TODO check 8
    # payload_data_rate = 20 * 0.2 * 8  # [Mb/s]

    prev_time = 0
    total_accumalated_data = 0
    backlog = 0
    backlog_history = [0]
    elevations = pass_info(lat, data_type="elevation")
    times = pass_info(lat, data_type="time")
    capacity = 0

    data_over_pass_list = np.zeros(91 - 6)  # hardcoded where does 7 and 91 come from
    for max_elev in tqdm(range(7, 91, 1), f"{antenna.name}", disable=disable_tqdm):
        data_over_pass_list[max_elev - 8] = (data_over_pass_antenna(max_elev, antenna, rotation_axes, delta, modulation))

    # data_over_pass_list = []
    #
    # # Calculate data per pass, about 6 seconds per pass
    # for max_elev in tqdm(range(7, 91, 1), f"{antenna.name}", disable=False):
    #     data_over_pass_list.append(data_over_pass_antenna(max_elev, antenna, rotation_axes, delta, modulation))

    time_data = [0]
    # for i in range(len(times)):
    #     if prev_time != 0:
    #         time_diff = times[i][0] - prev_time[1]
    #         time_data.append(time_diff.total_seconds() + time_data[-1])

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

            gathered_data = payload_data_rate * time_diff.total_seconds()
            total_accumalated_data += gathered_data
            backlog += gathered_data
            backlog_history.append(backlog)
            if int(elevations[i]) >= 7:
                sent_data = data_over_pass_list[int(elevations[i]) - 7]
            else:
                sent_data = 0
            capacity += sent_data
            gathered_data = payload_data_rate * pass_duration.total_seconds()
            total_accumalated_data += gathered_data
            backlog += gathered_data
            backlog -= sent_data  # calculates the max data sent over a certain pass

            if backlog < 0:
                backlog = 0
            backlog_history.append(backlog)

        prev_time = pass_time_list

    total_time = times[-1][1] - times[0][0]
    # print(total_time)
    # print(total_time.total_seconds())
    # print(time_data)

    time_data = [i / 86400 for i in time_data]
    plt.plot(time_data, backlog_history, label=f"Capacity {antenna.name}-{modulation.name}: {int(capacity) / 1000}Gbit")
    plt.ylabel("Data backlog [Mb]")
    plt.xlabel("Time passed [days]")
    plt.title(f"Total acc data = {int(total_accumalated_data) / 1000}[Gb]")


def plot_accumalated_data_all_mods(lat, antenna, rotation_axes, delta):
    with open("Modulations.csv") as file:
        lines = file.readlines()

        # Only iteration on modulations
        for i, line in enumerate(tqdm(lines)):
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

                # Call calculation for total modulation and Plot the data
                plot_accumalated_data(lat, antenna, rotation_axes, delta, mod, disable_tqdm=True)


# plot_gains()

# plot_all_data_vs_time(70, 2, 0.1)

# plot_single_data_vs_time(isoflux, 0, 0.1)
# plot_single_data_vs_time(low_gain_patch, 0, 0.1)
# plot_single_data_vs_time(high_gain_patch, 0, 0.1)
# plot_single_data_vs_time(comb_patch, 0, 0.1)
# plot_single_data_vs_time(reflect_array, 0, 0.1)

# plot_all_elev(70, 0)
# plot_all_elev(60, 0)
# plot_all_elev(50, 0)
# plot_all_elev(40, 0)

# for elevation in range(10, 91, 10):
#     plot_data_vs_elev(elevation, high_gain_patch, 2, f'elev = {elevation}')
# plt.legend()
# plt.show()


# plot_data_vs_elev(70, isoflux, 0, plot_label=f'{isoflux.name}')
# plot_data_vs_elev(70, low_gain_patch, 0, plot_label=f'{low_gain_patch.name}')
# plt.legend()
# plt.show()

# for antenna in antenna_list:
#     total_data(60, antenna, 0, 0.1)

# plot_accumalated_data_all_mods(78, isoflux, 0, 0.1)

plot_accumalated_data(0, isoflux, 0, 0.1, disable_tqdm=True)
# plot_accumalated_data(0, low_gain_patch, 0, 0.1)
# plot_accumalated_data(60, high_gain_patch, 2, 0.1)
# plot_accumalated_data(0, reflect_array, 0, 0.1)
plt.legend()
plt.show()

# print(vis_time_calc(66, 10, 415, 51.6))  # ISS pass, should be 6 min 36 sec, or 396 secs
# print(data_vs_time(70, isoflux, 0))