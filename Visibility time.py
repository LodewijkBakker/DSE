from tqdm import tqdm
import numpy as np
import os
import pickle

def vis_time_calc(max_elev, min_elev, orbit_alt=300, inclination=96.7):  # https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=805436
    earth_radius = 6378  # [km]
    orbit_radius = earth_radius + orbit_alt
    orbit_time = 2 * np.pi * np.sqrt(((orbit_radius * 1000) ** 3) / (6.67430e-11 * 5.9736e+24))

    a = 2 / (np.deg2rad(360 / orbit_time) - np.deg2rad(360 / 86400) * np.cos(np.deg2rad(inclination)))
    b = np.arccos(np.cos(np.arccos((earth_radius / orbit_radius) * np.cos(np.deg2rad(min_elev))) - np.deg2rad(min_elev)) /
                  np.cos(np.arccos((earth_radius / orbit_radius) * np.cos(np.deg2rad(max_elev))) - np.deg2rad(max_elev)))
    return a * b


delta = 0.1

full_times = []
for max_elev in tqdm(range(7, 91, 1)):
    times = len(np.arange(-max_elev + 5, max_elev - 5 + delta, delta)) * [0]
    prev_vis_time = -1

    for i in tqdm(np.arange(0, int((max_elev - 5) * (1 / delta)), 1), f"{max_elev} Visibility", disable=True):  # visibility time for each elev
        vis_time = vis_time_calc(max_elev, 5 + (i * delta))
        if i == 0:
            times[-1] = vis_time
            times[int(len(times) / 2)] = vis_time / 2
        else:
            times[i] = times[i - 1] + ((prev_vis_time - vis_time) / 2)
            times[-i - 1] = times[-i] - ((prev_vis_time - vis_time) / 2)

        prev_vis_time = vis_time

    full_times.append(times)

folder_loc = "C:\\Users\\TimLe\\PycharmProjects\\DSE"

if not os.path.exists(f'{folder_loc}\\Visibility_times\\'):
    os.mkdir(os.path.join(f'{folder_loc}\\Visibility_times\\'))

# for i in full_times:
#     print(i[-1])
with open(f'Visibility_times\\Visibility_times_{delta}.txt', 'wb') as f:
    f.truncate(0)
    pickle.dump(full_times, f)
