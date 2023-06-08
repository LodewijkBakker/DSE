import numpy as np


class Env:
    def __init__(self):
        self.sigma_boltzmann = 5.670374419e-8
        self.orbit_alt = (6371+300)*1e3
        self.mu_earth = 3.986004418e14
        self.v_orbit = ((self.mu_earth/self.orbit_alt)**0.5)
        # self.t_orbit = int(2*np.pi*self.orbit_alt/self.v_orbit)
        self.t_orbit = 5429
        self.f_orbit = 0.404
        self.albedo = 0.3