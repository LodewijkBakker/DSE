import numpy as np
import pandas as pd
import toml
from scipy.integrate import solve_ivp
from astropy.constants import G, M_earth, R_earth


class Satellite:
    def __init__(self, filename):
        # Dicktionary of all inputs
        self.dick = toml.load(filename)
        self.Pos_Sun = np.array(self.dick["Natural constants"]["Pos sun"])
        self.u_sun = np.array(self.dick["Natural constants"]["u_sun"])
        self.I_sp = self.dick["Satellite parameters"]["I_sp"]
        self.T = self.dick["Satellite parameters"]["Thrust"]
        self.R_min = self.dick["Natural constants"]["Rad Earth"] + self.dick["Orbital parameters"]["Min alt"]
        self.R_max = self.dick["Natural constants"]["Rad Earth"] + self.dick["Orbital parameters"]['Max alt']
        self.thrust = False
        self.dry_m = self.dick["Satellite parameters"]["dry mass"]
        self.m_dot = self.T / self.dick["Natural constants"]["g0"] / self.dick["Satellite parameters"]["I_sp"]
        self.A_solar = self.dick["Satellite parameters"]["A_solar"]
        self.output = pd.DataFrame(index=range(1), columns=range(1))
        self.mu = G.value * M_earth.value
        self.shadow = False
        self.counter = 0

    def angle_between(self, u1, u2):
        return np.arccos(np.clip(np.dot(u1, u2), -1.0, 1.0))

    def cart_2_kep(self, r_vec, v_vec, t):
        # 1
        h_bar = np.cross(r_vec, v_vec)
        h = np.linalg.norm(h_bar)
        # 2
        r = np.linalg.norm(r_vec)
        v = np.linalg.norm(v_vec)
        # 3
        E = 0.5 * (v ** 2) - self.mu / r
        # 4
        a = -self.mu / (2 * E)
        # 5
        e = np.sqrt(1 - (h ** 2) / (a * self.mu))
        # 6
        i = np.arccos(h_bar[2] / h)
        # 7
        omega_LAN = np.arctan2(h_bar[0], -h_bar[1])
        # 8
        # beware of division by zero here
        lat = np.arctan2(np.divide(r_vec[2], (np.sin(i))), \
                         (r_vec[0] * np.cos(omega_LAN) + r_vec[1] * np.sin(omega_LAN)))
        # 9
        p = a * (1 - e ** 2)
        nu = np.arctan2(np.sqrt(p / self.mu) * np.dot(r_vec, v_vec), p - r)
        # 10
        omega_AP = lat - nu
        # 11
        EA = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2))
        # 12
        n = np.sqrt(self.mu / (a ** 3))
        T = t - (1 / n) * (EA - e * np.sin(EA))
        return a, e, i, omega_AP, omega_LAN, T, EA

    def kep_2_cart(self, a, e, i, omega_AP, omega_LAN, T, EA, t):

        # 1
        n = np.sqrt(self.mu / (a ** 3))
        M = n * (t - T)
        # 2
        MA = EA - e * np.sin(EA)
        # 3
        #
        # ERROR WAS HERE
        # nu = 2*np.arctan(np.sqrt((1-e)/(1+e)) * np.tan(EA/2))
        nu = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(EA / 2))
        # 4
        r = a * (1 - e * np.cos(EA))
        # 5
        h = np.sqrt(self.mu * a * (1 - e ** 2))
        # 6
        Om = omega_LAN
        w = omega_AP

        X = r * (np.cos(Om) * np.cos(w + nu) - np.sin(Om) * np.sin(w + nu) * np.cos(i))
        Y = r * (np.sin(Om) * np.cos(w + nu) + np.cos(Om) * np.sin(w + nu) * np.cos(i))
        Z = r * (np.sin(i) * np.sin(w + nu))

        # 7
        p = a * (1 - e ** 2)

        V_X = (X * h * e / (r * p)) * np.sin(nu) - (h / r) * (np.cos(Om) * np.sin(w + nu) + \
                                                              np.sin(Om) * np.cos(w + nu) * np.cos(i))
        V_Y = (Y * h * e / (r * p)) * np.sin(nu) - (h / r) * (np.sin(Om) * np.sin(w + nu) - \
                                                              np.cos(Om) * np.cos(w + nu) * np.cos(i))
        V_Z = (Z * h * e / (r * p)) * np.sin(nu) + (h / r) * (np.cos(w + nu) * np.sin(i))

        return [X, Y, Z], [V_X, V_Y, V_Z]

    def circular_set_up(self, h, m_wet, time):
        self.Y0 = np.empty(7)
        R = self.dick["Natural constants"]["Rad Earth"] + h
        V = np.sqrt(self.dick["Natural constants"]["mu"] / R)

        u = np.array([1, 0, 0])
        u_bar = u / np.linalg.norm(u)
        self.Y0[:3] = u_bar * R

        v = np.array([0, 1, 1])
        v_bar = v / np.linalg.norm(v)
        self.Y0[3:6] = v_bar * V

        self.Y0[6] = m_wet

    def orbital_el_set_up(self, sim_length):
        # t,x,y,z,vx,vy,vz,Eclipse,Angle to Sun,Thrust
        self.Output = pd.DataFrame(index=range(sim_length), columns=["t","x","y","z","v_x","v_y","v_z","m","Eclipse","Angle to Sun","Thrust","f_aero"])
        self.counter=0
        self.time_range = [self.dick["Orbital parameters"]["t"], self.dick["Orbital parameters"]["t"]+sim_length]
        self.time_eval=np.arange(self.dick["Orbital parameters"]["t"], self.dick["Orbital parameters"]["t"]+sim_length)
        self.t0=self.dick["Orbital parameters"]["t"]
        a = self.dick["Orbital parameters"]["a"]
        R, V = self.kep_2_cart(a, self.dick["Orbital parameters"]["e"], self.dick["Orbital parameters"]["i"],
                               self.dick["Orbital parameters"]["omega_AP"],
                               self.dick["Orbital parameters"]["omega_LAN"], self.dick["Orbital parameters"]["T"],
                               self.dick["Orbital parameters"]["EA"], self.dick["Orbital parameters"]["t"])
        self.Y0=np.concatenate((R,V,np.array([self.dick["Satellite parameters"]["dry mass"]+self.dick["Satellite parameters"]["fuel mass"]])))
        print(self.Y0)
        print("Intial set up: Done\n")

    def dy_vector(self, t, Y):
        dy = np.empty(7)  # x,y,z, xdot,ydot,zdot ,m
        # Orbital parameters
        R = Y[:3]
        V = Y[3:6]
        m = Y[6]
        R_abs = np.linalg.norm(Y[:3])
        V_abs = np.linalg.norm(Y[3:6])
        # a, e, i, omega_AP, omega_LAN, T, EA = self.cart_2_kep(R, V, t)
        a= self.mu/(2/R_abs*self.mu-V_abs**2)
        # Consequences of position and orbital parameters
        rho = 2.30E-11  # Add function for the density
        if (a <= self.R_min):
            self.thrust = True
        elif (a >= self.R_max):
            self.thrust = False

        # Calculates cos of angle between sun vector an velocity
        sun_multiplier = abs(np.clip(np.dot(self.u_sun, V / V_abs), -1.0, 1.0))
        sh_ang = np.pi / 2 - self.angle_between(self.u_sun, R / R_abs)
        if (np.cos(sh_ang) * R_abs < self.dick["Natural constants"]["Rad Earth"] and sh_ang < np.pi / 2):  # alternatively sin(sh_ang)
            self.shadow = True
        else:
            self.shadow=False

        # Set up orbital forces
        Grav = -R * self.dick["Natural constants"]["mu"] / R_abs ** 3  # Acceleration due to gravity
        A = self.dick["Satellite parameters"]["A_extended"] + self.A_solar * sun_multiplier
        f_aero = 0.5 * rho * V_abs ** 2 * A * self.dick["Satellite parameters"]["Cd"]

        # Set up dy
        dy[:3] = V
        if (self.thrust):
            dy[3:6] = (self.T - f_aero) * V / V_abs / m + Grav
            dy[6] = -self.m_dot
        else:
            dy[3:6] = -f_aero * V / V_abs / m + Grav
            dy[6] = 0

        #print(f"t:{t} f_aero: {f_aero} Semi-major axis: {a}, Thrust: {self.thrust} Radius: {R_abs} Velocity: {V_abs}")

        self.Output.loc[round(t-self.t0)]=[round(t),R[0],R[1],R[2],V[0],V[1],V[2],m,int(self.shadow == True),sun_multiplier,int(self.shadow == False),f_aero]
        return dy

    def determine_char(self, Y):
        self.R = Y[:3]
        self.V = Y[3:6]
        self.m = Y[6]
        self.R_abs = np.linalg.norm(Y[:3])
        self.V_abs = np.linalg.norm(Y[3:6])
    # def orientation(self):
    #   self.theta==


if __name__ == "__main__":
    LAMP = Satellite("Inputs/Config.toml")
    # LAMP.circular_set_up(300e3,LAMP.dry_m+4)
    # print(LAMP.Y0)
    length=5431
    LAMP.orbital_el_set_up(length)
    sol=solve_ivp(LAMP.dy_vector,LAMP.time_range,LAMP.Y0, method="Radau", max_step=1 ,t_eval = LAMP.time_range, rtol = 1)
    LAMP.determine_char(sol.y[:,-1])
    print(f"Final result\n R: {LAMP.R_abs/1e3-LAMP.dick['Natural constants']['Rad Earth']/1e3}km V: {LAMP.V_abs/1e3}km")
    LAMP.Output.to_csv("Outputs/Test1")
    LAMP.Output["Delta v"]=LAMP.Output["f_aero"]/LAMP.Output["m"]

    print(5*365.25*24*60*60/length*LAMP.Output["Delta v"].sum())
    #Eclipse/Penumbra/E
    # t,x,y,z,vx,vy,vz,Eclipse,Angle to Sun,Thrust,f_aero
