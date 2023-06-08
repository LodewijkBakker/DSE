import numpy as np
import pandas as pd
import toml
from scipy.integrate import solve_ivp
from astropy.constants import G, M_earth, R_earth
import vg
from numba.experimental import jitclass
from numba import jit

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
        self.A_solar = self.dick["Satellite parameters"]["A_SP_side"]+self.dick["Satellite parameters"]["A_SP_back"]
        self.output = pd.DataFrame(index=range(1), columns=range(1))
        self.mu = G.value * M_earth.value
        self.mu_time=self.dick["Natural constants"]["mu"]
        self.shadow = False
        self.counter = 0
        self.u_z=np.array([0,0,1])
        self.a_earth=self.dick["Natural constants"]["a_earth"]
        self.b_earth=self.dick["Natural constants"]["b_earth"]
        self.Min_drag_config=self.dick["Satellite parameters"]["Min drag"]
        self.third_panel=self.dick["Satellite parameters"]["Third moving panel"]
        self.counter=0
        self.A_SP_side=self.dick["Satellite parameters"]["A_SP_side"]
        self.A_SP_back=self.dick["Satellite parameters"]["A_SP_back"]
        self.mean_sol_sens=np.zeros(2)
        self.cov=np.array([[1,0],[0,1]])
        self.mean_rot=np.zeros(2)

    def angle_between(self, u1, u2):
        #return np.arctan2(np.linalg.norm((np.cross(u1,u2))), np.dot(u1,u2))
        return vg.angle(u1,u2)*np.pi/180
        #return np.arccos(np.clip(np.dot(u1, u2), -1.0, 1.0))

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

    def orbit_time(self):
        return int(2*np.pi*np.sqrt(self.dick["Orbital parameters"]["a"]**3/self.mu))

    def orbital_el_set_up(self, sim_length):
        # t,x,y,z,vx,vy,vz,Eclipse,Angle to Sun,Thrust
        self.Output = pd.DataFrame(index=range(sim_length), columns=["t","x","y","z","v_x","v_y","v_z","m","SMA","Eclipse","Sun multiplier","SP angle","Thrust","f_aero","third panel","eff_f_side","eff_r_side","eff_r_back"])
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
        #print(self.Y0)
        #print("Intial set up: Done\n")

    # @staticmethod
    # @jit
    def dy_vector(self, t, Y):
        dy = np.empty(7)  # x,y,z, xdot,ydot,zdot ,m
        # Orbital parameters
        R = Y[:3]
        V = Y[3:6]
        m = Y[6]
        R_abs = np.linalg.norm(Y[:3])
        V_abs = np.linalg.norm(Y[3:6])
        u_R=R/R_abs
        u_V=V/V_abs
        #a, e, i, omega_AP, omega_LAN, T, EA = self.cart_2_kep(R, V, t)
        a= self.mu/(2/R_abs*self.mu-V_abs**2)

        # Consequences of position and orbital parameters
        z_ang=np.pi/2- self.angle_between(self.u_z,R/R_abs)
        self.height=R_abs-np.sqrt(((self.a_earth**2*np.cos(z_ang))**2+(self.b_earth**2*np.sin(z_ang))**2)/((self.a_earth*np.cos(z_ang))**2+(self.b_earth*np.sin(z_ang))**2))
        #print(a)
        rho = 2.30E-11*np.exp((300e3-self.height)/45.8e3)  # 4.39e-11 maximum

        if (a <= self.R_min):
            self.thrust = True
        elif (a >= self.R_max):
            self.thrust = False

        # Calculates cos of angle between sun vector an velocity
        sh_ang = np.pi / 2 - self.angle_between(self.u_sun, u_R)

        if (np.cos(sh_ang) * R_abs < self.dick["Natural constants"]["Rad Earth"] and np.dot(self.u_sun,
                                                                                            u_R) > 0):  # alternatively sin(sh_ang)
            self.shadow = True
            # Used to turn the satellie by 30 degree 5% of the time
            # if (int(t) % 100 <= 5):
            #     sun_multiplier= np.sin(np.pi/6)
            # else:
            sun_multiplier = abs(np.sin(self.dick["Satellite parameters"]["SP pointing accuracy"]))
            A = self.dick["Satellite parameters"]["A_extended"] + self.A_solar * sun_multiplier + \
                self.dick["Satellite parameters"]["SP Frontal Area"]

            self.check_third = False
        else:
            self.shadow = False
            if (self.Min_drag_config):
                # Used to turn the satellie by 30 degree 5% of the time
                # if (int(t) % 100 <= 5):
                #     sun_multiplier = np.sin(np.pi / 6)
                # else:
                sun_multiplier = np.sin(self.dick["Satellite parameters"]["SP pointing accuracy"])
                A = self.dick["Satellite parameters"]["A_extended"] + self.A_solar * sun_multiplier + \
                    self.dick["Satellite parameters"]["SP Frontal Area"]
                # if (self.third_panel and R[1] >= 0):
                #     self.check_third = True
                #     A+= abs(np.clip(np.dot(self.u_sun, V / V_abs), -1.0, 1.0))*self.dick["Satellite parameters"]["Third panel area"]
            else:
                sun_multiplier = abs(np.clip(np.dot(self.u_sun, V / V_abs), -1.0, 1.0))
                # print(sun_multiplier)
                A = self.dick["Satellite parameters"][
                        "A_extended"] + self.A_SP_side * sun_multiplier + self.A_SP_back * abs(
                    np.sin(self.dick["Satellite parameters"]["SP pointing accuracy"]))
                # print(A)
            self.check_third = False

        #Calculate efficiencies of solar panels
        cos_s_r=abs(np.clip(np.dot(self.u_sun, u_R), -1.0, 1.0))
        cos_s_v=abs(np.clip(np.dot(self.u_sun, u_V), -1.0, 1.0))
        cos_s_side=abs(np.clip(np.dot(self.u_sun, np.cross(u_V, u_R)),-1.0,1.0))
        rot_correction=np.random.normal(0, 1/3)*np.pi/180
        sol_correction=np.random.multivariate_normal(self.mean_sol_sens, self.cov)*5/3*np.pi/180
        #rot_correction=np.random.multivariate_normal(self.mean_rot, self.cov)*1/3*np.pi/180
        #Determine fixed side effeicieny

        if (self.shadow):
            eff_f_side = 0
            eff_r_side =0
            eff_r_back=0
        else:
            if (np.dot(self.u_sun, u_R) < 0):
                eff_f_side = cos_s_r
            else:
                eff_f_side = 0

            # eff_side = np.sqrt(
            #     1 - abs(np.clip(np.dot(self.u_sun, np.cross(V / V_abs, R / R_abs)), -1.0, 1.0)) ** 2)
            ang_side_inplane = np.pi / 2 - np.arccos(cos_s_side)+ sol_correction[0]
            ang_side_outofplane = sol_correction[1] + rot_correction
            eff_r_side = np.cos(ang_side_inplane) * np.cos(ang_side_outofplane)

            #eff_back = np.sqrt(1 - abs(np.clip(np.dot(self.u_sun, V / V_abs), -1.0, 1.0)) ** 2)
            ang_back_inplane=np.pi/2 - np.arccos(cos_s_v)+sol_correction[0]
            ang_back_outofplane=sol_correction[1]+rot_correction
            eff_r_back = np.cos(ang_back_inplane)*np.cos(ang_back_outofplane)

            #print(A)

        # Set up orbital forces
        Grav = -R * self.dick["Natural constants"]["mu"] / R_abs ** 3  # Acceleration due to gravity
        f_aero = 0.5 * rho * V_abs ** 2 * A * self.dick["Satellite parameters"]["Cd"]

        # Set up dy
        dy[:3] = V

        if (self.thrust):
            dy[3:6] = (self.T - f_aero) * V / V_abs / m + Grav
            dy[6] = -self.m_dot
        else:
            dy[3:6] = -f_aero * V / V_abs / m + Grav
            dy[6] = 0
        # print(f"t:{t} f_aero: {f_aero} Semi-major axis: {a}, Thrust: {self.thrust} Radius: {R_abs} Velocity: {V_abs}")
        self.Output.loc[round(t-self.t0)]=[round(t),R[0],R[1],R[2],V[0],V[1],V[2],m,a,int(self.shadow == True),sun_multiplier,np.pi/2-sh_ang,int(self.thrust == True),f_aero,int(self.check_third == True),eff_f_side,eff_r_side,eff_r_back]
        return dy

    def determine_char(self, Y):
        self.R = Y[:3]
        self.V = Y[3:6]
        self.m = Y[6]
        self.R_abs = np.linalg.norm(Y[:3])
        self.V_abs = np.linalg.norm(Y[3:6])
    # def orientation(self):
    #   self.theta==

    def final_properties(self, filename,sim_length, sol,row):
        print(f"writing to {filename}")
        file = open(filename, "w")
        file.write("---Setup---\n")
        file.write(f"Propulsion system: {row['Name']} \n")
        file.write(f"SP Area:{self.A_solar}m^2\n")
        file.write(f"Specific impulse: {row['Specific Impulse']}s\n")
        file.write(f"Thrust: {row['Thrust']}mN\n")
        file.write("---Results---\n")
        self.Output["Delta v"] = self.Output["f_aero"] / self.Output["m"]
        self.delta_v=5 * 365.25 * 24 * 60 * 60 / sim_length * self.Output['Delta v'].sum() * 1.2
        file.write(
            f"Total DeltaV requirement:{round(self.delta_v,2)}\n")
        fuel= (np.exp(self.delta_v/(self.I_sp*self.dick["Natural constants"]['g0']))-1)*self.dick['Satellite parameters']['dry mass']
        file.write(f"Fuel mass: {round(fuel,2)}\n")
        file.write(f"Maximum aerodynamic force: {round(self.Output['f_aero'].max()*1000,3)}mN\n")
        file.write(f"Average aerodynamic force: {round(self.Output['f_aero'].mean() * 1000, 3)}mN\n")
        file.write("---Debugging---\n")
        self.determine_char(sol.y[:, -1])
        file.write(
            f"Final result \nR: {self.R_abs / 1e3 - self.dick['Natural constants']['Rad Earth'] / 1e3}km V: {self.V_abs / 1e3}km\n")
        file.write(f"Eclipse fraction in orbit: {self.Output['Eclipse'].sum()/self.Output['Eclipse'].size}\n")
        file.write(f"Third panel full efficiency fraction: {self.Output['third panel'].sum()/self.Output['Eclipse'].size}\n")
        file.close()

def average_final_properties(filename,sat_bol,sat_eol,sim_length,sol_bol, sol_eol, row):
    print(f"writing to {filename}")
    file = open(filename, "w")

    file.write("---Setup---\n")
    file.write(f"Propulsion system: {row['Name']} \n")
    file.write(f"Peak power: {row['Power']}W\n")
    file.write(f"Solar panel setup {sat_bol.dick['Name']} \n")
    file.write(f"SP Area:{sat_bol.A_solar}m^2\n")
    file.write(f"Specific impulse: {row['Specific Impulse']}s\n")
    file.write(f"Thrust: {row['Thrust']}mN\n")
    file.write("---Results---\n")
    sat_bol.Output["Delta v"] = sat_bol.Output["f_aero"] / sat_bol.Output["m"]
    sat_eol.Output["Delta v"] = sat_eol.Output["f_aero"] / sat_eol.Output["m"]
    delta_v = 5 * 365.25 * 24 * 60 * 60 / sim_length * (sat_bol.Output['Delta v'].sum()+sat_eol.Output['Delta v'].sum())/2 * 1.2
    file.write(
        f"Design Total DeltaV requirement:{round(delta_v, 2)}s\n")
    file.write(f"Worst-worst case DeltaV requirement: {round(delta_v*1.908696)}s\n")

    fuel = (np.exp(delta_v / (sat_bol.I_sp * sat_bol.dick["Natural constants"]['g0'])) - 1) * \
           sat_bol.dick['Satellite parameters']['dry mass']
    file.write(f"Fuel mass: {round(fuel, 2)}kg\n")
    file.write(f"Fuel volume: {round(fuel/row['Density'], 2)}U\n")
    file.write(f"Very prelim. prop. system mass:{round(0.01+1.3*fuel+row['Mass']/1000, 2)}kg\n")
    file.write(f"Very prelim. prop. system volume: {round((0.01+0.3*fuel)/2.7+row['Volume']+fuel/row['Density'], 2)}U\n")
    file.write(f"Maximum aerodynamic force: {round(max(sat_bol.Output['f_aero'].max(),sat_eol.Output['f_aero'].max()) * 1000, 3)}mN\n")
    f_m=(sat_bol.Output['f_aero'].mean()+sat_eol.Output['f_aero'].mean())/2
    file.write(f"Average aerodynamic force: {round(f_m * 1000, 3)}mN\n")
    burn_time=round(f_m*sat_bol.orbit_time()*3/(sat_bol.T))
    file.write(f"Time thrusting every 3 orbits average overall: {burn_time}s\n")
    file.write(
        f"Time thrusting every 3 orbits maximum: {round(f_m*1.908696 * sat_bol.orbit_time() * 3 / (sat_bol.T))}s\n")
    file.write(f"Duty cycle: {f_m/(sat_bol.T)}")
    file.write(f"Duty cycle worst: {f_m*1.908696/ (sat_bol.T)}")
    file.write(f"Energy spent per 3 orbits(Avg burn*1.5):{burn_time*1.5*row['Power']/3600}Wh")
    # s = sat_bol.Output["Thrust"]
    # file.write(f"Longest continious thrust: {(~s).cumsum()[s].value_counts().max()}s\n")
    file.write("---Debugging---\n")
    sat_bol.determine_char(sol_bol.y[:, -1])
    file.write(
        f"Final result \nR: {sat_bol.R_abs / 1e3 - sat_bol.dick['Natural constants']['Rad Earth'] / 1e3}km V: {sat_bol.V_abs / 1e3}km\n")
    file.write(f"Eclipse fraction in orbit: {sat_bol.Output['Eclipse'].sum() / sat_bol.Output['Eclipse'].size}\n")
    file.write(
        f"Third panel full efficiency fraction: {sat_bol.Output['third panel'].sum() / sat_bol.Output['Eclipse'].size}\n")

    file.close()



if __name__ == "__main__":
    LAMP = Satellite("Inputs/Config-4.toml")
    # LAMP.circular_set_up(300e3,LAMP.dry_m+4)
    # print(LAMP.Y0)
    Orbital_period=LAMP.orbit_time()
    sim_length=int(1*Orbital_period)
    LAMP.orbital_el_set_up(sim_length)
    sol=solve_ivp(LAMP.dy_vector,LAMP.time_range,LAMP.Y0, method="Radau", max_step=1 ,t_eval = LAMP.time_range, rtol = 1)
    LAMP.determine_char(sol.y[:,-1])
    print(f"Final result \nR: {LAMP.R_abs/1e3-LAMP.dick['Natural constants']['Rad Earth']/1e3}km V: {LAMP.V_abs/1e3}km")
    LAMP.Output.to_csv("Outputs/Config-4.csv")
    LAMP.Output["Delta v"]=LAMP.Output["f_aero"]/LAMP.Output["m"]
    print(5*365.25*24*60*60/sim_length*LAMP.Output["Delta v"].sum()*1.2)
    LAMP.Output["Ones"] = 1
    ecl_frac=(LAMP.Output["Eclipse"]).sum()/LAMP.Output["Eclipse"].size
    print(ecl_frac)
    print(LAMP.Output["third panel"].sum()/LAMP.Output["Eclipse"].size)
    print(LAMP.Output["f_aero"].max())
    print("Height:",LAMP.height)
    print(LAMP.Output["Thrust"].sum()/LAMP.Output["Thrust"].size)
    #print(side_eff)
    if(LAMP.Min_drag_config):
        print(f"Eff side panel = {LAMP.Output.eff_side.sum()/(sim_length/2)}")
    else:
        print(f"Eff side panel = {LAMP.Output.eff_side.sum() / (sim_length* (1-ecl_frac))}")
    print(f"Eff back panel = {LAMP.Output.eff_back.sum()/(sim_length*(1-ecl_frac))}")

    #Eclipse/Penumbra/E
    # t,x,y,z,vx,vy,vz,Eclipse,Angle to Sun,Thrust,f_aero

    #M_struc=10+0.3Mp
    #Vol_struc=0.1*V_fuel
