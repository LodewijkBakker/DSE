import numpy as np
from scipy.integrate import solve_ivp
import toml

class Satellite:
    def __init__(self,filename,A_solar, I_sp, dry_m, Thrust, h_min, h_max):
        self.dick=toml.load(filename)
        self.Pos_Sun=np.array(self.dick["Natural Constants"]["Pos sun"])
        self.u_sun=np.array(self.dick["Natural Constants"]["u_sun"])
        self.I_sp = self.dick["Satellite parameters"]["I_sp"]
        self.T = Thrust
        self.R_min =self.dick["Natural constants"]["Rad Earth"] + h_min
        self.R_max = self.dick["Natural constants"]["Rad Earth"] + h_max
        self.thrust = False
        self.dry_m = dry_m
        self.m_dot = Thrust / self.dick["Natural constants"]["g0"] / I_sp
        self.A_solar=A_solar
    def circular_set_up(self,h,m_wet):
        self.Y0=np.empty(7)
        R=self.dick["Natural constants"]["Rad Earth"]+h
        V=np.sqrt(self.dick["Natural constants"]["mu"]/R)

        u=np.array([1,0,0])
        u_bar=u/np.linalg.norm(u)
        self.Y0[:3]=u_bar*R

        v=np.array([0,1,1])
        v_bar=v/np.linalg.norm(v)
        self.Y0[3:6]=v_bar*V

        self.Y0[6]=m_wet
    def angle_between(self, u1,u2):
        return np.arccos(np.clip(np.dot(u1, u2), -1.0, 1.0))
    def dy_vector(self, t, Y):
        dy = np.empty(7)  # x,y,z, xdot,ydot,zdot ,m

        # Orbital parameters
        R = Y[:3]
        V = Y[3:6]
        m = Y[6]
        R_abs = np.linalg.norm(Y[:3])
        V_abs = np.linalg.norm(Y[3:6])
        a = self.dick["Natural constants"]["mu"] / (2 * self.dick["Natural constants"]["mu"] / R_abs - V_abs**2)  # semi major axis


        #Consequences of position and orbital parameters
        rho = 2.30E-11 # Add function for the density
        if (a <= self.R_min):
            self.thrust = True
        elif (a >= self.R_max):
            self.thrust = False
        sun_multiplier=abs(np.clip(np.dot(self.u_sun, V/V_abs), -1.0, 1.0))


        # Set up orbital forces
        Grav = -R * self.dick["Natural constants"]["mu"] / R_abs ** 3  # Acceleration due to gravity
        A=self.dick["Natural constants"]["A extended"]+self.A_solar*sun_multiplier
        f_aero = 0.5 * rho * V_abs ** 2 * A * self.dick["Natural constants"]["Cd"]

        # Set up dy
        dy[:3] = V
        if (self.thrust):
            dy[3:6] = (self.T - f_aero) * V / V_abs / m + Grav
            dy[6] = -self.m_dot
        else:
            dy[3:6] = -f_aero * V / V_abs / m + Grav
            dy[6] = 0

        print(f"t:{t} f_aero: {f_aero} Semi-major axis: {a}, Thrust: {self.thrust} Radius: {R_abs} Velocity: {V_abs}")
        return dy

    def determine_char(self, Y):
        self.R = Y[:3]
        self.V = Y[3:6]
        self.m = Y[6]
        self.R_abs = np.linalg.norm(Y[:3])
        self.V_abs = np.linalg.norm(Y[3:6])
    #def orientation(self):
    #   self.theta==

if __name__ == "__main__":
    LAMP=Satellite("Inputs/Config.toml",1000,1.5,54,1e-3,295e3,305e3)
    LAMP.circular_set_up(300e3,LAMP.dry_m+4)
    print(LAMP.Y0)
    sol=solve_ivp(LAMP.dy_vector,[0,10000],LAMP.Y0, method="Radau")
    print(sol.y.shape)

    LAMP.determine_char(sol.y[:,-1])
    print(f"Final result\n R: {LAMP.R_abs/1e3-LAMP.dick['Natural constants']['Rad Earth']/1e3}km V: {LAMP.V_abs/1e3}km")

