######################################################################################################
# Description: This file contains the EPS class, which is used to simulate the power system of the   #
#              satellite, and the Solar Panel Sizing class.                                          #
######################################################################################################

# ----------------------------------- Imports -------------------------------------------------------#
import numpy as np
import toml 
import pandas as pd
import matplotlib.pyplot as plt
from colored import fg
green = fg('green')
white = fg('white')

class Solar_Array:
    def __init__(self,cst_file, sim_file):
        
        # Load constants file
        self.cst = toml.load(cst_file)

        # Load simulation file
        self.sim = pd.read_csv(sim_file)

        # Load orbital constants
        self.T_o = self.cst["Orbital"]["T_orbit"]                   
        f_eclipse = self.cst["Orbital"]["perc_eclipse"]
        self.S_flux = self.cst["Orbital"]["S_flux"]

        # Eclipse & Sun time in idealised cylindrical shadow
        self.T_e = self.T_o * f_eclipse
        self.T_s = self.T_o - self.T_e

        # Power Consumption
        self.P_prop = self.cst["Power Subsystems"]["Propulsion_P"]
        self.P_SM = self.cst["Power Subsystems"]["S&M_P"]
        self.P_ADCS = self.cst["Power Subsystems"]["ADCS_P"]
        self.P_TCS = self.cst["Power Subsystems"]["TCS_P"]
        self.P_EPS = self.cst["Power Subsystems"]["EPS_P"]
        self.P_SDH = self.cst["Power Subsystems"]["SDH_P"]
        self.P_COMMS = self.cst["Power Subsystems"]["COMMS_P"]
        self.P_pay = self.cst["Power Subsystems"]["Payload_P"]

        self.T_prop = self.cst["Power Subsystems"]["Peak Active Prop"]
        self.T_COMMS = self.cst["Power Subsystems"]["Peak Active COMMS"]
        self.T_ADCS = self.cst["Power Subsystems"]["Peak Active ADCS"]

        # Solar Panel Constants
        self.SA_eff = self.cst["Solar_efficiency"]["GaAs_ultra"]

        # Battery effciency
        self.B_eff = self.cst["Efficiencies Bat"]["n_LiIon"]

    
        print("\n------------- Configuration -------------------")
        print("Orbital Period = %.5f s" % self.T_o)
        print("Eclipse Time = %.5f s" % self.T_e)
        print("Sun Time = %.5f s" % self.T_s)
        print(r"Solar Flux = %lf W/m^2" % self.S_flux)
        print(" ----------------------------------------------")
        print("Propulsion Power = %.5f W" % self.P_prop[0][0])
        print("S&M Power = %.5f W" % self.P_SM[0][0])
        print("ADCS Power = %.5f W" % self.P_ADCS[0][0])
        print("TCS Power = %.5f W" % self.P_TCS[0][0])
        print("EPS Power = %.5f W" % self.P_EPS[0][0])
        print("SDH Power = %.5f W" % self.P_SDH[0][0])
        print("COMMS Power = %.5f W" % self.P_COMMS[0][0])
        print("Payload Power = %.5f W" % self.P_pay[0][0])
        print("----------------------------------------------")
        print("Peak Active Propulsion Power = %.5f W" % self.P_prop[0][1])
        print("Peak Active COMMS Power = %.5f W" % self.P_COMMS[0][1])
        print("Peak Active ADCS Power = %.5f W" % self.P_ADCS[0][1])
        print("Peak Active Propulsion Time = %.5f s" % self.T_prop)
        print("Peak Active COMMS Time = %.5f s" % self.T_COMMS)
        print("Peak Active ADCS Time = %.5f s" % self.T_ADCS)
        print("----------------------------------------------")
        print("Solar Array Efficiency = %.5f" % self.SA_eff)
        print("Battery Efficiency = %.5f" % self.B_eff)
        print("Angle Accuracy = %.5f" % self.cst["Sun Angle Accuracy"]["n_theta_fix"])
        print("----------------------------------------------")



        pass

    def E_out(self):
        P_nom = self.P_pay[0][0] + self.P_SM[0][0] + self.P_TCS[0][0] +  self.P_EPS[0][0] + self.P_SDH[0][0]
        E_peak = self.P_COMMS[0][1] * self.T_COMMS + self.P_prop[0][1] * self.T_prop + self.P_ADCS[0][1] * self.T_ADCS
        E_nom_f_peak = self.P_COMMS[0][0] * (self.T_o - self.T_COMMS) + self.P_prop[0][0] * (self.T_o - self.T_prop) + self.P_ADCS[0][0] * (self.T_o - self.T_ADCS)

        return P_nom * self.T_o + E_peak + E_nom_f_peak
    
    def Array_size(self):
        return self.E_out() / (self.S_flux * self.SA_eff * self.T_o/2 * self.B_eff * self.cst["Sun Angle Accuracy"]["n_theta_fix"])

    def plot_SA_efficiency_vs_Array_size(self):
        # Plot SA efficiency vs Array size
        SA_eff = np.linspace(0.1,0.4,100)
        Array_size = self.E_out() / (self.S_flux * SA_eff * self.T_s * self.B_eff * self.cst["Sun Angle Accuracy"]["n_theta_fix"])

        plt.plot(Array_size,SA_eff)
        plt.xlabel(r"Array Size ($m^2$)")
        plt.ylabel("Solar Array Efficiency")
        plt.title("Solar Array Efficiency vs Array Size")
        plt.show()
    
    def Power_SA(self):
        return self.Array_size() * self.S_flux * self.SA_eff * self.cst["Sun Angle Accuracy"]["n_theta_fix"]

    def Xtra_AR_from_SA(self):

        Area_1 = 0.63   # m^2
        E_SA_1 = self.T_o/2 * Area_1 * self.S_flux * self.SA_eff * self.B_eff * self.cst["Sun Angle Accuracy"]["n_theta_fix"]
        if self.E_out() > E_SA_1:
            xtr_array = (self.E_out() - E_SA_1) / ((self.S_flux * self.SA_eff * self.B_eff) * (self.T_o/2 * 0.99 + (self.T_o-self.T_e)/np.pi))
        else:
            xtr_array = 0
        return (green + "Extra Array Required = %.5f m^2 from the 0.63 m^2 basic" % xtr_array + white)


