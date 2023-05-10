######################################################################################################
# Description: This file contains the EPS class, which is used to simulate the power system of the   #
#              satellite, and the Solar Panel Sizing class.                                          #
######################################################################################################

# ----------------------------------- Imports -------------------------------------------------------#
import numpy as np
import toml 
import pandas as pd
import matplotlib.pyplot as plt

class Solar_Array:
    def __init__(self,cst_file):
        
        # Load constants file
        self.cst = toml.load(cst_file)

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

        # Solar Panel Constants
        self.SA_eff = self.cst["Solar_efficiency"]["GaAs_ultra"]

        # Battery effciency
        self.B_eff = self.cst["Efficiencies Bat"]["n_LiIon"]

        pass

    def E_out(self):

        P_nom = self.P_pay[0][0] + self.P_prop[0][0] + self.P_SM[0][0] + self.P_ADCS[0][0] + self.P_TCS[0][0] + \
                self.P_EPS[0][0] + self.P_SDH[0][0] + self.P_COMMS[0][0]
        E_peak = self.P_COMMS[0][1] * self.T_COMMS + self.P_prop[0][1] * self.T_prop

        return P_nom * self.T_o + E_peak 
    
    def Array_size(self):
        return self.E_out() / (self.S_flux * self.SA_eff * self.T_s * self.B_eff * self.cst["Sun Angle Accuracy"]["n_theta_rot"])

    def plot_SA_efficiency_vs_Array_size(self):
        # Plot SA efficiency vs Array size
        SA_eff = np.linspace(0.1,0.4,100)
        Array_size = self.E_out() / (self.S_flux * SA_eff * self.T_s * self.B_eff * self.cst["Sun Angle Accuracy"]["n_theta_rot"])

        plt.plot(Array_size,SA_eff)
        plt.xlabel(r"Array Size ($m^2$)")
        plt.ylabel("Solar Array Efficiency")
        plt.title("Solar Array Efficiency vs Array Size")
        plt.show()
    
    def Power_SA(self):
        return self.Array_size() * self.S_flux * self.SA_eff * self.cst["Sun Angle Accuracy"]["n_theta_rot"]
