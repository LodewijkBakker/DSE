####################################################################################################
# ---------------------------------- Solar Panel sizing class: ------------------------------------#
#  All the battery sizing and visualisaton functions are defined in this class.                    #
# -------------------------------------------------------------------------------------------------#
####################################################################################################

# ---------------------------------- Import libraries: --------------------------------------------#
import numpy as np
import toml
import matplotlib.pyplot as plt
from colored import fg
import json

from Battery_sizing import Battery_sizing

# ---------------------------------- Define spetial characters: -----------------------------------#
# --- Define Colors --- #
red = fg('red')
green = fg('green')
white = fg('white')

bat_s = Battery_sizing()

# --------------------------------- Solar Panel sizing class: -------------------------------------#

class Solar_panel_sizing:
    def __init__(self):

        # --- Define paths --- #
        path_c = './Inputs/Cst_Param.toml'

        # --- Unpack Parameters --- # 
        self.c_cst = toml.load(path_c)

        # Import orbital parameters:
        self.T_o = self.c_cst["Orbital"]["T_orbit"]
        self.T_e = self.T_o * self.c_cst["Orbital"]["perc_eclipse"]
        self.T_s = self.T_o - self.T_e

        # Import EPS parameters:
        self.n_pcdu = self.c_cst["EPS efficiency"]["PCDU"]
        self.n_eps_loss = self.c_cst["EPS efficiency"]["Loss"]

        # Import the power profile: --> 1st line - [Average, Peak] [W], 2nd line T_peak [s]
        self.ADCS_P = self.c_cst["Power Subsystems"]["ADCS_P"]
        self.T_ADCS = self.c_cst["Power Subsystems"]["Peak Active ADCS"]

        self.SDH_P = self.c_cst["Power Subsystems"]["SDH_P"]
        self.T_SDH = self.c_cst["Power Subsystems"]["Peak Active SDH"]

        self.COMMS_P = self.c_cst["Power Subsystems"]["COMMS_P"]
        self.T_COMMS = self.c_cst["Power Subsystems"]["Peak Active COMMS"]

        self.GNS_P = self.c_cst["Power Subsystems"]["GNS_P"]   
        self.T_GNS = self.c_cst["Power Subsystems"]["Peak Active GNS"]

        self.SM_P = self.c_cst["Power Subsystems"]["SM_P"]
        self.T_SM = self.c_cst["Power Subsystems"]["Peak Active SM"]

        self.TCS_P = self.c_cst["Power Subsystems"]["TCS_P"]
        self.T_TCS = self.c_cst["Power Subsystems"]["Peak Active TCS"]

        self.Payload_d = self.c_cst["Payload"]["Discontinuous"]
        self.T_pay_d = self.T_s

        self.Prop_tab = self.c_cst["Power Subsystems"]["Prop_P"] # [[Average, Peak],[Average, Peak],...,[Average, Peak]]
        self.T_prop_tab = self.c_cst["Power Subsystems"]["T_prop"] # [[ConfigI, ConfigII], [ConfigI, ConfigII], [ConfigI, ConfigII], ...], if above To -> use To
        
        self.EPS_P_tab = self.EPS_loss() + self.c_cst["Power Subsystems"]["EPS_P"] # [[ConfigI, ConfigII], [ConfigI, ConfigII], [ConfigI, ConfigII], ...] -> Peak = Average

    pass

    def EPS_loss(self):
        ''' Function that computes the EPS losses and adds them to the EPS power consumption, 
            based on the different types of propulsion & configurations. '''

        # --- Compute total energy per orbit without EPS or Prop: --- #
        E_p = self.ADCS_P[1] * self.T_ADCS + self.SDH_P[1] * self.T_SDH \
            + self.COMMS_P[1] * self.T_COMMS + self.GNS_P[1] * self.T_GNS \
            + self.SM_P[1] * self.T_SM + self.TCS_P[1] * self.T_TCS + self.Payload_d[1] * self.T_pay_d
        
        E_a = self.ADCS_P[0] * (self.T_o - self.T_ADCS) + self.SDH_P[0] * (self.T_o - self.T_SDH) \
            + self.COMMS_P[0] * (self.T_o - self.T_COMMS) + self.GNS_P[0] * (self.T_o - self.T_GNS) \
            + self.SM_P[0] * (self.T_o - self.T_SM) + self.TCS_P[0] * (self.T_o - self.T_TCS) \
            + self.Payload_d[0] * (self.T_o - self.T_pay_d)
        
        # --- Add total energy form Prop for the different configuration --- #
        P_avg = np.zeros((len(self.Prop_tab), len(self.Prop_tab[0]))) # [[ConfigI, ConfigII], [ConfigI, ConfigII], [ConfigI, ConfigII], ...]

        for i in range(len(P_avg)):
            # Configuration I
            P_avg[i][0] = (E_a + E_p + self.Prop_tab[i][1] * self.T_prop_tab[i][0] \
                        + (self.T_o - self.T_prop_tab[i][0]) * self.Prop_tab[i][0]) / self.T_o
            
            # Configuration II
            P_avg[i][1] = (E_a + E_p + self.Prop_tab[i][1] * self.T_prop_tab[i][1] \
                        + (self.T_o - self.T_prop_tab[i][1]) * self.Prop_tab[i][0]) / self.T_o

        EPS_Loss_tab = np.zeros((len(P_avg), len(P_avg[0]))) # [[ConfigI, ConfigII], [ConfigI, ConfigII], [ConfigI, ConfigII], ...]

        # --- Compute the average power loss for each configuration --- #
        for i in range(len(EPS_Loss_tab)):
            EPS_Loss_tab[i][0] = P_avg[i][0] / (self.n_pcdu * self.n_eps_loss) - P_avg[i][0]
            EPS_Loss_tab[i][1] = P_avg[i][1] / (self.n_pcdu * self.n_eps_loss) - P_avg[i][1]

        return EPS_Loss_tab
    
    def E_orbit(self):
        ''' Function that computes the entire energy of one orbital period'''

        E_p = self.ADCS_P[1] * self.T_ADCS + self.SDH_P[1] * self.T_SDH \
            + self.COMMS_P[1] * self.T_COMMS + self.GNS_P[1] * self.T_GNS \
            + self.SM_P[1] * self.T_SM + self.TCS_P[1] * self.T_TCS + self.Payload_d[1] * self.T_pay_d
        

pass