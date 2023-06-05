####################################################################################################
# ---------------------------------- Battery sizing class: ----------------------------------------#
#  All the battery sizing and visualisaton functions are defined in this class.                    #
# -------------------------------------------------------------------------------------------------#
####################################################################################################

# ---------------------------------- Import libraries: --------------------------------------------#
import numpy as np
import toml
import matplotlib.pyplot as plt
from colored import fg
import json

# ---------------------------------- Define spetial characters: -----------------------------------#
# --- Define Colors --- #
red = fg('red')
green = fg('green')
white = fg('white')

# ---------------------------------- Battery sizing class: ----------------------------------------#
class Battery_sizing:
    def __init__(self):
        # --- Define Paths: --- #
        self.path_capacity_per_engine = './Outputs/Batteries/Capacity_per_engine[Wh].toml'
        self.path_capacity_per_battery = './Outputs/Batteries/Capacity_per_battery[Wh].toml'
        self.path_mass_per_engine = './Outputs/Batteries/Mass_per_engine[kg].toml'
        self.path_mass_per_battery = './Outputs/Batteries/Mass_per_battery[kg].toml'
        self.path_volume_per_engine = './Outputs/Batteries/Volume_per_engine[U].toml'
        self.path_volume_per_battery = './Outputs/Batteries/Volume_per_battery[U].toml'
        
        # --- Unpack the battery sizing parameters: --- #
        path_b = './Inputs/Bat_Param.toml' # ----------------------------------->
        self.b_cst = toml.load(path_b)

        # Import battery parameters:
        self.n_liion = self.b_cst["LiIon"]["efficiency"]
        self.DOD_liion = self.b_cst["LiIon"]["DOD"]
        self.Esp_liion = self.b_cst["LiIon"]["E_sp"]
        self.rhosp_liion = self.b_cst["LiIon"]["rho_sp"]

        self.n_lipoly = self.b_cst["LiPoly"]["efficiency"]
        self.DOD_lipoly = self.b_cst["LiPoly"]["DOD"]
        self.Esp_lipoly = self.b_cst["LiPoly"]["E_sp"]
        self.rhosp_lipoly = self.b_cst["LiPoly"]["rho_sp"]

        self.n_nicd = self.b_cst["NiCd"]["efficiency"]
        self.DOD_nicd = self.b_cst["NiCd"]["DOD"]
        self.Esp_nicd = self.b_cst["NiCd"]["E_sp"]
        self.rhosp_nicd = self.b_cst["NiCd"]["rho_sp"]

        self.n_nih2 = self.b_cst["NiH2"]["efficiency"]
        self.DOD_nih2 = self.b_cst["NiH2"]["DOD"]
        self.Esp_nih2 = self.b_cst["NiH2"]["E_sp"]
        self.rhosp_nih2 = self.b_cst["NiH2"]["rho_sp"]

        self.n_life = self.b_cst["LiFe"]["efficiency"]
        self.DOD_life = self.b_cst["LiFe"]["DOD"]
        self.Esp_life = self.b_cst["LiFe"]["E_sp"]
        self.rhosp_life = self.b_cst["LiFe"]["rho_sp"]

        # ----------------------------------->

        path_c = './Inputs/Cst_Param.toml' # ----------------------------------->
        self.c_cst = toml.load(path_c)

        # Import orbital parameters:
        self.T_o = self.c_cst["Orbital"]["T_orbit"]
        self.T_e = self.T_o * self.c_cst["Orbital"]["perc_eclipse"]

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

        self.Pay = self.c_cst["Payload"]["Continuous"] # Constrianing battery scenario: Payload is always ON
        self.T_pay = self.T_o

        self.Prop_tab = self.c_cst["Power Subsystems"]["Prop_P"] # [[Average, Peak],[Average, Peak],...,[Average, Peak]]
        self.T_prop_tab = self.c_cst["Power Subsystems"]["T_prop"] # [[ConfigI, ConfigII], [ConfigI, ConfigII], [ConfigI, ConfigII], ...], if above To -> use To
        
        self.EPS_P_tab = self.EPS_loss() + self.c_cst["Power Subsystems"]["EPS_P"] # [[ConfigI, ConfigII], [ConfigI, ConfigII], [ConfigI, ConfigII], ...] -> Peak = Average
        # ----------------------------------->


        # --- Error messages: --- #
        # Check if the number of configurations is the same for Prop, T_prop, and EPS_P:
        if len(self.Prop_tab) != len(self.T_prop_tab) or len(self.Prop_tab) != len(self.EPS_P_tab):
            print(f"{red}Error: The number of configurations is not the same for Prop, T_prop, and EPS_P.{white}")
            exit()
        # Check if the peak times are all positieve and smaller or equal to the orbital period for all subsystems:
        for i in range(len(self.T_prop_tab)):
            for j in range(len(self.T_prop_tab[i])):
                if self.T_prop_tab[i][j] < 0 or self.T_prop_tab[i][j] > self.T_o:
                    print(f"{red}Error: The peak times must be positieve and smaller or equal to the orbital period check Cst_Param.toml.{white}")
                    exit()
        for i in range(len(self.EPS_P_tab)):
            for j in range(len(self.EPS_P_tab[i])):
                if self.EPS_P_tab[i][j] < 0 or self.EPS_P_tab[i][j] > self.T_o:
                    print(f"{red}Error: The peak times must be positieve and smaller or equal to the orbital period check Cst_Param.toml.{white}")
                    exit()
        if self.T_ADCS < 0 or self.T_ADCS > self.T_o or self.T_SDH < 0 or self.T_SDH > self.T_o or self.T_COMMS < 0 or self.T_COMMS > self.T_o or self.T_GNS < 0 or self.T_GNS > self.T_o or self.T_SM < 0 or self.T_SM > self.T_o or self.T_TCS < 0 or self.T_TCS > self.T_o or self.T_pay < 0 or self.T_pay > self.T_o:
            print(f"{red}Error: The peak times must be positieve and smaller or equal to the orbital period check Cst_Param.toml.{white}")
            exit()

        pass

    def EPS_loss(self):
        ''' Function that computes the EPS losses and adds them to the EPS power consumption, 
            based on the different types of propulsion & configurations. '''

        # --- Compute total energy per orbit without EPS or Prop: --- #
        E_p = self.ADCS_P[1] * self.T_ADCS + self.SDH_P[1] * self.T_SDH \
            + self.COMMS_P[1] * self.T_COMMS + self.GNS_P[1] * self.T_GNS \
            + self.SM_P[1] * self.T_SM + self.TCS_P[1] * self.T_TCS + self.Pay[1] * self.T_pay
        
        E_a = self.ADCS_P[0] * (self.T_o - self.T_ADCS) + self.SDH_P[0] * (self.T_o - self.T_SDH) \
            + self.COMMS_P[0] * (self.T_o - self.T_COMMS) + self.GNS_P[0] * (self.T_o - self.T_GNS) \
            + self.SM_P[0] * (self.T_o - self.T_SM) + self.TCS_P[0] * (self.T_o - self.T_TCS) \
            + self.Pay[0] * (self.T_o - self.T_pay)
        
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
    
    def E_eclipse(self):
        ''' Function that computes the energy during eclipse, based on the different types of propulsion & configurations. '''
        
        # --- Compute total energy per eclipse without EPS or Prop: --- #
        E_i = 0 # Energy during eclipse without EPS or Prop

        time_peak = [self.T_ADCS, self.T_SDH, self.T_COMMS, self.T_GNS, self.T_SM, self.T_TCS, self.T_pay]
        power_a = [self.ADCS_P[0], self.SDH_P[0], self.COMMS_P[0], self.GNS_P[0], self.SM_P[0], self.TCS_P[0], self.Pay[0]]
        power_p = [self.ADCS_P[1], self.SDH_P[1], self.COMMS_P[1], self.GNS_P[1], self.SM_P[1], self.TCS_P[1], self.Pay[1]]
        time_prop_peak = self.T_prop_tab

        for i in range(len(time_peak)): #P1
            if time_peak[i] > self.T_e or time_peak[i] == self.T_e:
                E_i += power_p[i] * self.T_e
            if time_peak[i] < self.T_e:
                E_i += power_p[i] * time_peak[i] + power_a[i] * (self.T_e - time_peak[i])
        
        # --- Add total energy form Prop & EPS for the different configuration --- #
        # Initialise E_eclipse_tab: [[ConfigI, ConfigII], [ConfigI, ConfigII], [ConfigI, ConfigII], ...]
        E_eclipse_tab = np.zeros((len(self.Prop_tab), len(self.Prop_tab[0])))

        for i in range(len(self.Prop_tab)): #P2
            # --- Config I --- #
            if time_prop_peak[i][0] > self.T_e or time_prop_peak[i][0] == self.T_e:
                E_eclipse_tab[i][0] = E_i + self.Prop_tab[i][1] * self.T_e + self.EPS_P_tab[i][0] * self.T_e
            if time_prop_peak[i][0] < self.T_e:
                E_eclipse_tab[i][0] = E_i + self.Prop_tab[i][1] * time_prop_peak[i][0] + self.EPS_P_tab[i][0] * self.T_e \
                                    + self.Prop_tab[i][0] * (self.T_e - time_prop_peak[i][0])
                
            # --- Config II --- #
            if time_prop_peak[i][1] >= self.T_e:
                E_eclipse_tab[i][1] = E_i + self.Prop_tab[i][1] * self.T_e + self.EPS_P_tab[i][1] * self.T_e
            if time_prop_peak[i][1] < self.T_e:
                E_eclipse_tab[i][1] = E_i + self.Prop_tab[i][1] * time_prop_peak[i][1] + self.EPS_P_tab[i][1] * self.T_e \
                                    + self.Prop_tab[i][0] * (self.T_e - time_prop_peak[i][1])
        return E_eclipse_tab

    def Bat_sizing(self):
        ''' Function that computes the battery sizing for the different scenarios & different batteries,
            the output is saved in ./Outputs/Batteries/ '''
        
        # --- Required energy for eclipse --- #
        E_eclipse_tab = self.E_eclipse()

        # --- Required capacity per battery and per configuration --- #
        # Initialise the Capacity per configuration fro the types of battery
        C_liion = np.zeros((len(self.Prop_tab), len(self.Prop_tab[0])))
        C_lipoly = np.zeros((len(self.Prop_tab), len(self.Prop_tab[0])))
        C_life = np.zeros((len(self.Prop_tab), len(self.Prop_tab[0])))
        C_nicd = np.zeros((len(self.Prop_tab), len(self.Prop_tab[0])))
        C_nih2 = np.zeros((len(self.Prop_tab), len(self.Prop_tab[0])))

        for i in range(len(E_eclipse_tab)):
            for j in range(len(E_eclipse_tab[0])):
                C_liion[i][j] = E_eclipse_tab[i][j] / (self.n_liion * self.DOD_liion * 3600)
                C_lipoly[i][j] = E_eclipse_tab[i][j] / (self.n_lipoly * self.DOD_lipoly * 3600)
                C_life[i][j] = E_eclipse_tab[i][j] / (self.n_life * self.DOD_life * 3600)
                C_nicd[i][j] = E_eclipse_tab[i][j] / (self.n_nicd * self.DOD_nicd * 3600)
                C_nih2[i][j] = E_eclipse_tab[i][j] / (self.n_nih2 * self.DOD_nih2 * 3600)
        
        # --- Organise Data --- #
        C_engine = {f'Engine_{i}':None for i in range(len(self.Prop_tab))} 
        C_bat = {'LiIon':None, 'LiPoly':None, 'LiFe':None, 'NiCd':None, 'NiH2':None}

        # --- Organise data per engine --- # 
    
        for i in range(len(self.Prop_tab)):
            C_engine[f'Engine_{i}'] = [{'LiIon':C_liion[i]}, {'LiPoly':C_lipoly[i]}, {'LiFe':C_life[i]}, {'NiCd':C_nicd[i]}, {'NiH2':C_nih2[i]}]
        
        # --- Organise data per battery type --- # 
        names = ['LiIon', 'LiPoly', 'LiFe', 'NiCd', 'NiH2']
        Esp = [self.Esp_liion, self.Esp_lipoly, self.Esp_life, self.Esp_nicd, self.Esp_nih2]
        rhosp = [self.rhosp_liion, self.rhosp_lipoly, self.rhosp_life, self.rhosp_nicd, self.rhosp_nih2]
        for i in range(len(names)):
            C_bat[names[i]] = [{f'Engine_{j}': C_engine[f'Engine_{j}'][i][names[i]]} for j in range(len(self.Prop_tab))]

        # --- Save data --- #
        # ---------------------------------------------- Capacity ---------------------------------------------- #
        # Save data per engine in txt file
        with open(self.path_capacity_per_engine, 'w') as f:
            f.write(f'# Capacity in [Wh] per engine for the different batteries\n')
            f.write(f'# Structure -> [Configuration I, Configuration II]\n\n')
            for i in range(len(self.Prop_tab)):
                f.write(f'["Engine_{i}"]\n')
                f.write(f'LiIon = {[float(C_engine[f"Engine_{i}"][0]["LiIon"][0]), float(C_engine[f"Engine_{i}"][0]["LiIon"][1])]}\n')
                f.write(f'LiPoly = {[float(C_engine[f"Engine_{i}"][1]["LiPoly"][0]), float(C_engine[f"Engine_{i}"][1]["LiPoly"][1])]}\n')
                f.write(f'LiFe = {[float(C_engine[f"Engine_{i}"][2]["LiFe"][0]), float(C_engine[f"Engine_{i}"][2]["LiFe"][1])]}\n')
                f.write(f'NiCd = {[float(C_engine[f"Engine_{i}"][3]["NiCd"][0]), float(C_engine[f"Engine_{i}"][3]["NiCd"][1])]}\n')
                f.write(f'NiH2 = {[float(C_engine[f"Engine_{i}"][4]["NiH2"][0]), float(C_engine[f"Engine_{i}"][4]["NiH2"][1])]}\n\n')

        # Save data per battery type in txt file
        with open(self.path_capacity_per_battery, 'w') as f:
            f.write(f'# Capacity in [Wh] per battery type for the different engines\n')
            f.write(f'# Structure -> [Configuration I, Configuration II]\n')
            for i in range(len(names)):
                f.write(f'\n["{names[i]}"]\n')
                for j in range(len(self.Prop_tab)):
                    f.write(f'Engine_{j} = {[float(C_bat[names[i]][j][f"Engine_{j}"][0]), float(C_bat[names[i]][j][f"Engine_{j}"][1])]}\n')
        
        # ---------------------------------------------- Mass ---------------------------------------------- #
        # Save data per engine in txt file
        with open(self.path_mass_per_engine, 'w') as f:
            f.write(f'# Mass in [kg] per engine for the different batteries\n')
            f.write(f'# Structure -> [Configuration I, Configuration II]\n\n')
            for i in range(len(self.Prop_tab)):
                f.write(f'["Engine_{i}"]\n')
                f.write(f'LiIon = {[float(C_engine[f"Engine_{i}"][0]["LiIon"][0])/Esp[0], float(C_engine[f"Engine_{i}"][0]["LiIon"][1])/Esp[0]]}\n')
                f.write(f'LiPoly = {[float(C_engine[f"Engine_{i}"][1]["LiPoly"][0])/Esp[1], float(C_engine[f"Engine_{i}"][1]["LiPoly"][1])/Esp[1]]}\n')
                f.write(f'LiFe = {[float(C_engine[f"Engine_{i}"][2]["LiFe"][0])/Esp[2], float(C_engine[f"Engine_{i}"][2]["LiFe"][1])/Esp[2]]}\n')
                f.write(f'NiCd = {[float(C_engine[f"Engine_{i}"][3]["NiCd"][0])/Esp[3], float(C_engine[f"Engine_{i}"][3]["NiCd"][1])/Esp[3]]}\n')
                f.write(f'NiH2 = {[float(C_engine[f"Engine_{i}"][4]["NiH2"][0])/Esp[4], float(C_engine[f"Engine_{i}"][4]["NiH2"][1])/Esp[4]]}\n\n')
        
        # Save data per battery type in txt file
        with open(self.path_mass_per_battery, 'w') as f:
            f.write(f'# Mass in [kg] per battery type for the different engines\n')
            f.write(f'# Structure -> [Configuration I, Configuration II]\n')
            for i in range(len(names)):
                f.write(f'\n["{names[i]}"]\n')
                for j in range(len(self.Prop_tab)):
                    f.write(f'Engine_{j} = {[float(C_bat[names[i]][j][f"Engine_{j}"][0])/Esp[i], float(C_bat[names[i]][j][f"Engine_{j}"][1])/Esp[i]]}\n')
        
        # ---------------------------------------------- Volume ---------------------------------------------- #
        # Save data per engine in txt file
        with open(self.path_volume_per_engine, 'w') as f:
            f.write(f'# Volume in [U] per engine for the different batteries\n')
            f.write(f'# Structure -> [Configuration I, Configuration II]\n\n')
            for i in range(len(self.Prop_tab)):
                f.write(f'["Engine_{i}"]\n')
                f.write(f'LiIon = {[float(C_engine[f"Engine_{i}"][0]["LiIon"][0])/(rhosp[0]), float(C_engine[f"Engine_{i}"][0]["LiIon"][1])/(rhosp[0])]}\n')
                f.write(f'LiPoly = {[float(C_engine[f"Engine_{i}"][1]["LiPoly"][0])/(rhosp[1]), float(C_engine[f"Engine_{i}"][1]["LiPoly"][1])/(rhosp[1])]}\n')
                f.write(f'LiFe = {[float(C_engine[f"Engine_{i}"][2]["LiFe"][0])/(rhosp[2]), float(C_engine[f"Engine_{i}"][2]["LiFe"][1])/(rhosp[2])]}\n')
                f.write(f'NiCd = {[float(C_engine[f"Engine_{i}"][3]["NiCd"][0])/(rhosp[3]), float(C_engine[f"Engine_{i}"][3]["NiCd"][1])/(rhosp[3])]}\n')
                f.write(f'NiH2 = {[float(C_engine[f"Engine_{i}"][4]["NiH2"][0])/(rhosp[4]), float(C_engine[f"Engine_{i}"][4]["NiH2"][1])/(rhosp[4])]}\n\n')
        
        # Save data per battery type in txt file
        with open(self.path_volume_per_battery, 'w') as f:
            f.write(f'# Volume in [U] per battery type for the different engines\n')
            f.write(f'# Structure -> [Configuration I, Configuration II]\n')
            for i in range(len(names)):
                f.write(f'\n["{names[i]}"]\n')
                for j in range(len(self.Prop_tab)):
                    f.write(f'Engine_{j} = {[float(C_bat[names[i]][j][f"Engine_{j}"][0])/(rhosp[i]), float(C_bat[names[i]][j][f"Engine_{j}"][1])/(rhosp[i])]}\n')

        pass

    pass
            
        


