######################################################################################################
#               Description: Program containing the Solar Array Sizing                               #
######################################################################################################

# ----------------------------------- Imports -------------------------------------------------------#
import numpy as np
import toml 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from colored import fg
green = fg('green')
white = fg('white')
red = fg('red')

class Solar_Array:
    def __init__(self, params_file, PARAM = True):
        
        # Unpack parameters file
        self.cst = toml.load(params_file)

        # Load orbital constants
        self.T_o = self.cst["Orbital"]["T_orbit"]                   
        f_eclipse = self.cst["Orbital"]["perc_eclipse"]
        self.S_flux = self.cst["Orbital"]["S_flux"]
        self.L = self.cst["Orbital"]["Mission_life"]

        # Eclipse & Sun time in idealised cylindrical shadow
        self.T_e = self.T_o * f_eclipse
        self.T_s = self.T_o - self.T_e

        # Power Consumptions
        self.P_prop = self.cst["Power Subsystems"]["Propulsion_P"]
        self.P_SM = self.cst["Power Subsystems"]["S&M_P"]
        self.P_ADCS = self.cst["Power Subsystems"]["ADCS_P"]
        self.P_TCS = self.cst["Power Subsystems"]["TCS_P"]
        self.P_EPS = self.cst["Power Subsystems"]["EPS_P"]
        self.P_SDH = self.cst["Power Subsystems"]["SDH_P"]
        self.P_COMMS = self.cst["Power Subsystems"]["COMMS_P"]
        self.P_pay = self.cst["Power Subsystems"]["Payload_P"]
        self.P_GNS = self.cst["Power Subsystems"]["GNS_P"]

        # Peak Power Consumption time
        self.T_prop = self.cst["Power Subsystems"]["Peak Active Prop"]
        self.T_COMMS = self.cst["Power Subsystems"]["Peak Active COMMS"]
        self.T_ADCS = self.cst["Power Subsystems"]["Peak Active ADCS"]
        self.T_TCS = self.cst["Power Subsystems"]["Peak Active TCS"]
        self.T_SDH = self.cst["Power Subsystems"]["Peak Active SDH"]
        self.T_pay = self.cst["Power Subsystems"]["Peak Active Pay"]
        self.T_SM = self.cst["Power Subsystems"]["Peak Active SM"]
        self.T_EPS = self.cst["Power Subsystems"]["Peak Active EPS"]
        self.T_GNS = self.cst["Power Subsystems"]["Peak Active GNS"]

        # Cell/Pannel parameters
        self.Id = self.cst["Cell efficiency"]["Inherent_degradation"]
        self.Cell_d = self.cst["Cell efficiency"]["cell_degradation"]
        self.Cell_D = (1-self.Cell_d) ** self.L

        if PARAM:
            print("------------- Configuration -------------------")
            print("Mission Life = %.5f years" % self.L)
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
            print("GNS Power = %.5f W" % self.P_GNS[0][0])
            print("----------------------------------------------")
            print("Peak Active Propulsion Power = %.5f W" % self.P_prop[0][1])
            print('Time Propulsion Peak Power = %.5f s' % self.T_prop)
            print("Peak Active COMMS Power = %.5f W" % self.P_COMMS[0][1])
            print('Time COMMS Peak Power = %.5f s' % self.T_COMMS)
            print("Peak Active ADCS Power = %.5f W" % self.P_ADCS[0][1])
            print('Time ADCS Peak Power = %.5f s' % self.T_ADCS)
            print("Peak Active TCS Power = %.5f W" % self.P_TCS[0][1])
            print('Time TCS Peak Power = %.5f s' % self.T_TCS)
            print("Peak Active EPS Power = %.5f W" % self.P_EPS[0][1])
            print('Time EPS Peak Power = %.5f s' % self.T_EPS)
            print("Peak Active SDH Power = %.5f W" % self.P_SDH[0][1])
            print('Time SDH Peak Power = %.5f s' % self.T_SDH)
            print("Peak Active Payload Power = %.5f W" % self.P_pay[0][1])
            print('Time Payload Peak Power = %.5f s' % self.T_pay)
            print("Peak Active S&M Power = %.5f W" % self.P_SM[0][1])
            print('Time S&M Peak Power = %.5f s' % self.T_SM)
            print("Peak Active GNS Power = %.5f W" % self.P_GNS[0][1])
            print('Time GNS Peak Power = %.5f s' % self.T_GNS)
            print("----------------------------------------------")
            print("Inherent Degradation = %.5f" % self.Id)
            print("Cell Degradation = %.5f" % self.Cell_d)
            print("Cell Degradation over Mission Life = %.5f" % self.Cell_D)
            print("----------------------------------------------\n")
        pass

    def E_out(self):
        ''' Function that computes the Energy per Orbit in [Ws]'''

        E_Peak = self.P_ADCS[0][1] * self.T_ADCS + self.P_COMMS[0][1] * self.T_COMMS + self.P_prop[0][1] * self.T_prop \
               + self.P_SDH[0][1] * self.T_SDH + self.P_pay[0][1] * self.T_pay + self.P_SM[0][1] * self.T_SM  \
               + self.P_TCS[0][1] * self.T_TCS + self.P_EPS[0][1] * self.T_EPS + self.P_GNS[0][1] * self.T_GNS
        
        E_nom = self.P_ADCS[0][0] * (self.T_o - self.T_ADCS) + self.P_COMMS[0][0] * (self.T_o - self.T_COMMS) + self.P_prop[0][0] * (self.T_o - self.T_prop) \
               + self.P_SDH[0][0] * (self.T_o - self.T_SDH) + self.P_pay[0][0] * (self.T_o - self.T_pay) + self.P_SM[0][0] * (self.T_o - self.T_SM)  \
               + self.P_TCS[0][0] * (self.T_o - self.T_TCS) + self.P_EPS[0][0] * (self.T_o - self.T_EPS) + self.P_GNS[0][0] * (self.T_o - self.T_GNS)
        return (E_nom + E_Peak)
        
    def P_avg(self, scenario = 0):
        ''' Function that computes the average power per orbit from E_out() that the solar array needs to provide in [W]'''
        if scenario == 0:
            # Exit code with error message
            print(f"{red}ERROR: Please specify a scenario for the average power computation{white}")
            exit()
        if scenario == 1:
            return self.E_out() / (self.T_o/2)
        if scenario == 2:
            return self.E_out() / (self.T_o/2)
        if scenario == 3 or scenario == 4:
            return self.E_out() / (self.T_s)
        if scenario == 5:
            return self.E_out() / (self.T_o)
    
    def disp_sc(self):
        print('--------------------------------------------------------------------------------------------------------------------')
        print('Scenario 1 description:')
        print('The power generation only occurs during half the orbit, the assumed incidence angle efficiency is 2/pi')
        print('--------------------------------------------------------------------------------------------------------------------')
        print('Scenario 2 description:')
        print('The power generation occurs during half the orbit, the assumed incidence angle is 23.5 degrees')
        print('--------------------------------------------------------------------------------------------------------------------')
        print('Scenario 3 description:')
        print('The power generation occurs during the entire sun orbit period, the assumed incidence efficiency is 0.99')
        print('--------------------------------------------------------------------------------------------------------------------')
        print('Scenario 4 description:')
        print('Scenario 3 with incidence angel of 23.5')
        print('--------------------------------------------------------------------------------------------------------------------')
        print('Scenario 5 description:')
        print('Constantly in sun with incidence efficiency of 0.99')
        print('--------------------------------------------------------------------------------------------------------------------')
        pass
    
    def SA_Area_I1(self, PLOT = False):
        ''' Solar Array area as a function of cell efficiency and in scenario with incidence angle 1:
            i.e. cos(theta) = 2/pi on average'''

        # Solar cell efficiency ranging from 16.8% to 32.2% in steps of 0.1%
        n_cell = np.arange(0.168, 0.322, 0.001)

        # Solar Array area
        SA_Area = self.P_avg(scenario=1) / (self.S_flux * n_cell * (2/np.pi) * self.Id * self.Cell_D)

        # Organise data 
        SA_Area_dict = {}
        SA_Area_dict[0] = {"cell efficiency": n_cell} 
        SA_Area_dict[1] = {"Area": SA_Area}
        SA_Area_dict[2] = {"units": f'm\u00b2'}

        # Plot
        if PLOT:
            # Efficiencies of specific cells
            n_SiS32 = 16.8
            n_XTE_SF = 32.2
            n_Z4J = 30.0
            n_TASC = 27.0

            plt.figure()
            plt.plot(n_cell*100, SA_Area, color='r', linewidth=2, linestyle='-')

            # plot vertical lines for the specific cells
            plt.axvline(x=n_SiS32, color='b', linestyle='--', linewidth=1, label='SiS32')
            plt.axvline(x=n_XTE_SF, color='g', linestyle='--', linewidth=1, label='XTE_SF')
            plt.axvline(x=n_Z4J, color='y', linestyle='--', linewidth=1, label='Z4J')
            plt.axvline(x=n_TASC, color='purple', linestyle='--', linewidth=1, label='TASC')
            
            # Compute the inetresction points of the vertical lines with the curve and plot a horizontal line
            x1 = np.interp(n_SiS32, n_cell*100, SA_Area)
            x2 = np.interp(n_XTE_SF, n_cell*100, SA_Area)
            x3 = np.interp(n_Z4J, n_cell*100, SA_Area)
            x4 = np.interp(n_TASC, n_cell*100, SA_Area)
            plt.axhline(y=x1, color='b', linestyle='--', linewidth=1)
            plt.axhline(y=x2, color='g', linestyle='--', linewidth=1)
            plt.axhline(y=x3, color='y', linestyle='--', linewidth=1)
            plt.axhline(y=x4, color='purple', linestyle='--', linewidth=1)

            plt.xlabel(r'Cell efficiency [%]')
            plt.ylabel(r'Solar Array Area [$m^2$]')
            plt.title('Solar Array Area vs. Cell Efficiency')
            plt.grid()
            plt.legend(loc='upper right')
            plt.title(r'Solar Array Area vs. Cell Efficiency, scenario 1')
            plt.savefig('./Outputs/SA_Area_I1.png')
            plt.show()
        
        return SA_Area_dict
    
    def SA_Area_I2(self, PLOT = False):
        ''' Solar Array area as a function of cell efficiency and in scenario with incidence angle 1: 
        i.e. cos(theta) = 2/pi on average'''

        # Solar cell efficiency ranging from 16.8% to 32.2% in steps of 0.1%
        n_cell = np.arange(0.168, 0.322, 0.001)

        # Solar Array area
        SA_Area = self.P_avg(scenario=2) / (self.S_flux * n_cell * (np.cos(np.radians(23.5))) * self.Id * self.Cell_D)

        # Organise data 
        SA_Area_dict = {}
        SA_Area_dict[0] = {"cell efficiency": n_cell} 
        SA_Area_dict[1] = {"Area": SA_Area}
        SA_Area_dict[2] = {"units": f'$m\u00b2$'}

        # Plot
        if PLOT:
            # Efficiencies of specific cells
            n_SiS32 = 16.8
            n_XTE_SF = 32.2
            n_Z4J = 30.0
            n_TASC = 27.0

            plt.figure()
            plt.plot(n_cell*100, SA_Area, color='r', linewidth=2, linestyle='-')

            # plot vertical lines for the specific cells
            plt.axvline(x=n_SiS32, color='b', linestyle='--', linewidth=1, label='SiS32')
            plt.axvline(x=n_XTE_SF, color='g', linestyle='--', linewidth=1, label='XTE_SF')
            plt.axvline(x=n_Z4J, color='y', linestyle='--', linewidth=1, label='Z4J')
            plt.axvline(x=n_TASC, color='purple', linestyle='--', linewidth=1, label='TASC')
            
            # Compute the inetresction points of the vertical lines with the curve and plot a horizontal line
            x1 = np.interp(n_SiS32, n_cell*100, SA_Area)
            x2 = np.interp(n_XTE_SF, n_cell*100, SA_Area)
            x3 = np.interp(n_Z4J, n_cell*100, SA_Area)
            x4 = np.interp(n_TASC, n_cell*100, SA_Area)
            plt.axhline(y=x1, color='b', linestyle='--', linewidth=1)
            plt.axhline(y=x2, color='g', linestyle='--', linewidth=1)
            plt.axhline(y=x3, color='y', linestyle='--', linewidth=1)
            plt.axhline(y=x4, color='purple', linestyle='--', linewidth=1)

            plt.xlabel(r'Cell efficiency [%]')
            plt.ylabel(r'Solar Array Area [$m^2$]')
            plt.title('Solar Array Area vs. Cell Efficiency')
            plt.grid()
            plt.legend(loc='upper right')
            plt.title(r'Solar Array Area vs. Cell Efficiency, scenario 2')
            plt.savefig('./Outputs/SA_Area_I2.png')
            plt.show()
        
        return SA_Area_dict

    def SA_Area_I3(self, PLOT = False):
            ''' Solar Array area as a function of cell efficiency and in scenario with incidence angle 1:
                i.e. cos(theta) = 2/pi on average'''

            # Solar cell efficiency ranging from 16.8% to 32.2% in steps of 0.1%
            n_cell = np.arange(0.168, 0.322, 0.001)

            # Solar Array area
            SA_Area = self.P_avg(scenario=3) / (self.S_flux * n_cell * (0.99) * self.Id * self.Cell_D)

            # Organise data 
            SA_Area_dict = {}
            SA_Area_dict[0] = {"cell efficiency": n_cell} 
            SA_Area_dict[1] = {"Area": SA_Area}
            SA_Area_dict[2] = {"units": f'm\u00b2'}

            # Plot
            if PLOT:
                # Efficiencies of specific cells
                n_SiS32 = 16.8
                n_XTE_SF = 32.2
                n_Z4J = 30.0
                n_TASC = 27.0

                plt.figure()
                plt.plot(n_cell*100, SA_Area, color='r', linewidth=2, linestyle='-')

                # plot vertical lines for the specific cells
                plt.axvline(x=n_SiS32, color='b', linestyle='--', linewidth=1, label='SiS32')
                plt.axvline(x=n_XTE_SF, color='g', linestyle='--', linewidth=1, label='XTE_SF')
                plt.axvline(x=n_Z4J, color='y', linestyle='--', linewidth=1, label='Z4J')
                plt.axvline(x=n_TASC, color='purple', linestyle='--', linewidth=1, label='TASC')

                # Compute the inetresction points of the vertical lines with the curve and plot a horizontal line
                x1 = np.interp(n_SiS32, n_cell*100, SA_Area)
                x2 = np.interp(n_XTE_SF, n_cell*100, SA_Area)
                x3 = np.interp(n_Z4J, n_cell*100, SA_Area)
                x4 = np.interp(n_TASC, n_cell*100, SA_Area)
                plt.axhline(y=x1, color='b', linestyle='--', linewidth=1)
                plt.axhline(y=x2, color='g', linestyle='--', linewidth=1)
                plt.axhline(y=x3, color='y', linestyle='--', linewidth=1)
                plt.axhline(y=x4, color='purple', linestyle='--', linewidth=1)

                plt.xlabel(r'Cell efficiency [%]')
                plt.ylabel(r'Solar Array Area [$m^2$]')
                plt.title('Solar Array Area vs. Cell Efficiency')
                plt.grid()
                plt.legend(loc='upper right')
                plt.title(r'Solar Array Area vs. Cell Efficiency, scenario 3')
                plt.savefig('./Outputs/SA_Area_I3.png')
                plt.show()

            return SA_Area_dict

    def SA_Area_I4(self, PLOT = False):
            ''' Solar Array area as a function of cell efficiency and in scenario with incidence angle 1:
                i.e. cos(theta) = 2/pi on average'''

            # Solar cell efficiency ranging from 16.8% to 32.2% in steps of 0.1%
            n_cell = np.arange(0.168, 0.322, 0.001)

            # Solar Array area
            SA_Area = self.P_avg(scenario=4) / (self.S_flux * n_cell * (np.cos(np.radians(23.5))) * self.Id * self.Cell_D)

            # Organise data 
            SA_Area_dict = {}
            SA_Area_dict[0] = {"cell efficiency": n_cell} 
            SA_Area_dict[1] = {"Area": SA_Area}
            SA_Area_dict[2] = {"units": f'm\u00b2'}

            # Plot
            if PLOT:
                # Efficiencies of specific cells
                n_SiS32 = 16.8
                n_XTE_SF = 32.2
                n_Z4J = 30.0
                n_TASC = 27.0

                plt.figure()
                plt.plot(n_cell*100, SA_Area, color='r', linewidth=2, linestyle='-')

                # plot vertical lines for the specific cells
                plt.axvline(x=n_SiS32, color='b', linestyle='--', linewidth=1, label='SiS32')
                plt.axvline(x=n_XTE_SF, color='g', linestyle='--', linewidth=1, label='XTE_SF')
                plt.axvline(x=n_Z4J, color='y', linestyle='--', linewidth=1, label='Z4J')
                plt.axvline(x=n_TASC, color='purple', linestyle='--', linewidth=1, label='TASC')

                # Compute the inetresction points of the vertical lines with the curve and plot a horizontal line
                x1 = np.interp(n_SiS32, n_cell*100, SA_Area)
                x2 = np.interp(n_XTE_SF, n_cell*100, SA_Area)
                x3 = np.interp(n_Z4J, n_cell*100, SA_Area)
                x4 = np.interp(n_TASC, n_cell*100, SA_Area)
                plt.axhline(y=x1, color='b', linestyle='--', linewidth=1)
                plt.axhline(y=x2, color='g', linestyle='--', linewidth=1)
                plt.axhline(y=x3, color='y', linestyle='--', linewidth=1)
                plt.axhline(y=x4, color='purple', linestyle='--', linewidth=1)

                plt.xlabel(r'Cell efficiency [%]')
                plt.ylabel(r'Solar Array Area [$m^2$]')
                plt.title('Solar Array Area vs. Cell Efficiency')
                plt.grid()
                plt.legend(loc='upper right')
                plt.title(r'Solar Array Area vs. Cell Efficiency, scenario 4')
                plt.savefig('./Outputs/SA_Area_I4.png')
                plt.show()

            return SA_Area_dict
    
    def SA_Area_I5(self, PLOT = False):
            ''' Solar Array area as a function of cell efficiency and in scenario with incidence angle 1:
                i.e. cos(theta) = 2/pi on average'''

            # Solar cell efficiency ranging from 16.8% to 32.2% in steps of 0.1%
            n_cell = np.arange(0.168, 0.322, 0.001)

            # Solar Array area
            SA_Area = self.P_avg(scenario=5) / (self.S_flux * n_cell * (0.99) * self.Id * self.Cell_D)

            # Organise data 
            SA_Area_dict = {}
            SA_Area_dict[0] = {"cell efficiency": n_cell} 
            SA_Area_dict[1] = {"Area": SA_Area}
            SA_Area_dict[2] = {"units": f'm\u00b2'}

            # Plot
            if PLOT:
                # Efficiencies of specific cells
                n_SiS32 = 16.8
                n_XTE_SF = 32.2
                n_Z4J = 30.0
                n_TASC = 27.0

                plt.figure()
                plt.plot(n_cell*100, SA_Area, color='r', linewidth=2, linestyle='-')

                # plot vertical lines for the specific cells
                plt.axvline(x=n_SiS32, color='b', linestyle='--', linewidth=1, label='SiS32')
                plt.axvline(x=n_XTE_SF, color='g', linestyle='--', linewidth=1, label='XTE_SF')
                plt.axvline(x=n_Z4J, color='y', linestyle='--', linewidth=1, label='Z4J')
                plt.axvline(x=n_TASC, color='purple', linestyle='--', linewidth=1, label='TASC')

                # Compute the inetresction points of the vertical lines with the curve and plot a horizontal line
                x1 = np.interp(n_SiS32, n_cell*100, SA_Area)
                x2 = np.interp(n_XTE_SF, n_cell*100, SA_Area)
                x3 = np.interp(n_Z4J, n_cell*100, SA_Area)
                x4 = np.interp(n_TASC, n_cell*100, SA_Area)
                plt.axhline(y=x1, color='b', linestyle='--', linewidth=1)
                plt.axhline(y=x2, color='g', linestyle='--', linewidth=1)
                plt.axhline(y=x3, color='y', linestyle='--', linewidth=1)
                plt.axhline(y=x4, color='purple', linestyle='--', linewidth=1)

                plt.xlabel(r'Cell efficiency [%]')
                plt.ylabel(r'Solar Array Area [$m^2$]')
                plt.title('Solar Array Area vs. Cell Efficiency')
                plt.grid()
                plt.legend(loc='upper right')
                plt.title(r'Solar Array Area vs. Cell Efficiency, scenario 4')
                plt.savefig('./Outputs/SA_Area_I5.png')
                plt.show()

            return SA_Area_dict
    

    def Bat_size_visual(self):
        ''' Function that sizes the battery based on eclipse'''
        time = [self.T_ADCS, self.T_pay, self.T_prop, self.T_SDH, self.T_COMMS, self.T_TCS, self.T_EPS, self.T_SM, self.T_GNS]
        P_nom = [self.P_ADCS[0][0], self.P_pay[0][0], self.P_prop[0][0], self.P_SDH[0][0], self.P_COMMS[0][0], self.P_TCS[0][0], self.P_EPS[0][0], self.P_SM[0][0], self.P_GNS[0][0]]
        P_peak = [self.P_ADCS[0][1], self.P_pay[0][1], self.P_prop[0][1], self.P_SDH[0][1], self.P_COMMS[0][1], self.P_TCS[0][1], self.P_EPS[0][1], self.P_SM[0][1], self.P_GNS[0][1]]
        E_peak = 0
        E_nom = 0
        for i in range(len(time)):
            if time[i] < self.T_e:
                E_peak += P_peak[i] * time[i]
                E_nom += P_nom[i] * (self.T_e - time[i])
            else:
                E_peak += P_peak[i] * self.T_e

        Capacity = E_peak + E_peak      # [Ws] Capacity
        Capacity_wh = Capacity / 3600      # [Wh] Capacity

        # Number of Cycles expected rounded to the upper integer
        N_c = np.ceil((self.L * (365*24*60*60) / self.T_o))        # [-] Number of cycles ~ 30,000

        # Assumed efficienncy of charging, discharging
        n_bat_Liion = 0.98
        n_bat_Nih2 = 0.85
        n_bat_NiCd = 0.80

        # Plotting Capacity as a function of DOD [%]
        DOD = np.linspace(10, 100, 90)

        C_Liion = Capacity / (DOD/100 * n_bat_Liion)
        C_Nih2 = Capacity / (DOD/100 * n_bat_Nih2)
        C_NiCd = Capacity / (DOD/100 * n_bat_NiCd)

        plt.figure()
        plt.plot(DOD, C_Liion/3600, color='r', linewidth=2, linestyle='-', label='Li-ion')
        plt.plot(DOD, C_Nih2/3600, color='b', linewidth=2, linestyle='-', label='NiH2')
        plt.plot(DOD, C_NiCd/3600, color='g', linewidth=2, linestyle='-', label='NiCd')

        # Ploting Vertical lines in relation to N_c, DOD_Liion = 20%, DOD_Nih2 = 65%, DOD_NiCd = 12% SMAD
        DOD_Liion = 20
        DOD_Nih2 = 80
        DOD_NiCd = 12

        plt.axvline(x=DOD_Liion, color='r', linestyle='--', linewidth=1)
        plt.axvline(x=DOD_Nih2, color='b', linestyle='--', linewidth=1)
        plt.axvline(x=DOD_NiCd, color='g', linestyle='--', linewidth=1)

        # Compute the inetresction points of the vertical lines with the curve and plot a horizontal line
        x1 = np.interp(DOD_Liion, DOD, C_Liion/3600)
        x2 = np.interp(DOD_Nih2, DOD, C_Nih2/3600)
        x3 = np.interp(DOD_NiCd, DOD, C_NiCd/3600)

        plt.axhline(y=x1, color='r', linestyle='--', linewidth=1)
        plt.axhline(y=x2, color='b', linestyle='--', linewidth=1)
        plt.axhline(y=x3, color='g', linestyle='--', linewidth=1)

        # Legend for the -- lines
        plt.plot([], [], color='black', linewidth=1, linestyle='--', label='DOD for ~30,000 cycles')

        plt.xlabel(r'DOD [%]')
        plt.ylabel(r'Capacity [Wh]')
        plt.title('Battery Capacity vs. DOD')
        plt.legend(loc='upper right')
        plt.savefig('./Outputs/Battery_capacity.png')
        plt.show()

        print('-------------------------------------------------------------------------------------')
        print('Battery parameters:')
        print(f'                  Capacity idealised (P_bat * t_dis): {Capacity} [Ws]')
        print(f'                  Capacity idealised (P_bat * t_dis): {Capacity_wh} [Wh]')
        print(f'                  Number of cycles: {N_c}')
        print(f'                  Capacity required Li-ion battery: {x1} [Wh]')
        print(f'                  Capacity required NiH2 battery: {x2} [Wh]')
        print(f'                  Capacity required NiCd battery: {x3} [Wh]')
        print('-------------------------------------------------------------------------------------')
        
        pass
    
