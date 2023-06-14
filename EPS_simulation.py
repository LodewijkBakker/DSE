####################################################################################################
# -------------------------------- EPS Exact Orbit Simulation -------------------------------------#
#            Program which simulates the behaviour of the EPS over 6 orbits & under 3LTAN          #
# -------------------------------------------------------------------------------------------------#
####################################################################################################

# ----------------------------------------- Imports -----------------------------------------------#
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from colored import fg
import os
import csv
import math
import pandas as pd

# ---------------------------------- Define spetial characters: -----------------------------------#
# --- Define Colors --- #
red = fg('red')
green = fg('green')
white = fg('white')

# ------------------------------------ EPS Simulation Class ---------------------------------------#

class EPS_Simulation_:
    def __init__(self):
        
        # --- Load Geommetry considerations, efficiency consideration EPS --- #
        self.A_back_panel = 220 * 350 / 1e6     # [m^2] for 1 fold
        self.A_side_panel = 200 * 356 / 1e6     # [m^2] for 1 fold
        self.A_top_panel_max = 329 * 359 / 1e6      # [m^2] for 1 fold
        self.A_top_panel_ant = (329-100-40) * (349-100-40) / 1e6      # [m^2] for 1 fold (- antenna size)
        self.n_pcdu = 0.96                                      # Pumpkin
        self.n_harness = 0.98                                   
        self.n_bat = 0.98                                       # Li-Ion NASA, SMAD
        self.n_cell = 0.32                                      # 32% NASA
        self.Esp_liion = 150                                    # [Wh/kg] NASA      
        self.rhosp_liion = 210                                  # [Wh/U] NASA
        # --- TESTED CONFIGURATION --- #
        self.num_back_folds = 3
        # self.num_side_folds = 2.5       # per side
        self.num_top_folds = 1          # per side

        # --- Load the Paths --- #
        self.path_dusk_dawn = './Inputs/SP_eff_Config-2_ltan18.csv'
        self.path_9h = './Inputs/SP_eff_Config-2_ltan9.csv'
        self.path_12h = './Inputs/SP_eff_Config-2_ltan12.csv'
        self.path_orbital_efficiency = './Inputs/Orbital_efficiencies.csv'

        # --- Load the Astrodynamical Parameters --- #
        self.Mission_life = 5                         # Mission Lifetime in years
        self.T_orbit = 5423                           # Orbital Period in seconds
        self.S_flux = 1322                            # Solar Flux in W/m^2 - COLD from THERMO SIMULATION

        self.Id = 0.88                                # SMAD

        self.Cell_d = 0.005                           # Cell Degradation per year 0.5% NASA
        
        self.Cell_D = (1-self.Cell_d) ** self.Mission_life   # Cell Degradation over Mission Lifetime

        ### --- Dusk-Dawn orbit (6h) --- ###
        eclipse_map_6h = np.array(self.get_col_from_csv(self.path_dusk_dawn, 13))
        n_incidence_6h_fside = np.array(self.get_col_from_csv(self.path_dusk_dawn, 19))
        n_incidence_6h_rback = np.array(self.get_col_from_csv(self.path_dusk_dawn, 21))
        n_incidence_6h_fback = np.array(self.get_col_from_csv(self.path_dusk_dawn, 22))
        n_incidence_6h_fside = np.where(n_incidence_6h_fside < 0, 0, n_incidence_6h_fside)
        n_incidence_6h_rback = np.where(n_incidence_6h_rback < 0, 0, n_incidence_6h_rback)
        n_incidence_6h_fback = np.where(n_incidence_6h_fback < 0, 0, n_incidence_6h_fback)
        temp = self.organise_orbit(eclipse_map_6h, eclipse_map_6h, n_incidence_6h_fside, n_incidence_6h_rback, n_incidence_6h_fback, length=self.T_orbit)
        self.eclipse_map_6h = temp[0]
        self.n_incidence_6h_fside = temp[1]
        self.n_incidence_6h_rback = temp[2]
        self.n_incidence_6h_fback = temp[3]
        ### --- 9h orbit --- ###
        eclipse_map_9h = np.array(self.get_col_from_csv(self.path_9h, 13))
        n_incidence_9h_fside = np.array(self.get_col_from_csv(self.path_9h, 19))
        n_incidence_9h_rback = np.array(self.get_col_from_csv(self.path_9h, 21))
        n_incidence_9h_fback = np.array(self.get_col_from_csv(self.path_9h, 22))
        n_incidence_9h_fside = np.where(n_incidence_9h_fside < 0, 0, n_incidence_9h_fside)
        n_incidence_9h_rback = np.where(n_incidence_9h_rback < 0, 0, n_incidence_9h_rback)
        n_incidence_9h_fback = np.where(n_incidence_9h_fback < 0, 0, n_incidence_9h_fback)
        temp = self.organise_orbit(eclipse_map_9h, eclipse_map_9h, n_incidence_9h_fside, n_incidence_9h_rback, n_incidence_9h_fback, length=self.T_orbit)
        self.eclipse_map_9h = temp[0]
        self.n_incidence_9h_fside = temp[1]
        self.n_incidence_9h_rback = temp[2]
        self.n_incidence_9h_fback = temp[3]
        ### --- 12h orbit --- ###
        eclipse_map_12h = np.array(self.get_col_from_csv(self.path_12h, 13))
        n_incidence_12h_fside = np.array(self.get_col_from_csv(self.path_12h, 19))
        n_incidence_12h_rback = np.array(self.get_col_from_csv(self.path_12h, 21))
        n_incidence_12h_fback = np.array(self.get_col_from_csv(self.path_12h, 22))
        n_incidence_12h_fside = np.where(n_incidence_12h_fside < 0, 0, n_incidence_12h_fside)
        n_incidence_12h_rback = np.where(n_incidence_12h_rback < 0, 0, n_incidence_12h_rback)
        n_incidence_12h_fback = np.where(n_incidence_12h_fback < 0, 0, n_incidence_12h_fback)
        temp = self.organise_orbit(eclipse_map_12h, eclipse_map_12h, n_incidence_12h_fside, n_incidence_12h_rback, n_incidence_12h_fback, length=self.T_orbit)
        self.eclipse_map_12h = temp[0]
        self.n_incidence_12h_fside = temp[1]
        self.n_incidence_12h_rback = temp[2]
        self.n_incidence_12h_fback = temp[3]
        
        # --- Base continuous power independant of scenario --- #
        # Base continuous power [W] per subsystem XXX_BC
        TCS_BC = 3.5
        Str_BC = 0
        GNS_BC = 0.2
        COMMS_BC = 2.1      # SEE COMMS - could make a video :D 
        SDH_BC = 4
        EPS_BC = 1          # only PCDU
        Prop_BC = 0         # Preheating and thermal maintenace via pipeheating and TCS
        ADCS_peakpower = max(self.ADCS_max_power_cont(POLES=True), self.ADCS_max_power_cont(POLES=False)) # [W] continuous max power

        self.base_c_power_all_sc = TCS_BC + Str_BC + GNS_BC + COMMS_BC \
                                 + SDH_BC + EPS_BC + Prop_BC + ADCS_peakpower                                  # [W]
        # Payload 
        self.Payload_cont = 10                          # [W]

        # --- Orbtial Powers --- #
        # Non continuous scenarios

        self.Prop_peak_power = 121                                   # [W] array of T_prop_peak values standby power 0W -> heating and thermal maintenance considered in thermal already


        ## Orbital Efficiency get data from CSV
        orb_ltan = np.array(self.get_col_from_csv(self.path_orbital_efficiency, 6))
        eclipse_ratio = np.array(self.get_col_from_csv(self.path_orbital_efficiency, 7))
        
        self.eclipse_ratio_6h = self.find_eclipse_idx(orb_ltan, eclipse_ratio, value=6)
        self.eclipse_ratio_9h = self.find_eclipse_idx(orb_ltan, eclipse_ratio, value=9)
        self.eclipse_ratio_12h = self.find_eclipse_idx(orb_ltan, eclipse_ratio, value=12)

        Pay_sun = 25                                                                                        # [W]
        self.Payload_sun_6h = Pay_sun * np.ones(int(self.T_orbit - self.eclipse_ratio_6h * self.T_orbit))   # [W]
        self.Payload_sun_9h = Pay_sun * np.ones(int(self.T_orbit - self.eclipse_ratio_9h * self.T_orbit))   # [W]
        self.Payload_sun_12h = Pay_sun * np.ones(int(self.T_orbit - self.eclipse_ratio_12h * self.T_orbit)) # [W]

        pass

    def Power_need(self, ltan, payload_cont=True):
        '''
            Function that outputs the power need of the satellite as a function of time over 6 orbits.

            Parameters:
                ltan (int)             -> local time of ascending node in hours (6, 9 or 12)
                payload_cont (bool)    -> True if payload is continuous, False if not

            Returns:
                Power need of the satellite as a function of time for 6 orbits [W]
            
            >>> Press 'q' to EXIT HELP mode <<<
        '''

        # --- Base continuous power independant of scenario --- #
        Base_Power_continuous = self.base_c_power_all_sc * np.ones(6*self.T_orbit) # [W]
        
        if payload_cont:
            Base_Power_continuous += self.Payload_cont * np.ones(6*self.T_orbit) # [W]
            if ltan == 6:
                Base_Power_continuous -= 5 * np.ones(6*self.T_orbit) # [W] COnstraining 12h to 5W continuously
        
        if not payload_cont:
            if ltan == 6:
                for i in range(6):
                    Base_Power_continuous[i*self.T_orbit:i*self.T_orbit+len(self.Payload_sun_6h)] += self.Payload_sun_6h
            elif ltan == 9:
                for i in range(6):
                    Base_Power_continuous[i*self.T_orbit:i*self.T_orbit+len(self.Payload_sun_9h)] += self.Payload_sun_9h
            elif ltan == 12:
                for i in range(6):
                    Base_Power_continuous[i*self.T_orbit:i*self.T_orbit+len(self.Payload_sun_12h)] += self.Payload_sun_12h

        # --- Propulsion system Contribution --- #
        Prop_power = self.Prop_peak_power * np.ones(845) # [W]
        mid_point_prop = math.floor(len(Prop_power)/2)
        # --- Inserting the propulsion into the timeline --- #
        
        if ltan == 6:
            # sun contribution
            mid_point_sun0 = math.floor(len(self.Payload_sun_6h)/2)
            temp_start = mid_point_sun0 - mid_point_prop
            if temp_start < 0:
                print(red + 'ERROR: Propulsion system starts before the satellite is in the sun during orbit 1' + white)
                print(red + '       Error in Power_need(ltan=6)' + white)
                exit()
            Base_Power_continuous[temp_start:temp_start+len(Prop_power)] += Prop_power

            # eclipse contribution
            mid_point_eclipse0 = math.floor(int(self.eclipse_ratio_6h * self.T_orbit)/2)
            temp_start = 4*self.T_orbit - mid_point_eclipse0 - mid_point_prop
            Base_Power_continuous[temp_start:temp_start+len(Prop_power)] += Prop_power
        
        elif ltan == 9:
            # sun contribution
            mid_point_sun0 = math.floor(len(self.Payload_sun_9h)/2)
            temp_start = mid_point_sun0 - mid_point_prop
            if temp_start < 0:
                print(red + 'ERROR: Propulsion system starts before the satellite is in the sun during orbit 1' + white)
                print(red + '       Error in Power_need(ltan=9)' + white)
                exit()
            Base_Power_continuous[temp_start:temp_start+len(Prop_power)] += Prop_power

            # eclipse contribution
            mid_point_eclipse0 = math.floor(int(self.eclipse_ratio_9h * self.T_orbit)/2)
            temp_start = 4*self.T_orbit - mid_point_eclipse0 - mid_point_prop
            Base_Power_continuous[temp_start:temp_start+len(Prop_power)] += Prop_power
        
        elif ltan == 12:
            # sun contribution
            mid_point_sun0 = math.floor(len(self.Payload_sun_12h)/2)
            temp_start = mid_point_sun0 - mid_point_prop
            if temp_start < 0:
                print(red + 'ERROR: Propulsion system starts before the satellite is in the sun during orbit 1' + white)
                print(red + '       Error in Power_need(ltan=12)' + white)
                exit()
            Base_Power_continuous[temp_start:temp_start+len(Prop_power)] += Prop_power

            # eclipse contribution
            mid_point_eclipse0 = math.floor(int(self.eclipse_ratio_12h * self.T_orbit)/2)
            temp_start = 4*self.T_orbit - mid_point_eclipse0 - mid_point_prop
            Base_Power_continuous[temp_start:temp_start+len(Prop_power)] += Prop_power

        # --- Insering EPS PCDU and harness losses --- #
        Base_Power_continuous = Base_Power_continuous / (self.n_pcdu * self.n_harness) # [W]
        return Base_Power_continuous

    def Power_intake(self, efficiency_side, efficiency_top, efficiency_back, num_side_folds):
        '''
            Function that computes the power intake based on the different configurations of hte solar panels.

            Parameters:
                efficiency_side (array) -> Efficiency of the solar panels on the side of the satellite as a function of time
                efficiency_top (array)  -> Efficiency of the solar panels on the top of the satellite as a function of time
                efficiency_back (array) -> Efficiency of the solar panels on the back of the satellite as a function of time
            
            Returns:
                Power intake of the solar panels as a function of time for 1 orbit [W]
            
            >>> Press 'q' to EXIT HELP mode <<<
        '''

        cst_value = self.S_flux * self.n_cell * self.Id * self.Cell_D
        a_side = 2 * self.A_side_panel * num_side_folds
        a_top = self.A_top_panel_ant * self.num_top_folds
        a_back = self.A_back_panel * self.num_back_folds

        P_in = cst_value * (a_side * efficiency_side + a_top * efficiency_top + a_back * efficiency_back) # [W]

        return P_in

    def organise_orbit(self, map, *args, length=5423):
        '''
            Function that reorders the output of one orbital period astodynamical parameters by begining at the end of eclipse.

            Parameters:
                *args (list)    -> Variable number of arrays to reorder
                map (array)     -> Array containing the booleans 0 (sun), 1 (eclipse), form which the arrays are going to get re-arranged.
                length          -> number of nodes (or seconds) in one orbital period

            Returns:
                all of the arrays in *args reorganised to [sun(l1), eclipse(l2)](length), l1+l2 = length
            
            >>> Press 'q' to EXIT HELP mode <<<
        '''
        # --- Get indexes of Eclipse i.e. map[i] = 1 --- #
        idx_eclipse = np.where((map>0.9) & (map<1.1))[0]

        # --- Get the indexes of the trantions between sun and eclipse --- #
        index_idx = []
        for i in range(len(idx_eclipse) -  1):
            if idx_eclipse[i+1] - idx_eclipse[i] != 1:
                index_idx.append(i)

        # --- Indentify start and stop for slicing --- #
        start = idx_eclipse[index_idx[0]] + 1
        stop = start + length

        # --- Slice the arrays --- #
        args = [arg[start:stop] for arg in args]

        return args

    def ADCS_max_power_cont(self, POLES=True):
        '''
            Function that calculates the base continuous power of the ADCS subsystem depending on LAMP's location over
            the poles.

            Parameters:
                POLES (boolean) -> True if LAMP is over the poles, False otherwise. Default: True
                                   

            Returns:
                ADCS_base_c_power (float)   -> Base continuous power of the ADCS subsystem [W]
            
            >>> Press 'q' to EXIT HELP mode <<<
        '''

        # --- Power per unit --- #
        # Sensors
        magn_meter = 0.75
        star_sensor = 1.4
        # Actuators
        magn_1_active = 1.8      # row
        magn_2_active = 1        # pitch
        magn_3_active = 0.5      # yaw

        rw_1_active = 4.5        # yaw, pitch
        rw_2_active = 0.65       # row
        

        magn_1_standby = 0.159
        magn_2_standby = 0
        magn_3_standby = 0

        rw_1_standby = 0.15
        rw_2_standby = 0.35
        
        # --- Quantities --- #
        magn_meter_quantity = 1
        star_sensor_quantity = 2

        magn_1 = 1
        magn_2 = 1
        magn_3 = 4

        rw_1 = 2
        rw_2 = 2

        # --- Power calculation --- #
        if POLES:
            # Formula considered over the poles to maintain NADIR:
            #        P = magn in row + rw in yaw & pitch & row + remaining in standby
            ADCS_base_c_power = magn_1_active * magn_1 \
                              + rw_1_active * rw_1 \
                              + magn_2_standby * magn_2 \
                              + magn_3_standby * magn_3 \
                              + rw_2_active * rw_2 \
                              + magn_meter_quantity * magn_meter \
                              + star_sensor_quantity * star_sensor
            
            return ADCS_base_c_power
        if not POLES:
            # Formula considered over the poles to maintain NADIR:
            #        P = magn in yaw & pitch + rw in yaw & pitch & row + remaining in standby
            ADCS_base_c_power = magn_2_active * magn_2 \
                              + magn_3_active * magn_3 \
                              + rw_2_active * rw_2 \
                              + magn_1_standby * magn_1 \
                              + rw_1_active * rw_1 \
                              + magn_meter_quantity * magn_meter \
                              + star_sensor_quantity * star_sensor
            
            return ADCS_base_c_power

    def ADCS_everything_on(self):
        '''
            Function that computes the maximum ADCS power consumption during payload operations, and propulsion.

            Parameters:
                None

            Returns:
                ADCS_max (float)    -> Maximum power consumption of the ADCS subsystem [W]

            >>> Press 'q' to EXIT HELP mode <<<
        '''

        # --- Power per unit --- #
        # Sensors
        magn_meter = 0.75
        star_sensor = 1.4
        # Actuators
        magn_1_active = 1.8      # row
        magn_2_active = 1        # pitch
        magn_3_active = 0.5      # yaw

        rw_1_active = 4.5        # yaw, pitch
        rw_2_active = 0.65       # row
        

        magn_1_standby = 0.159
        magn_2_standby = 0
        magn_3_standby = 0

        rw_1_standby = 0.15
        rw_2_standby = 0.35
        
        # --- Quantities --- #
        magn_meter_quantity = 1
        star_sensor_quantity = 2

        magn_1 = 1
        magn_2 = 1
        magn_3 = 4

        rw_1 = 2
        rw_2 = 2

        # --- Power calculation --- #
        P = magn_1_active * magn_1 \
          + rw_1_active * rw_1 \
          + magn_2_active * magn_2 \
          + magn_3_active * magn_3 \
          + rw_2_active * rw_2 \
          + magn_meter_quantity * magn_meter \
          + star_sensor_quantity * star_sensor
        
        return P
    
    def get_col_from_csv(self, path, col_number):
        '''
            Function that extracts a column from a csv file.

            Parameters:
                path (string)       -> Path to the csv file
                col_number (int)    -> Number of the column to extract

            Returns:
                col (array)          -> Array containing the values of the column
            
            >>> Press 'q' to EXIT HELP mode <<<
        '''
        column_values = []
        with open(path, 'r') as csvfile:
            reader = csv.reader(csvfile)
            next(reader)
            for row in reader:
                if col_number< len(row):
                    column_values.append(float(row[col_number]))
        return column_values

    def find_eclipse_idx(self, map, array, value):
        '''
            Function that finds the index of a value in the map and returns the value of the same index in the array.
            Used only for the eclipse fraction computations.

            Parameters:
                map (array)                  -> Array containing the ltan values
                array (array)                -> Array containing the eclipse fraction values
                value (float)                -> ltan of interest

            Returns:
                Eclipse_fraction (float)     -> Value of the eclipse fraction corresponding to the ltan of interest

            >>> Press 'q' to EXIT HELP mode <<<
        '''
        idx_interest = np.where((map < value+1e-3)&(map > value-1e-3))
        if len(idx_interest[0]) == 0:
            print(red + 'ERROR: ltan value not found in the map' + white)
            exit()
        return array[idx_interest[0][0]]
    
    
    def battery_size_computations(self, battery_array_cd, name, side_id,ltan, Pay_continuous=True):
        '''
            Function which takes in the battery charge discharge state array and computes feasability based on the charge, discharge states.
            Then computes the capacity, mass and volume. Also checks if the charge, discharge cycle work in more depth.

            Parameters:
                battery_array_cd (array)    -> Array containing the charge, discharge states of the battery
                name (string)               -> Name of the configuration analysed
                side_id (int)               -> number of side panels
                Pay_continuous (bool)       -> True if the payload is always on, False if the payload is on only during eclipse

            returns:
                array (array)               -> Array containing [Configuration, Capacity_prop, Mass_prop, Volume_prop, Capacity_remaining, Mass_remaining, Volume_remaining]
        '''

        # Check if feasable at all 
        charge = np.sum(positive for positive in battery_array_cd if positive > 0)
        discharge = np.sum(negative for negative in battery_array_cd if negative < 0)
        if (self.n_pcdu * charge) < abs(discharge):
            pass
        else:
            # Define the battery capacity used from the propulsion battery during sun operations only if in discharge
            midpoint_prop = math.floor((845)/2)
            if ltan == 6:
                mid_point_sun0 = math.floor(len(self.Payload_sun_6h)/2)
                temp_start = mid_point_sun0 - midpoint_prop
                C_needed_prop_sun = np.sum(battery_array_cd[temp_start:temp_start+845])
                if C_needed_prop_sun > 0:       # In case of a positive value it would indicate charginng
                    C_needed_prop_sun = 0
            elif ltan == 9:
                mid_point_sun0 = math.floor(len(self.Payload_sun_9h)/2)
                temp_start = mid_point_sun0 - midpoint_prop
                C_needed_prop_sun = np.sum(battery_array_cd[temp_start:temp_start+845])
                if C_needed_prop_sun > 0:       # In case of a positive value it would indicate charginng
                    C_needed_prop_sun = 0
            elif ltan == 12:
                mid_point_sun0 = math.floor(len(self.Payload_sun_12h)/2)
                temp_start = mid_point_sun0 - midpoint_prop
                C_needed_prop_sun = np.sum(battery_array_cd[temp_start:temp_start+845])
                if C_needed_prop_sun > 0:       # In case of a positive value it would indicate charginng
                    C_needed_prop_sun = 0

        
            # Compute the total needed capacity for operations of propulsion battery
            # Battery Li-Ion operates for 10 000 cycles - > 50% DOD: SOURCE - Analysis of On-Board Photovoltaics for a Battery Electric Bus and Their Impact on Battery Lifespan - Scientific Figure on ResearchGate. Available from: https://www.researchgate.net/figure/Depth-of-discharge-versus-cycle-life-of-the-lithium-ion-battery_fig4_318292540 [accessed 14 Jun, 2023]
            C_needed_prop = self.Prop_peak_power * 845 + abs(C_needed_prop_sun)     # [Ws]
            C_BOL_prop = C_needed_prop / (0.5 * self.n_bat)                         # [Ws]
            M_prop_bat = (C_BOL_prop / 3600) / self.Esp_liion                       # [kg] = [Ws / 3600] / [Wh/kg] = Wh / Wh/kg
            V_prop_bat = (C_BOL_prop / 3600) / self.rhosp_liion                     # [U] = [Ws / 3600] / [Wh/L] = Wh / Wh/U

            # Compute the remaining battrey power need for non-propulsion operations
            # Battery Li-Ion operates for 30 000 cycles - > 20% DOD: SOURCE - Analysis of On-Board Photovoltaics for a Battery Electric Bus and Their Impact on Battery Lifespan - Scientific Figure on ResearchGate. Available from: https://www.researchgate.net/figure/Depth-of-discharge-versus-cycle-life-of-the-lithium-ion-battery_fig4_318292540 [accessed 14 Jun, 2023]

            C_remainimg = (abs(discharge) - C_needed_prop)/6 # Per orbit
            C_BOL_remain = C_remainimg / (0.2 * self.n_bat)                         # [Ws]
            M_remain_bat = (C_BOL_remain / 3600) / self.Esp_liion                    # [kg] = [Ws / 3600] / [Wh/kg] = Wh / Wh/kg
            V_remain_bat = (C_BOL_remain / 3600) / self.rhosp_liion                  # [U] = [Ws / 3600] / [Wh/L] = Wh / Wh/U

            M_total = M_prop_bat + M_remain_bat
            V_total = V_prop_bat + V_remain_bat

            if (M_prop_bat + M_remain_bat) > 2 or (V_prop_bat + V_remain_bat) > 1.5:
                pass
            else:
                txt = name +'_Folds_per_side'+str(side_id) + 'Pay_continuous_' + str(Pay_continuous) + 'ltan_' + str(ltan) 
                return [txt, C_BOL_prop/3600, M_prop_bat, V_prop_bat, C_BOL_remain/3600, M_remain_bat, V_remain_bat, M_total, V_total]
                
    
    def EPS_simulation_run(self):
        '''
            All the results of the EPS simulation are stored in ./Outputs/EPS_simulation. This function runs the
            simulation and stores the results in the corresponding folder. All different possibilities are considered.
            Plots and csv files are generated. Power in is dependant of the configuration considered. Power needed by
            the Payload scenario considered.
            
            Configuration (str)    -> Configuration of the solar panels 'CX-XXx-XXx'
                                       >>> 'CI-Br-Sf_0' : configuration I, back panel rotates freely, side panels fixed at angle 0
                                       >>> 'CI-Br-Sf_i' : configuration I, back panel rotates freely, side panels fixed at angle i
                                       >>> 'CIII-Bfh_i-Sf_0' : configuration III, back panel fixed horizontally at angle i, side panels fixed at angle 0
                                       >>> 'CIII-Bfh_i-Sf_i' : configuration III, back panel fixed horizontally at angle i, side panels fixed at angle i
                                       >>> 'CIII-Bfv_i-Sf_0' : configuration III, back panel fixed vertically at angle i, side panels fixed at angle 0
                                       >>> 'CIII-Bfv_i-Sf_i' : configuration III, back panel fixed vertically at angle i, side panels fixed at angle i
                                       >>> 'CIV-Fr-Sf_i' : configuration IV, front panel rotate freely, side panels fixed at angle i
                                       >>> 'CIV-Ffh_i-Sf_i' : configuration IV, front panel fixed horizontally at angle i, side panels fixed at angle i

            Scenario (str)        -> Paylaod continuous/discontinuous

            LTAN (float)          -> LTAN of the orbit in hours (6, 9, 12)

            >>> Press 'q' to EXIT HELP mode <<<
        '''
        # --- Create the directories --- #
        if not os.path.exists('./Outputs/EPS_simulation'):
            os.makedirs('./Outputs/EPS_simulation')
        if not os.path.exists('./Outputs/EPS_simulation/Plots'):
            os.makedirs('./Outputs/EPS_simulation/Plots')
        if not os.path.exists('./Outputs/EPS_simulation/CSV'):
            os.makedirs('./Outputs/EPS_simulation/CSV')

        # --- Configuration data --- #

        configuration_name = ['CI-Br-Sf_0',
                              'CI-Br-Sf_i',
                              'CIII-Bfh_i-Sf_0',
                              'CIII-Bfh_i-Sf_i',
                              'CIII-Bfv_i-Sf_0',
                              'CIII-Bfv_i-Sf_i',
                              'CIV-Fr-Sf_i',
                              'CIV-Ffh_i-Sf_i']
        
        back_panel_efficiency_map_6h = [self.n_incidence_6h_rback,
                                        self.n_incidence_6h_rback,
                                        self.n_incidence_6h_fback,
                                        self.n_incidence_6h_fback,
                                        self.n_incidence_6h_fback,
                                        self.n_incidence_6h_fback,
                                        self.n_incidence_6h_rback,
                                        self.n_incidence_6h_fback]
        
        back_panel_efficiency_map_9h = [self.n_incidence_9h_rback,
                                        self.n_incidence_9h_rback,
                                        self.n_incidence_9h_fback,
                                        self.n_incidence_9h_fback,
                                        self.n_incidence_9h_fback,
                                        self.n_incidence_9h_fback,
                                        self.n_incidence_9h_rback,
                                        self.n_incidence_9h_fback]
        
        back_panel_efficiency_map_12h = [self.n_incidence_12h_rback,
                                         self.n_incidence_12h_rback,
                                         self.n_incidence_12h_fback,
                                         self.n_incidence_12h_fback,
                                         self.n_incidence_12h_fback,
                                         self.n_incidence_12h_fback,
                                         self.n_incidence_12h_rback,
                                         self.n_incidence_12h_fback]
        
        top_panel_efficiency_map_6h = [self.n_incidence_6h_fside, 
                                       self.n_incidence_6h_fside,
                                       self.n_incidence_6h_fside,
                                       self.n_incidence_6h_fside,
                                       self.n_incidence_6h_fside,
                                       self.n_incidence_6h_fside,
                                       self.n_incidence_6h_fside,
                                       self.n_incidence_6h_fside]
        
        top_panel_efficiency_map_9h = [self.n_incidence_9h_fside,
                                       self.n_incidence_9h_fside,
                                       self.n_incidence_9h_fside,
                                       self.n_incidence_9h_fside,
                                       self.n_incidence_9h_fside,
                                       self.n_incidence_9h_fside,
                                       self.n_incidence_9h_fside,
                                       self.n_incidence_9h_fside]
        
        top_panel_efficiency_map_12h = [self.n_incidence_12h_fside,
                                        self.n_incidence_12h_fside,
                                        self.n_incidence_12h_fside,
                                        self.n_incidence_12h_fside,
                                        self.n_incidence_12h_fside,
                                        self.n_incidence_12h_fside,
                                        self.n_incidence_12h_fside,
                                        self.n_incidence_12h_fside]
        
        side_panel_efficiency_map_6h = [self.n_incidence_6h_fside,
                                        self.n_incidence_6h_fback,
                                        self.n_incidence_6h_fside,
                                        self.n_incidence_6h_fback,
                                        self.n_incidence_6h_fside,
                                        self.n_incidence_6h_fback,
                                        self.n_incidence_6h_fback,
                                        self.n_incidence_6h_fback]
        
        side_panel_efficiency_map_9h = [self.n_incidence_9h_fside,
                                        self.n_incidence_9h_fback,
                                        self.n_incidence_9h_fside,
                                        self.n_incidence_9h_fback,
                                        self.n_incidence_9h_fside,
                                        self.n_incidence_9h_fback,
                                        self.n_incidence_9h_fback,
                                        self.n_incidence_9h_fback]
        
        side_panel_efficiency_map_12h = [self.n_incidence_12h_fside,
                                         self.n_incidence_12h_fback,
                                         self.n_incidence_12h_fside,
                                         self.n_incidence_12h_fback,
                                         self.n_incidence_12h_fside,
                                         self.n_incidence_12h_fback,
                                         self.n_incidence_12h_fback,
                                         self.n_incidence_12h_fback]
        
        

        # ------------------------------------------------------------------------------------------------------------- #
        # ================================ Battery Consideration Path ================================================= #

        if not os.path.exists('./Outputs/EPS_simulation/Battery_CSV'):
            os.makedirs('./Outputs/EPS_simulation/Battery_CSV')
        
        battery_array_csv = [['Config', 'C_BOL_prop', 'M_prop_bat', 'V_prop_bat', 'C_BOL_remain', 'M_remain_bat', 'V_remain_bat', 'M_total', 'V_total'],]

        # ------------------------------------------------------------------------------------------------------------- #
        # =============================== Continuous Payload Operations =============================================== #
        
        if not os.path.exists('./Outputs/EPS_simulation/CSV/Continuous'):
            os.makedirs('./Outputs/EPS_simulation/CSV/Continuous/6h/')
            os.makedirs('./Outputs/EPS_simulation/CSV/Continuous/9h/')
            os.makedirs('./Outputs/EPS_simulation/CSV/Continuous/12h/')
        if not os.path.exists('./Outputs/EPS_simulation/Plots/Continuous'):
            os.makedirs('./Outputs/EPS_simulation/Plots/Continuous/6h/')
            os.makedirs('./Outputs/EPS_simulation/Plots/Continuous/9h/')
            os.makedirs('./Outputs/EPS_simulation/Plots/Continuous/12h/')

        
        for i, config_name in tqdm(enumerate(configuration_name)):
            for j in [0.5, 1, 1.5, 2, 2.5]:
                # --- 6h orbit --- #
                power_consumption_6h = self.Power_need(ltan=6, payload_cont=True)
                power_intake_6h = self.Power_intake(side_panel_efficiency_map_6h[i], top_panel_efficiency_map_6h[i], back_panel_efficiency_map_6h[i],j)
                power_intake_6h = np.tile(power_intake_6h, 6)
                # Check len of power_intake_6h == len of power_consumption_6h
                if len(power_intake_6h) != len(power_consumption_6h):
                    print(red+'Error: power_intake_6h and power_consumption_6h have different length'+white)
                    exit()
                # Battery discharge, charge status
                battery_6h_contribution = power_intake_6h - power_consumption_6h
                temp_bat_6 = self.battery_size_computations(battery_array_cd=battery_6h_contribution, name=config_name, side_id= j, ltan=6, Pay_continuous=True)
                
                if temp_bat_6 != None:
                    battery_array_csv.append(temp_bat_6)

                # saving plot data: plotting power_consumption, power_intake, battery_6h_contribution.
                plt.plot(power_consumption_6h, label='Power Consumption')
                plt.plot(power_intake_6h, label='Power Generation')
                plt.plot(battery_6h_contribution, label='Battery Contribution (Charge/Discharge)')
                plt.legend()
                plt.savefig('./Outputs/EPS_simulation/Plots/Continuous/6h/'+config_name+'_Folds_per_side'+str(j)+'.png')
                plt.close()

                # saving csv data: power_consumption, power_intake, battery_6h_contribution. saved vertically
                df = pd.DataFrame({'Power Consumption': power_consumption_6h, 'Power Generation': power_intake_6h, 'Battery Contribution (Charge/Discharge)': battery_6h_contribution})
                df.to_csv('./Outputs/EPS_simulation/CSV/Continuous/6h/'+config_name+'_Folds_per_side'+str(j)+'.csv', index=False)

                # --- 9h orbit --- #
                power_consumption_9h = self.Power_need(ltan=9, payload_cont=True)
                power_intake_9h = self.Power_intake(side_panel_efficiency_map_9h[i], top_panel_efficiency_map_9h[i], back_panel_efficiency_map_9h[i], j)
                power_intake_9h = np.tile(power_intake_9h, 6)
                # Check len of power_intake_9h == len of power_consumption_9h
                if len(power_intake_9h) != len(power_consumption_9h):
                    print(red+'Error: power_intake_9h and power_consumption_9h have different length'+white)
                    exit()
                # Battery discharge, charge status
                battery_9h_contribution = power_intake_9h - power_consumption_9h
                temp_bat_9 = self.battery_size_computations(battery_array_cd=battery_9h_contribution, name=config_name, side_id= j, ltan=9, Pay_continuous=True)
                
                if temp_bat_9 != None:
                    battery_array_csv.append(temp_bat_9)

                # saving plot data: plotting power_consumption, power_intake, battery_9h_contribution.
                plt.plot(power_consumption_9h, label='Power Consumption')
                plt.plot(power_intake_9h, label='Power Generation')
                plt.plot(battery_9h_contribution, label='Battery Contribution (Charge/Discharge)')
                plt.legend()
                plt.savefig('./Outputs/EPS_simulation/Plots/Continuous/9h/'+config_name+'_Folds_per_side'+str(j)+'.png')
                plt.close()

                # saving csv data: power_consumption, power_intake, battery_9h_contribution. saved vertically
                df = pd.DataFrame({'Power Consumption': power_consumption_9h, 'Power Generation': power_intake_9h, 'Battery Contribution (Charge/Discharge)': battery_9h_contribution})
                df.to_csv('./Outputs/EPS_simulation/CSV/Continuous/9h/'+config_name+'_Folds_per_side'+str(j)+'.csv', index=False)

                # --- 12h orbit --- #
                power_consumption_12h = self.Power_need(ltan=12, payload_cont=True)
                power_intake_12h = self.Power_intake(side_panel_efficiency_map_12h[i], top_panel_efficiency_map_12h[i], back_panel_efficiency_map_12h[i], j)
                power_intake_12h = np.tile(power_intake_12h, 6)
                # Check len of power_intake_12h == len of power_consumption_12h
                if len(power_intake_12h) != len(power_consumption_12h):
                    print(red+'Error: power_intake_12h and power_consumption_12h have different length'+white)
                    exit()
                # Battery discharge, charge status
                battery_12h_contribution = power_intake_12h - power_consumption_12h
                temp_bat_12 = self.battery_size_computations(battery_array_cd=battery_12h_contribution, name=config_name, side_id= j, ltan=12, Pay_continuous=True)

                if temp_bat_12 != None:
                    battery_array_csv.append(temp_bat_12)

                # saving plot data: plotting power_consumption, power_intake, battery_12h_contribution.
                plt.plot(power_consumption_12h, label='Power Consumption')
                plt.plot(power_intake_12h, label='Power Generation')
                plt.plot(battery_12h_contribution, label='Battery Contribution (Charge/Discharge)')
                plt.legend()
                plt.savefig('./Outputs/EPS_simulation/Plots/Continuous/12h/'+config_name+'_Folds_per_side'+str(j)+'.png')
                plt.close()

                # saving csv data: power_consumption, power_intake, battery_12h_contribution. saved vertically
                df = pd.DataFrame({'Power Consumption': power_consumption_12h, 'Power Generation': power_intake_12h, 'Battery Contribution (Charge/Discharge)': battery_12h_contribution})
                df.to_csv('./Outputs/EPS_simulation/CSV/Continuous/12h/'+config_name+'_Folds_per_side'+str(j)+'.csv', index=False)

        # ------------------------------------------------------------------------------------------------------------- #
        # ============================== Discontinuous Payload Operations ============================================= #


        if not os.path.exists('./Outputs/EPS_simulation/Plots/Discontinuous'):
            os.makedirs('./Outputs/EPS_simulation/Plots/Discontinuous/6h/')
            os.makedirs('./Outputs/EPS_simulation/Plots/Discontinuous/9h/')
            os.makedirs('./Outputs/EPS_simulation/Plots/Discontinuous/12h/')
        if not os.path.exists('./Outputs/EPS_simulation/CSV/Discontinuous'):
            os.makedirs('./Outputs/EPS_simulation/CSV/Discontinuous/6h/')
            os.makedirs('./Outputs/EPS_simulation/CSV/Discontinuous/9h/')
            os.makedirs('./Outputs/EPS_simulation/CSV/Discontinuous/12h/')
        
        for i, config_name in tqdm(enumerate(configuration_name)):
            for j in [0.5, 1, 1.5, 2, 2.5]:
                # --- 6h orbit --- #
                power_consumption_6h = self.Power_need(ltan=6, payload_cont=False)
                power_intake_6h = self.Power_intake(side_panel_efficiency_map_6h[i], top_panel_efficiency_map_6h[i], back_panel_efficiency_map_6h[i], j)
                power_intake_6h = np.tile(power_intake_6h, 6)
                # Check len of power_intake_6h == len of power_consumption_6h
                if len(power_intake_6h) != len(power_consumption_6h):
                    print(red+'Error: power_intake_6h and power_consumption_6h have different length'+white)
                    exit()
                # Battery discharge, charge status
                battery_6h_contribution = power_intake_6h - power_consumption_6h
                temp_bat_6 = self.battery_size_computations(battery_array_cd=battery_6h_contribution, name=config_name, side_id= j, ltan=6, Pay_continuous=False)

                if temp_bat_6 != None:
                    battery_array_csv.append(temp_bat_6)

                # saving plot data: plotting power_consumption, power_intake, battery_6h_contribution.
                plt.plot(power_consumption_6h, label='Power Consumption')
                plt.plot(power_intake_6h, label='Power Generation')
                plt.plot(battery_6h_contribution, label='Battery Contribution (Charge/Discharge)')
                plt.legend()
                plt.savefig('./Outputs/EPS_simulation/Plots/Discontinuous/6h/'+config_name+'_Folds_per_side'+str(j)+'.png')
                plt.close()

                # saving csv data: power_consumption, power_intake, battery_6h_contribution. saved vertically
                df = pd.DataFrame({'Power Consumption': power_consumption_6h, 'Power Generation': power_intake_6h, 'Battery Contribution (Charge/Discharge)': battery_6h_contribution})
                df.to_csv('./Outputs/EPS_simulation/CSV/Discontinuous/6h/'+config_name+'_Folds_per_side'+str(j)+'.csv', index=False)

                # --- 9h orbit --- #
                power_consumption_9h = self.Power_need(ltan=9, payload_cont=False)
                power_intake_9h = self.Power_intake(side_panel_efficiency_map_9h[i], top_panel_efficiency_map_9h[i], back_panel_efficiency_map_9h[i], j)
                power_intake_9h = np.tile(power_intake_9h, 6)
                # Check len of power_intake_9h == len of power_consumption_9h
                if len(power_intake_9h) != len(power_consumption_9h):
                    print(red+'Error: power_intake_9h and power_consumption_9h have different length'+white)
                    exit()
                # Battery discharge, charge status
                battery_9h_contribution = power_intake_9h - power_consumption_9h
                temp_bat_9 = self.battery_size_computations(battery_array_cd=battery_9h_contribution, name=config_name, side_id= j, ltan=9, Pay_continuous=False)

                if temp_bat_9 != None:
                    battery_array_csv.append(temp_bat_9)


                # saving plot data: plotting power_consumption, power_intake, battery_9h_contribution.
                plt.plot(power_consumption_9h, label='Power Consumption')
                plt.plot(power_intake_9h, label='Power Generation')
                plt.plot(battery_9h_contribution, label='Battery Contribution (Charge/Discharge)')
                plt.legend()
                plt.savefig('./Outputs/EPS_simulation/Plots/Discontinuous/9h/'+config_name+'_Folds_per_side'+str(j)+'.png')
                plt.close()

                # saving csv data: power_consumption, power_intake, battery_9h_contribution. saved vertically
                df = pd.DataFrame({'Power Consumption': power_consumption_9h, 'Power Generation': power_intake_9h, 'Battery Contribution (Charge/Discharge)': battery_9h_contribution})
                df.to_csv('./Outputs/EPS_simulation/CSV/Discontinuous/9h/'+config_name+'_Folds_per_side'+str(j)+'.csv', index=False)

                # --- 12h orbit --- #
                power_consumption_12h = self.Power_need(ltan=12, payload_cont=False)
                power_intake_12h = self.Power_intake(side_panel_efficiency_map_12h[i], top_panel_efficiency_map_12h[i], back_panel_efficiency_map_12h[i], j)
                power_intake_12h = np.tile(power_intake_12h, 6)
                # Check len of power_intake_12h == len of power_consumption_12h
                if len(power_intake_12h) != len(power_consumption_12h):
                    print(red+'Error: power_intake_12h and power_consumption_12h have different length'+white)
                    exit()
                # Battery discharge, charge status
                battery_12h_contribution = power_intake_12h - power_consumption_12h
                temp_bat_12 = self.battery_size_computations(battery_array_cd=battery_12h_contribution, name=config_name, side_id= j, ltan=12, Pay_continuous=False)

                if temp_bat_12 != None:
                    battery_array_csv.append(temp_bat_12)

                # saving plot data: plotting power_consumption, power_intake, battery_12h_contribution.
                plt.plot(power_consumption_12h, label='Power Consumption')
                plt.plot(power_intake_12h, label='Power Generation')
                plt.plot(battery_12h_contribution, label='Battery Contribution (Charge/Discharge)')
                plt.legend()
                plt.savefig('./Outputs/EPS_simulation/Plots/Discontinuous/12h/'+config_name+'_Folds_per_side'+str(j)+'.png')
                plt.close()

                # saving csv data: power_consumption, power_intake, battery_12h_contribution. saved vertically
                df = pd.DataFrame({'Power Consumption': power_consumption_12h, 'Power Generation': power_intake_12h, 'Battery Contribution (Charge/Discharge)': battery_12h_contribution})
                df.to_csv('./Outputs/EPS_simulation/CSV/Discontinuous/12h/'+config_name+'_Folds_per_side'+str(j)+'.csv', index=False)

            # Save the battery array csv as each entry one row in the csv file and each [0][i] one column
            df = pd.DataFrame(battery_array_csv)
            df.to_csv('./Outputs/EPS_simulation/Battery_CSV/battery_sizings.csv', index=False, header=False)
            
        pass
    
    pass