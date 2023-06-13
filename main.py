####################################################################################################
# ---------------------------------- Main Running Module ------------------------------------------#
####################################################################################################

# ---------------------------------- Import Libraries ---------------------------------------------#
# from Battery_sizing import Battery_sizing
# from Solar_panel import Solar_panel_sizing
from EPS_simulation import EPS_Simulation_
# Bat_s = Battery_sizing()
# Bat_s.Bat_sizing()

# SP_s = Solar_panel_sizing()
# SP_s.Solar_Panel_size()

EPS_s = EPS_Simulation_()
EPS_s.EPS_simulation_run()






