####################################################################################################
# ---------------------------------- Main Running Module ------------------------------------------#
####################################################################################################

# ---------------------------------- Import Libraries ---------------------------------------------#
from Battery_sizing import Battery_sizing
from Solar_panel import Solar_panel_sizing

Bat_s = Battery_sizing()
Bat_s.Bat_sizing()

SP_s = Solar_panel_sizing()
SP_s.Solar_Panel_size()

