import matplotlib.pyplot as plt
import  numpy as np
import pandas as pd
from FIxed_back_panel_efficiency import max_eff_back_fixed
data=pd.read_csv("Outputs/Orbital_efficiencies.csv")
data=data.sort_values(["LTAN"],ascending=True)
print(data["LTAN"].to_numpy())
fig, ax1 = plt.subplots()
ax1.plot(data["LTAN"].to_numpy()[:-1],data["eff_fixed_side"].to_numpy()[:-1],label='Fixed side wings')
ax1.plot(data["LTAN"].to_numpy()[:-1],data["eff_rotatable_back"].to_numpy()[:-1],label="Rotatable back wing")
ax1.plot(data["LTAN"].to_numpy()[:-1],data["eff_rotatable_side"].to_numpy()[:-1],label="Rotatable side wing")
eff_f_back=max_eff_back_fixed()
ax1.plot(eff_f_back[:-1,0],eff_f_back[:-1,2],label="Fixed back wing")
ax2 = ax1.twinx()
ax2.plot(data["LTAN"].to_numpy()[:-1],data["Eclipse_fraction"].to_numpy()[:-1],color="black",label="Eclipse fraction")
ax1.legend()
ax1.set_ylim(0.1, 0.7)
ax2.set_ylim(0.1, 0.7)
ax1.set_xlabel("Local time of ascention [h]")
ax1.set_ylabel("Incidence angle efficiency")
ax2.set_ylabel("Eclipse fraction")
plt.title("Average cos (incidence angle) over an orbit")
plt.show()