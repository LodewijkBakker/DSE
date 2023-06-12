import matplotlib.pyplot as plt
import  numpy as np
import pandas as pd
from FIxed_back_panel_efficiency import max_eff_back_fixed
data=pd.read_csv("Outputs/Orbital_efficiencies.csv")
data=data.sort_values(["LTAN"],ascending=True)
print(data["LTAN"].to_numpy())
fig, ax1 = plt.subplots()
ln1=ax1.plot(data["LTAN"].to_numpy(),data["eff_fixed_side"].to_numpy(),label='Fixed side wings')
ln2=ax1.plot(data["LTAN"].to_numpy(),data["eff_rotatable_back"].to_numpy(),label="Rotatable back wing")
ln3=ax1.plot(data["LTAN"].to_numpy(),data["eff_rotatable_side"].to_numpy(),label="Rotatable side wing")
#eff_f_back=max_eff_back_fixed()
#data["eff_fixed_back"]=eff_f_back[:,2]
#data["Angle_max_eff"]=eff_f_back[:,1]
data.to_csv("Outputs/Orbital_efficiencies.csv")
ln4=ax1.plot(data["LTAN"].to_numpy(),data["eff_fixed_back"].to_numpy(),label="Fixed back wing")
ax2 = ax1.twinx()
ln5=ax2.plot(data["LTAN"].to_numpy(),data["Eclipse_fraction"].to_numpy(),color="black",label="Eclipse fraction")
#add legend
lns=ln1+ln2+ln3+ln4+ln5
labs = [l.get_label() for l in lns]
#ax1.legend(lns, labs, loc=0)

ax1.legend(loc="lower center")
ax2.legend(loc="upper center")
ax1.set_ylim(0.1, 0.75)
ax2.set_ylim(0.1, 0.75)
ax1.set_xlabel("Local time of ascention [h]")
ax1.set_ylabel("Incidence angle efficiency")
ax2.set_ylabel("Eclipse fraction")
plt.title("Average cos (incidence angle) over an orbit")
plt.show()