from Orbital_component import Satellite, average_final_properties
import numpy as np
import toml
import pandas as pd
from scipy.integrate import solve_ivp

#Used just to calculate the orbital period and extract the sim length from there
OP_sat= Satellite("Inputs/Config-2.toml")
Orbital_period=OP_sat.orbit_time()
print(Orbital_period)
sim_length=int(3*Orbital_period)
Save_df=pd.DataFrame(columns=["Config","LTAN","Eclipse_fraction","eff_fixed_side","eff_rotatable_side","eff_rotatable_back"])
#for conf in range(2):
conf=0
inp_file = toml.load(f"Inputs/Config-{(conf + 1) * 2}.toml")
#Equivalent to 0 to 90 (included )in steps of 7.5 if you divide it in two
for lan in range(0,180,15):
    if(lan<=90):
        ltan=12+int(lan/15)
    else:
        ltan=int(lan/15)
    input_file_name = f"Inputs/Config-{(conf + 1) * 2}_ltan{ltan}.toml"
    lan*=np.pi/180
    inp_file['Orbital parameters']['omega_LAN'] = lan
    with open(input_file_name, "w") as toml_file:
        toml.dump(inp_file, toml_file)
    Sat1 = Satellite(input_file_name)
    Sat1.orbital_el_set_up(sim_length)
    sol1 = solve_ivp(Sat1.dy_vector, Sat1.time_range, Sat1.Y0, method="Radau", max_step=1,
                     t_eval=Sat1.time_range, rtol=1)
    Sat1.Output.to_csv(f"Outputs/SP_eff_Config-{(conf + 1) * 2}_ltan{ltan}.csv")

    #Calculate efficiencies averaged over the whole orbit
    print(input_file_name)
    eff_f_side = Sat1.Output.eff_f_side.sum() / sim_length
    eff_r_side = Sat1.Output.eff_r_side.sum() / sim_length
    eff_r_back = Sat1.Output.eff_r_back.sum() / sim_length
    print(f"Conf={(conf+1)*2} LTAN={ltan}")
    print(f"Eff fixed side panel = {eff_f_side}")
    print(f"Eff rotatable side panel = {eff_r_side}")
    print(f"Eff rotatable back panel = {eff_r_back}")
    ecl_frac=Sat1.Output.Eclipse.sum()/sim_length
    Save_df.loc[len(Save_df)]= [(conf+1)*2,ltan,ecl_frac,eff_f_side,eff_r_side,eff_r_back]
    if ltan==18:
        Save_df.loc[len(Save_df)] = [(conf + 1) * 2, 6, ecl_frac, eff_f_side, eff_r_side, eff_r_back]
        Sat1.Output.to_csv(f"Outputs/SP_eff_Config-{(conf + 1) * 2}_ltan{6}.csv")
    print(Save_df)
# six=Save_df.loc[Save_df['LTAN'] == 18]
# six["LTAN"][0]=6
# print(six)
# Save_df.append(six)
Save_df1=Save_df.sort_values(["LTAN"],ascending=True)
print(Save_df1)
Save_df1.to_csv(f"Outputs/Orbital_efficiencies.csv")




