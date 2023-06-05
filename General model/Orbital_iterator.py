from Orbital_component import Satellite, average_final_properties
import numpy as np
import pandas as pd
import toml
from scipy.integrate import solve_ivp
from astropy.constants import G, M_earth, R_earth
import vg
import pandas as pd
LAMP = Satellite("Inputs/Config-2.toml")
Orbital_period=LAMP.orbit_time()
print(Orbital_period)
sim_length=int(3*Orbital_period)

df=pd.read_csv("Inputs/Propulsion_table4.csv")
for i,row in df.iterrows():
    for a in range(1,3):
            #Satellite BOL conditions
            inp_file=toml.load(f"Inputs/Config-{a*2}.toml")
            inp_file['Satellite parameters']['I_sp']=float(row['Specific Impulse'])
            inp_file['Satellite parameters']['Thrust'] = float(row['Thrust']/1000)
            inp_file['Satellite parameters']['fuel mass']=float((np.exp(1500/(row['Specific Impulse']*9.80665))-1)*inp_file["Satellite parameters"]["dry mass"])
            input_file_name = f"Inputs/Config-BOL-{a*2}.toml"
            with open(input_file_name, "w") as toml_file:
                toml.dump(inp_file, toml_file)
            Sat1=Satellite(input_file_name)
            Sat1.orbital_el_set_up(sim_length)
            sol1 = solve_ivp(Sat1.dy_vector, Sat1.time_range, Sat1.Y0, method="Radau", max_step=1,
                            t_eval=Sat1.time_range, rtol=1)
            Sat1.Output.to_csv(f"Outputs/MidtermResult-BOL-Config{a*2}.csv")
            Sat1.orbit_time()

            #Satellite at EOL conditions
            inp_file['Satellite parameters']['fuel mass'] = 0.2
            input_file_name = f"Inputs/Config-EOL-{a*2}.toml"
            with open(input_file_name, "w") as toml_file:
                toml.dump(inp_file, toml_file)
            Sat2 = Satellite(input_file_name)
            Sat2.orbital_el_set_up(sim_length)
            sol2 = solve_ivp(Sat2.dy_vector, Sat2.time_range, Sat2.Y0, method="Radau", max_step=1,
                            t_eval=Sat1.time_range, rtol=1)
            Sat2.Output.to_csv(f"Outputs/MidtermResult-EOL-Config{a*2}.csv")
            #Sat1.final_properties(f"Outputs/MidtermResult-Config{a*2}.txt",sim_length,sol2,row)
            average_final_properties(f"Outputs/MidtermResult-Config{a*2}.txt",Sat1,Sat2, sim_length,sol1,sol2,row)
