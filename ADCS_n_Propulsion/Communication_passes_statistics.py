import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
lat=90
filename=f"90d_lat{lat}"
with open("Orbital pass data/"+filename+".txt") as f:
    count=0
    passes=[]
    for line in f.readlines():
        count+=1
        if(count>=13 and line!="\n"):
            line= re.sub(' +', ' ',line)
            passes.append(line.split(" "))
data=np.empty(int(len(passes)/3))
for i in range(0,len(passes),3):
    print(passes[i + 1][4])
    data[int(i/3)]=float(passes[i+1][4])
weights = np.ones_like(data)/len(data)
plt.hist(data,bins=15,range=[5,90], weights=weights)
x_limits = plt.xlim()  # Get the x-limits of the histogram
x_center = (x_limits[1] - x_limits[0]) / 2
plt.title(f"Distribution of tracking passes over 90 days [Lat {lat}°]")
plt.text(x_center, -6, f'Total number of passes: {int(len(data))}',size=12, ha='center')
plt.savefig("Orbital_pass_graphs/hist_"+filename+"_normalized.png",dpi=400)
plt.show()

plt.hist(data,bins=15,range=[5,90])
plt.title(f"Distribution of tracking passes over 90 days [Lat {lat}°]")
x_limits = plt.xlim()  # Get the x-limits of the histogram
x_center = (x_limits[1] - x_limits[0]) / 2
plt.text(x_center, -160, f'Total number of passes: {int(len(data))}',size=12, ha='center')
# Adjust the plot layout
plt.tight_layout()
plt.savefig("Orbital_pass_graphs/hist_"+filename+".png",dpi=400)
plt.show()
