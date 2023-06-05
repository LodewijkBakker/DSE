import matplotlib.pyplot as plt
import numpy as np

def read_file(filename):
    a = []
    with open(f"Inputs/{filename}") as f:
        for line in f.readlines():
            line=line.strip("\n")
            l=line.split()
            if l[0]!="#" and l[0]!="#\n" :
                if l[0][0]!="#":
                    a.append(l[2:])
    array=np.array(a).astype(float)
    array[:, 0] -= array[0, 0]
    return array
R=6371.8
array1=read_file("oscar.oev")
perigee=array1[:,1]*(1-array1[:,2])-R
apogee=array1[:,1]*(1+array1[:,2])-R
plt.plot(array1[:,0],apogee,color="red", label="Apogee",linewidth=1.5)
plt.plot(array1[:,0],perigee, color="blue", label="Perigee", linewidth=1.5)
plt.title("DRAMA worst case deorbiting simulation")
# array2=read_file("oscar.oeb")
# plt.plot(array2[:,0],array2[:,1]-R, color="red", label="Apogee Nominal (Worst case)")
# array3=read_file("oscar.oew")
# #plt.plot(array3[:,0],array3[:,1]-R, color="blue")
plt.ylabel("Altitude [km]")
plt.xlabel("Days since last thrust")
plt.legend()
plt.savefig("Outputs/Deorbiting_graph.png")
plt.show()
