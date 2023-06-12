import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
def normalize_rows(p):
    sum_of_rows = p.sum(axis=1)
    return p / sum_of_rows[:, np.newaxis]
def rowise_dot(p,q):
    return np.sum(p*q,axis=1)

def max_eff_back_fixed():
    orb_max=np.empty((13,3))#0 is the LTAN 1 is the angle 2 is the eff
    for h in range(6,19):
        filename=f'Outputs/SP_eff_Config-2_ltan{h}.csv'
        orbit=pd.read_csv(filename)
        r=orbit[["x","y","z"]].to_numpy()
        r=normalize(r)
        v=orbit[["v_x","v_y","v_z"]].to_numpy()
        v=normalize(v)
        #sum=np.sum(r,1)
        u_sun=np.array([0.9171,0.0,0.399])

        day=1-orbit["Eclipse"].to_numpy()
        print(f"Eclipse_frac={1-np.mean(day)}")
        rxv=np.cross(r,v)
        ang_eff=np.empty(360)
        result=np.empty((r.shape[0],360))
        for i in range(0,360):
            #ang=10*np.pi/180
            ang=i*np.pi/180
            vdrxv=rowise_dot(v,rxv)
            rxvxv=np.cross(v,rxv)
            SP_norm=rxv*np.cos(ang)+rxvxv*np.sin(ang)+np.multiply(v, vdrxv[:, np.newaxis])*(1-np.cos(ang))
            eff=np.maximum(np.dot(SP_norm,u_sun),0)
            res=eff*day
            result[:,i]=res
            ang_eff[i]=np.mean(res)
            #print(f"For angle {i} the average eff is {np.mean(res)}")
        pos = np.argmax(ang_eff)
        orb_max[h - 6][0]=h
        orb_max[h-6][1]=pos
        orb_max[h-6][2]=ang_eff[pos]
        orbit["eff_f_back"]=result[:,pos]
        orbit.to_csv(filename)
        #print(f"For {h}h orbit the maximum eff is at angle {pos} is {ang_eff[pos]}")
    print(orb_max)
    return orb_max

if __name__ == '__main__':
    orb_max=max_eff_back_fixed()
    plt.plot(orb_max[:,0],orb_max[:,2])
    plt.show()

