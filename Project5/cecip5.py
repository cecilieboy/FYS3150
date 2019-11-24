#%%
import numpy as np 
from matplotlib import pyplot as plt
from tqdm import trange 
import pandas as pd
import random 
#%%
def g(b, S, I, N=400, a=4, c=0.5):
    return c * (N - S - I) - a * S * I / N

def h(b, S, I, N=400, a=4, c=0.5):
    return a * S * I / N - b * I

def r(b, S, I, N=400, a=4, c=0.5):
    return b * I - c * R



#%%
def runge_katta(stepsize, cutoff, b = [1,2,3,4], N=400, a=4, c=0.5):


    for i in range(len(b)):
        I_0 = 300
        S_0 = 100
        R_0 = 0
    
        S = [S_0]
        I = [I_0]
        R = [R_0]

        t = 0
        

        while t < cutoff: 
            k1_S = g(b[i], S[-1], I[-1])
            k1_I = h(b[i], S[-1], I[-1])
            k2_S = g(b[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I)
            k2_I = h(b[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I)
            k3_S = g(b[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k1_I)
            k3_I = h(b[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k1_I)
            k4_S = g(b[i], S[-1] + k3_S, I[-1] + k3_I)
            k4_I = h(b[i], S[-1] + k3_S, I[-1] + k3_I)

            S_val = g(b[i], S[-1], I[-1]) + 1/6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S) * stepsize
            I_val = h(b[i], S[-1], I[-1]) + 1/6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I) * stepsize

            S.append(S_val)
            I.append(I_val)
            R.append(N - S_val - I_val)


            Data = {'S':S, 'I':I, 'R': R}
            df = pd.DataFrame(Data)
            df.to_csv('Results/b_%i.csv'%(b[i]))

            t += 1 




runge_katta(1,100)

# %%
