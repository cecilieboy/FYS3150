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

#def r(b, S, I, N=400, a=4, c=0.5):
#    return b * I - c * R

#%%
def runge_katta(stepsize, numb_steps, b = [1,2,3,4], N=400, a=4, c=0.5):


    for i in trange(len(b)):
        I_0 = 100
        S_0 = 300
        R_0 = 0
    
        S = [S_0]
        I = [I_0]
        R = [R_0]

        t = stepsize

        cutoff = numb_steps * stepsize

        timearray = np.arange(0, cutoff + stepsize, stepsize)
        
        while t <= cutoff: 
            k1_S = g(b[i], S[-1], I[-1]) * stepsize
            k1_I = h(b[i], S[-1], I[-1]) * stepsize
            k2_S = g(b[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I) * stepsize
            k2_I = h(b[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I) * stepsize
            k3_S = g(b[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k1_I) * stepsize
            k3_I = h(b[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k1_I) * stepsize
            k4_S = g(b[i], S[-1] + k3_S, I[-1] + k3_I) * stepsize
            k4_I = h(b[i], S[-1] + k3_S, I[-1] + k3_I) * stepsize

            S_val = S[-1] + 1/6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S)
            I_val = I[-1] + 1/6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I)

            S.append(S_val)
            I.append(I_val)
            R.append(N - S_val - I_val)

            t += stepsize
        

        Data = {'time':timearray, 'S':S, 'I':I, 'R': R}
        df = pd.DataFrame(Data)
        df.to_csv('Results/b=%i.csv'%(b[i]))

  

#%%
runge_katta(0.5,1000)
#%%


# %%
