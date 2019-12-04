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
def a_runge_katta(stepsize, numb_steps, b = [1,2,3,4], N=400, a=4, c=0.5):


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

            k3_S = g(b[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I) * stepsize
            k3_I = h(b[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I) * stepsize

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
        df.to_csv('Results/RK_b=%i.csv'%(b[i]))


  
# %%
def S_to_I(b, S, I, R, stepsize, N=400, a=4, c=0.5):
    return a * S * I / N * stepsize

def I_to_R(b, S, I, R, stepsize, N=400, a=4, c=0.5):
    return b * I * stepsize 

def R_to_S(b, S, I, R, stepsize, N=400, a=4, c=0.5):
    return c * R * stepsize
#%%
def a_MC(numb_steps, b = [1,2,3,4], N=400, a=4, c=0.5):
    
    for i in trange (len(b)):
        stepsize = min([4 / (a * N), 1 / (b[i] * N), 1 / (c * N)])
        cutoff = numb_steps * stepsize
        timearray = np.arange(0, cutoff, stepsize)

        S = np.zeros(numb_steps)
        I = np.zeros(numb_steps)
        R = np.zeros(numb_steps)

        S[0] = 300
        I[0] = 100
        R[0] = 0

        for j in range(numb_steps - 1):
            psi_p = random.uniform(0,1)
            pir_p = random.uniform(0,1)
            prs_p = random.uniform(0,1)

            if S_to_I(b[i],S[j],I[j], R[j], stepsize) > psi_p:
                S[j + 1] = S[j] - 1
                I[j + 1] = I[j] + 1
            else:
                S[j + 1] = S[j]
                I[j + 1] = I[j]

            if I_to_R(b[i],S[j],I[j], R[j], stepsize) > pir_p:
                I[j + 1] = I[j + 1] - 1
                R[j + 1] = R[j] + 1
            else:
                I[j + 1] = I[j + 1]
                R[j + 1] = R[j]

            if R_to_S(b[i],S[j],I[j], R[j], stepsize) > pir_p:
                R[j + 1] = R[j + 1] - 1
                S[j + 1] = S[j + 1] + 1
            else:
                R[j + 1] = R[j + 1]
                S[j + 1] = S[j + 1]

        Data = {'time':timearray, 'S':S, 'I':I, 'R': R}
        df = pd.DataFrame(Data)
        df.to_csv('Results/MC_b=%i.csv'%(b[i]))
MC(100)

#%%
def gc(d_I, S, I, R, N=400, a=4, b = 1, c=0.5, e = 10, d = 3):
    return c * R - a * S * I / N - d * S + e * N

def hc(d_I, S, I, R, N=400, a=4, b = 1, c=0.5, e = 10, d = 3):
    return a * S * I / N - b * I - d * I - d_I * I 

def rc(d_I, S, I, R, N=400, a=4, b = 1, c=0.5, e = 10, d = 3):
    return b * I - c * R - d * R
#%%  

def c_runge_katta(stepsize, cutoff, d_I = [1,2,3,4,5,6], N=400, a=4, b =1, c=0.5, e = 10, d = 3):


    for i in trange(len(d_I)):
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
            k1_S = gc(d_I[i], S[-1], I[-1], R[-1]) * stepsize
            k1_I = hc(d_I[i], S[-1], I[-1], R[-1]) * stepsize
            k1_R = rc(d_I[i], S[-1], I[-1], R[-1]) * stepsize

            k2_S = gc(d_I[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, R[-1] + 1/2 * k1_R) * stepsize
            k2_I = hc(d_I[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, R[-1] + 1/2 * k1_R) * stepsize
            k2_R = rc(d_I[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, R[-1] + 1/2 * k1_R) * stepsize

            k3_S = gc(d_I[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k1_I, R[-1] + 1/2 * k2_R) * stepsize
            k3_I = hc(d_I[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k1_I, R[-1] + 1/2 * k2_R) * stepsize
            k3_R = rc(d_I[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k1_I, R[-1] + 1/2 * k2_R) * stepsize

            k4_S = gc(d_I[i], S[-1] + k3_S, I[-1] + k3_I, R[-1] + k3_R) * stepsize
            k4_I = hc(d_I[i], S[-1] + k3_S, I[-1] + k3_I, R[-1] + k3_R) * stepsize
            k4_R = rc(d_I[i], S[-1] + k3_S, I[-1] + k3_I, R[-1] + k3_R) * stepsize

            S_val = S[-1] + 1/6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S)
            I_val = I[-1] + 1/6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I)
            R_val = R[-1] + 1/6 * (k1_R + 2 * k2_R + 2 * k3_I + k4_R)

            S.append(S_val)
            I.append(I_val)
            R.append(R_val)

            t += stepsize
        

        Data = {'time':timearray, 'S':S, 'I':I, 'R': R}
        df = pd.DataFrame(Data)
        df.to_csv('Results/RKc)d_I=%i.csv'%(d_I[i]))

runge_katta(0.5, 40)



# %%
