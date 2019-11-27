#%%
import numpy as np 
from matplotlib import pyplot as plt
from tqdm import trange 
import pandas as pd
import random 
import seaborn as sns 
import matplotlib.pyplot as plt
#%%
def rhs_S(t, S, I, a= 4, b= 1, c= 0.5 , N=400):
    return c * (N - S - I) - a * S * I / N

def rhs_I(t, S, I, a= 4, b= 1, c= 0.5 , N=400):
    return a * S * I / N - b * I
#%%
def k_updates(RHS, X, RHS_kargs, stepsize):
    weights = np.array([1,2,2,1])
    temp = np.copy(X)
    print(temp)
    var_len = len(X)-1
    k = np.zeros(shape=(var_len, 4))
    for i in range(4):
        if i > 0:
            temp[0] = X[0] + stepsize/2
            temp[1:] = X[1:] + k[:, i-1]/2
        if i == 3:
            temp[0] = X[0] + stepsize
            temp[1:] = X[1:] + k[:, i-1]

        for j in range(var_len):
            k[j,i] = RHS[j](*temp, **RHS_kargs)

    return weights*k

def runge_kutta(stepsize, t_max, X_0, RHS, kargs_RHS, var_names =["S", "I"]):
    """
    flexible runge kutta for 
    X'(t) = RHS([t, x], **kargs)
    X is a vector i.e [S,I]

    Args:
        stepsize for Runge-Kutta
        t_max for Runge-Kutta
        X_0 is a list of initial values
        RHS a list of RHS functions rhs(t,x1,x2,...,x_n, **kargs_rhs),
        **kargs_RHS a list of keword arguments to each RHS entry
    Kargs:
        var_names: list of varibale names
    """
    num_steps = int(t_max / stepsize)
    X_t = np.zeros(shape = (len(X_0) +1, num_steps + 1))
    X_t[1:, 0]= X_0
    
    for i in range (1, num_steps+1):
        X_t[0, i] = i*stepsize
        k = k_updates(RHS, X_t[:,i -1 ], kargs_RHS, stepsize)
        X_t[1:, i] = X_t[1:, i - 1] + stepsize/6*np.sum(k, axis=1)
 
    return pd.DataFrame(X_t.T, columns=np.append(["time"], var_names)) 
#%%
df = runge_kutta(0.5, 50, [300, 100], [rhs_S, rhs_I], {} )
df["R"] = 400- df["S"] - df["I"]
plt.figure(figsize=(10,10))
for name in ["S", "I", "R"]:
    plt.plot(df["time"], df[name], label=name)
plt.legend(loc='best', fontsize = 28)
plt.xlabel("Time in a.u.", fontsize = 32)
plt.ylabel("Number", fontsize = 32)
plt.tick_params(size =24, labelsize=26)
plt.tight_layout()
plt.show()

#%%
def RHS_S(b, S, I, N=400, a=4, c=0.5):
    return c * (N - S - I) - a * S * I / N

def RHS_I(b, S, I, N=400, a=4, c=0.5):
    return a * S * I / N - b * I

#def r(b, S, I, N=400, a=4, c=0.5):
#    return b * I - c * R

def runge_kutta(stepsize, numb_steps, b = [1,2,3,4], N=400, a=4, c=0.5):


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
