#%%
import numpy as np 
from matplotlib import pyplot as plt
from tqdm import trange 
import pandas as pd
import random 
import seaborn as sns 
import matplotlib.pyplot as plt
#%%
def rhs_S(t, S, I, a= 4, b= 1, c= 0.5, A =0, omega =1, f = 0, N=400):
    a_t = max(0.5,A*np.cos(omega*t) +a)
    f_t = max(0,f*np.cos((omega-0.05)*t + 0.05))
    return c * (N - S - I) - a_t * S * I / N - f_t

def rhs_I(t, S, I, a= 4, b= 1, c= 0.5, A =0, omega =1,f = 0, N=400):
    a_t = max(A*np.cos(omega*t) +a,0-5)
    return a_t * S * I / N - b * I
#%%
def k_updates(RHS, X, RHS_kargs, stepsize):
    weights = np.array([1,2,2,1])
    temp = np.copy(X)
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
            k[j,i] = RHS[j](*temp, **RHS_kargs[j])

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
        new_X = X_t[1:, i - 1] + stepsize/6*np.sum(k, axis=1)
        #set 0 as lower limit
        X_t[1:, i] = np.where(new_X >0, new_X, 0)
 
    return pd.DataFrame(X_t.T, columns=np.append(["time"], var_names)) 
#%%
def save_plot(df, par):
    plt.figure(figsize=(10,10))
    for name in ["S", "I", "R"]:
        plt.plot(df["time"], df[name], label=name)
    plt.legend(loc='best', fontsize = 28)
    plt.xlabel("Time in a.u.", fontsize = 32)
    plt.ylabel("Number", fontsize = 32)
    plt.ylim(0,400)
    plt.tick_params(size =24, labelsize=26)
    plt.tight_layout()
    plt.savefig("./Results/exp_runge_%s_%i.pdf"%(par[0], int(par[1])))

def eval_model(par_to_vary, N = 400, X_0= [300, 100], t_max = 15, kargs = [{}, {}]):
    
    for p in par_to_vary[1]:
        kargs_all = [ {par_to_vary[0]: p} for i in range(len(X_0))]
        for i in range(len(X_0)):
            kargs_all[i].update(kargs[i])
        df = runge_kutta(0.01, t_max, X_0, [rhs_S, rhs_I], kargs_all)
        df["R"] = N- df["S"] - df["I"]
        save_plot(df, (par_to_vary[0], p))
"""
eval_model(('b', [1,2,3,4]))
eval_model(('A', [0.5, 1,2 ,4]), t_max= 25)
eval_model(('omega', [0.2, 0.7, 2, 4 ]), t_max= 25, kargs=[{'A':1}, {'A':1}])
eval_model(('f', [1, 5, 10, 20]))
"""
df = runge_kutta(0.01, 30, [300, 100], [rhs_S, rhs_I],
                [{'A':3, 'f':30, 'a':1, 'c':1.5}, {'A':3,'a':1, 'c':1.5}])
df["R"] = 400- df["S"] - df["I"]
save_plot(df, ('influenza_w_vac', 1))


# %%
