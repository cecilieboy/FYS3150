

#%%
import numpy as np
from scipy.special import legendre
import numpy.polynomial.legendre as leg
from numpy.linalg import inv

#%%

def returnroots (N):
    input = np.zeros(N+1)
    for i in range(N+1):
        input[i] = 0 
        input[N] = 1
    roots = leg.legroots(input)
    return roots 
#%%

def L_function (N):
    L = np.zeros((N, N))
    roots = returnroots(N)
    for i in range(N):
        input = np.zeros(N)
        for j in range(N):
            input[i] = 1
            L[j,i] = leg.legval(roots[j], input)  
    L_invers = inv(L)
    weights = 2 * L_invers[0,:]   
    return weights    

#%%
def function(N):
    x = returnroots(N)
    f = x**2
    return f
#%%
def integral(N):
    sum = 0
    f = function(N)
    weights = L_function(N)
    for i in range(N):
        value = weights[i] * f[i]
        sum += value 
    return sum
