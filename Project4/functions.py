#%%
import numpy as np
import random
from matplotlib import pyplot as plt
#%%
"""
defining all functions for Ising model,
expected input is the array of spin states state
"""
def E(state, J = 1):
    """
    calcualtes the Energy of a given state, with given J
    uses periodic boundary conditions
    """ 
    key = {0:1, 1:-1}
    sum_of_neighbours = 0
    shape = state.shape
    #looping backwards through array to use numpy [-1] index as last entry
    #no double counting this way!
    for i in reversed(range(shape[0])):
        for j in reversed(range(shape[1])):
            sum_of_neighbours += key[state[i,j]]*key[state[i, j-1]] + key[state[i,j]]*key[state[i-1,j]]
    return -J*sum_of_neighbours

#%%
def lattice(T,cutoff = 1000):
    t = 0
    Energies = []
    Magnetz = [] 
    init_lattice = np.random.randint(2,size=(2,2))
    
    E_init = E(init_lattice)
    
    while t < cutoff: 
        #print('====')
        #print('original Zustand:\n',init_lattice)
        #print('original Energie:',E_init)
        new_lattice = np.random.randint(2,size=(2,2))
        
        #print(new_lattice)
        E_new = E(new_lattice)
        #print(E_new)
        rnd_p = random.uniform(0,1)
        #print('rndp',rnd_p)
        #print('vergleich',np.exp(-T * (E_new-E_init)))
        if np.exp(-1/T * (E_new-E_init)) > rnd_p:
            init_lattice = new_lattice
            E_init = E_new
            #print('change')
        else: 
            pass
        t += 1
        Energies.append(E_init)
        
        #Magnetz.append()

    plt.figure()
    plt.plot(Energies)    
lattice(300)
#%%




# %%
