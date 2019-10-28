#%%
import numpy as np
import random
from matplotlib import pyplot as plt
#%%
"""
defining all functions for Ising model,
expected input is the array of spin states state
"""
import numpy as np 

def E(state, J = 1, b=0):
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
    return -J*sum_of_neighbours - b* np.where(state==1, 1, -1).sum()

def M(state,factor =1):
    """
    calculates the Magnetisation of one state for a given prefactor
    """
    spins = np.where(state == 1, 1, -1)
    return factor*np.sum(spins)

#%%
def lattice(T,cutoff = 1000, L =20):
    t = 0
    
    init_lattice = np.random.randint(2,size=(L,L))
    av_lattice = init_lattice
    E_init = E(init_lattice)
    Energies = [E_init]
    Magnetz = [M(init_lattice)] 
    while t < cutoff: 
        #print('====')
        #print('original Zustand:\n',init_lattice)
        #print('original Energie:',E_init)
        new_lattice = np.random.randint(2,size=(L,L))
        
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

        t += 1
        av_lattice += init_lattice
        Energies.append(E_init)
        Magnetz.append(M(init_lattice))

    plt.figure()
    plt.title('Energy')
    plt.plot(Energies)    
    plt.figure()
    plt.title('Magnet')
    plt.plot(Magnetz)
    plt.figure()
    c=plt.matshow(av_lattice/(cutoff+1))
    plt.colorbar(c)

lattice(1,cutoff=10000, L=20)
#%%
print(init_latice)



# %%

def mean_Energy(T, J):
    """
    claculates the mean energy for a given 
    temperature T in units J/kb
    pass