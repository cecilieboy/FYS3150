#%%
import numpy as np
import random
from matplotlib import pyplot as plt
#%%
"""
defining all functions for Ising model,
expected input is the array of spin states state 0, 1
interpret |0> -> -1, |1>-> 1
"""
import numpy as np 

def E(state, J = 1, b=0):
    """
    calcualtes the Energy of a given state, with given J
    uses periodic boundary conditions
    """ 
    state = np.where(state==1, 1, -1)
    sum_of_neighbours = 0
    shape = state.shape
    #looping backwards through array to use numpy [-1] index as last entry
    #no double counting this way!
    for i in reversed(range(shape[0])):
        for j in reversed(range(shape[1])):
            sum_of_neighbours += state[i,j]*state[i, j-1] + state[i,j]*state[i-1,j]
    return -J*sum_of_neighbours - b* np.where(state==1, 1, -1).sum()

def E_star(state, J = 1):
    """
    only even grids, 
    """
    state = np.where(state==1, 1, -1)
    sum_of_neighbours = 0
    shape = state.shape
    i = 0
    while i <shape[0]:
        j = i % 2
        while j < shape[1]:
            temp = state[i-1,j] +  state[i, j-1]
            if i == shape[0] - 1:
                pi = 0
            else:
                pi = i + 1
            if j == shape[1] - 1:
                pj = 0
            else:
                pj = j + 1
            temp += state[i, pj] + state [pi, j]
            sum_of_neighbours += state[i,j]*temp
            j += 2
        i += 1
    return -J*sum_of_neighbours

def M(state,factor =1):
    """
    calculates the Magnetisation of one state for a given prefactor
    """
    spins = np.where(state == 1, 1, -1)
    return factor*np.sum(spins)

from time import perf_counter
arr = np.random.randint(2, size =(20,20))
print(arr)
start = perf_counter()
e = E(arr)
te= perf_counter()- start 
print("Normal: E=%i, time = %.4f"%(e,te))
start = perf_counter()
e = E_star(arr)
te1= perf_counter()- start 
print("Star: E=%i, time = %.4f"%(e,te1))
print("time saving for 1e6 iter: ", (te-te1)*1e6)
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
        M_new = M(new_lattice)
        print(E_new)
        print(M_new)

        rnd_p = random.uniform(0,1)
        print('rndp',rnd_p)
        print('vergleich',np.exp(-T * (E_new-E_init)))
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



