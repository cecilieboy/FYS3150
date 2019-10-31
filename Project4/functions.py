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

def upper_neighbour_index(i,j, L):
    """
    returns the upper neighbour indices,
    which are appropriate for periodic boundary conditions
    """
    if i == L - 1:
        pi = 0
    else:
        pi = i + 1
    if j == L - 1:
        pj = 0
    else:
        pj = j + 1
    return pi, pj

def E_neighbourhood(i,j, state):
    shape = state.shape
    temp = state[i-1,j] +  state[i, j-1]
    pi, pj = upper_neighbour_index(i,j, shape[0])
    temp += state[i, pj] + state [pi, j]
    return state[i,j]*temp

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
            sum_of_neighbours += E_neighbourhood(i,j, state)
            j += 2
        i += 1
    return -J*sum_of_neighbours

def M(state,factor =1):
    """
    calculates the Magnetisation of one state for a given prefactor
    """
    spins = np.where(state == 1, 1, -1)
    return np.abs(factor*np.sum(spins))

#%%
def lattice(T,cutoff = 1000, L =10):
    t = 0
    #optimal Energy calculation method
    if L%2==0:
        energy_of_sate = E_star
    else:
        energy_of_sate = E

    init_lattice = np.random.randint(2,size=(L,L))
    av_lattice = np.copy(init_lattice)
    E_current = energy_of_sate(init_lattice)
    M_current = M(init_lattice)
    Energies = [E_current]
    Magnetz = [M_current] 
    

    
    while t < cutoff: 

        position_i = np.random.randint(L)
        position_j = np.random.randint(L)
        new_lattice = np.copy(init_lattice)

        new_lattice[position_i,position_j] = 1 - init_lattice[position_i,position_j]

        diff_E = E_neighbourhood(position_i,position_j,new_lattice) - E_neighbourhood(position_i,position_j,init_lattice)

        rnd_p = random.uniform(0,1)
        if np.exp(-1/T * diff_E) > rnd_p:
            M_current += new_lattice[position_i, position_j] - init_lattice[position_i, position_j]
            init_lattice = np.copy(new_lattice)
            E_current  += diff_E


        t += 1
        av_lattice += init_lattice
        Energies.append(E_current)
        Magnetz.append(M_current)

    plt.subplot(121)
    plt.title('Energy')
    plt.plot(Energies)    
    plt.subplot(122)
    plt.title('Magnet')
    plt.plot(Magnetz)
    plt.show()

    c=plt.matshow(av_lattice/(cutoff+1))
    plt.colorbar(c)
    plt.show()

lattice(1,cutoff=10000,L=20)
