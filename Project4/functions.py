"""
defining all functions for Ising model,
expected input is the array of spin states state
"""
import numpy as np 

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

def M(state,factor =1):
    """
    calculates the Magnetisation of one state for a given prefactor
    """
    spins = np.where(state == 1, state, -1)
    return factor*np.sum(spins)

def mean_Energy(T, J):
    """
    claculates the mean energy for a given 
    temperature T in units J/kb
