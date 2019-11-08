#%%
import numpy as np
import random
from matplotlib import pyplot as plt
import statistics
from tqdm import tqdm
#%%
"""
defining all functions for Ising model,
expected input is the array of spin states state 0, 1
interpret |0> -> -1, |1>-> 1
"""
import numpy as np 

def E(spins, J = 1, b=0):
    """
    calcualtes the Energy of a given state, with given J
    uses periodic boundary conditions
    """ 
    sum_of_neighbours = 0
    shape = spins.shape 
    #looping backwards through array to use numpy [-1] index as last entry
    #no double counting this way!
    for i in reversed(range(shape[0])):
        for j in reversed(range(shape[1])):
            sum_of_neighbours += spins[i,j]*spins[i, j-1] + spins[i,j]*spins[i-1,j]
    return -J*sum_of_neighbours - b* spins.sum()

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

def E_neighbourhood(i,j, spins, J=1):
    shape = spins.shape
    temp = spins[i-1,j] +  spins[i, j-1]
    pi, pj = upper_neighbour_index(i,j, shape[0])
    temp += spins[i, pj] + spins [pi, j]
    return - J*spins[i,j]*temp

def E_star(spins, J = 1):
    """
    only even grids, 
    """
    sum_of_neighbours = 0
    shape = spins.shape
    i = 0
    while i <shape[0]:
        j = i % 2
        while j < shape[1]:
            sum_of_neighbours += E_neighbourhood(i,j, spins, J)
            j += 2
        i += 1
    return sum_of_neighbours

def M(spins,factor =1):
    """
    calculates the Magnetisation of one state for a given prefactor
    """
    return factor*np.sum(spins)

#%%
def lattice(T,cutoff = 1000, L =2, plot = False):
    
    t = 0
    #optimal Energy calculation method
    if L%2==0:
        energy_of_sate = E_star
    else:
        energy_of_sate = E

    init_lattice = np.random.randint(2,size=(L,L))
    init_lattice = np.where(init_lattice==1, 1, -1)
    
    
    E_current = energy_of_sate(init_lattice)
    M_current = M(init_lattice)

    if plot:
        average = np.copy(init_lattice)
        Energies = [E_current]
        Magnetz = [M_current] 
    else:
        Energies_squared = []
        Magnetz_squared = []

    keys_diff_E = {diff_E:np.exp(-1/T * diff_E) for diff_E in [-8,-4,0,4,8]}
    
    while t < cutoff: 

        position_i = np.random.randint(L)
        position_j = np.random.randint(L)
        new_lattice = np.copy(init_lattice)

        new_lattice[position_i,position_j] = - init_lattice[position_i,position_j]

        diff_E = E_neighbourhood(position_i,position_j,new_lattice) - E_neighbourhood(position_i,position_j,init_lattice)
        rnd_p = random.uniform(0,1)
        if keys_diff_E[diff_E] > rnd_p:
            M_current += 2* new_lattice[position_i, position_j]
            init_lattice = np.copy(new_lattice)
            
            E_current  += diff_E

        t += 1

        if plot:
            Energies.append(E_current)
            Magnetz.append(M_current)
            average += np.copy(init_lattice)
        elif t > cutoff - 1e5:
            Energies.append(E_current)
            Magnetz.append(M_current)
        

    Energies = np.array(Energies) /L**2
    Magnetz = np.abs(Magnetz)   / L**2  
    if plot:
        
        plt.figure(figsize=(10,10))
        plt.title("T/k$_B$J = %.1f; L = %i"%(T,L), fontsize = 24)

        plt.subplot(121)
        plt.plot(Energies) 
        plt.xlabel('MC cycle', fontsize = 24)
        plt.ylabel('E/JL$^2$', fontsize = 24)  
        plt.xticks(fontsize=20, ticks=np.linspace(0, cutoff, 3))
        plt.yticks(fontsize=20) 

        plt.subplot(122)
        plt.plot(Magnetz)
        plt.xlabel('MC cycle', fontsize = 24)
        plt.ylabel('|M|/$\mu L^2$', fontsize = 24)
        plt.xticks(fontsize=20, ticks=np.linspace(0, cutoff, 3))
        plt.yticks(fontsize=20)
        
        plt.tight_layout()
        #plt.savefig("./Results/Random_Walk_L%i_T%i.pdf"%(L, 10*T))

        plt.figure(figsize=(10,10), clear = True)
        plt.title("T/k$_B$J = %.2f; L = %i"%(T,L), fontsize = 26)
        c = plt.imshow(average/cutoff, cmap='coolwarm')
        ax=plt.colorbar(c)
        ax.set_label("S", fontsize = 24)
        plt.xlabel('x a.u.', fontsize = 24)
        plt.ylabel('y a.u.', fontsize = 24)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.savefig("./Results/Average_Lattice_L%i_T%i.pdf"%(L, 10*T))

        return Energies[cutoff//2:], Magnetz[cutoff//2:]



    E_T = np.mean(Energies)
    cv_T = (np.mean(Energies**2)  - E_T**2)/T**2
    M_T = np.mean(Magnetz) #chi gives always right answer with this def of M_T
    chi_T = ( np.mean(Magnetz**2) - M_T**2)/T

    return T, E_T, np.var(Energies), cv_T, M_T, np.var(Magnetz), chi_T

#%%
def plot_lattice(L, cutoff=10**7, temp = [1.0, 2.4, 3.0]):
    """
    function to plot lattice random walks and distribution in stable tail of random walk
    """
    fig = np.array([plt.subplots(1, figsize=(10,10)) for i in range(2)])
    
    for T in tqdm(temp):
        stab_E, stab_M = lattice(T,cutoff=cutoff,L=L,plot = True)

        if T == 1.0:
            b = 5
        else:
            b = 15

        fig[0,1].hist(stab_E, bins = b, weights=np.ones(len(stab_E))/len(stab_E), density=False, alpha = 0.6, label="T/k$_B$J = %.1f"%T)
        fig[1,1].hist(stab_M, bins = b, weights=np.ones(len(stab_M))/len(stab_M), density=False, alpha = 0.6, label="T/k$_B$J = %.1f"%T)
       
    fig[0,1].set_xlabel('E/JL$^2$', fontsize = 24)  
    fig[0,1].set_ylabel('P(E)', fontsize = 24)  
    fig[0,1].tick_params(axis = 'both',labelsize = 20)
    fig[0,1].legend(loc='best', fontsize =22)

    fig[1,1].set_xlabel('|M|/$\mu L^2$', fontsize = 24)  
    fig[1,1].set_ylabel('P(|M|)', fontsize = 24)  
    fig[1,1].tick_params(axis = 'both',labelsize = 20)
    fig[1,1].legend(loc='best', fontsize =22)

    plt.figure(fig[0,0].number)
    plt.tight_layout()
    plt.savefig("Energies_%i.pdf"%L)

    plt.figure(fig[1,0].number)
    plt.tight_layout()
    plt.savefig("Magnet_%i.pdf"%L)

def anal_sol(T,kb = 1):
    '''''
    analytic solutions for L = 2. J is taken to be 1.
    '''''
    x = 8 / (kb * T) 
    part_func = 4 * (np.cosh(x) + 3)
    mean_e = (-8 * np.sinh(x)) / (np.cosh(x) + 3)
    cv = 1 / (8 * T) * x * ((8**2) * np.cosh(x) / (np.cosh(x) + 3) - ((8 * np.sinh(x)) / (np.cosh(x) + 3))**2)
    mean_abs_m = (8 * np.exp(x) + 16) / (4 * (np.cosh(x) + 3))
    chi = 1 / 8 * x * ((8 * (np.exp(x) + 1)) / (np.cosh(x) + 3) - ((2 * (np.exp(x) + 2)) / (np.cosh(x) + 3))**2)
    #chi_meanmsecond = (((2 * (np.exp(x) + 2)) / (np.cosh(x) + 3))**2)
    #chi_meanmsfirst = (8 * (np.exp(x) + 1)) / (np.cosh(x) + 3)


    return mean_e,cv,mean_abs_m,chi


#%%
def comp_AB():

    T = 1
    mean_e, cv, mean_abs_m, chi = anal_sol(T)
    #print('Anal.solutions:',mean_e,cv,mean_abs_m,chi)
    cutoff = np.asarray([10,100,1000,5000,10000])
    #print(len(cutoff))
    rel_err_e = np.zeros(len(cutoff))
    rel_err_m = np.zeros(len(cutoff))
    rel_err_cv = np.zeros(len(cutoff))
    rel_err_chi = np.zeros(len(cutoff))

    for i in range(len(cutoff)):
        #print(cutoff[i])
        en,mn,cvn,chin = repeat_calls(T,cutoff[i],L=2,plot = False,numb_run = 10,)
        #print(en,mn)
        rel_err_e[i] = np.abs((en - mean_e) / (mean_e))
        rel_err_m[i] = np.abs((mn - mean_abs_m) / (mean_abs_m)) 
        rel_err_cv[i] = np.abs((cvn - cv) / cv)
        rel_err_chi[i] = np.abs((chin - chi) / chi) 

    #print(rel_err_e)
    plt.subplot(1,2,1)
    
    plt.plot(cutoff,rel_err_e,'b')
    plt.plot(cutoff,rel_err_m,'r')

    plt.subplot(1,2,2)
    plt.plot(cutoff,rel_err_cv,'g')
    plt.plot(cutoff,rel_err_chi,'y')
    plt.show()



#%%
def repeat_calls(T=1,cutoff=10000,L=2,plot = False,numb_run = 4):
    repetition = 0
    E_run = []
    M_run = []
    cv_run = []
    chi_run = []
    while repetition < numb_run:
        _,E_T,_,cv_T,M_T,_,chi_T = lattice(T,cutoff,L,plot=False)
        #print(cv_T)
        E_run.append(E_T)
        M_run.append(M_T)
        cv_run.append(cv_T)
        chi_run.append(chi_T)
        repetition += 1
    tot_E_T = np.mean(E_run)
    tot_M_T = np.mean(M_run)
    tot_cv_T = np.mean(cv_T)
    tot_chi_T = np.mean(chi_T)

    return tot_E_T,tot_M_T,tot_cv_T,tot_chi_T
 


#%%



if __name__ =='__main__':
    plot_lattice(20, cutoff=10**6)
    comp_AB()
    repeat_calls()


