#%%
from mpi4py import MPI
import sys
import numpy as np
from time import perf_counter 
from functions import lattice

# %%

start = perf_counter()
comm = MPI.COMM_WORLD
worldSize = comm.Get_size()
rank = comm.Get_rank()
TaskMaster = 0

f = open('paralell4.txt','a+')

T = np.arange(0.1,5,0.05)
print(T)

L = int(sys.argv[1])
shape = T.shape

a = shape[0] // worldSize
N = 1

if rank == worldSize - 1: 

    for r in T[rank * a:]:
        for i in range(N):
            E,M = lattice(r)
            f.write(str(r) + '\t' + str(np.mean(E)) +'\t' + str(np.mean(M)) +  '\n')

else: 

    for r in T[rank * a : (rank + 1) * a]:
        for i in range(N):
            E,M = lattice(r)
            f.write(str(r) + '\t' + str(np.mean(E)) +'\t' + str(np.mean(M)) +  '\n')


    
f.close()
        


# %%
