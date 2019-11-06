#%%
from mpi4py import MPI
import sys
import numpy as np
from time import perf_counter 
from functions import lattice
import array
from tqdm import tqdm
# %%

#start = perf_counter()
comm = MPI.COMM_WORLD
worldSize = comm.Get_size() #(Mac: worldsize:2)
rank = comm.Get_rank()
TaskMaster = 0

L = int(sys.argv[1])
cut = 10
#       write only     | create file   | append at end
amode = MPI.MODE_WRONLY|MPI.MODE_CREATE|MPI.MODE_APPEND
f = MPI.File.Open(comm, 'paralell_L%i.txt'%L,amode)

#T = np.arange(2,2.305,0.05)
#print(T)

T = np.append(np.append(np.linspace(0.1, 2, 10 ),   #0.2 stepsize
                        np.linspace(2,2.5,50)),     #0.01 stepsize
                        np.linspace(2.5,4.5,10))    #0.2 stepsize

shape = T.shape

a = shape[0] // worldSize
N = 1



for r in tqdm(T[rank * a : (rank + 1) * a]):
    print('first cycle:',r)
    print('ranktimes_a:',rank * a)
    for i in range(N):
        #       T     E   var E  cv    M    var M  chi
        out = "%f \t %f \t %f \t %f \t %f \t %f \t %f \n" % lattice(r, cutoff = cut, L = L, plot = False)
        buf = bytearray()
        buf.extend(map(ord, out))
        f.Write_shared(buf)

comm.Barrier()

#for r in tqdm(T[shape[0] - (a * worldSize):shape[0]]):
#    print('second cycle:',r)
#    for i in range(N):
#        #       T     E   var E  cv    M    var M  chi
#        out = "%f \t %f \t %f \t %f \t %f \t %f \t %f \n" % lattice(r, cutoff = cut, L = L, plot = False)
#        buf = bytearray()
#        buf.extend(map(ord, out))
#        f.Write_shared(buf)

#comm.Barrier()

if rank < worldSize/2 : 
    for i in tqdm(range(N//2)):
        #       T       E     cv    M     chi
        out = "%f \t %f \t %f \t %f \t %f \n" % lattice(T[-2], cutoff = cut, L = L, plot = False)
        buf = bytearray()
        buf.extend(map(ord, out))
        f.Write_shared(buf)
else:
    for i in tqdm(range(N//2)):
        #       T       E     cv    M     chi
        out = "%f \t %f \t %f \t %f \t %f \n" % lattice(T[-1], cutoff = cut, L = L, plot = False)
        buf = bytearray()
        buf.extend(map(ord, out))
        f.Write_shared(buf)


comm.Barrier()

if rank == TaskMaster:    
    f.Close()
        


# %%