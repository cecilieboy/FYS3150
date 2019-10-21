import numpy as np
from mpi4py import MPI
import sys
import random
import time

from Project3 import integrand_radial

def Random(N):

    r1 = np.random.exponential(1,N)
    r2 = np.random.exponential(1,N)
    ct1, ct2 = 2 * np.random.rand(2 *n).reshape((2,n)) -1
    p1, p2 = 2 * np.pi * np.random.rand(2 * n).reshape((2,n)) - 1
    return r1, r2, ct1, ct2, p1, p2


def integral(n):   
    integral = np.sum(integrand_radial(*Random(n)))
    return (integral * (4 * np.pi)**2 / 4**5) 


comm = MPI.COMM_WORLD
Size = comm.Get_size()
rank = comm.Get_rank()

MC_exp = int(sys.argv[1])
MC_samp = int(sys.argv[2])
n = MC_samp // Size

    
final_int = np.zeros(MC_exp)
ints = np.zeros(MC_exp)

start = time.time()


for i in range(MC_exp):
    ints[i] = integral(n)/MC_samp


comm.Reduce(ints, final_int, op = MPI.SUM,root=0)

if rank == 0: 
    print('sum:',np.mean(final_int))
    f = open('paralell.txt','a+')
    f.write(str(Size) + '\t' +str(time.time() - start) + '\t' + str(np.mean(final_int)) +'\t' + str(np.var(final_int)) +  '\t' + str(MC_samp) + '\n')
    f.close()
<<<<<<< HEAD
    
    #print(Time)
=======
>>>>>>> f3439c940a95cef8383c8a75dd4baa1ec575ef35




#%%
