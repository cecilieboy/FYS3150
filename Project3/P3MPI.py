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

N = int(sys.argv[1])
n = N // Size
final_int = 0

start = time.time()





final_int += integral(n) / N 
final_int = comm.reduce(final_int,op = MPI.SUM,root=0)

if rank == 0: 
    print('sum:',final_int)
    f = open('time_array','a+')
    f.write(str(Size) + '\t' +str(time.time() - start) + '\t' + str(final_int) +  '\t' + str(N) + '\n')
    f.close()




#%%
