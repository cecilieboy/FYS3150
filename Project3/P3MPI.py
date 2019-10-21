#%%

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import legendre
from numpy.polynomial import legendre, laguerre
from scipy.special import roots_genlaguerre, eval_genlaguerre
from numpy.linalg import inv
from time import perf_counter
import pandas as pd
import random
from tqdm import trange
import time

#from Project3 import integrand_radial
from mpi4py import MPI
import sys

#%%
#radial_integration()

comm = MPI.COMM_WORLD
Size = comm.Get_size()
rank = comm.Get_rank()
N = int(sys.argv[1])
n = N // Size
final_int = 0

start = time.time()

def integrand_radial(r1, r2, ct1, ct2, p1,p2, gen_lag = False, cutoff = 10**-5):
    """
    returns the integrand in radial coordinates
    (r_i might be scaled for proper results, its integrated over cos theta)
    gen_lag-flag can be used wehn generalized Laguerre polynomials are used
    treatment of singual values is: usinig principle value + same curvature around r1 = r2
                                    --> return 0
    """
    if gen_lag:
        num = 1
    else:
        num = r1**2 * r2**2
    cos_b = ct1 * ct2 + np.sqrt((1 - ct1**2) * (1 - ct2**2)) * np.cos(p1 - p2)
    r12 = np.sqrt( r1**2 + r2**2 - 2 * r1 * r2 * cos_b )
    return np.where(r12 > cutoff,num / r12, 0)




def Random(N):

    r1 = np.random.exponential(1,N)
    r2 = np.random.exponential(1,N)
    ct1 = np.zeros(N)
    for i in range(N): ct1[i] = (random.random() * 2) - 1
    ct2 = np.zeros(N)
    for i in range(N): ct2[i] = (random.random() * 2) - 1
    p1 = np.zeros(N) 
    for i in range(N): p1[i] = (random.random()) * 2 * np.pi
    p2 = np.zeros(N) 
    for i in range(N): p2[i] = (random.random()) * 2 * np.pi

    return r1, r2, ct1, ct2, p1, p2


def integral(n):
    
    integral = np.sum(integrand_radial(*Random(n)))

    return (integral * (4 * np.pi)**2 / 4**5) 

final_int += integral(n) / N 


print(rank,final_int)


final_int = comm.reduce(final_int,op = MPI.SUM,root=0)

if rank == 0: 
    print()
    print('sum:',final_int)
    f = open('time_array','a+')
    f.write(str(time.time() - start) + '\t' + str(final_int) +  '\t' + str(N) + '\n')
    f.close()
    
    #print(Time)




#%%
