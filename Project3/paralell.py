from mpi4py import MPI
import sys
import numpy as np
from Project3 import GaussLaguerre, L_function, integrand_radial
from time import perf_counter

start = perf_counter()
comm = MPI.COMM_WORLD
worldSize = comm.Get_size()
rank = comm.Get_rank()
TaskMaster = 0

N = int( sys.argv[1])
transformation = np.pi ** 2 / 4**2
weights_r, r = GaussLaguerre(N)
x  = r / 4
weights_t, t = L_function(N)
weights_p, p = weights_t,  np.pi*( t + 1)
TaskMaster = 0

index = np.arange(0,N)

ix, iy ,iz = np.meshgrid(index, index, index)
ix = ix.flatten()
iy = iy.flatten()
iz = iz.flatten()


    
size = len(x) // worldSize
integral = np.zeros(worldSize, 'd')
print(rank, size)
comm.Barrier()
for i in ix[rank*size : (rank + 1)*size]:
    for j in iy[rank*size : (rank + 1)*size]:
        for k in iz[rank*size : (rank + 1)*size]:
            for l in range(N):
                for m in range(N):
                    for n in range(N):
                        
                        integral[rank] += (weights_r[i]* weights_r[j] * weights_t[k] * weights_t[l]* weights_p[m] * weights_p[n] * 
                                            integrand_radial(x[i], x[j], t[k], t[l], p[m], p[n], gen_lag=False))
    
if rank == TaskMaster:
    #calc rest
    for i in ix[worldSize*size : ]:
        for j in iy[worldSize*size : ]:
            for k in iz[worldSize*size :]:
                for l in range(N):
                    for m in range(N):
                        for n in range(N):
                            integral[TaskMaster] += (weights_r[i]* weights_r[j] * weights_t[k] * weights_t[l]* weights_p[m] * weights_p[n] * 
                                                integrand_radial(x[i], x[j], t[k], t[l], p[m], p[n], gen_lag=False))
comm.Barrier()
time = perf_counter() - start
if rank == TaskMaster:
    print(time, transformation * integral)
