

#%%
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import legendre
from numpy.polynomial import legendre, laguerre
from scipy.special import roots_genlaguerre, eval_genlaguerre
from numpy.linalg import inv
from time import perf_counter
import pandas as pd
#%%

def returnroots (N):
    coeff = np.zeros(N+1)
    coeff[N] = 1
    roots = legendre.legroots(coeff )
    return roots 
#%%

def L_function (N):
    L = np.zeros((N, N))
    roots = returnroots(N)
    for i in range(N):
        input = np.zeros(N)
        for j in range(N):
            input[i] = 1
            L[j,i] = legendre.legval(roots[j], input)  
    L_invers = inv(L)
    weights = 2 * L_invers[0,:]   
    return  weights  , roots  

#%%

def GaussLaguerre(N):
    """
    function to return roots and weights for Gauss-Laguerre integration
    with N mesh points
    no multiplication for weights for Laguerre!
    """
    #initialize coeffitents for Laguerre series
    coeff = np.zeros(N+1)
    #find root of N-th Laguerre polynomial
    coeff[N] = 1
    roots = laguerre.lagroots(coeff)
    #reset coeff.; initalize L matrix
    L = np.zeros((N,N), dtype = np.float64)
    coeff = np.zeros(N)
    for i in range(N):
        #fill j-th Laguerre poly with the i-th root.
        for j in range(N):
            coeff[j] = 1
            L[i, j] = laguerre.lagval(roots[i], coeff)
            coeff[j] = 0
    L_inv = inv(L)
    return   L_inv[0,:], roots


#%%
#defines function to be integrated
def function(N,x1,y1,z1,x2,y2,z2):
    if np.sqrt((x1 - x2)**2 + (y1 - y2)**2+ (z1 - z2)**2 ) < 0.000001:
        f = 0

    else:
        f = np.exp(-4 * (np.sqrt(x1**2 + y1**2 + z1**2) + np.sqrt(x2**2 + y2**2 + z2**2))) * 1 / (np.sqrt((x1 - x2)**2 + (y1 - y2)**2+ (z1 - z2)**2 ))

    return f
#%%
#function executing the integral
#summation over all possible combinations of chosen rootpoints for each variable
#And also sums over all possibla factors of weights
def integral(N):
    sum = 0

    x1 = y1 = z1 = x2 = y2 = z2 = returnroots(N)

    weights = L_function(N)
    
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    for m in range(N):
                        for n in range(N):
                            sum += weights[i] * weights[j] * weights[k] * weights[l] * weights[m] * weights[n] * function(N,x1[i],y1[j],z1[k],x2[l],y2[m],z2[n])
    return sum

integral(2)


#%%
#Function to be integrated with change of variables
def changedfunction(N,x1,y1,z1,x2,y2,z2,a,b): 

    r1 = (b - a) / 2 * x1 + (b + a) / 2
    s1 = (b - a) / 2 * y1 + (b + a) / 2
    t1 = (b - a) / 2 * z1 + (b + a) / 2
    
    r2 = (b - a) / 2 * x2 + (b + a) / 2
    s2 = (b - a) / 2 * y2 + (b + a) / 2
    t2 = (b - a) / 2 * z2 + (b + a) / 2

    if np.sqrt((r1 - r2)**2 + (s1 - s2)**2+ (t1 - t2)**2 ) < 0.00001:
        f = 0
    else:
        f = np.exp(-4 * np.sqrt(r1**2 + s1**2 + t1**2) + np.sqrt(r2**2 + s2**2 + t2**2)) * 1 / (np.sqrt((r1 - r2)**2 + (s1 - s2)**2+ (t1 - t2)**2 ))

    return f

#%%
#Integration function with borders a and b 
#Notice: a is lower border!!
def newintegral(N,a,b):
    sum = 0

    x1 = y1 = z1 = x2 = y2 = z2 = returnroots(N)

    weights = L_function(N)
    
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    for m in range(N):
                        for n in range(N):

                            value = weights[i] * weights[j] * weights[k] * weights[l] * weights[m] * weights[n] * changedfunction(N,x1[i],y1[j],z1[k],x2[l],y2[m],z2[n],a,b)

                            sum += value
    return sum * (b-a) / 2





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
    if r12 > cutoff:
        return num / r12
    else:
        return 0

def radial_integration(N):
    """
    integrate in radial coordinates: 0 to inf for r^2dr, -1 to 1 for dcos(theta), 0 to 2pi for dphi
    use laguerre polynomials for radial part
    -> def x = 4r <-> r = x/4; dr = dx/4
    use legandre polynomials for other; transformation of integralbounds 
    s -> pi*(p + 1); ds = pi dp
    """
    transformation = np.pi ** 2 / 4**2
    weights_r, r = GaussLaguerre(N)
    x  = r / 4
    weights_t, t = L_function(N)
    weights_p, p = weights_t,  np.pi*( t + 1)

    integral  = 0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):  
                    for m in range(N):
                        for n in range(N):
                            integral += (weights_r[i]* weights_r[j] * weights_t[k] * weights_t[l]* weights_p[m] * weights_p[n] * 
                                            integrand_radial(x[i], x[j], t[k], t[l], p[m], p[n], gen_lag=False))
    return transformation * integral


#%%

def norm(x, y, z):
    """
    vector norm
    """
    return np.sqrt(x**2 + y**2 + z**2)

def integrand(x1,x2,y1,y2,z1,z2, cutoff = 10**-5):
    """
    integrand function compatible with numpy arrays
    """
    r1 = norm(x1, y1, z1)
    r2 = norm(x2, y2, z2)
    r12 = norm(x1-x2, y1-y2, z1-z2)
    if np.min(r12) > cutoff:
        return np.exp(-4 * (r1 + r2)) / r12
    else:
        return 0
    
def brutforce_mc(samples, cutoff, cycles):
    """
    calculates MC integration with uniform sampling on the 
    interval [-cutoff, cutoff]^6 -> V (2*cutoff)^6
    The integral I is approximated cycels-times with 
    I = V <f>
    and the variance of the integral is
    Var(I) = V Var(f)/ cycels
    """
    result = np.zeros(cycles)
    var_theo = np.zeros(cycles)
    V =2**6* cutoff**6
    for i  in range(cycles):
        points = -cutoff + 2*cutoff * np.random.rand(samples * 6)
        points = points.reshape(6,samples)
        I = integrand(*points)
        result[i] = V*np.mean(I)
        var_theo[i] = V**2 * np.var(I)/samples 
    return np.mean(result), np.var(result), np.max(var_theo)

theory = 5*np.pi**2/16**2

N = np.arange(2,19, 2)
lenN = len(N)
toi_laguerre = pd.DataFrame(data=np.zeros(( lenN, 4)),columns=["N", "time", "I", "rel_err"], dtype=np.float64)
toi_laguerre["N"].iloc[:] = N

for i, n in enumerate(N):
    print(n)
    start = perf_counter()
    toi_laguerre["I"].iloc[i] = radial_integration(n)
    toi_laguerre["time"].iloc[i] = perf_counter() -start
    toi_laguerre["rel_err"].iloc[i] = np.abs( toi_laguerre["I"].iloc[i] - theory)/theory
toi_laguerre.to_csv("laguerre.csv")

samples = [10**i for i in range(3,6)]
cutoff = np.linspace(1,15,10)
cycles = [1, 10, 100]
index = 0
toi_mc = pd.DataFrame(data=np.zeros(( len(samples)*len(cutoff)*len(cycles), 8)),
            columns=["samples","cutoff", "cycles", "time", "I","var_mc","var_t", "rel_err"], dtype=np.float64)

for s in samples:
    for cu in cutoff:
        for cy in cycles:
            print (s, cu, cy)
            toi_mc["samples"].iloc[index] = s
            toi_mc["cutoff"].iloc[index] = cu
            toi_mc["cycles"].iloc[index] = cy

            start = perf_counter()
            I, var_mc , var_t = brutforce_mc(s, cu, cy)
            time = perf_counter() - start

            toi_mc["I"].iloc[index] = I
            toi_mc["var_mc"] = var_mc
            toi_mc["var_t"] = var_t
            toi_mc["rel_err"] = np.abs(I-theory)/theory
            index += 1
toi_mc.to_csv("brutforce_mc.csv")
#%%
