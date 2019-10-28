

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
#%%
 
def returnroots (N):
    """
    function returns roots for Gauss-Legendre integration
    """
    coeff = np.zeros(N+1)
    coeff[N] = 1
    roots = legendre.legroots(coeff )
    return roots 
#%%

def L_function (N):
    """
    function returns roots and weigths for Gauss-Legendre integration
    """
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
def function(N,x1,y1,z1,x2,y2,z2):
    """
    returns integrand in cartesian coordinates. For singular points: returns 0
    """
    if np.sqrt((x1 - x2)**2 + (y1 - y2)**2+ (z1 - z2)**2 ) == 0:
        f = 0
    else:
        f = np.exp(-4 * (np.sqrt(x1**2 + y1**2 + z1**2) + np.sqrt(x2**2 + y2**2 + z2**2))) * 1 / (np.sqrt((x1 - x2)**2 + (y1 - y2)**2+ (z1 - z2)**2 ))
    return f
#%%
#Integration function with borders a and b 
#Notice: a is lower border!!
def newintegral(N,a,b):
    """
    integrates in cartesian coordinates. arguments of function are the roots
    multiplied by (b - a)/2 where a is the lower integration border and b the upper. 
    These borders are chosen to be symmetric (a =-b)
    """
    summe = 0
    x1 = y1 = z1 = x2 = y2 = z2 = returnroots(N) * (b - a) / 2
    weights,_ = L_function(N)
    start = time.time()
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    for m in range(N):
                        for n in range(N):
                            value = value = weights[i] * weights[j] * weights[k] * weights[l] * weights[m] * weights[n] * function(N,x1[i],y1[j],z1[k],x2[l],y2[m],z2[n])
                            summe += value
    Time =  time.time() - start
    return summe * ((b-a) / 2)**6,Time

<<<<<<< HEAD
#%%
######################## Evaluation of integrals ######################
#N = np.arange(2,19,2)
#Integrals = np.zeros(len(N))
#theo_int = 5 * np.pi**2 / 16**2
#rel_error = np.zeros(len(N))
#Time = np.zeros(len(N))
#for i in trange(len(N)):
#    Integrals[i],Time[i] = newintegral(N[i],-5,5)
#    rel_error[i] = np.abs(Integrals[i] - theo_int) / theo_int
#print(Integrals)
#print(Time)
#print(rel_error)
#Data = {'N':N,'I':Integrals,'time': Time,  'rel_err': rel_error}
#data_legendre = pd.DataFrame(Data)   
#data_legendre.to_csv('Results/legendre_5.csv')
=======
>>>>>>> f3439c940a95cef8383c8a75dd4baa1ec575ef35

#%%


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
    V =2**6* cutoff**6
    for i in range(cycles):
        points = -cutoff + 2*cutoff * np.random.rand(samples * 6)
        points = points.reshape(6,samples)
        result[i] = V*np.mean(integrand(*points))
    return np.mean(result), np.var(result) 




#%%

def exp_mc(N):
    """
    integrates radial integrand by using importance sampling.
    defines the random variables scaled with the right borders for r, ct and p.
    """
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

    integral = 0 
    for i in range(N):
        integral += integrand_radial(r1[i],r2[i],ct1[i],ct2[i],p1[i],p2[i])  
    return (integral * (4 * np.pi)**2 / 4**5) / N
                        
#%%

def main():
<<<<<<< HEAD
    """
    carries out function:exp_mc for different number of Samples
    and for a given number of calls. return mean value and variance. 
    """
    nmb_samples = [10**2,10**3,10**4,10**5,10**6] 
    nmb_calls = 10
=======
>>>>>>> f3439c940a95cef8383c8a75dd4baa1ec575ef35
    theo_int = 5 * np.pi**2 / 16**2

    #################################################################################################################
    #Legendre Integration
    #################################################################################################################
    N = np.arange(2,19,1)
    Integrals = np.zeros(len(N))
    rel_error = np.zeros(len(N))
    Time = np.zeros(len(N))
    
    for i in trange(len(N)):
        Integrals[i],Time[i] = newintegral(N[i],-3,3)
        rel_error[i] = np.abs(Integrals[i] - theo_int) / theo_int
    
<<<<<<< HEAD
    for i in trange(len(nmb_samples)):
        start = time.time()
        for j in range(nmb_calls):
            results[i,j] = exp_mc(nmb_samples[i])
        Time[i] =  time.time() - start
    
    mean_int = np.mean(results,axis = 1)
    var_int = (np.var(results,axis = 1))
    rel_err = (np.abs(mean_int - theo_int) / theo_int)
 
    Data = {'mean of integrals':mean_int,'Samples(N)':nmb_samples,'log10_Samples(N)':np.log10(nmb_samples),'Time': Time, 'log10_Time':np.log10(Time), 'Variance': var_int,'log10_Var':np.log10(var_int), 'Relative Error': rel_err,'log10_rel_err':np.log10(rel_err)}
    data_exp_mc = pd.DataFrame(Data)
    data_exp_mc.to_csv('Results/main_results_%i.csv'%nmb_calls)

if __name__ == '__main__':
    main()

 #%%

    theory = 5*np.pi**2/16**2
=======
    Data = {'N':N,'I':Integrals,'time': Time,  'rel_err': rel_error} 
    data_legendre = pd.DataFrame(Data) 
    data_legendre.to_csv('Results/legendre.csv')
>>>>>>> f3439c940a95cef8383c8a75dd4baa1ec575ef35

    #################################################################################################################
    #Laguerre Integration
    #################################################################################################################
    N = np.arange(2,17, 2)
    lenN = len(N)
    toi_laguerre = pd.DataFrame(data=np.zeros(( lenN, 4)),columns=["N", "time", "I", "rel_err"], dtype=np.float64)
    toi_laguerre["N"].iloc[:] = N
    
    for i, n in enumerate(N):
        start = perf_counter()
        toi_laguerre["I"].iloc[i] = radial_integration(n)
        toi_laguerre["time"].iloc[i] = perf_counter() -start
        toi_laguerre["rel_err"].iloc[i] = np.abs( toi_laguerre["I"].iloc[i] - theo_int)/theo_int
    toi_laguerre.to_csv("Results/laguerre.csv")
    
    #################################################################################################################
    #Brute-force MC
    #################################################################################################################
    samples = [10**i for i in range(2,7)]
    cutoff = np.arange(1,11)
    cycles = [ 10, 100, 500]
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
                toi_mc["time"].iloc[index]= perf_counter() - start

                toi_mc["I"].iloc[index] = I
                toi_mc["var_mc"].iloc[index] = var_mc
                toi_mc["var_t"].iloc[index] = var_t
                toi_mc["rel_err"].iloc[index] = np.abs(I-theo_int)/theo_int
                index += 1
    toi_mc.to_csv("Results/brutforce_mc.csv")
    
    #################################################################################################################
    #Importance sampling
    #################################################################################################################
    nmb_samples = [10**2,10**3,10**4,10**5,10**6] 
    nmb_calls = 10
    theo_int = 5 * np.pi**2 / 16**2
    
    results = np.zeros((len(nmb_samples),nmb_calls))
    rel_err = np.zeros(len(nmb_samples))
    Time = np.zeros(len(nmb_samples))
    
    for i in trange(len(nmb_samples)):
        start = time.time()
        for j in range(nmb_calls):
            results[i,j] = exp_mc(nmb_samples[i])
        Time[i] =  time.time() - start
    
    mean_int = np.mean(results,axis = 1)
    var_int = (np.var(results,axis = 1))
    rel_err = (np.abs(mean_int - theo_int) / theo_int)
 
    Data = {'mean of integrals':mean_int,'Samples(N)':nmb_samples,'log10_Samples(N)':np.log10(nmb_samples),'Time': Time, 'Variance': var_int,'log10_Var':np.log10(var_int), 'Relative Error': rel_err,'log10_rel_err':np.log10(rel_err)}
    data_exp_mc = pd.DataFrame(Data)
    data_exp_mc.to_csv('Results/main_results_%i.csv'%nmb_calls)

if __name__ == '__main__':
    main()

#%%
#%%
