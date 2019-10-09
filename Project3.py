

#%%
import numpy as np
from scipy.special import legendre
from numpy.polynomial import legendre, laguerre
from scipy.special import roots_genlaguerre, eval_genlaguerre
from numpy.linalg import inv

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

def GaussLaguerre(N):
    """
    function to return roots and weights for Gauss-Laguerre integration
    with N mesh points
    """
    #initialize coeffitents for Laguerre series
    coeff = np.zeros(N+1)
    #find root of N-th Laguerre polynomial
    coeff[N] = 1
    roots = laguerre.lagroots(coeff)
    #reset coeff.; initalize L matrix
    L = np.zeros((N,N))
    coeff = np.zeros(N)
    for i in range(N):
        #fill j-th Laguerre poly with the i-th root.
        for j in range(N):
            coeff[j] = 1
            L[i, j] = laguerre.lagval(roots[i], coeff)
            coeff[j] = 0
    L_inv = inv(L)
    return   2* L_inv[0,:], roots


#%%
def function(N):
    x = returnroots(N)
    f = x**2
    return f
#%%
def integral(N):
    sum = 0
    f = function(N)
    weights = L_function(N)
    for i in range(N):
        value = weights[i] * f[i]
        sum += value 
    return sum

def integrand_radial(r1, r2, t1, t2, p1,p2, gen_lag = False):
    """
    returns the integrand in radial coordinates
    (r_i might be scaled for proper results)
    gen_lag-flag can be used wehn generalized Laguerre polynomials are used
    treatment of singual values is: usinig principle value + same curvature around r1 = r2
                                    --> return 0
    """
    if gen_lag:
        num = 1
    else:
        num = r1**2 * r2**2
    cos_b = np.cos(t1) * np.cos(t2) +np.sin(t1)*np.sin(t2) * np.cos(p1-p2)
    r12 = np.sqrt( r1**2 + r2**2 - 2 * r1*r2*cos_b )
    if r12 != 0:
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
    weights_p, p = weights_t, np.pi*(t + 1)

    integral  = 0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    for m in range(N):
                        for n in range(N):
                            integral += (weights_r[i]* weights_r[j] * weights_t[k] * weights_t[l]* weights_p[m] * weights_p[n] * 
                                            integrand_radial(x[i], x[j], t[k], t[l], p[m], p[n]))
    return transformation * integral

print(radial_integration(5), 5*np.pi**2/16)

#%%
