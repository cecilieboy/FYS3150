 

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt 

def rotation(A,l,k,v): 
# A er matricen der skal diagonaliseres
# l, k er positionerne for rotationsledende og positionerne i A der bliver minimeret
# v er en mastrix bestående af basisvektorene
    
    B = np.copy(A)    
    w = np.copy(v)
    tau = (A[l,l] - A[k,k]) / (2 * A[k,l])
#tau er defineret sådan så vi får en andengradsligning

    tplus  = 1.0/(tau + np.sqrt(1 + tau**2))
    tminus = 1.0/(tau - np.sqrt(1 + tau**2))
#Udfra betingelsen at b_kl = 0, får vi en andengradsligning med disse løsninger for t        

    if abs(tplus) < abs(tminus):
        t = tplus
    else:
        t = tminus

        
    c = 1.0 / (np.sqrt(1 + t**2))
    s = t * c
    for j in range(len(v)):

        w[j,k] =   c * v[j,k] - s * v[j,l]
        w[j,l] =   s * v[j,k] + c * v[j,l]
        # Her sker selve rotationen

        if j != k and j != l:
            B[j,k] = B[k,j] = A[j,k] * c - A[j,l] * s
            B[j,l] = B[l,j] = A[j,l] * c + A[j,k] * s

#Vi fandt ud af ved testing at, hvis vi bruger rotatione matricen, sker der kun noget med k'te og l'te række i vær emkel søjle

    B[k,k] = A[k,k] * c**2 - 2 * A[k,l] * c * s + A[l,l] * s**2
    B[l,l] = A[l,l] * c**2 + 2 * A[k,l] * c * s + A[k,k] * s**2
    B[k,l] = B[l,k] = 0


    return B,w

def diag_A(A, tol = 10**(-10), max_count =1000): 
    n = len(A)    
    v = np.eye(n,n)
    count = 0
    while count < max_count:    
        highestmatrixentry = -1  
        for i in range(len(A)-1):
            temp = np.argmax(np.abs(A[i, i+1:]))
            temp += i+1
            if np.abs(A[i,temp]) > highestmatrixentry:
                highestmatrixentry = np.abs(A[i,temp])
                l, k= i, temp
        if highestmatrixentry < tol:
            break
    #Finder positionen af den højeste matrice indgang, med et loop som går igennem alle indgangene og opdaterer den højeste indgangs vœrdi/position
        A,v = rotation(A,l,k,v)
        count += 1  
    #sort eigenvalues
    lam = np.diag(A)
    ind_sort = np.argsort(lam)
    vT = v.T
    lam_sort = np.zeros(lam.shape)
    vT_sort = np.zeros(vT.shape) 
    for i in range(len(lam)):
        ind = ind_sort[i]
        lam_sort[i] = lam[ind]
        vT_sort[i] = vT[ind]
    #returns sorted eigen values and sorted eigenvectors column-wise
    return lam_sort, vT_sort, count

def discretize_HO(integrationpoints,rho_min = 0, rho_max = 10):
    h = (rho_max -rho_min) / integrationpoints
    rho = np.array([rho_min +i *h for i in range(1,integrationpoints)])
    d = 2.0/h**2 * np.eye(integrationpoints -1) + np.diag(rho**2)
    e_upper = - 1.0/h**2 * np.eye(integrationpoints -1 , k=1)
    e_lower = - 1.0/h**2 * np.eye(integrationpoints -1 , k=-1)
    return d + e_lower + e_upper, rho
    
def charac_poly(a,b):
    n = len(a) + 1
    f, q = np.zeros((2,n))
    f[0] = 1
    for i in range(0, n -1):
        f[i + 1] = a[i] * f[i] 
        if i > 0:
            f[i+1] -= b[i-1]**2 * f[i - 1] 
        if abs(f[i]) < 10**(-15):
            q[i + 1] = a[i] - abs(b[i-1])*10**(-15)
        else:
            q[i +1] = f[i + 1]/f[i]
    return f[-1]

def stable_ver(a,b):
    n = len(a) 
    q = np.zeros(n + 1)
    q[0] = 1
    q[1] = a[0]
    for i in range(1, n):
        q[i+1] = a[i] - b[i-1]**2 / q[i - 1]
    return q[-1]

    
    
def bisect(a, b, low, up, max_count, tolerance):
    
    count = 0    
    while np.abs(up - low) > tolerance and count < max_count:
        c = 0.5*(low + up)
        fa = charac_poly(a - low, b)
        fb = charac_poly(a - up, b)
        fc = charac_poly(a - c, b)   
        if fa*fb > 0:
            raise ArithmeticError("Fail at ", count)   
            
        if fa*fc < 0:
            up = c 
        if fc*fb < 0:
            low = c
        count += 1 
    return up
      
def find_eigen(a,b, max_iter = 10**5, max_eigen = 3, tol = 10**-5):
    """
    UNSTABLE VERSION
    calculates the first max_eigen values of a tridiagonal, symetric matrix
    with a the diagonal and b the non diagonal entries
    a is expected to be sorted in ascending order
    max_iter is used for the maximum number of iterations
    returns the first max_eigen eigenvalues
    """
    n = len(a)
    lam = np.zeros(max_eigen + 1)
    #setup bounderis with Gash Goren
    bounds_low = np.array([a[i] - np.abs(b[i]) - np.abs(b[i+1]) for i in range(max_eigen + 1)])
    bounds_up = np.array([a[i] + np.abs(b[i]) + np.abs(b[i+1]) for i in range(max_eigen + 1)])
    #first iteration lam = a[0]
    bounds_up[0] = a[0]
    bounds_low[1] = a[0]
    for i in range(2,n):
        #print(i)
        blt = np.copy(bounds_low)
        but = np.copy(bounds_up)
        lt = np.copy(lam)
        for j in range(max_eigen + 1):
            if j < i:
                lam[j] = bisect(a[:i], b[:i - 1 ], blt[j], but[j], max_count=max_iter, tolerance=tol)
                #update boundaries with interlacing theorem
                bounds_up[j] = lam[j]
                if j + 1 < max_eigen + 1:
                    bounds_low[j + 1] = lam[j]

        if np.max(np.abs(lam-lt)) < tol:
            print('Converged')
            return lam[:-1]

    return lam[:-1]
    
"""
BEGIN OF STABel Implementation
def find_eigen(a,b, max_iter = 10**5, max_eigen = 3, tol = 10**-5):
    #
    #calculates the first max_eigen values of a tridiagonal, symetric matrix
    #with a the diagonal and b the non diagonal entries
    #a is expected to be sorted in ascending order
    #max_iter is used for the maximum number of iterations
    #returns the first max_eigen eigenvalues
    
    n = len(a)
    lam = np.zeros(max_eigen)
    #setup bounderis with Gash Goren
    bounds_low = a[0] - np.abs(b[0]) 
    bounds_up = a[0] + np.abs(b[0])
    for j in range(max_eigen):
        if j > 0 :
            c =a[j] - np.abs(b[j]) - abs(b[j+1])
            if c>bounds_up:
                bounds_low = c
            else: 
                bounds_low = bounds_up
            bounds_up = a[j] + abs(b[j]) + abs(b[j+1])  
            a = 0  
            q = 1        
            for i in range(1,n):
                q = a[i] -
                #update boundaries with interlacing theorem
                bounds_up = lam[j]
                if j + 1 < max_eigen + 1:
                    bounds_low[j + 1] = lam[j]

                if bounds_up- bounds_low < tol:
                    print('Converged')
                    return lam

    return lam
"""

def main():
    ## benchmark diag_A n = 5, 10, 15, 20 , 25, tol e-4,e-8, e-12
    N = [5* i for i in range(1,17)]
    tol = [10**(-1), 10**(-8), 10**(-16)]
    count = np.zeros((len(tol), len(N)))
    
    plt.figure(figsize=(10,10))
    for i, t in enumerate(tol):
        for j, n in enumerate(N):
            A = -2* np.eye(n) + np.eye(n,k=1) + np.eye(n,k=-1) 
            _,  _, count[i,j] = diag_A(A, tol= t, max_count = 10**9)
        plt.plot(N,count[i], label ='tol = %f' % t)
    plt.legend(loc ='best', fontsize = 24)
    plt.xlabel(r"Dimensions", fontsize = 24)
    plt.ylabel(r"# Rotation calls", fontsize = 24)
    plt.savefig('benchmark.pdf') 

   
    
    integration_points = [i*10 for i in np.arange(2 ,26, 4)]
    rho_max = [2, 5, 7.5, 10, 12.5]
    max_rel_err = np.zeros((len(rho_max),len(integration_points)))
    print("rel. err")
    for i, r in enumerate(rho_max):
        print("solution")
        A,rho = discretize_HO(150, rho_max=r)
        lam, u, _ = diag_A(A,tol=10**(-10), max_count=10**9)
        plt.figure(figsize=(10,10))
        for k in range(3):
            plt.plot(rho , u[k]**2, label = "$\lambda$ = %.5f" %lam[k] )
        plt.legend(loc = 'best', fontsize = 24)
        plt.xlabel(r"$\rho$", fontsize = 24)
        plt.ylabel(r"|u($\rho$)|$^2$", fontsize = 24)
        plt.savefig('solution_%i.pdf'%r)
        plt.close('all')
        for j, N in enumerate(integration_points):
            print(r,N)
            lam_theo = [3 + i*4 for i in range(N-1)]
            A, rho = discretize_HO(N, rho_max=r)
            lam, u, c = diag_A(A, tol=10**(-5), max_count= 10**6)
            max_rel_err[i,j] = np.max(np.abs(lam - lam_theo)/lam_theo)
            print("count: %i" %c)

    plt.figure(figsize=(10,10))
    for i in range(len(rho_max)):
        plt.plot(integration_points, max_rel_err[i], label=r"$\rho_{max}$ = " + str(rho_max[i]))
    plt.legend(loc = 'best', fontsize = 24)
    plt.xlabel(r"N", fontsize = 24)
    plt.ylabel(r"max$\frac{|\lambda-\lambda_{theo}|}{\lambda_{theo}}$", fontsize = 24)
    plt.savefig('error.pdf')
    
def test_eigen():
    A, rho = discretize_HO(100,0,10)
    a = np.diag(A)
    b = np.diag(A,k=1)
    lam = find_eigen(a,b, tol=10**(-8), max_iter= 10**6)
    print(lam)

if __name__ == '__main__':
    test_eigen()
    #main()

#%%
