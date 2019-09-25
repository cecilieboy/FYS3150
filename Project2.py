#!/usr/bin/env python
# coding: utf-8

# In[2]:

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
#Vi skal åbenbart altid vælge den mindre t?
        
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

   
    
    integration_points = [i*10 for i in np.arange(2 ,21, 4)]
    rho_max = [5, 7.5, 10, 12.5]
    max_rel_err = np.zeros((len(rho_max),len(integration_points)))
    print("rel. err")
    for i, r in enumerate(rho_max):
        print("solution")
        A,rho = discretize_HO(150, rho_max=r)
        lam, u, _ = diag_A(A,tol=10**(-10), max_count=10**12)
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
    
main()

#%%
