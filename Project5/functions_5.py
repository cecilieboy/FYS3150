#%%
import numpy as np 
from matplotlib import pyplot as plt
from tqdm import trange 
import pandas as pd
import random 
#%%
def g(b, S, I, N=400, a=4, c=0.5):
    return c * (N - S - I) - a * S * I / N

def h(b, S, I, N=400, a=4, c=0.5):
    return a * S * I / N - b * I

#def r(b, S, I, N=400, a=4, c=0.5):
#    return b * I - c * R

def a_runge_katta(stepsize, cutoff, b = [1,2,3,4], N=400, a=4, c=0.5):


    for i in trange(len(b)):
        I_0 = 100
        S_0 = 300
        R_0 = 0
    
        S = [S_0]
        I = [I_0]
        R = [R_0]

        t = stepsize

        numb_steps = cutoff / stepsize

        timearray = np.arange(0, cutoff + stepsize, stepsize)
        
        while t <= cutoff: 
            k1_S = g(b[i], S[-1], I[-1]) * stepsize
            k1_I = h(b[i], S[-1], I[-1]) * stepsize

            k2_S = g(b[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I) * stepsize
            k2_I = h(b[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I) * stepsize

            k3_S = g(b[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I) * stepsize
            k3_I = h(b[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I) * stepsize

            k4_S = g(b[i], S[-1] + k3_S, I[-1] + k3_I) * stepsize
            k4_I = h(b[i], S[-1] + k3_S, I[-1] + k3_I) * stepsize

            S_val = S[-1] + 1/6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S)
            I_val = I[-1] + 1/6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I)

            S.append(S_val)
            I.append(I_val)
            R.append(N - S_val - I_val)

            t += stepsize
        

        Data = {'time':timearray, 'S':S, 'I':I, 'R': R}
        df = pd.DataFrame(Data)
        df.to_csv('Results/a_RK_b=%i.csv'%(b[i]))

a_runge_katta(0.5,40)


  
# %%
def S_to_I(b, S, I, R, stepsize, N=400, a=4, c=0.5):
    return a * S * I / N * stepsize

def I_to_R(b, S, I, R, stepsize, N=400, a=4, c=0.5):
    return b * I * stepsize 

def R_to_S(b, S, I, R, stepsize, N=400, a=4, c=0.5):
    return c * R * stepsize

def a_MC(cutoff, b = [1,2,3,4], N=400, a=4, c=0.5):
    
    for i in trange (len(b)):
        stepsize = min([4 / (a * N), 1 / (b[i] * N), 1 / (c * N)])

        numb_steps = int(cutoff / stepsize)
  
        timearray = np.arange(0, cutoff, stepsize)

        S = np.zeros(numb_steps)
        I = np.zeros(numb_steps)
        R = np.zeros(numb_steps)

        S[0] = 300
        I[0] = 100
        R[0] = 0

        for j in range(numb_steps - 1):
            psi_p = random.uniform(0,1)
            pir_p = random.uniform(0,1)
            prs_p = random.uniform(0,1)

            if S_to_I(b[i],S[j],I[j], R[j], stepsize) > psi_p:
                S[j + 1] = S[j] - 1
                I[j + 1] = I[j] + 1
            else:
                S[j + 1] = S[j]
                I[j + 1] = I[j]

            if I_to_R(b[i],S[j],I[j], R[j], stepsize) > pir_p:
                I[j + 1] = I[j + 1] - 1
                R[j + 1] = R[j] + 1
            else:
                I[j + 1] = I[j + 1]
                R[j + 1] = R[j]

            if R_to_S(b[i],S[j],I[j], R[j], stepsize) > prs_p:
                R[j + 1] = R[j + 1] - 1
                S[j + 1] = S[j + 1] + 1
            else:
                R[j + 1] = R[j + 1]
                S[j + 1] = S[j + 1]

        Data = {'time':timearray, 'S':S, 'I':I, 'R': R}
        df = pd.DataFrame(Data)
        df.to_csv('Results/a_MC_b=%i.csv'%(b[i]))
#a_MC(40)

#%%
#def gc(d_I, S, I, R, a=4, b = 1, c=0.5, e = 0, d = 0):
#    return c * R - a * S * I / (S + I + R) - d * S + e * (S + I + R)

#def hc(d_I, S, I, R, a=4, b = 1, c=0.5, e = 4, d = 3):
#    return a * S * I / (S + I + R) - b * I - d * I - d_I * I 

#def rc(d_I, S, I, R, a=4, b = 1, c=0.5, e = 4, d = 3):
#    return b * I - c * R - d * R

#def c_runge_katta(stepsize, cutoff, d_I = [0], a=4, b =1, c=0.5, e = 0, d = 0):


#    for i in trange(len(d_I)):
#        I_0 = 100
#        S_0 = 300
#        R_0 = 0
#        N_0 = 400
    
#        S = [S_0]
#        I = [I_0]
#        R = [R_0]
#        N = [N_0]

#        t = stepsize

#        numb_steps = cutoff / stepsize

#        timearray = np.arange(0, cutoff + stepsize, stepsize)
        
#        while t <= cutoff: 
#            k1_S = gc(d_I[i], S[-1], I[-1], R[-1]) * stepsize
#            k1_I = hc(d_I[i], S[-1], I[-1], R[-1]) * stepsize
#            k1_R = rc(d_I[i], S[-1], I[-1], R[-1]) * stepsize

#            k2_S = gc(d_I[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, R[-1] + 1/2 * k1_R) * stepsize
#            k2_I = hc(d_I[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, R[-1] + 1/2 * k1_R) * stepsize
#            k2_R = rc(d_I[i], S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, R[-1] + 1/2 * k1_R) * stepsize

#            k3_S = gc(d_I[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I, R[-1] + 1/2 * k2_R) * stepsize
#            k3_I = hc(d_I[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I, R[-1] + 1/2 * k2_R) * stepsize
#            k3_R = rc(d_I[i], S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I, R[-1] + 1/2 * k2_R) * stepsize

#            k4_S = gc(d_I[i], S[-1] + k3_S, I[-1] + k3_I, R[-1] + k3_R) * stepsize
#            k4_I = hc(d_I[i], S[-1] + k3_S, I[-1] + k3_I, R[-1] + k3_R) * stepsize
#            k4_R = rc(d_I[i], S[-1] + k3_S, I[-1] + k3_I, R[-1] + k3_R) * stepsize

#            S_val = S[-1] + 1/6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S)
#            I_val = I[-1] + 1/6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I)
#            R_val = R[-1] + 1/6 * (k1_R + 2 * k2_R + 2 * k3_R + k4_R)

#            S.append(max(0, S_val))
#            I.append(max(0, I_val))
#            R.append(max(0, R_val))
#            N.append(S_val + I_val + R_val)

#            t += stepsize
        

#        Data = {'time':timearray, 'N':N, 'S':S, 'I':I, 'R': R}
#        df = pd.DataFrame(Data)
#        df.to_csv('Results/c_RK_d_I=%i.csv'%(d_I[i]))

#c_runge_katta(0.5,40)

#%%
def gc(S, I, N, e, d, d_I):
    return 0.5 * (N - S - I) - 4 * S * I / N - d * S + e * N

def hc(S, I, N, e, d, d_I):
    return 4 * S * I / N - 1 * I - d * I - d_I * I

def c_runge_katta(stepsize, cutoff, N, e=3, d=3, d_I = 0):
    
    I_0 = 100
    S_0 = 300
    R_0 = 0

    S = [S_0]
    I = [I_0]
    R = [R_0]

    t = stepsize

    numb_steps = cutoff / stepsize

    timearray = np.arange(0, cutoff + stepsize, stepsize)
    
    while t <= cutoff: 
        k1_S = gc(S[-1], I[-1], N, e, d, d_I) * stepsize
        k1_I = hc(S[-1], I[-1], N, e, d, d_I) * stepsize

        k2_S = gc(S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, N, e, d, d_I) * stepsize
        k2_I = hc(S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, N, e, d, d_I) * stepsize

        k3_S = gc(S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I, N, e, d, d_I) * stepsize
        k3_I = hc(S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I, N, e, d, d_I) * stepsize

        k4_S = gc(S[-1] + k3_S, I[-1] + k3_I, N, e, d, d_I) * stepsize
        k4_I = hc(S[-1] + k3_S, I[-1] + k3_I, N, e, d, d_I) * stepsize

        S_val = S[-1] + 1/6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S)
        I_val = I[-1] + 1/6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I)

        S.append(max(0, S_val))
        I.append(max(0, I_val))
        R.append(max(0, N - S_val - I_val))

        t += stepsize
    

    Data = {'time':timearray, 'S':S, 'I':I, 'R': R}
    df = pd.DataFrame(Data)
    df.to_csv('Results/Exercise C/c_RK_e=%.2f_d=%.2f_dI=%.2f.csv'%(e,d,d_I))



def call_c_RK():

    de_values = np.arange(0,3,0.1)
    for i in range(len(de_values)):
        c_runge_katta(0.5,50,N=400, e=de_values[i], d=de_values[i], d_I=0)

    di_values = np.arange(0,1,0.05)
    for j in range(len(di_values)):
        c_runge_katta(0.5,50,N=400, e=1, d=1 - di_values[j], d_I = di_values[j])

call_c_RK()

#%%
def S_to_I(S, I, R, stepsize, N, e, d, d_I):
    return 4 * S * I / N * stepsize

def I_to_R(S, I, R, stepsize, N, e, d, d_I):
    return 1 * I * stepsize 

def R_to_S(S, I, R, stepsize, N, e, d, d_I):
    return 0.5 * R * stepsize

def U_to_S(S, I, R, stepsize, N, e, d, d_I):
    return e * N * stepsize

def S_to_D(S, I, R, stepsize, N, e, d, d_I):
    return d * S * stepsize

def I_to_D(S, I, R, stepsize, N, e, d, d_I):
    return (d * I) * stepsize

#def I_to_DI(S, I, R, stepsize, N, e, d, d_I):
#    return (d_I * I) * stepsize

def R_to_D(S, I, R, stepsize, N, e, d, d_I):
    return d * R * stepsize

def c_MC(cutoff, N, e, d, d_I):
    
    stepsize = min([4 / (4 * N), 1 / (1 * N), 1 / (0.5 * N), 1 / (e * N), 1 / (d * N), 1 / ((d + d_I) * N)])

    numb_steps = int(cutoff / stepsize)

    timearray = np.arange(0, cutoff, stepsize)

    S = np.zeros(numb_steps)
    I = np.zeros(numb_steps)
    R = np.zeros(numb_steps)
    

    S[0] = 300
    I[0] = 100
    R[0] = 0
    change = 0
    examplechange = 0

    for j in range(numb_steps - 1):
        psi_p = random.uniform(0,1)
        pir_p = random.uniform(0,1)
        prs_p = random.uniform(0,1)
        #pus_p = random.uniform(0,1)
        
        psd_p = random.uniform(0,1)
        pid_p = random.uniform(0,1)
        #pidi_p = random.uniform(0,1)
        prd_p = random.uniform(0,1)

        #pus_p = psd_p + psi_p * pid_p + psi_p * pir_p * prd_p
        #pus_p = 1 - pid_p - psd_p - prd_p

        if U_to_S(S[j], I[j], R[j], stepsize, N, e, d, d_I) > pus_p:
            S[j + 1] = S[j] + 1
            change += 1
        else:
            S[j + 1] = S[j]

        if S_to_I(S[j], I[j], R[j], stepsize, N, e, d, d_I) > psi_p:
            S[j + 1] = S[j + 1] - 1
            I[j + 1] = I[j] + 1
        else:
            S[j + 1] = S[j + 1]
            I[j + 1] = I[j]

        if I_to_D(S[j], I[j], R[j], stepsize, N, e, d, d_I) > pid_p:
            I[j + 1] = I[j + 1] - 1
        else: 
            I[j + 1] = I[j + 1]

        #if I_to_DI(S[j], I[j], R[j], stepsize, N, e, d, d_I) > pidi_p:
        #    I[j + 1] = I[j + 1] - 1
        #else:
        #    I[j + 1] = I[j + 1]

        if I_to_R(S[j], I[j], R[j], stepsize, N, e, d, d_I) > pir_p:
            I[j + 1] = I[j + 1] - 1
            R[j + 1] = R[j] + 1            
            examplechange += 1
        else:
            I[j + 1] = I[j + 1]
            R[j + 1] = R[j]

        if R_to_D(S[j], I[j], R[j], stepsize, N, e, d, d_I) > prd_p:
            R[j + 1] = R[j + 1] - 1
        else:
            R[j + 1] = R[j + 1]

        if R_to_S(S[j], I[j], R[j], stepsize, N, e, d, d_I) > prs_p:
            R[j + 1] = R[j + 1] - 1
            S[j + 1] = S[j + 1] + 1
        else:
            R[j + 1] = R[j + 1]
            S[j + 1] = S[j + 1]

        if S_to_D(S[j], I[j], R[j], stepsize, N, e, d, d_I) > psd_p:
            S[j + 1] = S[j + 1] - 1 
        else:
            S[j + 1] = S[j + 1]

    print(change, numb_steps, examplechange)

    Data = {'time':timearray, 'S':S, 'I':I, 'R': R}
    df = pd.DataFrame(Data)
    df.to_csv('Results/cccccc_MC_e=%.2f_d=%.2f_di=%.2f.csv'%(e, d, d_I))

c_MC(10, 400, 1, 0.5, 0.5)


#%%
def g_N(S, I, R, N, a, b, c, e, d, d_I):
    return c * R - a * S * I / N - d * S + e * N
def h_N(S, I, R, N, a, b, c, e, d, d_I):
    return a * S * I / N - b * I - d * I - d_I * I
def r_N(S, I, R, N, a, b, c, e, d, d_I):
    return b * I - c * R - d * R
def n_N(S, I, R, N, a, b, c, e, d, d_I):
    return e * N - d * (S + I + R) - d_I * I


def jacobian(S, I, R, N, a, b, c, e, d, d_I):
    print(type(a))
    print(type(S))
    print(type(I))
    print(type(R))
    print(type(N))
    print(type(d_I))
    print(type(b))
    print(b)
    print(a * S / N - b - d - d_I)
    return [[-a * I / N - d, -a * S / N, c, a * S * I / N + e], \
        [a * I / N, a * S / N - b - d - d_I, 0, -a * S * I / N**2], \
        [0, b, -c - d, 0], \
        [-d, -d - d_I, -d, -e]]


def b_vector(S, I, R, N, a, b, c, e, d, d_I):
    return [-g_N(S, I, R, N, a, b, c, e, d, d_I), -h_N(S, I, R, N, a, b, c, e, d, d_I), -r_N(S, I, R, N, a, b, c, e, d, d_I), -n_N(S, I, R, N, a, b, c, e, d, d_I)]

def newton_method(cutoff, a, b, c, e, d, d_I):
    print(type(a))
    S = 1000
    I = 0.1
    R = 0.1
    N = 1000
    
    norm = 100
    
    while norm > cutoff:

        b_v = b_vector(S, I, R, N, a, b, c, e, d, d_I)
        print(type(S))
        print(type(a))
        J = jacobian(S, I, R, N, a, b, c, e, d, d_I)
        y = np.linalg.solve(J,b_v)
        
        S_p = np.copy(S)
        I_p = np.copy(I)
        R_p = np.copy(R)
        N_p = np.copy(N)
        
        S += y[0]
        I += y[1]
        R += y[2]
        N += y[3]
        
        norm = np.linalg.norm(np.array([S, I, R, N]) - np.array([S_p, I_p, R_p, N_p]))
        print(norm)
        
        
    return np.array([S, I, R, N])

newton_method(0.0000001, 4, 1, 0.5, 4, 3, 6)

#%%
