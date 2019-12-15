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

def a_runge_kutta(stepsize, cutoff, b = [1,2,3,4], N=400, a=4, c=0.5):
    '''''
    two folded runge kutta method for integrating functions g and h for a given set of b-values.
    '''''

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

#a_runge_kutta(0.5,40)
 
# %%
def S_to_I(b, S, I, R, stepsize, N=400, a=4, c=0.5):
    return a * S * I / N * stepsize

def I_to_R(b, S, I, R, stepsize, N=400, a=4, c=0.5):
    return b * I * stepsize 

def R_to_S(b, S, I, R, stepsize, N=400, a=4, c=0.5):
    return c * R * stepsize

def a_MC(cutoff, b = [1,2,3,4], N=400, a=4, c=0.5):
    '''''
    Monte Carlo method for a constant population using probabilites above.
    Evaluating one stepsize used, in order two let at max one person move
    compartment. Implementation for a given set of b-values.
    '''''
    
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

def gc(S, I, N, e, d, d_I):
    return 0.5 * (N - S - I) - 4 * S * I / N - d * S + e * N

def hc(S, I, N, e, d, d_I):
    return 4 * S * I / N - 1 * I - d * I - d_I * I

def nc(S, I, N, e, d, d_I):
    return e * N - d * I - d_I * I - d * S - d * (N - S - I)


def c_runge_kutta(stepsize, cutoff, e, d, d_I):
    '''''
    Runge Kutta method for three folded differential equations gc, hc and nc.
    A maximum value is always appended in order to avoid negative population.
    '''''
    I_0 = 100
    S_0 = 300
    R_0 = 0
    N_0 = 400

    S = [S_0]
    I = [I_0]
    R = [R_0]
    N = [N_0]

    t = stepsize

    numb_steps = cutoff / stepsize

    timearray = np.arange(0, cutoff + stepsize, stepsize)
    
    while t <= cutoff: 
        if N[-1] == 0:
            S.append(0)
            I.append(0)
            N.append(0)
            R.append(0)

        else:

            k1_S = gc(S[-1], I[-1], N[-1], e, d, d_I) * stepsize
            k1_I = hc(S[-1], I[-1], N[-1], e, d, d_I) * stepsize
            k1_N = nc(S[-1], I[-1], N[-1], e, d, d_I) * stepsize

            k2_S = gc(S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, N[-1] + 1/2 * k1_N, e, d, d_I) * stepsize
            k2_I = hc(S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, N[-1] + 1/2 * k1_N, e, d, d_I) * stepsize
            k2_N = nc(S[-1] + 1/2 * k1_S, I[-1] + 1/2 * k1_I, N[-1] + 1/2 * k1_N, e, d, d_I) * stepsize

            k3_S = gc(S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I, N[-1] + 1/2 * k2_N, e, d, d_I) * stepsize
            k3_I = hc(S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I, N[-1] + 1/2 * k2_N, e, d, d_I) * stepsize
            k3_N = nc(S[-1] + 1/2 * k2_S, I[-1] + 1/2 * k2_I, N[-1] + 1/2 * k2_N, e, d, d_I) * stepsize

            k4_S = gc(S[-1] + k3_S, I[-1] + k3_I, N[-1] + k3_N, e, d, d_I) * stepsize
            k4_I = hc(S[-1] + k3_S, I[-1] + k3_I, N[-1] + k3_N, e, d, d_I) * stepsize
            k4_N = nc(S[-1] + k3_S, I[-1] + k3_I, N[-1] + k3_N, e, d, d_I) * stepsize

            S_val = S[-1] + 1/6 * (k1_S + 2 * k2_S + 2 * k3_S + k4_S)
            I_val = I[-1] + 1/6 * (k1_I + 2 * k2_I + 2 * k3_I + k4_I)
            N_val = N[-1] + 1/6 * (k1_N + 2 * k2_N + 2 * k3_N + k4_N)

            S.append(max(0, S_val))
            I.append(max(0, I_val))
            N.append(max(0, N_val))
            R.append(max(0, N_val - S_val - I_val))

        t += stepsize

    Data = {'time':timearray, 'S':S, 'I':I, 'R': R, 'N':N}
    df = pd.DataFrame(Data)
    df.to_csv('Results/Alt Exercise C/c_RK_e=%.2f_d=%.2f_dI=%.2f.csv'%(e,d,d_I))
    return S, I, R, N

#c_runge_kutta(0.5,800,0.3,0.3,6.57)

#%%
def call_c_RK():
    '''''
    Calls c_RK for a set of e, d, and di. stable points are found and taken 
    to firstly be given by the last entrance of the time evalution.
    '''''

    e_values = np.append(np.arange(0,0.4,0.05),(np.append(np.arange(0.4,0.6,0.01), np.arange(0.6,0.6,0.05))))
    s_stable = []
    i_stable = []
    r_stable = []
    n_stable = []
    for i in range(len(e_values)):
        S, I, R, N = c_runge_katta(0.5,800, e=e_values[i], d=0.3, d_I=1.5)
        s_stable.append(np.mean(S[-1]))
        i_stable.append(np.mean(I[-1]))
        r_stable.append(np.mean(R[-1]))
        n_stable.append(np.mean(N[-1]))
    Data = {'e-values':e_values, 's*': s_stable, 'i*': i_stable, 'r*': r_stable, 'n*': n_stable}
    df = pd.DataFrame(Data)
    df.to_csv('Results/Alt Exercise C/stablepoints_e.csv')

    d_values = np.append(np.arange(0,0.2,0.05),(np.append(np.arange(0.2,0.4,0.01), np.arange(0.4,1,0.05))))
    s_stable = []
    i_stable = []
    r_stable = []
    n_stable = []
    for n in range(len(d_values)):
        S, I, R, N = c_runge_katta(0.5,800, e=0.5, d=d_values[n], d_I = 1.5)
        s_stable.append(np.mean(S[-1]))
        i_stable.append(np.mean(I[-1]))
        r_stable.append(np.mean(R[-1]))       
        n_stable.append(np.mean(N[-1]))
    Data = {'d-values':d_values, 's*': s_stable, 'i*': i_stable, 'r*': r_stable, 'n*': n_stable}
    df = pd.DataFrame(Data)
    df.to_csv('Results/Alt Exercise C/stablepoints_d.csv')


    di_values = np.arange(0,6.6,0.2)
    s_stable = []
    i_stable = []
    r_stable = []
    n_stable = []
    for j in range(len(di_values)):
        S, I, R, N = c_runge_katta(0.5,800, e=0.3, d=0.3, d_I = di_values[j])
        s_stable.append(np.mean(S[-1]))
        i_stable.append(np.mean(I[-1]))
        r_stable.append(np.mean(R[-1]))       
        n_stable.append(np.mean(N[-1]))
    Data = {'di-values':di_values, 's*': s_stable, 'i*': i_stable, 'r*': r_stable, 'n*': n_stable}
    df = pd.DataFrame(Data)
    df.to_csv('Results/Alt Exercise C/stablepoints_di.csv')

#call_c_RK()


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
    return ((d + d_I) * I) * stepsize

def R_to_D(S, I, R, stepsize, N, e, d, d_I):
    return d * R * stepsize

def c_MC(numb_steps, e, d, d_I):
    '''''
    Monte Carlo method implemented for vital dynamics as well. Therefore, 
    a new stepsize is defined for each iteration, due to varying N. 
    Population is saved not to become negative by letting it die out completely/
    become zero, as soon the given value is zero. Stable points of each random 
    walk are taken to be the mean of the last 70 entrances. 
    '''''
    
    #numb_steps = #int(cutoff / stepsize) 

    S = np.zeros(numb_steps)
    I = np.zeros(numb_steps)
    R = np.zeros(numb_steps)
    N = np.zeros(numb_steps)

    timevalue = 0
    time = np.zeros(numb_steps)
    
    S[0] = 300
    I[0] = 100
    R[0] = 0
    N[0] = 400 

    for j in range(numb_steps - 1):

        if N[j] == 0:
            S[j + 1] == 0
            I[j + 1] == 0
            R[j + 1] == 0
            N[j + 1] == 0

        else:
            stepsize = min([4 / (4 * N[j]), 1 / (1 * N[j]), 1 / (0.5 * N[j]), 1 / (e * N[j]), 1 / (d * N[j]), 1 / ((d + d_I) * N[j])])

            timevalue += stepsize
            time[j + 1] = timevalue

            psi_p = random.uniform(0,1)
            pir_p = random.uniform(0,1)
            prs_p = random.uniform(0,1)
            pus_p = random.uniform(0,1)
            psd_p = random.uniform(0,1)
            pid_p = random.uniform(0,1)
            prd_p = random.uniform(0,1)

            if U_to_S(S[j], I[j], R[j], stepsize, N[j], e, d, d_I) > pus_p:
                S[j + 1] = S[j] + 1
            else:
                S[j + 1] = S[j]

            if S_to_I(S[j], I[j], R[j], stepsize, N[j], e, d, d_I) > psi_p:
                S[j + 1] = S[j + 1] - 1
                I[j + 1] = I[j] + 1
            else:
                S[j + 1] = S[j + 1]
                I[j + 1] = I[j]


            if I_to_D(S[j], I[j], R[j], stepsize, N[j], e, d, d_I) > pid_p:
                I[j + 1] = I[j + 1] - 1
            else: 
                I[j + 1] = I[j + 1]


            if I_to_R(S[j], I[j], R[j], stepsize, N[j], e, d, d_I) > pir_p:
                I[j + 1] = I[j + 1] - 1
                R[j + 1] = R[j] + 1            
            else:
                I[j + 1] = I[j + 1]
                R[j + 1] = R[j]


            if R_to_D(S[j], I[j], R[j], stepsize, N[j], e, d, d_I) > prd_p:
                R[j + 1] = R[j + 1] - 1
            else:
                R[j + 1] = R[j + 1]

        
            if R_to_S(S[j], I[j], R[j], stepsize, N[j], e, d, d_I) > prs_p:
                R[j + 1] = R[j + 1] - 1
                S[j + 1] = S[j + 1] + 1
            else:
                R[j + 1] = R[j + 1]
                S[j + 1] = S[j + 1]

    
            if S_to_D(S[j], I[j], R[j], stepsize, N[j], e, d, d_I) > psd_p:
                S[j + 1] = S[j + 1] - 1 
            else:
                S[j + 1] = S[j + 1]
   
            N[j + 1] = S[j + 1] + I[j + 1] + R[j + 1]



    Data = {'time':time, 'S':S, 'I':I, 'R': R, 'N':N}
    df = pd.DataFrame(Data)
    df.to_csv('Results/MC Exercise C/ c_MC_e=%.2f_d=%.2f_di=%.2f.csv'%(e, d, d_I))

    #return np.mean(S[-70:]), np.mean(I[-70:]), np.mean(R[-70:]), np.mean(N[-70:])
    return S, I, R, N, time
#%%
def bootstrap(numb_rep, numb_steps, e, d, d_I):
    '''''
    Bootstrap for Monte Carlo method above. The stable points are found
    numb_rep number of times. Function return the mean and std of these stable points
    showing the precision/randomness of the random walks.
    '''''
    x = 0
    s_stables = []
    i_stables = []
    r_stables = []
    n_stables = []

    while x < numb_rep:
        s_star, i_star, r_star , n_star = c_MC(numb_steps, e, d, d_I)
        
        s_stables.append(s_star)
        i_stables.append(i_star)
        r_stables.append(r_star)
        n_stables.append(n_star)

        x += 1 

    Data = {'s_stables': s_stables, 'i_stables': i_stables, 'r_stables': r_stables, 'n_stables': n_stables}
    df = pd.DataFrame(Data)
    df.to_csv('Results/MC Exercise C/ stables_e=%.2f_d=%.2f_di=%.2f.csv'%(e, d, d_I))

    return np.mean(s_stables), np.mean(i_stables), np.mean(r_stables), np.mean(n_stables), \
        np.std(s_stables), np.std(i_stables), np.std(r_stables), np.std(n_stables)

#%%
def call_c_MC():
    '''''
    Function finds the total stable point for a given set of e, d, and di values.
    Number of repetitions is taken to be 30. And the random walks consist of 5000 steps.
    '''''
    e_values = np.append(np.arange(0,0.4,0.05),(np.append(np.arange(0.4,0.6,0.01), np.arange(0.6,0.6,0.05))))
    s_stable = []
    i_stable = []
    r_stable = []
    n_stable = []
    uncer_s_stable = []
    uncer_i_stable = []
    uncer_r_stable = []
    uncer_n_stable = []
    for i in trange(len(e_values)):
        S, I, R, N, u_s, u_i, u_r, u_n = bootstrap(30,5000,e_values[i], d = 0.3, d_I = 1.5)
        s_stable.append(S)
        i_stable.append(I)
        r_stable.append(R)
        n_stable.append(N)
        uncer_s_stable.append(u_s)
        uncer_i_stable.append(u_i)
        uncer_r_stable.append(u_r)
        uncer_n_stable.append(u_n)
    Data = {'e-values':e_values, 's*': s_stable, 'i*': i_stable, 'r*': r_stable, 'n*': n_stable, 'us*':uncer_s_stable, 'ui*': uncer_i_stable, 'ur*': uncer_r_stable, 'un*': uncer_n_stable }
    df = pd.DataFrame(Data)
    df.to_csv('Results/MC Exercise C/stablepoints_e.csv')

    d_values = np.append(np.arange(0,0.2,0.05),(np.append(np.arange(0.2,0.4,0.01), np.arange(0.4,1,0.05))))
    s_stable = []
    i_stable = []
    r_stable = []
    n_stable = []
    uncer_s_stable = []
    uncer_i_stable = []
    uncer_r_stable = []
    uncer_n_stable = []
    for n in trange(len(d_values)):
        S, I, R, N, u_s, u_i, u_r, u_n = bootstrap(30,5000, e=0.5, d=d_values[n], d_I = 1.5)
        s_stable.append(S)
        i_stable.append(I)
        r_stable.append(R)       
        n_stable.append(N)
        uncer_s_stable.append(u_s)
        uncer_i_stable.append(u_i)
        uncer_r_stable.append(u_r)
        uncer_n_stable.append(u_n)
    Data = {'d-values':d_values, 's*': s_stable, 'i*': i_stable, 'r*': r_stable, 'n*': n_stable, 'us*':uncer_s_stable, 'ui*': uncer_i_stable, 'ur*': uncer_r_stable, 'un*': uncer_n_stable }
    df = pd.DataFrame(Data)
    df.to_csv('Results/MC Exercise C/stablepoints_d.csv')


    di_values = np.arange(0,6.6,0.2) 
    s_stable = []
    i_stable = []
    r_stable = []
    n_stable = []
    uncer_s_stable = []
    uncer_i_stable = []
    uncer_r_stable = []
    uncer_n_stable = []
    for j in trange(len(di_values)):
        S, I, R, N, u_s, u_i, u_r, u_n = bootstrap(30,5000, e=0.3, d=0.3, d_I = di_values[j])
        s_stable.append(S)
        i_stable.append(I)
        r_stable.append(R)       
        n_stable.append(N)
        uncer_s_stable.append(u_s)
        uncer_i_stable.append(u_i)
        uncer_r_stable.append(u_r)
        uncer_n_stable.append(u_n)
    Data = {'di-values':di_values, 's*': s_stable, 'i*': i_stable, 'r*': r_stable, 'n*': n_stable, 'us*':uncer_s_stable, 'ui*': uncer_i_stable, 'ur*': uncer_r_stable, 'un*': uncer_n_stable }
    df = pd.DataFrame(Data)
    df.to_csv('Results/MC Exercise C/stablepoints_di.csv')


#%%
#newton method
def g_N(S, I, R, N, a, b, c, e, d, d_I):
    return c * R - a * S * I / N - d * S + e * N
def h_N(S, I, R, N, a, b, c, e, d, d_I):
    return a * S * I / N - b * I - d * I - d_I * I
def r_N(S, I, R, N, a, b, c, e, d, d_I):
    return b * I - c * R - d * R
def n_N(S, I, R, N, a, b, c, e, d, d_I):
    return e * N - d * (S + I + R) - d_I * I


def jacobian(S, I, R, N, a, b, c, e, d, d_I):
    return [[-a * I / N - d, -a * S / N, c, a * S * I / N + e], \
        [a * I / N, a * S / N - b - d - d_I, 0, -a * S * I / N**2], \
        [0, b, -c - d, 0], \
        [-d, -d - d_I, -d, -e]]


def b_vector(S, I, R, N, a, b, c, e, d, d_I):
    return [-g_N(S, I, R, N, a, b, c, e, d, d_I), -h_N(S, I, R, N, a, b, c, e, d, d_I), -r_N(S, I, R, N, a, b, c, e, d, d_I), -n_N(S, I, R, N, a, b, c, e, d, d_I)]

def newton_method(cutoff, a, b, c, e, d, d_I):
    S = 1000
    I = 0.1
    R = 0.1
    N = 1000
    
    norm = 100
    
    while norm > cutoff:

        b_v = b_vector(S, I, R, N, a, b, c, e, d, d_I)
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
        
        
    return np.array([S, I, R, N])

#newton_method(0.0000001, 4, 1, 0.5, 4, 3, 6)

#%%
