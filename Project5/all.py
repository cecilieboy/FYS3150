"""
script which can vary all parameters
"""
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from flex_MC import SIR_MC, S_to_I, S_to_R, I_to_R, R_to_S
from flex_runge_kutta import runge_kutta

def rhs_S(t, S, I, R, a= 4, b= 1, c= 0.5, d=0.1, dI=0, e=0.1, A =0, omega =1, f = 0):
    N = S + I + R
    a_t = A*np.cos(omega*t) +a
    return c * R - a_t * S * I / N - f +e*N -S*d

def rhs_I(t, S, I, R, a= 4, b= 1, c= 0.5, d=0.1, dI=0, e=0.1, A =0, omega =1, f = 0):
    N = S + I + R
    a_t = A*np.cos(omega*t) +a
    return a_t * S * I / N - b * I -d*I -dI*I

def rhs_R(t, S, I, R, a= 4, b= 1, c= 0.5, d=0.1, dI=0, e=0.1, A =0, omega =1, f = 0):
    return -c*R + b*I - d*R + f 

def birthS(timestep,t, S, I, R, e=0.1):
    return e*timestep*(S+I+R)
def deathS(timestep,t, S, I, R, d=0.1):
    return d*timestep*S
def deathI(timestep,t, S, I, R, d=0.1, dI=0):
    return (d +dI)*timestep*I
def deathR(timestep,t, S, I, R, d=0.1):
    return d*timestep*R

t_max = 50
X_0 = np.array([300,100,0])
a = 4
A = 1
omega = 1
b = 2.1
c = 0.5
d = 0.3
dI = 0.1*d
e = 1.1*d
f = 20
indices = np.array([0]+[i**4 for i in range(2,20)])
kargs = [{'a':a, 'A':A, 'omega':omega},#Sto I
         {'b':b},   # I to R
         {'c':c}, #R to S
         {'f':f},  #S to R
         {'e':e}, # birth S
         {"d":d}, #death S
         {"d":d, "dI":dI}, # death I
         {"d":d}]     # death R
probs = [S_to_I, I_to_R, R_to_S, S_to_R, birthS, deathS, deathI, deathR]
rules = [[(0,-1),(1,1)],#S to I
        [(1,-1), (2,1)],#I to R
        [(2,-1),(0,1)],#R to S
        [(0,-1), (2,1)],# S to R
        [(0,1)],#birth S
        [(0,-1)],#death S
        [(1,-1)], #death I
        [(2,-1)]] # death R       
rk_karga = {'a':a, 'A':A, 'omega':omega,'b':b,'c':c,'f':f, 'e':e, "d":d, "dI":dI}
plt.figure(figsize=(10,10))
s = 50
x_start =  np.reshape(np.random.rand(2*s), (s,2))

for i in tqdm(range(s)):
    SI_frac = int(400*x_start[i,0])  
    I_frac = int(SI_frac*x_start[i,1])
    X_0 = [SI_frac -I_frac, I_frac, 400-SI_frac]
                                             
    df_mc = SIR_MC(t_max,probs , kargs, rules , X_0  )
   

    sample = indices[indices<df_mc.shape[0]]
    X = [df_mc["S"].iloc[sample], df_mc["I"].iloc[sample]]
    x = [df_mc["S"], df_mc["I"]]
    args  = df_mc.iloc[sample].to_numpy().T
    dX_dt = [rhs_S(*args, **rk_karga), rhs_I(*args, **rk_karga)]
    plt.plot(*x, color='k', alpha=0.1)
    plt.quiver(*X,*dX_dt, np.sqrt(dX_dt[0]**2 + dX_dt[1]**2), angles='xy')

plt.tick_params(size =24, labelsize=26)
plt.xlabel("S", fontsize =28)
plt.ylabel("I", fontsize =28)
plt.tight_layout()
plt.savefig('./Results/MC/phase_space_all_mc.pdf')


df_rk = runge_kutta(0.01, t_max, X_0, [rhs_S, rhs_I, rhs_R], [rk_karga, rk_karga, rk_karga], var_names=["S", "I", "R"])
col = ["tab:blue","tab:orange","tab:green"]
f =plt.figure(figsize=(10,10))
for i, name in zip([0,1,2],["S", "I", "R"]):
    plt.plot(df_mc["time"], df_mc[name], label=name, color = col[i])
    plt.plot(df_rk["time"], df_rk[name], color = col[i], linestyle = '--')
plt.legend(loc='best', fontsize = 28)
plt.xlabel("Time in a.u.", fontsize = 32)
plt.ylabel("Number", fontsize = 32)
plt.ylim(0,1.1*np.max(df_mc[["S", "I","R"]].max()))
plt.tick_params(size =24, labelsize=26)
plt.tight_layout()
plt.savefig('./Results/MC/all.pdf')
plt.close()
