from flex_MC import SIR_MC
from flex_runge_kutta import runge_kutta
import numpy as np 
import matplotlib.pyplot as plt 

def rhs_S(t, S, I, R, a= 4, b= 1, c= 0.5, d=0.1, dI=0, e=0.1, A =0, omega =1, f = 0):
    N = S + I + R
    a_t = A*np.cos(omega*t) +a
    #f = np.where(t>15, f, 0)
    return c * R - a_t * S * I / N - f +e*N -S*d

def rhs_I(t, S, I, R, a= 4, b= 1, c= 0.5, d=0.1, dI=0, e=0.1, A =0, omega =1, f = 0):
    N = S + I + R
    a_t = A*np.cos(omega*t) +a
    return a_t * S * I / N - b * I -d*I -dI*I

def rhs_R(t, S, I, R, a= 4, b= 1, c= 0.5, d=0.1, dI=0, e=0.1, A =0, omega =1, f = 0):
    #f = np.where(t>15, f, 0)
    return -c*R + b*I - d*R + f 

def S_to_I(timestep,t, S, I, R, a= 4, A =0, omega =1, N=400):
    #dirty fix to have proper time in cos
    a_t = A*np.cos(omega*t) +a
    return a_t*S*I/N *timestep

def I_to_R(timestep,t, S, I, R, b= 1):
    return b * I * timestep

def R_to_S(timestep,t, S, I, R, c= 0.5):
    return c*R*timestep

def S_to_R(timestep,t, S, I, R, f = 1):
    #f = np.where(t>15, f, 0)
    return f *timestep
def birthS(timestep,t, S, I, R, e=0.1):
    return e*timestep*(S+I+R)
def deathS(timestep,t, S, I, R, d=0.1):
    return d*timestep*S
def deathI(timestep,t, S, I, R, d=0.1, dI=0):
    return (d +dI)*timestep*I
def deathR(timestep,t, S, I, R, d=0.1):
    return d*timestep*R

t_max = 50
X_0 = np.array([14950,50,0])

A = 0
omega = 1
b = 0.9
c = 0.01
d = 0.05
dI = 0.1
e = 0.065
f = 0#1800
a = 0.026#2* (b-c)/b/c 
print(a)
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


df_mc = SIR_MC(t_max,probs , kargs, rules , X_0  )
#df_rk = runge_kutta(0.01, t_max, X_0, [rhs_S, rhs_I, rhs_R], [rk_karga, rk_karga, rk_karga], var_names=["S", "I", "R"])

col = ["tab:blue","tab:orange","tab:green"]
f =plt.figure(figsize=(10,10))
for i, name in zip([0,1,2],["S", "I", "R"]):
    plt.plot(df_mc["time"], df_mc[name], label=name, color = col[i])
    #plt.plot(df_rk["time"], df_rk[name], color = col[i], linestyle='--')
plt.plot(df_mc["time"], df_mc["S"]+ df_mc["I"]+ df_mc["R"], color='k', linestyle='-.', label='N')
#plt.plot([15,15],[0,30001],linestyle='--', color='r')
plt.legend(loc='best', fontsize = 28)
plt.xlabel("Time in years", fontsize = 32)
plt.ylabel("Number in 10,000", fontsize = 32)
plt.ylim(1,30000)#np.max(df_mc[["S","I","R"]].max()))
plt.yscale('log')
plt.xlim(0,t_max)
plt.tick_params(size =24, labelsize=26)
plt.tight_layout()
plt.savefig('./Results/MC/compar.pdf')

#0.5 = (1-b/a)/(1+b/c) <--> (1 - 0.5*(1+b/c))/b = 1/a <--> (0.5-0.5b/c)/b = 0.5(1/b-1/c) = 1/a