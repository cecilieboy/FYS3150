from flex_MC import SIR_MC, S_to_I, I_to_R, R_to_S, S_to_R
from flex_runge_kutta import runge_kutta, rhs_I, rhs_S
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
from tqdm import tqdm
def var_b(b):
    figname ='_f%i' %(100*b)
    t_max = 50
    X_0 = np.array([300,100,0])
    kargs = [{'a':4, 'A':0, 'omega':0.1}, {'b':1}, {'c':0.5}, {'f':b}]
    df_mc = SIR_MC(t_max, [S_to_I, I_to_R, R_to_S], kargs, [[(0,-1),(1,1)],
                                                                        [(1,-1), (2,1)],
                                                                        [(2,-1),(0,1)],
                                                                        [(0,-1), (2,1)]], X_0  )

    df_rk = runge_kutta(0.01, t_max, X_0[:-1], [rhs_S, rhs_I], [{"b":b, "A":1},{"b":b, "A":1}])
    df_rk["R"] = X_0.sum() - df_rk["S"] - df_rk["I"]

    df_stable = df_mc[df_mc["time"]>20]
    dt = df_stable["time"].iloc[1]- df_stable["time"].iloc[0]
    freq = np.fft.fftfreq(len(df_stable["S"]), d= dt)
    freq = freq[freq >= 0]

    scaled_auto = 0
    for  name1 in ["S", "I", "R"]:
        fft1 = np.fft.rfft(df_stable[name1])
        fft2 = np.fft.rfft(df_stable[name1])
        cov = fft1*fft2
        cor = np.real(cov* np.conj(cov))
        try:
            scaled = cor[:]/np.max(cor[1:])
        except:
            scaled = np.zeros(cor.shape)
        scaled_auto += scaled/3
   
    return scaled_auto, freq

b = np.linspace(0,100,10)
b = np.sort(b)
bootstraps = 20
c, f = var_b(b[0])
mask_f = (f>0.01) & (f<1)
c = 0*c[mask_f]
f = f[mask_f]
t_max = 70
X_0 = np.array([300,100,0])
kargs = [{'a':4, 'A':0, 'omega':1}, {'b':2.2}, {'c':0.5}, {'f':b[0]}]

df = SIR_MC(t_max, [S_to_I, I_to_R, R_to_S, S_to_R], kargs, [[(0,-1),(1,1)],
                                                                        [(1,-1), (2,1)],
                                                                        [(2,-1),(0,1)],
                                                                        [(0,-1), (2,1)]], X_0 )
rk_stable = np.zeros((3,len(b)))
df["b"] = b[0]
for i in tqdm(range(0,len(b))):
    temp = 0
    kargs = [{'a':4, 'A':0, 'omega':1}, {'b':2.2}, {'c':0.5}, {'f':b[i]}]

    df_rk = runge_kutta(0.01, t_max,X_0[:-1], [rhs_S, rhs_I], [{"b":2.2, "A":0, "f":b[i]},{"b":2.2,"f":b[i], "A":0}])
    df_rk["R"] = 400 - df_rk["I"]- df_rk["S"]
    df_rk = df_rk[df_rk["time"]>50]
    df_rk.drop(columns=["time"], inplace=True)
    rk_stable[:, i] = df_rk.mean()

    for j in range(bootstraps):
        if (j == 0) & (i==0):
            continue 
        temp = SIR_MC(t_max, [S_to_I, I_to_R, R_to_S, S_to_R], kargs, [[(0,-1),(1,1)],
                                                                        [(1,-1), (2,1)],
                                                                        [(2,-1),(0,1)],
                                                                       [(0,-1), (2,1)]], X_0  )
        temp = temp[temp["time"]>50]
        temp["b"] = b[i]
        df = df.append(temp)
        #temp1, _ = var_b(b[i])
        #temp1 = np.nan_to_num(temp1[mask_f], nan=-0.1)
        #temp += temp1/bootstraps
    #c = np.vstack([c,temp])

#temp0 = df[df["f"]==b[-4]].copy()
#select stable times  
df_stab = df#[df["time"]>30]
#for histogram
#temp1 = df[df["f"]==b[0]].copy()
#temp2 = df[df["f"]==b[-1]].copy()

df_stab.to_csv("./Results/MC/b_var_os/measure.csv")
df_stab = df_stab.groupby("b").agg(['mean','std'])   

col = ["tab:blue", "tab:orange","tab:green"]
plt.figure(figsize=(10,10))
for i, name in enumerate(["S", "I", "R"]):
    plt.plot(b, rk_stable[i], linestyle ='--', color = col[i])
    plt.errorbar(b, df_stab[name,'mean'], yerr=df_stab[name,'std'], linestyle='', marker='o', color=col[i], label=name)
plt.legend(loc='best', fontsize = 28)
plt.xlabel("f in a.u.", fontsize = 32)
plt.ylabel("Number", fontsize = 32)
plt.ylim(0,400)
plt.tick_params(size =24, labelsize=26)
plt.tight_layout()
plt.savefig('./Results/MC/b_var_os/stable_b_os.pdf')

"""
ind=1
for t in [temp1, temp2]:
    plt.figure(figsize=(10,10))
    for i, name in enumerate(["S", "I", "R"]):
        plt.hist(t[name], bins = 50, density =True, label=name, alpha =0.5)
    plt.legend(loc='best', fontsize = 28)
    plt.xlabel("f in a.u.", fontsize = 32)
    plt.ylabel("Number", fontsize = 32)
    plt.tick_params(size =24, labelsize=26)
    plt.tight_layout()
    plt.savefig('./Results/MC/f_var/stable_f%i_dist_.pdf'%ind)
    ind+=2

df_rk = runge_kutta(0.01, t_max, X_0[:-1], [rhs_S, rhs_I], [{"f":b[-4], "A":0, "omega":0.1, "b":1},{"f":b[-4], "A":0, "omega":0.1, "b":1}])
df_rk["R"] = 400 - df_rk["S"] - df_rk["I"]
temp0 = temp0.groupby("time").agg(['mean', 'std'])
plt.figure(figsize=(10,10))
for i, name in enumerate(["S", "I", "R"]):
    plt.plot(df_rk["time"], df_rk[name], color = col[i], linestyle='--')
    plt.plot(temp0.index, temp0[name,'mean'], color = col[i], label = name)
    plt.fill_between(temp0.index, temp0[name, 'mean']-temp0[name, 'std'], temp0[name, 'mean']+temp0[name, 'std'], color= col[i], alpha=0.4)
plt.legend(titel='f = %i'%b[-4],loc='best', fontsize = 28)
plt.xlabel("t in a.u.", fontsize = 32)
plt.ylabel("Number", fontsize = 32)
plt.ylim(0,400)
plt.tick_params(size =24, labelsize=26)
plt.tight_layout()
plt.savefig('./Results/MC/f_var/av_f1.pdf')
"""
"""
np.savetxt("./Results/MC/b_var/spectrum.txt", c)

c = np.loadtxt("./Results/MC/b_var/spectrum.txt")
F, B = np.meshgrid(f[1:],b)


plt.figure(figsize=(10,10))
c =plt.contourf(F,B, c[1: ,1:],levels=200, vmin = -0.1, vmax = c.max())
cbar = plt.colorbar(c, label ="$\langle C/C_0\\rangle$")
cbar.ax.tick_params(labelsize =26, size =20)
cbar.ax.set_ylabel("$\langle C/C_0\\rangle$",fontsize =32)
cbar.set_ticks([-0.1, 0.3, 0.6,1])
plt.xlabel("$\omega$ in a.u.", fontsize = 32)
plt.tick_params(size =20, labelsize=26)
plt.ylabel('b in a.u.', fontsize=32)
plt.xlim(0.05,0.5)
plt.tight_layout()
for cs in c.collections:
    cs.set_edgecolor("face")
plt.savefig("./Results/MC/b_var/spectrum.pdf")
"""