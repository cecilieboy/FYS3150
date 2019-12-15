#%%
import numpy as np 
from matplotlib import pyplot as plt
from tqdm import trange 
import pandas as pd
import random 
fs = 32
ts = 24
#%%
b = [1,2,3,4]
for i in range(len(b)):
    df = pd.read_csv('Results/a_RK_b=%i.csv'%b[i],usecols=['time','S','I','R'])
    plt.figure()
    ax = plt.gca()
    df.plot.line(x='time',y='S',ax=ax)
    df.plot.line(x='time',y='I',ax=ax)
    df.plot.line(x='time',y='R',ax=ax)
    plt.ylabel('# Inhabitants')
    plt.yticks(np.arange(0,500,100))
    plt.xticks(np.arange(0,50,10))
    plt.savefig("./Results/runge_katta_b=%i.pdf"%b[i])


# %%
for i in range(len(b)):
    df = pd.read_csv('Results/a_MC_b=%i.csv'%b[i],usecols=['time','S','I','R'])
    plt.figure()
    ax = plt.gca()
    df.plot.line(x='time',y='S',ax=ax)
    df.plot.line(x='time',y='I',ax=ax)
    df.plot.line(x='time',y='R',ax=ax)
    plt.ylabel('# Inhabitants')
    plt.yticks(np.arange(0,500,100))
    plt.xticks(np.arange(0,50,10))
    plt.savefig("./Results/monte_carlo_b=%i.pdf"%b[i])

#%%
b = [1,2,3,4]
for i in range(len(b)):
    fs = 32
    ts = 24
    df = pd.read_csv('Results/a_RK_b=%i.csv'%b[i],usecols=['time','S','I','R'])
    plt.figure(figsize=(10,10))
    ax = plt.gca()
    df.plot.line(x='time',y='S', color = 'tab:blue', linestyle = '--', ax=ax, legend = False, fontsize = fs, label = '_nolegend_')
    df.plot.line(x='time',y='I', color = 'tab:orange', linestyle ='--', ax=ax, legend = False, fontsize = fs, label = '_nolegend_')
    df.plot.line(x='time',y='R', color = 'tab:green', linestyle = '--', ax=ax, legend = False, fontsize = fs, label = '_nolegend_')
    df = pd.read_csv('Results/a_MC_b=%i.csv'%b[i],usecols=['time','S','I','R'])
    df.plot.line(x='time',y='S',style = 'tab:blue', linewidth = 0.8, ax=ax, legend = True, fontsize = fs)
    df.plot.line(x='time',y='I', style = 'tab:orange', linewidth = 0.8, ax=ax, legend = True, fontsize = fs)
    df.plot.line(x='time',y='R', style = 'tab:green', linewidth = 0.8, ax=ax, legend = True, fontsize = fs)
    plt.legend(fontsize = ts)
    plt.yticks(np.arange(0,500,100), fontsize = ts)
    plt.xticks(np.arange(0,50,10), fontsize = ts)
    plt.xlabel('t in a.u.', fontsize = fs)
    plt.ylabel('Number', fontsize = fs)
    plt.tight_layout()
    plt.savefig("./Results/a_plots_b=%i.pdf"%b[i])
       
       

# %%

df = pd.read_csv('Results/Exercise C/stablepoints_d=e.csv',usecols=['de-values','s*','i*','r*'])
plt.figure(figsize=(10,10))
ax = plt.gca()
df.plot.line(x='de-values',y='s*',ax=ax)
df.plot.line(x='de-values',y='i*',ax=ax)
df.plot.line(x='de-values',y='r*',ax=ax)
plt.savefig("./Results/Exercise C/stablepoints_d=e.pdf")



df = pd.read_csv('Results/Exercise C/stablepoints_e=1.csv',usecols=['di-values','s*','i*','r*'])
plt.figure(figsize=(10,10))
ax = plt.gca()
df.plot.line(x='di-values',y='s*',ax=ax)
df.plot.line(x='di-values',y='i*',ax=ax)
df.plot.line(x='di-values',y='r*',ax=ax)
plt.savefig("./Results/Exercise C/stablepoints_e=1.pdf")


# %%
# RK e plot
df = pd.read_csv('Results/Alt Exercise C/stablepoints_e.csv', usecols=['e-values', 's*', 'i*', 'r*', 'n*'])
plt.figure(figsize=(10,10))
ax = plt.gca()
brug = df[16:20]
brug.plot.line(x='e-values',y='s*',ax=ax, label = 'S')
brug.plot.line(x='e-values',y='i*',ax=ax, label = 'I')
brug.plot.line(x='e-values',y='r*',ax=ax, label = 'R')
brug.plot.line(x='e-values',y='n*', style = 'k', ax=ax, label = 'N')
plt.xticks([0.477,0.48,0.49,0.5,0.51], ('0','0.48','0.49','0.50', '0.51'), fontsize = ts)
plt.yticks([0,60000,10.5**6],('0','309','inf'), fontsize = ts)
plt.xlabel('e in a.u.', fontsize = fs)
plt.ylabel('Number', fontsize = fs)
plt.legend(fontsize = ts)
plt.tight_layout()
plt.savefig("./Results/Alt Exercise C/stablepoints_e.pdf")
#%%
# MC e plot
df = pd.read_csv('Results/MC Exercise C/stablepoints_e.csv', usecols=['e-values', 's*', 'i*', 'r*', 'n*', 'us*', 'ui*', 'ur*', 'un*'])
plt.figure(figsize=(10,10))
ax = plt.gca()
df.plot.line(x='e-values',y='s*', yerr='us*', ax=ax, label = 'S')
df.plot.line(x='e-values',y='i*', yerr='ui*', ax=ax, label = 'I')
df.plot.line(x='e-values',y='r*', yerr='ur*', ax=ax, label = 'R')
df.plot.line(x='e-values',y='n*', color = 'k', yerr='un*', ax=ax, label = 'N')
plt.xticks([0,0.1,0.2,0.3,0.4,0.5], ('0','0.1','0.2','0.3', '0.4','0.5'), fontsize = ts)
plt.yticks([0,200,400,600,800,1000],fontsize = ts)
plt.xlabel('e in a.u.', fontsize = fs)
plt.ylabel('Number', fontsize = fs)
plt.legend(fontsize = ts)
plt.tight_layout()
plt.savefig("./Results/MC Exercise C/stablepoints_e.pdf")
# %%
# RK d plot
df = pd.read_csv('Results/Alt Exercise C/stablepoints_d.csv', usecols=['d-values', 's*', 'i*', 'r*', 'n*'])
plt.figure(figsize=(10,10))
ax = plt.gca()
brug = df[13:16]
brug.plot.line(x='d-values',y='s*',ax=ax, label = 'S')
brug.plot.line(x='d-values',y='i*',ax=ax, label = 'I')
brug.plot.line(x='d-values',y='r*',ax=ax, label = 'R')
brug.plot.line(x='d-values',y='n*', style = 'k', ax=ax, label = 'N')
plt.xticks([0.287,0.29,0.3,0.31], ('0.','0.29','0.30','0.31'), fontsize = fs)
plt.yticks([0,60000,10**6],('0','309','inf'), fontsize = ts)
plt.xlabel('d in a.u.', fontsize = fs)
plt.ylabel('Number', fontsize = fs)
plt.legend(fontsize = ts)
plt.tight_layout()
plt.savefig("./Results/Alt Exercise C/stablepoints_d.pdf")s

#%%
# MC d plot
df = pd.read_csv('Results/MC Exercise C/stablepoints_d.csv', usecols=['d-values', 's*', 'i*', 'r*', 'n*', 'us*', 'ui*', 'ur*', 'un*'])
plt.figure(figsize=(10,10))
ax = plt.gca()
df.plot.line(x='d-values',y='s*', yerr='us*', ax=ax, label = 'S')
df.plot.line(x='d-values',y='i*', yerr='ui*', ax=ax, label = 'I')
df.plot.line(x='d-values',y='r*', yerr='ur*', ax=ax, label = 'R')
df.plot.line(x='d-values',y='n*', color='k', yerr='un*', ax=ax, label = 'N')
plt.xticks([0,0.20,0.28,0.4,0.6,0.8], ('0','0.2','0.28','0.4','0.6','0.8'),fontsize = ts)
plt.yticks([0,200,400,600,800,1000,1200,1400],fontsize = ts)
plt.xlabel('d in a.u.', fontsize = fs)
plt.ylabel('Number', fontsize = fs)
plt.legend(fontsize = ts)
plt.tight_layout()
plt.savefig("./Results/MC Exercise C/stablepoints_d.pdf")

# %%
# RK di plot
df = pd.read_csv('Results/Alt Exercise C/stablepoints_di.csv', usecols=['di-values', 's*', 'i*', 'r*', 'n*'])
plt.figure(figsize=(10,10))
ax = plt.gca()
df.plot.line(x='di-values',y='s*',ax=ax, label = 'S')
df.plot.line(x='di-values',y='i*',ax=ax, label = 'I')
df.plot.line(x='di-values',y='r*',ax=ax, label = 'R')
df.plot.line(x='di-values',y='n*', style = 'k', ax=ax, label = 'N')
plt.xticks([0,1,2,2.8,4,5,6],('0','1','2','2.8', '4', '5', '6'), fontsize = ts)
plt.yticks([0,100,200,300,400], fontsize = ts)
plt.xlabel('di in a.u.', fontsize = fs)
plt.ylabel('Number', fontsize = fs)
plt.legend(fontsize = ts)
plt.tight_layout()
plt.savefig("./Results/Alt Exercise C/stablepoints_di.pdf")
#%%
# MC di plot
df = pd.read_csv('Results/MC Exercise C/stablepoints_di.csv', usecols=['di-values', 's*', 'i*', 'r*', 'n*', 'us*', 'ui*', 'ur*', 'un*'])
plt.figure(figsize=(10,10))
ax = plt.gca()
df.plot.line(x='di-values',y='s*', yerr='us*', ax=ax, label = 'S')
df.plot.line(x='di-values',y='i*', yerr='ui*', ax=ax, label = 'I')
df.plot.line(x='di-values',y='r*',yerr='ur*', ax=ax, label = 'R')
df.plot.line(x='di-values',y='n*', color = 'k', yerr='un*', ax=ax, label = 'N')
plt.xticks([0,1,2,3,4,5,6], fontsize = ts)
plt.yticks([0,100,200,300,400], fontsize = ts)
plt.xlabel('di in a.u.', fontsize = fs)
plt.ylabel('Number', fontsize = fs)
plt.legend(fontsize = ts)
plt.tight_layout()
plt.savefig("./Results/MC Exercise C/stablepoints_di.pdf")

# %%
plt.figure(figsize = (10,10))
plt.