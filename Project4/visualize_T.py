"""
script to visualize results of parallel computation
"""
#%%
import pandas as pd 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt 
import sys 

filename = "paralell_L"
curve_label = "L "
cols = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
marker = ["o", "v", "s", "D"]
y_label = [ "E J$^{-1}$", "$c_v$", "|M| $\mu^{-1}$", "$\chi$"]
figname = ["E.pdf","cv.pdf", "M.pdf", "chi.pdf"]

ax = np.array([plt.subplots(1, figsize =(10,10)) for i in range(4)])
L = [40, 60, 80, 100]

for j,l in enumerate(L):

    T, E, var_E, cv, M, var_M, chi = np.loadtxt(filename + str(l)+'.txt', unpack=True)

    Y = [ E, cv, M, chi]
    err = [np.sqrt(var_E), np.sqrt(var_M)]

    for i in range(4):
        if (i == 0) or (i==2):
            ax[i,1].errorbar(T, Y[i], yerr =err[i - int(0.5*i)], label = curve_label + str(l), linestyle='none', marker= marker[j], color=cols[j])
        else:
            ax[i,1].plot(T, Y[i], label = curve_label + str(l), linestyle='none', marker= marker[j], color=cols[j])

for i in range(4):
    plt.figure(ax[i,0].number)
    ax[i,1].legend(loc='best', fontsize =22)
    ax[i,1].set_xlabel("T $k_B^{-1}$ J$^{-1}$", fontsize=24)
    ax[i,1].set_ylabel(y_label[i]+' $L^{-2}$', fontsize=24)
    if i ==1:
        ax[i,1].set_ylim(0,0.0014)
    ax[i,1].tick_params(axis = 'both',labelsize = 20)
    plt.tight_layout()
    plt.savefig(figname[i])
