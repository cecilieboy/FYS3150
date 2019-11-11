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
marker = ["o", "v", "s"]#, "D"]
cols = ['tab:blue', 'tab:orange', 'tab:green']#, 'tab:red']
y_label = [ "E $J^{-1} \ L^{-2}$", "$c_v\ k_b^{-1}\ L^{-4}$", "|M| $\mu^{-1}\ L^{-2}$", "$\chi J\ \mu^{-2}\ L^{-4}$"]
Y = [ "E", "cv", "M", "chi"]
#L = [40, 60, 80, 100]
L = [60, 80, 100]
df = pd.DataFrame(columns=["T", "E", "var E", "cv", "M", "var M", "chi", "L"])
TC = pd.DataFrame(columns=["L", "from", "TC max", "TC min", "TC mean" , "TC var"])
TC_temp = {"L":-100, "from":'none', "TC max":-1, "TC min":-2 , "TC mean":-3, "TC var": -4}
ind = 0
for l in L:

    data = np.loadtxt(filename + str(l)+'.txt')
    
    temp = pd.DataFrame( data, columns=["T", "E", "var E", "cv", "M", "var M", "chi"])
    temp["L"] = l
    df = df.append(temp)
    
    temp = temp[ (temp["T"]>1.9) & (temp["T"]<2.6)]
    temp=temp.groupby("T").agg({"E":'mean', "var E": 'sum', "cv": ["max", "min", "mean", "var"], "M":"mean", "var M": "sum", "chi":["max", "min", "mean", "var"]})
    TC_temp["L"] = l
    for func in ["cv", "chi"]:
        TC_temp["from"] = func
        for agg in ["max", "min", "mean", "var"]:
            TC_temp["TC " + agg] = temp[func, agg].idxmax()
        TC = TC.append(pd.DataFrame(TC_temp, index = [ind]))
        ind += 1

    del temp

TC["err from min"] =  TC["TC min"]- TC["TC mean"]
TC["err from max"] =  TC["TC max"]- TC["TC mean"]
TC["err low"] = - np.sqrt(TC[["err from min", "err from max"]].min(axis=1)**2 +0.01**2)
TC["err up"] = np.sqrt(TC[["err from min", "err from max"]].max(axis=1) +0.01**2)
print(TC)
for i in range(4):
    plt.figure(figsize=(10,10))
    ax = sns.lineplot(x= "T", y= Y[i],hue = "L", style="L", data=df, palette=cols, markers=marker, ci ='sd',err_style='bars')
    for j in range(3):
        ax.lines[2*j].set_linestyle('none')
    plt.xlabel("T $k_B$ J$^{-1}$", fontsize=24)
    plt.ylabel(y_label[i], fontsize=24)
    plt.tick_params(axis = 'both',labelsize = 20)
    plt.tight_layout()
    plt.savefig(Y[i]+ ".pdf")

lab = ["$c_v$", "$\chi$"]
plt.figure(figsize=(10,10))
for i, filt in enumerate(["cv", "chi"]):
    temp = TC[TC["from"]== filt]

    as_err = [temp["TC mean"] -temp["TC min"], temp["TC max"]- temp["TC mean"]]
    #plt.errorbar(temp["L"], temp["TC mean"], yerr=0.01, color = cols[0], linestyle='none', solid_capstyle='projecting', capsize=5)
    #plt.errorbar(temp["L"], temp["TC mean"], yerr=as_err, color = cols[1], linestyle='none', solid_capstyle='projecting', capsize=5)
    plt.errorbar(temp["L"], temp["TC mean"], yerr=[-temp["err low"], temp["err up"]],
                                         color = cols[i], linestyle='none', marker=marker[i], barsabove=False, label =lab[i]
                                         , solid_capstyle='projecting', capsize=10)
    
plt.legend(loc='best', fontsize =22)
plt.ylabel("$T_c\ k_B$ J$^{-1}$", fontsize=24)
plt.xlabel("L", fontsize=24)
plt.tick_params(axis = 'both',labelsize = 20)
plt.tight_layout()
plt.savefig("TC.pdf")

# %%
