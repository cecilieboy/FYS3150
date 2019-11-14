"""
script to visualize results of parallel computation and find TC(L=inf)
"""
#%%
import pandas as pd 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt 
from scipy.stats import mode
from scipy.optimize import curve_fit
import sys 

def TCL(L, a, TCinf):
    return a/L + TCinf

fs =32
ts=24

filename = "paralell_L"
curve_label = "L "
marker = ["o", "v", "s", "D"]
cols = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
y_label = [ "E / $J^{-1} \ L^{-2}$", "$c_v\ /\ k_b^{-1}\ L^{-4}$", "|M| / $\mu\ L^{-2}$", "$\chi\ /\ J\ \mu^2\ L^{-4}$"]
Y = [ "E", "cv", "M", "chi"]
L = [40, 60, 80, 100]
df = pd.DataFrame(columns=["T", "E", "var E", "cv", "M", "var M", "chi", "L"])

def mode_app(x):
    return x.mode().mean()



central_value='median'
TC = pd.DataFrame(columns=["L", "from", "TC max", "TC min", "TC "+central_value , "TC var"])

TC_temp = {"L":-100, "from":'none', "TC max":-1, "TC min":-2 , "TC " + central_value:-3, "TC var": -4}
ind = 0

for l in L:

    data = np.loadtxt(filename + str(l)+'.txt')
    
    temp = pd.DataFrame( data, columns=["T", "E", "var E", "cv", "M", "var M", "chi"])
    temp["L"] = l
    df = df.append(temp)
    
    temp = temp[ (temp["T"]>1.9) & (temp["T"]<2.6)]
    if central_value == 'mode_app':
        temp=temp.groupby("T").agg({"E":'mean', "var E": 'sum', "cv": ["max", "min", mode_app, "var"],
                                    "M":"mean", "var M": "sum", "chi":["max", "min", mode_app, "var"]})
    else:
        temp=temp.groupby("T").agg({"E":'mean', "var E": 'sum', "cv": ["max", "min", central_value, "var"],
                                   "M":"mean", "var M": "sum", "chi":["max", "min", central_value, "var"]})
        
    TC_temp["L"] = l
    for func in ["cv", "chi"]:
        TC_temp["from"] = func
        for agg in ["max", "min", central_value, "var"]:
            TC_temp["TC " + agg] = temp[func, agg].idxmax()
        TC = TC.append(pd.DataFrame(TC_temp, index = [ind]))
        ind += 1

    del temp

TC["err from min"] =  TC["TC min"]- TC["TC " + central_value]
TC["err from max"] =  TC["TC max"]- TC["TC " + central_value]
#TC["err low"] = - np.sqrt(TC[["err from min", "err from max"]].min(axis=1)**2 +0.01**2)
#TC["err up"] = np.sqrt(TC[["err from min", "err from max"]].max(axis=1) +0.01**2)
TC["err low"] = -0.01 + np.where(TC["err from max"] * TC["err from min"] <= 0, TC[["err from min", "err from max"]].min(axis=1), 0 )
TC["err up"] =  TC[["err from min", "err from max"]].max(axis=1) +0.01
print(TC.to_latex())
for i in range(4):
    plt.figure(figsize=(10,10))
    ax = sns.lineplot(x= "T", y= Y[i],hue = "L", style="L", data=df, palette=cols, markers=marker,dashes = False, ci =68 ,err_style='bars')
    for j in range(4):
        ax.lines[2*j].set_linestyle('none')
    plt.legend(fontsize = fs-2)
    plt.xlabel("T $k_B$ J$^{-1}$", fontsize=fs)
    plt.ylabel(y_label[i], fontsize=fs)
    plt.tick_params(axis = 'both',labelsize = ts)
    plt.tight_layout()
    plt.savefig(Y[i]+ ".pdf")

lab = ["$c_v$", "$\chi$"]
plt.figure(figsize=(10,10))
for i, filt in enumerate(["cv", "chi"]):
    temp = TC[TC["from"]== filt]
    L= np.linspace(40,100,50)
    s = (temp["err up"]- temp["err low"])/2
    paropt, cov_opt = curve_fit( TCL,temp["L"], temp["TC " +central_value],p0=[1, 2.26], sigma = s, absolute_sigma=True)
    chisq = np.sum((temp["TC " + central_value ] - TCL(temp["L"],*paropt))**2/s)
    print(filt, paropt[1], np.sqrt(cov_opt[1,1]), chisq)
    plt.plot(L, TCL(L,*paropt), color= cols[i])
    #plt.errorbar(temp["L"], temp["TC " + central_value], yerr=np.sqrt(temp["TC var"] +0.01**2), color = cols[1], linestyle='none', solid_capstyle='projecting', capsize=5)
    plt.errorbar(temp["L"], temp["TC " + central_value], yerr=[-temp["err low"], temp["err up"]],
                                         color = cols[i], linestyle='none', marker=marker[i], barsabove=False, label =lab[i]
                                         , solid_capstyle='projecting', capsize=10)
    
plt.legend(loc='best', fontsize =fs-2)
plt.ylabel("$T_c\ k_B$ J$^{-1}$", fontsize=fs)
plt.xlabel("L", fontsize=fs)
plt.tick_params(axis = 'both',labelsize = ts)
plt.tight_layout()
plt.savefig("TC_%s.pdf" % central_value)

# %%
