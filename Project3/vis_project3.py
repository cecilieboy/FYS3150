#%%

import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt 
import numpy as np
#font size controles
SMALL_SIZE = 20
MEDIUM_SIZE = 24
BIGGER_SIZE = 28
sns.set_context("paper", rc={"font.size":MEDIUM_SIZE,"axes.titlesize":MEDIUM_SIZE,"axes.labelsize":MEDIUM_SIZE, 'legend.fontsize':MEDIUM_SIZE,
                 'xticks':SMALL_SIZE, 'yticks':SMALL_SIZE})

#%%
def gauss_quad_plots(df_legendre, df_laguerre):
    f =plt.figure(figsize=(10,10))
    plt.plot(df_legendre["N"].to_numpy(), df_legendre["time"].to_numpy(), label = 'Legendre')
    plt.plot(df_laguerre["N"].to_numpy(), df_laguerre["time"].to_numpy(), label = 'Laguerre')
    plt.xlabel("N", fontsize = 28)
    plt.ylabel("t in s", fontsize = 28)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(loc='best', fontsize = 28)
    plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.1)
    plt.savefig("Results/time.pdf")
    del f

    max_rel =max([ df_laguerre["rel_err"].max()])#, df_legendre["rel_err"].max])
    f =plt.figure(figsize=(10,10))
    plt.plot(df_legendre["N"].to_numpy(), df_legendre["rel_err"].to_numpy(), label = 'Legendre')
    plt.plot(df_laguerre["N"].to_numpy(), df_laguerre["rel_err"].to_numpy(), label = 'Laguerre')
    plt.xlabel("N", fontsize = 28)
    plt.ylabel("$\\frac{|I-I_T|}{I_T}$", fontsize = 28)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(loc='best', fontsize = 28)    
    plt.subplots_adjust(left=0.17, right=0.97, top=0.97, bottom=0.1)
    plt.savefig("Results/rel_err.pdf")
    del f



def plots_mc(df, name, err_max = 0.25, var_max = 0.05):
    #time plots
    f = plt.figure(figsize= (10,10))
    for cy in [10,100,500]:
        filters = df["cycles"] == cy
        temp = df[filters]

        temp = temp.groupby("samples").mean()
        plt.plot(temp.index.values, temp["time"], label = cy)
        del temp

    plt.xlabel("MC samples", fontsize = 28)
    plt.ylabel("t in s", fontsize = 28)
    plt.xticks(ticks=np.array([10**i for i in range(2,7)]), labels=["$10^{%i}$" %i for i in range(2,7)],fontsize=24)
    plt.yticks(fontsize=24)
    plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.1)
    plt.legend(title ='Validation cycles' , loc='best', fontsize = 28)
    plt.savefig("Results/"+ name+"_time.pdf")
    del f

    err_min = df["rel_err"].min()
    var_min = df["var_mc"].min()
    
    for cy in [10,100,500]:
        filters = (df["cycles"] == cy) 
        temp = df[filters]

        grid = temp.pivot(index = 'cutoff', columns='samples', values = 'rel_err')

        f =plt.figure(figsize=(10,10))
        sns.set(font_scale=2)
        g =sns.heatmap(grid, vmin =err_min,vmax= err_max, cmap ='RdYlGn_r', annot=True, fmt='.3f', cbar = True, cbar_kws={'label':"$\\frac{|I-I_T|}{I_T}$"}, 
                        label='big', xticklabels=["$10^{%i}$" %i for i in range(2,7)], yticklabels=[ str(i) for i in range(1,11)])
        plt.xlabel("MC Samples", fontsize =28)
        plt.ylabel("$\Lambda$", fontsize = 28)
        plt.subplots_adjust(left=0.1, right=0.92, top=0.97, bottom=0.1)
        plt.savefig("Results/"+name + "_err"+str(cy)+".pdf")
        del f, g

        grid = temp.pivot(index = 'cutoff', columns='samples', values = 'var_mc')

        f =plt.figure(figsize=(10,10))
        sns.set(font_scale=2)
        g =sns.heatmap(grid, vmin =var_min, vmax=var_max, cmap ='RdYlGn_r', annot=True, fmt='.3f', cbar = True, cbar_kws={'label':"Var(I)"},
                            xticklabels=["$10^{%i}$" %i for i in range(2,7)], yticklabels=[ str(i) for i in range(1,11)])
        plt.xlabel("MC Samples", fontsize =28)
        plt.ylabel("$\Lambda$", fontsize = 28)
        plt.subplots_adjust(left=0.1, right=0.95, top=0.97, bottom=0.1)
        plt.savefig("Results/"+name + "_var"+str(cy)+".pdf")
        del f, g

def parallel_plots():
    df_exp_mc_100 = pd.read_csv('Results/main_results_100.csv')
    df_mc = df_exp_mc_100[["Samples(N)","Time", "mean of integrals", "Variance"]]
    df_mc = df_mc.rename(columns={"Time":"time sum", "mean of integrals":"I mean", "Samples(N)":"mc" ,"Variance": "I var"})
    df_mc.insert(0, "threads" , ["original"  for i in df_mc.index.values], True)

    data = np.loadtxt("time_array.txt" )
    df  = pd.DataFrame(data= data, index=np.arange(0,len(data)), columns=["threads", "time", "I", "mc"])
    df = df.groupby(["threads", "mc"]).agg({"time":np.sum,
                                    "I":[np.mean, np.var]})
    df.columns = [' '.join(col).strip() for col in df.columns.values]
    df = df.reset_index()
    
    df["threads"].where(df["threads"] != 1 , 'vectorized', True)
    df = df.append(df_mc)
    
    plt.figure(figsize=(10,10))
    sns.lineplot(x="MC Samples", y="Time in s", hue = "Computation",
                 data=df.rename(columns={"time sum": "Time in s", "mc":"MC Samples", "threads":"Computation"}))
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig("Results/parallel.pdf")

parallel_plots()

df_laguerre = pd.read_csv("Results/laguerre.csv")
df_legendre = pd.read_csv("Results/legendre.csv")
df_brutforce_mc = pd.read_csv("Results/brutforce_mc.csv")
df_exp_mc_10 = pd.read_csv('Results/main_results_10.csv')
df_exp_mc_100 = pd.read_csv('Results/main_results_100.csv')
#df_exp_mc_200 = pd.read_csv('Results/main_results_200.csv')

gauss_quad_plots(df_legendre, df_laguerre)
plots_mc(df_brutforce_mc, "brutforce")

#%% 
#Visualization Importance sampling MC

#Time plot:
fig, ax = plt.subplots(figsize=(10, 10))
ax = plt.gca()
plot_time = pd.read_csv("Results/main_results_10.csv")
plot_time.plot(kind='line',x='log10_Samples(N)',y='Time',ax=ax,label = 10)
plot_time = pd.read_csv("Results/main_results_100.csv")
plot_time.plot(kind='line',x='log10_Samples(N)',y='Time', color='red', ax=ax,label = 100)
plot_time = pd.read_csv("Results/main_results_200.csv")
plot_time.plot(kind='line',x='log10_Samples(N)',y='Time', color='green', ax=ax,label = 200)
plt.xlabel("$log_{10}$(MC Samples)", fontsize =20)
plt.ylabel("t in s", fontsize = 20)
plt.legend(title = 'MC experiments', loc='best', fontsize = 20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.17, right=0.97, top=0.97, bottom=0.1)
plt.savefig("Results/time_ceci_mc.pdf")
plt.show()

#Rel error plot:
fig, ax = plt.subplots(figsize=(10, 10))
ax = plt.gca()
#ax.set_facecolor('white')
#fig.patch.set_facecolor('white')
#ax.patch.set_facecolor('white')
plot_rel_error = pd.read_csv("Results/main_results_10.csv")
plot_rel_error.plot(kind='line',x='log10_Samples(N)',y='Relative Error',ax=ax,label = 10)
plot_rel_error = pd.read_csv("Results/main_results_100.csv")
plot_rel_error.plot(kind='line',x='log10_Samples(N)',y='Relative Error', color='red', ax=ax,label = 100)
plot_rel_error = pd.read_csv("Results/main_results_200.csv")
plot_rel_error.plot(kind='line',x='log10_Samples(N)',y='Relative Error', color='green', ax=ax,label = 200)
plt.xlabel("$log_{10}$(MC Samples)", fontsize =20)
plt.ylabel("$\\frac{|I-I_T|}{I_T}$", fontsize = 28)
plt.legend(title = 'MC experiments', loc='best', fontsize = 20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.subplots_adjust(left=0.17, right=0.97, top=0.97, bottom=0.1)
plt.savefig("Results/Rel_error_ceci_mc.pdf")
plt.show()

#Variance plot:
fig, ax = plt.subplots(figsize=(10, 10))
ax = plt.gca()
plot_time = pd.read_csv("Results/main_results_10.csv")
plot_time.plot(kind='line',x='log10_Samples(N)',y='Variance',ax=ax,label = 10)
plot_time = pd.read_csv("Results/main_results_100.csv")
plot_time.plot(kind='line',x='log10_Samples(N)',y='Variance', color='red', ax=ax,label = 100)
plot_time = pd.read_csv("Results/main_results_200.csv")
plot_time.plot(kind='line',x='log10_Samples(N)',y='Variance', color='green', ax=ax,label = 200)
plt.xlabel("$log_{10}$(MC Samples)", fontsize =20)
plt.ylabel("Var(I)", fontsize = 20)
plt.legend(title = 'MC experiments', loc='best', fontsize = 20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("Results/var_ceci_mc.pdf")
plt.show()



#%%

#Time plot:
fig, ax = plt.subplots(figsize=(10, 10))
ax = plt.gca()
plot_time = pd.read_csv("Results/main_results_10.csv")
plot_time.plot(kind='line',x='Samples(N)',y='Time',ax=ax,label = 10)
plot_time = pd.read_csv("Results/main_results_100.csv")
plot_time.plot(kind='line',x='Samples(N)',y='Time', color='red', ax=ax,label = 100)
plot_time = pd.read_csv("Results/main_results_200.csv")
plot_time.plot(kind='line',x='Samples(N)',y='Time', color='green', ax=ax,label = 200)
plt.xlabel("$log_{10}$(MC Samples)", fontsize =20)
plt.ylabel("t in s", fontsize = 20)
plt.legend(title = 'MC experiments', loc='best', fontsize = 20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.17, right=0.97, top=0.97, bottom=0.1)
plt.savefig("Results/time_lin_ceci_mc.pdf")
plt.show()



#%%
