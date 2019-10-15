import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt 




def gauss_quad_plots(df_legendre, df_laguerre):
    f =plt.figure(figsize=(10,10))
    #plt.plot(df_legendre["N"].to_numpy(), df_legendre["time"].to_numpy(), label = 'Legendre')
    plt.plot(df_laguerre["N"].to_numpy(), df_laguerre["time"].to_numpy(), label = 'Laguerre')
    plt.xlabel("N", fontsize = 28)
    plt.ylabel("t in s", fontsize = 28)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(loc='best', fontsize = 28)
    plt.savefig("Results/time.pdf")
    del f

    f =plt.figure(figsize=(10,10))
    #plt.plot(df_legendre["N"].to_numpy(), df_legendre["rel_err"].to_numpy(), label = 'Legendre')
    plt.plot(df_laguerre["N"].to_numpy(), df_laguerre["rel_err"].to_numpy(), label = 'Laguerre')
    plt.xlabel("N", fontsize = 28)
    plt.ylabel("$\\frac{|I-I_T|}{I_T}$", fontsize = 28)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(loc='best', fontsize = 28)
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
        plt.savefig("Results/"+name + "_err"+str(cy)+".pdf")
        del f, g

        grid = temp.pivot(index = 'cutoff', columns='samples', values = 'var_mc')

        f =plt.figure(figsize=(10,10))
        sns.set(font_scale=2)
        g =sns.heatmap(grid, vmin =var_min, vmax=var_max, cmap ='RdYlGn_r', annot=True, fmt='.3f', cbar = True, cbar_kws={'label':"Var(I)"},
                            xticklabels=["$10^{%i}$" %i for i in range(2,7)], yticklabels=[ str(i) for i in range(1,11)])
        plt.xlabel("MC Samples", fontsize =28)
        plt.ylabel("$\Lambda$", fontsize = 28)
        plt.savefig("Results/"+name + "_var"+str(cy)+".pdf")
        del f, g

df_laguerre = pd.read_csv("Results/laguerre.csv")
#df_legendre = pd.read_csv("Results/legendre.csv")
df_brutforce_mc = pd.read_csv("Results/brutforce_mc.csv")

gauss_quad_plots(None, df_laguerre)
plots_mc(df_brutforce_mc, "brutforce")