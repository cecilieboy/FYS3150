"""
script to visualize results of parallel computation
"""
#%%
import pandas as pd 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt 
import sys 
#%%
filename = sys.argv[1]
L = int(sys.argv[2])

data = np.loadtxt(filename)

df = pd.DataFrame(data, columns=["T", "E", "cv", "M", "chi"])

Y = [ "E", "cv", "M", "chi"]
y_label = [ "E/J", "$c_v$", "M", "$\chi$"]

for i,y in enumerate(Y):
    plt.figure(figsize=(10,10))
    sns.lineplot(x=df["T"], y= df[y]/L**2)
    plt.xlabel("T/$k_B$ J", fontsize=24)
    plt.ylabel(y_label[i]+'/$L^2$', fontsize=24)
    plt.savefig("L_%i"%L+y+".pdf")
#%%
