#%%
import numpy as np 
from matplotlib import pyplot as plt
from tqdm import trange 
import pandas as pd
import random 
#%%
b = [1,2,3,4]
for i in range(len(b)):
    df = pd.read_csv('Results/b=%i.csv'%b[i],usecols=['time','S','I','R'])
    print(np.shape(df))
    print(df)
    plt.figure()
    ax = plt.gca()
    df.plot.line(x='time',y='S',ax=ax)
    df.plot.line(x='time',y='I',ax=ax)
    df.plot.line(x='time',y='R',ax=ax)
    plt.ylabel('# Inhabitants')
    plt.yticks(np.arange(0,500,100))
    plt.savefig("./Results/runge_katta_b=%i.pdf"%b[i])


# %%
