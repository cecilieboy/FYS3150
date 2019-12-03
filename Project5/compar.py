from flex_MC import SIR_MC, S_to_I, I_to_R, R_to_S
from flex_runge_kutta import runge_kutta, rhs_I, rhs_S
import numpy as np 
import matplotlib.pyplot as plt 

t_max = 50
X_0 = np.array([300,100,0])
kargs = [{'a':4, 'A':0, 'omega':1}, {'b':1}, {'c':0.5}, {'f':0}]
df_mc = SIR_MC(t_max, [S_to_I, I_to_R, R_to_S], kargs, [[(0,-1),(1,1)],
                                                                    [(1,-1), (2,1)],
                                                                    [(2,-1),(0,1)],
                                                                    [(0,-1), (2,1)]], X_0  )

df_rk = runge_kutta(0.01, t_max, X_0[:-1], [rhs_S, rhs_I], {})

col = ["tab:blue","tab:orange","tab:green"]
f =plt.figure(figsize=(10,10))
for i, name in zip([0,1,2],["S", "I", "R"]):
    plt.plot(df_mc["time"], df_mc[name], label=name, color = col[i])
    plt.plot(df_rk["time"], df_rk[name], color = col[i])
plt.legend(loc='best', fontsize = 28)
plt.xlabel("Time in a.u.", fontsize = 32)
plt.ylabel("Number", fontsize = 32)
plt.ylim(0,400)
plt.tick_params(size =24, labelsize=26)
plt.tight_layout()
plt.savefig('./Results/MC/compar.pdf')