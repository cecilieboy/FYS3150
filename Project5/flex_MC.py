from numpy import random 
import numpy as np 
import pandas as pd
from flex_runge_kutta import save_plot 



def S_to_I(timestep,t, S, I, R, a= 4, A =0, omega =1, N=400):
    #dirty fix to have proper time in cos
    a_t = A*np.cos(omega*t) +a
    return a_t*S*I/N *timestep

def I_to_R(timestep,t, S, I, R, b= 1):
    return b * I * timestep

def R_to_S(timestep,t, S, I, R, c= 0.5):
    return c*R*timestep

def S_to_R(timestep,t, S, I, R, f = 1):
    return f *timestep

def SIR_MC(cutoff, transition_rule, kargs_transition_rule, update_rule, start_condition,  group_names = ["S", "I", "R"]):
    """
    Markov-Chain evaluation of the SIR model with flexible groups G0=S, G1=I,... for number of models, transitions and update rules
    Args:
        cutoff: maximal         time value
        transition_rule         list of transition rule for each variable each rule is a function transition(stepsize, G1, G2,..,Gn, **kargs_transition_rule)
        kargs_transition_rule   list of dicts containing the kargs for each of  the transition_rules
        update_rule             list of list of tuples containing the index ind_g (stating at 0) specifying which group to update and the amount by which this groupe is updated
                                i.e. [[(0, -1), (1, 1)],...] descripbing the transition from S to I (group at 0 G0=S reduced by 1 and group at 1 G1=I is added 1)
        start_condition         list of initial population in each group, the total population is inferred from it
    
    Kargs:
        group_names             list of names for groups for pandas datafram


    """
    num_var = len(start_condition)
    N = sum(start_condition)

    estimated_w_N = [N for i in range(num_var)]
    stepsize = min([ 1. / transition_rule[i](1, 0, *estimated_w_N, **kargs_transition_rule[i]) for i in range(num_var)])

    num_steps = int(cutoff/ stepsize)
    counts = np.zeros((num_var + 1, num_steps))
    counts[1:, 0] = start_condition

    #iteration for time
    for i in range(1, num_steps):
       rands = random.rand(len(transition_rule))
       counts[0, i] = i*stepsize
       #initially set group poplulation to previous time step
       counts[1: , i] = counts[1:, i-1]
       
       #iteration for updates
       for j in range(len(transition_rule)):  

           if rands[j] < transition_rule[j](stepsize, *counts[:,i-1], **kargs_transition_rule[j]):
               #update all relevant groups
               for update_index, update in update_rule[j]:
                   to_update =  max(0,counts[1 + update_index, i] + update) #ensures count
                   counts[1 + update_index, i] = to_update

    return pd.DataFrame(counts.T, columns=np.append(["time"], group_names)) 

def eval_model(par_to_vary, N = 400, X_0= [300, 100, 0], t_max = 15, add_kargs = {}):
    kargs = [{'a':4, 'A':0, 'omega':1}, {'b':1}, {'c':0.5}, {'f':0}]
    for p in par_to_vary[1]:
        for i in range(len(kargs)):
            if par_to_vary[0] in kargs[i]:
                kargs[i][par_to_vary[0]] = p
            for key in add_kargs:
                if key in kargs[i]:
                    kargs[i][key] = add_kargs[key]

        df = SIR_MC(t_max, [S_to_I, I_to_R, R_to_S, S_to_R], kargs, [[(0,-1),(1,1)],
                                                                    [(1,-1), (2,1)],
                                                                    [(2,-1),(0,1)],
                                                                    [(0,-1), (2,1)]], X_0  )
        save_plot(df, (par_to_vary[0], p), path = "./Results/MC/")

if __name__ == '__main__':
    eval_model(('b', [1,2,3,4]), t_max= 50)
    eval_model(('A', [0.5, 1,2 ,4]), t_max= 50)
    eval_model(('omega', [0.2, 0.7, 2, 4 ]), t_max= 50, add_kargs={'A':1})
    eval_model(('f', [0.5, 1, 2, 4]), t_max= 50)
