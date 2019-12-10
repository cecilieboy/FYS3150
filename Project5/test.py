import pytest 
from flex_runge_kutta import runge_kutta, rhs_S, rhs_I
from flex_MC import SIR_MC, S_to_I, I_to_R, R_to_S
import numpy as np 

def test_runge_kutta():
    """
    testing runge kutta stable points
    for a=4, b=1, c=0.5 -> S* = 0.25*400, I* = 0.75/3*400
    """
    a=4
    b=1
    c=0.5
    t_max = 50
    X_0 = np.array([300,100,0])
    #kargs = [{'a':4, 'A':0, 'omega':1}, {'b':2.2}, {'c':0.5}, {'f':20}]
    df_rk = runge_kutta(0.01, t_max, X_0[:-1], [rhs_S, rhs_I], [{},{}])#no kargs because usage of default values
    mean = df_rk[df_rk["time"]>30].mean()
    
    assert np.abs(mean["S"]-400*b/a) == pytest.approx(0, abs=1e-8)
    assert np.abs(mean["I"]-400*(1-b/a)/(1+b/c)) == pytest.approx(0, abs=1e-8)

def test_MC():
    """
    testing runge kutta stable points
    for a=4, b=1, c=0.5 -> S* = 0.25*400, I* = 0.75/3*400
    """
    a=4
    b=1
    c=0.5
    t_max = 50
    X_0 = np.array([300,100,0])
    kargs = [{'a':4, 'A':0, 'omega':1}, {'b':1}, {'c':0.5},]
    stab = np.zeros((20,3))
    for i in range(20):
        df_mc = SIR_MC(t_max, [S_to_I, I_to_R, R_to_S], kargs, [[(0,-1),(1,1)],
                                                                [(1,-1), (2,1)],
                                                                [(2,-1),(0,1)]] , X_0  )
        stab[i] = df_mc[df_mc["time"]>30].drop(columns=["time"]).mean()
    mean = np.mean(stab,axis=0)
    std = np.std(stab, axis=0)
    assert np.abs(mean[0]-400*b/a) == pytest.approx(0, abs=std[0])
    assert np.abs(mean[1]-400*(1-b/a)/(1+b/c)) == pytest.approx(0, abs=std[1])
    assert np.abs(mean[2]-400*b/c*(1-b/a)/(1+b/c)) == pytest.approx(0, abs=std[2])