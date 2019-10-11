from Project3 import GaussLaguerre, L_function
from scipy.special import roots_laguerre, roots_legendre 
import numpy as np
import pytest

def test_gausslaguerre():
    """
    compar own implementation aginst scipy.special
    """
    N = 20
    w, r = GaussLaguerre(N)
    r_e, w_e = roots_laguerre(N)
    res_w = np.abs(w-w_e)
    res_r  = np.abs(r - r_e)
    assert np.max(res_w) == pytest.approx(0, abs=10**-5)
    assert np.max(res_r) == pytest.approx(0)
def test_gausslegendre():
    """
    compar own implementation aginst scipy.special
    """
    N = 20
    w, r = L_function(N)
    r_e, w_e = roots_legendre(N)
    res_w = np.abs(w-w_e)
    res_r  = np.abs(r - r_e)
    assert np.max(res_w) == pytest.approx(0)
    assert np.max(res_r) == pytest.approx(0)