from Project3 import GaussLaguerre
from scipy.special import roots_laguerre
import numpy as np
import pytest

def test_gausslaguerre():
    """
    compar own implementation aginst scipy.special,
    weights in scipy.special are not multiplyed with 2!
    """
    N = 10
    w, r = GaussLaguerre(N)
    r_e, w_e = roots_laguerre(N)
    res_w = np.abs(w/2-w_e)
    res_r  = np.abs(r - r_e)
    assert np.max(res_w) == pytest.approx(0)
    assert np.max(res_r) == pytest.approx(0)
