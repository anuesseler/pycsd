#-- global imports -------------------------------------------------------------
import numpy  as np
import pytest as pt

from numpy.linalg  import eigh
from numpy.testing import assert_almost_equal

#-- local imports --------------------------------------------------------------
from pycsd import cs_decomp

decimals  = 14

@pt.mark.parametrize('M, P, Q',
    [( 10,  5,  5),
     (100, 50, 50)])
def test_cs_decomp(M, P, Q):

    if P != Q:
        raise NotImplementedError

    #-- generate unitary matrix ------------------------------------------------
    H    = np.random.rand(M,M) + 1.9j * np.random.rand(M,M)
    H    = H + H.conj().T
    D, U = eigh(H)

    print("D = \n{}\n".format(np.array_str(D)))

    #-- csd routine ------------------------------------------------------------
    u1, u2, v1h, v2h, theta = cs_decomp(U, P, Q)


    #-- check correctness of decomposition -------------------------------------
    UD = np.vstack((np.hstack((u1,np.zeros((P, M-P)))),
                    np.hstack((np.zeros((M-P, P)),u2))))

    VDH = np.vstack((np.hstack((v1h,np.zeros((Q, M-Q)))),
                     np.hstack((np.zeros((M-Q, Q)),v2h))))

    C  = np.diag(np.cos(theta))
    S  = np.diag(np.sin(theta))
    CS = np.vstack((np.hstack((C, -S)), np.hstack((S,C))))

    assert_almost_equal(U, UD @ CS @ VDH, decimals)