#-- global imports -------------------------------------------------------------
import numpy as np

from numpy.linalg  import eigh
from numpy.testing import assert_almost_equal

#-- local imports --------------------------------------------------------------
from builds import zuncsd as csd

decimals  = 15
formatter = dict(precision=3, max_line_width=200, suppress_small=True)

def cs_decomp(  X, p, q, 
                comp_u1=True, comp_u2=True,
                comp_v1h=True, comp_v2h=True,
                sign_conv=''):
    """
        Computes cosine-sine (CS) decomposition of an M-by-M partitioned unitary
        matrix X, s.t.

                                         [  I  0  0 |  0  0  0 ]
                                         [  0  C  0 |  0 -S  0 ]
             [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**H
         X = [-----------] = [---------] [---------------------] [---------]   .
             [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
                                         [  0  S  0 |  0  C  0 ]
                                         [  0  0  I |  0  0  0 ]


        :param X:           square unitary to decompose
        :param p:           number of rows in upper left block
        :param q:           number of columns in the upper left block
        :param comp_u1:     compute unitary u1
        :param comp_u2:     compute unitary u2
        :param comp_v1h:    compute unitary v1^\dagger
        :param comp_v2h:    compute unitary v2^\dagger
        :param sign_conv:   'low'       minus sign is put in lower left block
                            otherwise:  minus sign is put in upper right block

        :returns u1:        (optional) upper left unitary block u1     
        :returns u2:        (optional) lower right unitary block u2
        :returns v1h:       (optional) upper left unitary block v1^\dagger
        :returns v2h:       (optional) lower right unitary block v2^\dagger
        :returns theta:     angles in radians to construct 
                            C = \diag(\cos(theta_1), \dots, \cos(theta_2))
                            S = \diag(\sin(theta_1), \dots, \sin(theta_2))

    """

    #-- split X and set lapack parameters --------------------------------------
    m = X.shape[0]

    # upper left
    x11   = X[:p,:q]
    ldx11 = p

    # upper right
    x12   = X[:p,q:]
    ldx12 = p

    # lower left
    x21   = X[p:,:q]
    ldx21 = m - p

    # lower right
    x22   = X[p:,q:]
    ldx22 = m - p

    ldu1  = p
    ldu2  = m - p

    ldv1t = q
    ldv2t = m - q

    #-- convert input parameters -----------------------------------------------
    jobu1  = 'Y' if comp_u1 else ''
    jobu2  = 'Y' if comp_u2 else ''
    jobv1t = 'Y' if comp_v1h else ''
    jobv2t = 'Y' if comp_v2h else ''
    signs  = 'O' if sign_conv == 'low' else ''

    #-- query work space dimensions --------------------------------------------
    lwork  = -1
    lrwork = -1

    _,_,_,_,_,_,_,_,_,work,rwork,_,info = \
    csd.zuncsd( m, p, q, 
                x11, ldx11, 
                x12, ldx12,
                x21, ldx21,
                x22, ldx22,
                ldu1, ldu2,
                ldv1t, ldv2t,
                lwork, lrwork,
                jobu1=jobu1,
                jobu2=jobu2,
                jobv1t=jobv1t,
                jobv2t=jobv2t,
                trans='',
                signs=signs)

    lwork  = work[0]
    lrwork = rwork[0]

    #-- compute csd ------------------------------------------------------------
    _,_,_,_,theta,u1,u2,v1h,v2h,_,_,_,info = \
    csd.zuncsd( m, p, q, 
                x11, ldx11, 
                x12, ldx12,
                x21, ldx21,
                x22, ldx22,
                ldu1, ldu2,
                ldv1t, ldv2t,
                lwork, lrwork,
                jobu1=jobu1,
                jobu2=jobu2,
                jobv1t=jobv1t,
                jobv2t=jobv2t,
                trans='',
                signs=signs)

    return u1, u2, v1h, v2h, theta

if __name__ == '__main__':
    

    #-- generate unitary matrix ------------------------------------------------
    num_modes = 4

    H = np.diag([1.0+0.0j] * (2*num_modes))\
      + np.diag([0.5+0.0j] * (2*num_modes -1), k=1)\
      + np.diag([0.5+0.0j] * (2*num_modes - 1), k=-1)
    D, U = eigh(H)


    #-- csd routine ------------------------------------------------------------
    u1, u2, v1h, v2h, theta = cs_decomp(U, num_modes, num_modes)


    #-- check correctness of decomposition -------------------------------------
    zero = np.zeros((num_modes, num_modes))
    UD = np.vstack((np.hstack((u1,zero)), np.hstack((zero,u2))))
    print("UD = \n{}\n".format(np.array_str(UD, **formatter)))

    VDH = np.vstack((np.hstack((v1h,zero)), np.hstack((zero,v2h))))
    print("VDH = \n{}\n".format(np.array_str(VDH, **formatter)))

    C  = np.diag(np.cos(theta))
    S  = np.diag(np.sin(theta))
    CS = np.vstack((np.hstack((C, -S)), np.hstack((S,C))))
    print("CS = \n{}\n".format(np.array_str(CS, **formatter)))

    print("U = \n{}\n".format(np.array_str(U, **formatter)))
    print("UD @ CS @ VDH = \n{}\n".format(np.array_str(UD @ CS @ VDH, **formatter)))

    assert_almost_equal(U, UD @ CS @ VDH, decimals)




