#-- imports --------------------------------------------------------------------
import cuncsd as csd
import numpy as np

from numpy.linalg import eigh
from numpy.testing import assert_almost_equal

formatter = dict(precision=3, max_line_width=200, suppress_small=True)

#-- generate unitary matrix ----------------------------------------------------
num_modes = 4

H = np.diag([1.0] * (2*num_modes)) + np.diag([0.5] * (2*num_modes -1), k=1) + np.diag([0.5] * (2*num_modes - 1), k=-1)
D, U = eigh(H)

print("U.shape = ", U.shape)

#-- print signature of cuncsd --------------------------------------------------
print(csd.__doc__)

#-- call csd -------------------------------------------------------------------
m = 2*num_modes
p = num_modes
q = num_modes
x11 = U[:num_modes, :num_modes]
ldx11 = x11.shape[0]
x12 = U[:num_modes, num_modes:]
ldx12 = x12.shape[0]
x21 = U[num_modes:, :num_modes]
ldx21 = x21.shape[0]
x22 = U[num_modes:, num_modes:]
ldx22 = x22.shape[0]

print("x11.shape: ", x11.shape)
print("x12.shape: ", x12.shape)
print("x21.shape: ", x21.shape)
print("x22.shape: ", x22.shape)

ldu1 = num_modes
ldu2 = num_modes
ldv1t = num_modes
ldv2t = num_modes

lwork  = -1
lrwork = -1

#-- query work space dimensions ------------------------------------------------
_,_,_,_,_,_,_,_,_,work,rwork,_,info = \
csd.cuncsd( m, p, q, 
            x11, ldx11, 
            x12, ldx12,
            x21, ldx21,
            x22, ldx22,
            ldu1, ldu2,
            ldv1t, ldv2t,
            lwork, lrwork,
            jobu1='Y',
            jobu2='Y',
            jobv1t='Y',
            jobv2t='Y',
            trans='T',
            signs='O',
            credit=0)

print("exit status: ", info)

print("work[0]: ", work[0])
print("rwork:[0] ", rwork[0])

lwork  = work[0]
lrwork = rwork[0]

#-- compute csd ----------------------------------------------------------------
x11,x12,x21,x22,theta,u1,u2,v1t,v2t,work,rwork,iwork,info = \
csd.cuncsd( m, p, q, 
            x11, ldx11, 
            x12, ldx12,
            x21, ldx21,
            x22, ldx22,
            ldu1, ldu2,
            ldv1t, ldv2t,
            lwork, lrwork,
            jobu1='Y',
            jobu2='Y',
            jobv1t='Y',
            jobv2t='Y',
            trans='T',
            signs='O',
            credit=0)

print("exit status: ", info)

print("work: ", work)
print("rwork: ", rwork)

# print("U1 = \n{}\n".format(np.array_str(u1, **formatter)))
# print("U2 = \n{}\n".format(np.array_str(u2, **formatter)))

#-- check unitarity ------------------------------------------------------------
decimals = 6
assert_almost_equal(u1 @ u1.conj().T, np.eye(num_modes), decimals)
assert_almost_equal(u2 @ u2.conj().T, np.eye(num_modes), decimals)

assert_almost_equal(v1t @ v1t.conj().T, np.eye(num_modes), decimals)
assert_almost_equal(v2t @ v2t.conj().T, np.eye(num_modes), decimals)

#-- print results --------------------------------------------------------------
print(np.array_str(theta, **formatter))