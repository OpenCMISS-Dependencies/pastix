from mpi4py import MPI
import argparse
import pypastix
import numpy as np
import re
from scipy.sparse import lil_matrix
from scipy.sparse import tril
from scipy.sparse import csc_matrix
from scipy.io import mminfo, mmread

from numpy import linalg as LA
# Build a simple n x n laplacian and solve
#   2  -1        X   1
#  -1   2 -1     X = 0
#      -1  2 -1  X   0
#         -1  2  X   1


# Argument parsing
parser = argparse.ArgumentParser(description=u'Solves Ax=b using PaStiX from pypastix with A an 1D laplacian matrix.')
parser.add_argument(u'--size', dest=u'n', type=int, default=4,
                    help='Size of th problem')
parser.add_argument(u'--type', dest=u'type', type=str, default=u'float',
                     help=u'arithemtic : [float, double, complex, dcomplex]')
parser.add_argument(u'--sym', dest=u'sym', action='store_true', default=False, help=u"Tells if the matrix is symmetric")
parser.add_argument(u'--her', dest=u'her', action='store_true', default=False, help=u"Tells if the matrix is hermitian")
parser.add_argument(u'--iparm', dest=u'iparm', nargs=2, metavar=('key', 'value'), action='append', help=u"PaStiX integer parameters", default=[])
parser.add_argument(u'--dparm', dest=u'dparm', nargs=2, metavar=('key', 'value'), action='append', help=u"PaStiX double parameters", default=[])
parser.add_argument(u"--matrix", dest=u'matrix', help=u"Path to the MatrixMarket file", default=None)

args = parser.parse_args()
if args.her:
    args.sym = True

n        = args.n
rhs_nbr  = 1


if args.type in [u"float", u"double"]:
    if args.type == u'float':
        dtype = np.float32
    else:
        dtype = np.double
    diag = 2.0
    subdiag = -1.0
    supdiag = -1.0
    rhs_mid=0.0
    rhs_1=1.0
    rhs_n=1.0
elif args.type in [u'complex', u'dcomplex'] :
    if args.type == u'complex':
        dtype = np.complex64
    else:
        dtype = np.complex128
    diag = 2.0
    subdiag = -1.0+2.j
    if args.her:
        supdiag = -1.0-2.j
        rhs_mid=0.0+ 0.j
        rhs_1=-1.0+3.j
    else:
        supdiag = -1.+2.j
        rhs_mid=-4.0+ 4.j
        rhs_1=3.0-1.j
    rhs_n=-1. + 3.j
else:
    parser.print_help()
    raise


#definition of the sparse matrix
if args.matrix:
    (n, m, nnz, form, arith, issym) = mminfo(args.matrix)
    n = int(n)
    A = mmread(args.matrix)
    rhs = [1] * n
    rhs = A*rhs
    dtype=np.dtype(A)

    if issym == u"symmetric":
        args.sym = True
    elif issym == u"hermitian":
        args.her = True
    else:
        args.sym = False
else:
    A = lil_matrix((n, n), dtype=dtype)
    A[0,0] =  diag
    A[1,0] =  subdiag
    rhs=np.ndarray(shape=(n), dtype=dtype)
    rhs[0]=rhs_1

    for i in range(1, n-1):
        A[i-1, i] = supdiag
        A[i, i]   = diag
        A[i+1, i] = subdiag
        rhs[i]=rhs_mid

    A[n-2, n-1] = supdiag
    A[n-1, n-1] = diag
    rhs[n-1]=rhs_n





# initialisation of the solver
comm = MPI.COMM_WORLD
nprocs = comm.Get_size()
myrank = comm.Get_rank()
if dtype == np.float32:
    pastix = pypastix.SPASTIX(comm)
elif dtype == np.float64:
    pastix = pypastix.DPASTIX(comm)
elif dtype == np.complex64:
    pastix = pypastix.CPASTIX(comm)
elif dtype == np.complex128:
    pastix = pypastix.ZPASTIX(comm)

pastix.setIparm(pastix.API.IPARM_START_TASK, pastix.API.API_TASK_ORDERING)
pastix.setIparm(pastix.API.IPARM_END_TASK, pastix.API.API_TASK_CLEAN)
pastix.setIparm(pastix.API.IPARM_VERBOSE, pastix.API.API_VERBOSE_UNBEARABLE)

if args.sym:
    if args.her:
        # LLT == LLH in complex
        pastix.setIparm(pastix.API.IPARM_SYM, pastix.API.API_SYM_HER)
        pastix.setIparm(pastix.API.IPARM_FACTORIZATION, pastix.API.API_FACT_LDLH)
    else:
        # LDLT != LDLH
        pastix.setIparm(pastix.API.IPARM_SYM, pastix.API.API_SYM_YES)
        pastix.setIparm(pastix.API.IPARM_FACTORIZATION, pastix.API.API_FACT_LDLT)
else:
    pastix.setIparm(pastix.API.IPARM_SYM, pastix.API.API_SYM_NO)
    pastix.setIparm(pastix.API.IPARM_FACTORIZATION, pastix.API.API_FACT_LU)
pastix.setIparm(pastix.API.IPARM_PRODUCE_STATS, pastix.API.API_YES)


for iparm in args.iparm:
    (keyword, value) = iparm
    if re.match(u"\d", value):
        value = int(value)
    else:
        value = eval(u"pastix.API."+value)

    if re.match(u"\d", keyword):
        keyword = int(keyword)
    else:
        keyword = eval(u"pastix.API."+keyword)

    pastix.setIparm(keyword, value)

for dparm in args.dparm:
    (keyword, value) = dparm
    if re.match(u"\d", keyword):
        keyword = int(keyword)
    else:
        keyword = eval(u"pastix.API."+keyword)
    value = float(value)
    pastix.setDparm(keyword, value)


if args.sym:
    cscA = tril(A, format='csc')
else:
    cscA = A.tocsc()
b=rhs.copy()
colptr = cscA.indptr
rows = cscA.indices
values = cscA.data

perm     = [x+1 for x in range(n)]
invp     = [x+1 for x in range(n)]
pastix.checkMatrix(n, pastix.integerType(colptr+1), pastix.integerType(rows+1), values, sym = pastix.getIparm(pastix.API.IPARM_SYM))
pastix.pastix(n, pastix.integerType(colptr+1),
              pastix.integerType(rows+1),
              values,
              np.array(perm,
                       dtype=pastix.integerType),
              np.array(invp,
                       dtype=pastix.integerType),
              rhs, rhs_nbr)

if myrank == 0:
    print u"Number of iterations : %d " % pastix.getIparm(pastix.API.IPARM_NBITER)
    print u"Relative error       : %g" % pastix.getDparm(pastix.API.DPARM_RELATIVE_ERROR)
    print u"Scaled residual      : %g" % pastix.getDparm(pastix.API.DPARM_SCALED_RESIDUAL)

    print u"||Ax-b||/||b|| = ", LA.norm(A*rhs-b, np.inf)/(LA.norm(b, np.inf))
    print u"max_i(|Ax-b|_i/(|A||x| + |b|)_i)", max(abs(A*rhs-b)/(abs(A)*abs(rhs)+abs(b)))

    if (n < 20):
        print "A", A.todense()
        print "solution : ", rhs
        print "Ax :", A*rhs
        print "b :", b
