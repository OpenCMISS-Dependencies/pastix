import pastix
from mpi4py import MPI
# Build a simple n x n laplacian and solve
#   2  -1        X   1
#  -1   2 -1     X = 0
#      -1  2 -1  X   0
#         -1  2  X   1

#print dir(pastix)

pdata    = 0
mpi_comm = MPI.COMM_WORLD
n        = 4
colptr   = pastix.new_pIntArray(n+1)
rows     = pastix.new_pIntArray(3*n-2)
values   = pastix.new_pFloatArray(3*n-2);
perm     = pastix.new_pIntArray(n)
invp     = pastix.new_pIntArray(n)
rhs      = pastix.new_pFloatArray(n)
rhs_nbr  = 1;
iparm    = pastix.new_pIntArray(64)
dparm    = pastix.new_doubleArray(64)

pastix.pIntArray_setitem(iparm, pastix.IPARM_MODIFY_PARAMETER, pastix.API_NO);
pdata = pastix.pastix(pdata, MPI.COMM_WORLD, 
                      n, colptr, rows, values, 
                      perm, invp, 
                      rhs , rhs_nbr, iparm, dparm)

pastix.pIntArray_setitem(iparm, pastix.IPARM_START_TASK, 
                         pastix.API_TASK_ORDERING);
pastix.pIntArray_setitem(iparm, pastix.IPARM_END_TASK,   
                         pastix.API_TASK_CLEAN);
pastix.pIntArray_setitem(iparm, pastix.IPARM_VERBOSE,    
                         pastix.API_VERBOSE_UNBEARABLE);
pastix.pIntArray_setitem(iparm, pastix.IPARM_SYM,    
                         pastix.API_SYM_NO);
pastix.pIntArray_setitem(iparm, pastix.IPARM_FACTORIZATION,    
                         pastix.API_FACT_LU);


pastix.pIntArray_setitem(colptr, 0, 1)
pastix.pIntArray_setitem(rows, 0, 1)
pastix.pIntArray_setitem(rows, 1, 2)
pastix.pFloatArray_setitem(values, 0, 2)
pastix.pFloatArray_setitem(values, 1, -1)
pastix.pFloatArray_setitem(rhs, 0, 1)

for i in range(1, n-1):
    pastix.pIntArray_setitem(colptr, i, i*3)
    pastix.pIntArray_setitem(rows, i*3-1, i)
    pastix.pIntArray_setitem(rows, i*3, i+1)
    pastix.pIntArray_setitem(rows, i*3+1, i+2)
    pastix.pFloatArray_setitem(values, i*3-1, -1)
    pastix.pFloatArray_setitem(values, i*3, 2)
    pastix.pFloatArray_setitem(values, i*3+1, -1)
    pastix.pFloatArray_setitem(rhs, i, 0)

pastix.pIntArray_setitem(colptr, n-1, 3*(n-1))
pastix.pIntArray_setitem(rows, (n-1)*3-1, n-1)
pastix.pIntArray_setitem(rows, (n-1)*3, n)
pastix.pFloatArray_setitem(values, (n-1)*3-1, -1)
pastix.pFloatArray_setitem(values, (n-1)*3, 2)
pastix.pFloatArray_setitem(rhs, n-1, 1)

pastix.pIntArray_setitem(colptr, n, 3*n-1)

pdata = pastix.pastix(pdata, mpi_comm, n, colptr, rows, values, perm, 
                      invp, rhs , rhs_nbr, iparm, dparm)

for i in range(0, n):
    print "rhs[%d] = %g" %(i,pastix.pFloatArray_getitem(rhs, i))

pastix.delete_pIntArray(colptr)
pastix.delete_pIntArray(rows)
pastix.delete_pFloatArray(values)
pastix.delete_pIntArray(perm)
pastix.delete_pIntArray(invp)
pastix.delete_pFloatArray(rhs)
pastix.delete_pIntArray(iparm)
pastix.delete_doubleArray(dparm)
