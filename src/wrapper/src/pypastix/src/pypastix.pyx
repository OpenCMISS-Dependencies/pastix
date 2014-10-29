#cython: embedsignature=True

__all__ = ['PASTIX']

# cython imports
from libc.stdlib cimport malloc, free
from libc.stdint cimport *
from posix.fcntl cimport O_RDONLY, O_WRONLY

from cpython cimport PyErr_SetFromErrno
from cpython.string cimport PyString_GET_SIZE, PyString_Check, PyString_AS_STRING

# python imports
import os
import numpy as np
cimport numpy as np
from collections import namedtuple
cimport cython
cdef extern from "mpi-compat.h": pass
cdef extern from "complex.h": pass

cimport mpi4py.MPI as MPI
# changed after 1.3
from mpi4py.mpi_c cimport *

#cimport mpi4py.MPI as MPI
#from mpi4py.mpi_c cimport *

# ...
global _nmatrices
_nmatrices = 0
# ...

# .................................
cdef extern from "mpi.h": pass
# .................................

# .................................

include "_pastix.pxd"
include "_murge.pxd"
include "_cscd_utils.pxd"
include "PASTIX.pxd"
include "MURGE.pxi"
include "MURGE.pyx"
include "Matrix.pyx"
