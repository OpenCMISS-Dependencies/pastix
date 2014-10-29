cdef extern from "pastix.h":
    ctypedef int32_t pastix_int_t
ctypedef np.int32_t np_pastix_int_t
ctypedef np.double_t  np_pastix_float_t

cdef extern from "murge.h":
    ctypedef double COEF
    ctypedef double REAL
    ctypedef int INTS
    ctypedef int INTL
ctypedef np.int32_t np_INTS
ctypedef np.int32_t np_INTL
ctypedef np.double_t  np_COEF
ctypedef np.double_t  np_REAL
