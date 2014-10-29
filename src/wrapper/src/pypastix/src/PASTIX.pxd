include "pastix_types.pxi"
class ENUM_API:
     MOD_MURGE            = 1200
     MOD_COMMON           = 600
     MOD_FAX              = 400
     MOD_KASS             = 1000
     MOD_BUBBLE           = 1100
     MOD_SCOTCH           = 300
     MOD_GRAPH            = 800
     MOD_SOPALIN          = 100
     MOD_UNKNOWN          = 0
     MOD_SYMBOL           = 900
     MOD_ORDER            = 500
     MOD_SI               = 700
     MOD_BLEND            = 200

     API_SYM_NO           = 1
     API_SYM_HER          = 2
     API_SYM_YES          = 0

     API_TASK_CLEAN       = 7
     API_TASK_REFINE      = 6
     API_TASK_ORDERING    = 1
     API_TASK_INIT        = 0
     API_TASK_SOLVE       = 5
     API_TASK_NUMFACT     = 4
     API_TASK_SYMBFACT    = 2
     API_TASK_ANALYSE     = 3

     API_IO_NO            = 0
     API_IO_LOAD          = 1
     API_IO_SAVE_GRAPH    = 8
     API_IO_LOAD_GRAPH    = 4
     API_IO_SAVE_CSC      = 32
     API_IO_SAVE          = 2
     API_IO_LOAD_CSC      = 16

     API_NO               = 0
     API_YES              = 1

     API_ORDER_METIS      = 1
     API_ORDER_SCOTCH     = 0
     API_ORDER_LOAD       = 3
     API_ORDER_PERSONAL   = 2

     API_TASK_BLEND       = 3
     API_TASK_FAX         = 2
     API_TASK_SOPALIN     = 4
     API_TASK_SCOTCH      = 1
     API_TASK_REFINEMENT  = 6
     API_TASK_UPDOWN      = 5

     API_RHS_0            = 3
     API_RHS_1            = 1
     API_RHS_B            = 0
     API_RHS_I            = 2

     API_VERBOSE_NO       = 1
     API_VERBOSE_CHATTERBOX = 3
     API_VERBOSE_UNBEARABLE = 4
     API_VERBOSE_YES      = 2
     API_VERBOSE_NOT      = 0

     DPARM_FACT_FLOPS     = 22
     DPARM_EPSILON_MAGN_CTRL = 10
     DPARM_MEM_MAX        = 2
     DPARM_ANALYZE_TIME   = 18
     DPARM_RELATIVE_ERROR = 6
     DPARM_PRED_FACT_TIME = 19
     DPARM_SOLV_TIME      = 21
     DPARM_RAFF_TIME      = 24
     DPARM_EPSILON_REFINEMENT = 5
     DPARM_FILL_IN        = 1
     DPARM_SCALED_RESIDUAL = 7
     DPARM_SIZE           = 64
     DPARM_FACT_TIME      = 20
     DPARM_SOLV_FLOPS     = 23

     API_BIND_TAB         = 2
     API_BIND_NO          = 0
     API_BIND_AUTO        = 1

     API_RAF_GRAD         = 1
     API_RAF_GMRES        = 0
     API_RAF_BCGSTAB      = 3
     API_RAF_PIVOT        = 2

     API_COMPLEXSINGLE    = 2
     API_REALSINGLE       = 0
     API_REALDOUBLE       = 1
     API_COMPLEXDOUBLE    = 3

     BAD_DEFINE_ERR       = 10
     FILE_ERR             = 9
     NOTIMPLEMENTED_ERR   = 4
     INTERNAL_ERR         = 7
     IO_ERR               = 12
     BADPARAMETER_ERR     = 8
     FLOAT_TYPE_ERR       = 14
     UNKNOWN_ERR          = 1
     INTEGER_TYPE_ERR     = 11
     ALLOC_ERR            = 2
     NO_ERR               = 0
     ASSERT_ERR           = 3
     MATRIX_ERR           = 13
     THREAD_ERR           = 6
     STEP_ORDER_ERR       = 15
     OUTOFMEMORY_ERR      = 5
     MPI_ERR              = 16

     IPARM_AMALGAMATION_LEVEL = 13
     IPARM_RHS_MAKING     = 38
     IPARM_TRANSPOSE_SOLVE = 65
     IPARM_ESP            = 43
     IPARM_NNZEROS_BLOCK_LOCAL = 31
     IPARM_ORDERING_SWITCH_LEVEL = 16
     IPARM_ORDERING_CMIN  = 17
     IPARM_ORDERING_CMAX  = 18
     IPARM_FACTORIZATION  = 30
     IPARM_PRODUCE_STATS  = 68
     IPARM_OOC_THREAD     = 48
     IPARM_ITERMAX        = 5
     IPARM_ORDERING       = 14
     IPARM_ERROR_NUMBER   = 63
     IPARM_BASEVAL        = 24
     IPARM_ESP_NBTASKS    = 55
     IPARM_STATIC_PIVOTING = 20
     IPARM_DOF_COST       = 57
     IPARM_FLOAT          = 61
     IPARM_ALLOCATED_TERMS = 23
     IPARM_START_TASK     = 1
     IPARM_FREE_CSCUSER   = 45
     IPARM_VERBOSE        = 3
     IPARM_GMRES_IM       = 44
     IPARM_CSCD_CORRECT   = 9
     IPARM_MIN_BLOCKSIZE  = 25
     IPARM_ORDERING_FRAT  = 19
     IPARM_ABS            = 42
     IPARM_ESP_THRESHOLD  = 56
     IPARM_NNZEROS        = 22
     IPARM_AUTOSPLIT_COMM = 60
     IPARM_GRAPHDIST      = 12
     IPARM_CPU_BY_NODE    = 32
     IPARM_IO_STRATEGY    = 37
     IPARM_TRACEFMT       = 11
     IPARM_REFINEMENT     = 39
     IPARM_MC64           = 7
     IPARM_INCOMPLETE     = 41
     IPARM_MATRIX_VERIFICATION = 6
     IPARM_OOC_ID         = 49
     IPARM_ONLY_RAFF      = 8
     IPARM_PID            = 62
     IPARM_SYM            = 40
     IPARM_OOC_LIMIT      = 47
     IPARM_SCHUR          = 27
     IPARM_METIS_PFACTOR  = 21
     IPARM_NB_SMP_NODE_USED = 50
     IPARM_MURGE_REFINEMENT = 58
     IPARM_THREAD_COMM_MODE = 51
     IPARM_RHSD_CHECK     = 29
     IPARM_END_TASK       = 2
     IPARM_ISOLATE_ZEROS  = 28
     IPARM_FREE_CSCPASTIX = 46
     IPARM_NB_THREAD_COMM = 52
     IPARM_CUDA_NBR       = 64
     IPARM_MAX_BLOCKSIZE  = 26
     IPARM_NBITER         = 10
     IPARM_THREAD_NBR     = 34
     IPARM_DISTRIBUTION_LEVEL = 35
     IPARM_INERTIA        = 54
     IPARM_DOF_NBR        = 4
     IPARM_MODIFY_PARAMETER = 0
     IPARM_STARPU_CTX_NBR = 67
     IPARM_BINDTHRD       = 33
     IPARM_STARPU_CTX_DEPTH = 66
     IPARM_FILL_MATRIX    = 53
     IPARM_LEVEL_OF_FILL  = 36
     IPARM_SIZE           = 128
     IPARM_DEFAULT_ORDERING = 15
     IPARM_STARPU         = 59

     API_TRACE_PICL       = 0
     API_TRACE_HUMREAD    = 2
     API_TRACE_PAJE       = 1
     API_TRACE_UNFORMATED = 3

     API_THREAD_MULTIPLE  = 1
     API_THREAD_FUNNELED  = 2
     API_THREAD_COMM_NBPROC = 16
     API_THREAD_COMM_ONE  = 4
     API_THREAD_COMM_DEFINED = 8

     API_CSC_PRESERVE     = 0
     API_CSC_FREE         = 1

     API_FACT_LLT         = 0
     API_FACT_LDLH        = 3
     API_FACT_LU          = 2
     API_FACT_LDLT        = 1

from libc.stdlib cimport malloc, free

cimport mpi4py.MPI as MPI
from mpi4py.mpi_c cimport *

# PaStiX Exception
class PastixException(Exception):
    pass

class PastixArithmeticException(PastixException):
    pass


cdef class ZPASTIX:
    API                  = ENUM_API

    cdef:

        cdef pastix_data_t  * pdata

        cdef MPI_Comm         mpi_comm

        cdef pastix_int_t   * iparm

        cdef double         * dparm

        cdef pastix_int_t     n

        cdef pastix_int_t   * colptr

        cdef pastix_int_t   * rows

        cdef double complex * values

        cdef pastix_int_t   * perm

        cdef pastix_int_t   * invp

        cdef double complex * rhs

        cdef pastix_int_t     rhs_nbr



    def __cinit__(self, MPI.Comm comm not None ):

        cdef MPI_Comm c_comm = comm.ob_mpi

        self.pdata = NULL

        self.mpi_comm = c_comm

        self.iparm = <pastix_int_t*>malloc(self.API.IPARM_SIZE*sizeof(pastix_int_t))

        self.dparm = <double*>malloc(self.API.DPARM_SIZE*sizeof(double))

        self.iparm[self.API.IPARM_MODIFY_PARAMETER] = self.API.API_NO;

        z_pastix(&self.pdata, self.mpi_comm,

                 0, NULL, NULL, NULL,

                 NULL, NULL, NULL, 1,

                 <pastix_int_t*>self.iparm, self.dparm);

        print self.dparm[self.API.DPARM_EPSILON_REFINEMENT]



    def __dealloc__(self):

        self.iparm[self.API.IPARM_START_TASK] = self.API.API_TASK_CLEAN;

        self.iparm[self.API.IPARM_END_TASK] = self.API.API_TASK_CLEAN;

        if (self.pdata != <pastix_data_t*>NULL):

            z_pastix(&self.pdata, self.mpi_comm,

                     0, NULL, NULL, NULL,

                     NULL, NULL, NULL, 1,

                     <pastix_int_t*>self.iparm, self.dparm);

        free(self.iparm)

        free(self.dparm)



    def setIparm(self, index, value):

        """ Set value of iparm[index] """

        self.iparm[index] = value

    def setDparm(self, index, value):

        """ Set value of dparm[index] """

        self.dparm[index] = value

    def getIparm(self, index):

        """ Get value of iparm[index] """

        return self.iparm[index]

    def getDparm(self, index):

        """ Get value of dparm[index] """

        return self.dparm[index]





    def pastix(self,

               ncol,

               np.ndarray[np_pastix_int_t,ndim=1] colptr,

               np.ndarray[np_pastix_int_t,ndim=1] rows,

               np.ndarray[np.complex128_t,  ndim=1] values,

               np.ndarray[np_pastix_int_t,ndim=1] perm,

               np.ndarray[np_pastix_int_t,ndim=1] invp,

               np.ndarray[np.complex128_t,  ndim=1] rhs,

               nrhs=1):

        '''

        Generic interface to [sdcz]_pastix.

        Arithmetic depends on the type of values (which should be same ad rhs).

        '''

        self.n = <pastix_int_t>ncol

        self.colptr = <pastix_int_t*>colptr.data

        self.rows = <pastix_int_t*>rows.data

        self.perm = <pastix_int_t*>perm.data

        self.invp = <pastix_int_t*>invp.data

        self.rhs_nbr = nrhs

        print "pouet", self.dparm[self.API.DPARM_EPSILON_REFINEMENT]



        z_pastix(&self.pdata, self.mpi_comm,

                 self.n, self.colptr, self.rows, <double complex*>values.data,

                 self.perm, self.invp,

                 <double complex*>rhs.data, self.rhs_nbr,

                 self.iparm, self.dparm)







    def dpastix(self,

                ncol,

                np.ndarray[np_pastix_int_t,ndim=1] colptr,

                np.ndarray[np_pastix_int_t,ndim=1] rows,

                np.ndarray[np.complex128_t,  ndim=1] values,

                np.ndarray[np_pastix_int_t,ndim=1] l2g,

                np.ndarray[np_pastix_int_t,ndim=1] perm,

                np.ndarray[np_pastix_int_t,ndim=1] invp,

                np.ndarray[np.complex128_t,  ndim=1] rhs,

                nrhs=1):

        '''

        Generic interface to [sdcz]_dpastix (Distributed interface to PaStiX).

        Arithmetic depends on the type of values (which should be same ad rhs).

        '''

        self.iparm[self.API.IPARM_GRAPHDIST] = self.API.API_YES

        self.n = <pastix_int_t>ncol

        self.colptr = <pastix_int_t*>colptr.data

        self.rows = <pastix_int_t*>rows.data

        self.perm = <pastix_int_t*>perm.data

        self.invp = <pastix_int_t*>invp.data

        self.rhs_nbr = nrhs

        z_dpastix(&self.pdata, self.mpi_comm,

                  self.n, self.colptr, self.rows, <double complex*>values.data,

                  <pastix_int_t*>l2g.data,

                  self.perm, self.invp,

                  <double complex*>rhs.data, self.rhs_nbr,

                  self.iparm, self.dparm)







    def checkMatrix(self,

                    n,

                    np.ndarray[np_pastix_int_t,ndim=1] colptr,

                    np.ndarray[np_pastix_int_t,ndim=1] rows,

                    np.ndarray[np.complex128_t,  ndim=1] values,

                    np.ndarray[np_pastix_int_t,ndim=1] loc2glob=None,

                    verbosity=None,

                    sym=None,

                    correction=None,

                    dofnbr=1):

        '''

        Generic interface to [sdcz]_pastix_checkMatrix()

        Arithmetic depends on the type of values.

        '''



        cdef pastix_int_t *   colptr_c = <pastix_int_t*>colptr.data

        cdef pastix_int_t *   rows_c   = <pastix_int_t*>rows.data

        cdef pastix_int_t *   l2g_c    = NULL

        cdef pastix_int_t **  l2g_ptr  = NULL

        cdef double complex * values_c = NULL



        cdef pastix_int_t ret = self.API.NO_ERR



        if sym == None:

            sym = self.API.API_SYM_NO

        if verbosity==None:

            verbosity = self.API.API_VERBOSE_NO

        if correction==None:

            correction = self.API.API_NO

        if not loc2glob == None:

            l2g_c = <pastix_int_t*>loc2glob.data

            l2g_ptr = &l2g_c

        values_c = <double complex *>values.data



        ret = z_pastix_checkMatrix(self.mpi_comm,

                                   <pastix_int_t>verbosity,

                                   <pastix_int_t>sym,

                                   <pastix_int_t>correction,

                                   <pastix_int_t>n,

                                   <pastix_int_t**>&colptr_c,

                                   <pastix_int_t**>&rows_c,

                                   <double complex**>&values_c,

                                   <pastix_int_t**>l2g_ptr,

                                   <pastix_int_t>dofnbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def bindThreads(self,np.ndarray[np_pastix_int_t,ndim=1] bindTab):

        cdef pastix_int_t nthreads_c

        (nthreads_py) = bindTab.shape

        nthreads_c = <pastix_int_t>nthreads_py

        z_pastix_bindThreads(self.pdata,

                             nthreads_c,

                             <pastix_int_t*>bindTab.data)



    def getLocalNodeLst(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist

        nodenbr = z_pastix_getLocalNodeNbr(&self.pdata)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)

        ret = z_pastix_getLocalNodeLst(&self.pdata,

                                       <pastix_int_t*>nodelist.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def getLocalUnknowLst(self):

        cdef pastix_int_t ret

        cdef pastix_int_t unknownnbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] unknownlist

        unknownnbr = z_pastix_getLocalUnknownNbr(&self.pdata)

        unknownlist = np.zeros(unknownnbr, dtype=self.integerType)

        ret = z_pastix_getLocalUnknownLst(&self.pdata,

                                          <pastix_int_t*>unknownlist.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



        return unknownlist





    def setSchurUnknowLst(self, np.ndarray[np_pastix_int_t,ndim=1] listSchur):

        cdef pastix_int_t ret

        cdef pastix_int_t n_c

        (n_py) = listSchur.shape

        n_c = <pastix_int_t>n_py



        ret = z_pastix_setSchurUnknownList(self.pdata,

                                           <pastix_int_t> n_c,

                                           <pastix_int_t *>listSchur.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getSchurLocalNodeList(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist



        ret = z_pastix_getSchurLocalNodeNbr(self.pdata, &nodenbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)



        ret = z_pastix_getSchurLocalNodeList(self.pdata,

                                             <pastix_int_t*>nodelist.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def getSchurLocalUnknownList(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist



        ret = z_pastix_getSchurLocalUnkownNbr(self.pdata, &nodenbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)



        ret = z_pastix_getSchurLocalUnknownList(self.pdata,

                                                <pastix_int_t*>nodelist.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def setSchurArray(self,

                      np.ndarray[np.complex128_t, ndim=1] schur):

        cdef pastix_int_t ret



        ret = z_pastix_setSchurArray(self.pdata,

                                     <double complex *>schur.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getSchur(self,

                 np.ndarray[np.complex128_t, ndim=1] schur):

        cdef pastix_int_t ret



        ret = z_pastix_getSchur(self.pdata,

                                <double complex *>schur.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getMemoryUsage(self):

        cdef long mem

        cdef long maxmem

        mem = z_pastix_getMemoryUsage()

        maxmem = z_pastix_getMaxMemoryUsage()

        return (mem, maxmem)



    def unscale(self, sym=None):

        if sym==None:

            sym = self.API.API_SYM_NO

        z_pastix_unscale(self.pdata,

                         <pastix_int_t> sym)



    property integerType:

            def __get__(self):

                cdef pastix_int_t integer

                integer = 1

                if sizeof(integer) == 4:

                    return np.int32

                elif sizeof(integer) == 8:

                    return np.int64

                else:

                    raise PastixException




cdef class CPASTIX:
    API = ENUM_API

    cdef:

        cdef pastix_data_t  * pdata

        cdef MPI_Comm         mpi_comm

        cdef pastix_int_t   * iparm

        cdef double         * dparm

        cdef pastix_int_t     n

        cdef pastix_int_t   * colptr

        cdef pastix_int_t   * rows

        cdef complex * values

        cdef pastix_int_t   * perm

        cdef pastix_int_t   * invp

        cdef complex * rhs

        cdef pastix_int_t     rhs_nbr



    def __cinit__(self, MPI.Comm comm not None ):

        cdef MPI_Comm c_comm = comm.ob_mpi

        self.pdata = NULL

        self.mpi_comm = c_comm

        self.iparm = <pastix_int_t*>malloc(self.API.IPARM_SIZE*sizeof(pastix_int_t))

        self.dparm = <double*>malloc(self.API.DPARM_SIZE*sizeof(double))

        self.iparm[self.API.IPARM_MODIFY_PARAMETER] = self.API.API_NO;

        c_pastix(&self.pdata, self.mpi_comm,

                 0, NULL, NULL, NULL,

                 NULL, NULL, NULL, 1,

                 <pastix_int_t*>self.iparm, self.dparm);

        print self.dparm[self.API.DPARM_EPSILON_REFINEMENT]



    def __dealloc__(self):

        self.iparm[self.API.IPARM_START_TASK] = self.API.API_TASK_CLEAN;

        self.iparm[self.API.IPARM_END_TASK] = self.API.API_TASK_CLEAN;

        if (self.pdata != <pastix_data_t*>NULL):

            c_pastix(&self.pdata, self.mpi_comm,

                     0, NULL, NULL, NULL,

                     NULL, NULL, NULL, 1,

                     <pastix_int_t*>self.iparm, self.dparm);

        free(self.iparm)

        free(self.dparm)



    def setIparm(self, index, value):

        """ Set value of iparm[index] """

        self.iparm[index] = value

    def setDparm(self, index, value):

        """ Set value of dparm[index] """

        self.dparm[index] = value

    def getIparm(self, index):

        """ Get value of iparm[index] """

        return self.iparm[index]

    def getDparm(self, index):

        """ Get value of dparm[index] """

        return self.dparm[index]





    def pastix(self,

               ncol,

               np.ndarray[np_pastix_int_t,ndim=1] colptr,

               np.ndarray[np_pastix_int_t,ndim=1] rows,

               np.ndarray[np.complex64_t,  ndim=1] values,

               np.ndarray[np_pastix_int_t,ndim=1] perm,

               np.ndarray[np_pastix_int_t,ndim=1] invp,

               np.ndarray[np.complex64_t,  ndim=1] rhs,

               nrhs=1):

        '''

        Generic interface to [sdcz]_pastix.

        Arithmetic depends on the type of values (which should be same ad rhs).

        '''

        self.n = <pastix_int_t>ncol

        self.colptr = <pastix_int_t*>colptr.data

        self.rows = <pastix_int_t*>rows.data

        self.perm = <pastix_int_t*>perm.data

        self.invp = <pastix_int_t*>invp.data

        self.rhs_nbr = nrhs

        print "pouet", self.dparm[self.API.DPARM_EPSILON_REFINEMENT]



        c_pastix(&self.pdata, self.mpi_comm,

                 self.n, self.colptr, self.rows, <complex*>values.data,

                 self.perm, self.invp,

                 <complex*>rhs.data, self.rhs_nbr,

                 self.iparm, self.dparm)







    def dpastix(self,

                ncol,

                np.ndarray[np_pastix_int_t,ndim=1] colptr,

                np.ndarray[np_pastix_int_t,ndim=1] rows,

                np.ndarray[np.complex64_t,  ndim=1] values,

                np.ndarray[np_pastix_int_t,ndim=1] l2g,

                np.ndarray[np_pastix_int_t,ndim=1] perm,

                np.ndarray[np_pastix_int_t,ndim=1] invp,

                np.ndarray[np.complex64_t,  ndim=1] rhs,

                nrhs=1):

        '''

        Generic interface to [sdcz]_dpastix (Distributed interface to PaStiX).

        Arithmetic depends on the type of values (which should be same ad rhs).

        '''

        self.iparm[self.API.IPARM_GRAPHDIST] = self.API.API_YES

        self.n = <pastix_int_t>ncol

        self.colptr = <pastix_int_t*>colptr.data

        self.rows = <pastix_int_t*>rows.data

        self.perm = <pastix_int_t*>perm.data

        self.invp = <pastix_int_t*>invp.data

        self.rhs_nbr = nrhs

        c_dpastix(&self.pdata, self.mpi_comm,

                  self.n, self.colptr, self.rows, <complex*>values.data,

                  <pastix_int_t*>l2g.data,

                  self.perm, self.invp,

                  <complex*>rhs.data, self.rhs_nbr,

                  self.iparm, self.dparm)







    def checkMatrix(self,

                    n,

                    np.ndarray[np_pastix_int_t,ndim=1] colptr,

                    np.ndarray[np_pastix_int_t,ndim=1] rows,

                    np.ndarray[np.complex64_t,  ndim=1] values,

                    np.ndarray[np_pastix_int_t,ndim=1] loc2glob=None,

                    verbosity=None,

                    sym=None,

                    correction=None,

                    dofnbr=1):

        '''

        Generic interface to [sdcz]_pastix_checkMatrix()

        Arithmetic depends on the type of values.

        '''



        cdef pastix_int_t *   colptr_c = <pastix_int_t*>colptr.data

        cdef pastix_int_t *   rows_c   = <pastix_int_t*>rows.data

        cdef pastix_int_t *   l2g_c    = NULL

        cdef pastix_int_t **  l2g_ptr  = NULL

        cdef complex * values_c = NULL



        cdef pastix_int_t ret = self.API.NO_ERR



        if sym == None:

            sym = self.API.API_SYM_NO

        if verbosity==None:

            verbosity = self.API.API_VERBOSE_NO

        if correction==None:

            correction = self.API.API_NO

        if not loc2glob == None:

            l2g_c = <pastix_int_t*>loc2glob.data

            l2g_ptr = &l2g_c

        values_c = <complex *>values.data



        ret = c_pastix_checkMatrix(self.mpi_comm,

                                   <pastix_int_t>verbosity,

                                   <pastix_int_t>sym,

                                   <pastix_int_t>correction,

                                   <pastix_int_t>n,

                                   <pastix_int_t**>&colptr_c,

                                   <pastix_int_t**>&rows_c,

                                   <complex**>&values_c,

                                   <pastix_int_t**>l2g_ptr,

                                   <pastix_int_t>dofnbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def bindThreads(self,np.ndarray[np_pastix_int_t,ndim=1] bindTab):

        cdef pastix_int_t nthreads_c

        (nthreads_py) = bindTab.shape

        nthreads_c = <pastix_int_t>nthreads_py

        c_pastix_bindThreads(self.pdata,

                             nthreads_c,

                             <pastix_int_t*>bindTab.data)



    def getLocalNodeLst(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist

        nodenbr = c_pastix_getLocalNodeNbr(&self.pdata)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)

        ret = c_pastix_getLocalNodeLst(&self.pdata,

                                       <pastix_int_t*>nodelist.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def getLocalUnknowLst(self):

        cdef pastix_int_t ret

        cdef pastix_int_t unknownnbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] unknownlist

        unknownnbr = c_pastix_getLocalUnknownNbr(&self.pdata)

        unknownlist = np.zeros(unknownnbr, dtype=self.integerType)

        ret = c_pastix_getLocalUnknownLst(&self.pdata,

                                          <pastix_int_t*>unknownlist.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



        return unknownlist





    def setSchurUnknowLst(self, np.ndarray[np_pastix_int_t,ndim=1] listSchur):

        cdef pastix_int_t ret

        cdef pastix_int_t n_c

        (n_py) = listSchur.shape

        n_c = <pastix_int_t>n_py



        ret = c_pastix_setSchurUnknownList(self.pdata,

                                           <pastix_int_t> n_c,

                                           <pastix_int_t *>listSchur.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getSchurLocalNodeList(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist



        ret = c_pastix_getSchurLocalNodeNbr(self.pdata, &nodenbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)



        ret = c_pastix_getSchurLocalNodeList(self.pdata,

                                             <pastix_int_t*>nodelist.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def getSchurLocalUnknownList(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist



        ret = c_pastix_getSchurLocalUnkownNbr(self.pdata, &nodenbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)



        ret = c_pastix_getSchurLocalUnknownList(self.pdata,

                                                <pastix_int_t*>nodelist.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def setSchurArray(self,

                      np.ndarray[np.complex64_t, ndim=1] schur):

        cdef pastix_int_t ret



        ret = c_pastix_setSchurArray(self.pdata,

                                     <complex *>schur.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getSchur(self,

                 np.ndarray[np.complex64_t, ndim=1] schur):

        cdef pastix_int_t ret



        ret = c_pastix_getSchur(self.pdata,

                                <complex *>schur.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getMemoryUsage(self):

        cdef long mem

        cdef long maxmem

        mem = c_pastix_getMemoryUsage()

        maxmem = c_pastix_getMaxMemoryUsage()

        return (mem, maxmem)



    def unscale(self, sym=None):

        if sym==None:

            sym = self.API.API_SYM_NO

        c_pastix_unscale(self.pdata,

                         <pastix_int_t> sym)



    property integerType:

            def __get__(self):

                cdef pastix_int_t integer

                integer = 1

                if sizeof(integer) == 4:

                    return np.int32

                elif sizeof(integer) == 8:

                    return np.int64

                else:

                    raise PastixException




cdef class DPASTIX:
    API = ENUM_API

    cdef:

        cdef pastix_data_t  * pdata

        cdef MPI_Comm         mpi_comm

        cdef pastix_int_t   * iparm

        cdef double         * dparm

        cdef pastix_int_t     n

        cdef pastix_int_t   * colptr

        cdef pastix_int_t   * rows

        cdef double * values

        cdef pastix_int_t   * perm

        cdef pastix_int_t   * invp

        cdef double * rhs

        cdef pastix_int_t     rhs_nbr



    def __cinit__(self, MPI.Comm comm not None ):

        cdef MPI_Comm c_comm = comm.ob_mpi

        self.pdata = NULL

        self.mpi_comm = c_comm

        self.iparm = <pastix_int_t*>malloc(self.API.IPARM_SIZE*sizeof(pastix_int_t))

        self.dparm = <double*>malloc(self.API.DPARM_SIZE*sizeof(double))

        self.iparm[self.API.IPARM_MODIFY_PARAMETER] = self.API.API_NO;

        d_pastix(&self.pdata, self.mpi_comm,

                 0, NULL, NULL, NULL,

                 NULL, NULL, NULL, 1,

                 <pastix_int_t*>self.iparm, self.dparm);

        print self.dparm[self.API.DPARM_EPSILON_REFINEMENT]



    def __dealloc__(self):

        self.iparm[self.API.IPARM_START_TASK] = self.API.API_TASK_CLEAN;

        self.iparm[self.API.IPARM_END_TASK] = self.API.API_TASK_CLEAN;

        if (self.pdata != <pastix_data_t*>NULL):

            d_pastix(&self.pdata, self.mpi_comm,

                     0, NULL, NULL, NULL,

                     NULL, NULL, NULL, 1,

                     <pastix_int_t*>self.iparm, self.dparm);

        free(self.iparm)

        free(self.dparm)



    def setIparm(self, index, value):

        """ Set value of iparm[index] """

        self.iparm[index] = value

    def setDparm(self, index, value):

        """ Set value of dparm[index] """

        self.dparm[index] = value

    def getIparm(self, index):

        """ Get value of iparm[index] """

        return self.iparm[index]

    def getDparm(self, index):

        """ Get value of dparm[index] """

        return self.dparm[index]





    def pastix(self,

               ncol,

               np.ndarray[np_pastix_int_t,ndim=1] colptr,

               np.ndarray[np_pastix_int_t,ndim=1] rows,

               np.ndarray[np.float64_t,  ndim=1] values,

               np.ndarray[np_pastix_int_t,ndim=1] perm,

               np.ndarray[np_pastix_int_t,ndim=1] invp,

               np.ndarray[np.float64_t,  ndim=1] rhs,

               nrhs=1):

        '''

        Generic interface to [sdcz]_pastix.

        Arithmetic depends on the type of values (which should be same ad rhs).

        '''

        self.n = <pastix_int_t>ncol

        self.colptr = <pastix_int_t*>colptr.data

        self.rows = <pastix_int_t*>rows.data

        self.perm = <pastix_int_t*>perm.data

        self.invp = <pastix_int_t*>invp.data

        self.rhs_nbr = nrhs

        print "pouet", self.dparm[self.API.DPARM_EPSILON_REFINEMENT]



        d_pastix(&self.pdata, self.mpi_comm,

                 self.n, self.colptr, self.rows, <double*>values.data,

                 self.perm, self.invp,

                 <double*>rhs.data, self.rhs_nbr,

                 self.iparm, self.dparm)







    def dpastix(self,

                ncol,

                np.ndarray[np_pastix_int_t,ndim=1] colptr,

                np.ndarray[np_pastix_int_t,ndim=1] rows,

                np.ndarray[np.float64_t,  ndim=1] values,

                np.ndarray[np_pastix_int_t,ndim=1] l2g,

                np.ndarray[np_pastix_int_t,ndim=1] perm,

                np.ndarray[np_pastix_int_t,ndim=1] invp,

                np.ndarray[np.float64_t,  ndim=1] rhs,

                nrhs=1):

        '''

        Generic interface to [sdcz]_dpastix (Distributed interface to PaStiX).

        Arithmetic depends on the type of values (which should be same ad rhs).

        '''

        self.iparm[self.API.IPARM_GRAPHDIST] = self.API.API_YES

        self.n = <pastix_int_t>ncol

        self.colptr = <pastix_int_t*>colptr.data

        self.rows = <pastix_int_t*>rows.data

        self.perm = <pastix_int_t*>perm.data

        self.invp = <pastix_int_t*>invp.data

        self.rhs_nbr = nrhs

        d_dpastix(&self.pdata, self.mpi_comm,

                  self.n, self.colptr, self.rows, <double*>values.data,

                  <pastix_int_t*>l2g.data,

                  self.perm, self.invp,

                  <double*>rhs.data, self.rhs_nbr,

                  self.iparm, self.dparm)







    def checkMatrix(self,

                    n,

                    np.ndarray[np_pastix_int_t,ndim=1] colptr,

                    np.ndarray[np_pastix_int_t,ndim=1] rows,

                    np.ndarray[np.float64_t,  ndim=1] values,

                    np.ndarray[np_pastix_int_t,ndim=1] loc2glob=None,

                    verbosity=None,

                    sym=None,

                    correction=None,

                    dofnbr=1):

        '''

        Generic interface to [sdcz]_pastix_checkMatrix()

        Arithmetic depends on the type of values.

        '''



        cdef pastix_int_t *   colptr_c = <pastix_int_t*>colptr.data

        cdef pastix_int_t *   rows_c   = <pastix_int_t*>rows.data

        cdef pastix_int_t *   l2g_c    = NULL

        cdef pastix_int_t **  l2g_ptr  = NULL

        cdef double * values_c = NULL



        cdef pastix_int_t ret = self.API.NO_ERR



        if sym == None:

            sym = self.API.API_SYM_NO

        if verbosity==None:

            verbosity = self.API.API_VERBOSE_NO

        if correction==None:

            correction = self.API.API_NO

        if not loc2glob == None:

            l2g_c = <pastix_int_t*>loc2glob.data

            l2g_ptr = &l2g_c

        values_c = <double *>values.data



        ret = d_pastix_checkMatrix(self.mpi_comm,

                                   <pastix_int_t>verbosity,

                                   <pastix_int_t>sym,

                                   <pastix_int_t>correction,

                                   <pastix_int_t>n,

                                   <pastix_int_t**>&colptr_c,

                                   <pastix_int_t**>&rows_c,

                                   <double**>&values_c,

                                   <pastix_int_t**>l2g_ptr,

                                   <pastix_int_t>dofnbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def bindThreads(self,np.ndarray[np_pastix_int_t,ndim=1] bindTab):

        cdef pastix_int_t nthreads_c

        (nthreads_py) = bindTab.shape

        nthreads_c = <pastix_int_t>nthreads_py

        d_pastix_bindThreads(self.pdata,

                             nthreads_c,

                             <pastix_int_t*>bindTab.data)



    def getLocalNodeLst(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist

        nodenbr = d_pastix_getLocalNodeNbr(&self.pdata)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)

        ret = d_pastix_getLocalNodeLst(&self.pdata,

                                       <pastix_int_t*>nodelist.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def getLocalUnknowLst(self):

        cdef pastix_int_t ret

        cdef pastix_int_t unknownnbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] unknownlist

        unknownnbr = d_pastix_getLocalUnknownNbr(&self.pdata)

        unknownlist = np.zeros(unknownnbr, dtype=self.integerType)

        ret = d_pastix_getLocalUnknownLst(&self.pdata,

                                          <pastix_int_t*>unknownlist.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



        return unknownlist





    def setSchurUnknowLst(self, np.ndarray[np_pastix_int_t,ndim=1] listSchur):

        cdef pastix_int_t ret

        cdef pastix_int_t n_c

        (n_py) = listSchur.shape

        n_c = <pastix_int_t>n_py



        ret = d_pastix_setSchurUnknownList(self.pdata,

                                           <pastix_int_t> n_c,

                                           <pastix_int_t *>listSchur.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getSchurLocalNodeList(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist



        ret = d_pastix_getSchurLocalNodeNbr(self.pdata, &nodenbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)



        ret = d_pastix_getSchurLocalNodeList(self.pdata,

                                             <pastix_int_t*>nodelist.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def getSchurLocalUnknownList(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist



        ret = d_pastix_getSchurLocalUnkownNbr(self.pdata, &nodenbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)



        ret = d_pastix_getSchurLocalUnknownList(self.pdata,

                                                <pastix_int_t*>nodelist.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def setSchurArray(self,

                      np.ndarray[np.float64_t, ndim=1] schur):

        cdef pastix_int_t ret



        ret = d_pastix_setSchurArray(self.pdata,

                                     <double *>schur.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getSchur(self,

                 np.ndarray[np.float64_t, ndim=1] schur):

        cdef pastix_int_t ret



        ret = d_pastix_getSchur(self.pdata,

                                <double *>schur.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getMemoryUsage(self):

        cdef long mem

        cdef long maxmem

        mem = d_pastix_getMemoryUsage()

        maxmem = d_pastix_getMaxMemoryUsage()

        return (mem, maxmem)



    def unscale(self, sym=None):

        if sym==None:

            sym = self.API.API_SYM_NO

        d_pastix_unscale(self.pdata,

                         <pastix_int_t> sym)



    property integerType:

            def __get__(self):

                cdef pastix_int_t integer

                integer = 1

                if sizeof(integer) == 4:

                    return np.int32

                elif sizeof(integer) == 8:

                    return np.int64

                else:

                    raise PastixException




cdef class SPASTIX:
    API = ENUM_API

    cdef:

        cdef pastix_data_t  * pdata

        cdef MPI_Comm         mpi_comm

        cdef pastix_int_t   * iparm

        cdef double         * dparm

        cdef pastix_int_t     n

        cdef pastix_int_t   * colptr

        cdef pastix_int_t   * rows

        cdef float * values

        cdef pastix_int_t   * perm

        cdef pastix_int_t   * invp

        cdef float * rhs

        cdef pastix_int_t     rhs_nbr



    def __cinit__(self, MPI.Comm comm not None ):

        cdef MPI_Comm c_comm = comm.ob_mpi

        self.pdata = NULL

        self.mpi_comm = c_comm

        self.iparm = <pastix_int_t*>malloc(self.API.IPARM_SIZE*sizeof(pastix_int_t))

        self.dparm = <double*>malloc(self.API.DPARM_SIZE*sizeof(double))

        self.iparm[self.API.IPARM_MODIFY_PARAMETER] = self.API.API_NO;

        s_pastix(&self.pdata, self.mpi_comm,

                 0, NULL, NULL, NULL,

                 NULL, NULL, NULL, 1,

                 <pastix_int_t*>self.iparm, self.dparm);

        print self.dparm[self.API.DPARM_EPSILON_REFINEMENT]



    def __dealloc__(self):

        self.iparm[self.API.IPARM_START_TASK] = self.API.API_TASK_CLEAN;

        self.iparm[self.API.IPARM_END_TASK] = self.API.API_TASK_CLEAN;

        if (self.pdata != <pastix_data_t*>NULL):

            s_pastix(&self.pdata, self.mpi_comm,

                     0, NULL, NULL, NULL,

                     NULL, NULL, NULL, 1,

                     <pastix_int_t*>self.iparm, self.dparm);

        free(self.iparm)

        free(self.dparm)



    def setIparm(self, index, value):

        """ Set value of iparm[index] """

        self.iparm[index] = value

    def setDparm(self, index, value):

        """ Set value of dparm[index] """

        self.dparm[index] = value

    def getIparm(self, index):

        """ Get value of iparm[index] """

        return self.iparm[index]

    def getDparm(self, index):

        """ Get value of dparm[index] """

        return self.dparm[index]





    def pastix(self,

               ncol,

               np.ndarray[np_pastix_int_t,ndim=1] colptr,

               np.ndarray[np_pastix_int_t,ndim=1] rows,

               np.ndarray[np.float32_t,  ndim=1] values,

               np.ndarray[np_pastix_int_t,ndim=1] perm,

               np.ndarray[np_pastix_int_t,ndim=1] invp,

               np.ndarray[np.float32_t,  ndim=1] rhs,

               nrhs=1):

        '''

        Generic interface to [sdcz]_pastix.

        Arithmetic depends on the type of values (which should be same ad rhs).

        '''

        self.n = <pastix_int_t>ncol

        self.colptr = <pastix_int_t*>colptr.data

        self.rows = <pastix_int_t*>rows.data

        self.perm = <pastix_int_t*>perm.data

        self.invp = <pastix_int_t*>invp.data

        self.rhs_nbr = nrhs

        print "pouet", self.dparm[self.API.DPARM_EPSILON_REFINEMENT]



        s_pastix(&self.pdata, self.mpi_comm,

                 self.n, self.colptr, self.rows, <float*>values.data,

                 self.perm, self.invp,

                 <float*>rhs.data, self.rhs_nbr,

                 self.iparm, self.dparm)







    def dpastix(self,

                ncol,

                np.ndarray[np_pastix_int_t,ndim=1] colptr,

                np.ndarray[np_pastix_int_t,ndim=1] rows,

                np.ndarray[np.float32_t,  ndim=1] values,

                np.ndarray[np_pastix_int_t,ndim=1] l2g,

                np.ndarray[np_pastix_int_t,ndim=1] perm,

                np.ndarray[np_pastix_int_t,ndim=1] invp,

                np.ndarray[np.float32_t,  ndim=1] rhs,

                nrhs=1):

        '''

        Generic interface to [sdcz]_dpastix (Distributed interface to PaStiX).

        Arithmetic depends on the type of values (which should be same ad rhs).

        '''

        self.iparm[self.API.IPARM_GRAPHDIST] = self.API.API_YES

        self.n = <pastix_int_t>ncol

        self.colptr = <pastix_int_t*>colptr.data

        self.rows = <pastix_int_t*>rows.data

        self.perm = <pastix_int_t*>perm.data

        self.invp = <pastix_int_t*>invp.data

        self.rhs_nbr = nrhs

        s_dpastix(&self.pdata, self.mpi_comm,

                  self.n, self.colptr, self.rows, <float*>values.data,

                  <pastix_int_t*>l2g.data,

                  self.perm, self.invp,

                  <float*>rhs.data, self.rhs_nbr,

                  self.iparm, self.dparm)







    def checkMatrix(self,

                    n,

                    np.ndarray[np_pastix_int_t,ndim=1] colptr,

                    np.ndarray[np_pastix_int_t,ndim=1] rows,

                    np.ndarray[np.float32_t,  ndim=1] values,

                    np.ndarray[np_pastix_int_t,ndim=1] loc2glob=None,

                    verbosity=None,

                    sym=None,

                    correction=None,

                    dofnbr=1):

        '''

        Generic interface to [sdcz]_pastix_checkMatrix()

        Arithmetic depends on the type of values.

        '''



        cdef pastix_int_t *   colptr_c = <pastix_int_t*>colptr.data

        cdef pastix_int_t *   rows_c   = <pastix_int_t*>rows.data

        cdef pastix_int_t *   l2g_c    = NULL

        cdef pastix_int_t **  l2g_ptr  = NULL

        cdef float * values_c = NULL



        cdef pastix_int_t ret = self.API.NO_ERR



        if sym == None:

            sym = self.API.API_SYM_NO

        if verbosity==None:

            verbosity = self.API.API_VERBOSE_NO

        if correction==None:

            correction = self.API.API_NO

        if not loc2glob == None:

            l2g_c = <pastix_int_t*>loc2glob.data

            l2g_ptr = &l2g_c

        values_c = <float *>values.data



        ret = s_pastix_checkMatrix(self.mpi_comm,

                                   <pastix_int_t>verbosity,

                                   <pastix_int_t>sym,

                                   <pastix_int_t>correction,

                                   <pastix_int_t>n,

                                   <pastix_int_t**>&colptr_c,

                                   <pastix_int_t**>&rows_c,

                                   <float**>&values_c,

                                   <pastix_int_t**>l2g_ptr,

                                   <pastix_int_t>dofnbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def bindThreads(self,np.ndarray[np_pastix_int_t,ndim=1] bindTab):

        cdef pastix_int_t nthreads_c

        (nthreads_py) = bindTab.shape

        nthreads_c = <pastix_int_t>nthreads_py

        s_pastix_bindThreads(self.pdata,

                             nthreads_c,

                             <pastix_int_t*>bindTab.data)



    def getLocalNodeLst(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist

        nodenbr = s_pastix_getLocalNodeNbr(&self.pdata)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)

        ret = s_pastix_getLocalNodeLst(&self.pdata,

                                       <pastix_int_t*>nodelist.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def getLocalUnknowLst(self):

        cdef pastix_int_t ret

        cdef pastix_int_t unknownnbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] unknownlist

        unknownnbr = s_pastix_getLocalUnknownNbr(&self.pdata)

        unknownlist = np.zeros(unknownnbr, dtype=self.integerType)

        ret = s_pastix_getLocalUnknownLst(&self.pdata,

                                          <pastix_int_t*>unknownlist.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



        return unknownlist





    def setSchurUnknowLst(self, np.ndarray[np_pastix_int_t,ndim=1] listSchur):

        cdef pastix_int_t ret

        cdef pastix_int_t n_c

        (n_py) = listSchur.shape

        n_c = <pastix_int_t>n_py



        ret = s_pastix_setSchurUnknownList(self.pdata,

                                           <pastix_int_t> n_c,

                                           <pastix_int_t *>listSchur.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getSchurLocalNodeList(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist



        ret = s_pastix_getSchurLocalNodeNbr(self.pdata, &nodenbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)



        ret = s_pastix_getSchurLocalNodeList(self.pdata,

                                             <pastix_int_t*>nodelist.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def getSchurLocalUnknownList(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist



        ret = s_pastix_getSchurLocalUnkownNbr(self.pdata, &nodenbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)



        ret = s_pastix_getSchurLocalUnknownList(self.pdata,

                                                <pastix_int_t*>nodelist.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def setSchurArray(self,

                      np.ndarray[np.float32_t, ndim=1] schur):

        cdef pastix_int_t ret



        ret = s_pastix_setSchurArray(self.pdata,

                                     <float *>schur.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getSchur(self,

                 np.ndarray[np.float32_t, ndim=1] schur):

        cdef pastix_int_t ret



        ret = s_pastix_getSchur(self.pdata,

                                <float *>schur.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getMemoryUsage(self):

        cdef long mem

        cdef long maxmem

        mem = s_pastix_getMemoryUsage()

        maxmem = s_pastix_getMaxMemoryUsage()

        return (mem, maxmem)



    def unscale(self, sym=None):

        if sym==None:

            sym = self.API.API_SYM_NO

        s_pastix_unscale(self.pdata,

                         <pastix_int_t> sym)



    property integerType:

            def __get__(self):

                cdef pastix_int_t integer

                integer = 1

                if sizeof(integer) == 4:

                    return np.int32

                elif sizeof(integer) == 8:

                    return np.int64

                else:

                    raise PastixException




cdef class PASTIX:
    API = ENUM_API

    cdef:

        cdef pastix_data_t  * pdata

        cdef MPI_Comm         mpi_comm

        cdef pastix_int_t   * iparm

        cdef double         * dparm

        cdef pastix_int_t     n

        cdef pastix_int_t   * colptr

        cdef pastix_int_t   * rows

        cdef pastix_float_t * values

        cdef pastix_int_t   * perm

        cdef pastix_int_t   * invp

        cdef pastix_float_t * rhs

        cdef pastix_int_t     rhs_nbr



    def __cinit__(self, MPI.Comm comm not None ):

        cdef MPI_Comm c_comm = comm.ob_mpi

        self.pdata = NULL

        self.mpi_comm = c_comm

        self.iparm = <pastix_int_t*>malloc(self.API.IPARM_SIZE*sizeof(pastix_int_t))

        self.dparm = <double*>malloc(self.API.DPARM_SIZE*sizeof(double))

        self.iparm[self.API.IPARM_MODIFY_PARAMETER] = self.API.API_NO;

        pastix(&self.pdata, self.mpi_comm,

                 0, NULL, NULL, NULL,

                 NULL, NULL, NULL, 1,

                 <pastix_int_t*>self.iparm, self.dparm);

        print self.dparm[self.API.DPARM_EPSILON_REFINEMENT]



    def __dealloc__(self):

        self.iparm[self.API.IPARM_START_TASK] = self.API.API_TASK_CLEAN;

        self.iparm[self.API.IPARM_END_TASK] = self.API.API_TASK_CLEAN;

        if (self.pdata != <pastix_data_t*>NULL):

            pastix(&self.pdata, self.mpi_comm,

                     0, NULL, NULL, NULL,

                     NULL, NULL, NULL, 1,

                     <pastix_int_t*>self.iparm, self.dparm);

        free(self.iparm)

        free(self.dparm)



    def setIparm(self, index, value):

        """ Set value of iparm[index] """

        self.iparm[index] = value

    def setDparm(self, index, value):

        """ Set value of dparm[index] """

        self.dparm[index] = value

    def getIparm(self, index):

        """ Get value of iparm[index] """

        return self.iparm[index]

    def getDparm(self, index):

        """ Get value of dparm[index] """

        return self.dparm[index]





    def pastix(self,

               ncol,

               np.ndarray[np_pastix_int_t,ndim=1] colptr,

               np.ndarray[np_pastix_int_t,ndim=1] rows,

               np.ndarray[np_pastix_float_t,  ndim=1] values,

               np.ndarray[np_pastix_int_t,ndim=1] perm,

               np.ndarray[np_pastix_int_t,ndim=1] invp,

               np.ndarray[np_pastix_float_t,  ndim=1] rhs,

               nrhs=1):

        '''

        Generic interface to [sdcz]_pastix.

        Arithmetic depends on the type of values (which should be same ad rhs).

        '''

        self.n = <pastix_int_t>ncol

        self.colptr = <pastix_int_t*>colptr.data

        self.rows = <pastix_int_t*>rows.data

        self.perm = <pastix_int_t*>perm.data

        self.invp = <pastix_int_t*>invp.data

        self.rhs_nbr = nrhs

        print "pouet", self.dparm[self.API.DPARM_EPSILON_REFINEMENT]



        pastix(&self.pdata, self.mpi_comm,

                 self.n, self.colptr, self.rows, <pastix_float_t*>values.data,

                 self.perm, self.invp,

                 <pastix_float_t*>rhs.data, self.rhs_nbr,

                 self.iparm, self.dparm)







    def dpastix(self,

                ncol,

                np.ndarray[np_pastix_int_t,ndim=1] colptr,

                np.ndarray[np_pastix_int_t,ndim=1] rows,

                np.ndarray[np_pastix_float_t,  ndim=1] values,

                np.ndarray[np_pastix_int_t,ndim=1] l2g,

                np.ndarray[np_pastix_int_t,ndim=1] perm,

                np.ndarray[np_pastix_int_t,ndim=1] invp,

                np.ndarray[np_pastix_float_t,  ndim=1] rhs,

                nrhs=1):

        '''

        Generic interface to [sdcz]_dpastix (Distributed interface to PaStiX).

        Arithmetic depends on the type of values (which should be same ad rhs).

        '''

        self.iparm[self.API.IPARM_GRAPHDIST] = self.API.API_YES

        self.n = <pastix_int_t>ncol

        self.colptr = <pastix_int_t*>colptr.data

        self.rows = <pastix_int_t*>rows.data

        self.perm = <pastix_int_t*>perm.data

        self.invp = <pastix_int_t*>invp.data

        self.rhs_nbr = nrhs

        dpastix(&self.pdata, self.mpi_comm,

                  self.n, self.colptr, self.rows, <pastix_float_t*>values.data,

                  <pastix_int_t*>l2g.data,

                  self.perm, self.invp,

                  <pastix_float_t*>rhs.data, self.rhs_nbr,

                  self.iparm, self.dparm)







    def checkMatrix(self,

                    n,

                    np.ndarray[np_pastix_int_t,ndim=1] colptr,

                    np.ndarray[np_pastix_int_t,ndim=1] rows,

                    np.ndarray[np_pastix_float_t,  ndim=1] values,

                    np.ndarray[np_pastix_int_t,ndim=1] loc2glob=None,

                    verbosity=None,

                    sym=None,

                    correction=None,

                    dofnbr=1):

        '''

        Generic interface to [sdcz]_pastix_checkMatrix()

        Arithmetic depends on the type of values.

        '''



        cdef pastix_int_t *   colptr_c = <pastix_int_t*>colptr.data

        cdef pastix_int_t *   rows_c   = <pastix_int_t*>rows.data

        cdef pastix_int_t *   l2g_c    = NULL

        cdef pastix_int_t **  l2g_ptr  = NULL

        cdef pastix_float_t * values_c = NULL



        cdef pastix_int_t ret = self.API.NO_ERR



        if sym == None:

            sym = self.API.API_SYM_NO

        if verbosity==None:

            verbosity = self.API.API_VERBOSE_NO

        if correction==None:

            correction = self.API.API_NO

        if not loc2glob == None:

            l2g_c = <pastix_int_t*>loc2glob.data

            l2g_ptr = &l2g_c

        values_c = <pastix_float_t *>values.data



        ret = pastix_checkMatrix(self.mpi_comm,

                                   <pastix_int_t>verbosity,

                                   <pastix_int_t>sym,

                                   <pastix_int_t>correction,

                                   <pastix_int_t>n,

                                   <pastix_int_t**>&colptr_c,

                                   <pastix_int_t**>&rows_c,

                                   <pastix_float_t**>&values_c,

                                   <pastix_int_t**>l2g_ptr,

                                   <pastix_int_t>dofnbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def bindThreads(self,np.ndarray[np_pastix_int_t,ndim=1] bindTab):

        cdef pastix_int_t nthreads_c

        (nthreads_py) = bindTab.shape

        nthreads_c = <pastix_int_t>nthreads_py

        pastix_bindThreads(self.pdata,

                             nthreads_c,

                             <pastix_int_t*>bindTab.data)



    def getLocalNodeLst(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist

        nodenbr = pastix_getLocalNodeNbr(&self.pdata)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)

        ret = pastix_getLocalNodeLst(&self.pdata,

                                       <pastix_int_t*>nodelist.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def getLocalUnknowLst(self):

        cdef pastix_int_t ret

        cdef pastix_int_t unknownnbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] unknownlist

        unknownnbr = pastix_getLocalUnknownNbr(&self.pdata)

        unknownlist = np.zeros(unknownnbr, dtype=self.integerType)

        ret = pastix_getLocalUnknownLst(&self.pdata,

                                          <pastix_int_t*>unknownlist.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



        return unknownlist





    def setSchurUnknowLst(self, np.ndarray[np_pastix_int_t,ndim=1] listSchur):

        cdef pastix_int_t ret

        cdef pastix_int_t n_c

        (n_py) = listSchur.shape

        n_c = <pastix_int_t>n_py



        ret = pastix_setSchurUnknownList(self.pdata,

                                           <pastix_int_t> n_c,

                                           <pastix_int_t *>listSchur.data)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getSchurLocalNodeList(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist



        ret = pastix_getSchurLocalNodeNbr(self.pdata, &nodenbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)



        ret = pastix_getSchurLocalNodeList(self.pdata,

                                             <pastix_int_t*>nodelist.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def getSchurLocalUnknownList(self):

        cdef pastix_int_t ret

        cdef pastix_int_t nodenbr

        cdef np.ndarray[np_pastix_int_t,ndim=1] nodelist



        ret = pastix_getSchurLocalUnkownNbr(self.pdata, &nodenbr)

        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        nodelist = np.zeros(nodenbr, dtype=self.integerType)



        ret = pastix_getSchurLocalUnknownList(self.pdata,

                                                <pastix_int_t*>nodelist.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)

        return nodelist



    def setSchurArray(self,

                      np.ndarray[np_pastix_float_t, ndim=1] schur):

        cdef pastix_int_t ret



        ret = pastix_setSchurArray(self.pdata,

                                     <pastix_float_t *>schur.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getSchur(self,

                 np.ndarray[np_pastix_float_t, ndim=1] schur):

        cdef pastix_int_t ret



        ret = pastix_getSchur(self.pdata,

                                <pastix_float_t *>schur.data)



        if not ret == self.API.NO_ERR:

            raise PastixException(u"Pastix returned %d" %ret)



    def getMemoryUsage(self):

        cdef long mem

        cdef long maxmem

        mem = pastix_getMemoryUsage()

        maxmem = pastix_getMaxMemoryUsage()

        return (mem, maxmem)



    def unscale(self, sym=None):

        if sym==None:

            sym = self.API.API_SYM_NO

        pastix_unscale(self.pdata,

                         <pastix_int_t> sym)



    property integerType:

            def __get__(self):

                cdef pastix_int_t integer

                integer = 1

                if sizeof(integer) == 4:

                    return np.int32

                elif sizeof(integer) == 8:

                    return np.int64

                else:

                    raise PastixException




