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

