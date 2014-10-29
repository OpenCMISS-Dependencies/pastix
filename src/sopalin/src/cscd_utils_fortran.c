/* Copyright 2008 BORDEAUX I UNIVERSITY & INRIA 
**
** This file is part of the PaStiX parallel sparse matrix package.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/*
  File: pastix_fortran.c

  Interface to the PaStiX API functions.

 */
#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif
#include "common_pastix.h"
#include "cscd_utils.h"


#if (defined X_ARCHpower_ibm_aix)
#define FORTRAN_CALL(nom) PASTIX_EXTERN_F(nom)
#else
#define FORTRAN_CALL(nom) PASTIX_EXTERN_F(nom ## _)
#endif
/*
  Struct: csc_data_

  Contains the new CSCD

*/
struct csc_data_ {
  PASTIX_INT     n;
  PASTIX_INT   * colptr;
  PASTIX_INT   * rows;
  PASTIX_FLOAT * values;
  PASTIX_FLOAT * rhs;
  PASTIX_INT     nrhs;
  PASTIX_INT   * perm;
  PASTIX_INT   * l2g;
  PASTIX_INT     dof;
};

/*
  Typedef: csc_data_t

  Type coresponding to the struct <csc_data_>
*/
typedef struct csc_data_ csc_data_t;

void FORTRAN_CALL(csc_dispatch_fortran)(csc_data_t ** csc_data,
                                        PASTIX_INT         *gN,
                                        PASTIX_INT         *gcolptr,
                                        PASTIX_INT         *grow,
                                        PASTIX_FLOAT       *gavals,
                                        PASTIX_FLOAT       *grhs,
                                        PASTIX_INT         *gperm,
                                        PASTIX_INT         *ginvp,
                                        int         *dispatch,
                                        PASTIX_INT         *newn,
                                        PASTIX_INT         *newnnz,
                                        MPI_Fint    *fortran_comm)
{
  MPI_Comm        pastix_comm;
  (void)fortran_comm;

  pastix_comm = MPI_Comm_f2c(*fortran_comm);
  MALLOC_INTERN(*csc_data, 1, struct csc_data_);
  (*csc_data)->n      = 0;
  (*csc_data)->colptr = NULL;
  (*csc_data)->rows   = NULL;
  (*csc_data)->values = NULL;
  (*csc_data)->rhs    = NULL;
  (*csc_data)->perm   = NULL;
  (*csc_data)->l2g    = NULL;

  csc_dispatch(*gN, gcolptr, grow, gavals, grhs, gperm, ginvp,
               &((*csc_data)->n), &((*csc_data)->colptr), &((*csc_data)->rows), &((*csc_data)->values),
               &((*csc_data)->rhs), &((*csc_data)->perm),
               &((*csc_data)->l2g), *dispatch, pastix_comm);

  *newn   = (*csc_data)->n;
  *newnnz = (*csc_data)->colptr[(*csc_data)->n]-1;

}

void FORTRAN_CALL(csc_dispatch_fortran_end)(csc_data_t ** csc_data,
                                            PASTIX_INT         *lcolptr,
                                            PASTIX_INT         *lrow,
                                            PASTIX_FLOAT       *lavals,
                                            PASTIX_FLOAT       *lrhs,
                                            PASTIX_FLOAT       *lperm,
                                            PASTIX_INT         *l2g)
{
  PASTIX_INT nnz = 0;
  if ((*csc_data)->colptr != NULL) {
    nnz = (*csc_data)->colptr[(*csc_data)->n]-1;
    memcpy(lcolptr, (*csc_data)->colptr, (1+(*csc_data)->n)*sizeof(PASTIX_INT));
    free((*csc_data)->colptr);
  }
  if ((*csc_data)->rows != NULL) {
    memcpy(lrow,    (*csc_data)->rows,   nnz*sizeof(PASTIX_INT));
    free((*csc_data)->rows);
  }
  if ((*csc_data)->values != NULL) {
    memcpy(lavals,  (*csc_data)->values,   nnz*sizeof(PASTIX_FLOAT));
    free((*csc_data)->values);
  }
  if ((*csc_data)->rhs != NULL) {
    memcpy(lrhs, (*csc_data)->rhs, (*csc_data)->n*sizeof(PASTIX_FLOAT));
    free((*csc_data)->rhs);
  }
  if ((*csc_data)->perm != NULL) {
    memcpy(lperm, (*csc_data)->perm, (*csc_data)->n*sizeof(PASTIX_INT));
    free((*csc_data)->perm);
  }
  if ((*csc_data)->l2g != NULL) {
    memcpy(l2g, (*csc_data)->l2g, (*csc_data)->n*sizeof(PASTIX_INT));
    free((*csc_data)->l2g);
  }
}

void FORTRAN_CALL(cscd_redispatch_fortran)(csc_data_t ** csc_data,
                                           PASTIX_INT         * n,
                                           PASTIX_INT         * ia,
                                           PASTIX_INT         * ja,
                                           PASTIX_FLOAT       * a,
                                           PASTIX_FLOAT       * rhs,
                                           PASTIX_INT         * nrhs,
                                           PASTIX_INT         * l2g,
                                           PASTIX_INT         * newn,
                                           PASTIX_INT         * newl2g,
                                           PASTIX_INT         * newnnz,
                                           MPI_Fint    * fortran_comm,
                                           PASTIX_INT         * dof,
                                           int         * ierr)
{

  MPI_Comm        pastix_comm;
  (void)newnnz; (void)fortran_comm;

  pastix_comm = MPI_Comm_f2c(*fortran_comm);
  MALLOC_INTERN(*csc_data, 1, struct csc_data_);
  (*csc_data)->n      = 0;
  (*csc_data)->colptr = NULL;
  (*csc_data)->rows   = NULL;
  (*csc_data)->values = NULL;
  (*csc_data)->rhs    = NULL;
  (*csc_data)->nrhs   = *nrhs;
  (*csc_data)->perm   = NULL;
  (*csc_data)->l2g    = NULL;
  (*csc_data)->dof    = *dof;

  *ierr = cscd_redispatch(*n,    ia, ja, a, rhs,  *nrhs, l2g,
                          *newn, &((*csc_data)->colptr), &((*csc_data)->rows), &((*csc_data)->values),  &((*csc_data)->rhs), newl2g,
                          pastix_comm, (*csc_data)->dof);
}

void FORTRAN_CALL(cscd_redispatch_fortran_end)(csc_data_t ** csc_data,
                                               PASTIX_INT         *dcolptr,
                                               PASTIX_INT         *drow,
                                               PASTIX_FLOAT       *davals,
                                               PASTIX_FLOAT       *drhs)
{
  PASTIX_INT nnz = 0;
  if ((*csc_data)->colptr != NULL) {
    nnz = (*csc_data)->colptr[(*csc_data)->n]-1;
    memcpy(dcolptr, (*csc_data)->colptr, (1+(*csc_data)->n)*sizeof(PASTIX_INT));
    free((*csc_data)->colptr);
  }
  if ((*csc_data)->rows != NULL) {
    memcpy(drow,    (*csc_data)->rows,   nnz*sizeof(PASTIX_INT));
    free((*csc_data)->rows);
  }
  if ((*csc_data)->values != NULL) {
    memcpy(davals,  (*csc_data)->values,
           (*csc_data)->dof*(*csc_data)->dof*nnz*sizeof(PASTIX_FLOAT));
    free((*csc_data)->values);
  }
  if ((*csc_data)->rhs != NULL) {
    memcpy(drhs, (*csc_data)->rhs, (*csc_data)->nrhs*(*csc_data)->n*sizeof(PASTIX_FLOAT));
    free((*csc_data)->rhs);
  }
  memFree_null(*csc_data);
}
