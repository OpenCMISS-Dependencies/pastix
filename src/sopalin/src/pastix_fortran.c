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
 * File: pastix_fortran.c
 *
 * Interface to the PaStiX API functions.
 *
 */
#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif
#include "common_pastix.h"
#include "pastix.h"
#include "cscd_utils.h"

#ifdef NO_TYPE_PREFIX
#define FORTRAN_NAME(nu,nl,pl,pc)   \
  void nl pl;                       \
  void nu pl                        \
  { nl pc; }                        \
  void nl ## _ pl                   \
  { nl pc; }                        \
  void nl ## __ pl                  \
  { nl pc; }
#else
#ifdef TYPE_COMPLEX
#ifdef PREC_DOUBLE
#define FORTRAN_NAME(nu,nl,pl,pc)         \
  void z_ ## nl pl;                       \
  void Z_ ## nu pl                        \
  { z_ ## nl pc; }                        \
  void z_ ## nl ## _ pl                   \
  { z_ ## nl pc; }                        \
  void z_ ## nl ## __ pl                  \
  { z_ ## nl pc; }
#else
#define FORTRAN_NAME(nu,nl,pl,pc)         \
  void c_ ## nl pl;                       \
  void C_ ## nu pl                        \
  { c_ ## nl pc; }                        \
  void c_ ## nl ## _ pl                   \
  { c_ ## nl pc; }                        \
  void c_ ## nl ## __ pl                  \
  { c_ ## nl pc; }
#endif
#else
#ifdef PREC_DOUBLE
#define FORTRAN_NAME(nu,nl,pl,pc)         \
  void d_ ## nl pl;                       \
  void D_ ## nu pl                        \
  { d_ ## nl pc; }                        \
  void d_ ## nl ## _ pl                   \
  { d_ ## nl pc; }                        \
  void d_ ## nl ## __ pl                  \
  { d_ ## nl pc; }
#else
#define FORTRAN_NAME(nu,nl,pl,pc)         \
  void s_ ## nl pl;                       \
  void S_ ## nu pl                        \
  { s_ ## nl pc; }                        \
  void s_ ## nl ## _ pl                   \
  { s_ ## nl pc; }                        \
  void s_ ## nl ## __ pl                  \
  { s_ ## nl pc; }
#endif
#endif
#endif

#define pastix_fortran                 PASTIX_EXTERN_F(pastix_fortran)
#define dpastix_fortran                PASTIX_EXTERN_F(dpastix_fortran)
#define pastix_fortran_setbindtab               \
  PASTIX_EXTERN_F(pastix_fortran_setbindtab)
#define pastix_fortran_checkmatrix              \
  PASTIX_EXTERN_F(pastix_fortran_checkmatrix)
#define pastix_fortran_checkmatrix_end            \
  PASTIX_EXTERN_F(pastix_fortran_checkmatrix_end)
#define pastix_fortran_getlocalnodenbr            \
  PASTIX_EXTERN_F(pastix_fortran_getlocalnodenbr)
#define pastix_fortran_getlocalnodelst            \
  PASTIX_EXTERN_F(pastix_fortran_getlocalnodelst)
#define cscd_addlocal_fortran          PASTIX_EXTERN_F(cscd_addlocal_fortran)
#define pastix_fortran_setschurunknownlist            \
  PASTIX_EXTERN_F(pastix_fortran_setschurunknownlist)
#define pastix_fortran_getschur        PASTIX_EXTERN_F(pastix_fortran_getschur)
#define pastix_fortran_getschurlocalnodenbr             \
  PASTIX_EXTERN_F(pastix_fortran_getschurlocalnodenbr)
#define pastix_fortran_getschurlocalunknownnbr            \
  PASTIX_EXTERN_F(pastix_fortran_getschurlocalunknownnbr)
#define pastix_fortran_getschurlocalnodelist            \
  PASTIX_EXTERN_F(pastix_fortran_getschurlocalnodelist)
#define pastix_fortran_getschurlocalunknownlist             \
  PASTIX_EXTERN_F(pastix_fortran_getschurlocalunknownlist)
#define pastix_fortran_setschurarray            \
  PASTIX_EXTERN_F(pastix_fortran_setschurarray)

/*
 * Struct: pastix_check_data_
 *
 * Contains the rows and values arrays.
 *
 */
struct pastix_check_data_ {
  PASTIX_INT     n;
  PASTIX_INT     nnz;
  PASTIX_INT   * colptr;
  PASTIX_INT   * rows;
  PASTIX_FLOAT * values;
};
/*
 * Typedef: pastix_check_data_t
 *
 * Type coresponding to the struct <pastix_check_data_>
 */
typedef struct pastix_check_data_ pastix_check_data_t;


/*
 * Function: pastix_fortran
 *
 * Fortran interface to the <pastix> function.
 *
 * Parameters:
 *   pastix_data  - Data used for a step by step execution.
 *   fortran_comm - MPI communicator which compute the resolution.
 *   n            - Size of the system.
 *   colptr       - Tabular containing the start of each column
 *                  in row and avals tabulars.
 *   row          - Tabular containing the row number for each element
 *                  sorted by column.
 *   avals        - Tabular containing the values of each element
 *                  sorted by column.
 *   perm         - Permutation tabular for the renumerotation of the unknowns.
 *   invp         - Reverse permutation tabular for the renumerotation
 *                  of the unknowns.
 *   b            - Right hand side vector(s).
 *   rhs          - Number of right hand side vector(s).
 *   iparm        - Integer parameters given to pastix.
 *   dparm        - Double parameters given to pâstix.
 */
#define pastix_fortran PASTIX_EXTERN_F(pastix_fortran)
void pastix_fortran ( void     *pastix_data,
                      MPI_Fint *fortran_comm,
                      PASTIX_INT      *n,
                      PASTIX_INT      *colptr,
                      PASTIX_INT      *row,
                      PASTIX_FLOAT    *avals,
                      PASTIX_INT      *perm,
                      PASTIX_INT      *invp,
                      PASTIX_FLOAT    *b,
                      PASTIX_INT      *rhs,
                      PASTIX_INT      *iparm,
                      double   *dparm )
{
  pastix_data_t **pastix_data2;
  MPI_Comm        pastix_comm;
  (void)fortran_comm;

  pastix_comm  = MPI_Comm_f2c(*fortran_comm);
  pastix_data2 = pastix_data;

#ifdef DEBUG_PASTIX_FORTRAN
  {
    int i;
    fprintf(stdout,"%ld %ld %ld %lf %ld %ld %lf %ld %ld %lf\n",
            (long)*n, (long)colptr[0],row[0],avals[0],perm[0],
            invp[0],b[0],*rhs,iparm[0],dparm[0]);
    for (i=0; i<10;i++)
      fprintf(stdout,"%lf\n",avals[i]);
  }
#endif

  pastix(pastix_data2, pastix_comm, *n, colptr, row, avals,
   perm, invp, b, *rhs, iparm, dparm);

}

FORTRAN_NAME(PASTIX_FORTRAN,
             pastix_fortran,
             ( void     *pastix_data,
               MPI_Fint *fortran_comm,
               PASTIX_INT      *n,
               PASTIX_INT      *colptr,
               PASTIX_INT      *row,
               PASTIX_FLOAT    *avals,
               PASTIX_INT      *perm,
               PASTIX_INT      *invp,
               PASTIX_FLOAT    *b,
               PASTIX_INT      *rhs,
               PASTIX_INT      *iparm,
               double   *dparm ),
             ( pastix_data,
               fortran_comm,
               n,
               colptr,
               row,
               avals,
               perm,
               invp,
               b,
               rhs,
               iparm,
               dparm ))
/*
 * Function: dpastix_fortran
 *
 * Fortran interface to the <dpastix> function.
 *
 * Parameters:
 *   pastix_data  - Data used for a step by step execution.
 *   fortran_comm - MPI communicator which compute the resolution.
 *   n            - Size of the system.
 *   colptr       - Tabular containing the start of each column
 *                  in row and avals tabulars.
 *   row          - Tabular containing the row number for each element
 *                  sorted by column.
 *   avals        - Tabular containing the values of each element
 *                  sorted by column.
 *   loc2glob     - Global column number of the local columns.
 *   perm         - Permutation tabular for the renumerotation of the unknowns.
 *   invp         - Reverse permutation tabular for the renumerotation
 *                  of the unknowns.
 *   b            - Right hand side vector(s).
 *   rhs          - Number of right hand side vector(s).
 *   iparm        - Integer parameters given to pastix.
 *   dparm        - Double parameters given to pâstix.
 */
#define dpastix_fortran PASTIX_EXTERN_F(dpastix_fortran)
void dpastix_fortran ( void     *pastix_data,
                       MPI_Fint *fortran_comm,
                       PASTIX_INT      *n,
                       PASTIX_INT      *colptr,
                       PASTIX_INT      *row,
                       PASTIX_FLOAT    *avals,
                       PASTIX_INT      *loc2glob,
                       PASTIX_INT      *perm,
                       PASTIX_INT      *invp,
                       PASTIX_FLOAT    *b,
                       PASTIX_INT      *rhs,
                       PASTIX_INT      *iparm,
                       double   *dparm )
{
  pastix_data_t **pastix_data2;
  MPI_Comm        pastix_comm;
  (void)fortran_comm;

  pastix_comm  = MPI_Comm_f2c(*fortran_comm);
  pastix_data2 = pastix_data;

#ifdef DEBUG_PASTIX_FORTRAN
  {
    int i;
    fprintf(stdout,"%ld %ld %ld %lf %ld %ld %lf %ld %ld %lf\n",
            (long)*n, (long)colptr[0],row[0],avals[0],perm[0],
            invp[0],b[0],*rhs,iparm[0],dparm[0]);
    for (i=0; i<10;i++)
      fprintf(stdout,"%lf\n",avals[i]);
  }
#endif

  dpastix(pastix_data2, pastix_comm, *n, colptr, row, avals,
    loc2glob, perm, invp, b, *rhs, iparm, dparm);
}

FORTRAN_NAME(DPASTIX_FORTRAN,
             dpastix_fortran,
             ( void     *pastix_data,
               MPI_Fint *fortran_comm,
               PASTIX_INT      *n,
               PASTIX_INT      *colptr,
               PASTIX_INT      *row,
               PASTIX_FLOAT    *avals,
               PASTIX_INT      *loc2glob,
               PASTIX_INT      *perm,
               PASTIX_INT      *invp,
               PASTIX_FLOAT    *b,
               PASTIX_INT      *rhs,
               PASTIX_INT      *iparm,
               double   *dparm ),
             ( pastix_data,
               fortran_comm,
               n,
               colptr,
               row,
               avals,
               loc2glob,
               perm,
               invp,
               b,
               rhs,
               iparm,
               dparm ))
/*
 * Function: pastix_fortran_bindthreads
 *
 * Set rules to bind thread on processors.
 *
 * Parameters:
 *   pastix_data - Structure to keep data between PaStiX calls
 *   thrdnbr     - Number of threads to use in pastix by MPI node.
 *   bindtab     - Tabular indicating which processor each thread
 *                 has to be binded to.
 */
#define pastix_fortran_bindthreads PASTIX_EXTERN_F(pastix_fortran_bindthreads)
void pastix_fortran_bindthreads ( void *pastix_data,
                                  PASTIX_INT  *thrdnbr,
                                  PASTIX_INT  *bindtab )
{
  pastix_bindThreads (*((pastix_data_t**)pastix_data), *thrdnbr, bindtab);
}

FORTRAN_NAME(PASTIX_FORTRAN_BINDTHREADS,
             pastix_fortran_bindthreads,
             ( void *pastix_data,
               PASTIX_INT  *thrdnbr,
               PASTIX_INT  *bindtab ),
             ( pastix_data,
               thrdnbr,
               bindtab ))
/*
 * Function: pastix_fortran_checkmatrix
 *
 * Check if the matrix is correct and can be given to PaStiX solver.
 *
 * Parameters:
 *   data_check   - Integer which will contain the <pastix_check_data_t>
 *                  structure adress (output)
 *   fortran_comm - PaStiX MPI communicator.
 *   verb         - Verbosity level.
 *   flagsym      - Flag indicating if the matrix is symetric
 *                  (API_SYM_YES or API_SYM_NO)
 *   flagcor      - Indicate if we permit the function to reallocate the matrix.
 *   n            - Size of the matrix.
 *   colptr       - Tabular containing the start of each column
 *                  in row and avals tabulars.
 *   row          - Tabular containing the row number for each element
 *                  sorted by column.
 *   avals        - Tabular containing the values of each element
 *                  sorted by column.
 *   loc2glob     - Global column number of the local columns,
 *                  -1 if not distributed.
 *   dof          - Number of degrees of freedom.
 */
#define pastix_fortran_checkmatrix PASTIX_EXTERN_F(pastix_fortran_checkmatrix)
void pastix_fortran_checkmatrix ( pastix_check_data_t **data_check,
                                  MPI_Fint             *fortran_comm,
                                  PASTIX_INT                  *verb,
                                  PASTIX_INT                  *flagsym,
                                  PASTIX_INT                  *flagcor,
                                  PASTIX_INT                  *n,
                                  PASTIX_INT                  *colptr,
                                  PASTIX_INT                  *row,
                                  PASTIX_FLOAT                *avals,
                                  PASTIX_INT                  *loc2glob,
                                  PASTIX_INT                  *dof)
{
  MPI_Comm        pastix_comm;
  PASTIX_INT            *tmprows;
  PASTIX_FLOAT          *tmpvalues;
  PASTIX_INT             nnz;
  (void)fortran_comm;

  pastix_comm = MPI_Comm_f2c(*fortran_comm);
  nnz = colptr[*n]-colptr[0];
  *data_check=NULL;
  if (*flagcor == API_YES)
    {
      MALLOC_INTERN(*data_check, 1, struct pastix_check_data_);
      MALLOC_EXTERN((*data_check)->rows,   nnz,        PASTIX_INT);
      MALLOC_EXTERN((*data_check)->values, nnz*(*dof), PASTIX_FLOAT);
      memcpy((*data_check)->rows,   row,   nnz*sizeof(PASTIX_INT));
      memcpy((*data_check)->values, avals, nnz*(*dof)*sizeof(PASTIX_FLOAT));
      tmprows   = (*data_check)->rows;
      tmpvalues = (*data_check)->values;
      (*data_check)->colptr=colptr;
      (*data_check)->n=*n;
    }
  else
    {
      tmprows   = row;
      tmpvalues = avals;
    }
  pastix_checkMatrix(pastix_comm, *verb, *flagsym, *flagcor,
                     *n, &colptr, &tmprows, &tmpvalues,
                     ((*loc2glob == -1)?NULL:(&loc2glob)), *dof);


  if (*flagcor == API_YES)
    {
      if (colptr[*n]-1 == nnz)
        {
          memcpy(row,   tmprows,   nnz*sizeof(PASTIX_INT));
          memcpy(avals, tmpvalues, nnz*(*dof)*sizeof(PASTIX_FLOAT));
          free(tmprows);
          free(tmpvalues);
          memFree_null(*data_check);
        }
      else
        {
          (*data_check)->nnz = colptr[*n]-1;
          (*data_check)->rows   = tmprows;
          (*data_check)->values = tmpvalues;

          if (*verb > API_VERBOSE_NOT)
            errorPrintW("Number of non zeros has changed,"
                        " don't forget to call pastix_fortran_checkmatrix_end");

        }
    }
}

FORTRAN_NAME( PASTIX_FORTRAN_CHECKMATRIX,
              pastix_fortran_checkmatrix,
              ( pastix_check_data_t **data_check,
                MPI_Fint             *fortran_comm,
                PASTIX_INT                  *verb,
                PASTIX_INT                  *flagsym,
                PASTIX_INT                  *flagcor,
                PASTIX_INT                  *n,
                PASTIX_INT                  *colptr,
                PASTIX_INT                  *row,
                PASTIX_FLOAT                *avals,
                PASTIX_INT                  *loc2glob,
                PASTIX_INT                  *dof ),
              ( data_check,
                fortran_comm,
                verb,
                flagsym,
                flagcor,
                n,
                colptr,
                row,
                avals,
                loc2glob,
                dof) )
/*
 * Function: pastix_fortran_checkmatrix_end
 *
 * Copy the corrected CSC in the user CSC.
 *
 * Parameters:
 *   data_check   - Integer which will contain the <pastix_check_data_t>
 *                  structure adress (output)
 *   verb         - Verbosity level.
 *   row          - Tabular into which to copy the new rows
 *   values       - Tabular into which to copy the new rows
 *   dof          - Number of degrees of freedom.
 */
#define pastix_fortran_checkmatrix_end PASTIX_EXTERN_F(pastix_fortran_checkmatrix_end)
void pastix_fortran_checkmatrix_end ( pastix_check_data_t **data_check,
                                      PASTIX_INT                  *verb,
                                      PASTIX_INT                  *row,
                                      PASTIX_FLOAT                *avals,
                                      PASTIX_INT                  *dof )
{
  (void)verb;
  /* pas de checkmatrix_end si pas de duplicates */
  ASSERT(*data_check!=NULL,MOD_SOPALIN);

  memcpy(row,   (*data_check)->rows,   (*data_check)->nnz*sizeof(PASTIX_INT));
  memcpy(avals, (*data_check)->values, (*data_check)->nnz*(*dof)*sizeof(PASTIX_FLOAT));
  free((*data_check)->rows);
  free((*data_check)->values);
  memFree_null(*data_check);

}
FORTRAN_NAME( PASTIX_FORTRAN_CHECKMATRIX_END,
              pastix_fortran_checkmatrix_end,
              ( pastix_check_data_t **data_check,
                PASTIX_INT                  *verb,
                PASTIX_INT                  *row,
                PASTIX_FLOAT                *avals,
                PASTIX_INT                  *dof ),
              ( data_check,
                verb,
                row,
                avals,
                dof ))
/*
 * Function: pastix_fortran_getlocalnodenbr
 *
 * Give the node number in the new distribution computed by blend.
 * Needs blend to be runned with pastix_data before.
 *
 * Parameters:
 *   pastix_data - Data used for a step by step execution.
 *   nodenbr     - pastix_int_t where to write node number.
 */
void pastix_fortran_getlocalnodenbr ( pastix_data_t ** pastix_data,
                                      PASTIX_INT            * nodenbr )
{
  *nodenbr = pastix_getLocalNodeNbr(pastix_data);
}

FORTRAN_NAME( PASTIX_FORTRAN_GETLOCALNODENBR,
              pastix_fortran_getlocalnodenbr,
              ( pastix_data_t ** pastix_data,
                PASTIX_INT            * nodenbr ),
              ( pastix_data,
                nodenbr ))
/*
 * Function: pastix_fortran_getlocalnodelst
 *
 * Fill in nodelst with the list of local nodes/clumns.
 * Needs nodelst to be allocated with nodenbr*sizeof(pastix_int_t),
 * where nodenbr has been computed by
 * <pastix_fortran_getLocalnodenbr>.
 *
 * Parameters:
 *   pastix_data - Data used for a step by step execution.
 *   nodelst     - An array where to write the list of local nodes/columns.
 */
void pastix_fortran_getlocalnodelst ( pastix_data_t **pastix_data,
                                     PASTIX_INT            *nodelst )
{
  pastix_getLocalNodeLst(pastix_data, nodelst);
}

FORTRAN_NAME( PASTIX_FORTRAN_GETLOCALNODELST,
              pastix_fortran_getlocalnodelst,
              ( pastix_data_t **pastix_data,
                PASTIX_INT            *nodelst ),
              ( pastix_data,
                nodelst ))
/*
 * Function: cscd_addlocal_fortran
 *
 * Interface to <cscd_addlocal>, for fortran calls.
 * Add second cscd to first cscd into third cscd (unallocated)
 * Only adds columns from the second CSCD which belongs to the first one.
 *
 * Parameters:
 *   n           - First cscd size
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   l2g         - local 2 global column numbers for first cscd
 *   addn        - CSCD to add size
 *   addia       - CSCD to add starting index of each column in *addja* and
 *                 *adda*
 *   addja       - Row of each element in second CSCD
 *   adda        - value of each cscd in second CSCD (can be NULL -> add 0)
 *   addl2g      - local 2 global column numbers for second cscd
 *   newn        - new cscd size (same as first)
 *   newia       - CSCD to add starting index of each column in *newja* and
 *                 *newwa*
 *   newja       - Row of each element in third CSCD
 *   newa        - value of each cscd in third CSCD
 *   OP          - Operation to manage common CSCD coefficients.
 *   dof         - Number of degrees of freedom.
 */
void cscd_addlocal_fortran ( PASTIX_INT   n,
                             PASTIX_INT *  ia,    PASTIX_INT *  ja,    PASTIX_FLOAT *  a,
                             PASTIX_INT * l2g,
                             PASTIX_INT   addn,
                             PASTIX_INT *  addia, PASTIX_INT *  addja, PASTIX_FLOAT *  adda,
                             PASTIX_INT * addl2g,
                             PASTIX_INT * newn,
                             PASTIX_INT ** newia, PASTIX_INT ** newja, PASTIX_FLOAT ** newa,
                             PASTIX_INT * OP, PASTIX_INT dof )
{
  cscd_addlocal(n   , ia   , ja   , a   ,  l2g,
                addn, addia, addja, adda,  addl2g,
                newn, newia, newja, newa,  (CSCD_OPERATIONS_t)(*OP), dof);
}
FORTRAN_NAME( CSCD_ADDLOCAL_FORTRAN,
              cscd_addlocal_fortran,
              ( PASTIX_INT   n,
                PASTIX_INT *  ia,    PASTIX_INT *  ja,    PASTIX_FLOAT *  a,
                PASTIX_INT * l2g,
                PASTIX_INT   addn,
                PASTIX_INT *  addia, PASTIX_INT *  addja, PASTIX_FLOAT *  adda,
                PASTIX_INT * addl2g,
                PASTIX_INT * newn,
                PASTIX_INT ** newia, PASTIX_INT ** newja, PASTIX_FLOAT ** newa,
                PASTIX_INT * OP, PASTIX_INT dof ),
              ( n, ia, ja, a, l2g,
                addn, addia, addja, adda, addl2g,
                newn, newia, newja, newa, OP, dof ))
/*
 * Function: pastix_fortran_setschurunknownlist
 *
 * Fortran interface to <pastix_setSchurUnknownList>
 *
 * Parameters:
 *   pastix_data - Data used for a step by step execution.
 *   n           - Number of unknowns.
 *   list        - List of unknowns.
 *   ierr        - error value.
 */
void pastix_fortran_setschurunknownlist (void * pastix_data,
                                         PASTIX_INT *n,
                                         PASTIX_INT *list,
                                         PASTIX_INT *ierr)
{
  pastix_data_t **pastix_data2;
  pastix_data2 = pastix_data;
  *ierr = pastix_setSchurUnknownList(*pastix_data2,
                                     *n,
                                     list);
}

FORTRAN_NAME( PASTIX_FORTRAN_SETSCHURUNKNOWNLIST,
              pastix_fortran_setschurunknownlist,
              (void * pastix_data,
               PASTIX_INT *n,
               PASTIX_INT *list,
               PASTIX_INT *ierr),
              (pastix_data,
               n,
               list,
               ierr))
/*
 * Function: pastix_getschur
 *
 * Fotran interface to <pastix_getSchur>
 *
 * Schur complement is a dense block in a
 * column scheme.
 *
 * Parameters:
 *   pastix_data - Data used for a step by step execution.
 *   schur       - Array to fill-in with Schur complement.
 *   ierr        - error value.
 */
void pastix_fortran_getschur ( void          * pastix_data,
                               PASTIX_FLOAT         * schur,
                               PASTIX_INT           * ierr )
{

  pastix_data_t **pastix_data2;
  pastix_data2 = pastix_data;
  *ierr = pastix_getSchur(*pastix_data2,
                          schur);
}

FORTRAN_NAME( PASTIX_FORTRAN_GETSCHUR,
              pastix_fortran_getschur,
              ( void          * pastix_data,
                PASTIX_FLOAT         * schur,
                PASTIX_INT           * ierr ),
              ( pastix_data,
                schur,
                ierr ))
/*
 * Function: pastix_fortran_getschurlocalnodenbr
 *
 * Fortran interface to <pastix_getSchurLocalNodeNbr>
 *
 * Compute the number of nodes in the local part of the Schur.
 *
 * Parameters:
 *   pastix_data - Common data structure for PaStiX calls.
 *   nodeNbr     - (out) Number of nodes in schur (local).
 *   ierr        - error value.
 *
 * Returns:
 *   NO_ERR      - For the moment
 *
 * TODO: Error management.
 */
void pastix_fortran_getschurlocalnodenbr ( void * pastix_data,
                                           PASTIX_INT * nodeNbr,
                                           PASTIX_INT * ierr )
{
  pastix_data_t **pastix_data2;
  pastix_data2 = pastix_data;
  *ierr = pastix_getSchurLocalNodeNbr(*pastix_data2,
                                      nodeNbr);

}

FORTRAN_NAME( PASTIX_FORTRAN_GETSCHURLOCALNODENBR,
              pastix_fortran_getschurlocalnodenbr,
              ( void * pastix_data,
                PASTIX_INT * nodeNbr,
                PASTIX_INT * ierr ),
              ( pastix_data,
                nodeNbr,
                ierr ))

/*
 * Function: pastix_fortran_getschurlocalunknownnbr
 *
 * Fortran interface to <pastix_getSchurLocalUnknownNbr>
 *
 * Compute the number of unknowns in the local part of the Schur.
 *
 * Parameters:
 *   pastix_data - Common data structure for PaStiX calls.
 *   unknownNbr  - (out) Number of nodes in schur (local).
 *
 * Returns:
 *   NO_ERR      - For the moment
 *
 * TODO: Error management.
 */
void pastix_fortran_getschurlocalunknownnbr ( void * pastix_data,
                                              PASTIX_INT  * unknownNbr,
                                              PASTIX_INT  * ierr )
{
  pastix_data_t **pastix_data2;
  pastix_data2 = pastix_data;
  *ierr = pastix_getSchurLocalUnkownNbr(*pastix_data2,
                                        unknownNbr);

}

FORTRAN_NAME( PASTIX_FORTRAN_GETSCHURLOCALUNKNOWNNBR,
              pastix_fortran_getschurlocalunknownnbr,
              ( void * pastix_data,
                PASTIX_INT  * unknownNbr,
                PASTIX_INT  * ierr ),
              ( pastix_data,
                unknownNbr,
                ierr ))
/*
 * Function: pastix_fortran_getschurlocalnodelist
 *
 * Fortran interface to <pastix_getSchurLocalNodeList>
 *
 * Compute the list of nodes in the local part of the Schur.
 *
 * Parameters:
 *   pastix_data - Common data structure for PaStiX calls.
 *   nodes       - (out) Nodes in schur (local).
 *   ierr        - error value.
 *
 * Returns:
 *   NO_ERR      - For the moment
 *
 * TODO: Error management.
 */
void pastix_fortran_getschurlocalnodelist ( void * pastix_data,
                                            PASTIX_INT  * nodes,
                                            PASTIX_INT  * ierr )
{
  pastix_data_t **pastix_data2;
  pastix_data2 = pastix_data;
  *ierr = pastix_getSchurLocalNodeList(*pastix_data2,
                                       nodes);
}

FORTRAN_NAME( PASTIX_FORTRAN_GETSCHURLOCALNODELIST,
              pastix_fortran_getschurlocalnodelist,
              ( void * pastix_data,
                PASTIX_INT  * nodes,
                PASTIX_INT  * ierr ),
              ( pastix_data,
                nodes,
                ierr ))
/*
 * Function: pastix_fortran_getschurlocalunkownlist
 *
 * Fortran interface to <pastix_getSchurLocalUnkownList>
 *
 * Compute the list of unknowns in the local part of the Schur.
 *
 * Parameters:
 *   pastix_data - Common data structure for PaStiX calls.
 *   unknowns    - (out) Unknowns in schur (local).
 *   ierr        - error value.
 *
 * Returns:
 *   NO_ERR      - For the moment
 *
 * TODO: Error management.
 */
void pastix_fortran_getschurlocalunknownlist (void * pastix_data,
                                              PASTIX_INT  * unknowns,
                                              PASTIX_INT  * ierr)
{
  pastix_data_t **pastix_data2;
  pastix_data2 = pastix_data;
  *ierr = pastix_getSchurLocalUnknownList(*pastix_data2,
                                          unknowns);
}

FORTRAN_NAME( PASTIX_FORTRAN_GETSCHURLOCALUNKNOWNLIST,
              pastix_fortran_getschurlocalunknownlist,
              (void * pastix_data,
               PASTIX_INT  * unknowns,
               PASTIX_INT  * ierr),
              (pastix_data,
               unknowns,
               ierr))
/*
 * Function: pastix_fortran_setschurarray
 *
 * Fortran interface to <pastix_setSchurArray>
 *
 * Give user memory area to store schur in PaStiX.
 *
 * Parameters:
 *   pastix_data - Common data structure for PaStiX calls.
 *   array       - Memory area to store the schur.
 *   ierr        - error value.
 *
 * Returns:
 *   NO_ERR      - For the moment
 *
 * TODO: Error management.
 */
void pastix_fortran_setschurarray ( void  * pastix_data,
                                    PASTIX_FLOAT * array,
                                    PASTIX_INT   * ierr )
{
  pastix_data_t **pastix_data2;
  pastix_data2 = pastix_data;
  *ierr = pastix_setSchurArray(*pastix_data2,
                               array);
}
FORTRAN_NAME( PASTIX_FORTRAN_SETSCHURARRAY,
              pastix_fortran_setschurarray,
              ( void  * pastix_data,
                PASTIX_FLOAT * array,
                PASTIX_INT   * ierr ),
              ( pastix_data,
                array,
                ierr ))
