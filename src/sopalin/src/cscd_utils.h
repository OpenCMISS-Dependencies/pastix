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
  File: cscd_utils.h

  Several operations on CSCD.

 */

#ifndef CSCD_UTILS_H
#define CSCD_UTILS_H

#ifndef _GLIBCXX_HAVE_COMPLEX_H
#  define _GLIBCXX_HAVE_COMPLEX_H 0
#endif

#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || defined __COMPLEX_H__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX || defined _STLP_template_complex  || defined _LIBCPP_COMPLEX )
#  define PASTIX_HAS_COMPLEX
#endif

#ifdef PASTIX_HAS_COMPLEX
#  ifdef   __cplusplus
#    define  COMPLEX  std::complex<float>
#    define  DCOMPLEX std::complex<double>
#  else /* not __cplusplus */
#    define  COMPLEX float complex
#    define  DCOMPLEX double complex
#  endif /* not __cplusplus */
#endif

/*
 * MULTIPLE_TYPE_DEFINE
 *
 * Automaticaly generate function for each PASTIX_FLOAT type.
 *
 * This macro is fitted for function not taking floating points arguments.
 */
#ifdef PASTIX_HAS_COMPLEX
#define MULTIPLE_TYPE_DEFINE(functype, funcname, funcargs) \
  functype s_ ## funcname funcargs;                        \
  functype d_ ## funcname funcargs;                        \
  functype c_ ## funcname funcargs;                        \
  functype z_ ## funcname funcargs;
#else
#define MULTIPLE_TYPE_DEFINE(functype, funcname, funcargs)  \
  functype s_ ## funcname funcargs;                         \
  functype d_ ## funcname funcargs;
#endif

/*
 * MULTIPLE_TYPE_DEFINE_F
 *
 * Automaticaly generate function for each PASTIX_FLOAT type.
 *
 * This macro is fitted for function taking floating points arguments.
 */
#ifdef PASTIX_HAS_COMPLEX
#define MULTIPLE_TYPE_DEFINE_F(functype,        \
                               funcname,        \
                               funcargs_s,			\
                               funcargs_d,			\
                               funcargs_c,			\
                               funcargs_z)			\
  functype s_ ## funcname funcargs_s;           \
  functype d_ ## funcname funcargs_d;           \
  functype c_ ## funcname funcargs_c;           \
  functype z_ ## funcname funcargs_z;
#else
#define MULTIPLE_TYPE_DEFINE_F(functype,        \
                               funcname,        \
                               funcargs_s,			\
                               funcargs_d,			\
                               funcargs_c,			\
                               funcargs_z)			\
  functype s_ ## funcname funcargs_s;           \
  functype d_ ## funcname funcargs_d;
#endif

/*
 * Enum: CSCD_OPERATIONS
 *
 * Operation when adding CSCD
 *
 * CSCD_ADD  - Add coefficient values.
 * CSCD_KEEP - Keep value from first CSCD.
 * CSCD_MAX  - Keep maximum of first and second CSCD.
 * CSCD_MIN  - Keep minimum of first and second CSCD.
 * CSCD_OVW  - Overwrite with second CSCD value.
 */
enum CSCD_OPERATIONS {
  CSCD_ADD,
  CSCD_KEEP,
  CSCD_MAX,
  CSCD_MIN,
  CSCD_OVW
};
typedef enum CSCD_OPERATIONS CSCD_OPERATIONS_t;

/*
 * Enum: CSC_DISPATCH_OP
 *
 * Operation when dispatching the csc into a cscd
 *
 * CSC_DISP_SIMPLE - Reparts linearly the columns over the proc.
 * CSC_DISP_CYCLIC - Reparts cyclicly the columns over the proc.
 *
 */
enum CSC_DISPATCH_OP {
  CSC_DISP_SIMPLE,
  CSC_DISP_CYCLIC
};
typedef enum CSC_DISPATCH_OP CSCD_DISPATCH_OP_t;

/* Section: Functions */
#if (defined PASTIX_FLOAT)
/*
 *  Function: csc_dispatch
 *
 *  Distribute a CSC to a CSCD
 *
 *  Parameters:
 *     gN                - global number of columns
 *     gcolptr           - global starting index of each column in grows ans gavals.
 *     grows             - global rows of each element.
 *     gavals            - global values of each element.
 *     gperm             - global permutation tabular.
 *     ginvp             - global reverse permutation tabular.
 *     lN                - local number of columns (output).
 *     lcolptr           - starting index of each local column (output).
 *     lrowptr           - row number of each local element (output).
 *     lavals            - values of each local element (output).
 *     lrhs              - local part of the right hand side (output).
 *     lperm             - local part of the permutation tabular (output).
 *     loc2glob          - global numbers of local columns (before permutation).
 *     dispatch          - choose how to dispatch the csc
 *     pastix_comm       - PaStiX MPI communicator.
 */
void csc_dispatch(PASTIX_INT  gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, PASTIX_FLOAT *  gavals,
                  PASTIX_FLOAT *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
                  PASTIX_INT *lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, PASTIX_FLOAT ** lavals,
                  PASTIX_FLOAT ** lrhs, PASTIX_INT ** lperm,
                  PASTIX_INT **loc2glob, int dispatch, MPI_Comm pastix_comm);
#endif
MULTIPLE_TYPE_DEFINE_F(void,
                       csc_dispatch,
                       (PASTIX_INT  gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, float *  gavals,
                        float *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
                        PASTIX_INT *lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, float ** lavals,
                        float ** lrhs, PASTIX_INT ** lperm,
                        PASTIX_INT **loc2glob, int dispatch, MPI_Comm pastix_comm),
                       (PASTIX_INT  gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, double *  gavals,
                        double *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
                        PASTIX_INT *lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, double ** lavals,
                        double ** lrhs, PASTIX_INT ** lperm,
                        PASTIX_INT **loc2glob, int dispatch, MPI_Comm pastix_comm),
                       (PASTIX_INT  gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, COMPLEX *  gavals,
                        COMPLEX *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
                        PASTIX_INT *lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, COMPLEX ** lavals,
                        COMPLEX ** lrhs, PASTIX_INT ** lperm,
                        PASTIX_INT **loc2glob, int dispatch, MPI_Comm pastix_comm),
                       (PASTIX_INT  gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, DCOMPLEX *  gavals,
                        DCOMPLEX *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
                        PASTIX_INT *lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, DCOMPLEX ** lavals,
                        DCOMPLEX ** lrhs, PASTIX_INT ** lperm,
                        PASTIX_INT **loc2glob, int dispatch, MPI_Comm pastix_comm))

#if (defined PASTIX_FLOAT)
/*
 * Function: csc_cyclic_distribution
 *
 * Distribute the CSC cyclicaly.
 *
 * Parameters:
 *   column      - column number to distribute
 *   columnnbr   - Number of colmuns.
 *   pastix_comm - PaStiX MPI communicator
 *
 * Return:
 *   owner of the column (column%commSize)
 */
PASTIX_INT csc_cyclic_distribution(PASTIX_INT column, PASTIX_INT columnnbr, MPI_Comm pastix_comm);
#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, csc_cyclic_distribution,
                     (PASTIX_INT column, PASTIX_INT columnnbr, MPI_Comm pastix_comm))

#if (defined PASTIX_FLOAT)
/*
 * Function: csc_simple_distribution
 *
 * Distribute the CSC.
 * First columns are for first proc and so on.
 *
 * Parameters:
 *   column      - column number to distribute
 *   columnnbr   - Number of colmuns.
 *   pastix_comm - PaStiX MPI communicator
 *
 * Return:
 *   owner of the column (column/commSize)
 */
PASTIX_INT csc_simple_distribution(PASTIX_INT column, PASTIX_INT columnnbr, MPI_Comm pastix_comm);

#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, csc_simple_distribution,
                     (PASTIX_INT column, PASTIX_INT columnnbr, MPI_Comm pastix_comm))

#if (defined PASTIX_FLOAT)
/*
 * Function: cscd_symgraph
 *
 * Check if the CSCD graph is symetric.
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - Starting index of each columns in *ja* and *a*
 *   ja          - Row of each element.
 *   a           - Values of each element.
 *   newn        - New number of local columns
 *   newia       - Starting index of each columns in *newja* and *newa*
 *   newja       - Row of each element.
 *   newa        - Values of each element.
 *   l2g         - global number of each local column.
 *   malloc_flag - flag to indicate if function call is intern to pastix or extern.
 */
int cscd_symgraph(PASTIX_INT      n, PASTIX_INT *     ia, PASTIX_INT *     ja, PASTIX_FLOAT *     a,
                  PASTIX_INT * newn, PASTIX_INT ** newia, PASTIX_INT ** newja, PASTIX_FLOAT ** newa,
                  PASTIX_INT *     l2g,  MPI_Comm comm);
#endif
MULTIPLE_TYPE_DEFINE(int, cscd_symgraph,
                     (PASTIX_INT      n, PASTIX_INT *     ia, PASTIX_INT *     ja, float *     a,
                      PASTIX_INT * newn, PASTIX_INT ** newia, PASTIX_INT ** newja, float ** newa,
                      PASTIX_INT *     l2g,  MPI_Comm comm))

#if (defined PASTIX_FLOAT)
/*
 * Function: cscd_addlocal
 *
 * Add second cscd to first cscd into third cscd (unallocated)
 *
 * Parameters:
 *   n           - First cscd size
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   l2g         - local 2 global column numbers for first cscd
 *   addn        - CSCD to add size
 *   addia       - CSCD to add starting index of each column in *addja* and *adda*
 *   addja       - Row of each element in second CSCD
 *   adda        - value of each cscd in second CSCD (can be NULL -> add 0)
 *   addl2g      - local 2 global column numbers for second cscd
 *   newn        - new cscd size (same as first)
 *   newia       - CSCD to add starting index of each column in *newja* and *newwa*
 *   newja       - Row of each element in third CSCD
 *   newa        - value of each cscd in third CSCD
 *   OP          - Operation to manage common CSCD coefficients.
 *   dof         - Number of degrees of freedom.
 */

int cscd_addlocal(PASTIX_INT   n   , PASTIX_INT *  ia   , PASTIX_INT *  ja   , PASTIX_FLOAT *  a   , PASTIX_INT * l2g,
      PASTIX_INT   addn, PASTIX_INT *  addia, PASTIX_INT *  addja, PASTIX_FLOAT *  adda, PASTIX_INT * addl2g,
      PASTIX_INT * newn, PASTIX_INT ** newia, PASTIX_INT ** newja, PASTIX_FLOAT ** newa, CSCD_OPERATIONS_t OP, int dof);
#endif
MULTIPLE_TYPE_DEFINE_F(int, cscd_addlocal,
                       (PASTIX_INT   n   , PASTIX_INT *  ia   , PASTIX_INT *  ja   , float *  a   , PASTIX_INT * l2g,
                        PASTIX_INT   addn, PASTIX_INT *  addia, PASTIX_INT *  addja, float *  adda, PASTIX_INT * addl2g,
                        PASTIX_INT * newn, PASTIX_INT ** newia, PASTIX_INT ** newja, float ** newa, CSCD_OPERATIONS_t OP, int dof),
                       (PASTIX_INT   n   , PASTIX_INT *  ia   , PASTIX_INT *  ja   , double *  a   , PASTIX_INT * l2g,
                        PASTIX_INT   addn, PASTIX_INT *  addia, PASTIX_INT *  addja, double *  adda, PASTIX_INT * addl2g,
                        PASTIX_INT * newn, PASTIX_INT ** newia, PASTIX_INT ** newja, double ** newa, CSCD_OPERATIONS_t OP, int dof),
                       (PASTIX_INT   n   , PASTIX_INT *  ia   , PASTIX_INT *  ja,
                        COMPLEX *  a   , PASTIX_INT * l2g,
                        PASTIX_INT   addn, PASTIX_INT *  addia, PASTIX_INT *  addja,
                        COMPLEX *  adda, PASTIX_INT * addl2g,
                        PASTIX_INT * newn, PASTIX_INT ** newia, PASTIX_INT ** newja,
                        COMPLEX ** newa, CSCD_OPERATIONS_t OP, int dof),
                       (PASTIX_INT   n   , PASTIX_INT *  ia   , PASTIX_INT *  ja,
                        DCOMPLEX *  a   , PASTIX_INT * l2g,
                        PASTIX_INT   addn, PASTIX_INT *  addia, PASTIX_INT *  addja,
                        DCOMPLEX *  adda, PASTIX_INT * addl2g,
                        PASTIX_INT * newn, PASTIX_INT ** newia, PASTIX_INT ** newja,
                        DCOMPLEX ** newa, CSCD_OPERATIONS_t OP, int dof))


#if (defined PASTIX_FLOAT)
/**
 *   Function: csc2cscd
 *
 *   Transform a csc to a cscd.
 *   Allocate the CSCD.
 *   If grhs == NULL forget right hand side part.
 *   If gperm == NULL forget permutation and reverse permutation part.
 *
 *   Parameters:
 *     gN       - global number of columns
 *     gcolptr  - global starting index of each column in grows ans gavals.
 *     grows    - global rows of each element.
 *     gavals   - global values of each element.
 *     gperm    - global permutation tabular.
 *     ginvp    - global reverse permutation tabular.
 *     lN       - local number of columns.
 *     lcolptr  - starting index of each local column.
 *     lrowptr  - row number of each local element.
 *     lavals   - values of each local element.
 *     lrhs     - local part of the right hand side (output).
 *     lperm    - local part of the permutation tabular (output).
 *     linvp    - local part of the reverse permutation tabular (output).
 *     loc2glob - global numbers of local columns (before permutation).
 */
void  csc2cscd(PASTIX_INT gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, PASTIX_FLOAT *  gavals,
               PASTIX_FLOAT *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
               PASTIX_INT lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, PASTIX_FLOAT ** lavals,
               PASTIX_FLOAT ** lrhs, PASTIX_INT ** lperm, PASTIX_INT ** linvp,
               PASTIX_INT *loc2glob);
#endif
MULTIPLE_TYPE_DEFINE_F(void,  csc2cscd,
                       (PASTIX_INT gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, float *  gavals,
                        float *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
                        PASTIX_INT lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, float ** lavals,
                        float ** lrhs, PASTIX_INT ** lperm, PASTIX_INT ** linvp,
                        PASTIX_INT *loc2glob),
                       (PASTIX_INT gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, double *  gavals,
                        double *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
                        PASTIX_INT lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, double ** lavals,
                        double ** lrhs, PASTIX_INT ** lperm, PASTIX_INT ** linvp,
                        PASTIX_INT *loc2glob),
                       (PASTIX_INT gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, COMPLEX *  gavals,
                        COMPLEX *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
                        PASTIX_INT lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, COMPLEX ** lavals,
                        COMPLEX ** lrhs, PASTIX_INT ** lperm, PASTIX_INT ** linvp,
                        PASTIX_INT *loc2glob),
                       (PASTIX_INT gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, DCOMPLEX *  gavals,
                        DCOMPLEX *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
                        PASTIX_INT lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, DCOMPLEX ** lavals,
                        DCOMPLEX ** lrhs, PASTIX_INT ** lperm, PASTIX_INT ** linvp,
                        PASTIX_INT *loc2glob))

#if (defined PASTIX_FLOAT)
/**
 *   Function: cscd2csc
 *
 *   Transform a cscd to a csc.
 *   colptr2, row2, avals2, rhs2, perm2, invp2 are allocated here.
 *
 *   Parameters:
 *      lN          - number of local column.
 *      lcolptr     - starting index of each local column in row and avals.
 *      lrow        _ row number of each local element.
 *      lavals      - values of each local element.
 *      lrhs        - local part of the right hand side.
 *      lperm       - local part of the permutation tabular.
 *      linvp       - local part of the reverse permutation tabular.
 *      gN          - global number of columns (output).
 *      gcolptr     - starting index of each column in row2 and avals2 (output).
 *      grow        - row number of each element (output).
 *      gavals      - values of each element (output).
 *      grhs        - global right hand side (output).
 *      gperm       - global permutation tabular (output).
 *      ginvp       - global reverse permutation tabular (output).
 *      loc2glob    - global number of each local column.
 *      pastix_comm - PaStiX MPI communicator.
 *      ndof        - Number of degree f freedom by node.
 *
 */

void  cscd2csc(PASTIX_INT  lN, PASTIX_INT *  lcolptr, PASTIX_INT * lrow, PASTIX_FLOAT * lavals,
               PASTIX_FLOAT * lrhs, PASTIX_INT * lperm, PASTIX_INT * linvp,
               PASTIX_INT *gN, PASTIX_INT ** gcolptr, PASTIX_INT **grow, PASTIX_FLOAT **gavals,
               PASTIX_FLOAT **grhs, PASTIX_INT **gperm, PASTIX_INT **ginvp,
               PASTIX_INT *loc2glob, MPI_Comm pastix_comm, PASTIX_INT ndof);
#endif
MULTIPLE_TYPE_DEFINE_F(void,  cscd2csc,
                       (PASTIX_INT  lN, PASTIX_INT *  lcolptr, PASTIX_INT * lrow, float * lavals,
                        float * lrhs, PASTIX_INT * lperm, PASTIX_INT * linvp,
                        PASTIX_INT *gN, PASTIX_INT ** gcolptr, PASTIX_INT **grow, float **gavals,
                        float **grhs, PASTIX_INT **gperm, PASTIX_INT **ginvp,
                        PASTIX_INT *loc2glob, MPI_Comm pastix_comm, PASTIX_INT ndof),
                       (PASTIX_INT  lN, PASTIX_INT *  lcolptr, PASTIX_INT * lrow, double * lavals,
                        double * lrhs, PASTIX_INT * lperm, PASTIX_INT * linvp,
                        PASTIX_INT *gN, PASTIX_INT ** gcolptr, PASTIX_INT **grow, double **gavals,
                        double **grhs, PASTIX_INT **gperm, PASTIX_INT **ginvp,
                        PASTIX_INT *loc2glob, MPI_Comm pastix_comm, PASTIX_INT ndof),
                       (PASTIX_INT  lN, PASTIX_INT *  lcolptr, PASTIX_INT * lrow, COMPLEX * lavals,
                        COMPLEX * lrhs, PASTIX_INT * lperm, PASTIX_INT * linvp,
                        PASTIX_INT *gN, PASTIX_INT ** gcolptr, PASTIX_INT **grow, COMPLEX **gavals,
                        COMPLEX **grhs, PASTIX_INT **gperm, PASTIX_INT **ginvp,
                        PASTIX_INT *loc2glob, MPI_Comm pastix_comm, PASTIX_INT ndof),
                       (PASTIX_INT  lN, PASTIX_INT *  lcolptr, PASTIX_INT * lrow, DCOMPLEX * lavals,
                        DCOMPLEX * lrhs, PASTIX_INT * lperm, PASTIX_INT * linvp,
                        PASTIX_INT *gN, PASTIX_INT ** gcolptr, PASTIX_INT **grow, DCOMPLEX **gavals,
                        DCOMPLEX **grhs, PASTIX_INT **gperm, PASTIX_INT **ginvp,
                        PASTIX_INT *loc2glob, MPI_Comm pastix_comm, PASTIX_INT ndof))

#if (defined PASTIX_FLOAT)
/*
 * Function: cscd_redispatch
 *
 * Redistribute the first cscd into a new one using *dl2g*.
 *
 * - gather all new loc2globs on all processors.
 * - allocate *dia*, *dja* and *da*.
 * - Create new CSC for each processor and send it.
 * - Merge all new CSC to the new local CSC with <cscd_addlocal_int>.
 *
 * If communicator size is one, check that n = dn and
 * l2g = dl2g and simply create a copy of the first cscd.
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   rhs         - right-hand-side member corresponding to the first CSCD
 *                 (can be NULL)
 *   nrhs        - number of right-hand-side.
 *   l2g         - local 2 global column numbers for first cscd
 *   dn          - Number of local columns
 *   dia         - New cscd starting index of each column in *ja* and *a*
 *   dja         - Row of each element in new CSCD
 *   da          - value of each cscd in new CSCD
 *   rhs         - right-hand-side member corresponding to the new CSCD
 *   dl2g        - local 2 global column numbers for new cscd
 *   comm        - MPI communicator
 *
 * Returns:
 *   EXIT_SUCCESS - If all goes well
 *   EXIT_FAILURE - If commsize = 1 and *n* != *dn* or *l2g* != *dl2g*.
 */
int cscd_redispatch(PASTIX_INT   n, PASTIX_INT *   ia, PASTIX_INT *   ja, PASTIX_FLOAT *   a,
                    PASTIX_FLOAT *  rhs,  PASTIX_INT nrhs, PASTIX_INT *   l2g,
                    PASTIX_INT  dn, PASTIX_INT ** dia, PASTIX_INT ** dja, PASTIX_FLOAT ** da,
                    PASTIX_FLOAT ** drhs,  PASTIX_INT *  dl2g,
                    MPI_Comm comm, PASTIX_INT dof);
#endif
MULTIPLE_TYPE_DEFINE_F(int, cscd_redispatch,
                       (PASTIX_INT   n, PASTIX_INT *   ia, PASTIX_INT *   ja, float *   a,
                        float *  rhs,  PASTIX_INT nrhs,   PASTIX_INT *   l2g,
                        PASTIX_INT  dn, PASTIX_INT ** dia, PASTIX_INT ** dja, float ** da,
                        float ** drhs,  PASTIX_INT *  dl2g,
                        MPI_Comm comm, PASTIX_INT dof),
                       (PASTIX_INT   n, PASTIX_INT *   ia, PASTIX_INT *   ja, double *   a,
                        double *  rhs,  PASTIX_INT nrhs,   PASTIX_INT *   l2g,
                        PASTIX_INT  dn, PASTIX_INT ** dia, PASTIX_INT ** dja, double ** da,
                        double ** drhs,  PASTIX_INT *  dl2g,
                        MPI_Comm comm, PASTIX_INT dof),
                       (PASTIX_INT   n, PASTIX_INT *   ia, PASTIX_INT *   ja, COMPLEX *   a,
                        COMPLEX *  rhs,  PASTIX_INT nrhs,   PASTIX_INT *   l2g,
                        PASTIX_INT  dn, PASTIX_INT ** dia, PASTIX_INT ** dja, COMPLEX ** da,
                        COMPLEX ** drhs,  PASTIX_INT *  dl2g,
                        MPI_Comm comm, PASTIX_INT dof),
                       (PASTIX_INT   n, PASTIX_INT *   ia, PASTIX_INT *   ja, DCOMPLEX *   a,
                        DCOMPLEX *  rhs,  PASTIX_INT nrhs,   PASTIX_INT *   l2g,
                        PASTIX_INT  dn, PASTIX_INT ** dia, PASTIX_INT ** dja, DCOMPLEX ** da,
                        DCOMPLEX ** drhs,  PASTIX_INT *  dl2g,
                        MPI_Comm comm, PASTIX_INT dof))

#if (defined PASTIX_FLOAT)
/*
 *  Function: cscd_save
 *
 *  save a distributed csc to disk.
 *  files are called $(filename) and $(filename)$(RANK)
 *  if filename is NULL then filename = cscd_matrix.
 *
 *  file filename contains the number of processors/files
 *  on first line. Then each line contain the name of each file
 *  (here $(filename)$(RANK)).
 *
 *
 *
 *  Parameters:
 *    n           - Number of local columns
 *    ia          - First cscd starting index of each column in *ja* and *a*
 *    ja          - Row of each element in first CSCD
 *    a           - value of each cscd in first CSCD (can be NULL)
 *    rhs         - Right hand side.
 *    l2g         - local 2 global column numbers for first cscd
 *    dof         - Number of degrees of freedom
 *    filename    - name of the files.
 *    comm        - MPI communicator
 */
int cscd_save(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT * a, PASTIX_FLOAT * rhs, PASTIX_INT* l2g,
              int dof, const char * filename, MPI_Comm comm);
#endif
MULTIPLE_TYPE_DEFINE_F(int, cscd_save,
                       (PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, float * a, float * rhs, PASTIX_INT* l2g,
                        int dof, const char * filename, MPI_Comm comm),
                       (PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, double * a, double * rhs, PASTIX_INT* l2g,
                        int dof, const char * filename, MPI_Comm comm),
                       (PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, COMPLEX * a, COMPLEX * rhs,
                        PASTIX_INT* l2g,  int dof,  const char * filename, MPI_Comm comm),
                       (PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, DCOMPLEX * a, DCOMPLEX * rhs,
                        PASTIX_INT* l2g,  int dof,  const char * filename, MPI_Comm comm))

#if (defined PASTIX_FLOAT)
/*
 *  Function: cscd_load
 *
 *  Loads a distributed csc from disk.
 *  if filename is NULL then filename = cscd_matrix.
 *
 *  Parameters:
 *    n           - Number of local columns
 *    ia          - First cscd starting index of each column in *ja* and *a*
 *    ja          - Row of each element in first CSCD
 *    a           - value of each cscd in first CSCD (can be NULL)
 *    rhs         - Right hand side.
 *    l2g         - local 2 global column numbers for first cscd
 *    filename    - name of the files.
 *    comm        - MPI communicator
 */
int cscd_load(PASTIX_INT *n, PASTIX_INT ** ia, PASTIX_INT ** ja, PASTIX_FLOAT ** a, PASTIX_FLOAT ** rhs, PASTIX_INT ** l2g,
              const char * filename, MPI_Comm mpi_comm);
#endif
MULTIPLE_TYPE_DEFINE_F(int, cscd_load,
                       (PASTIX_INT *n, PASTIX_INT ** ia, PASTIX_INT ** ja, float ** a, float ** rhs, PASTIX_INT ** l2g,
                        const char * filename, MPI_Comm mpi_comm),
                       (PASTIX_INT *n, PASTIX_INT ** ia, PASTIX_INT ** ja, double ** a, double ** rhs,
                        PASTIX_INT ** l2g, const char * filename, MPI_Comm mpi_comm),
                       (PASTIX_INT *n, PASTIX_INT ** ia, PASTIX_INT ** ja, COMPLEX ** a,
                        COMPLEX ** rhs, PASTIX_INT ** l2g, const char * filename,
                        MPI_Comm mpi_comm),
                       (PASTIX_INT *n, PASTIX_INT ** ia, PASTIX_INT ** ja, DCOMPLEX ** a,
                        DCOMPLEX ** rhs, PASTIX_INT ** l2g, const char * filename,
                        MPI_Comm mpi_comm))

#undef MULTIPLE_TYPE_DEFINE_F
#undef MULTIPLE_TYPE_DEFINE
#undef COMPLEX
#undef DCOMPLEX
#endif /* CSCD_UTILS_H */
