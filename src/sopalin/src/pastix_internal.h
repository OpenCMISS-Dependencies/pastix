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
 * File: pastix_internal.h
 *
 * Header for function internal to pastix.c
 *
 */
#ifndef PASTIX_INTERNAL_H
#define PASTIX_INTERNAL_H
/*
 * Function: pastix_fake_fillin_csc
 *
 * Fill in the internal csc based on the user csc and fill in the coeftab structure
 *
 * Parameters:
 * pastix_data - PaStiX data structure.
 * pastix_comm - PaStiX MPI communicator.
 * n           - Size of the matrix.
 * colptr      - starting index of each column in row and avals.
 * row         - row of each element of the matrix.
 * avals       - value of each element of the matrix.
 * b           - Right hand side.
 * nrhs        - Number of right-hand-sides.
 * loc2glob    - global number of local columns, NULL if not distributed.
 */
#define pastix_fake_fillin_csc PASTIX_PREFIX_F(pastix_fake_fillin_csc)
int pastix_fake_fillin_csc( pastix_data_t *pastix_data,
                            MPI_Comm       pastix_comm,
                            PASTIX_INT            n,
                            PASTIX_INT           *colptr,
                            PASTIX_INT           *row,
                            PASTIX_FLOAT         *avals,
                            PASTIX_FLOAT         *b,
                            PASTIX_INT            nrhs,
                            PASTIX_INT           *loc2glob);

/*
 *  Function: pastix_fillin_csc
 *
 *  Fill in the internal csc based on the user csc and fill in the coeftab structure
 *
 *  Parameters:
 *  pastix_data - PaStiX data structure.
 *  pastix_comm - PaStiX MPI communicator.
 *  n           - Size of the matrix.
 *  colptr      - starting index of each column in row and avals.
 *  row         - row of each element of the matrix.
 *  avals       - value of each element of the matrix.
 *  b           - Right hand side.
 *  nrhs        - Number of right-hand-sides.
 *  loc2glob    - global number of local columns, NULL if not distributed.
 */
#define pastix_fillin_csc PASTIX_PREFIX_F(pastix_fillin_csc)
int pastix_fillin_csc( pastix_data_t *pastix_data,
                       MPI_Comm       pastix_comm,
                       PASTIX_INT            n,
                       PASTIX_INT           *colptr,
                       PASTIX_INT           *row,
                       PASTIX_FLOAT         *avals,
                       PASTIX_FLOAT         *b,
                       PASTIX_INT            nrhs,
                       PASTIX_INT           *loc2glob);


/*
 * Function: pastix_checkMatrix_int
 *
 * Check the matrix :
 * - Renumbers in Fortran numerotation (base 1) if needed (base 0)
 * - Check that the matrix contains no doubles,  with flagcor == API_YES,
 *   correct it.
 * - Can scale the matrix if compiled with -DMC64 -DSCALING (untested)
 * - Checks the symetry of the graph in non symmetric mode.
 *   With non distributed matrices, with flagcor == API_YES,
 *   correct the matrix.
 * - sort the CSC.
 *
 * Parameters:
 *   pastix_comm - PaStiX MPI communicator
 *   verb        - Level of prints (API_VERBOSE_[NOT|NO|YES])
 *   flagsym     - Indicate if the given matrix is symetric
 *                 (API_SYM_YES or API_SYM_NO)
 *   flagcor     - Indicate if we permit the function to reallocate the matrix.
 *   n           - Number of local columns.
 *   colptr      - First element of each row in *row* and *avals*.
 *   row         - Row of each element of the matrix.
 *   avals       - Value of each element of the matrix.
 *   loc2glob    - Global column number of local columns
 *                 (NULL if not distributed).
 *   dof         - Number of degrees of freedom.
 *   flagalloc   - indicate if allocation on CSC uses internal malloc.
 */
PASTIX_INT pastix_checkMatrix_int(MPI_Comm pastix_comm,
                           PASTIX_INT      verb,
                           PASTIX_INT      flagsym,
                           PASTIX_INT      flagcor,
                           PASTIX_INT      n,
                           PASTIX_INT    **colptr,
                           PASTIX_INT    **row,
                           PASTIX_FLOAT  **avals,
                           PASTIX_INT    **loc2glob,
                           PASTIX_INT      dof,
                           PASTIX_INT      flagalloc);

/*
 * Function: pastix_welcome_print
 *
 * Will print welcome message, options and parameters.
 *
 * Parameters:
 * pastix_data - PaStiX data structure
 * colptr      - starting index of each column in the CSC.
 * n           - number of columns.
 *
 */
void pastix_welcome_print(pastix_data_t *pastix_data,
                          PASTIX_INT           *colptr,
                          PASTIX_INT            ln);

/*
 * Function: pastix_task_clean
 *
 * Cleaning task
 *
 * Parameters:
 *
 */
void pastix_task_clean(pastix_data_t **pastix_data,
                       MPI_Comm        pastix_comm);

/*
 * Function: pastix_task_init
 *
 * Allocate and fill-in pastix_data
 *
 * Parameters:
 *   pastix_data - structure to build
 *   pastix_comm - PaStiX MPI communicator
 *   iparm       - integer parameters, to fill-in pastix_data
 *   dparm       - floating parameters, to fill-in pastix_data
 */
void pastix_task_init(pastix_data_t **pastix_data,
                      MPI_Comm        pastix_comm,
                      PASTIX_INT            *iparm,
                      double         *dparm);

/*
 * Function: pastix_initParam
 *
 * sets default parameters for iparm and dparm
 *
 * Parameters:
 * iparm - tabular of IPARM_SIZE integer parameters.
 * dparm - tabular of DPARM_SIZE double parameters.
 */
void pastix_initParam(PASTIX_INT    *iparm,
                      double *dparm);




#endif /* not PASTIX_INTERNAL_H */
