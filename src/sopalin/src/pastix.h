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
 * Title: PaStiX interface
 *
 *   Parallel Sparse matriX solver interface
 *
 * Authors:
 *   Astrid  CASADEI  - astrid.casadei@inria.fr
 *   Mathieu FAVERGE  - faverge@labri.fr
 *   Xavier  LACOSTE  - xavier.lacoste@inria.fr
 *   Pierre  RAMET    - ramet@labri.fr
 *
 */

#ifndef _PASTIX_H_
#define _PASTIX_H_

#define PASTIX_MAJOR_VERSION   5
#define PASTIX_MEDIUM_VERSION  2
#define PASTIX_MINOR_VERSION   2
#define PASTIX_BUGFIX_VERSION  16

struct pastix_data_t;
typedef struct pastix_data_t pastix_data_t;

#ifndef _GLIBCXX_HAVE_COMPLEX_H
#  define _GLIBCXX_HAVE_COMPLEX_H 0
#endif

#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || defined __COMPLEX_H__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX || defined _STLP_template_complex || defined _LIBCPP_COMPLEX)
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

#if (defined PASTIX_FLOAT)

/*
 * Group: Main PaStiX functions
 */
/*
 * Function: pastix
 *
 * Computes steps of the resolution of Ax=b linear system,
 * using direct methods.
 *
 * The matrix is given in CSC format.
 *
 * Parameters:
 *   pastix_data - Data used for a step by step execution.
 *   pastix_comm - MPI communicator which compute the resolution.
 *   n           - Size of the system.
 *   colptr      - Tabular containing the start of each column in row
 *                 and avals tabulars.
 *   row         - Tabular containing the row number for each element
 *                 sorted by column.
 *   avals       - Tabular containing the values of each element
 *                 sorted by column.
 *   perm        - Permutation tabular for the renumerotation of the unknowns.
 *   invp        - Reverse permutation tabular for the renumerotation
 *                 of the unknowns.
 *   b           - Right hand side vector(s).
 *   rhs         - Number of right hand side vector(s).
 *   iparm       - Integer parameters given to pastix.
 *   dparm       - Double parameters given to pâstix.
 *
 * About: Example
 *
 *   from file <simple.c> :
 *
 *   > /\*******************************************\/
 *   > /\*    Check Matrix format                  *\/
 *   > /\*******************************************\/
 *   > /\*
 *   >  * Matrix needs :
 *   >  *    - to be in fortran numbering
 *   >  *    - to have only the lower triangular part in symmetric case
 *   >  *    - to have a graph with a symmetric structure in unsymmetric case
 *   >  *\/
 *   > mat_type = API_SYM_NO;
 *   > if (MTX_ISSYM(type)) mat_type = API_SYM_YES;
 *   > if (MTX_ISHER(type)) mat_type = API_SYM_HER;
 *   > pastix_checkMatrix( MPI_COMM_WORLD, verbosemode,
 *   >                     mat_sym,
 *   >                     API_YES,
 *   >                     ncol, &colptr, &rows, &values, NULL);
 *   >
 *   > /\*******************************************\/
 *   > /\* Initialize parameters to default values *\/
 *   > /\*******************************************\/
 *   > iparm[IPARM_MODIFY_PARAMETER] = API_NO;
 *   > pastix(&pastix_data, MPI_COMM_WORLD,
 *   >        ncol, colptr, rows, values,
 *   >        perm, invp, rhs, 1, iparm, dparm);
 *   >
 *   > /\*******************************************\/
 *   > /\*       Customize some parameters         *\/
 *   > /\*******************************************\/
 *   > iparm[IPARM_THREAD_NBR] = nbthread;
 *   > iparm[IPARM_SYM] = mat_type;
 *   > switch (mat_type)
 *   >   {
 *   >     case API_SYM_YES:
 *   >       iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
 *   >       break;
 *   >     case API_SYM_HER:
 *   >       iparm[IPARM_FACTORIZATION] = API_FACT_LDLH;
 *   >       break;
 *   >     default:
 *   >       iparm[IPARM_FACTORIZATION] = API_FACT_LU;
 *   >   }
 *   > iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
 *   > iparm[IPARM_END_TASK]            = API_TASK_CLEAN;
 *   >
 *   > /\*******************************************\/
 *   > /\*           Save the rhs                  *\/
 *   > /\*    (it will be replaced by solution)    *\/
 *   > /\*******************************************\/
 *   > rhssaved = malloc(ncol*sizeof(pastix_float_t));
 *   > memcpy(rhssaved, rhs, ncol*sizeof(pastix_float_t));
 *   >
 *   > /\*******************************************\/
 *   > /\*           Call pastix                   *\/
 *   > /\*******************************************\/
 *   > perm = malloc(ncol*sizeof(pastix_int_t));
 *   > invp = malloc(ncol*sizeof(pastix_int_t));
 *   >
 *   > pastix(&pastix_data, MPI_COMM_WORLD,
 *   >  ncol, colptr, rows, values,
 *   >  perm, invp, rhs, 1, iparm, dparm);
 */
void pastix(pastix_data_t **pastix_data, MPI_Comm pastix_comm,
            PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
            PASTIX_FLOAT *avals, PASTIX_INT *perm, PASTIX_INT *invp, PASTIX_FLOAT *b, PASTIX_INT rhs,
            PASTIX_INT *iparm, double *dparm);
#endif
MULTIPLE_TYPE_DEFINE_F(void, pastix,
                       ( pastix_data_t **pastix_data, MPI_Comm pastix_comm,
                         PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
                         float *avals, PASTIX_INT *perm, PASTIX_INT *invp, float *b, PASTIX_INT rhs,
                         PASTIX_INT *iparm, double *dparm),

                       ( pastix_data_t **pastix_data, MPI_Comm pastix_comm,
                         PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
                         double *avals, PASTIX_INT *perm, PASTIX_INT *invp,
                         double *b, PASTIX_INT rhs,
                         PASTIX_INT *iparm, double *dparm),

                       ( pastix_data_t **pastix_data, MPI_Comm pastix_comm,
                         PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
                         COMPLEX *avals, PASTIX_INT *perm, PASTIX_INT *invp,
                         COMPLEX *b, PASTIX_INT rhs, PASTIX_INT *iparm, double *dparm),

                       ( pastix_data_t **pastix_data, MPI_Comm pastix_comm,
                         PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
                         DCOMPLEX *avals, PASTIX_INT *perm, PASTIX_INT *invp,
                         DCOMPLEX *b, PASTIX_INT rhs, PASTIX_INT *iparm, double *dparm))

#if (defined PASTIX_FLOAT)
/*
 * Function: dpastix
 *
 *   Computes steps of the resolution of Ax=b linear system,
 *   using direct methods.
 *   Here the matrix is given distributed.
 *
 *   The matrix is given in CSCD format.
 *
 *   Parameters:
 *      pastix_data - Data used for a step by step execution.
 *      pastix_comm - MPI communicator which compute the resolution.
 *      n           - Size of the system.
 *      colptr      - Tabular containing the start of each column in
 *                    *row* and *avals* tabulars.
 *      row         - Tabular containing the row number for each element
 *                    sorted by column.
 *      avals       - Tabular containing the values of each element
 *                    sorted by column.
 *      loc2glob    - Global column number of the local columns.
 *      perm        - Permutation tabular for the renumerotation
 *                    of the unknowns.
 *      invp        - Reverse permutation tabular for the renumerotation
 *                    of the unknowns.
 *      b           - Right hand side vector(s).
 *      rhs         - Number of right hand side vector(s).
 *      iparm       - Integer parameters given to pastix.
 *      dparm       - Double parameters given to pâstix.
 *
 * About: Example
 *
 *   from file <simple_dist.c> :
 *
 *   > /\*******************************************\/
 *   > /\*    Check Matrix format                  *\/
 *   > /\*******************************************\/
 *   > /\*
 *   >  * Matrix needs :
 *   >  *    - to be in fortran numbering
 *   >  *    - to have only the lower triangular part in symmetric case
 *   >  *    - to have a graph with a symmetric structure in unsymmetric case
 *   >  *\/
 *   > pastix_checkMatrix( MPI_COMM_WORLD, verbosemode,
 *   >                     (MTX_ISSYM(type) ? API_SYM_YES : API_SYM_NO),
 *   >                     API_YES,
 *   >                     ncol, &colptr, &rows, &values, &loc2glob, 1);
 *   >
 *   > /\*******************************************\/
 *   > /\* Initialize parameters to default values *\/
 *   > /\*******************************************\/
 *   > iparm[IPARM_MODIFY_PARAMETER] = API_NO;
 *   > dpastix(&pastix_data, MPI_COMM_WORLD,
 *   >         ncol, colptr, rows, values, loc2glob,
 *   >         perm, invp, rhs, 1, iparm, dparm);
 *   >
 *   > /\*******************************************\/
 *   > /\*       Customize some parameters         *\/
 *   > /\*******************************************\/
 *   > iparm[IPARM_THREAD_NBR] = nbthread;
 *   > if (MTX_ISSYM(type))
 *   >   {
 *   >     iparm[IPARM_SYM]           = API_SYM_YES;
 *   >     iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
 *   >   }
 *   > else
 *   >   {
 *   >     iparm[IPARM_SYM]           = API_SYM_NO;
 *   >     iparm[IPARM_FACTORIZATION] = API_FACT_LU;
 *   >   }
 *   > iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
 *   >
 *   > iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
 *   > iparm[IPARM_END_TASK]            = API_TASK_BLEND;
 *   >
 *   > /\*******************************************\/
 *   > /\*           Call pastix                   *\/
 *   > /\*******************************************\/
 *   > perm = malloc(ncol*sizeof(pastix_int_t));
 *   > /\* No need to allocate invp in dpastix *\/
 *   >
 *   > dpastix(&pastix_data, MPI_COMM_WORLD,
 *   >         ncol, colptr, rows, NULL, loc2glob,
 *   >         perm, NULL, NULL, 1, iparm, dparm);
 *   >
 *   > /\* Redistributing the matrix *\/
 *   >
 *   > ncol2 = pastix_getLocalNodeNbr(&pastix_data);
 *   >
 *   > if (NULL == (loc2glob2 = malloc(ncol2 * sizeof(pastix_int_t))))
 *   >   {
 *   >     fprintf(stderr, "Malloc error\n");
 *   >     return EXIT_FAILURE;
 *   >   }
 *   >
 *   > pastix_getLocalNodeLst(&pastix_data, loc2glob2);
 *   >
 *   > if (EXIT_SUCCESS != cscd_redispatch(ncol,   colptr,   rows,
 *   >                                     values,    rhs,  loc2glob,
 *   >                                     ncol2, &colptr2, &rows2,
 *   >                                     &values2, &rhs2, loc2glob2,
 *   >                                     MPI_COMM_WORLD))
 *   >   return EXIT_FAILURE;
 *   >
 *   > free(colptr);
 *   > free(rows);
 *   > free(values);
 *   > free(rhs);
 *   > free(loc2glob);
 *   > free(perm);
 *   >
 *   > iparm[IPARM_START_TASK]          = API_TASK_NUMFACT;
 *   > iparm[IPARM_END_TASK]            = API_TASK_CLEAN;
 *   >
 *   >
 *   > dpastix(&pastix_data, MPI_COMM_WORLD,
 *   >         ncol2, colptr2, rows2, values2, loc2glob2,
 *   >         perm, invp, rhs2, 1, iparm, dparm);
 *   >
 */
void dpastix(pastix_data_t **pastix_data, MPI_Comm pastix_comm,
             PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
             PASTIX_FLOAT *avals, PASTIX_INT * loc2glob, PASTIX_INT *perm, PASTIX_INT *invp,
             PASTIX_FLOAT *b, PASTIX_INT rhs, PASTIX_INT *iparm, double *dparm);
#endif
MULTIPLE_TYPE_DEFINE_F(void, dpastix,
                       (pastix_data_t **pastix_data, MPI_Comm pastix_comm,
                        PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
                        float *avals, PASTIX_INT * loc2glob, PASTIX_INT *perm, PASTIX_INT *invp,
                        float *b, PASTIX_INT rhs, PASTIX_INT *iparm,
                        double *dparm),

                       (pastix_data_t **pastix_data, MPI_Comm pastix_comm,
                        PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
                        double *avals, PASTIX_INT * loc2glob, PASTIX_INT *perm, PASTIX_INT *invp,
                        double *b, PASTIX_INT rhs, PASTIX_INT *iparm,
                        double *dparm),

                       (pastix_data_t **pastix_data, MPI_Comm pastix_comm,
                        PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
                        COMPLEX *avals, PASTIX_INT * loc2glob, PASTIX_INT *perm, PASTIX_INT *invp,
                        COMPLEX *b, PASTIX_INT rhs, PASTIX_INT *iparm,
                        double *dparm),

                       (pastix_data_t **pastix_data, MPI_Comm pastix_comm,
                        PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
                        DCOMPLEX *avals, PASTIX_INT * loc2glob, PASTIX_INT *perm, PASTIX_INT *invp,
                        DCOMPLEX *b, PASTIX_INT rhs, PASTIX_INT *iparm,
                        double *dparm))
#if (defined PASTIX_FLOAT)
/*
 * Group: Thread functions
 */

/*
  Function: pastix_bindThreads

  Set bindtab in pastix_data, it gives for each thread the CPU to bind in to.
  bindtab follows this organisation :

  bindtab[threadnum] = cpu to set thread threadnum.

  Parameters:
    pastix_data - Structure de donnée pour l'utilisation step by step
    thrdnbr     - Nombre de threads / Taille du tableau
    bindtab     - Tableau de correspondance entre chaque thread et coeur de la machine
*/
void pastix_bindThreads(pastix_data_t *pastix_data, PASTIX_INT thrdnbr, PASTIX_INT *bindtab);
#endif
MULTIPLE_TYPE_DEFINE(void, pastix_bindThreads,
                     (pastix_data_t *pastix_data, PASTIX_INT thrdnbr, PASTIX_INT *bindtab))
#if (defined PASTIX_FLOAT)
/*
 * Group: Checking the matrix.
 */

/*
 * Function: pastix_checkMatrix
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
 */
PASTIX_INT pastix_checkMatrix(MPI_Comm pastix_comm, PASTIX_INT verb, PASTIX_INT flagsym, PASTIX_INT flagcor,
                       PASTIX_INT n, PASTIX_INT **colptr, PASTIX_INT **row, PASTIX_FLOAT **avals,
                       PASTIX_INT **loc2glob, PASTIX_INT dof);
#endif
MULTIPLE_TYPE_DEFINE_F(PASTIX_INT, pastix_checkMatrix,
                       (MPI_Comm pastix_comm, PASTIX_INT verb,
                        PASTIX_INT flagsym, PASTIX_INT flagcor,
                        PASTIX_INT n, PASTIX_INT **colptr, PASTIX_INT **row, float **avals,
                        PASTIX_INT **loc2glob, PASTIX_INT dof),

                       (MPI_Comm pastix_comm, PASTIX_INT verb,
                        PASTIX_INT flagsym, PASTIX_INT flagcor,
                        PASTIX_INT n, PASTIX_INT **colptr, PASTIX_INT **row, double **avals,
                        PASTIX_INT **loc2glob, PASTIX_INT dof),

                       (MPI_Comm pastix_comm, PASTIX_INT verb,
                        PASTIX_INT flagsym, PASTIX_INT flagcor,
                        PASTIX_INT n, PASTIX_INT **colptr, PASTIX_INT **row, COMPLEX **avals,
                        PASTIX_INT **loc2glob, PASTIX_INT dof),

                       (MPI_Comm pastix_comm, PASTIX_INT verb,
                        PASTIX_INT flagsym, PASTIX_INT flagcor,
                        PASTIX_INT n, PASTIX_INT **colptr, PASTIX_INT **row, DCOMPLEX **avals,
                        PASTIX_INT **loc2glob, PASTIX_INT dof))

#if (defined PASTIX_FLOAT)
/*
 * Group: Getting solver distribution.
 */

/*
  Function: pastix_getLocalNodeNbr

  Return the node number in the new distribution computed by analyze step.
  Needs analyze step to be runned with pastix_data before.

  Parameters:
    pastix_data - Data used for a step by step execution.

  Returns:
    Number of local nodes/columns in new distribution.
 */
PASTIX_INT pastix_getLocalNodeNbr(pastix_data_t ** pastix_data);
#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, pastix_getLocalNodeNbr,
                     (pastix_data_t ** pastix_data))

#if (defined PASTIX_FLOAT)
/*
  Function: pastix_getLocalNodeLst

  Fill in nodelst with the list of local nodes/columns.
  Needs nodelst to be allocated with nodenbr*sizeof(pastix_int_t),
  where nodenbr has been computed by <pastix_getLocalNodeNbr>.

  Parameters:
    pastix_data - Data used for a step by step execution.
    nodelst     - An array where to write the list of local nodes/columns.
 */
PASTIX_INT pastix_getLocalNodeLst(pastix_data_t ** pastix_data, PASTIX_INT * nodelst);
#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, pastix_getLocalNodeLst,
                     (pastix_data_t ** pastix_data, PASTIX_INT * nodelst))

#if (defined PASTIX_FLOAT)
/*
  Function: pastix_getLocalUnknownNbr

  Return the unknown number in the new distribution computed by analyze step.
  Needs analyze step to be runned with pastix_data before.

  Parameters:
    pastix_data - Data used for a step by step execution.

  Returns:
    Number of local unknowns/columns in new distribution.
 */
PASTIX_INT pastix_getLocalUnknownNbr(pastix_data_t ** pastix_data);
#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, pastix_getLocalUnknownNbr,
                     (pastix_data_t ** pastix_data))
#if (defined PASTIX_FLOAT)
/*
  Function: pastix_getLocalUnknownLst

  Fill in unknownlst with the list of local unknowns/clumns.
  Needs unknownlst to be allocated with unknownnbr*sizeof(pastix_int_t),
  where unknownnbr has been computed by <pastix_getLocalUnknownNbr>.

  Parameters:
    pastix_data - Data used for a step by step execution.
    unknownlst     - An array where to write the list of local unknowns/columns.
 */
PASTIX_INT pastix_getLocalUnknownLst(pastix_data_t ** pastix_data, PASTIX_INT * unknownlst);
#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, pastix_getLocalUnknownLst,
                     (pastix_data_t ** pastix_data, PASTIX_INT * unknownlst))

#if (defined PASTIX_FLOAT)


/*
 * Group: About the Schur complement.
 */


/*
  Function: pastix_setSchurUnknownList

  Set the list of unknowns to isolate at the end
  of the matrix via permutations.

  Has to be called if using IPARM_SCHUR = API_YES.

  Parameters:
    pastix_data - Data used for a step by step execution.
    n           - Number of unknowns.
    list        - List of unknowns.
*/
PASTIX_INT pastix_setSchurUnknownList(pastix_data_t * pastix_data,
             PASTIX_INT  n,
             PASTIX_INT *list);
#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, pastix_setSchurUnknownList,
                     (pastix_data_t * pastix_data,
                      PASTIX_INT  n,
                      PASTIX_INT *list))

#ifdef PASTIX_FLOAT
/*
  Function: pastix_getSchurLocalNodeNbr

  Compute the number of nodes in the local part of the Schur.

  Parameters:
    pastix_data - Common data structure for PaStiX calls.
    nodeNbr     - (out) Number of nodes in Schur (local).

  Returns:
    NO_ERR      - For the moment

  TODO: Error management.
*/
PASTIX_INT pastix_getSchurLocalNodeNbr(pastix_data_t * pastix_data, PASTIX_INT * nodeNbr);
#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, pastix_getSchurLocalNodeNbr,
                     (pastix_data_t * pastix_data, PASTIX_INT * nodeNbr))

#ifdef PASTIX_FLOAT
/*
  Function: pastix_getSchurLocalUnkownNbr

  Compute the number of unknowns in the local part of the Schur.

  Parameters:
    pastix_data - Common data structure for PaStiX calls.
    unknownNbr  - (out) Number of unknowns in Schur (local).

  Returns:
    NO_ERR      - For the moment

  TODO: Error management.
*/
PASTIX_INT pastix_getSchurLocalUnkownNbr(pastix_data_t * pastix_data,
                                  PASTIX_INT * unknownNbr);
#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, pastix_getSchurLocalUnkownNbr,
                     (pastix_data_t * pastix_data, PASTIX_INT * unknownNbr))

#ifdef PASTIX_FLOAT
/*
  Function: pastix_getSchurLocalNodeList

  Compute the list of nodes in the local part of the Schur.

  Parameters:
    pastix_data - Common data structure for PaStiX calls.
    nodes     - (out) Nodes in Schur (local).

  Returns:
    NO_ERR      - For the moment

  TODO: Error management.
*/
PASTIX_INT pastix_getSchurLocalNodeList(pastix_data_t * pastix_data, PASTIX_INT * nodes);
#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, pastix_getSchurLocalNodeList,
                     (pastix_data_t * pastix_data, PASTIX_INT * nodes))

#ifdef PASTIX_FLOAT
/*
  Function: pastix_getSchurLocalUnkownList

  Compute the list of unknowns in the local part of the Schur.

  Parameters:
    pastix_data - Common data structure for PaStiX calls.
    unknowns    - (out) Unknowns in Schur (local).

  Returns:
    NO_ERR      - For the moment

  TODO: Error management.
*/
PASTIX_INT pastix_getSchurLocalUnknownList(pastix_data_t * pastix_data,
                                    PASTIX_INT * unknowns);
#endif
MULTIPLE_TYPE_DEFINE(PASTIX_INT, pastix_getSchurLocalUnknownList,
                     (pastix_data_t * pastix_data, PASTIX_INT * unknowns))

#ifdef PASTIX_FLOAT
/*
  Function: pastix_setSchurArray

  Give user memory area to store Schur in PaStiX.

  Parameters:
    pastix_data - Common data structure for PaStiX calls.
    array       - Memory area to store the Schur.

  Returns:
    NO_ERR      - For the moment

  TODO: Error management.
*/
PASTIX_INT pastix_setSchurArray(pastix_data_t * pastix_data, PASTIX_FLOAT * array);
#endif
MULTIPLE_TYPE_DEFINE_F(PASTIX_INT, pastix_setSchurArray,
                       (pastix_data_t * pastix_data, float * array),
                       (pastix_data_t * pastix_data, double * array),
                       (pastix_data_t * pastix_data, COMPLEX * array),
                       (pastix_data_t * pastix_data, DCOMPLEX * array))

#if (defined PASTIX_FLOAT)
/*
  Function: pastix_getSchur

  Get the Schur complement from PaStiX.

  Schur complement is a dense block in a
  column scheme.

  Parameters:
    pastix_data - Data used for a step by step execution.
    schur - Array to fill-in with Schur complement.

*/
PASTIX_INT pastix_getSchur(pastix_data_t * pastix_data,
                    PASTIX_FLOAT * schur);

#endif
MULTIPLE_TYPE_DEFINE_F(PASTIX_INT, pastix_getSchur,
                       (pastix_data_t * pastix_data, float * schur),
                       (pastix_data_t * pastix_data, double * schur),
                       (pastix_data_t * pastix_data, COMPLEX * schur),
                       (pastix_data_t * pastix_data, DCOMPLEX* schur))

#if (defined PASTIX_FLOAT)


/*
 * Group: About parameters.
 */

/*
  Function: pastix_initParam

  Sets default parameters for iparm and dparm

  Parameters:
    iparm - tabular of IPARM_SIZE integer parameters.
    dparm - tabular of DPARM_SIZE double parameters.
*/
void pastix_initParam(PASTIX_INT    *iparm,
                      double *dparm);
#endif
MULTIPLE_TYPE_DEFINE(void, pastix_initParam, (PASTIX_INT    *iparm,
                                              double *dparm))
#if (defined PASTIX_FLOAT)
/*
 * Group: About Scaling.
 *
 * Working on it...
 */

/*
 * Function: pastix_unscale
 *
 * Unscale the factorized SolverMatrix.
 *
 * Unscale the factorized SolverMatrix in pastix_data,
 * according to pastix_data->scalerowtab,
 * pastix_data->iscalerowtab, pastix_data->scalecoltab and
 * pastix_data->iscalecoltab
 *
 * (Elements in pastix_data->iscaleXXXtab are supposed to be the inverse of
 * elements in pastix_data->scaleXXXtab).
 *
 * Parameters:
 *   pastix_data - Common data structure for PaStiX calls.
 *   sym         - API_YES if the matrix is symmetric
 *                   ( only pastix_data->scalerowtab and
`*                     pastix_data->iscalerowtab
 *                     are used in this case),
 *                 API_NO otherwise
 */
void pastix_unscale ( pastix_data_t *pastix_data, PASTIX_INT sym);
#endif /* PASTIX_FLOAT */
MULTIPLE_TYPE_DEFINE(void, pastix_unscale, ( pastix_data_t *pastix_data, PASTIX_INT sym))

unsigned long pastix_getMemoryUsage(void);
MULTIPLE_TYPE_DEFINE(unsigned long, pastix_getMemoryUsage, (void))
unsigned long pastix_getMaxMemoryUsage(void);
MULTIPLE_TYPE_DEFINE(unsigned long, pastix_getMaxMemoryUsage, (void))

/* clean up defines */
#undef MULTIPLE_TYPE_DEFINE
#undef MULTIPLE_TYPE_DEFINE_F
#undef COMPLEX
#undef DCOMPLEX
#undef PASTIX_HAS_COMPLEX
#endif /* not _PASTIX_H_ */
