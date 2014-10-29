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
/* File: reentrant.c

   A simple example :
   run two threads then run two instances of PaStiX in each.

*/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#ifndef FORCE_NOMPI
#include <mpi.h>
#endif
#include <complex.h>
/* to access functions from the libpastix, respect this order */
#include "pastix.h"
#include "read_matrix.h"
#include "get_options.h"

#ifdef FORCE_NOMPI
#define MPI_COMM_WORLD 0
#endif
/*
  Struct: solv_param

  Structure containing information to give
  to each thread.

  Contains:
    ncol        - number of columns
    colptr      - Indexes of first element of each column in row and values
    rows        - Row of each element of the matrix
    values      - Value of each element of the matrix
    rhs         - right-hand-side member
    type        - type of the matrix
    rhstype     - type of the right-hand-side member
    nbthread    - number of threads
    verbosemode - Verbose mode
    ordering    - ordering library.
 */
typedef struct solve_param {
  pastix_int_t    ncol;
  pastix_int_t   *colptr;
  pastix_int_t   *rows;
  pastix_float_t *values;
  pastix_float_t *rhs;
  char           *type;
  char           *rhstype;
  int             nbthread;
  int             verbosemode;
  int             ordering;
  int             incomplete;
  int             level_of_fill;
  int             amalgamation;
  int             ooc;
  int             argc;
  char          **argv;
} solve_param_t;

/*
  Function: solve_smp

  Thread routine to launch the solver

  Parameters:
    arg - a pointer to a <solve_param> structure.
*/
static void *solve_smp(void *arg)
{

  pastix_data_t    *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
  pastix_int_t      ncol;               /* Size of the matrix                                        */
  pastix_int_t     *colptr      = NULL; /* Indexes of first element of each column in row and values */
  pastix_int_t     *rows        = NULL; /* Row of each element of the matrix                         */
  pastix_float_t   *values      = NULL; /* Value of each element of the matrix                       */
  pastix_float_t   *rhs         = NULL; /* right hand side                                           */
  pastix_int_t      iparm[IPARM_SIZE];  /* integer parameters for pastix                             */
  double            dparm[DPARM_SIZE];  /* floating parameters for pastix                            */
  pastix_int_t     *perm        = NULL; /* Permutation tabular                                       */
  pastix_int_t     *invp        = NULL; /* Reverse permutation tabular                               */
  char             *type        = NULL; /* type of the matrix                                        */
  char             *rhstype     = NULL; /* type of the right-hand-side member                        */
  int               nbthread;           /* Number of thread wanted by user                           */
  int               verbosemode;        /* Level of verbose mode (0, 1, 2)                           */
  int               ordering;           /* ordering option value                                     */
  int               incomplete;         /* Indicate if we want to use incomplete factorisation       */
  int               level_of_fill;      /* Level of fill for incomplete factorisation                */
  int               amalgamation;       /* Level of amalgamation for Kass                             */
  int               ooc;                /* OOC limit (Mo/percent depending on compilation options)   */
  pastix_int_t      mat_type;

  solve_param_t param = *(solve_param_t *)arg;
  nbthread            = param.nbthread;
  colptr              = param.colptr;
  rows                = param.rows;
  values              = param.values;
  ncol                = param.ncol;
  type                = param.type;
  rhstype             = param.rhstype;
  rhs                 = param.rhs;
  verbosemode         = param.verbosemode;
  ordering            = param.ordering;
  incomplete          = param.incomplete;
  level_of_fill       = param.level_of_fill;
  amalgamation        = param.amalgamation;
  ooc                 = param.ooc;
  /*******************************************/
  /* Initialize parameters to default values */
  /*******************************************/
  iparm[IPARM_MODIFY_PARAMETER] = API_NO;
  pastix(&pastix_data, MPI_COMM_WORLD,
         ncol, colptr, rows, values,
         perm, invp, rhs, 1, iparm, dparm);

  /*******************************************/
  /*       Customize some parameters         */
  /*******************************************/
  iparm[IPARM_THREAD_NBR] = nbthread;

  mat_type = API_SYM_NO;
  if (MTX_ISSYM(type)) mat_type = API_SYM_YES;
  if (MTX_ISHER(type)) mat_type = API_SYM_HER;
  iparm[IPARM_SYM] = mat_type;
  switch (mat_type)
  {
  case API_SYM_YES:
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    break;
  case API_SYM_HER:
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLH;
    break;
  default:
    iparm[IPARM_FACTORIZATION] = API_FACT_LU;
  }
  iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
  iparm[IPARM_VERBOSE]             = verbosemode;
  iparm[IPARM_ORDERING]            = ordering;
  iparm[IPARM_INCOMPLETE]          = incomplete;
  iparm[IPARM_OOC_LIMIT]           = ooc;

  if (incomplete == 1)
    {
      dparm[DPARM_EPSILON_REFINEMENT] = 1e-7;
    }
  iparm[IPARM_LEVEL_OF_FILL]       = level_of_fill;
  iparm[IPARM_AMALGAMATION_LEVEL]  = amalgamation;
  if (EXIT_FAILURE == get_idparm(param.argc, param.argv,
                                 iparm,          dparm))
    return NULL;
  iparm[IPARM_RHS_MAKING]          = API_RHS_1;

  /*******************************************/
  /*           Call pastix                   */
  /*******************************************/
  perm = malloc(ncol*sizeof(pastix_int_t));
  invp = malloc(ncol*sizeof(pastix_int_t));

  pastix(&pastix_data, MPI_COMM_WORLD,
         ncol, colptr, rows, values,
         perm, invp, rhs, 1, iparm, dparm);

  free(colptr);
  free(rows);
  free(values);
  free(rhs);
  free(perm);
  free(invp);
  free(type);
  free(rhstype);

  return NULL;
}

/*
  Function: main

  see <reentrant.c>

 */
int main (int argc, char **argv)
{

#ifndef FORCE_NOMPI
  int               required;           /* MPI thread level required                                 */
  int               provided;           /* MPI thread level provided                                 */
  int               commSize;           /* MPI rank                                                  */
#endif
  driver_type_t    *driver_type;        /* Matrix driver(s) requested by user                        */
  char            **filename;           /* Filename(s) given by user                                 */
  int               nbmatrices;         /* Number of matrices given by user                          */
  int               nbthread;           /* Number of thread wanted by user                           */
  int               verbosemode;        /* Level of verbose mode (0, 1, 2)                           */
  int               ordering;           /* Ordering to use                                           */
  int               incomplete;         /* Indicate if we want to use incomplete factorisation       */
  int               level_of_fill;      /* Level of fill for incomplete factorisation                */
  int               amalgamation;       /* Level of amalgamation for Kass                             */
  pastix_int_t      ncolgen;            /* Ncol if matrix is generated                               */
  int               ooc;                /* OOC limit (Mo/percent depending on compilation options)   */
  pastix_int_t      mat_type;
  long              i;
  int               nbcallingthreads = 2;
  solve_param_t    *solve_param;
  pthread_t        *threads;


  /*******************************************/
  /*          MPI initialisation             */
  /*******************************************/
#ifndef FORCE_NOMPI
  required = MPI_THREAD_MULTIPLE;
  provided = -1;
  MPI_Init_thread(&argc, &argv, required, &provided);
  switch (provided)
    {
    case MPI_THREAD_SINGLE:
      printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
      break;
    case MPI_THREAD_FUNNELED:
      printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
      break;
    case MPI_THREAD_SERIALIZED:
      printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
      break;
    case MPI_THREAD_MULTIPLE:
      printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
      break;
    default:
      printf("MPI_Init_thread level = ???\n");
    }

  /*******************************************/
  /* Quit if communicator size > 1           */
  /*******************************************/
  /* Should work with duplicated MPI communicator.   */
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);
  if (commSize > 1)
    {
      fprintf(stderr, "ERROR : This example does not works with more than one MPI process\n");
      MPI_Finalize();
      return EXIT_FAILURE;
    }
#endif
  /*******************************************/
  /*    Get options from command line        */
  /*******************************************/
  if (EXIT_FAILURE == get_options(argc, argv,     &driver_type,
                                  &filename,      &nbmatrices,
                                  &nbthread,      &verbosemode,
                                  &ordering,      &incomplete,
                                  &level_of_fill, &amalgamation,
                                  &ooc,           &ncolgen))
    return EXIT_FAILURE;

  fprintf(stdout, "verbosemode %d\n",verbosemode);
  /*******************************************/
  /*      Check the number of matrices       */
  /*******************************************/
  if (nbmatrices != nbcallingthreads)
    {
      fprintf(stderr,"WARNING: should have %d matrix\n", nbcallingthreads);
      while ( nbmatrices < nbcallingthreads)
        {
          fprintf(stderr,"\tSetting matrice %d to the same as matrix 1\n", nbmatrices+1);
          driver_type[nbmatrices] = driver_type[0];
          if (filename[0] != NULL)
            {
              filename[nbmatrices]    = malloc((strlen(filename[0])+1)*sizeof(char));
              strcpy(filename[nbmatrices], filename[0]);
            }
          else
            {
              filename[nbmatrices] = NULL;
            }
          nbmatrices++;
        }
    }

  /*******************************************/
  /*    Set parameters for each thread       */
  /*******************************************/
  solve_param = (solve_param_t*) malloc(nbcallingthreads*sizeof(solve_param_t));
  threads     = (pthread_t*)     malloc(nbcallingthreads*sizeof(pthread_t));

  for (i = 0; i < nbcallingthreads; i++)
    {
      /*******************************************/
      /*         Read ieme Matrix                */
      /*******************************************/
      if (driver_type[i] == LAPLACIAN)
        solve_param[i].ncol = ncolgen;

      read_matrix(filename[i],
                  &(solve_param[i].ncol),
                  &(solve_param[i].colptr),
                  &(solve_param[i].rows),
                  &(solve_param[i].values),
                  &(solve_param[i].rhs),
                  &(solve_param[i].type),
                  &(solve_param[i].rhstype),
                  driver_type[i],
                  MPI_COMM_WORLD);

      solve_param[i].nbthread      = nbthread;
      solve_param[i].verbosemode   = verbosemode;
      solve_param[i].ordering      = ordering;
      solve_param[i].incomplete      = incomplete;
      solve_param[i].level_of_fill = level_of_fill;
      solve_param[i].amalgamation  = amalgamation;
      solve_param[i].ooc           = ooc;
      solve_param[i].argc          = argc;
      solve_param[i].argv          = argv;
      if (filename[i] != NULL)
        free(filename[i]);

      mat_type = API_SYM_NO;
      if (MTX_ISSYM(solve_param[i].type)) mat_type = API_SYM_YES;
      if (MTX_ISHER(solve_param[i].type)) mat_type = API_SYM_HER;
      /*******************************************/
      /*    Check Matrix format                  */
      /*******************************************/
      /*
       * Matrix needs :
       *    - to be in fortran numbering
       *    - to have only the lower triangular part in symmetric case
       *    - to have a graph with a symmetric structure in unsymmetric case
       */
      pastix_checkMatrix(MPI_COMM_WORLD,
                         API_VERBOSE_NO,
                         mat_type,
                         API_YES,
                         solve_param[i].ncol,
                         &(solve_param[i].colptr),
                         &(solve_param[i].rows),
                         &(solve_param[i].values),
                         NULL, 1);

      /*******************************************/
      /*    Launch instance of solver            */
      /*******************************************/
      pthread_create(&threads[i], NULL, solve_smp, (void *)&solve_param[i]);
   }

  /*******************************************/
  /*     Wait for the end of thread          */
  /*******************************************/
  for (i = 0; i < nbcallingthreads; i++)
    pthread_join(threads[i],(void**)NULL);

  free(filename);
  free(driver_type);
  free(threads);
  free(solve_param);

#ifndef FORCE_NOMPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
