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
/* File: mult-rhs-dist.c

   An example with a distributed matrix and multiple right-hand-side:
   read the matrix, check it is correct and correct it if needed,
   distribute it then run pastix in one call.

*/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifndef FORCE_NOMPI
#include <mpi.h>
#else
#define MPI_COMM_WORLD 0
#endif
#include <complex.h>
/* to access functions from the libpastix, respect this order */
#include "pastix.h"
#include "cscd_utils.h"
#include "read_matrix.h"
#include "get_options.h"
#include "utils.h"

int main (int argc, char **argv)
{

  pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
  pastix_int_t    ncol;               /* Number of local columns                                   */
  pastix_int_t   *colptr      = NULL; /* Indexes of first element of each column in row and values */
  pastix_int_t   *rows        = NULL; /* Row of each element of the matrix                         */
  pastix_int_t   *loc2glob    = NULL; /* Local to local column correspondance                      */
  pastix_float_t *values      = NULL; /* Value of each element of the matrix                       */
  pastix_int_t    ncol2;              /* Size of the matrix                                        */
  pastix_int_t   *colptr2     = NULL; /* Indexes of first element of each column in row and values */
  pastix_int_t   *rows2       = NULL; /* Row of each element of the matrix                         */
  pastix_int_t   *loc2glob2   = NULL; /* Local to local column correspondance                      */
  pastix_float_t *values2     = NULL; /* Value of each element of the matrix                       */
  pastix_float_t *rhs         = NULL; /* right-hand-side                                           */
  pastix_float_t *rhs2        = NULL; /* right-hand-side                                           */
  pastix_float_t *rhssaved    = NULL; /* right hand side (save)                                    */
  pastix_float_t *rhssaved_g  = NULL; /* right hand side (save, global)                            */
  pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                             */
  double          dparm[DPARM_SIZE];  /* floating parameters for pastix                            */
  pastix_int_t   *perm        = NULL; /* Permutation tabular                                       */
  pastix_int_t   *invp        = NULL; /* Reverse permutation tabular                               */
  char           *type        = NULL; /* type of the matrix                                        */
  char           *rhstype     = NULL; /* type of the right hand side                               */
#ifndef FORCE_NOMPI
  int             required;           /* MPI thread level required                                 */
  int             provided;           /* MPI thread level provided                                 */
#endif
  int             mpid;
  driver_type_t  *driver_type;        /* Matrix driver(s) requested by user                        */
  char          **filename;           /* Filename(s) given by user                                 */
  int             nbmatrices;         /* Number of matrices given by user                          */
  int             nbthread;           /* Number of thread wanted by user                           */
  int             verbosemode;        /* Level of verbose mode (0, 1, 2)                           */
  int             ordering;           /* Ordering to use                                           */
  int             incomplete;         /* Indicate if we want to use incomplete factorisation       */
  int             level_of_fill;      /* Level of fill for incomplete factorisation                */
  int             amalgamation;       /* Level of amalgamation for Kass                            */
  int             ooc;                /* OOC limit (Mo/percent depending on compilation options)   */
  pastix_int_t    mat_type;
  long            i;
  pastix_int_t    globn;
  pastix_int_t    nrhs = 2;
  /*******************************************/
  /*          MPI initialisation             */
  /*******************************************/
#ifndef FORCE_NOMPI
  required = MPI_THREAD_MULTIPLE;
  provided = -1;
  MPI_Init_thread(&argc, &argv, required, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
  if (mpid == 0)
    {
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
    }
#else
  mpid = 0;
#endif

  /*******************************************/
  /*    Get options from command line        */
  /*******************************************/
  if (EXIT_FAILURE == get_options(argc, argv,     &driver_type,
                                  &filename,      &nbmatrices,
                                  &nbthread,      &verbosemode,
                                  &ordering,      &incomplete,
                                  &level_of_fill, &amalgamation,
                                  &ooc,           &ncol))
    return EXIT_FAILURE;

  if (nbmatrices != 1)
    {
      /* Matrices for each iteration must have the same patern, this is why we only
         authorize one matrix in this exemple.
         But it could be used with several matrices with same patern and different values.
      */
      fprintf(stderr,"WARNING: should have only one matrix\n");
    }
  /*******************************************/
  /*      Read Matrice                       */
  /*******************************************/
  dread_matrix(filename[0], &ncol, &colptr, &rows, &loc2glob, &values,
               &rhs, &type, &rhstype, driver_type[0], MPI_COMM_WORLD);
  rhs = realloc(rhs, ncol*nrhs*sizeof(pastix_float_t));
  for (i = 1; i < nrhs; i++)
    memcpy(&(rhs[i*ncol]), rhs, ncol*sizeof(pastix_float_t));
  for (i = 0; i < nrhs; i++)
    {
      PRINT_RHS("RHS1", (&(rhs[i*ncol])), ncol, mpid, verbosemode);
    }

  mat_type = API_SYM_NO;
  if (MTX_ISSYM(type)) mat_type = API_SYM_YES;
  if (MTX_ISHER(type)) mat_type = API_SYM_HER;

  for (i = 0; i < nbmatrices; i++)
    if (filename[i] != NULL)
      free(filename[i]);
  free(filename);
  free(driver_type);

  /*******************************************/
  /*    Check Matrix format                  */
  /*******************************************/
  /*
   * Matrix needs :
   *    - to be in fortran numbering
   *    - to have only the lower triangular part in symmetric case
   *    - to have a graph with a symmetric structure in unsymmetric case
   */
  pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
                     mat_type,  API_YES,
                     ncol, &colptr, &rows, &values, &loc2glob, 1);

  /*******************************************/
  /* Initialize parameters to default values */
  /*******************************************/
  iparm[IPARM_MODIFY_PARAMETER] = API_NO;
  dpastix(&pastix_data, MPI_COMM_WORLD,
          0, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, 0, iparm, dparm);

  /*******************************************/
  /*       Customize some parameters         */
  /*******************************************/
  iparm[IPARM_THREAD_NBR] = nbthread;
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
  iparm[IPARM_RHS_MAKING]          = API_RHS_B;
  /* reread parameters to set IPARM/DPARM */
  if (EXIT_FAILURE == get_idparm(argc, argv,
                                 iparm,          dparm))
    return EXIT_FAILURE;

  iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
  iparm[IPARM_END_TASK]            = API_TASK_BLEND;


  /*******************************************/
  /*           Save the rhs                  */
  /*    (it will be replaced by solution)    */
  /*******************************************/
  rhssaved = malloc(ncol*sizeof(pastix_float_t));
  memcpy(rhssaved, rhs, ncol*sizeof(pastix_float_t));
#ifndef FORCE_NOMPI
  MPI_Allreduce(&ncol, &globn, 1, MPI_PASTIX_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  globn = ncol;
#endif
  rhssaved_g = malloc(globn*sizeof(pastix_float_t));
  memset(rhssaved_g, 0, globn*sizeof(pastix_float_t));
  for (i = 0; i < ncol; i++)
    {
      rhssaved_g[loc2glob[i]-1] = rhssaved[i];
    }

  free(rhssaved);
#ifndef FORCE_NOMPI
  {
    pastix_float_t * rhssaved_g_rcv;
    rhssaved_g_rcv = malloc(globn*sizeof(pastix_float_t));
    MPI_Allreduce(rhssaved_g, rhssaved_g_rcv, globn, MPI_PASTIX_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    free(rhssaved_g);
    rhssaved_g = rhssaved_g_rcv;
  }
#endif
  /*******************************************/
  /*           Call pastix                   */
  /*******************************************/
  perm = malloc(ncol*sizeof(pastix_int_t));
  /* No need to allocate invp in dpastix */

  dpastix(&pastix_data, MPI_COMM_WORLD,
          ncol, colptr, rows, NULL, loc2glob,
          perm, NULL, NULL, nrhs, iparm, dparm);

  /* Redistributing the matrix */

  ncol2 = pastix_getLocalNodeNbr(&pastix_data);

  if (ncol2 != 0)
    {
      if (NULL == (loc2glob2 = malloc(ncol2 * sizeof(pastix_int_t))))
        {
          fprintf(stderr, "%s:%d Malloc error\n", __FILE__, __LINE__);
          return EXIT_FAILURE;
        }
    }
  else
    {
      loc2glob2 = NULL;
    }

  pastix_getLocalNodeLst(&pastix_data, loc2glob2);

  if (EXIT_SUCCESS != cscd_redispatch(ncol,   colptr,   rows,  values,    rhs,  nrhs, loc2glob,
                                      ncol2, &colptr2, &rows2, &values2, &rhs2, loc2glob2,
                                      MPI_COMM_WORLD, iparm[IPARM_DOF_NBR]))
    return EXIT_FAILURE;


  free(colptr);
  free(rows);
  free(values);
  free(rhs);
  free(loc2glob);
  free(perm);

  iparm[IPARM_START_TASK]          = API_TASK_NUMFACT;
  iparm[IPARM_END_TASK]            = API_TASK_SOLVE;
  for (i = 0; i < nrhs; i++)
    {
      PRINT_RHS("RHS", (&(rhs2[i*ncol2])), ncol2, mpid, iparm[IPARM_VERBOSE]);
    }
  dpastix(&pastix_data, MPI_COMM_WORLD,
          ncol2, colptr2, rows2, values2, loc2glob2,
          perm, invp, rhs2, nrhs, iparm, dparm);
  for (i = 0; i < nrhs; i++)
    {
      PRINT_RHS("SOL", (&(rhs2[i*ncol2])), ncol2, mpid, iparm[IPARM_VERBOSE]);
      CHECK_DIST_SOL(colptr2, rows2, values2, (&(rhs2[i*ncol2])), ncol2, loc2glob2, globn, rhssaved_g);
    }
  iparm[IPARM_START_TASK]          = API_TASK_CLEAN;
  iparm[IPARM_END_TASK]            = API_TASK_CLEAN;
  dpastix(&pastix_data, MPI_COMM_WORLD,
          ncol2, colptr2, rows2, values2, loc2glob2,
          perm, invp, rhs2, nrhs, iparm, dparm);

  free(colptr2);
  free(rows2);
  free(values2);
  free(rhs2);
  free(rhssaved_g);
  free(type);
  free(rhstype);
  free(loc2glob2);
#ifndef FORCE_NOMPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
