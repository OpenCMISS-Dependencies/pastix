/* File: cppsimple.cpp
 *
 *  A simple example :
 *  read the matrix, check it is correct and correct it if needed,
 *  then run pastix in one call.
 *
 */

#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#ifndef FORCE_NOMPI
#include <mpi.h>
#endif
#include <complex>
#include <iostream>
/* to access functions from the libpastix, respect this order */
namespace PaStiX {
	extern "C" {
#include "pastix.h"
#include "read_matrix.h"
#include "get_options.h"
	}
}
#include "utils.hpp"
#ifdef FORCE_NOMPI
#define MPI_COMM_WORLD 0
#endif

int main (int argc, char **argv)
{
  PaStiX::pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
  PaStiX::pastix_int_t    ncol;               /* Size of the matrix                                        */
  PaStiX::pastix_int_t   *colptr      = NULL; /* Indexes of first element of each column in row and values */
  PaStiX::pastix_int_t   *rows        = NULL; /* Row of each element of the matrix                         */
  PaStiX::pastix_float_t *values      = NULL; /* Value of each element of the matrix                       */
  PaStiX::pastix_float_t *rhs         = NULL; /* right hand side                                           */
  PaStiX::pastix_float_t *rhssaved    = NULL; /* right hand side (save)                                    */
  PaStiX::pastix_float_t *ax          = NULL; /* A times X product                                         */
  PaStiX::pastix_int_t    iparm[PaStiX::IPARM_SIZE];  /* integer parameters for pastix                             */
  double                  dparm[PaStiX::DPARM_SIZE];  /* floating parameters for pastix                            */
  PaStiX::pastix_int_t   *perm        = NULL; /* Permutation tabular                                       */
  PaStiX::pastix_int_t   *invp        = NULL; /* Reverse permutation tabular                               */
  char                   *type        = NULL; /* type of the matrix                                        */
  char                   *rhstype     = NULL; /* type of the right hand side                               */
#ifndef FORCE_NOMPI
  int             required;           /* MPI thread level required                                 */
  int             provided;           /* MPI thread level provided                                 */
#endif
  int             mpid;
  PaStiX::driver_type_t  *driver_type;        /* Matrix driver(s) requested by user                        */
  char          **filename;           /* Filename(s) given by user                                 */
  int             nbmatrices;         /* Number of matrices given by user                          */
  int             nbthread;           /* Number of thread wanted by user                           */
  int             verbosemode;        /* Level of verbose mode (0, 1, 2)                           */
  int             ordering;           /* Ordering to use                                           */
  int             nbrhs = 1;
  int             incomplete;         /* Indicate if we want to use incomplete factorisation       */
  int             level_of_fill;      /* Level of fill for incomplete factorisation                */
  int             amalgamation;       /* Level of amalgamation for Kass                            */
  int             ooc;                /* OOC limit (Mo/percent depending on compilation options)   */
  PaStiX::pastix_int_t    mat_type;
  long            i;
  double norme1, norme2;

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
    /* Matrices for each iteration must have the same patern,
       this is why we only authorize one matrix in this exemple.
       But it could be used with several matrices with same patern
       and different values.
    */
    std::cerr << "WARNING: should have only one matrix" << std::endl;
  }

  /*******************************************/
  /*      Read Matrice                       */
  /*******************************************/
  PaStiX::read_matrix(filename[0], &ncol, &colptr, &rows, &values,
		      &rhs, &type, &rhstype, driver_type[0], MPI_COMM_WORLD);

  /* as it has been allocated in C with malloc, we need a free, not a delete [] */
  for (i = 0; i < nbmatrices; i++)
    if (filename[i] != NULL)
      free(filename[i]);
  free(filename);
  free(driver_type);

  mat_type = PaStiX::API_SYM_NO;
  if (MTX_ISSYM(type)) mat_type = PaStiX::API_SYM_YES;
  if (MTX_ISHER(type)) mat_type = PaStiX::API_SYM_HER;
  /*******************************************/
  /*    Check Matrix format                  */
  /*******************************************/
  /*
   * Matrix needs :
   *    - to be in fortran numbering
   *    - to have only the lower triangular part in symmetric case
   *    - to have a graph with a symmetric structure in unsymmetric case
   */
  if (PaStiX::NO_ERR != PaStiX::pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
						   mat_type,
						   PaStiX::API_YES,
						   ncol, &colptr, &rows, &values, NULL, 1))
    return 1;

  /*******************************************/
  /* Initialize parameters to default values */
  /*******************************************/
  iparm[PaStiX::IPARM_MODIFY_PARAMETER] = PaStiX::API_NO;
  PaStiX::pastix(&pastix_data, MPI_COMM_WORLD,
		 ncol, colptr, rows, values,
		 perm, invp, rhs, 1, iparm, dparm);

  /*******************************************/
  /*       Customize some parameters         */
  /*******************************************/
  iparm[PaStiX::IPARM_THREAD_NBR] = nbthread;
  iparm[PaStiX::IPARM_SYM] = mat_type;
  switch (mat_type)
  {
  case PaStiX::API_SYM_YES:
    iparm[PaStiX::IPARM_FACTORIZATION] = PaStiX::API_FACT_LDLT;
    break;
  case PaStiX::API_SYM_HER:
    iparm[PaStiX::IPARM_FACTORIZATION] = PaStiX::API_FACT_LDLH;
    break;
  default:
    iparm[PaStiX::IPARM_FACTORIZATION] = PaStiX::API_FACT_LU;
  }
  iparm[PaStiX::IPARM_MATRIX_VERIFICATION] = PaStiX::API_NO;
  iparm[PaStiX::IPARM_VERBOSE]             = verbosemode;
  iparm[PaStiX::IPARM_ORDERING]            = ordering;
  iparm[PaStiX::IPARM_INCOMPLETE]          = incomplete;
  iparm[PaStiX::IPARM_OOC_LIMIT]           = ooc;

  if (incomplete == 1)
    {
      dparm[PaStiX::DPARM_EPSILON_REFINEMENT] = 1e-7;
    }
  iparm[PaStiX::IPARM_LEVEL_OF_FILL]       = level_of_fill;
  iparm[PaStiX::IPARM_AMALGAMATION_LEVEL]  = amalgamation;
  iparm[PaStiX::IPARM_RHS_MAKING]          = PaStiX::API_RHS_B;
  /* read parameters to set IPARM/DPARM */
  if (EXIT_FAILURE == PaStiX::get_idparm(argc,   argv,
					 iparm,  dparm))
    return EXIT_FAILURE;
  iparm[PaStiX::IPARM_START_TASK]          = PaStiX::API_TASK_ORDERING;
  iparm[PaStiX::IPARM_END_TASK]            = PaStiX::API_TASK_CLEAN;

  /*******************************************/
  /*           Save the rhs                  */
  /*    (it will be replaced by solution)    */
  /*******************************************/
  rhssaved = new PaStiX::pastix_float_t[nbrhs*ncol];
  for(i=0; i<nbrhs; i++)
    memcpy(rhssaved+(i*ncol), rhs, ncol*sizeof(PaStiX::pastix_float_t));

  /*******************************************/
  /*           Call pastix                   */
  /*******************************************/
  perm = new PaStiX::pastix_int_t [ncol];
  invp = new PaStiX::pastix_int_t [ncol];

  PRINT_RHS("RHS", rhs, ncol, mpid, iparm[PaStiX::IPARM_VERBOSE]);

  PaStiX::pastix(&pastix_data, MPI_COMM_WORLD,
		 ncol, colptr, rows, values,
		 perm, invp, rhs, nbrhs, iparm, dparm);

  PRINT_RHS("SOL", rhs, ncol, mpid, iparm[PaStiX::IPARM_VERBOSE]);
  CHECK_SOL(rhs, rhssaved, ncol, mpid);

  /* as it has been allocated in C with malloc, we need a free, not a delete [] */
  free(colptr);
  free(rows);
  free(values);
  free(rhs);
  free(type);
  free(rhstype);
  delete [] perm;
  delete [] invp;
  delete [] rhssaved;
#ifndef FORCE_NOMPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
