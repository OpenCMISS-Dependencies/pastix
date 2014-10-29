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
  File: murge-product.c

  Example which build a distributed laplacian and compute a product
  using *Murge* interface.

  >  2  1  0  0 ... 0   1   4
  >  1  2  1  0 ... 0   2   8
  >   ....            x . = .
  >  0 ...    2  1  0   .   .
  >  0 ...    1  2  1   .   .
  >  0 ...    0  1  2   n   .

  Authors:
     X. LACOSTE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>
#include <math.h>
#ifdef TYPE_COMPLEX
#include <complex.h>
#endif /* TYPE_COMPLEX */
#include "murge.h"
#include "pastix.h"
#include <time.h>
#define BILLION  1000000000L;
#ifdef __APPLE__
#ifndef X_ARCHi686_mac
#define X_ARCHi686_mac
#endif
#endif
/* Section: Macros */

/*
   Macro: MIN

   Compute the minimum of two values.

   Parameters:
     x - The first value.
     y - The second value.
*/
#define MIN(x,y) (((x)<(y))?(x):(y))

#ifdef PREC_DOUBLE
#  define FLOAT_ARE_DIFFERENT(x, y) (sqrt(((x)-(y))*((x)-(y)))>1e-15)
#else
#  define FLOAT_ARE_DIFFERENT(x, y) (sqrtf(((x)-(y))*((x)-(y)))>1e-7)
#endif

/*
  Macro: MURGE_CALL

  Execute *call* fonction and test the return value.

  Parameters:
    call - The murge function to execute.
*/
#define MURGE_CALL(call)                                                \
  {                                                                     \
    INTS ret;                                                           \
    ret = call;                                                         \
    if (ret != MURGE_SUCCESS) {                                         \
      fprintf(stderr, "%s:%d error in murge call : %d\n",               \
              __FILE__, __LINE__, (int)ret);                            \
      exit(EXIT_FAILURE);                                               \
    }                                                                   \
  }

/* Section: Defines */

/* Define: DIAG_VALUE
   The value of the diagonal elements of the matrix
*/
#define DIAG_VALUE         2

/* Define: LOWER_DIAG_VALUE
   The value of the lower diagonal elements of the matrix.
*/
#define LOWER_DIAG_VALUE   1
/* Define: UPPER_DIAG_VALUE
   The value of the upper diagonal elements of the matrix.
*/
#define UPPER_DIAG_VALUE   1

/* Section: Functions */

/*
  Function: main

  Main function of <murge-product> example, solving laplacian.

  Initialize MURGE and default solver options,

  Set the distribution using <MURGE_ProductSetLocalNodeNbr>
  and <MURGE_ProductSetetLocalNodeList>.

  Fill the ditributed matrix, value by value using <MURGE_AssemblyBegin>,
  <MURGE_AssemblySetValue> and <MURGE_AssemblyEnd>.

  Set the right hand side member using <MURGE_SetLocalRHS>.

  Solve the problem and get the local solution using <MURGE_GetLocalSolution>.

*/
int main(int argc, char *argv[]) {
  COEF   *vector   = NULL;    /* Local vector to multiply                    */
  COEF   *ax       = NULL;    /* Local product                               */
  INTS   *loc2glob = NULL;    /* Local to global column numbers              */
  INTS    globalN;            /* Size  of the problem                        */
  INTL    lnnz;               /* Local number of non zeros                   */
  INTS    localn;             /* Local number of unkowns                     */
  INTL    start;              /* First local column                          */
  INTL    end;                /* Last local column                           */
  INTL    i,j;                /* iterator                                    */
  int     proc_id;            /* MPI process ID                              */
  int     nproc;              /* Number of MPI processes                     */
  INTS    id;                 /* Solver instance id                          */
  int     baseval  = 1;       /* base 0 or 1, used to define the matrix      */
  INTS    dof;                /* base 0 or 1, used to define the matrix      */
  INTS    thread_nbr = -1;
  INTS    sym;                /* indicate if the matrix pattern is symmetric */
  int     required;           /* MPI thread level required                   */
  int     provided;           /* MPI thread level provided                   */
  int     nb_prod = 10;
  COEF *  mat_elem_low = NULL, *mat_elem_up, * mat_elem_diag = NULL;
#ifndef X_ARCHi686_mac
  struct timespec t1, t2,tt1, tt2;
  double t, t_rcv;
#endif

  sym = MURGE_BOOLEAN_FALSE;

  /** Init MPI environment **/
  required = MPI_THREAD_MULTIPLE;
  provided = -1;
  MPI_Init_thread(&argc, &argv, required, &provided);
#ifndef X_ARCHi686_mac
  clock_gettime(CLOCK_REALTIME, &tt1);
#endif
  if (provided != required) {
    switch (provided) {
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
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  globalN = dof = 0;
  if (argc >= 3) {
    globalN = atoi(argv[1]);
    dof = atoi(argv[2]);
    thread_nbr = atoi(argv[3]);
  } else {
    if (proc_id == 0) {
      fprintf(stderr, "Usage: %s <size> <DofNumber> <thread_nbr>", argv[0]);
      return 1;
    }
  }

  mat_elem_low = malloc(dof*dof*sizeof(COEF));
  for (j = 0; j < dof*dof; j++)
    mat_elem_low[j] = LOWER_DIAG_VALUE;
  mat_elem_up = malloc(dof*dof*sizeof(COEF));
  for (j = 0; j < dof*dof; j++)
    mat_elem_up[j] = UPPER_DIAG_VALUE;
  mat_elem_diag = malloc(dof*dof*sizeof(COEF));
  for (j = 0; j < dof*dof; j++)
    mat_elem_diag[j] = DIAG_VALUE;

  /***************************************/
  /* Initialize MURGE for 1 problem      */
  /***************************************/
  MURGE_CALL(MURGE_Initialize(1));

  id = 0; /** id of the linear system **/

  /***************************************/
  /* Initialize Default solver options   */
  /***************************************/
  MURGE_CALL(MURGE_SetDefaultOptions(id, 0));

  /***************************************/
  /* Set the communicator                */
  /***************************************/
  MURGE_CALL(MURGE_SetCommunicator(id, MPI_COMM_WORLD));

  /**********************************/
  /* Read the matrix from file      */
  /**********************************/
  start  = proc_id*((globalN + 1)/nproc);
  end    = MIN(((proc_id+1)*((globalN + 1)/nproc))-1, globalN-1);

  if (proc_id == nproc-1)
    end = globalN -1;
  localn = end - start + 1;

  MURGE_CALL(MURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL, (INTS)baseval));
  MURGE_CALL(MURGE_SetOptionINT(id, MURGE_IPARAM_DOF, (INTS)dof));
  MURGE_CALL(MURGE_SetOptionINT(id, MURGE_IPARAM_SYM, sym));
  MURGE_CALL(MURGE_SetOptionINT(id, IPARM_VERBOSE, 2));
  MURGE_CALL(MURGE_SetOptionINT(id, IPARM_THREAD_NBR, thread_nbr));

  if (sym == MURGE_BOOLEAN_TRUE)
    lnnz = 2*localn;
  else
    lnnz = 3*localn;

  if (sym == MURGE_BOOLEAN_FALSE && start == 0)
    lnnz --;
  if (end == globalN-1)
    lnnz --;

  /***************************************************/
  /* ENTER THE GRAPH : PARALLEL INTERFACE            */
  /* Every processor deals with a part of the graph  */
  /***************************************************/
  loc2glob = malloc(localn*sizeof(INTS));
  for (i = 0; i < localn; i++) {
      loc2glob[i] = start+baseval+i;
    }

  /***************************************************/
  /*                                                 */
  /*           Set LOCAL Node LIST           */
  /***************************************************/

  MURGE_CALL(MURGE_ProductSetLocalNodeNbr (id, localn));
  MURGE_CALL(MURGE_ProductSetLocalNodeList(id, loc2glob));

  if (sym == 1)
    lnnz = 2*localn;
  else
    lnnz = 3*localn;
  for (i = 0; i < localn; i++)
    if (( sym == 0 && loc2glob[i] == baseval) ||
	loc2glob[i] == globalN - 1 + baseval)
      lnnz --;

  /***************************************************/
  /*                                                 */
  /* ENTER THE MATRIX COEFS : Distributed INTERFACE  */
  /***************************************************/
#ifndef X_ARCHi686_mac
  clock_gettime(CLOCK_REALTIME, &t1);
#endif
  MURGE_CALL(MURGE_AssemblyBegin(id, globalN, lnnz*dof*dof,
                                 MURGE_ASSEMBLY_ADD,
                                 MURGE_ASSEMBLY_ADD,
                                 MURGE_ASSEMBLY_RESPECT, sym));

  for (i = 0; i < localn; i++) {
    if ((sym == MURGE_BOOLEAN_FALSE) && (loc2glob[i] != baseval)) {
          MURGE_CALL(MURGE_AssemblySetNodeValues(id, loc2glob[i]-1,
                                                 loc2glob[i],
                                                 mat_elem_up ));
        }
      MURGE_CALL(MURGE_AssemblySetNodeValues(id, loc2glob[i],
                                             loc2glob[i], mat_elem_diag));

    if (loc2glob[i] != globalN + baseval - 1) {
          MURGE_CALL(MURGE_AssemblySetNodeValues(id, loc2glob[i]+1,
                                                 loc2glob[i],
                                                 mat_elem_low));
        }
    }
  MURGE_CALL(MURGE_AssemblyEnd(id));
#ifndef X_ARCHi686_mac
  clock_gettime(CLOCK_REALTIME, &t2);
  t = ( t2.tv_sec - t1.tv_sec )
    + (double)( t2.tv_nsec - t1.tv_nsec )
    / (double)BILLION;
  MPI_Reduce(&t,&t_rcv, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (proc_id == 0) {
      fprintf(stdout, "Assembly performed in %.3g s\n", t_rcv);
    }
#endif

  vector = (COEF*) malloc(localn*dof*sizeof(COEF));
  ax = (COEF *)malloc(sizeof(COEF)*globalN*dof);
  for (i=0;i<localn*dof;i++) {
      vector[i] = loc2glob[(i-i%dof)/dof];
    }
#ifndef X_ARCHi686_mac
  clock_gettime(CLOCK_REALTIME, &t1);
#endif
  for (i = 0; i < nb_prod; i++) {
    /****************************************************/
    /* Set the local rhs                               */
    /****************************************************/
    MURGE_CALL(MURGE_RHSReset(id));
    MURGE_CALL(MURGE_SetLocalRHS(id, vector, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW));

    /****************************************************/
    /* Get the global solution on processor 0           */
    /* Original ordering                                */
    /****************************************************/
    MURGE_CALL(MURGE_GetGlobalProduct(id, ax, 0));
  }
#ifndef X_ARCHi686_mac
  clock_gettime(CLOCK_REALTIME, &t2);
  t = ( t2.tv_sec - t1.tv_sec )
    + (double)( t2.tv_nsec - t1.tv_nsec )
    / (double)BILLION;
  MPI_Reduce(&t,&t_rcv, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  clock_gettime(CLOCK_REALTIME, &t1);
#endif
  if (proc_id == 0) {
      int OK = 1;
#ifndef X_ARCHi686_mac
      fprintf(stdout, "%d products performed in %.3g s\n", nb_prod, t_rcv);
#endif
    for (i=0;i<globalN;i++) {
      for (j = 0; j < dof; j++) {
              /*fprintf(stdout, "ax[%ld] = %lg\n", (long)i*dof+j, (double)ax[i*dof+j]);*/
        if (i == 0 || i == globalN-1) {
          if (i == 0) {
            if (FLOAT_ARE_DIFFERENT(ax[i*dof+j],
                                    dof * ((i+1)*DIAG_VALUE +
                                           (i+2)*UPPER_DIAG_VALUE))) {
              fprintf(stdout, "ax[%ld]  = %lg != %lg\n", (long)i*dof+j,
                      (double)(ax[i*dof+j]),
                      (double)(dof*((i+1)*DIAG_VALUE +
                                    (i+2)*UPPER_DIAG_VALUE)));
              OK = 0;
            }
          } else {
            if (FLOAT_ARE_DIFFERENT(ax[i*dof+j],
                                    dof * ((i+1)*DIAG_VALUE +
                                           (i)*LOWER_DIAG_VALUE))) {
              fprintf(stdout, "ax[%ld]  = %lg != %lg\n", (long)i*dof+j,
                      (double)(ax[i*dof+j]),
                      (double)(dof*((i+1)*DIAG_VALUE + (i)*LOWER_DIAG_VALUE)));
              OK = 0;
            }
          }
        } else {
          if (FLOAT_ARE_DIFFERENT(ax[i*dof+j],
                                  dof *((i)*LOWER_DIAG_VALUE + (i+1)*DIAG_VALUE +
                                        (i+2)*UPPER_DIAG_VALUE))) {
            fprintf(stdout, "ax[%ld] = %lg != %lg\n", (long)i*dof+j,
                    (double)(ax[i*dof+j]),
                    (double)(dof*((i)*LOWER_DIAG_VALUE +
                                  (i+1)*DIAG_VALUE + (i+2)* UPPER_DIAG_VALUE)));
            OK = 0;
          }
        }
      }
    }
    if (OK == 1) {
      fprintf(stdout, "Product computed correctly\n");
    }
  }
#ifndef X_ARCHi686_mac
  clock_gettime(CLOCK_REALTIME, &t2);
  t = ( t2.tv_sec - t1.tv_sec )
    + (double)( t2.tv_nsec - t1.tv_nsec )
    / (double)BILLION;
  MPI_Reduce(&t,&t_rcv, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (proc_id == 0) {
      fprintf(stdout, "Verification performed in %.3g s\n", t_rcv);
    }
#endif

  free(loc2glob);
  free(mat_elem_up);
  free(mat_elem_low);
  free(mat_elem_diag);
  loc2glob = NULL;
  /***************************************************/
  /* Free Solver internal structures for problem id  */
  /***************************************************/
  MURGE_CALL(MURGE_Clean(id));

  /************************************/
  /* Free Solver internal structures  */
  /************************************/
  MURGE_CALL(MURGE_Finalize());

  free(ax);
  free(vector);
#ifndef X_ARCHi686_mac
  clock_gettime(CLOCK_REALTIME, &tt2);
  t = ( tt2.tv_sec - tt1.tv_sec )
    + (double)( tt2.tv_nsec - tt1.tv_nsec )
    / (double)BILLION;
  MPI_Reduce(&t,&t_rcv, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (proc_id == 0) {
      fprintf(stdout, "Total time %.3g s\n", t_rcv);
    }
#endif

  /** End MPI **/
  MPI_Finalize();

  return 0;
}
