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
/**
 *  File: murge_sequence.c
 *
 *  Example that generate A Laplacian with multiple
 *  degrees of freedom and solves it.
 *
 *  That example only works with the -DDISTRIBUTED option.
 *
 *  The laplacian looks like :
 *  >    1       0                      ... 0
 *  > -2(n-1)  (n-1)  -2(n-1)   0       ... 0
 *  >    0    -2(n-1)  (n-1)  -2(n-1) 0 ... 0
 *  >              ....
 *  >              ....
 *  >    0 ...              -2(n-1)  n-1  -2(n-1)
 *  >    0 ...                        0     1
 *
 *  The Right-hand-side member is defined by
 *    $$RHS_i = - 4\pi^2 sin(2\pi x)$$
 *
 *  The solution should be $$X_i = sin(2\pi x)$$
 *
 *  The solution is stored in the file "result.plot" that can be ploted using
 *  > gnuplot> plot "./result.plot" u 2:3 t "reference", "./result.plot" \
 *             u 2:4 t "Result, first degree of freedom"
 *
 *
 *  Usage:
 *  > ./murge <size> <DofNbr>
 *
 *  Authors:
 *    Xavier LACOSTE - lacoste@labri.fr
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef TYPE_COMPLEX
#  include <complex.h>
#endif
#ifdef _OPENMP
#  include <omp.h>
#endif

#include "mpi.h"
#include "pastix.h"
#include "murge.h"

/**
 *  Function GetCoef
 *
 *  Computes the value for a given coefficient.
 *
 *  Parameters:
 *   val  - Value to set
 *   i    - Row of the coefficient.
 *   j    - Column of the coefficient.
 *   xmin - Minimum value of the interval.
 *   xmax - Maximum value of the interval.
 *   n    - Number of points in the interval.
 */
static
void GetCoef(COEF * matElem,INTS i, INTS j, double xmin, double xmax,
             INTS n, INTS dof) {
  double dx_1;
  COEF val;
  int ii, jj;

  dx_1 = (n-1.)/ (xmax - xmin);

  val = 0.;

  if (i==j) {
    if (i==1 || i == n) {
      /* Boundary Condition (Dirichlet) */
      val = 1;
    } else {
      /*  Interior diagonnal part */
      val = -2 * dx_1;
    }
  } else {
    val = dx_1;
  }
  for (ii = 0; ii <dof; ii++) {
    for (jj = 0; jj< dof; jj++) {
      if (ii == jj) {
        matElem[ii+jj*dof] = val;
      }
      else {
        matElem[ii+jj*dof] = 0;
      }
    }
  }
}

/**
 *  Function GetRhs
 *
 *  computes the value of a coefficient of the Right-hand-side member.
 *
 *  Parameters:
 *    val  - Value to set.
 *    i    - Index of the value.
 *    xmin - Minimum value of the interval.
 *    xmax - Maximum value of the interval.
 *    n    - Number of points in the interval.
 */
static
void GetRhs(COEF * val, INTS i, double xmin, double xmax, INTS n) {
  double dx , x,Pi;

  dx = (xmax - xmin) / (n-1.);
  x = xmin + (i-1)*dx;

  Pi =acos(-1.);

  *val = -4*Pi*Pi*sin(2*Pi*x);
  /* Boundary Condition (Dirichlet) */
  if (i == n || i==1) {
    *val = 0.;
  } else {
    *val = dx*(*val);
  }
}

/**
 *  Function: store
 *
 *  Write the solution into a file result.plot.
 *
 *  The file contains :
 *  > k x sin(2 \pi x) sol(k:k+dofnbr-1)
 *  Where k goes from 1 to n*dofnbr, dofnbr by dofnbr and x goes from
 *  xmin to xmax, with a step of (xmax - xmin) / (n-1).
 *
 *  Parameters:
 *    sol  - The solution of the problem
 *    xmin - The minimum value of the interval.
 *    xmax - The maximum value of the interval.
 *    dof  - The Number of degree of freedom.
 */
static
void store(COEF * sol, double xmin, double xmax, INTS n, INTS dof) {
  double x,dx,Pi2;
  INTL i,k,l;
  FILE * f;

  x = xmin;
  dx = (xmax - xmin) / (n-1.);
  Pi2 = 2.*acos(-1.);

  f = fopen("result.plot", "w");
  k = 1;
  for (i = 0; i< n; i++) {
    fprintf(f,"%ld\t%.15g\t%.15g",(long)k,x,sin(Pi2*x));
    for (l = 0; l < dof; l++)
      {
#ifdef TYPE_COMPLEX
        fprintf(f, "\t%.15g %.15g", creal(sol[i*dof+l]), cimag(sol[i*dof+l]));
#else
        fprintf(f, "\t%.15g", sol[i*dof+l]);
#endif
      }
    fprintf(f,"\n");
    k = k + dof;
    x = x + dx;
  }
  fclose(f);
}

int main( int argc, char ** argv) {

  INTS ierr;
  int me, NTasks, required, provided;
  INTS id, m;
  /* CSC Data */
  INTS n, dof;
  INTL nnzeros, edgenbr;
  COEF val;
  /* Local Data */
  INTS localnodenbr;
  INTS * nodelist;
  INTS new_localnodenbr;
  INTS * new_nodelist;
  INTS root;
  INTS base;
  COEF * lrhs;
  COEF * globrhs;
  COEF * globrhs_recv;
  COEF * globx;
  COEF * globprod;
  /* Other data */
  COEF * matElem;
  double prec, xmin, xmax, sum1, sum2;
  INTS i, j, k;
  INTS solver;
  INTS zero=0;
  INTS one=1;
  INTS nb_threads;
  INTS   id_seq;

  root = -1;
  base = 1;

  required=MPI_THREAD_MULTIPLE;
  MPI_Init_thread(&argc, &argv, required, &provided);

  MPI_Comm_size(MPI_COMM_WORLD, &NTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  n = dof = 0;
  if (argc >= 2) {
    n = atoi(argv[1]);
    dof = atoi(argv[2]);
  } else {
    if (me == 0) {
      fprintf(stderr, "Usage: %s <size> <DofNumber>\n", argv[0]);
      return 1;
    }
  }

  xmin = 0.0;
  xmax = 1.0;

  /* Starting MURGE*/
  ierr = MURGE_Initialize(one);
  if (ierr != MURGE_SUCCESS) {
    fprintf(stderr, "Error %ld in MURGE_Initialize\n", (long)ierr);
    return 1;
  }
  id = 0;

  /* Set Options */
  prec = 1e-7;
  /*
    Call MURGE_Get_Solver(solver)
  */
  solver = MURGE_SOLVER_PASTIX;




  if ( solver == MURGE_SOLVER_PASTIX ) {
    MURGE_SetDefaultOptions(id, zero);
    MURGE_SetOptionINT(id, IPARM_VERBOSE, API_VERBOSE_NO);
    MURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_YES);
    nb_threads = 1;
#ifdef _OPENMP
#pragma omp parallel shared(nb_threads)
    {
      nb_threads = omp_get_num_threads();
    }
#endif /* _OPENMP */

    if (me == 0) {
      fprintf(stdout, "Running on %ld threads and %d MPI Tasks\n",
              (long)nb_threads, NTasks);
    }
    MURGE_SetOptionINT(id, IPARM_THREAD_NBR, nb_threads);
  } else if (solver == MURGE_SOLVER_HIPS) {
#ifdef HIPS
    if ( method == 1 ) {
      MURGE_SetDefaultOptions(id, HIPS_ITERATIVE);
    } else {
      MURGE_SetDefaultOptions(id, HIPS_HYBRID);
      MURGE_SetOptionINT(id, HIPS_PARTITION_TYPE, zero);
      MURGE_SetOptionINT(id, HIPS_DOMSIZE, domsize);
    }
    MURGE_SetOptionINT(id, HIPS_SYMMETRIC, zero);
    MURGE_SetOptionINT(id, HIPS_LOCALLY, zero);
    MURGE_SetOptionINT(id, HIPS_ITMAX, itmax);
    MURGE_SetOptionINT(id, HIPS_KRYLOV_RESTART, restart);
    MURGE_SetOptionINT(id, HIPS_VERBOSE, verbose);
    MURGE_SetOptionINT(id, HIPS_DOMNBR, NTasks);
    MURGE_SetOptionINT(id, HIPS_CHECK_GRAPH, one);
    MURGE_SetOptionINT(id, HIPS_CHECK_MATRIX, one);
#endif
  }
  MURGE_SetOptionINT(id, MURGE_IPARAM_DOF, dof);
  MURGE_SetOptionINT(id, MURGE_IPARAM_SYM, MURGE_BOOLEAN_FALSE);
  MURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL, base);

  MURGE_SetOptionREAL(id, MURGE_RPARAM_EPSILON_ERROR, prec);
  /* Set the graph : all processor enter some edge of the
     graph that corresponds to non-zeros location in the matrix */

  /****************************************
   ** Enter the matrix non-zero pattern  **
   ** you can use any distribution       **
   ****************************************/

  /* this processor enters the A(myfirstrow:mylastrow, *)
     part of the matrix non-zero pattern */
  if (me == 0) {
    edgenbr = 3*n-4;

    MURGE_GraphBegin(id, n, edgenbr);

    /* Dirichlet boundary condition */
    MURGE_GraphEdge(id, one, one);
    MURGE_GraphEdge(id, n, n);

    /* Interior */
    for (i = 2; i < n; i++) {
      for (j = -1; j <= 1; j++) {
        MURGE_GraphEdge(id, i, i+j);
        /* if (j != 0) {
           MURGE_GraphEdge(id, j+i, i);
           } */
      }
    }
  } else {
    edgenbr = 0;
    MURGE_GraphBegin(id, n, edgenbr);
  }
  MURGE_GraphEnd(id);


  /*  Get Local nodes */
  MURGE_GetLocalNodeNbr(id, &localnodenbr);
  nodelist = (INTS*)malloc(localnodenbr*sizeof(INTS));
  MURGE_GetLocalNodeList(id, nodelist);

  new_localnodenbr = localnodenbr;
  new_nodelist     = nodelist;

  if (NTasks > 1)
    {
      /* move column 1 to next proc */
      int    found = 0;
      int    owner, owner_rcv, index;

      for (m = 0; m < localnodenbr; m++) {
        i = nodelist[m];
        if (i == 1)
          found = 1;
      }
      if (found == 1)
        owner = me;
      else
        owner = 0;
      MPI_Allreduce(&owner, &owner_rcv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (owner_rcv == me)
        {
          new_localnodenbr--;
          new_nodelist = (INTS*)malloc(new_localnodenbr*sizeof(INTS));
          index = 0;
          for (m = 0; m < localnodenbr; m++)
            {
              i = nodelist[m];
              if (i != 1)
                new_nodelist[index++] = nodelist[m];
            }
        }
      if (owner_rcv == (me+1)%NTasks)
        {
          index = 0;
          new_localnodenbr++;
          new_nodelist = (INTS*)malloc(new_localnodenbr*sizeof(INTS));
          new_nodelist[index++] = 1;
          for (m = 0; m < localnodenbr; m++)
            {
              new_nodelist[index++] = nodelist[m];
            }
        }
    }
  /* Compute the number of non-zeros; */
  nnzeros = 0;
  for (m = 0; m < new_localnodenbr; m++) {
    i = new_nodelist[m];

    if (i == 1 || i == n) {
      /*  Boundaries */
      nnzeros = nnzeros + 1;
    } else {
      /*  Interior */
      for (k = -1; k <= 1; k++) {
        nnzeros = nnzeros + 1;
      }
    }
  }
  /*  We are using dof so a non zero is in fact a block of size dof**2 */
  edgenbr = nnzeros;
  nnzeros = nnzeros * dof*dof;

  /* You can enter only coefficient (i,j) that are in A(nodelist, nodelist)
     on this processor */

  /* We enter the lower and upper triangular part of the matrix so sym = 0 */

  /* matElem is the identity matrix of size 'dof' stored by line */
  {
    INTS * ROWs  = malloc(edgenbr*sizeof(INTS));
    INTS * COLs  = malloc(edgenbr*sizeof(INTS));
    int    index = 0;

    for (m = 0; m < new_localnodenbr; m++) {
      i = new_nodelist[m];
      if ( i == 1 || i == n ) {
        /*  Boundaries */
        ROWs[index] = i;
        COLs[index] = i;
        index++;
      } else {
        for (k = -1; k <= 1; k++) {
          ROWs[index] = i+k;
          COLs[index] = i;
          index++;
        }
      }
    }
    ierr = MURGE_AssemblySetSequence(id, edgenbr, ROWs, COLs,
                                     MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW,
                                     MURGE_ASSEMBLY_FOOL, MURGE_BOOLEAN_TRUE,
                                     &id_seq);
    if (ierr != MURGE_SUCCESS) {
      fprintf(stderr, "Error %ld in MURGE_AssemblySetSequence\n", (long)ierr);
      exit(1);
    }
    free(ROWs);
    free(COLs);
  }

  {
    COEF * values = malloc(nnzeros*sizeof(COEF));
    int    index  = 0;
    matElem = (COEF*)malloc(dof*dof*sizeof(COEF));
    for (m = 0; m < new_localnodenbr; m++) {
      i = new_nodelist[m];
      if ( i == 1 || i == n ) {
        /*  Boundaries */
        GetCoef(matElem, i,i,xmin,xmax,n, dof);
        memcpy(&(values[dof*dof*index]),
               matElem, dof*dof*sizeof(COEF));
        index++;
      } else {
        for (k = -1; k <= 1; k++) {
          GetCoef(matElem,i+k,i,xmin,xmax,n, dof);
          memcpy(&(values[dof*dof*index]),
                 matElem, dof*dof*sizeof(COEF));
          index++;
        }
      }
    }

    free(matElem);
    ierr = MURGE_AssemblyUseSequence(id, id_seq, values);
    if (ierr != MURGE_SUCCESS) {
      fprintf(stderr, "Error %ld in MURGE_AssemblySetSequence\n", (long)ierr);
      exit(1);
    }
    free(values);
  }

  /* We build the rhs */
  lrhs = (COEF*)malloc(localnodenbr*dof*sizeof(COEF));
  globrhs = (COEF*)malloc(n*dof*sizeof(COEF));
  for (k = 0; k < n*dof; k++)
    globrhs[k] = 0.0;

  for (m = 0; m < localnodenbr; m++) {
    GetRhs(&val,nodelist[m],xmin,xmax,n);
    for (k = 0; k < dof; k++)
      globrhs[(nodelist[m]-1)*dof+k] = val;
    for (k = 0; k < dof; k++)
      lrhs[m*dof+k] = val;
  }

  globrhs_recv = (COEF*)malloc(n*dof*sizeof(COEF));
  MPI_Allreduce(globrhs, globrhs_recv, n*dof, MURGE_MPI_COEF,
                MPI_SUM, MPI_COMM_WORLD);
  free(globrhs);
  MURGE_SetLocalRHS(id, lrhs, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW);

  /* Get the global solution */
  globx = (COEF*)malloc(n*dof*sizeof(COEF));
  MURGE_GetGlobalSolution(id, globx, root);

  MURGE_SetGlobalRHS(id, globx, -one, MURGE_ASSEMBLY_OVW);
  globprod = (COEF*)malloc(n*dof*sizeof(COEF));
  MURGE_GetGlobalProduct(id, globprod, -one);
  sum1 = 0;
  sum2 = 0;
  for (k = 0; k < n*dof; k++) {
    sum1 += globprod[k]*globprod[k];
    sum2 += (globprod[k] - globrhs_recv[k])*(globprod[k]-globrhs_recv[k]);
  }
  fprintf(stdout, "||AX - B||/||AX||  : %.15g\n", sqrt(sum2/sum1));

  /* Store in a file */
  if (me == 0)
    store(globx,xmin,xmax,n,dof);

  {
    int iter;
    INTS  n_coefs = n/NTasks;
    INTS *coef_idx;
    COEF *coef_vals;

    if (me < n%NTasks)
      n_coefs++;

    fprintf(stdout, "Now using MURGE_SetRHS and MURGE_ASSEMBLY_FOOL\n");
    coef_idx  = malloc(n_coefs*sizeof(INTS));
    coef_vals = malloc(n_coefs*dof*sizeof(COEF));
    /* cyclic distribution of RHS */
    for (iter = 0; iter < n_coefs; iter++)
      {
        coef_idx[iter]  = me + iter*NTasks + 1; /* baseval == 1 */
        for (k = 0; k < dof; k++)
          coef_vals[iter*dof+k] = globrhs_recv[(me + iter*NTasks)*dof + k];

      }

    MURGE_SetRHS(id, n_coefs, coef_idx, coef_vals, MURGE_ASSEMBLY_OVW,
                 MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_FOOL);

    free(coef_vals);
    free(coef_idx);
    MURGE_GetGlobalSolution(id, globx, root);

    MURGE_SetGlobalRHS(id, globx, -one, MURGE_ASSEMBLY_OVW);
    MURGE_GetGlobalProduct(id, globprod, -one);
    sum1 = 0;
    sum2 = 0;
    for (k = 0; k < n*dof; k++) {
      sum1 += globprod[k]*globprod[k];
      sum2 += (globprod[k] - globrhs_recv[k])*(globprod[k]-globrhs_recv[k]);
    }
    fprintf(stdout, "||AX - B||/||AX||  : %.15g\n", sqrt(sum2/sum1));
  }

  /* I'm Free  */
  MURGE_Clean(id);
  MURGE_Finalize();
  MPI_Finalize();

  free(nodelist);
  free(lrhs);
  free(globx);
  free(globprod);
  free(globrhs_recv);
  return 0;
}
