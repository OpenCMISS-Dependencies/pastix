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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas.h>
#include <math.h>
#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#include "common_pastix.h"
#include "getra_cuda.h"
#include "sopalin_compute.h"
#include "sopalin_define.h"


#define STR(x) #x
#ifdef PREC_DOUBLE
#  define FLOAT_EPSILON DBL_EPSILON
#else
#  define FLOAT_EPSILON FLT_EPSILON
#endif

void usage(char * name)
{
  fprintf(stdout, "usage: %s <m> <lda>\n", name);
}

#define READ_INT(m, i) do {                     \
    char *       endptr;                        \
    m = strtol(argv[i], &endptr, 0);            \
    if (argv[i] == endptr)                      \
      {                                         \
        usage(argv[0]);                         \
        return 1;                               \
      }                                         \
  } while(0)

#define FILL(A,n) do {                          \
  int fill_i;                                   \
  for (fill_i = 0; fill_i < n; fill_i++)        \
    A[fill_i] = (PASTIX_FLOAT)(((double)rand())/       \
                        ((double)RAND_MAX));    \
  } while(0)


#define CUDA_CALL(x) do {                                           \
    cudaError_t CUAD_CALL_err;                                      \
    if (cudaSuccess != (CUAD_CALL_err = x))                         \
      {                                                             \
        errorPrint("%s %s (%s,%d)\n",                               \
                   cudaGetErrorString(CUAD_CALL_err),               \
#x, __FILE__,__LINE__);                          \
      }                                                             \
  } while(0)


#ifndef TYPE_COMPLEX
#  define conj(a)  (a)
#  define creal(a) (a)
#  define cimag(a) 0.0
#endif

#define COMPARE_RES(B1,B2) do {                                         \
    int cmpres_i;                                                       \
    int cmpres_p = 0;                                                   \
    double cmpres_norm =  0.0, cmpres_sum = 0.0;                        \
    double cmpres_maxdiff = 0.0;                                        \
    for (cmpres_i = 0; cmpres_i < lda*n; cmpres_i++)                    \
      {                                                                 \
        double cmpres_diff = (double)((B1[cmpres_i]-B2[cmpres_i])*      \
                                      conj(B1[cmpres_i]-B2[cmpres_i])); \
        cmpres_norm += cmpres_diff;                                     \
        cmpres_maxdiff = MAX(cmpres_maxdiff, cmpres_diff);              \
        cmpres_sum  += B1[cmpres_i]*conj(B1[cmpres_i]);                 \
        if (cmpres_p < 10 && sqrt(cmpres_diff) > 0.01)                  \
          {                                                             \
            fprintf(stdout,                                             \
                    "B[%d] %.10g %.10g != %.10g %.10g (%.10g)\n",       \
                    cmpres_i,                                           \
                    creal(B1[cmpres_i]), cimag(B1[cmpres_i]),           \
                    creal(B2[cmpres_i]), cimag(B2[cmpres_i]),           \
                    cmpres_diff);                                       \
            cmpres_p++;                                                 \
          }                                                             \
      }                                                                 \
    fprintf(stdout,                                                     \
            "sqrt( sum( (B1-B2)*conj(B1-B2) )/sum(B1*conj(B1)) ) "      \
            "= %e\n",                                                   \
            sqrt(cmpres_norm/cmpres_sum));                              \
    fprintf(stdout,                                                     \
            "sqrt( max( (B1-B2)*conj(B1-B2) ) )                  "      \
            "= %e\n",                                                   \
            sqrt(cmpres_maxdiff));                                      \
  } while(0)

#define PRINT_TIME(str, time, ops) do {                                 \
    fprintf(stdout,  str " : %.2g s, %.2g GFLOPS\n",                    \
            time, ops/(time*(1<<30)));                                  \
  } while (0)
#define COMPARE_TIME(str, time, ops, time_ref) do {                     \
    fprintf(stdout,  str " : %.2g s, %.2f GFLOPS, Acceleration %.2e\n", \
            time, ops/(time*(1<<30)), time_ref/time);                   \
  } while (0)

double compute_norme(PASTIX_FLOAT * A1, int lda, int ncol, int nrow)
{
  int i,j;
  double norm = 0.0;
  for (i = 0; i <ncol; i++)
    {
      for(j = 0; j < nrow; j++)
        {
          norm = MAX(norm, sqrt(A1[i*lda + j]*conj(A1[i*lda + j])));
        }
    }
  return norm;
}

#define DimTrans           API_CALL(DimTrans)
void DimTrans (PASTIX_FLOAT *A, PASTIX_INT stride, PASTIX_INT size, PASTIX_FLOAT *B);


/* #ifdef TYPE_COMPLEX */
/* #  ifdef PREC_DOUBLE */
/* #    define CU_FLOAT cuDoubleComplex */
/* #  else */
/* #    define CU_FLOAT cuComplex */
/* #  endif */
/* #else */
#  define CU_FLOAT PASTIX_FLOAT
/* #endif */

int
main(int argc, char ** argv)
{
  unsigned int  iseed = (unsigned int)time(NULL);
  int           n;
  int           lda;
  PASTIX_FLOAT        *A;
  PASTIX_FLOAT        *B;
  PASTIX_FLOAT        *B_save;
  PASTIX_FLOAT        *B_res;
  CU_FLOAT     *d_A;
  CU_FLOAT     *d_B;
  Clock         clk;
  Clock         clk_wt;
  PASTIX_FLOAT         alpha = 1.0;
  double        time_CPU;
  double        time_CUDA;
  double        time_CUDA_wt;
  int           ops = n*n;

  if (argc != 3)
    {
      usage(argv[0]);
      return 1;
    }

  READ_INT(n, 1);
  READ_INT(lda, 2);
  srand (iseed);

  MALLOC_INTERN(A,      n*lda, PASTIX_FLOAT);
  MALLOC_INTERN(B,      n*lda, PASTIX_FLOAT);
  MALLOC_INTERN(B_save, n*lda, PASTIX_FLOAT);
  MALLOC_INTERN(B_res,  n*lda, PASTIX_FLOAT);

  FILL(A, n*lda);
  FILL(B, n*lda);
  memcpy(B_save, B, n*lda*sizeof(PASTIX_FLOAT));

  clockInit(&(clk));
  clockStart(&(clk));
  DimTrans(A, lda, n, B);
  clockStop(&(clk));
  time_CPU = clockVal(&(clk));
  PRINT_TIME("GETRA on CPU", time_CPU, ops);

  clockInit(&(clk_wt));
  clockStart(&(clk_wt));
  CUDA_CALL(cudaMalloc((void*)&(d_A),
                       lda*n*sizeof(PASTIX_FLOAT)));
  CUDA_CALL(cudaMemcpy((void*)d_A, A,
                       lda*n*sizeof(PASTIX_FLOAT),
                       cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMalloc((void*)&(d_B),
                       lda*n*sizeof(PASTIX_FLOAT)));
  CUDA_CALL(cudaMemcpy((void*)d_B, B_save,
                       lda*n*sizeof(PASTIX_FLOAT),
                       cudaMemcpyHostToDevice));
  clockInit(&(clk));
  clockStart(&(clk));
  getra_cuda(d_A, lda,
             d_B, lda, n);
  clockStop(&(clk));

  CUDA_CALL(cudaMemcpy((void*)B_res, d_B,
                       lda*n*sizeof(PASTIX_FLOAT),
                       cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaFree(d_A));
  CUDA_CALL(cudaFree(d_B));
  clockStop(&(clk_wt));

  time_CUDA    = clockVal(&(clk));
  time_CUDA_wt = clockVal(&(clk_wt));

  COMPARE_TIME("GETRA on GPU",
               time_CUDA, ops, time_CPU);
  COMPARE_TIME("GETRA on GPU with transfer",
               time_CUDA_wt, ops, time_CPU);
  COMPARE_RES(B, B_res);

  memFree_null(A);
  memFree_null(B);
  memFree_null(B_save);
  memFree_null(B_res);

  return EXIT_SUCCESS;
}
