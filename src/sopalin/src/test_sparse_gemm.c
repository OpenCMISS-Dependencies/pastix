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
#ifdef WITH_MAGMABLAS
#include <magmablas.h>
#endif
#include <math.h>

#include "common_pastix.h"
#include "sparse_gemm.h"
#include "sopalin_compute.h"
#include "pastix_cuda_helper.h"

#if (CUDA_SM_VERSION >= 20)
#include "sparse_gemm_fermi.h"
#endif

#define TO_STR(x) STR(x)
#define STR(x) #x
#ifdef PREC_DOUBLE
#  define FLOAT_EPSILON DBL_EPSILON
#else
#  define FLOAT_EPSILON FLT_EPSILON
#endif

#define SAVE_MIN_TIME(time)                                             \
  do {                                                                  \
    if (run_idx == 0)                                                   \
      min_ ## time [run2_idx] = time;                                   \
    else                                                                \
      min_ ## time [run2_idx] = MIN(time, min_ ## time [run2_idx]);     \
  } while (0)

void usage(char * name)
{
  fprintf(stdout, "usage: %s <nrowsA> <ncolsA> <nrowsA11> [nruns_by_dist n_dist]\n", name);
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
    int fill_i;                                 \
    for (fill_i = 0; fill_i < n; fill_i++)      \
      A[fill_i] = (PASTIX_FLOAT)(((double)rand())/     \
                          ((double)RAND_MAX));  \
  } while(0)


#define CUDA_CALL(x) do {                               \
    cudaError_t CUAD_CALL_err;                          \
    if (cudaSuccess != (CUAD_CALL_err = x))             \
      {                                                 \
        errorPrint("%s %s (%s,%d)\n",                   \
                   cudaGetErrorString(CUAD_CALL_err),   \
                   #x, __FILE__,__LINE__);              \
      }                                                 \
  } while(0)


#define CUDA_SPARSE_GEMM(TRANSA, TRANSB,                                \
                         dimi, dimj, dima,                              \
                         alpha,                                         \
                         A,  stride_A,                                  \
                         B,  stride_B,                                  \
                         beta,                                          \
                         C, stride_C,                                   \
                         blocknbr, blocktab,                            \
                         fblocknbr, fblocktab)                          \
  do {                                                                  \
    magmablas_sparse_gemm_kernel_N_T_64_16_4_16_4((int)dimi,            \
                                                  (int)dimj,            \
                                                  (int)dima,            \
                                                  (float)alpha,         \
                                                  A,                    \
                                                  (int)stride_A,        \
                                                  B,                    \
                                                  (int)stride_B,        \
                                                  (float)beta,          \
                                                  C,                    \
                                                  (int)stride_C,        \
                                                  blocknbr,             \
                                                  blocktab,             \
                                                  fblocknbr,            \
                                                  fblocktab);           \
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
    for (cmpres_i = 0; cmpres_i < ldb*n; cmpres_i++)                    \
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
    fprintf(stdout, "%d.%d: norm2 = %e\n",                              \
            run2_idx, run_idx, sqrt(cmpres_norm/cmpres_sum));           \
    fprintf(stdout, "%d.%d: normM = %e\n",                              \
            run2_idx, run_idx, sqrt(cmpres_maxdiff));                   \
    fprintf(stdout, "%d.%d: normM/(sum(normA,At,C)*MAX(m,n)*EPS) = %e\n", \
            run2_idx, run_idx, sqrt(cmpres_maxdiff)/                    \
            ((norm_A+norm_At+norm_B)*n*FLOAT_EPSILON));                 \
  } while(0)

#define PRINT_TIME(str, time, ops) do {                         \
    fprintf(stdout,  "%d.%d: " str " : %.2g s, %.2g GFLOPS\n",  \
            run2_idx, run_idx, time, ops/(time*(1<<30)));       \
  } while (0)
#define COMPARE_TIME(str, time, ops, time_ref) do {                     \
    fprintf(stdout,                                                     \
            "%d.%d: " str " : %.2g s, %.2f GFLOPS, Acceleration %.2f\n", \
            run2_idx, run_idx, time, ops/(time*(1<<30)), time_ref/time); \
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

#ifdef TYPE_COMPLEX
#  ifdef WITH_MAGMABLAS
#    ifdef PREC_DOUBLE
#      define CUBLAS_GEMM magmablas_zgemm
#    else
#      define CUBLAS_GEMM magmablas_cgemm
#    endif
#  else
#    ifdef PREC_DOUBLE
#      define CUBLAS_GEMM cublasZgemm
#    else
#      define CUBLAS_GEMM cublasCgemm
#    endif
#  endif
#else
#  ifdef WITH_MAGMABLAS
#    ifdef PREC_DOUBLE
#      define CUBLAS_GEMM magmablas_dgemm
#    else
#      define CUBLAS_GEMM magmablas_sgemm
#    endif
#  else
#    ifdef PREC_DOUBLE
#      define CUBLAS_GEMM cublasDgemm
#    else
#      define CUBLAS_GEMM cublasSgemm
#    endif
#  endif
#endif
int
main(int argc, char ** argv)
{

  unsigned int  iseed      = (unsigned int)time(NULL);
  long          m,n,k,nrowsA11;
  int  *        blocks;
  int  *        fblocks;
  int  *        d_blocks;
  int  *        d_fblocks;
  int           before_size;
  int           after_size;
  int           nb_blocks_in;
  int           b;
  int           block_size;
  PASTIX_FLOAT        *A1, *D;
  int           lda;
  PASTIX_FLOAT        *B0, *B_ref, *B_res;
  PASTIX_FLOAT        *work;
  int           worksize;
  int           ldb;
  int           ldd;
  Clock         clk;
  Clock         clk_wt;

  PASTIX_FLOAT alpha = -1.0;
  PASTIX_FLOAT beta  = 1.0;

  CU_FLOAT cu_alpha;
  CU_FLOAT cu_beta;
  CU_FLOAT *d_A, *d_B, *d_D;

  double time_sparse_CPU,          *min_time_sparse_CPU;
  double time_sparse_GPU,          *min_time_sparse_GPU;
  double time_sparse_GPU_wt,       *min_time_sparse_GPU_wt;
#if (CUDA_SM_VERSION >= 20)
  double time_sparse_GPU_FERMI,    *min_time_sparse_GPU_FERMI;
  double time_sparse_GPU_FERMI_wt, *min_time_sparse_GPU_FERMI_wt;
  double time_sparse_GPU_FERMI_GEMDM,    *min_time_sparse_GPU_FERMI_GEMDM;
  double time_sparse_GPU_FERMI_GEMDM_wt, *min_time_sparse_GPU_FERMI_GEMDM_wt;
#endif
  double time_dense_CPU,           *min_time_dense_CPU;
  double time_dense_GPU,           *min_time_dense_GPU;
  double time_dense_GPU_wt,        *min_time_dense_GPU_wt;
  double ops;
  double norm_A;
  double norm_B;
  double norm_At;
  int    nruns = 1, nruns2 = 1, run_idx, run2_idx;

  cu_alpha = CU_FLOAT_INIT(creal(alpha), cimag(alpha));
  cu_beta  = CU_FLOAT_INIT(creal(beta),  cimag(beta));
  srand (iseed);

  if (argc != 4 && argc != 6)
    {
      usage(argv[0]);
      return 1;
    }

  READ_INT(m, 1);
  READ_INT(k, 2);
  READ_INT(nrowsA11, 3);
  assert(nrowsA11 <= m);
  if (argc == 6)
    {
      READ_INT(nruns, 4);
      READ_INT(nruns2, 5);
    }
  MALLOC_INTERN(min_time_sparse_CPU,    nruns2, double);
  MALLOC_INTERN(min_time_sparse_GPU,    nruns2, double);
  MALLOC_INTERN(min_time_sparse_GPU_wt, nruns2, double);
#if (CUDA_SM_VERSION >= 20)
  MALLOC_INTERN(min_time_sparse_GPU_FERMI,    nruns2, double);
  MALLOC_INTERN(min_time_sparse_GPU_FERMI_wt, nruns2, double);
  MALLOC_INTERN(min_time_sparse_GPU_FERMI_GEMDM,    nruns2, double);
  MALLOC_INTERN(min_time_sparse_GPU_FERMI_GEMDM_wt, nruns2, double);
#endif
  MALLOC_INTERN(min_time_dense_CPU,    nruns2, double);
  MALLOC_INTERN(min_time_dense_GPU,    nruns2, double);
  MALLOC_INTERN(min_time_dense_GPU_wt, nruns2, double);

  for (run2_idx = 0; run2_idx < nruns2; run2_idx++)
    {
      int           nb_blocks  = 0;
      int           size       = 0;
      int           last       = 0;
      int           nb_fblocks = 0;
      /* Build a sparse block column */
      /* The total size of all blocks must be m,
       so there are m blocks max. */
      MALLOC_INTERN(blocks, 2*m, int);
      do {
        block_size = (int)(((double)m*(double)rand())/
                           ((double)RAND_MAX)/10)+1;
        if (nb_blocks == 0)
          block_size =nrowsA11;

        if (size + block_size > m)
          block_size = m - block_size;

        blocks[2*nb_blocks]   = last + (int)(((double)m*(double)rand())/
                                             ((double)RAND_MAX)) +1;
        blocks[2*nb_blocks+1] = blocks[2*nb_blocks] + block_size - 1;
        if ( size + blocks[2*nb_blocks+1] - blocks[2*nb_blocks] + 1 > m)
          blocks[2*nb_blocks+1] = blocks[2*nb_blocks] + m - size - 1;
        fprintf(stdout, "block [%d, %d] ([%d, %d])\n",
                blocks[2*nb_blocks], blocks[2*nb_blocks+1],
                size, size + block_size - 1);
        size += blocks[2*nb_blocks+1] - blocks[2*nb_blocks] + 1;
        last = blocks[2*nb_blocks+1];
        nb_blocks++;
      } while(size < m);

      MALLOC_INTERN(fblocks, 2*m, int);
      b = 0;

      /* Build a facing sparse block column */
      size = 0;
      do {
        nb_blocks_in = (int)((nb_blocks*(double)rand())/
                             ((double)RAND_MAX));

        if (nb_fblocks == 0)
          before_size = blocks[0];
        else
          before_size = blocks[2*b]  - fblocks[2*nb_fblocks-1] - 1;
        before_size = (int)rint(((double)before_size*(double)rand())/
                                ((double)RAND_MAX));
        if (nb_blocks_in == 0 || b + nb_blocks_in >= nb_blocks)
          {
            nb_blocks_in = nb_blocks - b;
            after_size = 100;
          }
        else
          after_size = blocks[2*(b+nb_blocks_in)] -
            blocks[2*(b+nb_blocks_in-1)+1]-1;
        after_size = (int)rint(((double)after_size*(double)rand())/
                               ((double)RAND_MAX));

        fblocks[2*nb_fblocks]   = blocks[2*b] - before_size;
        block_size = blocks[2*(b+nb_blocks_in-1)+1] - blocks[2*b] + 1;
        fblocks[2*nb_fblocks+1] = fblocks[2*nb_fblocks] +
          before_size + block_size + after_size - 1;
        fprintf(stdout, "fblock [%d, %d] ([%d, %d]) %d %d %d %d\n",
                fblocks[2*nb_fblocks], fblocks[2*nb_fblocks+1],
                size , size + before_size + block_size + after_size - 1,
                nb_blocks_in, before_size, block_size, after_size);
	size += before_size + block_size + after_size; 
        b += nb_blocks_in;
        nb_fblocks++;
      } while(b != nb_blocks);

      /* check */
      {
        int fb = 0;
        for (b = 0; b < nb_blocks; b++)
          {
            while (!(blocks[2*b] >= fblocks[2*fb] &&
                     blocks[2*b+1] <= fblocks[2*fb+1]))
              {
                fb++;
                assert(fb < nb_fblocks);
              }
          }
      }

      /* allocate and fill the matrices */
      lda = m;
      ldd = lda;
      ldb = fblocks[2*(nb_fblocks-1)+1]-fblocks[0]+1;
      n = blocks[1] - blocks[0] +1;
      MALLOC_INTERN(A1, lda*k, PASTIX_FLOAT);
      MALLOC_INTERN(D,  k*ldd,     PASTIX_FLOAT);
      MALLOC_INTERN(B0, ldb*n, PASTIX_FLOAT);
      memset(B0, 0, ldb*n*sizeof(PASTIX_FLOAT));
      MALLOC_INTERN(B_ref, ldb*n, PASTIX_FLOAT);
      MALLOC_INTERN(B_res, ldb*n, PASTIX_FLOAT);
      MALLOC_INTERN(work, m*n, PASTIX_FLOAT);
      worksize = m*n;
      FILL(A1, lda*k);
      FILL(B0, ldb*n);
      norm_A  = compute_norme(A1, lda, k, lda);
      norm_At = compute_norme(A1, lda, k, blocks[1]-blocks[0]+1);
      norm_B  = compute_norme(B0, ldb, n, ldb);
      {
        int i,j;
        for (i = 0; i < ldd; i++)
	  for (j = 0; j < k; j++)
	    D[j*ldd + i] = 1.0*(i==j);
      }
      
      for (run_idx = 0; run_idx < nruns; run_idx++)
        {
          memcpy(B_ref, B0, ldb*n*sizeof(PASTIX_FLOAT));
          clockInit(&(clk));
          clockStart(&(clk));
          sparse_gemm_cpu("N", "C",
                          m, n, k,
                          alpha,
                          A1, lda,
                          A1, lda,
                          beta,
                          B_ref, ldb,
                          nb_blocks,  blocks,
                          nb_fblocks, fblocks,
                          work, worksize);
          clockStop(&(clk));
          time_sparse_CPU = clockVal(&(clk));
          ops = m*n*k*2;
          SAVE_MIN_TIME(time_sparse_CPU);
          PRINT_TIME("sparse GEMM on CPU", time_sparse_CPU, ops);





          memcpy(B_res, B0, ldb*n*sizeof(PASTIX_FLOAT));
          if (CUDA_SUCCESS != cuInit(0))
            {
              errorPrint("cuInit()");
              assert(0);
            }
          CUDA_CALL(cudaSetDevice(0));

          CUDA_CALL(cudaMalloc((void*)&(d_blocks),
                               2*nb_blocks*sizeof(int)));
          CUDA_CALL(cudaMemcpy((void*)d_blocks, blocks,
                               2*nb_blocks*sizeof(int),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_fblocks),
                               2*nb_fblocks*sizeof(int)));
          CUDA_CALL(cudaMemcpy((void*)d_fblocks, fblocks,
                               2*nb_fblocks*sizeof(int),
                               cudaMemcpyHostToDevice));
#if (CUDA_SM_VERSION >= 20 || !(defined PREC_DOUBLE && defined TYPE_COMPLEX))
          if (nruns <  10)
            fprintf(stdout, ">>> %s <<<\n",
                    TO_STR(PASTIX_PREFIX_F(magmablas_sparse_gemm_kernel_N_T_64_16_4_16_4)));

          clockInit(&(clk_wt));
          clockStart(&(clk_wt));
          CUDA_CALL(cudaMalloc((void*)&(d_A),
                               lda*k*sizeof(PASTIX_FLOAT)));
          CUDA_CALL(cudaMemcpy((void*)d_A, A1,
                               lda*k*sizeof(PASTIX_FLOAT),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_B),
                               ldb*n*sizeof(PASTIX_FLOAT)));
          CUDA_CALL(cudaMemcpy((void*)d_B, B0,
                               ldb*n*sizeof(PASTIX_FLOAT),
                               cudaMemcpyHostToDevice));

          CUDA_CALL(cudaThreadSynchronize());
          clockInit(&(clk));
          clockStart(&(clk));
          CUDA_SPARSE_GEMM("N", "C",
                           m, n, k,
                           alpha,
                           (PASTIX_FLOAT*)d_A, lda,
                           (PASTIX_FLOAT*)d_A, lda,
                           beta,
                           (PASTIX_FLOAT*)d_B, ldb,
                           nb_blocks,  d_blocks,
                           nb_fblocks, d_fblocks);
          CUDA_CALL(cudaThreadSynchronize());

          clockStop(&(clk));

          CUDA_CALL(cudaMemcpy((void*)B_res, d_B,
                               ldb*n*sizeof(PASTIX_FLOAT),
                               cudaMemcpyDeviceToHost));
          CUDA_CALL(cudaFree(d_A));
          CUDA_CALL(cudaFree(d_B));
          clockStop(&(clk_wt));

          time_sparse_GPU = clockVal(&(clk));
          SAVE_MIN_TIME(time_sparse_GPU);
          time_sparse_GPU_wt = clockVal(&(clk_wt));
          SAVE_MIN_TIME(time_sparse_GPU_wt);
          COMPARE_TIME("sparse GEMM on GPU (SM < 20)",
                       time_sparse_GPU, ops, time_sparse_CPU);
          COMPARE_TIME("sparse GEMM on GPU (SM < 20) with transfer",
                       time_sparse_GPU_wt, ops, time_sparse_CPU);
          COMPARE_RES(B_ref, B_res);

#endif




#if (CUDA_SM_VERSION >= 20)
          if (nruns <  10)
            fprintf(stdout, ">>> %s <<< (FERMI)\n",
                    TO_STR(GENERATE_SM_VERSION_NAME(gemm)));

          /* FERMI kernel */
          clockInit(&(clk_wt));
          clockStart(&(clk_wt));
          CUDA_CALL(cudaMalloc((void*)&(d_A),
                               lda*k*sizeof(PASTIX_FLOAT)));
          CUDA_CALL(cudaMemcpy((void*)d_A, A1,
                               lda*k*sizeof(PASTIX_FLOAT),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_B),
                               ldb*n*sizeof(PASTIX_FLOAT)));
          CUDA_CALL(cudaMemcpy((void*)d_B, B0,
                               ldb*n*sizeof(PASTIX_FLOAT),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaThreadSynchronize());
          clockInit(&(clk));
          clockStart(&(clk));
          GENERATE_SM_VERSION_NAME(gemm)('N',
                                         'T',
                                         m, n, k,
                                         cu_alpha,
                                         d_A, lda,
                                         d_A, lda,
                                         cu_beta,
                                         d_B, ldb,
                                         nb_blocks,  d_blocks,
                                         nb_fblocks, d_fblocks,
                                         0 );
          CUDA_CALL(cudaThreadSynchronize());
          clockStop(&(clk));
          CUDA_CALL(cudaMemcpy((void*)B_res, d_B,
                               ldb*n*sizeof(PASTIX_FLOAT),
                               cudaMemcpyDeviceToHost));
          CUDA_CALL(cudaFree(d_A));
          CUDA_CALL(cudaFree(d_B));
          clockStop(&(clk_wt));
          time_sparse_GPU_FERMI = clockVal(&(clk));
          SAVE_MIN_TIME(time_sparse_GPU_FERMI);
          time_sparse_GPU_FERMI_wt = clockVal(&(clk_wt));
          SAVE_MIN_TIME(time_sparse_GPU_FERMI_wt);
          COMPARE_TIME("sparse GEMM on GPU (Fermi)",
                       time_sparse_GPU_FERMI, ops, time_sparse_CPU);
          COMPARE_TIME("sparse GEMM on GPU (Fermi) with transfert",
                       time_sparse_GPU_FERMI_wt, ops, time_sparse_CPU);
          COMPARE_RES(B_ref,B_res);





          /* FERMI kernel */
          clockInit(&(clk_wt));
          clockStart(&(clk_wt));
          CUDA_CALL(cudaMalloc((void*)&(d_A),
                               lda*k*sizeof(PASTIX_FLOAT)));
          CUDA_CALL(cudaMemcpy((void*)d_A, A1,
                               lda*k*sizeof(PASTIX_FLOAT),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_D),
                               k*ldd*sizeof(PASTIX_FLOAT)));
          CUDA_CALL(cudaMemcpy(d_D, D,
                               k*ldd*sizeof(PASTIX_FLOAT),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_B),
                               ldb*n*sizeof(PASTIX_FLOAT)));
          CUDA_CALL(cudaMemcpy((void*)d_B, B0,
                               ldb*n*sizeof(PASTIX_FLOAT),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaThreadSynchronize());
          clockInit(&(clk));
          clockStart(&(clk));
          GENERATE_SM_VERSION_NAME(gemdm)('N',
                                          'T',
                                          m, n, k,
                                          cu_alpha,
                                          d_A, lda,
                                          d_D, ldd,
                                          d_A, lda,
                                          cu_beta,
                                          d_B, ldb,
                                          nb_blocks,  d_blocks,
                                          nb_fblocks, d_fblocks,
                                          0 );
	  CUDA_CALL(cudaThreadSynchronize());
          clockStop(&(clk));
          CUDA_CALL(cudaMemcpy((void*)B_res, d_B,
                               ldb*n*sizeof(PASTIX_FLOAT),
                               cudaMemcpyDeviceToHost));
          CUDA_CALL(cudaFree(d_A));
          CUDA_CALL(cudaFree(d_B));
          clockStop(&(clk_wt));
          time_sparse_GPU_FERMI_GEMDM = clockVal(&(clk));
          SAVE_MIN_TIME(time_sparse_GPU_FERMI_GEMDM);
          time_sparse_GPU_FERMI_GEMDM_wt = clockVal(&(clk_wt));
          SAVE_MIN_TIME(time_sparse_GPU_FERMI_GEMDM_wt);
          COMPARE_TIME("sparse GEMDM on GPU (Fermi)",
                       time_sparse_GPU_FERMI_GEMDM, ops, time_sparse_CPU);
          COMPARE_TIME("sparse GEMDM on GPU (Fermi) with transfert",
                       time_sparse_GPU_FERMI_GEMDM_wt, ops, time_sparse_CPU);
          COMPARE_RES(B_ref,B_res);
#endif




          CUDA_CALL(cudaFree(d_blocks));
          CUDA_CALL(cudaFree(d_fblocks));
          /* Just for timing, perform dense GEMM */

          clockInit(&(clk));
          clockStart(&(clk));
          SOPALIN_GEMM( "N", "C",
                        m, n, k,
                        alpha,
                        A1, lda,
                        A1, lda,
                        beta,
                        B_ref, ldb);
          clockStop(&(clk));
          time_dense_CPU = clockVal(&(clk));
          SAVE_MIN_TIME(time_dense_CPU);
          PRINT_TIME("dense GEMM on CPU", time_dense_CPU, ops);

          clockInit(&(clk_wt));
          clockStart(&(clk_wt));
          CUDA_CALL(cudaMalloc((void*)&(d_A),
                               lda*k*sizeof(PASTIX_FLOAT)));
          CUDA_CALL(cudaMemcpy((void*)d_A, A1,
                               lda*k*sizeof(PASTIX_FLOAT),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_B),
                               ldb*n*sizeof(PASTIX_FLOAT)));
          CUDA_CALL(cudaMemcpy((void*)d_B, B0,
                               ldb*n*sizeof(PASTIX_FLOAT),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaThreadSynchronize());

          clockInit(&(clk));
          clockStart(&(clk));
          CUDA_CALL(cudaThreadSynchronize());
          CUBLAS_GEMM('N', 'C', m, n, k, cu_alpha,
                      (CU_FLOAT*)d_A, lda, (CU_FLOAT*)d_A, lda,
                      cu_beta, (CU_FLOAT*)d_B, ldb);
          CUDA_CALL(cudaThreadSynchronize());
          clockStop(&(clk));


          CUDA_CALL(cudaMemcpy((void*)B_res, d_B,
                               ldb*n*sizeof(PASTIX_FLOAT),
                               cudaMemcpyDeviceToHost));
          CUDA_CALL(cudaFree(d_A));
          CUDA_CALL(cudaFree(d_B));
          clockStop(&(clk_wt));

          time_dense_GPU = clockVal(&(clk));
          SAVE_MIN_TIME(time_dense_GPU);
          time_dense_GPU_wt = clockVal(&(clk_wt));
          SAVE_MIN_TIME(time_dense_GPU_wt);
#ifdef WITH_MAGMABLAS
          COMPARE_TIME("dense magGEMM on GPU", time_dense_GPU, ops, time_dense_CPU);
          COMPARE_TIME("dense magGEMM on GPU with transfert", time_dense_GPU_wt, ops, time_dense_CPU);
#else
          COMPARE_TIME("dense cuGEMM on GPU", time_dense_GPU, ops, time_dense_CPU);
          COMPARE_TIME("dense cuGEMM on GPU with transfert", time_dense_GPU_wt, ops, time_dense_CPU);
#endif


        }
      PRINT_TIME("(min) sparse GEMM on CPU", min_time_sparse_CPU[run2_idx], ops);
      COMPARE_TIME("(min) sparse GEMM in GPU", min_time_sparse_GPU[run2_idx], ops, min_time_sparse_CPU[run2_idx]);
      COMPARE_TIME("(min) sparse GEMM in GPU with transfert", min_time_sparse_GPU_wt[run2_idx], ops, min_time_sparse_CPU[run2_idx]);
#if (CUDA_SM_VERSION >= 20)
      COMPARE_TIME("(min) sparse GEMM in GPU (FERMI)", min_time_sparse_GPU_FERMI[run2_idx], ops, min_time_sparse_CPU[run2_idx]);
      COMPARE_TIME("(min) sparse GEMM in GPU (FERMI) with transfert", min_time_sparse_GPU_FERMI_wt[run2_idx], ops, min_time_sparse_CPU[run2_idx]);
      COMPARE_TIME("(min) sparse GEMDM in GPU (FERMI)", min_time_sparse_GPU_FERMI_GEMDM[run2_idx], ops, min_time_sparse_CPU[run2_idx]);
      COMPARE_TIME("(min) sparse GEMDM in GPU (FERMI) with transfert", min_time_sparse_GPU_FERMI_GEMDM_wt[run2_idx], ops, min_time_sparse_CPU[run2_idx]);
#endif
      PRINT_TIME("(min) dense GEMM on CPU", min_time_dense_CPU[run2_idx], ops);
      COMPARE_TIME("(min) dense GEMM on GPU", min_time_dense_GPU[run2_idx], ops, min_time_dense_CPU[run2_idx]);
      COMPARE_TIME("(min) dense GEMM on GPU with transfert", min_time_dense_GPU_wt[run2_idx], ops, min_time_dense_CPU[run2_idx]);


      memFree_null(A1);
      memFree_null(B0);
      memFree_null(B_res);
      memFree_null(B_ref);
      memFree_null(blocks);
      memFree_null(fblocks);
    }
  memFree_null(min_time_sparse_CPU);
  memFree_null(min_time_sparse_GPU);
  memFree_null(min_time_sparse_GPU_wt);
#if (CUDA_SM_VERSION >= 20)
  memFree_null(min_time_sparse_GPU_FERMI);
  memFree_null(min_time_sparse_GPU_FERMI_wt);
#endif
  memFree_null(min_time_dense_CPU);
  memFree_null(min_time_dense_GPU);
  memFree_null(min_time_dense_GPU_wt);
  
  return EXIT_SUCCESS;
}
