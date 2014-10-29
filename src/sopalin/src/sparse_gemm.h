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
 * File sparse_gemm.h
 *
 * Header for sparse GEMM CUDA kernel.
 *
 * Author:
 *  Xavier Lacoste - xavier.lacoste@inria.fr
 */
#ifndef SPARSE_GEMM_H
#define SPARSE_GEMM_H

#define sparse_gemm_kernel_N_T_64_16_4_16_4             \
  PASTIX_PREFIX_F(sparse_gemm_kernel_N_T_64_16_4_16_4)
#define magmablas_sparse_gemm_kernel_N_T_64_16_4_16_4                   \
  PASTIX_PREFIX_F(magmablas_sparse_gemm_kernel_N_T_64_16_4_16_4)
#define sparse_gemm_cpu PASTIX_PREFIX_F(sparse_gemm_cpu)

#ifdef __cplusplus
extern "C"
__global__
#endif
/*
 * Function: sparse_gemm_kernel_N_T_64_16_4_16_4
 *
 * CUDA kernel to update C.
 *
 * Performs : $C \leftarrow \alpha A \times B^T + \beta B$
 *
 * Parameters:
 *   m          - Number of rows in *C*.
 *   n          - Number of columns in *C*.
 *   k          - Number of columns in *A*.
 *   alpha      - A coefficient .
 *   A          - $m \times k$ matrix.
 *   lda        - Leading dimension of *A*.
 *   B          - $n \times k$ matrix.
 *   ldb        - Leading dimension of *B*.
 *   beta       - A coefficient.
 *   C          - $m \times n$ matrix.
 *   ldc        - Leading dimension of *C*.
 *   blocknbr   - Number of blocks in *A*.
 *   blocktab   - Array containing first and last row of each block in *A*.
 *                blocktab[2i]   : first row of block i,
 *                blocktab[2i+1] : last row of block i.
 *   fblocknbr  - number of blocks in *C*.
 *   fblocktab  - Array containing first and last row of each block in *C*.
 */
void
sparse_gemm_kernel_N_T_64_16_4_16_4(int m, int n, int k,
                                    PASTIX_FLOAT alpha,
                                    const PASTIX_FLOAT *A, int lda,
                                    const PASTIX_FLOAT *B, int ldb,
                                    PASTIX_FLOAT beta,
                                    PASTIX_FLOAT       *c, int ldc,
                                    int blocknbr,  const int * blocktab,
                                    int fblocknbr, const int * fblocktab);
#ifdef __cplusplus
extern "C"
#endif

/*
 * Function: magmablas_sparse_gemm_kernel_N_T_64_16_4_16_4
 *
 * Interface to the CUDA kernel <sparse_gemm_kernel_N_T_64_16_4_16_4>.
 *
 * Parameters:
 *   m          - Number of rows in *C*.
 *   n          - Number of columns in *C*.
 *   k          - Number of columns in *A*.
 *   alpha      - A coefficient .
 *   A          - $m \times k$ matrix.
 *   lda        - Leading dimension of *A*.
 *   B          - $n \times k$ matrix.
 *   ldb        - Leading dimension of *B*.
 *   beta       - A coefficient.
 *   C          - $m \times n$ matrix.
 *   ldc        - Leading dimension of *C*.
 *   blocknbr   - Number of blocks in *A*.
 *   blocktab   - Array containing first and last row of each block in *A*.
 *                blocktab[2i]   : first row of block i,
 *                blocktab[2i+1] : last row of block i.
 *   fblocknbr  - number of blocks in *C*.
 *   fblocktab  - Array containing first and last row of each block in *C*.
 */
void
magmablas_sparse_gemm_kernel_N_T_64_16_4_16_4(int m, int n, int k,
                                              PASTIX_FLOAT alpha,
                                              const PASTIX_FLOAT *A, int lda,
                                              const PASTIX_FLOAT *B, int ldb,
                                              PASTIX_FLOAT beta,
                                              PASTIX_FLOAT       *c, int ldc,
                                              int blocknbr,  const int * blocktab,
                                              int fblocknbr, const int * fblocktab);


/*
 * Function: sparse_gemm_cpu
 *
 * Performs : $C \leftarrow \alpha A \times B^T + \beta B$
 *
 * Parameters:
 *   m          - Number of rows in *C*.
 *   n          - Number of columns in *C*.
 *   k          - Number of columns in *A*.
 *   alpha      - A coefficient .
 *   A          - $m \times k$ matrix.
 *   lda        - Leading dimension of *A*.
 *   B          - $n \times k$ matrix.
 *   ldb        - Leading dimension of *B*.
 *   beta       - A coefficient.
 *   C          - $m \times n$ matrix.
 *   ldc        - Leading dimension of *C*.
 *   blocknbr   - Number of blocks in *A*.
 *   blocktab   - Array containing first and last row of each block in *A*.
 *                blocktab[2i]   : first row of block i,
 *                blocktab[2i+1] : last row of block i.
 *   fblocknbr  - number of blocks in *C*.
 *   fblocktab  - Array containing first and last row of each block in *C*.
 *   work       - Array of PASTIX_FLOAT used for temporary storage.
 *   worksize   - Size of work.
 */
int
sparse_gemm_cpu( char * transa, char * transb,
                 int m, int n, int k,
                 PASTIX_FLOAT alpha,
                 const PASTIX_FLOAT * a, int lda,
                      const PASTIX_FLOAT * b, int ldb,
                 PASTIX_FLOAT beta,
                 PASTIX_FLOAT * c, unsigned int ldc,
                 int blocknbr,  const int * blocktab,
                 int fblocknbr, const int * fblocktab,
                 PASTIX_FLOAT *work, int worksize);

#ifdef WITH_STARPU
#  ifndef FORCE_NO_CUDA
#    if (CUDA_SM_VERSION >= 20)
#      include "sparse_gemm_fermi.h"
#    endif
#  endif
#endif

#endif /* SPARSE_GEMM_H */
