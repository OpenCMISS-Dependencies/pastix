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
#ifndef STARPU_KERNELS_H
#define STARPU_KERNELS_H

#define ARCH_CPU  0
#define ARCH_CUDA 1

#define MAT_zaxpy( m, n, alpha,                   \
                   A, lda,                        \
                   B, ldb )                       \
  {                                               \
    PASTIX_INT i,j;                                      \
    for (j=0; j<n; j++) {                         \
      for (i=0; i<m; i++) {                       \
        B[j*ldb+i] = B[j*ldb+i] - A[j*lda+i];     \
      }                                           \
    }                                             \
  }

#define MAT_zaxpyt( m, n, alpha,                  \
                    A, lda,                       \
                    B, ldb )                      \
  {                                               \
    PASTIX_INT i,j;                                      \
    for (j=0; j<n; j++) {                         \
      for (i=0; i<m; i++) {                       \
        B[j*ldb+i] = B[j*ldb+i] - A[j+i*lda];     \
      }                                           \
    }                                             \
}


#define getrfsp1d_starpu_common             PASTIX_PREFIX_F(getrfsp1d_starpu_common)
#define getrfsp1d_gemm_starpu_common        PASTIX_PREFIX_F(getrfsp1d_gemm_starpu_common)
#define getrfsp1d_sparse_gemm_starpu_common PASTIX_PREFIX_F(getrfsp1d_sparse_gemm_starpu_common)
#define getrfsp1d_starpu_cpu                PASTIX_PREFIX_F(getrfsp1d_starpu_cpu)
#define getrfsp1d_gemm_starpu_cpu           PASTIX_PREFIX_F(getrfsp1d_gemm_starpu_cpu)
#define getrfsp1d_sparse_gemm_starpu_cpu    PASTIX_PREFIX_F(getrfsp1d_sparse_gemm_starpu_cpu)
#define getrfsp1d_starpu_cuda               PASTIX_PREFIX_F(getrfsp1d_starpu_cuda)
#define getrfsp1d_gemm_starpu_cuda          PASTIX_PREFIX_F(getrfsp1d_gemm_starpu_cuda)
#define getrfsp1d_sparse_gemm_starpu_cuda   PASTIX_PREFIX_F(getrfsp1d_sparse_gemm_starpu_cuda)
void getrfsp1d_starpu_cpu(void * buffers[], void * _args);
void getrfsp1d_gemm_starpu_cpu(void * buffers[], void * _args);
void getrfsp1d_sparse_gemm_starpu_cpu(void * buffers[], void * _args);
void getrfsp1d_starpu_cuda(void * buffers[], void * _args);
void getrfsp1d_gemm_starpu_cuda(void * buffers[], void * _args);
void getrfsp1d_sparse_gemm_starpu_cuda(void * buffers[], void * _args);

#define potrfsp1d_starpu_common             PASTIX_PREFIX_F(potrfsp1d_starpu_common)
#define potrfsp1d_gemm_starpu_common        PASTIX_PREFIX_F(potrfsp1d_gemm_starpu_common)
#define potrfsp1d_sparse_gemm_starpu_common PASTIX_PREFIX_F(potrfsp1d_sparse_gemm_starpu_common)
#define potrfsp1d_starpu_cpu                PASTIX_PREFIX_F(potrfsp1d_starpu_cpu)
#define potrfsp1d_gemm_starpu_cpu           PASTIX_PREFIX_F(potrfsp1d_gemm_starpu_cpu)
#define potrfsp1d_sparse_gemm_starpu_cpu    PASTIX_PREFIX_F(potrfsp1d_sparse_gemm_starpu_cpu)
#define potrfsp1d_starpu_cuda               PASTIX_PREFIX_F(potrfsp1d_starpu_cuda)
#define potrfsp1d_gemm_starpu_cuda          PASTIX_PREFIX_F(potrfsp1d_gemm_starpu_cuda)
#define potrfsp1d_sparse_gemm_starpu_cuda   PASTIX_PREFIX_F(potrfsp1d_sparse_gemm_starpu_cuda)
void potrfsp1d_starpu_cpu(void * buffers[], void * _args);
void potrfsp1d_gemm_starpu_cpu(void * buffers[], void * _args);
void potrfsp1d_sparse_gemm_starpu_cpu(void * buffers[], void * _args);
void potrfsp1d_starpu_cuda(void * buffers[], void * _args);
void potrfsp1d_gemm_starpu_cuda(void * buffers[], void * _args);
void potrfsp1d_sparse_gemm_starpu_cuda(void * buffers[], void * _args);

#ifdef HERMITIAN
#  define hetrfsp1d_starpu_common      PASTIX_PREFIX_F(hetrfsp1d_starpu_common)
#  define hetrfsp1d_gemm_starpu_common PASTIX_PREFIX_F(hetrfsp1d_gemm_starpu_common)
#  define hetrfsp1d_starpu_cpu         PASTIX_PREFIX_F(hetrfsp1d_starpu_cpu)
#  define hetrfsp1d_gemm_starpu_cpu    PASTIX_PREFIX_F(hetrfsp1d_gemm_starpu_cpu)
#  define hetrfsp1d_starpu_cuda        PASTIX_PREFIX_F(hetrfsp1d_starpu_cuda)
#  define hetrfsp1d_gemm_starpu_cuda   PASTIX_PREFIX_F(hetrfsp1d_gemm_starpu_cuda)
#else
#  define hetrfsp1d_starpu_common      PASTIX_PREFIX_F(sytrfsp1d_starpu_common)
#  define hetrfsp1d_gemm_starpu_common PASTIX_PREFIX_F(sytrfsp1d_gemm_starpu_common)
#  define hetrfsp1d_starpu_cpu         PASTIX_PREFIX_F(sytrfsp1d_starpu_cpu)
#  define hetrfsp1d_gemm_starpu_cpu    PASTIX_PREFIX_F(sytrfsp1d_gemm_starpu_cpu)
#  define hetrfsp1d_starpu_cuda        PASTIX_PREFIX_F(sytrfsp1d_starpu_cuda)
#  define hetrfsp1d_gemm_starpu_cuda   PASTIX_PREFIX_F(sytrfsp1d_gemm_starpu_cuda)
#endif
void hetrfsp1d_starpu_cpu(void * buffers[], void * _args);
void hetrfsp1d_gemm_starpu_cpu(void * buffers[], void * _args);
void hetrfsp1d_starpu_cuda(void * buffers[], void * _args);
void hetrfsp1d_gemm_starpu_cuda(void * buffers[], void * _args);

#endif /* STARPU_KERNELS_H */
