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
#ifndef GEADD_CUDA_H
#define GEADD_CUDA_H
#define geadd_cuda_kernel PASTIX_PREFIX_F(geadd_cuda_kernel)
#define geadd_cuda        PASTIX_PREFIX_F(geadd_cuda)
/*
 * Function: geadd_cuda_kernel_nn
 *
 * kernel to compute b = alpha*a + beta*b.
 *
 * Parameters:
 *   m     - Number of rows in the matrices.
 *   n     - Number of columns in the matrices.
 *   alpha - Coefficient to multiply a.
 *   a     - Matrice to add.
 *   lda   - Leading dimension of a.
 *   beta  - Coefficient to multiply b.
 *   b     - Matrice to receive addition.
 *   ldb   - Leading dimension of b.
 */
#ifdef __cplusplus
extern "C"
__global__
#endif
void geadd_cuda_kernel_nn(int m, int n,
                          PASTIX_FLOAT alpha, PASTIX_FLOAT *a, int lda,
                          PASTIX_FLOAT beta,  PASTIX_FLOAT *b, int ldb);

/*
 * Function: geadd_cuda_kernel_nn
 *
 * kernel to compute b^T = alpha*a + beta*b^T.
 *
 * Parameters:
 *   m     - Number of rows in the matrices.
 *   n     - Number of columns in the matrices.
 *   alpha - Coefficient to multiply a.
 *   a     - Matrice to add.
 *   lda   - Leading dimension of a.
 *   beta  - Coefficient to multiply b.
 *   b     - Matrice to receive addition.
 *   ldb   - Leading dimension of b.
 */
#ifdef __cplusplus
extern "C"
__global__
#endif
void geadd_cuda_kernel_nt(int m, int n,
                          PASTIX_FLOAT alpha, PASTIX_FLOAT *a, int lda,
                          PASTIX_FLOAT beta,  PASTIX_FLOAT *b, int ldb);

/*
 * Function: geadd_cuda_kernel_nn
 *
 * kernel to compute b = alpha*a^T + beta*b.
 *
 * Parameters:
 *   m     - Number of rows in the matrices.
 *   n     - Number of columns in the matrices.
 *   alpha - Coefficient to multiply a.
 *   a     - Matrice to add.
 *   lda   - Leading dimension of a.
 *   beta  - Coefficient to multiply b.
 *   b     - Matrice to receive addition.
 *   ldb   - Leading dimension of b.
 */
#ifdef __cplusplus
extern "C"
__global__
#endif
void geadd_cuda_kernel_tn(int m, int n,
                          PASTIX_FLOAT alpha, PASTIX_FLOAT *a, int lda,
                          PASTIX_FLOAT beta,  PASTIX_FLOAT *b, int ldb);

/*
 * Function: geadd_cuda_kernel_nn
 *
 * kernel to compute b^T = alpha*a^T + beta*b^T.
 *
 * Parameters:
 *   m     - Number of rows in the matrices.
 *   n     - Number of columns in the matrices.
 *   alpha - Coefficient to multiply a.
 *   a     - Matrice to add.
 *   lda   - Leading dimension of a.
 *   beta  - Coefficient to multiply b.
 *   b     - Matrice to receive addition.
 *   ldb   - Leading dimension of b.
 */
#ifdef __cplusplus
extern "C"
__global__
#endif
void geadd_cuda_kernel_tt(int m, int n,
                          PASTIX_FLOAT alpha, PASTIX_FLOAT *a, int lda,
                          PASTIX_FLOAT beta,  PASTIX_FLOAT *b, int ldb);

/*
 * Function: geadd_cuda
 *
 * Interface to the kernel to compute b = alpha*a + beta*b.
 *
 * Parameters:
 *   transa - "T" if op(a) is a^T.
 *   transb - "T" if op(b) is b^T.
 *   m     - Number of rows in the matrices.
 *   n     - Number of columns in the matrices.
 *   alpha - Coefficient to multiply a.
 *   a     - Matrice to add.
 *   lda   - Leading dimension of a.
 *   beta  - Coefficient to multiply b.
 *   b     - Matrice to receive addition.
 *   ldb   - Leading dimension of b.
 */
#ifdef __cplusplus
extern "C"
#endif
void geadd_cuda(char * transa, char * transb,
                int m, int n,
                PASTIX_FLOAT alpha, PASTIX_FLOAT *a, int lda,
                PASTIX_FLOAT beta,  PASTIX_FLOAT *b, int ldb);

#endif /* GEADD_CUDA_H */
