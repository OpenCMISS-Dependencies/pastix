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
#ifndef GETRA_CUDA_H
#define GETRA_CUDA_H
#define getra_cuda_kernel PASTIX_PREFIX_F(getra_cuda_kernel)
#define getra_cuda        PASTIX_PREFIX_F(getra_cuda)
/*
 * Function: getra_cuda_kernel
 *
 * kernel to compute a = a ^T
 *
 * Parameters:
 *   a     - Matrix to transpose.
 *   lda   - Leading dimension of a.
 *   N     - Size of the matrix.
 */
void getra_cuda_kernel(PASTIX_FLOAT *A, int lda, int N);



/*
 * Function: getra_cuda
 *
 * Interface to the kernel to compute b = a ^T
 *
 * Parameters:
 *   A     - Matrix to transpose.
 *   lda   - Leading dimension of a.
 *   B     - Matrix to receive the transposed matrix.
 *   ldb   - Leading dimension of a.
 *   N     - Size of the matrix.
 */
#ifdef __cplusplus
extern "C"
#endif
void
getra_cuda(PASTIX_FLOAT *A, int lda, PASTIX_FLOAT *B, int ldb, int N);
#endif /* GETRA_CUDA_H */
