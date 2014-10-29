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
#ifndef SPARSE_GEMM_FERMI_H
#define SPARSE_GEMM_FERMI_H
#ifdef TYPE_COMPLEX
#  ifdef PREC_DOUBLE
#    define PRECISION_z
#  else
#    define PRECISION_c
#  endif
#else
#  ifdef PREC_DOUBLE
#    define PRECISION_d
#  else
#    define PRECISION_s
#  endif
#endif

#include <cuda.h>
#if defined(PRECISION_z) || defined(PRECISION_c)
#include <cuComplex.h>
#endif
#include <cuda_runtime.h>

#define magmablas_zgemm_fermi magmablas_zgemm
#define magmablas_zgemdm_fermi magmablas_zgemdm

////////////////////////////////////////////////////////////////////////////////
#ifdef PRECISION_z
#  define GENERATE_SM_VERSION_NAME_I(func, version) sparse_z##func##_SM##version
#endif
#ifdef PRECISION_c
#  define GENERATE_SM_VERSION_NAME_I(func, version) sparse_c##func##_SM##version
#  define cublasZgemm cublasCgemm
#  define cuDoubleComplex cuFloatComplex
#endif
#ifdef PRECISION_d
#  define GENERATE_SM_VERSION_NAME_I(func, version) sparse_d##func##_SM##version
#  define cublasZgemm cublasDgemm
#  define cuDoubleComplex double
#endif
#ifdef PRECISION_s
#  define GENERATE_SM_VERSION_NAME_I(func, version) sparse_s##func##_SM##version
#  define cublasZgemm cublasSgemm
#  define cuDoubleComplex float
#endif
#define GENERATE_SM_VERSION_KERNEL_NAME_I2(func, version) GENERATE_SM_VERSION_NAME_I(func, version)
#define GENERATE_SM_VERSION_KERNEL_NAME(func) GENERATE_SM_VERSION_NAME_I2(func, CUDA_SM_VERSION)


#define GENERATE_SM_VERSION_NAME_I2(func, version) GENERATE_SM_VERSION_NAME_I(func, version)
#define GENERATE_SM_VERSION_NAME(func) GENERATE_SM_VERSION_NAME_I2(func, CUDA_SM_VERSION)

#ifdef __cplusplus
extern "C"
#endif
void
GENERATE_SM_VERSION_NAME(gemm)( char TRANSA, char TRANSB,
                                int m , int n , int k ,
                                cuDoubleComplex alpha,
                                const cuDoubleComplex *d_A, int lda,
                                const cuDoubleComplex *d_B, int ldb,
                                cuDoubleComplex beta,
                                cuDoubleComplex *d_C, int ldc,
                                int blocknbr, const int *blocktab,
                                int fblocknbr, const int *fblocktab,
                                cudaStream_t stream );
#ifdef __cplusplus
extern "C"
#endif
void
GENERATE_SM_VERSION_NAME(gemdm)( char TRANSA, char TRANSB,
                                 int m , int n , int k ,
                                 cuDoubleComplex alpha,
                                 const cuDoubleComplex *d_A, int lda,
                                 const cuDoubleComplex *d_D, int ldd,
                                 const cuDoubleComplex *d_B, int ldb,
                                 cuDoubleComplex beta,
                                 cuDoubleComplex *d_C, int ldc,
                                 int blocknbr, const int *blocktab,
                                 int fblocknbr, const int *fblocktab,
                                 cudaStream_t stream );

#endif
