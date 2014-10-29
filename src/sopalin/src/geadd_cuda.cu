#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#include "redefine_functions.h"

#ifdef TYPE_COMPLEX
#  include <cuComplex.h>
#  ifdef PREC_DOUBLE
#    define PASTIX_FLOAT cuDoubleComplex
#    define MULT(a, b)  cuCmul(a,b)
#    define ADD(a, b)   cuCadd(a,b)
#  else
#    define PASTIX_FLOAT cuFloatComplex
#    define MULT(a,b)  cuCmulf(a,b)
#    define ADD(a, b)  cuCaddf(a,b)
#  endif
#else /* not TYPE_COMPLEX */
#  define MULT(a,b) ((a) * (b))
#  define ADD(a, b) ((a) + (b))
#  ifdef PREC_DOUBLE
#    define PASTIX_FLOAT double
#  else
#    define PASTIX_FLOAT float
#  endif
#endif

#include "geadd_cuda.h"

#define geadd_cuda_kernel_nn PASTIX_PREFIX_F(geadd_cuda_kernel_nn)
#define geadd_cuda_kernel_tn PASTIX_PREFIX_F(geadd_cuda_kernel_tn)
#define geadd_cuda_kernel_nt PASTIX_PREFIX_F(geadd_cuda_kernel_nt)
#define geadd_cuda_kernel_tt PASTIX_PREFIX_F(geadd_cuda_kernel_tt)

__global__ void geadd_cuda_kernel_nn(int m, int n,
				     PASTIX_FLOAT alpha, PASTIX_FLOAT *a, int lda,
				     PASTIX_FLOAT beta,  PASTIX_FLOAT *b, int ldb)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  int j=blockIdx.y*blockDim.y+threadIdx.y;
  int index_a =i+j*lda;
  int index_b =i+j*ldb;
  if ( i < m && j < n )
    b[index_b]= ADD( MULT(alpha,a[index_a]), MULT(beta,b[index_b]));
}

__global__ void geadd_cuda_kernel_nt(int m, int n,
				     PASTIX_FLOAT alpha, PASTIX_FLOAT *a, int lda,
				     PASTIX_FLOAT beta,  PASTIX_FLOAT *b, int ldb)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  int j=blockIdx.y*blockDim.y+threadIdx.y;
  int index_a =i+j*lda;
  int index_b =j+i*ldb;
  if ( i < m && j < n)
    b[index_b]= ADD( MULT(alpha,a[index_a]), MULT(beta,b[index_b]));
}



__global__ void geadd_cuda_kernel_tn(int m, int n,
				     PASTIX_FLOAT alpha, PASTIX_FLOAT *a, int lda,
				     PASTIX_FLOAT beta,  PASTIX_FLOAT *b, int ldb)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  int j=blockIdx.y*blockDim.y+threadIdx.y;
  int index_a =j+i*lda;
  int index_b =i+j*ldb;
  if ( i < m && j < n)
    b[index_b]= ADD( MULT(alpha,a[index_a]), MULT(beta,b[index_b]));
}

__global__ void geadd_cuda_kernel_tt(int m, int n,
				     PASTIX_FLOAT alpha, PASTIX_FLOAT *a, int lda,
				     PASTIX_FLOAT beta,  PASTIX_FLOAT *b, int ldb)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  int j=blockIdx.y*blockDim.y+threadIdx.y;
  int index_a =j+i*lda;
  int index_b =j+i*ldb;
  if ( i < m && j < n)
    b[index_b]= ADD( MULT(alpha,a[index_a]), MULT(beta,b[index_b]));
}
extern "C" void
geadd_cuda(char * transa, char * transb,
	   int m, int n,
	   PASTIX_FLOAT alpha, PASTIX_FLOAT *a, int lda,
	   PASTIX_FLOAT beta,  PASTIX_FLOAT *b, int ldb)
{
  dim3 threads( 16, 4 );
  dim3 grid (m/threads.x + (m%threads.x != 0), n/threads.y + (n%threads.y != 0));

  if (*transa == 'N') {
    if (*transb == 'N') {
      geadd_cuda_kernel_nn<<< grid, threads >>>(m, n,
						alpha, a, lda,
						beta,  b, ldb);
    }
    else {
      geadd_cuda_kernel_nt<<< grid, threads >>>(m, n,
						alpha, a, lda,
						beta,  b, ldb);
    }
  }
  else {
    if (*transb == 'N') {
      geadd_cuda_kernel_tn<<< grid, threads >>>(m, n,
						alpha, a, lda,
						beta,  b, ldb);
    }
    else {
      geadd_cuda_kernel_tt<<< grid, threads >>>(m, n,
						alpha, a, lda,
						beta,  b, ldb);
    }  }
}
