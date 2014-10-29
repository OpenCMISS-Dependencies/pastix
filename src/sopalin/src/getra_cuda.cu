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

#define TILE_DIM   16

#include "getra_cuda.h"

__global__ void getra_cuda_kernel(PASTIX_FLOAT *A, int lda, PASTIX_FLOAT * B, int ldb, int N)
{
  __shared__ PASTIX_FLOAT block[TILE_DIM][TILE_DIM+1];
	
  // read the matrix tile into shared memory
  unsigned int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  unsigned int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  if((xIndex < N) && (yIndex < N))
    {
      unsigned int index_in = yIndex * lda + xIndex;
      block[threadIdx.y][threadIdx.x] = A[index_in];
    }

  __syncthreads();

  // write the transposed matrix tile to global memory
  xIndex = blockIdx.y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx.x * TILE_DIM + threadIdx.y;
  if((xIndex < N) && (yIndex < N))
    {
      unsigned int index_out = yIndex * ldb + xIndex;
      B[index_out] = block[threadIdx.x][threadIdx.y];
    }
}



extern "C" void
getra_cuda(PASTIX_FLOAT *A, int lda, PASTIX_FLOAT * B, int ldb, int N)
{
  dim3 threads( TILE_DIM, TILE_DIM );
  dim3 grid (N/threads.x + (N%threads.x != 0), 1+N/threads.y + (N%threads.y != 0));

  getra_cuda_kernel<<< grid, threads >>>(A, lda, B, ldb, N);
}
