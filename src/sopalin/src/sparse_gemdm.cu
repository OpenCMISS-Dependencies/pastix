/*
 * File: sparse_gemdm.cu
 *
 * This file has been built from :
 * > -- MAGMA (version 1.1) --
 * >    Univ. of Tennessee, Knoxville
 * >    Univ. of California, Berkeley
 * >    Univ. of Colorado, Denver
 * >    November 2011
 *
 * It hads multi-arithmetic and works on sparse *A* and *C* matrices.
 *
 */


#ifdef STARPU_USE_DEPRECATED_API
#  undef STARPU_USE_DEPRECATED_API
#endif
#include <starpu.h>
#include <starpu_cuda.h>

#include "redefine_functions.h"

#define axpy PASTIX_PREFIX_F(axpy)
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

#include "sparse_gemdm.h"

/*
 * Function: axpy
 *
 * Perform : $C \leftarrow C + a \times B$
 *
 * Parameters:
 *   a - A coefficient.
 *   c - Vector of 16 <PASTIX_FLOAT>.
 *   b - Vector of 16 <PASTIX_FLOAT>.
 */
static __device__ void
axpy(PASTIX_FLOAT a, PASTIX_FLOAT *b, PASTIX_FLOAT *c) {
  c[0]  = ADD(c[0],  MULT(a, b[0]));
  c[1]  = ADD(c[1],  MULT(a, b[1]));
  c[2]  = ADD(c[2],  MULT(a, b[2]));
  c[3]  = ADD(c[3],  MULT(a, b[3]));
  c[4]  = ADD(c[4],  MULT(a, b[4]));
  c[5]  = ADD(c[5],  MULT(a, b[5]));
  c[6]  = ADD(c[6],  MULT(a, b[6]));
  c[7]  = ADD(c[7],  MULT(a, b[7]));
  c[8]  = ADD(c[8],  MULT(a, b[8]));
  c[9]  = ADD(c[9],  MULT(a, b[9]));
  c[10] = ADD(c[10], MULT(a, b[10]));
  c[11] = ADD(c[11], MULT(a, b[11]));
  c[12] = ADD(c[12], MULT(a, b[12]));
  c[13] = ADD(c[13], MULT(a, b[13]));
  c[14] = ADD(c[14], MULT(a, b[14]));
  c[15] = ADD(c[15], MULT(a, b[15]));
}

/*
 * Function: sparse_gemdm_kernel_N_T_64_16_4_16_4
 *
 * CUDA kernel to update C.
 *
 * Performs : $C \leftarrow \alpha A \times D \times B^T + \beta B$
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
extern "C" __global__ void
sparse_gemdm_kernel_N_T_64_16_4_16_4(int m, int n, int k,
				     PASTIX_FLOAT alpha,
				     const PASTIX_FLOAT *A, int lda,
				     const PASTIX_FLOAT *D, int ldd,
				     const PASTIX_FLOAT *B, int ldb,
				     PASTIX_FLOAT beta,
				     PASTIX_FLOAT       *C, int ldc,
				     int blocknbr,
				     const int * blocktab,
				     int fblocknbr,
				     const int * fblocktab)
{
  /*  -- MAGMA (version 1.1) --
      Univ. of Tennessee, Knoxville
      Univ. of California, Berkeley
      Univ. of Colorado, Denver
      November 2011

      Purpose:
      ========
      This routine computes
      C = alpha* A*B^T  + beta * C

      B is put into shared memory
      Parameters Used:
      blk_M=64 blk_N=16 blk_K=4 nthd_x=16 nthd_y=4

      This code should run for any matrix size.
      ===============================================================  */

  const int tx = threadIdx.x;
  const int ty = threadIdx.y;

  const int ibx = blockIdx.x * 64;
  const int iby = blockIdx.y *16;


  const int idt = ty * 16 + tx;

  if( iby + tx >=n )
    B+= iby+0;
  else
    B+= iby+tx;
  /*
    Taking care of boundary cases where K<4.
  */
  if( ty >=k )
    B+= __mul24( 0,ldb);
  else
    B+= __mul24( ty,ldb);

  if( ibx + idt >= m )
    A += ibx + 0 ;
  else
    A += ibx + idt;


  int s2=lda, s3=2*lda, s4=3*lda ;
  int t2=ldd, t3=2*ldd, t4=3*ldd ;

  switch (k){
  case 1:
    s2=0; s3=0; s4=0;
    t2=0; t3=0; t4=0;
    break ;
  case 2:
    s2=lda; s3=0; s4=0;
    t2=ldd; t3=0; t4=0;
    break ;
  case 3:
    s2=lda; s3=2*lda; s4=0;
    t2=ldd; t3=2*ldd; t4=0;
    break ;
  }

  {
#define FROWNUM(tab, b) tab[2*b]
#define LROWNUM(tab, b) tab[2*b+1]
#define BLOCKSIZE(tab, b) LROWNUM(tab, b) - FROWNUM(tab, b) + 1
    int idx_x = ibx +idt;
    int blocknum = 0, fblocknum = 0;
    size_t totalblocksize = 0;
    size_t blocksize = BLOCKSIZE(blocktab, blocknum);
    int rownum;
    int offset;
    while(totalblocksize + blocksize < idx_x + 1)
      {
        totalblocksize += blocksize;
        blocknum++;
        blocksize = BLOCKSIZE(blocktab, blocknum);
      }
    rownum = idx_x - totalblocksize + FROWNUM(blocktab, blocknum);
    offset = 0;
    while (LROWNUM(fblocktab, fblocknum) < rownum) {
      offset += BLOCKSIZE(fblocktab, fblocknum);
      fblocknum++;
    }
    offset += rownum - FROWNUM(fblocktab, fblocknum);

    C += offset + __mul24(iby,ldc);
#undef FROWNUM
#undef LROWNUM
  }



  PASTIX_FLOAT Ap[4]={A[0], A[s2], A[s3], A[s4]};
  PASTIX_FLOAT Dp[4]={D[0], D[t2], D[t3], D[t4]};
  PASTIX_FLOAT b=B[0];

  const PASTIX_FLOAT *Bend = B + ldb*(k-k%4);

  B+=4*ldb;
  A+=4*lda;

  __shared__ PASTIX_FLOAT Bb[4][16];

  PASTIX_FLOAT Cb[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  if(k>7)
    do {
      PASTIX_FLOAT Ab[4] = {MULT(Ap[0],Dp[0]), MULT(Ap[1],Dp[1]), MULT(Ap[2],Dp[2]), MULT(Ap[3],Dp[3])};

      Bb[ty][tx]=b;

      __syncthreads();

      Ap[0] = A[0];
      Ap[1] = A[s2];
      Ap[2] = A[s3];
      Ap[3] = A[s4];
      Dp[0] = D[0];
      Dp[1] = D[t2];
      Dp[2] = D[t3];
      Dp[3] = D[t4];

      b=B[0];


      axpy(Ab[0], &Bb[0][0], Cb);
      axpy(Ab[1], &Bb[1][0], Cb);
      axpy(Ab[2], &Bb[2][0], Cb);
      axpy(Ab[3], &Bb[3][0], Cb);

      A += 4*lda;
      D += 4*ldd;
      B += 4*ldb;

      __syncthreads();
    } while (B < Bend);

  if(k>3){

    Bb[ty][tx]=b;
    int k1 = k-k%4;

    if( (k1+ty) >=k)
      B-=4*ldb;
    else
      B-=0*ldb;

    if( (k1+0) >= k ) {s2=0;s3=0*lda;s4=0;A-=4*lda;} else
      if( (k1+1) >= k ) {s2=0;s3=0*lda;s4=0;A-=0*lda;} else
        if( (k1+2) >= k ) {s2=lda;s3=0*lda;s4=0;A-=0*lda;} else
          if( (k1+3) >= k ) {s2=lda;s3=2*lda;s4=0;A-=0*lda;}

    __syncthreads();


    b=B[0];

    axpy(MULT(Ap[0], Dp[0]), &Bb[0][0], Cb);        Ap[0] = A[0];
    axpy(MULT(Ap[1], Dp[1]), &Bb[1][0], Cb);        Ap[1] = A[s2];
    axpy(MULT(Ap[2], Dp[2]), &Bb[2][0], Cb);        Ap[2] = A[s3];
    axpy(MULT(Ap[3], Dp[3]), &Bb[3][0], Cb);        Ap[3] = A[s4];

  }

  k=k%4;

  if ( k!=0){

    __syncthreads();

    Bb[ty][tx]=b;

    __syncthreads();

    for(int i=0;i<k;i++){
      axpy(MULT(Ap[i], Dp[i]),&Bb[i][0], Cb);
    }
  }



  if( (iby+16)>=n) {
    lda = n-iby;
  }
  else{
    lda = 16;
  }

  if( (ibx+idt) >= m )
    lda = 0 ;
  else lda = lda ;



  switch(lda){
  case 16:

    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    C[6*ldc]  = ADD( MULT(alpha, Cb[6]),   MULT(beta, C[6*ldc]) );
    C[7*ldc]  = ADD( MULT(alpha, Cb[7]),   MULT(beta, C[7*ldc]) );
    C[8*ldc]  = ADD( MULT(alpha, Cb[8]),   MULT(beta, C[8*ldc]) );
    C[9*ldc]  = ADD( MULT(alpha, Cb[9]),   MULT(beta, C[9*ldc]) );
    C[10*ldc] = ADD( MULT(alpha, Cb[10]),  MULT(beta, C[10*ldc]) );
    C[11*ldc] = ADD( MULT(alpha, Cb[11]),  MULT(beta, C[11*ldc]) );
    C[12*ldc] = ADD( MULT(alpha, Cb[12]),  MULT(beta, C[12*ldc]) );
    C[13*ldc] = ADD( MULT(alpha, Cb[13]),  MULT(beta, C[13*ldc]) );
    C[14*ldc] = ADD( MULT(alpha, Cb[14]),  MULT(beta, C[14*ldc]) );
    C[15*ldc] = ADD( MULT(alpha, Cb[15]),  MULT(beta, C[15*ldc]) );

    break;
  case 15:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    C[6*ldc]  = ADD( MULT(alpha, Cb[6]),   MULT(beta, C[6*ldc]) );
    C[7*ldc]  = ADD( MULT(alpha, Cb[7]),   MULT(beta, C[7*ldc]) );
    C[8*ldc]  = ADD( MULT(alpha, Cb[8]),   MULT(beta, C[8*ldc]) );
    C[9*ldc]  = ADD( MULT(alpha, Cb[9]),   MULT(beta, C[9*ldc]) );
    C[10*ldc] = ADD( MULT(alpha, Cb[10]),  MULT(beta, C[10*ldc]) );
    C[11*ldc] = ADD( MULT(alpha, Cb[11]),  MULT(beta, C[11*ldc]) );
    C[12*ldc] = ADD( MULT(alpha, Cb[12]),  MULT(beta, C[12*ldc]) );
    C[13*ldc] = ADD( MULT(alpha, Cb[13]),  MULT(beta, C[13*ldc]) );
    C[14*ldc] = ADD( MULT(alpha, Cb[14]),  MULT(beta, C[14*ldc]) );
    break;
  case 14:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    C[6*ldc]  = ADD( MULT(alpha, Cb[6]),   MULT(beta, C[6*ldc]) );
    C[7*ldc]  = ADD( MULT(alpha, Cb[7]),   MULT(beta, C[7*ldc]) );
    C[8*ldc]  = ADD( MULT(alpha, Cb[8]),   MULT(beta, C[8*ldc]) );
    C[9*ldc]  = ADD( MULT(alpha, Cb[9]),   MULT(beta, C[9*ldc]) );
    C[10*ldc] = ADD( MULT(alpha, Cb[10]),  MULT(beta, C[10*ldc]) );
    C[11*ldc] = ADD( MULT(alpha, Cb[11]),  MULT(beta, C[11*ldc]) );
    C[12*ldc] = ADD( MULT(alpha, Cb[12]),  MULT(beta, C[12*ldc]) );
    C[13*ldc] = ADD( MULT(alpha, Cb[13]),  MULT(beta, C[13*ldc]) );
    break;
  case 13:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    C[6*ldc]  = ADD( MULT(alpha, Cb[6]),   MULT(beta, C[6*ldc]) );
    C[7*ldc]  = ADD( MULT(alpha, Cb[7]),   MULT(beta, C[7*ldc]) );
    C[8*ldc]  = ADD( MULT(alpha, Cb[8]),   MULT(beta, C[8*ldc]) );
    C[9*ldc]  = ADD( MULT(alpha, Cb[9]),   MULT(beta, C[9*ldc]) );
    C[10*ldc] = ADD( MULT(alpha, Cb[10]),  MULT(beta, C[10*ldc]) );
    C[11*ldc] = ADD( MULT(alpha, Cb[11]),  MULT(beta, C[11*ldc]) );
    C[12*ldc] = ADD( MULT(alpha, Cb[12]),  MULT(beta, C[12*ldc]) );
    break;
  case 12:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    C[6*ldc]  = ADD( MULT(alpha, Cb[6]),   MULT(beta, C[6*ldc]) );
    C[7*ldc]  = ADD( MULT(alpha, Cb[7]),   MULT(beta, C[7*ldc]) );
    C[8*ldc]  = ADD( MULT(alpha, Cb[8]),   MULT(beta, C[8*ldc]) );
    C[9*ldc]  = ADD( MULT(alpha, Cb[9]),   MULT(beta, C[9*ldc]) );
    C[10*ldc] = ADD( MULT(alpha, Cb[10]),  MULT(beta, C[10*ldc]) );
    C[11*ldc] = ADD( MULT(alpha, Cb[11]),  MULT(beta, C[11*ldc]) );
    break;
  case 11:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    C[6*ldc]  = ADD( MULT(alpha, Cb[6]),   MULT(beta, C[6*ldc]) );
    C[7*ldc]  = ADD( MULT(alpha, Cb[7]),   MULT(beta, C[7*ldc]) );
    C[8*ldc]  = ADD( MULT(alpha, Cb[8]),   MULT(beta, C[8*ldc]) );
    C[9*ldc]  = ADD( MULT(alpha, Cb[9]),   MULT(beta, C[9*ldc]) );
    C[10*ldc] = ADD( MULT(alpha, Cb[10]),  MULT(beta, C[10*ldc]) );
    break;
  case 10:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    C[6*ldc]  = ADD( MULT(alpha, Cb[6]),   MULT(beta, C[6*ldc]) );
    C[7*ldc]  = ADD( MULT(alpha, Cb[7]),   MULT(beta, C[7*ldc]) );
    C[8*ldc]  = ADD( MULT(alpha, Cb[8]),   MULT(beta, C[8*ldc]) );
    C[9*ldc]  = ADD( MULT(alpha, Cb[9]),   MULT(beta, C[9*ldc]) );
    break;
  case 9:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    C[6*ldc]  = ADD( MULT(alpha, Cb[6]),   MULT(beta, C[6*ldc]) );
    C[7*ldc]  = ADD( MULT(alpha, Cb[7]),   MULT(beta, C[7*ldc]) );
    C[8*ldc]  = ADD( MULT(alpha, Cb[8]),   MULT(beta, C[8*ldc]) );
    break;
  case 8:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    C[6*ldc]  = ADD( MULT(alpha, Cb[6]),   MULT(beta, C[6*ldc]) );
    C[7*ldc]  = ADD( MULT(alpha, Cb[7]),   MULT(beta, C[7*ldc]) );
    break;
  case 7:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    C[6*ldc]  = ADD( MULT(alpha, Cb[6]),   MULT(beta, C[6*ldc]) );
    break;
  case 6:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    C[5*ldc]  = ADD( MULT(alpha, Cb[5]),   MULT(beta, C[5*ldc]) );
    break;
  case 5:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    C[4*ldc]  = ADD( MULT(alpha, Cb[4]),   MULT(beta, C[4*ldc]) );
    break;
  case 4:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    C[3*ldc]  = ADD( MULT(alpha, Cb[3]),   MULT(beta, C[3*ldc]) );
    break;
  case 3:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    C[2*ldc]  = ADD( MULT(alpha, Cb[2]),   MULT(beta, C[2*ldc]) );
    break;
  case 2:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    C[1*ldc]  = ADD( MULT(alpha, Cb[1]),   MULT(beta, C[1*ldc]) );
    break;
  case 1:
    C[0]      = ADD( MULT(alpha, Cb[0]),   MULT(beta, C[0]) );
    break;
  case 0:
    break;
  }

}

/*
 * Function: magmablas_sparse_gemdm_kernel_N_T_64_16_4_16_4
 *
 * Interface to the CUDA kernel <sparse_gemdm_kernel_N_T_64_16_4_16_4>.
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
 *   blocknbr   - Number of blocks in *A*.
 *   blocktab   - Array containing first and last row of each block in *A*.
 *                blocktab[2i]   : first row of block i,
 *                blocktab[2i+1] : last row of block i.
 *   fblocknbr  - number of blocks in *C*.
 *   fblocktab  - Array containing first and last row of each block in *C*.
 */
extern "C" void
magmablas_sparse_gemdm_kernel_N_T_64_16_4_16_4(int m, int n, int k,
					       PASTIX_FLOAT alpha,
					       const PASTIX_FLOAT *A, int lda,
					       const PASTIX_FLOAT *D, int ldd,
					       const PASTIX_FLOAT *B, int ldb,
					       PASTIX_FLOAT beta,
					       PASTIX_FLOAT       *C, int ldc,
					       int blocknbr,
					       const int * blocktab,
					       int fblocknbr,
					       const int * fblocktab)
{
#if (defined PREC_DOUBLE && defined TYPE_COMPLEX)
  dim3 threads( 16, 4 );
  dim3 grid(m/64+(m%64!=0),n/16+(n%16!=0));
#else
  dim3 threads( 16, 4 );
  dim3 grid(m/64+(m%64!=0),n/16+(n%16!=0));
#endif
  sparse_gemdm_kernel_N_T_64_16_4_16_4
    <<< grid, threads, 0, 0 /*starpu_cuda_get_local_stream()*/ >>>(m,n,k,
							     alpha,
							     A, lda,
							     D, ldd,
							     B, ldb,
							     beta,
							     C, ldc,
							     blocknbr,
							     blocktab,
							     fblocknbr,
							     fblocktab);
}
