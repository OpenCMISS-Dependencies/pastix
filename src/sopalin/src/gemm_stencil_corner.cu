//on garde tous les tests sur M et N



 ///////////////////////////////////////////////////////////////////////////////////////////////////

// size of work for a thread
#define THR_M ( BLK_M / DIM_X )
#define THR_N ( BLK_N / DIM_Y )
 
///////////////////////////////////////////////////////////////////////////////////////////////////

#if   (version == trans_nn)
#define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(gemm_corner_nn)
#elif (version == trans_nt)
#define TRANS_B
#define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(gemm_corner_nt)
#elif (version == trans_nc)
#define TRANS_B 
#define CONJ_B
#define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(gemm_corner_nc)
#elif (version == trans_tn)
#define TRANS_A
#define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(gemm_corner_tn)
#elif (version == trans_tt)
#define TRANS_A
#define TRANS_B 
#define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(gemm_corner_tt)
#elif (version == trans_tc)
#define TRANS_A
#define TRANS_B
#define CONJ_B
#define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(gemm_corner_tc)
#elif (version == trans_cn)
#define TRANS_A
#define CONJ_A
#define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(gemm_corner_cn)
#elif (version == trans_ct)
#define TRANS_A
#define CONJ_A
#define TRANS_B
#define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(gemm_corner_ct)
#elif (version == trans_cc)
#define TRANS_A
#define CONJ_A
#define TRANS_B
#define CONJ_B
#define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(gemm_corner_cc)
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////

//#define MAX(a,b) (((a)>(b))?(a):(b))
//#define MIN(a,b) (((a)<(b))?(a):(b))

extern "C" __global__
void kernel_name (int M, int N, int K,
		  FloatingPoint_t alpha, 
		  const FloatingPoint_t *A, int LDA,
		  const FloatingPoint_t *B, int LDB, 
		  FloatingPoint_t beta, 
		  FloatingPoint_t *C, int LDC, 
		  int offsetA, int offsetB,
		  int blocknbr, const int *blocktab, 
		  int fblocknbr, const int *fblocktab)
{
  int offset[THR_M+1];

  int idx = threadIdx.x;  // thread's m dimension
  int idy = threadIdx.y;  // thread's n dimension
  
  int idt = DIM_X * idy + idx;    // thread's global number
  
  int idxA = idt % DIM_XA;    // idx within A
  int idyA = idt / DIM_XA;    // idy within A
  
  int idxB = idt % DIM_XB;    // idx within B
  int idyB = idt / DIM_XB;    // idy within B
  
  int blx = M/BLK_M; //blockIdx.x;   // block's m dimension
  int bly = N/BLK_N; //blockIdx.y;   // block's n dimension
  
  __shared__ FloatingPoint_t sA[BLK_K][BLK_M+1];      // +1 only required if A is transposed
  __shared__ FloatingPoint_t sB[BLK_N][BLK_K+1];      // +1 always required
  
  // Registers for the innermost loop
  FloatingPoint_t rC[THR_N][THR_M];
  FloatingPoint_t rA[THR_M];
  FloatingPoint_t rB[THR_N];
  
  
#ifdef TRANS_A
  const FloatingPoint_t *offs_dA = A + blx*BLK_M*LDA + idyA*LDA+idxA;
#else
  const FloatingPoint_t *offs_dA = A + blx*BLK_M     + idyA*LDA+idxA;
#endif
#ifdef TRANS_B
  const FloatingPoint_t *offs_dB = B + bly*BLK_N     + idyB*LDB+idxB;
#else
  const FloatingPoint_t *offs_dB = B + bly*BLK_N*LDB + idyB*LDB+idxB;
#endif

  
  int m, n, k, kk;
  int coordm,coordn;
  int kk_aux=-1*BLK_K;

  // Zero C
#pragma unroll
  for (n = 0; n < THR_N; n++)
#pragma unroll
    for (m = 0; m < THR_M; m++)
      rC[n][m] = make_FloatingPoint(0.0, 0.0);
  
  for (kk = 0; kk < K-BLK_K; kk += BLK_K)
    {
      kk_aux=kk;
      // Load A dev->shmem
	
#ifdef TRANS_A
#pragma unroll
	for (n = 0; n < BLK_M; n += DIM_YA){
#pragma unroll
	  for (m = 0; m < BLK_K; m += DIM_XA){
	    //coordm = m + idxA + kk;
	    coordn = n + blx*BLK_M + idyA;
	    if(/*coordm < K && */coordn < M)
	      sA[m+idxA][n+idyA] = fetch(A, m, n);// offs_dA[n*LDA+m]
	    else sA[m+idxA][n+idyA] = make_FloatingPoint(0.0,0.0);
	  }
	}
#else
#pragma unroll
    for (n = 0; n < BLK_K; n += DIM_YA){
#pragma unroll
      for (m = 0; m < BLK_M; m += DIM_XA){
	coordm = m + blx*BLK_M + idxA;
	//coordn = n + idyA + kk;
	if(coordm < M)// && coordn < K)
	  sA[n+idyA][m+idxA] = fetch(A, m, n);
	else sA[n+idyA][m+idxA] = make_FloatingPoint(0.0,0.0);
      }
    }
    
    
#endif
    
    // Load B dev->shmem
#ifdef TRANS_B
#pragma unroll
    for (n = 0; n < BLK_K; n += DIM_YB){
#pragma unroll
      for (m = 0; m < BLK_N; m += DIM_XB){
	coordm = m + bly*BLK_N + idxB;
	//coordn = n + idyB + kk;
	if(coordm < N)// && coordn < K)
	  sB[m+idxB][n+idyB] = fetch(B, m, n);
	else sB[m+idxB][n+idyB] = make_FloatingPoint(0.0,0.0);
      }
    }
#else
#pragma unroll
    for (n = 0; n < BLK_N; n += DIM_YB){
#pragma unroll
      for (m = 0; m < BLK_K; m += DIM_XB){
	//coordm = m + idxB + kk;
	coordn = n + bly*BLK_N + idyB;
	if(/*coordm < K && */coordn < N)
	  sB[n+idyB][m+idxB] = fetch(B, m, n);
	else sB[n+idyB][m+idxB] = make_FloatingPoint(0.0,0.0);
      }
    }
#endif
    
    __syncthreads();
    
    
    // Multiply
#pragma unroll
    for (k = 0; k < BLK_K; k++)
      {
	// Load A shmem->regs
#pragma unroll
	for (m = 0; m < THR_M; m++)
	  rA[m] = sA[k][m*DIM_X+idx];
	
	// Load B shmem->regs
#pragma unroll
	for (n = 0; n < THR_N; n++)
	  rB[n] = sB[n*DIM_Y+idy][k];
	
	// Compute
#pragma unroll
	for (n = 0; n < THR_N; n++)
#pragma unroll
	  for (m = 0; m < THR_M; m++)
	    
#ifdef CONJ_A
#ifdef CONJ_B
	    fma(conj(rA[m]), conj(rB[n]), rC[n][m]);
#else
	fma(conj(rA[m]), rB[n], rC[n][m]);
#endif
#else
#ifdef CONJ_B
	fma(rA[m], conj(rB[n]), rC[n][m]);
#else
	fma(rA[m], rB[n], rC[n][m]);
#endif
#endif
      }
    
    __syncthreads();
    
    //maj offset
#ifdef TRANS_A
    offs_dA += BLK_K;
#else
    offs_dA += BLK_K*LDA;
#endif
#ifdef TRANS_B
    offs_dB += BLK_K*LDB;
#else
    offs_dB += BLK_K;
#endif
    
    __syncthreads();
    
    }

  ///////////////////////////////////////////////////////////////////////////////////


  kk_aux+=BLK_K;
  // Load A dev->shmem
  
#ifdef TRANS_A
  #pragma unroll
  for (n = 0; n < BLK_M; n += DIM_YA){
    #pragma unroll
    for (m = 0; m < BLK_K; m += DIM_XA){
      coordm = m + idxA + kk_aux;
      coordn = n + blx*BLK_M + idyA;
      if(coordm < K && coordn < M)
	sA[m+idxA][n+idyA] = fetch(A, m, n);// offs_dA[n*LDA+m]
      else sA[m+idxA][n+idyA] = make_FloatingPoint(0.0,0.0);
    }
  }
#else
  #pragma unroll
  for (n = 0; n < BLK_K; n += DIM_YA){
    #pragma unroll
    for (m = 0; m < BLK_M; m += DIM_XA){
      coordm = m + blx*BLK_M + idxA;
      coordn = n + idyA + kk_aux;
      if(coordm < M && coordn < K)
	sA[n+idyA][m+idxA] = fetch(A, m, n);
      else sA[n+idyA][m+idxA] = make_FloatingPoint(0.0,0.0);
    }
  }
  
    
#endif
    
    // Load B dev->shmem
#ifdef TRANS_B
  #pragma unroll
  for (n = 0; n < BLK_K; n += DIM_YB){
    #pragma unroll
    for (m = 0; m < BLK_N; m += DIM_XB){
      coordm = m + bly*BLK_N + idxB;
      coordn = n + idyB + kk_aux;
      if(coordm < N && coordn < K)
	sB[m+idxB][n+idyB] = fetch(B, m, n);
      else sB[m+idxB][n+idyB] = make_FloatingPoint(0.0,0.0);
    }
  }
#else
  #pragma unroll
  for (n = 0; n < BLK_N; n += DIM_YB){
    #pragma unroll
    for (m = 0; m < BLK_K; m += DIM_XB){
      coordm = m + idxB + kk_aux;
      coordn = n + bly*BLK_N + idyB;
      if(coordm < K && coordn < N)
	sB[n+idyB][m+idxB] = fetch(B, m, n);
      else sB[n+idyB][m+idxB] = make_FloatingPoint(0.0,0.0);
    }
  }
#endif
    
    __syncthreads();
    
    
    // Multiply
#pragma unroll
    for (k = 0; k < BLK_K; k++)
      {
	// Load A shmem->regs
#pragma unroll
	for (m = 0; m < THR_M; m++)
	  rA[m] = sA[k][m*DIM_X+idx];
	
	// Load B shmem->regs
#pragma unroll
	for (n = 0; n < THR_N; n++)
	  rB[n] = sB[n*DIM_Y+idy][k];
	
	// Compute
#pragma unroll
	for (n = 0; n < THR_N; n++)
#pragma unroll
	  for (m = 0; m < THR_M; m++)
	    
#ifdef CONJ_A
#ifdef CONJ_B
	    fma(conj(rA[m]), conj(rB[n]), rC[n][m]);
#else
	fma(conj(rA[m]), rB[n], rC[n][m]);
#endif
#else
#ifdef CONJ_B
	fma(rA[m], conj(rB[n]), rC[n][m]);
#else
	fma(rA[m], rB[n], rC[n][m]);
#endif
#endif
      }
    
    __syncthreads();
    
    //maj offset
#ifdef TRANS_A
    offs_dA += BLK_K;
#else
    offs_dA += BLK_K*LDA;
#endif
#ifdef TRANS_B
    offs_dB += BLK_K*LDB;
#else
    offs_dB += BLK_K;
#endif
    
    __syncthreads();
    ////////////////////////////////////////////////////////////////////////

    {
#define FROWNUM(tab, b) tab[2*b]
#define LROWNUM(tab, b) tab[2*b+1]
#define BLOCKSIZE(tab, b) LROWNUM(tab, b) - FROWNUM(tab, b) + 1
      int blocknum = 0, fblocknum = 0;
      size_t totalblocksize = 0;
      size_t blocksize = BLOCKSIZE(blocktab, blocknum);
      int    rownum;
      
      offset[0] = 0;
      for (m = 0; m < THR_M; m++) {
	int coord_dCm = blx*BLK_M + m*DIM_X+idx;
	
	if (coord_dCm < M) {
	  
	  /*
	   * We should keep blocknum < blocknbr
	   */
	  while( totalblocksize + blocksize < coord_dCm + 1)
	    {
	      totalblocksize += blocksize;
	      blocknum++;
	      blocksize = BLOCKSIZE(blocktab, blocknum);
	    }
	  
	  /* Global row index */
	  rownum = coord_dCm - totalblocksize + FROWNUM(blocktab, blocknum);
	  
	  while (LROWNUM(fblocktab, fblocknum) < rownum) {
	    offset[m] += BLOCKSIZE(fblocktab, fblocknum);
	    fblocknum++;
	  }
	  offset[m+1] = offset[m];
	  offset[m] += rownum - FROWNUM(fblocktab, fblocknum);
	}
      }
      __syncthreads();
#undef FROWNUM
#undef LROWNUM
    }
    
    
    
    // Store C regs->dev
#pragma unroll
    for (n = 0; n < THR_N; n++) {
        int coord_dCn = bly*BLK_N + n*DIM_Y+idy;
#pragma unroll
        for (m = 0; m < THR_M; m++) {
	  int coord_dCm = blx*BLK_M + m*DIM_X+idx;
	  if (coord_dCm < M && coord_dCn < N) {
	    int offsC = coord_dCn*LDC + offset[m]; /*coord_dCm;*/
	    
	    FloatingPoint_t &regC = rC[n][m];
	    FloatingPoint_t &memC = C[offsC];
	    
	    memC = add(mul(alpha, regC), mul(beta, memC));
	    }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#undef TRANS_A
#undef TRANS_B
#undef CONJ_A
#undef CONJ_B
/*
#undef BLK_M
#undef BLK_N
#undef BLK_K

#undef DIM_X
#undef DIM_Y

#undef DIM_XA
#undef DIM_YA 

#undef DIM_XB
#undef DIM_YB
*/
#undef version

#undef THR_M
#undef THR_N

#undef kernel_name
