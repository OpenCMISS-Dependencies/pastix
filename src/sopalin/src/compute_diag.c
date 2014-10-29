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
 * File: compute_diag.c
 *
 * Computation functions.
 *
 * Authors:
 *   Pierre Ramet    : fev 2003
 *   Mathieu Faverge
 *   Xavier Lacoste
 */

/*
 * Constant: MAXSIZEOFBLOCKS
 *  Maximum size of blocks given to blas in factorization
 */
#define MAXSIZEOFBLOCKS   64    /*64 in LAPACK*/

#define PASTIX_getrf       API_CALL(PASTIX_getrf)
#define PASTIX_potrf       API_CALL(PASTIX_potrf)
#define PASTIX_sytrf       API_CALL(PASTIX_sytrf)
#define PASTIX_hetrf       API_CALL(PASTIX_hetrf)
#define PASTIX_getrf_block API_CALL(PASTIX_getrf_block)
#define PASTIX_potrf_block API_CALL(PASTIX_potrf_block)
#define PASTIX_sytrf_block API_CALL(PASTIX_sytrf_block)
#define PASTIX_hetrf_block API_CALL(PASTIX_hetrf_block)
#define DimTrans           API_CALL(DimTrans)
#define factor_diag        API_CALL(factor_diag)

void PASTIX_getrf       ( PASTIX_FLOAT *A, PASTIX_INT m, PASTIX_INT n, PASTIX_INT lda, PASTIX_INT *npvt, double crit );
void PASTIX_potrf       ( PASTIX_FLOAT *A, PASTIX_INT n,        PASTIX_INT lda, PASTIX_INT *npvt, double crit );
void PASTIX_sytrf       ( PASTIX_FLOAT *A, PASTIX_INT n,        PASTIX_INT lda, PASTIX_INT *npvt, double crit );
void PASTIX_hetrf       ( PASTIX_FLOAT *A, PASTIX_INT n,        PASTIX_INT lda, PASTIX_INT *npvt, double crit );
void PASTIX_getrf_block ( PASTIX_FLOAT *A, PASTIX_INT m, PASTIX_INT n, PASTIX_INT lda, PASTIX_INT *npvt, double crit );
void PASTIX_potrf_block ( PASTIX_FLOAT *A, PASTIX_INT n,        PASTIX_INT lda, PASTIX_INT *npvt, double crit );
void PASTIX_sytrf_block ( PASTIX_FLOAT *A, PASTIX_INT n,        PASTIX_INT lda, PASTIX_INT *npvt, double crit, PASTIX_FLOAT * tmp4   );
void PASTIX_hetrf_block ( PASTIX_FLOAT *A, PASTIX_INT n,        PASTIX_INT lda, PASTIX_INT *npvt, double crit, PASTIX_FLOAT * tmp4   );
void DimTrans           ( PASTIX_FLOAT *A, PASTIX_INT lda, PASTIX_INT size, PASTIX_FLOAT *B );


/****************************************************************************/
/* COMPUTE TASK DIAG                                                        */
/****************************************************************************/

#define TALIGNF (pa,n,lda,pb) {                 \
    PASTIX_INT i,l; PASTIX_FLOAT *pt,*p1,*p2;                 \
    p1= (pa); p2= (pb); l= (lda)+1;             \
    for (i=(n);i>0;i--) {                       \
      for (pt=p1+i;p1<pt;)                      \
        * (p2++)=*(p1++);                       \
      p1+=l-i; } }

#define TALIGNB (pa,n,lda,pb) {                 \
    PASTIX_INT i,l; PASTIX_FLOAT *pt,*p1,*p2;                 \
    p1= (pa); p2=(pb); l=(lda)+1;               \
    for (i=(n);i>0;i--) {                       \
      for (pt=p1+i;p1<pt;)                      \
        * (p1++)=*(p2++);                       \
      p1+=l-i; } }

/*
#define TALIGNF (pa,n,lda,pb) {PASTIX_INT i,j,k=0;\
  for (j=0;j<n;j++)\
  for (i=j;i<n;i++)\
  {\
  (pb)[k]=(pa)[j*lda+i];\
  k=k+1;\
  }}\

#define TALIGNB (pa,n,lda,pb) {PASTIX_INT i,j,k=(n*(n+1))/2-1;\
  for (j=n-1;j>=0;j--)\
  for (i=n-1;i>=j;i--)\
  {\
  (pa)[j*lda+i]=(pb)[k];\
  k=k-1;\
  }}\
*/

/*
 *  Function: PASTIX_potrf
 *
 *  Factorization LLt BLAS2 3 terms
 *
 * > A = LL^T
 *
 * Parameters:
 *    A       - Matrix to factorize
 *    n       - Size of A
 *    stride  - Stide between 2 columns of the matrix
 *    nbpivot - IN/OUT pivot number.
 *    critere - Pivoting threshold.
 *
 */
void
PASTIX_potrf (PASTIX_FLOAT * A, PASTIX_INT n, PASTIX_INT stride, PASTIX_INT *nbpivot, double critere)
{
  PASTIX_INT k;
  PASTIX_FLOAT *tmp,*tmp1;

  for (k=0;k<n;k++)
    {
      tmp=A+k* (stride+1);
#ifdef USE_CSC
      if (ABS_FLOAT(*tmp)<critere)
        {
          (*tmp) = (PASTIX_FLOAT)critere;
          (*nbpivot)++;
        }
#endif
#ifdef TYPE_COMPLEX
      *tmp = (PASTIX_FLOAT)csqrt(*tmp);
#else
      *tmp = (PASTIX_FLOAT)sqrt(*tmp);
      if (*tmp < 0)
        {
          errorPrint ("Negative diagonal term\n");
          EXIT (MOD_SOPALIN, INTERNAL_ERR);
        }
#endif
      tmp1=tmp+1;
      SOPALIN_SCAL (n-k-1,(fun/(*tmp)),tmp1,iun);
      SOPALIN_SYR ("L",n-k-1,-fun,tmp1,iun,tmp1+stride,stride);
    }
}

/*
 * Function: PASTIX_potrf_block
 *
 * Computes the block LL^T factorization of the
 * matrix A.
 *
 * > A = LL^T
 *
 * Parameters:
 *    A       - Matrix to factorize
 *    n       - Size of A
 *    stride  - Stide between 2 columns of the matrix
 *    nbpivot - IN/OUT pivot number.
 *    critere - Pivoting threshold.
 */
void
PASTIX_potrf_block (PASTIX_FLOAT * A, PASTIX_INT n, PASTIX_INT stride, PASTIX_INT *nbpivot,
                    double critere)
{
  PASTIX_FLOAT *tmp,*tmp1,*tmp2;
  PASTIX_INT    k, blocknbr, blocksize, matrixsize;

  blocknbr = (PASTIX_INT)ceil((double)n/(double)MAXSIZEOFBLOCKS);
  for (k=0; k<blocknbr; k++)
    {
      blocksize = MIN (MAXSIZEOFBLOCKS,n-k*MAXSIZEOFBLOCKS);
      tmp  = A+ (k*MAXSIZEOFBLOCKS)*(stride+1); /* Lk,k     */
      tmp1 = tmp+ blocksize;                   /* Lk+1,k   */
      tmp2 = tmp1 + stride* blocksize;         /* Lk+1,k+1 */

      /* Factorize the diagonal block Akk*/
      PASTIX_potrf ( tmp,blocksize, stride, nbpivot, critere );
      if ((k*MAXSIZEOFBLOCKS+blocksize) < n)
        {
          matrixsize = n- (k*MAXSIZEOFBLOCKS+blocksize);
          /* Compute the column Lk+1k */
          SOPALIN_TRSM ("R","L","T","N",
                        matrixsize,
                        blocksize,
                        fun, tmp,stride,
                        tmp1,stride);
          /* Update Ak+1k+1 = Ak+1k+1 - Lk+1k*Lk+1kT */
          SOPALIN_SYRK ("L","N",
                        matrixsize,blocksize,
                        -fun,tmp1,stride,
                        fun,tmp2,stride);
        }
    }
}


/*
 *  Function: PASTIX_sytrf
 *
 *  Factorization LDLt BLAS2 3 terms
 *
 *  In complex : A is symmetric.
 *
 * > A = LDL^T
 *
 * Parameters:
 *    A       - Matrix to factorize
 *    n       - Size of A
 *    stride  - Stide between 2 columns of the matrix
 *    nbpivot - IN/OUT pivot number.
 *    critere - Pivoting threshold.
 */
void
PASTIX_sytrf ( PASTIX_FLOAT * A, PASTIX_INT n, PASTIX_INT stride, PASTIX_INT *nbpivot, double critere)
{
  PASTIX_INT k;
  PASTIX_FLOAT *tmp,*tmp1;

  for (k=0;k<n;k++)
    {
      tmp=A+k* (stride+1);
#ifdef USE_CSC
      if (ABS_FLOAT(*tmp)<critere)
        {
          (*tmp) = (PASTIX_FLOAT)critere;
          (*nbpivot)++;
        }
#endif
      tmp1=tmp+1;
      SOPALIN_SCAL (n-k-1,(fun/(*tmp)),tmp1,iun);
      SOPALIN_SYR ("L",n-k-1,-(*tmp),tmp1,iun,tmp1+stride,stride);
    }
}

/*
 * Function: PASTIX_sytrf_block
 *
 * Computes the block LDL^T factorization of the
 * matrix A.
 *
 *  In complex : A is symmetric.
 *
 * > A = LDL^T
 *
 * Parameters:
 *    A       - Matrix to factorize
 *    n       - Size of A
 *    stride  - Stide between 2 columns of the matrix
 *    nbpivot - IN/OUT pivot number.
 *    critere - Pivoting threshold.
 */
void
PASTIX_sytrf_block ( PASTIX_FLOAT * A, PASTIX_INT n, PASTIX_INT stride, PASTIX_INT *nbpivot,
                     double critere, PASTIX_FLOAT * tmp4)
{
  PASTIX_INT k,blocknbr,blocksize,matrixsize,col;
  PASTIX_FLOAT *tmp,*tmp1,*tmp2;
  PASTIX_FLOAT alpha;

  blocknbr = (PASTIX_INT)ceil((double)n/(double)MAXSIZEOFBLOCKS);
  for (k=0;k<blocknbr;k++)
    {
      blocksize = MIN (MAXSIZEOFBLOCKS,n-k*MAXSIZEOFBLOCKS);
      tmp  = A+ (k*MAXSIZEOFBLOCKS)*(stride+1); /* Lk,k     */
      tmp1 = tmp+ blocksize;                   /* Lk+1,k   */
      tmp2 = tmp1 + stride* blocksize;         /* Lk+1,k+1 */

      /* Factorize the diagonal block Akk*/
      PASTIX_sytrf (tmp,blocksize, stride, nbpivot, critere);
      if ((k*MAXSIZEOFBLOCKS+blocksize) < n)
        {
          matrixsize = n- (k*MAXSIZEOFBLOCKS+blocksize);
          /* Compute the column Lk+1k */
          /** Compute Dk,k*Lk+1,k      */
          SOPALIN_TRSM ("R","L","T","U",
                        matrixsize,
                        blocksize,
                        fun, tmp,stride,
                        tmp1,stride);
          for (col = 0; col < blocksize; col++)
            {
              /** Copy Dk,k*Lk+1,k and compute Lk+1,k */
              SOPALIN_COPY (matrixsize, tmp1+col*stride, iun,
                            tmp4+col*matrixsize,iun);
              alpha = fun / *(tmp + col*(stride+1));
              SOPALIN_SCAL (matrixsize, alpha,
                            tmp1+col*stride,
                            iun);
            }
          /* Update Ak+1k+1 = Ak+1k+1 - Lk+1k*Dk,k*Lk+1kT */
          SOPALIN_GEMM ("N","T",matrixsize,matrixsize,
                        blocksize,
                        -fun,tmp4,matrixsize,
                        tmp1,stride,
                        fun,tmp2,stride);
        }
    }
}

/*
 *  Function: PASTIX_hetrf
 *
 *  Factorization LDLt BLAS2 3 terms.
 *
 *  In complex : A is hermitian.
 *
 * > A = LDL^T
 *
 * Parameters:
 *    A       - Matrix to factorize
 *    n       - Size of A
 *    stride  - Stide between 2 columns of the matrix
 *    nbpivot - IN/OUT pivot number.
 *    critere - Pivoting threshold.
 */
void
PASTIX_hetrf ( PASTIX_FLOAT * A, PASTIX_INT n, PASTIX_INT stride, PASTIX_INT *nbpivot, double critere)
{
  PASTIX_INT k;
  PASTIX_FLOAT *tmp,*tmp1;

  for (k=0;k<n;k++)
    {
      tmp=A+k* (stride+1);
#ifdef USE_CSC
      if (ABS_FLOAT(*tmp)<critere)
        {
          (*tmp) = (PASTIX_FLOAT)critere;
          (*nbpivot)++;
        }
#endif
      tmp1=tmp+1;
      SOPALIN_SCAL (n-k-1,(fun/(*tmp)),tmp1,iun);
      SOPALIN_HER ("L",n-k-1,-(*tmp),tmp1,iun,tmp1+stride,stride);
    }
}

/*
 * Function: PASTIX_hetrf_block
 *
 * Computes the block LDL^T factorization of the
 * matrix A.
 *
 *  In complex : A is hermitian.
 *
 * > A = LDL^T
 *
 * Parameters:
 *    A       - Matrix to factorize
 *    n       - Size of A
 *    stride  - Stide between 2 columns of the matrix
 *    nbpivot - IN/OUT pivot number.
 *    critere - Pivoting threshold.
 */
void
PASTIX_hetrf_block ( PASTIX_FLOAT * A, PASTIX_INT n, PASTIX_INT stride, PASTIX_INT *nbpivot,
                     double critere, PASTIX_FLOAT * tmp4)
{
  PASTIX_INT k,blocknbr,blocksize,matrixsize,col;
  PASTIX_FLOAT *tmp,*tmp1,*tmp2;
  PASTIX_FLOAT alpha;

  blocknbr = (PASTIX_INT)ceil((double)n/(double)MAXSIZEOFBLOCKS);
  for (k=0;k<blocknbr;k++)
    {
      blocksize = MIN (MAXSIZEOFBLOCKS,n-k*MAXSIZEOFBLOCKS);
      tmp  = A+ (k*MAXSIZEOFBLOCKS)*(stride+1); /* Lk,k     */
      tmp1 = tmp+ blocksize;                   /* Lk+1,k   */
      tmp2 = tmp1 + stride* blocksize;         /* Lk+1,k+1 */

      /* Factorize the diagonal block Akk*/
      PASTIX_hetrf (tmp,blocksize, stride, nbpivot, critere);
      if ((k*MAXSIZEOFBLOCKS+blocksize) < n)
        {
          matrixsize = n- (k*MAXSIZEOFBLOCKS+blocksize);
          /* Compute the column Lk+1k */
          /** Compute Dk,k*Lk+1,k      */
          SOPALIN_TRSM ("R","L","C","U",
                        matrixsize,
                        blocksize,
                        fun, tmp,stride,
                        tmp1,stride);
          for (col = 0; col < blocksize; col++)
            {
              /** Copy Dk,k*Lk+1,k and compute Lk+1,k */
              SOPALIN_COPY (matrixsize, tmp1+col*stride, iun,
                            tmp4+col*matrixsize,iun);
              alpha = fun / *(tmp + col*(stride+1));
              SOPALIN_SCAL (matrixsize, alpha,
                            tmp1+col*stride,
                            iun);
            }
          /* Update Ak+1k+1 = Ak+1k+1 - Lk+1k*Dk,k*Lk+1kT */
          SOPALIN_GEMM ("N","C",matrixsize,matrixsize,
                        blocksize,
                        -fun,tmp4,matrixsize,
                        tmp1,stride,
                        fun,tmp2,stride);
        }
    }
}
/*
 *  Function: PASTIX_getrf
 *
 *  LU Factorization of one (diagonal) block
 *  $A = LU$
 *
 *  For each column :
 *    - Divide the column by the diagonal element.
 *    - Substract the product of the subdiagonal part by
 *      the line after the diagonal element from the
 *      matrix under the diagonal element.
 *
 * Parameters:
 *    A       - Matrix to factorize
 *    m       - number of rows of the Matrix A
 *    n       - number of cols of the Matrix A
 *    stride  - Stide between 2 columns of the matrix
 *    nbpivot - IN/OUT pivot number.
 *    critere - Pivoting threshold.
 */
void
PASTIX_getrf ( PASTIX_FLOAT *A, PASTIX_INT m, PASTIX_INT n, PASTIX_INT stride, PASTIX_INT *nbpivot,
               double critere)
{
  PASTIX_INT j;
  PASTIX_FLOAT *tmp;
  PASTIX_FLOAT *tmp1;

  for (j=0; j<MIN(m,n); j++)
    {
      tmp = A + j* (stride+1); /* A[j][j] */
      tmp1 = tmp+1; /* A[j+1][j] */
#ifdef USE_CSC
      if (ABS_FLOAT(*tmp) < critere)
        {
          (*tmp) = (PASTIX_FLOAT)critere;
          (*nbpivot)++;
        }
#endif
      /* A[k][j] = A[k][j]/A[j] (j], k = j+1 .. n */
      SOPALIN_SCAL ((m-j-1), (fun/(*tmp)), tmp1, iun);
      if (j +1 < MIN(m,n))
        {
          /* A[k][l] = A[k][l] - A[k][j]*A[j,l] , k,l = j+1..n*/
          SOPALIN_GER ((m-j-1), (n-j-1), -fun, tmp1, iun,
                       tmp+stride, stride, tmp+stride+1, stride);
        }
    }

  /* Test sur la dernier valeur diag */
  tmp = A + (n-1)*(stride+1);
#ifdef USE_CSC
  if (ABS_FLOAT(*tmp) < critere)
    {
      (*tmp) = (PASTIX_FLOAT)critere;
      (*nbpivot)++;
    }
#endif
}

/*
 *  Function: PASTIX_getrf_block
 *
 *  Block LU Factorization of one (diagonal) big block
 *  > A = LU
 *
 * Parameters:
 *    A       - Matrix to factorize.
 *    rows    - Number of rows.
 *    cols    - Number of columns.
 *    stride  - Stide between 2 columns of the matrix.
 *    nbpivot - IN/OUT pivot number.
 *    critere - Pivoting threshold.
 */
void
PASTIX_getrf_block (PASTIX_FLOAT *A, PASTIX_INT rows, PASTIX_INT cols, PASTIX_INT stride, PASTIX_INT *nbpivot,
                    double critere)
{
  PASTIX_INT    k,blocknbr,blocksize,matsize;
  PASTIX_FLOAT *tmp,*tmp1,*tmp2,*tmp3;

  blocknbr = (PASTIX_INT)ceil((double)cols/(double)MAXSIZEOFBLOCKS);
  for (k=0;k<blocknbr;k++)
    {
      blocksize = MIN (MAXSIZEOFBLOCKS,cols-k*MAXSIZEOFBLOCKS);
      tmp  = A+ (k*MAXSIZEOFBLOCKS)*(stride+1);
      tmp1 = tmp + blocksize;            /* Lk+1,k   */
      tmp2 = tmp + stride*blocksize;     /* Uk,k+1   */
      tmp3 = tmp + (stride+1)*blocksize; /* Ak+1,k+1 */
      /* Factorize the diagonal block Akk*/
      PASTIX_getrf (tmp, (rows - k*MAXSIZEOFBLOCKS), blocksize, stride, nbpivot,
                    critere);
      if ((k*MAXSIZEOFBLOCKS+blocksize) < cols)
        {
          matsize = rows - k*MAXSIZEOFBLOCKS-blocksize;
          /* Compute the column Ukk+1 */
          SOPALIN_TRSM ("L","L","N","U",
                        blocksize,
                        matsize,
                        fun, tmp,stride,tmp2,stride);
          /* Update Ak+1,k+1 = Ak+1,k+1 - Lk+1,k*Uk,k+1 */
          SOPALIN_GEMM ("N","N",matsize, matsize,blocksize, -fun, tmp1, stride,
                        tmp2, stride, fun, tmp3, stride);

        }
    }

}

void
DimTrans (PASTIX_FLOAT *A, PASTIX_INT stride, PASTIX_INT size, PASTIX_FLOAT *B)
{
  PASTIX_INT i,j;

  for (i=0; i<size; i++)
    {
      for (j=0; j<size; j++)
        {
          B[i*stride+j] = A[j*stride+i];
        }
    }
}


/*
 * Factorization of diagonal block
 */
void factor_diag (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT c)
{

  PASTIX_INT    size,stride;
  PASTIX_FLOAT *ga = NULL;
#ifdef SOPALIN_LU
  PASTIX_FLOAT *gb = NULL;
#endif
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];

  /* check if diagonal column block */
  ASSERTDBG ( SYMB_FCOLNUM(c) == SYMB_FROWNUM(SYMB_BLOKNUM(c)),
              MOD_SOPALIN );

  /* Initialisation des pointeurs de blocs */
  ga = & (SOLV_COEFTAB(c)[ SOLV_COEFIND(SYMB_BLOKNUM(c))]);
#ifdef SOPALIN_LU
  gb = & (SOLV_UCOEFTAB(c)[SOLV_COEFIND(SYMB_BLOKNUM(c))]);
#endif
  size   = SYMB_LCOLNUM (c)-SYMB_FCOLNUM(c)+1;
  stride = SOLV_STRIDE (c);

#ifdef COMPUTE
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  PASTIX_getrf_block (ga, size, size, stride,
                      & (thread_data->nbpivot),
                      sopalin_data->critere);
  DimTrans (ga,stride,size,gb);
#    else /* SOPALIN_LU */
  /*SOPALIN_POF ("L",ga,stride,size);*/
  PASTIX_potrf_block (ga, size, stride,
                      & (thread_data->nbpivot),
                      sopalin_data->critere);

#    endif /* SOPALIN_LU */

#  else /* CHOL_SOPALIN */

  /* version avec PPF full storage
     gb=maxbloktab1[me];
     TALIGNF (ga,size,stride,gb);
     SOPALIN_PPF (gb,size,izero);
     TALIGNB (ga,size,stride,gb); */
#    ifdef HERMITIAN
  PASTIX_hetrf_block (ga, size, stride,
                      & (thread_data->nbpivot),
                      sopalin_data->critere,
                      thread_data->maxbloktab1);
#    else
  PASTIX_sytrf_block (ga, size, stride,
                      & (thread_data->nbpivot),
                      sopalin_data->critere,
                      thread_data->maxbloktab1);
#    endif
  /* version BLAS 3 */
  /* Copy diagonal for esp tasks */
  if (sopalin_data->sopar->iparm[IPARM_ESP])
    {
      stride++;
      SOPALIN_COPY (size, ga, stride,
                    & (SOLV_UCOEFTAB(c)[SOLV_COEFIND(SYMB_BLOKNUM(c))]), iun);
    }

#  endif /* CHOL_SOPALIN */
#endif /* COMPUTE */
}
