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
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @precisions normal z -> s d c

*/
#include <assert.h>
#ifdef CHOL_SOPALIN
#  undef CHOL_SOPALIN
#endif /* CHOL_SOPALIN */

#ifndef HERMITIAN
#  define HERMITIAN
#endif /* not HERMITIAN */

#include "common_pastix.h"
#include "sopalin_define.h"
#define inline static inline
#include "magma.h"
#undef inline

// === Define what BLAS to use ============================================
#define min                         MIN
#define max                         MAX
#ifdef TYPE_COMPLEX
#  ifdef PREC_DOUBLE
#    define magma_get_zhetrf_nb     magma_get_zgetrf_nb
#  else /* not PREC_DOUBLE */
#    define magma_zgemm             magma_cgemm
#    define magma_ztrsm             magma_ctrsm
#    define magma_zherk             magma_cherk
#    define magma_zgetmatrix        magma_cgetmatrix
#    define magma_zsetmatrix        magma_csetmatrix
#    define magma_zgetmatrix_async  magma_cgetmatrix_async
#    define magma_zsetmatrix_async  magma_csetmatrix_async
#    define magma_zmalloc_host      magma_cmalloc_host
#    define magma_get_zhetrf_nb     magma_get_cgetrf_nb
#    undef  MAGMA_Z_ONE
#    undef  MAGMA_Z_NEG_ONE
#    define MAGMA_Z_ONE             MAGMA_C_ONE
#    define MAGMA_Z_NEG_ONE         MAGMA_C_NEG_ONE
#  endif /* not PREC_DOUBLE */
#else /* not TYPE_COMPLEX */
#  ifdef PREC_DOUBLE
#    define magma_zgemm             magmablas_dgemm
#    define magma_ztrsm             magmablas_dtrsm
#    define magma_zherk             magmablas_dherk
#    define magma_zgetmatrix        magma_dgetmatrix
#    define magma_zsetmatrix        magma_dsetmatrix
#    define magma_zgetmatrix_async  magma_dgetmatrix_async
#    define magma_zsetmatrix_async  magma_dsetmatrix_async
#    define magma_zmalloc_host      magma_dmalloc_host
#    define magma_get_zhetrf_nb     magma_get_dgetrf_nb
#    undef  MAGMA_Z_ONE
#    undef  MAGMA_Z_NEG_ONE
#    define MAGMA_Z_ONE             MAGMA_D_ONE
#    define MAGMA_Z_NEG_ONE         MAGMA_D_NEG_ONE
#  else /* not PREC_DOUBLE */
#    define magma_zgemm             magmablas_sgemm
#    define magma_ztrsm             magmablas_strsm
#    define magma_zherk             magmablas_sherk
#    define magma_zgetmatrix        magma_sgetmatrix
#    define magma_zsetmatrix        magma_ssetmatrix
#    define magma_zgetmatrix_async  magma_sgetmatrix_async
#    define magma_zsetmatrix_async  magma_ssetmatrix_async
#    define magma_zmalloc_host      magma_smalloc_host
#    define magma_get_zhetrf_nb     magma_get_sgetrf_nb
#    undef  MAGMA_Z_ONE
#    undef  MAGMA_Z_NEG_ONE
#    define MAGMA_Z_ONE             MAGMA_S_ONE
#    define MAGMA_Z_NEG_ONE         MAGMA_S_NEG_ONE
#  endif /* not PREC_DOUBLE */
#endif /* not TYPE_COMPLEX */

#include "pastix_cuda_helper.h"
#include "zhetrf_stapiv_gpu.h"

#define cuDoubleComplex CU_FLOAT
#define uplo_lapackf77_lsame(a, b) ((a[0] == b[0])?1:0)
#define lapackf77_lsame(a, b) ((a[0] == b[0])?1:0)

#define PASTIX_hetrf_block API_CALL(PASTIX_hetrf_block)
void PASTIX_hetrf_block ( PASTIX_FLOAT *A, PASTIX_INT n, PASTIX_INT lda, PASTIX_INT *npvt, double crit, PASTIX_FLOAT * tmp4);

// === End defining what BLAS to use =======================================

#define dA(i, j)  (dA + (j)*ldda + (i))

magma_int_t
magma_zhetrf_stapiv_gpu(char uplo, magma_int_t n,
                        cuDoubleComplex *dA, magma_int_t ldda,
                        double criteria, PASTIX_INT * npiv, PASTIX_FLOAT * tmp4, magma_int_t *info)
{
/*  -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

    Purpose
    =======
    ZHETRF computes the Cholesky factorization of a complex Hermitian
    positive definite matrix dA.

    The factorization has the form
       dA = U**H * U,  if UPLO = 'U', or
       dA = L  * L**H,  if UPLO = 'L',
    where U is an upper triangular matrix and L is lower triangular.

    This is the block version of the algorithm, calling Level 3 BLAS.

    Arguments
    =========
    UPLO    (input) CHARACTER*1
            = 'U':  Upper triangle of dA is stored;
            = 'L':  Lower triangle of dA is stored.

    N       (input) INTEGER
            The order of the matrix dA.  N >= 0.

    dA      (input/output) COMPLEX_16 array on the GPU, dimension (LDDA,N)
            On entry, the Hermitian matrix dA.  If UPLO = 'U', the leading
            N-by-N upper triangular part of dA contains the upper
            triangular part of the matrix dA, and the strictly lower
            triangular part of dA is not referenced.  If UPLO = 'L', the
            leading N-by-N lower triangular part of dA contains the lower
            triangular part of the matrix dA, and the strictly upper
            triangular part of dA is not referenced.

            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization dA = U**H * U or dA = L * L**H.

    LDDA     (input) INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,N).
            To benefit from coalescent memory accesses LDDA must be
            dividable by 16.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
            > 0:  if INFO = i, the leading minor of order i is not
                  positive definite, and the factorization could not be
                  completed.
    =====================================================================   */


    magma_int_t     j, jb, nb;
    char            uplo_[2] = {uplo, 0};
    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
    cuDoubleComplex *work;
    double          d_one     =  1.0;
    double          d_neg_one = -1.0;
    long int        upper = uplo_lapackf77_lsame(uplo_, "U");

    *info = 0;
    if ( (! upper) && (! lapackf77_lsame(uplo_, "L")) ) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    nb = magma_get_zhetrf_nb(n);

    if (MAGMA_SUCCESS != magma_zmalloc_host( &work, nb*nb )) {
        *info = MAGMA_ERR_HOST_ALLOC;
        return *info;
    }

    static cudaStream_t stream[2];
    magma_queue_create( &stream[0] );
    magma_queue_create( &stream[1] );

    if ((nb <= 1) || (nb >= n)) {
        /*  Use unblocked code. */
        magma_zgetmatrix( n, n, dA, ldda, work, n );
        assert(!upper); /* PaStiX only works with lower */
        PASTIX_hetrf_block((PASTIX_FLOAT*)work, n, n,
                           npiv,
                           criteria, tmp4);
        magma_zsetmatrix( n, n, work, n, dA, ldda );
    } else {
        /* Use blocked code. */
        if (upper) {
          assert(0); /* PaStiX only works with lower */
      
            /* Compute the Cholesky factorization A = U'*U. */
            for (j=0; j<n; j+=nb) {

                /* Update and factorize the current diagonal block and test
                   for non-positive-definiteness. Computing MIN */
                jb = min(nb, (n-j));

                magma_zherk(MagmaUpper, MagmaConjTrans, jb, j,
                            d_neg_one, dA(0, j), ldda,
                            d_one,     dA(j, j), ldda);

                magma_zgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[1] );

                if ( (j+jb) < n) {
                    /* Compute the current block row. */
                    magma_zgemm(MagmaConjTrans, MagmaNoTrans,
                                jb, (n-j-jb), j,
                                c_neg_one, dA(0, j   ), ldda,
                                           dA(0, j+jb), ldda,
                                c_one,     dA(j, j+jb), ldda);
                }

                magma_queue_sync( stream[1] );

                /* lapackf77_zhetrf(MagmaUpperStr, &jb, work, &jb, info); */
                magma_zsetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[0] );
                if (*info != 0) {
                  *info = *info + j;
                  break;
                }

                if ( (j+jb) < n)
                    magma_ztrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit,
                                 jb, (n-j-jb),
                                 c_one, dA(j, j   ), ldda,
                                        dA(j, j+jb), ldda);
            }
        } else {
            //=========================================================
            // Compute the Cholesky factorization A = L*L'.
            for (j=0; j<n; j+=nb) {

                //  Update and factorize the current diagonal block and test
                //  for non-positive-definiteness. Computing MIN
                jb = min(nb, (n-j));

                magma_zherk(MagmaLower, MagmaNoTrans, jb, j,
                            d_neg_one, dA(j, 0), ldda,
                            d_one,     dA(j, j), ldda);

                magma_zgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[1] );

                if ( (j+jb) < n) {
                    magma_zgemm( MagmaNoTrans, MagmaConjTrans,
                                 (n-j-jb), jb, j,
                                 c_neg_one, dA(j+jb, 0), ldda,
                                            dA(j,    0), ldda,
                                 c_one,     dA(j+jb, j), ldda);
                }

                magma_queue_sync( stream[1] );
                PASTIX_hetrf_block((PASTIX_FLOAT*)work, jb, jb,
                                   npiv,
                                   criteria, tmp4);

                magma_zsetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[0] );
                if (*info != 0) {
                  *info = *info + j;
                  break;
                }

                if ( (j+jb) < n)
                    magma_ztrsm(MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit,
                                (n-j-jb), jb,
                                c_one, dA(j,    j), ldda,
                                       dA(j+jb, j), ldda);
            }

        }
    }

    magma_queue_destroy( stream[0] );
    magma_queue_destroy( stream[1] );
    magma_free_host( work );

    return *info;
} /* magma_zhetrf_gpu */
