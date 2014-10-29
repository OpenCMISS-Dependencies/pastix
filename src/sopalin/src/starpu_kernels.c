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
#ifdef WITH_STARPU

#ifdef WITH_MAGMABLAS
#include <magmablas.h>
#endif /* WITH_MAGMABLAS */

#ifdef STARPU_USE_DEPRECATED_API
#undef STARPU_USE_DEPRECATED_API
#endif
#include <starpu.h>
#include <starpu_cuda.h>
#include "common_pastix.h"
#include "starpu_kernels.h"
#include "starpu_submit_tasks.h"

#include "sparse_gemm.h"
#include "pastix_cuda_helper.h"

#include <inttypes.h>

#ifdef STARPU_USE_CUDA
#if ((!defined PREC_DOUBLE)  || (!(defined __CUDA_ARCH__) || __CUDA_ARCH__ >= 130))
#if !(defined PREC_DOUBLE && defined TYPE_COMPLEX && CUDA_SM_VERSION < 20)
#ifndef FORCE_NO_CUDA
#define STARPU_USE_CUDA_GEMM_FUNC
#endif
#endif
#endif
#endif



#if (!defined STARPU_USE_CUDA || defined FORCE_NO_CUDA)
#  ifdef  WITH_MAGMABLAS
#    undef  WITH_MAGMABLAS
#  endif /* WITH_MAGMABLAS */
#endif /* !STARPU_USE_CUDA || FORCE_NO_CUDA*/


#ifdef WITH_MAGMABLAS
#include "geadd_cuda.h"
#include "getra_cuda.h"
#include "zgetrf_stapiv_gpu.h"
#endif /* WITH_MAGMABLAS */

#define DimTrans                API_CALL(DimTrans)
#define kernel_trsm             API_CALL(kernel_trsm)
#define CORE_gemdm              API_CALL(CORE_gemdm)

#if !(defined STARPU_USE_CUDA_GEMM_FUNC)
#define CUDA_SPARSE_GEMM(TRANSA, TRANSB,                \
                         dimi, dimj, dima,              \
                         alpha,                         \
                         A,  stride_A,                  \
                         B,  stride_B,                  \
                         beta,                          \
                         C, stride_C,                   \
                         blocknbr,  blocktab,           \
                         fblocknbr, fblocktab)          \
  do {                                                  \
  } while(0)
#define CUDA_SPARSE_GEMDM(TRANSA, TRANSB,                \
                          dimi, dimj, dima,              \
                          alpha,                         \
                          A,  stride_A,                  \
                          D,  stride_D,                  \
                          B,  stride_B,                  \
                          beta,                          \
                          C, stride_C,                   \
                          blocknbr,  blocktab,           \
                          fblocknbr, fblocktab)          \
  do {                                                   \
    assert(0);                                           \
    /* avoid warnings */                                 \
    assert(fblocktab);                                   \
    assert(blocktab);                                    \
    assert(fblocknbr);                                   \
    assert(blocknbr);                                    \
    assert(beta);                                        \
    assert(alpha);                                       \
    assert(dima);                                        \
    assert(dimj);                                        \
    assert(dimi);                                        \
    assert(stridefc);                                    \
  } while(0)
#else
#if (CUDA_SM_VERSION >= 20)
#define CUDA_SPARSE_GEMM(TRANSA, TRANSB,                                \
                         dimi, dimj, dima,                              \
                         alpha,                                         \
                         A,  stride_A,                                  \
                         B,  stride_B,                                  \
                         beta,                                          \
                         C, stride_C,                                   \
                         blocknbr, blocktab,                            \
                         fblocknbr, fblocktab)                          \
  do {                                                                  \
    CU_FLOAT cu_alpha = CU_FLOAT_INIT(creal(alpha), cimag(alpha));      \
    CU_FLOAT cu_beta  = CU_FLOAT_INIT(creal(beta),  cimag(beta));       \
    GENERATE_SM_VERSION_NAME(gemm)(*TRANSA, *TRANSB,                    \
                                   (int)dimi, (int)dimj, (int)dima,     \
                                   cu_alpha,                            \
                                   (CU_FLOAT*)A, (int)stride_A,         \
                                   (CU_FLOAT*)B, (int)stride_B,         \
                                   cu_beta,                             \
                                   (CU_FLOAT*)C, (int)stride_C,         \
                                   blocknbr,  blocktab,                 \
                                   fblocknbr, fblocktab,                \
                                   starpu_cuda_get_local_stream());     \
  } while(0)
#define CUDA_SPARSE_GEMDM(TRANSA, TRANSB,                               \
                          dimi, dimj, dima,                             \
                          alpha,                                        \
                          A,  stride_A,                                 \
                          D, stride_D,                                  \
                          B,  stride_B,                                 \
                          beta,                                         \
                          C, stride_C,                                  \
                          blocknbr, blocktab,                           \
                          fblocknbr, fblocktab)                         \
    do {                                                                \
      CU_FLOAT cu_alpha = CU_FLOAT_INIT(creal(alpha), cimag(alpha));    \
      CU_FLOAT cu_beta  = CU_FLOAT_INIT(creal(beta),  cimag(beta));     \
      GENERATE_SM_VERSION_NAME(gemdm)(*TRANSA, *TRANSB,                 \
                                      (int)dimi, (int)dimj, (int)dima,  \
                                      cu_alpha,                         \
                                      (CU_FLOAT*)A, (int)stride_A,      \
                                      (CU_FLOAT*)D, (int)stride_D,      \
                                      (CU_FLOAT*)B, (int)stride_B,      \
                                      cu_beta,                          \
                                      (CU_FLOAT*)C, (int)stride_C,      \
                                      blocknbr,  blocktab,              \
                                      fblocknbr, fblocktab,             \
                                      starpu_cuda_get_local_stream());  \
    } while(0)
#else
#define CUDA_SPARSE_GEMM(TRANSA, TRANSB,                                \
                         dimi, dimj, dima,                              \
                         alpha,                                         \
                         A,  stride_A,                                  \
                         B,  stride_B,                                  \
                         beta,                                          \
                         C, stride_C,                                   \
                         blocknbr, blocktab,                            \
                         fblocknbr, fblocktab)                          \
  do {                                                                  \
    magmablas_sparse_gemm_kernel_N_T_64_16_4_16_4((int)dimi,            \
                                                  (int)dimj,            \
                                                  (int)dima,            \
                                                  (float)alpha,         \
                                                  A,                    \
                                                  (int)stride_A,        \
                                                  B,                    \
                                                  (int)stride_B,        \
                                                  (float)beta,          \
                                                  C,                    \
                                                  (int)stride_C,        \
                                                  blocknbr,             \
                                                  blocktab,             \
                                                  fblocknbr,            \
                                                  fblocktab);           \
  } while(0)
#define CUDA_SPARSE_GEMDM(TRANSA, TRANSB,                               \
                          dimi, dimj, dima,                             \
                          alpha,                                        \
                          A,  stride_A,                                 \
                          D,  stride_D,                                 \
                          B,  stride_B,                                 \
                          beta,                                         \
                          C, stride_C,                                  \
                          blocknbr, blocktab,                           \
                          fblocknbr, fblocktab)                         \
  do {                                                                  \
    magmablas_sparse_gemdm_kernel_N_T_64_16_4_16_4((int)dimi,           \
                                                   (int)dimj,           \
                                                   (int)dima,           \
                                                   (float)alpha,        \
                                                   A,                   \
                                                   (int)stride_A,       \
                                                   D,                   \
                                                   (int)stride_D,       \
                                                   B,                   \
                                                   (int)stride_B,       \
                                                   (float)beta,         \
                                                   C,                   \
                                                   (int)stride_C,       \
                                                   blocknbr,            \
                                                   blocktab,            \
                                                   fblocknbr,           \
                                                   fblocktab);          \
  } while(0)
#endif
#endif

#define SUBMIT_TRF_IF_NEEDED						\
  {									\
    PASTIX_INT STIN_t, STIN_j;							\
    PASTIX_INT STIN_fblknum      = SYMB_BLOKNUM(cblknum);			\
    PASTIX_INT STIN_lblknum      = SYMB_BLOKNUM(cblknum+1);			\
    PASTIX_INT STIN_n = bloknum - STIN_fblknum;				\
									\
    STIN_n = (STIN_n * (STIN_n - 1)) / 2;				\
    STIN_n = (bloknum - STIN_fblknum - 1) *				\
      (STIN_lblknum - STIN_fblknum) - STIN_n;				\
									\
    for (STIN_j = bloknum ; STIN_j < STIN_lblknum; STIN_j++)		\
      {									\
	STIN_t = SOLV_INDTAB[TASK_INDNUM(args->tasknum)+(STIN_n++)];	\
	if (STIN_t < 0)							\
	  {								\
	    TASK_CTRBCNT(-STIN_t)--;					\
	    if (TASK_CTRBCNT(-STIN_t) == 0)				\
	      {								\
		starpu_submit_one_trf(-STIN_t, sopalin_data);		\
	      }								\
	  }								\
      }									\
  }


/* #define SOPALIN_SPARSE_GEMM(TRANSA, TRANSB,                     \ */
/*                             dimi, dimj, dima,                   \ */
/*                             alpha,  A,  stride_A,               \ */
/*                             B,  stride_B,                       \ */
/*                             beta,  C, stride_C,                 \ */
/*                             nblocs, blocs_idx, facing_bloc_idx, \ */
/*                             wtmp, wtmpsize)                     \ */
/*   do {                                                          \ */
/*     sparse_gemm(TRANSA, TRANSB,                                 \ */
/*                 dimi, dimj, dima,                               \ */
/*                 alpha,  A,  stride_A,                           \ */
/*                 B,  stride_B,                                   \ */
/*                 beta,  C, stride_C,                             \ */
/*                 nblocs, blocs_idx, facing_bloc_idx,             \ */
/*                 wtmp, wtmpsize);                                \ */
/*   } while(0) */

#define SOPALIN_SPARSE_GEMM(TRANSA, TRANSB,                 \
                            dimi, dimj, dima,               \
                            alpha,  A,  stride_A,           \
                            B,  stride_B,                   \
                            beta,  C, stride_C,             \
                            blocknbr, blocktab,             \
                            fblocknbr, fblocktab,           \
                            wtmp, wtmpsize)                 \
  do {                                                      \
    assert(0);                                              \
    assert(wtmp); assert(wtmpsize);                         \
    /* Avoid warnings */                                    \
    assert(fblocktab);                                      \
    assert(blocktab );                                      \
    assert(fblocknbr);                                      \
    assert(blocknbr );                                      \
    assert(beta     );                                      \
    assert(alpha    );                                      \
    assert(dima     );                                      \
    assert(dimj     );                                      \
    assert(dimi     );                                      \
    assert(stridefc );                                      \
  } while(0)

#ifdef TYPE_COMPLEX
#  ifdef PREC_DOUBLE
#    define CUDA_TRSM(RL, UL, transa, transb, dimb, dima,        \
                      one_cuf,                                   \
                      A, stridea,                                \
                      B, strideb) do {                           \
      magma_ztrsm(RL, UL, transa, transb, dimb, dima, one_cuf,   \
                  (CU_FLOAT*)((void*)A), stridea,                \
                  (CU_FLOAT*)((void*)B), strideb);               \
    } while(0)
#  else /* not PREC_DOUBLE */
#    define CUDA_TRSM(RL, UL, transa, transb, dimb, dima,        \
                      one_cuf,                                   \
                      A, stridea,                                \
                      B, strideb) do {                           \
      magma_ctrsm(RL, UL, transa, transb, dimb, dima, one_cuf,   \
                  (CU_FLOAT*)((void*)A), stridea,                \
                  (CU_FLOAT*)((void*)B), strideb);               \
    } while(0)
#  endif /* not PREC_DOUBLE */
#else /* not TYPE_COMPLEX */
#  ifdef PREC_DOUBLE
#    define CUDA_TRSM(RL, UL, transa, transb, dimb, dima,        \
                      one_cuf,                                   \
                      A, stridea,                                \
                      B, strideb) do {                           \
      magma_dtrsm(RL, UL, transa, transb, dimb, dima, one_cuf,   \
                  (CU_FLOAT*)((void*)A), stridea,                \
                  (CU_FLOAT*)((void*)B), strideb);               \
    } while(0)
#  else /* not PREC_DOUBLE */
#    define CUDA_TRSM(RL, UL, transa, transb, dimb, dima,        \
                      one_cuf,                                   \
                      A, stridea,                                \
                      B, strideb) do {                           \
      magma_strsm(RL, UL, transa, transb, dimb, dima, one_cuf,   \
                  (CU_FLOAT*)((void*)A), stridea,                \
                  (CU_FLOAT*)((void*)B), strideb);               \
    } while(0)
#  endif /* not PREC_DOUBLE */
#endif /* not TYPE_COMPLEX */

#ifdef CHOL_SOPALIN
#ifdef SOPALIN_LU
/* General LU decomposition */
/*
 * Function: getrfsp1d_starpu_common
 *
 * Diagonal block factorization and column block update for LU decomposition.
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     0            - L column block
 *     1            - U column block
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
static inline
void getrfsp1d_starpu_common(void * buffers[], void * _args, int arch)
{
  starpu_trf_data_t * args         = (starpu_trf_data_t*)_args;
  Sopalin_Data_t    * sopalin_data = args->sopalin_data;
  SolverMatrix      * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT             * lDiag        = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT             * uDiag        = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
  PASTIX_INT                 stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                 cblknum      = args->cblknum;
  PASTIX_INT                 fblknum      = SYMB_BLOKNUM(cblknum);
  PASTIX_INT                 lblknum      = SYMB_BLOKNUM(cblknum+1);
  PASTIX_FLOAT             * lExtraDiag   = NULL;
  PASTIX_FLOAT             * uExtraDiag   = NULL;
  PASTIX_INT                 dima         = SYMB_LCOLNUM(cblknum) -
    SYMB_FCOLNUM(cblknum) + 1;
  PASTIX_INT                 dimb         = stride - dima;
  int                 me           = starpu_worker_get_id();
#ifdef WITH_MAGMABLAS
  CU_FLOAT            one_cuf      = CU_FLOAT_INIT(1.0, 0.0);
  magma_int_t         info;
#endif /* WITH_MAGMABLAS */

  /* check if diagonal column block */
  assert( SYMB_FCOLNUM(cblknum) == SYMB_FROWNUM(fblknum) );

  switch(arch) {
  case ARCH_CPU:
    /* Add U diagonal updates into L */
    SOPALIN_GEAM("T", "N", dima, dima, 1.0,
                 uDiag, stride,
                 lDiag, stride);

    /* Factorize diagonal block (two terms version with workspace) */
    PASTIX_getrf_block(lDiag, dima, dima, stride,
                       &(sopalin_data->thread_data[me]->nbpivot),
                       sopalin_data->critere);
    /* Transpose L_diag in U_diag Matrix */
    DimTrans(lDiag,stride, dima,uDiag);
    break;
#ifdef WITH_MAGMABLAS
  case ARCH_CUDA:
    geadd_cuda("T", "N", dima, dima,
               1.0, uDiag, stride,
               1.0, lDiag, stride);
    magma_zgetrf_stapiv_gpu(dima, dima, (CU_FLOAT*)lDiag, stride,
                            sopalin_data->critere,
                            &(sopalin_data->thread_data[me]->nbpivot), &info);
    getra_cuda(lDiag, stride, uDiag, stride, dima);
    break;
#endif /* WITH_MAGMABLAS */
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }

  /* if there is an extra-diagonal bloc in column block */
  if ( fblknum+1 < lblknum )
    {
      lExtraDiag = lDiag + SOLV_COEFIND(fblknum+1);
      uExtraDiag = uDiag + SOLV_COEFIND(fblknum+1);

      switch(arch) {
      case ARCH_CPU:
        kernel_trsm(dimb, dima,
                    lDiag,       uDiag,      stride,
                    lExtraDiag,  uExtraDiag, stride);
        break;
#ifdef WITH_MAGMABLAS
      case ARCH_CUDA:
        CUDA_TRSM('R', 'U', 'N', 'N', dimb, dima,
                  one_cuf, lDiag, stride, lExtraDiag, stride);
        CUDA_TRSM('R', 'U', 'N', 'U', dimb, dima,
                  one_cuf, uDiag, stride, uExtraDiag, stride);
        break;
#endif
      default:
        errorPrint("Unknown Architecture");
        assert(0);
        break;
      }
    }
#ifdef STARPU_SUBMIT_READY
  starpu_submit_bunch_of_gemm(args->tasknum, sopalin_data);
#endif
}



/*
 * Function: getrfsp1d_starpu_cpu
 *
 * Diagonal block factorization and column block update for LU decomposition.
 *
 * CPU function.
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     0            - L column block
 *     1            - U column block
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
void getrfsp1d_starpu_cpu(void * buffers[], void * _args)
{
  getrfsp1d_starpu_common(buffers, _args, ARCH_CPU);
}

#ifdef WITH_MAGMABLAS
/*
 * Function: getrfsp1d_starpu_cuda
 *
 * Diagonal block factorization and column block update for LU decomposition.
 *
 * CUDA function.
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     0            - L column block
 *     1            - U column block
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
void getrfsp1d_starpu_cuda(void * buffers[], void * _args)
{
  getrfsp1d_starpu_common(buffers, _args, ARCH_CUDA);
}
#endif /* WITH_MAGMABLAS */

/*
  Function: getrfsp1d_gemm_starpu_common

  General LU update of left block column facing current block.

  Common function for CPU and GPU.

  Update block by block.

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - U column block.
      3 - U facing column block.
      4 - Working memory area.
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
   arch     - indicate if the codelet is runned on CPU or CUDA node.
*/
static inline void
getrfsp1d_gemm_starpu_common(void * buffers[], void * _args, int arch)
{
  starpu_gemm_data_t * args         = (starpu_gemm_data_t*)_args;
  Sopalin_Data_t     * sopalin_data = args->sopalin_data;
  PASTIX_FLOAT              * L            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT              * Cl           = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
  PASTIX_FLOAT              * U            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[2]);
  PASTIX_FLOAT              * Cu           = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[3]);
  PASTIX_FLOAT              * work         = (PASTIX_FLOAT*)STARPU_VECTOR_GET_PTR(buffers[4]);
  PASTIX_INT                  stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                  stridefc     = STARPU_MATRIX_GET_LD(buffers[1]);
  PASTIX_INT                  cblknum      = args->cblknum;
  PASTIX_INT                  bloknum      = args->bloknum;
  PASTIX_INT                  fcblknum     = args->fcblknum;
  PASTIX_FLOAT              * wtmp;
  PASTIX_FLOAT              * C;
  PASTIX_FLOAT              * Aik;
  PASTIX_FLOAT              * Akj;
  PASTIX_FLOAT              * Aij;
  SolverMatrix       *datacode = sopalin_data->datacode;
  PASTIX_INT fblknum, lblknum, frownum;
  PASTIX_INT indblok;
  PASTIX_INT b, j;
  PASTIX_INT dimi, dimj, dima, dimb;

  fblknum = SYMB_BLOKNUM(cblknum);
  lblknum = SYMB_BLOKNUM(cblknum + 1);

  indblok = SOLV_COEFIND(bloknum);

  dimi = stride - indblok;
  dimj = SYMB_LROWNUM(bloknum) - SYMB_FROWNUM(bloknum) + 1;
  dima = STARPU_MATRIX_GET_NY(buffers[0]);

  /* Matrix A = Aik */
  Aik = L + indblok;
  Akj = U + indblok;
  wtmp = work;

  switch(arch) {
  case ARCH_CPU:
    /* fprintf(stdout, */
    /*         "L  %016"   PRIXPTR " U  %016"  PRIXPTR */
    /*         " Cl  %016" PRIXPTR " Cu  %016" PRIXPTR "\n", L, U, Cl, Cu); */

    SOPALIN_GEMM( "N", "T",
                  dimi, dimj, dima,
                  1.,  Aik,  stride,
                  Akj,  stride,
                  0.,  wtmp, dimi);
    break;
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }

  /*
   * Add contribution to facing cblk
   */
  b = SYMB_BLOKNUM( fcblknum );
  C = Cl + (SYMB_FROWNUM(bloknum) - SYMB_FCOLNUM(fcblknum)) * stridefc;
  /* for all following blocks in block column */
  for (j=bloknum; j<lblknum; j++) {
    frownum = SYMB_FROWNUM(j);

    /* Find facing bloknum */
    while (!BLOCK_ISFACING(j,b))
      {
        b++;
        assert( b < SYMB_BLOKNUM( fcblknum+1 ) );
      }

    Aij = C + SOLV_COEFIND(b) + frownum - SYMB_FROWNUM(b);
    dimb = SYMB_LROWNUM(j) - frownum + 1;

    switch(arch) {
    case ARCH_CPU:
      SOPALIN_GEAM("N", "N", dimb, dimj, -1.0,
                   wtmp, dimi,
                   Aij,  stridefc );
      break;
    default:
      errorPrint("Unknown Architecture");
      assert(0);
      break;
    }


    /* Displacement to next block */
    wtmp += dimb;
  }

  /*
   * Compute update on U
   */

  Aik = U + indblok;
  Akj = L + indblok;
  wtmp = work;
  switch(arch) {
  case ARCH_CPU:
    SOPALIN_GEMM( "N", "T",
                  dimi, dimj, dima,
                  1.,  Aik,  stride,
                  Akj,  stride,
                  0.,  wtmp, dimi  );
    break;
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }

  wtmp += SYMB_LROWNUM(bloknum) - SYMB_FROWNUM(bloknum) + 1;

  /*
   * Add contribution to facing cblk
   */
  b = SYMB_BLOKNUM( fcblknum );
  C = Cl + (SYMB_FROWNUM(bloknum) - SYMB_FCOLNUM(fcblknum));

  /* for all following blocks in block column */
  for (j=bloknum+1; j<lblknum; j++) {
    frownum = SYMB_FROWNUM(j);

    /* Find facing bloknum */
    /* WARNING: may not work for NAPA */
    if (
#ifdef NAPA_SOPALIN /* ILU(k) */
        !(((SYMB_FROWNUM(j)>=SYMB_FROWNUM(b)) &&
           (SYMB_LROWNUM(j)<=SYMB_LROWNUM(b))) ||
          ((SYMB_FROWNUM(j)<=SYMB_FROWNUM(b)) &&
           (SYMB_LROWNUM(j)>=SYMB_LROWNUM(b))) ||
          ((SYMB_FROWNUM(j)<=SYMB_FROWNUM(b)) &&
           (SYMB_LROWNUM(j)>=SYMB_FROWNUM(b))) ||
          ((SYMB_FROWNUM(j)<=SYMB_LROWNUM(b)) &&
           (SYMB_LROWNUM(j)>=SYMB_LROWNUM(b))))
#else
        !((SYMB_FROWNUM(j)>=SYMB_FROWNUM(b)) &&
          (SYMB_LROWNUM(j)<=SYMB_LROWNUM(b)))
#endif
        )
      break;


    Aij = C + (frownum - SYMB_FROWNUM(b))*stridefc;
    dimb = SYMB_LROWNUM(j) - frownum + 1;

    switch(arch) {
    case ARCH_CPU:
      SOPALIN_GEAM( "T", "N", dimj, dimb, -1.0,
                    wtmp, dimi,
                    Aij,  stridefc );
      break;
    default:
      errorPrint("Unknown Architecture");
      assert(0);
      break;
    }

    /* Displacement to next block */
    wtmp += dimb;
  }


  C = Cu + (SYMB_FROWNUM(bloknum) - SYMB_FCOLNUM(fcblknum)) * stridefc;

  /* Keep updating on U */
  for (; j<lblknum; j++) {
    frownum = SYMB_FROWNUM(j);

    /* Find facing bloknum */
    while (!BLOCK_ISFACING(j,b))
      {
        b++;
        assert( b < SYMB_BLOKNUM( fcblknum+1 ) );
      }

    dimb = SYMB_LROWNUM(j) - frownum + 1;
    Aij = C + SOLV_COEFIND(b) + frownum - SYMB_FROWNUM(b);
    switch(arch) {
    case ARCH_CPU:
      SOPALIN_GEAM("N", "N", dimb, dimj, -1.0,
                   wtmp, dimi,
                   Aij,  stridefc );
      break;
    default:
      errorPrint("Unknown Architecture");
      assert(0);
      break;
    }


    /* Displacement to next block */
    wtmp += dimb;
  }

#ifdef STARPU_SUBMIT_READY
  SUBMIT_TRF_IF_NEEDED;
#endif

}

/*
  Function: getrfsp1d_gemm_starpu_cpu

  General LU update of left block column facing current block.

  Update block by block.

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - U column block.
      3 - U facing column block.
      4 - Working memory area.
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
*/
void getrfsp1d_gemm_starpu_cpu(void * buffers[], void * _args)
{
  getrfsp1d_gemm_starpu_common(buffers, _args, ARCH_CPU);
}

#ifdef STARPU_USE_CUDA
/*
  Function: getrfsp1d_gemm_starpu_cuda

  General LU update of left block column facing current block.

  Update block by block.

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - U column block.
      3 - U facing column block.
      4 - Working memory area.
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
*/
void getrfsp1d_gemm_starpu_cuda(void * buffers[], void * _args)
{
  getrfsp1d_gemm_starpu_common(buffers, _args, ARCH_CUDA);
}
#endif


/*
  Function: getrfsp1d_sparse_gemm_starpu_common

  Sparse LU update of left block column facing current block.

  Common function for CPU and GPU.

  Update all the facing column block at once.

  TODO: Implement the CPU version

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - U column block.
      3 - U facing column block.
      4 - Working memory area.
      5 - blocktab (depending on -DSTARPU_BLOCKTAB_SELFCOPY).
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
      nblocs       - Number of blocks in current column block.
      d_blocktab   - pointers on blocktabs for CUDA nodes
                     (-DSTARPU_BLOCKTAB_SELFCOPY).
   arch     - indicate if the codelet is runned on CPU or CUDA node.
*/
static inline
void getrfsp1d_sparse_gemm_starpu_common(void * buffers[], void * _args,
                                         int arch)
{
  starpu_gemm_data_t * args         = (starpu_gemm_data_t*)_args;
  Sopalin_Data_t     * sopalin_data = args->sopalin_data;
  SolverMatrix       * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT              * L            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT              * Cl           = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
  PASTIX_FLOAT              * U            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[2]);
  PASTIX_FLOAT              * Cu           = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[3]);
  PASTIX_FLOAT              * work         = (PASTIX_FLOAT*)STARPU_VECTOR_GET_PTR(buffers[4]);
  size_t               worksize     = STARPU_VECTOR_GET_NX(buffers[4]);
  PASTIX_INT                  stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                  stridefc     = STARPU_MATRIX_GET_LD(buffers[1]);
#ifdef STARPU_BLOCKTAB_SELFCOPY
  int                  device;
  int                * all_blocktab;
#else
  int                * all_blocktab = (int*)STARPU_VECTOR_GET_PTR(buffers[5]);
#endif
  PASTIX_INT                  cblknum      = args->cblknum;
  PASTIX_INT                  bloknum      = args->bloknum;
  PASTIX_INT                  fcblknum     = args->fcblknum;
  PASTIX_INT                  indblok      = SOLV_COEFIND(bloknum);
  PASTIX_INT                  dimi         = stride - indblok;
  PASTIX_INT                  dimj         = BLOK_ROWNBR(bloknum);
  PASTIX_INT                  dima         = CBLK_COLNBR(cblknum);
  PASTIX_FLOAT                alpha        = -1.;
  PASTIX_FLOAT                beta         = 1.;
  PASTIX_FLOAT              * Aik;
  PASTIX_FLOAT              * Akj;
  PASTIX_FLOAT              * Aij;
  int                  blocknbr  = SYMB_BLOKNUM(cblknum+1) - bloknum;
  int                  fblocknbr = CBLK_BLOKNBR(fcblknum);
#ifdef STARPU_BLOCKTAB_SELFCOPY
  int                * blocktab;
  int                * fblocktab;
  if (arch == ARCH_CUDA) {
    cudaGetDevice(&device);
    all_blocktab=args->d_blocktab[device];
    blocktab   = &(all_blocktab[2*bloknum]);
    fblocktab  = &(all_blocktab[2*SYMB_BLOKNUM(fcblknum)]);
  }
#else
  int                * blocktab  = &(all_blocktab[2*bloknum]);
  int                * fblocktab = &(all_blocktab[2*SYMB_BLOKNUM(fcblknum)]);
#endif

  Aik = L + indblok;
  Akj = U + indblok;
  Aij = Cl + ( SYMB_FROWNUM(bloknum) -
               SYMB_FCOLNUM(fcblknum) )*SOLV_STRIDE(fcblknum);

  switch(arch) {
  case ARCH_CPU:
    SOPALIN_SPARSE_GEMM( "N", "T",
                         dimi, dimj, dima,
                         alpha,  Aik,  stride,
                         Akj,  stride,
                         beta,  Aij, stridefc,
                         blocknbr, blocktab, fblocknbr, fblocktab,
                         work, worksize);
    break;
#ifdef STARPU_USE_CUDA
  case ARCH_CUDA:
    CUDA_SPARSE_GEMM( "n", "t",
                      dimi, dimj, dima,
                      alpha,  Aik,  stride,
                      Akj,  stride,
                      beta,  Aij, stridefc,
                      blocknbr, blocktab, fblocknbr, fblocktab);
    break;
#endif
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }

  /*
   * Compute update on U
   */
  if ( blocknbr > 1 )
    {
      dimi = dimi - (SOLV_COEFIND(bloknum+1) - SOLV_COEFIND(bloknum));
      Aik = U + SOLV_COEFIND(bloknum+1);
      Akj = L + indblok;
      Aij = Cu + (SYMB_FROWNUM(bloknum) -
                  SYMB_FCOLNUM(fcblknum))*SOLV_STRIDE(fcblknum);
      switch(arch) {
      case ARCH_CPU:
        SOPALIN_SPARSE_GEMM( "N", "T",
                             dimi, dimj, dima,
                             alpha,  Aik,  stride,
                             Akj,  stride,
                             beta,  Aij, stridefc,
                             blocknbr-1, &(blocktab[2]),
                             fblocknbr, fblocktab,
                             work, worksize);
        break;
#ifdef STARPU_USE_CUDA
      case ARCH_CUDA:
        CUDA_SPARSE_GEMM( "n", "t",
                          dimi, dimj, dima,
                          alpha,  Aik,  stride,
                          Akj,  stride,
                          beta,  Aij, stridefc,
                          blocknbr-1, &(blocktab[2]),
                          fblocknbr, fblocktab);
        break;
#endif
      default:
        errorPrint("Unknown Architecture");
        assert(0);
        break;
      }
    }

#ifdef STARPU_SUBMIT_READY
  SUBMIT_TRF_IF_NEEDED;
#endif

}


/*
  Function: getrfsp1d_sparse_gemm_starpu_cpu

  Sparse LU update of left block column facing current block.

  Update all the facing column block at once.

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - U column block.
      3 - U facing column block.
      4 - Working memory area.
      5 - blocktab (depending on -DSTARPU_BLOCKTAB_SELFCOPY).
    _args   - codelet arguments
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
      nblocs       - Number of blocks in current column block.
      d_blocktab   - pointers on blocktabs for CUDA nodes
                     (-DSTARPU_BLOCKTAB_SELFCOPY).
*/
void getrfsp1d_sparse_gemm_starpu_cpu(void * buffers[], void * _args)
{
  getrfsp1d_gemm_starpu_common(buffers, _args, ARCH_CPU);
}

/*
  Function: getrfsp1d_sparse_gemm_starpu_cuda

  Sparse LU update of left block column facing current block.

  Update all the facing column block at once.

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - U column block.
      3 - U facing column block.
      4 - Working memory area.
      5 - blocktab (depending on -DSTARPU_BLOCKTAB_SELFCOPY).
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
      nblocs       - Number of blocks in current column block.
      d_blocktab   - pointers on blocktabs for CUDA nodes
                     (-DSTARPU_BLOCKTAB_SELFCOPY).
*/
void getrfsp1d_sparse_gemm_starpu_cuda(void * buffers[], void * _args)
{
  getrfsp1d_sparse_gemm_starpu_common(buffers, _args, ARCH_CUDA);
}

#else /* not SOPALIN_LU */
/* LLT decomposition */

/*
 * Function: potrfsp1d_starpu_common
 *
 * Diagonal block factorization and column block update for LLt decomposition.
 *
 * Parameters:
 *   buffers - Data handlers :
 *     0 - L column block
 *   _args   - Codelet arguments:
 *    sopalin_data - global PaStiX internal data.
 *    cblknum      - Current column block index.
 */
static inline void
potrfsp1d_starpu_common(void * buffers[], void * _args, int arch)
{
  starpu_trf_data_t * args         = (starpu_trf_data_t*)_args;
  Sopalin_Data_t    * sopalin_data = args->sopalin_data;
  SolverMatrix      * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT             * Diag         = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_INT                 stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                 cblknum      = args->cblknum;
  PASTIX_INT                 fblknum      = SYMB_BLOKNUM(cblknum);
  PASTIX_INT                 lblknum      = SYMB_BLOKNUM(cblknum+1);
  PASTIX_INT                 dima         = CBLK_COLNBR(cblknum);
  PASTIX_INT                 dimb         = stride - dima;
  PASTIX_FLOAT             * ExtraDiag    = NULL;
  int                 me           = starpu_worker_get_id();
#ifdef WITH_MAGMABLAS
  int                 info;
  CU_FLOAT            one_cuf      = CU_FLOAT_INIT(1.0, 0.0);
#endif /* WITH_MAGMABLAS */

  /* check if diagonal column block */
  assert( SYMB_FCOLNUM(cblknum) == SYMB_FROWNUM(fblknum) );

  /* Factorize diagonal block (two terms version with workspace) */
  switch(arch)
    {
    case ARCH_CPU:
      PASTIX_potrf_block(Diag, dima, stride,
                         &(sopalin_data->thread_data[me]->nbpivot),
                         sopalin_data->critere);
      break;
#ifdef WITH_MAGMABLAS
    case ARCH_CUDA:
      magma_zpotrf_stapiv_gpu("L", dima, (CU_FLOAT*)Diag, stride,
                              sopalin_data->critere,
                              &(sopalin_data->thread_data[me]->nbpivot), &info);
    break;
#endif /* WITH_MAGMABLAS */
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }
  /* if there is an extra-diagonal bloc in column block */
  if ( fblknum+1 < lblknum )
    {
      ExtraDiag    = Diag + SOLV_COEFIND(fblknum+1);
      switch(arch)
        {
        case ARCH_CPU:
          kernel_trsm(dimb, dima,
                      Diag,      stride,
                      ExtraDiag, stride);
          break;
#ifdef WITH_MAGMABLAS
        case ARCH_CUDA:
          CUDA_TRSM('R', 'L', 'T', 'N', dimb, dima,
                    one_cuf, Diag, stride, ExtraDiag, stride);
          break;
#endif /* WITH_MAGMABLAS */
        default:
          errorPrint("Unknown Architecture");
          assert(0);
          break;
        }
    }
#ifdef STARPU_SUBMIT_READY
  starpu_submit_bunch_of_gemm(args->tasknum, sopalin_data);
#endif
}

/*
 * Function: potrfsp1d_starpu_cpu
 *
 * Diagonal block factorization and column block update for LLt decomposition.
 *
 * Parameters:
 *   buffers - Data handlers :
 *     0 - L column block
 *   _args   - Codelet arguments:
 *    sopalin_data - global PaStiX internal data.
 *    cblknum      - Current column block index.
 */
void
potrfsp1d_starpu_cpu(void * buffers[], void * _args)
{
  potrfsp1d_starpu_common(buffers, _args, ARCH_CPU);
}

/*
 * Function: potrfsp1d_starpu_cuda
 *
 * Diagonal block factorization and column block update for LLt decomposition.
 *
 * Parameters:
 *   buffers - Data handlers :
 *     0 - L column block
 *   _args   - Codelet arguments:
 *    sopalin_data - global PaStiX internal data.
 *    cblknum      - Current column block index.
 */
void
potrfsp1d_starpu_cuda(void * buffers[], void * _args)
{
  potrfsp1d_starpu_common(buffers, _args, ARCH_CUDA);
}

/*
  Function: potrfsp1d_gemm_starpu_cpu

  General LLt update of left block column facing current block.

  Update block by block.

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - Working memory area.
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
*/
void potrfsp1d_gemm_starpu_cpu(void * buffers[], void * _args)
{
  starpu_gemm_data_t * args         = (starpu_gemm_data_t*)_args;
  Sopalin_Data_t     * sopalin_data = args->sopalin_data;
  SolverMatrix       * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT              * L            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT              * C            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
#ifndef LLT_NO_SCRATCH
  PASTIX_FLOAT              * work         = (PASTIX_FLOAT*)STARPU_VECTOR_GET_PTR(buffers[2]);
#else
  PASTIX_FLOAT              * Akj;
#endif
  PASTIX_INT                  stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                  stridefc     = STARPU_MATRIX_GET_LD(buffers[1]);
  PASTIX_INT                  cblknum      = args->cblknum;
  PASTIX_INT                  bloknum      = args->bloknum;
  PASTIX_INT                  fcblknum     = args->fcblknum;
  PASTIX_INT                  lblknum      = SYMB_BLOKNUM(cblknum+1);
  PASTIX_INT                  indblok      = SOLV_COEFIND(bloknum);
  PASTIX_INT                  dimi         = stride - indblok;
  PASTIX_INT                  dimj         = BLOK_ROWNBR(bloknum);
  PASTIX_INT                  dima         = CBLK_COLNBR(cblknum);

  PASTIX_FLOAT *Aik, *Aij;
  PASTIX_INT frownum;
  PASTIX_INT b, j;
  PASTIX_INT dimb;


  /* Matrix A = Aik */
  Aik = L + indblok;
#ifndef LLT_NO_SCRATCH
  /* Compute the contribution */
  SOPALIN_GEMM( "N", "C",
                dimi, dimj, dima,
                1.,  Aik,  stride,
                Aik,  stride,
                0.,  work, dimi);
#else
  Akj = Aij;
#endif
  /*
   * Add contribution to facing cblk
   */
  b = SYMB_BLOKNUM( fcblknum );
  C = C + (SYMB_FROWNUM(bloknum) - SYMB_FCOLNUM(fcblknum)) * stridefc;
  /* for all following blocks in block column */
  for (j=bloknum; j<lblknum; j++) {
    frownum = SYMB_FROWNUM(j);

    /* Find facing bloknum */
    while(!BLOCK_ISFACING(j,b))
      {
        b++;
        assert( b < SYMB_BLOKNUM( fcblknum+1 ) );
      }

    Aij = C + SOLV_COEFIND(b) + frownum - SYMB_FROWNUM(b);
    dimb = SYMB_LROWNUM(j) - frownum + 1;
#ifdef LLT_NO_SCRATCH
  /* Compute the contribution */
  SOPALIN_GEMM( "N", "C",
                dimi, dimj, dima,
                -1.0.,  Aik,  stride,
                Akj,  stride,
                1.0,  Aij, stridefc);
  Akj += dimb;
#else
    MAT_zaxpy( dimb, dimj, -1.0,
               work, dimi,
               Aij,  stridefc );

    /* Displacement to next block */
    work += dimb;
#endif
  }

#ifdef STARPU_SUBMIT_READY
  SUBMIT_TRF_IF_NEEDED;
#endif
}

/*
  Function: potrfsp1d_sparse_gemm_starpu_common

  Sparse LLt update of left block column facing current block.

  Common function for CPU and GPU.

  Update all the facing column block at once.

  TODO: Implement the CPU version

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - Working memory area.
      3 - blocktab (depending on -DSTARPU_BLOCKTAB_SELFCOPY).
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
      nblocs       - Number of blocks in current column block.
      d_blocktab   - pointers on blocktabs for CUDA nodes
                     (-DSTARPU_BLOCKTAB_SELFCOPY).
    arch    - indicate if the codelet is runned on CPU or CUDA node.
*/
static inline
void potrfsp1d_sparse_gemm_starpu_common(void * buffers[], void * _args,
                                         int arch)
{
  starpu_gemm_data_t * args         = (starpu_gemm_data_t*)_args;
  Sopalin_Data_t     * sopalin_data = args->sopalin_data;
  SolverMatrix       * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT              * L            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT              * C            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
  PASTIX_FLOAT              * work         = (PASTIX_FLOAT*)STARPU_VECTOR_GET_PTR(buffers[2]);
  size_t               worksize     = STARPU_VECTOR_GET_NX(buffers[2]);
  PASTIX_INT                  stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                  stridefc     = STARPU_MATRIX_GET_LD(buffers[1]);
#ifdef STARPU_BLOCKTAB_SELFCOPY
  int                  device;
  int                * all_blocktab;
#else
  int                * all_blocktab = (int*)STARPU_VECTOR_GET_PTR(buffers[3]);
#endif
  PASTIX_INT                  cblknum      = args->cblknum;
  PASTIX_INT                  bloknum      = args->bloknum;
  PASTIX_INT                  fcblknum     = args->fcblknum;
  PASTIX_INT                  indblok      = SOLV_COEFIND(bloknum);
  PASTIX_INT                  dimi         = stride - indblok;
  PASTIX_INT                  dimj         = BLOK_ROWNBR(bloknum);
  PASTIX_INT                  dima         = CBLK_COLNBR(cblknum);
  PASTIX_FLOAT                alpha = -1.;
  PASTIX_FLOAT                beta = 1.;
  PASTIX_FLOAT              * Aik;
  PASTIX_FLOAT              * Aij;
  int                  blocknbr     = SYMB_BLOKNUM(cblknum+1) - bloknum;
  int                  fblocknbr    = CBLK_BLOKNBR(fcblknum);
#ifdef STARPU_BLOCKTAB_SELFCOPY
  int                * blocktab;
  int                * fblocktab;
  cudaGetDevice(&device);
  all_blocktab=args->d_blocktab[device];
  blocktab   = &(all_blocktab[2*bloknum]);
  fblocktab  = &(all_blocktab[2*SYMB_BLOKNUM(fcblknum)]);
#else
  int                * blocktab     = &(all_blocktab[2*bloknum]);
  int                * fblocktab    = &(all_blocktab[2*SYMB_BLOKNUM(fcblknum)]);
#endif

  Aik = L + indblok;
  Aij = C + ( SYMB_FROWNUM(bloknum) -
              SYMB_FCOLNUM(fcblknum) )*SOLV_STRIDE(fcblknum);

  /* Compute the contribution */
  switch(arch) {
  case ARCH_CPU:
    SOPALIN_SPARSE_GEMM( "N", "C",
                         dimi, dimj, dima,
                         alpha,  Aik,  stride,
                         Aik,  stride,
                         beta,  Aij, stridefc,
                         blocknbr, blocktab, fblocknbr, fblocktab,
                         work, worksize);
    break;
#ifdef STARPU_USE_CUDA
  case ARCH_CUDA:
    CUDA_SPARSE_GEMM( "n", "c",
                      dimi, dimj, dima,
                      alpha,  Aik,  stride,
                      Aik,  stride,
                      beta,  Aij, stridefc,
                      blocknbr, blocktab, fblocknbr, fblocktab);
    break;
#endif
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }

#ifdef STARPU_SUBMIT_READY
  SUBMIT_TRF_IF_NEEDED;
#endif
}

/*
  Function: potrfsp1d_sparse_gemm_starpu_cpu

  Sparse LLt update of left block column facing current block.

  Update all the facing column block at once.

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - Working memory area.
      3 - blocktab (depending on -DSTARPU_BLOCKTAB_SELFCOPY).
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
      nblocs       - Number of blocks in current column block.
      d_blocktab   - pointers on blocktabs for CUDA nodes
                     (-DSTARPU_BLOCKTAB_SELFCOPY).
*/
void potrfsp1d_sparse_gemm_starpu_cpu(void * buffers[], void * _args)
{
  potrfsp1d_gemm_starpu_cpu(buffers, _args);
}

/*
  Function: potrfsp1d_sparse_gemm_starpu_cuda

  Sparse LLt update of left block column facing current block.

  Update all the facing column block at once.

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - Working memory area.
      3 - blocktab (depending on -DSTARPU_BLOCKTAB_SELFCOPY).
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
      nblocs       - Number of blocks in current column block.
      d_blocktab   - pointers on blocktabs for CUDA nodes
                     (-DSTARPU_BLOCKTAB_SELFCOPY).
*/
void potrfsp1d_sparse_gemm_starpu_cuda(void * buffers[], void * _args)
{
  potrfsp1d_sparse_gemm_starpu_common(buffers, _args, ARCH_CUDA);
}

#endif /* not SOPALIN_LU */
#else  /* not CHOL_SOPALIN */
/* LDLT Decomposition */
/*
 * Function: hetrfsp1d_starpu_common
 *
 * Diagonal block factorization and column block update for LDLt decomposition.
 *
 * Parameters:
 *   buffers - Data handlers :
 *     0 - L column block
 *     1 - Will receive L*DIAG(BLOCK_DIAG(L))
 *   _args   - Codelet arguments:
 *    sopalin_data - global PaStiX internal data.
 *    cblknum      - Current column block index.
 */
static inline
void hetrfsp1d_starpu_common(void * buffers[], void * _args, int arch)
{
  starpu_trf_data_t * args         = (starpu_trf_data_t*)_args;
  Sopalin_Data_t    * sopalin_data = args->sopalin_data;
  SolverMatrix      * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT             * Diag         = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT             * tmp4         = (PASTIX_FLOAT*)STARPU_VECTOR_GET_PTR(buffers[1]);
  PASTIX_INT                 stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                 cblknum      = args->cblknum;
  PASTIX_INT                 fblknum      = SYMB_BLOKNUM(cblknum);
  PASTIX_INT                 lblknum      = SYMB_BLOKNUM(cblknum+1);
  PASTIX_FLOAT             * ExtraDiag    = NULL;
  PASTIX_INT                 dima         = CBLK_COLNBR(cblknum);
  PASTIX_INT                 dimb         = stride - dima;
  int                 me           = starpu_worker_get_id();
#ifdef WITH_MAGMABLAS
  CU_FLOAT            one_cuf      = CU_FLOAT_INIT(1.0, 0.0);
  int                 info;
#endif /* WITH_MAGMABLAS */

  /* check if diagonal column block */
  assert( SYMB_FCOLNUM(cblknum) == SYMB_FROWNUM(fblknum) );

  /* Initialisation des pointeurs de blocs */

  /* Factorize diagonal block (two terms version with workspace) */
  switch(arch)
    {
    case ARCH_CPU:
#ifdef HERMITIAN
      PASTIX_hetrf_block(Diag, dima, stride,
                         &(sopalin_data->thread_data[me]->nbpivot),
                         sopalin_data->critere,
                         sopalin_data->thread_data[me]->maxbloktab1);
#else
      PASTIX_sytrf_block(Diag, dima, stride,
                         &(sopalin_data->thread_data[me]->nbpivot),
                         sopalin_data->critere,
                         sopalin_data->thread_data[me]->maxbloktab1);
#endif
      break;
#ifdef WITH_MAGMABLAS
    case ARCH_CUDA:
#ifdef HERMITIAN
      magma_zhetrf_stapiv_gpu('L', dima,
                              Diag, stride,
                              sopalin_data->critere,
                              &(sopalin_data->thread_data[me]->nbpivot),
                              sopalin_data->thread_data[me]->maxbloktab1,
                              info);
#else
      magma_zsytrf_stapiv_gpu('L', dima,
                              Diag, stride,
                              sopalin_data->critere,
                              &(sopalin_data->thread_data[me]->nbpivot),
                              sopalin_data->thread_data[me]->maxbloktab1,
                              info);
#endif
      break;
#endif /* WITH_MAGMABLAS */
    default:
      errorPrint("Unknown Architecture");
      assert(0);
      break;
    }

  /* if there is an extra-diagonal bloc in column block */
  if ( fblknum+1 < lblknum )
    {
      ExtraDiag = Diag + SOLV_COEFIND(fblknum+1);
      /* WARNING : Second "tmp4" should be a buffer to use for send. */
      switch(arch)
        {
        case ARCH_CPU:
          kernel_trsm(dimb, dima,
                      Diag,      stride,
                      ExtraDiag, tmp4, stride);
          break;
#ifdef WITH_MAGMABLAS
        case ARCH_CUDA:
#  ifdef HERMITIAN
          CUDA_TRSM('R', 'L', 'C', 'N', dimb, dima,
                    one_cuf, Diag, stride, ExtraDiag, stride);
#  else /* not HERMITIAN */
          CUDA_TRSM('R', 'L', 'T', 'N', dimb, dima,
                    one_cuf, Diag, stride, ExtraDiag, stride);
#  endif /* not HERMITIAN */
          /* TODO */
/*             for (k=0; k<n; k++) */
/*   { */
/*     PASTIX_FLOAT alpha; */
/* # ifdef COMPUTE */
/*     ASSERTDBG(dL[k+k*ldd] != 0., MOD_SOPALIN); */
/*     alpha = fun / dL[k+k*ldd]; */
/* # endif */
/*     SOPALIN_COPY(m, &(L[ k*ldl]), iun, */
/*             &(L2[k*m  ]), iun); */
/*     SOPALIN_SCAL(m, alpha, &(L[k*ldl]), iun); */
/*   } */

/*           break; */
#endif /* WITH_MAGMABLAS */
        default:
          errorPrint("Unknown Architecture");
          assert(0);
          break;
        }

    }
#ifdef STARPU_SUBMIT_READY
  starpu_submit_bunch_of_gemm(args->tasknum, sopalin_data);
#endif
}

/*
  Function: hetrfsp1d_starpu_cpu

  Diagonal block factorization and column block update for LDLt decomposition.

  Parameters:
    buffers - Data handlers :
      0 - L column block
      1 - Will receive L*DIAG(BLOCK_DIAG(L))
    _args   - Codelet arguments:
     sopalin_data - global PaStiX internal data.
     cblknum      - Current column block index.
*/
void hetrfsp1d_starpu_cpu(void * buffers[], void * _args)
{
  hetrfsp1d_starpu_common(buffers, _args, ARCH_CPU);
}

/*
  Function: hetrfsp1d_starpu_cuda

  Diagonal block factorization and column block update for LDLt decomposition.

  Parameters:
    buffers - Data handlers :
      0 - L column block
      1 - Will receive L*DIAG(BLOCK_DIAG(L))
    _args   - Codelet arguments:
     sopalin_data - global PaStiX internal data.
     cblknum      - Current column block index.
*/
void hetrfsp1d_starpu_cuda(void * buffers[], void * _args)
{
  hetrfsp1d_starpu_common(buffers, _args, ARCH_CUDA);
}

/*
  Function: hetrfsp1d_gemm_starpu_cpu

  General LDLt update of left block column facing current block.

  CPU function.

  Update block by block.

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - Working memory area.
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
*/
void hetrfsp1d_gemm_starpu_cpu(void * buffers[], void * _args)
{
  starpu_gemm_data_t * args         = (starpu_gemm_data_t*)_args;
  Sopalin_Data_t     * sopalin_data = args->sopalin_data;
  SolverMatrix       * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT              * L            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT              * C            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
  PASTIX_FLOAT              * work1        = (PASTIX_FLOAT*)STARPU_VECTOR_GET_PTR(buffers[2]);
  PASTIX_FLOAT              * work2        = work1 + SOLV_COEFMAX;
  PASTIX_INT                  cblknum      = args->cblknum;
  PASTIX_INT                  bloknum      = args->bloknum;
  PASTIX_INT                  fcblknum     = args->fcblknum;
  PASTIX_FLOAT *Aik, *Aij;
  PASTIX_INT fblknum, lblknum, frownum;
  PASTIX_INT stride, stridefc, indblok;
  PASTIX_INT b, j;
  PASTIX_INT dimi, dimj, dima, dimb;
  PASTIX_INT ldw = SOLV_COEFMAX;
  fblknum = SYMB_BLOKNUM(cblknum);
  lblknum = SYMB_BLOKNUM(cblknum + 1);

  indblok = SOLV_COEFIND(bloknum);
  stride  = SOLV_STRIDE(cblknum);

  dimi = stride - indblok;
  dimj = SYMB_LROWNUM(bloknum) - SYMB_FROWNUM(bloknum) + 1;
  dima = SYMB_LCOLNUM(cblknum) - SYMB_FCOLNUM(cblknum) + 1;

  /* Matrix A = Aik */
  Aik = L + indblok;

  /* Compute the contribution */
  CORE_gemdm( PastixNoTrans,
#ifdef HERMITIAN
              PastixConjTrans,
#else
              PastixTrans,
#endif
              dimi, dimj, dima,
              1.,  Aik,   stride,
              Aik,   stride,
              0.,  work1, dimi,
              L,     stride+1,
              work2, ldw );

  /*
   * Add contribution to facing cblk
   */
  b = SYMB_BLOKNUM( fcblknum );
  stridefc = SOLV_STRIDE(fcblknum);
  C = C + (SYMB_FROWNUM(bloknum) - SYMB_FCOLNUM(fcblknum)) * stridefc;

  /* for all following blocks in block column */
  for (j=bloknum; j<lblknum; j++) {
    frownum = SYMB_FROWNUM(j);

    /* Find facing bloknum */
    while (!BLOCK_ISFACING(j,b))
      {
        b++;
        assert( b < SYMB_BLOKNUM( fcblknum+1 ) );
      }


    Aij = C + SOLV_COEFIND(b) + frownum - SYMB_FROWNUM(b);
    dimb = SYMB_LROWNUM(j) - frownum + 1;

    MAT_zaxpy( dimb, dimj, -1.0,
               work1, dimi,
               Aij,   stridefc );

    /* Displacement to next block */
    work1 += dimb;
  }

#ifdef STARPU_SUBMIT_READY
  SUBMIT_TRF_IF_NEEDED;
#endif
}

/*
  Function: hetrfsp1d_gemm_starpu_cuda

  General LDLt update of left block column facing current block.

  Cuda function.

  Update block by block.

  Parameters:
    buffers - Data handlers :
      0 - L column block.
      1 - L facing column block.
      2 - blocktab (depending on -DSTARPU_BLOCKTAB_SELFCOPY).
    _args   - codelet arguments :
      sopalin_data - global PaStiX internal data.
      cblknum      - Current column block index.
      bloknum      - Current block index.
      fcblknum     - Facing column block index.
*/
void hetrfsp1d_gemm_starpu_cuda(void * buffers[], void * _args)
{
  starpu_gemm_data_t * args         = (starpu_gemm_data_t*)_args;
  Sopalin_Data_t     * sopalin_data = args->sopalin_data;
  SolverMatrix       * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT              * L            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT              * C            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
  PASTIX_INT                  stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                  stridefc     = STARPU_MATRIX_GET_LD(buffers[1]);
#ifdef STARPU_BLOCKTAB_SELFCOPY
  int                  device;
  int                * all_blocktab;
#else
  int                * all_blocktab = (int*)STARPU_VECTOR_GET_PTR(buffers[3]);
#endif
  PASTIX_INT                  cblknum      = args->cblknum;
  PASTIX_INT                  bloknum      = args->bloknum;
  PASTIX_INT                  fcblknum     = args->fcblknum;
  PASTIX_INT                  indblok      = SOLV_COEFIND(bloknum);
  PASTIX_INT                  dimi         = stride - indblok;
  PASTIX_INT                  dimj         = BLOK_ROWNBR(bloknum);
  PASTIX_INT                  dima         = CBLK_COLNBR(cblknum);
  PASTIX_FLOAT                alpha = -1.;
  PASTIX_FLOAT                beta = 1.;
  PASTIX_FLOAT              * Aik;
  PASTIX_FLOAT              * Aij;
  int                  blocknbr     = SYMB_BLOKNUM(cblknum+1) - bloknum;
  int                  fblocknbr    = CBLK_BLOKNBR(fcblknum);
#ifdef STARPU_BLOCKTAB_SELFCOPY
  int                * blocktab;
  int                * fblocktab;
  cudaGetDevice(&device);
  all_blocktab=args->d_blocktab[device];
  blocktab   = &(all_blocktab[2*bloknum]);
  fblocktab  = &(all_blocktab[2*SYMB_BLOKNUM(fcblknum)]);
#else
  int                * blocktab     = &(all_blocktab[2*bloknum]);
  int                * fblocktab    = &(all_blocktab[2*SYMB_BLOKNUM(fcblknum)]);
#endif

  /* Matrix A = Aik */
  Aik = L + indblok;
  Aij = C + (SYMB_FROWNUM(bloknum) -
             SYMB_FCOLNUM(fcblknum))*SOLV_STRIDE(fcblknum);

  /* Compute the contribution */
#ifdef HERMITIAN
  CUDA_SPARSE_GEMDM( "n",
                     "c",
                     dimi, dimj, dima,
                     alpha,
                     Aik,  stride,
                     L,    stride+1,
                     Aik,  stride,
                     beta,
                     Aij,  stridefc,
                     blocknbr, blocktab,
                     fblocknbr, fblocktab);
#else
  CUDA_SPARSE_GEMDM( "n",
                     "t",
                     dimi, dimj, dima,
                     alpha,
                     Aik,  stride,
                     L,    stride+1,
                     Aik,  stride,
                     beta,
                     Aij,  stridefc,
                     blocknbr, blocktab,
                     fblocknbr, fblocktab);
#endif

#ifdef STARPU_SUBMIT_READY
  SUBMIT_TRF_IF_NEEDED;
#endif
}
#endif /* not CHOL_SOPALIN */

#else  /* not WITH_STARPU */
/* ISO C forbids an empty source file */
#include "not_empty.h"
NOT_EMPTY( starpu_kernels )
#endif /* not WITH_STARPU */
