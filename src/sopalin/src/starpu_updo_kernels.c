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

#  ifdef WITH_MAGMABLAS
#    include <magmablas.h>
#  endif /* WITH_MAGMABLAS */


#  include <starpu.h>
#  include "common_pastix.h"
#  include "starpu_updo_kernels.h"
#  include <inttypes.h>



/*
 * Function: updo_trsm_starpu_common
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
static inline
void updo_trsm_starpu_common(void * buffers[], void * _args, int arch)
{
  starpu_updo_trsm_data_t * args         = (starpu_updo_trsm_data_t*)_args;
  Sopalin_Data_t          * sopalin_data = args->sopalin_data;
  SolverMatrix            * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT                   * L            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT                   * RHS          = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
  PASTIX_INT                       stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                       rhsnbr       = STARPU_MATRIX_GET_NY(buffers[1]);
  PASTIX_INT                       rhssze       = STARPU_MATRIX_GET_LD(buffers[1]);
  PASTIX_INT                       cblknum      = args->cblknum;
  char                    * transpose    = &(args->transpose);
  char                    * diag         = &(args->diag);
  PASTIX_INT                       colnbr       = CBLK_COLNBR(cblknum);
  PASTIX_FLOAT                     fun          = 1.0;

  ASSERTDBG(UPDOWN_SM2XNBR == rhsnbr, MOD_SOPALIN);
  ASSERTDBG(UPDOWN_SM2XSZE == rhssze, MOD_SOPALIN);
  switch(arch) {
  case ARCH_CPU:
    SOPALIN_TRSM("L","L",transpose,diag,colnbr,rhsnbr,fun,L,stride,RHS,rhssze);
    break;
  case ARCH_CUDA:
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }
}



/*
 * Function: updo_trsm_starpu_cpu
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
void updo_trsm_starpu_cpu(void * buffers[], void * _args)
{
  updo_trsm_starpu_common(buffers, _args, ARCH_CPU);
}


/*
 * Function: updo_down_gemm_starpu_common
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
static inline
void updo_down_gemm_starpu_common(void * buffers[], void * _args, int arch)
{
  starpu_updo_gemm_data_t * args         = (starpu_updo_gemm_data_t*)_args;
  Sopalin_Data_t          * sopalin_data = args->sopalin_data;
  SolverMatrix            * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT                   * L            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT                   * RHS          = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
  PASTIX_FLOAT                   * RHS2         = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[2]);
  PASTIX_INT                       stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                       rhsnbr       = STARPU_MATRIX_GET_NY(buffers[1]);
  PASTIX_INT                       rhssze       = STARPU_MATRIX_GET_LD(buffers[1]);
  PASTIX_INT                       cblknum      = args->cblknum;
  PASTIX_INT                       bloknum      = args->bloknum;
  char                    * transpose    = &(args->transpose);
  PASTIX_INT                       fcblknum     = SYMB_CBLKNUM(bloknum);
  PASTIX_INT                       colnbr       = CBLK_COLNBR(cblknum);
  PASTIX_INT                       rownbr       = BLOK_ROWNBR(bloknum);
  PASTIX_FLOAT                     fun          = 1.0;
  PASTIX_FLOAT                   * ga           = L + SOLV_COEFIND(bloknum);
  PASTIX_FLOAT                   * gc           = RHS2 +
    SYMB_FROWNUM(bloknum) - SYMB_FCOLNUM(fcblknum);

  ASSERTDBG(UPDOWN_SM2XNBR == rhsnbr, MOD_SOPALIN);
  ASSERTDBG(UPDOWN_SM2XSZE == rhssze, MOD_SOPALIN);

  switch(arch) {
  case ARCH_CPU:
    SOPALIN_GEMM(transpose,"N",rownbr,rhsnbr,colnbr,-fun,ga,stride,
                 RHS,rhssze,fun,gc, UPDOWN_SM2XSZE);
    break;
  case ARCH_CUDA:
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }
}



/*
 * Function: updo_down_gemm_starpu_cpu
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
void updo_down_gemm_starpu_cpu(void * buffers[], void * _args)
{
  updo_down_gemm_starpu_common(buffers, _args, ARCH_CPU);
}


/*
 * Function: updo_up_gemm_starpu_common
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
static inline
void updo_up_gemm_starpu_common(void * buffers[], void * _args, int arch)
{
  starpu_updo_gemm_data_t * args         = (starpu_updo_gemm_data_t*)_args;
  Sopalin_Data_t          * sopalin_data = args->sopalin_data;
  SolverMatrix            * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT                   * L            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT                   * RHS          = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
  PASTIX_FLOAT                   * RHS2         = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[2]);
  PASTIX_INT                       stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                       rhsnbr       = STARPU_MATRIX_GET_NY(buffers[1]);
  PASTIX_INT                       rhssze       = STARPU_MATRIX_GET_LD(buffers[1]);
  PASTIX_INT                       cblknum      = args->cblknum;
  PASTIX_INT                       bloknum      = args->bloknum;
  char                    * transpose    = &(args->transpose);
  PASTIX_INT                       fcblknum     = SYMB_CBLKNUM(bloknum);
  PASTIX_INT                       colnbr       = CBLK_COLNBR(cblknum);
  PASTIX_INT                       rownbr       = BLOK_ROWNBR(bloknum);
  PASTIX_FLOAT                     fun          = 1.0;
  PASTIX_FLOAT                   * ga           = L + SOLV_COEFIND(bloknum);
  PASTIX_FLOAT                   * gc           = RHS2 +
    SYMB_FROWNUM(bloknum) - SYMB_FCOLNUM(fcblknum);

  ASSERTDBG(UPDOWN_SM2XNBR == rhsnbr, MOD_SOPALIN);
  ASSERTDBG(UPDOWN_SM2XSZE == rhssze, MOD_SOPALIN);

  switch(arch) {
  case ARCH_CPU:
    SOPALIN_GEMM(transpose,"N",colnbr,rhsnbr,rownbr,-fun,ga,stride,
                 gc,rhssze,fun,RHS, UPDOWN_SM2XSZE);
    break;
  case ARCH_CUDA:
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }
}



/*
 * Function: updo_up_gemm_starpu_cpu
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
void updo_up_gemm_starpu_cpu(void * buffers[], void * _args)
{
  updo_up_gemm_starpu_common(buffers, _args, ARCH_CPU);
}

/*
 * Function: updo_diag_starpu_common
 *
 * Divide the right-hand-side(s) by the diagonal
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     0            - L column block
 *     1            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 *   arch       - Type of architecture : ARCH_CPU | ARCH_CUDA
 */
static inline
void updo_diag_starpu_common(void * buffers[], void * _args, int arch)
{
  starpu_updo_diag_data_t * args         = (starpu_updo_diag_data_t*)_args;
  Sopalin_Data_t          * sopalin_data = args->sopalin_data;
  SolverMatrix            * datacode     = sopalin_data->datacode;
  PASTIX_FLOAT                   * L            = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[0]);
  PASTIX_FLOAT                   * RHS          = (PASTIX_FLOAT*)STARPU_MATRIX_GET_PTR(buffers[1]);
  PASTIX_INT                       stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  PASTIX_INT                       rhsnbr       = STARPU_MATRIX_GET_NY(buffers[1]);
  PASTIX_INT                       rhssze       = STARPU_MATRIX_GET_LD(buffers[1]);
  PASTIX_INT                       cblknum      = args->cblknum;
  PASTIX_INT                       colnbr       = CBLK_COLNBR(cblknum);

  ASSERTDBG(UPDOWN_SM2XNBR == rhsnbr, MOD_SOPALIN);
  ASSERTDBG(UPDOWN_SM2XSZE == rhssze, MOD_SOPALIN);

  switch(arch) {
  case ARCH_CPU:
  {
    PASTIX_INT i, j;
    PASTIX_FLOAT * myRHS = RHS;
    for (j = 0; j < rhsnbr; j++)
      {
        for (i = 0; i < colnbr; i++)
          {
            myRHS[i] /= L[i*(stride+1)];
          }
        myRHS += rhssze;
      }
    break;
  }
  case ARCH_CUDA:
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }
}

/*
 * Function: updo_diag_starpu_cpu
 *
 * Divide the right-hand-side(s) by the diagonal.
 *
 * CPU interface.
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     0            - L column block
 *     1            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
void updo_diag_starpu_cpu(void * buffers[], void * args)
{
  updo_diag_starpu_common(buffers, args, ARCH_CPU);
}
#endif /* WITH_STARPU */
