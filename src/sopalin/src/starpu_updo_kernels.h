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
#ifndef STARPU_UPDO_KERNELS_H
#define STARPU_UPDO_KERNELS_H
#include "sopalin_define.h"
struct starpu_updo_trsm_data_ {
  PASTIX_INT              cblknum;
  Sopalin_Data_t * sopalin_data;
  char             transpose;
  char             diag;
};
typedef struct starpu_updo_trsm_data_ starpu_updo_trsm_data_t;

struct starpu_updo_gemm_data_ {
  PASTIX_INT              cblknum;
  PASTIX_INT              bloknum;
  Sopalin_Data_t * sopalin_data;
  char             transpose;
};
typedef struct starpu_updo_gemm_data_ starpu_updo_gemm_data_t;

struct starpu_updo_diag_trsm_data_ {
  PASTIX_INT              cblknum;
  Sopalin_Data_t * sopalin_data;
};
typedef struct starpu_updo_diag_trsm_data_ starpu_updo_diag_data_t;

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
#define updo_trsm_starpu_cpu API_CALL(updo_trsm_starpu_cpu)
void updo_trsm_starpu_cpu(void * buffers[], void * _args);
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
#define updo_down_gemm_starpu_cpu API_CALL(updo_down_gemm_starpu_cpu)
void updo_down_gemm_starpu_cpu(void * buffers[], void * _args);

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
#define updo_up_gemm_starpu_cpu API_CALL(updo_up_gemm_starpu_cpu)
void updo_up_gemm_starpu_cpu(void * buffers[], void * _args);

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
#define updo_diag_starpu_cpu API_CALL(updo_diag_starpu_cpu)
void updo_diag_starpu_cpu(void * buffers[], void * args);

#endif /* STARPU_UPDO_KERNELS_H */
