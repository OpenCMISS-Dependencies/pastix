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
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011
*/

#ifndef _CGEMDM_FERMI_DEFINE_RIGHT_H_
#define _CGEMDM_FERMI_DEFINE_RIGHT_H_

#define PRECISION_c

//#include "gemdm_stencil_defs.cu"
#include "gemdm_stencil_right.h"
///////////////////////////////////////////////////////////////////////////////////////////////////
// Common parameters

// size of work for a thread block
/* #define BLK_M_nn 64 */
/* #define BLK_N_nn 64 */
/* #define BLK_K_nn 16 */

/* #define BLK_M_nt 64 */
/* #define BLK_N_nt 64 */
/* #define BLK_K_nt 16 */

/* #define BLK_M_tt 64 */
/* #define BLK_N_tt 64 */
/* #define BLK_K_tt 16 */

/* #define BLK_M_tn 64 */
/* #define BLK_N_tn 64 */
/* #define BLK_K_tn 16 */

/* // size of work for a thread block */
/* #define BLK_M BLK_M_nn */
/* #define BLK_N BLK_N_nn  */
/* #define BLK_K BLK_K_nn */

/* // size of thread block for calculating C (innermost loop) */
/* #define DIM_X 16 */
/* #define DIM_Y 16 */

/* /////////////////////////////////////////////////////////////////////////////////////////////////// */
/* // */
/* //  NoTrans - NoTrans */
/* // */

/* // size of thread block for reading A (dev->regs->shmem) */
/* #define DIM_XA 32 */
/* #define DIM_YA  8 */
  
/* // size of thread block for reading B (dev->regs->shmem) */
/* #define DIM_XB 16 */
/* #define DIM_YB 16 */

/* #define version trans_nn */
/* #include "gemdm_stencil.cu" */
 
/* #undef DIM_XA */
/* #undef DIM_YA */

/* #undef DIM_XB */
/* #undef DIM_YB */

/* /////////////////////////////////////////////////////////////////////////////////////////////////// */
/* // */
/* //  NoTrans - Trans */
/* // */
 
/* // size of thread block for reading A (dev->regs->shmem) */
/* #define DIM_XA 16 */
/* #define DIM_YA 16 */

/* // size of thread block for reading B (dev->regs->shmem) */
/* #define DIM_XB 16 */
/* #define DIM_YB 16 */

/* #define version trans_nt */
/* #include "gemdm_stencil.cu" */

/* #define version trans_nc */
/* #include "gemdm_stencil.cu" */

/* #undef DIM_XA */
/* #undef DIM_YA */

/* #undef DIM_XB */
/* #undef DIM_YB */

/* /////////////////////////////////////////////////////////////////////////////////////////////////// */
/* // */
/* //  Trans - Trans */
/* // */
 
/* // size of thread block for reading A (dev->regs->shmem) */
/* #define DIM_XA 16 */
/* #define DIM_YA 16 */

/* // size of thread block for reading B (dev->regs->shmem) */
/* #define DIM_XB 32 */
/* #define DIM_YB  8 */

/* #define version trans_tt */
/* #include "gemdm_stencil.cu" */

/* #define version trans_cc */
/* #include "gemdm_stencil.cu" */

/* #define version trans_ct */
/* #include "gemdm_stencil.cu" */

/* #define version trans_tc */
/* #include "gemdm_stencil.cu" */

/* #undef DIM_XA */
/* #undef DIM_YA */

/* #undef DIM_XB */
/* #undef DIM_YB */

/* /////////////////////////////////////////////////////////////////////////////////////////////////// */
/* // */
/* //  Trans - NoTrans */
/* // */
 
/* // size of thread block for reading A (dev->regs->shmem) */
/* #define DIM_XA 16 */
/* #define DIM_YA 16 */

/* // size of thread block for reading B (dev->regs->shmem) */
/* #define DIM_XB 16 */
/* #define DIM_YB 16 */

/* #define version trans_tn */
/* #include "gemdm_stencil.cu" */

/* #define version trans_cn */
/* #include "gemdm_stencil.cu" */

#endif /* _CGEMDM_FERMI_DEFINE_RIGHT_H_ */
