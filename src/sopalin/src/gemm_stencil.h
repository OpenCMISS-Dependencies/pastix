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

#ifndef _GEMM_STENCIL_H_
#define _GEMM_STENCIL_H_

///////////////////////////////////////////////////////////////////////////////////////////////////
// Common parameters
///////////////////////////////////////////////////////////////////////////////////////////////////

//#define TEXTURE_1D

///////////////////////////////////////////////////////////////////////////////////////////////////

#define trans_nn 1
#define trans_nt 2
#define trans_nc 3

#define trans_tn 4
#define trans_tt 5
#define trans_tc 6

#define trans_cn 7
#define trans_ct 8
#define trans_cc 9

///////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(PRECISION_z)
    typedef cuDoubleComplex FloatingPoint_t;
#elif defined(PRECISION_c)
    typedef cuFloatComplex FloatingPoint_t;
#elif defined(PRECISION_d)
    typedef double FloatingPoint_t;
#elif defined(PRECISION_s)
    typedef float FloatingPoint_t;
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef TEXTURE_1D
#define fetch(A, m, n) tex_fetch(tex_ref_##A, coord_##A + n*LD##A+m)
#else
#define fetch(A, m, n) offs_d##A[n*LD##A+m]
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
#if defined(PRECISION_z)

#define conj(A)          cuConj(A)
#define add(A, B)        cuCadd(A, B)
#define mul(A, B)        cuCmul(A, B)
#define fma(A, B, C) C = cuCfma(A, B, C)
#define make_FloatingPoint(x, y) make_cuDoubleComplex(x, y);

#if defined(TEXTURE_1D)

static __device__
FloatingPoint_t tex_fetch(texture<int4> tex_ref, int coord)
{
    int4 v = tex1Dfetch(tex_ref, coord);
    return make_cuDoubleComplex(__hiloint2double(v.y, v.x), __hiloint2double(v.w, v.z));
}

texture<int4, 1, cudaReadModeElementType> tex_ref_A;
texture<int4, 1, cudaReadModeElementType> tex_ref_B;

#endif/* defined(TEXTURE_1D) */

///////////////////////////////////////////////////////////////////////////////////////////////////
#elif defined(PRECISION_c)

#define conj(A)          cuConjf(A)
#define add(A, B)        cuCaddf(A, B)
#define mul(A, B)        cuCmulf(A, B)
#define fma(A, B, C) C = cuCfmaf(A, B, C)
#define make_FloatingPoint(x, y) make_cuFloatComplex(x, y);

#if defined(TEXTURE_1D)

static __device__
FloatingPoint_t tex_fetch(texture<float2> tex_ref, int coord)
{
    return tex1Dfetch(tex_ref, coord);
}

texture<float2, 1, cudaReadModeElementType> tex_ref_A;
texture<float2, 1, cudaReadModeElementType> tex_ref_B;

#endif/* defined(TEXTURE_1D) */

///////////////////////////////////////////////////////////////////////////////////////////////////
#elif defined(PRECISION_d)

#define conj(A)           (A)
#define add(A, B)         (A+B)
#define mul(A, B)         (A*B)
#define fma(A, B, C) C += (A*B)
#define make_FloatingPoint(x, y) (x)

#if defined(TEXTURE_1D)

static __device__
FloatingPoint_t tex_fetch(texture<int2> tex_ref, int coord)
{
    int2 v = tex1Dfetch(tex_ref, coord);
    return __hiloint2double(v.y, v.x);
}

texture<int2, 1, cudaReadModeElementType> tex_ref_A;
texture<int2, 1, cudaReadModeElementType> tex_ref_B;

#endif/* defined(TEXTURE_1D) */

///////////////////////////////////////////////////////////////////////////////////////////////////
#elif defined(PRECISION_s)

#define conj(A)           (A)
#define add(A, B)         (A+B)
#define mul(A, B)         (A*B)
#define fma(A, B, C) C += (A*B)
#define make_FloatingPoint(x, y) (x)

#if defined(TEXTURE_1D)

static __device__
FloatingPoint_t tex_fetch(texture<float> tex_ref, int coord)
{
    return tex1Dfetch(tex_ref, coord);
}

texture<float, 1, cudaReadModeElementType> tex_ref_A;
texture<float, 1, cudaReadModeElementType> tex_ref_B;

#endif/* defined(TEXTURE_1D) */

///////////////////////////////////////////////////////////////////////////////////////////////////
#endif /* defined(PRECISION_x) */


///////////////////////////////////////////////////////////////////////////////////////////////////
//  Block sizes parameters
///////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(PRECISION_z)

#define DIM_X  8
#define DIM_Y  8

#define BLK_M_nn 24
#define BLK_N_nn 16
#define BLK_K_nn  8

#define DIM_XA_nn 8
#define DIM_YA_nn 8

#define DIM_XB_nn 8
#define DIM_YB_nn 8

///////////////////
#define BLK_M_nt 16
#define BLK_N_nt 24
#define BLK_K_nt  8

#define DIM_XA_nt 8
#define DIM_YA_nt 8

#define DIM_XB_nt 8
#define DIM_YB_nt 8

///////////////////
#define BLK_M_tt 16
#define BLK_N_tt 24
#define BLK_K_tt  8

#define DIM_XA_tt  4
#define DIM_YA_tt 16

#define DIM_XB_tt 8
#define DIM_YB_tt 8

///////////////////
#define BLK_M_tn 24
#define BLK_N_tn 16
#define BLK_K_tn  8

#define DIM_XA_tn 8
#define DIM_YA_tn 8

#define DIM_XB_tn 8
#define DIM_YB_tn 8

////////////////////////////////////////////////////////////////////////////
#elif defined(PRECISION_c)

#define DIM_X 16
#define DIM_Y 16

//////// NN ///////
#define BLK_M_nn 64
#define BLK_N_nn 64
#define BLK_K_nn 16

#define DIM_XA_nn 32
#define DIM_YA_nn  8

#define DIM_XB_nn 16
#define DIM_YB_nn 16

//////// NT ///////
#define BLK_M_nt 64
#define BLK_N_nt 64
#define BLK_K_nt 16

#define DIM_XA_nt 16
#define DIM_YA_nt 16

#define DIM_XB_nt 16
#define DIM_YB_nt 16

//////// TT ///////
#define BLK_M_tt 64
#define BLK_N_tt 64
#define BLK_K_tt 16

#define DIM_XA_tt 16
#define DIM_YA_tt 16

#define DIM_XB_tt 32
#define DIM_YB_tt  8

//////// TN ///////
#define BLK_M_tn 64
#define BLK_N_tn 64
#define BLK_K_tn 16

#define DIM_XA_tn 16
#define DIM_YA_tn 16

#define DIM_XB_tn 16
#define DIM_YB_tn 16

////////////////////////////////////////////////////////////////////////////
#elif defined(PRECISION_d)

#define DIM_X 16
#define DIM_Y 16

//////// NN ///////
#define BLK_M_nn 64
#define BLK_N_nn 64
#define BLK_K_nn 16

#define DIM_XA_nn 16
#define DIM_YA_nn 16

#define DIM_XB_nn 16
#define DIM_YB_nn 16

//////// NT ///////
#define BLK_M_nt 64
#define BLK_N_nt 64
#define BLK_K_nt 16

#define DIM_XA_nt 16
#define DIM_YA_nt 16

#define DIM_XB_nt 16
#define DIM_YB_nt 16

//////// TT ///////
#define BLK_M_tt 64
#define BLK_N_tt 64
#define BLK_K_tt 16

#define DIM_XA_tt 16
#define DIM_YA_tt 16

#define DIM_XB_tt 16
#define DIM_YB_tt 16

//////// TN ///////
#define BLK_M_tn 64
#define BLK_N_tn 64
#define BLK_K_tn 16

#define DIM_XA_tn 16
#define DIM_YA_tn 16

#define DIM_XB_tn 16
#define DIM_YB_tn 16

////////////////////////////////////////////////////////////////////////////
#elif defined(PRECISION_s)

#define DIM_X 16
#define DIM_Y 16

//////// NN ///////
#define BLK_M_nn 96
#define BLK_N_nn 96
#define BLK_K_nn 16

#define DIM_XA_nn 32
#define DIM_YA_nn  8

#define DIM_XB_nn  8
#define DIM_YB_nn 32

//////// NT ///////
#define BLK_M_nt 96
#define BLK_N_nt 96
#define BLK_K_nt 16

#define DIM_XA_nt 32
#define DIM_YA_nt  8

#define DIM_XB_nt 32
#define DIM_YB_nt  8

//////// TT ///////
#define BLK_M_tt 96
#define BLK_N_tt 96
#define BLK_K_tt 16

#define DIM_XA_tt 16
#define DIM_YA_tt 16

#define DIM_XB_tt 32
#define DIM_YB_tt  8

//////// TN ///////
#define BLK_M_tn 96
#define BLK_N_tn 96
#define BLK_K_tn 16

#define DIM_XA_tn 16
#define DIM_YA_tn 16

#define DIM_XB_tn 16
#define DIM_YB_tn 16

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  NoTrans - NoTrans
//

#define BLK_M  BLK_M_nn
#define BLK_N  BLK_N_nn
#define BLK_K  BLK_K_nn
#define DIM_XA DIM_XA_nn
#define DIM_YA DIM_YA_nn
#define DIM_XB DIM_XB_nn
#define DIM_YB DIM_YB_nn

#define version trans_nn
#include "gemm_stencil.cu"

#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  NoTrans - Trans
//

#define BLK_M  BLK_M_nt
#define BLK_N  BLK_N_nt
#define BLK_K  BLK_K_nt
#define DIM_XA DIM_XA_nt
#define DIM_YA DIM_YA_nt
#define DIM_XB DIM_XB_nt
#define DIM_YB DIM_YB_nt

#define version trans_nt
#include "gemm_stencil.cu"

#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#if defined(PRECISION_z) || defined(PRECISION_c)
#define BLK_M  BLK_M_nt
#define BLK_N  BLK_N_nt
#define BLK_K  BLK_K_nt
#define DIM_XA DIM_XA_nt
#define DIM_YA DIM_YA_nt
#define DIM_XB DIM_XB_nt
#define DIM_YB DIM_YB_nt

#define version trans_nc
#include "gemm_stencil.cu"
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Trans - Trans
//

#define BLK_M  BLK_M_tt
#define BLK_N  BLK_N_tt
#define BLK_K  BLK_K_tt
#define DIM_XA DIM_XA_tt
#define DIM_YA DIM_YA_tt
#define DIM_XB DIM_XB_tt
#define DIM_YB DIM_YB_tt

#define version trans_tt
#include "gemm_stencil.cu"
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#if defined(PRECISION_z) || defined(PRECISION_c)
#define BLK_M  BLK_M_tt
#define BLK_N  BLK_N_tt
#define BLK_K  BLK_K_tt
#define DIM_XA DIM_XA_tt
#define DIM_YA DIM_YA_tt
#define DIM_XB DIM_XB_tt
#define DIM_YB DIM_YB_tt

#define version trans_tc
#include "gemm_stencil.cu"
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#define BLK_M  BLK_M_tt
#define BLK_N  BLK_N_tt
#define BLK_K  BLK_K_tt
#define DIM_XA DIM_XA_tt
#define DIM_YA DIM_YA_tt
#define DIM_XB DIM_XB_tt
#define DIM_YB DIM_YB_tt

#define version trans_ct
#include "gemm_stencil.cu"
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#define BLK_M  BLK_M_tt
#define BLK_N  BLK_N_tt
#define BLK_K  BLK_K_tt
#define DIM_XA DIM_XA_tt
#define DIM_YA DIM_YA_tt
#define DIM_XB DIM_XB_tt
#define DIM_YB DIM_YB_tt

#define version trans_cc
#include "gemm_stencil.cu"
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Trans - NoTrans
//

#define BLK_M  BLK_M_tn
#define BLK_N  BLK_N_tn
#define BLK_K  BLK_K_tn
#define DIM_XA DIM_XA_tn
#define DIM_YA DIM_YA_tn
#define DIM_XB DIM_XB_tn
#define DIM_YB DIM_YB_tn

#define version trans_tn
#include "gemm_stencil.cu"
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#if defined(PRECISION_z) || defined(PRECISION_c)
#define BLK_M  BLK_M_tn
#define BLK_N  BLK_N_tn
#define BLK_K  BLK_K_tn
#define DIM_XA DIM_XA_tn
#define DIM_YA DIM_YA_tn
#define DIM_XB DIM_XB_tn
#define DIM_YB DIM_YB_tn

#define version trans_cn
#include "gemm_stencil.cu"
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////

#endif /* _GEMM_STENCIL_H_ */
