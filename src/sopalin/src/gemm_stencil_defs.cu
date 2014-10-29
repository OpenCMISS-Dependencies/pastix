
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

#ifdef TEXTURE_1D

#if defined(PRECISION_z)

#define conj(A)          cuConj(A)
#define add(A, B)        cuCadd(A, B)
#define mul(A, B)        cuCmul(A, B)
#define fma(A, B, C) C = cuCfma(A, B, C)
#define make_FloatingPoint(x, y) make_cuDoubleComplex(x, y);

static __device__
FloatingPoint_t tex_fetch(texture<int4> tex_ref, int coord)
{
    int4 v = tex1Dfetch(tex_ref, coord);
    return make_cuDoubleComplex(__hiloint2double(v.y, v.x), __hiloint2double(v.w, v.z));
}

texture<int4, 1, cudaReadModeElementType> tex_ref_A;
texture<int4, 1, cudaReadModeElementType> tex_ref_B;

///////////////////////////////////////////////////////////////////////////////////////////////////
#elif defined(PRECISION_c)

#define conj(A)          cuConjf(A)
#define add(A, B)        cuCaddf(A, B)
#define mul(A, B)        cuCmulf(A, B)
#define fma(A, B, C) C = cuCfmaf(A, B, C)
#define make_FloatingPoint(x, y) make_cuFloatComplex(x, y);

static __device__
FloatingPoint_t tex_fetch(texture<float2> tex_ref, int coord)
{
    return tex1Dfetch(tex_ref, coord);
}

texture<float2, 1, cudaReadModeElementType> tex_ref_A;
texture<float2, 1, cudaReadModeElementType> tex_ref_B;

///////////////////////////////////////////////////////////////////////////////////////////////////
#elif defined(PRECISION_d)

#define conj(A)           (A)
#define add(A, B)         (A+B)
#define mul(A, B)         (A*B)
#define fma(A, B, C) C += (A*B)
#define make_FloatingPoint(x, y) (x)

static __device__
FloatingPoint_t tex_fetch(texture<int2> tex_ref, int coord)
{
    int2 v = tex1Dfetch(tex_ref, coord);
    return __hiloint2double(v.y, v.x);
}

texture<int2, 1, cudaReadModeElementType> tex_ref_A;
texture<int2, 1, cudaReadModeElementType> tex_ref_B;

///////////////////////////////////////////////////////////////////////////////////////////////////
#elif defined(PRECISION_s)

#define conj(A)           (A)
#define add(A, B)         (A+B)
#define mul(A, B)         (A*B)
#define fma(A, B, C) C += (A*B)
#define make_FloatingPoint(x, y) (x)

static __device__
FloatingPoint_t tex_fetch(texture<float> tex_ref, int coord)
{
    return tex1Dfetch(tex_ref, coord);
}

texture<float, 1, cudaReadModeElementType> tex_ref_A;
texture<float, 1, cudaReadModeElementType> tex_ref_B;

#endif

#endif /* defined(TEXTURE_1D) */

///////////////////////////////////////////////////////////////////////////////////////////////////
