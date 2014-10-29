#define CS_VER 1    /* CSparse Version 1.2.0 */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Mar 6, 2006"    /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006"

typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
  PASTIX_INT nzmax ;    /* maximum number of entries */
  PASTIX_INT m ;    /* number of rows */
  PASTIX_INT n ;    /* number of columns */
  PASTIX_INT *p ;    /* column poINTers (size n+1) or col indices (size nzmax) */
  PASTIX_INT *i ;    /* row indices, size nzmax */
  PASTIX_FLOAT *x ;    /* numerical values, size nzmax */
  PASTIX_INT nz ;    /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

/* keep all large entries */
PASTIX_INT cs_droptol (cs *A, PASTIX_FLOAT tol) ;
/* keep all nonzero entries */
PASTIX_INT cs_dropzeros (cs *A) ;
/* C = alpha*A + beta*B */
cs *cs_add (const cs *A, const cs *B, PASTIX_FLOAT alpha, PASTIX_FLOAT beta) ;
/* removes duplicate entries from A */
PASTIX_INT cs_dupl (cs *A) ;
/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
PASTIX_INT cs_entry (cs *T, PASTIX_INT i, PASTIX_INT j, PASTIX_FLOAT x) ;
/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
PASTIX_INT cs_fkeep (cs *A, PASTIX_INT (*fkeep) (PASTIX_INT, PASTIX_INT, PASTIX_FLOAT, void *), void *other) ;
/* y = A*x+y */
PASTIX_INT cs_gaxpy (const cs *A, const PASTIX_FLOAT *x, PASTIX_FLOAT *y) ;
/* C = A*B */
cs *cs_multiply (const cs *A, const cs *B) ;
/* 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum */
PASTIX_FLOAT cs_norm (const cs *A) ;
/* C = A(P,Q) where P and Q are permutations of 0..m-1 and 0..n-1. */
cs *cs_permute (const cs *A, const PASTIX_INT *P, const PASTIX_INT *Q, PASTIX_INT values) ;
/* Pinv = P', or P = Pinv' */
PASTIX_INT *cs_pinv (const PASTIX_INT *P, PASTIX_INT n) ;
/* C = A' */
cs *cs_transpose (const cs *A, PASTIX_INT values) ;
/* C = compressed-column form of a triplet matrix T */
cs *cs_triplet (const cs *T) ;
/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
PASTIX_INT cs_scatter (const cs *A, PASTIX_INT j, PASTIX_FLOAT beta, PASTIX_INT *w, PASTIX_FLOAT *x, PASTIX_INT mark,
                cs *C, PASTIX_INT nz) ;
/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
PASTIX_INT cs_cumsum (PASTIX_INT *p, PASTIX_INT *c, PASTIX_INT n) ;

/* utilities */
/* wrapper for malloc */
void *cs_malloc (PASTIX_INT n, size_t size) ;
/* wrapper for calloc */
void *cs_calloc (PASTIX_INT n, size_t size) ;
/* wrapper for free */
void *cs_free (void *p) ;
/* wrapper for realloc */
void *cs_realloc (void *p, PASTIX_INT n, size_t size, PASTIX_INT *ok) ;
/* allocate a sparse matrix (triplet form or compressed-column form) */
cs *cs_spalloc (PASTIX_INT m, PASTIX_INT n, PASTIX_INT nzmax, PASTIX_INT values, PASTIX_INT triplet) ;
/* change the max # of entries sparse matrix */
PASTIX_INT cs_sprealloc (cs *A, PASTIX_INT nzmax) ;
/* free a sparse matrix */
cs *cs_spfree (cs *A) ;
/* free workspace and return a sparse matrix result */
cs *cs_done (cs *C, void *w, void *x, PASTIX_INT ok) ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (PASTIX_INT) size)
