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

static PASTIX_INT cs_tol (PASTIX_INT i, PASTIX_INT j, PASTIX_FLOAT aij, void *tol)
{
  return (fabs (aij) > *((PASTIX_FLOAT *) tol)) ;
}

PASTIX_INT cs_droptol (cs *A, PASTIX_FLOAT tol)
{
  return (cs_fkeep (A, &cs_tol, &tol)) ;    /* keep all large entries */
}

static PASTIX_INT cs_nonzero (PASTIX_INT i, PASTIX_INT j, PASTIX_FLOAT aij, void *other)
{
  return (aij != 0) ;
}

PASTIX_INT cs_dropzeros (cs *A)
{
  return (cs_fkeep (A, &cs_nonzero, NULL)) ;/* keep all nonzero entries */
}

/* C = alpha*A + beta*B */
cs *cs_add ( const cs *A, const cs *B, PASTIX_FLOAT alpha, PASTIX_FLOAT beta )
{
  PASTIX_INT p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values ;
  PASTIX_FLOAT *x, *Bx, *Cx ;
  cs *C ;
  if (!A || !B) return (NULL) ;/* check inputs */
  m = A->m ; anz = A->p [A->n] ;
  n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
  w = cs_calloc (m, sizeof (PASTIX_INT)) ;
  values = (A->x != NULL) && (Bx != NULL) ;
  x = values ? cs_malloc (m, sizeof (PASTIX_FLOAT)) : NULL ;
  C = cs_spalloc (m, n, anz + bnz, values, 0) ;
  if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (j = 0 ; j < n ; j++)
    {
      Cp [j] = nz ;/* column j of C starts here */
      nz = cs_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
      nz = cs_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
      if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
  Cp [n] = nz ;/* finalize the last column of C */
  cs_sprealloc (C, 0) ;/* remove extra space from C */
  return (cs_done (C, w, x, 1)) ;/* success; free workspace, return C */
}

/* removes duplicate entries from A */
PASTIX_INT cs_dupl (cs *A)
{
  PASTIX_INT i, j, p, q, nz = 0, n, m, *Ap, *Ai, *w ;
  PASTIX_FLOAT *Ax ;
  if (!A) return (0) ;/* check inputs */
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  w = cs_malloc (m, sizeof (PASTIX_INT)) ;/* get workspace */
  if (!w) return (0) ;/* out of memory */
  for (i = 0 ; i < m ; i++) w [i] = -1 ;/* row i not yet seen */
  for (j = 0 ; j < n ; j++)
    {
      q = nz ;/* column j will start at q */
      for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	  i = Ai [p] ;/* A(i,j) is nonzero */
	  if (w [i] >= q)
	    {
	      Ax [w [i]] += Ax [p] ;/* A(i,j) is a duplicate */
	    }
	  else
	    {
	      w [i] = nz ;/* record where row i occurs */
	      Ai [nz] = i ;/* keep A(i,j) */
	      Ax [nz++] = Ax [p] ;
	    }
	}
      Ap [j] = q ;/* record start of column j */
    }
  Ap [n] = nz ;/* finalize A */
  cs_free (w) ;/* free workspace */
  return (cs_sprealloc (A, 0)) ;/* remove extra space from A */
}

/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
PASTIX_INT cs_entry (cs *T, PASTIX_INT i, PASTIX_INT j, PASTIX_FLOAT x)
{
  if (!T || (T->nz >= T->nzmax && !cs_sprealloc (T, 2*(T->nzmax)))) return(0);
  if (T->x) T->x [T->nz] = x ;
  T->i [T->nz] = i ;
  T->p [T->nz++] = j ;
  T->m = CS_MAX (T->m, i+1) ;
  T->n = CS_MAX (T->n, j+1) ;
  return (1) ;
}

/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
PASTIX_INT cs_fkeep (cs *A, PASTIX_INT (*fkeep) (PASTIX_INT, PASTIX_INT, PASTIX_FLOAT, void *), void *other)
{
  PASTIX_INT baseval = 1 ;
  PASTIX_INT j, p, nz = 0, n, *Ap, *Ai ;
  PASTIX_FLOAT *Ax ;
  if (!A || !fkeep) return (-1) ;    /* check inputs */
  n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  for (j = 0 ; j < n ; j++)
    {
      p = Ap [j] - baseval ;    /* get current location of col j */
      Ap [j] = nz + baseval ;    /* record new location of col j */
      for ( ; p < Ap [j+1] - baseval ; p++)
	{
	  if (fkeep (Ai [p], j, Ax ? Ax [p] : 1, other))
	    {
	      if (Ax) Ax [nz] = Ax [p] ;  /* keep A(i,j) */
	      Ai [nz++] = Ai [p] ;
	    }
	  else printf("drop %ld,%ld\n",(long)j,(long)p);
	}
    }
  /* finalize A and return nnz(A) */
  return (Ap [n] = nz + baseval) ;
}

/* y = A*x+y */
PASTIX_INT cs_gaxpy (const cs *A, const PASTIX_FLOAT *x, PASTIX_FLOAT *y)
{
  PASTIX_INT p, j, n, *Ap, *Ai ;
  PASTIX_FLOAT *Ax ;
  if (!A || !x || !y) return (0) ;    /* check inputs */
  n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  for (j = 0 ; j < n ; j++)
    {
      for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	  y [Ai [p]] += Ax [p] * x [j] ;
	}
    }
  return (1) ;
}


/* C = A*B */
cs *cs_multiply (const cs *A, const cs *B)
{
  PASTIX_INT p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
  PASTIX_FLOAT *x, *Bx, *Cx ;
  cs *C ;
  if (!A || !B) return (NULL) ;/* check inputs */
  m = A->m ; anz = A->p [A->n] ;
  n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
  w = cs_calloc (m, sizeof (PASTIX_INT)) ;
  values = (A->x != NULL) && (Bx != NULL) ;
  x = values ? cs_malloc (m, sizeof (PASTIX_FLOAT)) : NULL ;
  C = cs_spalloc (m, n, anz + bnz, values, 0) ;
  if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
  Cp = C->p ;
  for (j = 0 ; j < n ; j++)
    {
      if (nz + m > C->nzmax && !cs_sprealloc (C, 2*(C->nzmax)+m))
	{
	  return (cs_done (C, w, x, 0)) ;/* out of memory */
	} 
      Ci = C->i ; Cx = C->x ;/* C may have been reallocated */
      Cp [j] = nz ;/* column j of C starts here */
      for (p = Bp [j] ; p < Bp [j+1] ; p++)
	{
	  nz = cs_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
	}
      if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
  Cp [n] = nz ;/* finalize the last column of C */
  cs_sprealloc (C, 0) ;/* remove extra space from C */
  return (cs_done (C, w, x, 1)) ;/* success; free workspace, return C */
}

/* 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum */
PASTIX_FLOAT cs_norm (const cs *A)
{
  PASTIX_INT p, j, n, *Ap ;
  PASTIX_FLOAT *Ax,  norm = 0, s ;
  if (!A || !A->x) return (-1) ;/* check inputs */
  n = A->n ; Ap = A->p ; Ax = A->x ;
  for (j = 0 ; j < n ; j++)
    {
      for (s = 0, p = Ap [j] ; p < Ap [j+1] ; p++) s += fabs (Ax [p]) ;
      norm = CS_MAX (norm, s) ;
    }
  return (norm) ;
}

/* C = A(P,Q) where P and Q are permutations of 0..m-1 and 0..n-1. */
cs *cs_permute (const cs *A, const PASTIX_INT *Pinv, const PASTIX_INT *Q, PASTIX_INT values)
{
  PASTIX_INT p, j, k, nz = 0, m, n, *Ap, *Ai, *Cp, *Ci ;
  PASTIX_FLOAT *Cx, *Ax ;
  cs *C ;
  if (!A) return (NULL) ;/* check inputs */
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  C = cs_spalloc (m, n, Ap [n], values && Ax != NULL, 0) ;
  if (!C) return (cs_done (C, NULL, NULL, 0)) ;   /* out of memory */
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (k = 0 ; k < n ; k++)
    {
      Cp [k] = nz ;/* column k of C is column Q[k] of A */
      j = Q ? (Q [k]) : k ;
      for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	  if (Cx) Cx [nz] = Ax [p] ;/* row i of A is row Pinv[i] of C */
	  Ci [nz++] = Pinv ? (Pinv [Ai [p]]) : Ai [p] ;
	}
    }
  Cp [n] = nz ;/* finalize the last column of C */
  return (cs_done (C, NULL, NULL, 1)) ;
}

/* Pinv = P', or P = Pinv' */
PASTIX_INT *cs_pinv (PASTIX_INT const *P, PASTIX_INT n)
{
  PASTIX_INT k, *Pinv ;
  if (!P) return (NULL) ;/* P = NULL denotes identity */
  Pinv = cs_malloc (n, sizeof (PASTIX_INT)) ;/* allocate resuult */
  if (!Pinv) return (NULL) ;/* out of memory */
  for (k = 0 ; k < n ; k++) Pinv [P [k]] = k ;/* invert the permutation */
  return (Pinv) ;/* return result */
}

/* C = A' */
cs *cs_transpose (const cs *A, PASTIX_INT values)
{
  PASTIX_INT p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w ;
  PASTIX_FLOAT *Cx, *Ax ;
  cs *C ;
  if (!A) return (NULL) ;
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  C = cs_spalloc (n, m, Ap [n], values && Ax, 0) ;   /* allocate result */
  w = cs_calloc (m, sizeof (PASTIX_INT)) ;
  if (!C || !w) return (cs_done (C, w, NULL, 0)) ;   /* out of memory */
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;   /* row counts */
  cs_cumsum (Cp, w, m) ;   /* row poINTers */
  for (j = 0 ; j < n ; j++)
    {
      for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	  Ci [q = w [Ai [p]]++] = j ;/* place A(i,j) as entry C(j,i) */
	  if (Cx) Cx [q] = Ax [p] ;
	}
    }
  return (cs_done (C, w, NULL, 1)) ;/* success; free w and return C */
}

/* C = compressed-column form of a triplet matrix T */
cs *cs_triplet (const cs *T)
{
  PASTIX_INT m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
  PASTIX_FLOAT *Cx, *Tx ;
  cs *C ;
  if (!T) return (NULL) ;/* check inputs */
  m = T->m ; n = T->n ; Ti = T->i ; Tj = T->p ; Tx = T->x ; nz = T->nz ;
  C = cs_spalloc (m, n, nz, Tx != NULL, 0) ;/* allocate result */
  w = cs_calloc (n, sizeof (PASTIX_INT)) ;/* get workspace */
  if (!C || !w) return (cs_done (C, w, NULL, 0)) ;/* out of memory */
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (k = 0 ; k < nz ; k++) w [Tj [k]]++ ;/* column counts */
  cs_cumsum (Cp, w, n) ;/* column poINTers */
  for (k = 0 ; k < nz ; k++)
    {
      Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
      if (Cx) Cx [p] = Tx [k] ;
    }
  return (cs_done (C, w, NULL, 1)) ;    /* success; free w and return C */
}

/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
PASTIX_INT cs_scatter (const cs *A, PASTIX_INT j, PASTIX_FLOAT beta, PASTIX_INT *w, PASTIX_FLOAT *x, PASTIX_INT mark,
		cs *C, PASTIX_INT nz)
{
  PASTIX_INT i, p, *Ap, *Ai, *Ci ;
  PASTIX_FLOAT *Ax ;
  if (!A || !w || !C) return (-1) ;/* ensure inputs are valid */
  Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
  for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
      i = Ai [p] ;/* A(i,j) is nonzero */
      if (w [i] < mark)
	{
	  w [i] = mark ;/* i is new entry in column j */
	  Ci [nz++] = i ;/* add i to pattern of C(:,j) */
	  if (x) x [i] = beta * Ax [p] ;/* x(i) = beta*A(i,j) */
	}
      else if (x) x [i] += beta * Ax [p] ;/* i exists in C(:,j) already */
    }
  return (nz) ;
}

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
PASTIX_INT cs_cumsum (PASTIX_INT *p, PASTIX_INT *c, PASTIX_INT n)
{
  PASTIX_INT i, nz = 0 ;
  if (!p || !c) return (-1) ;    /* check inputs */
  for (i = 0 ; i < n ; i++)
    {
      p [i] = nz ;
      nz += c [i] ;
      c [i] = p [i] ;
    }
  p [n] = nz ;
  return (nz) ;    /* return sum (c [0..n-1]) */
}

/* wrapper for malloc */
void *cs_malloc (PASTIX_INT n, size_t size)
{
  return (CS_OVERFLOW (n,size) ? NULL : malloc (CS_MAX (n,1) * size)) ;
}

/* wrapper for calloc */
void *cs_calloc (PASTIX_INT n, size_t size)
{
  return (CS_OVERFLOW (n,size) ? NULL : calloc (CS_MAX (n,1), size)) ;
}

/* wrapper for free */
void *cs_free (void *p)
{
  if (p) free (p) ;    /* free p if it is not already NULL */
  return (NULL) ;    /* return NULL to simplify the use of cs_free */
}

/* wrapper for realloc */
void *cs_realloc (void *p, PASTIX_INT n, size_t size, PASTIX_INT *ok)
{
  void *p2 ;
  *ok = !CS_OVERFLOW (n,size) ;    /* guard against PASTIX_INT overflow */
  if (!(*ok)) return (p) ;    /* p unchanged if n too large */
  p2 = realloc (p, CS_MAX (n,1) * size) ; /* realloc the block */
  *ok = (p2 != NULL) ;
  return ((*ok) ? p2 : p) ;    /* return original p if failure */
}

/* allocate a sparse matrix (triplet form or compressed-column form) */
cs *cs_spalloc (PASTIX_INT m, PASTIX_INT n, PASTIX_INT nzmax, PASTIX_INT values, PASTIX_INT triplet)
{
  cs *A = cs_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */
  if (!A) return (NULL) ;    /* out of memory */
  A->m = m ;    /* define dimensions and nzmax */
  A->n = n ;
  A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
  A->nz = triplet ? 0 : -1 ;    /* allocate triplet or comp.col */
  A->p = cs_malloc (triplet ? nzmax : n+1, sizeof (PASTIX_INT)) ;
  A->i = cs_malloc (nzmax, sizeof (PASTIX_INT)) ;
  A->x = values ? cs_malloc (nzmax, sizeof (PASTIX_FLOAT)) : NULL ;
  return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree (A) : A) ;
}

/* change the max # of entries sparse matrix */
PASTIX_INT cs_sprealloc (cs *A, PASTIX_INT nzmax)
{
  PASTIX_INT ok, oki, okj = 1, okx = 1 ;
  if (!A) return (0) ;
  nzmax = (nzmax <= 0) ? (A->p [A->n]) : nzmax ;
  A->i = cs_realloc (A->i, nzmax, sizeof (PASTIX_INT), &oki) ;
  if (A->nz >= 0) A->p = cs_realloc (A->p, nzmax, sizeof (PASTIX_INT), &okj) ;
  if (A->x) A->x = cs_realloc (A->x, nzmax, sizeof (PASTIX_FLOAT), &okx) ;
  ok = (oki && okj && okx) ;
  if (ok) A->nzmax = nzmax ;
  return (ok) ;
}

/* free a sparse matrix */
cs *cs_spfree (cs *A)
{
  if (!A) return (NULL) ;/* do nothing if A already NULL */
  cs_free (A->p) ;
  cs_free (A->i) ;
  cs_free (A->x) ;
  return (cs_free (A)) ;/* free the cs struct and return NULL */
}

/* free workspace and return a sparse matrix result */
cs *cs_done (cs *C, void *w, void *x, PASTIX_INT ok)
{
  cs_free (w) ;/* free workspace */
  cs_free (x) ;
  return (ok ? C : cs_spfree (C)) ;/* return result if OK, else free it */
}
