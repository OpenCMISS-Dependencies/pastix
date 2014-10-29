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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "common_pastix.h"
/* #include "symbol.h" */

#include "sparRow.h"


PASTIX_INT initCS(csptr amat, PASTIX_INT len)
{
/*---------------------------------------------------------------------- 
| Initialize SparRow structs.
|----------------------------------------------------------------------
| on entry: 
|========== 
| ( amat )  =  Pointer to a SparRow struct.
|     len   =  size of matrix
|
| On return:
|===========
|
|  amat->n
|      ->*nnzrow
|      ->**ja
|      ->**ma
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   amat->n = len;
   amat->inarow = 0;
   if(len > 0)
     {
       MALLOC_INTERN(amat->nnzrow, len, PASTIX_INT);
       MALLOC_INTERN(amat->ja,     len, PASTIX_INT *);
       MALLOC_INTERN(amat->ma,     len, double *);

       bzero( amat->nnzrow, sizeof(PASTIX_INT)*len);
       bzero( amat->ja, sizeof(PASTIX_INT *)*len);
       bzero( amat->ma, sizeof(double *)*len);
       amat->jatab = NULL;
       amat->matab = NULL;
     }
   else
     {
       amat->nnzrow = NULL;
       amat->ja = NULL;
       amat->ma = NULL;
       amat->jatab = NULL;
       amat->matab = NULL;
     }

   return 0;
}


PASTIX_INT cleanCS(csptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for SparRow structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SparRow struct.
|     len   =  size of matrix
|--------------------------------------------------------------------*/
   /*	*/
  PASTIX_INT i; 
  if (amat == NULL) return 0;
  /*if (amat->n < 1) return 0;*/

  if(amat->inarow == 0)
    {
      if(amat->n > 0) /** The sparse part has not been deallocated **/
	{
	  for (i=0; i<amat->n; i++) {
	    if (amat->nnzrow[i] > 0) {
	      if (amat->ma[i]) memFree_null(amat->ma[i]);
	      if (amat->ja[i]) memFree_null(amat->ja[i]);
	    }
	  }	
	}
    }
  else
    {
      if(amat->inarow != 2 && amat->jatab != NULL)
	{
	  memFree_null(amat->matab);
	  memFree_null(amat->jatab);
	}
      amat->inarow = 0;

    }

  if (amat->ma) memFree_null(amat->ma);
  if (amat->ja) memFree_null(amat->ja);
  if (amat->nnzrow) memFree_null(amat->nnzrow); 
     

  return 0;
}


PASTIX_INT CSnnz(csptr mat)
{
  /*-----------------------------------/
  / Return number of non zero entries  /
  / in a matrix in SparRow format      /
  /-----------------------------------*/
  PASTIX_INT nnz;
  PASTIX_INT i;

  nnz = 0;
  for(i=0;i<mat->n;i++)
    nnz += mat->nnzrow[i];
      
  return nnz;
}

PASTIX_INT CS_RowPerm(csptr mat, PASTIX_INT *perm)
{
/*----------------------------------------------------------------------
|
| This subroutine permutes the rows of a matrix in SparRow format. 
| rperm  computes B = P A  where P is a permutation matrix.  
| The permutation P is defined through the array perm: for each j, 
| perm[j] represents the destination row number of row number j. 
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat) = a matrix stored in SparRow format.
|
|
| on return:
| ----------
| (amat) = P A stored in SparRow format.
|
| INTeger value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
  PASTIX_INT **addj = NULL;
  PASTIX_INT * nnz  = NULL;
  PASTIX_INT i, size;
  double **addm = NULL;

   size = mat->n;
   if(size == 0)
     return 0;

   MALLOC_INTERN(addj, size, PASTIX_INT *);
   MALLOC_INTERN(addm, size, double *);
   MALLOC_INTERN(nnz,  size, PASTIX_INT);

   for (i=0; i<size; i++) {
      addj[perm[i]] = mat->ja[i];
      addm[perm[i]] = mat->ma[i];
      nnz[perm[i]] = mat->nnzrow[i];
   }
   for (i=0; i<size; i++) {
      mat->ja[i] = addj[i];
      mat->ma[i] = addm[i];
      mat->nnzrow[i] = nnz[i];
   }
   memFree(addj);
   memFree(addm);
   memFree(nnz);
   return 0;
}
/*------- end of rperm ------------------------------------------------- 
|---------------------------------------------------------------------*/

PASTIX_INT CS_ColPerm(csptr mat, PASTIX_INT *perm) 
{
/*----------------------------------------------------------------------
|
| This subroutine permutes the columns of a matrix in SparRow format.
| cperm computes B = A P, where P is a permutation matrix.
| that maps column j INTo column perm(j), i.e., on return 
| The permutation P is defined through the array perm: for each j, 
| perm[j] represents the destination column number of column number j. 
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (mat) = a matrix stored in SparRow format.
|
|
| on return:
| ----------
| (mat) = A P stored in SparRow format.
|
|---------------------------------------------------------------------*/
   PASTIX_INT i, j,  size, *aja;
   PASTIX_INT *newj = NULL;
   PASTIX_INT maxrow;
   size = mat->n;
   if(size == 0)
     return 0;
   
   maxrow = 0;
   for(i=0;i<mat->n;i++)
     if(mat->nnzrow[i]>maxrow)
       maxrow = mat->nnzrow[i];
   
   if(maxrow == 0)
     return 0;

   MALLOC_INTERN(newj, maxrow, PASTIX_INT);

   for (i=0; i<size; i++) {
      aja = mat->ja[i];
      for (j=0; j<mat->nnzrow[i]; j++)
	newj[j] = perm[aja[j]];
  
      for (j=0; j<mat->nnzrow[i]; j++)
	 aja[j] = newj[j];
      /*mat->ja[i] = aja;*/
   }
   memFree(newj);
   return 0;
}
/*------- end of cperm ------------------------------------------------- 
|---------------------------------------------------------------------*/
PASTIX_INT CS_Perm(csptr mat, PASTIX_INT *perm) 
{
/*----------------------------------------------------------------------
|
| This subroutine permutes the rows and columns of a matrix in 
| SparRow format.  dperm computes B = P^T A P, where P is a permutation 
| matrix.
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat) = a matrix stored in SparRow format.
|
|
| on return:
| ----------
| (amat) = P^T A P stored in SparRow format.
|
| INTeger value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
   if (CS_RowPerm(mat, perm)) return 1;
   if (CS_ColPerm(mat, perm)) return 1;
   return 0;
}

