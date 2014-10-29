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
/************************************************************/
/**                                                        **/
/**   NAME       : amalgamate.c                            **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 15/08/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "common_pastix.h"
/* #include "symbol.h" */
/* #include "queue.h" */
#include "sparRow.h"
/* #include "sort_row.h" */
#include "SF_Direct.h"


void UnionSet(PASTIX_INT *set1, PASTIX_INT n1, PASTIX_INT *set2, PASTIX_INT n2, PASTIX_INT *set, PASTIX_INT *n);

void SF_Direct(csptr A, PASTIX_INT cblknbr, PASTIX_INT *rangtab, PASTIX_INT *treetab, csptr P)
{
  /********************************************************/
  /* This function computes the direct factor nnz pattern */
  /* of a matrix A given the supernode partition          */
  /********************************************************/
  PASTIX_INT i,j,k;
  PASTIX_INT ind, nnznbr, father;
  PASTIX_INT *tmpj      = NULL;
  PASTIX_INT *tmp       = NULL;
  PASTIX_INT *tmp2      = NULL;
  PASTIX_INT *ja        = NULL;
  PASTIX_INT *node2cblk = NULL;

  MALLOC_INTERN(tmpj,      A->n, PASTIX_INT);
  MALLOC_INTERN(tmp,       A->n, PASTIX_INT);
  MALLOC_INTERN(node2cblk, A->n, PASTIX_INT);


  for(k=0;k<cblknbr;k++)
    for(i=rangtab[k];i<rangtab[k+1];i++)
      node2cblk[i] = k;

  /** Compute the nnz structure of each supernode in A **/
  for(k=0;k<cblknbr;k++)
    {
      ind = rangtab[k];


      /** Put the diagonal elements (A does not contains them) **/
      j = 0;
      for(i=rangtab[k];i<rangtab[k+1];i++)
	tmpj[j++] = i; 
      nnznbr = j;

      for(i=rangtab[k];i<rangtab[k+1];i++)
      {
	j = 0;
	while(j<A->nnzrow[i] && A->ja[i][j] <= i)
	  j++;

	/** Realise la fusion de 2 listes triees croissantes **/
	UnionSet(tmpj, nnznbr, A->ja[i]+j, A->nnzrow[i]-j, tmp, &ind);
	
	/** echange tmpj et le resultat de la fusion **/
	nnznbr = ind;
	tmp2 = tmpj;
        tmpj = tmp;
	tmp  = tmp2;
      }
#ifdef DEBUG_KASS
      ind = 0;
      for(j=rangtab[k]; j < rangtab[k+1];j++)
	ASSERT(tmpj[ind++] == j, MOD_KASS);
#endif

#ifdef DEBUG_KASS
      ASSERT(nnznbr > 0, MOD_KASS);
#endif
      P->nnzrow[k] = nnznbr;
      MALLOC_INTERN(P->ja[k], nnznbr, PASTIX_INT);
      memCpy(P->ja[k], tmpj, sizeof(PASTIX_INT)*nnznbr);
      P->ma[k]=NULL;
    }

 
  

  
  P->n = cblknbr;

  /** Compute the symbolic factorization **/
  for(k=0;k<cblknbr;k++)
    {
      /*father = treetab[k];*/
      i = 0;
      ja = P->ja[k];
      while(i < P->nnzrow[k] && node2cblk[ja[i]] <= k)
	i++;
      if(i<P->nnzrow[k])
	father = node2cblk[ja[i]];
      else
	father = -1;
      treetab[k] = father;

      if(father != k && father > 0)
	{
	  /*i = 0;
	  while(i < P->nnzrow[k] && P->ja[k][i] <= rangtab[father])
	  i++;
	  if(i == P->nnzrow[k])
	  continue;*/

	  UnionSet(P->ja[k]+i, P->nnzrow[k]-i, P->ja[father], P->nnzrow[father], tmpj, &nnznbr);

	  memFree(P->ja[father]);
	  MALLOC_INTERN(P->ja[father], nnznbr, PASTIX_INT);
	  memCpy(P->ja[father], tmpj, sizeof(PASTIX_INT)*nnznbr);
	  P->nnzrow[father] = nnznbr;
	}
    }


#ifdef DEBUG_KASS
  for(i=0;i<cblknbr;i++)
    {
      PASTIX_INT j;
      ASSERT(P->nnzrow[i] >= rangtab[i+1]-rangtab[i], MOD_KASS);
      k = 0;
      for(j=rangtab[i]; j < rangtab[i+1];j++)
	ASSERT(P->ja[i][k++] == j, MOD_KASS);
    }
  
  /** Check that all terms of A are in the pattern **/
  for(k=0;k<cblknbr;k++)
    {
      /** Put the diagonal elements (A does not contains them) **/
      for(i=rangtab[k];i<rangtab[k+1];i++)
	{
	  j = 0;
	  while(j<A->nnzrow[i] && A->ja[i][j] < i)
	    j++;

	  for(ind = j;ind < A->nnzrow[i];ind++)
	    assert(A->ja[i][ind] >= i);
	  for(ind = j+1;ind < A->nnzrow[i];ind++)
	    assert(A->ja[i][ind] > A->ja[i][ind-1]);  
	  
	  
	  UnionSet(P->ja[k], P->nnzrow[k], A->ja[i]+j, A->nnzrow[i]-j, tmp, &ind);
	  if(ind > P->nnzrow[k])
	    fprintf(stderr, "k=%ld [%ld %ld]  i=%ld ind %ld nnz %ld \n", 
		    (long)k, (long)rangtab[k], (long)rangtab[k+1], (long)i, (long)ind, (long)P->nnzrow[k]);

	  ASSERT(ind <= P->nnzrow[k], MOD_KASS);
	}
    }
#endif

  memFree(node2cblk);
  memFree(tmpj);
  memFree(tmp);
  
}


void UnionSet(PASTIX_INT *set1, PASTIX_INT n1, PASTIX_INT *set2, PASTIX_INT n2, PASTIX_INT *set, PASTIX_INT *n)
{
  /********************************************************/
  /* Compute the union of two sorted set                  */
  /* set must have a big enough size to contain the union */
  /* i.e n1+n2                                            */
  /********************************************************/
  PASTIX_INT ind, ind1, ind2;

  ind = 0;
  ind1 = 0;
  ind2 = 0;

#ifdef DEBUG_KASS
 {
   PASTIX_INT i;
   for(i=0;i<n1-1;i++)
     ASSERT(set1[i] < set1[i+1], MOD_KASS);
   for(i=0;i<n2-1;i++)
     ASSERT(set2[i] < set2[i+1], MOD_KASS);
 }
#endif

  while(ind1 < n1 && ind2 < n2)
    {
      if(set1[ind1] == set2[ind2])
	{
	  set[ind] = set1[ind1];
	  ind++;
	  ind1++;
	  ind2++;
	  continue;
	}
      
      if(set1[ind1] < set2[ind2])
	{
	  set[ind] = set1[ind1];
	  ind++;
	  ind1++;
	  continue;
	}
      
      if(set1[ind1] > set2[ind2])
	{
	  set[ind] = set2[ind2];
	  ind++;
	  ind2++;
	  continue;
	}
    }

  while(ind1 < n1)
    {
      set[ind] = set1[ind1];
      ind++;
      ind1++;
    }
  while(ind2 < n2)
    {
      set[ind] = set2[ind2];
      ind++;
      ind2++;
    }
#ifdef DEBUG_KASS
  ASSERT(ind <= ind1 +ind2, MOD_KASS);
  ASSERT(ind >= MAX(ind1, ind2), MOD_KASS);
#endif


  *n = ind;
}

