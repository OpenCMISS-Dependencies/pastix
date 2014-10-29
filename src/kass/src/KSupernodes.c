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
/**   NAME       : KSupernodes.c                           **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10/02/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "common_pastix.h"
/* #include "symbol.h" */
#include "queue.h"
#include "sparRow.h"
/* #include "sort_row.h" */
/* #include "SF_level.h" */
#include "KSupernodes.h"


PASTIX_INT  Merge_Cost(PASTIX_INT a, PASTIX_INT b, csptr P, PASTIX_INT *colweight);
PASTIX_INT  col_match(PASTIX_INT n1, PASTIX_INT *ja1, PASTIX_INT n2, PASTIX_INT *ja2);
void col_merge(PASTIX_INT n1, PASTIX_INT *ja1, PASTIX_INT n2, PASTIX_INT *ja2, PASTIX_INT *n, PASTIX_INT *ja );

void KSupernodes(csptr P, double rat, PASTIX_INT snodenbr, PASTIX_INT *snodetab, PASTIX_INT *cblknbr, PASTIX_INT **rangtab)
{
  /**********************************************************************/
  /* Find the supernodes in ILU(k) with a global fill-in tolerance of   */
  /* rat                                                                */
  /**********************************************************************/
  PASTIX_INT i, j, k, ind;
  PASTIX_INT *nnzadd    = NULL;
  PASTIX_INT *prevcol   = NULL;
  PASTIX_INT *colweight = NULL;
  double key;
  PASTIX_INT *tmpj      = NULL;

  Queue heap;
  long fillmax, fill;


  fillmax = (long)(CSnnz(P)*rat);
  queueInit(&heap, 2*P->n);
  MALLOC_INTERN(colweight, P->n,   PASTIX_INT);
  MALLOC_INTERN(prevcol,   P->n+1, PASTIX_INT);
  MALLOC_INTERN(nnzadd,    P->n,   PASTIX_INT);
  MALLOC_INTERN(tmpj,      P->n,   PASTIX_INT);
  
  for(k=0;k<snodenbr;k++)
    {
      prevcol[snodetab[k]] = -1;
      colweight[snodetab[k]] = 1;
      for(i=snodetab[k]+1;i<snodetab[k+1];i++)
	{
	  prevcol[i] = i-1; 
	  colweight[i] = 1;
	}
    }
  prevcol[P->n] = -1;

  /** Compute the fill to merge col i and col i+1 */
  for(k=0;k<snodenbr;k++)
    for(i=snodetab[k];i<snodetab[k+1]-1;i++)
      {
	/*ind = col_match(P->nnzrow[i], P->ja[i], P->nnzrow[i+1], P->ja[i+1]);
	  nnzadd[i] = colweight[i]*(P->nnzrow[i]-1-ind) +
  colweight[i+1]*(P->nnzrow[i+1]-ind);*/
	nnzadd[i] = Merge_Cost(i, i+1, P, colweight);
      }

  /** Merge all the columns so it doesn't add fill-in **/
  for(k=0;k<snodenbr;k++)
    for(i=snodetab[k];i<snodetab[k+1]-1;i++)
      {

#ifdef DEBUG_KASS
	ASSERT(colweight[i] > 0, MOD_KASS);
#endif
	ind = i;
	while(ind < snodetab[k+1]-1 && nnzadd[ind] == 0)
	  {
	    ind++;
	    colweight[i] += colweight[ind];
	    
	    /** This column does not exist any more **/
	    colweight[ind] = 0; 
	    prevcol[ind] = -1;
	    P->nnzrow[ind] = 0;
	    memFree_null(P->ja[ind]);
	    /*nnzadd[ind] = 0;*/
	  }
#ifdef DEBUG_KASS
	ASSERT(ind < snodetab[k+1], MOD_KASS);
#endif
	if(ind < snodetab[k+1]-1)
	  prevcol[ind+1] = i;

	i=ind;
      }

  /*#ifdef DEBUG_KASS*/
  ind = 0;
  for(k=0;k<snodenbr;k++)
    for(i=snodetab[k];i<snodetab[k+1]-1;i++)
      if(colweight[i] > 0)
	{
	  ind++;
#ifdef DEBUG_KASS
	  if(i+colweight[i]<snodetab[k+1])
	    {
	      ASSERT(prevcol[i+colweight[i]] == i, MOD_KASS);
	    }
	  else
	    {
	      ASSERT(prevcol[i+colweight[i]] == -1, MOD_KASS);
	    }
#endif
	}
  fprintf(stderr, "KASS: NUMBER OF CBLK IN THE NON AMALGAMATED SYMBOL MATRIX FOR ILUK = %ld \n", (long)ind);
    /*#endif*/



#ifdef DEBUG_KASS
  ind = 0;
  for(i=0;i<P->n;i++)
    ind += colweight[i];
  ASSERT(ind == P->n, MOD_KASS);
#endif
  
  /*** Recompute the nnzadd array for this partition ***/
   for(k=0;k<snodenbr;k++)
    for(i=snodetab[k];i<snodetab[k+1]-1;i++)
      {
	j = i+colweight[i];
	if(prevcol[j] >= 0)
	  {
#ifdef DEBUG_KASS
	    ASSERT(j < snodetab[k+1], MOD_KASS);
#endif
	    /*ind = col_match(P->nnzrow[i], P->ja[i], P->nnzrow[j], P->ja[j]);
	      nnzadd[i] = colweight[i]*(P->nnzrow[i]-1-ind) +
   colweight[j]*(P->nnzrow[j]-ind);*/
	    nnzadd[i] = Merge_Cost(i, j, P, colweight);
#ifdef DEBUG_KASS
	    ASSERT(nnzadd[i]>0, MOD_KASS);
#endif

	    /*** Add the merge fill in a heap ***/
#ifdef DEBUG_KASS
	    ASSERT(prevcol[j] == i, MOD_KASS);
#endif
	    queueAdd(&heap, i, (double) nnzadd[i]);
	  }

	i = j-1;
      }
   
   /*** Merge supernodes untill we reach the fillmax limit ****/
   fill = 0.0;
   while(queueSize(&heap)>0 && fill < fillmax)
     {
       i = queueGet2(&heap, &key, NULL);

       

       if(nnzadd[i] != key || colweight[i] <= 0)
	 continue;

       if(fill + nnzadd[i] > fillmax)
	 break;
       else
	 fill += nnzadd[i];

       j = i + colweight[i];
       if(prevcol[i+colweight[i]] == -1)
	 continue;



#ifdef DEBUG_KASS
       ASSERT(prevcol[j] == i, MOD_KASS);
#endif
       /*fprintf(stdout, "Merge col %ld and col %ld \n", i, j);*/

       /** Merge supernode i and supernode j **/
       col_merge(P->nnzrow[i], P->ja[i], P->nnzrow[j], P->ja[j], &ind, tmpj);
       P->nnzrow[i] = ind;
       P->ja[i] = (PASTIX_INT *)realloc(P->ja[i], sizeof(PASTIX_INT)*ind);
       memcpy(P->ja[i], tmpj, sizeof(PASTIX_INT)*ind);

       assert(P->nnzrow[j]>0);

       P->nnzrow[j] = 0;
       memFree(P->ja[j]);

#ifdef DEBUG_KASS
       ASSERT(colweight[j] > 0, MOD_KASS);
#endif
       colweight[i] += colweight[j];

#ifdef DEBUG_KASS
       for(k=0;k<colweight[i];k++)
	 assert(P->ja[i][k] == i+k);
#endif
	 


       /*** Suppress col j ***/
       colweight[j] = 0;
       nnzadd[j] = -1;
       prevcol[j] = -1;

       /** Recompute the nnzadd of coli  **/
       j = i+colweight[i];
       if(prevcol[j]>=0)
	 {
	   prevcol[j] = i;
	   /*ind = col_match(P->nnzrow[i], P->ja[i], P->nnzrow[j], P->ja[j]);
	     nnzadd[i] = colweight[i]*(P->nnzrow[i]-1-ind) +
       colweight[j]*(P->nnzrow[j]-ind);*/
	   nnzadd[i] = Merge_Cost(i, j, P, colweight);
	   /*** Add the merge fill in a heap ***/
#ifdef DEBUG_KASS
	   assert(i+colweight[i]<P->n);
	   assert(prevcol[i+colweight[i]] == i);
#endif
	   queueAdd(&heap, i, (double) nnzadd[i]);
	 }

       /** Recompute the nnzadd of previous cblk of col i  **/
       j = prevcol[i];

       if(j>=0)
	 {
	   /*ind = col_match(P->nnzrow[j], P->ja[j], P->nnzrow[i], P->ja[i]);
	     nnzadd[j] = colweight[j]*(P->nnzrow[j]-1-ind) +
	 colweight[i]*(P->nnzrow[i]-ind);*/
	   nnzadd[j] = Merge_Cost(j, i, P, colweight);
	   /*** Add the merge fill in a heap ***/
#ifdef DEBUG_KASS
	   assert(j+colweight[j] == i);
	   assert(j+colweight[j]<P->n);
	   assert(prevcol[j+colweight[j]] == j);
#endif
	   queueAdd(&heap, j, (double) nnzadd[j]);
	 }
     }

   /*** Compute the partition ****/
   ind = 0;
   for(i=0;i<P->n;i++)
     if(colweight[i]>0)
       ind++;
   
   *cblknbr = ind;
   MALLOC_INTERN(*rangtab, ind+1, PASTIX_INT);
   ind = 0;
   k = 0;
   while(ind<P->n)
     {
       (*rangtab)[k] = ind;
       k++;
       ind += colweight[ind];
     }

   (*rangtab)[(*cblknbr)] = ind;
#ifdef DEBUG_KASS
   for(k=0;k<*cblknbr;k++)
     assert(colweight[(*rangtab)[k]] == (*rangtab)[k+1]-(*rangtab)[k]);
#endif



   /**  NEW SINCE amalgamate.c OIMBE **/
   /** Compress P ***/
   ind = 0;
   while(ind < P->n && P->nnzrow[ind] >0 )
     {
#ifdef DEBUG_KASS
       ASSERT(colweight[ind] > 0, MOD_KASS);
#endif
       ind++;
     }
   for(i=ind;i<P->n;i++)
     if(P->nnzrow[i] > 0)
       {
#ifdef DEBUG_KASS
	 ASSERT(ind < i, MOD_KASS);
	 ASSERT(colweight[i] > 0, MOD_KASS);
#endif
	 P->nnzrow[ind] = P->nnzrow[i];
	 P->ja[ind] = P->ja[i];
	 P->nnzrow[i] = 0;
	 P->ja[i] = NULL;
	 ind++;
       }
   P->n = *cblknbr;
   /******/



   memFree(tmpj);
   memFree(colweight);
   memFree(nnzadd);
   memFree(prevcol);
   queueExit(&heap);
}



PASTIX_INT col_match(PASTIX_INT n1, PASTIX_INT *ja1, PASTIX_INT n2, PASTIX_INT *ja2)
{
  PASTIX_INT i1, i2;
  PASTIX_INT matched;
  i1 = 0;
  i2 = 0;
  matched = 0;
  while(i1 < n1 && i2 < n2)
    {
      if(ja1[i1] < ja2[i2])
	{
	  i1++;
	  continue;
	}

      if(ja1[i1] > ja2[i2])
	{
	  i2++;
	  continue;
	}
      
      matched++;
      i1++;
      i2++;
    }
#ifdef DEBUG_KASS
  assert(matched <= n2);
  assert(matched <= n1);
#endif
  
  return matched;
}


void col_merge(PASTIX_INT n1, PASTIX_INT *ja1, PASTIX_INT n2, PASTIX_INT *ja2, PASTIX_INT *n, PASTIX_INT *ja )
{
  PASTIX_INT i1, i2, i;
  
  i1 = 0;
  i2 = 0;
  i = 0;
  while(i1 < n1 && i2 < n2)
    {
      if(ja1[i1] < ja2[i2])
	{
	  ja[i] = ja1[i1];
	  i1++;
	  i++;
	  continue;
	}

      if(ja1[i1] > ja2[i2])
	{
	  ja[i] = ja2[i2];
	  i2++;
	  i++;
	  continue;
	}
      
      ja[i] = ja1[i1];
      i++;
      i1++;
      i2++;
    }

#ifdef DEBUG_KASS
  assert(i1 == n1 || i2 == n2);
#endif

  for(;i1<n1;i1++)
    ja[i++] = ja1[i1];

  for(;i2<n2;i2++)
    ja[i++] = ja2[i2];

  *n = i;

#ifdef DEBUG_KASS
  assert(i >= n1 && i >= n2);
  assert(i <= n1+n2);
#endif

}



PASTIX_INT Merge_Cost(PASTIX_INT a, PASTIX_INT b, csptr P, PASTIX_INT *colweight)
{
  PASTIX_INT i1, n1, i2, n2;
  PASTIX_INT *ja1, *ja2;
  PASTIX_INT cost;




  ja1 = P->ja[a];
  ja2 = P->ja[b];


  n1 = P->nnzrow[a];
  n2 = P->nnzrow[b];

  i1 = i2 = 0;
  /** The diagonal elements of row a does not create fill-in **/ 
  while(i1 < n1 && ja1[i1] < ja2[0])
    i1++;


  /*fprintf(stderr, "MERGECOST %ld (%ld) + %ld (%ld)  i1 = %ld \n", a, n1, b, n2, i1);*/

  cost = 0;
  while(i1 < n1 && i2 < n2)
    {
      if(ja1[i1] < ja2[i2])
	{
	  cost += colweight[b];
	  /*fprintf(stderr, "TOTO cost1 %ld \n", cost);*/
	  i1++;
	  continue;
	}
      if(ja1[i1] > ja2[i2])
	{
	  cost += colweight[a];
	  /*fprintf(stderr, "TOTO cost2 %ld \n", cost);*/
	  i2++;
	  continue;
	}

      /** ja1[i1] == ja2[i2] **/
      i1++;
      i2++;
    }

  while(i1 < n1)
    {
      cost += colweight[b];
      /*fprintf(stderr, "TOTO costR1 %ld \n", cost);*/
      i1++;
      continue;
    }

  while(i2 < n2)
    {
      cost += colweight[a];
      /*fprintf(stderr, "TOTO costR2 %ld \n", cost);*/
      i2++;
      continue;
    }

  return cost;
}
