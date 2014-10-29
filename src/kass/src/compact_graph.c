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
/**   NAME       : compact_graph.c                         **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10/08/2006      **/
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
#include "compact_graph.h"

PASTIX_INT is_col_match(PASTIX_INT n, PASTIX_INT *ja1, PASTIX_INT *ja2);


void compact_graph(csptr P, PASTIX_INT *snodenbr, PASTIX_INT *snodetab, PASTIX_INT *iperm)
{
  /************************************************/
  /* This function compress a graph: all column   */
  /* of the graph that have the same pattern are  */
  /* merged                                       */
  /* The compressed graph is compute in place (P) */
  /* the set of unknowns that are in vertex i are */
  /* nodelist[snodetab[i]:snodetab[i+1][          */
  /* On entry/return:                             */
  /* P: a symmetric graph (with the diagonal !)   */
  /* iperm : the permutation that was             */
  /* On return only                               */
  /* snodetab, iperm                              */
  /* vwgt: the vertex weight (number of unknowns) */
  /************************************************/
  PASTIX_INT i, j, k;
  PASTIX_INT *colflag = NULL;
  PASTIX_INT cnt;
  PASTIX_INT *ja;

#ifdef DEBUG_KASS
  for(i=0;i<P->n;i++)
    {
      k = 0;
      /** Check the row are sorted **/
      for(j=1;j<P->nnzrow[i];j++)
        assert(P->ja[i][j]>P->ja[i][j-1]);

      /** Check the diagonal is here **/
      for(j=0;j<P->nnzrow[i];j++)
        if(P->ja[i][j] == i)
          k = 1;
      assert(k>0);
    }
#endif

  MALLOC_INTERN(colflag, P->n, PASTIX_INT);

  bzero(colflag, sizeof(PASTIX_INT)*P->n);
  cnt = 0;
  for(k=0;k<P->n;k++)
    {
      if(colflag[k] != 0)
        continue;
      /** Find the set of column that can be merged
          in the set of columns that have a non-zero
          in the row k (NB the graph is symmetric !) **/
      ja = P->ja[k];
      for(i=0;i<P->nnzrow[k];k++)
        {
          j = ja[i];
          if(j != i && colflag[j] == 0 && P->nnzrow[j] == P->nnzrow[i])
            if(is_col_match(P->nnzrow[i], ja, P->ja[j]) != 0)
              {
                colflag[j] = 1;
                cnt++;
              }
        }
    }
  fprintf(stderr, "Number of compression %ld cblk = %ld \n", (long)cnt, (long)(P->n-cnt));


  memFree(colflag);

}


PASTIX_INT is_col_match(PASTIX_INT n, PASTIX_INT *ja1, PASTIX_INT *ja2)
{
  PASTIX_INT i;
  for(i=0;i<n;i++)
    if(ja1[i] != ja2[i])
      return 0;


  return 1;
}
