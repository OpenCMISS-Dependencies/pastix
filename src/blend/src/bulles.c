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
#include <math.h>
#include <assert.h>

#include "common_pastix.h"
#include "queue.h"
#include "bulles.h"

#define NBLEAVES btree->leavesnbr
#define NBBUBMAX btree->nodemax

/* 
   Function: Bubble_InitTree
   
   Initialize a tree of bubbles
   
   Parameters:
     btree        - Pointer to bubble tree to allocate and initialize.
     nbleaves     - Number of leaves in the tree.

   Returns:

*/
void Bubble_InitTree(BubbleTree * btree, int nbleaves)
{
  
  int i,n;
  int nblevel=2;
  int nbbubblesmax;

  n=nbleaves;
  while(n >>= 1)
    nblevel++;

  nbbubblesmax = (1 << nblevel)-1;
  print_debug(DBG_BUBBLES, "Bulles : nbleaves = %d, nblevel = %d et nbbubblesmax = %d\n", 
	      nbleaves, nblevel, nbbubblesmax); 

  btree->leavesnbr = nbleaves;
  btree->nodenbr   = nbleaves;
  btree->nodemax   = nbbubblesmax;
  MALLOC_INTERN(btree->sonstab, nbbubblesmax, int);
  MALLOC_INTERN(btree->nodetab, nbbubblesmax, BubbleTreeNode);
  
  for(i=0; i<NBLEAVES; i++)
    {
      btree->sonstab[i] = -1;
      btree->nodetab[i].fathnum   = -1;
      btree->nodetab[i].sonsnbr   = 0;
      btree->nodetab[i].fsonnum   = 0;
      btree->nodetab[i].fcandnum  = i;
      btree->nodetab[i].lcandnum  = i;
      btree->nodetab[i].costlevel = 0.0;
      btree->nodetab[i].treelevel = 0;
      btree->nodetab[i].priomin   = INTVALMAX;
      btree->nodetab[i].priomax   = 0;
      btree->nodetab[i].treelevel = 0;
      MALLOC_INTERN(btree->nodetab[i].taskheap, 1, Queue);
      queueInit(btree->nodetab[i].taskheap, 1000);
    }

  for(i=NBLEAVES; i<NBBUBMAX; i++)
    {
      btree->sonstab[i] = -1;
      btree->nodetab[i].fathnum   = -1;
      btree->nodetab[i].sonsnbr   = 0;
      btree->nodetab[i].fsonnum   = 0;
      btree->nodetab[i].fcandnum  = -1;
      btree->nodetab[i].lcandnum  = -1;
      btree->nodetab[i].costlevel = 0.0;
      btree->nodetab[i].treelevel = 0;
      btree->nodetab[i].priomin   = INTVALMAX;
      btree->nodetab[i].priomax   = 0;
      MALLOC_INTERN(btree->nodetab[i].taskheap, 1, Queue);
      queueInit(btree->nodetab[i].taskheap, 100);
    }
}

/* 
   Function: Bubble_Free
   
   Free the BubbleTree structure.
   
   Parameters:
     btree        - Pointer to the bubble tree to be free, NULL is returned.

   Returns:

*/
void Bubble_Free(BubbleTree *btree)
{
  int i;
  
  for (i=0; i<NBBUBMAX; i++)
    {
      queueExit(btree->nodetab[i].taskheap);
      memFree_null(btree->nodetab[i].taskheap);
    }

  memFree_null(btree->nodetab);
  memFree_null(btree->sonstab);
}

/* 
   Function: Bubble_Add
   
   Add a task in the BubbleTree structure defined by the set of candidates 
   processors (first and last) and by the cost and the deep of this task in the 
   dependence's tree.
   
   Parameters:
     btree        - Pointer to Bubble tree in which we need to add some task.
     fcandnum     - Indice of the first candidate for the task to add.
     lcandnum     - Indice of the last candidate for the task to add.
     costlevel    - Cost of the path from the root to this task.
     treelevel    - Deep of the task in the task's tree.

   Returns:
     The indice of the bubble corresponding to the set of candidates.
     Exit if the search or the add of the bubble fails.

*/
int Bubble_Add(BubbleTree * btree, 
	       PASTIX_INT          fcandnum, 
	       PASTIX_INT          lcandnum, 
	       double       costlevel, 
	       PASTIX_INT          treelevel){

  int i;

  /* It's a leaf */
  if( fcandnum == lcandnum )
    {
#ifdef BUBBLE_DEBUG
      errorPrintW("Bubble_Add : Ajout bulle processeur %d, costlevel %f, treelevel %d",
	      fcandnum, costlevel, treelevel);
#endif
      btree->nodetab[fcandnum].costlevel = MIN(btree->nodetab[fcandnum].costlevel, -costlevel);
      btree->nodetab[fcandnum].treelevel = MIN(btree->nodetab[fcandnum].treelevel, -treelevel);
      return fcandnum;
    }

  for (i=NBLEAVES; i<NBBUBMAX; i++)
    {
      /* The node already exists */
      if ((btree->nodetab[i].fcandnum == fcandnum) && (btree->nodetab[i].lcandnum == lcandnum))
	{
	  btree->nodetab[i].costlevel = MIN(btree->nodetab[i].costlevel, -costlevel);
	  btree->nodetab[i].treelevel = MIN(btree->nodetab[i].treelevel, -treelevel);
	  return i;
	}
      /* The node doesn't exist, it is created */
      else if (btree->nodetab[i].fcandnum == -1)
	{
	  btree->nodetab[i].fcandnum = fcandnum; 
	  btree->nodetab[i].lcandnum = lcandnum;
	  btree->nodetab[i].costlevel = -costlevel;
	  btree->nodetab[i].treelevel = -treelevel;
	  btree->nodenbr++;
	  return i;
	}
    }

  /* Don't never be here */
  errorPrint("La bulle [%ld %ld] n'a pas pu etre ajoutee\n", 
	     (long)fcandnum, (long)lcandnum);
  ASSERT(0==1, MOD_BUBBLE);
  return -1;
} 


/* 
   Function: Bubble_buildtree
   
   Construct the tree by filling nodetab and sonstab tables thanks to 
   the bubbles added whith Bubble_Add. This function must be call only 
   one time for one bubble tree. And only after the add of all bubbles.
   
   Parameters:
     btree        - Pointer to Bubble tree
     verbose      - If verbose, the tree is printed

   Returns:
     void
*/
void Bubble_BuildTree(const BubbleTree * btree){

  int i, j;
  int bubsize;
  int fprocnum, lprocnum, fathnum;
  int maxlevel    = 0;
  int *sonsnbrtmp = NULL;

  for (i=0; i<btree->nodenbr; i++)
    {
      bubsize  = NBLEAVES;
      fprocnum = btree->nodetab[i].fcandnum;
      lprocnum = btree->nodetab[i].lcandnum;

      for(j=NBLEAVES; j<btree->nodenbr; j++)
	{
	  /* The father node is the smallest bubble which contains it */
	  if ((i != j) && (fprocnum >= btree->nodetab[j].fcandnum)
	      && (lprocnum <= btree->nodetab[j].lcandnum)
	      && (bubsize > (btree->nodetab[j].lcandnum - btree->nodetab[j].fcandnum)))
	    {
	      bubsize  = btree->nodetab[j].lcandnum - btree->nodetab[j].fcandnum;
	      btree->nodetab[i].fathnum = j;
	    }
	}

      if (btree->nodetab[i].fathnum != -1)
	btree->nodetab[btree->nodetab[i].fathnum].sonsnbr++;      
    }
  
  for(i=NBLEAVES; i<btree->nodenbr; i++)
    maxlevel = MAX(maxlevel, btree->nodetab[i].treelevel);
  for(i=NBLEAVES; i<btree->nodenbr; i++)
    btree->nodetab[i].treelevel = maxlevel + 1 - btree->nodetab[i].treelevel;

  /* Compute the number of son and indice of first son */
  MALLOC_INTERN(sonsnbrtmp, btree->nodenbr, int);
  sonsnbrtmp[0] = 0;
  for (i=1; i<btree->nodenbr; i++)
    {
      btree->nodetab[i].fsonnum = btree->nodetab[i-1].fsonnum + btree->nodetab[i-1].sonsnbr;
      sonsnbrtmp[i] = 0;
    }
      
  /* Fill in sonstab */
  for (i=0; i<btree->nodenbr; i++)
    {
      fathnum = btree->nodetab[i].fathnum;
      if (fathnum != -1)
	{
	  btree->sonstab[btree->nodetab[fathnum].fsonnum + sonsnbrtmp[fathnum]] = i;
	  sonsnbrtmp[fathnum]++;
	  ASSERTDBG(sonsnbrtmp[fathnum] <= btree->nodetab[fathnum].sonsnbr, MOD_BLEND);
	}
    }
  memFree_null(sonsnbrtmp);
  
  return;
}


/* 
   Function: Bubble_print
   
   Generate a dot file for graphviz with the bubble tree.
   
   Parameters:
     btree        - Pointer to Bubble tree
     bcost        - Cost in time of all tasks included in each bubble

   Returns:
     void
*/
void Bubble_Print(const BubbleTree * btree,
		  const double     * bcost,
                  double totalcost, 
		  FILE *out)
{
  PASTIX_INT i;
  PASTIX_INT fprocnum, lprocnum;
  double percent;
  double *subtrees_costs = (double*)malloc( btree->nodenbr * sizeof(double) );
  memset( subtrees_costs, 0, btree->nodenbr * sizeof(double) );
  

  fprintf(out,"digraph G {\n"
	  "\tcolor=white\n"
	  "\trankdir=BT;\n");

  for (i=0; i<btree->nodenbr; i++) {
      int father = btree->nodetab[i].fathnum;
      subtrees_costs[i] += bcost[i];

      while ( father != -1 ) {
          subtrees_costs[father] += bcost[i];
          father = btree->nodetab[father].fathnum;
      }
  }

  for (i=0; i<btree->nodenbr; i++)
    {
      fprocnum = btree->nodetab[i].fcandnum;
      lprocnum = btree->nodetab[i].lcandnum;
      
      percent = subtrees_costs[i] / totalcost;

      if (btree->nodetab[i].fathnum == -1)
	fprintf(out, "\t%ld [style=filled, label=\"%ld-%ld\\n#tsk: %ld\\n"
                "Local Time: %lf, %.2lf\\nSubtree Time: %lf, %.2lf\",color=\"/set39/%d\"]\n",
		(long)i, (long)fprocnum, (long)lprocnum,
		(long)queueSize(btree->nodetab[i].taskheap), 
                bcost[i], bcost[i] / totalcost,
                subtrees_costs[i],
                subtrees_costs[i] / totalcost,
                ( percent > ( (double)( lprocnum-fprocnum+1 ) / (double)(btree->leavesnbr) ) ) ? 4 : 7 );
      else
	fprintf(out, "\t%ld -> %ld\n\t%ld [style=filled, label=\"%ld-%ld\\n#tsk : %ld\\nLocal Time: %lf, %.2lf\\nSubtree Time: %lf, %.2lf\",color=\"/set39/%d\"]\n",
		(long)i, (long)btree->nodetab[i].fathnum, 
		(long)i, (long)fprocnum, (long)lprocnum,
		(long)queueSize(btree->nodetab[i].taskheap), 
                bcost[i], bcost[i] / totalcost,
                subtrees_costs[i],
                subtrees_costs[i] / totalcost,
                ( percent > ( (double)( lprocnum-fprocnum+1 ) / (double)(btree->leavesnbr) ) ) ? 4 : 7 );
    }
      
  fprintf(out, "}\n");

  free(subtrees_costs);
}
