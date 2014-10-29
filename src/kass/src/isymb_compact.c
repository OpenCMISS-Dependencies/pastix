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
#include "common_pastix.h"
/* #include "param_comm.h" */
#ifdef WITH_SCOTCH
#ifdef DISTRIBUTED
#include <mpi.h>
#include "ptscotch.h"
#else
#include "scotch.h"
#endif
#endif
#include "graph.h"
#include "cost.h"
#include "ftgt.h"
#include "symbol.h"
#include "solver.h"
#include "dof.h"
#include "elimin.h"
#include "extrastruct.h"
#include "perf.h"
#include "symbol_cost.h"
#include "queue.h"
#include "extendVector.h"
#include "cand.h"
#include "param_blend.h"
#include "blendctrl.h"
#include "splitfunc.h"
#include "costfunc.h"
#include "isymb_compact.h"

#define FATHER(k)   (symbmtx->bloktab[symbmtx->cblktab[k].bloknum+1].cblknum);
#define MERGED -13
#define END_COMPACT -14

double symbol_nnz(SymbolMatrix *symbmtx);

void isymb_compact(SymbolMatrix *symbmtx, double g_ratio)
{
  /*******************************************************************************/
  /*  This function merges some column-block in the symbol                       */
  /* matrix to decrease the factorization time.                                  */
  /* The function is allowed to create new non-zero term in the incomplete       */
  /* block matrix but it a to respect the                                        */
  /*g_ratio = (NNZ in the new symbol matrix) / (NNZ in the initial symbol matrix */
  /*******************************************************************************/
  PASTIX_INT i, j, k;
  PASTIX_INT root;
  EliminTree *T;
  MALLOC_INTERN(T, 1, EliminTree);
  treeInit(T);

  /** Build the elimination tree from the symbolic partition **/
  etreeBuild(T, symbmtx);
  root = ROOT(T); /** The root may become another node if a son of a root is merged with the root **/

  /*** Find a new partition of supernode by merging some supernodes with their fathers ***/
  subtree_compact(root, 0, T, symbmtx);

  

  etreeExit(T);
}

double subtree_compact(PASTIX_INT node, double gmax, EliminTree *T, SymbolMatrix *symbmtx) 
{
  double gnode;
  PASTIX_INT i, max_son;
  PASTIX_INT gi_max;

  do
    {
      gnode = compute_gain(node, T, symbmtx);
      if(gnode > gmax)
	gmax = gnode;
      
      gi_max = 0;
      for(i=0;i<T->nodetab[node].sonsnbr;i++)
	{
	  PASTIX_INT son;
	  son = TSON(T, node, i);
	  if(son == MERGED)
	    continue; /** This son does not exist anymore **/
	  
	  gi = subtree_compact(son, gmax, T, symbmtx);
	  if(gi == END_COMPACT)
	    return END_COMPACT;

	  if(gi >= gi_max)
	    {
	      max_son = son;
	      gi_max = gi;
	    }
	}
      if(gi_max >= gmax && gi_max > 0)
	{
	  nnz_add += merge_node(max_son, node, T, symbmtx);
	  if(nnz_add > nnz_add_limit)
	    return END_COMPACT;
	}

    }while(gi_max > gmax);
  
  return gnode;
}



double symbol_nnz(SymbolMatrix *symbmtx)
{
  recursive_sum(0, symbmtx->cblknbr-1, nnz, symbmtx, dofptr); 
}





