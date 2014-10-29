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
#include <strings.h>
#include <assert.h>
#include <math.h>
#include "common_pastix.h"
#include "dof.h"
#include "cost.h"
#include "ftgt.h"
#include "symbol.h"
#include "elimin.h"
#include "perf.h"
#include "cand.h"
#include "queue.h"
#include "bulles.h"
#include "extendVector.h"
#include "param_blend.h"
#include "blendctrl.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "symbol_cost.h"
#include "extrastruct.h"
#include "splitfunc.h"
#include "simu.h"
#include "costfunc.h"


#define BLEND_CHOLESKY /** LLt version **/

void   subtreeSetNullCost    (PASTIX_INT, const BlendCtrl * ctrl, const SymbolMatrix *, const SimuCtrl *,  PASTIX_INT);
double cblkComputeCost2DLocal(PASTIX_INT, const BlendCtrl * ctrl, const SymbolMatrix *, const Dof *, const SimuCtrl *);

/*+ Compute cost time  for each cblk and contribution of blok in the matrix +*/
void costMatrixBuild(CostMatrix *costmtx, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
  PASTIX_INT i;
  MALLOC_INTERN(costmtx->cblktab, symbmtx->cblknbr, CostCblk);
  MALLOC_INTERN(costmtx->bloktab, symbmtx->bloknbr, CostBlok);
  
  for(i=0;i<symbmtx->cblknbr;i++)
    cblkComputeCost(i, costmtx, symbmtx, dofptr);
}


void costMatrixCorrect(CostMatrix *costmtx, const SymbolMatrix *symbmtx, Cand * candtab, const Dof * dofptr)
{
  PASTIX_INT i;
  for(i=0;i<symbmtx->cblknbr;i++)
    if(candtab[i].distrib == D1)
      cblkComputeCost(i, costmtx, symbmtx, dofptr);
    else
      cblkComputeCost2D(i, costmtx, symbmtx, dofptr);
}


/*+ Summ the subtree cost node ; do not recompute node cost +*/  
double subtreeUpdateCost(PASTIX_INT rootnum, CostMatrix *costmtx, const EliminTree *etree)
{
  PASTIX_INT i;
  costmtx->cblktab[rootnum].subtree = costmtx->cblktab[rootnum].total;
  for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
    costmtx->cblktab[rootnum].subtree += subtreeUpdateCost(TSON(etree, rootnum, i), costmtx, etree);
  return costmtx->cblktab[rootnum].subtree;
}


/*+ Summ the subtree cost local node ; do not recompute node cost +*/  
double subtreeUpdateCostLocal(PASTIX_INT rootnum, const BlendCtrl * ctrl, const SymbolMatrix *symbmtx, 
			      const SimuCtrl *simuctrl, const Dof * dofptr,  PASTIX_INT clustnum)
{

  CostMatrix *costmtx = ctrl->costmtx;
  EliminTree *etree   = ctrl->etree;
  Cand       *candtab = ctrl->candtab;
  PASTIX_INT         i;

  /* Update cost of local task 1D */
  if (candtab[rootnum].distrib == D1)
    {
      /* Subtree is not local */
      if (candtab[rootnum].cluster != clustnum)
	{
	  costmtx->cblktab[rootnum].total = 0.0;
	}
      /* If not, we don't touch the cost */
    }
  /* 2D */
  else
    {
      /* Update cost */
      costmtx->cblktab[rootnum].total = cblkComputeCost2DLocal(rootnum, ctrl, symbmtx, dofptr, simuctrl);
    }
  
  costmtx->cblktab[rootnum].subtree = costmtx->cblktab[rootnum].total;

  if ((candtab[rootnum].fccandnum <= clustnum) &&
      (candtab[rootnum].lccandnum >= clustnum))
    {
      for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
	costmtx->cblktab[rootnum].subtree += subtreeUpdateCostLocal(TSON(etree, rootnum, i), ctrl, 
								    symbmtx, simuctrl, dofptr, clustnum);
    }
#ifdef DEBUG_BLEND
  else
    {
      for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
	subtreeSetNullCost(TSON(etree, rootnum, i), ctrl, 
			   symbmtx, simuctrl, clustnum);
    }
#endif

  /* Sort the sons by decreasing order */
  { 
      PASTIX_INT son, i, sonsnbr;
      double cumul_cost, soncost;
      Queue *queue_tree;

      sonsnbr = etree->nodetab[rootnum].sonsnbr;

      MALLOC_INTERN(queue_tree, 1, Queue);
      queueInit(queue_tree, sonsnbr);
      for(i=0;i<sonsnbr;i++)
      {
          son = TSON(etree, rootnum, i);
              
          /** Cost in the current subtree to be mapped **/
          cumul_cost = -ctrl->costmtx->cblktab[son].subtree;
          
          /* Cost of the root node in the subtree */
          soncost    = -ctrl->costmtx->cblktab[son].total;
          
          queueAdd2(queue_tree, son, cumul_cost, soncost);
      }

      for(i=0;i<sonsnbr;i++)
      {
          TSON(etree, rootnum, i) = queueGet(queue_tree);
      }
      queueExit(queue_tree);
      memFree(queue_tree);


      for(i=1;i<sonsnbr;i++)
      {
          assert( ctrl->costmtx->cblktab[TSON(etree, rootnum, i)].subtree 
                  <= ctrl->costmtx->cblktab[TSON(etree, rootnum, i-1)].subtree );
      }
  }


  return costmtx->cblktab[rootnum].subtree;
}

void subtreeSetNullCost(PASTIX_INT rootnum, const BlendCtrl * ctrl, 
			const SymbolMatrix *symbmtx, const SimuCtrl *simuctrl, 
			PASTIX_INT clustnum)
{
  CostMatrix *costmtx = ctrl->costmtx;
  EliminTree *etree   = ctrl->etree;
  PASTIX_INT         i;
  
  ASSERT(ctrl->candtab[rootnum].cluster != clustnum, MOD_BLEND);
  ASSERT(ctrl->proc2clust[simuctrl->blprtab[symbmtx->cblktab[rootnum].bloknum]] != clustnum, MOD_BLEND);
  
  costmtx->cblktab[rootnum].total   = 0.0;
  costmtx->cblktab[rootnum].subtree = 0.0;
  for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
    subtreeSetNullCost(TSON(etree, rootnum, i), ctrl, symbmtx, simuctrl, clustnum);
  
  return;
}

double cblkComputeCost2D(PASTIX_INT cblknum, CostMatrix *costmtx, const SymbolMatrix *symbptr, const Dof * dofptr)
{
  PASTIX_INT i, j;
  PASTIX_INT L, h, g;
  double cost = 0.0;

  L    = (symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum + 1);
  L   *= (dofptr)->noddval;
  cost = DIAGCost(L);
  for(i=symbptr->cblktab[cblknum].bloknum+1;i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
      h = symbptr->bloktab[i].lrownum - symbptr->bloktab[i].frownum + 1;
      h *= (dofptr)->noddval;
      cost += E1Cost(L, h);
      for(j=i;j<symbptr->cblktab[cblknum+1].bloknum;j++)
	{

	  g = symbptr->bloktab[j].lrownum - symbptr->bloktab[j].frownum + 1;
	  g *= (dofptr)->noddval;
	  cost += E2Cost(L, h, g);
#ifdef DEBUG_BLEND
	  ASSERT(L > 0,MOD_BLEND);
	  ASSERT(h > 0,MOD_BLEND);
	  ASSERT(g > 0,MOD_BLEND);
#endif
	}
    }
#ifdef DEBUG_BLEND
/*  ASSERT(cost >= 0,MOD_BLEND);*/
#endif
  costmtx->cblktab[cblknum].total = cost;
  return cost;
}
	
double cblkComputeCost2DLocal(PASTIX_INT cblknum, const BlendCtrl * ctrl, const SymbolMatrix *symbptr, 
                              const Dof * dofptr, const SimuCtrl *simuctrl)
{
  double      cost = 0.0;
  PASTIX_INT         i, j;
  PASTIX_INT         L, h, g;

  L  = (symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum + 1);
  L *= (dofptr)->noddval;
  
  /*  if (simuctrl->bloktab[symbptr->cblktab[cblknum].bloknum].tasknum != -1)*/
  if (ctrl->proc2clust[simuctrl->blprtab[symbptr->cblktab[cblknum].bloknum]] == ctrl->clustnum)
    cost = DIAGCost(L);
  
  for(i=symbptr->cblktab[cblknum].bloknum+1;i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
      if (ctrl->proc2clust[simuctrl->blprtab[i]] != ctrl->clustnum)
        continue;

      h     = symbptr->bloktab[i].lrownum - symbptr->bloktab[i].frownum + 1;
      h    *= (dofptr)->noddval;
      cost += E1Cost(L, h);

      for(j=i; j<symbptr->cblktab[cblknum+1].bloknum; j++)
        {
          g     = symbptr->bloktab[j].lrownum - symbptr->bloktab[j].frownum + 1;
          g    *= (dofptr)->noddval;
          cost += E2Cost(L, h, g);
        }
    }
  return cost;
}  

/*+ Compute cost of the cblk, return total cost +*/

/** Assure that cblkComputeCost and cblkCost compute the same things !!!! **/
double cblkComputeCost(PASTIX_INT cblknum, CostMatrix *costmtx, const SymbolMatrix *symbmtx, const Dof * dofptr)
{
    PASTIX_INT l, h, g;
    PASTIX_INT k;
#ifndef DOF_CONSTANT
    PASTIX_INT i;
#endif
    
    /** we need the height of cblk non empty lines  and the broadness 
      of the cbl to compute the local compute cost **/
#ifdef DOF_CONSTANT
    l = (symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1)*(dofptr)->noddval;
    /*l = (symbmtx->bloktab[symbmtx->cblktab[cblknum].bloknum].lrownum - symbmtx->bloktab[symbmtx->cblktab[cblknum].bloknum].frownum+ 1)*(dofptr)->noddval;*/
#else
    for(i=symbmtx->cblktab[cblknum].fcolnum;i<=symbmtx->cblktab[cblknum].lcolnum;i++)
      l+= noddDlt(dofptr, i);
#endif

    g = 0;
    for(k=symbmtx->cblktab[cblknum].bloknum;k<symbmtx->cblktab[cblknum+1].bloknum;k++)
      {
#ifdef  DOF_CONSTANT
	g += (symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1)*(dofptr)->noddval;
#else
	for(i=symbmtx->bloktab[k].frownum;i<=symbmtx->bloktab[k].lrownum;i++)
	  g+= noddDlt(dofptr, i);
#endif
      }


    costmtx->bloktab[symbmtx->cblktab[cblknum].bloknum].linenbr = g;
    
    /** retrieve diag height so let g be the odb non empty lines height **/
    g -= l;

    /** compute the local compute cost **/
    if(l!=0)
      {
	costmtx->cblktab[cblknum].compute = computeCost(l, g);
      }
    else
      costmtx->cblktab[cblknum].compute = 0;
    /** compute for each odb its contribution compute cost and add cost **/
    costmtx->cblktab[cblknum].send = 0;
    for(k=symbmtx->cblktab[cblknum].bloknum+1;k<symbmtx->cblktab[cblknum+1].bloknum;k++)
      {
#ifdef  DOF_CONSTANT
	h = (symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1)*(dofptr)->noddval;
#endif

	/* g is the odb lines number above this odb (odb lines include)*/
	costmtx->bloktab[k].linenbr     = g;
	/*if(l!=0 && h != 0 && g != 0)*/
	costmtx->bloktab[k].contrib     = contribCompCost(l, h, g);
	/*else
	  costmtx->bloktab[k].contrib     = 0;*/
	/*if(h != 0 && g != 0)*/
	costmtx->bloktab[k].contrib    += contribAddCost(h, g);

	costmtx->cblktab[cblknum].send += costmtx->bloktab[k].contrib;
	g -= h;
      }
    costmtx->cblktab[cblknum].total = costmtx->cblktab[cblknum].compute 
                                          + costmtx->cblktab[cblknum].send;

#ifdef DEBUG_BLEND
   {
      PASTIX_INT stride=0;
      double cost2=0;
     
      for(k=symbmtx->cblktab[cblknum].bloknum; k<symbmtx->cblktab[cblknum+1].bloknum;k++)
	stride += symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1;
      ASSERT( costmtx->bloktab[symbmtx->cblktab[cblknum].bloknum].linenbr == stride * dofptr->noddval,MOD_BLEND);
    
      ASSERT(costmtx->cblktab[cblknum].total > 0,MOD_BLEND);
      cost2 = cblkCost(symbmtx->cblktab[cblknum+1].bloknum -  symbmtx->cblktab[cblknum].bloknum, 
		       &(symbmtx->bloktab[symbmtx->cblktab[cblknum].bloknum]), 
		       dofptr);
      /* Values should be equals but we accept the machine computational error */
      ASSERT(costmtx->cblktab[cblknum].total - cost2 < 10e-15, 
	     MOD_BLEND);

    }
#endif


    return costmtx->cblktab[cblknum].total;
}


/*****************************************************************************************
 *      There are the cost functions of the compute phase of factorization algorithm     *
 *****************************************************************************************/
#ifndef BLEND_CHOLESKY
double computeCost(PASTIX_INT L, PASTIX_INT g_total)
{
  double total = 0;
  total =(double)(L*PERF_COPY(L)+ PERF_PPF(L) + PERF_TRSM(L, g_total) + L*PERF_SCAL(g_total)
		 + L*PERF_COPY(g_total)); 
  return (total>0)?total:0;
}
#else
double computeCost(PASTIX_INT L, PASTIX_INT g_total)
{
  double total = 0;
  total =(double)(PERF_POF(L) + PERF_TRSM(L, g_total)) ;
  return (total>0)?total:0;
}
#endif


double contribCompCost(PASTIX_INT L, PASTIX_INT h, PASTIX_INT g)
{
  double total = 0;
#ifdef DEBUG_BLEND
  ASSERT(L>0,MOD_BLEND);
  ASSERT(h>=0,MOD_BLEND);
  ASSERT(g>=0,MOD_BLEND);
#endif
  total = (double)(PERF_GEMM(g,h,L));
#ifdef DEBUG_BLEND
  /*
  if(total>1)
    return 0.99;*/
#endif
  
  return (total>0)?total:0;
}

double contribAddCost(PASTIX_INT h, PASTIX_INT g)
{
  double total = 0;
#ifdef DEBUG_BLEND
  ASSERT(h>=0,MOD_BLEND);
  ASSERT(g>0,MOD_BLEND);
#endif
  
  total = (double)(PERF_GEAM(g, h));

  return (total>0)?total:0;
}


double costFtgtSend(PASTIX_INT clustsrc, PASTIX_INT sync_comm_nbr, FanInTarget *ftgt, BlendCtrl *ctrl, const Dof * dofptr)
{
  PASTIX_INT ddl_coefnbr = 0;
  PASTIX_INT ddl_delta   = 0;

  if(clustsrc == ctrl->proc2clust[ftgt->infotab[FTGT_PROCDST]])
    return 0.0;


#ifdef DEBUG_BLEND
  ASSERT(clustsrc >= 0,MOD_BLEND);
  /*if(procsrc < 0)
    fprintf(stdout, "Procsrc %ld procdest %ld \n", (long)procsrc, (long)ftgt->infotab[FTGT_PROCDST]);*/
#endif
#ifdef DOF_CONSTANT
  ddl_delta   = (ftgt->infotab[FTGT_LCOLNUM]-ftgt->infotab[FTGT_FCOLNUM]+1)*dofptr->noddval;
/*  ddl_coefnbr = (ftgt->indtab[ftgt->infotab[FTGT_BLOKNBR]]*(dofptr->noddval))*ddl_delta;*/
  ddl_coefnbr = (ftgt->infotab[FTGT_LROWNUM] -  ftgt->infotab[FTGT_FROWNUM]+1)*ddl_delta*dofptr->noddval;
#else
  /** Oimbe Pas implemente **/
  fprintf(stderr, "costFtgtSend not implemented for the case dof non constant \n");
  EXIT(MOD_BLEND,NOTIMPLEMENTED_ERR);
#endif
#ifdef DEBUG_BLEND
  ASSERT(ddl_coefnbr > 0,MOD_BLEND);
#endif
  perfcluster2(clustsrc, ctrl->proc2clust[ftgt->infotab[FTGT_PROCDST]], sync_comm_nbr, ctrl->perfptr, ctrl);

  return (ctrl->perfptr->startup 
	  + ctrl->perfptr->bandwidth * (ddl_coefnbr*sizeof(double) + MAXINFO * sizeof(PASTIX_INT)));

}

double costFtgtAdd(FanInTarget *ftgt, const Dof * dofptr)
{
   PASTIX_INT ddl_delta   = 0;
   PASTIX_INT ddl_stride  = 0;
#ifdef DOF_CONSTANT
  ddl_delta   = (ftgt->infotab[FTGT_LCOLNUM]-ftgt->infotab[FTGT_FCOLNUM]+1)*dofptr->noddval;
  /*ddl_stride  = ftgt->indtab[ftgt->infotab[FTGT_BLOKNBR]]*(dofptr->noddval);*/
   ddl_stride = (ftgt->infotab[FTGT_LROWNUM] -  ftgt->infotab[FTGT_FROWNUM]+1)*dofptr->noddval;
#else
  /** Oimbe Pas implemente **/
    EXIT(MOD_BLEND,NOTIMPLEMENTED_ERR);
#endif
#ifdef DEBUG_BLEND
  ASSERT(ddl_stride>0,MOD_BLEND);
  ASSERT( ddl_delta>0,MOD_BLEND);
#endif
  return contribAddCost(ddl_stride, ddl_delta);
}

/**********************************************/
/*     Pour le 2D                             */
/**********************************************/
#ifndef BLEND_CHOLESKY
double DIAGCost(PASTIX_INT L)
{
  return (double)(L*PERF_COPY(L)+ PERF_PPF(L));
}
double E1Cost(PASTIX_INT L, PASTIX_INT g)
{
  return (double)(PERF_TRSM(L, g) + L*PERF_SCAL(g)
		 + L*PERF_COPY(g));
}
double E2Cost(PASTIX_INT L, PASTIX_INT h, PASTIX_INT g)
{
  return (double)(PERF_GEMM(g,h,L) + PERF_GEAM(g, h));
}
#else
double DIAGCost(PASTIX_INT L)
{
  return (double)(PERF_POF(L));
}
double E1Cost(PASTIX_INT L, PASTIX_INT g)
{
  return (double)(PERF_TRSM(L, g));
}
double E2Cost(PASTIX_INT L, PASTIX_INT h, PASTIX_INT g)
{
  return (double)(PERF_GEMM(g,h,L));
}
#endif
/*****************************************************************************************
 *                  END of cost functions                                                *
 *****************************************************************************************/


double cblkMaxCost(PASTIX_INT cblknbr, const CostMatrix *costmtx)
{
    PASTIX_INT i;
    double maxcost;
    maxcost = 0;
    for(i=0;i< cblknbr;i++)
	if(costmtx->cblktab[i].total > maxcost)
	    maxcost = costmtx->cblktab[i].total;
    return maxcost;
}



double totalCost(PASTIX_INT cblknbr, const CostMatrix *costmtx)
{
    PASTIX_INT i;
    double total=0;
    for(i=0;i<cblknbr;i++)
	total += costmtx->cblktab[i].total;
    return total;
}
    

double memorySpaceCost(const SolverMatrix *solvmtx)
{
  double space=0;
  space += solverSpaceCost(solvmtx);
  return space;
}


double solverSpaceCost(const SolverMatrix *solvmtx)
{
  double space=0;
  /*PASTIX_INT i;*/
  space += sizeof(double)*(solvmtx->coefnbr);

  space += sizeof(SolverCblk)*(solvmtx->cblknbr);
  space += sizeof(SolverBlok)*(solvmtx->bloknbr);
  space += sizeof(FanInTarget)*(solvmtx->ftgtnbr);
  space += sizeof(PASTIX_INT)*MAXINFO*(solvmtx->ftgtnbr);
/*  for(i=0;i<solvmtx->ftgtnbr;i++)
    space += sizeof(double) * solvmtx->ftgttab[i].indtab[solvmtx->ftgttab[i].infotab[FTGT_BLOKNBR]]
             * (solvmtx->ftgttab[i].infotab[FTGT_LCOLNUM] -solvmtx->ftgttab[i].infotab[FTGT_FCOLNUM] + 1) ;
  */
  return space;
}

double symbolSpaceCost(const SymbolMatrix *symbmtx)
{
  double space=0;
  space += sizeof(SymbolCblk)*(symbmtx->cblknbr+1);
  space += sizeof(SymbolBlok)*(symbmtx->bloknbr+1);
  return space;
}


void printSolverInfo(FILE *out, const SolverMatrix * solvmtx, const SymbolMatrix * symbmtx, const Dof * const dofptr)
{
  PASTIX_INT procnum = solvmtx->clustnum;
  double totalspace = memorySpaceCost(solvmtx);
  fprintf(out,   " %ld : Number of operations             : %g \n", (long)procnum,
	  recursive_sum(0, symbmtx->cblknbr-1, crout_blok, symbmtx, dofptr));
  fprintf(out,   " %ld : Number of Column Bloks           : %ld \n", (long)procnum, (long)symbmtx->cblknbr);
  fprintf(out,   " %ld : Number of Bloks                  : %ld \n", (long)procnum, (long)symbmtx->bloknbr);
  fprintf(out,   " %ld : Number of Non Null Coeff         : %ld --> %g Octs \n",
	  (long)procnum, (long)solvmtx->coefnbr, (double)solvmtx->coefnbr*sizeof(double));
  fprintf(out,   " %ld : ExtraStructure Memory Space      : %g Octs \n", (long)procnum, totalspace - (double)solvmtx->coefnbr*sizeof(double));
  fprintf(out,   " %ld : Total Memory space               : %g  Octs\n", (long)procnum, totalspace); 
}




