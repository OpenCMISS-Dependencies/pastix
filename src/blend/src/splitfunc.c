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
#include "symbol.h"
#include "cost.h"
#include "extrastruct.h"
#include "dof.h"
#include "param_blend.h"
#include "elimin.h"
#include "cand.h"
#include "queue.h"
#include "extendVector.h"
#include "bulles.h"
#include "blendctrl.h"
#include "ftgt.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "simu.h"
#include "costfunc.h"
/* #include "assembly.h" */
/* #include "eliminfunc.h" */
/* #include "splitpart.h" */
/* #include "assert.h" */
#include "splitfunc.h"

PASTIX_INT splitSeqCblk2D( PASTIX_INT cblknum, PASTIX_INT procnbr, const SymbolMatrix * symbptr, const ExtraSymbolMatrix * extrasymbptr, const Dof * dofptr, const BlendCtrl * ctrl,
                    PASTIX_INT (*P)(PASTIX_INT, PASTIX_INT, const SymbolMatrix *,const ExtraSymbolMatrix *, const Dof *, PASTIX_INT, const BlendCtrl *) )
{
  PASTIX_INT a, c;
  a = 1;
  c = symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum + 1;

  /*(!P(cblknum, procnbr, symbptr, extrasymbptr, dofptr, 1, ctrl))
    return 1;*/
  if( P(cblknum, procnbr, symbptr, extrasymbptr, dofptr, c, ctrl))
    return c;

  while(a<c-1)
    {
      if(P(cblknum, procnbr, symbptr, extrasymbptr, dofptr, (a+c)/2, ctrl))
        a = (a+c)/2;
      else
        c = (a+c)/2;
    }
#ifdef DEBUG_BLEND
  ASSERT(!P(cblknum, procnbr, symbptr, extrasymbptr, dofptr, c, ctrl),MOD_BLEND);
  ASSERT(a>=1 && a<= symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum + 1,MOD_BLEND);
#endif

  return a;
}

#if 0
/** Basic criterion based on broadness of cblk **/
PASTIX_INT P1D(PASTIX_INT cblknum, PASTIX_INT procnbr, const SymbolMatrix *symbptr, const ExtraSymbolMatrix * const extrasymbptr, const Dof * dofptr, PASTIX_INT nseq, const BlendCtrl * ctrl)
{
  PASTIX_INT delta = symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum + 1;
  if( dofptr->noddval*delta/nseq >= ctrl->option->blcolmin )
    return 1;
  return 0;
}

/** Criterion based on overhead du to the splitting **/
PASTIX_INT P2D(PASTIX_INT cblknum, PASTIX_INT procnbr, const SymbolMatrix * symbptr, const ExtraSymbolMatrix * const extrasymbptr, const Dof * dofptr, PASTIX_INT nseq, const BlendCtrl * ctrl)
{
  PASTIX_INT bloknbr;
  PASTIX_INT sbloknbr;
  SymbolBlok *bloktab  = NULL;
  SymbolBlok *sbloktab = NULL;
  PASTIX_INT        *indtab  = NULL;
  double cost1;
  double cost2;


  PASTIX_INT k;

  /* Number of blocks */
  bloknbr = cblkNbr(cblknum, symbptr, extrasymbptr);

  MALLOC_INTERN(bloktab, bloknbr, SymbolBlok);
  build_cblk(cblknum, symbptr, extrasymbptr, bloktab);

  /* Number of splitted blocks */
  sbloknbr = (nseq*(nseq+1))/2 + (bloknbr-1)*nseq;


  /* Cost of the non splitted cblk */
  cost1 = ctrl->costmtx->cblktab[cblknum].total;
#ifdef DEBUG_BLEND
  if(cost1 != cblkCost(bloknbr, bloktab, dofptr))
    fprintf(stderr, "cblknum %ld ; nseq %ld bloknbr %ld costmtx %g mycost %g \n", (long)cblknum,  (long)nseq, (long)bloknbr, cost1, cblkCost(bloknbr, bloktab, dofptr));
/*  ASSERT(cost1 == cblkCost(bloknbr, bloktab, dofptr),MOD_BLEND);*/
#endif

  /** give the structure of splitted cblk (it is not repercuted on the symbol matrix) **/
  MALLOC_INTERN(indtab,   nseq+1,   PASTIX_INT);
  MALLOC_INTERN(sbloktab, sbloknbr, SymbolBlok);

  virtualSplit(nseq, bloknbr, bloktab, indtab, sbloktab);

  /* Cost of splitted cblk */
  cost2 = 0;
  for(k=0;k<nseq;k++)
    {
      cost2 += cblkCost(bloknbr+nseq-k-1, &(sbloktab[indtab[k]]), dofptr);
    }

  memFree(bloktab);
  memFree(indtab);
  memFree(sbloktab);
  fprintf(stdout, "cblknum %ld NSEQ %ld COST2 : %g COST1 : %g Cost2/Cost1 %g\n",  (long)cblknum,(long)nseq, cost2, cost1, cost2/cost1);
  if(cost1==0)
    return 0;

  procnbr = MIN(procnbr, nseq/4);

  if(cost2/cost1 < 0.5*procnbr)
    return 1;
  return 0;

}
#endif

/** OIMBE TO DO: le calcul total du cout splitter peut etre bp plus rapide en groupant
  les calculs sur les odb dans les cblk decouper **/
void virtualSplit(PASTIX_INT nseq, PASTIX_INT bloknbr, const SymbolBlok * src_bloktab,
                         PASTIX_INT *indtab, SymbolBlok * dest_bloktab)
{
  PASTIX_INT k;
  PASTIX_INT i;
  PASTIX_INT step;

  /** Index of splitted cblk in dest_bloktab **/
  /* + a virtual extra end cblk */

  indtab[0] = 0;
  for(k=1;k<=nseq;k++)
    indtab[k] = indtab[k-1] + bloknbr + nseq-k;

  step = (src_bloktab[0].lrownum - src_bloktab[0].frownum + 1)/nseq;

  /** Compute the splitted cblk **/
  for(k=0;k<nseq;k++)
    {
      /* We just have to compute diagonal splitting */
      for(i=0;i < nseq-k-1;i++)
        {
          dest_bloktab[i+indtab[k]].frownum = src_bloktab[0].frownum + (i+k)*step;
          dest_bloktab[i+indtab[k]].lrownum = src_bloktab[0].frownum + (i+k+1)*step-1;
        }
      dest_bloktab[i+indtab[k]].frownum = src_bloktab[0].frownum + (i+k)*step;
      dest_bloktab[i+indtab[k]].lrownum = src_bloktab[0].lrownum;
      i++;
      /* The odbs remain the same */
      if(bloknbr>1)
        memCpy(&(dest_bloktab[i+indtab[k]]), &(src_bloktab[1]), sizeof(SymbolBlok)*(bloknbr-1));
    }
  /*fprintf(stdout, "CBlok Orig \n");
  for(i=0;i<bloknbr;i++)
    printf("O [%ld %ld]\n",src_bloktab[i].frownum,src_bloktab[i].lrownum );

  for(k=0;k<nseq;k++)
    {
      printf("CBlok %ld \n", k);
      for(i=indtab[k];i<indtab[k+1];i++)
        printf("S [%ld %ld]\n",dest_bloktab[i].frownum,dest_bloktab[i].lrownum );
    }*/

}

/** Assure that cblkComputeCost and cblkCost compute the same things !!!! **/
double cblkCost(PASTIX_INT bloknbr, const SymbolBlok * bloktab, const Dof * dofptr)
{
    PASTIX_INT l, h, g;
    PASTIX_INT k;
    double total_cost = 0;
    double compute_cost = 0;
    double send_cost    = 0;
    double contrib_cost = 0;
    /** we need the height of cblk non empty lines  and the broadness
      of the cblk to compute the local compute cost **/
#ifdef DOF_CONSTANT
    l = (bloktab[0].lrownum -bloktab[0].frownum + 1)*(dofptr)->noddval;
#endif

    g = 0;
    for(k=0;k<bloknbr;k++)
      {
#ifdef  DOF_CONSTANT
        g += (bloktab[k].lrownum - bloktab[k].frownum + 1)*(dofptr)->noddval;
#endif
      }

    /** retrieve diag height so let g be the odb non empty lines height **/
    g -= l;

    /** compute the local compute cost **/
    if(l!=0)
        compute_cost += computeCost(l, g);
    else
      compute_cost = 0;

    /** compute for each odb its contribution compute cost and add cost **/
    for(k=1;k<bloknbr;k++)
      {
#ifdef  DOF_CONSTANT
        h = (bloktab[k].lrownum - bloktab[k].frownum + 1)*(dofptr)->noddval;
#endif
        /* g is the odb lines number above this odb (odb lines include)*/
        /*if(l!=0 && h != 0 && g != 0)*/
        contrib_cost     = contribCompCost(l, h, g);
        /*if(h != 0 && g != 0)*/
        contrib_cost     += contribAddCost(h, g);
        send_cost += contrib_cost;
        g -= h;
      }
    total_cost = compute_cost + send_cost;
    return total_cost;
}

PASTIX_INT cblkNbr(PASTIX_INT cblknum,  const SymbolMatrix * symbptr, const ExtraSymbolMatrix * extrasymbptr)
{
  PASTIX_INT bloknbr = 0;
  PASTIX_INT i;
  for(i=symbptr->cblktab[cblknum].bloknum;i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
      if(extrasymbptr->sptblnb[i]>1)
        bloknbr += extrasymbptr->sptblnb[i];
      else
        bloknbr++;
    }
  return bloknbr;
}

void build_cblk(PASTIX_INT cblknum, const SymbolMatrix * symbptr, const ExtraSymbolMatrix * extrasymbptr, SymbolBlok *bloktab)
{
  PASTIX_INT i;
  PASTIX_INT blokcur = 0;
  PASTIX_INT sptnbr;
  for(i=symbptr->cblktab[cblknum].bloknum;i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
      if(extrasymbptr->sptblok[i]>=0)
        {
          sptnbr = extrasymbptr->sptblnb[i];
#ifdef DEBUG_BLEND
          ASSERT(sptnbr>0,MOD_BLEND);
#endif
          memCpy(&(bloktab[blokcur]), &(extrasymbptr->bloktab[extrasymbptr->sptblok[i]]), sptnbr*sizeof(SymbolBlok));
          blokcur += sptnbr;
        }
      else
        {
          bloktab[blokcur].frownum = symbptr->bloktab[i].frownum;
          bloktab[blokcur].lrownum = symbptr->bloktab[i].lrownum;
          blokcur++;
        }
    }

}
