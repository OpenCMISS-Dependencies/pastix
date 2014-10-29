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
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common_pastix.h"
#include "dof.h"
#include "cost.h"
#include "symbol.h"
#include "elimin.h"
#include "extrastruct.h"
#include "extendVector.h"
#include "cand.h"
#include "param_blend.h"
#include "queue.h"
#include "bulles.h"
#include "blendctrl.h"
#include "ftgt.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "simu.h"
#include "costfunc.h"
#include "partbuild.h"
#include "eliminfunc.h"
#include "splitpart.h"

#define CROSS_TOLERANCE 0.1

/*
  Function: splitPart

  Repartitioning of the initial symbolic factorization
  and processing of candidate processors group for
  each colum bloc

  Parameters:
    symbmtx - Symbolic matrix.
    ctrl    -
    dofptr  -
*/
void splitPart(SymbolMatrix *symbmtx,
               BlendCtrl    *ctrl,
               const Dof    *dofptr)
{
    PASTIX_INT i;
    ExtraSymbolMatrix *extrasymb;
    ExtraCostMatrix   *extracost;

    MALLOC_INTERN(extrasymb, 1, ExtraSymbolMatrix);
    extrasymbolInit(extrasymb);

    MALLOC_INTERN(extracost, 1, ExtraCostMatrix);
    extracostInit(extracost);

    /* initialize candtab */
    MALLOC_INTERN(ctrl->candtab, symbmtx->cblknbr, Cand);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
            ctrl->candtab[i].costlevel = 0.0;
            ctrl->candtab[i].treelevel = 0;
            ctrl->candtab[i].fcandnum  = -1;
            ctrl->candtab[i].lcandnum  = -1;
            ctrl->candtab[i].fccandnum = -1;
            ctrl->candtab[i].lccandnum = -1;
            ctrl->candtab[i].distrib   = D1;
            ctrl->candtab[i].cluster   = -1;
        }

    /** set tree level of each cblk **/
    /** OIMBE le faire apres partbuild **/
    setTreeLevel(ctrl->candtab, ctrl->etree);
    if(ctrl->option->costlevel)
      setTreeCostLevel(ctrl->candtab, ctrl->etree, ctrl->costmtx);

    /* initialize spt[tab] */
    MALLOC_INTERN(extrasymb->sptcblk,      symbmtx->cblknbr, PASTIX_INT);
    MALLOC_INTERN(extrasymb->sptcbnb,      symbmtx->cblknbr, PASTIX_INT);
    MALLOC_INTERN(extrasymb->sptblok,      symbmtx->bloknbr, PASTIX_INT);
    MALLOC_INTERN(extrasymb->sptblnb,      symbmtx->bloknbr, PASTIX_INT);
    MALLOC_INTERN(extrasymb->subtreeblnbr, symbmtx->cblknbr, PASTIX_INT);

    /* set spt* array to -1 ,
       positive or null value means cblk/blok has been splitted
       the value is, in this case the index in the extra- cblk/blok -tab */
    for(i=0;i<symbmtx->cblknbr;i++)
        {
            extrasymb->sptcblk[i] = -1;
            extrasymb->sptcbnb[i] =  1;
        }

    for(i=0;i<symbmtx->bloknbr;i++)
        extrasymb->sptblok[i] = -1;
    bzero(extrasymb->sptblnb, symbmtx->bloknbr * sizeof(PASTIX_INT));

    extrasymb->curcblk = 0;
    extrasymb->curblok = 0;

    /* We choose an arbitrary size for initial allocation of bloktab and cblktab */
    MALLOC_INTERN(extrasymb->cblktab, symbmtx->cblknbr/3 + 1, SymbolCblk);
    extrasymb->sizcblk = symbmtx->cblknbr/3 + 1;
    MALLOC_INTERN(extrasymb->bloktab, symbmtx->bloknbr/3 + 1, SymbolBlok);
    extrasymb->sizblok = symbmtx->bloknbr/3 + 1;

    MALLOC_INTERN(extracost->cblktab, symbmtx->cblknbr/3 + 1, CostCblk);
    MALLOC_INTERN(extracost->bloktab, symbmtx->bloknbr/3 + 1, CostBlok);


    if((ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO) &&
       (ctrl->option->leader == ctrl->clustnum))
      fprintf(stdout, "   Using proportionnal mapping \n");

    if(ctrl->option->candcorrect)
      {
        propMappTreeNoSplit(symbmtx, ctrl, dofptr);
        setDistribType(symbmtx->cblknbr, symbmtx, ctrl->candtab, ctrl->option->level2D);
        costMatrixCorrect(ctrl->costmtx, symbmtx, ctrl->candtab, dofptr);
        subtreeUpdateCost(ROOT(ctrl->etree), ctrl->costmtx, ctrl->etree);
      }

    /* Repartition symbolic matrix using proportionnal mapping method */
    propMappTree(symbmtx, extrasymb, extracost, ctrl, dofptr);

    ASSERTDBG(check_candidat(symbmtx, ctrl) >= 0,MOD_BLEND);

    /** set the distribution type (1D or 2D) for cblk **/
    if(ctrl->option->autolevel)
      {
        /*setSubtreeBlokNbr(ROOT(ctrl->etree), ctrl->etree, symbmtx, extrasymb, ctrl->option->blcolmin);
          fprintf(stdout, "Total blok %ld \n", extrasymb->subtreeblnbr[ROOT(ctrl->etree)]);*/
        setSubtreeDistribType(symbmtx, ctrl->costmtx, ROOT(ctrl->etree), ctrl, D2);
      }
    else
      setDistribType(symbmtx->cblknbr, symbmtx, ctrl->candtab, ctrl->option->level2D);

    /** For a dense end block **/
    if(ctrl->option->dense_endblock)
      ctrl->candtab[symbmtx->cblknbr-1].distrib = DENSE;

    /* Rebuild the symbolic matrix */
    partBuild(symbmtx, extrasymb, ctrl->costmtx, extracost, ctrl, dofptr);

    ASSERTDBG(check_candidat(symbmtx, ctrl) >= 0,MOD_BLEND);

    if(ctrl->option->candcorrect)
      {
        treeExit(ctrl->etree);
        MALLOC_INTERN(ctrl->etree, 1, EliminTree);
        treeInit(ctrl->etree);
        eliminTreeBuild(symbmtx, ctrl);

        costMatrixCorrect(ctrl->costmtx, symbmtx, ctrl->candtab, dofptr);
        subtreeUpdateCost(ROOT(ctrl->etree), ctrl->costmtx, ctrl->etree);
        propMappTreeNoSplit(symbmtx, ctrl, dofptr);
        setDistribType(symbmtx->cblknbr, symbmtx, ctrl->candtab, ctrl->option->level2D);
      }

    /**************************/
    /** Reset the tree level **/
    /**************************/
    treeExit(ctrl->etree);
    MALLOC_INTERN(ctrl->etree, 1, EliminTree);
    treeInit(ctrl->etree);
    eliminTreeBuild(symbmtx, ctrl);
    setTreeLevel(ctrl->candtab, ctrl->etree);
    if(ctrl->option->costlevel)
      setTreeCostLevel(ctrl->candtab, ctrl->etree, ctrl->costmtx);

    extrasymbolExit(extrasymb);
    extracostExit(extracost);

#ifdef SEB
    {
      SCOTCH_Graph    graphptr;
      SCOTCH_Mapping  mappptr;
      SCOTCH_Arch     archptr;
      SCOTCH_Num     *parttab;
      SCOTCH_Strat    stratptr;
      FILE           *archfile;
      PASTIX_INT             vertnbr;

      /* Construction du graph scotch */
      tree2graph(ctrl->etree, ctrl->costmtx, &graphptr);

      SCOTCH_graphData(&graphptr, NULL, &vertnbr, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
      ASSERT((ctrl->etree->nodenbr == vertnbr), MOD_BLEND);

      MALLOC_INTERN(parttab, ctrl->etree->nodenbr, PASTIX_INT);

      /* Init strat */
      SCOTCH_stratInit(&stratptr);    /* Strategy for sub-architecture computation */
      SCOTCH_stratGraphMap (&stratptr, "b{job=t,map=t,poli=S,strat=m{type=h,vert=80,low=h{pass=10}f{bal=0.0005,move=80},asc=b{bnd=d{pass=40}f{bal=0.005,move=80},org=f{bal=0.005,move=80}}}|m{type=h,vert=80,low=h{pass=10}f{bal=0.0005,move=80},asc=b{bnd=d{pass=100}f{bal=0.005,move=80},org=f{bal=0.005,move=80}}}}");

      /* Arch init */
      SCOTCH_archInit(&archptr);
      archfile = fopen("arch.tgt", "r");
      SCOTCH_archLoad(&archptr, archfile);

      /* Appel a ta fonction */
      SCOTCH_graphMapInit(&graphptr, &mappptr, &archptr, parttab);

      SCOTCH_graphMapCompute(&graphptr, &mappptr, &stratptr);

      /* Ecrasement des candidats */
      SetCandtab(ctrl->candtab, parttab, ctrl->etree->nodenbr);

    }
#endif

    /* printTree(stderr, ctrl->etree, ROOT(ctrl->etree)); */

    /*------------------------------------------------/
    /    Set the cluster candidat according to the    /
    /    processor candidat                           /
    /------------------------------------------------*/
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        ctrl->candtab[i].fccandnum =   ctrl->proc2clust[ctrl->candtab[i].fcandnum] ;
        ctrl->candtab[i].lccandnum =   ctrl->proc2clust[ctrl->candtab[i].lcandnum] ;
      }


#ifdef DEBUG_BLEND
    if(ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, " New Cost of the Matrix %g \n", totalCost(symbmtx->cblknbr, ctrl->costmtx));
#endif
}

void printTree(FILE *stream, const EliminTree *etree, PASTIX_INT rootnum)
{
  int i, sonsnbr;
  PASTIX_INT son;

  sonsnbr = etree->nodetab[rootnum].sonsnbr;

  fprintf(stream, "Rootnum %ld %d\n", (long)rootnum, sonsnbr);
  for(i=0;i<sonsnbr;i++)
    fprintf(stream,"       (%4ld)\n",  (long)TSON(etree, rootnum, i));
  /*     fprintf(stdout, "\n"); */

  for(i=0;i<sonsnbr;i++)
    {
      son = TSON(etree, rootnum, i);
      if (etree->nodetab[son].sonsnbr)
        printTree(stream, etree, son);
    }
}

void setTreeLevel(Cand *candtab, const EliminTree *etree)
{
  candtab[ROOT(etree)].treelevel = -1;
  setSubtreeLevel(ROOT(etree), candtab, etree);
}

void setSubtreeLevel(PASTIX_INT rootnum, Cand *candtab, const EliminTree *etree)
{
  PASTIX_INT i;
  for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
    {
      candtab[TSON(etree, rootnum, i)].treelevel = candtab[rootnum].treelevel -1;
      setSubtreeLevel(TSON(etree, rootnum, i), candtab, etree);
    }
}





void setTreeCostLevel(Cand *candtab, const EliminTree *etree, const CostMatrix *costmtx)
{
  candtab[ROOT(etree)].costlevel = -1.0;
  setSubtreeCostLevel(ROOT(etree), candtab, etree, costmtx);
}

void setSubtreeCostLevel(PASTIX_INT rootnum, Cand *candtab, const EliminTree *etree, const CostMatrix *costmtx)
{
  PASTIX_INT i;
  for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
    {
      candtab[TSON(etree, rootnum, i)].costlevel = candtab[rootnum].costlevel - costmtx->cblktab[rootnum].total;
      setSubtreeCostLevel(TSON(etree, rootnum, i), candtab, etree, costmtx);
    }
}


void setDistribType(const PASTIX_INT cblknbr, SymbolMatrix *symbptr, Cand *candtab, const PASTIX_INT level2D)
{
  PASTIX_INT i;
  (void)symbptr;

  for(i=0;i<cblknbr;i++)
    {
      if(candtab[i].fcandnum != candtab[i].lcandnum)
        {
          if( (-candtab[i].treelevel <= level2D))
            candtab[i].distrib = D2;
          else
            candtab[i].distrib = D1;
        }
      else
        candtab[i].distrib = D1;
    }
}

PASTIX_INT PERF_limit_2D(PASTIX_INT procnbr, PASTIX_INT bloklimit)
{
  return bloklimit*procnbr;
}

void setSubtreeDistribType(const SymbolMatrix *symbptr, const CostMatrix *costmtx, PASTIX_INT rootnum, const BlendCtrl *ctrl, PASTIX_INT distrib_type)
{
  PASTIX_INT i;

  if(distrib_type == D1)
    {
      ctrl->candtab[rootnum].distrib = D1;
      for(i=0;i<ctrl->etree->nodetab[rootnum].sonsnbr;i++)
        setSubtreeDistribType(symbptr, costmtx, TSON(ctrl->etree, rootnum, i), ctrl, D1);
    }
  else
    {
      PASTIX_INT candnbr; /* , bloknbr; */

      candnbr = ctrl->candtab[rootnum].lcandnum - ctrl->candtab[rootnum].fcandnum+1;

      if(candnbr > (PASTIX_INT)ctrl->option->ratiolimit)
         {
          ctrl->candtab[rootnum].distrib = D2;
          for(i=0;i<ctrl->etree->nodetab[rootnum].sonsnbr;i++)
            setSubtreeDistribType(symbptr, costmtx, TSON(ctrl->etree, rootnum, i),  ctrl, D2);
        }
      else
        {
          if (ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
            fprintf(stdout, "level %ld procnbr %ld subtreework %g \n",
                    (long)-ctrl->candtab[rootnum].treelevel, (long)candnbr,
                    costmtx->cblktab[rootnum].subtree);

          ctrl->candtab[rootnum].distrib = D1;
          for(i=0;i<ctrl->etree->nodetab[rootnum].sonsnbr;i++)
            setSubtreeDistribType(symbptr, costmtx, TSON(ctrl->etree, rootnum, i),  ctrl, D1);
        }
    }
}
/*
  Function: splitOnProcs

  Parameters:
    symbmtx    - Symbolic matrix
    extrasymb  -
    extracost  -
    ctrl       -
    dofptr     -
    cblknum    -
    procnbr    -
*/
void splitOnProcs(SymbolMatrix      *symbmtx,
                  ExtraSymbolMatrix *extrasymb,
                  ExtraCostMatrix   *extracost,
                  BlendCtrl         *ctrl,
                  const Dof         *dofptr,
                  PASTIX_INT                cblknum,
                  PASTIX_INT                procnbr)
{
    PASTIX_INT i;
    PASTIX_INT blas_min_col;
    PASTIX_INT blas_max_col;
    PASTIX_INT pas;
    PASTIX_INT *seq;
    PASTIX_INT nseq;


    /* if only one proc : no need to split */
    /*if(procnbr == 1)
      return;*/


    /* Compute minimun broadness for splitting this cblk */
    if(procnbr > ctrl->option->ratiolimit)
      {
        blas_min_col = ctrl->option->blblokmin;
        blas_max_col = ctrl->option->blblokmax;
      }
    else
      {
        blas_min_col = ctrl->option->blcolmin;
        blas_max_col = ctrl->option->blcolmax;
      }


#ifdef DOF_CONSTANT
    blas_min_col /= dofptr->noddval;
    if(blas_min_col == 0)
      blas_min_col = 1;
#endif

    /* number of split, for instance
       we choose to split at the maximum */
#ifdef DOF_CONSTANT

    if(procnbr == 1)
      {
        /*** Need to split big supernode because
             the diagonal block factorization is written
             in BLAS1 (due to the pivoting in LDLt and LU) ***/
        pas = symbmtx->cblktab[cblknum].lcolnum -
          symbmtx->cblktab[cblknum].fcolnum + 1;

        /*blas_max_col = 100;*/
        /*
         * if the column block size is small enough there is no need to
         * split it.
         */
        if(pas <= blas_max_col)
          return;
        else
          pas = blas_max_col;

        nseq = (symbmtx->cblktab[cblknum].lcolnum -
                symbmtx->cblktab[cblknum].fcolnum + 1)/pas;

      }
    else
      {
        PASTIX_INT abs = ctrl->option->abs;
        if(procnbr > ctrl->option->ratiolimit)
          {
            abs *= 2; /* Increase abs for 2D */
          }

        /***  If option adaptative block size is set then compute the size of a column block ***/
        if(abs > 0)
          {
            pas = (symbmtx->cblktab[cblknum].lcolnum -
                   symbmtx->cblktab[cblknum].fcolnum + 1)/(abs * procnbr);

            pas = MAX(pas, blas_min_col);
            pas = MIN(pas, blas_max_col);

            nseq = (symbmtx->cblktab[cblknum].lcolnum -
                    symbmtx->cblktab[cblknum].fcolnum + 1)/pas;
          }
        else
          nseq =(int)ceil((symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1.0)/ blas_min_col);
      }

    /*nseq = splitSeqCblk2D(cblknum, procnbr, symbmtx, extrasymb, dofptr, ctrl, P2D);*/
#endif

   /** No parallelism available above 4 splitted cblk **/
    if(nseq < 4)
      return;

#ifdef SMART_CBLK_SPLIT
    {
      int old_nseq = nseq;
      smart_cblk_split(ctrl,
                       symbmtx,
                       cblknum,
                       procnbr,
                       blas_min_col,
                       blas_max_col,
                       &nseq,
                       &seq);
      /** No parallelism available above 4 splitted cblk **/
      if(nseq < 4)
        return;
    }
#else /* SMART_CBLK_SPLIT */
    pas = (int)((symbmtx->cblktab[cblknum].lcolnum -
                 symbmtx->cblktab[cblknum].fcolnum+1.0)/nseq);
    /* extend intvec to contain seq */
    extendint_ToSize(2*nseq, ctrl->intvec);
    seq = ctrl->intvec->inttab;

    for(i=0;i<nseq-1;i++)
        {
            seq[2*i]   = symbmtx->cblktab[cblknum].fcolnum + pas*i;
            seq[2*i+1] = symbmtx->cblktab[cblknum].fcolnum + pas*(i+1)-1;
        }
    seq[2*(nseq-1)] =  symbmtx->cblktab[cblknum].fcolnum + pas*(nseq-1);
    seq[2*nseq-1]   =  symbmtx->cblktab[cblknum].lcolnum;

    ASSERTDBG(seq[2*(nseq-1)]<= seq[2*nseq-1],MOD_BLEND);

#endif /* SMART_CBLK_SPLIT */
    splitCblk(symbmtx, extrasymb, extracost, ctrl, dofptr, cblknum, nseq, seq);
}


/*
  Function: splitCblk

  Split a column bloc in nseq column blocs according to the seq array.

  Update symbolic matrix and extra symbolic matrix.

  Update cost of the node descendant of the splitted column bloc
  in the elimination tree.

  Parameters:
    symbmtx   - Symbol matrix.
    extrasymb - Extra symbol matrix.
    extracost -
    ctrl      - Blend control structure.
    dofptr    - Structure for degree of freedom,
                ! does not work anymore !.
    cblknum   - Column bloc to split.
    nseq      - Number of part of the split
    *seq      - Splitting indexes array.
*/
void  splitCblk(SymbolMatrix      *symbmtx,
                ExtraSymbolMatrix *extrasymb,
                ExtraCostMatrix   *extracost,
                BlendCtrl         *ctrl,
                const Dof         *dofptr,
                PASTIX_INT         cblknum,
                PASTIX_INT         nseq,
                PASTIX_INT        *seq)
{
    PASTIX_INT i, j, s;
    PASTIX_INT bloknbr      = 0;
    PASTIX_INT splitbloknbr = 0;
    PASTIX_INT saved_nseq   = 0;

    /**no need to split **/
    if(nseq == 1)
      return;

    /* ASSERT that the cblk can be splitted according to seq array */
    ASSERTDBG(symbmtx->cblktab[cblknum].fcolnum == seq[0],MOD_BLEND);
    ASSERTDBG(symbmtx->cblktab[cblknum].lcolnum == seq[2*nseq-1],MOD_BLEND);

    /* number of blok in the cbl to be split (diag included)  */
    bloknbr = 0;
    for(j=symbmtx->cblktab[cblknum].bloknum;j<symbmtx->cblktab[cblknum+1].bloknum;j++)
        if(extrasymb->sptblok[j]>=0)
          bloknbr += extrasymb->sptblnb[j];
        else
          bloknbr++;

    /*
     * XL: For Schur complement we keep the last column block.as full
     */
    if (ctrl->option->iparm[IPARM_SCHUR] == API_YES &&
         symbmtx->cblktab[cblknum].lcolnum == ctrl->option->n -1) {
      PASTIX_INT tmp;
      saved_nseq = nseq;
      tmp = seq[1];
      seq[1] = seq[(nseq-1)*2+1];
      seq[2*(nseq-1)+1] = tmp;
      nseq = 1;
    }

    /** mark the cblk to be splitted **/
    /** NB odb of this cblk don't need to be marked
     *  because they won't be splitted any more in our 
     *  top-down tree strategy **/
    extrasymb->sptcblk[cblknum] = extrasymb->curcblk;
    extrasymb->sptcbnb[cblknum] = nseq;

    /** now create the new cblk and associated blok **/

    for(i=0;i<nseq;i++)
      {
        extrasymb->cblktab[extrasymb->curcblk].fcolnum = seq[2*i];
        extrasymb->cblktab[extrasymb->curcblk].lcolnum = seq[2*i+1];
        extrasymb->cblktab[extrasymb->curcblk].bloknum = extrasymb->curblok;

        /* create diag blok */
        extrasymb->bloktab[extrasymb->curblok].frownum = seq[2*i];
        extrasymb->bloktab[extrasymb->curblok].lrownum = seq[2*i+1];
        extrasymb->bloktab[extrasymb->curblok].cblknum = extrasymb->curcblk;
        extra_inc_blok(extrasymb, extracost);
        splitbloknbr++;
        /* create odb due to the splitting of the diag blok */
        for(j=i+1;j<nseq;j++)
          {
            extrasymb->bloktab[extrasymb->curblok].frownum = seq[2*j];
            extrasymb->bloktab[extrasymb->curblok].lrownum = seq[2*j+1];
            extrasymb->bloktab[extrasymb->curblok].cblknum = extrasymb->sptcblk[cblknum] + j;
            extra_inc_blok(extrasymb, extracost);
            splitbloknbr++;
          }
        /* create other odb */
        /* We have to test if some of them have been splitted before */
        for(j=symbmtx->cblktab[cblknum].bloknum+1;j<symbmtx->cblktab[cblknum+1].bloknum;j++)
          {
            /* this odb hasn't been splitted */
            if(extrasymb->sptblok[j]<0)
              {
                extrasymb->bloktab[extrasymb->curblok].frownum = symbmtx->bloktab[j].frownum;
                extrasymb->bloktab[extrasymb->curblok].lrownum = symbmtx->bloktab[j].lrownum;
                extrasymb->bloktab[extrasymb->curblok].cblknum = symbmtx->bloktab[j].cblknum;
                extra_inc_blok(extrasymb, extracost);
                splitbloknbr++;
              }
            /* this odb has been splitted before (its facing diag blok has been splitted) */
            else
              {
                for(s=extrasymb->sptblok[j];s<extrasymb->sptblok[j]+extrasymb->sptblnb[j];s++)
                  {
                    extrasymb->bloktab[extrasymb->curblok].frownum = extrasymb->bloktab[s].frownum;
                    extrasymb->bloktab[extrasymb->curblok].lrownum = extrasymb->bloktab[s].lrownum;
                    extrasymb->bloktab[extrasymb->curblok].cblknum = extrasymb->bloktab[s].cblknum;
                    extra_inc_blok(extrasymb, extracost);
                    splitbloknbr++;
                  }
              }
          }
        extra_inc_cblk(extrasymb, extracost);
      }
    /* update extracblk and extrablok */
    extrasymb->addcblk += nseq-1;
    extrasymb->addblok += splitbloknbr - bloknbr;

    if (ctrl->option->iparm[IPARM_SCHUR] == API_YES &&
        symbmtx->cblktab[cblknum].lcolnum == ctrl->option->n -1) {
      PASTIX_INT tmp;
      nseq = saved_nseq;
      tmp = seq[1];
      seq[1] = seq[(nseq-1)*2+1];
      seq[2*(nseq-1)+1] = tmp;
    }

    /** we have to add an extra cblk to extrasymb because of side effect **/
    /* NB extrasymb->cblktab have at least one extra-allocated cells
     *    so don't worry about allocated memory size */
    ASSERTDBG(extrasymb->sizcblk > extrasymb->curcblk,MOD_BLEND);
    extrasymb->cblktab[extrasymb->curcblk].bloknum = extrasymb->curblok;


    /** Now we're going to split odb that hit our diag blok **/
    /** We have to mark them as splitted bloks because they may be split again later
      if their cblk owner is splitted **/
    /** NB odb blok that hit the diag (those that we're about to split )
      can have been splitted before because of our top-down tree strategy **/
    for(i=0;i<ctrl->egraph->verttab[cblknum].innbr;i++)
      {
        PASTIX_INT sptbloknbr; /* number of splitted bloks resulting */
        PASTIX_INT bloknum;
        bloknum = ctrl->egraph->inbltab[ctrl->egraph->verttab[cblknum].innum+i];

        ASSERTDBG(symbmtx->bloktab[bloknum].cblknum == cblknum,MOD_BLEND);

        sptbloknbr = 0;
        /* mark this blok as splitted */
        extrasymb->sptblok[bloknum] = extrasymb->curblok;
        for(j=0;j<nseq;j++)
          {
            /* there are six possible cases of intersection
             beetween the blok and the seq we consider ,
             among these six cases, only 4 of them are
             not empty intersection */

            /* empty intersections */
            if(symbmtx->bloktab[bloknum].frownum > seq[2*j+1])
              continue;

            if(symbmtx->bloktab[bloknum].lrownum < seq[2*j])
              /* in this case there will no more splitted blok to create */
              break;

            /* not empty intersections */
            if((symbmtx->bloktab[bloknum].frownum >= seq[2*j])
               && (symbmtx->bloktab[bloknum].lrownum >= seq[2*j+1]))
              {
                extrasymb->bloktab[extrasymb->curblok].frownum = symbmtx->bloktab[bloknum].frownum;
                extrasymb->bloktab[extrasymb->curblok].lrownum = seq[2*j+1];
                goto endloop;
              }
            if((symbmtx->bloktab[bloknum].frownum <= seq[2*j])
               && (symbmtx->bloktab[bloknum].lrownum >= seq[2*j+1]))
              {
                extrasymb->bloktab[extrasymb->curblok].frownum = seq[2*j];
                extrasymb->bloktab[extrasymb->curblok].lrownum = seq[2*j+1];
                goto endloop;
              }
            if((symbmtx->bloktab[bloknum].frownum <= seq[2*j])
               && (symbmtx->bloktab[bloknum].lrownum <= seq[2*j+1]))
              {
                extrasymb->bloktab[extrasymb->curblok].frownum = seq[2*j];
                extrasymb->bloktab[extrasymb->curblok].lrownum = symbmtx->bloktab[bloknum].lrownum;
                goto endloop;
              }
            if((symbmtx->bloktab[bloknum].frownum >= seq[2*j])
               && (symbmtx->bloktab[bloknum].lrownum <= seq[2*j+1]))
              {
                extrasymb->bloktab[extrasymb->curblok].frownum = symbmtx->bloktab[bloknum].frownum;
                extrasymb->bloktab[extrasymb->curblok].lrownum = symbmtx->bloktab[bloknum].lrownum;
                goto endloop;
              }
            endloop:
            if (ctrl->option->iparm[IPARM_SCHUR] == API_YES &&
                symbmtx->cblktab[cblknum].lcolnum == ctrl->option->n -1)
              extrasymb->bloktab[extrasymb->curblok].cblknum = extrasymb->sptcblk[cblknum];
            else
              extrasymb->bloktab[extrasymb->curblok].cblknum = extrasymb->sptcblk[cblknum]+j;
            sptbloknbr++;
            extra_inc_blok(extrasymb, extracost);
          }
        extrasymb->sptblnb[bloknum] = sptbloknbr;
        extrasymb->addblok += sptbloknbr-1;

        /** update cost of the cblk owning the splitted bloks **/
        blokUpdateCost(bloknum, ctrl->egraph->ownetab[bloknum], ctrl->costmtx, extracost, symbmtx, extrasymb, ctrl, dofptr);
      }
}

/*
  Function: propMappTree

  Repartition symbolic matrix using proportionnal mapping method

  Parameters:
    symbmtx   - Symbolic matrix
    extrasymb -
    extracost -
    ctrl      -
    dofptr    -
*/
void propMappTree(SymbolMatrix      *symbmtx,
                  ExtraSymbolMatrix *extrasymb,
                  ExtraCostMatrix   *extracost,
                  BlendCtrl         *ctrl,
                  const Dof         *dofptr)
{
  double isocost;
  PASTIX_INT pr;
  double *cost_remain = NULL;

  MALLOC_INTERN(cost_remain, ctrl->procnbr, double);

  isocost = ctrl->costmtx->cblktab[ROOT(ctrl->etree)].subtree / ctrl->procnbr;

  /* compute cost to get in the elimination tree for each processor */
  for(pr=0;pr<ctrl->procnbr;pr++)
    cost_remain[pr] = isocost;

  if (!ctrl->option->nocrossproc)
    propMappSubtree(symbmtx, extrasymb, extracost, ctrl, dofptr,
                    ROOT(ctrl->etree), 0, ctrl->procnbr-1,
                    NOCLUSTER, cost_remain);
  else
    propMappSubtreeNC(symbmtx, extrasymb, extracost, ctrl, dofptr,
                      ROOT(ctrl->etree), 0, ctrl->procnbr-1,
                      NOCLUSTER, cost_remain);
}

/*
  Function: propMappSubtree


  Parameters:
    symbmtx     - Symbolic matrix.
    extrasymb   -
    extracost   -
    ctrl        -
    dofptr      -
    rootnum     -
    fcandnum    -
    lcandnum    -
    cluster     -
    cost_remain -


*/
void propMappSubtree(SymbolMatrix      *symbmtx,
                     ExtraSymbolMatrix *extrasymb,
                     ExtraCostMatrix   *extracost,
                     BlendCtrl         *ctrl,
                     const Dof         *dofptr,
                     PASTIX_INT                rootnum,
                     PASTIX_INT                fcandnum,
                     PASTIX_INT                lcandnum,
                     PASTIX_INT                cluster,
                     double            *cost_remain)
{
  PASTIX_INT p;
  PASTIX_INT candnbr;
  PASTIX_INT fcand = 0;
  PASTIX_INT lcand = 0;
  PASTIX_INT i;
  PASTIX_INT sonsnbr;
  int nbthrdbyproc;
  double isocost;
  double aspt_cost;
  double cumul_cost;
  double *sub_cost_remain = NULL;
  double epsilon;

  /* initialize nbthrdbyproc */
  if (ctrl->cudanbr > 0)
    nbthrdbyproc = MAX(ctrl->thrdnbr / (ctrl->procnbr*2*ctrl->cudanbr), 1);
  else
    nbthrdbyproc = MAX(ctrl->thrdnbr / (ctrl->procnbr), 1);
  candnbr      = lcandnum - fcandnum + 1;

  /* Si il n'y a qu'un candidat, tous le sous-arbre lui appartient */
  if (candnbr == 1)
    {
      memFree_null(cost_remain);
      propMappSubtreeOn1P(symbmtx, extrasymb, extracost, ctrl, dofptr,
                          rootnum, fcandnum, lcandnum, cluster);
      return;
    }

  /* Set the cand group for this tree node */
  if(ctrl->option->allcand)
    {
      ctrl->candtab[rootnum].fcandnum = 0;
      ctrl->candtab[rootnum].lcandnum = ctrl->procnbr - 1;
      ctrl->candtab[rootnum].cluster  = cluster;
    }
  else
    {
      ctrl->candtab[rootnum].fcandnum = fcandnum;
      ctrl->candtab[rootnum].lcandnum = lcandnum;
      ctrl->candtab[rootnum].cluster  = cluster;

    }

  /* this treenode is a leave split it and return */
  if(ctrl->etree->nodetab[rootnum].sonsnbr == 0)
    {
      splitOnProcs(symbmtx, extrasymb, extracost, ctrl,
                   dofptr, rootnum, candnbr*nbthrdbyproc);
      memFree_null(cost_remain);
      return;
    }

  /* work that each processor is intended to get from this treenode */
  isocost = ctrl->costmtx->cblktab[rootnum].total / candnbr;
  for(p=0;p<candnbr;p++)
    cost_remain[p] -= isocost;

  /* split the treenode */
  splitOnProcs(symbmtx, extrasymb, extracost, ctrl, dofptr,
               rootnum, candnbr*nbthrdbyproc);

  /* Correct the subtree cost */
  subtreeUpdateCost(rootnum, ctrl->costmtx, ctrl->etree);

  /* Cost remaining in the descendance of the treenode after split*/
  aspt_cost = ctrl->costmtx->cblktab[rootnum].subtree - ctrl->costmtx->cblktab[rootnum].total;

  /* if the proc cand has reached its cost to get
     forget it (only the first and last proc in this group
     can have reached their cost to get */
  if(cost_remain[0] <= 0)
    fcand = 1;
  else
    fcand = 0;
  if (cost_remain[candnbr-1] <= 0)
    candnbr--;

  /** ASSERT that the sum of cost_remain in used proc equals to
    after split cost **/
  cumul_cost = 0;
  for(i=fcand;i<candnbr;i++)
    cumul_cost += cost_remain[i];
  for(i=fcand;i<candnbr;i++)
    cost_remain[i] *= aspt_cost/cumul_cost;

  /** Compute the minimun participation rate of a candidat processor**/
  epsilon = CROSS_TOLERANCE*cumul_cost/candnbr ;

#ifdef DEBUG_BLEND_SPLIT
  {
    fprintf(stdout, "Rootnum %ld fproc %ld lproc %ld [ ", (long)rootnum, (long)fcandnum, (long)lcandnum);
    for(p=0;p<candnbr;p++)
      fprintf(stdout, "%g ", cost_remain[p]);
    fprintf(stdout, " ]\n");
    fprintf(stdout, " Sons ");
    for(i=0;i<ctrl->etree->nodetab[rootnum].sonsnbr;i++)
      fprintf(stdout, " [%ld, %g] ", (long)i, ctrl->costmtx->cblktab[TSON(ctrl->etree, rootnum, i)].subtree);
    fprintf(stdout, "\n");
  }
#endif

  lcand = fcand;

  /* compute the cand group for each proc */
  sonsnbr = ctrl->etree->nodetab[rootnum].sonsnbr;

  for(i=0;i<sonsnbr;i++)
    {

      /** Cost in the current subtree to be mapped **/
      cumul_cost = ctrl->costmtx->cblktab[TSON(ctrl->etree, rootnum, i)].subtree;


      /** Are we to take the last processor used in the previous sons ?
          (== crossing processor)**/
      if(!ctrl->option->nocrossproc)
        {
          /* if crossing proc is fed , forget it */
          if( (cost_remain[lcand] <= epsilon) && (lcand < candnbr-1) )
            {
              fcand = lcand + 1;
            }
          else
            fcand = lcand;
        }
      else
        {
          /* if last proc can compute all the subtree */
          if (cost_remain[fcand] > (cumul_cost - epsilon))
            {
              fcand = lcand ;
              MALLOC_INTERN(sub_cost_remain, lcand-fcand+1, double);
              sub_cost_remain[0] = cost_remain[fcand];
              propMappSubtree(symbmtx, extrasymb, extracost, ctrl, dofptr, TSON(ctrl->etree, rootnum, i),
                              fcandnum+fcand, fcandnum+lcand, cluster, sub_cost_remain);
              cost_remain[fcand] -= cumul_cost;
              continue;
            }
          else if(lcand < candnbr - 1)
            fcand = lcand+1;
        }

      /* Update subtree cost and remain cost of fcand proc */
      cumul_cost -= cost_remain[fcand];

      if(cost_remain[fcand]<=0)
        print_debug(DBG_BUBBLESPLIT, "candnbr %ld fcand %ld cost_remain[fcand] %g epsilon %g \n",
                    (long)candnbr, (long)fcand, cost_remain[fcand], epsilon);

      lcand = fcand;
      while ((cumul_cost > epsilon) && (lcand < candnbr - 1))
        {
          if (cost_remain[lcand+1]<=0)
            break;
          lcand++;
          cumul_cost -= cost_remain[lcand];
          ASSERTDBG(cost_remain[lcand]>0,MOD_BLEND);
        }

      MALLOC_INTERN(sub_cost_remain, lcand-fcand+1, double);

      /* compute the cost to get for each proc cand for this subtree */
      for(p=0;p<lcand-fcand;p++)
        sub_cost_remain[p] = cost_remain[fcand+p];
      sub_cost_remain[lcand-fcand] = cost_remain[lcand] + cumul_cost; /* cumul_cost <= 0 */

      if(!ctrl->option->nocrossproc)
        {
          /* retrieve cost for the crossing proc */
          cost_remain[lcand] = -cumul_cost;
        }
      else
        {
          if(lcand > fcand)
            if( (-cumul_cost)> (cost_remain[lcand]/2.0))
              lcand--;
        }

      /* go on to subtree */
      propMappSubtree(symbmtx, extrasymb, extracost, ctrl, dofptr, TSON(ctrl->etree, rootnum, i),
                      fcandnum+fcand, fcandnum+lcand, cluster, sub_cost_remain);

    }

  memFree_null(cost_remain);
  return;
}

/*
  Function: propMappSubtreeNC

  Parameters:
    symbmtx     - Symbolic matrix.
    extrasymb   -
    extracost   -
    ctrl        -
    dofptr      -
    rootnum     -
    fcandnum    -
    lcandnum    -
    cluster     -
    cost_remain -

*/
void propMappSubtreeNC(SymbolMatrix      *symbmtx,
                       ExtraSymbolMatrix *extrasymb,
                       ExtraCostMatrix   *extracost,
                       BlendCtrl         *ctrl,
                       const Dof         *dofptr,
                       PASTIX_INT                rootnum,
                       PASTIX_INT                fcandnum,
                       PASTIX_INT                lcandnum,
                       PASTIX_INT                cluster,
                       double            *cost_remain)
{
  PASTIX_INT p;
  PASTIX_INT candnbr;
  PASTIX_INT fcand = 0;
  PASTIX_INT lcand = 0;
  PASTIX_INT i;
  PASTIX_INT sonsnbr;
  double isocost;
  double bspt_cost;
  double aspt_cost;
  double cumul_cost;
  double *sub_cost_remain = NULL;
  double epsilon;
  double soncost;
  Queue *queue_tree;
  Queue *queue_proc;

  candnbr = lcandnum - fcandnum + 1;

  /* Si il n'y a qu'un candidat, tous le sous-arbre lui appartient */
  if(candnbr == 1)
    {
      memFree(cost_remain);
      propMappSubtreeOn1P(symbmtx, extrasymb, extracost, ctrl, dofptr,
                          rootnum, fcandnum, lcandnum, cluster);
      return;
    }

  /* Set the cand group for this tree node */
  if(ctrl->option->allcand)
    {
      ctrl->candtab[rootnum].fcandnum = 0;
      /*ctrl->candtab[rootnum].lcandnum = ctrl->clustnbr - 1;*/
      ctrl->candtab[rootnum].lcandnum = ctrl->procnbr - 1;
      ctrl->candtab[rootnum].cluster  = cluster;
    }
  else
    {
      ctrl->candtab[rootnum].fcandnum = fcandnum;
      ctrl->candtab[rootnum].lcandnum = lcandnum;
      ctrl->candtab[rootnum].cluster  = cluster;
    }

  /* this treenode is a leave split it and return */
  if(ctrl->etree->nodetab[rootnum].sonsnbr == 0)
    {
      splitOnProcs(symbmtx, extrasymb, extracost, ctrl, dofptr, rootnum, lcandnum-fcandnum+1);
      memFree(cost_remain);
      return;
    }

  /* work that each processor is intended to get from this treenode */
  isocost = ctrl->costmtx->cblktab[rootnum].total / candnbr;
  for(p=0;p<candnbr;p++)
    cost_remain[p] -= isocost;

  /* Cost remaining in the descendance of the treenode before split*/
  bspt_cost = ctrl->costmtx->cblktab[rootnum].subtree - ctrl->costmtx->cblktab[rootnum].total;

  /* split the treenode */
  splitOnProcs(symbmtx, extrasymb, extracost, ctrl, dofptr, rootnum, lcandnum-fcandnum+1);

  /* Correct the subtree cost */
  subtreeUpdateCost(rootnum, ctrl->costmtx, ctrl->etree);

  /* Cost remaining in the descendance of the treenode after split*/
  aspt_cost = ctrl->costmtx->cblktab[rootnum].subtree - ctrl->costmtx->cblktab[rootnum].total;

  if(cost_remain[0] <= 0)
    fcand = 1;
  else
    fcand = 0;
  if(cost_remain[candnbr-1] <= 0)
    candnbr--;


  /** ASSERT that the sum of cost_remain in used proc equals to
    after split cost **/
  cumul_cost = 0;
  for(i=fcand;i<candnbr;i++)
    cumul_cost += cost_remain[i];
  for(i=fcand;i<candnbr;i++)
    cost_remain[i] *= aspt_cost/cumul_cost;

  /** Compute the minimun participation rate of a candidat processor**/
  epsilon = CROSS_TOLERANCE*cumul_cost/candnbr ;


#ifdef DEBUG_BLEND_SPLIT
  {
    fprintf(stdout, "Rootnum %ld fproc %ld lproc %ld [ ", (long)rootnum, (long)fcandnum, (long)lcandnum);
    for(p=0;p<candnbr;p++)
      fprintf(stdout, "%g ", cost_remain[p]);
    fprintf(stdout, " ]\n");
    fprintf(stdout, " Sons ");
    for(i=0;i<ctrl->etree->nodetab[rootnum].sonsnbr;i++)
      fprintf(stdout, " [%ld, %g] ", i, ctrl->costmtx->cblktab[TSON(ctrl->etree, rootnum, i)].subtree);
    fprintf(stdout, "\n");
  }
#endif

  lcand = fcand;

  /* compute the cand group for each proc */
  bspt_cost = 0.0;
  sonsnbr = ctrl->etree->nodetab[rootnum].sonsnbr;

  /* allocate queue of proc and subtree */
  MALLOC_INTERN(queue_tree, 1, Queue);
  MALLOC_INTERN(queue_proc, 1, Queue);
  queueInit(queue_tree, sonsnbr);
  queueInit(queue_proc, candnbr);

  /* Fill queue subtree order by cost descending */
  for(i=0;i<sonsnbr;i++)
    {

      /** Cost in the current subtree to be mapped **/
      cumul_cost = -ctrl->costmtx->cblktab[TSON(ctrl->etree, rootnum, i)].subtree;

      /* Cost of the root node in the subtree */
      soncost    = -ctrl->costmtx->cblktab[TSON(ctrl->etree, rootnum, i)].total;

      queueAdd2(queue_tree, i, cumul_cost, (PASTIX_INT)soncost);
      bspt_cost += cumul_cost;
    }

  /* Proportionnal mapping of the subtree on cand */
  /* Premiere repartition des sous-arbres de cout superieur a un proc */
  while (queueSize(queue_tree) > 0)
    {
      i = queueGet(queue_tree);

      /** Cost in the current subtree to be mapped **/
      cumul_cost = ctrl->costmtx->cblktab[TSON(ctrl->etree, rootnum, i)].subtree;

      if (cumul_cost < cost_remain[fcand]){
        /* No more split subtree */
        cost_remain[fcand] -= cumul_cost;
        propMappSubtreeOn1P(symbmtx, extrasymb, extracost, ctrl, dofptr,
                            TSON(ctrl->etree, rootnum, i), fcandnum+fcand, fcandnum+fcand, cluster);
        break;
      }

      /* Update subtree cost and remain cost of fcand proc */
      cumul_cost -= cost_remain[fcand];

#ifdef DEBUG_BLEND
      if(cost_remain[fcand]<=0)
        fprintf(stderr, "candnbr %ld fcand %ld cost_remain[fcand] %g epsilon %g \n",
                (long)candnbr, (long)fcand, cost_remain[fcand], epsilon);
      /*ASSERT(cost_remain[fcand]>0,MOD_BLEND);*/
#endif

      lcand = fcand;
      while ((cumul_cost > ((lcand - fcand +1)*epsilon)) &&
             (lcand < candnbr - 1))
        {
          lcand++;
          cumul_cost -= cost_remain[lcand];
          ASSERTDBG(cost_remain[lcand]>0,MOD_BLEND);
        }

      MALLOC_INTERN(sub_cost_remain, lcand-fcand+1, double);
      if (cumul_cost > 0)
        {
          for(p=0;p<lcand-fcand+1;p++)
            {
              sub_cost_remain[p]   = cost_remain[fcand+p] + (cumul_cost / (lcand - fcand + 1));
              cost_remain[fcand+p] = - cumul_cost / (lcand - fcand + 1);
            }
        }
      else
        {
          /* compute the cost to get for each proc cand for this subtree */
          for(p=0;p<lcand-fcand;p++)
            {
              sub_cost_remain[p] = cost_remain[fcand+p];
              cost_remain[fcand+p] = 0.0;
            }
          sub_cost_remain[lcand-fcand] = cost_remain[lcand] + cumul_cost; /* cumul_cost <= 0 */
          cost_remain[fcand+p] = - cumul_cost ;
        }

      /* go on to subtree */
      propMappSubtreeNC(symbmtx, extrasymb, extracost, ctrl, dofptr, TSON(ctrl->etree, rootnum, i),
                        fcandnum+fcand, fcandnum+lcand, cluster, sub_cost_remain);

      if (lcand < candnbr - 1)
        fcand = lcand+1;
      else
        break;
    }


  /* Fill queue proc order by remain cost descending */
  for (i=0; i<candnbr; i++)
    queueAdd(queue_proc, i, -cost_remain[i]);

  while (queueSize(queue_tree) > 0)
    {
      i = queueGet(queue_tree);
      fcand = queueGet(queue_proc);

      /** Cost in the current subtree to be mapped **/
      cumul_cost = ctrl->costmtx->cblktab[TSON(ctrl->etree, rootnum, i)].subtree;

      propMappSubtreeOn1P(symbmtx, extrasymb, extracost, ctrl, dofptr,
                          TSON(ctrl->etree, rootnum, i), fcandnum+fcand, fcandnum+fcand, cluster);

      cost_remain[fcand] -= cumul_cost;
      queueAdd(queue_proc, fcand, -cost_remain[fcand]);
    }

  queueExit(queue_tree);
  queueExit(queue_proc);
  memFree(queue_tree);
  memFree(queue_proc);
  memFree(cost_remain);
  return;
}


void propMappSubtreeOn1P(SymbolMatrix *symbmtx, ExtraSymbolMatrix *extrasymb, ExtraCostMatrix *extracost,
                         BlendCtrl *ctrl, const Dof * dofptr,
                         PASTIX_INT rootnum, PASTIX_INT fcandnum, PASTIX_INT lcandnum, PASTIX_INT cluster)
{
  PASTIX_INT i;
  PASTIX_INT sonsnbr;

  ctrl->candtab[rootnum].fcandnum = fcandnum;
  ctrl->candtab[rootnum].lcandnum = lcandnum;
  ctrl->candtab[rootnum].cluster  = cluster;

  /* split the treenode */
  splitOnProcs(symbmtx, extrasymb, extracost, ctrl, dofptr, rootnum, 1);

  /* Correct the subtree cost */
  subtreeUpdateCost(rootnum, ctrl->costmtx, ctrl->etree);

  sonsnbr = ctrl->etree->nodetab[rootnum].sonsnbr;

  /* Fill queue subtree order by cost descending */
  for(i=0;i<sonsnbr;i++)
    propMappSubtreeOn1P(symbmtx, extrasymb, extracost, ctrl, dofptr,
                        TSON(ctrl->etree, rootnum, i), fcandnum, lcandnum, cluster);

  return;
}

void subtreeSetCand(PASTIX_INT procnum, PASTIX_INT rootnum, BlendCtrl *ctrl, double rcost)
{
    PASTIX_INT i;

    ctrl->candtab[rootnum].fcandnum = procnum;
    ctrl->candtab[rootnum].lcandnum = procnum;
    ctrl->candtab[rootnum].cluster  = NOCLUSTER;
    for(i=0;i<ctrl->etree->nodetab[rootnum].sonsnbr;i++)
        subtreeSetCand(procnum, TSON(ctrl->etree, rootnum, i), ctrl, rcost);
}

double maxProcCost(double *proc_cost, PASTIX_INT procnbr)
{
    double maxcost = 0;
    PASTIX_INT p;
    for(p=0;p<procnbr;p++)
        if(proc_cost[p]>maxcost)
            maxcost = proc_cost[p];
    return maxcost;
}

void propMappTreeNoSplit(SymbolMatrix *symbmtx, BlendCtrl *ctrl, const Dof *dofptr)
{
  double *cost_remain = NULL;
  double isocost;
  PASTIX_INT p;
  MALLOC_INTERN(cost_remain, ctrl->clustnbr, double);
  isocost = ctrl->costmtx->cblktab[ROOT(ctrl->etree)].subtree / ctrl->clustnbr;

  /* compute cost to get in the elimination tree for each processor */
  for(p=0;p<ctrl->clustnbr;p++)
    cost_remain[p] = isocost;

  propMappSubtreeNoSplit(symbmtx, ctrl, dofptr, ROOT(ctrl->etree), 0, ctrl->clustnbr-1, cost_remain);

}

void propMappSubtreeNoSplit(SymbolMatrix *symbmtx, BlendCtrl *ctrl, const Dof * dofptr,
                            PASTIX_INT rootnum, PASTIX_INT fcandnum, PASTIX_INT lcandnum, double *cost_remain)
{
  PASTIX_INT p;
  PASTIX_INT procnbr;
  PASTIX_INT fcand = 0;
  PASTIX_INT lcand = 0;
  PASTIX_INT i;
  double isocost;
  double aspt_cost;
  double cumul_cost;
  double *sub_cost_remain = NULL;
  double epsilon;

  procnbr = lcandnum - fcandnum + 1;

  if(procnbr == 1)
    {
      /* only one proc Cand, all the subtree is devoted
         to this processeur */
      /*fprintf(stdout, "Mapp Root %ld on  proc %ld \n", (long) rootnum, fcandnum);*/
      memFree(cost_remain);
      subtreeSetCand(fcandnum, rootnum, ctrl, ctrl->costmtx->cblktab[rootnum].subtree);
      return;
    }

  /* Set the cand group for this tree node */
  if(ctrl->option->allcand)
    {
      ctrl->candtab[rootnum].fcandnum = 0;
      ctrl->candtab[rootnum].lcandnum = ctrl->clustnbr - 1;
    }
  else
    {
      /*if(ctrl->candtab[rootnum].fcandnum != fcandnum
         || ctrl->candtab[rootnum].lcandnum != lcandnum)
        fprintf(stdout, "[%ld %ld] --> [%ld %ld]\n", ctrl->candtab[rootnum].fcandnum, ctrl->candtab[rootnum].lcandnum, fcandnum, lcandnum);*/
      ctrl->candtab[rootnum].fcandnum = fcandnum;
      ctrl->candtab[rootnum].lcandnum = lcandnum;
    }
  /* this treenode is a leave split it and return */
  if(ctrl->etree->nodetab[rootnum].sonsnbr == 0)
    {
      /*splitOnProcs(symbmtx, extrasymb, extracost, ctrl, dofptr, rootnum, lcandnum-fcandnum+1);*/
      memFree(cost_remain);
      return;
    }

  /* work that each processor is intended to get from this treenode */
  /** Considere splitted cblk as the original one **/
  isocost = 0;
  isocost = ctrl->costmtx->cblktab[rootnum].total;
  while(ctrl->etree->nodetab[rootnum].sonsnbr == 1)
    {
      rootnum = TSON(ctrl->etree, rootnum, 0);
      isocost += ctrl->costmtx->cblktab[rootnum].total;
      if(ctrl->option->allcand)
        {
          ctrl->candtab[rootnum].fcandnum = 0;
          ctrl->candtab[rootnum].lcandnum = ctrl->clustnbr - 1;
        }
      else
        {
          ctrl->candtab[rootnum].fcandnum = fcandnum;
          ctrl->candtab[rootnum].lcandnum = lcandnum;
        }

    }
  isocost /=procnbr;



  for(p=0;p<procnbr;p++)
    cost_remain[p] -= isocost;


  /* split the treenode */
  /*splitOnProcs(symbmtx, extrasymb, extracost, ctrl, dofptr, rootnum, lcandnum-fcandnum+1);*/


  /* Correct the subtree cost */
  /*subtreeUpdateCost(rootnum, ctrl->costmtx, ctrl->etree);*/

  /* Cost remaining in the descendance of the treenode after split*/
  aspt_cost = ctrl->costmtx->cblktab[rootnum].subtree - ctrl->costmtx->cblktab[rootnum].total;


  /* if the proc cand has reached its cost to get
     forget it (only the first and last proc in this group
     can have reached their cost to get */
  if(cost_remain[0] <= 0)
    fcand = 1;
  else
    fcand = 0;
  if(cost_remain[procnbr-1] <= 0)
    procnbr--;

  /** ASSERT that the sum of cost_remain in used proc equals to
    after split cost **/
  /* OIMBE on doit avoir cumul_cost == aspt_cost ?? car on n'a pas redecoupe */
  cumul_cost = 0;
  for(i=fcand;i<procnbr;i++)
    cumul_cost += cost_remain[i];

  /*fprintf(stdout, "Cumul_cost %g aspt_cost %g \n", cumul_cost, aspt_cost);*/
  for(i=fcand;i<procnbr;i++)
    cost_remain[i] *= aspt_cost/cumul_cost;



  /** Compute the minimun participation rate of a candidat processor**/
  epsilon = CROSS_TOLERANCE*cumul_cost/procnbr ;

  /* For debug */
  /*{
    fprintf(stdout, "Rootnum %ld fproc %ld lproc %ld [ ", (long)rootnum, (long)fcandnum, (long)lcandnum);
    for(p=0;p<procnbr;p++)
      fprintf(stdout, "%g ", cost_remain[p]);
    fprintf(stdout, " ]\n");
    fprintf(stdout, " Sons ");
    for(i=0;i<ctrl->etree->nodetab[rootnum].sonsnbr;i++)
      fprintf(stdout, " [%ld, %g] ", i, ctrl->costmtx->cblktab[TSON(ctrl->etree, rootnum, i)].subtree);
    fprintf(stdout, "\n");
    fprintf(stdout, "Espilon %g \n", epsilon);
  }*/

  lcand = fcand;
  /* compute the cand group for each proc */
  for(i=0;i<ctrl->etree->nodetab[rootnum].sonsnbr;i++)
    {
      /** Cost in the current subtree to be mapped **/
      cumul_cost = ctrl->costmtx->cblktab[TSON(ctrl->etree, rootnum, i)].subtree;


      /** Are we to take the last processor used in the previous sons ?
        (== crossing processor)**/
      /*if(!ctrl->option->nocrossproc)
        {*/
      /* if crossing proc is fed , forget it */
      if( (cost_remain[lcand] <= epsilon) && (lcand < procnbr-1) )
        {
          /*fprintf(stdout, "yoyo \n");*/
          fcand = lcand + 1;
        }
      else
        fcand = lcand;

      /*else
        if(lcand < procnbr)
              fcand = lcand+1;*/

      cumul_cost -= cost_remain[fcand];
#ifdef DEBUG_BLEND
      if(cost_remain[fcand]<=0)
        fprintf(stderr, "procnbr %ld fcand %ld cost_remain[fcand] %g epsilon %g \n",
                (long)procnbr, (long)fcand, cost_remain[fcand], epsilon);
      /*ASSERT(cost_remain[fcand]>0,MOD_BLEND);*/
#endif
      lcand = fcand;

      while((cumul_cost > epsilon) && (lcand <procnbr-1))
        {
          if(cost_remain[lcand+1]<=0)
            break;
          lcand++;
          cumul_cost -= cost_remain[lcand];
          ASSERTDBG(cost_remain[lcand]>0,MOD_BLEND);
        }
      MALLOC_INTERN(sub_cost_remain, lcand-fcand+1, double);

      /* compute the cost to get for each proc cand for this subtree */
      for(p=0;p<lcand-fcand;p++)
        sub_cost_remain[p] = cost_remain[fcand+p];
      sub_cost_remain[lcand-fcand] = cost_remain[lcand] + cumul_cost; /* cumul_cost <= 0 */

      if(!ctrl->option->nocrossproc)
        {
          /* retrieve cost for the crossing proc */
          cost_remain[lcand] = -cumul_cost;

        }
      else
        {
          if(lcand > fcand)
            if( (-cumul_cost)> (cost_remain[lcand]/2.0))
              lcand--;
        }

      /*fprintf(stdout, "Cand proc [%ld %ld] for son %ld\n", fcandnum+fcand, fcandnum+lcand, i);*/
      /* go on to subtree */
      propMappSubtreeNoSplit(symbmtx, ctrl, dofptr, TSON(ctrl->etree, rootnum, i),
                             fcandnum+fcand, fcandnum+lcand, sub_cost_remain);



    }
  memFree(cost_remain);
  return;
}

/*+ Recompute cost of cblk which some odb have been splitted, return new total cost - old total cost +*/
double blokUpdateCost(PASTIX_INT bloknum, PASTIX_INT cblknum, CostMatrix *costmtx, ExtraCostMatrix *extracost, const SymbolMatrix *symbmtx, const ExtraSymbolMatrix *extrasymb, BlendCtrl *ctrl, const Dof * dofptr)
{
    PASTIX_INT L, h, g;
    PASTIX_INT s;
    double oldcost, newcost;
#ifndef DOF_CONSTANT
    PASTIX_INT i;
#endif

    ASSERTDBG(extrasymb->sptblok[bloknum] >= 0,MOD_BLEND);

    /** we need the height of the odb and the broadness
      of the cbl to compute the local compute cost **/

#ifdef DOF_CONSTANT
    L = (symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1)*(dofptr)->noddval;
#else
    for(i=symbmtx->cblktab[cblknum].fcolnum;i<=symbmtx->cblktab[cblknum].lcolnum;i++)
      L+= noddDlt(dofptr, i);
#endif
    /** no need to recompute the local compute cost because odb lines number
      is not changed **/


    /** recompute for each splitted odb its contribution compute cost and add cost **/
    oldcost = 0;
    newcost = 0;


    oldcost = costmtx->bloktab[bloknum].contrib; /* cost of blok not splitted */
    g = costmtx->bloktab[bloknum].linenbr;

    /* now compute the new cost of each resulting blok of the spliting */
    for(s = extrasymb->sptblok[bloknum];s<extrasymb->sptblok[bloknum]+extrasymb->sptblnb[bloknum];s++)
        {

#ifdef  DOF_CONSTANT
        h = (extrasymb->bloktab[s].lrownum - extrasymb->bloktab[s].frownum + 1)*(dofptr)->noddval;
#else
        for(i=extrasymb->bloktab[s].frownum;i<=extrasymb->bloktab[s].lrownum;i++)
          h+= noddDlt(dofptr, i);
#endif
            extracost->bloktab[s].linenbr = g;
            extracost->bloktab[s].contrib =  contribCompCost(L, h, g) + contribAddCost(h, g);
            newcost += extracost->bloktab[s].contrib;
            g -= h;
        }


    costmtx->cblktab[cblknum].send  += newcost - oldcost;

    if(ctrl->candtab[cblknum].distrib == D1)
      costmtx->cblktab[cblknum].total += newcost - oldcost;

    return newcost - oldcost;
}


PASTIX_INT countBlok(PASTIX_INT cblknum, SymbolMatrix *symbptr, PASTIX_INT blcolmin)
{
  PASTIX_INT i;
  PASTIX_INT bloknbr;
  double delta;
  double stride = 0;
  delta = (double)(symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum+1);
  delta = ceil(delta/blcolmin);

  for(i=symbptr->cblktab[cblknum].bloknum+1;i<symbptr->cblktab[cblknum+1].bloknum;i++)
    stride += symbptr->bloktab[i].lrownum-symbptr->bloktab[i].frownum + 1;
  stride = ceil(stride/blcolmin);
  /*fprintf(stdout, "delta %g stride %g blcolmin %ld \n", delta, stride, blcolmin); */
  bloknbr = 0;
  bloknbr += (PASTIX_INT) (((delta + 1)*delta)/2);
  bloknbr += (PASTIX_INT) (stride*delta);

  return bloknbr;
}

PASTIX_INT setSubtreeBlokNbr(PASTIX_INT rootnum, const EliminTree *etree, SymbolMatrix *symbptr, ExtraSymbolMatrix *extrasymb, PASTIX_INT blcolmin)
{
  PASTIX_INT i;
  extrasymb->subtreeblnbr[rootnum] =  countBlok(rootnum, symbptr, blcolmin);
  /*fprintf(stdout, "Rootnum %ld bloknbr %ld \n", rootnum, extrasymb->blnbtab[rootnum]);*/
  for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
     extrasymb->subtreeblnbr[rootnum] += setSubtreeBlokNbr(TSON(etree, rootnum, i), etree, symbptr, extrasymb, blcolmin);
  return extrasymb->subtreeblnbr[rootnum];
}



PASTIX_INT check_candidat(SymbolMatrix *symbmtx, BlendCtrl *ctrl)
{
  PASTIX_INT i, j;
  PASTIX_INT facecblknum;
  for(i=0;i<symbmtx->cblknbr;i++)
    {
      for(j= symbmtx->cblktab[i].bloknum; j < symbmtx->cblktab[i+1].bloknum;j++)
        {
          facecblknum = symbmtx->bloktab[j].cblknum;
          if(ctrl->candtab[i].fcandnum < ctrl->candtab[facecblknum].fcandnum
             || ctrl->candtab[i].lcandnum > ctrl->candtab[facecblknum].lcandnum)
            {
              errorPrint("bad processor candidat sets : cblk %ld candidat =[%ld %ld] father %ld candidat = [%ld %ld].",
                         (long)i, (long)ctrl->candtab[i].fcandnum, (long)ctrl->candtab[i].lcandnum,
                         (long)facecblknum, (long)ctrl->candtab[facecblknum].fcandnum,
                         (long)ctrl->candtab[facecblknum].lcandnum);
              return -1;
            }
        }
    }
  return 1;
}
