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
#include "ftgt.h"
#include "symbol.h"
#include "extrastruct.h"
#include "elimin.h"
#include "extendVector.h"
#include "cand.h"
#include "queue.h"
#include "bulles.h"
#include "param_blend.h"
#include "blendctrl.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "simu.h"
#include "costfunc.h"
#include "partbuild.h"


/*+
  Make a new SymbolMatrix and CostMatrix from the former ones
  and the Extra ones (that contains splitted bloks)
  +*/
/* OIMBE : retravailler partBuild avec des memcpy */

void partBuild(SymbolMatrix *symbmtx, ExtraSymbolMatrix *extrasymb,
               CostMatrix *costmtx, ExtraCostMatrix *extracost,
               BlendCtrl *ctrl, const Dof * dofptr)
{
    PASTIX_INT i, j, k, s;
    PASTIX_INT curbloknum;
    PASTIX_INT sptbloknum;
    PASTIX_INT *newnum      = NULL;
    PASTIX_INT *extranewnum = NULL;
    SymbolMatrix *tmp;
    CostMatrix   *tmp2;
    Cand         *tmp3;
    Cand         *candtab;
    PASTIX_INT    facing_splitted_cnt = 0;

    /* No splitted cblk: partition remains the same */
    if(extrasymb->curcblk == 0)
        return;

    MALLOC_INTERN(tmp,  1, SymbolMatrix);
    MALLOC_INTERN(tmp2, 1, CostMatrix);
    symbolInit(tmp);
    costInit(tmp2);

    tmp->baseval  = symbmtx->baseval;
    tmp->cblknbr  = symbmtx->cblknbr;
    tmp->bloknbr  = symbmtx->bloknbr;
    tmp->cblktab  = symbmtx->cblktab;
    tmp->bloktab  = symbmtx->bloktab;

    tmp2->cblktab = costmtx->cblktab;
    tmp2->bloktab = costmtx->bloktab;

    tmp3  = ctrl->candtab;

    symbmtx->cblknbr += extrasymb->addcblk;
    symbmtx->bloknbr += extrasymb->addblok;

    if (ctrl->option->iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
      {
        fprintf(stdout, "Number of column blocks created by splitting : %d\n", (int)(extrasymb->addcblk));
        fprintf(stdout, "Number of blocks creating by splitting       : %d\n", (int)(extrasymb->addblok));
      }
    /** allocate space for the new matrix **/
    MALLOC_INTERN(symbmtx->cblktab, symbmtx->cblknbr+1, SymbolCblk);
    MALLOC_INTERN(symbmtx->bloktab, symbmtx->bloknbr,   SymbolBlok);
    MALLOC_INTERN(costmtx->cblktab, symbmtx->cblknbr+1, CostCblk);
    MALLOC_INTERN(costmtx->bloktab, symbmtx->bloknbr,   CostBlok);
    MALLOC_INTERN(ctrl->candtab,    symbmtx->cblknbr,   Cand);

    candtab = ctrl->candtab;

    /** we use sptcbnb to get new num of cblk in the new symbolic matrix **/
    MALLOC_INTERN(newnum, tmp->cblknbr+1, PASTIX_INT);
    memcpy(&(newnum[1]), extrasymb->sptcbnb, tmp->cblknbr * sizeof(PASTIX_INT));
    newnum[0] = 0;
    for(i=1;i < tmp->cblknbr+1;i++)
        newnum[i] += newnum[i-1];
    /** newnum[i] is the new decomp number of the first splitted cblk from former
        cblk number i**/
#ifdef DEBUB_BLEND
    for(i=0;i<tmp->cblknbr;i++)
        ASSERT((newnum[i]>=0) && (newnum[i]<symbmtx->cblknbr),MOD_BLEND);
#endif


    /** Now, we use sptcblk and newind to get the new decomp number of all cblk owned
        by the extra symbolic matrix **/
    MALLOC_INTERN(extranewnum, extrasymb->curcblk, PASTIX_INT);
    for(i=0;i<tmp->cblknbr;i++)
        if(extrasymb->sptcblk[i]>=0)
            for(j=0;j<extrasymb->sptcbnb[i];j++)
                extranewnum[extrasymb->sptcblk[i]+j] = newnum[i]+j;
    /** extranewnum[i] is the decomp number of the cblk[i] of the extra symbolic
        matrix */

#ifdef DEBUB_BLEND
    for(i=0;i<extrasymb->curcblk;i++)
        ASSERT((extranewnum[i]>=0) && (extranewnum[i]<symbmtx->cblknbr),MOD_BLEND);
#endif



    /** fill the new symbolic matrix resulting from splitting of the former one **/
    curbloknum = 0;
    for(i=0;i<tmp->cblknbr;i++)
        {
            if(extrasymb->sptcblk[i] < 0) /* not a splitted cblk */
                {
                    symbmtx->cblktab[newnum[i]].fcolnum = tmp->cblktab[i].fcolnum;
                    symbmtx->cblktab[newnum[i]].lcolnum = tmp->cblktab[i].lcolnum;
                    symbmtx->cblktab[newnum[i]].bloknum = curbloknum;
#ifdef STARPU_GET_TASK_CTX
                    symbmtx->cblktab[newnum[i]].ctx     = tmp->cblktab[i].ctx;
#endif

                    /* no need to copy subtree cost, we'll have to update it */
                    costmtx->cblktab[newnum[i]].total     = tmp2->cblktab[i].total;
                    costmtx->cblktab[newnum[i]].compute   = tmp2->cblktab[i].compute;
                    costmtx->cblktab[newnum[i]].send      = tmp2->cblktab[i].send;

                    candtab[newnum[i]].treelevel = tmp3[i].treelevel;
                    candtab[newnum[i]].costlevel = tmp3[i].costlevel;

                    candtab[newnum[i]].fcandnum  = tmp3[i].fcandnum;
                    candtab[newnum[i]].lcandnum  = tmp3[i].lcandnum;
                    candtab[newnum[i]].cluster   = tmp3[i].cluster;
                    candtab[newnum[i]].distrib   = tmp3[i].distrib;


                    for(j=tmp->cblktab[i].bloknum; j<tmp->cblktab[i+1].bloknum;j++)
                        {
                            if(extrasymb->sptblok[j] < 0) /* not a splitted blok
                                                             so its facing diag is not splitted */
                                /* NB: even if a blok is not splitted while its facing cblk is
                                   splitted , it's considered as splitted */
                                {
                                    symbmtx->bloktab[curbloknum].frownum = tmp->bloktab[j].frownum;
                                    symbmtx->bloktab[curbloknum].lrownum = tmp->bloktab[j].lrownum;
                                    symbmtx->bloktab[curbloknum].cblknum = newnum[tmp->bloktab[j].cblknum];
#ifdef DEBUG_BLEND
                                    ASSERT((newnum[tmp->bloktab[j].cblknum] >=0)&&(newnum[tmp->bloktab[j].cblknum] <= symbmtx->cblknbr),MOD_BLEND);
#endif
                                    costmtx->bloktab[curbloknum].contrib = tmp2->bloktab[j].contrib;
                                    costmtx->bloktab[curbloknum].linenbr = tmp2->bloktab[j].linenbr;
                                    curbloknum++;
                                }
                            else      /* splitted blok in a non splitted cblk
                                         -> the facing diagblok is splitted */
                                {
                                    facing_splitted_cnt += extrasymb->sptblnb[j]-1;
                                    for(k=extrasymb->sptblok[j];k < extrasymb->sptblok[j]+extrasymb->sptblnb[j];k++)
                                        {
                                            symbmtx->bloktab[curbloknum].frownum = extrasymb->bloktab[k].frownum;
                                            symbmtx->bloktab[curbloknum].lrownum = extrasymb->bloktab[k].lrownum;
                                            symbmtx->bloktab[curbloknum].cblknum = extranewnum[extrasymb->bloktab[k].cblknum];
#ifdef DEBUG_BLEND
                                            ASSERT((extranewnum[extrasymb->bloktab[k].cblknum] >=0)&&(extranewnum[extrasymb->bloktab[k].cblknum] <= symbmtx->cblknbr),MOD_BLEND);
#endif
                                            costmtx->bloktab[curbloknum].contrib = extracost->bloktab[k].contrib;
                                            costmtx->bloktab[curbloknum].linenbr = extracost->bloktab[k].linenbr;

                                            curbloknum++;
                                        }
                                }
                        }
                }
            else    /* splitted cblk */
                {
                    for(j=extrasymb->sptcblk[i]; j < extrasymb->sptcblk[i]+extrasymb->sptcbnb[i];j++)
                        {
                            symbmtx->cblktab[extranewnum[j]].fcolnum = extrasymb->cblktab[j].fcolnum;
                            symbmtx->cblktab[extranewnum[j]].lcolnum = extrasymb->cblktab[j].lcolnum;
#ifdef STARPU_GET_TASK_CTX
                            symbmtx->cblktab[extranewnum[j]].ctx = extrasymb->cblktab[j].ctx;
#endif
                            symbmtx->cblktab[extranewnum[j]].bloknum = curbloknum;

                            candtab[extranewnum[j]].treelevel = tmp3[i].treelevel;
                            candtab[extranewnum[j]].costlevel = tmp3[i].costlevel;
                            candtab[extranewnum[j]].fcandnum  = tmp3[i].fcandnum;
                            candtab[extranewnum[j]].lcandnum  = tmp3[i].lcandnum;
                            candtab[extranewnum[j]].cluster   = tmp3[i].cluster;
                            candtab[extranewnum[j]].distrib   = tmp3[i].distrib;

                            /** treat blok created by splitting of the diag blok **/
                            for(k=extrasymb->cblktab[j].bloknum;
                                k < (extrasymb->cblktab[j].bloknum +extrasymb->sptcblk[i] +extrasymb->sptcbnb[i]-j);k++)
                                {
                                    symbmtx->bloktab[curbloknum].frownum = extrasymb->bloktab[k].frownum;
                                    symbmtx->bloktab[curbloknum].lrownum = extrasymb->bloktab[k].lrownum;
                                    symbmtx->bloktab[curbloknum].cblknum = extranewnum[extrasymb->bloktab[k].cblknum];
#ifdef DEBUG_BLEND
                                    ASSERT(symbmtx->bloktab[curbloknum].frownum >= extrasymb->cblktab[extrasymb->bloktab[k].cblknum].fcolnum,MOD_BLEND);
                                    ASSERT(symbmtx->bloktab[curbloknum].lrownum <= extrasymb->cblktab[extrasymb->bloktab[k].cblknum].lcolnum,MOD_BLEND);
#endif
                                    curbloknum++;
                                }

                            sptbloknum = k;
                            for(k=tmp->cblktab[i].bloknum+1;k<tmp->cblktab[i+1].bloknum;k++)
                                {
                                    if(extrasymb->sptblok[k]<0)
                                        {
                                            symbmtx->bloktab[curbloknum].frownum = tmp->bloktab[k].frownum;
                                            symbmtx->bloktab[curbloknum].lrownum = tmp->bloktab[k].lrownum;
                                            symbmtx->bloktab[curbloknum].cblknum = newnum[extrasymb->bloktab[sptbloknum].cblknum];
                                            sptbloknum++;
                                            curbloknum++;
                                        }
                                    else
                                        {
                                            facing_splitted_cnt += extrasymb->sptblnb[k]-1;
                                            for(s=extrasymb->sptblok[k];s<extrasymb->sptblok[k]+extrasymb->sptblnb[k];s++)
                                                {
                                                    symbmtx->bloktab[curbloknum].frownum = extrasymb->bloktab[s].frownum;
                                                    symbmtx->bloktab[curbloknum].lrownum = extrasymb->bloktab[s].lrownum;
                                                    symbmtx->bloktab[curbloknum].cblknum = extranewnum[extrasymb->bloktab[s].cblknum];
                                                    sptbloknum++;
                                                    curbloknum++;
                                                }
                                        }
                                }

                        }

                }

        }

#ifdef DEBUG_BLEND
    ASSERT(curbloknum == symbmtx->bloknbr,MOD_BLEND);
#endif
    /* virtual cblk to avoid side effect in the loops on cblk bloks */
    symbmtx->cblktab[symbmtx->cblknbr].fcolnum = symbmtx->cblktab[symbmtx->cblknbr-1].lcolnum+1;
    symbmtx->cblktab[symbmtx->cblknbr].lcolnum = symbmtx->cblktab[symbmtx->cblknbr-1].lcolnum+1;
    symbmtx->cblktab[symbmtx->cblknbr].bloknum = curbloknum;

    costmtx->cblktab[symbmtx->cblknbr].total   = 0;
    costmtx->cblktab[symbmtx->cblknbr].compute = 0;
    costmtx->cblktab[symbmtx->cblknbr].send    = 0;


    /** we have to compute the cost of a splitted cblk **/
    for(i=0;i<tmp->cblknbr;i++)
        {
            if(extrasymb->sptcblk[i] >= 0) /* not a splitted cblk */
                for(j=extrasymb->sptcblk[i]; j < extrasymb->sptcblk[i]+extrasymb->sptcbnb[i];j++)
                    cblkComputeCost(extranewnum[j], costmtx, symbmtx, dofptr);
        }


    if (ctrl->option->iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
      {
        double        block_height_sum = 0.0;
        double        cblk_width_sum = 0.0;
        for (j = 0; j < symbmtx->cblknbr; j++)
          {
            cblk_width_sum += (double)(symbmtx->cblktab[j].lcolnum - symbmtx->cblktab[j].fcolnum + 1);

            for (i = symbmtx->cblktab[j].bloknum+1; i < symbmtx->cblktab[j+1].bloknum; i++)
              {
                block_height_sum += (double)(symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum + 1);
              }
          }
        fprintf(stdout, "Average cblk size : %g\n", cblk_width_sum/symbmtx->cblknbr);
        fprintf(stdout, "Average extra diagonal block height : %g\n", block_height_sum/(symbmtx->bloknbr-symbmtx->cblknbr));
        fprintf(stdout, "Number of blocks created due to facing block splitting : %d\n", (int)facing_splitted_cnt);
      }

    /** Free memory **/
    memFree_null(newnum);
    memFree_null(extranewnum);
    symbolExit(tmp);
    memFree_null(tmp);
    costExit(tmp2);
    memFree_null(tmp3);
}
