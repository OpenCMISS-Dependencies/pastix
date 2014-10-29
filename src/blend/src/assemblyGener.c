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
#include <assert.h>

#include "common_pastix.h"
#include "dof.h"
#include "extendVector.h"
#include "elimin.h"
#include "cost.h"
/* #include "ftgt.h" */
#include "symbol.h"
/* #include "updown.h" */
/* #include "solver.h" */
#include "cand.h"
#include "queue.h"
#include "bulles.h"
/* #include "simu.h" */
#include "param_blend.h"
#include "blendctrl.h"
#include "assembly.h"
#include "assemblyGener.h"

void assemblyGener(PASTIX_INT clustnum, Assembly1D *assemb1D, Assembly2D *assemb2D,
                   PASTIX_INT clustnbr, const SymbolMatrix *symbmtx, const PASTIX_INT *blprtab,
                   BlendCtrl *ctrl, const Dof * const dofptr)
/**************************************************************************/
/*  Function to gener assembly structure                                  */
/* IN SMP procnbr == clustnbr ?? (modifier blprtab dans ce fichier)       */
/* ATTENTION blprtab --> blclusttab                                       */
/**************************************************************************/
{
    PASTIX_INT i, j, p;
    PASTIX_INT *localnbr       = NULL;
    PASTIX_INT *proc2rownbr    = NULL;
    PASTIX_INT pr;
    PASTIX_INT maxrow;
    double nlcoefnbr; /** Non local coeff nbr in 2D distribution / assem 1D distribution**/
    double nlbloknbr; /** Non local block nbr in 2D distribution / assem 1D distribution**/
    double lcoefnbr; /** Local coeff nbr in 2D distribution / assem 1D distribution**/
    double lbloknbr; /** Local block nbr in 2D distribution / assem 1D distribution**/
    PASTIX_INT *cbprtab        = NULL;    /** cblk (in global num) to processor owner in 1D distribution**/

    PASTIX_INT *bloklocalnum1D = NULL; /** Global to local 1D blocknum **/
    PASTIX_INT *cblklocalnum1D = NULL; /** Global to local 2D cblknum if cblk is local in the 2D distribution **/
    PASTIX_INT *cblklocalnum2D = NULL; /** Global to local 2D cblknum if cblk is local in the 2D distribution **/
    PASTIX_INT *bloklocalnum2D = NULL; /** Global to local 2D blocknum **/
    PASTIX_INT clustid;
    PASTIX_INT bloknbr1D;
    PASTIX_INT cblknbr1D;
    PASTIX_INT bloknbr2D;
    PASTIX_INT *bloknum        = NULL;
    PASTIX_INT *cblknum        = NULL;
    PASTIX_INT *blcltab        = NULL; /** blcltab[i] == Cluster owner for block i **/
    PASTIX_INT delta;


    MALLOC_INTERN(localnbr, clustnbr, PASTIX_INT);
    bzero(localnbr, clustnbr*sizeof(PASTIX_INT));

    /** Convert blprtab in blcltab **/
    MALLOC_INTERN(blcltab, symbmtx->bloknbr, PASTIX_INT);
    for(i=0;i<symbmtx->bloknbr;i++)
      {
#ifdef DEBUG_BLEND
        ASSERT(blprtab[i]>=0,MOD_BLEND);
        ASSERT(blprtab[i]<ctrl->procnbr,MOD_BLEND);
#endif
        blcltab[i] = ctrl->proc2clust[blprtab[i]];
        /*fprintf(stderr, "%ld : blprtab %ld  blcltab %ld \n", (long)i, (long)blprtab[i], (long)blcltab[i]);*/
      }

    /** Fill cbprtab **/
    MALLOC_INTERN(cbprtab, symbmtx->cblknbr, PASTIX_INT);

    /***************************************************************************************************************/
    /*  if a cblk is mapped with a 1D distribution the processor owner is the same than the one decided in the     */
    /*  factorization distribution                                                                                 */
    /*  if a cblk is mapped with a 2D distribution then the processor owner in assembly is the processor that      */
    /* owns the biggest area of block in the cblk                                                                  */
    /***************************************************************************************************************/
    MALLOC_INTERN(proc2rownbr, clustnbr, PASTIX_INT);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        if(ctrl->candtab[i].distrib == D1)
          /*cbprtab[i] = ctrl->proc2clust[blprtab[symbmtx->cblktab[i].bloknum]];*/
          cbprtab[i] = blcltab[symbmtx->cblktab[i].bloknum];
        else
          {
            bzero(proc2rownbr, sizeof(PASTIX_INT)*clustnbr);
            for(j=symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
              /*proc2rownbr[blprtab[j]] += symbmtx->bloktab[j].lrownum - symbmtx->bloktab[j].frownum + 1;*/
              proc2rownbr[blcltab[j]] += symbmtx->bloktab[j].lrownum - symbmtx->bloktab[j].frownum + 1;

            /** Find the processor that has the largest number of row in this cblk **/
            maxrow = -1;
            pr = -1;
            for(p=0;p<clustnbr;p++)
              if(proc2rownbr[p] > maxrow)
                {
                  maxrow = proc2rownbr[p];
                  pr = p;
                }
            cbprtab[i] = pr;
          }
      }
    memFree(proc2rownbr);


    /** Compute the local numbering for cblk, blok, for each proc **/
#ifdef DEBUG_M
    ASSERT(clustnbr > 0,MOD_BLEND);
#endif
    MALLOC_INTERN(bloknum, clustnbr, PASTIX_INT);
    MALLOC_INTERN(cblknum, clustnbr, PASTIX_INT);

    /** GLOBAL TO LOCAL ORDERING 1D **/
    MALLOC_INTERN(cblklocalnum1D, symbmtx->cblknbr, PASTIX_INT);
    MALLOC_INTERN(bloklocalnum1D, symbmtx->bloknbr, PASTIX_INT);

    bzero(bloknum, sizeof(PASTIX_INT)*clustnbr);
    bzero(cblknum, sizeof(PASTIX_INT)*clustnbr);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        clustid = cbprtab[i];
        cblklocalnum1D[i] = cblknum[clustid];
        cblknum[clustid]++;
        for(j = symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
          {
            bloklocalnum1D[j] = bloknum[clustid];
            bloknum[clustid]++;
          }
      }
    bloknbr1D = bloknum[clustnum];
    cblknbr1D = cblknum[clustnum];

    /** GLOBAL TO LOCAL ORDERING 2D **/
    MALLOC_INTERN(cblklocalnum2D, symbmtx->cblknbr, PASTIX_INT);
    MALLOC_INTERN(bloklocalnum2D, symbmtx->bloknbr, PASTIX_INT);


    bzero(bloknum, sizeof(PASTIX_INT)*clustnbr);
    bzero(cblknum, sizeof(PASTIX_INT)*clustnbr);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        PASTIX_INT flag;
        flag = 0;

        for(j = symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
          {
            /*clustid = blprtab[j];*/
            clustid = blcltab[j];
            if(clustid == clustnum)
              flag = 1; /** This is a local cblk in the 2D distribution on this processor **/

            bloklocalnum2D[j] = bloknum[clustid];
            bloknum[clustid]++;
          }
        if(flag == 1)
          {
            cblklocalnum2D[i] = cblknum[clustnum];
            cblknum[clustnum]++;
          }
        else
          cblklocalnum2D[i] = -1;
      }
    bloknbr2D = bloknum[clustnum];

    memFree(cblknum);
    memFree(bloknum);


    /** Fill assemb1D->blprtab:  ATTENTION for the 1D distribution blprtab means cbprtab **/
    MALLOC_INTERN(assemb1D->blprtab, symbmtx->cblknbr, PASTIX_INT);
    memCpy(assemb1D->blprtab, cbprtab, sizeof(PASTIX_INT)*symbmtx->cblknbr);


    /** Fill nocbtab **/
    MALLOC_INTERN(assemb1D->nocbtab, symbmtx->nodenbr, PASTIX_INT);
    for(i=0;i<symbmtx->cblknbr;i++)
        for(j=symbmtx->cblktab[i].fcolnum;j<=symbmtx->cblktab[i].lcolnum;j++)
            assemb1D->nocbtab[j] = i;
#ifdef DEBUG_BLEND
    ASSERT(assemb1D->nocbtab[0] == 0,MOD_BLEND);
    ASSERT(assemb1D->nocbtab[symbmtx->nodenbr-1] == symbmtx->cblknbr-1,MOD_BLEND);
#endif

    /** Fill rnumtab **/
    MALLOC_INTERN(assemb1D->rnumtab, symbmtx->cblknbr, PASTIX_INT);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        assemb1D->rnumtab[i] = localnbr[cbprtab[i]];
        localnbr[cbprtab[i]]++;
      }
    memFree(localnbr);

    if(ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
      {
        /*** Estimated amount of non local block ***/
        nlcoefnbr = 0;
        nlbloknbr = 0;
        lcoefnbr = 0;
        lbloknbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
          {
            if(ctrl->candtab[i].distrib != D1)
              {
                delta = symbmtx->cblktab[i].lcolnum - symbmtx->cblktab[i].lcolnum + 1;
                for(j=symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
                  {

                    /*if(blprtab[j] != cbprtab[i])*/
                    if(blcltab[j] != cbprtab[i])
                      {
                        nlcoefnbr += (symbmtx->bloktab[j].lrownum - symbmtx->bloktab[j].frownum + 1)*delta;
                        nlbloknbr++;
                      }
                    else
                      {
                        lcoefnbr += (symbmtx->bloktab[j].lrownum - symbmtx->bloktab[j].frownum + 1)*delta;
                        lbloknbr++;
                      }
                  }
              }

          }

        if(lbloknbr>0)
          {
            fprintf(stdout, "In assembly: in 2D distributed column block:  %g percent coef are not local \n", nlcoefnbr*100/lcoefnbr);
            fprintf(stdout, "In assembly: in 2D distributed column block:  %g percent block are not local \n", nlbloknbr*100/lbloknbr);
          }
      }


    /**************************************************************************************************/
    /*** Generation of the assembly structure 2D                                                   ****/
    /**************************************************************************************************/
    MALLOC_INTERN(assemb2D->blok2proc_tab, bloknbr1D, PASTIX_INT);  /*+ local block i in 1D --> processor owner in 2D distribution +*/
    MALLOC_INTERN(assemb2D->blok2cblk_tab, bloknbr2D, PASTIX_INT);  /*+ local block i in 2D --> local cblk on the same processor in
                                                                                         the 2D distribution  +*/
    MALLOC_INTERN(assemb2D->blok2blok_tab, bloknbr1D, PASTIX_INT);  /*+ local block i in 1D --> local block i in the 2D distribution +*/

    for(i=0;i<symbmtx->cblknbr;i++)
      {
        if(cbprtab[i] == clustnum)
          for(j = symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
            {
              /*assemb2D->blok2proc_tab[bloklocalnum1D[j]] = blprtab[j]; */
              assemb2D->blok2proc_tab[bloklocalnum1D[j]] = blcltab[j];

              assemb2D->blok2blok_tab[bloklocalnum1D[j]] = bloklocalnum2D[j];
            }

        for(j = symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
          /*if(blprtab[j] == clustnum)*/
          if(blcltab[j] == clustnum)
            assemb2D->blok2cblk_tab[bloklocalnum2D[j]] = cblklocalnum2D[i];

      }
    /* Generation of the local symbol matrix distributed by column block for the assembly phase      */
    MALLOC_INTERN(assemb2D->symbmtx, 1, SymbolMatrix);
    symbolInit(assemb2D->symbmtx);
    symbolGener(clustnum, cblklocalnum1D, bloknbr1D, cblknbr1D, cbprtab, symbmtx, assemb2D->symbmtx, dofptr);

    memFree(cbprtab);
    memFree(blcltab);
    memFree(cblklocalnum1D);
    memFree(cblklocalnum2D);
    memFree(bloklocalnum1D);
    memFree(bloklocalnum2D);
}

void symbolGener(PASTIX_INT clustnum,
                 const PASTIX_INT *cblklocalnum1D,
                 PASTIX_INT bloknbr1D,
                 PASTIX_INT cblknbr1D,
                 const PASTIX_INT *cbprtab,
                 const SymbolMatrix *symbmtx,
                 SymbolMatrix *symb1D,
                 const Dof * const dofptr)
{
  PASTIX_INT i, j;
  PASTIX_INT cblknum, bloknum;
  PASTIX_INT nodenbr;

  symb1D->cblknbr = cblknbr1D;
  symb1D->bloknbr = bloknbr1D;

  /*************************/
  /**   Fill symb1D       **/
  /*************************/
  /* Allocations */
  MALLOC_INTERN(symb1D->cblktab, symb1D->cblknbr+1, SymbolCblk);
  MALLOC_INTERN(symb1D->bloktab, symb1D->bloknbr, SymbolBlok);

  cblknum = 0;
  bloknum = 0;
  nodenbr = 0;
  for(i=0;i<symbmtx->cblknbr;i++)
    {
      if(cbprtab[i] == clustnum)
        {
          symb1D->cblktab[cblknum].fcolnum = symbmtx->cblktab[i].fcolnum * dofptr->noddval;
          symb1D->cblktab[cblknum].lcolnum = symbmtx->cblktab[i].lcolnum * dofptr->noddval + dofptr->noddval-1;
          symb1D->cblktab[cblknum].bloknum = bloknum;
          nodenbr += symbmtx->cblktab[i].lcolnum - symbmtx->cblktab[i].fcolnum + 1;
          cblknum++;

          for(j=symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
            {
              symb1D->bloktab[bloknum].frownum = symbmtx->bloktab[j].frownum * dofptr->noddval;
              symb1D->bloktab[bloknum].lrownum = symbmtx->bloktab[j].lrownum * dofptr->noddval + dofptr->noddval-1;
              symb1D->bloktab[bloknum].cblknum = cblklocalnum1D[symbmtx->bloktab[j].cblknum];
              bloknum ++;
            }
        }
    }

  symb1D->nodenbr = nodenbr;
#ifdef DEBUG_BLEND
  ASSERT(symb1D->cblknbr == cblknum,MOD_BLEND);
  if(symb1D->bloknbr != bloknum)
    fprintf(stderr, "bloknbr %ld bloknum %ld \n", (long)symb1D->bloknbr, (long)bloknum);
  ASSERT(symb1D->bloknbr == bloknum,MOD_BLEND);
#endif

  /*  virtual cblk to avoid side effect in the loops on cblk bloks */
  symb1D->cblktab[cblknum].fcolnum = symb1D->cblktab[cblknum-1].lcolnum+1;
  symb1D->cblktab[cblknum].lcolnum = symb1D->cblktab[cblknum-1].lcolnum+1;
  symb1D->cblktab[cblknum].bloknum = bloknum;

}





#ifdef OLD_ASSEMBLY
/** AssemblyGener have to be used before solverMatrxGen (no expansion in ddl) **/
void assemblyGener(Assembly1D *assemb, PASTIX_INT procnbr, const SymbolMatrix *symbmtx, const PASTIX_INT *cbprtab)
/**************************************************************************/
/*  Function that was used in 1D to gener assmebly structure              */
/**************************************************************************/
{
    PASTIX_INT i, j;
    PASTIX_INT *localnbr;

    MALLOC_INTERN(localnbr, procnbr, PASTIX_INT);
    bzero(localnbr, procnbr*sizeof(PASTIX_INT));


    /** Fill blprtab **/
    /** IN 1D BL MEANS CBL **/
    MALLOC_INTERN(assemb->blprtab, symbmtx->cblknbr, PASTIX_INT);
    for(i=0;i<symbmtx->cblknbr;i++)
      assemb->blprtab[i] = cbprtab[i];

    /** Fill nocbtab **/
    MALLOC_INTERN(assemb->nocbtab, symbmtx->nodenbr, PASTIX_INT);
    for(i=0;i<symbmtx->cblknbr;i++)
        for(j=symbmtx->cblktab[i].fcolnum;j<=symbmtx->cblktab[i].lcolnum;j++)
            assemb->nocbtab[j] = i;

    /** Fill rnumtab **/
    MALLOC_INTERN(assemb->rnumtab, symbmtx->cblknbr, PASTIX_INT);
    for(i=0;i<symbmtx->cblknbr;i++)
      {
        assemb->rnumtab[i] = localnbr[cbprtab[i]];
        localnbr[cbprtab[i]]++;
      }
    memFree(localnbr);

}
#endif
