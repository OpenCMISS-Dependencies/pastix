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

/*
**  The defines and includes.
*/

#define SOLVER_IO

#include "common_pastix.h"
#include "symbol.h"
#include "ftgt.h"
#include "queue.h"
#include "bulles.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "solverRealloc.h"
#include "solver_io.h"



/*****************************************/
/* The solver matrix handling routines.  */
/*                                       */
/*****************************************/

/*+ This routine reads the given
*** solver matrix structure from
*** the given stream. The solver matrix
*** structure must have been
*** previously initialized.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

PASTIX_INT solverLoad(SolverMatrix *solvptr, FILE *stream)
{
    PASTIX_INT i,j;
    PASTIX_INT clustnbr, clustnum;
    SolverCblk   *cblkptr;
    SolverCblk   *cblktnd;
    SolverBlok   *blokptr;
    SolverBlok   *bloktnd;
    FanInTarget  *ftgtptr;
    FanInTarget  *ftgttnd;
    BlockTarget  *btagptr;
    BlockTarget  *btagtnd;
    BlockCoeff   *bcofptr;
    BlockCoeff   *bcofnd;
    Task         *taskptr;
    Task         *tasknd;

    PASTIX_INT                 versval;
    PASTIX_INT                 baseval;
    PASTIX_INT                 nodenbr;
    PASTIX_INT                 cblknbr;
    PASTIX_INT                 cblknum;
    PASTIX_INT                 bloknbr;
    PASTIX_INT                 bloknum;

    solverInit(solvptr);

    /** Load the symbol matrix **/
    if ((intLoad (stream, &versval) +               /* Read header */
         intLoad (stream, &cblknbr) +
         intLoad (stream, &bloknbr) +
         intLoad (stream, &nodenbr) +
         intLoad (stream, &baseval) != 5) ||
        (versval < 0)                     ||        /* Version should be 0 or 1 */
        (versval > 1)                     ||
        (bloknbr < cblknbr)               ||
        (nodenbr < cblknbr)) {
      errorPrint ("solverLoad: bad input (1)");
      return     (1);
    }
    MALLOC_INTERN(solvptr->cblktab, cblknbr + 1, SolverCblk);
    MALLOC_INTERN(solvptr->bloktab, bloknbr,     SolverBlok);
    if (solvptr->cblktab == NULL || solvptr->bloktab == NULL) {
      errorPrint ("solverLoad: out of memory");
      solverExit (solvptr);
      solverInit (solvptr);
      return     (1);
    }
    solvptr->baseval = baseval;
    solvptr->cblknbr = cblknbr;
    solvptr->bloknbr = bloknbr;
    solvptr->nodenbr = nodenbr;

    for (cblknum = 0; cblknum < cblknbr; cblknum ++) {
      if ((intLoad (stream, &solvptr->cblktab[cblknum].fcolnum) + /* Read column blocks */
           intLoad (stream, &solvptr->cblktab[cblknum].lcolnum) +
           intLoad (stream, &solvptr->cblktab[cblknum].bloknum) != 3) ||
          (solvptr->cblktab[cblknum].fcolnum > solvptr->cblktab[cblknum].lcolnum)) {
        errorPrint ("solverLoad: bad input (2)");
        /* solverExit (solvptr); */
        /* solverInit (solvptr); */
        return     (1);
      }
    }
    solvptr->cblktab[cblknbr].fcolnum =             /* Set last column block */
      solvptr->cblktab[cblknbr].lcolnum = nodenbr + baseval;
    solvptr->cblktab[cblknbr].bloknum = bloknbr + baseval;

    for (bloknum = 0; bloknum < bloknbr; bloknum ++) {
      if ((intLoad (stream, &solvptr->bloktab[bloknum].frownum) + /* Read column blocks */
           intLoad (stream, &solvptr->bloktab[bloknum].lrownum) +
           intLoad (stream, &solvptr->bloktab[bloknum].cblknum) != 3) ||
          (solvptr->bloktab[bloknum].frownum > solvptr->bloktab[bloknum].lrownum)) {
        errorPrint ("solverLoad: bad input (3)");
        solverExit (solvptr);
        solverInit (solvptr);
        return     (1);
    }

    solvptr->bloktab[bloknum].levfval = 0;        /* Assume version 0 */
    if ((versval > 0) &&
        ((intLoad (stream, &solvptr->bloktab[bloknum].levfval) != 1) ||
         (solvptr->bloktab[bloknum].levfval < 0))) {
      errorPrint ("solverLoad: bad input (4)");
      solverExit (solvptr);
      solverInit (solvptr);
      return     (1);
    }
  }


    if(  intLoad (stream, &solvptr->coefnbr) +
         intLoad (stream, &solvptr->ftgtnbr) +
         intLoad (stream, &solvptr->coefmax) +
         intLoad (stream, &solvptr->bpftmax) +
         intLoad (stream, &solvptr->cpftmax) +
         intLoad (stream, &solvptr->nbftmax) +
         intLoad (stream, &solvptr->arftmax) +
         intLoad (stream, &clustnum) +
         intLoad (stream, &clustnbr) +
         intLoad (stream, &solvptr->btagnbr) +
         intLoad (stream, &solvptr->bcofnbr) +
         intLoad (stream, &solvptr->indnbr) +
         intLoad (stream, &solvptr->tasknbr) +
         intLoad (stream, &solvptr->procnbr) +
         intLoad (stream, &solvptr->thrdnbr) +
         intLoad (stream, &solvptr->gridldim) + intLoad (stream, &solvptr->gridcdim)
         != 13)
        {
            errorPrint ("solverLoad: bad input (1)");
            return     (1);
        }

    solvptr->clustnbr = (PASTIX_INT)clustnbr;
    solvptr->clustnum = (PASTIX_INT)clustnum;


    if (((solvptr->cblktab = (SolverCblk *) memAlloc ((solvptr->cblknbr + 1) * sizeof (SolverCblk)))  == NULL) ||
        ((solvptr->bloktab = (SolverBlok *) memAlloc (solvptr->bloknbr       * sizeof (SolverBlok)))  == NULL) ||
        ((solvptr->ftgttab = (FanInTarget *)memAlloc (solvptr->ftgtnbr               * sizeof (FanInTarget))) == NULL) ||
        ((solvptr->btagtab = (BlockTarget *)memAlloc (solvptr->btagnbr               * sizeof (BlockTarget))) == NULL) ||
        ((solvptr->bcoftab = (BlockCoeff *) memAlloc (solvptr->bcofnbr               * sizeof (BlockCoeff)))  == NULL) ||
        ((solvptr->indtab  = (PASTIX_INT *)        memAlloc (solvptr->indnbr                * sizeof (PASTIX_INT)))         == NULL) ||
        ((solvptr->tasktab = (Task *)       memAlloc ((solvptr->tasknbr+1)           * sizeof (Task)))        == NULL) ||
        ((solvptr->ttsknbr = (PASTIX_INT *)        memAlloc ((solvptr->thrdnbr)             * sizeof (PASTIX_INT)))         == NULL) ||
        ((solvptr->ttsktab = (PASTIX_INT **)       memAlloc ((solvptr->thrdnbr)             * sizeof (PASTIX_INT *)))       == NULL)
        ) {
        errorPrint ("solverLoad: out of memory (1)");
        if (solvptr->cblktab != NULL)
          memFree_null (solvptr->cblktab);
        if (solvptr->bloktab != NULL)
          memFree_null (solvptr->bloktab);
        if (solvptr->ftgttab != NULL)
          memFree_null (solvptr->ftgttab);
        if (solvptr->btagtab != NULL)
          memFree_null (solvptr->btagtab);
        if (solvptr->bcoftab != NULL)
          memFree_null (solvptr->bcoftab);
        if (solvptr->indtab != NULL)
          memFree_null (solvptr->indtab);
        if (solvptr->tasktab != NULL)
          memFree_null (solvptr->tasktab);
        return (1);
    }

    for (cblkptr = solvptr->cblktab,                /* Read column block data */
           cblktnd = cblkptr + solvptr->cblknbr;
         cblkptr < cblktnd; cblkptr ++){
      if (intLoad (stream, &cblkptr->stride) +
          intLoad (stream, &cblkptr->color) +
          intLoad (stream, &cblkptr->procdiag) +
          intLoad (stream, &cblkptr->cblkdiag) != 4)
          {
            errorPrint ("solverlLoad: bad input (2)");
            solverExit (solvptr);
            return     (1);
          }
#ifdef SOLVER_DEBUG
      /*intLoad (stream, &cblkptr->color);*/
#endif

      }

        for (blokptr = solvptr->bloktab,                /* Read block data */
             bloktnd = blokptr + solvptr->bloknbr;
             blokptr < bloktnd; blokptr ++) {
          if (intLoad (stream, &blokptr->coefind) != 1)
            {
              errorPrint ("solverLoad: bad input (3)");
              solverExit (solvptr);
              return     (1);
            }

        }


    for (ftgtptr = solvptr->ftgttab,                /* Read fan in target data */
           ftgttnd = ftgtptr + solvptr->ftgtnbr;
         ftgtptr < ftgttnd; ftgtptr ++)
      {
        for(i=0;i<MAXINFO;i++)
          intLoad (stream, &(ftgtptr->infotab[i]));
        ftgtptr->coeftab = NULL;
      }

    for (bcofptr = solvptr->bcoftab,                /* Read BlockCoeff data */
           bcofnd = bcofptr + solvptr->bcofnbr;
         (bcofptr < bcofnd); bcofptr ++)
      {
        for(i=0;i<BCOFINFO;i++)
          intLoad(stream, &(bcofptr->infotab[i]));
        bcofptr->coeftab = NULL;
        intLoad(stream, &(bcofptr->sendcnt));
      }

    for (btagptr = solvptr->btagtab,                /* Read BlockTarget data */
           btagtnd = btagptr + solvptr->btagnbr;
         (btagptr < btagtnd); btagptr ++)
      {
        PASTIX_INT bcofind;
        for(i=0;i<BTAGINFO;i++)
          intLoad(stream, &btagptr->infotab[i]);

        intLoad(stream, &bcofind);
        /*fprintf(stdout, "%d\t", bcofind);*/
        if(bcofind < 0)
          btagptr->bcofptr = NULL;
        else
          {
#ifdef DEBUG_BLEND
            ASSERT(bcofind < solvptr->bcofnbr,MOD_BLEND);
#endif
            btagptr->bcofptr = &(solvptr->bcoftab[bcofind]);
          }
      }



   for(i=0;i<solvptr->indnbr;i++)                   /** Read indtab **/
     intLoad(stream, &(solvptr->indtab[i]));


   for (taskptr = solvptr->tasktab,                /** Read Task data **/
        tasknd = taskptr + solvptr->tasknbr +1;
        (taskptr < tasknd); taskptr ++)
     {
       intLoad(stream, &(taskptr->taskid));
       intLoad(stream, &(taskptr->prionum));
       intLoad(stream, &(taskptr->cblknum));
       intLoad(stream, &(taskptr->bloknum));
       {
         /* volatile pb alpha */
         PASTIX_INT temp;
         intLoad(stream, &temp);
         taskptr->ftgtcnt = temp;
         intLoad(stream, &temp);
         taskptr->ctrbcnt = temp;
       }
       intLoad(stream, &(taskptr->indnum));
       intLoad(stream, &(taskptr->tasknext));
       taskptr->btagptr = NULL;
     }

   for(i=0;i<solvptr->thrdnbr;i++)                 /** Read task by thread data **/
     {
       intLoad(stream, &(solvptr->ttsknbr[i]));
       MALLOC_INTERN(solvptr->ttsktab[i], solvptr->ttsknbr[i], PASTIX_INT);
       if (solvptr->ttsktab[i] == NULL)
         {
           errorPrint ("solverLoad: out of memory (1)");
           return 1;
         }
       for (j=0;j<solvptr->ttsknbr[i];j++)
         {
           intLoad(stream, &(solvptr->ttsktab[i][j]));
         }
     }

   for(i=0;i<solvptr->procnbr;i++)                 /** Read proc -> cluster **/
     {
       intLoad(stream, &(solvptr->proc2clust[i]));
     }

   if(  intLoad (stream, &(solvptr->updovct.sm2xmax)) + /** Read updown **/
        intLoad (stream, &(solvptr->updovct.sm2xsze)) +
        intLoad (stream, &(solvptr->updovct.sm2xnbr)) +
        intLoad (stream, &(solvptr->updovct.gcblk2listnbr)) +
        intLoad (stream, &(solvptr->updovct.listptrnbr)) +
        intLoad (stream, &(solvptr->updovct.listnbr)) +
        intLoad (stream, &(solvptr->updovct.loc2globnbr)) +
        intLoad (stream, &(solvptr->updovct.gcblknbr)) +
        intLoad (stream, &(solvptr->updovct.gnodenbr))
        != 9)
     {
       errorPrint ("solverLoad: bad input (1)");
       return     (1);
     }

   MALLOC_INTERN(solvptr->updovct.cblktab,    solvptr->cblknbr,       UpDownCblk);
   MALLOC_INTERN(solvptr->updovct.gcblk2list, solvptr->updovct.gcblk2listnbr, PASTIX_INT);
   MALLOC_INTERN(solvptr->updovct.listptr,    solvptr->updovct.listptrnbr,    PASTIX_INT);
   MALLOC_INTERN(solvptr->updovct.listcblk,   solvptr->updovct.listnbr,       PASTIX_INT);
   MALLOC_INTERN(solvptr->updovct.listblok,   solvptr->updovct.listnbr,       PASTIX_INT);
   MALLOC_INTERN(solvptr->updovct.loc2glob,   solvptr->updovct.loc2globnbr,   PASTIX_INT);
   MALLOC_INTERN(solvptr->updovct.lblk2gcblk, solvptr->bloknbr,       PASTIX_INT);

   for (i=0;i<solvptr->cblknbr;i++)
     {
       intLoad(stream, &(solvptr->updovct.cblktab[i].sm2xind));
       intLoad(stream, &(solvptr->updovct.cblktab[i].browprocnbr));
       intLoad(stream, &(solvptr->updovct.cblktab[i].msgnbr));
       {
         PASTIX_INT msgcnt = solvptr->updovct.cblktab[i].msgcnt;
         intLoad(stream, &msgcnt);
       }
       intLoad(stream, &(solvptr->updovct.cblktab[i].ctrbnbr));
       {
         PASTIX_INT ctrbcnt = solvptr->updovct.cblktab[i].ctrbcnt;
         intLoad(stream, &ctrbcnt);
       }

       MALLOC_INTERN(solvptr->updovct.cblktab[i].browproctab,
                     solvptr->updovct.cblktab[i].browprocnbr,
                     PASTIX_INT);
       MALLOC_INTERN(solvptr->updovct.cblktab[i].browcblktab,
                     solvptr->updovct.cblktab[i].browprocnbr,
                     PASTIX_INT);

       for (j=0;j<solvptr->updovct.cblktab[i].browprocnbr;j++)
         {
           intLoad(stream, &(solvptr->updovct.cblktab[i].browproctab[j]));
           intLoad(stream, &(solvptr->updovct.cblktab[i].browcblktab[j]));
         }
     }

   for (i=0;i<solvptr->updovct.gcblk2listnbr;i++)
     {
       intLoad(stream, &(solvptr->updovct.gcblk2list[i]));
     }

   for (i=0;i<solvptr->updovct.listptrnbr;i++)
     {
       intLoad(stream, &(solvptr->updovct.listptr[i]));
     }

   for (i=0;i<solvptr->updovct.listnbr;i++)
     {
       intLoad(stream, &(solvptr->updovct.listcblk[i]));
       intLoad(stream, &(solvptr->updovct.listblok[i]));
     }

   for (i=0;i<solvptr->updovct.loc2globnbr;i++)
     {
       intLoad(stream, &(solvptr->updovct.loc2glob[i]));
     }

   for (i=0;i<solvptr->bloknbr;i++)
     {
       intLoad(stream, &(solvptr->updovct.lblk2gcblk[i]));
     }

   return (0);
}


PASTIX_INT solverSave(const SolverMatrix * solvptr, FILE *stream)
{
   PASTIX_INT i;
   PASTIX_INT j;
   SolverCblk   *cblkptr;
   SolverCblk   *cblktnd;
   SolverBlok   *blokptr;
   SolverBlok   *bloktnd;
   FanInTarget  *ftgtptr;
   FanInTarget  *ftgttnd;
   BlockTarget  *btagptr;
   BlockTarget  *btagtnd;
   BlockCoeff   *bcofptr;
   BlockCoeff   *bcofnd;
   Task         *taskptr;
   Task         *tasknd;

   PASTIX_INT          o;


   /** Save the solver matrix **/
   {
     const SolverCblk *  cblktnd;
     const SolverCblk *  cblkptr;
     const SolverBlok *  bloktnd;
     const SolverBlok *  blokptr;
     int                 o;

     o = (fprintf (stream, "1\n%ld\t%ld\t%ld\t%ld\n", /* Write file header */
                   (long) solvptr->cblknbr,
                   (long) solvptr->bloknbr,
                   (long) solvptr->nodenbr,
                   (long) solvptr->baseval) == EOF);
     for (cblkptr = solvptr->cblktab, cblktnd = cblkptr + solvptr->cblknbr;
          (cblkptr < cblktnd) && (o == 0); cblkptr ++) {
       o = (fprintf (stream, "%ld\t%ld\t%ld\n",
                     (long) cblkptr->fcolnum,
                     (long) cblkptr->lcolnum,
                     (long) cblkptr->bloknum) == EOF);
     }
     for (blokptr = solvptr->bloktab, bloktnd = blokptr + solvptr->bloknbr;
          (blokptr < bloktnd) && (o == 0); blokptr ++) {
       o = (fprintf (stream, "%ld\t%ld\t%ld\t%ld\n",
                     (long) blokptr->frownum,
                     (long) blokptr->lrownum,
                     (long) blokptr->cblknum,
                     (long) blokptr->levfval) == EOF);
     }

   }



   /*fprintf(stream, "File header\n");*/
   o = (fprintf (stream, "\n%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", /* Write file header */
                 (long) solvptr->coefnbr,
                 (long) solvptr->ftgtnbr,
                 (long) solvptr->coefmax,
                 (long) solvptr->bpftmax,
                 (long) solvptr->cpftmax,
                 (long) solvptr->nbftmax,
                 (long) solvptr->arftmax,
                 (long) solvptr->clustnum,
                 (long) solvptr->clustnbr,
                 (long) solvptr->btagnbr,
                 (long) solvptr->bcofnbr,
                 (long) solvptr->indnbr,
                 (long) solvptr->tasknbr,
                 (long) solvptr->procnbr,
                 (long) solvptr->thrdnbr,
                 (long) solvptr->gridldim,
                 (long) solvptr->gridcdim
                 ) == EOF);
   /* write cblk data */
   /*fprintf(stream, "cblk data\n");*/
   for (cblkptr = solvptr->cblktab, cblktnd = cblkptr + solvptr->cblknbr;
        (cblkptr < cblktnd) && (o == 0); cblkptr ++)
     {
       o = (fprintf (stream, "%ld\t%ld\t%ld\t%ld\n",
                     (long) cblkptr->stride,
                     (long) cblkptr->color,
                     (long) cblkptr->procdiag,
                     (long) cblkptr->cblkdiag) == EOF);
#ifdef SOLVER_DEBUG
       /*fprintf(stream, "%ld\n", (long)cblkptr->color);*/
#endif
     }


   /* write blok data */
   /*fprintf(stream, "blok data\n");*/
   for (blokptr = solvptr->bloktab, bloktnd = blokptr + solvptr->bloknbr;
        (blokptr < bloktnd) && (o == 0); blokptr ++) {
     o = (fprintf (stream, "%ld\n",(long) blokptr->coefind) == EOF);

   }

   /*fprintf(stream, "fan in target data\n");*/
   for (ftgtptr = solvptr->ftgttab,                /* Write fan in target data */
        ftgttnd = ftgtptr + solvptr->ftgtnbr;
        (ftgtptr < ftgttnd) && (o==0); ftgtptr ++) {
     for(i=0;i<MAXINFO;i++)
       o = (fprintf(stream, "%ld\t", (long)ftgtptr->infotab[i]) == EOF);
     fprintf(stream, "\n");
     fprintf(stream, "\n");
   }

   /*fprintf(stream, "blockcoeff data\n");*/
   for (bcofptr = solvptr->bcoftab,                /* Write BlockCoeff data !!BEFORE BLOCKTARGET FOR LOADING!!  */
          bcofnd = bcofptr + solvptr->bcofnbr;
        (bcofptr < bcofnd) && (o==0); bcofptr ++) {
     for(i=0;i<BCOFINFO;i++)
       o = (fprintf(stream, "%ld\t", (long)bcofptr->infotab[i]) == EOF);

     o = (fprintf(stream, "%ld\t", (long)bcofptr->sendcnt) == EOF);

     fprintf(stream, "\n");
     fprintf(stream, "\n");
   }


   /*fprintf(stream, "blocktarget data\n");*/
   for (btagptr = solvptr->btagtab,                /* Write BlockTarget data */
        btagtnd = btagptr + solvptr->btagnbr;
        (btagptr < btagtnd) && (o==0); btagptr ++) {
     for(i=0;i<BTAGINFO;i++)
       o = (fprintf(stream, "%ld\t", (long)btagptr->infotab[i]) == EOF);
     if(btagptr->bcofptr == NULL)
       {
         fprintf(stream, "-1\t");
         /*fprintf(stdout, "-1\t");*/
       }
     else
       {
         fprintf(stream, "%ld\t", (long)( ((long)btagptr->bcofptr-(long)solvptr->bcoftab)/sizeof(BlockCoeff) ));
         /*fprintf(stdout, "%ld\t", (long)( ((long)btagptr->bcofptr-(long)solvptr->bcoftab)/sizeof(BlockCoeff) ));*/

       }


     fprintf(stream, "\n");
     fprintf(stream, "\n");
   }


   /*fprintf(stream, "indtab data\n");*/
   for(i=0;i<solvptr->indnbr;i++)                   /** Write indtab **/
     fprintf(stream, "%ld\t", (long)solvptr->indtab[i]);
   fprintf(stream, "\n");
   fprintf(stream, "\n");

   /*fprintf(stream, "task data\n");*/
   for (taskptr = solvptr->tasktab,                /* Write Task data */
        tasknd = taskptr + solvptr->tasknbr+1;
        (taskptr < tasknd) && (o==0); taskptr ++)
     {
       fprintf(stream, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",
               (long)taskptr->taskid, (long)taskptr->prionum, (long)taskptr->cblknum, (long)taskptr->bloknum,
               (long)taskptr->ftgtcnt, (long)taskptr->ctrbcnt, (long)taskptr->indnum, (long)taskptr->tasknext);
       fprintf(stream, "\n");
       fprintf(stream, "\n");
     }

   /*fprintf(stream, "ttsktab\n");*/
   for (i=0; i<solvptr->thrdnbr; i++) /* Write ttsktab */
     {
       fprintf(stream, "%ld\n", (long)solvptr->ttsknbr[i]);
       for (j=0; j<solvptr->ttsknbr[i]; j++)
         {
           fprintf(stream, "%ld\n", (long)solvptr->ttsktab[i][j]);
         }
     }

   /*fprintf(stream, "proc2clust\n");*/
   for (i=0; i<solvptr->procnbr; i++) /* Write proc2clust */
     {
       fprintf(stream, "%ld\n", (long)solvptr->proc2clust[i]);
     }

   /*fprintf(stream, "updo\n");*/
   fprintf(stream, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",
           (long) solvptr->updovct.sm2xmax, (long) solvptr->updovct.sm2xsze, (long) solvptr->updovct.sm2xnbr,
           (long) solvptr->updovct.gcblk2listnbr, (long) solvptr->updovct.listptrnbr, (long) solvptr->updovct.listnbr,
           (long) solvptr->updovct.loc2globnbr, (long) solvptr->updovct.gcblknbr, (long) solvptr->updovct.gnodenbr);

   /*fprintf(stream, "updown cblk\n");*/
   for (i=0; i<solvptr->cblknbr; i++)
     {
       fprintf(stream, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",
               (long) solvptr->updovct.cblktab[i].sm2xind, (long) solvptr->updovct.cblktab[i].browprocnbr, (long) solvptr->updovct.cblktab[i].msgnbr,
               (long) solvptr->updovct.cblktab[i].msgcnt, (long) solvptr->updovct.cblktab[i].ctrbnbr, (long) solvptr->updovct.cblktab[i].ctrbcnt);

       for (j=0; j<solvptr->updovct.cblktab[i].browprocnbr; j++)
         {
           fprintf(stream, "%ld\t%ld\t\n",
                   (long) solvptr->updovct.cblktab[i].browproctab[j],
                   (long) solvptr->updovct.cblktab[i].browcblktab[j]);
         }
     }

   /*fprintf(stream, "updown gcblk2list\n");*/
   for (i=0; i<solvptr->updovct.gcblk2listnbr; i++)
     {
       fprintf(stream, "%ld\n", (long) solvptr->updovct.gcblk2list[i]);
     }

   /*fprintf(stream, "updown listptr\n");*/
   for (i=0; i<solvptr->updovct.listptrnbr; i++)
     {
       fprintf(stream, "%ld\n", (long) solvptr->updovct.listptr[i]);
     }

   /*fprintf(stream, "updown listcblk & listblok\n");*/
   for (i=0; i<solvptr->updovct.listnbr; i++)
     {
       fprintf(stream, "%ld\t%ld\n", (long) solvptr->updovct.listcblk[i], (long) solvptr->updovct.listblok[i]);
     }

   /*fprintf(stream, "updown loc2globcblk\n");*/
   for (i=0; i<solvptr->updovct.loc2globnbr; i++)
     {
       fprintf(stream, "%ld\n", (long) solvptr->updovct.loc2glob[i]);
     }

   /*fprintf(stream, "updown lblk2gcblk\n");*/
   for (i=0; i<solvptr->bloknbr; i++)
     {
       fprintf(stream, "%ld\n", (long) solvptr->updovct.lblk2gcblk[i]);
     }


   return o;
 }
