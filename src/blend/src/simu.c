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
#include <strings.h>
#include <math.h>

#include "common_pastix.h"
#include "ftgt.h"
#include "symbol.h"
#include "cand.h"
#include "queue.h"
#include "extendVector.h"
#include "simu.h"

#define TIMEBASE 10.0

PASTIX_INT simuInit(SimuCtrl *simuctrl, SymbolMatrix *symbptr, PASTIX_INT clustnbr, PASTIX_INT procnbr,
             PASTIX_INT cblknbr, PASTIX_INT bloknbr, Cand *candtab)
{
    PASTIX_INT i, j;
    PASTIX_INT p;
    PASTIX_INT ftgtcur;
    PASTIX_INT candnbr;
    PASTIX_INT step;

    simuctrl->cblknbr  = cblknbr;
    simuctrl->ftgtprio = 0;
    simuctrl->tasktab  = NULL;
    simuctrl->ftgtcnt  = 0;

    /** Processor initialisation **/
    MALLOC_INTERN(simuctrl->proctab, procnbr, SimuProc);
    for(i=0;i<procnbr;i++)
        {
          timerSet(TIMER(i), 0.0); /** for paragraph numeric tolerance **/
          simuctrl->proctab[i].prionum   = 0;
          MALLOC_INTERN(simuctrl->proctab[i].taskheap,  1, Queue);
          MALLOC_INTERN(simuctrl->proctab[i].taskheap2, 1, Queue);
          queueInit(simuctrl->proctab[i].taskheap,  100);
          queueInit(simuctrl->proctab[i].taskheap2, 100);

          MALLOC_INTERN(simuctrl->proctab[i].tasktab, 1, ExtendVectorINT);
          extendint_Init(simuctrl->proctab[i].tasktab, bloknbr/procnbr + 1);
        }

    /** Cluster initialization **/
    MALLOC_INTERN(simuctrl->clustab, clustnbr, SimuCluster);
    step = procnbr / clustnbr;
    for(i=0;i<clustnbr;i++)
      {
        simuctrl->clustab[i].fprocnum = i*step;
        simuctrl->clustab[i].lprocnum = simuctrl->clustab[i].fprocnum + step - 1;
        MALLOC_INTERN(simuctrl->clustab[i].ftgtsend, clustnbr, ExtendVectorINT);
        simuctrl->clustab[i].prionum  = 0;
        for(p=0;p<clustnbr;p++)
          extendint_Init(&(simuctrl->clustab[i].ftgtsend[p]), cblknbr/(2*clustnbr)+1);
      }
    simuctrl->clustab[clustnbr-1].lprocnum = procnbr-1;

    MALLOC_INTERN(simuctrl->ownetab, cblknbr, PASTIX_INT);
    MALLOC_INTERN(simuctrl->blprtab, bloknbr, PASTIX_INT);

    /* affect a negative value to cblk not mapped */
    for(i=0;i<cblknbr;i++)
      simuctrl->ownetab[i] = -1;
    for(i=0;i<bloknbr;i++)
      simuctrl->blprtab[i] = -1;


    MALLOC_INTERN(simuctrl->cblktab, cblknbr+1, SimuCblk);
    MALLOC_INTERN(simuctrl->bloktab, bloknbr+1, SimuBlok);
    ftgtcur = 0;

    for(i=0;i<cblknbr;i++)
      {
        candnbr = candtab[i].lccandnum - candtab[i].fccandnum + 1;
        simuctrl->cblktab[i].ctrbcnt = 0;

        for(j=symbptr->cblktab[i].bloknum;j<symbptr->cblktab[i+1].bloknum;j++)
          {
            simuctrl->bloktab[j].ftgtnum = ftgtcur;
            simuctrl->bloktab[j].tasknum = -1;
            simuctrl->bloktab[j].ctrbcnt = 0;
            /*if(candnbr > 1)*/
            ftgtcur += candnbr;
          }
      }
    /* one extracblk for avoiding side effect */
    simuctrl->bloktab[bloknbr].ftgtnum = ftgtcur;
    simuctrl->ftgtnbr = ftgtcur;

    if(simuctrl->ftgtnbr > 0)
        {
          /** Allocate and Initialize the timer for the reception of each ftgt on a candidate cluster **/
          MALLOC_INTERN(simuctrl->ftgttimetab, simuctrl->ftgtnbr, SimuTimer);
          for(i=0;i<simuctrl->ftgtnbr;i++)
            timerSet(&(simuctrl->ftgttimetab[i]), 0.0);

          MALLOC_INTERN(simuctrl->ftgttab, ftgtcur, SimuFtgt);
          for(i=0;i<simuctrl->ftgtnbr;i++)
            {
              simuctrl->ftgttab[i].clustnum = -1;
              timerSet(&(simuctrl->ftgttab[i].timerecv), 0.0);
              simuctrl->ftgttab[i].costsend = 0.0;
              simuctrl->ftgttab[i].costadd  = 0.0;
              bzero(simuctrl->ftgttab[i].ftgt.infotab,MAXINFO*sizeof(PASTIX_INT));
              simuctrl->ftgttab[i].ftgt.infotab[FTGT_FCOLNUM] = INTVALMAX;
              simuctrl->ftgttab[i].ftgt.infotab[FTGT_FROWNUM] = INTVALMAX;
              simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBNBR] = 0;
              simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBCNT] = 0;
            }
        }
    else
        {
            simuctrl->ftgttab     = NULL;
            simuctrl->tasktimetab = NULL;
            simuctrl->ftgttimetab = NULL;
        }

    return 1;
}


PASTIX_INT simuRealloc(SimuCtrl *simuctrl, PASTIX_INT procnbr, PASTIX_INT thrdlocnbr)
{
    PASTIX_INT i;

    /* Free processor structure */
    for(i=0;i<procnbr;i++)
      {
        queueExit     (simuctrl->proctab[i].taskheap);
        memFree_null  (simuctrl->proctab[i].taskheap);
        queueExit     (simuctrl->proctab[i].taskheap2);
        memFree_null  (simuctrl->proctab[i].taskheap2);
        extendint_Exit(simuctrl->proctab[i].tasktab);
        memFree_null  (simuctrl->proctab[i].tasktab);
      }
    memFree_null(simuctrl->proctab);

    /** Initialisation for local thread **/
    MALLOC_INTERN(simuctrl->proctab, thrdlocnbr, SimuProc);
    for(i=0;i<thrdlocnbr;i++)
        {
          MALLOC_INTERN(simuctrl->proctab[i].tasktab, 1, ExtendVectorINT);
          /* On initialise pas les vecteur d'entier, car il nous faudrait la structure d'arbre de bulles */
        }

    return 1;
}

void simuExit(SimuCtrl *simuctrl, PASTIX_INT clustnbr, PASTIX_INT procnbr, PASTIX_INT thrdlocnbr)
{
    PASTIX_INT i,j;
    (void)thrdlocnbr; (void)procnbr;

#ifndef PASTIX_DYNSCHED
    for(i=0;i<procnbr;i++)
      {
        queueExit(simuctrl->proctab[i].taskheap);
        memFree_null(simuctrl->proctab[i].taskheap);
        queueExit(simuctrl->proctab[i].taskheap2);
        memFree_null(simuctrl->proctab[i].taskheap2);
        extendint_Exit(simuctrl->proctab[i].tasktab);
        memFree_null(simuctrl->proctab[i].tasktab);
      }
#else
    for(i=0;i<thrdlocnbr;i++)
      {
        extendint_Exit(simuctrl->proctab[i].tasktab);
        memFree_null(simuctrl->proctab[i].tasktab);
      }
#endif

    for(i=0;i<clustnbr;i++)
      {
        for(j=0;j<clustnbr;j++)
          extendint_Exit(&(simuctrl->clustab[i].ftgtsend[j]));
        memFree_null(simuctrl->clustab[i].ftgtsend);
      }

    if(simuctrl->ftgttab != NULL)
        {
            memFree_null(simuctrl->ftgttab);
            memFree_null(simuctrl->ftgttimetab);
        }
    /* memFree_null(simuctrl->tasktimetab); */
    memFree_null(simuctrl->tasktab);
    memFree_null(simuctrl->proctab);
    memFree_null(simuctrl->clustab);
    memFree_null(simuctrl->ownetab);
    memFree_null(simuctrl->blprtab);
    memFree_null(simuctrl->cblktab);
    memFree_null(simuctrl->bloktab);
    memFree_null(simuctrl);
}


PASTIX_INT compTimer(SimuTimer *t1, SimuTimer *t2)
{
    /** Return 1 if t1 < t2 **/
    /** 0 in other cases **/
    if(t1->s < t2->s)
        return 1;
    /*if(t1->s == t2->s && t1->ms < t2->ms)
      return 1;*/
    return 0;
}

void timerAdd(SimuTimer *timer, double t)
{
    timer->s += t;

    /*timer->s += floor(t);
    timer->ms += (t-floor(t))*TIMEBASE;
    if(timer->ms >= TIMEBASE)
        {
            timer->s += 1;
            timer->ms = timer->ms - TIMEBASE;
            }*/
}

double timerVal(SimuTimer *t)
{
/*#ifdef DEBUG_BLEND
    ASSERT(t->ms < TIMEBASE,MOD_BLEND);
#endif*/
    return (t->s /* + t->ms/TIMEBASE */);
}

void timerSet(SimuTimer *timer, double t)
{
    timer->s = t;
    /*timer->s = floor(t);
      timer->ms = (t-floor(t))*TIMEBASE;*/
}


/** OIMBE pour optimisation faire un timerMax pour distrib **/
/** OIMBE pour optimisation faire un timerAffect(timer) pour distrib **/
