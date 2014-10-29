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

#include "common_pastix.h"
#include "symbol.h"
#include "ftgt.h"
#include "cand.h"
#include "queue.h"
#include "extendVector.h"
#include "simu.h"
#include "dof.h"
#include "elimin.h"
#include "cost.h"
#include "bulles.h"
#include "param_blend.h"
#include "blendctrl.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "costfunc.h"
/* #include "extrastruct.h" */
/* #include "param_comm.h" */
/* #include "perf.h" */
/* #include "assert.h" */
#include "task.h"

void taskBuild(SimuCtrl *simuctrl, SymbolMatrix *symbptr, Cand *candtab,
               const Dof * dofptr, EliminGraph *egraph, BlendCtrl *ctrl)
{
  PASTIX_INT i, j, k;
  PASTIX_INT tasknbr = 0;
  PASTIX_INT odb_nbr;
  PASTIX_INT L, h, g;
  PASTIX_INT firstE1task, firstE2task;
  (void)egraph;

  SimuTask *task = NULL;
  /** Count number of task **/
  for(i=0;i<symbptr->cblknbr;i++)
    {
      switch(candtab[i].distrib)
        {
          case D1 :
            /* Task 1D */
            tasknbr++;
            break;
          case D2 :
            /* number of extra diagonal blocks in the cblk */
            odb_nbr = symbptr->cblktab[i+1].bloknum - symbptr->cblktab[i].bloknum-1;
            tasknbr += 1 + odb_nbr + (odb_nbr*(odb_nbr+1))/2;
            break;
          case DENSE :
            fprintf(stderr, "DENSE cblk \n");
            tasknbr++;
            goto tnext;
          default:
            fprintf(stderr, "Task No %ld has wrong type \n", (long)i);
            EXIT(MOD_BLEND,INTERNAL_ERR);
        }
    }
 tnext:
  simuctrl->tasknbr = tasknbr;
  /*fprintf(stderr, "TASKNBR %ld \n" , tasknbr);
    EXIT(MOD_BLEND,INTERNAL_ERR);*/
  MALLOC_INTERN(simuctrl->tasktab, tasknbr, SimuTask);
#ifdef DEBUG_BLEND
  ASSERT(simuctrl->tasktab != NULL,MOD_BLEND);
#endif
  task = simuctrl->tasktab;
  tasknbr = 0;
  for(i=0;i<symbptr->cblknbr;i++)
    {
      switch(candtab[i].distrib)
        {
        case D1:
          task->taskid   = COMP_1D;
          task->prionum  = -1;
          task->cblknum  = i;
          task->bloknum  = symbptr->cblktab[i].bloknum;
          task->bloknum2 = -1;
          task->ctrbcnt  = 0;
          task->ftgtcnt  = 0;
          task->facebloknum = -1;
          task->cost     = -1;
          timerSet(&(task->time), 0.0);
          task->mesglen  = 0.0;
          task->tasknext = -1;
          for(j=symbptr->cblktab[i].bloknum;j<symbptr->cblktab[i+1].bloknum;j++)
            simuctrl->bloktab[j].tasknum = tasknbr;
          tasknbr++;
          task++;
          break;


        case D2:
          L  = symbptr->cblktab[i].lcolnum - symbptr->cblktab[i].fcolnum + 1;
          L *= (dofptr)->noddval;

          /** Factorize diagonal block task**/
          task->taskid   = DIAG;
          task->prionum  = -1;
          task->cblknum  = i;
          task->bloknum  = symbptr->cblktab[i].bloknum;
          task->bloknum2 = -1;
          task->ftgtcnt  = 0;
          task->ctrbcnt  = 0;
          task->facebloknum = -1;
          /* set task indice in bloktab */
          simuctrl->bloktab[task->bloknum].tasknum = tasknbr;
          task->cost = (double)DIAGCost(L);
          /*fprintf(stderr, "DIAGCost %f dost %f cost %f \n", DIAGCost(L), dost, task->cost);*/
          /*fprintf(stderr, "cost %f copy %g PPF %g \n", dost, PERF_COPY(L),  PERF_PPF(L) );*/

          timerSet(&(task->time), 0.0);
          /* OIMBE il faut prendre en compte toute la taille du message (pareil pour les fan in targets */
          /*task->timesend =TIME_STARTUP + TIME_BANDWIDTH * L*L*sizeof(double);*/
          task->mesglen = L*L*sizeof(double);
          task->tasknext = -1;
          tasknbr++;
          task++;
          firstE1task = tasknbr;
          for(j=symbptr->cblktab[i].bloknum+1;j<symbptr->cblktab[i+1].bloknum;j++)
            {
              h  = symbptr->bloktab[j].lrownum - symbptr->bloktab[j].frownum + 1;
              h *= (dofptr)->noddval;
              /** Solxve odb task **/
              task->taskid   = E1;
              task->prionum  = -1;
              task->cblknum  = i;
              task->bloknum  = j;
              task->bloknum2 = symbptr->cblktab[i].bloknum;
              task->ftgtcnt  = 0;
              task->ctrbcnt  = 0;
              task->facebloknum = -1;
              if(j<symbptr->cblktab[i+1].bloknum - 1)
                task->tasknext = tasknbr + symbptr->cblktab[i+1].bloknum-j+1;
              else
                {
                  /*task->tasknext = -1;*/
                  /** Now we make a cycle **/
                  task->tasknext = firstE1task;
                }

              /* set task indice in bloktab */
              simuctrl->bloktab[j].tasknum = tasknbr;
              task->cost = E1Cost(L, h);
              timerSet(&(task->time), 0.0);
              /*task->timesend = TIME_STARTUP + TIME_BANDWIDTH * L*h * sizeof(double);*/
              task->mesglen = L*h * sizeof(double);
              tasknbr++;
              task++;
              firstE2task = tasknbr;
              {
                PASTIX_INT lastbloknum;
                lastbloknum = 0;
                for(k=j;k<symbptr->cblktab[i+1].bloknum;k++)
                  {
                    g = symbptr->bloktab[k].lrownum - symbptr->bloktab[k].frownum + 1;
                    g *= (dofptr)->noddval;
                    /** Compute contrib task **/
                    task->taskid = E2;
                    task->prionum = -1;
                    task->cblknum = i;
                    task->bloknum = k;
                    task->ftgtcnt = 0;
                    task->ctrbcnt = 0;
                    task->bloknum2 = j;
                    task->facebloknum = getFaceBlockE2(lastbloknum, j, k, symbptr, ctrl->option->ricar);
                    lastbloknum = task->facebloknum;
                    task->cost = E2Cost(L, h, g);
                    timerSet(&(task->time), 0.0);
                    task->mesglen =  0.0;
                    if(k<symbptr->cblktab[i+1].bloknum - 1)
                      task->tasknext = tasknbr+1;
                    else
                      {
                        /*task->tasknext = -1;*/
                        /** We chain the E2 in a cycle now **/
                        task->tasknext = firstE2task;
                      }
                    tasknbr++;
                    task++;
                  }
              }
            }
          break;

        case DENSE:
          task->taskid = DRUNK;
          task->prionum = -1;
          task->cblknum = i;
          task->bloknum = symbptr->cblktab[i].bloknum;
          task->bloknum2 = (symbptr->nodenbr - symbptr->cblktab[i].fcolnum)*dofptr->noddval;
          task->ctrbcnt = 0;
          task->facebloknum = -1;
          task->cost = -1;
          timerSet(&(task->time), 0.0);
          task->mesglen = 0.0;
          task->tasknext = -1;
          for(j=symbptr->cblktab[i].bloknum;j<symbptr->bloknbr;j++)
            simuctrl->bloktab[j].tasknum = tasknbr;
          tasknbr++;
          task++;
          goto tend;
        default:
          fprintf(stderr, "Task No %ld has wrong type \n", (long)i);
          EXIT(MOD_BLEND,INTERNAL_ERR);
        }
    }
 tend:
#ifdef DEBUG_BLEND
  ASSERT(simuctrl->tasknbr == tasknbr,MOD_BLEND);
  for(i=0;i<simuctrl->tasknbr;i++)
    if(simuctrl->tasktab[i].tasknext != -1)
      if(simuctrl->tasktab[simuctrl->tasktab[i].tasknext].taskid != -1)
        ASSERT(simuctrl->tasktab[simuctrl->tasktab[i].tasknext].taskid == simuctrl->tasktab[i].taskid,MOD_BLEND);
#endif
;
}

/** Get face block for task E2 **/
PASTIX_INT getFaceBlockE2(PASTIX_INT startsearch, PASTIX_INT bloksrc, PASTIX_INT bloknum, const SymbolMatrix *symbptr, int ricar)
{
  PASTIX_INT i;

  if(startsearch < symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum].bloknum)
    startsearch = symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum].bloknum;

  ASSERT(startsearch < symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum+1].bloknum,MOD_BLEND);


  /*#ifndef NAPA*/
  if(ricar == 0)
    {
      /*for(i=symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum].bloknum;i<symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum+1].bloknum;i++)*/
      for(i=startsearch;i<symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum+1].bloknum;i++)
        if(symbptr->bloktab[i].lrownum >= symbptr->bloktab[bloknum].frownum)
          break;

#ifdef DEBUG_BLEND
      ASSERT(symbptr->bloktab[i].frownum <= symbptr->bloktab[bloknum].frownum
             && symbptr->bloktab[i].lrownum >= symbptr->bloktab[bloknum].lrownum,MOD_BLEND);
#endif

      return i;
    }
    /*#else*/


  /*for(i=symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum].bloknum;i<symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum+1].bloknum;i++)*/
  for(i=startsearch;i<symbptr->cblktab[symbptr->bloktab[bloksrc].cblknum+1].bloknum;i++)
    {
      if( (symbptr->bloktab[bloknum].frownum >= symbptr->bloktab[i].frownum &&  symbptr->bloktab[bloknum].frownum <= symbptr->bloktab[i].lrownum)
         || (symbptr->bloktab[bloknum].lrownum >= symbptr->bloktab[i].frownum &&  symbptr->bloktab[bloknum].lrownum <= symbptr->bloktab[i].lrownum)
         || (symbptr->bloktab[bloknum].frownum <= symbptr->bloktab[i].frownum &&  symbptr->bloktab[bloknum].lrownum >= symbptr->bloktab[i].lrownum)
         )
        return i;  /** We found the first block that matches **/
      if(symbptr->bloktab[bloknum].lrownum < symbptr->bloktab[i].frownum)
        {
          return -1;
        }

    }
  return -1;
  /*#endif*/



}

double taskSendCost(SimuTask *taskptr, const PASTIX_INT clustsrc, const PASTIX_INT clustdst, BlendCtrl *ctrl)
{
#ifdef DEBUG_BLEND
  ASSERT(clustsrc>=0 && clustsrc < ctrl->clustnbr,MOD_BLEND);
  ASSERT(clustdst>=0 && clustdst < ctrl->clustnbr,MOD_BLEND);

#endif
  if(clustsrc == clustdst)
    return 0.0;

  perfcluster2(clustsrc, clustdst, ctrl->candtab[taskptr->cblknum].lccandnum-ctrl->candtab[taskptr->cblknum].fccandnum+1,ctrl->perfptr, ctrl);

#ifdef DEBUG_BLEND
  ASSERT(taskptr->taskid != COMP_1D,MOD_BLEND);
#endif
  return (ctrl->perfptr->startup + ctrl->perfptr->bandwidth * taskptr->mesglen);
}
