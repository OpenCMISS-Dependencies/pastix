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
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "common_pastix.h"
#include "symbol.h"
#include "extendVector.h"
#include "queue.h"
#include "ftgt.h"
#include "cand.h"
#include "simu.h"
#include "param_blend.h"
#include "elimin.h"
#include "cost.h"
#include "bulles.h"
#include "blendctrl.h"
#include "dof.h"
#include "task.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "costfunc.h"
#include "trace.h"
#include "distribPart.h"

/* #include "param_comm.h" */
/* #include "extrastruct.h" */
/* #include "assembly.h" */
/* #include "perf.h" */

/*#define TC         TIME_STARTUP + TIME_BANDWIDTH*4096
#define ALPHA(ts)  (ts>=TC?1.0:0.4)*/
#define FACE

/** OIMBE a t on  tjs besoins de egraph->ownetab ???? **/

/******** VERIFIER clustdst et clustsrc dans costFtgtSend et taskSendCost *******/
/******** VERIFIER clustdst et clustsrc dans costFtgtSend et taskSendCost *******/
/******** VERIFIER clustdst et clustsrc dans costFtgtSend et taskSendCost *******/
/******** VERIFIER clustdst et clustsrc dans costFtgtSend et taskSendCost *******/
/******** VERIFIER clustdst et clustsrc dans costFtgtSend et taskSendCost *******/

void distribPart(SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl, const Dof * dofptr)
{

    PASTIX_INT             i, j, b;
    PASTIX_INT             cblknum, bloknum;
    /*PASTIX_INT             c;*/
    PASTIX_INT             pr;
    Queue           q;

    /** this queue is used to sort ftgt in receive order in cblkTimeFtgtAdd **/
    queueInit(&q, ctrl->procnbr);

    for(i=0;i<symbptr->cblknbr;i++)
      simuctrl->cblktab[i].ctrbcnt = ctrl->egraph->verttab[i].innbr;

    /* OIMBE attention les ctrbcnt des cblk en COMP1D sont recalculee dans computeBlockCtrbNbr */
    /** Compute number of contributions for blocks **/
    computeBlockCtrbNbr(simuctrl, symbptr, ctrl);

    if(ctrl->option->tracegen)
      /** Init tracefile **/
      for(pr=0;pr < ctrl->procnbr; pr++)
        trace_start(ctrl->tracefile, timerVal(TIMER(pr)),
                    pr/(ctrl->procnbr/ctrl->clustnbr), /* Numero de proc */
                    pr%(ctrl->procnbr/ctrl->clustnbr) /* Numero de thread */
                    );

    /*----------------------------------------------------------------------------------------*/
    /** Put the task with no contribution in the ready task heaps of their candidat processor */
    /*----------------------------------------------------------------------------------------*/
    for(i=0;i<symbptr->cblknbr;i++)
      {
        if(simuctrl->cblktab[i].ctrbcnt == 0) /* if cblk is a leave */
          {
            ASSERTDBG(ctrl->candtab[i].treelevel < 0,MOD_BLEND);
            if(ctrl->option->costlevel)
              ASSERTDBG(ctrl->candtab[i].costlevel < 0,MOD_BLEND);
            ASSERTDBG(simuctrl->tasktab[simuctrl->bloktab[symbptr->cblktab[i].bloknum].tasknum].taskid == COMP_1D
                   || simuctrl->tasktab[simuctrl->bloktab[symbptr->cblktab[i].bloknum].tasknum].taskid == DIAG,MOD_BLEND);
            if(simuctrl->tasktab[simuctrl->bloktab[symbptr->cblktab[i].bloknum].tasknum].taskid == COMP_1D)
              ASSERTDBG(ctrl->candtab[i].distrib == D1,MOD_BLEND);

            for(pr=ctrl->candtab[i].fcandnum;pr<=ctrl->candtab[i].lcandnum;pr++)
              putInReadyQueue(pr, simuctrl->bloktab[symbptr->cblktab[i].bloknum].tasknum,
                              symbptr, simuctrl, ctrl);
          }
      }

    /*--------------------------------*/
    /** Simules and mapps the tasks  **/
    /*--------------------------------*/
     while(1)
        {
          /** Get the next earlier task index and the processor on which it is mapped **/
          i = getNextTaskNextProc(simuctrl, ctrl, &pr);

          if(i==-1)
            break;

          bloknum = simuctrl->tasktab[i].bloknum;
          cblknum = simuctrl->tasktab[i].cblknum;

          /** Compute the time at which each proc cand
            will have added its ftgt and received block target if the task is mapped on **/
          if(simuctrl->tasktab[i].taskid != E2)
            {

              if(simuctrl->tasktab[i].taskid != COMP_1D)
                ASSERTDBG(simuctrl->blprtab[bloknum]<0,MOD_BLEND)
              else
                ASSERTDBG(simuctrl->ownetab[cblknum]<0,MOD_BLEND);

              /***** Set processor owners *****/
              switch(simuctrl->tasktab[i].taskid)
                {
                case DIAG:
                case E1:
                  simuctrl->blprtab[bloknum] = pr;
                  break;
                case COMP_1D:
                  simuctrl->ownetab[cblknum] = pr;
                  for(j=symbptr->cblktab[cblknum].bloknum;
                      j<symbptr->cblktab[cblknum+1].bloknum;j++)
                    {
                      simuctrl->blprtab[j] = pr;
                    }
                  break;

                default:
                  fprintf(stderr, "Task has wrong type \n");
                  EXIT(MOD_BLEND,INTERNAL_ERR);
                }
            }
          else
            ASSERTDBG(simuctrl->blprtab[bloknum] >= 0,MOD_BLEND);

          simuctrl->tasktab[i].prionum = simuctrl->clustab[ctrl->proc2clust[pr]].prionum;
          simuctrl->clustab[ctrl->proc2clust[pr]].prionum++;

          /* Ajout de la tache a la file du proc pour version standard */
          extendint_Add(simuctrl->proctab[pr].tasktab, i);
#if defined(TRACE_SOPALIN) || defined(PASTIX_DYNSCHED)
          ctrl->candtab[simuctrl->tasktab[i].cblknum].cand = pr;
#endif

          /* Sauvegarde du processus MPI devant executer la tache pour version MARCEL */
          if ((simuctrl->tasktab[i].taskid == DIAG)
              || (simuctrl->tasktab[i].taskid == COMP_1D))
            {
              ASSERTDBG(simuctrl->tasktab[i].cblknum < symbptr->cblknbr, MOD_BLEND);
              ctrl->candtab[simuctrl->tasktab[i].cblknum].cluster = ctrl->proc2clust[pr];
            }

          /*-------------------------------------------------------------/
          /   UPDATE TIMER OF THE PROC ON WHICH IS MAPPED THE TASK       /
          /   TO THE DATE THE PROC WILL BEGIN THE TASK INNER COMPUTATION /
          /-------------------------------------------------------------*/
          if(simuctrl->tasktab[i].taskid == E2
             || ctrl->candtab[simuctrl->tasktab[i].cblknum].fccandnum == ctrl->candtab[simuctrl->tasktab[i].cblknum].lccandnum)
            /** Time do not depend on the reception of a ftgt **/
            timerSet(TIMER(pr), MAX(timerVal(TIMER(pr)), timerVal(&(simuctrl->tasktab[i].time))));
          else
            /** Time depends on the reception of a ftgt **/
            timerSet(TIMER(pr),
                     MAX(timerVal(TIMER(pr)),
                         timerVal(&(simuctrl->ftgttimetab[CLUST2INDEX(bloknum, ctrl->proc2clust[pr])]))));

          /*------------------------------------------------------------------------/
          /  Fill some fanintarget info (task of type E2 does not have any ftgt)    /
          /------------------------------------------------------------------------*/
          if(simuctrl->tasktab[i].taskid != E2)
            {
              switch(ctrl->candtab[cblknum].distrib)
                {
                case D2:

                  /** Task DIAG or E1 **/
                  if(simuctrl->bloktab[bloknum].ftgtnum< simuctrl->bloktab[bloknum+1].ftgtnum)
                    {
                      for(j=simuctrl->bloktab[bloknum].ftgtnum;j<simuctrl->bloktab[bloknum+1].ftgtnum;j++)
                        {
                          if((simuctrl->ftgttab[j].ftgt.infotab[FTGT_CTRBNBR] >0)
                             && (j != CLUST2INDEX(bloknum, ctrl->proc2clust[pr])))
                            {
                              /*simuctrl->ftgttab[j].procnum = PROC(j,bloknum);*/
                              simuctrl->ftgttab[j].clustnum = INDEX2CLUST(j, bloknum);
                              /*simuctrl->ftgttab[j].ftgt.infotab[FTGT_PRIONUM] = simuctrl->ftgtprio;*/
                              simuctrl->ftgttab[j].ftgt.infotab[FTGT_PRIONUM] = simuctrl->tasktab[i].prionum;
                              simuctrl->ftgttab[j].ftgt.infotab[FTGT_PROCDST] = pr;
                              simuctrl->ftgttab[j].ftgt.infotab[FTGT_BLOKDST] = bloknum;
                              simuctrl->ftgttab[j].ftgt.infotab[FTGT_TASKDST] = simuctrl->bloktab[bloknum].tasknum;
#ifdef OOC_FTGT
                              simuctrl->ftgttab[j].ftgt.infotab[FTGT_GCBKDST] = simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].cblknum;
#endif
                              extendint_Add(&(simuctrl->clustab[INDEX2CLUST(j,bloknum)].ftgtsend[ctrl->proc2clust[pr]]), j);
                              /*simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].ftgtcnt += simuctrl->ftgttab[j].ftgt.infotab[FTGT_CTRBNBR];*/

                              simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].ftgtcnt++;

                              if (ctrl->proc2clust[pr] == ctrl->clustnum)
                                simuctrl->ftgtcnt++;
                            }
                        }

                      simuctrl->ftgtprio++;
                    }
                  else
                    ASSERTDBG(ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum,MOD_BLEND);

                  break;

                case D1:
                  if(simuctrl->bloktab[bloknum].ftgtnum< simuctrl->bloktab[bloknum+1].ftgtnum)
                    {
                      /** Task COMP_1D with several cand cluster **/
                      for(b=bloknum;b<symbptr->cblktab[cblknum+1].bloknum;b++)
                        {
                          for(j=simuctrl->bloktab[b].ftgtnum;j<simuctrl->bloktab[b+1].ftgtnum;j++)
                            {
                              if((simuctrl->ftgttab[j].ftgt.infotab[FTGT_CTRBNBR] >0)
                                 && (j != CLUST2INDEX(b, ctrl->proc2clust[pr])))
                                {
                                  /*simuctrl->ftgttab[j].procnum = PROC(j,b);*/
                                  simuctrl->ftgttab[j].clustnum = INDEX2CLUST(j, b);
                                  /*simuctrl->ftgttab[j].ftgt.infotab[FTGT_PRIONUM] = simuctrl->ftgtprio;*/
                                  simuctrl->ftgttab[j].ftgt.infotab[FTGT_PRIONUM] = simuctrl->tasktab[i].prionum;
                                  simuctrl->ftgttab[j].ftgt.infotab[FTGT_PROCDST] = pr;
                                  simuctrl->ftgttab[j].ftgt.infotab[FTGT_BLOKDST] = b;
                                  simuctrl->ftgttab[j].ftgt.infotab[FTGT_TASKDST] = simuctrl->bloktab[bloknum].tasknum;
#ifdef OOC_FTGT
                                  simuctrl->ftgttab[j].ftgt.infotab[FTGT_GCBKDST] = simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].cblknum;
#endif
                                  extendint_Add(&(simuctrl->clustab[INDEX2CLUST(j,b)].ftgtsend[ctrl->proc2clust[pr]]), j);

                                  /*simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].ftgtcnt += simuctrl->ftgttab[j].ftgt.infotab[FTGT_CTRBNBR];*/
                                  simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].ftgtcnt++;

                                  if (ctrl->proc2clust[pr] == ctrl->clustnum)
                                    simuctrl->ftgtcnt++;
                                }
                            }

                        }
                      simuctrl->ftgtprio++;

                    }
                  else
                    ASSERTDBG(ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum,MOD_BLEND);

                  break;
                default:
                  EXIT(MOD_BLEND,INTERNAL_ERR);
                }
            }

          /* Simule the computing of the task */
          switch(simuctrl->tasktab[i].taskid)
            {

            case COMP_1D:
              taskExec_COMP1D(i, symbptr, simuctrl, ctrl, dofptr);
              break;
            case DIAG:
              /*fprintf(stderr, "SMP Version not yet in 2D \n");
              EXIT(MOD_BLEND,NOTIMPLEMENTED_ERR);*/
              taskExec_DIAG(i, symbptr, simuctrl, ctrl, dofptr);
              break;
            case E1:
              /*fprintf(stderr, "SMP Version not yet in 2D \n");
              EXIT(MOD_BLEND,NOTIMPLEMENTED_ERR);*/
              taskExec_E1(i, symbptr, simuctrl, ctrl, dofptr);
              break;
            case E2:
              /*fprintf(stderr, "SMP Version not yet in 2D \n");
              EXIT(MOD_BLEND,NOTIMPLEMENTED_ERR);*/
              taskExec_E2(i, symbptr, simuctrl, ctrl, dofptr);
              break;
            default:
              errorPrint("task has no type.");
              EXIT(MOD_BLEND,INTERNAL_ERR);
            }

          /*fprintf(stderr, "TIMER %f \n", timerVal(TIMER(pr)));*/
          queueReorder(i, symbptr, simuctrl, ctrl);

        }



    if(ctrl->option->tracegen)
      {
        for(pr=0; pr<ctrl->procnbr; pr++)
          trace_finish(ctrl->tracefile, timerVal(TIMER(pr)),
                       pr/(ctrl->procnbr/ctrl->clustnbr), /* Numero de proc */
                       pr%(ctrl->procnbr/ctrl->clustnbr) /* Numero de thread */
                       );
      }


    double maxtime = 0;
    for(pr=0;pr<ctrl->procnbr;pr++)
      {
        if(timerVal(TIMER(pr)) > maxtime)
          maxtime = timerVal(TIMER(pr));
      }
    set_dparm(ctrl->option->dparm, DPARM_PRED_FACT_TIME, maxtime);

    if((ctrl->clustnum == 0) &&
       (ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO))
      {
        FILE *out;
        OUT_OPENFILEINDIR(ctrl->option->iparm, out, "taskrepartition", "w");

        /* fprintf(out," \n TIME PREDICTED FOR FACTORIZATION: \n");  */
        for(pr=0;pr<ctrl->procnbr;pr++)
          {
            /* fprintf(stdout," Proc %ld : %g / %ld tasks\n", (long)pr, timerVal(TIMER(pr)),  */
            /*              (long)extendint_Size(simuctrl->proctab[pr].tasktab)); */
            fprintf(out, "%ld %g %ld\n", (long)pr, timerVal(TIMER(pr)),
                    (long)extendint_Size(simuctrl->proctab[pr].tasktab));
          }
        /* fprintf(stdout," \n *** SOLVER RESOLUTION PREDICTED IN %g ***\n",maxtime);  */

        OUT_CLOSEFILEINDIR(out);
      }

#ifdef DEBUG_BLEND

    for(i=0;i<simuctrl->cblknbr;i++)
      if(ctrl->candtab[i].distrib == D1)
        if(simuctrl->ownetab[i] < 0) /** Check valid for 1D distribution only **/
          fprintf(stderr, "CBLK %ld has no processor \n", (long)i);

    for(i=0;i<symbptr->bloknbr;i++)
      if(!(simuctrl->blprtab[i]>=0))
        {
          fprintf(stderr, "BLOCK %ld has no processor \n", (long)i);
          fprintf(stdout, "blprtab [ %ld ] = %ld type %ld \n", (long)i,
                  (long)simuctrl->blprtab[i],
                  (long)simuctrl->tasktab[simuctrl->bloktab[i].tasknum].taskid);
          EXIT(MOD_BLEND,INTERNAL_ERR);
        }
    for(i=0;i<symbptr->bloknbr;i++)
      ASSERTDBG(simuctrl->blprtab[i]>=0,MOD_BLEND);
#endif
    /** Free allocated memory **/
    queueExit(&q);

}

void taskExec_DIAG(PASTIX_INT tasknum, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl, const Dof * dofptr)
{
  PASTIX_INT i;
  PASTIX_INT pr;
  PASTIX_INT cblknum, bloknum;
  PASTIX_INT procnum;
  bloknum = simuctrl->tasktab[tasknum].bloknum;
  cblknum = simuctrl->tasktab[tasknum].cblknum;
  procnum = simuctrl->blprtab[bloknum];

  if(ctrl->option->tracegen)
    {
      trace_begin_task2(ctrl->tracefile, timerVal(TIMER(procnum)),
                        procnum/(ctrl->procnbr/ctrl->clustnbr), /* Numero de proc */
                        procnum%(ctrl->procnbr/ctrl->clustnbr), /* Numero de thread */
                        STATE_DIAG,
                        tasknum,
                        ctrl->candtab[cblknum].fcandnum,
                        ctrl->candtab[cblknum].lcandnum,
                        procnum
                       );
    }

  /** Add the time to compute the task to the proc timer **/
  timerAdd(TIMER(procnum), simuctrl->tasktab[tasknum].cost);


  /** Create tasks E1 whose all contributions are computed **/
  for(i=symbptr->cblktab[cblknum].bloknum+1; i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
      /** Set timer of task E1 to the date at which  the DIAG task is finished **/
      timerSet(&(simuctrl->tasktab[simuctrl->bloktab[i].tasknum].time),
               MAX(timerVal(&(simuctrl->tasktab[simuctrl->bloktab[i].tasknum].time)), timerVal(TIMER(procnum))));

      simuctrl->bloktab[i].ctrbcnt--;
      if(simuctrl->bloktab[i].ctrbcnt == 0)
        {

          if(ctrl->candtab[cblknum].fccandnum != ctrl->candtab[cblknum].lccandnum)
            {
              /* Compute time received task for this task E1 of the cblk */
              computeTaskReceiveTime(simuctrl->bloktab[i].tasknum, symbptr, simuctrl, ctrl, dofptr);
            }

          /* Put the task E1 in all cand processors queue */
          for(pr = ctrl->candtab[cblknum].fcandnum; pr<=ctrl->candtab[cblknum].lcandnum;pr++)
            {
              putInReadyQueue(pr, simuctrl->bloktab[i].tasknum, symbptr, simuctrl, ctrl);
            }

          ASSERTDBG(simuctrl->bloktab[i].tasknum < simuctrl->tasknbr,MOD_BLEND);
        }
    }

  if(ctrl->option->tracegen)
    {
      trace_end_task(ctrl->tracefile, timerVal(TIMER(procnum)),
                     procnum/(ctrl->procnbr/ctrl->clustnbr), /* Numero de proc */
                     procnum%(ctrl->procnbr/ctrl->clustnbr), /* Numero de thread */
                     1,
                     STATE_DIAG,
                     tasknum
                     );
    }

}

void taskExec_E1(PASTIX_INT tasknum, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl, const Dof * dofptr)
{
 PASTIX_INT i;
 PASTIX_INT cblknum, bloknum;
 PASTIX_INT procnum;
 (void)dofptr;

 bloknum = simuctrl->tasktab[tasknum].bloknum;
 cblknum = simuctrl->tasktab[tasknum].cblknum;
 procnum = simuctrl->blprtab[bloknum];

 if(ctrl->option->tracegen)
   {
      trace_begin_task2(ctrl->tracefile, timerVal(TIMER(procnum)),
                        procnum/(ctrl->procnbr/ctrl->clustnbr), /* Numero de proc */
                        procnum%(ctrl->procnbr/ctrl->clustnbr), /* Numero de thread */
                        STATE_E1,
                        tasknum,
                        ctrl->candtab[cblknum].fcandnum,
                        ctrl->candtab[cblknum].lcandnum,
                        procnum
                       );
   }

 /** Add the time to compute the task to the proc timer **/
 timerAdd(TIMER(procnum), simuctrl->tasktab[tasknum].cost);


 /*---------------------------------------------------------------------/
 / Update timer for task E2 on the same block than the current task     /
 / and put them in their ready task heap                                /
 /---------------------------------------------------------------------*/
 for(i=symbptr->cblktab[cblknum].bloknum+1;i<=bloknum;i++)
   {
     ASSERTDBG(simuctrl->tasktab[simuctrl->bloktab[i].tasknum + bloknum - i+1].taskid == E2,MOD_BLEND);
     ASSERTDBG(simuctrl->tasktab[simuctrl->bloktab[i].tasknum + bloknum - i+1].bloknum == bloknum,MOD_BLEND);
     ASSERTDBG(simuctrl->tasktab[simuctrl->bloktab[i].tasknum + bloknum - i+1].bloknum2 == i,MOD_BLEND);

     if(simuctrl->blprtab[i] < 0) /** Task E2 of a block not yet mapped **/
       {
         ASSERTDBG(timerVal(&(simuctrl->tasktab[simuctrl->bloktab[i].tasknum + bloknum - i+1].time)) == 0,MOD_BLEND);
         timerSet(&(simuctrl->tasktab[simuctrl->bloktab[i].tasknum + bloknum - i+1].time), timerVal(TIMER(procnum)));
       }
     else
       {
         timerSet(&(simuctrl->tasktab[simuctrl->bloktab[i].tasknum + bloknum - i+1].time),
                  MAX(timerVal(&(simuctrl->tasktab[simuctrl->bloktab[i].tasknum + bloknum - i+1].time))
                      + taskSendCost( &(simuctrl->tasktab[simuctrl->bloktab[i].tasknum+1]), ctrl->proc2clust[simuctrl->blprtab[i]], ctrl->proc2clust[procnum], ctrl ),
                      timerVal(TIMER(procnum))));
         /** add this task E2 to processor ready queue **/
         putInReadyQueue(procnum, simuctrl->bloktab[i].tasknum+bloknum -i+1, symbptr, simuctrl, ctrl);
       }
   }


 /*----------------------------------------------------------------/
 / tasks E2 for block below whose task E1 are processed            /
 / must have their TIMER UPDATED and put in the READY TASK HEAPs   /
 /----------------------------------------------------------------*/
 for(i=bloknum+1; i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
      if(simuctrl->blprtab[i] < 0)
        timerSet(&(simuctrl->tasktab[tasknum + i-bloknum+1].time), timerVal(TIMER(procnum)));
      else
        {
          timerSet(&(simuctrl->tasktab[tasknum + i-bloknum+1].time), MAX(timerVal(&(simuctrl->tasktab[tasknum + i-bloknum+1].time)),timerVal(TIMER(procnum))+ taskSendCost( &(simuctrl->tasktab[tasknum]), ctrl->proc2clust[procnum], ctrl->proc2clust[simuctrl->blprtab[i]], ctrl)) );

          putInReadyQueue(simuctrl->blprtab[i], tasknum + i-bloknum+1, symbptr, simuctrl, ctrl);
        }
      ASSERTDBG(tasknum + i-bloknum + 1< simuctrl->tasknbr,MOD_BLEND);
      ASSERTDBG(simuctrl->tasktab[tasknum + i-bloknum+1].taskid == E2,MOD_BLEND);
      ASSERTDBG(simuctrl->tasktab[tasknum + i-bloknum+1].bloknum == i,MOD_BLEND);
      ASSERTDBG(simuctrl->tasktab[tasknum + i-bloknum+1].bloknum2 == bloknum,MOD_BLEND);
    }

 if(ctrl->option->tracegen)
   {
     trace_end_task(ctrl->tracefile, timerVal(TIMER(procnum)),
                    procnum/(ctrl->procnbr/ctrl->clustnbr), /* Numero de proc */
                    procnum%(ctrl->procnbr/ctrl->clustnbr), /* Numero de thread */
                    1,
                    STATE_E1,
                    tasknum
                    );
   }
}

void taskExec_E2(PASTIX_INT tasknum, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl, const Dof * dofptr)
{
  PASTIX_INT bloknum;
  PASTIX_INT facebloknum;
  PASTIX_INT facetasknum;
  PASTIX_INT ftgtnum;
  PASTIX_INT procnum;
  /*PASTIX_INT clustnum;*/
  PASTIX_INT local;  /** if == 1 facing block is local **/
  PASTIX_INT pr;

  bloknum = simuctrl->tasktab[tasknum].bloknum;
  procnum = simuctrl->blprtab[bloknum];

#if defined(TRACE_SOPALIN)
  if(ctrl->option->tracegen)
    {
      PASTIX_INT cblknum = simuctrl->tasktab[tasknum].cblknum;
      trace_begin_task2(ctrl->tracefile, timerVal(TIMER(procnum)),
                        procnum/(ctrl->procnbr/ctrl->clustnbr), /* Numero de proc */
                        procnum%(ctrl->procnbr/ctrl->clustnbr), /* Numero de thread */
                        STATE_E2,
                        tasknum,
                        ctrl->candtab[cblknum].fcandnum,
                        ctrl->candtab[cblknum].lcandnum,
                        procnum
                       );
    }
#endif

  timerAdd(TIMER(simuctrl->blprtab[bloknum]), simuctrl->tasktab[tasknum].cost);

  /** get face block **/
  facebloknum = simuctrl->tasktab[tasknum].facebloknum;
  facetasknum = simuctrl->bloktab[facebloknum].tasknum;

  ASSERTDBG(simuctrl->tasktab[facetasknum].taskid != E2,MOD_BLEND);
  ASSERTDBG(simuctrl->tasktab[facetasknum].bloknum == facebloknum,MOD_BLEND);

  /** Is the facing task local or not **/
  if(ctrl->candtab[simuctrl->tasktab[facetasknum].cblknum].fccandnum
     ==  ctrl->candtab[simuctrl->tasktab[facetasknum].cblknum].lccandnum)
    local = 1;
  else
    local = 0;


  if(!local)
    {
      /************************************************/
      /** Update fan in target timer of facing block **/
      /************************************************/
      /*ftgtnum = INDEX(facebloknum, procnum); */
      ftgtnum = CLUST2INDEX(facebloknum, ctrl->proc2clust[procnum]);

      timerSet(&(simuctrl->ftgttimetab[ftgtnum]), MAX(timerVal(TIMER(procnum)), timerVal(&(simuctrl->ftgttimetab[ftgtnum]))));
      /** Update FanInTarget structure **/
      updateFtgtStruct(bloknum, simuctrl->tasktab[tasknum].bloknum2, ftgtnum, symbptr, simuctrl, ctrl);
    }
  else
    {
      timerSet(&(simuctrl->tasktab[facetasknum].time), MAX(timerVal(TIMER(procnum)), timerVal(&(simuctrl->tasktab[facetasknum].time))));
    }
  /**********************************************************************************/
  /** Update timer of facing block candidat  proc in the same cluster than procnum **/
  /**********************************************************************************/
  /** Update contribution counter of the facing block **/
  simuctrl->bloktab[facebloknum].ctrbcnt--;

  if(simuctrl->bloktab[facebloknum].ctrbcnt == 0)
    {
      if(!local)
        computeTaskReceiveTime(facetasknum, symbptr, simuctrl, ctrl, dofptr);

      /* The task DIAG or E1 can be put in task heaps of candidat processor */
      for(pr=ctrl->candtab[simuctrl->tasktab[facetasknum].cblknum].fcandnum;
          pr<=ctrl->candtab[simuctrl->tasktab[facetasknum].cblknum].lcandnum; pr++)
        {
          putInReadyQueue(pr, facetasknum, symbptr, simuctrl, ctrl);
        }

      ASSERTDBG(facetasknum < simuctrl->tasknbr,MOD_BLEND);
    }

  if(ctrl->option->tracegen)
    {
      trace_end_task(ctrl->tracefile, timerVal(TIMER(procnum)),
                     procnum/(ctrl->procnbr/ctrl->clustnbr), /* Numero de proc */
                     procnum%(ctrl->procnbr/ctrl->clustnbr), /* Numero de thread */
                     1,
                     STATE_E2,
                     tasknum
                     );
    }
}

void taskExec_COMP1D(PASTIX_INT tasknum, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl, const Dof * dofptr)
{
    PASTIX_INT          i, j;
    PASTIX_INT          pr;
    PASTIX_INT          cblknum;
    PASTIX_INT          facebloknum;
    PASTIX_INT          facecblknum;
    PASTIX_INT          local;
    PASTIX_INT          facetasknum;
    PASTIX_INT          ftgtnum;
    PASTIX_INT          procnum;
    /*PASTIX_INT          clustnum;*/
    SimuProc     *proc;
    CostMatrix   *costmtx;

    cblknum = simuctrl->tasktab[tasknum].cblknum; /* in case of COMP1D bloknum in a SimuTask struct means cblknum */
    procnum = simuctrl->ownetab[cblknum];
    proc    = &(simuctrl->proctab[procnum]);
    costmtx = ctrl->costmtx;

#ifdef DEBUG_BLEND
    if (procnum < ctrl->candtab[cblknum].fcandnum || procnum > ctrl->candtab[cblknum].lcandnum)
      fprintf(stderr, "procnum : %ld, fcandnum : %ld, lcandnum : %ld\n",
              (long)procnum, (long)ctrl->candtab[cblknum].fcandnum, (long)ctrl->candtab[cblknum].lcandnum);
    ASSERT(procnum >= ctrl->candtab[cblknum].fcandnum && procnum <= ctrl->candtab[cblknum].lcandnum,MOD_BLEND);
#endif

    if(ctrl->option->tracegen)
      {
        trace_begin_task2(ctrl->tracefile, timerVal(TIMER(procnum)),
                          procnum/(ctrl->procnbr/ctrl->clustnbr), /* Numero de proc */
                          procnum%(ctrl->procnbr/ctrl->clustnbr), /* Numero de thread */
                          STATE_COMP1D,
                          tasknum,
                          ctrl->candtab[cblknum].fcandnum,
                          ctrl->candtab[cblknum].lcandnum,
                          procnum
                          );
      }

    /** Add time for factorizatoin of the diag blok and repercution on the off diag bloks **/
    timerAdd(&(proc->timer), costmtx->cblktab[cblknum].compute);

    for(i=symbptr->cblktab[cblknum].bloknum+1;i<symbptr->cblktab[cblknum+1].bloknum;i++)
        {
            /** Add time for compute of the contrib due to this odb **/
            /** OIMBE: pour l'instant je considere que les contribs sont calculees
                en un bloc **/
            timerAdd(&(proc->timer), costmtx->bloktab[i].contrib);

            facecblknum = symbptr->bloktab[i].cblknum;

            /*if(simuctrl->tasktab[simuctrl->bloktab[symbptr->cblktab[facecblknum].bloknum].tasknum].taskid == COMP_1D
               && (ctrl->candtab[facecblknum].fccandnum ==  ctrl->candtab[facecblknum].lccandnum) )*/
            if(ctrl->candtab[facecblknum].fccandnum ==  ctrl->candtab[facecblknum].lccandnum)
              local = 1; /** Facing task is COMP1D and is on the local cluster subtree **/
            else
              local = 0;

            if(!local || ctrl->candtab[facecblknum].distrib != D1 )
              {
                /*facebloknum = symbptr->cblktab[facecblknum].bloknum;*/
                facebloknum = 0;

                for(j=i;j<symbptr->cblktab[cblknum+1].bloknum;j++)
                  {
                    /* OIMBE trop couteux !! */
                    facebloknum = getFaceBlockE2(facebloknum, i, j, symbptr, ctrl->option->ricar);

#ifdef DEBUG_M
                    if(ctrl->option->ricar == 0)
                      ASSERT(facebloknum >= 0,MOD_BLEND);
#endif
                    /*#ifdef NAPA*/
                    if(facebloknum >= 0)
                      {
                        /*#endif*/
                        if(ctrl->candtab[facecblknum].distrib == D1)
                          simuctrl->cblktab[facecblknum].ctrbcnt--;
                        else
                          simuctrl->bloktab[facebloknum].ctrbcnt--;

                        if(!local)
                          {
                            ftgtnum = CLUST2INDEX(facebloknum, ctrl->proc2clust[procnum]);
                            updateFtgtStruct(j, i, ftgtnum, symbptr, simuctrl, ctrl);

                            /** Update timer ready for receiver of the ftgt **/

                            if(ctrl->candtab[facecblknum].distrib == D2)
                              {
                                /*timerSet(&(simuctrl->ftgttimetab[ftgtnum]), timerVal(&(proc->timer)));*/
                              }
                            else
                              {
                                ftgtnum = CLUST2INDEX(symbptr->cblktab[facecblknum].bloknum, ctrl->proc2clust[procnum]);
                                /*timerSet(&(simuctrl->ftgttimetab[ftgtnum]) , MAX( timerVal(&(simuctrl->ftgttimetab[ftgtnum])) ,timerVal(&(proc->timer))));*/
                              }
                            timerSet(&(simuctrl->ftgttimetab[ftgtnum]) , MAX( timerVal(&(simuctrl->ftgttimetab[ftgtnum])) ,
                                                                              timerVal(&(proc->timer))));
                          }
                        else
                          {
                            /*** LOCAL task ***/
                            if(ctrl->candtab[facecblknum].distrib == D2)
                              {
                                facetasknum = simuctrl->bloktab[facebloknum].tasknum;
                                timerSet(&(simuctrl->tasktab[facetasknum].time), MAX(timerVal(&(proc->timer)),
                                                                                     timerVal(&(simuctrl->tasktab[facetasknum].time))));
                              }
                            else
                              {
                                facetasknum = simuctrl->bloktab[symbptr->cblktab[facecblknum].bloknum].tasknum;
                                timerSet(&(simuctrl->tasktab[facetasknum].time), MAX(timerVal(&(proc->timer)),
                                                                                     timerVal(&(simuctrl->tasktab[facetasknum].time))));
                              }
                          }

                        if( (ctrl->candtab[facecblknum].distrib == D2 && simuctrl->bloktab[facebloknum].ctrbcnt == 0)
                            || (ctrl->candtab[facecblknum].distrib == D1 && simuctrl->cblktab[facecblknum].ctrbcnt == 0) )
                          {
                            facetasknum = simuctrl->bloktab[facebloknum].tasknum;
                            ASSERTDBG(simuctrl->tasktab[facetasknum].taskid != E2,MOD_BLEND);

                            if(!local)
                              computeTaskReceiveTime(facetasknum, symbptr, simuctrl, ctrl, dofptr);

                            for(pr=ctrl->candtab[facecblknum].fcandnum;
                                pr<=ctrl->candtab[facecblknum].lcandnum; pr++)
                              {
                                putInReadyQueue(pr, facetasknum, symbptr, simuctrl, ctrl);
                              }
                            ASSERTDBG(facetasknum < simuctrl->tasknbr,MOD_BLEND);
                          }
                        /*#ifdef NAPA*/
                      }
                    /*#endif*/
                  }
              }
            else
              {
                /** The facing task is local COMP_1D**/
                /*#ifdef NAPA*/
                if(ctrl->option->ricar == 1)
                  {
                    /*facebloknum = symbptr->cblktab[facecblknum].bloknum;*/
                    facebloknum = 0;
                    for(j=i;j<symbptr->cblktab[cblknum+1].bloknum;j++)
                      {
                        /* OIMBE trop couteux ON PEUT FAIRE MIEUX EN PARCOURANT EN DESCENDANT!! */
                        facebloknum = getFaceBlockE2(facebloknum, i, j, symbptr, ctrl->option->ricar);
                        if(facebloknum>=0)
                          simuctrl->cblktab[facecblknum].ctrbcnt--;
                      }
                    /*#else*/
                  }
                else
                  simuctrl->cblktab[facecblknum].ctrbcnt -= symbptr->cblktab[cblknum+1].bloknum - i;            /** A REVOIR POUR NAPA **/
                /*#endif*/

                ASSERTDBG(simuctrl->cblktab[facecblknum].ctrbcnt >= 0,MOD_BLEND);

                if(simuctrl->cblktab[facecblknum].ctrbcnt == 0)
                  {
                    facetasknum = simuctrl->bloktab[symbptr->cblktab[facecblknum].bloknum].tasknum;
                    /** Update timer of the task (owned by the diag block**/
                    ASSERTDBG(ctrl->candtab[facecblknum].fccandnum == ctrl->candtab[facecblknum].lccandnum,MOD_BLEND);

                    /* NB: The facing cblk is local !! */
                    /*timerSet(&(simuctrl->tasktab[facetasknum].time), timerVal(&(proc->timer)));*/
                    timerSet(&(simuctrl->tasktab[facetasknum].time), MAX(timerVal(&(proc->timer)),
                                                                         timerVal(&(simuctrl->tasktab[facetasknum].time))) );
#ifdef DEBUG_BLEND
                    ASSERTDBG(simuctrl->tasktab[facetasknum].taskid == COMP_1D,MOD_BLEND);
                    /*#ifdef NAPA*/
                    if(ctrl->option->ricar == 1)
                      ASSERTDBG(ctrl->candtab[facecblknum].lccandnum == ctrl->candtab[facecblknum].fccandnum,MOD_BLEND);
                    /*#endif*/
                    if(ctrl->proc2clust[procnum] != ctrl->candtab[facecblknum].fccandnum)
                      {
                        fprintf(stderr, "clustnum %ld  face proc cand %ld \n", (long)ctrl->proc2clust[procnum], (long)ctrl->candtab[facecblknum].fccandnum);
                        fprintf(stderr, "%ld candidat [%ld %ld] => %ld candidat [%ld %ld ]\n", (long)cblknum, (long)ctrl->candtab[cblknum].fccandnum, (long)ctrl->candtab[cblknum].lccandnum,
                                 (long)facecblknum, (long)ctrl->candtab[facecblknum].fccandnum, (long)ctrl->candtab[facecblknum].lccandnum);
                      }
                    ASSERTDBG(ctrl->proc2clust[procnum] == ctrl->candtab[facecblknum].fccandnum,MOD_BLEND);
                    ASSERTDBG(ctrl->proc2clust[procnum] == ctrl->candtab[facecblknum].lccandnum,MOD_BLEND);
                    ASSERTDBG(facetasknum<simuctrl->tasknbr,MOD_BLEND);
#endif

                    /** Put the task in the ready heap of its local candidat processor **/
                    for(pr = ctrl->candtab[facecblknum].fcandnum; pr <= ctrl->candtab[facecblknum].lcandnum;pr++)
                      putInReadyQueue(pr, facetasknum, symbptr, simuctrl, ctrl);
                  }
              }
          }

    if(ctrl->option->tracegen)
      {
        trace_end_task(ctrl->tracefile, timerVal(TIMER(procnum)),
                       procnum/(ctrl->procnbr/ctrl->clustnbr), /* Numero de proc */
                       procnum%(ctrl->procnbr/ctrl->clustnbr), /* Numero de thread */
                       1,
                       STATE_COMP1D,
                       tasknum
                       );
      }
}

void computeTaskReceiveTime(const PASTIX_INT tasknum, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl, const Dof * dofptr)
{
  /*-------------------------------------------------------------------------/
  / Compute the time the cblk would have RECEIVED and ADDED                  /
  / all its contributions if it was mapped on a given cand CLUSTER           /
  / !! These times not include add time for fan in target !!                 /
  /--------------------------------------------------------------------------*/

  PASTIX_INT i, j;
  double lftgttime = 0;
  double sftgttime = 0;
  PASTIX_INT   lftgtnum  = -1;
  PASTIX_INT   cblknum;
  PASTIX_INT   bloknum;
  PASTIX_INT   clustsrc;
  PASTIX_INT   clustdst;

  bloknum = simuctrl->tasktab[tasknum].bloknum;
  cblknum = simuctrl->tasktab[tasknum].cblknum;

  /* no fan_in_target-> no need treatment */
  if(ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum)
    return;

  ASSERTDBG(simuctrl->tasktab[tasknum].taskid != E2,MOD_BLEND);

  /*------------------------------------------------------------------------------------------------/
  / Compute the cblk on proc timer that is time the cblk would have received                        /
  / all its contributions if it was mapped on a given cand processor                                /
  / These times INCLUDE add time for fan in target !!                                               /
  /------------------------------------------------------------------------------------------------*/

  /** Compute receive time (time at which a non-local processor should received the target **/
  /* find the latest ftgt receive time and the second latest*/
  for(i=simuctrl->bloktab[bloknum].ftgtnum; i<simuctrl->bloktab[bloknum+1].ftgtnum;i++)
    {
      /* Source of this ftgt */
      clustdst = INDEX2CLUST(i,bloknum);

      /** Compute cost of the ftgt **/
      if(ctrl->candtab[cblknum].distrib == D2)
        {

          /** Task DIAG or E1 **/
          if(simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBNBR]>0)
            {
              simuctrl->ftgttab[i].costadd = costFtgtAdd(&(simuctrl->ftgttab[i].ftgt), dofptr);
              simuctrl->ftgttab[i].costsend = costFtgtSend(clustdst, ctrl->candtab[cblknum].lccandnum-ctrl->candtab[cblknum].fccandnum+1 ,&(simuctrl->ftgttab[i].ftgt), ctrl, dofptr);
            }
        }
      else
        /** Task COMP_1D with several cand proc **/
        /** The information about ftgt costs are in the ftgt of the diagonal block;
            this loop sums the cost of all the ftgt received by the blocks in this column block **/
        if(simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBNBR]>0)
          for(j=bloknum;j<symbptr->cblktab[cblknum+1].bloknum;j++)
            {
              if(simuctrl->ftgttab[simuctrl->bloktab[j].ftgtnum + i-simuctrl->bloktab[bloknum].ftgtnum].ftgt.infotab[FTGT_CTRBNBR]>0)
                {
                  simuctrl->ftgttab[i].costadd +=
                    costFtgtAdd(&(simuctrl->ftgttab[CLUST2INDEX(j, clustdst)].ftgt), dofptr);

                  simuctrl->ftgttab[i].costsend +=
                    costFtgtSend(clustdst, ctrl->candtab[cblknum].lccandnum-ctrl->candtab[cblknum].fccandnum+1,
                                 &(simuctrl->ftgttab[CLUST2INDEX(j, clustdst)].ftgt), ctrl, dofptr);
                }
            }
#ifdef DEBUG_BLEND
      if(!(simuctrl->ftgttab[i].costsend >= 0.0))
        errorPrint("ftgt %ld costsend %f", (long)i, simuctrl->ftgttab[i].costsend);
      if(!(simuctrl->ftgttab[i].costadd >= 0.0))
        errorPrint("ftgt %ld costadd %f", (long)i, simuctrl->ftgttab[i].costadd);

      ASSERTDBG(simuctrl->ftgttab[i].costsend >= 0.0,MOD_BLEND);
      ASSERTDBG(simuctrl->ftgttab[i].costadd >= 0.0,MOD_BLEND);
#endif

      /** ftgttab[].timerecv is the time this ftgt will be receive **/
      timerSet(&(simuctrl->ftgttab[i].timerecv), timerVal(&(simuctrl->ftgttimetab[i])) + simuctrl->ftgttab[i].costsend + simuctrl->ftgttab[i].costadd);

      /** If this ftgt the last reveived or the second last received ?? **/
      if(timerVal(&(simuctrl->ftgttab[i].timerecv)) > lftgttime)
        {
          /*lftgttime = simuctrl->ftgttab[i].timerecv;*/
          lftgttime = timerVal(&(simuctrl->ftgttab[i].timerecv));
          lftgtnum  = i;
        }
      else
        if(timerVal(&(simuctrl->ftgttab[i].timerecv)) > sftgttime)
          sftgttime = timerVal(&(simuctrl->ftgttab[i].timerecv));
    }


  /*------------------------------------------------------/
  / Put in ftgttimetab[] the date at which the cluster    /
  / would have received and add all the ftgt if the task  /
  /  was mapped on it                                     /
  /------------------------------------------------------*/
  for(i=simuctrl->bloktab[bloknum].ftgtnum; i<simuctrl->bloktab[bloknum+1].ftgtnum;i++)
    {
      if(i != lftgtnum)
        timerSet(&(simuctrl->ftgttimetab[i]), lftgttime);
      else
        timerSet(&(simuctrl->ftgttimetab[i]), MAX(timerVal(&(simuctrl->ftgttimetab[i])), sftgttime));

      /** Now ftgttimetab[?] is the time the cluster CLUST would have received
        all its contributions (ftgt) if the block was mapped on **/
      if(simuctrl->tasktab[tasknum].taskid == E1)
        {
          clustsrc = ctrl->proc2clust[simuctrl->blprtab[simuctrl->tasktab[tasknum].bloknum2]];
          clustdst = INDEX2CLUST(i, bloknum);
          /** The diagonal block have to be locally present for the task to be ready **/
          timerSet(&(simuctrl->ftgttimetab[i]),MAX(timerVal(&(simuctrl->ftgttimetab[i])), timerVal(&(simuctrl->tasktab[tasknum].time)) + taskSendCost(&(simuctrl->tasktab[simuctrl->bloktab[simuctrl->tasktab[tasknum].bloknum2].tasknum]), clustsrc, clustdst, ctrl)));
        }
    }
}


PASTIX_INT comp_int(const PASTIX_INT * a, const PASTIX_INT * b)
{
    if(a[0]>b[0])
        return 1;
    if(a[0]<b[0])
        return -1;

    /*a == b*/
    return 0;
}

/* OIMBE utiliser un tas pour getNextProc --> cout log(P) au lieu de P */
PASTIX_INT getNextProc(SimuProc *proctab, PASTIX_INT procnbr)
{
    double min;
    PASTIX_INT pr;
    PASTIX_INT procnum = -1;

    min = (double)INTVALMAX;
    for(pr=0;pr<procnbr;pr++)
        if((timerVal(&(proctab[pr].timer)) < min) && ( (queueSize(proctab[pr].taskheap)>0) || (queueSize(proctab[pr].taskheap2)>0) ))
            {
                min = timerVal(&(proctab[pr].timer));
                procnum = pr;
            }
    return procnum;
}

PASTIX_INT getTaskUnmapped(Queue *q1, Queue *q2, SimuCtrl *simuctrl)
{
    PASTIX_INT next = -1;
    while(queueSize(q2)>0)
        {
            next = queueGet(q2);
            switch(simuctrl->tasktab[next].taskid)
              {
              case COMP_1D:
                if(simuctrl->blprtab[simuctrl->tasktab[next].bloknum]<0)
                  {
                    ASSERTDBG(simuctrl->ownetab[simuctrl->tasktab[next].cblknum]<0,MOD_BLEND);
                    goto end;
                  }
                break;
              case E2:
                ASSERTDBG(simuctrl->blprtab[simuctrl->tasktab[next].bloknum]>=0,MOD_BLEND);
                goto end;
              case DIAG:
                if(simuctrl->blprtab[simuctrl->tasktab[next].bloknum]<0)
                  goto end;
                break;
              case E1:
                if(simuctrl->blprtab[simuctrl->tasktab[next].bloknum]<0)
                  goto end;
                break;
              default:
                errorPrint("q2, taskid %ld", (long)simuctrl->tasktab[next].taskid);
                errorPrint("in %s:%d getTaskUnmapped",__FILE__,__LINE__);
                EXIT(MOD_BLEND,INTERNAL_ERR);
              }
        }
    while(queueSize(q1)>0)
        {
            next = queueGet(q1);
            switch(simuctrl->tasktab[next].taskid)
              {
              case COMP_1D:
                if(simuctrl->blprtab[simuctrl->tasktab[next].bloknum]<0)
                {
                  ASSERTDBG(simuctrl->ownetab[simuctrl->tasktab[next].cblknum]<0,MOD_BLEND);
                  goto end;
                }
                break;
              case E2:
                ASSERTDBG(simuctrl->blprtab[simuctrl->tasktab[next].bloknum]>=0,MOD_BLEND);
                goto end;
              case DIAG:
                if(simuctrl->blprtab[simuctrl->tasktab[next].bloknum]<0)
                  goto end;
                break;
              case E1:
                if(simuctrl->blprtab[simuctrl->tasktab[next].bloknum]<0)
                  goto end;
                break;
              default:
                errorPrint("q1, taskid %ld", (long)simuctrl->tasktab[next].taskid);
                errorPrint("in %s:%d getTaskUnmapped",__FILE__,__LINE__);
                EXIT(MOD_BLEND,INTERNAL_ERR);
              }
        }
    /** no unmapped task found **/
    RETURN(-1,MOD_BLEND,INTERNAL_ERR);

end:
    return next;
}


PASTIX_INT getNextTaskNextProc(SimuCtrl *simuctrl, BlendCtrl *ctrl, PASTIX_INT *procnumptr)
{
  /*----------------------------------------------------------------------------------------------------/
  /  Get the next task and the next proc in order that they are the first that can compute something    /
  /  On return : earlier task index                                                                     /
  /----------------------------------------------------------------------------------------------------*/

  PASTIX_INT p;
  PASTIX_INT procnum = -1;
  PASTIX_INT tasknum;
  double earlytimeready = INTVALMAX;
  double earlyproctimer = INTVALMAX;
  double timeready;
  PASTIX_INT earlytask = -1;

  /** Find the earlier task in the processor heaps **/
  for(p=0;p<ctrl->procnbr;p++)
    {
      tasknum = -1;
      /** First we search the earlier task in the set of task whose ready date is < proc timer **/
      while(queueSize(simuctrl->proctab[p].taskheap2)>0)
        {
          tasknum = queueRead(simuctrl->proctab[p].taskheap2);
          if( (simuctrl->tasktab[tasknum].taskid != E2) && (simuctrl->blprtab[simuctrl->tasktab[tasknum].bloknum]>=0) )
            {
              /** This task have to be remove from the heap (already mapped) **/
              queueGet(simuctrl->proctab[p].taskheap2);
              tasknum = -1;
            }
          else
            break;
        }
      /** We found no task which ready date is < proc timer so we search one that minimizes ready date - proc-timer **/
      if(tasknum == -1)
        {
          while(queueSize(simuctrl->proctab[p].taskheap)>0)
            {
              tasknum = queueRead(simuctrl->proctab[p].taskheap);
              if( (simuctrl->tasktab[tasknum].taskid != E2) && (simuctrl->blprtab[simuctrl->tasktab[tasknum].bloknum]>=0) )
                {
                  /** This task have to be remove from the heap (already mapped) **/
                  queueGet(simuctrl->proctab[p].taskheap);
                  tasknum = -1;
                }
              else
                break;
            }
        }

      if(tasknum != -1)
        {
          timeready = MAX(timerVal(TIMER(p)), timerVal(&(simuctrl->ftgttimetab[CLUST2INDEX(simuctrl->tasktab[tasknum].bloknum, ctrl->proc2clust[p])])));

          /** We prevent to distribute on the same processor set when all time are equals **/
          if((timeready == earlytimeready) && (timerVal(TIMER(p)) < earlyproctimer))
            {
              procnum = p;
              earlyproctimer = timerVal(TIMER(p));
              earlytask = tasknum;
              earlytimeready = timeready;
            }

          if(timeready < earlytimeready)
            {
              procnum  = p;
              earlytask = tasknum;
              earlytimeready = timeready;
            }
        }
    }
  if(procnum != -1)
    {
      if(queueSize(simuctrl->proctab[procnum].taskheap2)>0)
        {ASSERT(earlytask == queueGet(simuctrl->proctab[procnum].taskheap2),MOD_BLEND);}
      else
        ASSERT(earlytask == queueGet(simuctrl->proctab[procnum].taskheap),MOD_BLEND);
    }
  *procnumptr = procnum;
  return earlytask;
}

void queueReorder(PASTIX_INT t, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl)
{
  PASTIX_INT tasknum;
  PASTIX_INT cblknum;
  PASTIX_INT facingcblk;
  PASTIX_INT procnum;

  /*for(procnum = ctrl->candtab[simuctrl->tasktab[t].cblknum].fcandnum;procnum<=ctrl->candtab[simuctrl->tasktab[t].cblknum].lcandnum;procnum++)*/
      {
        procnum = simuctrl->blprtab[simuctrl->tasktab[t].bloknum];
        while(queueSize(simuctrl->proctab[procnum].taskheap)>0)
          {
            tasknum = queueRead(simuctrl->proctab[procnum].taskheap);
            cblknum = simuctrl->tasktab[tasknum].cblknum;

#ifdef FACE
            if(simuctrl->tasktab[tasknum].taskid != COMP_1D)
              {
                if(simuctrl->tasktab[tasknum].taskid == E2)
                  facingcblk = symbptr->bloktab[simuctrl->tasktab[tasknum].bloknum2].cblknum;
                else
                  facingcblk = symbptr->bloktab[simuctrl->tasktab[tasknum].bloknum].cblknum;
              }
            else
              facingcblk = cblknum;
#endif
            /*if(ctrl->candtab[cblknum].fcandnum == ctrl->candtab[cblknum].lcandnum)
            {
                tasknum = queueGet(simuctrl->proctab[procnum].taskheap);
#ifndef FACE
                queueAdd2(simuctrl->proctab[procnum].taskheap2, tasknum, (double)ctrl->candtab[cblknum].treelevel, simuctrl->tasktab[tasknum].bloknum);
#else
                queueAdd2(simuctrl->proctab[procnum].taskheap2, tasknum, (double)ctrl->candtab[cblknum].treelevel, simuctrl->tasktab[tasknum].bloknum);
#endif
              }
            else*/
              if(!compTimer(TIMER(procnum), &(simuctrl->ftgttimetab[CLUST2INDEX(simuctrl->tasktab[tasknum].bloknum, ctrl->proc2clust[procnum])])))
                {
                  tasknum = queueGet(simuctrl->proctab[procnum].taskheap);
                  cblknum = simuctrl->tasktab[tasknum].cblknum;
#ifdef FACE
                  queueAdd2(simuctrl->proctab[procnum].taskheap2, tasknum,(double)ctrl->candtab[facingcblk].treelevel, simuctrl->tasktab[tasknum].facebloknum );
#else
                  queueAdd2(simuctrl->proctab[procnum].taskheap2, tasknum,(double)ctrl->candtab[cblknum].treelevel, simuctrl->tasktab[tasknum].bloknum );
#endif

                }
              else
                break;
          }
      }
}


void computeBlockCtrbNbr(SimuCtrl *simuctrl, SymbolMatrix *symbptr, BlendCtrl *ctrl)
{
  PASTIX_INT i, j, k;
  PASTIX_INT facebloknum;
  PASTIX_INT cursor;
  SimuTask *task;
  for(i=0;i<symbptr->cblknbr;i++)
    {
      if(simuctrl->tasktab[simuctrl->bloktab[symbptr->cblktab[i].bloknum].tasknum].taskid !=COMP_1D)
        for(j=symbptr->cblktab[i].bloknum+1;j<symbptr->cblktab[i+1].bloknum;j++)
          {
            /** Add one for DIAG contribution **/
            simuctrl->bloktab[j].ctrbcnt++;

            /** Add contribution due to E2 **/
            cursor = simuctrl->bloktab[j].tasknum+1;
            for(k=j;k<symbptr->cblktab[i+1].bloknum;k++)
              {
                ASSERTDBG(simuctrl->tasktab[cursor].taskid == E2,MOD_BLEND);
                simuctrl->bloktab[simuctrl->tasktab[cursor].facebloknum].ctrbcnt++;
                cursor++;
            }
          }
      else
        {
          /** 1D cblk computed **/
          for(j=symbptr->cblktab[i].bloknum+1;j<symbptr->cblktab[i+1].bloknum;j++)
            {
              facebloknum = 0;
              /** Add contribution due to E2 **/
              for(k=j;k<symbptr->cblktab[i+1].bloknum;k++)
                {
                  /** we don't care if facing task is 1D or 2D **/
                  facebloknum = getFaceBlockE2(facebloknum, j, k, symbptr, ctrl->option->ricar);

                  /*#ifdef NAPA*/
                  if(facebloknum >= 0)
                    /*#endif*/
                    simuctrl->bloktab[facebloknum].ctrbcnt++;
                }
            }
        }
    }

  /** Set up the task ctrbcnt and cblkcnt **/
  task = simuctrl->tasktab;
  for(i=0;i<simuctrl->tasknbr;i++)
    {
      switch(task->taskid)
        {
        case COMP_1D:
          task->ctrbcnt = 0;
          for(j=symbptr->cblktab[task->cblknum].bloknum;j<symbptr->cblktab[task->cblknum+1].bloknum;j++)
            task->ctrbcnt += simuctrl->bloktab[j].ctrbcnt;/*simuctrl->cblktab[task->cblknum].ctrbcnt;*/

          simuctrl->cblktab[task->cblknum].ctrbcnt = task->ctrbcnt;
          break;
        case DIAG:
          task->ctrbcnt = simuctrl->bloktab[task->bloknum].ctrbcnt;
          break;
        case E1:
          task->ctrbcnt = simuctrl->bloktab[task->bloknum].ctrbcnt-1;
          break;
        case E2:
          task->ctrbcnt = 0;
          break;
        default:
          errorPrint("task has no type.");
          EXIT(MOD_BLEND,INTERNAL_ERR);
        }
      task++;
    }

  if(ctrl->option->forceE2)
    for(i=0;i<simuctrl->tasknbr;i++)
      if(simuctrl->tasktab[i].taskid == E2)
        for(j=simuctrl->tasktab[i].bloknum+1;j<symbptr->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum;j++)
          simuctrl->bloktab[j].ctrbcnt++;

}

void updateFtgtStruct(PASTIX_INT bloknum, PASTIX_INT bloknum2, PASTIX_INT ftgtnum, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl)
{
  SymbolBlok * blokptr;
  SymbolBlok * blokptr2;
  (void)ctrl;

  blokptr  = &(symbptr->bloktab[bloknum]);
  blokptr2 = &(symbptr->bloktab[bloknum2]);
  simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_CTRBNBR]++;
  if(blokptr2->frownum < simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FCOLNUM])
    simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FCOLNUM] = blokptr2->frownum;
  if(blokptr2->lrownum > simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LCOLNUM])
    simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LCOLNUM] = blokptr2->lrownum;
  if(blokptr->frownum < simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FROWNUM])
    simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FROWNUM] = blokptr->frownum;
  if(blokptr->lrownum > simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LROWNUM])
    simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LROWNUM] = blokptr->lrownum;

  ASSERTDBG(simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LCOLNUM] - simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FCOLNUM]+1 > 0,MOD_BLEND);
  ASSERTDBG(simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LROWNUM] - simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FROWNUM]+1 > 0,MOD_BLEND);

}



void putInReadyQueue(PASTIX_INT procnum, PASTIX_INT tasknum, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl)
{
  /*---------------------------------------------------------/
  / This function according to the ready date of a task      /
  / put this task on the ready queue of a processor          /
  / NOTE: when the ready date of a task is inferior to the   /
  / proc timer then he task is ordered according to its      /
  / priorities in the elimination tree                       /
  /---------------------------------------------------------*/
  double ready_date = 0.0;
  PASTIX_INT bloknum;
  PASTIX_INT cblknum;
  PASTIX_INT facingcblk;

  bloknum = simuctrl->tasktab[tasknum].bloknum;
  cblknum = simuctrl->tasktab[tasknum].cblknum;
  /** Get the ready date of the task on the processor passed in parameter **/
  if(simuctrl->tasktab[tasknum].taskid == E2 || ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum)
    ready_date = timerVal(&(simuctrl->tasktab[tasknum].time));
  else
    ready_date = timerVal(&(simuctrl->ftgttimetab[CLUST2INDEX(bloknum,ctrl->proc2clust[procnum])]));

#ifdef FACE
  if(simuctrl->tasktab[tasknum].taskid != COMP_1D)
    {
      if(simuctrl->tasktab[tasknum].taskid == E2)
        facingcblk = symbptr->bloktab[simuctrl->tasktab[tasknum].bloknum2].cblknum;
      else
        facingcblk = symbptr->bloktab[simuctrl->tasktab[tasknum].bloknum].cblknum;
    }
  else
    facingcblk = cblknum;

  if(ready_date > timerVal(TIMER(procnum)))
    queueAdd2(simuctrl->proctab[procnum].taskheap, tasknum, ready_date, ctrl->candtab[facingcblk].treelevel);
  else
    queueAdd2(simuctrl->proctab[procnum].taskheap2, tasknum, (double)(ctrl->candtab[facingcblk].treelevel), bloknum);
#else
  if(ready_date > timerVal(TIMER(procnum)))
    queueAdd2(simuctrl->proctab[procnum].taskheap, tasknum, ready_date, ctrl->candtab[cblknum].treelevel);
  else
    queueAdd2(simuctrl->proctab[procnum].taskheap2, tasknum, (double)(ctrl->candtab[cblknum].treelevel), bloknum);
#endif
}
