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
#include <sys/stat.h>
#include <unistd.h>

#include "common_pastix.h"
#include "elimin.h"
#include "cost.h"
#include "extendVector.h"
#include "cand.h"
#include "queue.h"
#include "bulles.h"
#include "param_blend.h"
/* #include "param_comm.h" */
/* #include "symbol.h" */
/* #include "ftgt.h" */
/* #include "simu.h" */
/* #include "assert.h" */
#include "blendctrl.h"
#include "perf.h"


/*void perfcluster(PASTIX_INT procsrc, PASTIX_INT procdst, netperf *np, BlendCtrl *ctrl)
{
#ifdef DEBUG_BLEND
  ASSERT(procsrc>=0 && procsrc < ctrl->procnbr,MOD_BLEND);
  ASSERT(procdst>=0 && procdst < ctrl->procnbr,MOD_BLEND);
#endif

  if(procsrc == procdst)
    {
      np->startup   = 0;
      np->bandwidth = 0;
      return;
    }
    
  if(CLUSTNUM(procsrc) == CLUSTNUM(procdst))
  {
      np->startup   = TIME_STARTUP_1;
      np->bandwidth = TIME_BANDWIDTH_1;
      return;
    }
  else
    {
      np->startup   = CLUSTER_STARTUP_1;
      np->bandwidth = CLUSTER_BANDWIDTH_1;
      return;
    }
}
*/

void perfcluster2(PASTIX_INT clustsrc, PASTIX_INT clustdst, PASTIX_INT sync_comm_nbr, netperf *np, BlendCtrl *ctrl)
{
#ifdef DEBUG_BLEND
  ASSERT(clustsrc>=0 && clustsrc < ctrl->clustnbr,MOD_BLEND);
  ASSERT(clustdst>=0 && clustdst < ctrl->clustnbr,MOD_BLEND);
  ASSERT(sync_comm_nbr>0 && sync_comm_nbr <= ctrl->clustnbr,MOD_BLEND);
#endif


  if(clustsrc == clustdst)
    {
      np->startup   = 0;
      np->bandwidth = 0;
      return;
    }
  
      
  if(SMPNUM(clustsrc) == SMPNUM(clustdst))
    {
      
      /*fprintf(stderr, "MPI_SHARED for %ld \n", (long)sync_comm_nbr);*/
      if(sync_comm_nbr <=2)
	{
	  np->startup   = SHARED_STARTUP_1;
	  np->bandwidth = SHARED_BANDWIDTH_1;
	  return;
	}
      if(sync_comm_nbr<=4)
	{
	  np->startup   = SHARED_STARTUP_2;
	  np->bandwidth = SHARED_BANDWIDTH_2;
	  return;
	}
      if(sync_comm_nbr<=8)
	{
	  np->startup   = SHARED_STARTUP_4;
	  np->bandwidth = SHARED_BANDWIDTH_4;
	  return;
	}
      if(sync_comm_nbr > 8)
	{
	  /*fprintf(stdout, "intra %ld extra %ld\n", intra, extra);*/
	  np->startup   = SHARED_STARTUP_8;
	  np->bandwidth = SHARED_BANDWIDTH_8;
	  return;
	}
      
    }
  else
    {
      /*fprintf(stderr, "MPI for %ld \n", (long)sync_comm_nbr);*/
      /*      extra++;*/
      if(sync_comm_nbr<=2)
	{
	  np->startup   = CLUSTER_STARTUP_1;
	  np->bandwidth = CLUSTER_BANDWIDTH_1;
	  return;
	}
      if(sync_comm_nbr<=4)
	{
	  np->startup   = CLUSTER_STARTUP_2;
	  np->bandwidth = CLUSTER_BANDWIDTH_2;
	  return;
	}
      if(sync_comm_nbr<=8)
	{
	  np->startup   = CLUSTER_STARTUP_4;
	  np->bandwidth = CLUSTER_BANDWIDTH_4;
	  return;
	}
      if(sync_comm_nbr > 8)
	{
	  /*fprintf(stdout, "intra %ld extra %ld\n", intra, extra);*/
	  np->startup   = CLUSTER_STARTUP_8;
	  np->bandwidth = CLUSTER_BANDWIDTH_8;
	  return;
	}
    }
  
}

PASTIX_INT blendCtrlInit(BlendCtrl *ctrl,
                  PASTIX_INT clustnbr,
                  PASTIX_INT thrdlocnbr,
                  PASTIX_INT cudanbr,
                  PASTIX_INT clustnum,
                  BlendParam *param)
{
    PASTIX_INT i;

    ctrl->option  = param;
    MALLOC_INTERN(ctrl->perfptr, 1, netperf);
    
    /* Nombre et numéro de processus MPI */
    ctrl->clustnbr   = clustnbr;
    ctrl->clustnum   = clustnum;
    ctrl->cudanbr    = cudanbr;

#ifdef PASTIX_DYNSCHED
    /* Le nombre de coeur par cpu est est donné par iparm[IPARM_CPU_BY_NODE] 
       et a defaut par sysconf(_SC_NPROCESSORS_ONLN)                          */
    if (ctrl->option->iparm[IPARM_CPU_BY_NODE] != 0)
      ctrl->option->procnbr = ctrl->option->iparm[IPARM_CPU_BY_NODE];

    /* Calcul du nombre de cpu total a notre disposition */
    if (ctrl->option->smpnbr)
      ctrl->procnbr = ctrl->option->procnbr * ctrl->option->smpnbr; 
    else
      ctrl->procnbr = ctrl->option->procnbr * clustnbr; 

    /* Cas ou on a moins de proc que de processus MPI */
    if (ctrl->clustnbr > ctrl->procnbr)
      {
	errorPrintW("blenctrlinit: plus de processus MPI que de processeurs disponible");
 	ctrl->procnbr = ctrl->clustnbr;
/* 	ctrl->option->procnbr = (int)ceil((double)ctrl->procnbr / (double)ctrl->option->smpnbr); */
      }
    
    /* Nombre de processeurs utilisés pour un processus MPI */
    /* et nombre de threads demandés pour un processus MPI  */
    ctrl->proclocnbr = ctrl->procnbr / ctrl->clustnbr;
    ctrl->thrdlocnbr = thrdlocnbr;
#else

    ctrl->proclocnbr = thrdlocnbr;
    ctrl->thrdlocnbr = thrdlocnbr;
    ctrl->procnbr    = ctrl->proclocnbr * ctrl->clustnbr;
#endif
 
    ctrl->thrdnbr = ctrl->thrdlocnbr * clustnbr;
    ctrl->bublnbr = ctrl->thrdlocnbr;

    /* Tableau d'affectation de processeur par processus MPI */
    MALLOC_INTERN(ctrl->proc2clust, ctrl->procnbr, PASTIX_INT);
    for(i=0;i<ctrl->procnbr;i++)
      ctrl->proc2clust[i] = CLUSTNUM(i);

    ctrl->egraph  = NULL;
    ctrl->etree   = NULL;
    ctrl->costmtx = NULL;
    ctrl->candtab = NULL;     
    MALLOC_INTERN(ctrl->lheap, 1, Queue); 
    queueInit(ctrl->lheap, 1000);
    MALLOC_INTERN(ctrl->intvec, 1, ExtendVectorINT);
    MALLOC_INTERN(ctrl->intvec2, 1, ExtendVectorINT);
    if(param->tracegen)
      OUT_OPENFILEINDIR(param->iparm, ctrl->tracefile, param->trace_filename, "w");

#ifdef PASTIX_DYNSCHED
    MALLOC_INTERN(ctrl->btree, 1, BubbleTree);
#endif
    return ( (extendint_Init(ctrl->intvec, 10) != NULL) && (extendint_Init(ctrl->intvec2, 10) != NULL));
}


void blendCtrlExit(BlendCtrl *ctrl)
{
  if(ctrl->option->tracegen)
    OUT_CLOSEFILEINDIR(ctrl->tracefile);
  if(ctrl->perfptr)
    memFree_null(ctrl->perfptr);
  queueExit(ctrl->lheap);
  memFree_null(ctrl->lheap);
  extendint_Exit(ctrl->intvec);
  memFree_null(ctrl->intvec);
  extendint_Exit(ctrl->intvec2);
  memFree_null(ctrl->intvec2);

  if(ctrl->proc2clust)
    memFree_null(ctrl->proc2clust);
  if(ctrl->candtab)
    memFree_null(ctrl->candtab);

  memFree_null(ctrl);
}

