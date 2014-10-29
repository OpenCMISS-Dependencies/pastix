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
#include <assert.h>

#include "common_pastix.h"
#include "ftgt.h"
#include "symbol.h"
#include "queue.h"
#include "bulles.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "cand.h"
#include "extendVector.h"
#include "simu.h"
#include "dof.h"
#include "elimin.h"
#include "cost.h"
/* #include "extrastruct.h" */
/* #include "param_comm.h" */
#include "param_blend.h"
#include "blendctrl.h"
/* #include "perf.h" */
/* #include "costfunc.h" */
#include "task.h"
#include "solver_check.h"

/*#define DEBUG_PRIO*/

void solverCheck(SolverMatrix *solvmtx)
{
    PASTIX_INT i, j, k = 0;
    PASTIX_INT cblknum, bloknum, ftgtnum;
    PASTIX_INT indnum, tasknum, btagnum;
/*     PASTIX_INT stride; */
    Task *task = NULL;
    PASTIX_INT  *sendcnt = NULL;
    PASTIX_INT total;
    

    BlockCoeff *bcofptr = NULL;


    /** Check the task **/
    for(i=0;i<solvmtx->tasknbr;i++)
      {
	cblknum = solvmtx->tasktab[i].cblknum;
	bloknum = solvmtx->tasktab[i].bloknum;
	indnum  = solvmtx->tasktab[i].indnum;
	ASSERT(cblknum < solvmtx->cblknbr,MOD_BLEND);
	ASSERT(bloknum < solvmtx->bloknbr,MOD_BLEND);
	/*ASSERT(solvmtx->tasktab[i].btagptr == NULL,MOD_BLEND);*/
	if(indnum >= solvmtx->indnbr)
	  printf("tasknbr %ld Tasknum %ld type %ld indnum %ld indnbr %ld \n", (long)solvmtx->tasknbr, (long)i, (long)solvmtx->tasktab[i].taskid, (long)indnum, (long)solvmtx->indnbr);
	/** OIMBE ce test foire (ligne precedente) si il y a du 1D jusqu'au bout !! mais a priori on s'en fout **/
	/*if(solvmtx->tasktab[i].taskid != 0 && solvmtx->cblktab[solvmtx->tasktab[i].cblknum] 
	  ASSERT(indnum < solvmtx->indnbr,MOD_BLEND);*/
	if(indnum >= solvmtx->indnbr)
	  printf("cblknbr %ld cblknum %ld indnum %ld indnbr %ld \n", (long)solvmtx->cblknbr, (long)solvmtx->tasktab[i].cblknum, (long)indnum, (long)solvmtx->indnbr);
	ASSERT(solvmtx->tasktab[i].taskid >= 0,MOD_BLEND);
	ASSERT(solvmtx->tasktab[i].prionum >= 0,MOD_BLEND);
	switch(solvmtx->tasktab[i].taskid)
	  {
	    case COMP_1D:
	      ASSERT(solvmtx->tasktab[i].tasknext == -1,MOD_BLEND);
	      ASSERT(bloknum  == solvmtx->cblktab[cblknum].bloknum,MOD_BLEND);
	      for(j = bloknum+1;j<solvmtx->cblktab[cblknum+1].bloknum;j++)
		{
		  for(k=j;k<solvmtx->cblktab[cblknum+1].bloknum;k++)
		    {
		      /*#ifdef NAPA*/
		      if(solvmtx->indtab[indnum] > solvmtx->ftgtnbr) /** No ftgt **/
			{
			  indnum++;
			  continue;
			}
		      /*#endif*/
		      if(solvmtx->indtab[indnum] < 0)
			{
			  tasknum = -solvmtx->indtab[indnum];
			  switch(solvmtx->tasktab[tasknum].taskid)
			      {
			      case COMP_1D:
				{
				  PASTIX_INT facebloknum, facecblknum;
				  facecblknum = solvmtx->bloktab[j].cblknum;
				  ASSERT(facecblknum >= 0,MOD_BLEND);
				  ASSERT(facecblknum == solvmtx->tasktab[tasknum].cblknum,MOD_BLEND);
                                  facebloknum = solvmtx->cblktab[facecblknum].bloknum;

                                  while ( ! (    ( solvmtx->bloktab[k].frownum >= solvmtx->bloktab[facebloknum].frownum &&
                                                   solvmtx->bloktab[k].frownum <= solvmtx->bloktab[facebloknum].lrownum)
                                              || ( solvmtx->bloktab[k].lrownum >= solvmtx->bloktab[facebloknum].frownum &&
                                                   solvmtx->bloktab[k].lrownum <= solvmtx->bloktab[facebloknum].lrownum)
                                              || ( solvmtx->bloktab[k].frownum <= solvmtx->bloktab[facebloknum].frownum &&
                                                   solvmtx->bloktab[k].lrownum >= solvmtx->bloktab[facebloknum].lrownum)))
                                    facebloknum++;

				  ASSERT(solvmtx->bloktab[k].frownum >= solvmtx->bloktab[facebloknum].frownum,MOD_BLEND);
				  ASSERT(solvmtx->bloktab[k].lrownum <= solvmtx->bloktab[facebloknum].lrownum,MOD_BLEND);
				  ASSERT(solvmtx->bloktab[j].frownum >= solvmtx->cblktab[facecblknum].fcolnum,MOD_BLEND);
				  ASSERT(solvmtx->bloktab[j].lrownum <= solvmtx->cblktab[facecblknum].lcolnum,MOD_BLEND);
				}
				break;
			      case DIAG:
			      case E1:
				{
				  PASTIX_INT facecblknum;
				  facecblknum = solvmtx->bloktab[j].cblknum;
				  ASSERT(facecblknum == solvmtx->tasktab[tasknum].cblknum,MOD_BLEND);
				  ASSERT(solvmtx->bloktab[k].frownum >= solvmtx->bloktab[solvmtx->tasktab[tasknum].bloknum].frownum,MOD_BLEND);
				  ASSERT(solvmtx->bloktab[k].lrownum <= solvmtx->bloktab[solvmtx->tasktab[tasknum].bloknum].lrownum,MOD_BLEND);
				  ASSERT(solvmtx->bloktab[j].frownum >= solvmtx->cblktab[solvmtx->tasktab[tasknum].cblknum].fcolnum,MOD_BLEND);
				  ASSERT(solvmtx->bloktab[j].lrownum <= solvmtx->cblktab[solvmtx->tasktab[tasknum].cblknum].lcolnum,MOD_BLEND);
				  break;
				}
			      case DRUNK:
				break;
			      default:
				  errorPrint("tasknum %ld tasknbr %ld taskid %ld",
					     (long)tasknum, (long)solvmtx->tasknbr, 
					     (long)solvmtx->tasktab[tasknum].taskid);
				  EXIT(MOD_BLEND,INTERNAL_ERR);
			      }
			}
		      else
			{
			  ftgtnum = solvmtx->indtab[indnum];


			  ASSERT(solvmtx->bloktab[j].frownum >= solvmtx->ftgttab[ftgtnum].infotab[FTGT_FCOLNUM],MOD_BLEND);
			  ASSERT(solvmtx->bloktab[j].lrownum <= solvmtx->ftgttab[ftgtnum].infotab[FTGT_LCOLNUM],MOD_BLEND);
			  ASSERT(solvmtx->bloktab[k].frownum >= solvmtx->ftgttab[ftgtnum].infotab[FTGT_FROWNUM],MOD_BLEND);
			  ASSERT(solvmtx->bloktab[k].lrownum <= solvmtx->ftgttab[ftgtnum].infotab[FTGT_LROWNUM],MOD_BLEND);
			  
			  /*ASSERT( (solvmtx->ftgttab[ftgtnum].infotab[FTGT_LCOLNUM] - solvmtx->ftgttab[ftgtnum].infotab[FTGT_FCOLNUM] + 1)* (solvmtx->ftgttab[ftgtnum].infotab[FTGT_LROWNUM] - solvmtx->ftgttab[ftgtnum].infotab[FTGT_FROWNUM] + 1) <= solvmtx->cpftmax,MOD_BLEND);*/

			  solvmtx->ftgttab[ftgtnum].infotab[FTGT_CTRBCNT]--;
			  
#ifdef DEBUG_PRIO
			  if( solvmtx->ftgttab[ftgtnum].infotab[FTGT_CTRBCNT] == 0)
			    if(solvmtx->ftgttab[ftgtnum].infotab[FTGT_PRIONUM] != solvmtx->tasktab[i].prionum)
			      fprintf(stdout, "Task1D %ld FTGT %ld  taskprio %ld ftgtprio %ld \n", (long)i, (long)ftgtnum, (long)solvmtx->tasktab[i].prionum, (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_PRIONUM]);
#endif
			  /*fprintf(stdout ," [ %ld %ld ] [%ld %ld ] \n", (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_FCOLNUM],
				  (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_LCOLNUM],
				  (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_FROWNUM],
				  (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_LROWNUM]);*/
			}
		      indnum++;
		    }
		}
	      break;
	  case DIAG:
	      ASSERT(solvmtx->tasktab[i].tasknext == -1,MOD_BLEND);
	      ASSERT(bloknum  == solvmtx->cblktab[cblknum].bloknum,MOD_BLEND);
	      btagnum = solvmtx->indtab[indnum];
	      if(btagnum == -1)
		break;
	      bcofptr = solvmtx->btagtab[btagnum].bcofptr;
	      ASSERT(solvmtx->bloktab[bloknum].frownum == bcofptr->infotab[BCOF_FROWNUM],MOD_BLEND);
	      ASSERT(solvmtx->bloktab[bloknum].lrownum == bcofptr->infotab[BCOF_LROWNUM],MOD_BLEND);
	      ASSERT(solvmtx->cblktab[cblknum].fcolnum == bcofptr->infotab[BCOF_FCOLNUM],MOD_BLEND);
	      ASSERT(solvmtx->cblktab[cblknum].lcolnum == bcofptr->infotab[BCOF_LCOLNUM],MOD_BLEND);
	      break;
	  case E1:
	      btagnum = solvmtx->indtab[indnum];
	      if(btagnum >=0)
		{
		  bcofptr = solvmtx->btagtab[btagnum].bcofptr;
		  ASSERT(solvmtx->bloktab[bloknum].frownum == bcofptr->infotab[BCOF_FROWNUM],MOD_BLEND);
		  ASSERT(solvmtx->bloktab[bloknum].lrownum == bcofptr->infotab[BCOF_LROWNUM],MOD_BLEND);
		  ASSERT(solvmtx->cblktab[cblknum].fcolnum == bcofptr->infotab[BCOF_FCOLNUM],MOD_BLEND);
		  ASSERT(solvmtx->cblktab[cblknum].lcolnum == bcofptr->infotab[BCOF_LCOLNUM],MOD_BLEND);
		  
		  /** Verify the chain **/
		  task = &(solvmtx->tasktab[i]);
		  bloknum = 0;
		  while(task->tasknext != i)
		    {
		      task = &(solvmtx->tasktab[task->tasknext]);
		      bloknum++;
		      ASSERT(task->taskid == E1,MOD_BLEND);
		    }
		  
		}
	      else
		{
		  fprintf(stderr, " Fucking matrix \n");
		  EXIT(MOD_BLEND,BADPARAMETER_ERR);
		}
	    break;
	  case E2:
	    if(solvmtx->indtab[indnum] < 0)
	      {
		tasknum = -solvmtx->indtab[indnum];
	      }
	    else
	      {
		ftgtnum = solvmtx->indtab[indnum];

		ASSERT(solvmtx->bloktab[bloknum].frownum >= solvmtx->ftgttab[ftgtnum].infotab[FTGT_FROWNUM],MOD_BLEND);
		ASSERT(solvmtx->bloktab[bloknum].lrownum <= solvmtx->ftgttab[ftgtnum].infotab[FTGT_LROWNUM],MOD_BLEND);
		solvmtx->ftgttab[ftgtnum].infotab[FTGT_CTRBCNT]--;
#ifdef DEBUG_PRIO
		if(solvmtx->ftgttab[ftgtnum].infotab[FTGT_CTRBCNT] == 0)
		  if(solvmtx->ftgttab[ftgtnum].infotab[FTGT_PRIONUM] != solvmtx->tasktab[i].prionum)
		    fprintf(stdout, "Task2D %ld FTGT %ld  taskprio %ld ftgtprio %ld \n", (long)i, 
			    (long)ftgtnum, (long)solvmtx->tasktab[i].prionum, (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_PRIONUM]);
#endif
		/*if(solvmtx->ftgttab[ftgtnum].infotab[FTGT_TASKDST] == -DRUNK)
		  fprintf(stderr, "Task drunk ctrbcnt %ld \n", (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_CTRBCNT]);*/
		/*fprintf(stdout ,"E2 [ %ld %ld ] [%ld %ld ] \n", (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_FCOLNUM],
				  (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_LCOLNUM],
				  (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_FROWNUM],
				  (long)solvmtx->ftgttab[ftgtnum].infotab[FTGT_LROWNUM]);*/
	      }
	    /** Verify the chain **/
	      task = &(solvmtx->tasktab[i]);
	      while(task->tasknext != i)
		{
		  task = &(solvmtx->tasktab[task->tasknext]);
		  ASSERT(task->taskid == E2,MOD_BLEND);
		}


	    break;
	  case DRUNK:
	      /*fprintf(stdout, "Task DRUNK CHECK \n");*/
	      ASSERT(i == solvmtx->tasknbr-1,MOD_BLEND);
	      break;
	  default:
	    fprintf(stderr, "solver_check: The task %ld has no type \n", (long)i);
	    EXIT(MOD_BLEND,INTERNAL_ERR);
	  }
      }
    for(i=0;i<solvmtx->ftgtnbr;i++)
      {
	ASSERT(solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR]>0,MOD_BLEND);
	ASSERT(solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT]==0,MOD_BLEND);
	/*if(solvmtx->ftgttab[i].infotab[FTGT_TASKDST] == -DRUNK)
	  fprintf(stderr, "Ftgt %ld drunk \n", (long)i);*/
      }

    /**** Check the btagptr *****/
    for(i=0;i<solvmtx->btagnbr;i++)
      {
	if(solvmtx->proc2clust[solvmtx->btagtab[i].infotab[BTAG_PROCDST]] == solvmtx->clustnum)
	  {
	    ASSERT(solvmtx->btagtab[i].infotab[BTAG_TASKDST] < solvmtx->tasknbr,MOD_BLEND);
	    ASSERT(solvmtx->tasktab[solvmtx->btagtab[i].infotab[BTAG_TASKDST] ].btagptr 
		   == &(solvmtx->btagtab[i]),MOD_BLEND);
	  }
      }


    /**To check the send counter of ftgt and btag **/
    if (solvmtx->bcofnbr != 0)
      {
	MALLOC_INTERN(sendcnt, solvmtx->bcofnbr, PASTIX_INT);
      }
    for(i=0;i<solvmtx->bcofnbr;i++)
      sendcnt[i] = solvmtx->bcoftab[i].sendcnt;

    for(i=0;i<solvmtx->btagnbr;i++)
      solvmtx->btagtab[i].bcofptr->sendcnt--;


    /** Reset the bcof->sendcnt **/
    for(i=0;i<solvmtx->bcofnbr;i++)
     {
       ASSERT(solvmtx->bcoftab[i].sendcnt == 0,MOD_BLEND);
       /*if(solvmtx->bcoftab[i].sendcnt != 0)
	 fprintf(stdout, "bcof %ld sendcntinit %ld sendcnt %ld \n", (long)i, (long)sendcnt[i], (long)solvmtx->bcoftab[i].sendcnt);*/
       solvmtx->bcoftab[i].sendcnt = sendcnt[i];
     }
    memFree(sendcnt);
    
    /** Reset the ftgt ctrbcnt **/
    for(i=0;i<solvmtx->ftgtnbr;i++)
      {
	ASSERT(solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR]>0,MOD_BLEND);
	ASSERT(solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT]==0,MOD_BLEND);
	solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT] = solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR];
      }

    /** Test the task partition on the thread of the cluster **/
    total = 0;
    for(i=0;i<solvmtx->bublnbr;i++){
      printf("i = %d, ttsknbr = %d\n", i, solvmtx->ttsknbr[i]);
      total += solvmtx->ttsknbr[i];
    }
    if(total != solvmtx->tasknbr)
      fprintf(stderr, " total %ld tasknbr %ld \n", (long)total, (long)solvmtx->tasknbr);

    ASSERT(total == solvmtx->tasknbr,MOD_BLEND);

#ifdef DEBUG_BLEND
    {
      PASTIX_INT * flag;
      MALLOC_INTERN(flag, solvmtx->tasknbr, PASTIX_INT);
      bzero(flag, sizeof(PASTIX_INT)*solvmtx->tasknbr);
      
      for(i=0;i<solvmtx->bublnbr;i++)
	for(j=0;j<solvmtx->ttsknbr[i];j++)
	  {
	    if(flag[solvmtx->ttsktab[i][j]] != 0)
	      fprintf(stderr, "flag %ld thread %ld task %ld already on another thread \n", (long)flag[solvmtx->ttsktab[i][j]], (long)i, (long)solvmtx->ttsktab[i][j]);
	    flag[solvmtx->ttsktab[i][j]]++;
	    
	  }
      
      for(i=0;i<solvmtx->tasknbr;i++)
	ASSERT(flag[i] == 1,MOD_BLEND);
      
      memFree(flag);
    }
#endif

#ifdef PASTIX_DYNSCHED
    {
      int k, bubnbr = solvmtx->bublnbr;
      
      for(k = 0; k<bubnbr; k++)
        {
          int father = BFATHER( solvmtx->btree, k );
          if ( (father != -1) && 
               (  solvmtx->btree->nodetab[k].priomax >  solvmtx->btree->nodetab[father].priomin ) )
            fprintf(stderr, "We have a problem of task distribution\n"
                    " Bubble[%d] priorities (%ld,%ld) intersect with bubble[%d](%ld,%ld)\n",
                    k,      (long)solvmtx->btree->nodetab[k].priomin,      (long)solvmtx->btree->nodetab[k].priomax,
                    father, (long)solvmtx->btree->nodetab[father].priomin, (long)solvmtx->btree->nodetab[father].priomax) ;
        }
    }
#endif

}
