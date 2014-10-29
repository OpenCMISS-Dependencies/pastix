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
#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#ifdef X_ARCHi686_pc_linux
#  include <sched.h>
#  ifdef X_ARCHi686_mac
#    include <mach/thread_act.h>
#    include <mach/mach_init.h>
#  endif
#endif
#ifdef FORCE_NOMPI
#  include "nompi.h"
#else
#  include <mpi.h>
#endif

#include "common_pastix.h"
#include "out.h"
#ifdef PASTIX_EZTRACE
#  include "pastix_eztrace.h"
#else
#  include "trace.h"
#endif
#include "sopalin_define.h"
#include "symbol.h"
#include "ftgt.h"
#include "csc.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
#include "sopalin_thread.h"
#include "sopalin_acces.h"
#include "sopalin_time.h"
#include "stack.h"
#include "sopalin3d.h"
#include "coefinit.h"
#include "sopalin_init.h"
#include "ooc.h"

#define print_onempi(fmt, ...) if( SOLV_PROCNUM == 0 )           fprintf(stdout, fmt, __VA_ARGS__)
#define print_one(fmt, ...)    if( me == 0 && SOLV_PROCNUM == 0) fprintf(stdout, fmt, __VA_ARGS__)
#define print_all(fmt, ...)    fprintf(stdout, fmt, ##__VA_ARGS__)
#define print_error(...)

/*********************************/
/*
  Function: sopalin_init

  Allocate fields of the sopalin_data structures and initialize them.
  This function is mono-thread and must be called just by one thread.

  Parameters:
    sopalin_data - Sopalin_data structure to initialize
    datacode     - SolverMatrix structure (common data)
    sopaparam    - sopalin parameters.
    fact         - Boolean for factorisation step or not

  Returns:
    void
 */
/*********************************/
void sopalin_init(Sopalin_Data_t *sopalin_data,
                  SolverMatrix   *m,
                  SopalinParam   *sopaparam,
                  int             fact)
{
  SolverMatrix *datacode;
  PASTIX_INT           i, j, task;

  /* Initialize global var */
  if (m != NULL)
    {
      sopalin_data->datacode  = m;
      sopalin_data->sopar     = sopaparam;

      datacode = sopalin_data->datacode;

      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        print_onempi("%s", OUT2_SOP_BINITG);

      sopalin_data->thread_data       = NULL;
      sopalin_data->fanintgtsendqueue = NULL;
      sopalin_data->blocktgtsendqueue = NULL;
      sopalin_data->taskmark          = NULL;
      /* Trace */
#ifdef TRACE_SOPALIN
      sopalin_data->tracefile   = NULL;
      sopalin_data->timestamp   = 0.0;
#endif

      /* Alocation FTgt */
#if (defined COMPUTE_ALLOC) || (defined STATS_SOPALIN)
      sopalin_data->current_alloc = 0;
#endif
#ifdef ALLOC_FTGT
      sopalin_data->max_alloc   = 0;
      sopalin_data->alloc_init  = 0;
#  ifdef STATS_SOPALIN
      pthread_mutex_init(&(sopalin_data->mutex_alloc),NULL);
#  endif
#endif

      /* Variables nécessaire à la version CSC */
#ifdef USE_CSC
      sopalin_data->critere       = 0;
      sopalin_data->stop          = 0;
      sopalin_data->berr          = 0;
      sopalin_data->lberr         = 0;
      sopalin_data->raffnbr       = 0;
      sopalin_data->count_iter    = 0;
      sopalin_data->flag_gmres    = 1;
#endif

#ifdef SMP_SOPALIN
      sopalin_data->mutex_task  = NULL;
      sopalin_data->cond_task   = NULL;
      sopalin_data->mutex_fanin = NULL;
      sopalin_data->cond_fanin  = NULL;
      sopalin_data->mutex_blok  = NULL;
      sopalin_data->mutex_queue_fanin = NULL;
      sopalin_data->mutex_queue_block = NULL;
#endif

#ifdef STORAGE
      sopalin_data->grhs          = NULL;
      sopalin_data->flagtab       = NULL;
      sopalin_data->mutex_flagtab = NULL;
      sopalin_data->cond_flagtab  = NULL;
#endif

#ifdef PASTIX_DYNSCHED
      sopalin_data->tasktab_mutex  = NULL;
      sopalin_data->tasktab_cond   = NULL;
      sopalin_data->tasktab_indice = NULL;
#endif

      sopalin_data->common_flt = NULL;
      sopalin_data->common_dbl = NULL;
      sopalin_data->ptr_csc    = NULL;
      /*
       * Allocation du tableau de structure de thread
       */
      {
        int threadnbr = SOLV_THRDNBR;
#ifdef WITH_STARPU
        if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
          threadnbr += sopalin_data->sopar->iparm[IPARM_CUDA_NBR];
        else
#endif
          threadnbr += sopaparam->nbthrdcomm;

        MALLOC_INTERN(sopalin_data->thread_data,
                      threadnbr,
                      Thread_Data_t*);
      for (i=0;i<threadnbr;i++)
        sopalin_data->thread_data[i] = NULL;
      }


      /*
       * Allocation des données nécessaires à la version multi-thread
       */
#ifdef SMP_SOPALIN
      /* Tableau de mutex pour la protection des tâches et des comms */
      MALLOC_INTERN(sopalin_data->cond_task,
                    SOLV_TASKNBR,
                    pthread_cond_t);
      MALLOC_INTERN(sopalin_data->mutex_task,
                    SOLV_TASKNBR,
                    pthread_mutex_t);
      if (SOLV_FTGTNBR != 0)
        MALLOC_INTERN(sopalin_data->mutex_fanin,
                      SOLV_FTGTNBR,
                      pthread_mutex_t);

      /* Initialisation des mutex et des conditions */
      for (i=0;i<SOLV_TASKNBR;i++)
        {
          pthread_mutex_init(&(sopalin_data->mutex_task[i]),NULL);
          pthread_cond_init(&(sopalin_data->cond_task[i]), NULL);
        }
      for (i=0;i<SOLV_FTGTNBR;i++)
        pthread_mutex_init(&(sopalin_data->mutex_fanin[i]),NULL);

      pthread_mutex_init(&(sopalin_data->mutex_raff),NULL);
      pthread_cond_init(&(sopalin_data->cond_raff),NULL);

      sopalin_data->barrier.instance = 0;
      sopalin_data->barrier.blocked_threads  = 0;
      pthread_mutex_init(&(sopalin_data->barrier.sync_lock), NULL);
      pthread_cond_init(&(sopalin_data->barrier.sync_cond), NULL);
#endif /* SMP_SOPALIN */

      /*
       * Données pour les threads de communications dédiés
       */
      if (THREAD_COMM_ON)
        {
          pthread_mutex_init(&(sopalin_data->mutex_comm), NULL);
          pthread_cond_init(&(sopalin_data->cond_comm), NULL);

          sopalin_data->step_comm = COMMSTEP_INIT;
          print_debug(DBG_THCOMM, "%s:%d INIT\n", __FILE__, __LINE__);

          if (THREAD_FUNNELED_ON)
            {
              MALLOC_INTERN(sopalin_data->sendqueue, 1, Queue);
              if (sopaparam->iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
                queueInit(sopalin_data->sendqueue, MAX(MAX(UPDOWN_UPMSGNBR,UPDOWN_DOWNMSGNBR),\
                                                       SOLV_FTGTNBR+SOLV_BTAGNBR));
              else
                queueInit(sopalin_data->sendqueue, SOLV_FTGTNBR+SOLV_BTAGNBR);
            }
        }

      /*
       * Allocation pour l'ordonnancement dynamique
       */
#ifdef PASTIX_DYNSCHED
      MALLOC_INTERN(sopalin_data->tasktab_mutex,
                    SOLV_BUBLNBR,
                    pthread_mutex_t);
      MALLOC_INTERN(sopalin_data->tasktab_cond,
                    SOLV_BUBLNBR,
                    pthread_cond_t);
      MALLOC_INTERN(sopalin_data->tasktab_indice,
                    SOLV_BUBLNBR,
                    PASTIX_INT);

      for (i=0; i<SOLV_BUBLNBR; i++)
        {
          pthread_mutex_init(&(sopalin_data->tasktab_mutex[i]), NULL);
          pthread_cond_init(&(sopalin_data->tasktab_cond[i]), NULL);
          sopalin_data->tasktab_indice[i] = 0;
        }

      /* Allocation du tableau de file de tâches */
      MALLOC_INTERN(sopalin_data->taskqueue,
                    sopalin_data->datacode->btree->nodenbr,
                    Queue);
      for(i=0; i<sopalin_data->datacode->btree->nodenbr; i++)
        sopalin_data->taskqueue[i].size = 0;
#endif

#ifdef OOC
      ooc_init(sopalin_data, sopalin_data->sopar->iparm[IPARM_OOC_LIMIT]);
#endif

#ifndef WITH_HWLOC
#  ifdef PASTIX_GET_SCHED_AFFINITY
      {
	int k,len, mpisize, nbproc;
	char procname[MPI_MAX_PROCESSOR_NAME];
	int color, key;
	MPI_Comm intra_node_comm;
	int intra_node_rank, intra_node_size;
	cpu_set_t mask;

	MPI_Get_processor_name(procname,&len);
	MPI_Comm_rank(sopalin_data->sopar->pastix_comm, &(key));
	color = 0;
	for (i = 0; i < len; i++)
	  color = color*256*sizeof(char) + procname[i];
	MPI_Comm_split(sopalin_data->sopar->pastix_comm, color, key, &intra_node_comm);
	MPI_Comm_rank(intra_node_comm, &intra_node_rank);
	MPI_Comm_size(intra_node_comm, &intra_node_size);

	CPU_ZERO(&mask);
	if (sched_getaffinity(0, sizeof(mask), &mask) < 0)
	  {
	    perror("sched_getaffinity");
	  }

	sopalin_data->ncore_avail = 0;
#    ifdef MARCEL
        nbproc = marcel_nbvps();
#    else
        nbproc = sysconf(_SC_NPROCESSORS_ONLN);
#    endif

	for (i = 0; i < nbproc; i++)
	  {
	    int cpu;
	    cpu = CPU_ISSET(i, &mask);
	    if (cpu)
	      sopalin_data->ncore_avail++;
	  }
	MALLOC_INTERN(sopalin_data->allowed_cpus, SOLV_THRDNBR, int);
	for (i = 0, j= 0, k = 0; i < nbproc; i++)
	  {
	    int cpu;
	    cpu = CPU_ISSET(i, &mask);
	    if (cpu)
	      {
		if (k<(intra_node_rank*SOLV_THRDNBR)%sopalin_data->ncore_avail)
		  {
		    k++;
		    continue;
		  }
		sopalin_data->allowed_cpus[j++] = i;
		fprintf(stdout, "%d : core %d  local (%d)\n", SOLV_PROCNUM, i, j);
	      }
	    if (j == MIN(SOLV_THRDNBR, sopalin_data->ncore_avail)) break;
	  }
	if (j == sopalin_data->ncore_avail)
	  {
	    for (i = sopalin_data->ncore_avail; i < SOLV_THRDNBR; i++)
	      {
		k = ((i - i%sopalin_data->ncore_avail)/sopalin_data->ncore_avail)-1;
		sopalin_data->allowed_cpus[k] = sopalin_data->allowed_cpus[i%sopalin_data->ncore_avail];
		fprintf(stdout, "%d : core %d  local (%d)\n", SOLV_PROCNUM, sopalin_data->allowed_cpus[i%sopalin_data->ncore_avail], k);
	      }

	  }
      }
#  endif /* PASTIX_GET_SCHED_AFFINITY */
#endif /* not WITH_HWLOC */
    }

  /* Fin initialisation Commune */
  datacode = sopalin_data->datacode;

  /*
   * Specifique a la factorisation
   */
  if (fact)
    {
      /*
       * Allocation des données nécessaires à la version monothread
       */
#ifndef SMP_SOPALIN
      /* Allocation de la file de tâches */
      queueInit(&(sopalin_data->taskqueue),SOLV_TASKNBR);
      for (i=0;i<SOLV_TASKNBR;i++)
        queueAdd(&(sopalin_data->taskqueue),i,((double)TASK_PRIONUM(i)));
#else

      /* Fanin and Block */
      if (SOLV_FTGTNBR != 0)
        MALLOC_INTERN(sopalin_data->cond_fanin,
                      SOLV_FTGTNBR,
                      pthread_cond_t);
      MALLOC_INTERN(sopalin_data->mutex_blok,
                    SYMB_BLOKNBR,
                    pthread_mutex_t);
      MALLOC_INTERN(sopalin_data->mutex_queue_fanin,
                    SOLV_PROCNBR,
                    pthread_mutex_t);
      MALLOC_INTERN(sopalin_data->mutex_queue_block,
                    SOLV_PROCNBR,
                    pthread_mutex_t);

      for (i=0;i<SOLV_FTGTNBR;i++)
        pthread_cond_init(&(sopalin_data->cond_fanin[i]),NULL);

      for (i=0;i<SYMB_BLOKNBR;i++)
        pthread_mutex_init(&(sopalin_data->mutex_blok[i]),NULL);
      for (i=0;i<SOLV_PROCNBR;i++)
        {
          pthread_mutex_init(&(sopalin_data->mutex_queue_fanin[i]),NULL);
          pthread_mutex_init(&(sopalin_data->mutex_queue_block[i]),NULL);
        }
#endif /* SOPALIN_SMP */

      /* Allocation des structures necessaires aux communications */
      if (SOLV_PROCNBR > 1)
        {
          /* Fanin */
          MALLOC_INTERN(sopalin_data->fanintgtsendqueue,
                        SOLV_PROCNBR,
                        Queue);

          /* Block */
          MALLOC_INTERN(sopalin_data->blocktgtsendqueue,
                        SOLV_PROCNBR,
                        Queue);

          for (i=0;i<SOLV_PROCNBR;i++)
            queueInit(&(sopalin_data->fanintgtsendqueue[i]),SOLV_FTGTNBR);
          for (i=0;i<SOLV_PROCNBR;i++)
            queueInit(&(sopalin_data->blocktgtsendqueue[i]),SOLV_BTAGNBR);


#ifndef COMM_REORDER
          for (i=0;i<SOLV_FTGTNBR;i++)
            queueAdd2(&(sopalin_data->fanintgtsendqueue[FANIN_PROCDST(i)]), i,
                      ((double)FANIN_PRIONUM(i)), i);
          for (i=0;i<SOLV_BTAGNBR;i++)
            queueAdd(&(sopalin_data->blocktgtsendqueue[BTAG_PROCDST(i)]), i,
                     ((double)BTAG_PRIONUM(i)));
#endif
        }

      /* mark task first btag */
      MALLOC_INTERN(sopalin_data->taskmark,
                    SOLV_TASKNBR, PASTIX_INT);

      for (i=0;i<SOLV_TASKNBR;i++)
        sopalin_data->taskmark[i]=-1;
#ifndef PASTIX_DYNSCHED
      for (i=0;i<SOLV_TASKNBR;i++)
        if (sopalin_data->taskmark[i]==-1)
          {
            PASTIX_INT mintask,minprio;
            sopalin_data->taskmark[i] = 0;
            task   = i;
            mintask = i;
            minprio = TASK_PRIONUM(i);
            while ((TASK_TASKNEXT(task) != i) && (TASK_TASKNEXT(task) != -1))
              {
                task = TASK_TASKNEXT(task);
                sopalin_data->taskmark[task] = 0;
                if (TASK_PRIONUM(task)<minprio)
                  {
                    minprio = TASK_PRIONUM(task);
                    mintask = task;
                  }
              }
            sopalin_data->taskmark[mintask] = 1;
          }
#endif
    }
  /*
   * Specifique au solve
   */
  else
    {
#ifdef STORAGE
      /* Allocate flagtab structures */
      MALLOC_INTERN(sopalin_data->mutex_flagtab,
                    UPDOWN_GCBLKNBR,
                    pthread_mutex_t);

      for (i=0;i<UPDOWN_GCBLKNBR;i++)
        pthread_mutex_init(&(sopalin_data->mutex_flagtab[i]),NULL);

      /* Données non-nécessaire à la facto et que l'on doit donc obligatoirement initialiser */
      MALLOC_INTERN(sopalin_data->cond_flagtab,
                    UPDOWN_GCBLKNBR,
                    pthread_cond_t);

      for (i=0;i<UPDOWN_GCBLKNBR;i++)
        pthread_cond_init(&(sopalin_data->cond_flagtab[i]), NULL);

      MALLOC_INTERN(sopalin_data->grhs,    UPDOWN_GNODENBR*UPDOWN_SM2XNBR, PASTIX_FLOAT);
      MALLOC_INTERN(sopalin_data->flagtab, UPDOWN_GCBLKNBR,                PASTIX_INT);

#endif /* STORAGE */

      /* Structures pour le rafinement */
      MALLOC_INTERN(sopalin_data->common_flt, 2*SOLV_THRDNBR, PASTIX_FLOAT);
      MALLOC_INTERN(sopalin_data->common_dbl, 2*SOLV_THRDNBR, double);
      MALLOC_INTERN(sopalin_data->ptr_csc,    SOLV_THRDNBR,   void *);

      for(i = 0; i < SOLV_THRDNBR; i++)
        {
          sopalin_data->ptr_csc[i]        = NULL;
          sopalin_data->common_flt[2*i]   = 0.0;
          sopalin_data->common_flt[2*i+1] = 0.0;
          sopalin_data->common_dbl[2*i]   = 0.0;
          sopalin_data->common_dbl[2*i+1] = 0.0;
        }

#ifndef OOC
      /* SYMB_CBLKNUM(x)=indice de la fanin en negatif */
      for (task=0;task<SOLV_TASKNBR;task++)
        {
          PASTIX_INT c,n,t;
          c = TASK_CBLKNUM(task);
          n = 0;
          for (i=SYMB_BLOKNUM(c)+1;i<SYMB_BLOKNUM(c+1);i++)
            {
              for (j=i;j<SYMB_BLOKNUM(c+1);j++)
                if ((t=SOLV_INDTAB[TASK_INDNUM(task)+(n++)])>=0)
                  {
                    if (i==j) SYMB_CBLKNUM(i)=-t;
                  }
            }
        }
#endif
    }

  /* Initialisation du timestamp pour les traces */
#ifdef TRACE_SOPALIN
  if (sopalin_data->timestamp == 0)
    {

      CALL_MPI MPI_Barrier(PASTIX_COMM);
      TEST_MPI("MPI_Barrier");
      sopalin_data->timestamp = clockGet();

      if (SOLV_PROCNBR < 100 )
        {
          char filename[12];
          sprintf(filename,"traceGen.%02d", (int)SOLV_PROCNUM);
          filename[11] = '\0';
          OUT_OPENFILEINDIR(sopaparam->iparm, sopalin_data->tracefile, filename, "w");
        }
      else
        {
          errorPrint("Trace impossible car trop de proc.");
          EXIT(MOD_SOPALIN, BADPARAMETER_ERR);
        }

      memAllocTrace(sopalin_data->tracefile, sopalin_data->timestamp, SOLV_PROCNUM);

      trace_start(sopalin_data->tracefile, SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, -1);
    }
#endif

  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_onempi("%s", OUT2_SOP_EINITG);
}

/*********************************/
/*
  Function: sopalin_clean

  Clean the fields in sopalin_data structures.
  If step = 1, fields used only for factorization are cleaned.
  If step = 2, all fields are freed.
  This function is mono-thread and must be called just by one thread.

  Parameters:
    sopalin_data - Sopalin_data structure to clean
    step         - Cleaning step (1 or 2, partial or complete)

  Returns:
    void
 */
/*********************************/
void sopalin_clean(Sopalin_Data_t *sopalin_data, int step)
{
  PASTIX_INT i;
  SolverMatrix *datacode = sopalin_data->datacode;

  if (step == 1)
    {
      if (sopalin_data->fanintgtsendqueue != NULL)
        {
          for (i=0;i<SOLV_PROCNBR;i++)
            queueExit(&(sopalin_data->fanintgtsendqueue[i]));
          memFree_null(sopalin_data->fanintgtsendqueue);
        }
      if (sopalin_data->blocktgtsendqueue != NULL)
        {
          for (i=0;i<SOLV_PROCNBR;i++)
            queueExit(&(sopalin_data->blocktgtsendqueue[i]));
          memFree_null(sopalin_data->blocktgtsendqueue);
        }

#ifndef SMP_SOPALIN
      queueExit(&(sopalin_data->taskqueue));
#else
      if (sopalin_data->mutex_blok != NULL)
        {
          for (i=0;i<SYMB_BLOKNBR;i++)
            pthread_mutex_destroy(&(sopalin_data->mutex_blok[i]));
          memFree_null(sopalin_data->mutex_blok);
        }
      if (sopalin_data->mutex_queue_fanin != NULL)
        {
          for (i=0;i<SOLV_PROCNBR;i++)
            pthread_mutex_destroy(&(sopalin_data->mutex_queue_fanin[i]));
          memFree_null(sopalin_data->mutex_queue_fanin);
        }
      if (sopalin_data->mutex_queue_block != NULL)
        {
          for (i=0;i<SOLV_PROCNBR;i++)
            pthread_mutex_destroy(&(sopalin_data->mutex_queue_block[i]));
          memFree_null(sopalin_data->mutex_queue_block);
        }

      pthread_mutex_destroy(&(sopalin_data->mutex_raff));
      pthread_cond_destroy(&(sopalin_data->cond_raff));
#  ifdef STATS_SOPALIN
      pthread_mutex_destroy(&(sopalin_data->mutex_alloc));
#  endif
#endif

      if (sopalin_data->taskmark != NULL)
        memFree_null(sopalin_data->taskmark);
    }

  if (step == 2)
    {
#ifdef OOC
      ooc_exit(sopalin_data);
#endif

#ifdef SMP_SOPALIN
#  ifdef TRYLOCK
      for (i=1; i<SOLV_THRDNBR; i++)
        {
          sopalin_data->thread_data[0]->ptbusy += sopalin_data->thread_data[i]->ptbusy;
          sopalin_data->thread_data[0]->ptfree += sopalin_data->thread_data[i]->ptfree;
          sopalin_data->thread_data[0]->ptwait += sopalin_data->thread_data[i]->ptwait;
        }

      printf("MUTEX : busy %ld , %f\n", sopalin_data->thread_data[0]->ptbusy,
             (double)sopalin_data->thread_data[0]->ptbusy / (double)(sopalin_data->thread_data[0]->ptbusy+
                                                                     sopalin_data->thread_data[0]->ptfree) *100.);
      printf("MUTEX : free %ld , %f\n", sopalin_data->thread_data[0]->ptfree,
             (double)sopalin_data->thread_data[0]->ptfree / (double)(sopalin_data->thread_data[0]->ptbusy+
                                                                     sopalin_data->thread_data[0]->ptfree) *100.);
      printf("MUTEX : wait %ld\n", sopalin_data->thread_data[0]->ptwait);
#  endif

      if (sopalin_data->mutex_task != NULL)
        {
          for (i=0;i<SOLV_TASKNBR;i++)
            {
              pthread_mutex_destroy(&(sopalin_data->mutex_task[i]));
              pthread_cond_destroy(&(sopalin_data->cond_task[i]));
            }
          memFree_null(sopalin_data->mutex_task);
          memFree_null(sopalin_data->cond_task);
        }
      if (sopalin_data->mutex_fanin != NULL)
        {
          for (i=0;i<SOLV_FTGTNBR;i++)
            {
              pthread_mutex_destroy(&(sopalin_data->mutex_fanin[i]));
            }
          memFree_null(sopalin_data->mutex_fanin);
        }
      if (sopalin_data->cond_fanin != NULL)
        {
          for (i=0;i<SOLV_FTGTNBR;i++)
            {
              pthread_cond_destroy(&(sopalin_data->cond_fanin[i]));
            }
          memFree_null(sopalin_data->cond_fanin);
        }

      pthread_mutex_destroy(&(sopalin_data->barrier.sync_lock));
      pthread_cond_destroy(&(sopalin_data->barrier.sync_cond));
#endif /* SMP_SOPALIN */

      if (THREAD_COMM_ON)
        {
          pthread_mutex_destroy(&(sopalin_data->mutex_comm));
          pthread_cond_destroy(&(sopalin_data->cond_comm));
          if (THREAD_FUNNELED_ON)
            {
              queueExit(sopalin_data->sendqueue);
              memFree_null(sopalin_data->sendqueue);
            }
        }

#if (defined PASTIX_DYNSCHED)
      for (i=0; i<SOLV_BUBLNBR; i++)
        {
          pthread_mutex_destroy(&(sopalin_data->tasktab_mutex[i]));
          pthread_cond_destroy(&(sopalin_data->tasktab_cond[i]));
          queueExit(&(sopalin_data->taskqueue[i]));
        }
      memFree_null(sopalin_data->taskqueue);

      if (sopalin_data->tasktab_mutex != NULL)
        memFree_null(sopalin_data->tasktab_mutex);
      if (sopalin_data->tasktab_cond != NULL)
        memFree_null(sopalin_data->tasktab_cond);
      if (sopalin_data->tasktab_indice != NULL)
        memFree_null(sopalin_data->tasktab_indice);
#endif

      if (sopalin_data->thread_data != NULL)
        {
          int threadnbr = SOLV_THRDNBR;
#ifdef WITH_STARPU
          if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
            threadnbr += sopalin_data->sopar->iparm[IPARM_CUDA_NBR];
          else
#endif
            threadnbr += sopalin_data->sopar->nbthrdcomm;

          for(i=0; i<threadnbr; i++)
            if (sopalin_data->thread_data[i] != NULL)
	      {
#if defined(PASTIX_DYNSCHED) && !defined(PASTIX_DYNSCHED_WITH_TREE)
                if (NO_ERR != tabtravel_deinit(sopalin_data->thread_data[i])) {
                  errorPrint("tabtravel_deinit (%s:%d)", __FILE__, __LINE__);
                }
#endif
		  memFree_null(sopalin_data->thread_data[i]);
	      }
          memFree_null(sopalin_data->thread_data);
        }

#ifdef STORAGE
      if (sopalin_data->mutex_flagtab != NULL)
        {
          for (i=0;i<UPDOWN_GCBLKNBR;i++)
            pthread_mutex_destroy(&(sopalin_data->mutex_flagtab[i]));
          memFree_null(sopalin_data->mutex_flagtab);
        }

      if (sopalin_data->grhs != NULL)
        memFree_null(sopalin_data->grhs);
      if (sopalin_data->flagtab != NULL)
        memFree_null(sopalin_data->flagtab);

      /* Allocate pthread_cond structures */
      if (sopalin_data->cond_flagtab != NULL)
        {
          for (i=0;i<UPDOWN_GCBLKNBR;i++)
            pthread_cond_destroy(&(sopalin_data->cond_flagtab[i]));
          memFree_null(sopalin_data->cond_flagtab);
        }

#endif
#if (!defined WITH_HWLOC && defined PASTIX_GET_SCHED_AFFINITY )
      memFree_null(sopalin_data->allowed_cpus);
#endif
      if (sopalin_data->common_flt != NULL)
        memFree_null(sopalin_data->common_flt);
      if (sopalin_data->common_dbl != NULL)
        memFree_null(sopalin_data->common_dbl);
      if (sopalin_data->ptr_csc != NULL)
        memFree_null(sopalin_data->ptr_csc);

#ifdef TRACE_SOPALIN
      trace_finish(sopalin_data->tracefile, SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, -1);
      memAllocUntrace();
      OUT_CLOSEFILEINDIR(sopalin_data->tracefile);
#endif
    }
}

/******************************************************************/
/*                                                                */
/*          Init and clean thread_data structure                  */
/*                                                                */
/******************************************************************/
#if (defined PASTIX_DYNSCHED && !(defined PASTIX_DYNSCHED_WITH_TREE))
static inline
int tabtravel_init(Sopalin_Data_t * sopalin_data,
                   Thread_Data_t  * thread_data,
                   int              me) {
  PASTIX_INT *visited = NULL;
  PASTIX_INT position = 0;
  PASTIX_INT father, son, i, j, bubnum = me;
  faststack_t stack;
  SolverMatrix  *datacode = sopalin_data->datacode;


  if (thread_data->tabtravel != NULL)
    return NO_ERR;

  MALLOC_INTERN(stack.tab,              datacode->btree->nodenbr+1, PASTIX_INT);
  MALLOC_INTERN(visited,                datacode->btree->nodenbr,   PASTIX_INT);
  MALLOC_INTERN(thread_data->tabtravel, datacode->thrdnbr,   PASTIX_INT);

  memset( thread_data->tabtravel, 0, datacode->thrdnbr * sizeof(PASTIX_INT) );

  FASTSTACK_INIT(stack);

  for(i=0; i<datacode->btree->nodenbr; i++) {
    stack.tab[i+1] = 0;
    visited[i] = 0;
  }

  thread_data->tabtravel[position] = bubnum;
  visited[bubnum] = 1;
  position++;

  do {
    /* Add father */
    father = BFATHER(datacode->btree, bubnum);
    if ( (father != -1)  &&
         (visited[father] == 0) ) {
      /*thread_data->tabtravel[position] = father;*/
      /*position++;*/
      visited[father] = 1;
      FASTSTACK_ADD(stack, father);
    }

    /* Add sons */
    for (j=0; j<BNBSON(datacode->btree, bubnum); j++) {
      son = BSON(datacode->btree, bubnum, j);
      if ( visited[son] == 0 ) {
        visited[son] = 1;

        if ( son < SOLV_THRDNBR ) {
          thread_data->tabtravel[position] = son;
          position++;
        }
        FASTSTACK_ADD(stack, son);
      }
    }

    FASTSTACK_TOP(stack, bubnum);

  } while ( bubnum != -1 );

  memFree_null(visited);
  memFree_null(stack.tab);

  return NO_ERR;
}

static inline
int tabtravel_deinit(Thread_Data_t * thread_data) {
  if (thread_data->tabtravel != NULL)
    memFree_null(thread_data->tabtravel);
  return NO_ERR;
}
#endif /* (PASTIX_DYNSCHED && !(defined PASTIX_DYNSCHED_WITH_TREE)) */

/*********************************/
/*
  Function: sopalin_init_smp

  Allocate and initialize thread_data
  This function is mono-thread and must be called by each thread.

  Parameters:
    sopalin_data - Sopalin_data structure
    me           - Thread indice
    fact         - Boolean for factorisation step or not

  Returns:
    void
 */
/*********************************/

void sopalin_init_smp(Sopalin_Data_t *sopalin_data, PASTIX_INT me, int fact, int init)
{
  Thread_Data_t *thread_data;
  SolverMatrix  *datacode = sopalin_data->datacode;
  SopalinParam  *sopar    = sopalin_data->sopar;
  PASTIX_INT            i;

  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_one("%s", OUT2_SOP_BINITL);

  thread_data = sopalin_data->thread_data[me];

  if (thread_data == NULL)
    {
      MALLOC_INTERN(thread_data, 1, Thread_Data_t);
#if defined(PASTIX_DYNSCHED) && !defined(PASTIX_DYNSCHED_WITH_TREE)
      thread_data->tabtravel = NULL;
#endif

      /* Initialisation des variables de la structure */
      thread_data->nbpivot   = 0;
      thread_data->flag_bind = sopalin_data->sopar->iparm[IPARM_BINDTHRD];
#ifdef TRYLOCK
      thread_data->ptbusy = 0;
      thread_data->ptfree = 0;
      thread_data->ptwait = 0;
#endif
#ifdef PASTIX_DYNSCHED
      thread_data->esp = 0;
#endif

      /* On associe les threads a un proc dans la version SMP */
#ifdef SMP_SOPALIN

#  if (defined X_ARCHalpha_compaq_osf1) || ( defined X_ARCHpower_ibm_aix && defined TEST_IRECV)
      thread_data->flag_bind = 0; /* Dans ces cas la on ne bind pas */
#  elif (defined X_ARCHpower_ibm_aix)
      if (THREAD_COMM_ON)
        thread_data->flag_bind = 0; /* Dans ces cas la on ne bind pas */
#  endif

#  ifdef STARPU_INIT_SMP
      if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
	thread_data->flag_bind = 0;
#  endif /* STARPU_INIT_SMP */
      {
        int nbproc;

#  ifdef MARCEL
        nbproc = marcel_nbvps();
#  else
        nbproc = sysconf(_SC_NPROCESSORS_ONLN);
#  endif
        if (INIT_COMPUTE & init)
          {
            if ((thread_data->flag_bind)
                && (SOLV_THRDNBR <= nbproc))
              {
                int cpu;

                if (thread_data->flag_bind == 1)
                  {
                    /* Calcul du proc sur lequel binder le thread */
#  if ( defined WITH_HWLOC || !defined PASTIX_GET_SCHED_AFFINITY )
		    cpu = (me + SOLV_PROCNUM * SOLV_THRDNBR)%nbproc;
#  else
		    cpu = sopalin_data->allowed_cpus[me];
#  endif
                  }
                else
                  {
                    if (sopalin_data->sopar->bindtab == NULL)
                      {
                        errorPrint("Bindtab not defined.");
                        EXIT(MOD_SOPALIN,BADPARAMETER_ERR);
                      }
                    cpu = sopalin_data->sopar->bindtab[me];
                  }

                cpu = sopalin_bindthread(cpu);

                if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
                  print_one("%s", OUT2_SOP_BIND);
                if (
#  ifdef WITH_STARPU
                  sopalin_data->sopar->iparm[IPARM_STARPU] == API_NO &&
#  endif /* WITH_STARPU */
                  sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO
                    )
                  {
                    int il;
                    static volatile int *tab = NULL;

                    MONOTHREAD_BEGIN;
                    MALLOC_INTERN(tab,  SOLV_THRDNBR, int);
                    for (il=0; il < SOLV_THRDNBR; il++)
                      tab[il] = 0;
                    MONOTHREAD_END;


                    SYNCHRO_X_THREAD(SOLV_THRDNBR, sopalin_data->barrier);
                    tab[me] = cpu;
                    SYNCHRO_X_THREAD(SOLV_THRDNBR, sopalin_data->barrier);

                    MONOTHREAD_BEGIN;
                    int *tab2 = NULL;
                    if (SOLV_PROCNUM == 0)
                      MALLOC_INTERN(tab2,  SOLV_THRDNBR*SOLV_PROCNBR, int);

                    MPI_Gather ( (void*)tab,  SOLV_THRDNBR, MPI_INT,
                                 tab2,        SOLV_THRDNBR, MPI_INT,
                                 0, PASTIX_COMM );

                    if (SOLV_PROCNUM == 0){
                      FILE *out;
                      int  jl;
                      OUT_OPENFILEINDIR(sopar->iparm, out, "threadbinding.txt", "w");
                      fprintf(out, "# Proc Thread Core\n");
                      for(il=0; il < SOLV_PROCNBR; il++)
                        for(jl=0; jl < SOLV_THRDNBR; jl++)
                          fprintf(out, "%ld %ld %ld\n", (long)il, (long)jl, (long)tab2[jl]);
                      OUT_CLOSEFILEINDIR(out);

                      memFree_null(tab2);
                    }
                    memFree_null(tab);
                    MONOTHREAD_END;
                  }
              }
            else
              {
                if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
                  print_one("%s", OUT2_SOP_NOTBIND);
              }
          }
      }
#endif /* SMP_SOPALIN */

#ifdef TRACE_SOPALIN
      thread_data->tracefile   = NULL;
      thread_data->traceactive = 0;
      thread_data->traceid     = 0;
#endif
    }


  thread_data->maxbloktab1 = NULL;
  thread_data->maxbloktab2 = NULL;
#ifdef COMPACT_SMX
  MALLOC_INTERN(thread_data->maxbloktab2, SOLV_COEFMAX, PASTIX_FLOAT);
#endif

  /* Données et buffers d'envoi */
  thread_data->gtabsize = NULL;
  thread_data->gtaboffs = NULL;
  thread_data->gtabtype = NULL;

  thread_data->send_block_requests     = NULL;
  thread_data->send_block_target       = NULL;
  thread_data->send_fanin_requests     = NULL;
  thread_data->send_fanin_target       = NULL;
  thread_data->send_fanin_target_extra = NULL;
#ifdef NO_MPI_TYPE
  thread_data->send_fanin_buffer      = NULL;
  thread_data->send_fanin_buffer_size = NULL;
  thread_data->send_block_buffer      = NULL;
  thread_data->send_block_buffer_size = NULL;
#else
  thread_data->send_fanin_mpitypes     = NULL;
  thread_data->send_fanin_infotab      = NULL;
#endif /* NO_MPI_TYPE */

  /* Données pour les réceptions */
  thread_data->recv_buffer = NULL;
#ifdef TEST_IRECV
  thread_data->recv_fanin_request = NULL;
  thread_data->recv_block_request = NULL;
  thread_data->recv_fanin_buffer  = NULL;
  thread_data->recv_block_buffer  = NULL;
#endif
  thread_data->maxsrequest_fanin = 0;
  thread_data->maxsrequest_block = 0;
#ifndef FORCE_NOMPI
  thread_data->srteststatus  = NULL;
  thread_data->srtestindices = NULL;
#endif


  /*
   * Donnees specifiques a la factorisation
   */
  if (fact)
    {

      if (INIT_COMPUTE & init)
        {
          /* work block initialization */
          MALLOC_INTERN(thread_data->maxbloktab1, SOLV_COEFMAX, PASTIX_FLOAT);
          memset(thread_data->maxbloktab1, 0, SOLV_COEFMAX*sizeof(PASTIX_FLOAT));
          MALLOC_INTERN(thread_data->maxbloktab2, SOLV_COEFMAX, PASTIX_FLOAT);
          memset(thread_data->maxbloktab2, 0, SOLV_COEFMAX*sizeof(PASTIX_FLOAT));
        }

      /* Donnes pour les envois */
      if ((INIT_SEND & init) && (SOLV_PROCNBR > 1))
        {
          MALLOC_INTERN(thread_data->send_block_requests,     MAX_S_REQUESTS, MPI_Request);
          MALLOC_INTERN(thread_data->send_block_target,       MAX_S_REQUESTS, PASTIX_INT);
          MALLOC_INTERN(thread_data->send_fanin_requests,     MAX_S_REQUESTS, MPI_Request);
          MALLOC_INTERN(thread_data->send_fanin_target,       MAX_S_REQUESTS, PASTIX_INT);
          MALLOC_INTERN(thread_data->send_fanin_target_extra, MAX_S_REQUESTS, PASTIX_INT*);
#if (!defined NO_MPI_TYPE)
          MALLOC_INTERN(thread_data->send_fanin_mpitypes,     MAX_S_REQUESTS,
                        MPI_Datatype);
          MALLOC_INTERN(thread_data->send_fanin_infotab,      MAX_S_REQUESTS, PASTIX_INT*);
#endif

          if (PACKMAX!=0)
            {
              MALLOC_INTERN(thread_data->gtabsize, 2*PACKMAX, int);
#ifndef NO_MPI_TYPE
              MALLOC_INTERN(thread_data->gtaboffs, 2*PACKMAX, MPI_Aint);
              MALLOC_INTERN(thread_data->gtabtype, 2*PACKMAX, MPI_Datatype);
#else
              MALLOC_INTERN(thread_data->gtaboffs, 2*PACKMAX, void *);
              MALLOC_INTERN(thread_data->gtabtype, 2*PACKMAX, int);

              MALLOC_INTERN(thread_data->send_fanin_buffer,      MAX_S_REQUESTS, void *);
              MALLOC_INTERN(thread_data->send_fanin_buffer_size, MAX_S_REQUESTS, PASTIX_INT);
              MALLOC_INTERN(thread_data->send_block_buffer,      MAX_S_REQUESTS, void *);
              MALLOC_INTERN(thread_data->send_block_buffer_size, MAX_S_REQUESTS, PASTIX_INT);
#endif /* NO_MPI_TYPE */
            }

          for (i=0;i<MAX_S_REQUESTS;i++)
            {
              thread_data->send_block_requests[i]     = MPI_REQUEST_NULL;
              thread_data->send_block_target[i]       = 0;
              thread_data->send_fanin_requests[i]     = MPI_REQUEST_NULL;
              thread_data->send_fanin_target[i]       = 0;
#if (!defined NO_MPI_TYPE)
              thread_data->send_fanin_infotab[i]      = NULL;
#endif
              thread_data->send_fanin_target_extra[i] = NULL;
            }
        }

      /* Donnees pour les receptions */
      if ((INIT_RECV & init) && (SOLV_PROCNBR > 1))
        {
#ifndef TEST_IRECV
          MALLOC_INTERN(thread_data->recv_buffer,
                        MAX(PACKMAX*(sizeof(PASTIX_INT)*MAXINFO)+
                            /* XL: apparently area is missing with dof
			     *     so I multiply by DOF...*/
			    sopalin_data->sopar->iparm[IPARM_DOF_NBR]*
			    PACKAREA*sizeof(PASTIX_FLOAT),
                            sizeof(PASTIX_INT)*(BTAGINFO+BCOFINFO)+
                            sizeof(PASTIX_FLOAT)*SOLV_BPFTMAX),
                        char);
#else

          if (MAX_R_REQUESTS)
            {
              PASTIX_INT sizef = PACKMAX*(sizeof(PASTIX_INT)*MAXINFO)+PACKAREA*sizeof(PASTIX_FLOAT);
              PASTIX_INT sizeb = sizeof(PASTIX_INT)*(BTAGINFO+BCOFINFO)+sizeof(PASTIX_FLOAT)*SOLV_BPFTMAX;

              MALLOC_INTERN(thread_data->recv_fanin_request, MAX_R_REQUESTS, MPI_Request);
              MALLOC_INTERN(thread_data->recv_block_request, MAX_R_REQUESTS, MPI_Request);
              MALLOC_INTERN(thread_data->recv_fanin_buffer,  MAX_R_REQUESTS, void*);
              MALLOC_INTERN(thread_data->recv_block_buffer,  MAX_R_REQUESTS, void*);

              for (i=0;i<MAX_R_REQUESTS;i++)
                {
                  MALLOC_INTERN(thread_data->recv_fanin_buffer[i], sizef, char);
                  MALLOC_INTERN(thread_data->recv_block_buffer[i], sizeb, char);
                }
            }
#endif /* TEST_IRECV  */
        }

      if (INIT_COMPUTE & init)
        {
#ifdef WITH_STARPU
	  if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_NO ||
	      me < SOLV_THRDNBR)
#endif
#if (defined PASTIX_DYNSCHED)
          {
            PASTIX_INT bubnum  = me;
            PASTIX_INT bubnum2 = me;

            while (bubnum != -1)
              {
    /* Allocation de la file de tâches */
                queueInit(&(sopalin_data->taskqueue[bubnum]), datacode->ttsknbr[bubnum]);
                {
                  for (i=0;i<datacode->ttsknbr[bubnum];i++)
                    {
                      PASTIX_INT task = datacode->ttsktab[bubnum][i];
                      if ((!TASK_CTRBCNT(task)) &&
                          ((TASK_TASKID(task) == COMP_1D) || (TASK_TASKID(task) == DIAG)))
                        {
#  if (DBG_PASTIX_DYNSCHED > 0)
                          ASSERTDBG(sopalin_data->taskmark[task] == -1, MOD_SOPALIN);
#  endif
                          sopalin_data->taskmark[task]++;
                          queueAdd(&(sopalin_data->taskqueue[bubnum]), task, ((double)TASK_PRIONUM(task)));
                        }
                      if (task > (SOLV_TASKNBR-1))
                        errorPrint("Pb task trop grand\n");
                    }
                }

                bubnum = BFATHER(datacode->btree, bubnum2);
                if ((bubnum != -1) &&
                    (datacode->btree->nodetab[bubnum].fcandnum != datacode->btree->nodetab[bubnum2].fcandnum))
                  bubnum = -1;
                bubnum2 = bubnum;
              }
#  if !defined(PASTIX_DYNSCHED_WITH_TREE)
            /*
             * Create the ordered list where to steal jobs
             */
            if (NO_ERR != tabtravel_init(sopalin_data, thread_data, me)) {
              errorPrint("tabtravel_init (%s:%d)", __FILE__, __LINE__);
            }
#  endif
          }
#endif /* PASTIX_DYNSCHED */

          /* Initialisation de la matrice dans le cas d'allocation locale a chaque thread */
#ifdef SMP_SOPALIN
          CoefMatrix_Allocate(sopar, datacode, &(sopalin_data->barrier.sync_lock), sopar->factotype, me);
#else
          CoefMatrix_Allocate(sopar, datacode, NULL, sopar->factotype, me);
#endif
          /* Need to synchronize thread to be sure allocation is done in 2D */
#ifdef STARPU_INIT_SMP
	  if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_NO)
#endif /* STARPU_INIT_SMP */
	    SYNCHRO_X_THREAD(SOLV_THRDNBR+((sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)?sopalin_data->sopar->iparm[IPARM_CUDA_NBR]:0), sopalin_data->barrier);
/* #define OOC_NOCOEFINIT */
#ifndef OOC_NOCOEFINIT
          CoefMatrix_Init(datacode, &(sopalin_data->barrier), me,
                          sopar->iparm, &(sopar->transcsc), sopalin_data);
#endif
          ooc_set_step(sopalin_data, API_TASK_NUMFACT);
        }
    }

  /*
   * Données spécifiques au solve (les buffers n'ont pas les memes dimensions)
   */
  else
    {
      int maxsrequest = MAX_S_REQUESTS;
      /* if (INIT_COMPUTE & init) */
      /*  { */
      /*    thread_data->maxsrequest_fanin = 1; */
      /*  } */
#if (defined PASTIX_DYNSCHED && !(defined PASTIX_DYNSCHED_WITH_TREE))
      if (INIT_COMPUTE & init) {
        if (NO_ERR != tabtravel_init(sopalin_data, thread_data, me)) {
          errorPrint("tabtravel_init (%s:%d)", __FILE__, __LINE__);
        }
      }
#endif /* (defined PASTIX_DYNSCHED && !(defined PASTIX_DYNSCHED_WITH_TREE)) */

      if ((INIT_SEND & init) && (SOLV_PROCNBR > 1))
        {
          MALLOC_INTERN(thread_data->send_fanin_requests, maxsrequest, MPI_Request);
          for (i=0;i<maxsrequest; i++)
            thread_data->send_fanin_requests[i] = MPI_REQUEST_NULL;
          MALLOC_INTERN(thread_data->gtabsize, UPDOWN_SM2XNBR+1, int);
#ifdef NO_MPI_TYPE
          MALLOC_INTERN(thread_data->send_fanin_buffer, maxsrequest, void*);
          thread_data->send_fanin_buffer[0] = NULL;

          MALLOC_INTERN(thread_data->gtaboffs, UPDOWN_SM2XNBR+1, void *);
          MALLOC_INTERN(thread_data->gtabtype, UPDOWN_SM2XNBR+1, int);
#else /* NO_MPI_TYPE */
          MALLOC_INTERN(thread_data->send_fanin_mpitypes, maxsrequest, MPI_Datatype);
          MALLOC_INTERN(thread_data->send_fanin_infotab, maxsrequest, PASTIX_INT*);
          memset(thread_data->send_fanin_infotab, 0, maxsrequest*sizeof(PASTIX_INT*));
          MALLOC_INTERN(thread_data->send_fanin_target, maxsrequest,      PASTIX_INT);
          MALLOC_INTERN(thread_data->gtaboffs,          UPDOWN_SM2XNBR+1, MPI_Aint);
          MALLOC_INTERN(thread_data->gtabtype,          UPDOWN_SM2XNBR+1, MPI_Datatype);
#endif /* NO_MPI_TYPE */
        }
    }

#ifdef TRACE_SOPALIN
  if ((thread_data->tracefile == NULL)
      && (SOLV_PROCNUM < 100 && SOLV_THRDNBR < 100))
    {
      char filename[18];
      char *mode;

      if (fact
          || thread_data->traceactive)
        {
          sprintf(filename,"traceFacto.%02d.%02d", (int)SOLV_PROCNUM, (int)me);
          filename[17] = '\0';
          if (thread_data->traceactive)
            {
              mode = "a";
            }
          else
            {
              thread_data->traceactive = 1;
              mode = "w";
            }
        }
      else
        {
          sprintf(filename,"traceSolve.%02d.%02d", (int)SOLV_PROCNUM, (int)me);
          filename[17] = '\0';
          if (thread_data->traceactive)
            {
              mode = "a";
            }
          else
            {
              thread_data->traceactive = 1;
              mode = "w";
            }
        }
      OUT_OPENFILEINDIR(sopar->iparm, thread_data->tracefile, filename, mode);
      trace_start(thread_data->tracefile, SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me);
    }
  else
    {
      errorPrint("Trace impossible car trop de proc ou de thread.");
      EXIT(MOD_SOPALIN, BADPARAMETER_ERR);
    }
#endif
#ifdef OOC
  ooc_defreeze(sopalin_data);
#endif

  sopalin_data->thread_data[me] = thread_data;

  return;
}


/*********************************/
/*
 * Function: sopalin_clean_smp
 *
 * Clean the thread_data structure.
 * This function is mono-thread and must be called just each thread.
 *
 * Parameters:
 *   sopalin_data - Sopalin_data structure
 *   thread_data  - thread_data structure to clean
 *
 * Returns:
 *   void
 */
/*********************************/
void sopalin_clean_smp(Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#ifdef TEST_IRECV
  PASTIX_INT i;
#endif

#ifdef OOC
  ooc_freeze(sopalin_data);
#endif
  if (thread_data->maxbloktab1 != NULL)
    memFree_null(thread_data->maxbloktab1);
  if (thread_data->maxbloktab2 != NULL)
    memFree_null(thread_data->maxbloktab2);

  if (thread_data->gtabsize != NULL)
    memFree_null(thread_data->gtabsize);
  if (thread_data->gtaboffs != NULL)
    memFree_null(thread_data->gtaboffs);
  if (thread_data->gtabtype != NULL)
    memFree_null(thread_data->gtabtype);

  if (SOLV_PROCNBR > 1)
    {
      memFree_null(thread_data->send_fanin_requests);
#ifdef NO_MPI_TYPE
      memFree_null(thread_data->send_fanin_buffer);
#else
      memFree_null(thread_data->send_fanin_mpitypes);
      memFree_null(thread_data->send_fanin_infotab);
#endif

      /* Uniquement en facto ou updo/sans NO_MPI_TYPE */
      if (thread_data->send_fanin_target != NULL)
  memFree_null(thread_data->send_fanin_target);

      /* Uniquement en facto */
      if (thread_data->send_fanin_target_extra != NULL)
  {
    memFree_null(thread_data->send_fanin_target_extra);
    memFree_null(thread_data->send_block_requests);
    memFree_null(thread_data->send_block_target);
#ifdef NO_MPI_TYPE
    memFree_null(thread_data->send_fanin_buffer_size);
    memFree_null(thread_data->send_block_buffer);
    memFree_null(thread_data->send_block_buffer_size);
#endif /* NO_MPI_TYPE */
  }
    }

  if (SOLV_PROCNBR > 1)
    {
#ifndef TEST_IRECV
      if (thread_data->recv_buffer != NULL)
  memFree_null(thread_data->recv_buffer);
#else
      if (MAX_R_REQUESTS
    && (thread_data->recv_fanin_buffer != NULL))
  {
    for (i=0;i<MAX_R_REQUESTS;i++){
      memFree_null(thread_data->recv_fanin_buffer[i]);
      memFree_null(thread_data->recv_block_buffer[i]);
    }
    memFree_null(thread_data->recv_fanin_buffer);
    memFree_null(thread_data->recv_block_buffer);
    memFree_null(thread_data->recv_fanin_request);
    memFree_null(thread_data->recv_block_request);
  }
#endif /* TEST_IRECV  */
    }

#ifndef FORCE_NOMPI
  if (thread_data->srteststatus != NULL)
    {
      memFree_null(thread_data->srteststatus);
      memFree_null(thread_data->srtestindices);
    }
#endif

#ifdef TRACE_SOPALIN
  if (thread_data->tracefile != NULL)
    {
      trace_finish(thread_data->tracefile, SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me);
      OUT_CLOSEFILEINDIR(thread_data->tracefile);
      thread_data->tracefile = NULL;
    }
#endif

  return;
}

/****************************************************************************/
/*                  RESTORE/BACKUP compteurs                                */
/*                  (Fonctions Mono-thread)                                 */
/****************************************************************************/
/*********************************/
/*
   Function: sopalin_backup

   Backup solver matrix counters.
   (Fonction Mono-thread)

   Parameters:
     datacode - SolverMatrix structure (common data)
     b        - Backup_t structure

   Returns:
     void

 */
/*********************************/
void sopalin_backup(SolverMatrix *datacode, Backup *b)
{
  PASTIX_INT i;

  print_debug(DBG_SOPALIN_ALLOC, "Backup Solver...\n");

  b->cpftmax = SOLV_CPFTMAX;
  b->arftmax = PACKAREA;
  b->nbftmax = PACKMAX;

  if (SOLV_TASKNBR)
    {
      MALLOC_INTERN(b->task_ctrbcnt, SOLV_TASKNBR, PASTIX_INT);
      MALLOC_INTERN(b->task_ftgtcnt, SOLV_TASKNBR, PASTIX_INT);
    }
  for (i=0;i<SOLV_TASKNBR;i++)
    {
      b->task_ctrbcnt[i] = TASK_CTRBCNT(i);
      b->task_ftgtcnt[i] = TASK_FTGTCNT(i);
    }

  if (SOLV_FTGTNBR)
    {
      MALLOC_INTERN(b->fanin_ctrbnbr, SOLV_FTGTNBR, PASTIX_INT);
      MALLOC_INTERN(b->fanin_prionum, SOLV_FTGTNBR, PASTIX_INT);
    }
  for (i=0;i<SOLV_FTGTNBR;i++)
    {
      b->fanin_ctrbnbr[i] = FANIN_CTRBNBR(i);
      b->fanin_prionum[i] = FANIN_PRIONUM(i);
    }

  if (SOLV_BTAGNBR)
    MALLOC_INTERN(b->btagtaskcnt, SOLV_BTAGNBR, PASTIX_INT);
  for (i=0;i<SOLV_BTAGNBR;i++)
    {
      b->btagtaskcnt[i] = BTAG_TASKCNT(i);
    }

  if (SOLV_BCOFNBR)
    MALLOC_INTERN(b->bcofsendcnt, SOLV_BCOFNBR, PASTIX_INT);
  for (i=0;i<SOLV_BCOFNBR;i++)
    {
      b->bcofsendcnt[i] = SOLV_BCOFTAB[i].sendcnt;
    }

  if (SYMB_BLOKNBR)
    MALLOC_INTERN(b->symbol_cblknum, SYMB_BLOKNBR, PASTIX_INT);
  for (i=0;i<SYMB_BLOKNBR;i++)
    b->symbol_cblknum[i] = SYMB_CBLKNUM(i);

  b->symbol_nodenbr = SYMB_NODENBR;
}

/*********************************/
/*
   Function: sopalin_restore

   Restore solver matrix counters.
   (Fonction Mono-thread)

   Parameters:
     datacode - SolverMatrix structure (common data)
     b        - Backup structure

   Returns
     void
 */
/*********************************/
void sopalin_restore(SolverMatrix *datacode, Backup *b)
{
  PASTIX_INT i;

  print_debug(DBG_SOPALIN_ALLOC, "Restore Solver...\n");

  SOLV_CPFTMAX = b->cpftmax;
  PACKAREA     = b->arftmax;
  PACKMAX      = b->nbftmax;

  for (i=0;i<SOLV_TASKNBR;i++)
    {
      TASK_CTRBCNT(i) = b->task_ctrbcnt[i];
      TASK_FTGTCNT(i) = b->task_ftgtcnt[i];
    }
  if (SOLV_TASKNBR)
    {
      memFree_null(b->task_ctrbcnt);
      memFree_null(b->task_ftgtcnt);
    }
  else
    {
      b->task_ctrbcnt = NULL;
      b->task_ftgtcnt = NULL;
    }
  for (i=0;i<SOLV_FTGTNBR;i++)
    {
      FANIN_CTRBNBR(i) = b->fanin_ctrbnbr[i];
      FANIN_CTRBCNT(i) = b->fanin_ctrbnbr[i];
      FANIN_PRIONUM(i) = b->fanin_prionum[i];
    }
  if (SOLV_FTGTNBR)
    {
      memFree_null(b->fanin_ctrbnbr);
      memFree_null(b->fanin_prionum);
    }
  else
    {
      b->fanin_ctrbnbr = NULL;
      b->fanin_prionum = NULL;
    }

  for (i=0;i<SOLV_BTAGNBR;i++)
    {
      BTAG_TASKCNT(i) = b->btagtaskcnt[i];
    }
  if (SOLV_BTAGNBR)
    memFree_null(b->btagtaskcnt);
  else
    b->btagtaskcnt = NULL;

  for (i=0;i<SOLV_BCOFNBR;i++)
    {
      SOLV_BCOFTAB[i].sendcnt = b->bcofsendcnt[i];
    }
  if (SOLV_BCOFNBR)
    memFree_null(b->bcofsendcnt);
  else
    b->bcofsendcnt = NULL;

  for (i=0;i<SYMB_BLOKNBR;i++)
    SYMB_CBLKNUM(i) = b->symbol_cblknum[i];
  if (SYMB_BLOKNBR)
    memFree_null(b->symbol_cblknum);
  else
    b->symbol_cblknum = NULL;

  SYMB_NODENBR = b->symbol_nodenbr;

  /* Restauration des pointeurs sur les btags pour step-by-step */
  {
    PASTIX_INT task;

    for (i=0; i<datacode->tasknbr; i++)
      datacode->tasktab[i].btagptr = NULL;
    for (i=0; i<datacode->btagnbr; i++)
      {
  if (datacode->proc2clust[datacode->btagtab[i].infotab[BTAG_PROCDST]] == datacode->clustnum)
    {
      task = datacode->btagtab[i].infotab[BTAG_TASKDST];
      datacode->tasktab[task].btagptr = &(datacode->btagtab[i]);

      while (datacode->tasktab[task].tasknext != datacode->btagtab[i].infotab[BTAG_TASKDST])
        {
    task = datacode->tasktab[task].tasknext;
    datacode->tasktab[task].btagptr = &(datacode->btagtab[i]);
        }
    }
      }
  }
}

/*********************************/
/*
   Function: solve_backup

   Backup solver matrix fanin contrib counters.
   (Fonction Mono-thread)

   Parameters:
     datacode - SolverMatrix structure (common data)
     b        - BackupSolve_t structure

   Returns:
     void

 */
/*********************************/
void solve_backup(SolverMatrix *datacode, BackupSolve_t *b)
{
  PASTIX_INT i;

  print_debug(DBG_UPDO,"Backup Solver...\n");

  if (SOLV_FTGTNBR)
    {
      MALLOC_INTERN(b->fanin_ctrbnbr, SOLV_FTGTNBR, PASTIX_INT);
    }
  for (i=0;i<SOLV_FTGTNBR;i++)
    b->fanin_ctrbnbr[i] = FANIN_CTRBNBR(i);

  MALLOC_INTERN(b->symbol_cblknum, SYMB_BLOKNBR, PASTIX_INT);
  for (i=0;i<SYMB_BLOKNBR;i++)
    b->symbol_cblknum[i] = SYMB_CBLKNUM(i);

}

/*********************************/
/*
   Function: solve_restore

   Restore solver matrix fanin contrib counters.
   (Fonction Mono-thread)

   Parameters:
     datacode - SolverMatrix structure (common data)
     b        - BackupSolve_t structure

   Returns
     void
 */
/*********************************/
void solve_restore(SolverMatrix *datacode, BackupSolve_t *b)
{
  PASTIX_INT i;

  print_debug(DBG_UPDO, "Restore Solver...\n");

  for (i=0;i<SOLV_FTGTNBR;i++)
    {
      FANIN_CTRBNBR(i) = b->fanin_ctrbnbr[i];
      FANIN_CTRBCNT(i) = b->fanin_ctrbnbr[i];
    }
  if (SOLV_FTGTNBR)
      memFree_null(b->fanin_ctrbnbr);

  for (i=0;i<SYMB_BLOKNBR;i++)
    SYMB_CBLKNUM(i) = b->symbol_cblknum[i];
  memFree_null(b->symbol_cblknum);

}
