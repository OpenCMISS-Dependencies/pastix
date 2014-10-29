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
   file: ooc.c
   Out Of core option for PaStiX

   author: Xavier Lacoste
   date: April 2008
*/

#if (defined OOC) || (defined OOC_FTGT)
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#include "common_pastix.h"
#include "out.h"
#include "ftgt.h"
#include "symbol.h"
#include "queue.h"
#include "csc.h"
#include "updown.h"
#include "bulles.h"
#include "solver.h"
#include "sopalin_define.h"
#include "sopalin_thread.h"
#include "sopalin3d.h"
#include "sopalin_acces.h"
#ifdef OOC_NOCOEFINIT
#include "csc_intern_solve.h"
#endif
/* #include "pastix.h" */
#include "ooc.h"

#ifndef OOC_DIR
#define OOC_DIR                "/tmp/pastix"
#endif

#define STR_SIZE               60
#define OOC_ME                 (me-SOLV_THRDNBR-COMM_THRNBR)
#define COMM_THRNBR            sopalin_data->sopar->nbthrdcomm
#define HITCBLKTAB(cblknum)    ooc->hitcblktab[cblknum].count
#define HACKEDCBLK(cblknum)    ooc->hackedcblk[cblknum].count
#define HASBEENSAVED(cblknum)  ooc->hasbeensaved[cblknum].count
#define HITFTGTTAB(ftgtnum)    ooc->hitftgttab[ftgtnum].count
#define HACKEDFTGT(ftgtnum)    ooc->hackedftgt[ftgtnum].count

#define MARGE_LOAD             2
#define HIGH                   0
#define MEDIUM                 1*MAX((SOLV_PROCNBR-1),1)
#define LOW                    1*MAX((SOLV_PROCNBR-1),1)

#ifdef DEBUG_OOC
#define CBLKWAIT(c) ooc->compteur_debug[c*4]++
#define CBLKSAVE(c) ooc->compteur_debug[c*4+1]++
#define CBLKREAD(c) ooc->compteur_debug[c*4+2]++
#define CBLKWRTE(c) ooc->compteur_debug[c*4+3]++
#else
#define CBLKWAIT(c)
#define CBLKSAVE(c)
#define CBLKREAD(c)
#define CBLKWRTE(c)
#endif

#define PRINT_STATS(string){                                            \
    if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT &&	\
        ooc->stop == 0)                                                 \
      {                                                                 \
        fprintf(stdout, OOC_IN_STEP, (int)SOLV_PROCNUM, string);        \
        fprintf(stdout, OOC_WRITTEN,                                    \
                (int)SOLV_PROCNUM,                                      \
                MEMORY_WRITE(ooc->ooc_thread[ooc_threadnum]->written),	\
                MEMORY_UNIT_WRITE(ooc->ooc_thread[ooc_threadnum]->written), \
                MEMORY_WRITE(ooc->allocated),                           \
                MEMORY_UNIT_WRITE(ooc->allocated));                     \
        fprintf(stdout,OOC_READ,                                        \
                (int)SOLV_PROCNUM,                                      \
                MEMORY_WRITE(ooc->ooc_thread[ooc_threadnum]->read),     \
                MEMORY_UNIT_WRITE(ooc->ooc_thread[ooc_threadnum]->read)); \
        fprintf(stdout,OOC_ALLOCATED,                                   \
                (int)SOLV_PROCNUM,                                      \
                MEMORY_WRITE(memAllocGetCurrent()),                     \
                MEMORY_UNIT_WRITE(memAllocGetCurrent()));               \
        fprintf(stdout,OOC_MAX_ALLOCATED,                               \
                (int)SOLV_PROCNUM,                                      \
                MEMORY_WRITE(memAllocGetMax()),                         \
                MEMORY_UNIT_WRITE(memAllocGetMax()));                   \
      }                                                                 \
  }

#ifdef OOC_PERCENT_COEFTAB

#define NOT_ENOUGH_MEM {                                                \
    min_mem += ooc->allocated;                                          \
    errorPrint("%s:%d: OOC limit (%d %%) too small : need greater than %d %%", \
               __FILE__,__LINE__,                                       \
               ooc->global_limit, min_mem*100/ooc->coeftabsize );       \
    EXIT(MOD_SOPALIN,OUTOFMEMORY_ERR);                                  \
  }
#else
#define NOT_ENOUGH_MEM {                                                \
    min_mem += memAllocGetCurrent();                                    \
    errorPrint("%s:%d: OOC limit (%lg %s) too small : need greater than %lg %s", \
               __FILE__,__LINE__,                                       \
               MEMORY_WRITE(ooc->global_limit),                         \
               MEMORY_UNIT_WRITE(ooc->global_limit),                    \
               MEMORY_WRITE(min_mem),                                   \
               MEMORY_UNIT_WRITE(min_mem));                             \
    EXIT(MOD_SOPALIN,OUTOFMEMORY_ERR);                                  \
  }

#endif

#define CHECK_HACK_QUEUE {                                              \
    PASTIX_INT cblk;                                                           \
    MUTEX_LOCK(&ooc->mutex_toload);                                     \
    while (queueSize(&ooc->toload) > 0)                                 \
      {                                                                 \
        cblk = queueGet(&ooc->toload);                                  \
        MUTEX_UNLOCK(&ooc->mutex_toload);                               \
        if (cblk >= 0)                                                  \
          ooc_prefetch(sopalin_data, cblk , ooc_threadnum, HIGH);       \
        else                                                            \
          ooc_prefetch_ftgt(sopalin_data, -(cblk+1),                    \
                            ooc_threadnum, HIGH);                       \
        MUTEX_LOCK(&ooc->mutex_toload);                                 \
      }                                                                 \
    MUTEX_UNLOCK(&ooc->mutex_toload);                                   \
  }

#define FACTO_FTGTSIZE(t)                                   \
  ((FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1)                    \
   * (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1)                  \
   * ((sopalin_data->sopar->factotype == API_FACT_LU)?2:1))

#define SOLVE_FTGTSIZE(t)                       \
  (UPDOWN_SM2XNBR                               \
   * (FANIN_LCOLNUM(t) - FANIN_FCOLNUM(t) + 1))

#define FTGTSIZE(t)                             \
  ((ooc->current_iterooc == 0) ?                \
   FACTO_FTGTSIZE(t) : SOLVE_FTGTSIZE(t))

#define CBLKSIZE(cblk)                              \
  (SOLV_STRIDE(cblk)                                \
   * (SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1))


#ifdef OOC_CLOCK
#define OOC_CLOCK_INIT  {clockInit(&(clk));clockStart(&(clk));}
#define OOC_CLOCK_STOP  {clockStop(&(clk));}
#define OOC_CLOCK_GET   ooc_clock_get(sopalin_data, &(clk), me);
#else
#define OOC_CLOCK_INIT
#define OOC_CLOCK_STOP
#define OOC_CLOCK_GET
#endif

#ifdef OOC_FTGT
#define MUTEX_LOCK_ALLOC    MUTEX_LOCK(&(ooc->mutex_allocated))
#define MUTEX_UNLOCK_ALLOC  MUTEX_UNLOCK(&(ooc->mutex_allocated))
#else
#define MUTEX_LOCK_ALLOC
#define MUTEX_UNLOCK_ALLOC
#endif


/*
   section: Out-Of-core structures
 */

/*
   struct: ooc_count_

   Structure containing a volatile integer, used for access count in ooc.
 */
struct ooc_count_ {
  /* variable: count
     volatile int to store shared counter */
  volatile int count;
};
/*
   struct: ooc_thread_

   ooc structure for each ooc thread
 */
struct ooc_thread_ {
  unsigned long         read;                 /*+ read total size read by thread                             +*/
  unsigned long         written;              /*+ total size written by thread                               +*/
#ifdef  OOC_CLOCK
  double                time_reading;
  double                time_writting;
#endif
};

/*
   struct: ooc_

   ooc structure
 */
struct ooc_ {
  PASTIX_INT                   id;                   /*+ id for saving files                                          +*/
  int                   stop;                 /*+ indicated when the thread should die                         +*/
  int                   threadnbr;            /*+ ooc thrdnbr                                                  +*/
  struct ooc_thread_ ** ooc_thread;           /*+ struct local for each thread                                 +*/
  unsigned long         global_limit;         /*+ Limit ask by user                                            +*/
  PASTIX_INT                   allocated;            /*+ OOC allocated memory                                         +*/
  pthread_mutex_t      *mutex_cblk;           /*+ Mutex for each cblk                                          +*/
  pthread_cond_t       *cond_cblk;            /*+ Cond for each cblk                                           +*/
  int                   max_cblk_size;        /*+ Max size of a cblk                                           +*/
  struct ooc_count_    *hitcblktab;           /*+ Number of access predicted for one thread and column block   +*/
  struct ooc_count_    *hackedcblk;           /*+ Number of hacked accessed except for receiving column block  +*/
  struct ooc_count_    *thrdUsingCblk;        /*+ Count for each cblk number of thread using it                +*/
  struct ooc_count_    *hitftgttab;           /*+ count for ftgt                                               +*/
  struct ooc_count_    *hackedftgt;           /*+ Number of hacked accessed for ftgt                           +*/
  struct ooc_count_    *hasbeensaved;         /*+ indicate if cblk has been saved after factorisation finished +*/
  PASTIX_INT                  *tasktab;              /*+ Ordered task tabular                                         +*/
  Queue                 queue;                /*+ queue of cblk to free                                        +*/
  pthread_mutex_t       mutex_queue;          /*+ mutex on the queue of cblk to free                           +*/
  pthread_cond_t        cond_queue;           /*+ con for the queue of cblk to free                            +*/
  Queue                 toload;               /*+ Queue of block column not taken into account in loops        +*/
  pthread_mutex_t       mutex_toload;         /*+ mutex on the queue of cblk to load                           +*/
  pthread_cond_t        cond_toload;          /*+ cond on the queue of cblk to load                            +*/
  PASTIX_INT                   frozen;               /*+ Used to stop loadings temporatly                             +*/
  pthread_mutex_t       mutex_frozen;         /*+ Mutex on frozen                                              +*/
  pthread_cond_t        cond_frozen;          /*+ Cond on frozen                                               +*/
  int                   current_step;         /*+ Current computing step for computing threads                 +*/
  int                   current_itercomp;
  int                   current_taskooc;
  pthread_mutex_t       mutex_current_task;
  pthread_cond_t        cond_current_task;
  int                   current_iterooc;
#ifdef OOC_DETECT_DEADLOCKS
  volatile int          threads_waiting;      /*+ Number of thread waiting for queue to be fill-in             +*/
  pthread_mutex_t       mutex_thrd_wting;     /*+ mutex corresponding to threads_waiting                       +*/
  PASTIX_INT                  *blocks_waited_for;    /*+ Blocks computing threads are waiting for                     +*/
  PASTIX_INT                   lastloaded;           /*+ Last block loaded                                            +*/
#endif
  int                   threads_receiving;    /*+ Number of thread waiting for a contribution                  +*/
  pthread_mutex_t       mutex_thrd_rcving;    /*+ mutex corresponding to threads_receiving                     +*/
  pthread_mutex_t       mutex_step;           /*+ mutex on the current_step variable                           +*/
  pthread_cond_t        cond_step;            /*+ cond for the current_step variable                           +*/
  PASTIX_INT                   cblkloaded;           /*+ Number of blocks loaded                                      +*/
  PASTIX_INT                   ftgtloaded;           /*+ Number of blocks loaded                                      +*/
#ifdef OOC_FTGT
  pthread_mutex_t       mutex_allocated;      /*+ mutex on allocation                                          +*/
  PASTIX_INT                  *fanin_ctrbcnt;        /*+ Ordered task tabular                                         +*/
#endif
#ifdef OOC_PERCENT_COEFTAB
  PASTIX_INT                   coeftabsize;          /*+ size of the coeftab                                          +*/
#endif
#ifdef OOC_CLOCK
  double                *time_waiting;        /*+ time spent waiting for a block                               +*/
#endif
#ifdef DEBUG_OOC
  int                   opened;               /*+ number of file opened                                        +*/
  pthread_mutex_t       opened_mutex;         /*+ mutex on opened                                              +*/
  int                  *compteur_debug;
#endif
};

/*
   section: Out-Of-core functions
 */
/* internal functions declarations */
struct ooc_thread_ * get_ooc_thread      (Sopalin_Data_t * sopalin_data, PASTIX_INT cblknum);
int                  ooc_init_thread     (Sopalin_Data_t * sopalin_data, int threadnum);
int                  ooc_clock_get       (Sopalin_Data_t * sopalin_data, Clock * clk, int computing_thread);
int                  ooc_do_save_coef    (Sopalin_Data_t * sopalin_data, PASTIX_INT cblknum, int me);
int                  ooc_do_load_coef    (Sopalin_Data_t * sopalin_data, PASTIX_INT cblknum, int me);
int                  ooc_allocate        (Sopalin_Data_t * sopalin_data, PASTIX_INT cblknum, int me);
#ifdef OOC_FTGT
int                  ooc_do_load_ftgt    (Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum, int me);
int                  ooc_remove_ftgt     (Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum, int me);
int                  ooc_allocate_ftgt   (Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum, int me);
int                  ooc_do_save_ftgt    (Sopalin_Data_t * sopalin_data, PASTIX_INT indbloc, int me);
#else
#define              ooc_do_load_ftgt(sopalin_data, ftgtnum, me) 1
#define              ooc_do_save_ftgt(sopalin_data, indbloc, me) 1
#define              ooc_remove_ftgt(sopalin_data, ftgtnum, me)
#define              ooc_allocate_ftgt(sopalin_data, ftgtnum, me);
#endif
void                 reduceMem           (Sopalin_Data_t * sopalin_data, PASTIX_INT size_wanted, int me);
void                 cblkNextAccess      (Sopalin_Data_t * sopalin_data, PASTIX_INT tasknum, PASTIX_INT cblknum,
            int computing_thread, PASTIX_INT *key1, PASTIX_INT *key2);
PASTIX_INT                  cblkAndContribSize  (Sopalin_Data_t * sopalin_data, PASTIX_INT task);
PASTIX_INT                  updo_init_mem       (Sopalin_Data_t * sopalin_data );
PASTIX_INT                  updo_init_smp_mem   (Sopalin_Data_t * sopalin_data );
PASTIX_INT                  sopalin_init_smp_mem(Sopalin_Data_t * sopalin_data);
PASTIX_INT                  raff_mem            (Sopalin_Data_t * sopalin_data);
int                  recursive_mkdir     (char * path, mode_t mode);


/*
   Function: ooc_prefetch

   Prefetch one column bloc

   Parameters:
     sopalin_data - sopalin_data structure (common data)
     cblknum      - column block number
     me           - ooc thread number
                    (should be on as only one is possible)
     priority     - priority of the load (HIGH or LOW)

   Returns:
     EXIT_SUCCESS - if the bloc has been successfully loaded
     EXIT_SUCCESS - if factorisation and/or solve is finished
                    (depending of what is asked)
     Exit program - if load wasn't succesfull.
*/
int ooc_prefetch(Sopalin_Data_t * sopalin_data, PASTIX_INT cblknum,
                 int me, int priority)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  PASTIX_INT            limitsze;
  int            ret;

  print_debug(DBG_OOC_TRACE_V2,"->ooc_prefetch %ld\n", (long)cblknum);
#ifdef OOC_DETECT_DEADLOCKS
  ooc->lastloaded = cblknum;
#endif
  limitsze = CBLKSIZE(cblknum)*sizeof(PASTIX_FLOAT) + (priority*ooc->max_cblk_size);

  pthread_mutex_lock(&(ooc->mutex_frozen));
  while ((ooc->frozen == SOLV_THRDNBR) && (ooc->stop != 1))
    {
      pthread_cond_wait(&(ooc->cond_frozen),
      &(ooc->mutex_frozen));
    }
  pthread_mutex_unlock(&(ooc->mutex_frozen));

  if (ooc->stop == 1)
    return EXIT_SUCCESS;

  /** Add one access to THE current CBLK**/
  MUTEX_LOCK(&(ooc->mutex_cblk[cblknum]));
  if (HACKEDCBLK(cblknum) == 0)
    HITCBLKTAB(cblknum)++;
  else
    HACKEDCBLK(cblknum)--;

  if ((SOLV_COEFTAB(cblknum) != NULL)
      || (HITCBLKTAB(cblknum) <= 0))
    {
      MUTEX_UNLOCK(&(ooc->mutex_cblk[cblknum]));
      return EXIT_SUCCESS;
    }

  while ((SOLV_COEFTAB(cblknum) == NULL)
         && (ooc->stop != 1)
#ifdef OOC_PERCENT_COEFTAB
         && (ooc->allocated +
             limitsze > ooc->global_limit * ooc->coeftabsize / 100)
#else
         && ((memAllocGetCurrent() + limitsze) > ooc->global_limit)
#endif
         )
    {

      MUTEX_UNLOCK(&(ooc->mutex_cblk[cblknum]));

      ASSERTDBG((sopalin_data->sopar->factotype != API_FACT_LU) ||
                (SOLV_UCOEFTAB(cblknum) == NULL),
                MOD_SOPALIN);

      reduceMem(sopalin_data, limitsze, me);

      MUTEX_LOCK(&(ooc->mutex_cblk[cblknum]));
    }


  /* If compute is ended */
  if (ooc->stop == 1)
    {
      MUTEX_UNLOCK(&(ooc->mutex_cblk[cblknum]));
      return EXIT_SUCCESS;
    }

  /* allocate the cblk and read it if file already exist*/
  ret = ooc_do_load_coef(sopalin_data, cblknum, me);
  if (ret != EXIT_SUCCESS)
    {
      errorPrint("Couldn't load cblk %d : %s\n", cblknum,
                 ((ret == EXIT_FAILURE_OUT_OF_MEMORY)?"out of memory":
                  ((ret == EXIT_FAILURE_CBLK_NOT_NULL)?"cblk not null":
                   ((ret == EXIT_FAILURE_FILE_OPENING)?"opening file":
                    "unknown error"))));
      EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

  MUTEX_UNLOCK(&(ooc->mutex_cblk[cblknum]));
  pthread_cond_broadcast(&(ooc->cond_cblk[cblknum]));

  print_debug(DBG_OOC_TRACE_V2,"<-ooc_prefetch %ld\n", (long)cblknum);

  return EXIT_SUCCESS;
}


/*
   Function: ooc_prefetch_ftgt

   Prefetch one column fanin target buffer

   Parameters:
     sopalin_data - sopalin_data structure (common data)
     ftgtnum      - fanin target block number
     me           - ooc thread number
                    (should be on as only one is possible)

   Returns:
     EXIT_SUCCESS - if the bloc has been successfully loaded
     EXIT_SUCCESS - if factorisation and/or solve is finished
                    (depending of what is asked)
     Exit program - if load wasn't succesfull.
*/
int ooc_prefetch_ftgt(Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum,
                      int me, int priority)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  PASTIX_INT            ftgtsize = FTGTSIZE(ftgtnum)*sizeof(PASTIX_FLOAT);
  PASTIX_INT            limitsze;
  int            ret;

  print_debug(DBG_OOC_TRACE_V2,"->ooc_prefetch_ftgt %ld\n", (long)ftgtnum);
  limitsze = ftgtsize + (priority*ooc->max_cblk_size);

  pthread_mutex_lock(&(ooc->mutex_frozen));
  while (ooc->frozen == SOLV_THRDNBR)
    {
      pthread_cond_wait(&(ooc->cond_frozen),
                        &(ooc->mutex_frozen));
    }
  pthread_mutex_unlock(&(ooc->mutex_frozen));

  /** Add one access to THE current CBLK**/
  MUTEX_LOCK(&(sopalin_data->mutex_fanin[ftgtnum]));

  if (HACKEDFTGT(ftgtnum) == 0)
    HITFTGTTAB(ftgtnum)++;
  else
    HACKEDFTGT(ftgtnum)--;

  if (FANIN_COEFTAB(ftgtnum) != NULL)
    {
      MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[ftgtnum]));
      return EXIT_SUCCESS;
    }

  /* if the bloc is NULL and the current access hasn't been done by
     computation thread */
  while ( (FANIN_COEFTAB(ftgtnum) == NULL)
          && (HITFTGTTAB(ftgtnum) > 0)
          && (ooc->stop != 1)
#ifdef OOC_PERCENT_COEFTAB
          && (ooc->allocated + limitsze >
              ooc->global_limit * ooc->coeftabsize / 100)
#else
          && (memAllocGetCurrent() + limitsze > ooc->global_limit)
#endif
          )
    {

      MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[ftgtnum]));

      reduceMem(sopalin_data, limitsze, me);

      MUTEX_LOCK(&(sopalin_data->mutex_fanin[ftgtnum]));

    }

  /* If compute is ended */
  if (ooc->stop == 1)
    {
      MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[ftgtnum]));
      return EXIT_SUCCESS;
    }

  /* allocate the cblk and read it if file already exist*/
  ret = ooc_do_load_ftgt(sopalin_data, ftgtnum, me);
  if (ret != EXIT_SUCCESS)
    {
      errorPrint("Couldn't load ftgt %d : %s\n",ftgtnum,
                 ((ret == EXIT_FAILURE_OUT_OF_MEMORY)?"out of memory":
                  ((ret == EXIT_FAILURE_CBLK_NOT_NULL)?"cblk not null":
                   ((ret == EXIT_FAILURE_FILE_OPENING)?"opening file":
                    "unknown error"))));
      EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }
  MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[ftgtnum]));
  pthread_cond_broadcast(&(sopalin_data->cond_fanin[ftgtnum]));

  print_debug(DBG_OOC_TRACE_V2,"<-ooc_prefetch %ld\n", (long)ftgtnum);

  return EXIT_SUCCESS;
}

/*
   Function: ooc_coefMatrixInit

   Follows the  CoefMatrix_Init algorithm
   to prefetch column blocks in correct order.

   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     ooc_threadnum - ooc thread number
                     (should be on as only one is possible)

   Returns:
     EXIT_SUCCESS

   See Also:
     <CoefMatrix_Init>
*/
int ooc_coefMatrixInit(Sopalin_Data_t * sopalin_data, int ooc_threadnum)
{
  SolverMatrix *datacode = sopalin_data->datacode;
  ooc_t        *ooc      = sopalin_data->ooc;
  PASTIX_INT           task;
  PASTIX_INT           cblknum;

  print_debug(DBG_OOC_TRACE_V1,"->ooc_coefMatrixInit\n");


  for(task = SOLV_TASKNBR-1; task > -1; task--)
    {
      cblknum = TASK_CBLKNUM(ooc->tasktab[task]);
      ooc_prefetch(sopalin_data, cblknum, ooc_threadnum, LOW);
    }

  PRINT_STATS("COEFINIT");
  ooc->ooc_thread[ooc_threadnum]->written = 0;
  ooc->ooc_thread[ooc_threadnum]->read    = 0;
  print_debug(DBG_OOC_TRACE_V1,"<-ooc_coefMatrixInit\n");
  return EXIT_SUCCESS;
}

/*
   Function: ooc_sopalin

   Follows the sopalin_smp algorithm
   to prefetch column blocks and
   fanin target buffers in correct order


   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     ooc_threadnum - ooc thread number
                     (should be on as only one is possible)

   Returns:
     EXIT_SUCCESS

   See Also:
     <sopalin_smp>
*/
int ooc_sopalin(Sopalin_Data_t * sopalin_data, int ooc_threadnum)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  PASTIX_INT i;
  PASTIX_INT task1;
  PASTIX_INT task;
  PASTIX_INT cblknum;
  PASTIX_INT bloknum;
  PASTIX_INT n;
  int ret = EXIT_SUCCESS;

  print_debug(DBG_OOC_TRACE_V1,"->ooc_sopalin\n");

  for (i=0;i<SOLV_FTGTNBR;i++)
    FANIN_COEFTAB(i) = NULL;

  for (i = 0; i < SOLV_TASKNBR; i++)
    {
      /* Hack to load non predictable blocks */
      CHECK_HACK_QUEUE;

      task1   = ooc->tasktab[i];
      MUTEX_LOCK(&(ooc->mutex_current_task));
      ooc->current_taskooc = task1;
      MUTEX_UNLOCK(&(ooc->mutex_current_task));
      pthread_cond_broadcast(&(ooc->cond_current_task));
      cblknum = TASK_CBLKNUM(task1);
      ooc_prefetch(sopalin_data, cblknum, ooc_threadnum, LOW);
      switch(TASK_TASKID(task1))
        {
        case COMP_1D :
          n = 0;
          for(bloknum=SYMB_BLOKNUM(cblknum)+1;
              bloknum < SYMB_BLOKNUM(cblknum+1);
              bloknum++)
            {
              /* Hack to load non predictable blocks */
              CHECK_HACK_QUEUE

                task = SOLV_INDTAB[TASK_INDNUM(task1)+(n)];
              /* Local task */
              if (task < 0)
                {
                  ooc_prefetch(sopalin_data, TASK_CBLKNUM(-task),
                               ooc_threadnum, LOW);
                  n += SYMB_BLOKNUM(cblknum+1) - bloknum;
                }
              /* Remote task => Fanin to send */
              else
                {
#ifdef OOC_FTGT
                  PASTIX_INT bloknum2;

      /* Fanin target */
                  for (bloknum2=bloknum;
                       bloknum2<SYMB_BLOKNUM(cblknum+1);
                       bloknum2++)
                    {
                      task = SOLV_INDTAB[TASK_INDNUM(task1)+(n++)];
                      /* ftgt */
                      if (task < SOLV_FTGTNBR)
                        {
                          ooc_prefetch_ftgt(sopalin_data, task,
                                            ooc_threadnum, LOW);
                          ooc->fanin_ctrbcnt[task]--;
                          if (ooc->fanin_ctrbcnt[task] == 0)
                            ooc_prefetch_ftgt(sopalin_data, task,
                                              ooc_threadnum, LOW);
                        }
                    }
#else
                  n += SYMB_BLOKNUM(cblknum+1) - bloknum;
#endif
                }
            }
    break;
        case E2 :

          /* Hack to load non predictable blocks */
          CHECK_HACK_QUEUE;

          task = SOLV_INDTAB[TASK_INDNUM(task1)];
          if (task < 0)
            {
              ooc_prefetch(sopalin_data,TASK_CBLKNUM(-task),ooc_threadnum,LOW);
            }
#ifdef OOC_FTGT
          else
            {
              /* ftgt */
              if (task < SOLV_FTGTNBR)
                {
                  ooc_prefetch_ftgt(sopalin_data, task, ooc_threadnum, LOW);
                }
            }
#endif
          break;
        case DIAG :
          break;
        case E1 :
          break;
        default:
          errorPrint("Taskid unknown for task %ld\n", (long)task1);
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }
    }

  PRINT_STATS("OOC_SOPALIN");
  if (ooc->ooc_thread[ooc_threadnum]->written == 0)
    ret = EXIT_SUCCESS_ALL_LOADED;

  ooc->ooc_thread[ooc_threadnum]->written = 0;
  ooc->ooc_thread[ooc_threadnum]->read    = 0;
  print_debug(DBG_OOC_TRACE_V1,"<-ooc_sopalin\n");
  return ret;
}

/*
   Function: ooc_updo

   Follows the up_down_smp algorithm
   to prefetch column blocks and
   fanin target buffers in correct order


   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     ooc_threadnum - ooc thread number
                     (should be on as only one is possible)

   Returns:
     EXIT_SUCCESS

   See Also:
     <up_down_smp>
*/
int ooc_updo(Sopalin_Data_t * sopalin_data, int ooc_threadnum)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  PASTIX_INT task;
  PASTIX_INT cblknum;
  int ret = EXIT_SUCCESS;
  PASTIX_INT i;

  print_debug(DBG_OOC_TRACE_V1,"->ooc_updo\n");

  /* SYMB_CBLKNUM(x)=indice de la fanin en negatif */
  if (ooc->current_iterooc == 1)
    {
      for (task=0;task<SOLV_TASKNBR;task++)
        {
          PASTIX_INT j,c,n,t;
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
    }

  MUTEX_LOCK(&(ooc->mutex_step));
  while (((ooc->current_iterooc - ooc->current_itercomp) > 1) &&
         (ooc->stop !=1))
    pthread_cond_wait(&(ooc->cond_step), &(ooc->mutex_step));
  MUTEX_UNLOCK(&(ooc->mutex_step));

  /* down */
  for (i = 0; i < SOLV_TASKNBR; i++)
    {
      /* Hack to load non predictable blocks */
      /* Utile pour les bloc non chargé de la facto */
      CHECK_HACK_QUEUE;

      task = ooc->tasktab[i];

      if ((TASK_TASKID(task)!=COMP_1D) &&
    (TASK_TASKID(task)!=DIAG)) continue;

      cblknum = TASK_CBLKNUM(task);

      ooc_prefetch(sopalin_data, cblknum, ooc_threadnum, LOW);

#ifdef OOC_FTGT
      {
        PASTIX_INT bloknum;
        for (bloknum=SYMB_BLOKNUM(cblknum)+1;
             bloknum<SYMB_BLOKNUM(cblknum+1);
             bloknum++)
          {
            if (SYMB_CBLKNUM(bloknum)<=0)
              {
                ooc_prefetch_ftgt(sopalin_data, SOLV_FTGTIND(bloknum),
                                  ooc_threadnum, LOW);
              }
          }
      }
#endif
    }

  MUTEX_LOCK(&(ooc->mutex_step));
  while (((ooc->current_iterooc - ooc->current_itercomp) > 0) &&
         (ooc->stop != 1))
    {
      COND_TIMEWAIT(&(ooc->cond_step), &(ooc->mutex_step));
      CHECK_HACK_QUEUE;
    }
  MUTEX_UNLOCK(&(ooc->mutex_step));

  /* up */
  for (i = SOLV_TASKNBR-1; i>-1; i--)
    {
      task = ooc->tasktab[i];

      if ((TASK_TASKID(task)!=COMP_1D) &&
          (TASK_TASKID(task)!=DIAG)) continue;

      cblknum = TASK_CBLKNUM(task);

      ooc_prefetch(sopalin_data, cblknum, ooc_threadnum, HIGH);
    }

  if (ooc->ooc_thread[ooc_threadnum]->written == 0)
    ret = EXIT_SUCCESS_ALL_LOADED;

  PRINT_STATS("OOC_UPDO");
  ooc->ooc_thread[ooc_threadnum]->written = 0;
  ooc->ooc_thread[ooc_threadnum]->read    = 0;
  print_debug(DBG_OOC_TRACE_V1,"<-ooc_updo\n");

  return ret;
}

/*
   Function: ooc_thread

   ooc thread main routine.
   The ooc thread will have to prefetch all column blocks and
   fanin target buffer in correct order to make computing threads
   abble to compute.

   Parameters:
     arg - address of one sopthread_data_t containing the
           thread number and the sopalin_data address.

   Returns:
     EXIT_SUCCESS

*/
void * ooc_thread(void * arg)
{
  sopthread_data_t   * argument     = (sopthread_data_t * ) arg;
  Sopalin_Data_t     * sopalin_data = (Sopalin_Data_t *)argument->data;
  int                  me           = argument->me;
#ifndef FORCE_NOSMP
  SolverMatrix       * datacode     = sopalin_data->datacode;
#endif
  ooc_t              * ooc          = sopalin_data->ooc;
  struct ooc_thread_ * ooc_thread;
#ifdef TRACE
  double             * memtrace;
  Queue                lastaccess;
#endif
  int                  stop;

  me = OOC_ME;
  print_debug(DBG_OOC_TRACE_V1,"[%d]->ooc_thread\n", (int)me);

  ooc_thread = ooc->ooc_thread[me];
  stop       = 0;

  /* start the threads */
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    {
#ifndef OOC_NOCOEFINIT
      ooc_coefMatrixInit(sopalin_data, me);
#endif
      ooc_sopalin(sopalin_data, me);
    }
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] >= API_TASK_SOLVE)
    {
      while (ooc->stop != 1)
        {
          ooc->current_iterooc++;
          ooc_updo(sopalin_data, me);
        }
    }

  print_debug(DBG_OOC_TRACE_V1,"[%d]<-ooc_thread\n",(int)me);

  return EXIT_SUCCESS;
}


/*
   Function: ooc_init

   init ooc structure

   Parameters:
     sopalin_data - sopalin_data structure (common data)
     limit        - asked memory limit

   Returns:
     EXIT_SUCCESS - if init successfull
     EXIT         - if not enough memory.

*/
int ooc_init(Sopalin_Data_t * sopalin_data, PASTIX_INT limit)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc;
  PASTIX_INT            i,j,me, cblknum, task;
  PASTIX_INT            total_size;
  PASTIX_INT            secondsize;
  char           str[STR_SIZE];
  Queue          taskqueue;

  print_debug(DBG_OOC_TRACE_V1,"->ooc_init\n");

  /* allocating ooc structure */
  MALLOC_INTERN(sopalin_data->ooc, 1, struct ooc_);

  ooc = sopalin_data->ooc;
  ooc->id                = -1;
  ooc->stop              = 0;
  ooc->threadnbr         = sopalin_data->sopar->iparm[IPARM_OOC_THREAD];
  ooc->ooc_thread        = 0;
  ooc->allocated         = 0;
  ooc->mutex_cblk        = NULL;
  ooc->cond_cblk         = NULL;
  ooc->max_cblk_size     = 0;
  ooc->hitcblktab        = NULL;
  ooc->hackedcblk        = NULL;
  ooc->thrdUsingCblk     = NULL;
  ooc->hitftgttab        = NULL;
  ooc->hackedftgt        = NULL;
  ooc->tasktab           = NULL;
  ooc->hasbeensaved      = NULL;
  ooc->current_itercomp  = 0;
  ooc->current_iterooc   = 0;
#ifdef OOC_DETECT_DEADLOCKS
  ooc->threads_waiting   = 0;
  ooc->lastloaded        = SYMB_CBLKNBR+1;
#endif
  ooc->threads_receiving = 0;
#ifdef OOC_CLOCK
  ooc->time_waiting      = 0;
#endif
#ifdef DEBUG_OOC
  ooc->opened            = 0;
  ooc->compteur_debug    = NULL;
#endif

  /* Allocation and fill in of the ordered task tab used for ooc prefetches */
  MALLOC_INTERN(ooc->tasktab, SOLV_TASKNBR, PASTIX_INT);
  /* fill-in */
  queueInit(&taskqueue, SOLV_TASKNBR);
  for (task = 0; task < SOLV_TASKNBR; task++)
    {
      queueAdd(&taskqueue, task, (double)TASK_PRIONUM(task));
    }
  i = 0;
  while (queueSize(&taskqueue) > 0)
    {
      ooc->tasktab[i] = queueGet(&taskqueue);
      i++;
    }
  queueExit(&taskqueue);

  /* First ooc step is coefinit */
  ooc->current_step = OOCSTEP_COEFINIT;
  pthread_mutex_init(&(ooc->mutex_step),NULL);
  pthread_cond_init(&(ooc->cond_step),NULL);

  /* Frozen is used to freeze ooc activities during
     Cleans and inits
  */
  ooc->frozen = 0;
  pthread_mutex_init(&ooc->mutex_frozen,NULL);
  pthread_cond_init(&ooc->cond_frozen,NULL);

  ooc->current_taskooc   = -1;
  pthread_mutex_init(&ooc->mutex_current_task,NULL);
  pthread_cond_init(&ooc->cond_current_task,NULL);

  /* Number of bloc loaded,
     used to decide if there is not enopugh memory
  */
  ooc->cblkloaded = 0;
  ooc->ftgtloaded = 0;

  /* Mutex on threads_waiting */
#ifdef OOC_DETECT_DEADLOCKS
  pthread_mutex_init(&(ooc->mutex_thrd_wting), NULL);
  MALLOC_INTERN(ooc->blocks_waited_for, SOLV_THRDNBR, PASTIX_INT);
  for (i = 0; i < SOLV_THRDNBR; i++)
    ooc->blocks_waited_for[i] = SYMB_CBLKNBR+1;
#endif
  pthread_mutex_init(&(ooc->mutex_thrd_rcving), NULL);

  /* set id for ooc session
     And creates ooc directories
  */
  if (sopalin_data->sopar->iparm[IPARM_OOC_ID] == -1)
    {
      srand(time(0) + getpid());
      ooc->id = rand();
      sopalin_data->sopar->iparm[IPARM_OOC_ID] = ooc->id;

      /* creating ooc directories */
      sprintf(str,"%s/pastix_coef_%d",OOC_DIR,(int)ooc->id);
      if (EXIT_SUCCESS != recursive_mkdir(str,
                                          S_IRUSR | S_IWUSR | S_IRGRP |
                                          S_IXUSR | S_IWGRP | S_IROTH |
                                          S_IWOTH))
        {
          perror("mkdir 1");
          EXIT(MOD_SOPALIN,UNKNOWN_ERR);
        }
      if (sopalin_data->sopar->iparm[IPARM_FACTORIZATION] == API_FACT_LU)
        {
          sprintf(str,"%s/pastix_ucoef_%d",OOC_DIR,(int)ooc->id);
          if (-1 == mkdir(str,
                          S_IRUSR | S_IWUSR  | S_IXUSR|
                          S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH))
            {
              perror("mkdir 2");
              EXIT(MOD_SOPALIN,UNKNOWN_ERR);
            }
        }
#ifdef OOC_FTGT
      sprintf(str,"%s/pastix_ftgt_%d",OOC_DIR,(int)ooc->id);
      if (-1 == mkdir(str,
                      S_IRUSR | S_IWUSR | S_IRGRP |
                      S_IXUSR | S_IWGRP | S_IROTH | S_IWOTH))
        {
          perror("mkdir 3");
          EXIT(MOD_SOPALIN,UNKNOWN_ERR);
        }
#endif
    }
  else
    {
      ooc->id = sopalin_data->sopar->iparm[IPARM_OOC_ID];
    }

  /* Init queue for not predictable blocs */
  queueInit(&ooc->toload, SYMB_CBLKNBR/10);
  pthread_mutex_init(&ooc->mutex_toload,NULL);
  pthread_cond_init(&ooc->cond_toload,NULL);

  /* set globalk limit */
  ooc->global_limit = limit;
#ifndef OOC_PERCENT_COEFTAB
  ooc->global_limit <<= 20;
#endif
  /* Mutex on hicblktab and thrUsingCblk */
  MALLOC_INTERN(ooc->mutex_cblk, SYMB_CBLKNBR, pthread_mutex_t);
  MALLOC_INTERN(ooc->cond_cblk,  SYMB_CBLKNBR, pthread_cond_t);
  for (i = 0; i < SYMB_CBLKNBR; i++)
    {
      pthread_mutex_init(&(ooc->mutex_cblk[i]), NULL);
      pthread_cond_init(&(ooc->cond_cblk[i]), NULL);
    }

  /* Saving Queue */
  pthread_mutex_init(&ooc->mutex_queue,NULL);
  pthread_cond_init (&ooc->cond_queue,NULL);
  queueInit(&ooc->queue, SYMB_CBLKNBR);


  /* Clock */
#ifdef OOC_CLOCK
  MALLOC_INTERN(ooc->time_waiting, SOLV_THRDNBR, double);
  for (i = 0; i < SOLV_THRDNBR; i++)
    ooc->time_waiting[i] = 0;
#endif

  /*set max cblk size */
  total_size = 0;
  {
    PASTIX_INT maxsize = 0;
    secondsize = 0;
    for (cblknum = 0 ; cblknum < SYMB_CBLKNBR; cblknum++)
    {
      PASTIX_INT cblksize = CBLKSIZE(cblknum);
      if ( cblksize > maxsize )
        {
          secondsize = maxsize;
          maxsize    = cblksize;
        }
      total_size += cblksize;
    }
    total_size *= sizeof(PASTIX_FLOAT);
    ooc->max_cblk_size = maxsize * sizeof(PASTIX_FLOAT);
  }

  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    fprintf(stdout, "Total coeftab size : %.3g %s\n",
            MEMORY_WRITE(total_size), MEMORY_UNIT_WRITE(total_size));

#ifdef OOC_PERCENT_COEFTAB
  ooc->coeftabsize = total_size;
#ifdef OOC_FTGT
  for (i = 0; i < SOLV_FTGTNBR; i ++)
    ooc->coeftabsize += FACTO_FTGTSIZE(i)*sizeof(PASTIX_FLOAT);
#endif
#endif


  if (sopalin_data->sopar->factotype == API_FACT_LU)
    {
      ooc->max_cblk_size *= 2;
      secondsize         *= 2;
    }

  /* hit tabs */
  MALLOC_INTERN(ooc->hitcblktab,    SYMB_CBLKNBR, struct ooc_count_);
  MALLOC_INTERN(ooc->hackedcblk,    SYMB_CBLKNBR, struct ooc_count_);
  MALLOC_INTERN(ooc->thrdUsingCblk, SYMB_CBLKNBR, struct ooc_count_);
  MALLOC_INTERN(ooc->hitftgttab,    SOLV_FTGTNBR, struct ooc_count_);
  MALLOC_INTERN(ooc->hackedftgt,    SOLV_FTGTNBR, struct ooc_count_);
  MALLOC_INTERN(ooc->hasbeensaved,  SYMB_CBLKNBR, struct ooc_count_);

  for (j = 0; j < SYMB_CBLKNBR; j++)
    {
      ooc->thrdUsingCblk[j].count = 0;
      ooc->hitcblktab[j].count    = 0;
      ooc->hackedcblk[j].count    = 0;
      ooc->hasbeensaved[j].count  = 0;
    }

  for (j = 0; j < SOLV_FTGTNBR; j++)
    {
      ooc->hitftgttab[j].count = 0;
      ooc->hackedftgt[j].count = 0;
    }

  /* structure ooc_thread */
  MALLOC_INTERN(ooc->ooc_thread, ooc->threadnbr, struct ooc_thread_*);

  for (i = 0; i < ooc->threadnbr; i++)
    ooc->ooc_thread[i] = NULL;
  for ( me = SOLV_THRDNBR+COMM_THRNBR;
        me < SOLV_THRDNBR+COMM_THRNBR+ooc->threadnbr;
        me++ )
    ooc_init_thread(sopalin_data,me);

#ifdef OOC_FTGT
  pthread_mutex_init(&(ooc->mutex_allocated), NULL);

  /* Copy of fanin ctrbcnt */
  MALLOC_INTERN(ooc->fanin_ctrbcnt, SOLV_FTGTNBR, PASTIX_INT);
  for (i=0;i<SOLV_FTGTNBR;i++)
    ooc->fanin_ctrbcnt[i] = FANIN_CTRBCNT(i);
#endif

  /* Debug Part */
#ifdef DEBUG_OOC
  pthread_mutex_init(&ooc->opened_mutex, NULL);

  MALLOC_INTERN(ooc->compteur_debug, 4*SYMB_CBLKNBR, int);
  for (i = 0; i < 4*SYMB_CBLKNBR; i++)
    ooc->compteur_debug[i] = 0;

#endif

  /* Check if we have enough memory */
  if (
#ifdef OOC_PERCENT_COEFTAB
      ooc->global_limit*ooc->coeftabsize / 100 <
      (ooc->max_cblk_size + (LOW)*secondsize) + secondsize
#else
      ooc->global_limit  < memAllocGetCurrent() +
      (ooc->max_cblk_size + (LOW)*secondsize) + secondsize
#endif
      )
    {
      PASTIX_INT min_mem = (ooc->max_cblk_size + (LOW)*secondsize) *
        SOLV_THRDNBR + secondsize;

      NOT_ENOUGH_MEM;
    }
  else
    {
      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        {
#ifdef OOC_PERCENT_COEFTAB
          double limit_in_octet =  (double)ooc->global_limit *
            (double)ooc->coeftabsize /100.0;
          fprintf(stdout,OOC_MEM_LIM_PERCENT,
                  (int)ooc->global_limit,
                  (double)MEMORY_WRITE(limit_in_octet),
                  MEMORY_UNIT_WRITE(limit_in_octet));
#else
          fprintf(stdout,OOC_MEM_LIM,
                  (double)MEMORY_WRITE(ooc->global_limit),
                  MEMORY_UNIT_WRITE(ooc->global_limit));

#endif
        }
    }

  print_debug(DBG_OOC_TRACE_V1,"<-ooc_init\n");


  return EXIT_SUCCESS;
}

/*
   Function: ooc_exit

   clean ooc structure

   Parameters:
     sopalin_data - sopalin_data structure (common data)

   Returns:
     EXIT_SUCCESS
*/
int ooc_exit (Sopalin_Data_t * sopalin_data)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  PASTIX_INT            i;
  int            me;

  print_debug(DBG_OOC_TRACE_V1,"->ooc_exit\n");

  for (i = 0; i < SYMB_CBLKNBR; i++)
    ASSERTDBG((HITCBLKTAB(i) == 0) || (HITCBLKTAB(i) == 1), MOD_SOPALIN);

#ifdef DEBUG_OOC
  for (i = 0; i < SYMB_CBLKNBR; i++)
    fprintf(stderr, "%4d %4d %4d %4d %4d\n", i, CBLKWAIT(i), CBLKSAVE(i),
            CBLKREAD(i), CBLKWRTE(i));
  memFree_null(ooc->compteur_debug);

  pthread_mutex_destroy(&ooc->opened_mutex);
#endif

#ifdef OOC_DETECT_DEADLOCKS
  pthread_mutex_destroy(&(ooc->mutex_thrd_wting));
  memFree_null(ooc->blocks_waited_for);
#endif
  pthread_mutex_destroy(&(ooc->mutex_thrd_rcving));
#ifdef OOC_FTGT
  pthread_mutex_destroy(&(ooc->mutex_allocated));
  memFree_null(ooc->fanin_ctrbcnt);
#endif

#ifdef OOC_CLOCK
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    {
      for (i = 0; i < SOLV_THRDNBR; i++)
        {
          fprintf(stdout,"Total time waiting     : %.3g s\n",
                  ooc->time_waiting[i]);
        }
      for (i = 0; i < ooc->threadnbr; i++)
        {
          fprintf(stdout,"Total time writting    : %.3g s\n",
                  ooc->ooc_thread[i]->time_writting);
    fprintf(stdout,"Total time reading     : %.3g s\n",
            ooc->ooc_thread[i]->time_writting);
        }
    }
  memFree_null(ooc->time_waiting);
#endif

  for ( me = SOLV_THRDNBR+COMM_THRNBR;
        me < SOLV_THRDNBR+COMM_THRNBR +ooc->threadnbr;
        me++)
    memFree_null(ooc->ooc_thread[OOC_ME]);

  memFree_null(ooc->tasktab);

  queueExit(&ooc->queue);
  queueExit(&ooc->toload);
  pthread_mutex_destroy(&ooc->mutex_queue);
  pthread_cond_destroy(&ooc->cond_queue);
  pthread_mutex_destroy(&ooc->mutex_toload);
  pthread_cond_destroy(&ooc->cond_toload);
  pthread_mutex_destroy(&ooc->mutex_frozen);
  pthread_cond_destroy(&ooc->cond_frozen);

  for (i = 0; i < SYMB_CBLKNBR; i++)
    {
      pthread_mutex_destroy(&(ooc->mutex_cblk[i]));
      pthread_cond_destroy(&(ooc->cond_cblk[i]));
    }
  memFree_null(ooc->mutex_cblk);
  memFree_null(ooc->cond_cblk);
  memFree_null(ooc->thrdUsingCblk);
  memFree_null(ooc->hitcblktab);
  memFree_null(ooc->hitftgttab);
  memFree_null(ooc->hackedcblk);
  memFree_null(ooc->hackedftgt);
  memFree_null(ooc->hasbeensaved);

  if (NULL != ooc->ooc_thread)
    memFree_null(ooc->ooc_thread);

  if (NULL != ooc)
    memFree_null(ooc);


  print_debug(DBG_OOC_TRACE_V1,"<-ooc_exit\n");

  return EXIT_SUCCESS;
}

/*
  Function: ooc_init_thread

  init ooc thread structure

  Parameters:
     sopalin_data - sopalin_data structure (common data)
     me           - ooc thread number
                    (should be on as only one is possible)

   Returns:
     EXIT_SUCCESS
*/
int ooc_init_thread (Sopalin_Data_t * sopalin_data, int me)
{
#ifndef FORCE_NOSMP
  SolverMatrix       * datacode = sopalin_data->datacode;
#endif
  ooc_t              * ooc      = sopalin_data->ooc;
  struct ooc_thread_ * ooc_thread ;

  me = OOC_ME;
  print_debug(DBG_OOC_TRACE_V1,"[%d]->ooc_init_thread\n",(int)me);

  MALLOC_INTERN(ooc->ooc_thread[me], 1, struct ooc_thread_);

  ooc_thread = ooc->ooc_thread[me];
  ooc->ooc_thread[me]->read    = 0;
  ooc->ooc_thread[me]->written = 0;
#ifdef OOC_CLOCK
  ooc->ooc_thread[me]->time_writting = 0.;
  ooc->ooc_thread[me]->time_reading  = 0.;
#endif

  print_debug(DBG_OOC_TRACE_V1,"[%d]<-ooc_init_thread\n",(int)me);

  return EXIT_SUCCESS;
}

/*
   Function: ooc_save_coef

   Had column block to the save queue.
   Called by computing thread when
   work on a cblknum has been done.

   Parameters:
     sopalin_data - sopalin_data structure (common data)
     cblknum      - column block number.
     task         - task the computing thread is computing.
     step         - computation step (sopalin, up-down...).
     me           - computing thread number.
   Returns:
     EXIT_SUCCESS
 */
int ooc_save_coef(Sopalin_Data_t * sopalin_data, PASTIX_INT tasknum,
                  PASTIX_INT cblknum, int me)
{
  ooc_t        * ooc      = sopalin_data->ooc;
  PASTIX_INT key1;
  PASTIX_INT key2      = 0;
  int save      = 1;

  print_debug(DBG_OOC_TRACE_V2,"[%d]->ooc_save_coef %d\n",
        (int) me, (int)cblknum);

  MUTEX_LOCK(&(ooc->mutex_cblk[cblknum]));
  CBLKSAVE(cblknum);

  HITCBLKTAB(cblknum) --;
  ooc->thrdUsingCblk[cblknum].count--;

  if (ooc->thrdUsingCblk[cblknum].count < 0)
    {
      fprintf(stdout,"Et merde !\n");
      while(1);
    }

  if ((HITCBLKTAB(cblknum) > 0) || (ooc->thrdUsingCblk[cblknum].count > 0))
    save = 0;

  MUTEX_UNLOCK(&(ooc->mutex_cblk[cblknum]));

  if (save)
    {
      cblkNextAccess(sopalin_data, tasknum, cblknum, me, &key1, &key2);

      print_debug(DBG_OOC_PREDICTION, "NEXT %4ld %4ld %4ld %4ld\n",
                  (long)tasknum, (long)cblknum, (long)key1, (long)key2);

      MUTEX_LOCK(&(ooc->mutex_queue));
      queueAdd2(&(ooc->queue), cblknum, (double)(-key1), (-key2));
      MUTEX_UNLOCK(&(ooc->mutex_queue));
      pthread_cond_signal(&(ooc->cond_queue));
    }
  else
    print_debug(DBG_OOC_PREDICTION, "NEXT %4ld %4ld NOTSAVE\n",
                (long)tasknum, (long)cblknum);


  print_debug(DBG_OOC_TRACE_V2,"[%d]<-ooc_save_coef %d\n",
        (int) me, (int)cblknum);

  return EXIT_SUCCESS;
}

/*
   Function: ooc_do_save_coef

   Save one block column on disk.
   Call to function must be protected by mutex_cblk[cblknum].

   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     cblknum       - column block number
     ooc_threadnum - ooc thread number
                     (should be on as only one is possible)
   Returns:
     EXIT_SUCCESS                    - saved
     EXIT_FAILURE_SAVING_NULL_BUFFER - buffer is null
     EXIT_FAILURE_FILE_OPENING       - error opening file
     EXIT_FAILURE                    - error writing
 */
int ooc_do_save_coef (Sopalin_Data_t * sopalin_data, PASTIX_INT cblknum,
                      int ooc_threadnum)
{
  SolverMatrix       * datacode   = sopalin_data->datacode;
  ooc_t              * ooc        = sopalin_data->ooc;
  char                 str[STR_SIZE];
  int                  fd;
  PASTIX_INT                  coefnbr;
  PASTIX_INT                  ret;
  PASTIX_INT                  written;
  PASTIX_INT                  towrite;
  PASTIX_INT                  saved;
#ifdef OOC_CLOCK
  Clock                 clk;
#endif

  print_debug(DBG_OOC_TRACE_V2,
              "->ooc_do_save_coef %d (%d blocks %d allocated)\n",
              (int)cblknum, (int)ooc->cblkloaded, (int)ooc->allocated);

  saved = HASBEENSAVED(cblknum);
  if (SOLV_COEFTAB(cblknum) == NULL)
    {
      return EXIT_FAILURE_SAVING_NULL_BUFFER;
    }

  coefnbr = CBLKSIZE(cblknum);
  towrite = coefnbr*sizeof(PASTIX_FLOAT);
  /* in updown, no need to write coef on disk */
  if (!(ooc->current_step > OOCSTEP_SOPALIN) ||
      !(HASBEENSAVED(cblknum) == 1))
    {
      if ((ooc->current_step > OOCSTEP_SOPALIN) ||
          ((ooc->current_step == OOCSTEP_SOPALIN) && (saved == -1)))
        saved = 1;
      /* opening the file */
      sprintf(str, "%s/pastix_coef_%d/%d", OOC_DIR, (int)ooc->id, (int)cblknum);
      if (-1 == (fd = open(str,O_WRONLY|O_CREAT|O_SYNC,S_IRUSR |
                           S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH)))
        {
          errorPrint("Error opening file %s for write", str);
          perror("open (write)");
          exit(-1);
          return EXIT_FAILURE_FILE_OPENING;
        }

#ifdef DEBUG_OOC
      MUTEX_LOCK(&ooc->opened_mutex);
      ooc->opened++;
      MUTEX_UNLOCK(&ooc->opened_mutex);
#endif

      /* only save if not in use */
      /* A voir pour mettre des assert ou tester la valeur de retour */
      if (ooc->thrdUsingCblk[cblknum].count > 0)
        {
          return EXIT_FAILURE_CBLK_USED;
        }

      CBLKWRTE(cblknum);

      /* writing the file */
      written = 0;

      OOC_CLOCK_INIT;
      while ( written < towrite )
        {
          if (0 == (ret = write(fd,SOLV_COEFTAB(cblknum)+written,
                                towrite-written)))
            {
              errorPrint("Coulnd't write to disk");
              return EXIT_FAILURE;
            }
          written += ret;
        }
      OOC_CLOCK_STOP;

#ifdef OOC_CLOCK
      ooc->ooc_thread[ooc_threadnum]->time_writting += clockVal(&clk);
#endif

      if (EXIT_SUCCESS != close(fd))
        perror("close");

      ooc->ooc_thread[ooc_threadnum]->written += towrite;

#ifdef DEBUG_OOC
      MUTEX_LOCK(&ooc->opened_mutex);
      ooc->opened--;
      MUTEX_UNLOCK(&ooc->opened_mutex);
#endif
    }
  /* freeing saved coeftab */
  MUTEX_LOCK_ALLOC;
  memFree_null(SOLV_COEFTAB(cblknum));
  ooc->allocated -= towrite;
  MUTEX_UNLOCK_ALLOC;

  /* LU , work on ucoeftab*/
  if (sopalin_data->sopar->factotype == API_FACT_LU)
    {
      if (SOLV_UCOEFTAB(cblknum) == NULL)
        {
          errorPrint("Couldn't save non allocated buffer, %d",cblknum);
          return EXIT_FAILURE_SAVING_NULL_BUFFER;
        }
      if (!(ooc->current_step > OOCSTEP_SOPALIN) ||
          !(HASBEENSAVED(cblknum) == 1))
        {
          if ((ooc->current_step > OOCSTEP_SOPALIN) ||
              ((ooc->current_step == OOCSTEP_SOPALIN) && (saved == -1)))
            saved = 1;

    /* opening the file */
          sprintf(str,"%s/pastix_ucoef_%d/%d", OOC_DIR,
                  (int)ooc->id, (int)cblknum);

          if (-1 == (fd = open(str,O_WRONLY|O_CREAT|O_SYNC,
                               S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP |
                               S_IROTH | S_IWOTH)))
            {
              errorPrint("Error opening file %s for write", str);
              perror("open (write)");
              exit(-1);
              return EXIT_FAILURE_FILE_OPENING;
            }

#ifdef DEBUG_OOC
          MUTEX_LOCK(&ooc->opened_mutex);
          ooc->opened++;
          MUTEX_UNLOCK(&ooc->opened_mutex);
#endif

    /* writing the file */
    written = 0;
    OOC_CLOCK_INIT;
    while ( written < towrite )
      {
        if (0 == (ret = write(fd, SOLV_UCOEFTAB(cblknum)+written,
                              towrite-written)))
          {
            errorPrint("Coulnd't write to disk");
            return EXIT_FAILURE;
          }
        written += ret;
      }
    ooc->ooc_thread[ooc_threadnum]->written += towrite;
    OOC_CLOCK_STOP;
#ifdef OOC_CLOCK
    ooc->ooc_thread[ooc_threadnum]->time_writting += clockVal(&clk);
#endif
    if (EXIT_SUCCESS != close(fd))
      perror("close");
#ifdef DEBUG_OOC
    MUTEX_LOCK(&ooc->opened_mutex);
    ooc->opened--;
    MUTEX_UNLOCK(&ooc->opened_mutex);
#endif
  }
      /* freeing saved coeftab */
      MUTEX_LOCK_ALLOC;
      memFree_null(SOLV_UCOEFTAB(cblknum));
      ooc->allocated -= towrite;
      MUTEX_UNLOCK_ALLOC;
    }

  HASBEENSAVED(cblknum) = saved;

  ooc->cblkloaded--;
  print_debug(DBG_OOC_TRACE_V2,
              "<-ooc_do_save_coef %d %d\n", (int)cblknum, (int)ooc->allocated);

  return EXIT_SUCCESS;
}

#ifdef OOC_FTGT
/*
   Function: ooc_save_ftgt

   Had fanin target buffer to the save queue.
   Called by computing thread when
   work on a fanin target buffer
   has been done.

   Parameters:
     sopalin_data - sopalin_data structure (common data)
     ftgtnum      - fanin target buffer number
     me           - computing thread number

   Returns:
     EXIT_SUCCESS
 */
int ooc_save_ftgt(Sopalin_Data_t * sopalin_data, PASTIX_INT tasknum, PASTIX_INT ftgtnum,
                  int me)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  PASTIX_INT            key1, key2;

  print_debug(DBG_OOC_TRACE_V2,"->ooc_save_ftgt %d\n",(int)ftgtnum);

  HITFTGTTAB(ftgtnum) --;

  if (!(HITFTGTTAB(ftgtnum) > 0) &&
      (ooc->fanin_ctrbcnt[ftgtnum] > 0) &&
      (FANIN_CTRBCNT(ftgtnum) > 0))
    {
      cblkNextAccess(sopalin_data, tasknum, -(ftgtnum+1), me, &key1, &key2);

      MUTEX_LOCK(&(ooc->mutex_queue));
      queueAdd2(&(ooc->queue), -(ftgtnum+1), -key1, -key2);
      MUTEX_UNLOCK(&(ooc->mutex_queue));
      pthread_cond_signal(&(ooc->cond_queue));
    }
  else
    print_debug(DBG_OOC_PREDICTION, "NEXT %4ld %4ld NOTSAVE\n", (long)tasknum,
                (long)ftgtnum);

  print_debug(DBG_OOC_TRACE_V2,"<-ooc_save_coef\n");

  return EXIT_SUCCESS;
}

/*
   Function: ooc_remove_ftgt

   Delete fanin target buffer from disk.

   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     ftgtnum       - fanin buffer buffer number
     me            - ooc thread number
                     (should be on as only one is possible)
   Returns:
     EXIT_SUCCESS                    - deleted or file doesn't exist
     EXIT_FAILURE                    - couldn't delete
 */
int ooc_remove_ftgt (Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum, int me)
{
  ooc_t     * ooc      = sopalin_data->ooc;
  char        str[STR_SIZE];
  struct stat stFileInfo;

  print_debug(DBG_OOC_TRACE_V2,
              "->ooc_remove_ftgt %d (%d blocks %d allocated)\n",
              (int)ftgtnum, (int)ooc->cblkloaded, (int)ooc->allocated);

  /* removing the file */
  sprintf(str,"%s/pastix_ftgt_%d/%ld",OOC_DIR,(int)ooc->id,(long)ftgtnum);
  if (-1 == stat(str,&stFileInfo))
    {
      /* the file doesn't exist */
      return EXIT_SUCCESS;
    }
  if (-1 == remove(str))
    {
      errorPrint("Error opening file %s for write", str);
      perror("remove");
      exit(EXIT_FAILURE);
    }

  print_debug(DBG_OOC_TRACE_V2,"<-ooc_remove_ftgt %d %d\n",
        (int)ftgtnum, (int)ooc->allocated);

  return EXIT_SUCCESS;
}


/*
   Function: ooc_reset_ftgt

   Set all fanin target buffer to 0 or delete it from disk.

   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     ftgtnum       - fanin buffer buffer number
     me            - ooc thread number
                     (should be on as only one is possible)
   Returns:
     EXIT_SUCCESS                    - done succesfully
 */
int ooc_reset_ftgt (Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum, int me)
{
  SolverMatrix  *datacode = sopalin_data->datacode;
  ooc_t         *ooc      = sopalin_data->ooc;
  PASTIX_INT            j, oldsize, newsize, newsize2;

  /* Warning: on considere qu'on est forcement en memory_usage donc
   * l'octet précédent le pointeur contient la taille allouee */
  oldsize  = (PASTIX_INT)*((double*)FANIN_COEFTAB(ftgtnum)-1);

  newsize2 = SOLVE_FTGTSIZE(ftgtnum);
  newsize  = SOLVE_FTGTSIZE(ftgtnum)*sizeof(PASTIX_FLOAT);

  ASSERTDBG(((oldsize == FACTO_FTGTSIZE(ftgtnum)*sizeof(PASTIX_FLOAT)) ||
             (oldsize == SOLVE_FTGTSIZE(ftgtnum)*sizeof(PASTIX_FLOAT))), MOD_SOPALIN);

  print_debug(DBG_OOC_TRACE_V2,
              "->ooc_reset_ftgt %d (old %d new %d)(%d blocks %d allocated)\n",
              (int)ftgtnum, (int)oldsize, (int)newsize,
              (int)ooc->cblkloaded, (int)ooc->allocated);

  MUTEX_LOCK(&(sopalin_data->mutex_fanin[ftgtnum]));
  HITFTGTTAB(ftgtnum)--;

  ASSERTDBG(NULL != FANIN_COEFTAB(ftgtnum), MOD_SOPALIN);

  /* If hitcblktab > 0 need to realloc ftgt for next use */
  if (HITFTGTTAB(ftgtnum) > 0)
    {
      ASSERT(newsize <= oldsize, MOD_SOPALIN);
      if (newsize != oldsize)
        {
          /* Free Fanin target */
          MUTEX_LOCK_ALLOC;
          memFree_null(FANIN_COEFTAB(ftgtnum));
          ooc->allocated -= oldsize;

          MALLOC_INTERN(FANIN_COEFTAB(ftgtnum), newsize2, PASTIX_FLOAT);
          ooc->allocated += newsize;
          MUTEX_UNLOCK_ALLOC;
        }

      for (j=0; j<newsize2; j++)
        FANIN_COEFTAB(ftgtnum)[j] = 0.0;
    }
  else
    {
      MUTEX_LOCK_ALLOC;
      memFree_null(FANIN_COEFTAB(ftgtnum));
      FANIN_COEFTAB(ftgtnum) = NULL;
      ooc->allocated -= oldsize;
      MUTEX_UNLOCK_ALLOC;

      /* Remove file if ftgt is no more in memory */
      ooc_remove_ftgt(sopalin_data, ftgtnum, me);
    }
  MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[ftgtnum]));

  print_debug(DBG_OOC_TRACE_V2,"<-ooc_reset_ftgt %d (%d blocks %d allocated)\n",
              (int)ftgtnum, (int)ooc->cblkloaded, (int)ooc->allocated);

  return EXIT_SUCCESS;
}

/*
   Function: ooc_do_save_ftgt

   Save one fanin target buffer on disk.

   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     ftgtnum       - fanin target buffer number
     ooc_threadnum - ooc thread number
                     (should be on as only one is possible)
   Returns:
     EXIT_SUCCESS                    - saved
     EXIT_FAILURE_SAVING_NULL_BUFFER - buffer is null
     EXIT_FAILURE_FILE_OPENING       - error opening file
     EXIT_FAILURE                    - error writing
 */
int ooc_do_save_ftgt (Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum,
                      int ooc_threadnum)
{
  SolverMatrix       * datacode = sopalin_data->datacode;
  ooc_t              * ooc      = sopalin_data->ooc;
  char                 str[STR_SIZE];
  int                  fd;
  PASTIX_INT                  ret;
  PASTIX_INT                  written;
  PASTIX_INT                  towrite;

  print_debug(DBG_OOC_TRACE_V2,
              "->ooc_do_save_ftgt %d (%d blocks %d allocated)\n",
              (int)ftgtnum, (int)ooc->cblkloaded, (int)ooc->allocated);

  if (FANIN_COEFTAB(ftgtnum) == NULL)
    {
      return EXIT_FAILURE_SAVING_NULL_BUFFER;
    }
  /* opening the file */
  sprintf(str,"%s/pastix_ftgt_%d/%d",OOC_DIR,(int)ooc->id, (int)ftgtnum);
  if (-1 == (fd = open(str,O_WRONLY|O_CREAT|O_SYNC,
                       S_IRUSR | S_IWUSR | S_IRGRP |
                       S_IWGRP | S_IROTH | S_IWOTH)))
    {
      errorPrint("Error opening file %s for write", str);
      perror("open (write)");
      exit(-1);
      return EXIT_FAILURE_FILE_OPENING;
    }

#ifdef DEBUG_OOC
  MUTEX_LOCK(&ooc->opened_mutex);
  ooc->opened++;
  MUTEX_UNLOCK(&ooc->opened_mutex);
#endif

  /* writing the file : write only during facto */
  towrite = FACTO_FTGTSIZE(ftgtnum)*sizeof(PASTIX_FLOAT);
  written = 0;
  while ( written < towrite )
    {
      if (0 == (ret = write(fd,FANIN_COEFTAB(ftgtnum)+written,towrite-written)))
        {
          errorPrint("Coulnd't write to disk");
          return EXIT_FAILURE;
        }
      written += ret;
    }

  if (EXIT_SUCCESS != close(fd))
    perror("close");
  ooc->ooc_thread[ooc_threadnum]->written += towrite;
#ifdef DEBUG_OOC
  MUTEX_LOCK(&ooc->opened_mutex);
  ooc->opened--;
  MUTEX_UNLOCK(&ooc->opened_mutex);
#endif
  /* freeing saved coeftab */
  MUTEX_LOCK_ALLOC;
  memFree_null(FANIN_COEFTAB(ftgtnum));
  FANIN_COEFTAB(ftgtnum) = NULL;
  ooc->allocated -= towrite;
  MUTEX_UNLOCK_ALLOC;

  ooc->ftgtloaded--;

  print_debug(DBG_OOC_TRACE_V2,"<-ooc_do_save_ftgt %d %d\n",
              (int)ftgtnum, (int)ooc->allocated);

  return EXIT_SUCCESS;
}

/*
   Function: ooc_do_load_ftgt

   Allocate fanin target buffer and read it from disk if it exists.

   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     ftgtnum       - fanin target buffer number
     me            - ooc thread number
                     (should be on as only one is possible)
   Returns:
     EXIT_SUCCESS                    - loaded
     EXIT_FAILURE_OUT_OF_MEMORY      - No more memory
     EXIT_FAILURE_CBLK_NOT_NULL      - block column not *NULL*
     EXIT_FAILURE_FILE_OPENING       - error opening file
     EXIT_FAILURE                    - error writing
 */
int ooc_do_load_ftgt (Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum, int me)
{
  SolverMatrix       * datacode   = sopalin_data->datacode;
  ooc_t              * ooc      = sopalin_data->ooc;
  struct stat          stFileInfo;
  char                 str[STR_SIZE];
  int                  fd;
  PASTIX_INT                  toread;
  PASTIX_INT                  rd = 0;
  PASTIX_INT                  ret;
  PASTIX_INT                  j;
  PASTIX_INT                  coefnbr;

  print_debug(DBG_OOC_TRACE_V2,"->ooc_do_load_ftgt %d (%d allocated)\n",
              (int)ftgtnum, (int)ooc->allocated );

#ifdef DEBUG_OOC2
  if (EBUSY != pthread_mutex_trylock(&(sopalin_data->mutex_fanin[ftgtnum])))
    {
      errorPrint("Attention mutex pas verouille\n");
      while(1);
    }
#endif

  /* wait until we have enough memory */
  ret = EXIT_FAILURE_OUT_OF_MEMORY;

  if (EXIT_FAILURE_OUT_OF_MEMORY == (ret = ooc_allocate_ftgt(sopalin_data,
                                                             ftgtnum, me)))
    {
      errorPrint("%s:%d no more memory, shouldn't appear\n",__FILE__,__LINE__);
    }

  ooc->ftgtloaded ++;
  switch(ret)
    {
    case EXIT_FAILURE_OUT_OF_MEMORY:
      return EXIT_FAILURE_OUT_OF_MEMORY;
      break;
    case EXIT_FAILURE_CBLK_NOT_NULL:
      return EXIT_FAILURE_CBLK_NOT_NULL;
      break;
    case EXIT_SUCCESS :
      coefnbr = FTGTSIZE(ftgtnum);
      toread  = coefnbr*sizeof(PASTIX_FLOAT);

      /* reading the file */
      sprintf(str,"%s/pastix_ftgt_%d/%d",OOC_DIR,(int)ooc->id,(int)ftgtnum);
      if (-1 == stat(str,&stFileInfo))
        {
          /* the file doesn't exist, fill it with 0 */
          for (j = 0; j < coefnbr; j++)
            {
              FANIN_COEFTAB(ftgtnum)[j] = ZERO;
            }
          return EXIT_SUCCESS;
        }
      if (-1 == (fd = open(str,O_RDONLY)))
        {
          errorPrint("Error opening file %s for read",str);
          perror("open");
          EXIT(MOD_SOPALIN,0);
          return EXIT_FAILURE_FILE_OPENING;
        }

      rd = 0;
      while (rd < toread)
        {
          ret = read(fd,FANIN_COEFTAB(ftgtnum)+rd,toread - rd);
          switch (ret)
            {
            case 0:
              errorPrint("%s:%d File truncated",__FILE__,__LINE__);
              if (EXIT_SUCCESS != close(fd))
                perror("close");
              EXIT(MOD_SOPALIN,UNKNOWN_ERR);
              break;
            case -1:
              errorPrint("%s:%d reading file : %s",
                         __FILE__,__LINE__,strerror(errno));
              if (EXIT_SUCCESS != close(fd))
                perror("close");
              EXIT(MOD_SOPALIN,UNKNOWN_ERR);
              break;
            default:
              rd += ret;
            }
          ooc->ooc_thread[me]->read += toread;

        }
      if (EXIT_SUCCESS != close(fd))
        perror("close");

      print_debug(DBG_OOC_TRACE_V2,"<-ooc_do_load_ftgt %d %d\n",
                  (int)ftgtnum, (int)ooc->allocated);

      return EXIT_SUCCESS;
      break;
    }
  return EXIT_SUCCESS;
}


/*
   Function: ooc_wait_for_ftgt

   Check if fanin target buffer is NULL.
   If yes wait untill it is loaded.
   If not continue

   Parameters:
     sopalin_data - sopalin_data structure (common data)
     cblk         - fanin target number
     me           - computing thread

   Returns:
     EXIT_SUCCESS
 */
int ooc_wait_for_ftgt(Sopalin_Data_t * sopalin_data, PASTIX_INT ftgt, int me)
{
  SolverMatrix *datacode = sopalin_data->datacode;
  ooc_t        * ooc     = sopalin_data->ooc;
  int            hacked  = 0;
#ifdef OOC_DETECT_DEADLOCKS
  int            wasnull = 0;
#endif /* OOC_DETECT_DEADLOCKS */
#ifdef OOC_CLOCK
  Clock         clk;
#endif

  print_debug(DBG_OOC_TRACE_V2, "[%d]-> ooc_wait_for_ftgt %d\n",
              (int)me, (int)ftgt);

  OOC_CLOCK_INIT;
#ifdef OOC_DETECT_DEADLOCKS
  if (FANIN_COEFTAB(ftgt) == NULL)
    {
      print_debug(DBG_OOC_TRACE_V1,"waiting for ftgt... (%d)\n", (int)ftgt);
      wasnull = 1;
      MUTEX_LOCK(&(ooc->mutex_thrd_wting));
      ooc->blocks_waited_for[me] = - (ftgt+1);
      ooc->threads_waiting++;
      MUTEX_UNLOCK(&(ooc->mutex_thrd_wting));
    }
#endif /* OOC_DETECT_DEADLOCKS */
  while (FANIN_COEFTAB(ftgt) == NULL)
    {
      COND_TIMEWAIT(&(sopalin_data->cond_fanin[ftgt]),
                    &(sopalin_data->mutex_fanin[ftgt]));

      if (FANIN_COEFTAB(ftgt) == NULL)
        {
          if (SOLV_THRDNBR > 0)
            {
              MUTEX_LOCK(&(ooc->mutex_thrd_rcving));
              if (ooc->threads_receiving == SOLV_THRDNBR -1)
                {
                  ASSERTDBG(hacked == 0, MOD_SOPALIN);
                  hacked = 1;
                  /*      fprintf(stderr,"Hacking FTGT %d\n",ftgt); */
                  ooc_hack_load(sopalin_data, -(ftgt+1),  me);
                  HACKEDFTGT(ftgt)++;
                }
              MUTEX_UNLOCK(&(ooc->mutex_thrd_rcving));
            }
        }

    }
#ifdef OOC_DETECT_DEADLOCKS
  if (wasnull == 1)
    {
      print_debug(DBG_OOC_TRACE_V1,"done ftgt (%d)n", (int)ftgt);

      MUTEX_LOCK(&(ooc->mutex_thrd_wting));
      ooc->blocks_waited_for[me] = SYMB_CBLKNBR+1;
      ooc->threads_waiting--;
      MUTEX_UNLOCK(&(ooc->mutex_thrd_wting));
    }
#endif /* OOC_DETECT_DEADLOCKS */
  print_debug(DBG_OOC_TRACE_V2, "<- ooc_wait_for_ftgt %d\n",(int)ftgt);
  OOC_CLOCK_STOP;
  OOC_CLOCK_GET;


  return EXIT_SUCCESS;

}

/*
   Function: ooc_allocate_ftgt

   Allocate fanin target buffer.

   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     ftgtnum       - fanin target buffer number
     me            - ooc thread number
                     (should be on as only one is possible)

   Returns:
     EXIT_SUCCESS               - Allocated
     EXIT_FAILURE_OUT_OF_MEMORY - No more memory
     EXIT_FAILURE_CBLK_NOT_NULL - fanin target buffer not *NULL*
 */
int ooc_allocate_ftgt(Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum, int me)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  PASTIX_INT            size;
  print_debug(DBG_OOC_TRACE_V2,"%d:%d ->ooc_allocate_ftgt %d \n",
              (int)SOLV_PROCNUM, (int)me, (int)ftgtnum);

  /* ftgt */
  if (FANIN_COEFTAB(ftgtnum) != NULL)
    {
      return EXIT_FAILURE_CBLK_NOT_NULL;
    }

  size = FTGTSIZE(ftgtnum)*sizeof(PASTIX_FLOAT);

  MUTEX_LOCK_ALLOC;
  MALLOC_INTERN(FANIN_COEFTAB(ftgtnum), FTGTSIZE(ftgtnum), PASTIX_FLOAT);
  /* allocated is already protected above */
  ooc->allocated += size;
  MUTEX_UNLOCK_ALLOC;
  /*   pthread_cond_broadcast(&ooc->cond_allocated); */

  print_debug(DBG_OOC_TRACE_V2,"<-ooc_allocate_ftgt\n");

  return EXIT_SUCCESS;
}

#endif /* OOC_FTGT */

/*
   Function: ooc_do_load_coef

   Allocate block column and read it from disk if it exists.

   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     cblknum       - column block number
     me            - ooc thread number
                     (should be on as only one is possible)
   Returns:
     EXIT_SUCCESS                    - loaded
     EXIT_FAILURE_OUT_OF_MEMORY      - No more memory
     EXIT_FAILURE_CBLK_NOT_NULL      - block column not *NULL*
     EXIT_FAILURE_FILE_OPENING       - error opening file
     EXIT_FAILURE                    - error writing
 */
int ooc_do_load_coef (Sopalin_Data_t * sopalin_data, PASTIX_INT cblknum, int me)
{
  SolverMatrix       * datacode   = sopalin_data->datacode;
  ooc_t              * ooc        = sopalin_data->ooc;
  struct stat          stFileInfo;
  char                 str[STR_SIZE];
  int                  fd;
  PASTIX_INT                  coefnbr;
  PASTIX_INT                  toread;
  PASTIX_INT                  rd = 0;
  PASTIX_INT                  ret;
  PASTIX_INT                  j;
#ifdef OOC_CLOCK
  Clock clk;
#endif

  print_debug(DBG_OOC_TRACE_V2,
              "->ooc_do_load_coef %d (%d blocs, %d allocated)\n",
              (int)cblknum, (int)ooc->cblkloaded,
              (int)ooc->allocated );

#ifdef DEBUG_OOC2
  if (EBUSY != pthread_mutex_trylock(&(ooc->mutex_cblk[cblknum])))
    {
      errorPrint("Attention mutex pas verouille\n");
      while(1);
    }
#endif

  /* wait until we have enough memory */
  ret = EXIT_FAILURE_OUT_OF_MEMORY;

  if (EXIT_FAILURE_OUT_OF_MEMORY == (ret = ooc_allocate(sopalin_data,
                                                        cblknum, me)))
    {
      errorPrint("%s:%d no more memory, shouldn't appear\n",__FILE__,__LINE__);
    }
  ooc->cblkloaded ++;
  switch(ret)
    {
    case EXIT_FAILURE_OUT_OF_MEMORY:
      return EXIT_FAILURE_OUT_OF_MEMORY;
      break;
    case EXIT_FAILURE_CBLK_NOT_NULL:
      return EXIT_FAILURE_CBLK_NOT_NULL;
      break;
    case EXIT_SUCCESS :

      coefnbr = CBLKSIZE(cblknum);
      toread  = coefnbr * sizeof(PASTIX_FLOAT);

      /* reading the file */
      sprintf(str,"%s/pastix_coef_%d/%d",OOC_DIR,(int)ooc->id,(int)cblknum);
      if (-1 == stat(str,&stFileInfo))
        {
          /* the file doesn't exist, fill it with 0 */
          for (j = 0; j < coefnbr; j++)
            {
              SOLV_COEFTAB(cblknum)[j] = ZERO;
            }
          if (sopalin_data->sopar->factotype == API_FACT_LU )
            for (j = 0; j < coefnbr; j++)
              {
                SOLV_UCOEFTAB(cblknum)[j] = ZERO;
              }
#ifdef OOC_NOCOEFINIT
          Csc2solv_cblk(&(datacode->cscmtx), datacode,
                        sopalin_data->sopar->transcsc, cblknum);
#endif
          return EXIT_SUCCESS;
        }
      if (-1 == (fd = open(str,O_RDONLY)))
        {
          errorPrint("Error opening file %s for read",str);
          perror("open");
          EXIT(MOD_SOPALIN,0);
          return EXIT_FAILURE_FILE_OPENING;
        }
#ifdef DEBUG_OOC
      MUTEX_LOCK(&ooc->opened_mutex);
      ooc->opened++;
      MUTEX_UNLOCK(&ooc->opened_mutex);
#endif
      rd = 0;
      OOC_CLOCK_INIT;
      while (rd < toread)
        {
          ret = read(fd,SOLV_COEFTAB(cblknum)+rd,toread - rd);
          switch (ret)
            {
            case 0:
              errorPrint("%s:%d File truncated",__FILE__,__LINE__);
              if (EXIT_SUCCESS != close(fd))
                perror("close");
              EXIT(MOD_SOPALIN,UNKNOWN_ERR);
              break;
            case -1:
              errorPrint("%s:%d reading file : %s",__FILE__,__LINE__,
                         strerror(errno));
              if (EXIT_SUCCESS != close(fd))
                perror("close");
              EXIT(MOD_SOPALIN,UNKNOWN_ERR);
              break;
            default:
              rd += ret;
              break;
            }
          ooc->ooc_thread[me]->read += toread;
        }
      if (EXIT_SUCCESS != close(fd))
        perror("close");
      OOC_CLOCK_STOP;
      CBLKREAD(cblknum);
#ifdef OOC_CLOCK
      ooc->ooc_thread[me]->time_reading += clockVal(&clk);
#endif
#ifdef DEBUG_OOC
      MUTEX_LOCK(&ooc->opened_mutex);
      ooc->opened--;
      MUTEX_UNLOCK(&ooc->opened_mutex);
#endif

      /* read the file for LU */
      if (sopalin_data->sopar->factotype == API_FACT_LU)
        {
          sprintf(str,"%s/pastix_ucoef_%d/%d",OOC_DIR,
                  (int)ooc->id,(int)cblknum);
          if (-1 == (fd = open(str,O_RDONLY)))
            {
              errorPrint("Error opening file %s for read",str);
              perror("");
              EXIT(MOD_SOPALIN,0);
              return EXIT_FAILURE_FILE_OPENING;
            }
#ifdef DEBUG_OOC
          MUTEX_LOCK(&ooc->opened_mutex);
          ooc->opened++;
          MUTEX_UNLOCK(&ooc->opened_mutex);
#endif
          rd = 0;
          OOC_CLOCK_INIT;
          while (rd < toread)
            {
              if (0 == (ret = read(fd,SOLV_UCOEFTAB(cblknum)+rd,toread - rd)))
                {
                  errorPrint("File truncated");
                  return EXIT_FAILURE_FILE_TRUNCATED;
                }
              rd += ret;
            }
          ooc->ooc_thread[me]->read += toread;
          if (EXIT_SUCCESS != close(fd))
            perror("close");
          OOC_CLOCK_STOP;
#ifdef OOC_CLOCK
          ooc->ooc_thread[me]->time_reading += clockVal(&clk);
#endif
#ifdef DEBUG_OOC
          MUTEX_LOCK(&ooc->opened_mutex);
          ooc->opened--;
          MUTEX_UNLOCK(&ooc->opened_mutex);
#endif

        }

      print_debug(DBG_OOC_TRACE_V2,"<-ooc_do_load_coef %d %d\n",
                  (int)cblknum, (int)ooc->allocated);

      return EXIT_SUCCESS;
      break;
    }
  return EXIT_SUCCESS;
}




/*
   Function: ooc_allocate

   Allocate column block.

   Parameters:
     sopalin_data  - sopalin_data structure (common data)
     cblknum       - column block number
     me            - ooc thread number
                     (should be on as only one is possible)

   Returns:
     EXIT_SUCCESS               - Allocated
     EXIT_FAILURE_OUT_OF_MEMORY - No more memory
     EXIT_FAILURE_CBLK_NOT_NULL - Column block not *NULL*
 */
int ooc_allocate(Sopalin_Data_t * sopalin_data, PASTIX_INT cblknum, int me)
{
  SolverMatrix *datacode = sopalin_data->datacode;
  ooc_t        *ooc      = sopalin_data->ooc;
  PASTIX_INT           coefnbr  = CBLKSIZE(cblknum);
  PASTIX_INT           size     = coefnbr*sizeof(PASTIX_FLOAT);

  print_debug(DBG_OOC_TRACE_V2,"->ooc_allocate %d (%d blocks)\n",
              (int)cblknum, (int)ooc->cblkloaded);

  if (SOLV_COEFTAB(cblknum) != NULL)
    return EXIT_FAILURE_CBLK_NOT_NULL;

  MUTEX_LOCK_ALLOC;
  MALLOC_INTERN(SOLV_COEFTAB(cblknum), coefnbr, PASTIX_FLOAT);

  /* WARNING : allocated*/
  ooc->allocated += size;
  MUTEX_UNLOCK_ALLOC;

  if (sopalin_data->sopar->factotype == API_FACT_LU)
    {
      /* allocated is already protected above */
      if (SOLV_UCOEFTAB(cblknum) != NULL)
        {
          errorPrint("%s:%d SOLV_UCOEFTAB(%d) not null",
                     __FILE__,__LINE__,cblknum);
          EXIT(MOD_SOPALIN,UNKNOWN_ERR);
        }

      MUTEX_LOCK_ALLOC;
      MALLOC_INTERN(SOLV_UCOEFTAB(cblknum), coefnbr, PASTIX_FLOAT);

      ooc->allocated += size;
      MUTEX_UNLOCK_ALLOC;
    }

  print_debug(DBG_OOC_TRACE_V2,"<-ooc_allocate\n");

  return EXIT_SUCCESS;
}


/*
   Function: ooc_wait_for_cblk

   Check if column block is NULL.
   If yes wait untill it is loaded.
   If not continue

   Parameters:
     sopalin_data - sopalin_data structure (common data)
     cblk         - column block number
     me           - computing thread

   Returns:
     EXIT_SUCCESS
 */
int ooc_wait_for_cblk(Sopalin_Data_t * sopalin_data, PASTIX_INT cblk, int me)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  int            hacked   = 0;
#ifdef OOC_DETECT_DEADLOCKS
  int            wasnull  = 0;
#endif
#ifdef OOC_CLOCK
  Clock clk;
#endif

  print_debug(DBG_OOC_TRACE_V2, "[%d]-> ooc_wait_for_cblk %d\n",
        (int)me, (int)cblk);


  OOC_CLOCK_INIT;

  MUTEX_LOCK(&(ooc->mutex_cblk[cblk]));
  ooc->thrdUsingCblk[cblk].count++;

  CBLKWAIT(cblk);
#ifdef OOC_DETECT_DEADLOCKS
  if ((SOLV_COEFTAB(cblk) == NULL)
      || ((sopalin_data->sopar->factotype == API_FACT_LU)
          && (SOLV_UCOEFTAB(cblk) == NULL)))
    {
      wasnull = 1;
      print_debug(DBG_OOC_TRACE_V1,"waiting... (%d)\n", (int)cblk);
      MUTEX_LOCK(&(ooc->mutex_thrd_wting));
      ooc->blocks_waited_for[me] = cblk;
      ooc->threads_waiting++;
      MUTEX_UNLOCK(&(ooc->mutex_thrd_wting));
    }
#endif /* OOC_DETECT_DEADLOCKS */
  /* Thread Wait for cblk */
  while ((SOLV_COEFTAB(cblk) == NULL)
         || ((sopalin_data->sopar->factotype == API_FACT_LU)
             && (SOLV_UCOEFTAB(cblk) == NULL)))
    {

      COND_TIMEWAIT(&(ooc->cond_cblk[cblk]), &(ooc->mutex_cblk[cblk]));
      if (SOLV_COEFTAB(cblk) == NULL)
        {
          if (SOLV_THRDNBR > 0)
            {
              MUTEX_LOCK(&(ooc->mutex_thrd_rcving));
              if (ooc->threads_receiving == SOLV_THRDNBR -1)
                {
                  ASSERTDBG(hacked == 0, MOD_SOPALIN);
                  hacked = 1;
                  /*      fprintf(stderr,"Hacking cblk %d\n",cblk); */
                  ooc_hack_load(sopalin_data, cblk,  me);
                  HACKEDCBLK(cblk)++;
                }
              MUTEX_UNLOCK(&(ooc->mutex_thrd_rcving));
            }
        }
    }

#ifdef OOC_DETECT_DEADLOCKS
  if (wasnull == 1)
    {
      print_debug(DBG_OOC_TRACE_V1,"done\n");
      MUTEX_LOCK(&(ooc->mutex_thrd_wting));
      ooc->blocks_waited_for[me] = SYMB_CBLKNBR+1;
      ooc->threads_waiting--;
      MUTEX_UNLOCK(&(ooc->mutex_thrd_wting));
    }
#endif /* OOC_DETECT_DEADLOCKS */
  MUTEX_UNLOCK(&(ooc->mutex_cblk[cblk]));

  OOC_CLOCK_STOP;
  OOC_CLOCK_GET;

  print_debug(DBG_OOC_TRACE_V2, "[%d]<- ooc_wait_for_cblk %d\n",
        (int)me, (int)cblk);

  return EXIT_SUCCESS;
}

/*
   Function: ooc_stop_thread

   Stops all ooc thread.

   Parameters:
     sopalin_data - sopalin_data structure (common data)

   Returns:
     EXIT_SUCCESS
 */
int ooc_stop_thread(Sopalin_Data_t * sopalin_data)
{
  ooc_t * ooc = sopalin_data->ooc;
  ooc->stop = 1;

  pthread_cond_signal(&(ooc->cond_toload));
  pthread_cond_signal(&(ooc->cond_frozen));
  pthread_cond_signal(&(ooc->cond_queue));
  pthread_cond_signal(&(ooc->cond_step));
  /*   for (i = 0 ; i < ooc->threadnbr; i++) */
  /*     pthread_cond_broadcast(&ooc->ooc_thread[i]->cond_queue); */
  /*   OOC_SYNCHRO_THREAD; */
  return EXIT_SUCCESS;
}


/*
   Function: cblkNextAccess

   Compute the number of the next access to a column block or fanin target.

   Parameters:
     sopalin_data - sopalin_data structure (common data)
     cblk         - column block or fanin target number
     comp_thread  - computing thread saving cblk.

   Returns:
     Next access number.
*/
void cblkNextAccess(Sopalin_Data_t *sopalin_data, PASTIX_INT task, PASTIX_INT cblk,
                    int comp_thread,
                    PASTIX_INT * key1, PASTIX_INT *key2)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  PASTIX_INT n;
  PASTIX_INT cblknum;
  PASTIX_INT bloknum;
  PASTIX_INT gcblk;
  int step;

  step  = ooc->current_step;
  *key1 = 0;
  *key2 = 0;

  /* On passe dans le if que si on est dans OOCSTEP_SOPALIN
   * et qu'il s'agit d'une fanin                             */
  if (cblk < 0)
    {
      if (step != OOCSTEP_SOPALIN)
  {
    errorPrint("Ftgt can't be saved outside factorization");
    EXIT(MOD_SOPALIN, INTERNAL_ERR);
  }

      cblk  = -(cblk+1);
      gcblk = UPDOWN_GCBLK2LIST(FANIN_GCBKDST(cblk));

      /* Next acces in next contribution */
      cblknum = TASK_CBLKNUM(task);
      n       = UPDOWN_LISTPTR(gcblk);

      /* Search for first cblknum in list of block which contribute to cblk */
      while ((UPDOWN_LISTCBLK(n) != cblknum)
             && (n < UPDOWN_LISTPTR(gcblk+1)))
        n++;

      /* Search for last cblknum in list */
      while((UPDOWN_LISTCBLK(n) == cblknum)
            && (n < UPDOWN_LISTPTR(gcblk+1)))
        n++;

      /* Search for first extra-block facing with ftgt */
      while( (n < UPDOWN_LISTPTR(gcblk+1) ) &&
             (!((SYMB_FROWNUM(UPDOWN_LISTBLOK(n)) >= FANIN_FROWNUM(cblk)) &&
                (SYMB_LROWNUM(UPDOWN_LISTBLOK(n)) <= FANIN_LROWNUM(cblk)))))
        n++;

      if ( n < UPDOWN_LISTPTR(gcblk+1) )
        {
          *key1 = TASK_PRIONUM(UPDOWN_LISTCBLK(n));
          *key2 = UPDOWN_LISTBLOK(n);

          print_debug(DBG_OOC_FTGT,
                      "NEXT %4ld %4ld %4ld %4ld\n", (long)cblk, (long)task,
                      (long)UPDOWN_LISTCBLK(n), (long)UPDOWN_LISTBLOK(n));
        }
      else
        {
          *key1 = 1;

          print_debug(DBG_OOC_FTGT, "NEXT %4ld %4ld %4ld %4ld\n",
                      (long)cblk, (long)task,
                      (long)*key1, (long)*key2);
        }
      return;
    }

  /* Cas d'un coeftab */
  gcblk = UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(cblk));

  /* Hack for recpetion */
  if (task == -1)
    step = OOCSTEP_COEFINIT;

  switch ( step )
    {
    case OOCSTEP_COEFINIT:
      /* Next access in factorization step (warn to first cblk) */
      if (gcblk == -1)
        {
          *key1 = TASK_PRIONUM(cblk);
        }
      else
        {
          *key1 = TASK_PRIONUM(UPDOWN_LISTCBLK(UPDOWN_LISTPTR(gcblk)));
          *key2 = UPDOWN_LISTBLOK(UPDOWN_LISTPTR(gcblk));
        }
      break;

    case OOCSTEP_SOPALIN:
       /* Next access in Down step */
      cblknum = TASK_CBLKNUM(task);
      if ((gcblk == -1) || (cblknum == cblk))
        {
          HASBEENSAVED(cblknum) = -1;
          *key1 = (SYMB_CBLKNBR + TASK_PRIONUM(cblk));
          return;
        }

      /* Next acces in next contribution */
      bloknum = SYMB_BLOKNUM(cblknum)+1;
      n       = UPDOWN_LISTPTR(gcblk);

      /* Search for first cblknum occurence in the list of block which
         contribute to cblk */
      while ((n < UPDOWN_LISTPTR(gcblk+1) &&
              (UPDOWN_LISTCBLK(n) != cblknum)))
        n++;

      /* Search for last cblknum occurence in the list */
      while((n < UPDOWN_LISTPTR(gcblk+1) &&
             (UPDOWN_LISTCBLK(n) == cblknum)))
        n++;

      /* Search for first extra-block facing with cblk */
      while( (n < UPDOWN_LISTPTR(gcblk+1) ) &&
             (!((SYMB_FROWNUM(UPDOWN_LISTBLOK(n)) >= SYMB_FCOLNUM(cblk)) &&
                (SYMB_LROWNUM(UPDOWN_LISTBLOK(n)) <= SYMB_LCOLNUM(cblk)))))
        n++;

      if ( n == (UPDOWN_LISTPTR(gcblk+1)))
        {
          *key1 = TASK_PRIONUM(cblk);
        }
      else if ( n < UPDOWN_LISTPTR(gcblk+1) )
        {
          *key1 = TASK_PRIONUM(UPDOWN_LISTCBLK(n));
          *key2 = UPDOWN_LISTBLOK(n);
        }
      else
        {
          errorPrint("%s:%d cblk %d not found in list of cblk using %d \n",
                     __FILE__,__LINE__,cblknum,cblk);
          EXIT(INTERNAL_ERR, MOD_SOPALIN);
        }
      break;

    case OOCSTEP_DOWN:
      *key1 = (SYMB_CBLKNBR + (2 * (SYMB_CBLKNBR - TASK_PRIONUM(task))));
      break;

    case OOCSTEP_UP:
      *key1 = (SYMB_CBLKNBR + (2 * TASK_PRIONUM(task)));
      break;
    }

  return;
}

/*
   Function: reduceMem

   Save some column block or fanin target buffers until more than
   *size_wanted* memory is available.

   Parameters:
     sopalin_data - sopalin_data structure (common data)
     size_wanted  - minimal memory needed
     me           - ooc thread number
                     (should be on as only one is possible)

   Returns:
     void
 */
void reduceMem(Sopalin_Data_t * sopalin_data, PASTIX_INT size_wanted, int me)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  PASTIX_INT            indbloc;
#ifdef OOC_DETECT_DEADLOCKS
  PASTIX_INT            i;
#endif

  print_debug(DBG_OOC_TRACE_V2,"->reduceMem %d %d %d\n",
              (int)ooc->allocated, (int)size_wanted,
              (int)queueSize(&ooc->queue));

  MUTEX_LOCK(&(ooc->mutex_queue));


  while ((
#ifdef OOC_PERCENT_COEFTAB
          ooc->allocated + size_wanted >
          ooc->global_limit * ooc->coeftabsize / 100
#else
          (memAllocGetCurrent() + size_wanted) > ooc->global_limit
#endif
          )
         && (ooc->stop != 1))

    {
      /* If there is something to save */
      if (queueSize(&(ooc->queue)) > 0)
        {
          indbloc = queueGet(&(ooc->queue));
          MUTEX_UNLOCK(&(ooc->mutex_queue));

          /* Free a cblk */
          if( indbloc >= 0 )
            {
              MUTEX_LOCK(&ooc->mutex_cblk[indbloc]);
              if ((SOLV_COEFTAB(indbloc) != NULL)
                  && (!(HITCBLKTAB(indbloc) > 0))
                  && (!(ooc->thrdUsingCblk[indbloc].count > 0)))
                {
                  print_debug(DBG_OOC_SAVE, "SAVE %ld %ld %ld\n",
                              (long)indbloc, (long)size_wanted,
                              (long)memAllocGetCurrent());
                  ooc_do_save_coef(sopalin_data, indbloc, me);
                }
              MUTEX_UNLOCK(&ooc->mutex_cblk[indbloc]);
            }
#ifdef OOC_FTGT
          /* Free a ftgt */
          else
            {
              indbloc = -(indbloc + 1);

              ASSERT(me < SOLV_THRDNBR,MOD_SOPALIN);
              ASSERT(indbloc < SOLV_FTGTNBR,MOD_SOPALIN);

              MUTEX_LOCK(&sopalin_data->mutex_fanin[indbloc]);
              if ((FANIN_COEFTAB(indbloc) != NULL)
                  && (!(HITFTGTTAB(indbloc) > 0))
                  && (FANIN_CTRBCNT(indbloc) > 0)
                  && (ooc->fanin_ctrbcnt[indbloc] > 0))
                {
                  print_debug(DBG_OOC_FTGT, "SAVE %ld %ld %ld\n",
                              (long)indbloc, (long)size_wanted,
                              (long)memAllocGetCurrent());
                  ooc_do_save_ftgt(sopalin_data, indbloc, me);
                }
              MUTEX_UNLOCK(&sopalin_data->mutex_fanin[indbloc]);
            }
#endif

          MUTEX_LOCK(&(ooc->mutex_queue));
        }
      else
        {
          /* Test if limit isn't too small  */
#ifdef OOC_DETECT_DEADLOCKS
          MUTEX_LOCK(&(ooc->mutex_thrd_wting));
          if (ooc->threads_waiting == SOLV_THRDNBR)
            {
              for (i = 0; i <  SOLV_THRDNBR; i++)
                {
                  if (ooc->blocks_waited_for[i] != -1 &&
                      ( ooc->blocks_waited_for[i] == SYMB_CBLKNBR+1 ||
                        SOLV_COEFTAB(ooc->blocks_waited_for[i]) != NULL))
                    {
                      goto NoDeadLock;
                    }
                }
              MUTEX_UNLOCK(&(ooc->mutex_thrd_wting));

              PASTIX_INT min_mem = size_wanted;

              MUTEX_UNLOCK(&(ooc->mutex_queue));
#ifdef DEBUG_OOC
              {
                PASTIX_INT i;
                fprintf(stderr, "CBLK chargés  : ");
                for (i=0; i< SYMB_CBLKNBR; i++)
                  {
                    if (SOLV_COEFTAB(i) != NULL) fprintf(stderr, "%d ", i);
                  }
                fprintf(stderr, "\nFTGT chargées : ");
                for (i=0; i< SOLV_FTGTNBR; i++)
                  {
                    if (FANIN_COEFTAB(i) != NULL)
                      fprintf(stderr, "(%d %d %d) ", i,
                              FANIN_CTRBCNT(i), HITFTGTTAB(i));
                  }
                fprintf(stderr, "\n");
                queuePrint(&(ooc->queue));
              }
#endif
              NOT_ENOUGH_MEM;
              EXIT(MOD_SOPALIN,OUTOFMEMORY_ERR);
            }
          else
#endif /* OOC_DETECT_DEADLOCKS */
            {
#ifdef  OOC_DETECT_DEADLOCKS
            NoDeadLock:
              MUTEX_UNLOCK(&(ooc->mutex_thrd_wting));
#endif /* OOC_DETECT_DEADLOCKS */
              /* Check the Hack Queue */
              MUTEX_LOCK(&ooc->mutex_toload);
              while (queueSize(&ooc->toload) > 0)
                ooc_prefetch(sopalin_data, queueGet(&ooc->toload),
                             0, HIGH); /* WARNING : 0... */
              MUTEX_UNLOCK(&ooc->mutex_toload);
              /* Wait for something to save */
              COND_TIMEWAIT(&(ooc->cond_queue), &(ooc->mutex_queue));
            }
        }
    }
  MUTEX_UNLOCK(&(ooc->mutex_queue));

  print_debug(DBG_OOC_TRACE_V2,"<-reduceMem %d\n",
              (int)ooc->allocated);

}

#ifdef OOC_CLOCK
/*
   Function: ooc_clock_get

   Adds time in Clock *clk* to the total time waiting for OOC to load data.
   Defined only if *OOC_CLOCK* is defined.

   Parameters:
     sopalin_data     - sopalin_data structure (common data)
     clk              - clock from wich we get time
     computing_thread - Computing thread waitiing
   Returns:
     EXIT_SUCCESS

*/
int ooc_clock_get(Sopalin_Data_t * sopalin_data, Clock * clk,
                  int computing_thread)
{
  ooc_t        * ooc      = sopalin_data->ooc;

  ooc->time_waiting[computing_thread] += clockVal(clk);
  return EXIT_SUCCESS;
}
#endif

/*
   Function: recursive_mkdir

   Create directory and subdirectory indicated by *path*

   Parameters:
     path - path to the directory to create
     mode - mode for the created directory
            (Ex: S_IRUSR | S_IWUSR | S_IRGRP | S_IXUSR | S_IWGRP |
                 S_IROTH | S_IWOTH)
   Returns:
     EXIT_SUCCESS

*/
int recursive_mkdir(char * path, mode_t mode)
{
  char *copy, *p;

  MALLOC_INTERN(p, strlen(path)+1, char);
  copy = p;
  sprintf(p,"%s",path);

  do {
    p = strchr (p + 1, '/');
    if (p)
      *p = '\0';
    if (access (copy, F_OK) == -1) {
      if (mkdir (copy, mode ) == -1) {
        errorPrint("mkdir %s failed: %s", path, strerror(errno));
        return EXIT_FAILURE;
      }
    }
    if (p)
      *p = '/';
  } while (p);
  memFree_null(copy);
  return EXIT_SUCCESS;
}

/*
  Function: ooc_hack_load

  Add a bloc column to priority load queue.

  Parameters:
    sopalin_data - Global PaStiX structure.
    cblknum      - column bloc to load.
    me           - Computing thread asking for cblk to be loaded.

*/
int ooc_hack_load(Sopalin_Data_t *sopalin_data, PASTIX_INT cblknum, int me)
{
  ooc_t        * ooc      = sopalin_data->ooc;

  print_debug(DBG_OOC_TRACE_V1,"Hacking loading queue %d\n", (int)cblknum);
  MUTEX_LOCK(&ooc->mutex_toload);
  queueAdd(&ooc->toload,cblknum,1);
  MUTEX_UNLOCK(&ooc->mutex_toload);
  pthread_cond_signal(&ooc->cond_toload);
  return EXIT_SUCCESS;
}

int ooc_set_step(Sopalin_Data_t *sopalin_data, int step)
{
  ooc_t        * ooc      = sopalin_data->ooc;

  MUTEX_LOCK(&(ooc->mutex_step));
  ooc->current_step = step;
  if (step == OOCSTEP_DOWN)
    ooc->current_itercomp++;
  MUTEX_UNLOCK(&(ooc->mutex_step));
  pthread_cond_signal(&(ooc->cond_step));
  return EXIT_SUCCESS;
}

int ooc_freeze(Sopalin_Data_t *sopalin_data)
{
  ooc_t        * ooc      = sopalin_data->ooc;

/*   fprintf(stdout, "Freezing OOC %ld\n",memAllocGetCurrent()); */
  MUTEX_LOCK(&(ooc->mutex_frozen));
  ooc->frozen++;
  MUTEX_UNLOCK(&(ooc->mutex_frozen));
  pthread_cond_signal(&ooc->cond_frozen);
  return EXIT_SUCCESS;
}

int ooc_defreeze(Sopalin_Data_t *sopalin_data)
{
  ooc_t        * ooc      = sopalin_data->ooc;

/*   fprintf(stdout, "DeFreezing OOC %ld\n",memAllocGetCurrent()); */
  MUTEX_LOCK(&(ooc->mutex_frozen));
  ooc->frozen--;
  MUTEX_UNLOCK(&(ooc->mutex_frozen));
  pthread_cond_signal(&(ooc->cond_frozen));
  return EXIT_SUCCESS;
}

void ooc_receiving(Sopalin_Data_t * sopalin_data)
{
  ooc_t        * ooc      = sopalin_data->ooc;
  MUTEX_LOCK(&(ooc->mutex_thrd_rcving));
  ooc->threads_receiving ++;
  MUTEX_UNLOCK(&(ooc->mutex_thrd_rcving));
}
void ooc_received(Sopalin_Data_t * sopalin_data)
{
  ooc_t        * ooc      = sopalin_data->ooc;
  MUTEX_LOCK(&(ooc->mutex_thrd_rcving));
  ooc->threads_receiving --;
  MUTEX_UNLOCK(&(ooc->mutex_thrd_rcving));
}

void ooc_wait_task(Sopalin_Data_t * sopalin_data, PASTIX_INT task, int me)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  ooc_t        * ooc      = sopalin_data->ooc;
  MUTEX_LOCK(&(ooc->mutex_current_task));

  print_debug(DBG_OOC_TRACE_V2, "[%d]-> ooc_wait_task %d\n",
              (int)me, (int)task);

/*   fprintf(stdout, "%ld %ld %ld\n",(long)task,  */
/*    (long)ooc->current_taskooc,  */
/*    (long)datacode->tasknbr); */
#ifdef OOC_DETECT_DEADLOCKS
  MUTEX_LOCK(&(ooc->mutex_thrd_wting));
  ooc->threads_waiting++;
  ooc->blocks_waited_for[me] = -1;
  MUTEX_UNLOCK(&(ooc->mutex_thrd_wting));
#endif /* OOC_DETECT_DEADLOCKS */

  while (ooc->current_taskooc == -1)
    COND_WAIT(&(ooc->cond_current_task),&(ooc->mutex_current_task));

  ASSERT(task < datacode->tasknbr, MOD_SOPALIN);
  ASSERT(task >= 0, MOD_SOPALIN);
  ASSERT(ooc->current_taskooc < datacode->tasknbr, MOD_SOPALIN);
  ASSERT(ooc->current_taskooc >= 0, MOD_SOPALIN);
  while (TASK_PRIONUM(task) > TASK_PRIONUM(ooc->current_taskooc))
    COND_WAIT(&(ooc->cond_current_task),&(ooc->mutex_current_task));
  MUTEX_UNLOCK(&(ooc->mutex_current_task));

#ifdef OOC_DETECT_DEADLOCKS
  MUTEX_LOCK(&(ooc->mutex_thrd_wting));
  ooc->threads_waiting--;
  ooc->blocks_waited_for[me] = SYMB_CBLKNBR+1;
  MUTEX_UNLOCK(&(ooc->mutex_thrd_wting));
#endif /* OOC_DETECT_DEADLOCKS */

  print_debug(DBG_OOC_TRACE_V2, "[%d]<- ooc_wait_task %d\n",
        (int)me, (int)task);
}
#else
/* ISO C forbids an empty source file */
#include "not_empty.h"
NOT_EMPTY(ooc)
#endif
