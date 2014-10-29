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
  File: sopalin_thread.h

  Structures, function declarations and macros used for thread management.
 */
#ifndef SOPALIN_THREAD_H
#define SOPALIN_THREAD_H

#include <sys/time.h>
#include <time.h>


/*
  Struct: sopthread_data

  Structure conotaining the thread number and a pointer to data given to thread in parameters.
 */
typedef struct sopthread_data {
  int    me;                   /*+ Thread number                    +*/
  void * data;                 /*+ Data given to thread as argument +*/
} sopthread_data_t;

/*
 * Function: sopalin_launch_thread
 *
 * Launch all PaStiX threads on wanted functions.
 *
 * Parameters:
 *   procnum       - MPI process rank
 *   procnbr       - Number of MPI process
 *   ptr           - Pointer the the bubble structure to use (if MARCEL version)
 *   verbose       - Verbosity level
 *   calc_thrdnbr  - Number of computing threads.
 *   calc_routine  - Computing function.
 *   calc_data     - Parameters for computing function.
 *   comm_thrdnbr  - Number of communicating threads.
 *   comm_routine  - communication function.
 *   comm_data     - Parameters for communication function.
 *   ooc_thrdnbr   - Number of out-of-core threads.
 *   ooc_routine   - Out-of-core function.
 *   ooc_data      - Parameters for *ooc_routine*.
 */
void sopalin_launch_thread(void *sopalin_data,
                           PASTIX_INT procnum, PASTIX_INT procnbr, void *ptr, PASTIX_INT verbose,
			   PASTIX_INT calc_thrdnbr, void * (*calc_routine)(void *), void *calc_data,
			   PASTIX_INT comm_thrdnbr, void * (*comm_routine)(void *), void *comm_data,
			   PASTIX_INT ooc_thrdnbr,  void * (*ooc_routine) (void *), void *ooc_data);

/*
 * Function: sopalin_launch_comm
 *
 * Launch communication threads
 *
 * Parameters
 *   nbthrdcomm    - Number of threads to launch.
 *   comm_routine  - Communication function.
 *   data          - Data for communication function.
 */
void sopalin_launch_comm(int nbthrdcomm, void * (*comm_routine)(void *), void *data);


/*
 * Function: sopalin_bindthread
 *
 * Bind threads onto processors.
 *
 * Parameters:
 *   cpu - Processor to bind to.
 */
PASTIX_INT  sopalin_bindthread(PASTIX_INT);

/* Version SMP */
#ifndef FORCE_NOSMP

#  define MONOTHREAD_BEGIN if(me==0){
#  define MONOTHREAD_END   }

/*
 * Struct: sopthread_barrier
 *
 * Computing threads synchronisation barrier.
 */
typedef struct sopthread_barrier {
  int volatile    instance;         /*+ ID of the barrier                +*/
  int volatile    blocked_threads;  /*+ Number of threads in the barrier +*/
  pthread_mutex_t sync_lock;        /*+ mutex for the barrier            +*/
  pthread_cond_t  sync_cond;        /*+ cond for the barrier             +*/
} sopthread_barrier_t;

/*
 *  Macro: SYNCHRO_X_THREAD
 *
 *  Synchronize *nbthread* threads.
 *
 *  Parameters:
 *    nbthread - Number of threads to synchronize.
 *    barrier  - sopthread_barrier structure associated with the synchronisation.
 *
 */
#  define SYNCHRO_X_THREAD(nbthread, barrier)                           \
  {                                                                     \
    int instance;                                                       \
    pthread_mutex_lock(&((barrier).sync_lock));                         \
    instance = (barrier).instance;                                      \
    (barrier).blocked_threads++;                                        \
    if ((barrier).blocked_threads == (nbthread))                        \
      {                                                                 \
        (barrier).blocked_threads = 0;                                  \
        (barrier).instance++;                                           \
        pthread_cond_broadcast(&((barrier).sync_cond));                 \
      }                                                                 \
    while (instance == (barrier).instance)                              \
      {                                                                 \
        pthread_cond_wait(&((barrier).sync_cond),                       \
                          &((barrier).sync_lock));                      \
      }                                                                 \
    pthread_mutex_unlock(&((barrier).sync_lock));                       \
  }

/*
 * Définition de MUTEX_LOCK et COND_WAIT (avec ou sans compteurs)
 */
#  ifdef TRYLOCK
/*
PASTIX_INT *ptbusy,*ptfree;
PASTIX_INT *ptwait;
*/
#    define MUTEX_LOCK(x)   if (pthread_mutex_trylock(x)) {     \
    thread_data->ptbusy;pthread_mutex_lock(x);}                 \
  else thread_data->ptfree++
#    define COND_WAIT(x,y)  pthread_cond_wait(x,y); thread_data->ptwait++

/* Cond de 5ms */
#    define COND_TIMEWAIT(x,y) {                        \
    struct timeval  now;				\
    struct timespec timeout;				\
    gettimeofday(&now, NULL);				\
    timeout.tv_sec  = now.tv_sec;			\
    timeout.tv_nsec = now.tv_usec * 1000 + 5 * 1000000;	\
    pthread_cond_timedwait(x,y,&timeout);		\
    thread_data->ptwait++;				\
  }

#  else

#    define MUTEX_LOCK(x)      pthread_mutex_lock(x)
#    define COND_WAIT(x,y)     pthread_cond_wait(x,y)
/* Cond de 5ms */
#    define COND_TIMEWAIT(x,y) {                        \
    struct timeval  now;				\
    struct timespec timeout;				\
    gettimeofday(&now, NULL);				\
    timeout.tv_sec  = now.tv_sec;			\
    timeout.tv_nsec = now.tv_usec * 1000 + 5 * 1000000;	\
    pthread_cond_timedwait(x,y,&timeout);		\
  }

#  endif

#  define MUTEX_UNLOCK(x) pthread_mutex_unlock(x)

/* #if (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) */
/* #define SEM_TIMEWAIT(x, err) {						\ */
/*     struct timeval  now;						\ */
/*     struct timespec timeout;						\ */
/*     gettimeofday(&now, NULL);						\ */
/*     timeout.tv_sec  = now.tv_sec;					\ */
/*     timeout.tv_nsec = now.tv_usec * 1000 + 5 * 1000000;			\ */
/*     err = sem_timedwait(x,&timeout);					\ */
/*   } */
/* #else */
/* #define SEM_TIMEDWAIT(x, err) {err = sem_wait(x);} */
/* #endif */

/* Version Non-SMP */
#else /* FORCE_NOSMP */

#  define pthread_mutex_lock(x)
#  define pthread_mutex_unlock(x)

#  define MUTEX_LOCK(x)    {}
#  define MUTEX_UNLOCK(x)  {}

#  define pthread_cond_signal(x)
#  define pthread_cond_broadcast(x)
#  define COND_WAIT(x,y)
#  define COND_TIMEWAIT(x,y)

#  define SYNCHRO_X_THREAD(nbthrd, barrier)
#  define MONOTHREAD_BEGIN
#  define MONOTHREAD_END

#  define sopthread_barrier_t int
#endif /* FORCE_NOSMP */

/* SMP Macros required even with FORCE_NOSMP */
#define SYNCHRO_THREAD  SYNCHRO_X_THREAD(SOLV_THRDNBR, sopalin_data->barrier)
#define MAXTHRDS        SOLV_THRDNBR

#endif /* SOPALIN_THREAD_H */

