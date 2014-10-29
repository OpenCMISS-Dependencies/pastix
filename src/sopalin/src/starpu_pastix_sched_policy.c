#if (defined WITH_STARPU && defined STARPU_PASTIX_SCHED)
/* StarPU --- Runtime system for heterogeneous multicore architectures.
 *
 * Copyright (C) 2010-2012  Université de Bordeaux 1
 * Copyright (C) 2010, 2011, 2012  Centre National de la Recherche Scientifique
 * Copyright (C) 2011, 2012  INRIA
 *
 * StarPU is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * StarPU is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License in COPYING.LGPL for more details.
 */

/* PaStiX Work stealing policy */

#include <pthread.h>
#ifdef STARPU_USE_DEPRECATED_API
#  undef STARPU_USE_DEPRECATED_API
#endif
#include <starpu.h>
#include <starpu_task_list.h>
#include <starpu_scheduler.h>
#include "common_pastix.h"
#include "queue.h"
#include "sopalin_thread.h"
#ifdef STARPU_USE_FXT
#  include <fxt/fxt.h>
#  include <fxt/fut.h>
#  define _STARPU_FUT_WORK_STEALING        0x5115
#  define _STARPU_TRACE_WORK_STEALING(empty_q, victim_q)		\
  FUT_DO_PROBE2(_STARPU_FUT_WORK_STEALING, empty_q, victim_q)
#else
#  define _STARPU_TRACE_WORK_STEALING(empty_q, victim_q)  do {} while(0);
#endif
struct _pastix_work_stealing_data
{
  struct starpu_task_list ** task_lists;
  int * nprocessed;
  int * njobs;
  unsigned rr_worker;
  /* keep track of the work performed from the beginning of the algorithm to make
   * better decisions about which queue to select when stealing or deferring work
   */
  unsigned performed_total;
  unsigned last_pop_worker;
  unsigned last_push_worker;
};

#ifdef USE_OVERLOAD

/**
 * Minimum number of task we wait for being processed before we start assuming
 * on which worker the computation would be faster.
 */
static int calibration_value = 0;

#endif /* USE_OVERLOAD */


/**
 * Return a worker from which a task can be stolen.
 * Selecting a worker is done in a round-robin fashion, unless
 * the worker previously selected doesn't own any task,
 * then we return the first non-empty worker.
 */
static unsigned select_victim_round_robin(unsigned sched_ctx_id)
{
  struct _pastix_work_stealing_data *pws = (struct _pastix_work_stealing_data*)starpu_sched_ctx_get_policy_data(sched_ctx_id);
  unsigned worker   = pws->last_pop_worker;
  unsigned nworkers = starpu_sched_ctx_get_nworkers(sched_ctx_id);

  /* If the worker's queue is empty, let's try
   * the next ones */
  while (!pws->njobs[worker])
    {
      worker = (worker + 1) % nworkers;
      if (worker == pws->last_pop_worker)
        {
          /* We got back to the first worker,
           * don't go in infinite loop */
          break;
        }
    }

  pws->last_pop_worker = (worker + 1) % nworkers;

  return worker;
}

/**
 * Return a worker to whom add a task.
 * Selecting a worker is done in a round-robin fashion.
 */
static unsigned select_worker_round_robin(unsigned sched_ctx_id)
{
  struct _pastix_work_stealing_data *pws = (struct _pastix_work_stealing_data*)starpu_sched_ctx_get_policy_data(sched_ctx_id);
  unsigned worker = pws->last_push_worker;
  unsigned nworkers = starpu_sched_ctx_get_nworkers(sched_ctx_id);

  pws->last_push_worker = (pws->last_push_worker + 1) % nworkers;

  return worker;
}

#ifdef USE_OVERLOAD

/**
 * Return a ratio helpful to determine whether a worker is suitable to steal
 * tasks from or to put some tasks in its queue.
 *
 * \return	a ratio with a positive or negative value, describing the current state of the worker :
 * 		a smaller value implies a faster worker with an relatively emptier queue : more suitable to put tasks in
 * 		a bigger value implies a slower worker with an reletively more replete queue : more suitable to steal tasks from
 */
static float overload_metric(unsigned sched_ctx_id, unsigned id)
{
  struct _pastix_work_stealing_data *pws = (struct _pastix_work_stealing_data*)starpu_sched_ctx_get_policy_data(sched_ctx_id);
  float execution_ratio = 0.0f;
  float current_ratio = 0.0f;

  int nprocessed = pws->nprocessed[id];
  unsigned njobs = pws->njobs[id];

  /* Did we get enough information ? */
  if (performed_total > 0 && nprocessed > 0)
    {
      /* How fast or slow is the worker compared to the other workers */
      execution_ratio = (float) nprocessed / performed_total;
      /* How replete is its queue */
      current_ratio = (float) njobs / nprocessed;
    }
  else
    {
      return 0.0f;
    }

  return (current_ratio - execution_ratio);
}

/**
 * Return the most suitable worker from which a task can be stolen.
 * The number of previously processed tasks, total and local,
 * and the number of tasks currently awaiting to be processed
 * by the tasks are taken into account to select the most suitable
 * worker to steal task from.
 */
static unsigned select_victim_overload(unsigned sched_ctx_id)
{
  unsigned worker;
  float  worker_ratio;
  unsigned best_worker = 0;
  float best_ratio = FLT_MIN;

  /* Don't try to play smart until we get
   * enough informations. */
  if (performed_total < calibration_value)
    return select_victim_round_robin(sched_ctx_id);

  struct starpu_sched_ctx_worker_collection *workers = starpu_sched_ctx_get_worker_collection(sched_ctx_id);

  struct starpu_iterator it;
  if(workers->init_iterator)
    workers->init_iterator(workers, &it);

  while(workers->has_next(workers, &it))
    {
      worker = workers->get_next(workers, &it);
      worker_ratio = overload_metric(sched_ctx_id, worker);

      if (worker_ratio > best_ratio)
        {
          best_worker = worker;
          best_ratio = worker_ratio;
        }
    }

  return best_worker;
}

/**
 * Return the most suitable worker to whom add a task.
 * The number of previously processed tasks, total and local,
 * and the number of tasks currently awaiting to be processed
 * by the tasks are taken into account to select the most suitable
 * worker to add a task to.
 */
static unsigned select_worker_overload(unsigned sched_ctx_id)
{
  unsigned worker;
  float  worker_ratio;
  unsigned best_worker = 0;
  float best_ratio = FLT_MAX;

  /* Don't try to play smart until we get
   * enough informations. */
  if (performed_total < calibration_value)
    return select_worker_round_robin(sched_ctx_id);

  struct starpu_sched_ctx_worker_collection *workers = starpu_sched_ctx_get_worker_collection(sched_ctx_id);

  struct starpu_iterator it;
  if(workers->init_iterator)
    workers->init_iterator(workers, &it);

  while(workers->has_next(workers, &it))
    {
      worker = workers->get_next(workers, &it);

      worker_ratio = overload_metric(sched_ctx_id, worker);

      if (worker_ratio < best_ratio)
        {
          best_worker = worker;
          best_ratio = worker_ratio;
        }
    }

  return best_worker;
}

#endif /* USE_OVERLOAD */


/**
 * Return a worker from which a task can be stolen.
 * This is a phony function used to call the right
 * function depending on the value of USE_OVERLOAD.
 */
static inline unsigned select_victim(unsigned sched_ctx_id)
{
#ifdef USE_OVERLOAD
  return select_victim_overload(sched_ctx_id);
#else
  return select_victim_round_robin(sched_ctx_id);
#endif /* USE_OVERLOAD */
}

/**
 * Return a worker from which a task can be stolen.
 * This is a phony function used to call the right
 * function depending on the value of USE_OVERLOAD.
 */
static inline unsigned select_worker(unsigned sched_ctx_id)
{
#ifdef USE_OVERLOAD
  return select_worker_overload(sched_ctx_id);
#else
  return select_worker_round_robin(sched_ctx_id);
#endif /* USE_OVERLOAD */
}


#ifdef STARPU_DEVEL
#warning TODO rewrite ... this will not scale at all now
#endif
static struct starpu_task *pws_pop_task(unsigned sched_ctx_id)
{
  struct _pastix_work_stealing_data *pws = (struct _pastix_work_stealing_data*)starpu_sched_ctx_get_policy_data(sched_ctx_id);

  struct starpu_task *task;
  struct starpu_task_list *q;

  int workerid = starpu_worker_get_id();

  STARPU_ASSERT(workerid != -1);

  q = pws->task_lists[workerid];

  task = starpu_task_list_pop_back(q);
  if (task)
    {
      /* there was a local task */
      pws->performed_total++;
      pws->nprocessed[workerid]++;
      pws->njobs[workerid]--;
      return task;
    }

  pthread_mutex_t *worker_sched_mutex;
  pthread_cond_t  *worker_sched_cond;
  starpu_worker_get_sched_condition(workerid, &worker_sched_mutex, &worker_sched_cond);
  MUTEX_UNLOCK(worker_sched_mutex);


  /* we need to steal someone's job */
  unsigned victim = select_victim(sched_ctx_id);

  pthread_mutex_t *victim_sched_mutex;
  pthread_cond_t *victim_sched_cond;

  starpu_worker_get_sched_condition(victim, &victim_sched_mutex, &victim_sched_cond);
  MUTEX_LOCK(victim_sched_mutex);
  struct starpu_task_list *victimq = pws->task_lists[victim];

  task = starpu_task_list_pop_back(victimq);
  if (task)
    {
      _STARPU_TRACE_WORK_STEALING(q, workerid);
      pws->performed_total++;

      /* Beware : we have to increase the number of processed tasks of
       * the stealer, not the victim ! */
      pws->nprocessed[workerid]++;
      pws->njobs[victim]--;
    }

  MUTEX_UNLOCK(victim_sched_mutex);

  MUTEX_LOCK(worker_sched_mutex);
  if(!task)
    {
      task = starpu_task_list_pop_back(q);
      if (task)
        {
          /* there was a local task */
          pws->performed_total++;
          pws->nprocessed[workerid]++;
          pws->njobs[workerid]--;
          return task;
        }
    }

  return task;
}

int pws_push_task(struct starpu_task *task)
{
  unsigned sched_ctx_id = task->sched_ctx;
  struct _pastix_work_stealing_data *pws = (struct _pastix_work_stealing_data*)starpu_sched_ctx_get_policy_data(sched_ctx_id);

  struct starpu_task_list *list;

  int workerid = task->workerid;
  if (starpu_worker_get_id() != -1)
    workerid = starpu_worker_get_id();

  int ret_val = -1;

  unsigned worker = 0;
  struct starpu_worker_collection *workers = starpu_sched_ctx_get_worker_collection(sched_ctx_id);
  struct starpu_sched_ctx_iterator it;
  if(workers->init_iterator)
    workers->init_iterator(workers, &it);

  while(workers->has_next(workers, &it))
    {
      worker = workers->get_next(workers, &it);
      pthread_mutex_t *sched_mutex;
      pthread_cond_t *sched_cond;
      starpu_worker_get_sched_condition(worker, &sched_mutex, &sched_cond);
      MUTEX_LOCK(sched_mutex);
    }


  /* If the current thread is not a worker but
   * the main thread (-1), we find the better one to
   * put task on its queue */
  if (workerid == -1)
    workerid = select_worker(sched_ctx_id);

  list = pws->task_lists[workerid];

/* #ifdef HAVE_AYUDAME_H */
/*   if (AYU_event) */
/*     { */
/*       int id = workerid; */
/*       AYU_event(AYU_ADDTASKTOQUEUE, j->job_id, &id); */
/*     } */
/* #endif */

  starpu_task_list_push_back(list, task);
  pws->njobs[workerid]++;

  while(workers->has_next(workers, &it))
    {
      worker = workers->get_next(workers, &it);
      pthread_mutex_t *sched_mutex;
      pthread_cond_t *sched_cond;
      starpu_worker_get_sched_condition(worker, &sched_mutex, &sched_cond);
      pthread_cond_signal(sched_cond);
      MUTEX_UNLOCK(sched_mutex);
    }

  return 0;
}

static void pws_add_workers(unsigned sched_ctx_id, int *workerids,unsigned nworkers)
{
  struct _pastix_work_stealing_data *pws = (struct _pastix_work_stealing_data*)starpu_sched_ctx_get_policy_data(sched_ctx_id);

  unsigned i;
  int workerid;

  for (i = 0; i < nworkers; i++)
    {
      workerid = workerids[i];
      MALLOC_INTERN(pws->task_lists[workerid], 1, struct starpu_task_list);
      starpu_task_list_init(pws->task_lists[workerid]);
      /**
       * The first PWS_POP_TASK will increase NPROCESSED though no task was actually performed yet,
       * we need to initialize it at -1.
       */
      pws->nprocessed[workerid] = -1;
      pws->njobs[workerid] = 0;
    }
}

static void pws_remove_workers(unsigned sched_ctx_id, int *workerids, unsigned nworkers)
{
  struct _pastix_work_stealing_data *pws = (struct _pastix_work_stealing_data*)starpu_sched_ctx_get_policy_data(sched_ctx_id);

  unsigned i;
  int workerid;

  for (i = 0; i < nworkers; i++)
    {
      workerid = workerids[i];
      memFree_null(pws->task_lists[workerid]);
    }
}

static void initialize_pws_policy(unsigned sched_ctx_id)
{
  starpu_sched_ctx_create_worker_collection(sched_ctx_id, STARPU_WORKER_LIST);

  struct _pastix_work_stealing_data *pws;
  MALLOC_INTERN(pws, 1, struct _pastix_work_stealing_data);
  starpu_sched_ctx_set_policy_data(sched_ctx_id, (void*)pws);

  pws->last_pop_worker = 0;
  pws->last_push_worker = 0;

  /**
   * The first PWS_POP_TASK will increase PERFORMED_TOTAL though no task was actually performed yet,
   * we need to initialize it at -1.
   */
  pws->performed_total = -1;

  MALLOC_INTERN(pws->task_lists, STARPU_NMAXWORKERS, struct starpu_task_list *);
  MALLOC_INTERN(pws->nprocessed, STARPU_NMAXWORKERS, int);
  MALLOC_INTERN(pws->njobs,      STARPU_NMAXWORKERS, int);
}

static void deinit_pws_policy(unsigned sched_ctx_id)
{
  struct _pastix_work_stealing_data *pws = (struct _pastix_work_stealing_data*)starpu_sched_ctx_get_policy_data(sched_ctx_id);

  memFree_null(pws->task_lists);
  memFree_null(pws->nprocessed);
  memFree_null(pws->njobs);
  memFree_null(pws);
  starpu_sched_ctx_delete_worker_collection(sched_ctx_id);
}

struct starpu_sched_policy starpu_pastix_sched_policy =
{
  .init_sched = initialize_pws_policy,
  .deinit_sched = deinit_pws_policy,
  .add_workers = pws_add_workers,
  .remove_workers = pws_remove_workers,
  .push_task = pws_push_task,
  .pop_task = pws_pop_task,
  .pre_exec_hook = NULL,
  .post_exec_hook = NULL,
  .pop_every_task = NULL,
  .policy_name = "pws",
  .policy_description = "work stealing"
};
#endif /* WITH_STARPU */
