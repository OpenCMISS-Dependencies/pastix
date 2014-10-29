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
 * Updown step written using StarPU.
 *
 */
#ifdef WITH_STARPU
#  ifdef STARPU_USE_DEPRECATED_API
#    undef STARPU_USE_DEPRECATED_API
#  endif
#  include <starpu.h>
#  ifdef FORCE_NOMPI
#    include "nompi.h"
#  else
#    include <mpi.h>
#  endif
#  include "common_pastix.h"
#  include "symbol.h"
#  include "ftgt.h"
#  include "csc.h"
#  include "updown.h"
#  include "queue.h"
#  include "bulles.h"
#  include "solver.h"
#  include "sopalin_thread.h"
#  include "sopalin_define.h"
#  include "sopalin3d.h"
#  include "sopalin_acces.h"
#  include "starpu_updo.h"
#  include "starpu_updo_kernels.h"
#  define USE_TASK_DEP
int starpu_register_sm2x(Sopalin_Data_t       * sopalin_data,
                         starpu_data_handle_t * SM2X_handles)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  int itercblk;
  for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++)
    {

      starpu_matrix_data_register(&(SM2X_handles[itercblk]), 0,
                                  (uintptr_t)&(UPDOWN_SM2XTAB[
                                                 UPDOWN_SM2XIND(itercblk)]),
                                  (uint32_t)UPDOWN_SM2XSZE,
                                  CBLK_COLNBR(itercblk),
                                  (uint32_t)UPDOWN_SM2XNBR,
                                  sizeof(PASTIX_FLOAT));
    }
  return NO_ERR;
}

#  define updown_TRSM_model API_CALL(updown_TRSM_model)
static struct starpu_perfmodel updo_TRSM_model = {
  .type = STARPU_HISTORY_BASED,
  .symbol = "updo_TRSM",
};

#  define updown_GEMM_model API_CALL(updown_GEMM_model)
static struct starpu_perfmodel updo_GEMM_model = {
  .type = STARPU_HISTORY_BASED,
  .symbol = "updo_GEMM",
};

#  define updown_DIAG_model API_CALL(updown_DIAG_model)
static struct starpu_perfmodel updo_DIAG_model = {
  .type = STARPU_HISTORY_BASED,
  .symbol = "updo_DIAG",
};

#  define updo_trsm_cl API_CALL(updo_trsm_cl)
static struct starpu_codelet updo_trsm_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = updo_trsm_starpu_cpu,
  .model = &updo_TRSM_model,
  .nbuffers = 2,
  .modes = {
    STARPU_R,
    STARPU_RW}
};

#  define updo_up_gemm_cl API_CALL(updo_up_gemm_cl)
static struct starpu_codelet updo_up_gemm_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = updo_up_gemm_starpu_cpu,
  .model = &updo_GEMM_model,
  .nbuffers = 3,
  .modes = {
    STARPU_R,
    STARPU_R,
    STARPU_RW
  }
};

#  define updo_down_gemm_cl API_CALL(updo_down_gemm_cl)
static struct starpu_codelet updo_down_gemm_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = updo_down_gemm_starpu_cpu,
  .model = &updo_GEMM_model,
  .nbuffers = 3,
  .modes = {
    STARPU_R,
    STARPU_R,
    STARPU_RW
  }
};

#  define updo_diag_cl API_CALL(updo_diag_cl)
struct starpu_codelet updo_diag_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = updo_diag_starpu_cpu,
  .model = &updo_DIAG_model,
  .nbuffers = 2,
  .modes = {
    STARPU_R,
    STARPU_RW}
};

#  define DOWN 0
#  define UP   1

static inline
int starpu_submit_up_or_down(Sopalin_Data_t * sopalin_data,
                             starpu_data_handle_t * cblk_handles,
                             starpu_data_handle_t * SM2X_handles,
                             starpu_updo_trsm_data_t * trsm_args,
                             starpu_updo_gemm_data_t * gemm_args,
                             Queue cblreadyqueue,
                             struct starpu_task ** tasktab,
                             int DOWN_OR_UP,
                             int * sched_ctxs)
{
  PASTIX_INT ii;
  SolverMatrix * datacode = sopalin_data->datacode;
#  ifdef USE_TASK_DEP
  PASTIX_INT firstupdotask = 0;
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    firstupdotask = SYMB_BLOKNBR;
  if (DOWN_OR_UP == UP)
    {
      firstupdotask += SYMB_BLOKNBR;
#    ifndef CHOL_SOPALIN
      firstupdotask += SYMB_CBLKNBR;
#    endif  /* CHOL_SOPALIN */
    }
#  endif /* USE_TASK_DEP */
  for (ii=0;ii<SYMB_CBLKNBR;ii++)
    {
      int itercblk;
      int iterblok;
      itercblk = queueGet(&cblreadyqueue);
      /* Task bloc itercblk */
      {
        int ret;
        struct starpu_task * task_trsm;
        task_trsm = starpu_task_create();
        trsm_args[itercblk].cblknum      = itercblk;
        trsm_args[itercblk].sopalin_data = sopalin_data;
        if (DOWN_OR_UP == DOWN)
          {
            trsm_args[itercblk].transpose    = 'N';
#  if (defined CHOL_SOPALIN) && (!defined SOPALIN_LU)
            trsm_args[itercblk].diag         = 'N';
#  else
            trsm_args[itercblk].diag         = 'U';
#  endif

          }
        else
          {
#  ifdef CHOL_SOPALIN
            trsm_args[itercblk].diag         = 'N';
            trsm_args[itercblk].transpose    = 'T';
#  else
            trsm_args[itercblk].diag         = 'U';
#    ifdef HERMITIAN
            trsm_args[itercblk].transpose    = 'C';
#    else
            trsm_args[itercblk].transpose    = 'T';
#    endif
#  endif
          }

        task_trsm->cl         = &updo_trsm_cl;
        task_trsm->cl_arg     = &(trsm_args[itercblk]);
        task_trsm->handles[0] = cblk_handles[itercblk];
        task_trsm->handles[1] = SM2X_handles[itercblk];
#ifdef STARPU_PASTIX_SCHED
        task_trsm->workerid = TASK_THREADID(itercblk);
        task_trsm->priority = 0;//TASK_PRIONUM(itercblk);
#endif

#  ifdef USE_TASK_DEP
        task_trsm->destroy = 0;
        tasktab[firstupdotask + SYMB_BLOKNUM(itercblk)] = task_trsm;
        if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
          {
            struct starpu_task * deps[1];
            deps[0] = tasktab[SYMB_BLOKNUM(itercblk)];
            starpu_task_declare_deps_array(task_trsm, 1, deps);
          }
#  endif /* USE_TASK_DEP */
#  ifdef STARPU_GET_TASK_CTX
        ret = starpu_task_submit_to_ctx(task_trsm,
                                        sched_ctxs[0]);
#  else
        ret = starpu_task_submit(task_trsm);
#  endif
        STARPU_ASSERT(!ret);
      }
      if (DOWN_OR_UP == DOWN)
        {
          for (iterblok = SYMB_BLOKNUM(itercblk)+1;
               iterblok < SYMB_BLOKNUM(itercblk+1);
               iterblok++)
            {
              int ret;
              struct starpu_task * task_gemm;
              task_gemm = starpu_task_create();
              gemm_args[iterblok].cblknum = itercblk;
              gemm_args[iterblok].bloknum = iterblok;
              gemm_args[iterblok].sopalin_data = sopalin_data;
              gemm_args[iterblok].transpose = 'N';
              task_gemm->cl = &updo_down_gemm_cl;
              task_gemm->cl_arg = &(gemm_args[iterblok]);
              task_gemm->handles[0] = cblk_handles[itercblk];
              task_gemm->handles[1] = SM2X_handles[itercblk];
              task_gemm->handles[2] = SM2X_handles[SYMB_CBLKNUM(iterblok)];
#ifdef STARPU_PASTIX_SCHED
              task_gemm->workerid = TASK_THREADID(itercblk);
              task_gemm->priority = 0;//TASK_PRIONUM(itercblk);
#endif
#  ifdef USE_TASK_DEP
              task_gemm->destroy = 0;
              tasktab[firstupdotask + iterblok] = task_gemm;
              if (sopalin_data->sopar->iparm[IPARM_START_TASK] <=
                  API_TASK_NUMFACT)
                {
                  struct starpu_task * deps[1];
                  deps[0] = tasktab[SYMB_BLOKNUM(itercblk)];
                  starpu_task_declare_deps_array(task_gemm, 1, deps);
                }
#  endif /* USE_TASK_DEP */
#  ifdef STARPU_GET_TASK_CTX
              ret = starpu_task_submit_to_ctx(task_gemm,
                                              sched_ctxs[0]);
#  else
              ret = starpu_task_submit(task_gemm);
#  endif
              STARPU_ASSERT(!ret);
            }
        }
      else
        {
          PASTIX_INT proc;
          for (proc=UPDOWN_BROWPROCNBR(itercblk)-1;proc>=0;proc--)
            {
              if (UPDOWN_BROWPROCTAB(itercblk)[proc] != SOLV_PROCNUM)
                {
                  assert(0);
                }
              else
                {
                  PASTIX_INT count;
                  PASTIX_INT listptridx = UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(itercblk));
                  for (count=UPDOWN_LISTPTR(listptridx);
                       count<UPDOWN_LISTPTR(listptridx+1);
                       count++)
                    {
                      PASTIX_INT cblk;
                      int ret;
                      struct starpu_task * task_gemm;

                      task_gemm = starpu_task_create();
                      cblk      = UPDOWN_LISTCBLK(count);
                      iterblok  = UPDOWN_LISTBLOK(count);

                      ASSERTDBG((SYMB_FROWNUM(iterblok)>=
                                 SYMB_FCOLNUM(itercblk)) &&
                                (SYMB_LROWNUM(iterblok)<=
                                 SYMB_LCOLNUM(itercblk)),
                                MOD_SOPALIN);

                      gemm_args[iterblok].cblknum = cblk;

                      gemm_args[iterblok].bloknum = iterblok;
                      gemm_args[iterblok].sopalin_data = sopalin_data;

#  ifdef HERMITIAN
                      gemm_args[iterblok].transpose = 'C';
#  else  /* not HERMITIAN */
                      gemm_args[iterblok].transpose = 'T';
#  endif /* not HERMITIAN */
                      task_gemm->cl = &updo_up_gemm_cl;
                      task_gemm->cl_arg = &(gemm_args[iterblok]);
                      task_gemm->handles[0] = cblk_handles[cblk];
                      task_gemm->handles[1] = SM2X_handles[cblk];
                      task_gemm->handles[2] = SM2X_handles[
                        SYMB_CBLKNUM(iterblok)];
#ifdef STARPU_PASTIX_SCHED
                      task_gemm->workerid = TASK_THREADID(itercblk);
                      task_gemm->priority = 0;//TASK_PRIONUM(itercblk);
#endif
#  ifdef USE_TASK_DEP
                      task_gemm->destroy = 0;
                      tasktab[firstupdotask + iterblok] = task_gemm;

                      if (sopalin_data->sopar->iparm[IPARM_START_TASK] <=
                          API_TASK_NUMFACT)
                        {
                          struct starpu_task * deps[1];
                          deps[0] = tasktab[SYMB_BLOKNUM(cblk)];
                          starpu_task_declare_deps_array(task_gemm, 1, deps);
                        }
#  endif /* USE_TASK_DEP */
#  ifdef STARPU_GET_TASK_CTX
                      ret = starpu_task_submit_to_ctx(task_gemm,
                                                      sched_ctxs[0]);
#  else
                      ret = starpu_task_submit(task_gemm);
#  endif
                      STARPU_ASSERT(!ret);
                    }
                }


            }
        }
    }
  return NO_ERR;
}

int starpu_submit_updown(Sopalin_Data_t * sopalin_data,
                         starpu_data_handle_t * L_handles,
                         starpu_data_handle_t * U_handles,
                         starpu_data_handle_t * SM2X_handles,
                         struct starpu_task  ** tasktab,
                         int                  * sched_ctxs)
{
  SolverMatrix * datacode = sopalin_data->datacode;
  Queue cblreadyqueue;
  starpu_updo_trsm_data_t * trsm_args;
  starpu_updo_gemm_data_t * gemm_args;
  starpu_updo_diag_data_t * diag_args;
  PASTIX_INT   itercblk;

  MALLOC_INTERN(trsm_args,    2*SYMB_CBLKNBR, starpu_updo_trsm_data_t);
  MALLOC_INTERN(gemm_args,    2*SYMB_BLOKNBR, starpu_updo_gemm_data_t);
  MALLOC_INTERN(diag_args,    SYMB_CBLKNBR, starpu_updo_diag_data_t);

  queueInit(&cblreadyqueue,SYMB_CBLKNBR);
  for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++)
    {
      queueAdd(&cblreadyqueue,itercblk,(double)(TASK_PRIONUM(itercblk)));
    }

  /* Down step */
  starpu_submit_up_or_down(sopalin_data,
                           L_handles,
                           SM2X_handles,
                           trsm_args,
                           gemm_args,
                           cblreadyqueue, tasktab, DOWN, sched_ctxs);

  /* Diag step */
#  ifndef CHOL_SOPALIN
#    ifdef USE_TASK_DEP
  {
    PASTIX_INT firstupdotask;
    firstupdotask = SYMB_BLOKNBR;
    if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
      firstupdotask += SYMB_BLOKNBR;
#    endif
    for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++)
      {
        int ret;
        struct starpu_task * task_diag;
        diag_args[itercblk].cblknum = itercblk;
        diag_args[itercblk].sopalin_data = sopalin_data;
        task_diag = starpu_task_create();
        task_diag->cl = &updo_diag_cl;
        task_diag->cl_arg = &(diag_args[itercblk]);
        task_diag->handles[0] = L_handles[itercblk];
        task_diag->handles[1] = SM2X_handles[itercblk];
#ifdef STARPU_PASTIX_SCHED
        task_diag->workerid = TASK_THREADID(itercblk);
        task_diag->priority = 0;//TASK_PRIONUM(itercblk);
#endif

#    ifdef USE_TASK_DEP
        task_diag->destroy = 0;
        tasktab[firstupdotask + itercblk] = task_diag;

        if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
          {
            struct starpu_task * deps[1];
            deps[0] = tasktab[SYMB_BLOKNUM(itercblk)];
            starpu_task_declare_deps_array(task_diag, 1, deps);
          }
#    endif /* USE_TASK_DEP */
#  ifdef STARPU_GET_TASK_CTX
        ret = starpu_task_submit_to_ctx(task_diag,
                                        sched_ctxs[0]);
#  else
        ret = starpu_task_submit(task_diag);
#  endif
        STARPU_ASSERT(!ret);
      }
  }
#  endif /* not CHOL_SOPALIN */
  /* Up step */
  queueInit(&cblreadyqueue,SYMB_CBLKNBR);
  for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++)
    {
      queueAdd(&cblreadyqueue,itercblk,-(double)(TASK_PRIONUM(itercblk)));
    }

  starpu_submit_up_or_down(sopalin_data,
#  ifdef SOPALIN_LU
                           U_handles,
#  else /* not SOPALIN_LU */
                           L_handles,
#  endif /* not SOPALIN_LU */
                           SM2X_handles,
                           trsm_args+SYMB_CBLKNBR,
                           gemm_args+SYMB_BLOKNBR,
                           cblreadyqueue, tasktab, UP,
                           sched_ctxs);
  return NO_ERR;
}
#else
/* ISO C forbids an empty source file */
#  include "not_empty.h"
NOT_EMPTY(starpu_updo)
#endif /* WITH_STARPU */
