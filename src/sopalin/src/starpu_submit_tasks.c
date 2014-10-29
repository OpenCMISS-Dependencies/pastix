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
#ifdef WITH_STARPU
#  ifndef FORCE_NO_CUDA
#    include <cuda.h>
#  endif
#  ifdef STARPU_USE_DEPRECATED_API
#    undef STARPU_USE_DEPRECATED_API
#  endif
#  include <starpu.h>
#  include <starpu_profiling.h>
#  include <pthread.h>
#  include <string.h>
#  ifdef FORCE_NOMPI
#    include "nompi.h"
#  else
#    include <mpi.h>
#  endif
#  include "common_pastix.h"
#  include "out.h"
#  include "sopalin_define.h"
#  include "sopalin_acces.h"
#  include "symbol.h"
#  include "ftgt.h"
#  include "csc.h"
#  include "updown.h"
#  include "queue.h"
#  include "bulles.h"
#  include "solver.h"
#  include "sopalin_thread.h"
#  include "sopalin_time.h"
#  include "sopalin3d.h"
#  include "starpu_kernels.h"
#  include "../../perf/src/perf.h"
#  include "starpu_submit_tasks.h"
#  include "starpu_updo.h"
#  include "starpu_pastix_sched_policy.h"
#  include "sopalin_init.h"
#define dump_all API_CALL(dump_all)
void  dump_all                 (SolverMatrix *, CscMatrix * cscmtx, int);

#  ifdef STARPU_USE_CUDA
#    if ((!defined PREC_DOUBLE) || (!(defined __CUDA_ARCH__) || __CUDA_ARCH__ >= 130))
#      if !(defined PREC_DOUBLE && defined TYPE_COMPLEX && CUDA_SM_VERSION < 20)
#        ifndef FORCE_NO_CUDA
#          define STARPU_USE_CUDA_GEMM_FUNC
#        endif
#      endif
#    endif
#  endif

#  ifdef TYPE_COMPLEX
#    ifdef PREC_DOUBLE
#      define PREFIX  "Z"
#    else
#      define PREFIX  "C"
#    endif
#  else
#    ifdef PREC_DOUBLE
#      define PREFIX  "D"
#    else
#      define PREFIX  "S"
#    endif
#  endif

#  define CUDA_CALL(x) do                                       \
    {                                                           \
      if (cudaSuccess != x)                                     \
        {                                                       \
          errorPrint("%s (%s,%d)\n",#x, __FILE__,__LINE__);     \
          assert(0);                                            \
        }                                                       \
    } while(0)


#  define USE_TASK_DEP

#  if (STARPU_MAJOR_VERSION < 1)
#    error "PaStiX requires STARPU >= 1"
#  endif /* STARPU_MAJOR_VERSION */

#  if (STARPU_MAJOR_VERSION == 1 &&  STARPU_MINOR_VERSION == 0)
/* 1.1 : starpu_perf_archtype => starpu_perfmodel_archtype
 *       starpu_profiling_worker_info => starpu_profiling_worker_info
 */
#    define starpu_perfmodel_archtype        starpu_perf_archtype
#    define starpu_profiling_worker_info     starpu_worker_profiling_info
#    define starpu_profiling_task_info       starpu_task_profiling_info
#    define starpu_profiling_worker_get_info starpu_worker_get_profiling_info
#    define starpu_profiling_set_id          starpu_set_profiling_id
#  endif

#  define ARCH_CPU  0
#  define ARCH_CUDA 1

#define starpu_id API_CALL(starpu_id)
static int starpu_id = 0;

static size_t trf_size(struct starpu_task *task,
                       enum starpu_perfmodel_archtype arch,
                       unsigned nimpl)
{
  starpu_trf_data_t * args         = (starpu_trf_data_t*)task->cl_arg;
  Sopalin_Data_t    * sopalin_data = args->sopalin_data;
  SolverMatrix      * datacode     = sopalin_data->datacode;
  PASTIX_INT                 stride       = STARPU_MATRIX_GET_LD(task->handles[0]);
  size_t              dima         = CBLK_COLNBR(args->cblknum);
  return OPS_PPF(dima) + OPS_TRSM(dima, stride);

}

static size_t gemm_size(struct starpu_task *task,
                        enum starpu_perfmodel_archtype arch,
                        unsigned nimpl)
{
  starpu_gemm_data_t         * args         = (starpu_gemm_data_t*)task->cl_arg;
  Sopalin_Data_t             * sopalin_data = args->sopalin_data;
  PASTIX_INT                   bloknum      = args->bloknum;
  SolverMatrix               * datacode     = sopalin_data->datacode;
  PASTIX_INT                   stride       = SOLV_STRIDE(args->cblknum);
  PASTIX_INT                   indblok      = SOLV_COEFIND(bloknum);
  PASTIX_INT                   dimi         = stride - indblok;
  PASTIX_INT                   dimj         = BLOK_ROWNBR(bloknum);
  PASTIX_INT                   dima         = CBLK_COLNBR(args->cblknum);

  return OPS_GEMM(dimi,dimj,dima);
}

static struct starpu_perfmodel GEMM_model = {
  .type = STARPU_HISTORY_BASED,
  .per_arch =
  {
    [STARPU_CPU_DEFAULT][0]  = { .size_base = gemm_size },
    [STARPU_CUDA_DEFAULT][0] = { .size_base = gemm_size }
  },
  .symbol = PREFIX "GEMM"
};

#  define starpu_partition_data                  API_CALL(starpu_partition_data)


#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
/* LU */
static struct starpu_perfmodel GETRF_TRSM_model = {
  .type = STARPU_REGRESSION_BASED,
  .symbol = PREFIX "GETRF_TRSM",
  .per_arch =
  {
    [STARPU_CPU_DEFAULT][0]  = { .size_base = trf_size },
    [STARPU_CUDA_DEFAULT][0] = { .size_base = trf_size }
  }
};

#      define getrfsp1d_cl             PASTIX_PREFIX_F(getrfsp1d_cl)
#      define getrfsp1d_gemm_cl        PASTIX_PREFIX_F(getrfsp1d_gemm_cl)
#      define getrfsp1d_sparse_gemm_cl PASTIX_PREFIX_F(getrfsp1d_sparse_gemm_cl)
#      define getrfsp1d_sparse_gemm_cpu_cl PASTIX_PREFIX_F(getrfsp1d_sparse_gemm_cpu_cl)

struct starpu_codelet getrfsp1d_cl =
{
  .where = STARPU_CPU
#ifndef FORCE_NO_CUDA
#  ifdef WITH_MAGMABLAS
  |STARPU_CUDA
#  endif /* WITH_MAGMABLAS */
#endif /* not FORCE_NO_CUDA */
  ,
  .cpu_funcs[0] = getrfsp1d_starpu_cpu,
#ifndef FORCE_NO_CUDA
#  ifdef STARPU_USE_CUDA
#    ifdef WITH_MAGMABLAS
  .cuda_funcs[0] = getrfsp1d_starpu_cuda,
#    endif /* WITH_MAGMABLAS */
#  endif /* STARPU_USE_CUDA */
#endif /* not FORCE_NO_CUDA */
  .model = &GETRF_TRSM_model,
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_RW}
};
struct starpu_codelet getrfsp1d_sparse_gemm_cpu_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = getrfsp1d_sparse_gemm_starpu_cpu,
  .model = &GEMM_model,
#      if (defined STARPU_BLOCKTAB_SELFCOPY)
  .nbuffers = 5,
#      else
  .nbuffers = 6,
#      endif
  .modes = {STARPU_R,
            STARPU_RW,
            STARPU_R,
            STARPU_RW,
            STARPU_SCRATCH,
#      if !(defined STARPU_BLOCKTAB_SELFCOPY)
            STARPU_R
#      endif
  }
};

struct starpu_codelet getrfsp1d_sparse_gemm_cl =
{
  .where = STARPU_CPU
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  |STARPU_CUDA
#      endif
  ,
  .cpu_funcs[0] = getrfsp1d_sparse_gemm_starpu_cpu,
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  .cuda_funcs[0] = getrfsp1d_sparse_gemm_starpu_cuda,
#      endif
  .model = &GEMM_model,
#      if (defined STARPU_BLOCKTAB_SELFCOPY)
  .nbuffers = 5,
#      else
  .nbuffers = 6,
#      endif
  .modes = {STARPU_R,
            STARPU_RW,
            STARPU_R,
            STARPU_RW,
            STARPU_SCRATCH,
#      if !(defined STARPU_BLOCKTAB_SELFCOPY)
            STARPU_R
#      endif
  }
};
#      define cl_trf       getrfsp1d_cl
#      define cl_gemm      getrfsp1d_sparse_gemm_cl
#      define cl_cpu_gemm      getrfsp1d_sparse_gemm_cpu_cl
#    else /* SOPALIN_LU */
/* LLT */
static struct starpu_perfmodel POTRF_TRSM_model = {
  .type = STARPU_REGRESSION_BASED,
  .symbol = PREFIX "POTRF_TRSM",
  .per_arch =
  {
    [STARPU_CPU_DEFAULT][0]  = { .size_base = trf_size },
    [STARPU_CUDA_DEFAULT][0] = { .size_base = trf_size }
  }

};
#      define potrfsp1d_cl             PASTIX_PREFIX_F(potrfsp1d_cl)
#      define potrfsp1d_gemm_cl        PASTIX_PREFIX_F(potrfsp1d_gemm_cl)
#      define potrfsp1d_sparse_gemm_cl PASTIX_PREFIX_F(potrfsp1d_sparse_gemm_cl)
#      define potrfsp1d_sparse_gemm_cpu_cl PASTIX_PREFIX_F(potrfsp1d_sparse_gemm_cpu_cl)
struct starpu_codelet potrfsp1d_cl =
{
  .where = STARPU_CPU
#ifndef FORCE_NO_CUDA
#  ifdef WITH_MAGMABLAS
  |STARPU_CUDA
#  endif
#endif /* not FORCE_NO_CUDA */
  ,
  .cpu_funcs[0] = potrfsp1d_starpu_cpu,
#ifndef FORCE_NO_CUDA
#  ifdef WITH_MAGMABLAS
  .cuda_funcs[0] = getrfsp1d_starpu_cuda,
#  endif /* WITH_MAGMABLAS */
#endif /* not FORCE_NO_CUDA */
  .model = &POTRF_TRSM_model,
  .nbuffers = 1,
  .modes = {
    STARPU_RW}
};

struct starpu_codelet potrfsp1d_sparse_gemm_cl =
{
  .where = STARPU_CUDA
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  | STARPU_CUDA
#      endif
  ,
  .cpu_funcs[0] = potrfsp1d_sparse_gemm_starpu_cpu,
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  .cuda_funcs[0] = potrfsp1d_sparse_gemm_starpu_cuda,
#      endif
  .model = &GEMM_model,
#      if (defined STARPU_BLOCKTAB_SELFCOPY)
  .nbuffers = 3,
#      else
  .nbuffers = 4,
#      endif

  .modes = {STARPU_R,
            STARPU_RW,
            STARPU_SCRATCH,
#      if !(defined STARPU_BLOCKTAB_SELFCOPY)
            STARPU_R
#      endif
  }
};

struct starpu_codelet potrfsp1d_sparse_gemm_cpu_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = potrfsp1d_sparse_gemm_starpu_cpu,
  .model = &GEMM_model,
#      if (defined STARPU_BLOCKTAB_SELFCOPY)
  .nbuffers = 3,
#      else
  .nbuffers = 4,
#      endif

  .modes = {STARPU_R,
            STARPU_RW,
            STARPU_SCRATCH,
#      if !(defined STARPU_BLOCKTAB_SELFCOPY)
            STARPU_R
#      endif
  }
};
#      define cl_trf       potrfsp1d_cl
#      define cl_gemm      potrfsp1d_sparse_gemm_cl
#      define cl_cpu_gemm      potrfsp1d_sparse_gemm_cpu_cl
#    endif /* SOPALIN_LU */
#  else  /* CHOL_SOPALIN */
/* LDLT */
static struct starpu_perfmodel HETRF_TRSM_model = {
  .type = STARPU_REGRESSION_BASED,
  .symbol = PREFIX "HETRF_TRSM",
  .per_arch =
  {
    [STARPU_CPU_DEFAULT][0]  = { .size_base = trf_size },
    [STARPU_CUDA_DEFAULT][0] = { .size_base = trf_size }
  }

};
#    ifdef HERMITIAN
#      define hetrfsp1d_cl             PASTIX_PREFIX_F(hetrfsp1d_cl)
#      define hetrfsp1d_sparse_gemm_cl PASTIX_PREFIX_F(hetrfsp1d_sparse_gemm_cl)
#      define hetrfsp1d_sparse_gemm_cpu_cl PASTIX_PREFIX_F(hetrfsp1d_sparse_gemm_cpu_cl)
#    else
#      define hetrfsp1d_cl             PASTIX_PREFIX_F(sytrfsp1d_cl)
#      define hetrfsp1d_sparse_gemm_cl PASTIX_PREFIX_F(sytrfsp1d_sparse_gemm_cl)
#      define hetrfsp1d_sparse_gemm_cpu_cl PASTIX_PREFIX_F(sytrfsp1d_sparse_gemm_cpu_cl)
#    endif
struct starpu_codelet hetrfsp1d_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = hetrfsp1d_starpu_cpu,
  .model = &HETRF_TRSM_model,
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_SCRATCH}
};

struct starpu_codelet hetrfsp1d_sparse_gemm_cl =
{
  .where = STARPU_CPU
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  | STARPU_CUDA
#      endif
,
  .cpu_funcs[0] = hetrfsp1d_gemm_starpu_cpu,
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  .cuda_funcs[0] = hetrfsp1d_gemm_starpu_cuda,
#      endif
  .model = &GEMM_model,
#    if (defined STARPU_BLOCKTAB_SELFCOPY)
  .nbuffers = 3,
#    else
  .nbuffers = 4,
#    endif
  .modes = { STARPU_R,
             STARPU_RW,
             STARPU_SCRATCH,
#    if !(defined STARPU_BLOCKTAB_SELFCOPY)
             STARPU_R
#    endif
  }
};
struct starpu_codelet hetrfsp1d_sparse_gemm_cpu_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = hetrfsp1d_gemm_starpu_cpu,
  .model = &GEMM_model,
#    if (defined STARPU_BLOCKTAB_SELFCOPY)
  .nbuffers = 3,
#    else
  .nbuffers = 4,
#    endif
  .modes = { STARPU_R,
             STARPU_RW,
             STARPU_SCRATCH,
#    if !(defined STARPU_BLOCKTAB_SELFCOPY)
             STARPU_R
#    endif
  }
};
#    define cl_trf       hetrfsp1d_cl
#    define cl_gemm      hetrfsp1d_sparse_gemm_cl
#    define cl_cpu_gemm  hetrfsp1d_sparse_gemm_cpu_cl
#  endif /* CHOL_SOPALIN */



/*
 Function: starpu_data_partition

 Initialize column blocks handlers.

 Parameters:
 sopalin_data - PaStiX global data structure.
 L_handle     - Handles for L column blocks.
 U_handle     - Handles for U column blocks.
 */
static void starpu_partition_data(Sopalin_Data_t * sopalin_data,
                                  starpu_data_handle_t * L_handle,
                                  starpu_data_handle_t * U_handle)
{
  SolverMatrix       * datacode         = sopalin_data->datacode;
  PASTIX_INT itercblk;

  for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++)
    {
      starpu_matrix_data_register(&(L_handle[itercblk]), 0,
                                  (uintptr_t)SOLV_COEFTAB(itercblk),
                                  (uint32_t)SOLV_STRIDE(itercblk),
                                  (uint32_t)SOLV_STRIDE(itercblk),
                                  CBLK_COLNBR(itercblk),
                                  sizeof(PASTIX_FLOAT));
#  if defined USE_TASK_DEP
      starpu_data_set_sequential_consistency_flag(L_handle[itercblk], 0);
#  endif

#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
      starpu_matrix_data_register(&(U_handle[itercblk]), 0,
                                  (uintptr_t)SOLV_UCOEFTAB(itercblk),
                                  (uint32_t)SOLV_STRIDE(itercblk),
                                  (uint32_t)SOLV_STRIDE(itercblk),
                                  CBLK_COLNBR(itercblk),
                                  sizeof(PASTIX_FLOAT));
#      if (defined USE_TASK_DEP)
      starpu_data_set_sequential_consistency_flag(U_handle[itercblk], 0);
#      endif
#    endif
#  endif
    }
}

/*
 * Function: starpu_init_smp
 *
 * Initialize thread data structure for factorization when using StarPU.
 */
#  define starpu_init_smp API_CALL(starpu_init_smp)
void*
starpu_init_smp (void * arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  int init;
  init = INIT_COMPUTE;
  if (THREAD_FUNNELED_OFF)
    {
      init = init | INIT_SEND;
      if (THREAD_COMM_OFF)
        {
          init = init | INIT_RECV;
        }
    }
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    {
      sopalin_init_smp(sopalin_data, argument->me, API_YES, init);
    }
  else
    {
      sopalin_init_smp(sopalin_data, argument->me, API_NO, init);
    }
  return NULL;
}

#define starpu_init_kernel API_CALL(starpu_init_kernel)
void starpu_init_kernel(void * buffers[], void * _args)
{
  starpu_init_smp(_args);
}
#define sopalin_init_cl API_CALL(sopalin_init_cl)
struct starpu_codelet sopalin_init_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = starpu_init_kernel,
  .nbuffers = 0,
  .modes = {}
};

/*
 Struct: sopthread_data

 Structure conotaining the thread number and a pointer to data given to thread in parameters.
 */
struct starpu_loop_data_ {
  int                    me;
  starpu_data_handle_t * L_handle;
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
  starpu_data_handle_t * U_handle;
#  endif
  starpu_data_handle_t   WORK_handle;
#  if (defined STARPU_BLOCKTAB_SELFCOPY)
  int                 ** d_blocktab;
#  else
  starpu_data_handle_t   blocktab_handle;
#  endif
  starpu_trf_data_t    * trf_args;
  starpu_gemm_data_t   * gemm_args;
  struct starpu_task  ** starpu_tasktab;
  Sopalin_Data_t       * sopalin_data;                 /*+ Data given to thread as argument +*/
  int                    ctx;
  int                    ctx_nbr;
  int                    thread_per_ctx;
  int                  * sched_ctxs;
  int                    first, last;
  int                    ndiag_submitted;
  int                    ngemm_submitted;
#  if (defined USE_TASK_DEP)
  int                  *  deps_tab_size;
  struct starpu_task  *** deps;
#  endif
  pthread_cond_t        * cond_end_facto;
  pthread_mutex_t       * mutex_end_facto;
  int                   * cpu_workerids;
  int                     ncpus;
  int                   * gpu_workerids;
  int                     ngpus;
  int                   * gpu_gemm_count;
};


int
starpu_submit_one_trf(PASTIX_INT itertask, Sopalin_Data_t * sopalin_data)
{
  SolverMatrix        *datacode          = sopalin_data->datacode;
  starpu_loop_data_t  *starpu_loop_data  = (starpu_loop_data_t *)sopalin_data->starpu_loop_data;
  starpu_trf_data_t   *trf_args          = starpu_loop_data->trf_args;
  struct starpu_task **starpu_tasktab    = starpu_loop_data->starpu_tasktab;
  starpu_gemm_data_t  *gemm_args         = starpu_loop_data->gemm_args;
  int                  me                = starpu_loop_data->me;
  int                 *sched_ctxs        = starpu_loop_data->sched_ctxs;

  PASTIX_INT itercblk = TASK_CBLKNUM(itertask);
  PASTIX_INT ret;
  PASTIX_INT iterbloc;
  PASTIX_INT handle_idx;
  struct starpu_task *task_diag;
  int workerid;
#  ifdef STARPU_CONTEXT
  PASTIX_INT threadid = TASK_THREADID(itertask);
  PASTIX_INT my_ctx;
  if (threadid > SOLV_THRDNBR)
    my_ctx = 0;
  else
    my_ctx = 1+threadid/starpu_loop_data->thread_per_ctx;
#  endif
  workerid = starpu_worker_get_id();
  if (workerid == -1) workerid = 0;

  starpu_loop_data->ndiag_submitted++;
  task_diag = starpu_task_create();
  /* We compute diagonal factorisation and TRSM */
  trf_args[itercblk].cblknum      = itercblk;
#  ifdef STARPU_SUBMIT_READY
  trf_args[itercblk].tasknum      = itertask;
#  endif
  trf_args[itercblk].sopalin_data = sopalin_data;

  task_diag->cl = &cl_trf;
  task_diag->cl_arg = &(trf_args[itercblk]);
  task_diag->destroy = 0;
#  ifdef STARPU_PASTIX_SCHED
  if (TASK_THREADID(itertask) < starpu_loop_data->ncpus) {
    task_diag->workerid = starpu_loop_data->cpu_workerids[TASK_THREADID(itertask)];
  } else {
    task_diag->workerid = starpu_loop_data->gpu_workerids[TASK_THREADID(itertask) -
                                                          starpu_loop_data->ncpus];
  }
  task_diag->priority = TASK_PRIONUM(itertask);
#  endif
  handle_idx = 0;

  task_diag->handles[handle_idx++] = starpu_loop_data->L_handle[itercblk];
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  /* LU */
  task_diag->handles[handle_idx++] = starpu_loop_data->U_handle[itercblk];
#    endif
#  else /* CHOL_SOPALIN */
  /* LDLT */
  task_diag->handles[handle_idx++] = starpu_loop_data->WORK_handle;
#  endif

  ASSERT(handle_idx == cl_trf.nbuffers, MOD_SOPALIN);

#  if (defined USE_TASK_DEP)
  {
    PASTIX_INT gcblk2list = UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB( itercblk ));
    PASTIX_INT browk      = ( gcblk2list != -1 )?
      UPDOWN_LISTPTR( gcblk2list    ):-1;
    PASTIX_INT browk1     = ( gcblk2list != -1 )?
      UPDOWN_LISTPTR( gcblk2list + 1):-1;
    PASTIX_INT ndeps, iter;
    ndeps = browk1-browk;

    starpu_tasktab[SYMB_BLOKNUM(itercblk)] = task_diag;
    if (ndeps > starpu_loop_data->deps_tab_size[workerid])
      {
        memFree_null(starpu_loop_data->deps[workerid]);
        starpu_loop_data->deps_tab_size[workerid] = ndeps;
        MALLOC_INTERN(starpu_loop_data->deps[workerid], ndeps, struct starpu_task *);
      }
    for (iter = 0; iter < ndeps; iter++)
      {
        starpu_loop_data->deps[workerid][iter] = starpu_tasktab[UPDOWN_LISTBLOK(browk+iter)];
      }
    if (ndeps != 0)
      starpu_task_declare_deps_array(task_diag, ndeps, starpu_loop_data->deps[workerid]);
  }
#  endif /* USE_TASK_DEP */
#  ifdef STARPU_CONTEXT
  ret = starpu_task_submit_to_ctx(task_diag,
                                  sched_ctxs[my_ctx]);
#  else
  ret = starpu_task_submit(task_diag);
#  endif
#  ifdef STARPU_SUBMIT_READY
  if (itercblk == SYMB_CBLKNBR-1)
  if (SYMB_BLOKNUM(itercblk) == SYMB_BLOKNBR - 1)
    {
      int rc;
      fprintf(stdout, "SIGNAL\n");
      rc = pthread_mutex_lock(starpu_loop_data->mutex_end_facto);
      if (rc) {
        perror("pthread_mutex_lock");
        exit(1);
      }

      rc = pthread_cond_signal(starpu_loop_data->cond_end_facto);
      if (rc) {
        pthread_mutex_unlock(starpu_loop_data->mutex_end_facto);
        perror("pthread_cond_signal");
        exit(1);
      }
      rc = pthread_mutex_unlock(starpu_loop_data->mutex_end_facto);
      if (rc) {
        perror("pthread_mutex_unlock");
        exit(1);
      }

    }
#  endif
  STARPU_ASSERT(!ret);
  return NO_ERR;
}


int
starpu_submit_bunch_of_gemm (PASTIX_INT itertask, Sopalin_Data_t * sopalin_data)
{
  SolverMatrix        *datacode          = sopalin_data->datacode;
  starpu_loop_data_t  *starpu_loop_data  = (starpu_loop_data_t *)sopalin_data->starpu_loop_data;
  starpu_trf_data_t   *trf_args          = starpu_loop_data->trf_args;
  struct starpu_task **starpu_tasktab    = starpu_loop_data->starpu_tasktab;
  starpu_gemm_data_t  *gemm_args         = starpu_loop_data->gemm_args;
  int                  me                = starpu_loop_data->me;
  int                 *sched_ctxs        = starpu_loop_data->sched_ctxs;

  PASTIX_INT itercblk = TASK_CBLKNUM(itertask);
  PASTIX_INT ret;
  PASTIX_INT iterbloc;
  PASTIX_INT handle_idx;
  int workerid;

#  ifdef STARPU_CONTEXT
  PASTIX_INT threadid = TASK_THREADID(itertask);
  PASTIX_INT my_ctx;
  if (threadid > SOLV_THRDNBR)
    my_ctx = 0;
  else
    my_ctx = 1+threadid/starpu_loop_data->thread_per_ctx;
#  endif

  workerid = starpu_worker_get_id();
  if (workerid == -1) workerid = 0;

  for (iterbloc = SYMB_BLOKNUM(itercblk)+1;
       iterbloc < SYMB_BLOKNUM(itercblk+1);
       iterbloc ++)
    {
      struct starpu_task * task_gemm;
      PASTIX_INT blocnbr;
      starpu_loop_data->ngemm_submitted++;
      task_gemm = starpu_task_create();

      blocnbr = SYMB_BLOKNUM(itercblk+1) - iterbloc;
      /* We compute GEMM */
      gemm_args[iterbloc].cblknum      = itercblk;
#  ifdef STARPU_SUBMIT_READY
      gemm_args[iterbloc].tasknum      = itertask;
#  endif
      gemm_args[iterbloc].bloknum      = iterbloc;
      gemm_args[iterbloc].nblocs       = blocnbr;
      gemm_args[iterbloc].fcblknum     = SYMB_CBLKNUM(iterbloc);
      gemm_args[iterbloc].sopalin_data = sopalin_data;
#  if (defined STARPU_BLOCKTAB_SELFCOPY)
      gemm_args[iterbloc].d_blocktab   = starpu_loop_data->d_blocktab;
#  endif

      if ( starpu_loop_data->ngpus == 0 ||
           SOLV_COLOR(SYMB_CBLKNUM(iterbloc)) < 0) {
        task_gemm->cl = &cl_cpu_gemm;
      } else {
        task_gemm->cl = &cl_gemm;
        task_gemm->workerid = starpu_loop_data->gpu_workerids[SOLV_COLOR(SYMB_CBLKNUM(iterbloc))];
        task_gemm->execute_on_a_specific_worker=1;
        starpu_loop_data->gpu_gemm_count[SOLV_COLOR(SYMB_CBLKNUM(iterbloc))]++;
      }
      task_gemm->cl_arg = &(gemm_args[iterbloc]);
      task_gemm->destroy = 0;
/* #ifdef STARPU_PASTIX_SCHED */
/*       fprintf(stdout, "TASK_THREADID(itertask) %d\n", TASK_THREADID(itertask)); */
/*       if (TASK_THREADID(itertask) < starpu_loop_data->ncpus) { */
/*         task_gemm->workerid = starpu_loop_data->cpu_workerids[TASK_THREADID(itertask)]; */
/*       } else { */
/*         task_gemm->workerid = starpu_loop_data->gpu_workerids[TASK_THREADID(itertask) - */
/*                                                               starpu_loop_data->ncpus]; */
/*         task_gemm->execute_on_a_specific_worker=1; */
/*         starpu_loop_data->gpu_gemm_count[TASK_THREADID(itertask) - */
/*                                          starpu_loop_data->ncpus]++; */
/*       } */
/*       task_gemm->priority = TASK_PRIONUM(itertask); */
/* #endif */

      handle_idx = 0;

      task_gemm->handles[handle_idx++] = starpu_loop_data->L_handle[itercblk];
      task_gemm->handles[handle_idx++] = starpu_loop_data->L_handle[SYMB_CBLKNUM(iterbloc)];
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
      task_gemm->handles[handle_idx++] = starpu_loop_data->U_handle[itercblk];
      task_gemm->handles[handle_idx++] = starpu_loop_data->U_handle[SYMB_CBLKNUM(iterbloc)];
#    endif
#  endif /* CHOL_SOPALIN */
      task_gemm->handles[handle_idx++] = starpu_loop_data->WORK_handle;
#  ifndef STARPU_BLOCKTAB_SELFCOPY
      task_gemm->handles[handle_idx++] = starpu_loop_data->blocktab_handle;
#  endif /* STARPU_BLOCKTAB_SELFCOPY */

      ASSERT(handle_idx == cl_gemm.nbuffers, MOD_SOPALIN);

#  ifdef USE_TASK_DEP
      starpu_tasktab[iterbloc] = task_gemm;
      starpu_loop_data->deps[workerid][0] = starpu_loop_data->starpu_tasktab[SYMB_BLOKNUM(itercblk)];
      starpu_task_declare_deps_array(task_gemm, 1, starpu_loop_data->deps[workerid]);
#  endif /* USE_TASK_DEP */
      /* fprintf(stdout, "GEMM size %d %d %d %d  %d %d\n", */
      /*         task_gemm->cl->model->per_arch[STARPU_CPU_DEFAULT][0].size_base(task_gemm, */
      /*                                                                         STARPU_CPU, */
      /*                                                                         1), */
              /* SOLV_STRIDE(itercblk)-SOLV_COEFIND(iterbloc), */
              /* BLOK_ROWNBR(iterbloc), */
              /* CBLK_COLNBR(itercblk), */
              /* STARPU_MATRIX_GET_LD(task_gemm->handles[0]), */
              /* STARPU_MATRIX_GET_NY(task_gemm->handles[0])); */
#  ifdef STARPU_CONTEXT
      ret = starpu_task_submit_to_ctx(task_gemm,
                                          sched_ctxs[my_ctx]);
#  else
      ret = starpu_task_submit(task_gemm);
#  endif
      if (ret == -ENODEV) {
        fprintf(stderr, "No worker may execute this task (%d, %d)\n",
                task_gemm->execute_on_a_specific_worker, task_gemm->workerid);
      }
      STARPU_ASSERT(!ret);
    }

  return NO_ERR;
};
/*
 * Function: starpu_submit_loop
 *
 * Submit the tasks.
 */
#  define starpu_submit_loop API_CALL(starpu_submit_loop)
void*
starpu_submit_loop (void * arg)
{
  starpu_loop_data_t  *starpu_loop_data  = (starpu_loop_data_t*)(arg);
  Sopalin_Data_t      *sopalin_data      = (Sopalin_Data_t *)(starpu_loop_data->sopalin_data);
  SolverMatrix        *datacode          = sopalin_data->datacode;
  starpu_trf_data_t   *trf_args          = starpu_loop_data->trf_args;
  struct starpu_task **starpu_tasktab    = starpu_loop_data->starpu_tasktab;
  starpu_gemm_data_t  *gemm_args         = starpu_loop_data->gemm_args;
  int                  me                = starpu_loop_data->me;
  int                 *sched_ctxs        = starpu_loop_data->sched_ctxs;
  PASTIX_INT itertask;
  PASTIX_INT n_cblks = 0, n_tasks = 0;
  char      *prefetch;
  if ((prefetch = getenv("PASTIX_STARPU_PREFETCH_ON_NODE")) && !strcmp(prefetch, "1")) {
    /* Prefetch data on GPUs */
    PASTIX_INT   iterworker;
    PASTIX_INT * memory_nodes;
    fprintf(stdout, "Prefetching data on GPUs %s\n", prefetch);
    MALLOC_INTERN(memory_nodes, starpu_loop_data->ngpus, PASTIX_INT);
    for (iterworker = 0; iterworker < starpu_loop_data->ngpus; iterworker++) {
      memory_nodes[iterworker] = starpu_worker_get_memory_node(starpu_loop_data->gpu_workerids[iterworker]);
    }
    for (itertask=0;itertask<SOLV_TASKNBR;itertask++) {
      PASTIX_INT itercblk = TASK_CBLKNUM(itertask);
      if (starpu_loop_data->ngpus > 0 &&
          SOLV_COLOR(itercblk) >= 0) {
        PASTIX_INT workerid = starpu_loop_data->gpu_workerids[SOLV_COLOR(itercblk)];
        PASTIX_INT node = memory_nodes[workerid];
        starpu_data_prefetch_on_node(starpu_loop_data->L_handle[itercblk],
                                     node, 1);
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
        starpu_data_prefetch_on_node(starpu_loop_data->U_handle[itercblk],
                                     node, 1);
#    endif
#  endif
      }
    }
    memFree_null(memory_nodes);
  }
#  if (defined USE_TASK_DEP)
  {
    int i;
    int nworkers = 1;
#    ifdef STARPU_SUBMIT_READY
    nworkers = starpu_loop_data->last;
#    endif
    MALLOC_INTERN(starpu_loop_data->deps, nworkers, struct starpu_task **);
    MALLOC_INTERN(starpu_loop_data->deps_tab_size, nworkers, int);
    for (i = 0; i < nworkers; i++)
      {
        starpu_loop_data->deps_tab_size[i] = 1024;

        MALLOC_INTERN(starpu_loop_data->deps[i],
                      starpu_loop_data->deps_tab_size[i], struct starpu_task *);
      }
  }
#  endif

  /* For all column blocks we add a diag+trsm task.
   For all bloc in column block, we add a gemm task.
   */
  for (itertask=0;itertask<SOLV_TASKNBR;itertask++)
    {
      PASTIX_INT itercblk = TASK_CBLKNUM(itertask);
#  ifdef STARPU_SUBMIT_READY
      if (TASK_CTRBCNT((itertask))) continue;
#  endif
      starpu_submit_one_trf(itertask, sopalin_data);
#  ifndef STARPU_SUBMIT_READY
      starpu_submit_bunch_of_gemm(itertask, sopalin_data);
#  endif
    }

  /* Wait for all tasks to be finished */
/*   fprintf(stdout, "%ld cblks %ld tasks in context %d\n", */
/*           (long int)n_cblks, (long int)n_tasks, (int)me); */
  return NULL;
}

/*
 * Funciton starpu_clean_smp
 *
 * Clean thread data structures when using starpu.
 */
#  define starpu_clean_smp API_CALL(starpu_clean_smp)
void*
starpu_clean_smp (void * arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  sopalin_clean_smp ( sopalin_data, argument->me );
  return NULL;
}

/*
 Function: starpu_submit_tasks

 Submit tasks to perform the decomposition of the matrix.

 Parameters:
 sopalin_data - PaStiX global data structure.

 Returns:
 NO_ERR
 */
int
starpu_submit_tasks(Sopalin_Data_t  * sopalin_data) {
  SolverMatrix         * datacode         = sopalin_data->datacode;
  Thread_Data_t        * thread_data;
  starpu_trf_data_t    * trf_args;
  starpu_gemm_data_t   * gemm_args;
  starpu_data_handle_t * L_handle;
  starpu_data_handle_t * SM2X_handles = NULL;
#  ifdef USE_TASK_DEP
  PASTIX_INT                    task_number;
  struct starpu_task  ** starpu_tasktab;
#  endif
  starpu_data_handle_t * U_handle = NULL;
  starpu_data_handle_t   WORK_handle;
  PASTIX_INT itertask;
  int * blocktab;
#  ifdef STARPU_BLOCKTAB_SELFCOPY
  int ** d_blocktab;
#  else
  starpu_data_handle_t   blocktab_handle;
#  endif

  struct starpu_conf     conf;
  int                    ret;
  int                    cuda_nbr = sopalin_data->sopar->iparm[IPARM_CUDA_NBR];
  unsigned int * sched_ctxs = NULL;

#  ifdef STARPU_CONTEXT
  int * devices;
  int iter;
  PASTIX_INT thread_per_ctx;
#  endif

#ifndef STARPU_INIT_SMP
  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR+cuda_nbr, starpu_init_smp, sopalin_data,
                        0, NULL, NULL,
                        0, NULL, NULL);
#endif /* not STARPU_INIT_SMP ¨*/

#ifdef STARPU_PASTIX_SCHED
  {
    int k, bubnbr = datacode->bublnbr;
    int priomin = INT_MAX, priomax = INT_MIN;

    for(k = 0; k<bubnbr; k++)
      {
        priomin = MIN(priomin, datacode->btree->nodetab[k].priomin);
        priomax = MAX(priomax, datacode->btree->nodetab[k].priomax);
      }
    starpu_sched_set_min_priority(priomin);
    starpu_sched_set_max_priority(priomax);
  }
#endif

  starpu_conf_init(&conf);
  /* 1 GB */
#if (STARPU_MAJOR_VERSION > 1 || (STARPU_MAJOR_VERSION == 1 && STARPU_MINOR_VERSION >= 1))
  conf.trace_buffer_size = 1<<30;
#else
#error
#endif
  if (NULL != conf.sched_policy_name)
    {
      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        fprintf(stdout, OUT_STARPU_TP, conf.sched_policy_name);
    }
  else
    {
#ifdef STARPU_PASTIX_SCHED
      conf.sched_policy_name = NULL;
      conf.sched_policy = &starpu_pastix_sched_policy;
#else
      conf.sched_policy_name = "dmda";
      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        fprintf(stdout, OUT_STARPU_STP, conf.sched_policy_name);
#endif
    }
  conf.ncpus = SOLV_THRDNBR;
  conf.ncuda = cuda_nbr;
  conf.nopencl = 0;

  //starpu_profiling_set_id(starpu_id++);
  if (0 != (ret = starpu_init(&conf)))
    {
      errorPrint("Error %d initializing StarPU\n", ret);
    }

#  ifdef STARPU_CONTEXT
  MALLOC_INTERN(devices, SOLV_THRDNBR+cuda_nbr, int);
  starpu_worker_get_ids_by_type(STARPU_CPU_WORKER, devices, SOLV_THRDNBR);

#    ifdef STARPU_USE_CUDA
  starpu_worker_get_ids_by_type(STARPU_CUDA_WORKER, devices+SOLV_THRDNBR, cuda_nbr);
#    endif
  /*create contexts however you want*/
  fprintf(stdout, "creating %d contexts \n", (int)sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR]);
  thread_per_ctx = SOLV_THRDNBR/(sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR]-1);
  if (SOLV_THRDNBR%(sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR]-1))
    thread_per_ctx++;

  MALLOC_INTERN(sched_ctxs, sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR], unsigned);
  sched_ctxs[0] = starpu_sched_ctx_create("dmda", devices+SOLV_THRDNBR, cuda_nbr, "ctx_0");
  for (iter = 1; iter < sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR]; iter++)
    {
      char string[128];
      int nthreads = thread_per_ctx;
      if (iter == sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR]-1 &&
          SOLV_THRDNBR%(thread_per_ctx) != 0)
        nthreads = SOLV_THRDNBR%(thread_per_ctx);
      sprintf(string, "ctx_%d", iter);
      fprintf(stdout, "creating %s contexts with %d cores %d\n", string, nthreads, thread_per_ctx);

      sched_ctxs[iter] = starpu_sched_ctx_create("dmda",
                                                 devices+(iter-1)*thread_per_ctx,
                                                 nthreads, string);
      starpu_sched_ctx_set_inheritor(sched_ctxs[iter], sched_ctxs[0]);
    }
#endif

#  ifdef STARPU_PROFILING
  if ((ret = starpu_profiling_status_set(STARPU_PROFILING_ENABLE) < 0))
    {
      errorPrint("Error %d in starpu_profiling_enable\n", ret);
    }
#  endif


#  ifdef STARPU_INIT_SMP
  {
    int threadid;
    sopthread_data_t * init_arg;
    MALLOC_INTERN(init_arg, SOLV_THRDNBR+cuda_nbr, sopthread_data_t);
    for (threadid = 0; threadid < SOLV_THRDNBR+cuda_nbr; threadid++)
      {
        struct starpu_task * task_init;
#    ifdef STARPU_CONTEXT
        PASTIX_INT my_ctx;
        my_ctx = 1+threadid/thread_per_ctx;
#    endif
        task_init = starpu_task_create();
        /* We compute GEMM */
        init_arg[threadid].me   = threadid;
        init_arg[threadid].data = sopalin_data;
        task_init->cl = &sopalin_init_cl;
        task_init->cl_arg = &(init_arg[threadid]);
#    ifdef STARPU_PASTIX_SCHED
        task_init->workerid = threadid;
        task_init->priority = 0; /* No priority needed as the task needs
                                  * to be runned first due to dependancies */
#    endif

#    ifdef STARPU_CONTEXT
        ret = starpu_task_submit_to_ctx(task_init,
                                        sched_ctxs[my_ctx]);
#    else
        ret = starpu_task_submit(task_init);
#    endif
        STARPU_ASSERT(!ret);
      }
    /* wait for end of init */
    starpu_task_wait_for_all();
    memFree_null(init_arg);
  }
#  endif

  {
    /* build blocktab */
    int iterblock;
    MALLOC_INTERN(blocktab, SYMB_BLOKNBR*2, int);
    if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
      fprintf(stdout, "sizeof blocktab : %d integers\n",
              (int)(2*SYMB_BLOKNBR));

    for (iterblock = 0; iterblock < SYMB_BLOKNBR; iterblock++)
      {
        blocktab[2*iterblock]   = SYMB_FROWNUM(iterblock);
        blocktab[2*iterblock+1] = SYMB_LROWNUM(iterblock);
      }
#  ifdef STARPU_BLOCKTAB_SELFCOPY
    {
      int ndevices;
      int device_id;
      CUDA_CALL(cudaGetDeviceCount(&ndevices));
      MALLOC_INTERN(d_blocktab, ndevices, int*);
      for (device_id = 0; device_id < ndevices; device_id++)
        {
          cudaSetDevice(device_id);
          CUDA_CALL(cudaMalloc((void*)&(d_blocktab[device_id]),
                               2*SYMB_BLOKNBR*sizeof(int)));
          CUDA_CALL(cudaMemcpy((void*)d_blocktab[device_id], blocktab,
                               2*SYMB_BLOKNBR*sizeof(int),
                               cudaMemcpyHostToDevice));
        }
    }
#  else
    starpu_vector_data_register(&blocktab_handle, 0,
                                (uintptr_t)blocktab, 2*SYMB_BLOKNBR,
                                sizeof(int));
    starpu_data_set_sequential_consistency_flag(blocktab_handle, 0);
#  endif
  }

#  ifdef PASTIX_DUMP_FACTO
  dump_all(datacode, sopalin_data->sopar->cscmtx,
           ((datacode->updovct.sm2xtab!=NULL)?
            (DUMP_CSC | DUMP_SOLV | DUMP_SMB):(DUMP_CSC | DUMP_SOLV)));
#  endif

  thread_data = sopalin_data->thread_data[0];
  sopalin_data->sopar->diagchange = 0;
  SOPALIN_CLOCK_INIT;
#  ifdef STARPU_USE_CUDA
  /* starpu_helper_cublas_init(); */
#  endif


#  ifdef CHOL_SOPALIN
  starpu_vector_data_register(&WORK_handle, -1, (uintptr_t)NULL, SOLV_COEFMAX,
                              sizeof(PASTIX_FLOAT));
#  else
  starpu_vector_data_register(&WORK_handle, -1, (uintptr_t)NULL, 2*SOLV_COEFMAX,
                              sizeof(PASTIX_FLOAT));
#  endif

  MALLOC_INTERN(trf_args,    SYMB_CBLKNBR, starpu_trf_data_t);

  MALLOC_INTERN(L_handle,    SYMB_CBLKNBR, starpu_data_handle_t);

#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  MALLOC_INTERN(U_handle,    SYMB_CBLKNBR, starpu_data_handle_t);
#    endif
#  endif

  MALLOC_INTERN(gemm_args,      SYMB_BLOKNBR, starpu_gemm_data_t);


  {
    int itercblk;
    int max_cblksize = 0;
    int max_cblkcolnbr = 0;
    for (itercblk = 0; itercblk < SYMB_CBLKNBR; itercblk++)
      {
        max_cblksize   = MAX(max_cblksize,
                             CBLK_COLNBR(itercblk)*SOLV_STRIDE(itercblk));
        max_cblkcolnbr = MAX(max_cblkcolnbr,
                             CBLK_COLNBR(itercblk));
      }

    if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
      fprintf(stdout, "Maximum cblk size %d, maximu cblk colnbr %d\n",
              max_cblksize, max_cblkcolnbr);
  }

#  ifdef USE_TASK_DEP
  task_number = 0;
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    task_number = SYMB_BLOKNBR;
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    task_number += 2*SYMB_BLOKNBR+SYMB_CBLKNBR;
  MALLOC_INTERN(starpu_tasktab, task_number, struct starpu_task *);
  {
    PASTIX_INT iter;
    for (iter = 0; iter < task_number; iter++)
      starpu_tasktab[iter] = NULL;
  }
#  endif
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
  starpu_partition_data(sopalin_data, L_handle, U_handle);
#  else
  starpu_partition_data(sopalin_data, L_handle, NULL);
#  endif
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    {
      if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_SOLVE)
        if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
          errorPrintW("Raffinement not available with StarPU,"
                      " only performing solve\n");
      MALLOC_INTERN(SM2X_handles, SYMB_CBLKNBR, starpu_data_handle_t);
      starpu_register_sm2x(sopalin_data, SM2X_handles);
    }
  SOPALIN_CLOCK_STOP;
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    fprintf(stdout,"----- Time after data registration %lf s\n",
            SOPALIN_CLOCK_GET);

  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    {
      starpu_loop_data_t * starpu_loop_data;
      pthread_cond_t  cond_end_facto;
      pthread_mutex_t mutex_end_facto;

      pthread_cond_init(&cond_end_facto, NULL);
      pthread_mutex_init(&mutex_end_facto, NULL);

      MALLOC_INTERN(starpu_loop_data, 1, starpu_loop_data_t);
      sopalin_data->starpu_loop_data = starpu_loop_data;
      starpu_loop_data->me               = 0;
      starpu_loop_data->trf_args         = trf_args;
      starpu_loop_data->gemm_args        = gemm_args;
      starpu_loop_data->starpu_tasktab   = starpu_tasktab;
      starpu_loop_data->L_handle         = L_handle;
#    if (defined CHOL_SOPALIN && defined SOPALIN_LU)
      starpu_loop_data->U_handle         = U_handle;
#    endif
      starpu_loop_data->WORK_handle      = WORK_handle;
#    ifdef STARPU_BLOCKTAB_SELFCOPY
      starpu_loop_data->d_blocktab       = d_blocktab;
#    else
      starpu_loop_data->blocktab_handle  = blocktab_handle;
#    endif
      starpu_loop_data->sopalin_data     = sopalin_data;
      starpu_loop_data->ctx_nbr          = 1;
      starpu_loop_data->first            = 0;
      starpu_loop_data->last             = SOLV_THRDNBR;
#    ifdef STARPU_CONTEXT
      starpu_loop_data->thread_per_ctx   = thread_per_ctx;
      starpu_loop_data->sched_ctxs       = sched_ctxs;
#    endif
      starpu_loop_data->ndiag_submitted = 0;
      starpu_loop_data->ngemm_submitted = 0;
      starpu_loop_data->cond_end_facto  = &cond_end_facto;
      starpu_loop_data->mutex_end_facto = &mutex_end_facto;
      MALLOC_INTERN(starpu_loop_data->cpu_workerids, SOLV_THRDNBR, int);
      starpu_loop_data->ncpus = SOLV_THRDNBR;
      MALLOC_INTERN(starpu_loop_data->gpu_workerids, cuda_nbr, int);
      MALLOC_INTERN(starpu_loop_data->gpu_gemm_count, cuda_nbr, int);
      memset(starpu_loop_data->gpu_gemm_count, 0, cuda_nbr*sizeof(int));
      starpu_loop_data->ngpus = cuda_nbr;
      starpu_worker_get_ids_by_type(STARPU_CPU_WORKER,
				    starpu_loop_data->cpu_workerids,
				    SOLV_THRDNBR);
      starpu_worker_get_ids_by_type(STARPU_CUDA_WORKER,
				    starpu_loop_data->gpu_workerids,
				    cuda_nbr);
      fprintf(stdout, "cuda_nbr %d\n", cuda_nbr);
      {
        int j;
        for (j = 0; j < cuda_nbr; j++)
          fprintf(stdout, "cuda_id %d\n", starpu_loop_data->gpu_workerids[j]);
      }
      starpu_submit_loop (starpu_loop_data);
    }
  else
    sopalin_data->starpu_loop_data = NULL;

  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    starpu_submit_updown(sopalin_data, L_handle, U_handle, SM2X_handles,
                         starpu_tasktab, (int*)sched_ctxs);
#    ifdef STARPU_SUBMIT_READY
  if (sopalin_data->starpu_loop_data != NULL)
    {
      while(starpu_tasktab[SYMB_BLOKNBR-1] == NULL)
	{
	  COND_WAIT(sopalin_data->starpu_loop_data->cond_end_facto,
		    sopalin_data->starpu_loop_data->mutex_end_facto);
	}
      fprintf(stdout, "SUBMITED GEMM %d DIAG %d\n",
	      sopalin_data->starpu_loop_data->ngemm_submitted,
	      sopalin_data->starpu_loop_data->ndiag_submitted);

      pthread_mutex_unlock(sopalin_data->starpu_loop_data->mutex_end_facto);

      pthread_cond_destroy(sopalin_data->starpu_loop_data->cond_end_facto);
      pthread_mutex_destroy(sopalin_data->starpu_loop_data->mutex_end_facto);
    }
#    endif
  starpu_task_wait_for_all();
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    {
      int i;
      int nworkers = 1;
#  if ( defined STARPU_SUBMIT_READY)
      nworkers = SOLV_THRDNBR;
#  endif
      for (i =0; i < nworkers; i++)
	if (sopalin_data->starpu_loop_data->deps_tab_size[i] > 0)
	  memFree_null(sopalin_data->starpu_loop_data->deps[i]);
      memFree_null(sopalin_data->starpu_loop_data->deps_tab_size);
      memFree_null(sopalin_data->starpu_loop_data->deps);
      memFree_null(sopalin_data->starpu_loop_data->cpu_workerids);
      memFree_null(sopalin_data->starpu_loop_data->gpu_workerids);
      for (i = 0; i < sopalin_data->sopar->iparm[IPARM_CUDA_NBR]; i++)
        fprintf(stdout, "%d GEMMs forced on GPU %d\n",
                sopalin_data->starpu_loop_data->gpu_gemm_count[i], i);
      memFree_null(sopalin_data->starpu_loop_data->gpu_gemm_count);
      memFree_null(sopalin_data->starpu_loop_data);
    }
  SOPALIN_CLOCK_STOP;
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    fprintf(stdout,"----- submission and wait for all %lf s (%d tasks)\n",
            SOPALIN_CLOCK_GET, SYMB_BLOKNBR);


  /* Unregister buffers and leave starpu */
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    {
#  if(defined STARPU_PROFILING && defined USE_TASK_DEP)
      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        {
          int worker;
          /* Display the occupancy of all workers during the test */

          for (worker = 0; worker < starpu_worker_get_count(); worker++)
            {
              struct starpu_profiling_worker_info worker_info;
              int ret = starpu_profiling_worker_get_info(worker, &worker_info);
              STARPU_ASSERT(!ret);

              double total_time     = starpu_timing_timespec_to_us(&worker_info.total_time);
              double executing_time = starpu_timing_timespec_to_us(&worker_info.executing_time);
              double sleeping_time  = starpu_timing_timespec_to_us(&worker_info.sleeping_time);

              float executing_ratio = 100.0*executing_time/total_time;
              float sleeping_ratio  = 100.0*sleeping_time/total_time;

              char workername[128];

              double delay_sum[2]  = {0.0, 0.0};
              double length_sum[2] = {0.0, 0.0};
              unsigned int cnt[2]  = {0,   0};

              for (itertask=0;itertask<SOLV_TASKNBR;itertask++)
                {
                  PASTIX_INT itercblk = TASK_CBLKNUM(itertask);
                  PASTIX_INT iterbloc;
                  struct starpu_profiling_task_info *info;
                  info = starpu_tasktab[SYMB_BLOKNUM(itercblk)]->profiling_info;
                  if (info->workerid == worker)
                    {
                      /* How much time did it take before the task started ? */
                      delay_sum[0] += starpu_timing_timespec_delay_us(&info->submit_time,
                                                                      &info->start_time);

                      /* How long was the task execution ? */
                      length_sum[0] += starpu_timing_timespec_delay_us(&info->start_time,
                                                                       &info->end_time);
                      cnt[0]++;
                    }

                  for (iterbloc = SYMB_BLOKNUM(itercblk)+1;
                       iterbloc < SYMB_BLOKNUM(itercblk+1);
                       iterbloc ++)
                    {
                      info = starpu_tasktab[iterbloc]->profiling_info;
                      if (info->workerid == worker)
                        {
                          /* How much time did it take before the task started ? */
                          delay_sum[1] += starpu_timing_timespec_delay_us(&info->submit_time,
                                                                          &info->start_time);

                          /* How long was the task execution ? */
                          length_sum[1] += starpu_timing_timespec_delay_us(&info->start_time,
                                                                           &info->end_time);
                          cnt[1]++;

                        }
                    }
                }
              starpu_worker_get_name(worker, workername, 128);
              fprintf(stdout, "Worker %s:\n", workername);
              if (cnt[0] != 0)
                {
                  fprintf(stdout, "Avg. delay on XXTRF : %2.2lf us, %d tasks\n",
                          (delay_sum[0])/cnt[0], cnt[0]);
                  fprintf(stdout, "Avg. length on XXTRF : %2.2lf us\n",
                          (length_sum[0])/cnt[0]);
                }
              if (cnt[1] != 0)
                {
                  fprintf(stdout, "Avg. delay on XXMM : %2.2lf us, %d tasks\n",
                          (delay_sum[1])/cnt[1], cnt[1]);
                  fprintf(stdout, "Avg. length on XXMM : %2.2lf us\n",
                          (length_sum[1])/cnt[1]);
                }


              fprintf(stdout, "\ttotal time : %.2lf ms\n", total_time*1e-3);
              fprintf(stdout, "\texec time  : %.2lf ms (%.2f %%)\n", executing_time*1e-3, executing_ratio);
              fprintf(stdout, "\tblocked time  : %.2lf ms (%.2f %%)\n", sleeping_time*1e-3, sleeping_ratio);
            }
        }
#  endif /* USE_TASK_DEP && STARPU_PROFILING*/
    }
  for (itertask=0;itertask<SOLV_TASKNBR;itertask++)
    {
      PASTIX_INT itercblk = TASK_CBLKNUM(itertask);
      PASTIX_INT iterbloc;
#  ifdef USE_TASK_DEP
      for (iterbloc = SYMB_BLOKNUM(itercblk);
           iterbloc < SYMB_BLOKNUM(itercblk+1);
           iterbloc ++)
        {
          PASTIX_INT first_task = 0;
          if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
            {
              starpu_task_destroy(starpu_tasktab[iterbloc]);
              first_task = SYMB_BLOKNBR;
            }
          if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
            {
              starpu_task_destroy(starpu_tasktab[first_task+iterbloc]);
              first_task += SYMB_BLOKNBR;
#ifndef CHOL_SOPALIN
              first_task += SYMB_CBLKNBR;
#endif /* not CHOL_SOPALIN */
              starpu_task_destroy(starpu_tasktab[first_task+iterbloc]);
            }
        }
#    ifndef CHOL_SOPALIN
      {
        if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
          {
            PASTIX_INT first_task = SYMB_BLOKNBR;
            if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
              first_task += SYMB_BLOKNBR;
            starpu_task_destroy(starpu_tasktab[first_task+itercblk]);
          }
      }
#    endif /* CHOL_SOPALIN */
#  endif
      starpu_data_unregister(L_handle[itercblk]);
      if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
        {
          starpu_data_unregister(SM2X_handles[itercblk]);
        }
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
      starpu_data_unregister(U_handle[itercblk]);
#    endif
#  endif
    }
  starpu_data_unregister(WORK_handle);
#  ifdef USE_TASK_DEP
  memFree_null(starpu_tasktab);
#  endif

#  ifdef STARPU_BLOCKTAB_SELFCOPY
  {
    int ndevices;
    int device_id;
    CUDA_CALL(cudaGetDeviceCount(&ndevices));
    for (device_id = 0; device_id < ndevices; device_id++)
      {
        cudaSetDevice(device_id);
        CUDA_CALL(cudaFree(d_blocktab[device_id]));
      }
    memFree_null(d_blocktab);
  }
#  else
  starpu_data_unregister(blocktab_handle);
#  endif /* STARPU_BLOCKTAB_SELFCOPY */
  memFree_null(blocktab);

  /* Reduction on pivot number */
  sopalin_data->sopar->diagchange = 0;
  {
    PASTIX_INT me;
    for (me = 0; me < SOLV_THRDNBR+cuda_nbr; me++)
      {
        sopalin_data->sopar->diagchange += sopalin_data->thread_data[me]->nbpivot;
      }
  }
#  ifdef STARPU_USE_CUDA
  /*starpu_helper_cublas_shutdown();*/
#  endif

  SOPALIN_CLOCK_STOP;
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    fprintf(stdout,"----- sopalin time %lf\n",
            SOPALIN_CLOCK_GET);
  sopalin_data->sopar->dparm[DPARM_FACT_TIME] = SOPALIN_CLOCK_GET;
#  ifdef PASTIX_DUMP_FACTO
  dump_all(datacode, sopalin_data->sopar->cscmtx,
           ((datacode->updovct.sm2xtab!=NULL)?
            (DUMP_CSC | DUMP_SOLV | DUMP_SMB):(DUMP_CSC | DUMP_SOLV)));
#  endif
  memFree_null(trf_args);
  memFree_null(gemm_args);
  memFree_null(L_handle);
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    memFree_null(SM2X_handles);
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  memFree_null(U_handle);
#    endif /* CHOL_SOPALIN */
#  endif /* SOPALIN_LU   */

  starpu_shutdown();
  sopalin_clean(sopalin_data, 1);

  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR, starpu_clean_smp, sopalin_data,
                        0, NULL, NULL,
                        0, NULL, NULL);
  return NO_ERR;
}
#else
/* ISO C forbids an empty source file */
#  include "not_empty.h"
NOT_EMPTY(starpu_submit_tasks)
#endif
