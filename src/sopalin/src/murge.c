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
 * File: murge.c
 *
 * This file implements <Murge> interface.
 *
 * About: Authors
 *   Mathieu Faverge - faverge@labri.fr
 *   Xavier Lacoste  - xavier.lacoste@inria.fr
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#ifndef FORCE_NOSMP
#  include <pthread.h>
#else
#  ifdef MURGE_THREADSAFE
#    undef MURGE_THREADSAFE
#  endif
#endif
#ifdef FORCE_NOMPI
#  include "nompi.h"
#else /* not FORCE_NOMPI */
#  include <mpi.h>
#endif /* not FORCE_NOMPI */
#ifdef WITH_SEM_BARRIER
#  include <semaphore.h>
#endif
#include "common_pastix.h"
#include "tools.h"
#include "sopalin_define.h"

#ifdef WITH_SCOTCH
#  ifdef    DISTRIBUTED
#    include "ptscotch.h"
#  else
#    include "scotch.h"
#  endif /* DISTRIBUTED */
#endif /* WITH_SCOTCH */

#include "ftgt.h"
#include "symbol.h"
#include "csc.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
#include "order.h"
#include "sopalin_thread.h"
#include "stack.h"
#include "sopalin3d.h"
#include "pastixstr.h"
#include "pastix.h"
#include "cscd_utils.h"
#include "cscd_utils_intern.h"
#include "csc_intern_compute.h"
#include "csc_intern_updown.h"
#include "sopalin_init.h"
#include "sopalin_compute.h"

#ifdef ZMURGE
#  include "zmurge.h"
#  include "zmurge_pastix.h"
#endif
#ifdef CMURGE
#  include "cmurge.h"
#  include "cmurge_pastix.h"
#endif
#ifdef DMURGE
#  include "dmurge.h"
#  include "dmurge_pastix.h"
#endif
#ifdef SMURGE
#  include "smurge.h"
#  include "smurge_pastix.h"
#endif
#if !(defined ZMURGE || defined DMURGE || defined CMURGE || defined SMURGE)
#  include "murge.h"
#  include "murge_pastix.h"
#endif
#include "murge_defines.h"

// ALWAYS_INLINE ---------------------------------------------//
// Macro to use in place of 'inline' to force a function to be inline
#define ALWAYS_INLINE
#if !defined(ALWAYS_INLINE)
#  if defined(_MSC_VER)
#    define ALWAYS_INLINE __forceinline
#  elif defined(__GNUC__) && __GNUC__ > 3
// Clang also defines __GNUC__ (as 4)
#    define ALWAYS_INLINE __attribute__ ((__always_inline__))
#  else
#    define ALWAYS_INLINE
#  endif
#endif
/* #define ALWAYS_INLINE __attribute__((always_inline)) */

/* internal function declarations */
static inline
INTS MURGE_GraphSetEdge_ (INTS id, INTS ROW, INTS COL) ALWAYS_INLINE;

static inline
INTS MURGE_AssemblySetValue_    (INTS id, INTS ROW, INTS COL, COEF value) ALWAYS_INLINE;
static inline
INTS MURGE_AssemblySetNodeValues_ (INTS id, INTS ROW, INTS COL, COEF *values) ALWAYS_INLINE;
static inline
INTS MURGE_AssemblySetBlockValues_(INTS id, INTS nROW, INTS *ROWlist,
                                   INTS nCOL, INTS *COLlist, COEF *values) ALWAYS_INLINE;

static inline
INTS MURGE_GraphSetEdge_ (INTS id, INTS ROW, INTS COL) ALWAYS_INLINE;


/******************************************************************************/
/***                           Section: Structures                          ***/
/******************************************************************************/

/*
 * Structure: ijv_
 *
 * Structure to represente coefficients.
 *
 * Contains:
 *   i     - row
 *   j     - column
 *   v     - pointer to the value array (can be several degree of freedom)
 *   owner - process which own the coefficient.
 */
struct ijv_ {
  PASTIX_INT    i;
  PASTIX_INT    j;
  PASTIX_FLOAT* v;
#ifdef MURGE_FOLLOW_INPUT_ORDER
  PASTIX_INT    idx;
#endif
  int    owner;
};

/*
 * Typedef: ijv_t
 *
 * Alias to structure <ijv_>.
 */
typedef struct ijv_ ijv_t;

/*
 * struct: murge_seq_t
 *
 * Structure used to store assembly sequence.
 *
 * Contains:
 *   indexes              - Order sequences of indexes to be set.
 *   coefnbr              - Number of entries in the sequence.
 *   recv_nbr             - Number of entries to receive from each processor.
 *   recv_indexes         - Indexes of entries which will be received.
 *   mode                 - Local entries or communicating mode
 *   fusion_local_entries - Operation to perform when a coefficient appear twice
 *   fusion_dist_entries  - Operation to perform when a coefficient appear
 *                            twice, given by two processors.
 *   ijv_size             - size of the required array to store
 *                            not local entries.
 *   nodes                - 0 entries are entered value by value,
 *                          1 entries are entries node by node.
 *   next                 - Next entry in the list of sequences.
 */
typedef struct murge_seq_ murge_seq_t;
struct murge_seq_ {
  INTL        * indexes;
  INTL          coefnbr;
  INTL        * recv_nbr;
  INTL       ** recv_indexes;
  INTS          mode;
  PASTIX_FLOAT       (*fusion_local_entries)(PASTIX_FLOAT , PASTIX_FLOAT);
  PASTIX_FLOAT       (*fusion_dist_entries)(PASTIX_FLOAT , PASTIX_FLOAT);
  INTL          ijv_size;
  INTS          nodes;
  murge_seq_t * next;
  PASTIX_INT           ID;
};

typedef struct variable_csc_s * variable_csc_t;


/*
 * Typedef: murge_data_t
 *
 * alias to structure <murge_data_>.
 */
typedef struct murge_data_ murge_data_t;

/*
 * struct: murge_product_data_
 *
 * thread_id - thread identification number
 * solver    - ponter to <murge_data_t> structure.
 * t_prod    - local part of the product computation.
 * all_prods - pointer to local products of each threads.
 * ret       - return value.
 */
struct murge_product_data_
{
  INTS           thread_id;
  murge_data_t * solver;
  COEF         * t_prod;
  COEF        ** all_prods;
  INTS           ret;
  INTS           counter;
};

/*
 * typedef: murge_product_data_t
 *
 * alias to struct <murge_product_data_>
 */
typedef struct murge_product_data_ murge_product_data_t;

struct murge_thread_data_ {
  murge_product_data_t * pdata;
  murge_data_t         * solver;
};

/*
 * typedef: murge_thread_data_t
 *
 * alias to struct <murge_thread_data_>
 */
typedef struct murge_thread_data_ murge_thread_data_t;

typedef enum murge_thread_state_ {
  MURGE_THREAD_WAIT,
  MURGE_THREAD_PRODUCT,
  MURGE_THREAD_END
} murge_thread_state_t;


/*
 * struct: murge_data_t
 *
 * Structure used to store murge data
 *
 * Contains:
 *   pastix_data    - Pointer to the <pastix_data_t> associated to
 *                    the solver instance
 *   n              - Number of local column indicated by murge user
 *   N              - Number of global column indicated by murge user
 *   colptr         - Colptr in murge's user CSCd
 *   rows           - Rows in murge's user CSCd
 *   values         - Values in murge's user CSCd
 *   l2g            - Local to global column number in murge's user CSCd
 *   g2l            - Global to local column number in murge's user CSCd
 *   perm           - Permtab for murge's user
 *   b              - Right-hand-side member(s) given by murge's user
 *   nrhs           - Number of right-hand-side member(s) given by murge's user
 *   cnt            - Iterator for number of entered edges
 *   edgenbr        - Number of edges
 *   state          - State of the solver
 *   mode           - Local entries or communicating mode
 *   op             - Operation to perform when a coefficient appear twice
 *   op2            - Operation to perform when a coefficient appear twice,
 *                    given by two processors.
 *   sym            - Indicate if we have to check that the matrix is symmetric
 *   sequences      - List of saved sequences.
 *   malloc_size    - Current allocated memory.
 *   malloc_maxsize - Maximum allocated memory.
 *   threadnbr      - Number of thread launched.
 *   threads        - Array of pthread_t.
 *   threads_data   - data associated to each thread.
 *   thread_state   - flag to control the threads.
 *   barrier        - barrier used for the threads.
 */
struct murge_data_ {
  pastix_data_t          *pastix_data;
  PASTIX_INT              n;
  PASTIX_INT              N;
  PASTIX_INT             *colptr;
  PASTIX_INT             *rows;
  PASTIX_FLOAT           *values;
  PASTIX_INT             *l2g;
  PASTIX_INT             *g2l;
  PASTIX_INT             *perm;
#ifdef CENTRALISED
  PASTIX_INT             *invp;
#endif
  PASTIX_FLOAT           *b;
  PASTIX_INT              nrhs;
  variable_csc_t         vcsc;
#ifdef MURGE_THREADSAFE
  pthread_mutex_t         mutex_tmpmatrix;
#endif
  PASTIX_INT              cnt;
  PASTIX_INT              cnt_zero;
  PASTIX_INT              cnt_node;
  PASTIX_INT              edgenbr;
  PASTIX_INT              coefnbr;
  PASTIX_INT              nodenbr;
  int                     state;
  int                     mode;
  int                     op;
  int                     op2;
  int                     sym;
  int                     dynamic;
  murge_seq_t            *sequences;
  PASTIX_INT              seq_ID;
  int                     ndump;
  char                   *dropmask;
  char                   *droprows;
  char                   *dropcols;
  int64_t                 malloc_size;
  int64_t                 malloc_maxsize;
  int                     threadnbr;
  pthread_t              *threads;
  murge_thread_data_t    *threads_data;
  murge_thread_state_t    threads_state;
  sopthread_barrier_t     barrier;
  pthread_mutex_t         mutex_state;
  pthread_cond_t          cond_state;
};

/******************************************************************************/
/***                           Section: Global variables                    ***/
/******************************************************************************/


/*
 * Variables: Global variables
 *
 * idnbr   - Number of solvers instances.
 * solvers - Murge solver instances array (<murge_data_t>).
 */
INTS           idnbr   = 0;
murge_data_t **solvers = NULL;
#include "variable_csc.c"


/******************************************************************************/
/***                           Section: Functions                           ***/
/******************************************************************************/


/*******************************************************************************
 * Group: Auxilary functions
 */

#ifdef CENTRALISED
#  define ALLOC_INVP                                                      \
  if (NULL == solvers[id]->invp)                                        \
    MURGE_MEMALLOC(solvers[id]->invp, solvers[id]->n, PASTIX_INT);
#  define INVP solvers[id]->invp
#else
#  define ALLOC_INVP {};
#  define INVP NULL
#endif


#define product_thread PASTIX_PREFIX_F(product_thread)
void* product_thread(void * data) {
  INTS first, last;
  INTS gfirst, glast;
  INTS nnz, nnz_per_thread;
  INTS itercol;
  INTS iterrows;
  INTS dof, baseval;
  INTS row;
  COEF * mat = NULL;
  murge_product_data_t * pdata = (murge_product_data_t *)data;
  INTS id = &(pdata->solver) - solvers;
  INTS threadnbr = pdata->solver->pastix_data->iparm[IPARM_THREAD_NBR];

  baseval = pdata->solver->pastix_data->iparm[IPARM_BASEVAL];
  dof = pdata->solver->pastix_data->iparm[IPARM_DOF_NBR];
  pdata->ret = MURGE_SUCCESS;

  /* Initialize buffer */
  if (NULL == pdata->t_prod) {
    MURGE_MEMALLOC_RET(pdata->t_prod, pdata->solver->N*dof, COEF, pdata->ret);
    pdata->all_prods[pdata->thread_id] = pdata->t_prod;
  }
  memset(pdata->t_prod, 0, pdata->solver->N*dof*sizeof(COEF));

  if (pdata->solver->n > 0)
    nnz = pdata->solver->colptr[pdata->solver->n]-baseval;
  else
    nnz = 0;
  nnz_per_thread = nnz/threadnbr;
  gfirst = pdata->thread_id * nnz_per_thread;
  glast = (pdata->thread_id +1)* nnz_per_thread;

  if (pdata->thread_id == threadnbr-1) {
    glast = nnz;
  }
  for (itercol = 0; itercol < pdata->solver->n; itercol++) {
    if (pdata->solver->colptr[itercol+1]-baseval > gfirst ||
        pdata->solver->colptr[itercol]-baseval  < glast) {
      first = MAX(gfirst, pdata->solver->colptr[itercol]-baseval);
      last  = MIN(glast,  pdata->solver->colptr[itercol+1]-baseval);
      for (iterrows = first;
           iterrows < last;
           iterrows++) {
        row = pdata->solver->rows[iterrows]-baseval;
        mat = &(pdata->solver->values[iterrows*dof*dof]);
        SOPALIN_GEMV("N", dof, dof, 1.0, mat, dof,
                     &(pdata->solver->b[itercol*dof]), 1, 1.0,
                     &(pdata->t_prod[row*dof]), 1);
      }

    }
  }

  /* reduction */
  SYNCHRO_X_THREAD(threadnbr, pdata->solver->barrier);
  {
    int i, j;
    for (i = threadnbr/2; i > 0; i/=2) {
      if (pdata->thread_id < i) {
        /* fuse t_prod from thread_id with t_prod from thread_id+i */
        for (j = 0; j < pdata->solver->N*dof; j++) {
          pdata->all_prods[pdata->thread_id][j] += pdata->all_prods[pdata->thread_id+i][j];
        }
      }
      SYNCHRO_X_THREAD(threadnbr, pdata->solver->barrier);
    }

    if (threadnbr%2 && pdata->thread_id == 0 && threadnbr>1) {
      /* fuse t_prod from thread 0 and t_prod from threadnbr-1 */
      for (j = 0; j < pdata->solver->N*dof; j++) {
        pdata->all_prods[0][j] += pdata->all_prods[threadnbr-1][j];
      }
    }
  }

  if (pdata->thread_id == 0) {
    pdata->solver->threads_state = MURGE_THREAD_WAIT;
    pthread_cond_broadcast(&(pdata->solver->cond_state));
  }
  SYNCHRO_X_THREAD(threadnbr, pdata->solver->barrier);
  return 0;
}

#define control_thread PASTIX_PREFIX_F(control_thread)
void* control_thread(void * data) {
  murge_thread_data_t * tdata= (murge_thread_data_t*)data;
#ifndef FORCE_NOSMP
  pthread_mutex_lock(&(tdata->solver->mutex_state));
#endif
  while (tdata->solver->threads_state != MURGE_THREAD_END) {
#ifndef FORCE_NOSMP
    pthread_cond_wait(&(tdata->solver->cond_state),
                      &(tdata->solver->mutex_state));
#endif
    switch(tdata->solver->threads_state) {
    case MURGE_THREAD_END:
      pthread_mutex_unlock(&(tdata->solver->mutex_state));
      return NULL;
      break;
    case MURGE_THREAD_WAIT:
      break;
    case MURGE_THREAD_PRODUCT:
      pthread_mutex_unlock(&(tdata->solver->mutex_state));
      product_thread(tdata->pdata);
      break;
    default:
      errorPrint("Undefined state : %d", tdata->solver->threads_state);
      break;
    }
  }
  return NULL;
}

#define start_threads PASTIX_PREFIX_F(start_threads)
INTS start_threads(INTS id) {
  INTS iter;
#ifndef FORCE_NOSMP
  int ret;
#endif
  COEF ** all_prods;
  solvers[id]->threadnbr = solvers[id]->pastix_data->iparm[IPARM_THREAD_NBR];
  MURGE_MEMALLOC(solvers[id]->threads_data,
                 solvers[id]->threadnbr,
                 murge_thread_data_t);
  MURGE_MEMALLOC(solvers[id]->threads,
                 solvers[id]->threadnbr,
                 pthread_t);
  MURGE_MEMALLOC(all_prods, solvers[id]->threadnbr, COEF*);
  solvers[id]->threads_state = MURGE_THREAD_WAIT;
  for (iter = 0; iter < solvers[id]->threadnbr; iter++) {
    MURGE_MEMALLOC(solvers[id]->threads_data[iter].pdata,
                   solvers[id]->threadnbr,
                   murge_product_data_t);
    solvers[id]->threads_data[iter].solver = solvers[id];
    {
      /* Initialize pdata */
      murge_product_data_t * pdata = solvers[id]->threads_data[iter].pdata;
      pdata->thread_id = iter;
      pdata->solver = solvers[id];
      pdata->t_prod = NULL;
      pdata->all_prods = all_prods;
      pdata->counter = 0;
    }
  }
  for (iter = 0; iter < solvers[id]->threadnbr; iter++) {
#ifndef FORCE_NOSMP
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    ret = pthread_create(&(solvers[id]->threads[iter]), &attr,
                         control_thread, (void *)&(solvers[id]->threads_data[iter]));
    if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
#endif
  }
  return MURGE_SUCCESS;
}

#ifndef FORCE_NOSMP
#define stop_threads PASTIX_PREFIX_F(stop_threads)
INTS stop_threads(INTS id) {
  int ret;
  INTS iter;
  pthread_mutex_lock(&(solvers[id]->mutex_state));
  solvers[id]->threads_state = MURGE_THREAD_END;
  pthread_mutex_unlock(&(solvers[id]->mutex_state));
  pthread_cond_broadcast(&(solvers[id]->cond_state));
  MURGE_FREE(solvers[id]->threads_data[0].pdata->all_prods);
  for (iter = 0; iter < solvers[id]->threadnbr; iter++) {
    ret = pthread_join(solvers[id]->threads[iter],(void**)NULL);
    if (ret) {errorPrint("thread join."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    MURGE_FREE(solvers[id]->threads_data[iter].pdata->t_prod);
    MURGE_FREE(solvers[id]->threads_data[iter].pdata);
  }
  MURGE_FREE(solvers[id]->threads);
  MURGE_FREE(solvers[id]->threads_data);
  solvers[id]->threadnbr = 0;
  return MURGE_SUCCESS;
}
#endif
/*
 * Function: check_preprocessing
 *
 * Checks if preprocessing (blend) has been called.
 *
 * If it hasn't, it will allocate permutation tabular
 * and call preprocessing step.
 *
 * After calling preprocessing, it will set local number of column
 * and local to global column number tabular to their new values.
 *
 * Colptr and rows will be destroyed because it is obsolete,
 * and state will be set to indicate that preprocessing has been performed.
 *
 * Parameters:
 *   id - Solver instance ID we want to check
 *
 * Returns:
 *   MURGE_ERR_ALLOCATE - If any allocation error occurs.
 */
static inline
int check_preprocessing(int id) {
  pastix_data_t   *pastix_data = solvers[id]->pastix_data;
  PASTIX_INT             *iparm       = pastix_data->iparm;
#ifdef MURGE_TIME
  Clock            clock;
#endif

  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK) &&
      !(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_BLEND_OK))) {
    /* Si il n'est pas fait on effectue le pretraitement */
    if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_SYMB_OK)))
      iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
    else
      iparm[IPARM_START_TASK]  = API_TASK_BLEND;
    iparm[IPARM_END_TASK]    = API_TASK_BLEND;

    CLOCK_INIT;

    if (iparm[IPARM_MATRIX_VERIFICATION] == API_YES) {
      PASTIX_INT err;
#ifdef MEMORY_USAGE
      size_t new_size, old_size = PTR_MEMSIZE(solvers[id]->rows);
#endif
      if (NO_ERR !=
          ((err =
            pastix_checkMatrix_int(pastix_data->pastix_comm,
                                   iparm[IPARM_VERBOSE],
                                   iparm[IPARM_SYM],
                                   API_YES,
                                   solvers[id]->n,
                                   &(solvers[id]->colptr),
                                   &(solvers[id]->rows),
                                   NULL,
                                   &(solvers[id]->l2g),
                                   iparm[IPARM_DOF_NBR],
                                   API_YES)))) {
        errorPrint("pastix_checkMatrix : err %ld\n", (long)err);
        return MURGE_ERR_PARAMETER;
      }
#ifdef MEMORY_USAGE
      new_size = PTR_MEMSIZE(solvers[id]->rows);
      MURGE_TRACE_MALLOC(new_size-old_size, char);
#endif
    }

    pastix_welcome_print(solvers[id]->pastix_data,
                         solvers[id]->colptr,
                         solvers[id]->n);

    if (solvers[id]->perm == NULL) {
      MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, PASTIX_INT);
      memset(solvers[id]->perm, 0, solvers[id]->n*sizeof(PASTIX_INT));
    }

    ALLOC_INVP;

    DPASTIX(&(pastix_data),
            pastix_data->pastix_comm,
            solvers[id]->n,
            solvers[id]->colptr,
            solvers[id]->rows,
            solvers[id]->values,
            solvers[id]->l2g,
            solvers[id]->perm,
            INVP,
            solvers[id]->b,
            solvers[id]->nrhs,
            iparm,
            pastix_data->dparm);
    CLOCK_PRINT("Preprocessing in PaStiX");
#ifdef MURGE_INSERT_DIRECTLY
    /* We build the new CSC which will receive the coefficients */
    {
      PASTIX_INT newn;
      PASTIX_INT * newl2g;
      PASTIX_INT * newcolptr;
      PASTIX_INT * newrows;

      newn = pastix_getLocalNodeNbr(&(pastix_data));
      MURGE_FREE(solvers[id]->perm);
      MURGE_MEMALLOC(solvers[id]->perm, newn, PASTIX_INT);
      memset(solvers[id]->perm, 0, newn*sizeof(PASTIX_INT));
      MURGE_MEMALLOC(newl2g,            newn, PASTIX_INT);
      pastix_getLocalNodeLst(&(solvers[id]->pastix_data),
                             newl2g);
      CLOCK_PRINT("pastix_getLocalNodeLst");

      cscd_redispatch_int(solvers[id]->n,  solvers[id]->colptr,
                          solvers[id]->rows,
                          NULL, NULL, 0,  solvers[id]->l2g,
                          newn,           &newcolptr,           &newrows,
                          NULL, NULL, newl2g, API_YES,
                          solvers[id]->pastix_data->pastix_comm,
                          iparm[IPARM_DOF_NBR]);
      MURGE_TRACE_MALLOC(PTR_MEMSIZE(newcolptr), char);
      MURGE_TRACE_MALLOC(PTR_MEMSIZE(newcolrows), char);
      CLOCK_PRINT("cscd_redispatch_int");

      MURGE_FREE(solvers[id]->l2g);
      if (solvers[id]->g2l)
        MURGE_FREE(solvers[id]->g2l);
      MURGE_FREE(solvers[id]->colptr);
      MURGE_FREE(solvers[id]->rows);

      solvers[id]->n      = newn;
      solvers[id]->N      = -1;
      solvers[id]->l2g    = newl2g;
      solvers[id]->colptr = newcolptr;
      solvers[id]->rows   = newrows;
      solvers[id]->values = NULL;
    }
#else /* not MURGE_INSERT_DIRECTLY */
    /* On corrige n et l2g */

    solvers[id]->n = pastix_getLocalNodeNbr(&(pastix_data));
    MURGE_FREE(solvers[id]->l2g);
    MURGE_FREE(solvers[id]->g2l);
    MURGE_FREE(solvers[id]->perm);
    MURGE_FREE(solvers[id]->colptr);
    MURGE_FREE(solvers[id]->rows);
    MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, PASTIX_INT);
    memset(solvers[id]->perm, 0, solvers[id]->n*sizeof(PASTIX_INT));
    MURGE_MEMALLOC(solvers[id]->l2g, solvers[id]->n, PASTIX_INT);
    pastix_getLocalNodeLst(&(solvers[id]->pastix_data),
                           solvers[id]->l2g);
    solvers[id]->N=-1;

#endif /* not MURGE_INSERT_DIRECTLY */
    /*
     * Building global to local column number array
     */
    cscd_build_g2l(solvers[id]->n,
                   solvers[id]->l2g,
                   solvers[id]->pastix_data->pastix_comm,
                   &solvers[id]->N,
                   &solvers[id]->g2l);
    MURGE_TRACE_MALLOC(PTR_MEMSIZE(solvers[id]->g2l), char);
    MURGE_STATE_TRUE(solvers[id]->state, MURGE_BLEND_OK);
    MURGE_STATE_TRUE(solvers[id]->state, MURGE_SYMB_OK);
    CLOCK_PRINT("cscd_build_g2l");
  }
  return MURGE_SUCCESS;
}

static inline
int check_fact(INTS id) {
  pastix_data_t   *pastix_data = solvers[id]->pastix_data;
  PASTIX_INT      *iparm       = pastix_data->iparm;

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_VALUES_OK))) {
    errorPrint("Need to set values before.");
    return MURGE_ERR_ORDER;
  }
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_RHS_OK))) {
    errorPrint("Need to set right-hand-side member before.");
    return MURGE_ERR_ORDER;
  }
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_FACTO_OK))) {
    if (iparm[IPARM_MATRIX_VERIFICATION] == API_YES) {
      PASTIX_INT err;
#ifdef MEMORY_USAGE
      size_t new_size;
      size_t old_size = PTR_MEMSIZE(solvers[id]->rows) +
	PTR_MEMSIZE(solvers[id]->values);
#endif
      if (NO_ERR !=
	  ((err =
	    pastix_checkMatrix_int(pastix_data->pastix_comm,
				   iparm[IPARM_VERBOSE],
				   iparm[IPARM_SYM],
				   API_YES,
				   solvers[id]->n,
				   &(solvers[id]->colptr),
				   &(solvers[id]->rows),
				   &(solvers[id]->values),
				   &(solvers[id]->l2g),
				   iparm[IPARM_DOF_NBR],
				   API_YES)))) {
	errorPrint("pastix_checkMatrix : err %ld\n", (long)err);
	return MURGE_ERR_PARAMETER;
      }
#ifdef MEMORY_USAGE
      new_size = PTR_MEMSIZE(solvers[id]->rows) +
	PTR_MEMSIZE(solvers[id]->values);
      MURGE_TRACE_MALLOC(new_size-old_size, char);
#endif
    }

    if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_BLEND_OK))) {
      if (!(MURGE_STATE_ISTRUE(solvers[id]->state,
                               MURGE_SYMB_OK)))
        iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
      else
        iparm[IPARM_START_TASK]  = API_TASK_BLEND;
      if (solvers[id]->perm == NULL) {
        MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, PASTIX_INT);
        memset(solvers[id]->perm, 0, solvers[id]->n*sizeof(PASTIX_INT));
      }
      pastix_welcome_print(solvers[id]->pastix_data,
                           solvers[id]->colptr,
                           solvers[id]->n);

      ALLOC_INVP;
      /* Call preprocessing separatly if required so that we can delete CSC
       * before Factorization call
       */
      iparm[IPARM_END_TASK]    = API_TASK_BLEND;
      DPASTIX(&pastix_data,
	      pastix_data->pastix_comm,
	      solvers[id]->n,
	      solvers[id]->colptr,
	      solvers[id]->rows,
	      solvers[id]->values,
	      solvers[id]->l2g,
	      solvers[id]->perm,
	      NULL,
	      solvers[id]->b,
	      solvers[id]->nrhs,
	      pastix_data->iparm,
	      pastix_data->dparm);
    }

    /* If not performed yet we
     * - fill the internal CSC
     * - delete murge CSC and
     * - perform factorization
     */
    PASTIX_FILLIN_CSC(solvers[id]->pastix_data,
		      solvers[id]->pastix_data->pastix_comm,
		      solvers[id]->n,
		      solvers[id]->colptr,
		      solvers[id]->rows,
		      solvers[id]->values,
		      solvers[id]->b,
		      solvers[id]->nrhs,
		      solvers[id]->l2g);
    solvers[id]->pastix_data->cscInternFilled = API_YES;
    if (solvers[id]->pastix_data->iparm[IPARM_FREE_CSCUSER] == API_CSC_FREE) {
      MURGE_FREE(solvers[id]->colptr);
      MURGE_FREE(solvers[id]->rows);
      MURGE_FREE(solvers[id]->values);
    }
    solvers[id]->pastix_data->iparm[IPARM_START_TASK]  =
      API_TASK_NUMFACT;
  } else {
    /* Else, solve is enough */
    if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_SOLVE_DONE))) {
      solvers[id]->pastix_data->iparm[IPARM_START_TASK]  =
        API_TASK_SOLVE;
    } else {
      solvers[id]->pastix_data->iparm[IPARM_START_TASK]  =
        API_TASK_REFINE;
    }
  }

  if (iparm[IPARM_ONLY_RAFF] ==  API_YES) {
    /* When only raff we need to call solve to set RHS */
    iparm[IPARM_END_TASK] = API_TASK_SOLVE;
    DPASTIX(&pastix_data,
            pastix_data->pastix_comm,
            solvers[id]->n,
            solvers[id]->colptr,
            solvers[id]->rows,
            solvers[id]->values,
            solvers[id]->l2g,
            solvers[id]->perm,
            NULL,
            solvers[id]->b,
            solvers[id]->nrhs,
            pastix_data->iparm,
            pastix_data->dparm);
  }
  if (iparm[IPARM_MURGE_REFINEMENT] == API_YES) {
    iparm[IPARM_END_TASK] = API_TASK_REFINE;
    if (solvers[id]->pastix_data->cscInternFilled == API_NO) {
      errorPrint("Trying to refine without internal CSC\n"
                 "\t You need to keep PaStiX internal CSC using"
                 " IPARM_MURGE_MAY_REFINE = API_YES.");
      return MURGE_ERR_PARAMETER;
    }
  } else {
    /* No need to IPARM_FREE_CSCUSER as it is done when facto is performed,
     * after fillin internal CSC, not possible to free it inside PaStiX as it
     * was allocated with memAlloc
     */
    if (iparm[IPARM_MURGE_MAY_REFINE] == API_NO)
      iparm[IPARM_FREE_CSCPASTIX] = API_CSC_FREE;
    iparm[IPARM_END_TASK] = API_TASK_SOLVE;
  }
  DPASTIX(&pastix_data,
          pastix_data->pastix_comm,
          solvers[id]->n,
          solvers[id]->colptr,
          solvers[id]->rows,
          solvers[id]->values,
          solvers[id]->l2g,
          solvers[id]->perm,
          NULL,
          solvers[id]->b,
          solvers[id]->nrhs,
          pastix_data->iparm,
          pastix_data->dparm);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_SYMB_OK);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_SOLVE_DONE);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_REFINE_DONE);
  if (iparm[IPARM_MURGE_REFINEMENT] == API_NO &&
      iparm[IPARM_MURGE_MAY_REFINE] == API_NO) {
    solvers[id]->pastix_data->cscInternFilled = API_NO;
  }
  return MURGE_SUCCESS;
}

static inline
void murge_pastix_fillin_csc(INTS            id,
                             pastix_data_t * data,
                             MPI_Comm        comm,
                             PASTIX_INT      n,
                             PASTIX_INT    * colptr,
                             PASTIX_INT    * rows,
                             PASTIX_FLOAT  * values,
                             PASTIX_FLOAT  * b,
                             PASTIX_INT      nrhs,
                             PASTIX_INT    * l2g) {
  PASTIX_INT tmpN;
  PASTIX_INT   *tmpcolptr = NULL, *tmprows = NULL;
  PASTIX_INT   *tmpperm = NULL, *tmpinvp = NULL;
  PASTIX_FLOAT *tmpvalues = NULL, *tmprhs = NULL;
  INTS has_values = API_NO;
  INTS has_rhs = API_NO;

  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_VALUES_OK))
    has_values = API_YES;
  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_RHS_OK))
    has_rhs = API_NO; /* We always enter NULL ptr here... */


  if (data->procnum == 0 && data->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    errorPrintW("To get an optimal MURGE installation"
                " please build PaStiX with -DDISTRIBUTED\n");

  cscd2csc_int(n, colptr, rows, values,
               b, NULL, NULL,
               &tmpN, &tmpcolptr, &tmprows,
               (has_values==API_NO)?NULL:&tmpvalues,
               (has_rhs==API_NO)?NULL:&tmprhs, NULL, NULL,
               l2g, comm, data->iparm[IPARM_DOF_NBR], API_YES);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpcolptr), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprows), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpvalues), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpperm), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpinvp), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprhs), char);
  pastix_fillin_csc(data,
                    comm,
                    tmpN,
                    tmpcolptr,
                    tmprows,
                    tmpvalues,
                    tmprhs,
                    nrhs,
                    NULL);
  if (NULL != tmpcolptr) MURGE_FREE(tmpcolptr);
  if (NULL != tmprows) MURGE_FREE(tmprows);
  if (NULL != values) MURGE_FREE(tmpvalues);
  if (NULL != tmpperm) MURGE_FREE(tmpperm);
  if (NULL != tmpinvp) MURGE_FREE(tmpinvp);
  if (NULL != b) MURGE_FREE(tmprhs);
}

static inline
void murge_dpastix(INTS             id,
                   pastix_data_t ** data,
                   MPI_Comm         comm,
                   PASTIX_INT       n,
                   PASTIX_INT     * colptr,
                   PASTIX_INT     * rows,
                   PASTIX_FLOAT   * values,
                   PASTIX_INT     * l2g,
                   PASTIX_INT     * perm,
                   PASTIX_INT     * invp,
                   PASTIX_FLOAT   * b,
                   PASTIX_INT       nrhs,
                   PASTIX_INT     * iparm,
                   double         * dparm) {
  PASTIX_INT tmpN;
  PASTIX_INT   *tmpcolptr = NULL, *tmprows = NULL;
  PASTIX_INT   *tmpperm = NULL, *tmpinvp = NULL;
  PASTIX_FLOAT *tmpvalues = NULL, *tmprhs = NULL;
  PASTIX_INT save_veri = iparm[IPARM_MATRIX_VERIFICATION];
  PASTIX_INT save_free = iparm[IPARM_FREE_CSCUSER];
  INTS i;
  INTS has_values = API_NO;
  INTS has_rhs = API_NO;
  INTS has_perm = API_NO;

  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_VALUES_OK))
    has_values = API_YES;
  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_RHS_OK))
    has_rhs = API_YES;
  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK))
    has_perm = API_YES;

  if ((*data)->procnum == 0 && (*data)->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    errorPrintW("To get an optimal MURGE installation"
                " please build PaStiX with -DDISTRIBUTED\n");

  cscd2csc_int(n, colptr, rows, values,
               b, perm, NULL,
               &tmpN, &tmpcolptr, &tmprows,
               (has_values==API_NO)?NULL:&tmpvalues,
               (has_rhs==API_NO)?NULL:&tmprhs,
               (has_perm==API_NO)?NULL:&tmpperm, &tmpinvp,
               l2g, comm, (*data)->iparm[IPARM_DOF_NBR], API_YES);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpcolptr), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprows), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpvalues), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpperm), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpinvp), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprhs), char);

  iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
  iparm[IPARM_FREE_CSCUSER] = API_CSC_PRESERVE;
  pastix(data, comm, tmpN, tmpcolptr, tmprows, tmpvalues,
         tmpperm, tmpinvp, tmprhs, nrhs, iparm, dparm);
  iparm[IPARM_MATRIX_VERIFICATION] = save_veri;
  iparm[IPARM_FREE_CSCUSER] = save_free;

  if (has_perm == API_YES)
    for (i = 0; i < n;i++)
      perm[i] = tmpperm[l2g[i]-1];
  if (has_rhs == API_YES)
    for (i = 0; i < n;i++)
      b[i] = tmprhs[l2g[i]-1];
  if (NULL != tmpcolptr) MURGE_FREE(tmpcolptr);
  if (NULL != tmprows) MURGE_FREE(tmprows);
  if (NULL != values) MURGE_FREE(tmpvalues);
  if (NULL != tmpperm) MURGE_FREE(tmpperm);
  if (NULL != tmpinvp) MURGE_FREE(tmpinvp);
  if (NULL != b) MURGE_FREE(tmprhs);
}


/*******************************************************************************
 * Group: Solver setup functions
 */

/*
 * Function: MURGE_GetSolver
 *
 * returns MURGE_SOLVER_PASTIX
 */
INTS MURGE_GetSolver(INTS * solver_id) {
  *solver_id = MURGE_SOLVER_PASTIX;
  return MURGE_SUCCESS;
}

static inline
INTS _MURGE_InitId(INTS murge_id) {
  
  solvers[murge_id] = (murge_data_t*)malloc(sizeof(murge_data_t));
  solvers[murge_id]->n           = 0;
  solvers[murge_id]->N           = 0;
  solvers[murge_id]->colptr      = NULL;
  solvers[murge_id]->rows        = NULL;
  solvers[murge_id]->values      = NULL;
  solvers[murge_id]->l2g         = NULL;
  solvers[murge_id]->g2l         = NULL;
  solvers[murge_id]->perm        = NULL;
#ifdef CENTRALISED
  solvers[murge_id]->invp        = NULL;
#endif
  solvers[murge_id]->b           = NULL;
  solvers[murge_id]->nrhs        = 1;
  solvers[murge_id]->state       = MURGE_INIT_OK;
  solvers[murge_id]->pastix_data = NULL;
  solvers[murge_id]->sequences   = NULL;
  solvers[murge_id]->seq_ID      = 0;
  solvers[murge_id]->ndump       = 0;
  solvers[murge_id]->dropmask    = NULL;
  solvers[murge_id]->droprows    = NULL;
  solvers[murge_id]->dropcols    = NULL;
  solvers[murge_id]->malloc_size = 0;
  solvers[murge_id]->malloc_maxsize = 0;
  solvers[murge_id]->threadnbr   = 0;
  solvers[murge_id]->threads_state = 0;
#ifndef FORCE_NOSMP
  solvers[murge_id]->barrier.instance         = 0;
  solvers[murge_id]->barrier.blocked_threads  = 0;
  pthread_mutex_init(&(solvers[murge_id]->barrier.sync_lock), NULL);
  pthread_cond_init(&(solvers[murge_id]->barrier.sync_cond), NULL);
#endif
  pastix_task_init(&(solvers[murge_id]->pastix_data), MPI_COMM_WORLD, NULL, NULL);
  return MURGE_SUCCESS;
}


/*
 * Function: MURGE_Initialize
 *
 * Allocate the instance arrays which will keeps intern data for all
 * solver instances.
 *
 * If user is creating several threads calling the solver, this function
 * has to be called before creating threads to insure solver is thread safe.
 *
 * Parameters:
 *   idnbr - Maximum number of solver instances that will be
 *           launched.
 *
 * Returns:
 *   MURGE_SUCCESS      - If function runned successfully.
 *   MURGE_ERR_ALLOCATE - If for some reason, allocation was not
 *                        successfull.
 */
INTS MURGE_Initialize(INTS id_nbr) {
  INTS i;

  print_debug(DBG_MURGE, ">> Murge_Initialize\n");

  if (sizeof(COEF) != sizeof(PASTIX_FLOAT)) {
    errorPrint("Incompatible coefficient type\n");
    return MURGE_ERR_PARAMETER;
  }

  if ( (solvers != NULL) ) {
    errorPrint("MURGE_Initialize has been already called");
    return MURGE_ERR_ORDER;
  }

  idnbr = id_nbr;

  solvers = (murge_data_t**)malloc(idnbr*sizeof(murge_data_t*));

  for (i=0; i< idnbr; i++) {
    solvers[i] = NULL;
    solvers[i] = (murge_data_t*)malloc(sizeof(murge_data_t));

    solvers[i]->n           = 0;
    solvers[i]->N           = 0;
    solvers[i]->colptr      = NULL;
    solvers[i]->rows        = NULL;
    solvers[i]->values      = NULL;
    solvers[i]->l2g         = NULL;
    solvers[i]->g2l         = NULL;
    solvers[i]->perm        = NULL;
#ifdef CENTRALISED
    solvers[i]->invp        = NULL;
#endif
    solvers[i]->b           = NULL;
    solvers[i]->nrhs        = 1;
    solvers[i]->state       = MURGE_INIT_OK;
    solvers[i]->pastix_data = NULL;
    solvers[i]->sequences   = NULL;
    solvers[i]->seq_ID      = 0;
    solvers[i]->ndump       = 0;
    solvers[i]->dropmask    = NULL;
    solvers[i]->droprows    = NULL;
    solvers[i]->dropcols    = NULL;
    solvers[i]->malloc_size = 0;
    solvers[i]->malloc_maxsize = 0;
    solvers[i]->threadnbr   = 0;
    solvers[i]->threads_state = 0;
#ifndef FORCE_NOSMP
    solvers[i]->barrier.instance         = 0;
    solvers[i]->barrier.blocked_threads  = 0;
    pthread_mutex_init(&(solvers[i]->barrier.sync_lock), NULL);
    pthread_cond_init(&(solvers[i]->barrier.sync_cond), NULL);
#endif
    pastix_task_init(&(solvers[i]->pastix_data), MPI_COMM_WORLD, NULL, NULL);
  }

  print_debug(DBG_MURGE, "<< Murge_Initialize\n");
  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_SetDefaultOptions
 *
 * Create a solver instance if not created yet.
 *
 * Sets default options, for solver instance number *id*.
 *
 * The default option set correspond to *stratnum* strategy ID,
 * depending on the solver.
 *
 * Needs <MURGE_Initialize> to be called before
 * to allocate solver instances array.
 *
 * Parameters:
 *   id       - Solver instance identification number.
 *   stratnum - Strategy for the default option Set.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_Initialize> was not called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *stratnum* is not valid.
 *   MURGE_ERR_ALLOCATE  - If couldn't create solver instance.
 */
INTS MURGE_SetDefaultOptions(INTS id, INTS stratnum) {
  print_debug(DBG_MURGE, ">> MURGE_SetDefaultOptions\n");
  CHECK_SOLVER_ID(id);

  if (solvers[id]->pastix_data->iparm == NULL) {
    MURGE_MEMALLOC_EXT(solvers[id]->pastix_data->iparm,
		       IPARM_SIZE, PASTIX_INT);
    MURGE_MEMALLOC_EXT(solvers[id]->pastix_data->dparm,
		       DPARM_SIZE, double);
#ifdef MURGE_THREADSAFE
    pthread_mutex_init(&solvers[id]->mutex_tmpmatrix, NULL);
#endif
  }

  pastix_initParam(solvers[id]->pastix_data->iparm,
                   solvers[id]->pastix_data->dparm);

  solvers[id]->pastix_data->iparm[IPARM_PID] =
    solvers[id]->pastix_data->pastix_id;
  solvers[id]->pastix_data->iparm[IPARM_RHS_MAKING] = API_RHS_B;
  solvers[id]->pastix_data->iparm[IPARM_FREE_CSCUSER] = API_YES;
  return MURGE_SUCCESS;
}


/*
 * Function: MURGE_SetOptionINT
 *
 * Sets integer option, indicated by *number*, to *value* for the
 * solver instance number *id*.
 *
 * Needs <MURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   number - Identification of the integer parameter.
 *   value  - value to set the parameter to.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_SetDefaultOption> was not
 *                         called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
 *                         *number* or *value* are not valid.
 *
 */
INTS MURGE_SetOptionINT (INTS id, INTS number, INTS value) {
  PASTIX_INT murge_param[64];
  PASTIX_INT * iparm = NULL;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  iparm = solvers[id]->pastix_data->iparm;

  murge_param[MURGE_IPARAM_BASEVAL       - 1024] =  IPARM_BASEVAL;
  murge_param[MURGE_IPARAM_DOF           - 1024] =  IPARM_DOF_NBR;
  murge_param[MURGE_IPARAM_SYM           - 1024] =  IPARM_SYM;

  if (number >= 1024) {
    if (number == MURGE_IPARAM_SYM) {
      if (value == MURGE_BOOLEAN_TRUE) {
        value = API_SYM_YES;
        if (iparm[IPARM_FACTORIZATION] == API_FACT_LU) {
          iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
        }
      } else {
        if (value == MURGE_BOOLEAN_FALSE) {
          value = API_SYM_NO;
          if (iparm[IPARM_FACTORIZATION] != API_FACT_LU) {
            iparm[IPARM_FACTORIZATION] = API_FACT_LU;
          }
        } else {
          errorPrint("Invalid value");
          return MURGE_ERR_PARAMETER;
        }
      }
    }
    number = murge_param[number-1024];
  }

  if (!( number < IPARM_SIZE )) {
    errorPrint("option number '%d' is too big", number);
    return MURGE_ERR_PARAMETER;
  }

  if (number < 0) {
    errorPrint("option number '%d' is negative", number);
    return MURGE_ERR_PARAMETER;
  }

  /* TODO : Est-ce qu'on ajoute des tests sur les valeurs rentrées ???? */
  iparm[number] = (PASTIX_INT)value;
  return MURGE_SUCCESS;
}

/*
 Function: MURGE_SetOptionREAL

 Sets real option, indicated by *number*, to *value* for the
 solver instance number *id*.

 Needs <MURGE_SetDefaultOption> to be called before to initiate
 solver instance data.

 Parameters:
 id     - Solver instance identification number.
 number - Identification of the integer parameter.
 value  - value to set the parameter to.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_ORDER     - If <MURGE_SetDefaultOption> was not
 called before.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
 *number* or *value* are not valid.

 */
INTS MURGE_SetOptionREAL(INTS id, INTS number, REAL value) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  PASTIX_INT murge_param[64];

  murge_param[MURGE_RPARAM_EPSILON_ERROR- 1024] =  DPARM_EPSILON_REFINEMENT;

  if (number >= 1024) {
    number = murge_param[number-1024];
  }

  if (!( number < DPARM_SIZE )) {
    errorPrint("number is too big");
    return MURGE_ERR_PARAMETER;
  }

  if (number < 0) {
    errorPrint("number is negative");
    return MURGE_ERR_PARAMETER;
  }

  /* TODO : Est-ce qu'on ajoute des tests sur les valeurs rentrées ???? */
  solvers[id]->pastix_data->dparm[number] = (double)value;
  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_SetCommunicator
 *
 * Sets MPI communicator for the given solver instance.
 *
 * Needs <MURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Musn't be called before <MURGE_SAVE>, <MURGE_LOAD>,
 * <MURGE_GetLocalNodeNbr> nor <MURGE_GetLocalUnknownNbr>
 * because the solver as to be runned with the same MPI
 * communicator all along.
 *
 * If this function is not called, MPI communicator will be
 * *MPI_COMM_WORLD*.
 *
 * This function may not exist if the solver
 * has been compiled without MPI.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   mpicomm - MPI communicator to affect the solver to.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_SetDefaultOption> was not
 *                         called before or if it is called after
 *                         the solver starts its computing tasks.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
 *                         *number* or *value* are not valid.
 */
INTS MURGE_SetCommunicator(INTS id, MPI_Comm mpicomm) {
  CHECK_SOLVER_ID(id);
  solvers[id]->pastix_data->pastix_comm = mpicomm;
  solvers[id]->pastix_data->inter_node_comm = mpicomm;
  solvers[id]->pastix_data->intra_node_comm = MPI_COMM_SELF;
  MPI_Comm_size((solvers[id]->pastix_data)->inter_node_comm,
                &((solvers[id]->pastix_data)->inter_node_procnbr));
  MPI_Comm_rank((solvers[id]->pastix_data)->inter_node_comm,
                &((solvers[id]->pastix_data)->inter_node_procnum));
  MPI_Comm_size((solvers[id]->pastix_data)->intra_node_comm,
                &((solvers[id]->pastix_data)->intra_node_procnbr));
  MPI_Comm_rank((solvers[id]->pastix_data)->intra_node_comm,
                &((solvers[id]->pastix_data)->intra_node_procnum));

  MPI_Comm_size(mpicomm, &(solvers[id]->pastix_data->procnbr));
  MPI_Comm_rank(mpicomm, &(solvers[id]->pastix_data->procnum));
  return MURGE_SUCCESS;
}


/*******************************************************************************
 * Group: I/O functions
 */

/*
 * Function: MURGE_Save
 *
 * Runs preprocessing step, if not done yet, and save the result to disk,
 * into *directory*, so that it can be resume using <MURGE_Load>.
 *
 * Needs <MURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   directory - Path to the directory where to save the solver step.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_SetDefaultOption> was not
 *                         called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 *   MURGE_ERR_IO        - If file(s) couldn't be writen.
 */
INTS MURGE_Save(INTS id, char* directory) {
  char * dest    = NULL;
  char * src     = NULL;
  int    procnum;
  FILE * stream  = NULL;

  CHECK_SOLVER_ID(id);
  procnum = (int)solvers[id]->pastix_data->procnum;

  if (solvers[id]->pastix_data->iparm == NULL) {
    errorPrint("You need to call MURGE_SetDefaultOptions before");
    return MURGE_ERR_ORDER;
  }
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("You need to set graph before");
    return MURGE_ERR_ORDER;
  }

  solvers[id]->pastix_data->iparm[IPARM_IO_STRATEGY] = API_IO_SAVE;
  solvers[id]->pastix_data->iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
  solvers[id]->pastix_data->iparm[IPARM_END_TASK]    = API_TASK_SYMBFACT;

  if (NULL == solvers[id]->perm) {
    MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, PASTIX_INT);
    memset(solvers[id]->perm, 0, solvers[id]->n*sizeof(PASTIX_INT));
  }
  pastix_welcome_print(solvers[id]->pastix_data,
                       solvers[id]->colptr,
                       solvers[id]->n);
  DPASTIX(&(solvers[id]->pastix_data),
          solvers[id]->pastix_data->pastix_comm,
          solvers[id]->n,
          solvers[id]->colptr,
          solvers[id]->rows,
          solvers[id]->values,
          solvers[id]->l2g,
          solvers[id]->perm,
          NULL,
          solvers[id]->b,
          solvers[id]->nrhs,
          solvers[id]->pastix_data->iparm,
          solvers[id]->pastix_data->dparm);

  MURGE_MEMALLOC(dest, (strlen(directory)+20), char);
  MURGE_MEMALLOC(src, 20, char);

  if (procnum == 0) {
    sprintf(dest,"%s/ordergen", directory);
    RENAME("ordergen", dest);
    sprintf(dest,"%s/symbgen", directory);
    RENAME("symbgen", dest);

    sprintf(dest,"%s/murge.save", directory);
    PASTIX_FOPEN(stream, dest, "w");
    fclose(stream);

  }

  MURGE_FREE(src);
  MURGE_FREE(dest);

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_BLEND_OK);
  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_Load
 *
 * Loads preprocessing result from disk, into *directory*,
 * where it had been saved by <MURGE_Save>.
 *
 * If preprocessing data was already computed or loaded, it will
 * be overwriten.
 *
 * Needs <MURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   directory - Path to the directory where to load the solver
 *               preprocessing data.
 *
 * In Fortran, *STR_LEN* is the length of the string directory.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_SetDefaultOption> was not
 *                         called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 *   MURGE_ERR_IO        - If file(s) couldn't be read.
 */
INTS MURGE_Load(INTS id, char* directory) {
  char * src     = NULL;
  int    procnum;

  CHECK_SOLVER_ID(id);
  procnum = (int)solvers[id]->pastix_data->procnum;

  if (solvers[id]->pastix_data->iparm == NULL) {
    errorPrint("You need to call MURGE_SetDefaultOptions before");
    return MURGE_ERR_ORDER;
  }

  solvers[id]->pastix_data->iparm[IPARM_IO_STRATEGY] = API_IO_LOAD;
  solvers[id]->pastix_data->iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
  solvers[id]->pastix_data->iparm[IPARM_END_TASK]    = API_TASK_BLEND;

  pastix_welcome_print(solvers[id]->pastix_data,
                       solvers[id]->colptr,
                       solvers[id]->n);

  if (NULL == solvers[id]->perm) {
    MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, PASTIX_INT);
    memset(solvers[id]->perm, 0, solvers[id]->n*sizeof(PASTIX_INT));
  }
  MURGE_MEMALLOC(src, (strlen(directory)+20), char);
  if (procnum == 0) {
    sprintf(src,"%s/ordergen", directory);
    LINK(src, "ordername");
    sprintf(src,"%s/symbgen", directory);
    LINK(src, "symbname");
  }
  MPI_Barrier(solvers[id]->pastix_data->pastix_comm);
  MURGE_FREE(src);

  DPASTIX(&(solvers[id]->pastix_data),
          solvers[id]->pastix_data->pastix_comm,
          solvers[id]->n,
          solvers[id]->colptr,
          solvers[id]->rows,
          solvers[id]->values,
          solvers[id]->l2g,
          solvers[id]->perm,
          NULL,
          solvers[id]->b,
          solvers[id]->nrhs,
          solvers[id]->pastix_data->iparm,
          solvers[id]->pastix_data->dparm);
  solvers[id]->N = solvers[id]->pastix_data->ordemesh.rangtab[solvers[id]->pastix_data->ordemesh.cblknbr];
  if (procnum == 0) {
    UNLINK("ordername");
    UNLINK("symbname");
  }


  MURGE_STATE_TRUE(solvers[id]->state, MURGE_SYMB_OK);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_GRAPH_OK);
  return MURGE_SUCCESS;
}


/*******************************************************************************
 * Group: Getting solver's distribution
 */

/*
 * Function: MURGE_GetLocalNodeNbr
 *
 * Computes preprocessing step, if not done, and the number of
 * Nodes in the new ditribution of the matrix.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   nodenbr   - *INTS* where to store number of nodes.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *nodenbr* is *NULL* (can occur in C).
 */
INTS MURGE_GetLocalNodeNbr    (INTS id, INTS *nodenbr){
#ifdef MURGE_TIME
  Clock clock;
#endif
  CLOCK_INIT;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);

  *nodenbr = (INTS)pastix_getLocalNodeNbr(&(solvers[id]->pastix_data));
  CLOCK_PRINT("CALL pastix_getLocalNodeNbr");

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_NODENBR_OK);
  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_GetLocalNodeList
 *
 * Computes the local node list, corresponding to
 * the new distribution, after preprocessing.
 *
 * *nodelist* array has to be allocated before calling
 * this function.
 *
 * As it's result determines the size of *nodelist*
 * array, <MURGE_GetLocalNodeNbr> should be run before it.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   nodelist  - Array where to store the list of local nodes.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - if <MURGE_GetLocalNodeNbr> has not been called
 *                         before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *nodelist* is *NULL* (can occur in C).
 */
INTS MURGE_GetLocalNodeList   (INTS id, INTS *nodelist) {
  int ret;
  PASTIX_INT nodenbr = 0;
  PASTIX_INT i;
  PASTIX_INT * intern_nodelist;
#ifdef MURGE_TIME
  Clock clock;
#endif
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_NODENBR_OK))) {
    errorPrint("You need to call MURGE_GetLocalNodeNbr before");
    return MURGE_ERR_ORDER;
  }
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_NODELST_OK);

  nodenbr = pastix_getLocalNodeNbr(&(solvers[id]->pastix_data));
  if (sizeof(PASTIX_INT) != sizeof(INTS)) {
    MURGE_MEMALLOC(intern_nodelist, nodenbr, PASTIX_INT);
  }
  else
    {
      intern_nodelist = (PASTIX_INT*)nodelist;
    }

    CLOCK_INIT;
    if (EXIT_SUCCESS != ( ret =
                        pastix_getLocalNodeLst(&(solvers[id]->pastix_data),
                                               intern_nodelist)))
    return MURGE_ERR_SOLVER;
    CLOCK_PRINT("CALL pastix_getLocalNodeLst");

  if (sizeof(PASTIX_INT) != sizeof(INTS)) {
    for (i = 0; i < nodenbr; i++) {
      nodelist[i] = (INTS) intern_nodelist[i];
    }
    MURGE_FREE(intern_nodelist);
  }
  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
    for (i = 0; i < nodenbr; i++)
      nodelist[i] -= 1;
  }

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_NODELST_OK);
  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_GetLocalUnkownNbr
 *
 * Computes preprocessing step, if not done, and the number of
 * Unkowns in the new ditribution of the matrix.
 *
 * Parameters:
 *   id            - Solver instance identification number.
 *   unkownnbr     - *INTS* where to store number of unkowns.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *unkownnbr* is *NULL* (can occur in C).
 */
INTS MURGE_GetLocalUnknownNbr (INTS id, INTS *unkownnbr) {
#ifdef MURGE_TIME
  Clock clock;
#endif
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CLOCK_INIT;
  *unkownnbr = (INTS)pastix_getLocalUnknownNbr(&(solvers[id]->pastix_data));
  CLOCK_PRINT("CALL pastix_getLocalUnknownNbr");
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_NODENBR_OK);
  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_GetLocalUnkownList
 *
 * Computes the local unkown list, corresponding to
 * the new distribution, after preprocessing.
 *
 * *unkownlist* array has to be allocated before calling
 * this function.
 *
 * As it's result determines the size of *unkownlist*
 * array, <MURGE_GetLocalUnkownNbr> should be run before it.
 *
 * Parameters:
 *   id          - Solver instance identification number.
 *   unkownlist  - Array where to store the list of local unkowns.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - if <MURGE_GetLocalUnkownNbr> has not been called
 *                         before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *unkownlist* is *NULL* (can occur in C).
 */
INTS MURGE_GetLocalUnknownList(INTS id, INTS *unkownlist){
  int ret;
  PASTIX_INT nodenbr = 0;
  PASTIX_INT i;
  PASTIX_INT * intern_nodelist = NULL;
#ifdef MURGE_TIME
  Clock clock;
#endif
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_NODENBR_OK))) {
    errorPrint("You need to call MURGE_GetLocalNodeNbr before");
    return MURGE_ERR_ORDER;
  }

  nodenbr = pastix_getLocalUnknownNbr(&(solvers[id]->pastix_data));
  if (sizeof(PASTIX_INT) != sizeof(INTS)) {
    MURGE_MEMALLOC(intern_nodelist, nodenbr, PASTIX_INT);
  }
  else
    {
      intern_nodelist = (PASTIX_INT*)unkownlist;
    }
  CLOCK_INIT;
  if (EXIT_SUCCESS != ( ret =
                        pastix_getLocalUnknownLst(&(solvers[id]->pastix_data),
                                                  intern_nodelist)))
    return MURGE_ERR_SOLVER;
  CLOCK_PRINT("CALL pastix_getLocalUnknownLst");

  if (sizeof(PASTIX_INT) != sizeof(INTS)) {
    for (i = 0; i < nodenbr; i++) {
      unkownlist[i] = (INTS) intern_nodelist[i];
    }
    MURGE_FREE(intern_nodelist);
  }
  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
    for (i = 0; i < nodenbr; i++)
      unkownlist[i] -= 1;
  }
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_NODELST_OK);
  return MURGE_SUCCESS;
}



/*******************************************************************************
 * Group: Graph setup functions
 */

/*
 * Function: MURGE_GraphBegin
 *
 * - Allocate temporary structure which will contain graph entries.
 * - Set number of unkowns in the graph.
 * - Set the number of entries that are expected in this building session.
 * - Reset the number of entries for this build session.
 * - Set all states except MURGE_GRAPH_BUILD to FALSE (graph, values, blend,
 * nodelst, nodenbr, facto)
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   N       - Number of unkowns.
 *   edgenbr - Number of edges in this building session.
 *             If edgenbr is negative, PaStiX will perform dynamic
 *             reallocation of the array, with the first allocation of
 *             size -edgenbr.
 *
 * Returns:
 *   MURGE_ERR_ORDER     - MURGE_GraphBegin has already been called, or if
 *                         *solvers* or *solvers[id]* are not allocated,
 *                         or if *iparm* or *dparm* are not allocated.
 *   MURGE_ERR_PARAMETER - If *id* is not in correct range.
 *   MURGE_SUCCESS       - Otherwise.
 */
INTS MURGE_GraphBegin(INTS id, INTS N, INTL edgenbr) {
  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD)) {
    errorPrint("MURGE_GraphBegin has been called before");
    return MURGE_ERR_ORDER;
  }

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#endif
  if (edgenbr < 0) {
    edgenbr = -edgenbr;
    solvers[id]->dynamic = API_YES;
  }
  else {
    solvers[id]->dynamic = API_NO;
  }
#ifdef HASH_MTX
  MURGE_MEMALLOC(solvers[id]->hashgraphtab, edgenbr, hash_graph_entry_t);
  solvers[id]->hashgraphtab_size = edgenbr;
#endif
  vcsc_init(&solvers[id]->vcsc, N, edgenbr, 0, id);

  solvers[id]->N        = N;
  solvers[id]->edgenbr  = edgenbr;
  solvers[id]->cnt      = 0;
  solvers[id]->cnt_zero = 0;

#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODELST_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_FACTO_OK);

  MURGE_STATE_TRUE(solvers[id]->state,  MURGE_GRAPH_BUILD);
  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_GraphSetEdge
 *
 * - Check that the number of entries has not been reach for
 * this session.
 * - Increments ROW and COL if baseval is set to 0.
 * - Checks that ROW and COL ranges are corrects.
 * - Adds an entry to the temporary ijv structure.
 *
 * Parameters:
 *   id  - Solver instance identification number.
 *   ROW - Row of the entry.
 *   COL - Column of the entry.
 *
 * Return:
 *   MURGE_ERR_ORDER     - if we are not in a graph building session, or if
 *                         two many edges have been entered, or if
 *                         *solvers* or *solvers[id]* are not allocated,
 *                         or if *iparm* or *dparm* are not allocated.
 *   MURGE_ERR_PARAMETER - *ROW* or *COL* are out of range or if *id* is not
 *                         in correct range.
 *   MURGE_SUCCESS       - Otherwise
 */
INTS MURGE_GraphSetEdge (INTS id, INTS ROW, INTS COL) {
  return MURGE_GraphSetEdge_ (id, ROW, COL);
}
INTS MURGE_GraphEdge (INTS id, INTS ROW, INTS COL) {
  return MURGE_GraphSetEdge_ (id, ROW, COL);
}
static inline
INTS MURGE_GraphSetEdge_ (INTS id, INTS ROW, INTS COL) {

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD))) {
    errorPrint("Need to call MURGE_GraphBegin first");
    return MURGE_ERR_ORDER;
  }
  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
    COL += 1;
    ROW += 1;
  }

  if (solvers[id]->dropmask != NULL && ROW != COL)
    if (solvers[id]->dropmask[(COL-1)] && solvers[id]->dropmask[(ROW-1)]) {
      solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }
  if (solvers[id]->droprows != NULL && ROW != COL)
    if (solvers[id]->droprows[(ROW-1)]) {
      solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }
  if (solvers[id]->dropcols != NULL && ROW != COL)
    if (solvers[id]->dropcols[(COL-1)]) {
      solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }

  if (ROW < 1 || COL < 1 || ROW > solvers[id]->N || COL > solvers[id]->N) {
    errorPrint("ROW %ld or COL %ld is out of range [%ld-%ld]",
               (long)ROW, (long)COL, (long)1, (long)solvers[id]->N);
    return MURGE_ERR_PARAMETER;
  }

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#endif
  if (solvers[id]->cnt + solvers[id]->cnt_zero + 1 > solvers[id]->edgenbr) {
    if (solvers[id]->dynamic == API_NO) {
      errorPrint("Too many edges in graph description (%ld > %ld)",
                 solvers[id]->cnt + solvers[id]->cnt_zero + 1, solvers[id]->edgenbr);
#ifdef MURGE_THREADSAFE
      pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif

      return MURGE_ERR_ORDER;
    }
  }
  COEF * VAL = NULL;
  vcsc_add_node(solvers[id]->vcsc, COL, ROW, VAL, NULL, id);
  solvers[id]->cnt++;
#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
  return MURGE_SUCCESS;
}

INTS MURGE_GraphSetBlockEdges(INTS id, INTS nROW, INTS *ROWlist,
                              INTS nCOL, INTS *COLlist) {
  INTS r, c, ret;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  for (r = 0; r < nROW; r++)
    for (c = 0; c < nCOL; c++)
      if (MURGE_SUCCESS !=
          (ret = MURGE_GraphSetEdge_(id, ROWlist[r], COLlist[c])))
        return ret;
  return MURGE_SUCCESS;
}


/*
 * Function: MURGE_GraphEnd
 *
 * - Sort temporary IJV structure with cols as key.
 * - Distribute columns onto processors.
 * (first column on first proc and so on...)
 * - Build a distributed CSC that will be given to PaStiX.
 *
 * TODO:
 * - In the case of a triangular matrix, count each extra-diagonal twice.
 * - Use initial distribution to compute column distribution,
 * in order to reduce communications.
 *
 * Parameters :
 *   id  - Solver instance identification number.
 *
 * Returns:
 *   MURGE_ERR_ORDER     - if we are not in a graph building session, or if
 *                         all edges have not been entered, or if
 *                         *solvers* or *solvers[id]* are not allocated,
 *                         or if *iparm* or *dparm* are not allocated.
 *   MURGE_ERR_PARAMETER - *ROW* or *COL* are out of range or if *id* is not
 *                         in correct range.
 *   MURGE_SUCCESS       - Otherwise
 */
INTS MURGE_GraphEnd  (INTS id) {
  pastix_data_t   *pastix_data = solvers[id]->pastix_data;
#ifdef MURGE_TIME
  Clock            clock;
#endif
  CLOCK_INIT;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

#ifdef CENTRALISED
  pastix_data->iparm[IPARM_GRAPHDIST] = API_NO;
#endif

  /*
   * Checking that the function is called at the right time
   */
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD))) {
    errorPrint("Need to call MURGE_GraphBegin first");
    return MURGE_ERR_ORDER;
  }

  if (solvers[id]->cnt_zero != 0) {
    if (solvers[id]->pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
      fprintf(stdout,
              "%ld (%.2g %%) entries were skipped on proc %ld\n",
              (long)solvers[id]->cnt_zero,
              (double)((double)100.0*((double)solvers[id]->cnt_zero)/
                       ((double)(solvers[id]->cnt + solvers[id]->cnt_zero))),
              (long)solvers[id]->pastix_data->procnum);
    }
    else
      {
        if (solvers[id]->pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NO) {
          PASTIX_INT nz_glob;
          PASTIX_INT zeros_glob;
          PASTIX_INT nz = solvers[id]->colptr[solvers[id]->n]-1;
          MPI_Reduce( &(solvers[id]->cnt_zero), &zeros_glob,
                      1, COMM_INT,
                      MPI_SUM, 0, solvers[id]->pastix_data->pastix_comm);
          MPI_Reduce( &nz, &nz_glob,
                      1, COMM_INT,
                      MPI_SUM, 0, solvers[id]->pastix_data->pastix_comm);
          if (solvers[id]->pastix_data->procnum == 0) {
            fprintf(stdout,
                    "%ld entries were skipped"
                    " (from %ld (%.3g%%))\n",
                    (long int)zeros_glob,
                    (long int)nz_glob,
                    100.0*((double)(zeros_glob)/
                           ((double)(nz_glob))));
          }

        }
      }
  }

  if (solvers[id]->dynamic == API_NO &&
      solvers[id]->cnt + solvers[id]->cnt_zero < solvers[id]->edgenbr) {
    errorPrint("Missing edges entry, expected %ld, entered %ld",
               (long)solvers[id]->edgenbr,
               (long)(solvers[id]->cnt + solvers[id]->cnt_zero));
    return MURGE_ERR_ORDER;
  }

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#endif
  if (solvers[id]->cnt + solvers[id]->cnt_zero > solvers[id]->edgenbr) {
    errorPrint("Too many edges entry, expected %ld, entered %ld (%ld dropped)",
               (long)solvers[id]->edgenbr,
               (long)(solvers[id]->cnt + solvers[id]->cnt_zero),
               (long)solvers[id]->cnt_zero);
#ifdef MURGE_THREADSAFE
    pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
    return MURGE_ERR_ORDER;
  }
  vcsc_to_cscd(solvers[id]->vcsc, pastix_data->pastix_comm,
               &(solvers[id]->n), &(solvers[id]->colptr),
               &(solvers[id]->rows), NULL,
               &(solvers[id]->l2g), &(solvers[id]->g2l), NULL, 0, id);
  vcsc_destroy(solvers[id]->vcsc, id);

  MURGE_DUMP_GRAPH;

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_GRAPH_BUILD);

#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
  CLOCK_PRINT("MURGE_GraphEnd");

  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_GraphGlobalCSR
 *
 * Enter the adjency graph in a Column Sparse Row form.
 *
 *
 * If the matrix is symmetric, calls <MURGE_GraphGlobalCSC>
 * else uses <MURGE_GraphBegin>, <MURGE_GraphEdge>,
 * <MURGE_GraphEnd> sequence.
 *
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   N      - Number of rows in the CSR.
 *   rowptr - Indexes of each row in COLS array.
 *   root   - Rank of the processor owning the CSR (-1 for all processors)
 */
INTS MURGE_GraphGlobalCSR(INTS id, INTS N, INTL *rowptr, INTS *COLS, INTS root) {
  PASTIX_INT  iter;
  PASTIX_INT  iter2;
  INTS ret;
  int  baseval;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if ( MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD) ) {
    errorPrint("Do not call MURGE_GraphBegin before");
    return MURGE_ERR_ORDER;
  }

  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  /* Si on a un graph symetrique autant faire un MURGE_GraphGlobalCSC */
  if (solvers[id]->pastix_data->iparm[IPARM_SYM] == API_SYM_YES)
    return MURGE_GraphGlobalCSC(id, N, rowptr, COLS, root);


  if (solvers[id]->pastix_data->procnum == root || root == -1) {
    if (MURGE_SUCCESS != (ret = MURGE_GraphBegin(id, N, rowptr[N]- baseval)))
      return ret;

    for (iter = 0; iter < N; iter++) {
      for (iter2 = rowptr[iter]; iter2 < rowptr[iter+1]; iter2++) {
        ret =
          MURGE_GraphEdge(id,
                          iter + baseval,
                          COLS[iter2 - baseval]);
        if (MURGE_SUCCESS != ret)
          return ret;
      }
    }
    if (MURGE_SUCCESS != (ret = MURGE_GraphEnd(id)))
      return ret;
  }
  else
    {
      if (MURGE_SUCCESS != (ret = MURGE_GraphBegin(id, N, 0)))
        return ret;
      if (MURGE_SUCCESS != (ret = MURGE_GraphEnd(id)))
        return ret;
    }

  MURGE_STATE_TRUE (solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODELST_OK);

  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_GraphGlobalCSC
 *
 * Distribute the CSC on the processors and use it for PaStiX calls.
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   N      - Number of columns in the CSR.
 *   colptr - Indexes of each columns in ROWS array.
 *   root   - Rank of the processor owning the CSR (-1 for all processors)
 */
INTS MURGE_GraphGlobalCSC(INTS id, INTS N, INTL *colptr, INTS *ROWS, INTS root) {
  PASTIX_INT          globaledgenbr;
  PASTIX_INT          averageedgenbr;
  PASTIX_INT          localedgenbr;
  PASTIX_INT          firstcol;
  PASTIX_INT          lastcol         = 0;
  PASTIX_INT          iter;
  PASTIX_INT          procnum;
  PASTIX_INT          firstlast[2];
  INTS        *tmpj            = NULL;
  MPI_Request *requests_fl     = NULL;
  MPI_Request *requests_colptr = NULL;
  MPI_Request *requests_rows   = NULL;
  MPI_Status   status;
  int          baseval;             /* User baseval               */

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if ( MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD) ) {
    errorPrint("Do not call MURGE_GraphBegin before");
    return MURGE_ERR_ORDER;
  }

  solvers[id]->N = N;
  baseval        = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];

  if (root == solvers[id]->pastix_data->procnum || root == -1) {
    globaledgenbr   = colptr[N] - baseval;
    averageedgenbr  = globaledgenbr / solvers[id]->pastix_data->procnbr;
    firstcol        = 0;
    if (root != -1) {
      MURGE_MEMALLOC(requests_fl,
                     solvers[id]->pastix_data->procnbr,
                     MPI_Request);
      MURGE_MEMALLOC(requests_colptr,
                     solvers[id]->pastix_data->procnbr,
                     MPI_Request);
      MURGE_MEMALLOC(requests_rows,
                     solvers[id]->pastix_data->procnbr,
                     MPI_Request);
    }
    /* Pour chaque prrocesseur,
     on attribue au processeur un certain nombre de colonnes et donc
     d'arrêtes

     Si le processeur est local on construit le loc2glob
     On construit le colptr
     On copie rows

     Sinon, on envoi les numéros de première et dernière colonnes et
     les morceaux de colptr et rows correspondants.

     */
    for (procnum = 0; procnum <  solvers[id]->pastix_data->procnbr; procnum++) {
      while (lastcol < N - 1&&
             colptr[lastcol+1] - colptr[firstcol]< averageedgenbr)
        lastcol++;
      localedgenbr =  colptr[lastcol+1] - colptr[firstcol];

      if (procnum == solvers[id]->pastix_data->procnum) {
        solvers[id]->n = lastcol-firstcol+1;

        MURGE_MEMALLOC(solvers[id]->l2g, solvers[id]->n, PASTIX_INT);

        for (iter = 0; iter < solvers[id]->n; iter++) {
          solvers[id]->l2g[iter] = firstcol+iter+1;
        }

        MURGE_MEMALLOC(solvers[id]->colptr, solvers[id]->n+1, PASTIX_INT);

        for (iter = 0; iter < solvers[id]->n+1; iter++) {
          solvers[id]->colptr[iter] =
            colptr[firstcol+iter] - colptr[firstcol]+1;
        }

        MURGE_MEMALLOC(solvers[id]->rows, localedgenbr, PASTIX_INT);

        for (iter = 0; iter < localedgenbr; iter++)
          solvers[id]->rows[iter] = ROWS[colptr[firstcol]+iter-1];

      }
      else
        {
          if (root != -1) {
            firstlast[0] = firstcol;
            firstlast[1] = lastcol;


            MPI_Isend(firstlast,
                      2, COMM_INT, procnum,
                      TAG_FL, solvers[id]->pastix_data->pastix_comm,
                      &requests_fl[procnum]);

            MPI_Isend(&colptr[firstcol],
                      lastcol-firstcol+2,
                      COMM_INT, procnum,
                      TAG_COL, solvers[id]->pastix_data->pastix_comm,
                      &requests_colptr[procnum]);

            MPI_Isend(&ROWS[colptr[firstcol]],
                      localedgenbr*sizeof(INTS),
                      MPI_BYTE, procnum,
                      TAG_ROW, solvers[id]->pastix_data->pastix_comm,
                      &requests_rows[procnum]);

          }
        }
      firstcol = lastcol + 1;
    }
    if (root != -1) {
      for (procnum = 0;
           procnum < solvers[id]->pastix_data->procnbr;
           procnum++) {
        if (procnum != solvers[id]->pastix_data->procnum) {
          MPI_Wait(&requests_fl[procnum], &status);
          MPI_Wait(&requests_colptr[procnum], &status);
          MPI_Wait(&requests_rows[procnum], &status);
        }
      }
    }
    MURGE_FREE(requests_rows);
    MURGE_FREE(requests_colptr);
    MURGE_FREE(requests_fl);
  }
  else
    {
      /* Si on est pas le processeur racine

       On recoit les numeros de première et dernière colonnes
       On en déduit le loca2glob
       On recoit les parties locales de colptr et rows
       On construit le colptr local et rows local.
       */
      MPI_Recv(firstlast,
               2, COMM_INT, root,
               TAG_FL, solvers[id]->pastix_data->pastix_comm,
               &status);
      firstcol = firstlast[0];
      lastcol  = firstlast[1];

      solvers[id]->n = lastcol-firstcol+1;

      MURGE_MEMALLOC(solvers[id]->l2g, lastcol-firstcol+1, PASTIX_INT);

      for (iter = 0; iter < lastcol-firstcol+1; iter++) {
        solvers[id]->l2g[iter] = firstcol+iter;
      }

      MURGE_MEMALLOC(solvers[id]->colptr, lastcol-firstcol+2, PASTIX_INT);

      MPI_Recv(solvers[id]->colptr,
               lastcol-firstcol+2, COMM_INT, root,
               TAG_COL, solvers[id]->pastix_data->pastix_comm,
               &status);


      for (iter = 0; iter < lastcol-firstcol+2; iter++) {
        solvers[id]->colptr[lastcol-firstcol+1 - iter] -=
          solvers[id]->colptr[0];
      }

      localedgenbr = solvers[id]->colptr[lastcol-firstcol+1]-1;

      MURGE_MEMALLOC(tmpj, localedgenbr, INTS);

      MPI_Recv(tmpj,
               localedgenbr*sizeof(INTS), MPI_BYTE, root,
               TAG_ROW, solvers[id]->pastix_data->pastix_comm,
               &status);

      if (sizeof(INTS) == sizeof(PASTIX_INT)) {
        solvers[id]->rows = (PASTIX_INT*)tmpj;
      }
      else
        {
          MURGE_MEMALLOC(solvers[id]->rows, localedgenbr, PASTIX_INT);

          for (iter = 0; iter < localedgenbr; iter++)
            solvers[id]->rows[iter] = (PASTIX_INT)tmpj[iter];
          MURGE_FREE(tmpj);
        }
    }

  MURGE_DUMP_GRAPH;

  MURGE_STATE_TRUE (solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODELST_OK);

  return MURGE_SUCCESS;
}
/*
 * Function: MURGE_GraphGlobalIJV
 *
 * Distribute the graph on the processors, compress the columns
 * array and use the built CSCd to call PaStiX.
 *
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   N      - Number of columns in the CSR.
 *   NNZ    - Number of non-zeros in the matrix.
 *   ROWS   - Rows array.
 *   COLS   - Columns array.
 *   root   - Rank of the processor owning the CSR (-1 for all processors)
 */
INTS MURGE_GraphGlobalIJV(INTS id, INTS N, INTL NNZ, INTS *ROWS,
                          INTS *COLS, INTS root) {
  PASTIX_INT       *localn     = NULL;   /* Number of local column on each proc  */
  PASTIX_INT       *localedges = NULL;   /* Number of local edges on each proc   */
  PASTIX_INT        lnnz;                /* Local number of edges                */
  PASTIX_INT        sizes[2];            /* Array to send n and nnz to each proc */
  INTS      *tmprows  = NULL;     /* Temporary local rows tabular         */
  INTS      *tmpcols  = NULL;     /* Temporary local columns tabular      */
  PASTIX_INT       *sizecols = NULL;     /* Number of rows in each column        */
  PASTIX_INT        totaledgenbr;        /* Total number of edges                */
  PASTIX_INT        avredgenbr;          /* Average number of edges              */
  int        baseval_int     = 1; /* Internal baseval, always 1           */
  int        baseval;             /* User baseval                         */
  PASTIX_INT        iter, iter2, index;  /* Iterators                            */
  PASTIX_INT       *coldist  = NULL;     /* Owner of each column                 */
  int        procnum;             /* Processor number iterator            */
  MPI_Status status;              /* MPI status */
  INTS       currentcol;
  void      *sortptr[2];

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if ( MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD) ) {
    errorPrint("Do not call MURGE_GraphBegin before");
    return MURGE_ERR_ORDER;
  }


  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];


  if ((solvers[id]->pastix_data->procnum == root) || (root == -1)) {
    solvers[id]->N = N;
  }
  if (root != -1) {
    MPI_Bcast(&solvers[id]->N, 1, COMM_INT, root,
              solvers[id]->pastix_data->pastix_comm);
    MPI_Bcast(&NNZ, sizeof(INTL), MPI_BYTE, root,
              solvers[id]->pastix_data->pastix_comm);
  }
  if ((solvers[id]->pastix_data->procnum == root) || (root == -1)) {
    /* Sort Col and Rows with first key, col */
    sortptr[0] = COLS;
    sortptr[1] = ROWS;
#ifdef INTSSIZE64
    qsort2IntAsc(sortptr, NNZ);
#else
    qsort2SmallIntAsc(sortptr, NNZ);
#endif
    /* decide how to distribute the graph */
    MURGE_MEMALLOC(sizecols, solvers[id]->N, PASTIX_INT);
    memset(sizecols, 0, solvers[id]->N*sizeof(PASTIX_INT));
    totaledgenbr = 0;
    /* Count how long is each column */
    for (iter = 0; iter < NNZ; iter ++) {
      /* TODO: Dans le cas ou la matrice donnée ne contient que la partie
       triangulaire, il faudra  compter 2 fois les elements
       non diagonaux.
       */
      sizecols[COLS[iter] - 1]++;
      totaledgenbr++;
    }

    avredgenbr = totaledgenbr/solvers[id]->pastix_data->procnbr;

    MURGE_MEMALLOC(coldist, solvers[id]->N, PASTIX_INT);

    for (iter = 0; iter < solvers[id]->N; iter++)
      coldist[iter] = solvers[id]->pastix_data->procnbr - 1;

    procnum    = 0;
    iter       = 0;

    MURGE_MEMALLOC(localedges, solvers[id]->pastix_data->procnbr, PASTIX_INT);
    MURGE_MEMALLOC(localn,     solvers[id]->pastix_data->procnbr, PASTIX_INT);
    memset(localedges, 0, solvers[id]->pastix_data->procnbr*sizeof(PASTIX_INT));
    memset(localn, 0, solvers[id]->pastix_data->procnbr*sizeof(PASTIX_INT));

    while (iter < solvers[id]->N) {
      localedges[procnum] = 0;
      while (iter < solvers[id]->N &&
             (localedges[procnum] < avredgenbr||
              (procnum == solvers[id]->pastix_data->procnbr -1))) {
        coldist[iter] = procnum;
        localn[procnum]++;
        localedges[procnum] +=  sizecols[iter];
        iter ++;
      }
      procnum++;
    }

    MURGE_FREE(coldist);

    /* Send data to each processor */

    for (index = 0, procnum = 0;
         procnum < solvers[id]->pastix_data->procnbr;
         procnum++) {
      if (procnum != solvers[id]->pastix_data->procnum) {
        if (root != -1) {
          sizes[0] = localn[procnum];
          sizes[1] = localedges[procnum];

          /* envoi du nombre de non zeros */
          MPI_Send(sizes,
                   2, COMM_INT, procnum,
                   TAG_SIZE, solvers[id]->pastix_data->pastix_comm);
          /* envoi des lignes */
          MPI_Send(&(ROWS[index]),
                   localedges[procnum]*sizeof(INTS), MPI_BYTE, procnum,
                   TAG_ROW, solvers[id]->pastix_data->pastix_comm);
          /* envoi des colonnes */
          MPI_Send(&(COLS[index]),
                   localedges[procnum]*sizeof(INTS), MPI_BYTE, procnum,
                   TAG_COL, solvers[id]->pastix_data->pastix_comm);
        }
      }
      else
        {
          tmprows = &(ROWS[index]);
          tmpcols = &(COLS[index]);
        }
      index += localedges[procnum];
    }

    solvers[id]->n = localn[solvers[id]->pastix_data->procnum];
    lnnz = localedges[solvers[id]->pastix_data->procnum];
    MURGE_FREE(localn);
    MURGE_FREE(localedges);
  }
  else
    {
      MPI_Recv(sizes,
               2, COMM_INT, root,
               TAG_SIZE, solvers[id]->pastix_data->pastix_comm,
               &status);
      solvers[id]->n = sizes[0];
      lnnz           = sizes[1];

      MURGE_MEMALLOC(tmprows, lnnz, INTS);
      MPI_Recv(tmprows,
               lnnz*sizeof(INTS), MPI_BYTE, root,
               TAG_ROW, solvers[id]->pastix_data->pastix_comm,
               &status);
      MURGE_MEMALLOC(tmpcols, lnnz, INTS);
      MPI_Recv(tmpcols,
               lnnz*sizeof(INTS), MPI_BYTE, root,
               TAG_COL, solvers[id]->pastix_data->pastix_comm,
               &status);
    }

  MURGE_MEMALLOC(solvers[id]->colptr, solvers[id]->n+1, PASTIX_INT);
  MURGE_MEMALLOC(solvers[id]->l2g,    solvers[id]->n  , PASTIX_INT);
  /* convert tmpcols/tmprows to CSCd */
  iter2=baseval_int;
  for (iter=0; iter<(solvers[id]->n); iter++) {
    solvers[id]->colptr[iter] = iter2;
    solvers[id]->l2g[iter]    = tmpcols[iter2-baseval_int]
      - baseval + baseval_int;
    currentcol = tmpcols[iter2-baseval_int];
    while (((iter2-baseval) < lnnz) &&
           ((tmpcols[iter2-baseval_int]) == (currentcol))) {
      iter2++;
    }
  }

  if ((solvers[id]->pastix_data->procnum != root) && (root != -1))
    MURGE_FREE(tmpcols);

  solvers[id]->colptr[solvers[id]->n] = iter2;

  if (iter2 != lnnz+baseval) {
    errorPrint("Mauvais nombre d'arrête");
    return MURGE_ERR_PARAMETER;
  }


  MURGE_MEMALLOC(solvers[id]->rows, lnnz, PASTIX_INT);

  for (iter=0; iter<lnnz; iter++) {
    solvers[id]->rows[iter] = tmprows[iter];
  }
  if (!((solvers[id]->pastix_data->procnum == root) || (root == -1)))
    MURGE_FREE(tmprows);

  MURGE_STATE_TRUE (solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODELST_OK);
  return MURGE_SUCCESS;
}

INTS MURGE_SetOrdering(INTS id, INTS * permutation) {
  INTS i;
  print_debug(DBG_MURGE, ">> MURGE_SetOrdering\n");
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Build graph before calling MURGE_SetOrdering");
    return MURGE_ERR_ORDER;
  }
  if (solvers[id]->l2g == NULL) {
    errorPrint("Local to global array is not set");
    return MURGE_ERR_PARAMETER;
  }
  if (permutation == NULL) {
    errorPrint("NULL parameter");
    return MURGE_ERR_PARAMETER;
  }
  MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, PASTIX_INT);
  memset(solvers[id]->perm, 0, solvers[id]->n*sizeof(PASTIX_INT));
  for (i = 0; i < solvers[id]->n; i++)
    solvers[id]->perm[i] = permutation[solvers[id]->l2g[i]-1];

  solvers[id]->pastix_data->iparm[IPARM_ORDERING] = API_ORDER_PERSONAL;
  solvers[id]->pastix_data->iparm[IPARM_LEVEL_OF_FILL] = -1;

  print_debug(DBG_MURGE, "<< MURGE_SetOrdering\n");
  return MURGE_SUCCESS;
}

/*******************************************************************************
 * Group: Matrix assembly functions
 */

/* Function: MURGE_AssemblySetSequence
 *
 * Create a sequence of entries to build a matrix and store it for being reused.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   coefnbr - Number of entries.
 *   ROWs    - List of rows in the sequence.
 *   COLs    - List of columns in the sequence.
 *   op      - Operation to perform for coefficient which appear
 *             several tim (see <MURGE_ASSEMBLY_OP>).
 *   op2     - Operation to perform when a coefficient is set by
 *             two different processors (see <MURGE_ASSEMBLY_OP>).
 *   mode    - Indicates if user ensure he will respect solvers distribution
 *             (see <MURGE_ASSEMBLY_MODE>).
 *   nodes   - 0 entries are entered value by value,
 *             1 entries are entries node by node.
 *   id_seq  - Sequence ID.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *op*, *mode*, *sym*, or *coefnbr* are not valid.
 */
int MURGE_AssemblySetSequence (INTS id, INTL coefnbr, INTS * ROWs, INTS * COLs,
                               INTS op, INTS op2, INTS mode, INTS nodes,
                               INTS * id_seq) {
  murge_seq_t * sequence     = NULL;
  INTL         iter;
  ijv_t      **send_ijv      = NULL;
  INTL        *send_ijv_size = NULL;
  INTL        *send_nbr      = NULL;
  INTS         dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Graph has to be built before");
    return MURGE_ERR_ORDER;
  }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);


  MURGE_MEMALLOC(sequence, 1, murge_seq_t);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (solvers[id]->sequences == NULL) {
    solvers[id]->sequences = sequence;
  }
  else
    {
      murge_seq_t *last_sequence = solvers[id]->sequences;
      while(sequence->next != NULL) {
        last_sequence = sequence->next;
      }
      last_sequence->next = sequence;
    }

  sequence->next         = NULL;
  sequence->indexes      = NULL;
  sequence->recv_indexes = NULL;
  sequence->recv_nbr     = NULL;
  CHOOSE_FUNC(sequence->fusion_local_entries, op);
  CHOOSE_FUNC(sequence->fusion_dist_entries, op);
  sequence->mode  = mode;
  sequence->nodes = (nodes==MURGE_BOOLEAN_FALSE)?API_NO:API_YES;
  sequence->ID = solvers[id]->seq_ID++;
  *id_seq = sequence->ID;

  sequence->coefnbr = coefnbr;
  MURGE_MEMALLOC(sequence->indexes, coefnbr, INTL);
  if (sequence->mode == MURGE_ASSEMBLY_FOOL) {
    MURGE_MEMALLOC(send_nbr, solvers[id]->pastix_data->procnbr, INTL);
    MURGE_MEMALLOC(sequence->recv_nbr,
                   solvers[id]->pastix_data->procnbr,
                   INTL);
    MURGE_MEMALLOC(sequence->recv_indexes,
                   solvers[id]->pastix_data->procnbr,
                   INTL*);
    MURGE_MEMALLOC(send_ijv, solvers[id]->pastix_data->procnbr, ijv_t*);
    MURGE_MEMALLOC(send_ijv_size, solvers[id]->pastix_data->procnbr, INTL);

    for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++) {
      send_nbr[iter] = 0;
      send_ijv_size[iter] = 1+coefnbr/10;
      MURGE_MEMALLOC(send_ijv[iter], send_ijv_size[iter], ijv_t);
      sequence->recv_indexes[iter] = NULL;
    }
  }

  if (solvers[id]->colptr == NULL) {
    /* Need to build a CSC */
    INTS inc = (solvers[id]->pastix_data->iparm[IPARM_BASEVAL]==0)?1:0;

    vcsc_init(&solvers[id]->vcsc, solvers[id]->N, coefnbr/dof, 0, id);
    for (iter = 0; iter < coefnbr; iter++) {
      COEF * VAL = NULL;
      vcsc_add_node(solvers[id]->vcsc, COLs[iter]+inc, ROWs[iter]+inc, VAL, NULL, id);
    }
    vcsc_to_cscd(solvers[id]->vcsc, solvers[id]->pastix_data->pastix_comm,
                 &(solvers[id]->n), &(solvers[id]->colptr),
                 &(solvers[id]->rows), NULL,
                 &(solvers[id]->l2g), &(solvers[id]->g2l), NULL, 0, id);
  }

  coefnbr = sequence->coefnbr;

  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_ONLY_PROD)) {
    /* TODO: CHECK That all entries exists in CSC and if not insert it*/
  }

  for (iter = 0; iter < coefnbr; iter++) {
    INTL iter2;
    INTS inc = (solvers[id]->pastix_data->iparm[IPARM_BASEVAL]==0)?1:0;
    INTS col = COLs[iter]+inc; /* 1 based */
    INTS row = ROWs[iter]+inc; /* 1 based */
    INTS node_col; /* 0 based */
    INTS node_row; /* 0 based */
    INTS in_node_col; /* 0 based */
    INTS in_node_row; /* 0 based */
    INTS node_col_loc; /* 1 based */

    if (dof > 1 && sequence->nodes == API_NO) {
      node_col     = (col-1 - (col-1)%dof)/dof;
      in_node_col  = (col-1)%dof;
      node_row     = (row-1 - (row-1)%dof)/dof;
      in_node_row  = (row-1)%dof;
    }
    else
      {
        node_col     = col-1;
        in_node_col  = 0;
        node_row     = row-1;
        in_node_row  = 0;
      }

    node_col_loc = solvers[id]->g2l[node_col];
    if ( node_col_loc > 0 ) {
      node_col_loc--;
      /* Entry is local */
      for (iter2 = solvers[id]->colptr[node_col_loc]-1;
           iter2 < solvers[id]->colptr[node_col_loc+1]-1;
           iter2++) {
        if (solvers[id]->rows[iter2] == row)
          break;
      }
      if (solvers[id]->colptr[node_col_loc+1]-1 == iter2) {
        /* Entry not found in CSC */
        errorPrint("ROW (%ld:%ld) not found in COL (%d:%d) %d",
                   (long)row, node_row, (long)col, node_col,node_col_loc);
        return MURGE_ERR_PARAMETER;
      }
      else
        {
          if (iter2*dof*dof + in_node_col*dof+in_node_row >
              dof*dof*(solvers[id]->colptr[solvers[id]->n]-1)) {
            return MURGE_ERR_PARAMETER;
          }

          sequence->indexes[iter] = iter2*dof*dof + in_node_col*dof+in_node_row;
        }
    } else {
      /* Entry not local */
      if (sequence->mode == MURGE_ASSEMBLY_RESPECT) {
        errorPrint("COL (%d) is not local (row %d, owner %d)",
                   (long)col+1, (long)row+1, -node_col_loc);
        return MURGE_ERR_PARAMETER;
      } else {
        int owner = -node_col_loc;
        sequence->indexes[iter] = -(owner+1);

        /* => send buffer */
        if (send_nbr[owner] == send_ijv_size[owner]) {
          send_ijv_size[owner] = (INTL)(1.5*send_ijv_size[owner] + 1);
          MURGE_REALLOC(send_ijv[owner], send_ijv_size[owner], ijv_t);
        }
        send_ijv[owner][send_nbr[owner]].i = row;
        send_ijv[owner][send_nbr[owner]].j = col;
        send_nbr[owner]++;
      }
    }
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
    if (node_col == (MURGE_TRACE_COL-1)/dof && node_row == (MURGE_TRACE_ROW-1)/dof)
      fprintf(stdout, "%d, %d index => %d, iter %d node_col_loc %d send_nbr %d\n", node_col, node_row, sequence->indexes[iter], iter, node_col_loc, (node_col_loc>0)?-1:send_nbr[-node_col_loc] );
#endif
  }

  if (sequence->mode != MURGE_ASSEMBLY_RESPECT) {
    MPI_Request * requests;
    ijv_t       * recv_ijv;
    int           size;
    int           lastsize;
    PASTIX_INT           iter_coef;

    MPI_Alltoall(send_nbr,           1, MPI_INTL,
                 sequence->recv_nbr, 1, MPI_INTL,
                 solvers[id]->pastix_data->pastix_comm);

    MURGE_MEMALLOC(requests, solvers[id]->pastix_data->procnbr, MPI_Request);
    for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++) {
      if (send_nbr[iter] > 0)
        MPI_Isend(send_ijv[iter], send_nbr[iter]*sizeof(ijv_t), MPI_BYTE,
                  iter, TAG_IJV, solvers[id]->pastix_data->pastix_comm,
                  &(requests[iter]));
    }

    lastsize = sequence->recv_nbr[0];
    MURGE_MEMALLOC(recv_ijv, lastsize, ijv_t);
    for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++) {
      MPI_Status status;
      size = sequence->recv_nbr[iter];
      MURGE_MEMALLOC(sequence->recv_indexes[iter], size, INTL);
      if (lastsize < size) {
        MURGE_REALLOC(recv_ijv, size, ijv_t);
        lastsize = size;
      }
      if (size > 0)
        MPI_Recv(recv_ijv, size*sizeof(ijv_t), MPI_BYTE,
                 iter, TAG_IJV, solvers[id]->pastix_data->pastix_comm,
                 &status);

      for (iter_coef = 0; iter_coef < size; iter_coef++) {
        INTL iter2;
        INTS col = recv_ijv[iter_coef].j;
        INTS row = recv_ijv[iter_coef].i;
        INTS node_col;
        INTS node_row;
        INTS in_node_col;
        INTS in_node_row;
        INTS node_col_loc;


        if (dof > 1 && sequence->nodes == API_NO) {
          node_col = (col-1-(col-1)%dof)/dof;
          node_row = (row-1-(row-1)%dof)/dof;
          in_node_row = (row-1)%dof;
          in_node_col = (col-1)%col;
        }
        else
          {
            node_col = col-1;
            node_row = row-1;
            in_node_row = 0;
            in_node_col = 0;
          }
        node_col_loc = solvers[id]->g2l[node_col];

        if ( node_col_loc > 0 ) {
          /* Entry is local */
          for (iter2 = solvers[id]->colptr[node_col_loc-1]-1;
               iter2 < solvers[id]->colptr[node_col_loc]-1;
               iter2++) {
            if (solvers[id]->rows[iter2] == node_row+1)
              break;
          }
          if (solvers[id]->colptr[node_col_loc]-1 == iter2) {
            /* Entry not found in CSC */
            errorPrint("ROW (%ld) not found in COL (%d) %d",
                       (long)row, (long)col, node_col_loc);

            return MURGE_ERR_PARAMETER;
          }
          else
            {
              sequence->recv_indexes[iter][iter_coef] =
                iter2*dof*dof + in_node_col*dof+in_node_row;
            }
        }
        else
          {
            /* Entry not local */
            errorPrint("%s:%d COL (%d) is not local (row %d, owner %d)",
                       __FILE__, __LINE__, (long)col, row, -node_col_loc);
            return MURGE_ERR_PARAMETER;
          }
      }
    }

    MURGE_FREE(recv_ijv);

    for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++) {
      MPI_Status status;

      if (send_nbr[iter] > 0)
        MPI_Wait(&(requests[iter]), &status);
      MURGE_FREE(send_ijv[iter]);
    }
    MURGE_FREE(send_ijv);
    MURGE_FREE(send_ijv_size);
    MURGE_FREE(requests);
    MURGE_FREE(send_nbr);
  }

  return MURGE_SUCCESS;
}

static inline
int try_complete_reception(INTS    id,
                           murge_seq_t * sequence,
                           int     src,
                           int     send_nbr,
                           int     blocksize,
                           PASTIX_FLOAT * recv_data,
                           PASTIX_INT   * recv_cnt) {
  MPI_Status TCR_status;
  int        TCR_flag = 1;

  while (TCR_flag) {
    /* Do we have something to receive ? */
    PASTIX_INT        TCR_iter;
    MPI_Iprobe(src, TAG_VAL,
               solvers[id]->pastix_data->pastix_comm,
               &TCR_flag, &TCR_status);
    if (TCR_flag) {
      /* Receive and add it */
      int     TCR_src = TCR_status.MPI_SOURCE;
      int     TCR_tag = TCR_status.MPI_TAG;
      PASTIX_FLOAT * TCR_vals_ptr;
      MPI_Recv(recv_data, send_nbr*blocksize, COMM_FLOAT,
               TCR_src, TCR_tag,
               solvers[id]->pastix_data->pastix_comm, &TCR_status);
      TCR_vals_ptr=recv_data;
      for (TCR_iter = 0; TCR_iter < send_nbr; TCR_iter++) {
        PASTIX_INT     TCR_index =
          sequence->recv_indexes[TCR_src][recv_cnt[TCR_src]];
        PASTIX_FLOAT * TCR_node  = &(solvers[id]->values[TCR_index]);
        PASTIX_INT     TCR_iterdof;
        recv_cnt[TCR_src]++;
        for (TCR_iterdof = 0; TCR_iterdof < blocksize;
             TCR_iterdof++, TCR_node++, TCR_vals_ptr++)
          *TCR_node= sequence->fusion_local_entries(*TCR_node,
                                                    *TCR_vals_ptr);
      }
    }
  }
  return MURGE_SUCCESS;
}

/*
 * MURGE_AssemblyUseSequence
 *
 * Assembly the matrix using a stored sequence.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   id_seq  - Sequence ID.
 *   values  - Values to insert in the CSC.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *id_seq* or *values* are not valid.
 */
INTS MURGE_AssemblyUseSequence(INTS id, INTS id_seq, COEF * values) {
  murge_seq_t * sequence;
  INTL          iter;
  INTS          dof;
  PASTIX_FLOAT       * send_array = NULL;
  PASTIX_INT         * send_cnt   = NULL;
  MPI_Request * send_reqs  = NULL;
  MPI_Request * send_reqs2 = NULL;
  PASTIX_FLOAT       * recv_data  = NULL;
  PASTIX_INT         * recv_cnt   = NULL;
  int           send_nbr   = -1;
  int           blocksize;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  sequence = solvers[id]->sequences;
  while(sequence != NULL && sequence->ID != id_seq)
    sequence = sequence->next;

  if (sequence == NULL) {
    errorPrint("Sequence %d not found", id_seq);
    sequence = solvers[id]->sequences;
    return MURGE_ERR_PARAMETER;
  }

  if (values == NULL) {
    errorPrint("NULL value Pointer");
    return MURGE_ERR_PARAMETER;
  }

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Graph has to be built before");
    return MURGE_ERR_ORDER;
  }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (solvers[id]->values == NULL) {
    MURGE_MEMALLOC(solvers[id]->values,
                   (solvers[id]->colptr[solvers[id]->n]-1)*dof*dof,
                   PASTIX_FLOAT);
    memset(solvers[id]->values, 0,
           (solvers[id]->colptr[solvers[id]->n]-1)*dof*dof*sizeof(PASTIX_FLOAT));
  }

  if (sequence->nodes == API_YES) {
    blocksize = dof*dof;
  }
  else
    {
      blocksize = 1;
    }

  if (sequence->mode != MURGE_ASSEMBLY_RESPECT) {
    PASTIX_INT coefnbr, recv_coefnbr;
    coefnbr = sequence->coefnbr;
    MPI_Allreduce(&coefnbr, &recv_coefnbr, 1, COMM_INT, MPI_SUM,
                  solvers[id]->pastix_data->pastix_comm);
    send_nbr = recv_coefnbr/100 + 1;

    MURGE_MEMALLOC(send_reqs, solvers[id]->pastix_data->procnbr, MPI_Request);
    MURGE_MEMALLOC(send_reqs2, solvers[id]->pastix_data->procnbr, MPI_Request);
    MURGE_MEMALLOC(send_array,
                   blocksize*send_nbr*solvers[id]->pastix_data->procnbr, PASTIX_FLOAT);
    MURGE_MEMALLOC(recv_data, blocksize*send_nbr, PASTIX_FLOAT);
    MURGE_MEMALLOC(recv_cnt,  solvers[id]->pastix_data->procnbr, PASTIX_INT);
    MURGE_MEMALLOC(send_cnt, solvers[id]->pastix_data->procnbr, PASTIX_INT);
    for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++) {
      recv_cnt[iter] = 0;
      send_cnt[iter] = 0;
    }
  }

  for (iter = 0; iter < sequence->coefnbr; iter++) {
    INTL index;
    PASTIX_FLOAT*node;
    INTS iterdof;

    if (sequence->mode != MURGE_ASSEMBLY_RESPECT) {
      try_complete_reception(id, sequence,
                             MPI_ANY_SOURCE, send_nbr, blocksize,
                             recv_data, recv_cnt);
    }

    index = sequence->indexes[iter];

    if (index>=0) {
      if (index > dof*dof*(solvers[id]->colptr[solvers[id]->n]-1)) {
        return -1;
      }
      node = &(solvers[id]->values[index]);
      for (iterdof = 0; iterdof < blocksize; iterdof++, node++, values++)
        *node= sequence->fusion_local_entries(*node, *values);
    } else {
      if (sequence->mode == MURGE_ASSEMBLY_RESPECT) {
        errorPrint("Non local entry incompatible with MURGE_ASSEMBLY_RESPECT");
        return MURGE_ERR_PARAMETER;
      }

      if (send_cnt[-index-1] == send_nbr) {
        MPI_Status status;
        int flag = 0;
        while(!flag) {
          MPI_Test(&(send_reqs[-index-1]), &flag, &status);
          if (!flag) {
            try_complete_reception(id, sequence,
                                   MPI_ANY_SOURCE, send_nbr, blocksize,
                                   recv_data, recv_cnt);
          } else {
            send_cnt[-index-1] = 0;
          }
        }
      }

      /* Prepare to send, if send_buff is full send */
      node = &(send_array[blocksize*(send_cnt[-index-1]+send_nbr*(-index-1))]);
      for (iterdof = 0; iterdof < blocksize; iterdof++, node++, values++)
        *node= *values;

      send_cnt[-index-1]++;

      if (send_cnt[-index-1] == send_nbr) {
        MPI_Isend(&(send_array[send_nbr*(-index-1)*blocksize]),
                  send_cnt[-index-1]*blocksize, COMM_FLOAT,
                  -index-1, TAG_VAL,
                  solvers[id]->pastix_data->pastix_comm,
                  &(send_reqs[-index-1]));
      }
    }
  }

  if (sequence->mode != MURGE_ASSEMBLY_RESPECT) {
    PASTIX_INT * done, done_cnt;
    MURGE_MEMALLOC(done, solvers[id]->pastix_data->procnbr, PASTIX_INT);
    for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
      done[iter] = API_NO;
    done_cnt = 0;
    for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++) {
      /* Wait last sends */
      if (send_cnt[iter] == send_nbr) {
        MPI_Status status;
        MPI_Wait(&(send_reqs[iter]), &status);
        send_cnt[iter] = 0;
      }

      /* Send last entries */
      MPI_Isend(&(send_cnt[iter]), 1, COMM_INT,
                iter, TAG_SIZE,
                solvers[id]->pastix_data->pastix_comm,
                &(send_reqs[iter]));
      MPI_Isend(&(send_array[blocksize*send_nbr*(iter)]),
                send_cnt[iter]*blocksize, COMM_FLOAT,
                iter, TAG_VAL2,
                solvers[id]->pastix_data->pastix_comm,
                &(send_reqs2[iter]));
    }

    while (done_cnt < solvers[id]->pastix_data->procnbr) {
      for (iter =0; iter < solvers[id]->pastix_data->procnbr; iter++) {
        if (done[iter] == API_NO) {
          try_complete_reception(id, sequence,
                                 iter, send_nbr, blocksize,
                                 recv_data, recv_cnt);

          /* recv last count /entries */
          MPI_Status status;
          PASTIX_INT        cnt;
          PASTIX_INT        iter2;
          PASTIX_FLOAT     *myvals_ptr;
          int       flag;
          /* Receive and add it */
          MPI_Iprobe(iter, TAG_SIZE,
                     solvers[id]->pastix_data->pastix_comm, &flag, &status);
          if (flag) {
            MPI_Recv(&cnt, 1, COMM_INT,
                     iter, TAG_SIZE,
                     solvers[id]->pastix_data->pastix_comm, &status);
            MPI_Recv(recv_data, cnt*blocksize, COMM_FLOAT,
                     iter, TAG_VAL2,
                     solvers[id]->pastix_data->pastix_comm, &status);
            myvals_ptr = recv_data;

            for (iter2 = 0; iter2 < cnt; iter2++) {
              INTS iterdof;
              PASTIX_FLOAT*node;
              PASTIX_INT index;

              index = sequence->recv_indexes[iter][recv_cnt[iter]];
              recv_cnt[status.MPI_SOURCE]++;

              node  = &(solvers[id]->values[index]);
              for (iterdof = 0; iterdof < dof*dof;
                   iterdof++, node++, myvals_ptr++)
                *node= sequence->fusion_local_entries(*node, *myvals_ptr);
            }
            done[iter] = API_YES;
            done_cnt++;
          }
        }
      }
    }
    for (iter =0; iter < solvers[id]->pastix_data->procnbr; iter++) {
      MPI_Status status;
      MPI_Wait(&(send_reqs[iter]), &(status));
      MPI_Wait(&(send_reqs2[iter]), &(status));
    }

    MURGE_FREE(done);
    MURGE_FREE(send_reqs);
    MURGE_FREE(send_reqs2);
    MURGE_FREE(send_array);
    MURGE_FREE(recv_data);
    MURGE_FREE(recv_cnt);
    MURGE_FREE(send_cnt);
  }
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_DUMP_MATRIX;
  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_AssemblyDeleteSequence
 *
 * Destroy an assembly sequence
 *
 *   id      - Solver instance identification number.
 *   id_seq  - Sequence ID.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *id_seq* is not valid.
 */
INTS MURGE_AssemblyDeleteSequence(INTS id, INTS id_seq) {
  murge_seq_t * sequence;
  murge_seq_t * psequence;
  PASTIX_INT iter;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  psequence = NULL;
  sequence  = solvers[id]->sequences;
  while(sequence != NULL && sequence->ID != id_seq) {
    psequence = sequence;
    sequence  = sequence->next;
  }

  if (sequence == NULL) {
    errorPrint("Sequence %d not found", id_seq);
    return MURGE_ERR_PARAMETER;
  }

  if (psequence != NULL) {
    psequence->next = sequence->next;
  }
  else
    {
      solvers[id]->sequences = sequence->next;
    }
  MURGE_FREE(sequence->indexes);
  MURGE_FREE(sequence->recv_nbr);
  if (NULL != sequence->recv_indexes) {
    for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
      MURGE_FREE(sequence->recv_indexes[iter]);
    MURGE_FREE(sequence->recv_indexes);
  }
  MURGE_FREE(sequence);
  return MURGE_SUCCESS;
}
/*
 * Function: MURGE_AssemblyBegin
 *
 * Check that preprocessing has been performed, if not performs it.
 *
 * Allocate ijv structure which will be used to store I,J,v[dof*dof].
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   N       - Global number of nodes in the graph.
 *   coefnbr - Number of coeficients the calling MPI node will add
 *             (in term of coefficients, not node).
 *   op      - Operation to perform for coefficient which appear
 *             several time (see <MURGE_ASSEMBLY_OP>).
 *   op2     - Operation to perform when a coefficient is set by
 *             two different processors (see <MURGE_ASSEMBLY_OP>).
 *   mode    - Indicates if user ensure he will respect solvers distribution
 *             (see <MURGE_ASSEMBLY_MODE>).
 *   sym     - Indicates if user will give coefficient in a symmetric way
 *             (ie: only triangullar part) or not.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *op*, *mode*, *sym*, or *coefnbr* are not valid.
 */
INTS MURGE_AssemblyBegin(INTS id, INTS N, INTL coefnbr, INTS op,
                         INTS op2, INTS mode, INTS sym) {
  int dof;

  print_debug(DBG_MURGE, ">> MURGE_AssemblyBegin\n");
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];


  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK)) {
    CHECK_PREPROCESSING(id);

    CHECK_L2G(id);
  }
  else {
    if (mode == MURGE_ASSEMBLY_RESPECT) {
      errorPrint("Can't use MURGE_ASSEMBLY_RESPECT if you didn't build the graph before\n");
    }
  }

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#endif
  if (coefnbr < 0) {
    solvers[id]->dynamic = API_YES;
    coefnbr = -coefnbr;
  }
  else {
#ifdef MURGE_INSERT_DIRECTLY
    if (solvers[id]->colptr != NULL) {
      solvers[id]->dynamic = API_YES;
      coefnbr = 1+coefnbr/1000;
    }
    else
#endif /* MURGE_INSERT_DIRECTLY */
      {
        solvers[id]->dynamic = API_NO;
      }
  }
#ifdef MURGE_INSERT_DIRECTLY
  /* allocate values */
  if (solvers[id]->colptr != NULL && solvers[id]->values == NULL) {
    MURGE_MEMALLOC(solvers[id]->values, dof*dof*(solvers[id]->colptr[solvers[id]->n]-1), PASTIX_FLOAT);
    memset(solvers[id]->values, 0, (dof*dof*(solvers[id]->colptr[solvers[id]->n]-1))*sizeof(PASTIX_FLOAT));
  }
#endif /* MURGE_INSERT_DIRECTLY */
  solvers[id]->N           = N;
  solvers[id]->coefnbr     = coefnbr;
  solvers[id]->nodenbr     = coefnbr/(dof*dof);

  vcsc_init(&solvers[id]->vcsc, N, coefnbr/dof*dof, dof, id);

  solvers[id]->edgenbr  = coefnbr;
  solvers[id]->cnt      = 0;
  solvers[id]->cnt_zero = 0;
  solvers[id]->cnt_node = 0;
  solvers[id]->mode     = mode;
  solvers[id]->op       = op;
  solvers[id]->op2      = op2;
  solvers[id]->sym      = sym;


#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif

  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_FACTO_OK);

  MURGE_STATE_TRUE(solvers[id]->state,  MURGE_MATR_BUILD);
  MURGE_MEMORY_USAGE_PRINT("MURGE_AssemblyBegin");
  return MURGE_SUCCESS;
}

/*
 * Function:MURGE_AssemblySetValue
 *
 * Check that we are in an assembly section.
 *
 * Check that the number of coefficient entered will not
 * overpass the number of coefficient waited.
 *
 * Adds ROW, COL and value into the solvers[id]->tmpijv structure.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   ROW     - Global row number of the coefficient.
 *   COL     - Global column number of the coefficient.
 *   value   - value of the coefficient.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If we are not in an assembly section.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 */
INTS MURGE_AssemblySetValue     (INTS id, INTS ROW, INTS COL, COEF value) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_MATR_BUILD))) {
    errorPrint("Need to call MURGE_AssemblyBegin first");
    return MURGE_ERR_ORDER;
  }
  return MURGE_AssemblySetValue_(id, ROW, COL, value);
}

static inline
INTS MURGE_AssemblySetValue_    (INTS id, INTS ROW, INTS COL, COEF value) {
  int dof;


  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];



  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
    COL += 1;
    ROW += 1;
  }

  if (solvers[id]->dropmask != NULL && ROW != COL) {
    if (solvers[id]->dropmask[(COL-1)/dof] && solvers[id]->dropmask[(ROW-1)/dof]) {
      solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }
  }

  if (solvers[id]->droprows != NULL && ROW != COL) {
    if (solvers[id]->droprows[(ROW-1)/dof]) {
      solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }
  }

  if (solvers[id]->dropcols != NULL && ROW != COL) {
    if (solvers[id]->dropcols[(COL-1)/dof]) {
      solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }
  }

  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_ONLY_PROD) && value == 0.0) {
      solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }

  if ((COL < 1) || (ROW < 1) ||
      (ROW > (((solvers[id]->N)*dof))) ||
      (COL > (((solvers[id]->N)*dof)))) {
    errorPrint("COL (%ld) or ROW (%ld) is out of range [1-%ld]",
               (long)COL, (long)ROW, (long)(solvers[id]->N*dof));
    return MURGE_ERR_PARAMETER;
  }

#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
  if (COL == MURGE_TRACE_COL && MURGE_TRACE_ROW == ROW)
    fprintf(stdout, "Setting A(%d,%d) <- %g\n",
            MURGE_TRACE_ROW, MURGE_TRACE_COL, value);
#endif

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&(solvers[id]->mutex_tmpmatrix));
#endif

  if ((solvers[id]->dynamic == API_NO) &&
      (solvers[id]->cnt + solvers[id]->cnt_zero +
       (solvers[id]->cnt_node )*dof*dof + 1 > solvers[id]->coefnbr)) {
    errorPrint("Too many coef added in matrix building session (%ld > %ld)",
               (long)(solvers[id]->cnt + solvers[id]->cnt_zero +
                      (solvers[id]->cnt_node )*dof*dof + 1),
               (long)solvers[id]->coefnbr);
#ifdef MURGE_THREADSAFE
    pthread_mutex_unlock(&(solvers[id]->mutex_tmpmatrix));
#endif
    return MURGE_ERR_ORDER;
  }

  {
    PASTIX_INT node_col_glob;
    PASTIX_INT node_col_loc;


    node_col_glob = (COL-1 - (COL-1)%dof)/dof;

    if (solvers[id]->g2l != NULL) {
      node_col_loc = solvers[id]->g2l[node_col_glob];
    }
    else {
      /* If we didn't entered the graph we cannot decide if a column is local or not */
      node_col_loc = -1;
    }
#ifdef MURGE_INSERT_DIRECTLY
    if ( solvers[id]->colptr != NULL && node_col_loc > 0 ) {
      PASTIX_FLOAT (*func)(PASTIX_FLOAT , PASTIX_FLOAT);
      PASTIX_INT node_row_glob;
      PASTIX_INT node_idx;
      int coef_idx;

      CHOOSE_FUNC(func, solvers[id]->op);

      node_row_glob = (ROW-1 - (ROW-1)%dof)/dof;

      node_col_loc--;
      /* The column is local we add it into the local CSC */
      for (node_idx = solvers[id]->colptr[node_col_loc]-1;
           node_idx < solvers[id]->colptr[node_col_loc+1]-1;
           node_idx++)
        if (solvers[id]->rows[node_idx]-1 == node_row_glob) break;

      if (node_idx == solvers[id]->colptr[node_col_loc+1]-1) {
        if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_ONLY_PROD)) {
          /* we will add it later */
          vcsc_add(solvers[id]->vcsc, COL, ROW, value, func, id);
          solvers[id]->cnt++;
        }
        else
          {
            errorPrint("ROW (%ld) not found in COL (%d)",
                       (long)ROW, (long)COL);
#  ifdef MURGE_THREADSAFE
            pthread_mutex_unlock(&(solvers[id]->mutex_tmpmatrix));
#  endif
            return MURGE_ERR_PARAMETER;
          }
      }
      else
        {
          /* we found the correct node */
          coef_idx = node_idx*dof*dof + ((COL-1)%dof)*dof + (ROW-1)%dof;
          solvers[id]->values[coef_idx] =
            func(solvers[id]->values[coef_idx], value);
        }
    }
    else
#endif /* MURGE_INSERT_DIRECTLY */
      {
        /* The column has to be sent to the correct CPU */
        PASTIX_FLOAT (*func)(PASTIX_FLOAT , PASTIX_FLOAT);
        if (solvers[id]->g2l != NULL) {
          if (node_col_loc <= 0) {
            if (solvers[id]->mode == MURGE_ASSEMBLY_RESPECT) {
              errorPrint("Column %d is not local", COL);
              return MURGE_ERR_PARAMETER;
            }
          }
        }

        CHOOSE_FUNC(func, solvers[id]->op);
        vcsc_add(solvers[id]->vcsc, COL, ROW, value, func, id);

        solvers[id]->cnt++;
      }
  }
#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&(solvers[id]->mutex_tmpmatrix));
#endif
  return MURGE_SUCCESS;
}


/*
 * Function:MURGE_AssemblySetNodeValues
 *
 * Check that we are in an assembly section.
 *
 * Check that the number of coefficient entered will not
 * overpass the number of coefficient waited.
 *
 * Adds ROW, COL and value into the solvers[id]->tmpijv structure.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   ROW     - Global row number of the coefficient.
 *   COL     - Global column number of the coefficient.
 *   values  - value of the coefficient.(dof^2 matrix)
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If we are not in an assembly section.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 */
INTS MURGE_AssemblySetNodeValues (INTS id, INTS ROW, INTS COL, COEF *values) {
  return MURGE_AssemblySetNodeValues_(id, ROW, COL, values);
}

static inline
INTS MURGE_AssemblySetNodeValues_ (INTS id, INTS ROW, INTS COL, COEF *values) {
  int dof;

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (dof == 1)
    return MURGE_AssemblySetValue_(id, ROW, COL, *values);

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_MATR_BUILD))) {
    errorPrint("Need to call MURGE_AssemblyBegin first");
    return MURGE_ERR_ORDER;
  }
#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#endif
  if (solvers[id]->dynamic == API_NO &&
      ( (solvers[id]->cnt_node+1)*dof*dof +
        solvers[id]->cnt + solvers[id]->cnt_zero > solvers[id]->coefnbr)) {
    errorPrint("Too many coef added in matrix building session (%ld > %ld)",
               (long)((solvers[id]->cnt_node+1)*dof*dof + solvers[id]->cnt),
               (long)(solvers[id]->coefnbr));
#ifdef MURGE_THREADSAFE
    pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
    return MURGE_ERR_ORDER;
  }


  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
    ROW += 1;
    COL += 1;
  }


#if (defined MURGE_TRACE_ROW && defined MURGE_TRACE_COL)
  if ((COL-1)*dof <= MURGE_TRACE_COL-1 &&
      (COL)*dof > MURGE_TRACE_COL-1 &&
      (ROW-1)*dof <= MURGE_TRACE_ROW-1 &&
      ROW*dof > MURGE_TRACE_ROW-1) {

    fprintf(stdout, "Setting A(%d*%d+%d,%d*%d+%d) <- %g\n",
            ROW-1, dof, (MURGE_TRACE_ROW-1-(ROW-1)*dof)+1,
            COL-1, dof, (MURGE_TRACE_COL-1-(COL-1)*dof)+1,
            values[(MURGE_TRACE_COL-1-(COL-1)*dof)*dof +
                   (MURGE_TRACE_ROW-1 - (ROW-1)*dof)]);
  }
#endif

  if (solvers[id]->dropmask != NULL && ROW != COL)
    if (solvers[id]->dropmask[(COL-1)] && solvers[id]->dropmask[(ROW-1)]) {
      solvers[id]->cnt_zero+=dof*dof;
      return MURGE_SUCCESS;
    }
  if (solvers[id]->droprows != NULL && ROW != COL)
    if (solvers[id]->droprows[(ROW-1)]) {
      solvers[id]->cnt_zero+=dof*dof;
      return MURGE_SUCCESS;
    }
  if (solvers[id]->dropcols != NULL && ROW != COL)
    if (solvers[id]->dropcols[(COL-1)]) {
      solvers[id]->cnt_zero+=dof*dof;
      return MURGE_SUCCESS;
    }

  if ((COL < 1) || (ROW < 1) ||
      (ROW > ((solvers[id]->N))) ||
      (COL > ((solvers[id]->N)))) {
    errorPrint("COL (%ld) or ROW (%ld) is out of range [1-%ld]",
               (long)COL, (long)ROW, (long)(solvers[id]->N));
    return MURGE_ERR_PARAMETER;
  }
  {
    PASTIX_INT node_col_loc;
    if (solvers[id]->g2l != NULL) {
      node_col_loc = solvers[id]->g2l[COL-1];
    }
    else {
      /* If we didn't entered the graph we cannot decide if a column is local or not */
      node_col_loc = -1;
    }

#ifdef MURGE_INSERT_DIRECTLY
    if ( solvers[id]->colptr != NULL && node_col_loc > 0 ) {
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
      if ((COL-1)*dof <= MURGE_TRACE_COL-1 &&
          (COL)*dof > MURGE_TRACE_COL-1 &&
          (ROW-1)*dof <= MURGE_TRACE_ROW-1 &&
          ROW*dof > MURGE_TRACE_ROW-1) {
        fprintf(stdout, "The TRACED column is local\n");
      }
#endif
      PASTIX_FLOAT (*func)(PASTIX_FLOAT , PASTIX_FLOAT);
      PASTIX_INT node_idx;
      int coef_idx;

      CHOOSE_FUNC(func, solvers[id]->op);

      node_col_loc--;
      /* The column is local we add it into the local CSC */
      for (node_idx = solvers[id]->colptr[node_col_loc]-1;
           node_idx < solvers[id]->colptr[node_col_loc+1]-1;
           node_idx++)
        if (solvers[id]->rows[node_idx] == ROW) break;

      if (node_idx == solvers[id]->colptr[node_col_loc+1]-1) {
        if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_ONLY_PROD)) {
          vcsc_add_node(solvers[id]->vcsc, COL, ROW, values, func, id);
          solvers[id]->cnt_node++;
        }
        else
          {
            errorPrint("ROW (%ld) not found in COL (%d)",
                       (long)ROW, (long)COL);
#  ifdef MURGE_THREADSAFE
            pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#  endif
            return MURGE_ERR_PARAMETER;
          }
      }
      /* we found the correct node */
      for ( coef_idx = 0;
            coef_idx < dof*dof;
            coef_idx++)
        solvers[id]->values[node_idx*dof*dof + coef_idx] =
          func(solvers[id]->values[node_idx*dof*dof+coef_idx],
               values[coef_idx]);
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
      if ((COL-1)*dof <= MURGE_TRACE_COL-1 &&
          (COL)*dof > MURGE_TRACE_COL-1 &&
          (ROW-1)*dof <= MURGE_TRACE_ROW-1 &&
          ROW*dof > MURGE_TRACE_ROW-1) {
        fprintf(stdout, "Setting A(%d*%d+%d,%d*%d+%d) = %g\n",
                ROW-1, dof, (MURGE_TRACE_ROW-1-(ROW-1)*dof)+1,
                COL-1, dof, (MURGE_TRACE_COL-1-(COL-1)*dof)+1,
                solvers[id]->values[node_idx*dof*dof +
                                    (MURGE_TRACE_COL-1-(COL-1)*dof)*dof +
                                    (MURGE_TRACE_ROW-1 - (ROW-1)*dof)]);
      }
#endif
    } else
#endif /* MURGE_INSERT_DIRECTLY */
      {
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
        if ((COL-1)*dof <= MURGE_TRACE_COL-1 &&
            (COL)*dof > MURGE_TRACE_COL-1 &&
            (ROW-1)*dof <= MURGE_TRACE_ROW-1 &&
            ROW*dof > MURGE_TRACE_ROW-1) {
          fprintf(stdout, "The TRACED column is not local (node_col_loc %d)\n", node_col_loc);
        }
#endif

        PASTIX_FLOAT (*func)(PASTIX_FLOAT , PASTIX_FLOAT);
        if (solvers[id]->g2l != NULL) {
          if (node_col_loc <= 0) {
            if (solvers[id]->mode == MURGE_ASSEMBLY_RESPECT) {
              errorPrint("Column %d is not local", COL);
              return MURGE_ERR_PARAMETER;
            }
          }
        }
        CHOOSE_FUNC(func, solvers[id]->op);
        vcsc_add_node(solvers[id]->vcsc, COL, ROW, values, func, id);
        solvers[id]->cnt_node++;
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
        if ((COL-1)*dof <= MURGE_TRACE_COL-1 &&
            (COL)*dof > MURGE_TRACE_COL-1 &&
            (ROW-1)*dof <= MURGE_TRACE_ROW-1 &&
            ROW*dof > MURGE_TRACE_ROW-1) {
          int index = 0;
          while (index < solvers[id]->vcsc.colsizes[COL-1] &&
                 solvers[id]->vcsc.rows[COL-1][index] !=  ROW)
            index++;

          fprintf(stdout, "Setting vcsc[%d-1][%d*%d*%d + %d *%d + %d] = %g (vcsc.rows[%d] %d)\n",
                  COL, index, dof, dof, (MURGE_TRACE_COL-1-(COL-1)*dof), dof, 
                  (MURGE_TRACE_ROW-1 - (ROW-1)*dof),
                  solvers[id]->vcsc.values[COL-1][index*dof*dof + 
                                                  (MURGE_TRACE_COL-1-(COL-1)*dof)*dof +
                                                  (MURGE_TRACE_ROW-1 - (ROW-1)*dof)],
                  index,
                  solvers[id]->vcsc.rows[COL-1][index]);
        }
#endif
      }
  }
#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
  return MURGE_SUCCESS;
}

INTS MURGE_AssemblySetBlockValues(INTS id, INTS nROW, INTS *ROWlist,
                                  INTS nCOL, INTS *COLlist, COEF *values) {
  return MURGE_AssemblySetBlockValues_(id, nROW, ROWlist, nCOL, COLlist, values);
}
static inline
INTS MURGE_AssemblySetBlockValues_(INTS id, INTS nROW, INTS *ROWlist,
                                   INTS nCOL, INTS *COLlist, COEF *values) {
  PASTIX_INT  iter;
  PASTIX_INT  iter2;
  INTS ret;
  PASTIX_INT  iterv;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_MATR_BUILD))) {
    errorPrint("Need to call MURGE_GraphBegin first");
    return MURGE_ERR_ORDER;
  }

  iterv = 0;
  for (iter = 0; iter < nCOL; iter ++) {
    for (iter2 = 0; iter2 < nROW; iter2++) {
      if (solvers[id]->pastix_data->iparm[IPARM_DOF_NBR] == 1) {
        ret = MURGE_AssemblySetValue_(id,
                                      ROWlist[iter2],
                                      COLlist[iter],
                                      values[iterv]);
        if (MURGE_SUCCESS != ret)
          return ret;
        iterv++;
      }
      else
        {
          if (solvers[id]->pastix_data->iparm[IPARM_DOF_NBR] < 1)
            return MURGE_ERR_PARAMETER;
        ret = MURGE_AssemblySetNodeValues_(id,
                                           ROWlist[iter2],
                                           COLlist[iter] ,
                                           &(values[iterv]));
        if (MURGE_SUCCESS != ret)
          return ret;

        iterv+=solvers[id]->pastix_data->iparm[IPARM_DOF_NBR]*
          solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
      }
    }
  }

  return MURGE_SUCCESS;
}

INTS MURGE_AssemblySetListOfBlockValues(INTS id, INTS nBlocks,
                                        INTS nROW, INTS *ROWlist,
                                        INTS nCOL, INTS *COLlist,
                                        COEF *values) {
  INTS i, ierr;
  int dof;

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  for (i = 0; i < nBlocks; i++) {
    ierr = MURGE_AssemblySetBlockValues_(id,
                                         nROW, ROWlist,
                                         nCOL, COLlist,
                                         values);
    if (ierr != MURGE_SUCCESS)
      return ierr;
    ROWlist+=nROW;
    COLlist+=nCOL;
    values+=nROW*nCOL*dof*dof;
  }
  return MURGE_SUCCESS;
}
/*
 Function: MURGE_AssemblyEnd

 We have on each proc a part of the matrix in
 two structure, one containing nodes to add
 to the CSCd the other containing simple values.

 We send all data to his owner:
 - We sort our data structures (IJV structures)
 using the "owner" attribute.
 - We send non local data to other processors.

 We merge all data in the node structure.
 - We receive Data and merge node structure with simple
 values one.
 - We look for each coef in node structure, if present we modify the node, if
 not, we search in the CSCd and directly modify it. Else we construct
 a new node and add it.

 We Add this structure to the local CSCd.

 */

INTS MURGE_AssemblyEnd(INTS id) {
  int          dof;
  PASTIX_FLOAT (*func)(PASTIX_FLOAT , PASTIX_FLOAT);
#ifdef CENTRALISED
  PASTIX_INT         *total_nodelist;
#endif
#ifdef MURGE_TIME
  Clock       clock;
#endif
  CLOCK_INIT;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  print_debug(DBG_MURGE, ">> MURGE_AssemblyEnd\n");

  /* Check that we are in a matrix building session */
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_MATR_BUILD))) {
    errorPrint("Need to call MURGE_AssemblyBegin first");
    return MURGE_ERR_ORDER;
  }

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  /* check that we entered the correct number of values */
  if (solvers[id]->dynamic == API_NO &&
      solvers[id]->cnt_node*dof*dof + solvers[id]->cnt + solvers[id]->cnt_zero != solvers[id]->edgenbr) {
    errorPrint("Wrong number of entries  (%ld != %ld) ",
               (long)solvers[id]->cnt_node*dof*dof + solvers[id]->cnt + solvers[id]->cnt_zero ,
               (long)solvers[id]->edgenbr);
    return MURGE_ERR_ORDER;
  }

  /* some information about skipped zeros entries */
  if (solvers[id]->cnt_zero != 0) {
    if (solvers[id]->pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
      fprintf(stdout,
              "%ld (%.2g %%) zero entries were skipped on proc %ld\n",
              (long)solvers[id]->cnt_zero,
              (double)((double)100.0*((double)solvers[id]->cnt_zero)/
                       ((double)(solvers[id]->cnt_node*dof*dof +
                                 solvers[id]->cnt + solvers[id]->cnt_zero))),
              (long)solvers[id]->pastix_data->procnum);
    }
    else
      {
        if (solvers[id]->pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NO) {
          PASTIX_INT nz_glob;
          PASTIX_INT zeros_glob;
          PASTIX_INT nz = (solvers[id]->colptr[solvers[id]->n]-1)*dof*dof;
          MPI_Reduce( &(solvers[id]->cnt_zero), &zeros_glob,
                      1, COMM_INT,
                      MPI_SUM, 0, solvers[id]->pastix_data->pastix_comm);
          MPI_Reduce( &nz, &nz_glob,
                      1, COMM_INT,
                      MPI_SUM, 0, solvers[id]->pastix_data->pastix_comm);
          if (solvers[id]->pastix_data->procnum == 0) {
            fprintf(stdout,
                    "%ld zero entries were skipped"
                    " (from %ld (%.3g%%))\n",
                    (long int)zeros_glob,
                    (long int)nz_glob,
                    100.0*((double)(zeros_glob)/
                           ((double)(nz_glob))));
          }

        }
      }
  }

#ifdef CENTRALISED
  MURGE_MEMALLOC(total_nodelist, solvers[id]->N, PASTIX_INT);
  for (iter = 0; iter < solvers[id]->N; iter++)
    total_nodelist[iter] = iter+1;
#endif

  CHOOSE_FUNC(func, solvers[id]->op);

  dof     = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];


  /* If the CSCd is existing it will fill it with new values, else it will create it.
   *
   * NB : VCSC has to fit into CSCd
   */
  vcsc_to_cscd(solvers[id]->vcsc, solvers[id]->pastix_data->pastix_comm,
               &(solvers[id]->n), &(solvers[id]->colptr),
               &(solvers[id]->rows), &(solvers[id]->values),
               &(solvers[id]->l2g), &(solvers[id]->g2l), func,
               MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_ONLY_PROD), id);
  vcsc_destroy(solvers[id]->vcsc, id);

#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_VALUES_OK);

  MURGE_STATE_FALSE(solvers[id]->state, MURGE_MATR_BUILD);
  MURGE_DUMP_MATRIX;
  CLOCK_PRINT("MURGE_AssemblyEnd");
  return MURGE_SUCCESS;
}

INTS MURGE_MatrixReset(INTS id){
  PASTIX_INT nbcoef;
  PASTIX_INT dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_ONLY_PROD) &&
      solvers[id]->colptr) {
    MURGE_FREE(solvers[id]->colptr);
    MURGE_FREE(solvers[id]->rows);
    MURGE_FREE(solvers[id]->values);
  }
  else {
    if (solvers[id]->values != NULL) {
      nbcoef = solvers[id]->colptr[solvers[id]->n]-1;
      dof    = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

      memset(solvers[id]->values, 0, nbcoef*dof*dof*sizeof(PASTIX_FLOAT));
    }
    MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  }
  return MURGE_SUCCESS;
}

INTS MURGE_MatrixGlobalCSR(INTS id, INTS N, INTL *rowptr, INTS *COLS,
                           COEF *values, INTS root, INTS op, INTS sym) {
  PASTIX_INT  dof;
  PASTIX_INT  coefnbr;
  PASTIX_INT  iter;
  PASTIX_INT  iter2;
  INTS ret;
  int baseval;

  CHECK_SOLVER_ID(id);
  if (solvers[id]->pastix_data->procnum == 0) {
    errorPrintW("MURGE_MatrixGlobalCSR is not optimal with PaStiX, try using MURGE_MatrixGlobalCSC instead");
  }
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Graph has to be built before");
    return MURGE_ERR_ORDER;
  }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);
  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  dof     = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  /* Si on a une matrice symetrique autant faire un MURGE_MatrixGlobalCSC */
  if (solvers[id]->pastix_data->iparm[IPARM_SYM] == API_SYM_YES)
    return  MURGE_MatrixGlobalCSC(id, N, rowptr, COLS, values, root, op, sym);

  coefnbr = rowptr[N] - baseval;

  if (solvers[id]->pastix_data->procnum == root ||
      (root == -1 && solvers[id]->pastix_data->procnum == 0)) {
    ret = MURGE_AssemblyBegin(id, N, coefnbr,  op, op,MURGE_ASSEMBLY_FOOL,
                              solvers[id]->pastix_data->iparm[IPARM_SYM]);
    if (MURGE_SUCCESS != ret)
      return ret;
    for (iter = 0; iter < N; iter++) {
      for (iter2 = rowptr[iter]; iter2 < rowptr[iter+1]; iter2++) {
        if (dof == 1) {

          ret = MURGE_AssemblySetValue_(id,
                                        iter+baseval,
                                        COLS[iter2-baseval],
                                        values[iter2-baseval]);
          if (MURGE_SUCCESS != ret)
            return ret;
        }
        else
          {
            ret = MURGE_AssemblySetNodeValues_(id,
                                               iter+baseval,
                                               COLS[iter2-baseval],
                                               &(values[(iter2-baseval)*
                                                        dof*dof]));
            if (MURGE_SUCCESS != ret)
              return ret;
          }
      }
    }
  }
  else
    {
      ret = MURGE_AssemblyBegin(id, N, 0, op, op,
                                MURGE_ASSEMBLY_FOOL,
                                solvers[id]->pastix_data->iparm[IPARM_SYM]);
      if (MURGE_SUCCESS != ret)
        return ret;
    }

  if (MURGE_SUCCESS != (ret = MURGE_AssemblyEnd(id)))
    return ret;

  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  return MURGE_SUCCESS;
}
/*
 Function: MURGE_MatrixGlobalCSC

 Give a CSC on one processor to PaStiX.
 */
INTS MURGE_MatrixGlobalCSC(INTS id, INTS N, INTL *COLPTR, INTS *ROWS,
                           COEF *values, INTS root, INTS op, INTS sym) {
  PASTIX_INT          *l2g = NULL;
  PASTIX_INT           procnum;
  PASTIX_INT           localn;
  PASTIX_INT          *tmpcolptr;
  PASTIX_INT          *tmprows;
  PASTIX_FLOAT        *tmpvalues;
  PASTIX_INT          *tmpcolptr2;
  PASTIX_INT          *tmprows2;
  PASTIX_INT           tmpn;
  MPI_Status    status;
  PASTIX_INT           iter;
  int           dof;
  PASTIX_FLOAT (*func)(PASTIX_FLOAT , PASTIX_FLOAT);

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);


  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Graph has to be built before");
    return MURGE_ERR_ORDER;
  }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);

  CHOOSE_FUNC(func, op);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  /* Si tout le monde est racine */
  if (root == -1 || root == solvers[id]->pastix_data->procnum) {

    if (sizeof(INTS) != sizeof(PASTIX_INT)) {
      MURGE_MEMALLOC(tmprows2, (COLPTR[N]-1), PASTIX_INT);
      for (iter = 0; iter <  COLPTR[N]-1; iter++) {
        tmprows2[iter] = (PASTIX_INT)ROWS[iter];
      }
    }
    else
      {
        tmprows2 = (PASTIX_INT*)ROWS;
      }

    if (sizeof(INTL) != sizeof(PASTIX_INT)) {
      MURGE_MEMALLOC(tmpcolptr2, N+1, PASTIX_INT);
      for (iter = 0; iter <  N+1; iter++) {
        tmpcolptr2[iter] = (PASTIX_INT)COLPTR[iter];
      }
    }
    else
      {
        tmpcolptr2 = (PASTIX_INT*)COLPTR;
      }


    if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
      tmprows2 --;
      values   --;
    }

    MURGE_MEMALLOC(l2g, N, PASTIX_INT);

    for (iter = 0; iter < N; iter++)
      l2g[iter] = iter+1;

    /* colptr must be allocated for cscd_addlocal_int */
    if (NULL == solvers[id]->colptr) {
      MURGE_MEMALLOC(solvers[id]->colptr, solvers[id]->n+1, PASTIX_INT);
      for (iter = 0; iter < solvers[id]->n+1; iter++) {
        solvers[id]->colptr[iter] = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
      }
    }

    cscd_addlocal_int(solvers[id]->n,
                      solvers[id]->colptr,
                      solvers[id]->rows,
                      solvers[id]->values,
                      solvers[id]->l2g,
                      N, tmpcolptr2, tmprows2, (PASTIX_FLOAT*)values, l2g,
                      &tmpn, &tmpcolptr, &tmprows, &tmpvalues, func,
                      dof,  API_YES);
    MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpcolptr), char);
    MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprows), char);
    MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpvalues), char);
    if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
      tmprows2 ++;
      values   ++;
    }

    if (root == -1 && sizeof(INTS) != sizeof(PASTIX_INT)) {
      MURGE_FREE(tmprows2);
    }
    if (root == -1 && sizeof(INTL) != sizeof(PASTIX_INT)) {
      MURGE_FREE(tmpcolptr2);
    }
    MURGE_FREE(solvers[id]->colptr);
    MURGE_FREE(solvers[id]->rows);
    MURGE_FREE(solvers[id]->values);

    solvers[id]->colptr = tmpcolptr;
    solvers[id]->rows   = tmprows;
    solvers[id]->values = tmpvalues;

  }
  if (root != -1) {
    /* si on est le processeur racine */
    if (root == solvers[id]->pastix_data->procnum) {

      /* Pour chaque processeur, on calcule
       la CSCD a ajouter puis on l'envoi.
       */
      for (procnum = 0; procnum < solvers[id]->pastix_data->procnbr; procnum++) {
        if (procnum != solvers[id]->pastix_data->procnum) {
          MPI_Recv(&localn, 1,  COMM_INT, procnum, TAG_SIZE,
                   solvers[id]->pastix_data->pastix_comm, &status);
          MPI_Recv(l2g, localn, COMM_INT, procnum, TAG_L2G,
                   solvers[id]->pastix_data->pastix_comm, &status);

          MURGE_MEMALLOC(tmpcolptr, localn + 1, PASTIX_INT);

          for (iter = 0; iter < localn+1; iter++) {
            tmpcolptr[iter] = 1;
          }
          if (sizeof(INTS) != sizeof(PASTIX_INT)) {
            MURGE_MEMALLOC(tmprows2, (COLPTR[N]-1), PASTIX_INT);
            for (iter = 0; iter <  COLPTR[N]-1; iter++) {
              tmprows2[iter] = (PASTIX_INT)ROWS[iter];
            }
          }
          else
            {
              tmprows2 = (PASTIX_INT*)ROWS;
            }

          if (sizeof(INTL) != sizeof(PASTIX_INT)) {
            MURGE_MEMALLOC(tmpcolptr2, N+1, PASTIX_INT);
            for (iter = 0; iter <  N+1; iter++) {
              tmpcolptr2[iter] = (PASTIX_INT)COLPTR[iter];
            }
          }
          else
            {
              tmpcolptr2 = (PASTIX_INT*)ROWS;
            }

          if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
            tmprows2 --;
            values   --;
          }
          cscd_addlocal_int(localn,
                            tmpcolptr,
                            NULL,
                            NULL,
                            l2g,
                            N, tmpcolptr2, tmprows2, (PASTIX_FLOAT*)values, l2g,
                            &localn,
                            &tmpcolptr,
                            &tmprows,
                            &tmpvalues, func, dof, API_YES);
	  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpcolptr), char);
	  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprows), char);
	  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpvalues), char);

          if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
            tmprows2 ++;
            values   ++;
          }

          /* On envoi tmpcolptr, tmprows, tmpvalues */
          MPI_Send(tmpcolptr,
                   localn+1, COMM_INT, procnum,
                   TAG_COL, solvers[id]->pastix_data->pastix_comm);
          MURGE_FREE(tmpcolptr);

          MPI_Send(tmprows,
                   tmpcolptr[localn - 1], COMM_INT, procnum,
                   TAG_ROW, solvers[id]->pastix_data->pastix_comm);
          MURGE_FREE(tmprows);

          MPI_Send(tmpvalues,
                   tmpcolptr[localn - 1], COMM_FLOAT, procnum,
                   TAG_VAL, solvers[id]->pastix_data->pastix_comm);
          MURGE_FREE(tmpvalues);

          if (sizeof(INTS) != sizeof(PASTIX_INT)) {
            MURGE_FREE(tmprows2);
          }
          if (sizeof(INTL) != sizeof(PASTIX_INT)) {
            MURGE_FREE(tmpcolptr2);
          }
        }
        else
          {
            /* La CSCd local a déjà été traitée */
          }

      }
      if (sizeof(INTS) != sizeof(PASTIX_INT)) {
        MURGE_FREE(tmprows2);
      }
      if (sizeof(INTL) != sizeof(PASTIX_INT)) {
        MURGE_FREE(tmpcolptr2);
      }
    }
    else
      {
        /* Si on est pas la racine, on recoit de la racine la CSCd a ajouter
         et on l'ajoute
         */

        MPI_Send(&solvers[id]->n,
                 1, COMM_INT, root,
                 TAG_SIZE, solvers[id]->pastix_data->pastix_comm);
        localn = solvers[id]->n;
        MPI_Send(solvers[id]->l2g,
                 solvers[id]->n,
                 COMM_INT, root,
                 TAG_L2G, solvers[id]->pastix_data->pastix_comm);

        MPI_Recv(solvers[id]->colptr,
                 localn+1, COMM_INT, root,
                 TAG_COL, solvers[id]->pastix_data->pastix_comm, &status);

        MPI_Recv(solvers[id]->rows,
                 tmpcolptr[localn - 1], COMM_INT, root,
                 TAG_ROW, solvers[id]->pastix_data->pastix_comm, &status);

        MPI_Recv(solvers[id]->values,
                 tmpcolptr[localn - 1], COMM_FLOAT, root,
                 TAG_VAL, solvers[id]->pastix_data->pastix_comm, &status);


      }
  }

  MURGE_DUMP_MATRIX;


  MURGE_STATE_TRUE(solvers[id]->state, MURGE_VALUES_OK);

  return MURGE_SUCCESS;
}

/*
 Function: MURGE_MatrixGlobalIJV

 Add the given global Compress Sparse Column matrix to the matrix.

 Parameters:
 id      - Solver instance identification number.
 N       - Number of edges.
 NNZ     - Number of non zeros.
 ROWS    - Global row number array.
 COLS    - Global column number array.
 values  - values array.
 root    - Root processor for MPI communications.
 op      - Operation to perform if a coefficient appear twice
 (see <MURGE_ASSEMBLY_OP>).
 sym     - Indicates if user will give coefficient in a symmetric way
 (ie: only triangullar part) or not.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range,
 if *root*, *op*, *ROWS* or *COLS* are not valid.

 Fortran interface:
 >
 > SUBROUTINE MURGE_MATRIXGLOBALIJV(ID, N, NNZ, ROWS, COLS, VALUES, &
 >                                & ROOT, OP, SYM, IERROR)
 >   INTS,               INTENT(IN)  :: ID, ROOT, OP, SYM, N
 >   INTL,               INTENT(IN)  :: NNZ
 >   INTS, DIMENSION(0), INTENT(IN)  :: ROWS, COLS
 >   COEF, DIMENSION(0), INTENT(IN)  :: VALUES
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_MATRIXGLOBALIJV
 */
INTS MURGE_MatrixGlobalIJV(INTS id, INTS N, INTL NNZ, INTS *ROWS, INTS *COLS,
                           COEF *values, INTS root, INTS op, INTS sym) {
  PASTIX_INT        iter;                /* Iterators                            */
  INTS       ret;                 /* Return value                         */
  int        dof;                 /* Number of degree of freedom          */

  CHECK_SOLVER_ID(id);
  if (solvers[id]->pastix_data->procnum == 0)
    errorPrintW("MURGE_MatrixGlobalIJV is not optimal with PaStiX, try using MURGE_MatrixGlobalCSC instead");
  CHECK_SOLVER_PARAM(id);


  if ( MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD) ) {
    errorPrint("Do not call MURGE_GraphBegin before");
    return MURGE_ERR_ORDER;
  }

  if ( MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_MATR_BUILD) ) {
    errorPrint("Do not call MURGE_AssemblyBegin before");
    return MURGE_ERR_ORDER;
  }
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Graph has to be built before");
    return MURGE_ERR_ORDER;
  }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);
  dof     = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (solvers[id]->pastix_data->procnum == root ||
      (root == -1 && solvers[id]->pastix_data->procnum == 0)) {
    ret = MURGE_AssemblyBegin(id, N, NNZ, op, op, MURGE_ASSEMBLY_FOOL,
                              solvers[id]->pastix_data->iparm[IPARM_SYM]);
    if (MURGE_SUCCESS != ret)
      return ret;

    for (iter = 0; iter < NNZ; iter++) {
      if (dof == 1) {

        ret = MURGE_AssemblySetValue_(id,
                                      ROWS[iter],
                                      COLS[iter],
                                      values[iter]);
        if (MURGE_SUCCESS != ret)
          return ret;
      }
      else
        {
          ret = MURGE_AssemblySetNodeValues_(id,
                                             ROWS[iter],
                                             COLS[iter],
                                             &(values[iter*dof*dof]));
          if (MURGE_SUCCESS != ret)
            return ret;
        }
    }
  }
  else
    {
      ret = MURGE_AssemblyBegin(id, N, 0, op, op,
                                MURGE_ASSEMBLY_FOOL,
                                solvers[id]->pastix_data->iparm[IPARM_SYM]);
      if (MURGE_SUCCESS != ret)
        return ret;
    }

  if (MURGE_SUCCESS != (ret = MURGE_AssemblyEnd(id)))
    return ret;

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_VALUES_OK);
  return MURGE_SUCCESS;
}


/*******************************************************************************
 * Group: Filling the right-hand-side member
 */


INTS MURGE_SetGlobalRHS(INTS id, COEF *b, INTS root, INTS op) {
  PASTIX_INT        iter;
  PASTIX_INT        procnum;
  PASTIX_INT        localn;
  PASTIX_INT        lastn = 0;
  PASTIX_INT       *l2g   = NULL;
  PASTIX_FLOAT     *tmpb  = NULL;
  PASTIX_FLOAT    (*func)(PASTIX_FLOAT , PASTIX_FLOAT);
  MPI_Status status;
  int        dof;
  int        iterdof;
  COEF      *b_recv;
  int        allocated = API_NO;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);

  if (root < -1 || root >= solvers[id]->pastix_data->procnbr) {
    errorPrint("Invalid root value");
    return MURGE_ERR_PARAMETER;
  }
  CHECK_L2G(id);

  CHOOSE_FUNC(func, op);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (root != -1) {
    /* Broadcast and then run algorithm for root = -1 */
    /* TODO : fix mistake and remove BCast            */
    MURGE_MEMALLOC(b_recv, solvers[id]->N*dof, PASTIX_FLOAT);
    if (root ==  solvers[id]->pastix_data->procnum) {
      memcpy(b_recv, b, solvers[id]->N*dof*sizeof(PASTIX_FLOAT));
    }
    MPI_Bcast(b_recv, solvers[id]->N*dof, COMM_FLOAT, root,
              solvers[id]->pastix_data->pastix_comm);
    root = -1;
    allocated = API_YES;
  }
  else {
    b_recv = b;
  }

  if (root == -1) {
    /* Si tous le monde à la racine */
    if (NULL == solvers[id]->b) {
      MURGE_MEMALLOC(solvers[id]->b, solvers[id]->n*dof, PASTIX_FLOAT);
      memset(solvers[id]->b, 0, solvers[id]->n*dof*sizeof(PASTIX_FLOAT));
    }
    for (iter = 0; iter < solvers[id]->n; iter++) {
      for (iterdof = 0; iterdof < dof; iterdof++) {

        solvers[id]->b[iter*dof +iterdof] =
          func(solvers[id]->b[iter*dof +iterdof],
               b_recv[(solvers[id]->l2g[iter]-1)*dof+iterdof]);
      }
    }
  }
  else
    {
      /* Sinon, on recupère sur la racine tous les loc2globs
       On construit et on envoi tous les right-hand-side locaux
       */
      if (root == solvers[id]->pastix_data->procnum) {
        lastn = 0;
        for (procnum = 0;
             procnum < solvers[id]->pastix_data->procnbr;
             procnum++) {
          if (procnum != solvers[id]->pastix_data->procnum) {

            MPI_Recv(&localn, 1, COMM_INT, procnum, TAG_SIZE,
                     solvers[id]->pastix_data->pastix_comm, &status);

            if (lastn < localn) {
              if (NULL != l2g)
                MURGE_FREE(l2g);
              MURGE_MEMALLOC(l2g, localn, PASTIX_INT);
              if (tmpb != NULL)
                MURGE_FREE(tmpb);
              MURGE_MEMALLOC(tmpb, localn*dof, PASTIX_FLOAT);
              lastn = localn;
            }
            MPI_Recv(l2g, localn, COMM_INT, procnum, TAG_L2G,
                     solvers[id]->pastix_data->pastix_comm, &status);

            for (iter = 0; iter < localn; iter++) {
              for (iterdof = 0; iterdof < dof; iterdof++) {
                tmpb[iter*dof+iterdof]= b[(l2g[iter]-1)*dof+iterdof];
              }
            }
            MPI_Send(tmpb,
                     localn*dof, COMM_FLOAT, procnum,
                     TAG_VAL, solvers[id]->pastix_data->pastix_comm);
          }
        }
        if (NULL != tmpb)
          MURGE_FREE(tmpb);
        if (NULL != l2g)
          MURGE_FREE(l2g);

        /* Le processeur racine construit son bout de RHS */
        if (NULL == solvers[id]->b) {
          MURGE_MEMALLOC(solvers[id]->b, solvers[id]->n*dof, PASTIX_FLOAT);
          memset(solvers[id]->b, 0, solvers[id]->n*dof*sizeof(PASTIX_FLOAT));
        }

        for (iter = 0; iter < solvers[id]->n; iter++) {
          for (iterdof = 0; iterdof < dof; iterdof++) {
            solvers[id]->b[iter*dof+iterdof] =
              func(solvers[id]->b[iter*dof+iterdof],
                   b[(solvers[id]->l2g[iter]-1)*dof+iterdof]);
          }
        }
      }
      else
        {
          /* Sur les procs non racine on recoit simplement le RHS a ajouter */
          MPI_Send(&solvers[id]->n,
                   1, COMM_INT, root,
                   TAG_SIZE, solvers[id]->pastix_data->pastix_comm);
          MPI_Send(solvers[id]->l2g,
                   solvers[id]->n,
                   COMM_INT, root,
                   TAG_L2G, solvers[id]->pastix_data->pastix_comm);

          MURGE_MEMALLOC(tmpb, solvers[id]->n, PASTIX_FLOAT);
          if (NULL == solvers[id]->b) {
            solvers[id]->b = tmpb;
          }

          MPI_Recv(tmpb,
                   solvers[id]->n*dof, COMM_FLOAT, root,
                   TAG_VAL, solvers[id]->pastix_data->pastix_comm, &status);

          if (tmpb != solvers[id]->b) {
            for (iter = 0; iter < solvers[id]->n*dof; iter++) {
              solvers[id]->b[iter] =
                func(solvers[id]->b[iter],
                     tmpb[iter]);
            }
            MURGE_FREE(tmpb);
          }
        }
    }
  if (allocated == API_YES) MURGE_FREE(b_recv);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_RHS_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_SOLVE_DONE);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_REFINE_DONE);
  return MURGE_SUCCESS;
}

INTS MURGE_SetLocalRHS (INTS id, COEF *b, INTS op, INTS op2) {
  PASTIX_INT        iter;
  PASTIX_FLOAT    (*func)(PASTIX_FLOAT , PASTIX_FLOAT) = NULL;
  int        dof;
#ifdef CENTRALISED
  PASTIX_INT        nodenbr;
  PASTIX_INT       *intern_nodelist;
  PASTIX_FLOAT     *tmpb;
#endif

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CHECK_L2G(id);

  CHOOSE_FUNC(func, op);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

#ifdef CENTRALISED
  nodenbr = pastix_getLocalNodeNbr(&(solvers[id]->pastix_data));
  MURGE_MEMALLOC(intern_nodelist, nodenbr, PASTIX_INT);
  if (NO_ERR != ( pastix_getLocalNodeLst(&(solvers[id]->pastix_data),
                                         intern_nodelist)))
    return MURGE_ERR_SOLVER;

  if (NULL == solvers[id]->b) {
    MURGE_MEMALLOC(solvers[id]->b, solvers[id]->N*dof, PASTIX_FLOAT);
    memset(solvers[id]->b, 0, solvers[id]->N*dof*sizeof(PASTIX_FLOAT));
  }
  MURGE_MEMALLOC(tmpb, solvers[id]->N*dof /* solvers[id]->n*dof */, PASTIX_FLOAT);
  memset(tmpb, 0, solvers[id]->N*dof*sizeof(PASTIX_FLOAT));

  for (iter = 0; iter < nodenbr*dof; iter++) {
    tmpb[(intern_nodelist[(iter-iter%dof)/dof]-1)*dof+iter%dof] =
      func(tmpb[(intern_nodelist[(iter-iter%dof)/dof]-1)*dof+iter%dof],
           b[iter]);
  }
  MPI_Allreduce(tmpb, solvers[id]->b, solvers[id]->N*dof, COMM_FLOAT, MPI_SUM,
                solvers[id]->pastix_data->pastix_comm);

  MURGE_FREE(tmpb);
  MURGE_FREE(intern_nodelist);

#else /* CENTRALISED */
  if (NULL == solvers[id]->b) {
    MURGE_MEMALLOC(solvers[id]->b, solvers[id]->n*dof , PASTIX_FLOAT);
    memset(solvers[id]->b, 0, solvers[id]->n*dof*sizeof(PASTIX_FLOAT));
  }

  for (iter = 0; iter < solvers[id]->n*dof; iter++) {
    solvers[id]->b[iter] = func(solvers[id]->b[iter],
                                b[iter]);
  }
#endif /* CENTRALISED */

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_RHS_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_SOLVE_DONE);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_REFINE_DONE);
  return MURGE_SUCCESS;
}

static inline
INTS MURGE_SetRHS_local_ (INTS id, INTS n, INTS *coefsidx, COEF *b, INTS op) {
  PASTIX_INT             iter;
  PASTIX_FLOAT         (*func)(PASTIX_FLOAT , PASTIX_FLOAT) = NULL;
  PASTIX_INT             index;
  int                    baseval;
  int                    dof;
  int                    iterdof;
  pastix_data_t *        pastix_data = solvers[id]->pastix_data;

  baseval = pastix_data->iparm[IPARM_BASEVAL];
  dof = pastix_data->iparm[IPARM_DOF_NBR];
  CHOOSE_FUNC(func, op);

  for (iter = 0; iter < n; iter++) {
    index = coefsidx[iter]- baseval;
    index = solvers[id]->g2l[index] - 1;

    for (iterdof = 0; iterdof < dof; iterdof++) {
      solvers[id]->b[index*dof+iterdof] =
        func(solvers[id]->b[index*dof+iterdof],
             b[iter*dof+iterdof]);
    }
  }
  return MURGE_SUCCESS;
}

INTS MURGE_SetRHS      (INTS id, INTS n, INTS *coefsidx, COEF *b, INTS op,
                        INTS op2, INTS mode) {
  PASTIX_INT             iter;
  PASTIX_FLOAT         (*func)(PASTIX_FLOAT , PASTIX_FLOAT) = NULL;
  int             baseval;
  int             dof;
  int             iterdof;
  pastix_data_t * pastix_data = solvers[id]->pastix_data;
  INTS            procnbr     = pastix_data->procnbr;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CHECK_L2G(id);

  baseval = pastix_data->iparm[IPARM_BASEVAL];
  dof = pastix_data->iparm[IPARM_DOF_NBR];
  CHOOSE_FUNC(func, op);

  for (iter = 0; iter < n; iter++)
    if (coefsidx[iter]-baseval >= solvers[id]->N) {
      errorPrint("coefsidx[%ld] (%ld) is out of range [1-%ld]",
                 (long)iter, (long)coefsidx[iter]-baseval, (long)(solvers[id]->N));
      return MURGE_ERR_PARAMETER;

    }

  if (NULL == solvers[id]->b) {
    MURGE_MEMALLOC(solvers[id]->b, solvers[id]->n*dof , PASTIX_FLOAT);
    memset(solvers[id]->b, 0, solvers[id]->n*dof*sizeof(PASTIX_FLOAT));
  }

  if (mode == MURGE_ASSEMBLY_FOOL) {
    INTS            coefs_rcv_size;
    INTS         *  coefnbr;
    INTS         ** coefs_idx;
    COEF         ** coefs_vals;
    MPI_Request  *  request_cfnbr;
    MPI_Request  *  request_cfidx;
    MPI_Request  *  request_cfvals;

    MURGE_MEMALLOC(coefnbr, procnbr, INTS);
    for (iter = 0; iter <procnbr; iter++)
      coefnbr[iter] = 0;

    /* Count the entries to send to each processor */
    for (iter = 0; iter <n; iter++) {
      PASTIX_INT procnum;
      if (solvers[id]->g2l[coefsidx[iter]-baseval] > 0)
        procnum = pastix_data->procnum;
      else
        procnum = -solvers[id]->g2l[coefsidx[iter]-baseval];

      coefnbr[procnum]++;
    }
    MURGE_MEMALLOC(coefs_idx,  procnbr, INTS*);
    MURGE_MEMALLOC(coefs_vals, procnbr, COEF*);

    for (iter = 0; iter < procnbr; iter++) {
      MURGE_MEMALLOC(coefs_idx[iter],  coefnbr[iter],     INTS);
      MURGE_MEMALLOC(coefs_vals[iter], coefnbr[iter]*dof, COEF);
      coefnbr[iter] = 0;
    }
    /* Prepare the arrays to send to each processors */
    for (iter = 0; iter <n; iter++) {
      PASTIX_INT procnum;
      if (solvers[id]->g2l[coefsidx[iter]-baseval] > 0)
        procnum = pastix_data->procnum;
      else
        procnum = -solvers[id]->g2l[coefsidx[iter]-baseval];

      coefs_idx[procnum][coefnbr[procnum]] = coefsidx[iter];
      for (iterdof = 0; iterdof < dof; iterdof++) {
        coefs_vals[procnum][coefnbr[procnum]*dof+iterdof] =
          b[iter*dof+iterdof];
      }

      coefnbr[procnum]++;
    }

    MURGE_MEMALLOC(request_cfnbr,  procnbr, MPI_Request);
    MURGE_MEMALLOC(request_cfidx,  procnbr, MPI_Request);
    MURGE_MEMALLOC(request_cfvals, procnbr, MPI_Request);

    /* Send the data to the processors */
    for (iter = 0; iter < procnbr; iter++) {
      if (iter == solvers[id]->pastix_data->procnum) continue;
      MPI_Isend(&(coefnbr[iter]),    1,                 MPI_INTS,
                iter, TAG_SIZE, pastix_data->pastix_comm,
                &(request_cfnbr[iter]));
      if (coefnbr[iter] > 0) {
        MPI_Isend(coefs_idx[iter],     coefnbr[iter],     MPI_INTS,
                  iter, TAG_ROW, pastix_data->pastix_comm,
                  &(request_cfidx[iter]));
        MPI_Isend(coefs_vals[iter],    coefnbr[iter]*dof, MURGE_MPI_COEF,
                  iter, TAG_VAL, pastix_data->pastix_comm,
                  &(request_cfvals[iter]));
      }
    }

    /* receive the data and run MPI_SetRHS with MURGE_ASSEMBLY_RESPECT */
    {
      coefs_rcv_size = 0;
      INTS       * coefs_idx_rcv  = NULL;
      COEF       * coefs_vals_rcv = NULL;

      for (iter = 0; iter < procnbr; iter++) {
        INTS         coefnbr_rcv;
        MPI_Status   status;
        if (iter == solvers[id]->pastix_data->procnum) {
          MURGE_SetRHS_local_(id, coefnbr[iter], coefs_idx[iter], coefs_vals[iter],
                              op);
        } else {
          MPI_Recv(&coefnbr_rcv, 1, MPI_INTS, iter, TAG_SIZE,
                   pastix_data->pastix_comm, &status);

          if (coefnbr_rcv > 0) {
            if (coefnbr_rcv > coefs_rcv_size) {
              if (coefs_rcv_size != 0) {
                MURGE_FREE(coefs_idx_rcv);
                MURGE_FREE(coefs_vals_rcv);
              }
              MURGE_MEMALLOC(coefs_idx_rcv,  coefnbr_rcv,         INTS);
              MURGE_MEMALLOC(coefs_vals_rcv, coefnbr_rcv*dof,     COEF);
              coefs_rcv_size = coefnbr_rcv;
            }
            MPI_Recv(coefs_idx_rcv, coefnbr_rcv,     MPI_INTS,
                     iter, TAG_ROW, pastix_data->pastix_comm, &status);
            MPI_Recv(coefs_vals_rcv,coefnbr_rcv*dof, MURGE_MPI_COEF,
                     iter, TAG_VAL, pastix_data->pastix_comm, &status);

            MURGE_SetRHS_local_(id, coefnbr_rcv, coefs_idx_rcv, coefs_vals_rcv, op2);
          }
        }
        if (iter == procnbr-1 && coefs_rcv_size != 0) {
          MURGE_FREE(coefs_idx_rcv);
          MURGE_FREE(coefs_vals_rcv);
        }
      }
    }
    /* Now we clean it all */
    for (iter = 0; iter < procnbr; iter++) {
      MPI_Status status;
      if (iter == solvers[id]->pastix_data->procnum) continue;
      MPI_Wait(&(request_cfnbr[iter]), &status);
      if (coefnbr[iter] > 0) {
        MPI_Wait(&(request_cfidx[iter]), &status);
        MPI_Wait(&(request_cfvals[iter]), &status);
      }
    }
    MURGE_FREE(request_cfnbr);
    MURGE_FREE(request_cfidx);
    MURGE_FREE(request_cfvals);

    for (iter = 0; iter <procnbr; iter++) {
      MURGE_FREE(coefs_idx[iter]);
      MURGE_FREE(coefs_vals[iter]);
    }
    MURGE_FREE(coefs_idx);
    MURGE_FREE(coefs_vals);
    MURGE_FREE(coefnbr);
  }
  else
    {
      if (!MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_NODELST_OK)) {
        errorPrint("Assembly can be respected only if user asked for node list");
        return MURGE_ERR_ORDER;
      }
      MURGE_SetRHS_local_(id, n, coefsidx, b, op);
    }

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_RHS_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_SOLVE_DONE);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_REFINE_DONE);
  return MURGE_SUCCESS;
}

INTS MURGE_RHSReset(INTS id){
  PASTIX_INT iter, dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CHECK_L2G(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (NULL == solvers[id]->b)
    MURGE_MEMALLOC(solvers[id]->b, solvers[id]->n*dof, PASTIX_FLOAT);

  for (iter = 0; iter < solvers[id]->n*dof; iter++)
    solvers[id]->b[iter] =0;
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_RHS_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_SOLVE_DONE);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_REFINE_DONE);
  return MURGE_SUCCESS;
}


/*******************************************************************************
 * Group: Getting the solution
 */
static inline
INTS MURGE_GetGlobalSolution_(INTS id, COEF *x, INTS root) {
  int    dof;
#ifndef CENTRALISED
  PASTIX_FLOAT *tmpx = NULL;
  PASTIX_INT    iter;
  int    iterdof;
#endif

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if ((root == solvers[id]->pastix_data->procnum ||
       root == -1)
      && NULL == x)
    return MURGE_ERR_PARAMETER;
#ifdef CENTRALISED
  if ((root == solvers[id]->pastix_data->procnum ||
       root == -1))
    memcpy(x, solvers[id]->b, solvers[id]->N*dof*sizeof(PASTIX_FLOAT));
#else
  MURGE_MEMALLOC(tmpx, solvers[id]->N*dof, PASTIX_FLOAT);

  for (iter = 0; iter < solvers[id]->N*dof; iter ++)
    tmpx[iter] = 0;
  for (iter = 0; iter < solvers[id]->n; iter ++) {
    for (iterdof = 0; iterdof < dof; iterdof++) {
      tmpx[(solvers[id]->l2g[iter]-1)*dof+iterdof] =
        solvers[id]->b[iter*dof+iterdof];
    }
  }

  if (root == -1) {
    MPI_Allreduce(tmpx, x, solvers[id]->N*dof, COMM_FLOAT, COMM_SUM,
                  solvers[id]->pastix_data->pastix_comm);
  }
  else
    {
      MPI_Reduce(tmpx, x, solvers[id]->N*dof, COMM_FLOAT, COMM_SUM, root,
                 solvers[id]->pastix_data->pastix_comm);
    }

  MURGE_FREE(tmpx);
#endif
  return MURGE_SUCCESS;
}

INTS MURGE_GetGlobalSolution(INTS id, COEF *x, INTS root) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);
  CHECK_FACT(id);
  MURGE_GetGlobalSolution_(id, x, root);
  return MURGE_SUCCESS;
}

INTS MURGE_GetLocalSolution (INTS id, COEF *x) {
  PASTIX_INT    iter;
  int    dof;
#ifdef CENTRALISED
  PASTIX_INT    nodenbr;
  PASTIX_INT   *intern_nodelist;
  int    iterdof;
#endif
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);
  CHECK_FACT(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (NULL == x) {
    return MURGE_ERR_PARAMETER;
  }
#ifdef CENTRALISED
  nodenbr = pastix_getLocalNodeNbr(&(solvers[id]->pastix_data));
  MURGE_MEMALLOC(intern_nodelist, nodenbr, PASTIX_INT);
  if (NO_ERR != ( pastix_getLocalNodeLst(&(solvers[id]->pastix_data),
                                         intern_nodelist)))
    return MURGE_ERR_SOLVER;

  for (iter = 0; iter < nodenbr; iter ++) {
    for (iterdof = 0; iterdof < dof; iterdof++) {
      x[iter*dof+iterdof] = solvers[id]->b[(intern_nodelist[iter]-1)*
                                           dof+iterdof];
    }
  }
  MURGE_FREE(intern_nodelist);
#else
  for (iter = 0; iter < solvers[id]->n*dof; iter ++) {
    x[iter] = solvers[id]->b[iter];
  }
#endif
  return MURGE_SUCCESS;
}

INTS MURGE_GetSolution      (INTS id, INTS n, INTS *coefsidx, COEF *x,
                             INTS mode) {
  PASTIX_INT    iter;
  COEF  *tmpx = NULL;
  int    dof;
  int    iterdof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);
  CHECK_FACT(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  MURGE_MEMALLOC(tmpx, solvers[id]->N*dof, COEF);
  MURGE_GetGlobalSolution_(id, tmpx, -1);
  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  for (iter = 0; iter < n; iter ++) {
    for (iterdof = 0; iterdof < dof; iterdof++) {
      x[iter*dof+iterdof] = tmpx[(coefsidx[iter]-
                                  solvers[id]->pastix_data->iparm[IPARM_BASEVAL])*
                                 dof+iterdof];
    }
  }
  MURGE_FREE(tmpx);


  return MURGE_SUCCESS;
}

/*******************************************************************************
 * Group: Cleaning up this mess
 */

INTS MURGE_Clean(INTS id){
  PASTIX_INT    * iparm;
  double * dparm;
#ifdef MEMORY_USAGE
  int rank, verb;
  MPI_Comm comm;
  rank =  solvers[id]->pastix_data->procnum;
  comm =  solvers[id]->pastix_data->pastix_comm;
  verb =  iparm[IPARM_VERBOSE];
#endif
  iparm = solvers[id]->pastix_data->iparm;
  dparm = solvers[id]->pastix_data->dparm;
  if (NULL != solvers[id]->colptr)
    MURGE_FREE(solvers[id]->colptr);
  if (NULL != solvers[id]->rows)
    MURGE_FREE(solvers[id]->rows);
  if (NULL != solvers[id]->values)
    MURGE_FREE(solvers[id]->values);
  if (NULL != solvers[id]->b)
    MURGE_FREE(solvers[id]->b);
  if (NULL != solvers[id]->l2g)
    MURGE_FREE(solvers[id]->l2g);
  if (NULL != solvers[id]->g2l)
    MURGE_FREE(solvers[id]->g2l);
  if (NULL != solvers[id]->perm)
    MURGE_FREE(solvers[id]->perm);
#ifdef CENTRALISED
  if (NULL != solvers[id]->invp)
    MURGE_FREE(solvers[id]->invp);
#endif
  if (NULL != solvers[id]->dropmask)
    MURGE_FREE(solvers[id]->dropmask);
  if (NULL != solvers[id]->dropcols)
    MURGE_FREE(solvers[id]->dropcols);
  if (NULL != solvers[id]->droprows)
    MURGE_FREE(solvers[id]->droprows);

  while (NULL != solvers[id]->sequences) {
    MURGE_AssemblyDeleteSequence(id, solvers[id]->sequences->ID);
  }
  if (NULL != solvers[id]->pastix_data)
    pastix_task_clean(&(solvers[id]->pastix_data),
                      solvers[id]->pastix_data->pastix_comm);

#ifndef FORCE_NOSMP
  if (solvers[id]->threadnbr)
    stop_threads(id);
  pthread_mutex_destroy(&(solvers[id]->barrier.sync_lock));
  pthread_cond_destroy(&(solvers[id]->barrier.sync_cond));
#endif

#ifdef MURGE_THREADSAFE
  pthread_mutex_destroy(&solvers[id]->mutex_tmpmatrix);
#endif
  MURGE_FREE_EXT(iparm, IPARM_SIZE*sizeof(PASTIX_INT));
  MURGE_FREE_EXT(dparm, DPARM_SIZE*sizeof(double));
  MURGE_MEMORY_USAGE_PRINT2("MURGE_Clean", verb, rank, comm);
  free(solvers[id]); solvers[id] = NULL;
  /* Hack for clean restart */
  _MURGE_InitId(id);
  return MURGE_SUCCESS;
}


INTS MURGE_Finalize(){
  PASTIX_INT i;

  for (i=0; i< idnbr; i++) {
      MURGE_Clean(i);
    if (NULL != solvers[i]->pastix_data)
      pastix_task_clean(&(solvers[i]->pastix_data),
                        solvers[i]->pastix_data->pastix_comm);
    free(solvers[i]); solvers[i] = NULL;
  }

  free(solvers); solvers = NULL;

  return MURGE_SUCCESS;
}

INTS MURGE_GetInfoINT(INTS id,  INTS metric, INTL * value) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  PASTIX_INT murge_param[1];

  murge_param[MURGE_IINFO_NNZ  - 1024] =  IPARM_NNZEROS;

  if (metric >= 1024) {
    metric = murge_param[metric-1024];
  }

  if (!( metric < IPARM_SIZE )) {
    errorPrint("metric is too big");
    return MURGE_ERR_PARAMETER;
  }

  if (metric < 0) {
    errorPrint("metric is negative");
    return MURGE_ERR_PARAMETER;
  }

  /* TODO : Est-ce qu'on ajoute des tests sur les valeurs rentrées ???? */
  *value = solvers[id]->pastix_data->iparm[metric];

  return MURGE_SUCCESS;
}

INTS MURGE_GetInfoREAL(INTS id,  INTS metric, REAL * value) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  PASTIX_INT murge_param[1];

  murge_param[MURGE_RPARAM_EPSILON_ERROR  - 1024] =  DPARM_RELATIVE_ERROR;

  if (metric >= 1024) {
    metric = murge_param[metric-1024];
  }

  if (!( metric < IPARM_SIZE )) {
    errorPrint("metric is too big");
    return MURGE_ERR_PARAMETER;
  }

  if (metric < 0) {
    errorPrint("metric is negative");
    return MURGE_ERR_PARAMETER;
  }

  *value = (REAL)solvers[id]->pastix_data->dparm[metric];
  return MURGE_SUCCESS;
}
/*
 Function: MURGE_PrintError

 Print the error message corresponding to ierror
 Parameters:
 error_number  - Error identification number.

 Returns:
 MURGE_ERR_PARAMETER - If ierror does not match an error number
 MURGE_SUCCESS       - If function runned successfully.

 Fortran interface:
 >
 > SUBROUTINE MURGE_PRINTERROR(ERROR_NUMBER, IERROR)
 >   INTS, INTENT(IN)  :: IERROR
 >   INTS, INTENT(OUT) :: ERROR_NUMBER
 > END SUBROUTINE MURGE_PRINTERROR
 */
INTS MURGE_PrintError(INTS error_number) {
  return MURGE_SUCCESS;
}

/*
 Function: MURGE_ExitOnError

 Print the error message corresponding to ierror.
 If the ierr is not MURGE_SUCCESS then the program is stopped.

 Parameters:
 ierror         - Error identification number.

 Returns:
 MURGE_SUCCESS   - If function runned successfully,
 stop the program otherwise.

 Fortran interface:
 >
 > SUBROUTINE MURGE_EXITONERROR(ERROR_NUMBER, IERROR)
 >   INTS, INTENT(IN)  :: IERROR
 >   INTS, INTENT(OUT) :: ERROR_NUMBER
 > END SUBROUTINE MURGE_EXITONERROR
 */
INTS MURGE_ExitOnError(INTS error_number) {
  if  (error_number == MURGE_SUCCESS)
    return MURGE_SUCCESS;
  else
    exit(1);
}


/*
 Group: Scaling
 */

/*
 Function: MURGE_GetGlobalNorm

 Compute the global norm array following a norm rule.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 norm    - Array of size global column number*dof which will contain
 the norm values
 root    - Indicates which processor will have the norm array
 at the end of the call, -1 for all.
 rule    - Rule to follow to build norm array, see <MURGE_NORM_RULES>

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_GETGLOBALNORM(ID, NORM, ROOT, RULE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, ROOT, RULE
 >   REAL, DIMENSION(0), INTENT(OUT) :: NORM
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_GETGLOBALNORM
 */
INTS MURGE_GetGlobalNorm(INTS id, REAL *norm, INTS root, INTS rule) {
  INTS itercol;       /* Each column*/
  PASTIX_INT  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS scal_idx;    /* Column number of the value */
  PASTIX_INT  value_idx;     /* Index of the value */
  REAL*local_norm = NULL;
  INTS dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  print_debug(DBG_MURGE, ">> MURGE_GetGlobalNorm\n");
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_NORM_RULE(rule);


  MURGE_MEMALLOC(local_norm, solvers[id]->N*dof, REAL);

  for(itercol = 0; itercol <  solvers[id]->N*dof; itercol++) {
    local_norm[itercol] = 0;
  }
  for(itercol = 0; itercol <  solvers[id]->n; itercol++) {
    for (iterrow = solvers[id]->colptr[itercol]-1;
         iterrow < solvers[id]->colptr[itercol+1]-1;
         iterrow++) {
      for (iterdof_col = 0;
           iterdof_col < dof;
           iterdof_col++) {
        for (iterdof_row = 0;
             iterdof_row < dof;
             iterdof_row++) {
          value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;
          switch(rule) {
          case MURGE_NORM_MAX_COL:
            scal_idx = iterdof_col + (solvers[id]->l2g[itercol]-1) * dof;
            local_norm[scal_idx] =
              MAX(local_norm[scal_idx],
                  ABS_FLOAT(solvers[id]->values[value_idx]));
            break;

          case MURGE_NORM_MAX_ROW:
            scal_idx = iterdof_row + (solvers[id]->rows[iterrow]-1) * dof;
            local_norm[scal_idx] =
              MAX(local_norm[scal_idx],
                  ABS_FLOAT(solvers[id]->values[value_idx]));
            break;

          case MURGE_NORM_2_COL:
            scal_idx = iterdof_col + (solvers[id]->l2g[itercol]-1) * dof;
            local_norm[scal_idx] = local_norm[scal_idx] +
              (REAL)(solvers[id]->values[value_idx]*
                     CONJ_FLOAT(solvers[id]->values[value_idx]));
            break;

          case MURGE_NORM_2_ROW:
            scal_idx = iterdof_row + (solvers[id]->rows[iterrow]-1) * dof;
            local_norm[scal_idx] = local_norm[scal_idx] +
              (REAL)(solvers[id]->values[value_idx]*
                     CONJ_FLOAT(solvers[id]->values[value_idx]));
            break;

          default:
            errorPrint("Rule not implemented");
            return MURGE_ERR_NOT_IMPLEMENTED;
          }
        }
      }
    }
  }


  if (rule == MURGE_NORM_2_COL ||
      rule == MURGE_NORM_2_ROW) {
    fprintf(stderr, "Reduce on norm\n");
    MPI_Allreduce(local_norm,norm,solvers[id]->N*dof,
                  MURGE_MPI_REAL,
                  MPI_SUM,
                  solvers[id]->pastix_data->pastix_comm);
    for(itercol = 0; itercol <  solvers[id]->N*dof; itercol++) {
      local_norm[itercol] = (REAL)sqrt((double)local_norm[itercol]);
    }
  }
  else
    {
      fprintf(stderr, "Reduce on norm 2\n");
      MPI_Allreduce(local_norm,norm,solvers[id]->N*dof,
                    MURGE_MPI_REAL,
                    MPI_MAX,
                    solvers[id]->pastix_data->pastix_comm);
    }
  return MURGE_SUCCESS;
}

/*
 Function: MURGE_GetLocalNorm

 Compute the local norm array following a norm rule.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 norm    - Array of size local column number*dof which will contain
 the solution
 rule    - Rule to follow to build norm array, see <MURGE_NORM_RULES>

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_GETLOCALNORM(ID, NORM, RULE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, RULE
 >   REAL, DIMENSION(0), INTENT(OUT) :: NORM
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_GETLOCALNORM
 */
INTS MURGE_GetLocalNorm(INTS id, REAL *norm, INTS rule){
  INTS itercol;       /* Each column*/
  PASTIX_INT  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS column_dof;    /* Column number of the value */
  PASTIX_INT  value_idx;     /* Index of the value */
  INTS dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_NORM_RULE(rule);
  CHECK_L2G(id);

  if (rule == MURGE_NORM_MAX_ROW || rule == MURGE_NORM_2_ROW) {
    errorPrint("PaStiX uses column distribution, local norm can't be a norm on rows");
    return MURGE_ERR_PARAMETER;
  }

  for(itercol = 0; itercol <  solvers[id]->n*dof; itercol++)
    norm[itercol] = 0;
  for(itercol = 0; itercol <  solvers[id]->n; itercol++) {
    for (iterrow = solvers[id]->colptr[itercol]-1;
         iterrow < solvers[id]->colptr[itercol+1]-1;
         iterrow++) {
      for (iterdof_col = 0;
           iterdof_col < dof;
           iterdof_col++) {
        for (iterdof_row = 0;
             iterdof_row < dof;
             iterdof_row++) {
          column_dof = iterdof_col + itercol * dof;
          value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;

          if (rule == MURGE_NORM_2_COL) {
            norm[column_dof] = norm[column_dof] +
              (REAL)(solvers[id]->values[value_idx]*
                     CONJ_FLOAT(solvers[id]->values[value_idx]));
          }
          else
            {
              norm[column_dof] =
                MAX(norm[column_dof],
                    ABS_FLOAT(solvers[id]->values[value_idx]));
            }
        }
      }
    }
  }

  for(itercol = 0; itercol <  solvers[id]->n*dof; itercol++)
    norm[itercol] = (REAL)sqrt(norm[itercol]);

  return MURGE_SUCCESS;
}


/*
 Function: MURGE_GetNorm

 Compute the indicated part of the norm array
 following a norm rule.

 Must be performed after assembly step.


 Parameters:
 id       - Solver instance identification number.
 n        - Number of coefficients user wants to get norm of.
 coefsidx - List of the coefficients user wants to get norm of.
 norm     - Array of size dof*n which will contain
 the solution.
 rule     - Rule to follow to build norm array, see <MURGE_NORM_RULES>
 mode     - Indicates if the user is sure to respect the distribution.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_GETNORM(ID, N, COEFSIDX, NORM, RULE, MODE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, MODE, N, RULE
 >   INTS, DIMENSION(0), INTENT(IN)  :: COEFSIDX
 >   COEF, DIMENSION(0), INTENT(OUT) :: NORM
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_GETNORM
 */
INTS MURGE_GetNorm(INTS id,  INTS n, INTS *coefsidx, REAL *norm, INTS rule, INTS mode){
  errorPrint("Not yet implemented");
  return MURGE_ERR_NOT_IMPLEMENTED;

}


/*
 Function: MURGE_ApplyGlobalScaling

 Apply scaling to local unknowns.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 scal    - Scaling user wants to apply.
 sc_mode - Indicate if the scaling is applied on rows or on columns.
 root    - Indicates which processor that posses the scaling array,
 -1 for all.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_APPLYGLOBALSCALING(ID, SCAL, SC_MODE, ROOT, IERROR)
 >   INTS,               INTENT(IN)  :: ID, ROOT, SC_MODE
 >   REAL, DIMENSION(0), INTENT(OUT) :: SCAL
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_APPLYGLOBALSCALING

 */
INTS MURGE_ApplyGlobalScaling(INTS id, REAL *scal, INTS root, INTS sc_mode){
  INTS itercol;       /* Each column*/
  PASTIX_INT  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS scal_idx;      /* Scaling array index */
  PASTIX_INT  value_idx;     /* Index of the value */
  INTS dof;
  REAL *scaling = NULL;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_SCALING_MODE(sc_mode);
  CHECK_L2G(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (root == -1) {
    scaling = scal;
  }
  else
    {
      if (root != (solvers[id]->pastix_data)->procnum) {
        MURGE_MEMALLOC(scaling, solvers[id]->N*dof, REAL);
      }
      else
        {
          scaling = scal;
        }
      MPI_Bcast( scaling, solvers[id]->N*dof,
                 MURGE_MPI_REAL,
                 root,
                 solvers[id]->pastix_data->pastix_comm);
    }

  for(itercol = 0; itercol <  solvers[id]->n; itercol++) {
    for (iterrow = solvers[id]->colptr[itercol]-1;
         iterrow < solvers[id]->colptr[itercol+1]-1;
         iterrow++) {
      for (iterdof_col = 0;
           iterdof_col < dof;
           iterdof_col++) {
        for (iterdof_row = 0;
             iterdof_row < dof;
             iterdof_row++) {
          if (sc_mode == MURGE_SCAL_COL) {
            scal_idx = iterdof_col + (solvers[id]->l2g[itercol]-1) * dof;
          }
          else
            {
              scal_idx = iterdof_row + (solvers[id]->rows[iterrow]-1) * dof;
            }
          value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;
          solvers[id]->values[value_idx] =
            solvers[id]->values[value_idx] /scaling[scal_idx];
        }
      }
    }
  }
  if (root != -1 && root != (solvers[id]->pastix_data)->procnum)
    MURGE_FREE(scaling);

  return MURGE_SUCCESS;
}

/*
 Function: MURGE_ApplyLocalScaling

 Apply the local scaling array on the matrix.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 scal    - Array of size local column number*dof which will contain
 the solution.
 sc_mode - Indicate if the scaling is applied on rows or on columns.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_APPLYLOCALSCALING(ID, SCAL, SC_MODE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, SC_MODE
 >   REAL, DIMENSION(0), INTENT(OUT) :: SCAL
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_APPLYLOCALSCALING
 */
INTS MURGE_ApplyLocalScaling(INTS id, REAL *scal, INTS sc_mode){
  INTS itercol;       /* Each column*/
  PASTIX_INT  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS scal_idx;      /* Index in scaling array */
  PASTIX_INT  value_idx;     /* Index of the value */
  INTS dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_SCALING_MODE(sc_mode);
  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (sc_mode == MURGE_SCAL_ROW) {
    /*
     * Building global to local column number array
     */
    CHECK_L2G(id);
  }
  for(itercol = 0; itercol <  solvers[id]->n; itercol++) {
    for (iterrow = solvers[id]->colptr[itercol]-1;
         iterrow < solvers[id]->colptr[itercol+1]-1;
         iterrow++) {
      for (iterdof_col = 0;
           iterdof_col < dof;
           iterdof_col++) {
        for (iterdof_row = 0;
             iterdof_row < dof;
             iterdof_row++) {

          if (sc_mode == MURGE_SCAL_COL) {
            scal_idx = iterdof_col + itercol * dof;
          }
          else
            {
              scal_idx = iterdof_row + (solvers[id]->g2l[solvers[id]->rows[iterrow]-1]-1)*dof;
            }
          value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;
          solvers[id]->values[value_idx] =
            solvers[id]->values[value_idx] /scal[scal_idx];
        }
      }
    }
  }
  return MURGE_SUCCESS;
}

/*
 Function: MURGE_ApplyScaling

 Apply the scaling array on the indicated part of the matrix

 Must be performed after assembly step.


 Parameters:
 id       - Solver instance identification number.
 n        - Number of coefficients user wants to scale.
 coefsidx - List of the coefficients user wants to scale.
 scal     - Array of size dof*n which will contain
 the solution.
 sc_mode  - Indicate if the scaling is applied on rows or on columns.
 mode     - Indicates if the user is sure to respect the distribution.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_APPLYSCALING(ID, N, COEFSIDX, SCAL, SC_MODE, MODE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, SC_MODE, MODE, N
 >   INTS, DIMENSION(0), INTENT(IN)  :: COEFSIDX
 >   COEF, DIMENSION(0), INTENT(OUT) :: SCAL
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_APPLYSCALING
 */
INTS MURGE_ApplyScaling(INTS id,  INTS n, INTS *coefsidx, REAL *scal,
                        INTS sc_mode, INTS mode){
  INTS itercol;       /* Each column*/
  INTS iterdof_col;   /* each dof on column */
  REAL *scaling = NULL;
  INTS dof;
  INTS baseval;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_SCALING_MODE(sc_mode);
  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (mode == MURGE_ASSEMBLY_RESPECT) {
    /*
     * Building global to local column number array
     */
    CHECK_L2G(id);

    MURGE_MEMALLOC(scaling, solvers[id]->n*dof, REAL);
    for (itercol = 0; itercol < solvers[id]->n*dof; itercol++)
      scaling[itercol] = 1.0;
    for (itercol = 0; itercol < n; itercol++)
      for (iterdof_col =0; iterdof_col < dof; iterdof_col++)
        scaling[(solvers[id]->g2l[coefsidx[itercol]-baseval]-1)*dof+iterdof_col] = scal[itercol*dof+iterdof_col];
    MURGE_ApplyLocalScaling(id, scaling, sc_mode);
    MURGE_FREE(scaling);
  }
  else
    {
      REAL * scaling_recv = NULL;
      MURGE_MEMALLOC(scaling, solvers[id]->N*dof, REAL);
      for (itercol = 0; itercol < solvers[id]->N*dof; itercol++)
        scaling[itercol] = 0.0;
      for (itercol = 0; itercol < n; itercol++)
        for (iterdof_col =0; iterdof_col < dof; iterdof_col++)
          scaling[(coefsidx[itercol]-baseval)*dof+iterdof_col] = scal[itercol*dof+iterdof_col];

      MURGE_MEMALLOC(scaling_recv, solvers[id]->N*dof, REAL);
      MPI_Allreduce(scaling, scaling_recv,
                    solvers[id]->N*dof,
                    MURGE_MPI_REAL,
                    MPI_SUM,
                    solvers[id]->pastix_data->pastix_comm);
      MURGE_FREE(scaling);
      for (itercol = 0; itercol < solvers[id]->N*dof; itercol++)
        if (scaling_recv[itercol] == 0.0)
          scaling_recv[itercol] = 1.0;

      for (itercol = 0; itercol < n; itercol++) {
        for (iterdof_col =0; iterdof_col < dof; iterdof_col++) {
          if (scaling[(coefsidx[itercol]-baseval)*dof+iterdof_col] != scal[itercol*dof+iterdof_col]) {
            errorPrint("Multiple entries for the same scaling entry");
            return MURGE_ERR_PARAMETER;
          }
        }
      }


      MURGE_ApplyGlobalScaling(id, scaling_recv, sc_mode, -1);
      MURGE_FREE(scaling_recv);
    }
  return MURGE_SUCCESS;
}





/******************************************************************************
 * Group: Specific PaStiX functions.                                          *
 ******************************************************************************/

/******************************************************************************
 * Function: MURGE_Analyze                                                    *
 *                                                                            *
 * Perform matrix analyze:                                                    *
 *   - Compute a new ordering of the unknows                                  *
 *   - Compute the symbolic factorisation of the matrix                       *
 *   - Distribute column blocks and computation on processors                 *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned succesfuly.                     *
 *   MURGE_ERR_ORDER     - If function the graph is not built.                *
 *   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_Analyze(INTS id){
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);

  return MURGE_SUCCESS;
}

/******************************************************************************
 * Function: MURGE_Factorize                                                  *
 *                                                                            *
 * Perform matrix factorization.                                              *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned succesfuly.                     *
 *   MURGE_ERR_ORDER     - If function the graph is not built.                *
 *   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_Factorize(INTS id){
  pastix_data_t   *pastix_data = solvers[id]->pastix_data;
  PASTIX_INT             *iparm       = pastix_data->iparm;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_VALUES_OK))) {
    errorPrint("Need to set values before.");
    return MURGE_ERR_ORDER;
  }
  if (iparm[IPARM_ONLY_RAFF] ==  API_YES) {
    errorPrint("MURGE_Factorize is not compatible with IPARM_ONLY_RAFF == API_YES\n");
    return MURGE_ERR_PARAMETER;
  }

  /* 
   * - fill the internal CSC
   * - delete murge CSC and
   * - perform factorization
   */
  PASTIX_FILLIN_CSC(solvers[id]->pastix_data,
		    solvers[id]->pastix_data->pastix_comm,
		    solvers[id]->n,
		    solvers[id]->colptr,
		    solvers[id]->rows,
		    solvers[id]->values,
		    NULL,
		    0,
		    solvers[id]->l2g);
  solvers[id]->pastix_data->cscInternFilled = API_YES;

  MURGE_FREE(solvers[id]->colptr);
  MURGE_FREE(solvers[id]->rows);
  MURGE_FREE(solvers[id]->values);

  pastix_data->iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
  iparm[IPARM_END_TASK]                = API_TASK_NUMFACT;

  DPASTIX(&pastix_data,
          pastix_data->pastix_comm,
          solvers[id]->n,
          solvers[id]->colptr,
          solvers[id]->rows,
          solvers[id]->values,
          solvers[id]->l2g,
          solvers[id]->perm,
          NULL,
          NULL,
          solvers[id]->nrhs,
          pastix_data->iparm,
          pastix_data->dparm);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_FACTO_OK);

  return MURGE_SUCCESS;
}

/******************************************************************************
 * Function: MURGE_ForceNoFacto                                               *
 *                                                                            *
 * Prevent Murge from running factorisation even if matrix has changed.       *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 * Returns:                                                                   *
 *   MURGE_SUCCESS                                                            *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_ForceNoFacto(INTS id) {
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_FACTO_OK);
  return MURGE_SUCCESS;
}

/******************************************************************************
 * Function: MURGE_ProductSetLocalNodeNbr                                     *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   n  - Number of local nodes.                                              *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_ProductSetLocalNodeNbr (INTS id, INTS n) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  solvers[id]->n = n;
  MPI_Allreduce(&solvers[id]->n,
                &solvers[id]->N, 1, COMM_INT,
                MPI_SUM, solvers[id]->pastix_data->pastix_comm);
  return MURGE_SUCCESS;
}


/******************************************************************************
 * Function: MURGE_ProductSetGlobalNodeNbr                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   N  - Number of global nodes.                                             *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_ProductSetGlobalNodeNbr (INTS id, INTS N) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_NODELST_OK)) {
    errorPrint("%s must be called before MURGE_ProductSetLocalNodeList",
               __FUNCTION__);
    return MURGE_ERR_ORDER;
  }
  solvers[id]->N = N;
  return MURGE_SUCCESS;
}
/******************************************************************************
 * Function: MURGE_ProductSetLocalNodeList                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id  - Solver instance identification number.                             *
 *   l2g - Local to global node numbers.                                      *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_ProductSetLocalNodeList (INTS id, INTS * l2g) {
  INTS i;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  MURGE_MEMALLOC(solvers[id]->l2g, solvers[id]->n, PASTIX_INT);
  for (i = 0; i < solvers[id]->n; i++) {
    solvers[id]->l2g[i] = l2g[i];
  }

  cscd_build_g2l(solvers[id]->n,
                 solvers[id]->l2g,
                 solvers[id]->pastix_data->pastix_comm,
                 &solvers[id]->N,
                 &solvers[id]->g2l);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(solvers[id]->g2l), char);

  /* No need to work on the graph nor factorize
   when we only perform product */
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_ONLY_PROD);
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_NODELST_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);


  return MURGE_SUCCESS;
}

/******************************************************************************
 * Function: MURGE_GetLocalProduct                                            *
 *                                                                            *
 * Perform the product A * X.                                                 *
 *                                                                            *
 * The vector must have been given trough <MURGE_SetLocalRHS> or              *
 * <MURGE_SetGlobalRHS>.                                                      *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   x  - Array in which the local part of the product will be stored.        *
 * Returns:                                                                   *
 *   MURGE_ERR_ORDER  - If values have not been set.                          *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_GetLocalProduct (INTS id, COEF *x) {
  COEF * glob_prod;
  INTS ierr, iter;
  MURGE_MEMALLOC(glob_prod, solvers[id]->N, COEF);

  if (MURGE_SUCCESS != (ierr = MURGE_GetGlobalProduct(id, glob_prod, -1)))
    return ierr;

  for (iter  = 0; iter < solvers[id]->n; iter++) {
    x[iter] = glob_prod[solvers[id]->l2g[iter]];
  }
  return MURGE_SUCCESS;
}


/******************************************************************************
 * Function: MURGE_GetGlobalProduct                                           *
 *                                                                            *
 * Perform the product A * X.                                                 *
 *                                                                            *
 * The vector must have been given trough <MURGE_SetLocalRHS> or              *
 * <MURGE_SetGlobalRHS>.                                                      *
 *                                                                            *
 * Parameters:                                                                *
 *   id   - Solver instance identification number.                            *
 *   x    - Array in which the product will be stored.                        *
 *   root - Rank of the process which will own the product at end of call,    *
 *          use -1 for all processes.                                         *
 * Returns:                                                                   *
 *   MURGE_ERR_ORDER  - If values have not been set.                          *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_GetGlobalProduct (INTS id, COEF *x, INTS root) {
  INTS dof;

  COEF * my_prod = NULL;
#ifdef MURGE_TIME
  Clock            clock;
#endif

  CLOCK_INIT;
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_VALUES_OK))) {
    errorPrint("Need to set values before.");
    return MURGE_ERR_ORDER;
  }

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);
  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (solvers[id]->pastix_data->iparm[IPARM_SYM] == API_SYM_YES) {
    errorPrint("Product only available with unsymmetric matrices.");
    return MURGE_ERR_NOT_IMPLEMENTED;
  }

#ifdef MURGE_PRODUCT_CHECK_ZEROS
  {
    INTS itercol, iterrows, row;
    COEF max, sum, norm, critere;
    PASTIX_INT cnt, cnt2, cnt_sum, cnt2_sum, nz_glob;
    COEF *values = solvers[id]->values;
    INTS baseval, iter, iter2;
    critere = solvers[id]->pastix_data->dparm[DPARM_EPSILON_MAGN_CTRL];
    baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];

    if (critere < 0.0) {
      critere = -critere;
    }
    else
      {
        max = 0;
        for (itercol = 0; itercol < solvers[id]->n; itercol++) {
          for (iter = 0; iter < dof; iter++) {
            sum = 0;
            for (iterrows = solvers[id]->colptr[itercol]-baseval;
                 iterrows < solvers[id]->colptr[itercol+1]-baseval;
                 iterrows++) {
              row = solvers[id]->rows[iterrows]-baseval;
              for (iter2 = 0; iter2 < dof; iter2++) {
                PASTIX_INT idx = iterrows*dof*dof+iter*dof+iter2;

                sum = sum + ABS_FLOAT(values[idx]);
              }

            }
            max = MAX(max,sum);
          }
        }

        MPI_Allreduce(&max, &norm, 1, MURGE_MPI_COEF, MPI_MAX,
                      solvers[id]->pastix_data->pastix_comm);

        critere = norm*sqrt(critere);
      }
    cnt = 0;
    cnt2 = 0;
    for (itercol = 0; itercol < solvers[id]->n; itercol++) {
      for (iter = 0; iter < dof; iter++) {
        sum = 0;
        for (iterrows = solvers[id]->colptr[itercol]-baseval;
             iterrows < solvers[id]->colptr[itercol+1]-baseval;
             iterrows++) {
          row = solvers[id]->rows[iterrows]-baseval;
          for (iter2 = 0; iter2 < dof; iter2++) {
            PASTIX_INT idx = iterrows*dof*dof+iter*dof+iter2;
            if (ABS_FLOAT(values[idx]) < critere)
              cnt = cnt + 1;
            if (values[idx] ==  0.0)
              cnt2 = cnt2 + 1;
          }
        }
        max = MAX(max,sum);
      }
    }
    cnt_sum = 0;
    MPI_Reduce(&cnt, &cnt_sum, 1, COMM_INT, MPI_SUM, 0,
               solvers[id]->pastix_data->pastix_comm);
    MPI_Reduce(&cnt2, &cnt2_sum, 1, COMM_INT, MPI_SUM, 0,
               solvers[id]->pastix_data->pastix_comm);
    cnt = solvers[id]->colptr[solvers[id]->n]-1;
    MPI_Reduce(&cnt, &nz_glob, 1, COMM_INT, MPI_SUM, 0,
               solvers[id]->pastix_data->pastix_comm);
    nz_glob = nz_glob *dof*dof;
    if ((solvers[id]->pastix_data)->procnum == 0) {
      fprintf(stdout, "%d zeros in matrix from %d (%.3lg %%) critere :"
              " %.20lg\n",
              cnt_sum, nz_glob,
              (100.0*(double)(cnt_sum)/(double)(nz_glob)), critere);
      fprintf(stdout, "%d real zeros in matrix from %d (%.3lg %%)\n",
              cnt2_sum, nz_glob,
              (100.0*(double)(cnt2_sum)/(double)(nz_glob)));
    }
  }
#endif

#ifndef FORCE_NOSMP
  if (solvers[id]->threadnbr != 0 &&
      solvers[id]->threadnbr != solvers[id]->pastix_data->iparm[IPARM_THREAD_NBR])
    stop_threads(id);
  if (solvers[id]->threadnbr == 0)
    start_threads(id);
  pthread_mutex_lock(&(solvers[id]->mutex_state));
  solvers[id]->threads_state = MURGE_THREAD_PRODUCT;
  pthread_mutex_unlock(&(solvers[id]->mutex_state));
  pthread_cond_broadcast(&(solvers[id]->cond_state));
  pthread_mutex_lock(&(solvers[id]->mutex_state));
  while (solvers[id]->threads_state != MURGE_THREAD_WAIT) {
    pthread_cond_wait(&(solvers[id]->cond_state),
                      &(solvers[id]->mutex_state));
  }

  if (solvers[id]->threads_data[0].pdata->ret != MURGE_SUCCESS)
    return solvers[id]->threads_data[0].pdata->ret;
#else
  product_thread(solvers[id]->threads_data[0].pdata);
#endif
  my_prod = solvers[id]->threads_data[0].pdata->t_prod;

  if (root == -1) {
    MPI_Allreduce(my_prod, x, solvers[id]->N*dof, MURGE_MPI_COEF, MPI_SUM,
                  solvers[id]->pastix_data->pastix_comm);
  }
  else {
    MPI_Reduce(my_prod, x, solvers[id]->N*dof, MURGE_MPI_COEF, MPI_SUM,
               root, solvers[id]->pastix_data->pastix_comm);
  }
  CLOCK_PRINT("MURGE_GetGlobalProduct");

  return MURGE_SUCCESS;
}

/*
 WARNING: NEEDS TO BE CHECKED !
 */
INTS MURGE_SetLocalNodeList   (INTS id, INTS nodenbr, INTS *nodelist) {
  PASTIX_INT i;
  solvers[id]->n = nodenbr;
  /* On detruit colptr et rows, la distribution a changé */
  MURGE_FREE(solvers[id]->colptr);
  MURGE_FREE(solvers[id]->rows);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_GRAPH_OK)
    MURGE_STATE_TRUE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_SYMB_OK);
  MURGE_MEMALLOC(solvers[id]->l2g, nodenbr, PASTIX_INT);
  for (i = 0; i < nodenbr; i++) {
    solvers[id]->l2g[i] = nodelist[i];
  }
  return MURGE_SUCCESS;
}

INTS MURGE_GetCommRank(INTS id, int * rank) {
  CHECK_SOLVER_ID(id);
  *rank = (solvers[id]->pastix_data)->procnum;
  return MURGE_SUCCESS;

}

INTS MURGE_GetCommSize(INTS id, int * size) {
  CHECK_SOLVER_ID(id);
  *size = (solvers[id]->pastix_data)->procnbr;
  return MURGE_SUCCESS;
}

INTS MURGE_GetOptionINT(INTS id, INTS index, INTS * value) {
  PASTIX_INT murge_param[64];
  PASTIX_INT * iparm = NULL;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  iparm = solvers[id]->pastix_data->iparm;

  murge_param[MURGE_IPARAM_BASEVAL       - 1024] =  IPARM_BASEVAL;
  murge_param[MURGE_IPARAM_DOF           - 1024] =  IPARM_DOF_NBR;
  murge_param[MURGE_IPARAM_SYM           - 1024] =  IPARM_SYM;

  if (index == MURGE_IPARAM_SYM) {
    if (iparm[IPARM_SYM] == API_SYM_YES)
      *value = MURGE_BOOLEAN_TRUE;
    else
      *value = MURGE_BOOLEAN_FALSE;
  }
  else
    {
      if (index >= 1024)
        index = murge_param[index-1024];

      *value = iparm[index];
    }
  return MURGE_SUCCESS;
}


INTS MURGE_GetComm(INTS id, MPI_Comm * comm) {
  CHECK_SOLVER_ID(id);
  * comm = solvers[id]->pastix_data->pastix_comm;
  return MURGE_SUCCESS;
}

INTS MURGE_SetDropNodes(INTS id, INTS nodenbr, INTS * dropmask) {
  INTS i;
  CHECK_SOLVER_ID(id);
  MURGE_MEMALLOC(solvers[id]->dropmask, nodenbr, char);
  for (i = 0; i < nodenbr; i++) {
    solvers[id]->dropmask[i] = (char)dropmask[i];
  }
  return MURGE_SUCCESS;
}


INTS MURGE_SetDropRows(INTS id, INTS nodenbr, INTS * droprows) {
  INTS i;
  CHECK_SOLVER_ID(id);
  MURGE_MEMALLOC(solvers[id]->droprows, nodenbr, char);
  for (i = 0; i < nodenbr; i++) {
    solvers[id]->droprows[i] = (char)droprows[i];
  }
  return MURGE_SUCCESS;
}

INTS MURGE_SetDropCols(INTS id, INTS nodenbr, INTS * dropcols) {
  INTS i;
  CHECK_SOLVER_ID(id);
  MURGE_MEMALLOC(solvers[id]->dropcols, nodenbr, char);
  for (i = 0; i < nodenbr; i++) {
    solvers[id]->dropcols[i] = (char)dropcols[i];
  }
  return MURGE_SUCCESS;
}

INTS MURGE_ColGetNonZerosNbr(INTS id, INTS COL, INTS * nnzNbr) {
  INTS lcol;
  INTS base = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  CHECK_SOLVER_ID(id);
  CHECK_PREPROCESSING(id);
  COL = COL - base;
  if (COL < 0 || COL > solvers[id]->N) {
    errorPrint("invalid column index");
    return MURGE_ERR_PARAMETER;
  }
  lcol = solvers[id]->g2l[COL];
  if (lcol > 0)
    *nnzNbr = (INTS)(solvers[id]->colptr[lcol] - solvers[id]->colptr[lcol-1]);
  else
    *nnzNbr = 0;
  return MURGE_SUCCESS;
}

INTS MURGE_ColGetNonZerosIdx(INTS id, INTS COL, INTS * indexes) {
  INTS lcol;
  INTS base = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  CHECK_SOLVER_ID(id);
  CHECK_PREPROCESSING(id);
  COL = COL - base;
  if (COL < 1 || COL > solvers[id]->N) {
    errorPrint("invalid column index");
    return MURGE_ERR_PARAMETER;
  }
  lcol = solvers[id]->g2l[COL];
  if (lcol > 0) {
    INTS i;
    for (i = solvers[id]->colptr[lcol-1]-1; i < solvers[id]->colptr[lcol]-1; i++)
      indexes[i- solvers[id]->colptr[lcol-1]-1] = solvers[id]->rows[i];
  }
  return MURGE_SUCCESS;
}
