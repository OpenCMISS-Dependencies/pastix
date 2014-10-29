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
 * File: murge_defines.h
 *
 *  This file define defines, macros and external funtion definition
 *  used to build <Murge> interface.
 *
 * About: Authors
 *   Xavier Lacoste  - xavier.lacoste@inria.fr
 */
#ifndef MURGE_DEFINE_H
#define MURGE_DEFINE_H
#include "pastix_internal.h"
#include "cscd_utils_intern.h"

/******************************************************************************/
/***                           Section: Defines                             ***/
/******************************************************************************/

/*
 * Defines: MPI Tags
 *
 * Tags for MPI communication.
 *
 * TAG_SIZE - To send size of a buffer.
 * TAG_ROW  - To send rows.
 * TAG_COL  - To send columns.
 * TAG_FL   - To send an interval (first last).
 * TAG_IJV  - To send <ijv_t> array.
 * TAG_L2G  - To send local to global column numbers.
 * TAG_VAL  - To send values.
 */
#define TAG_SIZE  1
#define TAG_ROW   2
#define TAG_COL   3
#define TAG_FL    4
#define TAG_IJV   5
#define TAG_L2G   6
#define TAG_VAL   7
#define TAG_SIZE2 8
#define TAG_VAL2  9
#define TAG_IJV2  10

/*
 * Defines: State masks
 *
 * Bit masks used to define state variable.
 *
 *   MURGE_INIT_OK     - If initialisation step has been called.
 *   MURGE_GRAPH_OK    - If graph of non zeros has been built.
 *   MURGE_GRAPH_BUILD - If we are in a graph building session.
 *   MURGE_VALUES_OK   - If Values of the matrix have been set.
 *   MURGE_BLEND_OK    - If preprocessing has been performed.
 *   MURGE_MATR_BUILD  - If we are in a matrix building session.
 *   MURGE_FACTO_OK    - If Factorization has been computed.
 *   MURGE_NODENBR_OK  - If node number has been given to user.
 *   MURGE_NODELST_OK  - If node list has been given to user.
 *   MURGE_RHS_OK      - If Right hand side has been set by user.
 *   MURGE_SYMB_OK     - If Symbolic factorization has been performed.
 *   MURGE_ONLY_PROD   - If we only compute producte
 */
#define MURGE_INIT_OK     1
#define MURGE_GRAPH_OK    2
#define MURGE_GRAPH_BUILD 4
#define MURGE_VALUES_OK   8
#define MURGE_BLEND_OK    16
#define MURGE_MATR_BUILD  32
#define MURGE_FACTO_OK    64
#define MURGE_NODENBR_OK  128
#define MURGE_NODELST_OK  256
#define MURGE_RHS_OK      512
#define MURGE_SYMB_OK     1024
#define MURGE_ONLY_PROD   2048
#define MURGE_SOLVE_DONE  4096
#define MURGE_REFINE_DONE 8192

/*
 * Defines: outputs strings
 *
 * OUT_TIMER        - For timing prints.
 * OUT_MEMORY       - Current memory usage in murge.
 * OUT_MAXMEMORY    - Maximum memory usage in murge.
 * OUT_CALL_DPASTIX - When calling dpastix().
 * OUT_CALL_PASTIX  - When calling pastix().
 * OUT_ERR_ALLOC    - Allocation error message.
 */
#define OUT_TIMER        " > MURGE < %-40s computed in          %8.3g s  (id=%d)\n"
#define OUT_MEMORY       " > MURGE < %-40s memory used          %8.3g %-2s (id=%d)\n"
#define OUT_MAXMEMORY    " > MURGE < %-40s maximum memory used  %8.3g %-2s (id=%d)\n"
#define OUT_CALL_DPASTIX            "CALL distributed PaStiX [%d-%d]"
#define OUT_CALL_PASTIX             "CALL centralised PaStiX [%d-%d]"
#define OUT_ERR_ALLOC    "%s:%d Memory allocation error"

/******************************************************************************/
/***                           Section: Macros                              ***/
/******************************************************************************/

#ifdef MURGE_THREADSAFE
#  define MEMORY_LOCK   pthread_mutex_lock(&solvers[id]->mutex_memory);
#  define MEMORY_UNLOCK pthread_mutex_unlock(&solvers[id]->mutex_memory);
#  define MATRIX_LOCK   pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#  define MATRIX_UNLOCK pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#else
#  define MEMORY_LOCK
#  define MEMORY_UNLOCK
#  define MATRIX_LOCK
#  define MATRIX_UNLOCK
#endif

#ifdef MEMORY_USAGE
#  define MURGE_TRACE_MALLOC(size, type)                \
  do {                                                  \
    MEMORY_LOCK;                                        \
    solvers[id]->malloc_size += (size)*sizeof(type);    \
    if (solvers[id]->malloc_size >                      \
        solvers[id]->malloc_maxsize)                    \
      solvers[id]->malloc_maxsize =                     \
        solvers[id]->malloc_size;                       \
    MEMORY_UNLOCK;                                      \
  } while(0)
#  define MURGE_TRACE_FREE(ptr)                                         \
  do {                                                                  \
    if (ptr != NULL) {                                                  \
      double * tf_memptr = (double*)(ptr);                              \
      uint64_t tf_size;                                                 \
      tf_memptr --;                                                     \
      tf_size = (uint64_t) tf_memptr[0];                                \
      MEMORY_LOCK;                                                      \
      solvers[id]->malloc_size -= tf_size;                              \
      MEMORY_UNLOCK;                                                    \
    }                                                                   \
  } while(0)
#  define MURGE_TRACE_FREE_EXT(size)                                    \
  do {                                                                  \
    MEMORY_LOCK;                                                        \
    solvers[id]->malloc_size -= (uint64_t)size;                         \
    MEMORY_UNLOCK;                                                      \
  } while(0)

#  define MURGE_MEMORY_USAGE_PRINT(string)                              \
  MURGE_MEMORY_USAGE_PRINT2(string,                                     \
                            solvers[id]->pastix_data->iparm[IPARM_VERBOSE], \
                            solvers[id]->pastix_data->procnum,          \
                            solvers[id]->pastix_data->pastix_comm)
#  define MURGE_MEMORY_USAGE_PRINT2(string, verb, rank, comm)           \
  do {                                                                  \
    if (verb > API_VERBOSE_NO) {                                        \
      int64_t mem, maxmem;                                              \
      MPI_Reduce(&(solvers[id]->malloc_size),                           \
                 &mem,                                                  \
                 1, MPI_INTEGER8, MPI_MAX, 0,                           \
                 comm);                                                 \
      MPI_Reduce(&(solvers[id]->malloc_maxsize),                        \
                 &maxmem,                                               \
                 1, MPI_INTEGER8, MPI_MAX, 0,                           \
                 comm);                                                 \
      if (rank == 0) {                                                  \
        fprintf(stdout,                                                 \
                OUT_MEMORY,                                             \
                string,                                                 \
                MEMORY_WRITE(solvers[id]->malloc_size),                 \
                MEMORY_UNIT_WRITE(solvers[id]->malloc_size), id);       \
        fprintf(stdout,                                                 \
                OUT_MAXMEMORY,                                          \
                string,                                                 \
                MEMORY_WRITE(solvers[id]->malloc_maxsize),              \
                MEMORY_UNIT_WRITE(solvers[id]->malloc_maxsize), id);    \
      }                                                                 \
    }                                                                   \
  } while(0)
#else
#  define MURGE_TRACE_MALLOC(type, size)
#  define MURGE_TRACE_FREE(ptr)
#  define MURGE_TRACE_FREE_EXT(size)
#  define MURGE_MEMORY_USAGE_PRINT(string)
#  define MURGE_MEMORY_USAGE_PRINT2(string, verb, rank, comm)
#endif
/*
 * Macro: MURGE_MEMALLOC
 *
 * Allocate a space of size *size* x sizeof(*type*)
 * at the adress indicated by ptr.
 *
 * Parameters:
 *   ptr   - address where to allocate.
 *   size  - Number of elements to allocate.
 *   types - Type of the elements to allocate.
 *
 * Returns:
 * MURGE_ERR_ALLOCATE - If allocation fails.
 */
#define MURGE_MEMALLOC(ptr, size, type)                                 \
  do {                                                                  \
    if ((size)*sizeof(type) == 0) {                                     \
      ptr = NULL;                                                       \
    } else {                                                            \
      if (NULL == ((ptr) = (type *) memAlloc((size) * sizeof(type)))) { \
        errorPrint(OUT_ERR_ALLOC,__FILE__, __LINE__); \
        return MURGE_ERR_ALLOCATE;                                      \
      }                                                                 \
      PRINT_ALLOC(ptr, ((size)*sizeof(type)), __FILE__, __LINE__);      \
    }                                                                   \
  } while(0)

/*
 * Macro: MURGE_MEMALLOC_RET
 *
 * Allocate a space of size *size* x sizeof(*type*)
 * at the adress indicated by ptr.
 *
 * Parameters:
 *   ptr   - address where to allocate.
 *   size  - Number of elements to allocate.
 *   types - Type of the elements to allocate.
 *   ret   - return value.
 *
 * Returns:
 * MURGE_ERR_ALLOCATE - If allocation fails.
 */
#define MURGE_MEMALLOC_RET(ptr, size, type, ret)                        \
  do {                                                                  \
    if (size*sizeof(type) == 0) {                                       \
      ptr = NULL;                                                       \
    } else {                                                            \
      if (NULL == ((ptr) = (type *) memAlloc((size) * sizeof(type)))) { \
        errorPrint(OUT_ERR_ALLOC,__FILE__, __LINE__);                   \
        ret = MURGE_ERR_ALLOCATE;                                       \
      }                                                                 \
      PRINT_ALLOC(ptr, ((size)*sizeof(type)), __FILE__, __LINE__);      \
    }                                                                   \
  } while(0)

/*
 * Macro: MURGE_MEMALLOC_EXT
 *
 * Allocate a space of size *size* x sizeof(*type*)
 * at the adress indicated by ptr.
 *
 * Use classical malloc to avoid counting int with memAlloc.
 *
 * Parameters:
 *   ptr   - address where to allocate.
 *   size  - Number of elements to allocate.
 *   types - Type of the elements to allocate.
 *
 * Returns:
 * MURGE_ERR_ALLOCATE - If allocation fails.
 */
#define MURGE_MEMALLOC_EXT(ptr, size, type)                             \
    do {                                                                \
      if (size*sizeof(type) == 0) {                                     \
        ptr = NULL;                                                     \
      } else {                                                          \
        if (NULL == ((ptr) = (type *) malloc((size) * sizeof(type)))) { \
          errorPrint(OUT_ERR_ALLOC);                                    \
          return MURGE_ERR_ALLOCATE;                                    \
        }                                                               \
        MURGE_TRACE_MALLOC(size, type);                                 \
        PRINT_ALLOC(ptr, ((size)*sizeof(type)), __FILE__, __LINE__);    \
      }                                                                 \
    } while(0)

/*
 * Macro: MURGE_REALLOC
 *
 * Reallocate a space of size *size* x sizeof(*type*)
 * over the adress indicated by ptr.
 *
 *
 * Parameters:
 *   ptr   - address where to allocate.
 *   size  - Number of elements to allocate.
 *   types - Type of the elements to allocate.
 *
 * Returns:
 * MURGE_ERR_ALLOCATE - If allocation fails.
 */
#define MURGE_REALLOC(ptr, size, type)                                  \
  do {                                                                  \
    if (size*sizeof(type) == 0) {                                       \
      ptr = NULL;                                                       \
    } else {                                                            \
      if (NULL == ((ptr) = (type *) memRealloc(ptr,                     \
                                               (size) *                 \
                                               sizeof(type)))) {        \
        errorPrint(OUT_ERR_ALLOC,__FILE__, __LINE__);                   \
        return MURGE_ERR_ALLOCATE;                                      \
      }                                                                 \
    }                                                                   \
  } while(0)
/*
 * Macro: MURGE_FREE
 *
 * Free a memory area and nullify the pointer.
 *
 * Parameter:
 *   ptr - Pointer to the memory area.
 * Returns:
 *   nothing
 */
#define MURGE_FREE(ptr)                                                 \
  do {                                                                  \
    MURGE_TRACE_FREE(ptr);                                              \
    memFree_null(ptr);                                                  \
  } while(0)
/*
 * Macro: MURGE_FREE_EXT
 *
 * Free a memory area allocated with classical malloc
 * and nullify the pointer.
 *
 * Parameter:
 *   ptr - Pointer to the memory area.
 * Returns:
 *   nothing
 */
#define MURGE_FREE_EXT(ptr, size)                                       \
  do {                                                                  \
    PRINT_DEALLOC_EXT(ptr, size,                                        \
                      __FILE__, __LINE__);                              \
    MURGE_TRACE_FREE_EXT(size);                                         \
    free(ptr); ptr = NULL;                                              \
  } while(0)
/*
 * Macro: MURGE_STATE_FALSE
 *
 * Set *states* bit corresponding to *mask* to 0.
 *
 * Parameters:
 *   state - state variable
 *   mask  - information we want to set.
 */
#define MURGE_STATE_FALSE(state, mask)          \
  {                                             \
    state |= ( mask );                          \
    state ^= ( mask );                          \
  }


/*
 * Macro: MURGE_STATE_TRUE
 *
 * Set *states* bit corresponding to *mask* to 1.
 *
 * Parameters:
 *   state - state variable
 *   mask  - information we want to set.
 */
#define MURGE_STATE_TRUE(state, mask)           \
  {                                             \
    state |= (mask);                            \
  }

/*
 * Macro: MURGE_STATE_ISTRUE
 *
 * Check if *states* bit corresponding to *mask* is set to 1.
 *
 * Parameters:
 *   state - state variable
 *   mask  - information we want to test.
 *
 * Returns:
 *   true  - if *mask* bit is set to 1
 *   false - otherwise
 */
#define MURGE_STATE_ISTRUE(state, mask) ( state == ( (state) | (mask) ) )

/*
 * Macro: CHECK_SOLVER_ID
 *
 * Checks if solvers structure has been correctly set for instance *id*
 *
 * It checks that *id* value is in correct range and that *solvers*
 * and solvers[id] have been allocated.
 *
 * Parameters:
 *   id - Solver instance ID we want to check.
 *
 * Returns:
 *   MURGE_ERR_PARAMETER - If *id* is not in correct range.
 *   MURGE_ERR_ORDER     - If *solvers* or *solvers[id]* are not allocated.
 */
#ifdef CENTRALISED
#  define INIT_CENTRALISED                      \
  solvers[id]->invp        = NULL;
#else
#  define INIT_CENTRALISED {}
#endif
#define CHECK_SOLVER_ID(id)                                             \
  {                                                                     \
    if ( (idnbr > 0) && ((id < 0) || (id >= idnbr)))                    \
      {                                                                 \
        errorPrint("Id is not in solvers array range");                 \
        return MURGE_ERR_PARAMETER;                                     \
      }                                                                 \
    if (solvers == NULL)                                                \
      {                                                                 \
        errorPrint("You need to call MURGE_Initialize before");         \
        return MURGE_ERR_ORDER;                                         \
      }                                                                 \
    if (solvers[id] == NULL)                                            \
      {                                                                 \
        MURGE_MEMALLOC(solvers[id], 1, murge_data_t);                   \
                                                                        \
        solvers[id]->n           = 0;                                   \
        solvers[id]->N           = 0;                                   \
        solvers[id]->colptr      = NULL;                                \
        solvers[id]->rows        = NULL;                                \
        solvers[id]->values      = NULL;                                \
        solvers[id]->l2g         = NULL;                                \
        solvers[id]->perm        = NULL;                                \
        INIT_CENTRALISED;                                               \
        solvers[id]->b           = NULL;                                \
        solvers[id]->nrhs        = 1;                                   \
        solvers[id]->state       = MURGE_INIT_OK;                       \
        solvers[id]->pastix_data = NULL;                                \
                                                                        \
        pastix_task_init(&(solvers[id]->pastix_data), MPI_COMM_WORLD,   \
                         NULL, NULL);                                   \
                                                                        \
      }                                                                 \
  }


/*
 * Macro: CHECK_SOLVER_PARAM
 *
 * Checks if parameters have been set once for solver instance *id*.
 *
 * Checks if *iparm* or *dparm* are allocated.
 *
 * Parameters:
 *   id - Solver instance ID we want to check.
 *
 * Returns:
 *   MURGE_ERR_ORDER - If *iparm* or *dparm* are not allocated.
 */
#define CHECK_SOLVER_PARAM(id)                                          \
  {                                                                     \
    if (solvers[id]->pastix_data->iparm == NULL ||                      \
        solvers[id]->pastix_data->dparm == NULL)                        \
      {                                                                 \
        errorPrint("You need to call MURGE_SetDefaultOptions before");  \
        return MURGE_ERR_ORDER;                                         \
      }                                                                 \
  }

/*
 * Macro: CHECK_PREPROCESSING
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

#define CHECK_PREPROCESSING(id)                 \
  {                                             \
    int err = check_preprocessing(id);          \
    if (err != MURGE_SUCCESS)                   \
      return err;                               \
  }

/*
 * Macro: CHECK_FACT
 *
 * Checks if matrix values have been set.
 *
 * Checks if factorization has been performed.
 *
 * If not, it will call for it and set state.
 *
 * Parameters:
 *   id - Solver instance ID we want to check
 *
 * Returns:
 *   MURGE_ERR_ORDER - If values or right-hand-side member have not been set.
 *
 */
#define CHECK_FACT(id) do {                        \
    INTS CF_err;                                   \
    CF_err = check_fact(id);                       \
    if (CF_err != MURGE_SUCCESS)                   \
      return CF_err;                               \
  } while(0)

/*
 * Macro: CHECK_L2G
 *
 * Checks if local to global array has been allocated.
 *
 * If not, it will correct number of local columns and set local to global
 * column array.
 *
 * Parameters:
 *   id - Solver instance ID we want to check
 *
 * Returns:
 *   MURGE_ERR_ALLOCATE - If any allocation error occurs.
 *   MURGE_ERR_SOLVER   - If local to global array setting fails.
 */
#define CHECK_L2G(id)                                                   \
  do {                                                                  \
    if (solvers[id]->l2g == NULL) {                                     \
      solvers[id]->n =                                                  \
        pastix_getLocalNodeNbr(&(solvers[id]->pastix_data));            \
      if (solvers[id]->n != 0) {                                        \
        MURGE_MEMALLOC(solvers[id]->l2g, solvers[id]->n, PASTIX_INT);   \
        if (EXIT_SUCCESS !=                                             \
            (pastix_getLocalNodeLst(&(solvers[id]->pastix_data),        \
                                    solvers[id]->l2g)))                 \
          return MURGE_ERR_SOLVER;                                      \
                                                                        \
        solvers[id]->N=-1;                                              \
        cscd_build_g2l(solvers[id]->n,                                  \
                       solvers[id]->l2g,                                \
                       solvers[id]->pastix_data->pastix_comm,           \
                       &solvers[id]->N,                                 \
                       &solvers[id]->g2l);                              \
        MURGE_TRACE_MALLOC(PTR_MEMSIZE(solvers[id]->g2l), char);        \
      }                                                                 \
    }                                                                   \
  } while(0)
/*
 * Macro: CHECK_SCALING_MODE
 *
 * Check that the given mode exists.
 *
 * Parameters:
 *   mode - The given scaling mode.
 * Returns:
 *   MURGE_ERR_PARAMETER - If mode is not defined.
 */
#define CHECK_SCALING_MODE(mode)                        \
  if (mode != MURGE_SCAL_COL && mode != MURGE_SCAL_ROW) \
    {                                                   \
      errorPrint("Invalid scaling mode");               \
      return MURGE_ERR_PARAMETER;                       \
    }
/*
 * Macro: CHECK_NORM_RULE
 *
 * Check that the given norm rule exists.
 *
 * Parameters:
 *   rule - The given norm rule.
 * Returns:
 *   MURGE_ERR_PARAMETER - If rule is not defined.
 */
#define CHECK_NORM_RULE(rule)                   \
  if (rule != MURGE_NORM_MAX_COL &&             \
      rule != MURGE_NORM_2_COL   &&             \
      rule != MURGE_NORM_MAX_ROW &&             \
      rule != MURGE_NORM_2_ROW)                 \
    {                                           \
      errorPrint("Invalid  norm rule");         \
      return MURGE_ERR_PARAMETER;               \
    }
/*
 * Macro: CHOOSE_FUNC
 *
 * Will set *func* to the good function, depending to *op*.
 *
 * Parameters:
 *   op    - Operation flag.
 *   func  - Function pointer to set.
 *
 * Returns:
 *   MURGE_ERR_PARAMETER - If *op* doesn't exists.
 */
#ifdef TYPE_COMPLEX
#  define CHOOSE_FUNC(func, op)                                 \
  {                                                             \
    /* Choix de l'operation a effectuer pour les doublons */    \
    switch(op)                                                  \
      {                                                         \
      case MURGE_ASSEMBLY_ADD:                                  \
        func = &add_two_floats;                                 \
        break;                                                  \
      case MURGE_ASSEMBLY_OVW:                                  \
        func = &keep_last;                                      \
        break;                                                  \
      default:                                                  \
        func = NULL;                                            \
        return MURGE_ERR_PARAMETER;                             \
        break;                                                  \
      }                                                         \
  }
#else
#  define CHOOSE_FUNC(func, op)                                 \
  {                                                             \
    /* Choix de l'operation a effectuer pour les doublons */    \
    switch(op)                                                  \
      {                                                         \
      case MURGE_ASSEMBLY_ADD:                                  \
        func = &add_two_floats;                                 \
        break;                                                  \
      case MURGE_ASSEMBLY_OVW:                                  \
        func = &keep_last;                                      \
        break;                                                  \
      case MURGE_ASSEMBLY_MAX:                                  \
        func = &get_max;                                        \
        break;                                                  \
      case MURGE_ASSEMBLY_MIN:                                  \
        func = &get_min;                                        \
        break;                                                  \
      default:                                                  \
        func = NULL;                                            \
        return MURGE_ERR_PARAMETER;                             \
        break;                                                  \
      }                                                         \
  }
#endif

/*
 * Macro: UNLINK
 *
 * Suppress a file from the disk.
 *
 * Parameters:
 *   file - The file to suppress
 * Returns:
 *   MURGE_ERR_IO - If an error occur.
 */
#define UNLINK(file)                                    \
  if (0 != unlink(file))                                \
    {                                                   \
      errorPrint("could not unlink %s\n", file);        \
      return MURGE_ERR_IO;                              \
    }                                                   \

/*
 * Macro: LINK
 *
 * Create a symbolink link.
 *
 * Parameters:
 *   src  - file to link.
 *   dest - link path.
 * Returns:
 *   MURGE_ERR_IO - If an error occur.
 */
#define LINK(src, dest)                                                 \
  if (0 != symlink(src,dest))                                           \
    {                                                                   \
      if (errno == EEXIST)                                              \
        {                                                               \
          UNLINK(dest);                                                 \
          if (0 != symlink(src,dest))                                   \
            {                                                           \
              errorPrint("Could not link %s to %s\n", dest, src);       \
              return MURGE_ERR_IO;                                      \
            }                                                           \
        }                                                               \
      else                                                              \
        {                                                               \
          errorPrint("Could not link %s to %s\n", dest, src);           \
          return MURGE_ERR_IO;                                          \
        }                                                               \
    }
/*
 * Macro: RENAME
 *
 * Move a file on disk.
 *
 * Parameters:
 *   src  - File to move.
 *   dest - New path.
 * Returns:
 *   MURGE_ERR_IO - If an error occur.
 */
#define RENAME(src, dest)                                       \
  if (0 != rename(src, dest))                                   \
    {                                                           \
      errorPrint("couldnt rename %s into %s", src, dest);       \
      return MURGE_ERR_IO;                                      \
    }


/*
 * Macros: Time macros
 *
 * CLOCK_INIT - Start a clok
 * CLOCK_STOP - Save clock time
 * CLOCK_GET  - Get saved time (double value)
 */
#ifdef MURGE_TIME
#  define CLOCK_INIT {clockInit(&clock);clockStart(&clock);}
#  define CLOCK_STOP {clockStop(&clock);}
#  define CLOCK_GET  clockVal(&clock)
#  define CLOCK_PRINT(str) do {                                         \
    CLOCK_STOP;                                                         \
    if (solvers[id]->pastix_data->iparm[IPARM_VERBOSE] >                \
        API_VERBOSE_NO) {                                               \
      MURGE_MEMORY_USAGE_PRINT(str);                                    \
      fprintf(stdout, OUT_TIMER, str, (double)CLOCK_GET, id);           \
    }                                                                   \
  } while (0)

#  define MYCLOCK_INIT {clockInit(&myclock);clockStart(&myclock);}
#  define MYCLOCK_STOP {clockStop(&myclock);}
#  define MYCLOCK_GET  clockVal(&myclock)
#  define MYCLOCK_PRINT(str) do {                                       \
    MYCLOCK_STOP;                                                       \
    if (solvers[id]->pastix_data->iparm[IPARM_VERBOSE] >                \
        API_VERBOSE_NO)							\
      fprintf(stdout, OUT_TIMER, str, (double)MYCLOCK_GET, id);         \
  } while (0)

#else
#  define CLOCK_INIT
#  define CLOCK_STOP
#  define CLOCK_GET  0
#  define CLOCK_PRINT(str) do { } while (0)
#endif

#ifdef MURGE_DUMP_CSC
#  define MURGE_DUMP_GRAPH do {                                 \
    char name[256];                                             \
    int id_comm = solvers[id]->pastix_data->pastix_comm;        \
    MPI_Bcast(&id_comm, 1, MPI_INT, 0,                          \
              solvers[id]->pastix_data->pastix_comm);           \
    sprintf(name,"Graph_%ld_%ld_",                              \
            (long)id_comm,                                      \
            (long)solvers[id]->ndump);                          \
    solvers[id]->ndump++;                                       \
    cscd_save(solvers[id]->n,                                   \
              solvers[id]->colptr,                              \
              solvers[id]->rows,                                \
              NULL,                                             \
              NULL,                                             \
              solvers[id]->l2g,                                 \
              solvers[id]->pastix_data->iparm[IPARM_DOF_NBR],   \
              name,                                             \
              solvers[id]->pastix_data->pastix_comm);           \
  } while(0)

#  define MURGE_DUMP_MATRIX do {                                \
    char name[256];                                             \
    int id_comm = solvers[id]->pastix_data->pastix_comm;        \
    MPI_Bcast(&id_comm, 1, MPI_INT, 0,                          \
              solvers[id]->pastix_data->pastix_comm);           \
    sprintf(name,"Matrix_%ld_%ld_",                             \
            (long)id_comm,                                      \
            (long)solvers[id]->ndump);                          \
    solvers[id]->ndump++;                                       \
    cscd_save(solvers[id]->n,                                   \
              solvers[id]->colptr,                              \
              solvers[id]->rows,                                \
              solvers[id]->values,                              \
              solvers[id]->b,                                   \
              solvers[id]->l2g,                                 \
              solvers[id]->pastix_data->iparm[IPARM_DOF_NBR],   \
              name,                                             \
              solvers[id]->pastix_data->pastix_comm);           \
  } while(0)

#else  /* not MURGE_DUMP_CSC */
#  define MURGE_DUMP_GRAPH do {} while (0)
#  define MURGE_DUMP_MATRIX do {} while (0)
#endif /* not MURGE_DUMP_CSC */

#define MPI_INTL (sizeof(long) == sizeof(INTL))?MPI_LONG:MPI_INT
#define MPI_INTS (sizeof(long) == sizeof(INTS))?MPI_LONG:MPI_INT

#ifdef DISTRIBUTED
#  ifdef MURGE_TIME
#    define DPASTIX(data, comm,                                         \
                    n, colptr, rows, values, l2g, perm, invp,           \
                    b, nrhs, iparm, dparm )                             \
  do {                                                                  \
    PASTIX_INT save_veri = iparm[IPARM_MATRIX_VERIFICATION];            \
    PASTIX_INT save_free = iparm[IPARM_FREE_CSCUSER];                   \
    Clock myclock;                                                      \
    int stask = iparm[IPARM_START_TASK], etask = iparm[IPARM_END_TASK]; \
    char string[256];                                                   \
    MYCLOCK_INIT;                                                       \
    iparm[IPARM_MATRIX_VERIFICATION] = API_NO;                          \
    iparm[IPARM_FREE_CSCUSER] = API_CSC_PRESERVE;                       \
    dpastix(data, comm,                                                 \
            n, colptr, rows, values, l2g, perm, invp,                   \
            b, nrhs, iparm, dparm );                                    \
    sprintf(string, OUT_CALL_DPASTIX, stask, etask);                    \
    MYCLOCK_PRINT(string);                                              \
    iparm[IPARM_MATRIX_VERIFICATION] = save_veri;                       \
    iparm[IPARM_FREE_CSCUSER] = save_free;                              \
  } while (0)
#  else /* MURGE_TIME */
#    define DPASTIX(data, comm,                                         \
                    n, colptr, rows, values, l2g, perm, invp,           \
                    b, nrhs, iparm, dparm )                             \
  do {                                                                  \
    PASTIX_INT save_veri = iparm[IPARM_MATRIX_VERIFICATION];            \
    PASTIX_INT save_free = iparm[IPARM_FREE_CSCUSER];                   \
    iparm[IPARM_MATRIX_VERIFICATION] = API_NO;                          \
    iparm[IPARM_FREE_CSCUSER] = API_CSC_PRESERVE;                       \
    dpastix(data, comm,                                                 \
            n, colptr, rows, values, l2g, perm, invp,                   \
            b, nrhs, iparm, dparm );                                    \
    iparm[IPARM_MATRIX_VERIFICATION] = save_veri;                       \
    iparm[IPARM_FREE_CSCUSER] = save_free;                              \
  } while (0)
#  endif /* MURGE_TIME */
#  define PASTIX_FILLIN_CSC(data, comm,                                 \
                            n, colptr, rows, values, b, nrhs, l2g)      \
    do {                                                                \
      PASTIX_INT save_free = iparm[IPARM_FREE_CSCUSER];                 \
      iparm[IPARM_FREE_CSCUSER] = API_CSC_PRESERVE;                     \
      pastix_fillin_csc(data, comm,                                     \
                        n, colptr, rows, values, b, nrhs, l2g);         \
      iparm[IPARM_FREE_CSCUSER] = save_free;                            \
    } while(0)
#else /* DISTRIBUTED */
static inline
void murge_dpastix(INTS id,
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
                   double         * dparm);
#  ifdef MURGE_TIME
#    define DPASTIX(data, comm, n, colptr, rows, values,                \
                    l2g, perm, invp, b, nrhs, iparm, dparm )            \
  do {                                                                  \
    Clock myclock;                                                      \
    int stask = iparm[IPARM_START_TASK], etask = iparm[IPARM_END_TASK]; \
    char string[256];                                                   \
    MYCLOCK_INIT;                                                       \
    murge_dpastix(id, data, comm, n, colptr, rows, values,              \
                  l2g, perm, invp, b, nrhs, iparm, dparm );             \
    sprintf(string, OUT_CALL_PASTIX, stask, etask);                     \
    MYCLOCK_PRINT(string);                                              \
  } while(0)
#  else
#    define DPASTIX(data, comm, n, colptr, rows, values,                \
                    l2g, perm, invp, b, nrhs, iparm, dparm )            \
      do {                                                              \
        murge_dpastix(id, data, comm, n, colptr, rows, values,          \
                      l2g, perm, invp, b, nrhs, iparm, dparm );         \
      } while(0)
#  endif
static inline
void murge_pastix_fillin_csc(INTS id,
                             pastix_data_t * data,
                             MPI_Comm        comm,
                             PASTIX_INT      n,
                             PASTIX_INT    * colptr,
                             PASTIX_INT    * rows,
                             PASTIX_FLOAT  * values,
                             PASTIX_FLOAT  * b,
                             PASTIX_INT      nrhs,
                             PASTIX_INT    * l2g);
#  define PASTIX_FILLIN_CSC(data, comm,                                 \
                            n, colptr, rows, values, b, nrhs, l2g)      \
  do {                                                                  \
    PASTIX_INT save_free = iparm[IPARM_FREE_CSCUSER];                   \
    iparm[IPARM_FREE_CSCUSER] = API_CSC_PRESERVE;                       \
    murge_pastix_fillin_csc(id, data, comm,                             \
                            n, colptr, rows, values, b, nrhs, l2g);     \
    iparm[IPARM_FREE_CSCUSER] = save_free;                              \
  } while(0)

#define HASH_FIND_INTL(head,findint,out)        \
  HASH_FIND(hh,head,findint,sizeof(INTL),out)
#define HASH_ADD_INTL(head,intfield,add)        \
  HASH_ADD(hh,head,intfield,sizeof(INTL),add)
#define HASH_REPLACE_INTL(head,intfield,add,replaced)    \
  HASH_REPLACE(hh,head,intfield,sizeof(INTL),add,replaced)

#  define HASH_MATRIX_ADD_NODE(id, ROW, COL, VALUES, op) \
  do {                                                   \
    struct hash_mtx_entry_s * s;                         \
    struct hash_mtx_entry_s add;                         \
    INTS dof = solvers[id]->iparm[IPARM_DOF_NBR];        \
    add.key = (COL-1)*solvers[id]->N+ROW-1;              \
    add.val = val;                                       \
    HASH_FIND_INTL(solvers[id]->hashmtx, add.key, s);    \
    if (s != NULL)                                       \
      {                                                  \
        INTS i;                                          \
        /* add using given operation */                  \
        for ( i = 0; i < dof*dof; i++)                   \
          s.val[i] = op(s.val[i], add.val[i]);           \
      }                                                  \
    else                                                 \
      {                                                  \
        HASH_ADD_INTL(solvers[id]->hashmtx, add);        \
      }                                                  \
  } while(0)

#  define HASH_MATRIX_ADD(id, ROW, COL, VALUES, op)      \
  do {                                                   \
    struct hash_mtx_entry_s * s;                         \
    struct hash_mtx_entry_s add;                         \
    INTS dof = solvers[id]->iparm[IPARM_DOF_NBR];        \
    INTS NODE_ROW = (ROW-1)/dof;                         \
    INTS NODE_COL = (COL-1)/dof;                         \
    add.key = NODE_COL*solvers[id]->N+NODE_ROW;          \
    add.val = val;                                       \
    HASH_FIND_INTL(solvers[id]->hashmtx, add.key, s);    \
    if (s != NULL) {                                     \
        INTS i = (COL-1)%dof*dof+(ROW-1)%dof;            \
        /* add using given operation */                  \
        s.val[i] = op(s.val[i], add.val[i]);             \
    }                                                    \
    else {                                               \
      HASH_ADD_INTL(solvers[id]->hashmtx, add);          \
    }                                                    \
  } while(0)

#endif /* DISTRIBUTED */
#endif /* MURGE_DEFINE_H */
