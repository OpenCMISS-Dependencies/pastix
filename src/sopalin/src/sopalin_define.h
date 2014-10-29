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
 * file: sopalin_define.h
 *
 * define maccros for sopalin
 */
#ifndef SOPALIN_DEFINE_H
#define SOPALIN_DEFINE_H
/*
 * Enum: TAG_COMM
 *
 * MPI tags
 *
 *   TAG_FANIN - Tags for fanin buffers
 *   TAG_BLOCK - Tags for block buffers
 *   TAG_DOWN  - Tags for down step buffers
 *   TAG_UP    - Tags for up step buffers
 */
typedef enum TAG_COMM {
  TAG_FANIN = 1,
  TAG_BLOCK = 2,
  TAG_DOWN  = 3,
  TAG_UP    = 4,
  TAG_END   = 5
} TagComm_t;

/*
 * Enum: COMMSTEP
 *
 * Communication steps.
 *
 *  COMMSTEP_END       - All communications are ended
 *  COMMSTEP_INIT      - Signal the return of communication thread in
 *                       initial state
 *  COMMSTEP_FACTO     - Ask to begin the factorization step
 *  COMMSTEP_FACTOEND  - Signal the end of communications in factorization step
 *  COMMSTEP_DOWN      - Ask to begin the down step
 *  COMMSTEP_UP        - Ask to begin the up step
 *  COMMSTEP_UPDOEND   - Signal the end of different steps during
 *                       backward/forward substitution
 *  COMMSTEP_ALLREDUCE - Ask for an MPI_ALLREDUCE in Funneled mode
 *  COMMSTEP_REDUCE    - Ask for an MPI_REDUCE in Funneled mode
 */
typedef enum COMMSTEP {
  COMMSTEP_END       = 0,
  COMMSTEP_INIT      = 1,
  COMMSTEP_FACTO     = 2,
  COMMSTEP_FACTOEND  = 3,
  COMMSTEP_DOWN      = 4,
  COMMSTEP_UP        = 5,
  COMMSTEP_UPDOEND   = 6,
  COMMSTEP_ALLREDUCE = 7,
  COMMSTEP_REDUCE    = 8
} CommStep_t;

#if (defined TRACE_SOPALIN) || (defined PASTIX_DUMP_SOLV_COMM)
#  define UPDOWN_SIZETAB    8 /* 3 -> 4 pour alignement */
#else
#  define UPDOWN_SIZETAB    4 /* 3 -> 4 pour alignement */
#endif
#define MAX_SMX_REQUESTS 32 /* < MAX_S_REQUESTS */

#define DUMP_CSC  0x1
#define DUMP_SOLV 0x2
#define DUMP_SMB  0x4

/*
 * Defines: Init parameters
 *
 * Values used to init coefficcients
 *
 *   ZERO - 0.0
 *   UN   - 1.0
 *   DEUX - 2.0
 */
#ifdef CPLX
#  define ZERO (0.0+0.0*I)
#  define UN   (1.0+0.0*I)
#  define DEUX (2.0+0.0*I)
#else
#  define ZERO 0.0
#  define UN   1.0
#  define DEUX 2.0
#endif

/**************************************/
/*           OPTIONS DE COMPIL        */
/**************************************/

/*
 * Defines: Communications
 *
 * EXACT_TAG     - Do not use MPI_ANY_TAG
 * EXACT_THREAD  - The tag used is the receiving thread,
 *                 otherwise the task number ( limited by the number of
 *                 tag allowed by MPI)
 */
/* #define EXACT_TAG */
#define EXACT_THREAD

/*
 Define: COMPUTE
 Effective computaiton or not.
 */
#define COMPUTE

/*
 Define: STATS_SOPALIN
 Statistics about memory used.
 */
/* #define STATS_SOPALIN */

/*
  Define: TRYLOCK
  Compute the number of pthread_mutex_lock et cond_wait "missed".
*/

/* #define TRYLOCK */

/*
 Define: USE_CSC
 TODO : fill-in
 */
#define USE_CSC

/*
 Define: MULT_SMX
 multi-right-hand-side method
 */
/*#define MULT_SMX*/

/* TODO : a completer */
/*#define DEADCODE*/

/* TODO : a completer */
/*#define GRAD_CONJ*/

/* TODO : a completer */
#define GMRES

/* TODO : a completer */
#define PRECOND

/* TODO : (desactiver avec ILU_EXTERN de main_sopalin3d.c) */
#define FACTORIZE

/* Pour facto incomplete, blocs ne correspondent pas tout a fait,
 * il faut calculer la zone de contribution */
#define NAPA_SOPALIN

/* Pour ajouter le support des threads dans la factorisation
 * et la descente remontée.
 * Desactive si FORCE_NOSMP
 */
#define SMP_SOPALIN

/*
  Define: SMP_RAFF

  Add multi-threading in reffinement steps.

 Deactivated if FORCE_NOSMP ou NOSMP_RAFF is defined.
 */
#define SMP_RAFF


/* TODO : a completer */
/*#define TRACE_PICL*/

/* Utilisation d'envoi (par thread) non bloquant */
#define TEST_ISEND

/* Utilisation de réception (par thread) non bloquante */
/*#define TEST_IRECV*/
#ifdef PASTIX_DYNSCHED
#  undef TEST_IRECV
#endif
/* TODO : a completer */
/*#define RECV_FANIN_OR_BLOCK*/

/* Forcer la consommation des réception en attente */
/*#define FORCE_CONSO*/
#ifdef PASTIX_DYNSCHED
#  undef FORCE_CONSO
#endif
/* Allocation des fanin TODO : a completer */
#define ALLOC_FTGT

/* TODO : a completer */
/*#define HPM_SOPALIN*/

/*
 * STORAGE : stocke le blco-colonne pour faire un blas global,
 * plutot qu'un blas par bloc.
 * STORAGE va couter cher en multi-seconds-membres ?
 */
#define STORAGE

/*
 * DEP_SMX = utilisation d'une file de job prêt. (efficient ? (no))
 * WARNING DEP_SMX don't work
 */
/*#define DEP_SMX*/

/*#define UPDO_DEADCODE*/

/*#define TRACE_SOPALIN*/
/*
 #define COMPUTE_ALLOC
 */

/*
 * If you want to activate the old way of stealing task by following the bubble tree
 */
/* #define PASTIX_DYNSCHED_WITH_TREE */

/**************************************/
/*        Flag dependencies           */
/**************************************/

/* Suppression des options liées au SMP */
#if (defined FORCE_NOSMP)
#  undef SMP_SOPALIN
#  undef SMP_RAFF
#  undef TRYLOCK
#  undef PASTIX_DYNSCHED
#endif

#ifdef NOSMP_RAFF
#  undef SMP_RAFF
#endif

/* Suppression des options liées à MPI */
#if (defined FORCE_NOMPI)
#  undef PASTIX_UPDO_ISEND
#endif


#if (!(defined SOLVER_UPDOWN))
#  define SOLVER_UPDOWN
#endif

/*#undef SOLVER_UPDOWN*/

/*
 * La version a bulles fonctionne uniquement avec l'option thread comm
 * en raison des flags utilisés dans les comms et du nombre de threads
 * inconnus
 */
#if (defined PASTIX_DYNSCHED) && (!(defined COMM_REORDER))
#  define COMM_REORDER
#endif

/* Ne peut pas fonctionner a cause de la répartition
   des taches dans les superieurs */
#if (defined PASTIX_DYNSCHED) && (defined SMP_RAFF)
#  undef SMP_RAFF
#endif

/* COMPUTE_ALLOC need ALLOC_FTGT */
#if (defined COMPUTE_ALLOC) && (!(defined ALLOC_FTGT))
#  define ALLOC_FTGT
#endif

/*
 * EXACT_TAG not compatible with RECV_FANIN_OR_BLOCK or FORCE_CONSO
 * TODO : cf si il faut aussi les desactiver avec EXACT_THREAD
 */
#if (defined EXACT_TAG)
#  undef RECV_FANIN_OR_BLOCK
#  undef FORCE_CONSO
#endif

/* SOLVER_UPDOWN need ALLOC_FTGT   */
#if (defined SOLVER_UPDOWN) && (!(defined ALLOC_FTGT))
#  define ALLOC_FTGT
#endif

/* SMP_SOPALIN with not THREAD_COMM needs EXACT_TAG or EXACT_THREAD ??? */
/* NEED CHECK => THREAD_COMM does not exists anymore ! */
#if (defined SMP_SOPALIN) && (!(defined THREAD_COMM))
#  if   (!(defined EXACT_TAG)) && (!(defined EXACT_THREAD))
#    define EXACT_THREAD
#  endif
#endif

#if (defined HERMITIAN && defined CHOL_SOPALIN)
#error "HERMITIAN and LLT/LU are not compatible"
#endif

/* SOPALIN_LU need CHOL_SOPALIN    */
#if (defined SOPALIN_LU) && (!(defined CHOL_SOPALIN))
#  define CHOL_SOPALIN
#endif


/* SMP_SOPALIN not compatible with DEP_SMX */
#if (defined SMP_SOPALIN) && (defined DEP_SMX)
#  undef DEP_SMX
#endif

/* DEP_SMX not compatible with EXACT_TAG */
#if (defined DEP_SMX) && (defined EXACT_TAG)
#  undef EXACT_TAG
#endif

#if (defined TEST_IRECV) && (!(defined FORCE_CONSO))
#  define FORCE_CONSO
#endif

#if (defined STATS_SOPALIN) && (!(defined ALLOC_FTGT))
#  undef STATS_SOPALIN
#endif


#if (defined TRACE_SOPALIN) && (defined PASTIX_DUMP_SOLV_COMM)
#  error "TRACE_SOPALIN and PASTIX_DUMP_SOLV_COMM can't be both activated"
#endif

#ifdef OOC
#  ifdef PASTIX_ESC /* PASTIX_ESC not compatible with OOC */
#    undef PASTIX_ESC
#  endif
#  ifdef PASTIX_DYNSCHED
#    undef PASTIX_DYNSCHED
#  endif
#endif /* OOC */

#ifdef OOC_FTGT
#  ifndef FORCE_CONSO /* OOC works better with FORCE_CONSO */
#    define FORCE_CONSO
#  endif
#  ifndef COMM_REORDER /* OOC needs COMM_REORDER to remove FTGT as soon as possible */
#    define COMM_REORDER
#  endif
#  ifndef ALLOC_FTGT   /* OOC needs ALLOC_FTGT to have pointers on each FTGT */
#    define ALLOC_FTGT
#  endif
#endif /* OOC_FTGT */

#if (defined PASTIX_ISEND) && (!defined PASTIX_UPDO_ISEND)
#  define PASTIX_UPDO_ISEND
#endif

#if (defined PASTIX_DYNSCHED && !defined PASTIX_UPDO_ISEND)
#  define PASTIX_UPDO_ISEND
#endif

/* Storage est trop couteux en multi-rhs */
#if (defined MULT_SMX) && (defined STORAGE)
#  undef STORAGE
#endif

/* On re-active storage en PASTIX_DYNSCHED qui implique thread comm meme en multi-rhs */
#if (defined PASTIX_DYNSCHED) && !(defined STORAGE)
#  define STORAGE
#endif

#ifdef PASTIX_THREAD_COMM
#  undef TEST_IRECV
#endif

#ifdef PASTIX_FUNNELED
#  ifndef PASTIX_THREAD_COMM
#    define PASTIX_THREAD_COMM
#  endif
#  ifndef COMM_REORDER
#    define COMM_REORDER
#  endif
#endif

#ifdef PASTIX_THREAD_COMM
#  ifndef PASTIX_UPDO_ISEND
#    define PASTIX_UPDO_ISEND
#  endif
#  ifndef STORAGE
#    define STORAGE
#  endif
#  ifdef FORCE_CONSO
#    undef FORCE_CONSO
#  endif
#  ifdef FORCE_NOSMP
#    error "FORCE_NOSMP is not compatible with PASTIX_THREAD_COMM nor PASTIX_FUNNELED"
#  endif
#  ifdef FORCE_NOMPI
#    error "FORCE_NOMPI is not compatible with PASTIX_THREAD_COMM nor PASTIX_FUNNELED"
#  endif
#  ifdef OOC
#    error "OOC is not compatible with PASTIX_THREAD_COMM nor PASTIX_FUNNELED"
#  endif
#endif


/*
 * Protection des appels MPI
 */

#ifdef VERIF_MPI
extern int err_mpi;
#  define MPI_PRINT_ERR(x,y) {\
                char s[MPI_MAX_ERROR_STRING]; int l;\
    MPI_Error_string(y,s,&l);\
    fprintf(stderr, "error in %s (%s)(line=%d,file=%s)\n",x,s,__LINE__,__FILE__); \
    EXIT(MOD_SOPALIN,UNKNOWN_ERR);}
#  define CALL_MPI    err_mpi=
#  define TEST_MPI(x) if (err_mpi!=MPI_SUCCESS) MPI_PRINT_ERR(x,err_mpi)
#else
#  define CALL_MPI
#  define TEST_MPI(x)
#endif


#ifdef USE_CSC
/* #define CSC_SOPALIN_HACK */
#endif

#if (defined STARPU_CONTEXT && defined STARPU_PASTIX_SCHED)
#  error "STARPU_CONTEXT and STARPU_PASTIX_SCHED ar not compatible"
#endif
/*
 macro: API_CALL
 Add prefixe (U (LU),L (LLt) or D (LDLt)) to sopalin functions.
 */
#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU
#    define API_CALL(nom) PASTIX_PREFIX_F(ge_ ## nom)
#  else
#    define API_CALL(nom) PASTIX_PREFIX_F(po_ ## nom)
#  endif
#else
#  ifdef HERMITIAN
#    define API_CALL(nom) PASTIX_PREFIX_F(he_ ## nom)
#  else
#    define API_CALL(nom) PASTIX_PREFIX_F(sy_ ## nom)
#  endif
#endif

/************************************************************/
/*   Statistiques to compute memory overhead                */
/************************************************************/
#ifdef STATS_SOPALIN
#  define STATS_ADD(siz)                                        \
  MUTEX_LOCK(&(sopalin_data->mutex_alloc));                     \
  sopalin_data->current_alloc += siz;                           \
  if ( sopalin_data->current_alloc > sopalin_data->max_alloc)   \
    sopalin_data->max_alloc = sopalin_data->current_alloc;      \
  MUTEX_UNLOCK(&(sopalin_data->mutex_alloc));

#  define STATS_SUB(siz)                                        \
  MUTEX_LOCK(&(sopalin_data->mutex_alloc));                     \
  sopalin_data->current_alloc -= siz;                           \
  MUTEX_UNLOCK(&(sopalin_data->mutex_alloc));
#else
#  define STATS_ADD(siz)
#  define STATS_SUB(siz)
#endif

/************************************************************/
/*   Redefinition de MPI_Allreduce pour thread_funneled     */
/************************************************************/

/*
 Macro: MyMPI_Allreduce

 Redefine MPI_Allreduce for Funneled mode.

  Parameters:
    lsendbuf  - Send buffer.
    lrecvbuf  - Receiving buffer.
    lcount    - Number of elements.
    ldatatype - Type of elements.
    lop       - Reduction operation.
    lcomm     - MPI communicator.
 */
#define MyMPI_Allreduce(lsendbuf, lrecvbuf, lcount, ldatatype, lop, lcomm) \
  {                                                                     \
    if (sopalin_data->sopar->iparm[IPARM_THREAD_COMM_MODE] &            \
        API_THREAD_FUNNELED)                                            \
      {                                                                 \
        Pastix_Allreduce_t *allreduce = &(sopalin_data->allreduce);     \
        allreduce->sendbuf  = lsendbuf;                                 \
        allreduce->recvbuf  = lrecvbuf;                                 \
        allreduce->count    = lcount;                                   \
        allreduce->datatype = ldatatype;                                \
        allreduce->op       = lop;                                      \
        /* Signal d'un reduce a faire */                                \
        MUTEX_LOCK(&(sopalin_data->mutex_comm));                        \
        sopalin_data->step_comm = COMMSTEP_ALLREDUCE;                   \
        print_debug(DBG_THCOMM, "%s:%d ALLREDUCE\n",                    \
                    __FILE__, __LINE__);                                \
        MUTEX_UNLOCK(&(sopalin_data->mutex_comm));                      \
        pthread_cond_broadcast(&(sopalin_data->cond_comm));             \
        /* Attente de la fin du reduce */                               \
        MUTEX_LOCK(&(sopalin_data->mutex_comm));                        \
        while(sopalin_data->step_comm != COMMSTEP_INIT)                 \
          COND_WAIT(&(sopalin_data->cond_comm),                         \
                    &(sopalin_data->mutex_comm));                       \
        MUTEX_UNLOCK(&(sopalin_data->mutex_comm));                      \
      }                                                                 \
    else                                                                \
      {                                                                 \
        MPI_Allreduce(lsendbuf, lrecvbuf, lcount, ldatatype, lop,       \
                      lcomm);                                           \
      }                                                                 \
  }
/*
  Macro: MyMPI_Reduce

  Redefine MPI_Reduce for Funneled mode.

  Parameters:
    lsendbuf  - Send buffer.
    lrecvbuf  - Receiving buffer.
    lcount    - Number of elements.
    ldatatype - Type of elements.
    lop       - Reduction operation.
    lroot     - Root of the operation.
    lcomm     - MPI communicator.
 */
#define MyMPI_Reduce(lsendbuf, lrecvbuf, lcount, ldatatype, lop, lroot, lcomm) \
  {                                                                     \
    if (sopalin_data->sopar->iparm[IPARM_THREAD_COMM_MODE] &            \
        API_THREAD_FUNNELED)                                            \
      {                                                                 \
        Pastix_Allreduce_t *reduce = &(sopalin_data->allreduce);        \
        reduce->sendbuf  = lsendbuf;                                    \
        reduce->recvbuf  = lrecvbuf;                                    \
        reduce->count    = lcount;                                      \
        reduce->datatype = ldatatype;                                   \
        reduce->op       = lop;                                         \
        /* Signal d'un reduce a faire */                                \
        MUTEX_LOCK(&(sopalin_data->mutex_comm));                        \
        sopalin_data->step_comm = COMMSTEP_REDUCE;                      \
        print_debug(DBG_THCOMM, "%s:%d REDUCE\n", __FILE__, __LINE__);  \
        MUTEX_UNLOCK(&(sopalin_data->mutex_comm));                      \
        pthread_cond_broadcast(&(sopalin_data->cond_comm));             \
        /* Attente de la fin du reduce */                               \
        MUTEX_LOCK(&(sopalin_data->mutex_comm));                        \
        while(sopalin_data->step_comm != COMMSTEP_INIT)                 \
          COND_WAIT(&(sopalin_data->cond_comm),                         \
                    &(sopalin_data->mutex_comm));                       \
        MUTEX_UNLOCK(&(sopalin_data->mutex_comm));                      \
      }                                                                 \
    else                                                                \
      {                                                                 \
        MPI_Reduce(lsendbuf, lrecvbuf, lcount, ldatatype, lop,          \
                   lroot, lcomm);                                       \
      }                                                                 \
  }
#endif
