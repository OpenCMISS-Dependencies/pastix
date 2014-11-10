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
 * File: pastix.c
 *
 * PaStiX external functions implementations.
 *
 * Authors:
 *   Mathieu FAVERGE  - faverge@labri.fr
 *   Xavier  LACOSTE  - lacoste@labri.fr
 *   Pierre  RAMET    - ramet@labri.fr
 */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef MARCEL
#  include <pthread.h>
#endif
#ifdef WITH_SEM_BARRIER
#  include <semaphore.h>
#  include <fcntl.h>
#  include <sys/ipc.h>
#  include <sys/shm.h>
#endif
#include <sys/stat.h>

#ifdef FORCE_NOMPI
#  include "nompi.h"
#else
#  include <mpi.h>
#endif

#ifdef METIS
#  include "metis.h" /* must be before common_pastix */
#  ifdef ASSERT
#    undef ASSERT
#  endif
#endif

#include "common_pastix.h"
#include "tools.h"
#include "sopalin_define.h"
#ifdef WITH_STARPU
#include "compute_context_nbr.h"
#endif

#ifdef WITH_SCOTCH
#  ifdef    DISTRIBUTED
#    include <ptscotch.h>
#  else
#    include <scotch.h>
#  endif /* DISTRIBUTED */
#endif /* WITH_SCOTCH */


#include "dof.h"
#include "ftgt.h"
#include "symbol.h"
#include "csc.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
#include "assembly.h"
#include "param_blend.h"
#include "order.h"
#include "fax.h"
#include "kass.h"
#include "blend.h"
#include "solverRealloc.h"
#include "sopalin_thread.h"
#include "stack.h"
#include "sopalin3d.h"
#include "sopalin_init.h"
#include "sopalin_option.h"
#include "csc_intern_updown.h"
#include "csc_intern_build.h"
#include "coefinit.h"
#include "out.h"
#include "pastix.h"
#include "pastix_internal.h"
#include "pastixstr.h"

#include "csc_utils.h"
#include "cscd_utils.h"
#include "cscd_utils_intern.h"
#include "bordi.h"
#include "sopalin_acces.h"
#include "scaling.c"
#include "perf.h"

/*******************************************************************************
 * Section: Defines and Macros
 */

/*
 * Macro: FORTRAN_CALL
 *
 * Call a fortran function.
 *
 * Parameters:
 * name - Fortran function name.
 *
 */
#ifndef USE_NOFORTRAN
#  if (defined X_ARCHpower_ibm_aix)
#    define FORTRAN_CALL(name) name
#  else
#    define FORTRAN_CALL(name) name ## _
#  endif
#else
#  define FORTRAN_CALL(name)
#endif

#undef THREAD_FUNNELED_ON
#undef THREAD_FUNNELED_OFF
#define THREAD_FUNNELED_ON (                                    \
    sopar->iparm[IPARM_THREAD_COMM_MODE] &                      \
    API_THREAD_FUNNELED)
#define THREAD_FUNNELED_OFF (!THREAD_FUNNELED_ON)

#undef THREAD_COMM_ON
#undef THREAD_COMM_OFF
#define THREAD_COMM_ON  (                                        \
    sopar->iparm[IPARM_THREAD_COMM_MODE] &                       \
    ( API_THREAD_FUNNELED|API_THREAD_COMM_ONE|                   \
      API_THREAD_COMM_DEFINED|API_THREAD_COMM_NBPROC ) )
#define THREAD_COMM_OFF (!THREAD_COMM_ON)

#define OPEN_SEM(sem, name, value) do {                              \
    sem = sem_open(name, O_CREAT, S_IRUSR | S_IWUSR, value);        \
    if (sem == SEM_FAILED)                                           \
      {                                                              \
        errorPrint("sem_open failed\n");                             \
        perror("sem_open");                                          \
      }                                                              \
  } while(0)



/*
 * Defines: pastix.c defines
 *
 *   PASTIX_LOG            - If defined, start and end of this file
 *                           functions will be printed on stdout.
 *   COMPUTE               - If not defined, PaStiX will not use user's
 *                           coefficient.
 *   FORGET_PARTITION      - If defined, PaStiX won't use Scotch partition.
 *   DUMP_SYMBOLMATRIX     - Write the symbol matrix in a postscript format.
 *   STR_SIZE              - The default size of a string.
 *   TAG_RHS               - MPI tag used to communicate right-hand-side member.
 *   SCOTCH_STRAT_DIRECT   - Default Scotch strategy for the direct solver.
 *   SCOTCH_STRAT_INCOMP   - Default Scotch strategy for the incomplete solver.
 *   SCOTCH_STRAT_PERSO    - Parametrisable Scotch strategy for the direct
 *                           solver, can be set using several <IPARM_ACCESS>.
 *   PTSCOTCH_STRAT_DIRECT - Default PT-Scotch strategy for the direct solver.
 *   PTSCOTCH_STRAT_INCOMP - Default PT-Scotch strategy for the incomplete
 *                           solver.
 *   PTSCOTCH_STRAT_PERSO  - Parametrisable PT-Scotch strategy for the direct
 *                           solver, can be set using several <IPARM_ACCESS>.
 */
/*#define PASTIX_LOG*/
#define COMPUTE
/* #define FORGET_PARTITION  */
/* #define DUMP_SYMBOLMATRIX */

#define STR_SIZE               60
#define TAG_RHS                7
#define TAG_RHS2               8
/* define pour l'affichage */
#define SCOTCH_STRAT_DIRECT                                             \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{rat=0.8,"                                  \
  ""                        "vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{rat=0.8,"                                  \
  ""                        "vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.0},"                       \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{rat=0.8,"                                 \
  ""                         "vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                        "m{rat=0.8,"                                \
  ""                          "vert=100,"                               \
  ""                          "low=h{pass=10},"                         \
  ""                          "asc=f{bal=0.2}};,"                       \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"

#define SCOTCH_STRAT_INCOMP                                             \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.08},"                        \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                       "m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""            "ose=g}}"
#define SCOTCH_STRAT_PERSO                                              \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>%ld)?m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=%ld,cmax=%ld,frat=%f},"                           \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>%ld)?(m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                       "m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=%ld,cmax=%ld,frat=%f},"                           \
  ""      "ose=g}}"

#define PTSCOTCH_STRAT_DIRECT                                           \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{rat=0.8,"                                 \
  ""                        "vert=100,"                                \
  ""                        "low=h{pass=10},"                          \
  ""                        "asc=f{bal=0.2}}|"                         \
  ""                      "m{rat=0.8,"                                 \
  ""                        "vert=100,"                                \
  ""                        "low=h{pass=10},"                          \
  ""                        "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.0},"                         \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{type=h,"                                  \
  ""                         "rat=0.8,"                                 \
  ""                         "vert=100000,"                             \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=08.2}})|"                       \
  ""                       "m{type=h,"                                  \
  ""                         "rat=0.8,"                                 \
  ""                         "vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"
#define PTSCOTCH_STRAT_INCOMP                                           \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.08},"                        \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                       "m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"
#define PTSCOTCH_STRAT_PERSO  "c{rat=0.7,cpr=n{sep=/(vert>%ld)?m{vert=100,low=h{pass=10},asc=f{bal=0.2}}|m{vert=100,low=h{pass=10},asc=f{bal=0.2}};,ole=f{cmin=%ld,cmax=%ld,frat=%f},ose=g},unc=n{sep=/(vert>%ld)?(m{vert=100,low=h{pass=10},asc=f{bal=0.2}})|m{vert=100,low=h{pass=10},asc=f{bal=0.2}};,ole=f{cmin=%ld,cmax=%ld,frat=%f},ose=g}}"


/*******************************************************************************
 *  Section: Macros
 */

/*
  macro: print_onempi

  Print a string using processor 0.
  Uses printf syntax.

  Parameters:
  fmt - Format string (see printf manual).
  ... - Arguments depending on the format string.
*/
#define print_onempi(fmt, ...) if(procnum == 0) fprintf(stdout, fmt, __VA_ARGS__)

/*******************************************************************************
 * Section: Functions
 */

/*
 * Function: pastix_initParam
 *
 * sets default parameters for iparm and dparm
 *
 * Parameters:
 * iparm - tabular of IPARM_SIZE integer parameters.
 * dparm - tabular of DPARM_SIZE double parameters.
 */
void pastix_initParam(PASTIX_INT    *iparm,
                      double *dparm)
{
#ifdef PASTIX_LOG
  fprintf(stdout, "-> api_init\n");
#endif

  PASTIX_INT i;
  for (i=0; i<IPARM_SIZE; i++)
    iparm[i] = 0;

  for (i=0; i<DPARM_SIZE; i++)
    dparm[i] = 0;

  iparm[IPARM_MODIFY_PARAMETER]      = API_YES;             /* Indicate if parameters have been set by user         */
  iparm[IPARM_START_TASK]            = API_TASK_ORDERING;   /* Indicate the first step to execute (see PaStiX steps)*/
  iparm[IPARM_END_TASK]              = API_TASK_CLEAN;      /* Indicate the last step to execute (see PaStiX steps) */
  iparm[IPARM_VERBOSE]               = API_VERBOSE_NO;      /* Verbose mode (see Verbose modes)                     */
  iparm[IPARM_DOF_NBR]               = 1;                   /* Degree of freedom per node                           */
  iparm[IPARM_DOF_COST]              = 0;                   /* Degree of freedom for cost computation
                                                               (If different from IPARM_DOF_NBR) */
  iparm[IPARM_ITERMAX]               = 250;                 /* Maximum iteration number for refinement              */
  iparm[IPARM_MATRIX_VERIFICATION]   = API_YES;             /* Check the input matrix                               */
  iparm[IPARM_MC64]                  = 0;                   /* MC64 operation <pastix.h> IGNORE                     */
  iparm[IPARM_ONLY_RAFF]             = API_NO;              /* Refinement only                                      */
  iparm[IPARM_TRACEFMT]              = API_TRACE_PICL;      /* Trace format (see Trace modes)                       */
  iparm[IPARM_GRAPHDIST]             = API_YES;             /* Specify if the given graph is distributed or not     */
  iparm[IPARM_AMALGAMATION_LEVEL]    = 5;                   /* Amalgamation level                                   */
  iparm[IPARM_ORDERING]              = API_ORDER_SCOTCH;    /* Choose ordering                                      */
  iparm[IPARM_DEFAULT_ORDERING]      = API_YES;             /* Use default ordering parameters with scotch or metis */
  iparm[IPARM_ORDERING_SWITCH_LEVEL] = 120;                 /* Ordering switch level    (see Scotch User's Guide)   */
  iparm[IPARM_ORDERING_CMIN]         = 0;                   /* Ordering cmin parameter  (see Scotch User's Guide)   */
  iparm[IPARM_ORDERING_CMAX]         = 100000;              /* Ordering cmax parameter  (see Scotch User's Guide)   */
  iparm[IPARM_ORDERING_FRAT]         = 8;                   /* Ordering frat parameter  (see Scotch User's Guide)   */
  iparm[IPARM_STATIC_PIVOTING]       = 0;                   /* number of control of diagonal magnitude              */
  iparm[IPARM_METIS_PFACTOR]         = 0;                   /* Metis pfactor                                        */
  iparm[IPARM_NNZEROS]               = 0;                   /* memory space for coefficients                        */
  iparm[IPARM_ALLOCATED_TERMS]       = 0;                   /* number of non zero in factorized sparse matrix       */
  iparm[IPARM_MIN_BLOCKSIZE]         = 60;                  /* min blocksize                                        */
  iparm[IPARM_MAX_BLOCKSIZE]         = 120;                 /* max blocksize                                        */
  iparm[IPARM_SCHUR]                 = API_NO;              /* Schur mode */
  iparm[IPARM_ISOLATE_ZEROS]         = API_NO;              /* Isolate null diagonal terms at the end of the matrix */
  iparm[IPARM_FACTORIZATION]         = API_FACT_LDLT;       /* LdLt     */
  iparm[IPARM_CPU_BY_NODE]           = 0;                   /* cpu/node */
  iparm[IPARM_BINDTHRD]              = API_BIND_AUTO;       /* Default binding method */
  iparm[IPARM_THREAD_NBR]            = 1;                   /* thread/mpi */
  iparm[IPARM_CUDA_NBR]              = 0;                   /* CUDA devices */
  iparm[IPARM_DISTRIBUTION_LEVEL]    = 0;                   /* 1d / 2d */
  iparm[IPARM_LEVEL_OF_FILL]         = 1;                   /* level of fill */
  iparm[IPARM_IO_STRATEGY]           = API_IO_NO;           /* I/O */
  iparm[IPARM_RHS_MAKING]            = API_RHS_B;           /* generate rhs */
  iparm[IPARM_REFINEMENT]            = API_RAF_GMRES;       /* gmres */
  iparm[IPARM_SYM]                   = API_SYM_YES;         /* Symmetric */
  iparm[IPARM_INCOMPLETE]            = API_NO;              /* direct */
  iparm[IPARM_ABS]                   = 1;                   /* ABS level to 1 */
  iparm[IPARM_ESP]                   = API_NO;              /* no esp */
#ifdef OOC
  iparm[IPARM_GMRES_IM]              = 1;                   /* gmres_im */
  iparm[IPARM_ITERMAX]               = 1;
#else
  iparm[IPARM_GMRES_IM]              = 25;                  /* gmres_im */
#endif
  iparm[IPARM_FREE_CSCUSER]          = API_CSC_PRESERVE;    /* Free user csc after coefinit */
  iparm[IPARM_FREE_CSCPASTIX]        = API_CSC_PRESERVE;    /* Free internal csc after coefinit */
  iparm[IPARM_OOC_LIMIT]             = 2000;                /* memory limit */
  iparm[IPARM_OOC_THREAD]            = 1;                   /* ooc thrdnbr */
  iparm[IPARM_OOC_ID]                = -1;                  /* Out of core run ID */
  iparm[IPARM_NB_SMP_NODE_USED]      = 0;                   /* Nb SMP node used (0 for 1 per MPI process) */
  iparm[IPARM_MURGE_REFINEMENT]      = API_YES;
  iparm[IPARM_MURGE_MAY_REFINE]      = API_NO;
  iparm[IPARM_TRANSPOSE_SOLVE]       = API_NO;
  /* Mode pour thread_comm :
     0 -> inutilis�
     1 -> 1 seul
     2 -> iparm[IPARM_NB_THREAD_COMM]
     3 -> Nbproc(iparm[IPARM_THREAD_NBR]))
  */
#ifdef PASTIX_FUNNELED
  iparm[IPARM_THREAD_COMM_MODE]      = API_THREAD_FUNNELED;
#else
  iparm[IPARM_THREAD_COMM_MODE]      = API_THREAD_MULTIPLE;
#endif
  iparm[IPARM_NB_THREAD_COMM]        = 1;                   /* Nb thread quand iparm[IPARM_THREAD_COMM_MODE] == API_THCOMM_DEFINED */
  iparm[IPARM_FILL_MATRIX]           = API_NO;              /* fill matrix */
  iparm[IPARM_INERTIA]               = -1;
  iparm[IPARM_ESP_NBTASKS]           = -1;
  iparm[IPARM_ESP_THRESHOLD]         = 16384;               /* Taille de bloc minimale pour passer en esp (2**14) = 128 * 128 */
  iparm[IPARM_ERROR_NUMBER]          = NO_ERR;
  iparm[IPARM_RHSD_CHECK]            = API_YES;
  iparm[IPARM_STARPU]                = API_NO;
  iparm[IPARM_AUTOSPLIT_COMM]        = API_NO;
#ifdef TYPE_COMPLEX
#  ifdef PREC_DOUBLE
  iparm[IPARM_FLOAT]                 = API_COMPLEXDOUBLE;
#  else
  iparm[IPARM_FLOAT]                 = API_COMPLEXSINGLE;
#  endif
#else
#  ifdef PREC_DOUBLE
  iparm[IPARM_FLOAT]                 = API_REALDOUBLE;
#  else
  iparm[IPARM_FLOAT]                 = API_REALSINGLE;
#  endif
#endif
  iparm[IPARM_STARPU_CTX_DEPTH]      = 3;
  iparm[IPARM_STARPU_CTX_NBR]        = -1;
  iparm[IPARM_PRODUCE_STATS]         = API_NO;
#ifdef PREC_DOUBLE
  dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
#else
  dparm[DPARM_EPSILON_REFINEMENT] = 1e-6;
#endif
  dparm[DPARM_RELATIVE_ERROR]     = -1;
  dparm[DPARM_SCALED_RESIDUAL]    = -1;
  dparm[DPARM_EPSILON_MAGN_CTRL]  = 1e-31;
  dparm[DPARM_FACT_FLOPS]         = 0;
  dparm[DPARM_SOLV_FLOPS]         = 0;
  dparm[63]                       = 0;

#ifdef PASTIX_LOG
  fprintf(stdout, "<- api_init\n");
#endif
}


/*
  Function: redispatch_rhs

  Redistribute right-hand-side member from *l2g*
  to *newl2g*

  Parameters:
  n        - Size of first right-hand-side member
  rhs      - Right-hand-side member
  l2g      - local to global column numbers
  newn     - New right-hand-side member size
  newrhs   - New right-hand-side member
  newl2g   - New local to global column numbers
  commSize - MPI communicator size
  commRank - MPI rank
  comm     - MPI communicator
*/
int redispatch_rhs(PASTIX_INT      n,
                   PASTIX_FLOAT   *rhs,
                   PASTIX_INT      nrhs,
                   PASTIX_INT     *l2g,
                   PASTIX_INT      newn,
                   PASTIX_FLOAT   *newrhs,
                   PASTIX_INT     *newl2g,
                   int      commSize,
                   int      commRank,
                   MPI_Comm comm,
                   PASTIX_INT dof)
{
  PASTIX_INT           i,j, irhs;
#ifndef FORCE_NOMPI
  MPI_Status    status;
#endif
  PASTIX_INT           gN = -1;
  PASTIX_INT          *g2l;
  PASTIX_INT         **count;
  MPI_Request * requests = NULL;
  PASTIX_INT   ** toSendIdx;
  PASTIX_FLOAT ** toSendValues;

  cscd_build_g2l(newn,
                 newl2g,
                 comm,
                 &gN,
                 &g2l);

  /* allocate counter */
  MALLOC_INTERN(count, commSize, PASTIX_INT*);
  for (i = 0; i < commSize; i++)
    {
      MALLOC_INTERN(count[i], commSize, PASTIX_INT);
      for (j = 0; j < commSize; j++)
        count[i][j] = 0;
    }

  /* count how many entries we have to send
     to each processor
     complexity : n
  */
  for (i = 0; i < n; i++)
    {
      PASTIX_INT globalidx = l2g[i];
      PASTIX_INT localidx  = g2l[globalidx-1];
      if ( localidx > 0) {
        count[commRank][commRank] ++;
      }
      else {
        count[commRank][-localidx]++;
      }
    }

  /* Broadcast counters */
  for (i = 0; i < commSize; i++)
    MPI_Bcast(count[i], commSize, COMM_INT, i, comm);

  MALLOC_INTERN(toSendIdx, commSize, PASTIX_INT *);
  MALLOC_INTERN(toSendValues, commSize, PASTIX_FLOAT *);
  for (irhs = 0; irhs < nrhs; irhs++)
    {
      /* Send data */
      for (i = 0; i < commSize; i++)
        {
          if (i != commRank) {
            MALLOC_INTERN(toSendIdx[i], count[commRank][i], PASTIX_INT);
            MALLOC_INTERN(toSendValues[i], count[commRank][i]*dof, PASTIX_FLOAT);
            memset(toSendIdx[i], 0, count[commRank][i]*sizeof(PASTIX_INT));
            memset(toSendValues[i], 0, count[commRank][i]*sizeof(PASTIX_FLOAT)*dof);
          }
        }

      for (i = 0; i < commSize; i++)
        count[commRank][i] = 0;

      for (i = 0; i < n; i++)
        {
          PASTIX_INT globalidx = l2g[i];
          PASTIX_INT localidx  = g2l[globalidx-1];
          if ( localidx > 0) {
            PASTIX_INT d;
            for (d = 0; d < dof; d++)
              newrhs[dof*(irhs*newn + localidx-1)+d] = rhs[dof*(irhs*n+i)+d];
          }
          else {
            PASTIX_INT d;
            toSendIdx[-localidx][count[commRank][-localidx]] = globalidx;
            for (d = 0; d < dof; d++)
              toSendValues[-localidx][dof*count[commRank][-localidx]+d] =
                rhs[dof*(irhs*n+i)+d];
            count[commRank][-localidx]++;
          }
        }

      MALLOC_INTERN(requests, commSize, MPI_Request);

      for (i = 0; i < commSize; i++)
        {
          if (commRank != i)
            {
              MPI_Isend(toSendIdx[i], count[commRank][i],
                        COMM_INT, i, TAG_RHS, comm, &requests[i]);
              MPI_Isend(toSendValues[i], dof*count[commRank][i],
                        COMM_FLOAT, i, TAG_RHS2, comm, &requests[i]);
            }
        }

      for (i = 0; i < commSize; i++)
        {
          if (commRank != i)
            {
              PASTIX_INT   * tmpIdx;
              PASTIX_FLOAT * tmpValues;
              MALLOC_INTERN(tmpIdx, count[i][commRank], PASTIX_INT);
              MALLOC_INTERN(tmpValues, count[i][commRank]*dof, PASTIX_FLOAT);

              MPI_Recv(tmpIdx, count[i][commRank], COMM_INT,
                       i, TAG_RHS, comm, &status);
              MPI_Recv(tmpValues, dof*count[i][commRank], COMM_FLOAT,
                       i, TAG_RHS2, comm, &status);
              for (j= 0; j < count[i][commRank]; j++)
                {
                  PASTIX_INT d;
                  for (d = 0; d < dof; d++)
                    newrhs[dof*(irhs*newn+g2l[tmpIdx[j]-1]-1)+d] = tmpValues[dof*j+d];
                }
              memFree_null(tmpIdx);
              memFree_null(tmpValues);
            }
        }

      for (i = 0; i < commSize; i++)
        {
          if (i != commRank)
            {
              MPI_Wait(&requests[i],&status);
              memFree_null(toSendIdx[i]);
              memFree_null(toSendValues[i]);
            }
        }
      memFree_null(requests);

    }
  memFree_null(toSendIdx);
  memFree_null(toSendValues);

  for (i = 0; i < commSize; i++)
    memFree_null(count[i]);
  memFree_null(count);
  memFree_null(g2l);
  return NO_ERR;
}

/*
  Function: buildUpDoVect

  Build UpDownVector from user B vector or
  computes its to obtain $$ \forall i X[i] = 1 $$ or $$ \forall i X[i] = i $$
  as the solution. (depending on iparm)

  Parameters:
  pastix_data - PaStiX global data structure.
  loc2glob2   - Global  column number of local columns.
  b           - User right-hand-side member.
  pastix_comm - MPI communicator.
*/
int buildUpdoVect(pastix_data_t *pastix_data,
                  PASTIX_INT           *loc2glob,
                  PASTIX_FLOAT         *b,
                  MPI_Comm       pastix_comm)
{
  PASTIX_INT           * iparm    = pastix_data->iparm;
  SolverMatrix  * solvmatr = &(pastix_data->solvmatr);
  Order         * ordemesh = &(pastix_data->ordemesh);
  PASTIX_INT             procnum  = pastix_data->procnum;
  PASTIX_INT           * invp     = ordemesh->peritab;
  (void)loc2glob;

  /* Rhs taking from b */
  if (iparm[IPARM_RHS_MAKING]==API_RHS_B)
    {
      /* Using b */
      if (solvmatr->updovct.sm2xsze > 0 &&  b==NULL)
        {
          errorPrint("b must be allocated.");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }
      else
        {
          /* Permuter b avec la permutation inverse */
          if (iparm[IPARM_GRAPHDIST] == API_NO )
            {
              CscUpdownRhs(&(solvmatr->updovct),
                           solvmatr,
                           b,
                           invp,
                           (int)iparm[IPARM_DOF_NBR]);
            }
#ifdef DISTRIBUTED
          else
            {
              CscdUpdownRhs(&(solvmatr->updovct),
                            solvmatr,
                            b,
                            invp,
                            pastix_data->glob2loc,
                            pastix_data->n,
                            (int)iparm[IPARM_DOF_NBR]);
            }
#endif /* DISTRIBUTED */
        }
    }
  /* Generate rhs */
  else
    {
      if (procnum == 0)
        {
          if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
            {
              if (iparm[IPARM_RHS_MAKING]==API_RHS_1)
                fprintf(stdout,GEN_RHS_1);
              else
                fprintf(stdout,GEN_RHS_I);
            }
        }

      /* In updo step, if we only want to use reffinement.
         Set first solution.
      */
      if ((iparm[IPARM_ONLY_RAFF] == API_YES) && (iparm[IPARM_START_TASK] < API_TASK_REFINE))
        {
          fprintf(stdout,GEN_SOL_0);
          Csc2updown_X0(&(solvmatr->updovct),
                        solvmatr,
                        API_RHS_0,
                        pastix_comm);
        }
      else
        {
          Csc2updown(&(pastix_data->cscmtx),
                     &(solvmatr->updovct),
                     solvmatr,
                     iparm[IPARM_RHS_MAKING],
                     pastix_comm);
        }
    }
  return NO_ERR;
}

/*
  Function: global2localrhs

  Converts global right hand side to local right hand side.

  Parameters:
  lN       - local number of columns.
  lrhs     - local right hand side.
  grhs     - global right hand side.
  loc2glob - global index of local columns.
*/

void global2localrhs(PASTIX_INT     lN,
                     PASTIX_FLOAT *lrhs,
                     PASTIX_FLOAT *grhs,
                     PASTIX_INT   *loc2glob)
{
  PASTIX_INT   i;

  for (i = 0; i < lN; i++)
    lrhs[i] = grhs [loc2glob[i]-1];
}

/*
  Function: global2localperm

  Converts global permutation (resp. reverse permutation) tabular to local
  permutation (resp. reverse permutation) tabular.

  Parameters:
  lN       - local number of columns.
  lperm    - local permutation tabular.
  gperm    - global permutation tabular.
  loc2glob - global index of local columns.
*/
void global2localperm(PASTIX_INT  lN,
                      PASTIX_INT *lperm,
                      PASTIX_INT *gperm,
                      PASTIX_INT *loc2glob)
{
  PASTIX_INT   i;
  for (i = 0; i < lN; i++)
    {
      lperm[i] = gperm[loc2glob[i]-1];
    }
}

/**
   Function: sizeofsolver

   Computes the size in memory of the SolverMatrix.

   Parameters:
   solvptr - address of the SolverMatrix

   Returns:
   SolverMatrix size.
*/
PASTIX_INT sizeofsolver(SolverMatrix *solvptr, PASTIX_INT *iparm)
{
  PASTIX_INT result=sizeof(SolverMatrix);
  PASTIX_INT iter;

  /* symbol */
  result += solvptr->cblknbr*sizeof(SymbolCblk);
  result += solvptr->bloknbr*sizeof(SymbolBlok);
  result += solvptr->cblknbr*sizeof(SolverCblk);
  result += solvptr->bloknbr*sizeof(SolverBlok);

  /* fanin target */
  result += solvptr->ftgtnbr*sizeof(FanInTarget);
  result += solvptr->btagnbr*sizeof(BlockTarget);
  result += solvptr->bcofnbr*sizeof(BlockTarget);
  result += solvptr->indnbr *sizeof(PASTIX_INT);

  /* task */
  result += solvptr->tasknbr*sizeof(Task);
  result += solvptr->thrdnbr*sizeof(PASTIX_INT);
  for (iter=0; iter<solvptr->thrdnbr; iter++)
    {
      result += solvptr->ttsknbr[iter]*sizeof(PASTIX_INT);
    }
  result += solvptr->procnbr*sizeof(PASTIX_INT);

  /* Pour l'instant uniquement si on est en 1d */
  if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
    {
      /* Updown */
      result += solvptr->cblknbr      *sizeof(UpDownCblk);
      for (iter=0; iter<solvptr->cblknbr; iter++)
        {
          result += 2*solvptr->updovct.cblktab[iter].browprocnbr*sizeof(PASTIX_INT);
        }
      result += solvptr->updovct.gcblk2listnbr*sizeof(PASTIX_INT);
      result += solvptr->updovct.listptrnbr   *sizeof(PASTIX_INT);
      result += 2*solvptr->updovct.listnbr    *sizeof(PASTIX_INT);
      result += solvptr->updovct.loc2globnbr  *sizeof(PASTIX_INT);
      result += solvptr->bloknbr              *sizeof(PASTIX_INT);
    }

  return result;
}

/*
 * Function: pastix_task_init
 *
 * Allocate and fill-in pastix_data
 *
 * Parameters:
 *   pastix_data - structure to build
 *   pastix_comm - PaStiX MPI communicator
 *   iparm       - integer parameters, to fill-in pastix_data
 *   dparm       - floating parameters, to fill-in pastix_data
 */
void pastix_task_init(pastix_data_t **pastix_data,
                      MPI_Comm        pastix_comm,
                      PASTIX_INT     *iparm,
                      double         *dparm)
{

  /* Allocate pastix_data structure when we enter PaStiX for the first time.
   */
  MALLOC_INTERN(*pastix_data, 1, struct pastix_data_t);

  /* Initialisation des champs de la structure */
  (*pastix_data)->n                = -1;
  (*pastix_data)->bmalcolrow       = 0;
#ifdef WITH_SCOTCH
  (*pastix_data)->malgrf           = 0;
#endif
  (*pastix_data)->malord           = 0;
  (*pastix_data)->malcsc           = 0;
  (*pastix_data)->malsmx           = 0;
  (*pastix_data)->malslv           = 0;
  (*pastix_data)->malcof           = 0;
  (*pastix_data)->iparm            = iparm;
  (*pastix_data)->dparm            = dparm;
  (*pastix_data)->pastix_comm      = pastix_comm;
  if (iparm != NULL && iparm[IPARM_AUTOSPLIT_COMM] == API_YES)
    {
      int i,len, mpisize;
      char procname[MPI_MAX_PROCESSOR_NAME];
      int color, key;

      MPI_Get_processor_name(procname,&len);
      MPI_Comm_rank(pastix_comm, &(key));
      color = 0;
      for (i = 0; i < len; i++)
        color = color*256*sizeof(char) + procname[i];
      MPI_Comm_split(pastix_comm, color, key, &(*pastix_data)->intra_node_comm);
      MPI_Comm_rank((*pastix_data)->intra_node_comm, &color);
      MPI_Comm_rank((*pastix_data)->intra_node_comm, &(key));
      MPI_Comm_split(pastix_comm, color, key, &(*pastix_data)->inter_node_comm);
      MPI_Comm_size((*pastix_data)->intra_node_comm, &mpisize);
      iparm[IPARM_THREAD_NBR] = mpisize;
    }
  else
    {
      (*pastix_data)->inter_node_comm      = pastix_comm;
      (*pastix_data)->intra_node_comm      = MPI_COMM_SELF;
    }

  (*pastix_data)->sopar.bindtab    = NULL;
  (*pastix_data)->sopar.b          = NULL;
  (*pastix_data)->sopar.transcsc   = NULL;
  (*pastix_data)->sopar.stopthrd   = API_NO;
  (*pastix_data)->bindtab          = NULL;
  (*pastix_data)->col2             = NULL;
  (*pastix_data)->row2             = NULL;

#ifdef DISTRIBUTED
  (*pastix_data)->loc2glob2        = NULL;
  (*pastix_data)->malrhsd_int      = API_NO;
  (*pastix_data)->l2g_int          = NULL;
  (*pastix_data)->mal_l2g_int      = API_NO;
  (*pastix_data)->glob2loc         = NULL;
  (*pastix_data)->PTS_permtab      = NULL;
  (*pastix_data)->PTS_peritab      = NULL;
#endif
  (*pastix_data)->schur_tab        = NULL;
  (*pastix_data)->schur_tab_set    = API_NO;
  (*pastix_data)->listschur        = NULL;

  (*pastix_data)->solvmatr.updovct.cblktab    = NULL;
  (*pastix_data)->solvmatr.updovct.lblk2gcblk = NULL;
  (*pastix_data)->solvmatr.updovct.listblok   = NULL;
  (*pastix_data)->solvmatr.updovct.listcblk   = NULL;
  (*pastix_data)->solvmatr.updovct.gcblk2list = NULL;
  (*pastix_data)->solvmatr.updovct.loc2glob   = NULL;
  (*pastix_data)->solvmatr.updovct.cblktab    = NULL;
  (*pastix_data)->solvmatr.updovct.listptr    = NULL;

  (*pastix_data)->scaling  = API_NO;
  (*pastix_data)->scalerowtab = NULL;
  (*pastix_data)->iscalerowtab = NULL;
  (*pastix_data)->scalecoltab = NULL;
  (*pastix_data)->iscalecoltab = NULL;

  (*pastix_data)->cscInternFilled = API_NO;
  /* R�cup�ration du nombre de proc */
  MPI_Comm_size(pastix_comm, &((*pastix_data)->procnbr));
  MPI_Comm_rank(pastix_comm, &((*pastix_data)->procnum));
  MPI_Comm_size((*pastix_data)->inter_node_comm, &((*pastix_data)->inter_node_procnbr));
  MPI_Comm_rank((*pastix_data)->inter_node_comm, &((*pastix_data)->inter_node_procnum));
  MPI_Comm_size((*pastix_data)->intra_node_comm, &((*pastix_data)->intra_node_procnbr));
  MPI_Comm_rank((*pastix_data)->intra_node_comm, &((*pastix_data)->intra_node_procnum));
  if ((*pastix_data)->procnum == 0)
    {
      (*pastix_data)->pastix_id = getpid();
    }
  MPI_Bcast(&((*pastix_data)->pastix_id), 1, COMM_INT, 0, pastix_comm);

#ifdef WITH_SEM_BARRIER
  if ((*pastix_data)->intra_node_procnbr > 1)
    {
      char sem_name[256];
      sprintf(sem_name, "/pastix_%d", (*pastix_data)->pastix_id);
      OPEN_SEM((*pastix_data)->sem_barrier, sem_name, 0);
    }
#endif

  if (iparm != NULL)
    {
      if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        {
          fprintf(stdout, "AUTOSPLIT_COMM : global rank : %d,"
                  " inter node rank %d,"
                  " intra node rank %d, threads %d\n",
                  (int)((*pastix_data)->procnum),
                  (int)((*pastix_data)->inter_node_procnum),
                  (int)((*pastix_data)->intra_node_procnum),
                  (int)iparm[IPARM_THREAD_NBR]);
        }

      iparm[IPARM_PID] = (*pastix_data)->pastix_id;
    }
  /* DIRTY Initialization for Scotch */
  srand(1);

  /* Environement variables */
  /* On Mac set VECLIB_MAXIMUM_THREADS if not setted */
  setenv("VECLIB_MAXIMUM_THREADS", "1", 0);
}

#ifdef MEMORY_USAGE
/*
  Function: pastix_print_memory_usage

  print memory usage during pastix.

  Parameters:
  iparm       - integer par�ameters.
  pastix_comm - PaStiX MPI communicator
*/
void pastix_print_memory_usage(PASTIX_INT      *iparm,
                               MPI_Comm  pastix_comm)
{
  unsigned long smem[2], rmem[2];
  int           procnum;

  MPI_Comm_rank(pastix_comm,&procnum);
  smem[0] = memAllocGetMax();
  smem[1] = memAllocGetCurrent();
  MPI_Reduce(smem,rmem,2,MPI_LONG,MPI_MAX,0,pastix_comm);
  iparm[DPARM_MEM_MAX] = rmem[0];
  if (procnum == 0)
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
      {
        fprintf(stdout,MAX_MEM_AF_CL,
                MEMORY_WRITE(rmem[0]),
                MEMORY_UNIT_WRITE(rmem[0]));

        fprintf(stdout,MEM_USED_AF_CL,
                MEMORY_WRITE(rmem[1]),
                MEMORY_UNIT_WRITE(rmem[1]));
      }
}
#else
#  ifdef pastix_print_memory_usage
#    undef pastix_print_memory_usage
#    define pastix_print_memory_usage(pastix_data,pastix_comm)
#  endif
#endif /* MEMORY_USAGE */

/*
 * Function: pastix_welcome_print
 *
 * Will print welcome message, options and parameters.
 *
 * Parameters:
 * pastix_data - PaStiX data structure
 * colptr      - starting index of each column in the CSC.
 * n           - number of columns.
 *
 */
void pastix_welcome_print(pastix_data_t *pastix_data,
                          PASTIX_INT           *colptr,
                          PASTIX_INT            ln)
{
  PASTIX_INT    * iparm = pastix_data->iparm;
  double * dparm = pastix_data->dparm;
  PASTIX_INT      gN    = 0;
  PASTIX_INT      lNnz  = 0;
  PASTIX_INT      gNnz  = 0;

  if (colptr != NULL)
    lNnz = colptr[ln]-colptr[0];

  if (iparm[IPARM_GRAPHDIST] == API_YES)
    {
      gN = 0;
      MPI_Allreduce(&ln, &gN, 1, COMM_INT, MPI_SUM,
                    pastix_data->pastix_comm);
      MPI_Allreduce(&lNnz, &gNnz, 1, COMM_INT, MPI_SUM,
                 pastix_data->pastix_comm);
    }
  else
    {
      gN   = ln;
      gNnz = lNnz;
    }
  dparm[DPARM_FILL_IN]       = 1.0/(double)gNnz;

  if ((iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) &&
      (pastix_data->procnum == 0 ))
    {
      fprintf(stdout, OUT_ENTETE_LINE1);
      fprintf(stdout, OUT_ENTETE_LINE2);
      fprintf(stdout, OUT_ENTETE_LINE3);

      /* TODO : en distribu� afficher pour chaque proc... */
      if ( gN != 0)
        {
          fprintf(stdout, OUT_MATRIX_SIZE, (long)gN, (long)gN);
          if (gNnz != 0) fprintf(stdout, OUT_NNZ, (long)gNnz);
        }
      sopalin_option();
    }
}

#ifdef WITH_SCOTCH
/*
  Function: pastix_order_save

  Save ordering structures on disk.

  Parameters:
  ordemesh - Scotch ordering structure to save.
  grafmesh - Scotch Graph structure to save.
  ncol     - Number of column in the CSC
  colptr   - starting index of each column in row
  rows     - row of each element.
  values   - value of each element.
  strategy - IO strategy.

*/
int pastix_order_save(Order        * ordemesh,
                      SCOTCH_Graph * grafmesh,
                      int            procnum,
                      PASTIX_INT            ncol,
                      PASTIX_INT          * colptr,
                      PASTIX_INT          * rows,
                      PASTIX_INT            strategy)
{
  FILE             * stream;
  int                retval     = NO_ERR;
#  ifndef WITH_SCOTCH
  errorPrint("Saving strategy needs to compile PaStiX with -DWITH_SCOTCH");
  retval = BADPARAMETER_ERR;

#  else
  if (procnum == 0)
    {
      PASTIX_FOPEN(stream, "ordergen","w");
      if (orderSave (ordemesh, stream) != 0)
        {
          errorPrint ("cannot save order");
          retval = INTERNAL_ERR;
        }
      fclose(stream);
      if (PASTIX_MASK_ISTRUE(strategy, API_IO_SAVE_GRAPH))
        {
          PASTIX_FOPEN(stream, "graphgen","w");
          if (SCOTCH_graphSave (grafmesh, stream) != 0)
            {
              errorPrint ("cannot save graph");
              retval = INTERNAL_ERR;
            }
          fclose(stream);
        }
      if (PASTIX_MASK_ISTRUE(strategy, API_IO_SAVE_CSC))
        {
          PASTIX_FOPEN(stream, "cscgen","w");
          retval = csc_save(ncol, colptr, rows, NULL, 1, stream);
          fclose(stream);
        }
    }
#  endif
  return retval;
}


/*
  Function: pastix_order_load
  Load ordering structures from disk.

  Parameters:
  ordemesh - Scotch ordering structure to save.
  grafmesh - Scotch Graph structure to save.
  ncol     - Number of column in the CSC
  colptr   - starting index of each column in row
  rows     - row of each element.
  values   - value of each element.
  startegy - IO strategy.
  comm     - MPI communicator.


*/
int pastix_order_load(Order        *  ordemesh,
                      SCOTCH_Graph *  grafmesh,
                      int             procnum,
                      PASTIX_INT          *  ncol,
                      PASTIX_INT          ** colptr,
                      PASTIX_INT          ** rows,
                      PASTIX_INT             strategy,
                      MPI_Comm        comm)
{
  FILE             * stream;
  int                retval     = NO_ERR;
  int                dof;
  (void)comm;

#  ifndef WITH_SCOTCH
  errorPrint("Loading strategy needs to compile PaStiX with -DWITH_SCOTCH");
  retval = BADPARAMETER_ERR;
  break;
#  else

  /* Load scotch result */

  if (PASTIX_MASK_ISTRUE(strategy, API_IO_LOAD_GRAPH))
    {
      PASTIX_FOPEN(stream, "graphname","r");
      if (SCOTCH_graphLoad(grafmesh, stream, 0, 0) != 0) {
        errorPrint ("test: cannot load mesh");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
      }
      fclose (stream);
    }
  PASTIX_FOPEN(stream, "ordername", "r");
  if (orderLoad(ordemesh, stream) != 0)
    {
      errorPrint("test: cannot load order");
      EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }
  fclose(stream);
  if (PASTIX_MASK_ISTRUE(strategy, API_IO_LOAD_CSC))
    {
      if (procnum == 0)
        {
          PASTIX_FOPEN(stream, "cscname","r");
          retval = csc_load(ncol, colptr, rows, NULL, &dof, stream);
          fclose(stream);
        }
      MPI_Bcast(ncol, 1, COMM_INT, 0, comm);
      if (procnum != 0)
        {
          MALLOC_INTERN((*colptr), *ncol+1, PASTIX_INT);
        }
      MPI_Bcast(*colptr, *ncol+1, COMM_INT, 0, comm);
      if  (procnum != 0)
        {
          MALLOC_INTERN(*rows, (*colptr)[*ncol]-(*colptr)[0], PASTIX_INT);
        }
      MPI_Bcast(*rows, (*colptr)[*ncol]-(*colptr)[0], COMM_INT, 0, comm);
    }
#  endif
  return retval;
}
#endif
/*
  Function: pastix_order_prepare_csc

  Create a copy of user's CSC and prepare it for ordering step.

  Symmetrize the graph and removes diagonal coefficients.

  Parameters:
  pastix_data - PaStiX internal data structure
  n           - Number of column in the CSC.
  colptr      - Start of each column in *rows* array.
  rows        - Row number of each non zeros.
*/
int pastix_order_prepare_csc(pastix_data_t * pastix_data,
                             PASTIX_INT             n,
                             PASTIX_INT           * colptr,
                             PASTIX_INT           * rows)
{
  PASTIX_INT * iparm;
  int procnum;

  iparm   = pastix_data->iparm;
  procnum = (int)pastix_data->procnum;
  /* Allocate a copy of col, row */
  pastix_data->bmalcolrow = 1;
  pastix_data->n2   = n;
  MALLOC_INTERN(pastix_data->col2, (n+1),         PASTIX_INT);
  MALLOC_INTERN(pastix_data->row2, (colptr[n]-1), PASTIX_INT);
  memcpy((void*) pastix_data->col2,(void*)colptr,        (n+1)*sizeof(PASTIX_INT));
  memcpy((void*) pastix_data->row2,(void*)rows,  (colptr[n]-1)*sizeof(PASTIX_INT));

  /* Symmetrize the graph */
  if (iparm[IPARM_SYM] == API_SYM_YES || iparm[IPARM_SYM] == API_SYM_HER)
    {
      PASTIX_INT tmpn;
      PASTIX_INT *tmpcol;
      PASTIX_INT *tmprow;

      if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        print_onempi("%s", OUT_SYMGRAPH);
      csc_symgraph_int(n,     pastix_data->col2,   pastix_data->row2,   NULL,
                       &tmpn, &tmpcol, &tmprow, NULL, API_YES);

      memFree_null((pastix_data->col2));
      memFree_null((pastix_data->row2));
      pastix_data->col2 = tmpcol;
      pastix_data->row2 = tmprow;
    }

  /* Remove diagonal coefficient */
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
    print_onempi("%s", OUT_NODIAG);
  csc_noDiag(1, n, pastix_data->col2, pastix_data->row2, NULL);

  return NO_ERR;
}
/*
  Function: pastix_task_scotch

  Execute ordering task, with a centralised graph.

  Free *col2*  and *row2* entries of pastix_data if <pastix_task_scotch>
  has already been called.

  Set *col2*, *row2* and *n2* to a copy of user's CSC.

  Symmetrize this CSC.

  Remove diagonal elements from it.

  Clean last oredering if it exists.
  Depending on *IPARM_ORDERING* :
  - Calls Scotch ordering,
  - Calls Metis ordering,
  - Uses user oredering,
  - Loads oredering stored on disk in a Scotch format.

  Can save computed ordering on disk.

  returns compuited ordering into user arrays.

  Parameters:
  pastix_data - PaStiX data structure.
  pastix_comm - PaStiX MPI communicator.
  n           - size of the matrix/number of vertices.
  colptr      - starting index of each column in row
  row         - row of each element
  perm        - permutation tabular
  invp        - reverse permutation tabular
*/
int pastix_task_scotch(pastix_data_t **pastix_data,
                       MPI_Comm        pastix_comm,
                       PASTIX_INT             n,
                       PASTIX_INT            *colptr,
                       PASTIX_INT            *row,
                       PASTIX_INT            *perm,
                       PASTIX_INT            *invp)
{
  PASTIX_INT              * iparm       = (*pastix_data)->iparm;
  PASTIX_INT             ** col2;
  PASTIX_INT             ** row2;
#ifdef WITH_SCOTCH
  SCOTCH_Graph     * grafmesh;
  SCOTCH_Strat       stratdat;
  char               strat[550];
  PASTIX_INT               *colptr_schur  = NULL;
  PASTIX_INT               *rows_schur    = NULL;
  PASTIX_INT               *perm_schur    = NULL;
  PASTIX_INT               *revperm_schur = NULL;
#endif
  Order            * ordemesh;
  Clock              timer1;
  PASTIX_INT                procnum;
  PASTIX_INT                iter;
  int                retval     = NO_ERR;
  int                retval_rcv;

#ifdef WITH_SCOTCH
  grafmesh  = &((*pastix_data)->grafmesh);
#endif
  ordemesh  = &((*pastix_data)->ordemesh);
  procnum   =   (*pastix_data)->procnum;

  col2      = &((*pastix_data)->col2);
  row2      = &((*pastix_data)->row2);

  print_debug(DBG_STEP,"-> pastix_task_scotch\n");
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_onempi("%s", OUT_STEP_ORDER);

  if ((*pastix_data)->bmalcolrow == 1)
    {
      if ((*col2)      != NULL) memFree_null(*col2);
      if ((*row2)      != NULL) memFree_null(*row2);
      (*pastix_data)->bmalcolrow = 0;
    }

  /* Clean ordering if it exists */
  if ((*pastix_data)->malord)
    {
      orderExit(ordemesh);
      (*pastix_data)->malord=0;
    }

  /* Prepare a copy of user's CSC */
  if (!(PASTIX_MASK_ISTRUE(iparm[IPARM_ORDERING], API_ORDER_LOAD)))
    pastix_order_prepare_csc(*pastix_data, n, colptr, row);


  if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
    print_onempi("%s", OUT_ORDERINIT);
  orderInit(ordemesh);
  if (iparm[IPARM_ORDERING] != API_ORDER_LOAD )
    {
      MALLOC_INTERN(ordemesh->permtab, n,     PASTIX_INT);
      MALLOC_INTERN(ordemesh->peritab, n,     PASTIX_INT);
      MALLOC_INTERN(ordemesh->rangtab, n + 1, PASTIX_INT);
    }
  (*pastix_data)->malord=1;

  clockInit(&timer1);
  clockStart(&timer1);

  switch (iparm[IPARM_ORDERING])
    {
      /*
       * Scotch Ordering
       */
    case API_ORDER_SCOTCH:
#ifndef WITH_SCOTCH
      errorPrint("Scotch ordering needs to compile PaStiX with -DWITH_SCOTCH");
      retval = BADPARAMETER_ERR;
      break;
#else
      {
        int ret;

        if (sizeof(PASTIX_INT) != sizeof(SCOTCH_Num))
          {
            errorPrint("Inconsistent integer type\n");
            retval = INTEGER_TYPE_ERR;
            break;
          }

        /* On nettoie grafmesh et col2/row2 si ils sont d�ja allou�s */
        if ((*pastix_data)->malgrf == 1)
          {
            SCOTCH_graphExit(grafmesh);
            (*pastix_data)->malgrf = 0;
          }

        /* construction du graphe */
        print_debug(DBG_SCOTCH, "> SCOTCH_graphInit <\n");
        SCOTCH_graphInit(grafmesh);

        print_debug(DBG_SCOTCH, "> SCOTCH_graphBuild <\n");
        if (iparm[IPARM_SCHUR] == API_YES ||
            iparm[IPARM_ISOLATE_ZEROS] == API_YES)
          {

            MALLOC_INTERN(colptr_schur, n+1, PASTIX_INT);
            MALLOC_INTERN(rows_schur, (*col2)[n]-1, PASTIX_INT);
            memcpy(colptr_schur, *col2, (n+1)*sizeof(PASTIX_INT));
            memcpy(rows_schur,   *row2, ((*col2)[n] - 1)*sizeof(PASTIX_INT));
            MALLOC_INTERN(perm_schur, n, PASTIX_INT);
            MALLOC_INTERN(revperm_schur, n, PASTIX_INT);
            CSC_isolate(n,
                        colptr_schur,
                        rows_schur,
                        (*pastix_data)->nschur,
                        (*pastix_data)->listschur,
                        perm_schur,
                        revperm_schur);

            memFree_null(revperm_schur);

            if (SCOTCH_graphBuild(grafmesh,                                   /* Graph to build     */
                                  1,                                          /* baseval            */
                                  n-(*pastix_data)->nschur,                   /* Number of vertices */
                                  colptr_schur,                               /* Vertex array       */
                                  NULL,
                                  NULL,                                       /* Array of vertex weights (DOFs) */
                                  NULL,
                                  (colptr_schur[n-(*pastix_data)->nschur]-1), /* Number of arcs     */
                                  rows_schur,                                 /* Edge array         */
                                  NULL))
              {
                errorPrint("pastix : graphBuildGraph");
                EXIT(MOD_SOPALIN,INTERNAL_ERR);
              }
          }
        else
          {
            if (SCOTCH_graphBuild(grafmesh,       /* Graph to build     */
                                  1,              /* baseval            */
                                  n,              /* Number of vertices */
                                  *col2,          /* Vertex array       */
                                  NULL,
                                  NULL,           /* Array of vertex weights (DOFs) */
                                  NULL,
                                  ((*col2)[n]-1), /* Number of arcs     */
                                  *row2,          /* Edge array         */
                                  NULL))
              {
                errorPrint("pastix : graphBuildGraph");
                EXIT(MOD_SOPALIN,INTERNAL_ERR);
              }
          }
        (*pastix_data)->malgrf=1;

        print_debug(DBG_SCOTCH, "> SCOTCH_graphCheck <\n");
        if (SCOTCH_graphCheck(grafmesh))
          {
            errorPrint("pastix : graphCheck");
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
          }

        SCOTCH_stratInit(&stratdat);
        SCOTCH_graphBase(grafmesh, 0);

        if (iparm[IPARM_DEFAULT_ORDERING] == API_YES) /* default ordering */
          {
            if (iparm[IPARM_INCOMPLETE] == API_NO)
              {
                if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
                  print_onempi("%s", "Scotch direct strategy\n");
                sprintf(strat, SCOTCH_STRAT_DIRECT);
              }
            else
              {
                if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
                  print_onempi("%s", "Scotch incomplete strategy\n");
                sprintf(strat, SCOTCH_STRAT_INCOMP);
              }
          }
        else /* personal ordering */
          {
            sprintf(strat, SCOTCH_STRAT_PERSO,
                    (long) iparm[IPARM_ORDERING_SWITCH_LEVEL],
                    (long) iparm[IPARM_ORDERING_CMIN],
                    (long) iparm[IPARM_ORDERING_CMAX],
                    ((float)iparm[IPARM_ORDERING_FRAT])/100,
                    (long) iparm[IPARM_ORDERING_SWITCH_LEVEL],
                    (long) iparm[IPARM_ORDERING_CMIN],
                    (long) iparm[IPARM_ORDERING_CMAX],
                    ((float)iparm[IPARM_ORDERING_FRAT])/100);
            if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
              print_onempi("Scotch personal strategy |%s|\n", strat);
          }

        ret = SCOTCH_stratGraphOrder (&stratdat, strat);
        if (ret == 0)
          {
            if (iparm[IPARM_SCHUR] == API_YES ||
                iparm[IPARM_ISOLATE_ZEROS] == API_YES)
              {
                PASTIX_INT *tmpperm = NULL;
                /* Compute graph ordering */
                ret =
                  SCOTCH_graphOrderList(grafmesh,
                                        (SCOTCH_Num)(n-(*pastix_data)->nschur),
                                        (SCOTCH_Num *) NULL,
                                        &stratdat,
                                        (SCOTCH_Num *)  ordemesh->permtab,
                                        (SCOTCH_Num *)  ordemesh->peritab,
                                        (SCOTCH_Num *) &ordemesh->cblknbr,
                                        (SCOTCH_Num *)  ordemesh->rangtab,
                                        NULL);

                /* Ajouter la permutation du Schur au permtab/peritab
                   Et un bloc au rangtab.
                */
                memFree_null(colptr_schur);
                memFree_null(rows_schur);
                ordemesh->rangtab[ordemesh->cblknbr+1] = n;
                ordemesh->cblknbr++;

                for(iter = n-(*pastix_data)->nschur; iter < n; iter++)
                  ordemesh->permtab[iter] = iter;

                for(iter = 0; iter < n; iter++)
                  {
                    ASSERT(ordemesh->permtab[iter] < n, MOD_SOPALIN);
                    ASSERT(ordemesh->permtab[iter] > -1, MOD_SOPALIN);
                    ASSERT(perm_schur[iter] < n, MOD_SOPALIN);
                    ASSERT(perm_schur[iter] > -1, MOD_SOPALIN);
                  }
                MALLOC_INTERN(tmpperm, n, PASTIX_INT);
                for(iter = 0; iter < n; iter++)
                  tmpperm[iter] = ordemesh->permtab[perm_schur[iter]];
                memcpy(ordemesh->permtab, tmpperm, n*sizeof(PASTIX_INT));
                memFree_null(tmpperm);
                memFree_null(perm_schur);
                for(iter = 0; iter < n; iter++)
                  ordemesh->peritab[ordemesh->permtab[iter]] = iter;
                for(iter = 0; iter < n; iter++)
                  {
                    ASSERT(ordemesh->peritab[iter] < n, MOD_SOPALIN);
                    ASSERT(ordemesh->peritab[iter] > -1, MOD_SOPALIN);
                  }
                /* rebuild graph for fax */
                if ((*pastix_data)->malgrf == 1)
                  {
                    SCOTCH_graphExit(grafmesh);
                    (*pastix_data)->malgrf = 0;
                  }
                if (SCOTCH_graphBuild(grafmesh,       /* Graph to build     */
                                      1,              /* baseval            */
                                      n,              /* Number of vertices */
                                      *col2,          /* Vertex array       */
                                      NULL,
                                      NULL,           /* Array of vertex weights (DOFs) */
                                      NULL,
                                      ((*col2)[n]-1), /* Number of arcs     */
                                      *row2,          /* Edge array         */
                                      NULL))
                  {
                    errorPrint("pastix : graphBuildGraph");
                    EXIT(MOD_SOPALIN,INTERNAL_ERR);
                  }
                (*pastix_data)->malgrf = 1;
                SCOTCH_graphBase(grafmesh, 0);
              }
            else
              {
                /* Compute graph ordering */
                ret = SCOTCH_graphOrderList (grafmesh,
                                             (SCOTCH_Num)   n,
                                             (SCOTCH_Num *) NULL,
                                             &stratdat,
                                             (SCOTCH_Num *)  ordemesh->permtab,
                                             (SCOTCH_Num *)  ordemesh->peritab,
                                             (SCOTCH_Num *) &ordemesh->cblknbr,
                                             (SCOTCH_Num *)  ordemesh->rangtab,
                                             NULL);
              }
          }
        SCOTCH_stratExit (&stratdat);
        if (ret != 0) {           /* If something failed in Scotch */
          orderExit (ordemesh);    /* Free ordering arrays          */
          orderInit (ordemesh);
          retval = INTERNAL_ERR;
          break;
        }

        /* Redimensionnement de rangtab a cblknbr */
        ordemesh->rangtab =
          (PASTIX_INT *) memRealloc (ordemesh->rangtab,
                              (ordemesh->cblknbr + 1)*sizeof (PASTIX_INT));

#  ifdef FORGET_PARTITION
        {
          memFree_null(ordemesh->rangtab);
          MALLOC_INTERN(ordemesh->rangtab, n+1, PASTIX_INT);
          for (iter=0;iter<n+1;iter++)
            ordemesh->rangtab[iter] = iter;
          ordemesh->cblknbr = n;
        }
#  endif
      }
#endif /* WITH_SCOTCH */
      break;


      /*
       *  METIS ordering
       */
    case API_ORDER_METIS:
#ifndef METIS
      errorPrint("Metis ordering needs to compile PaStiX with -DMETIS");
      retval = BADPARAMETER_ERR;
      break;
#else /* METIS */
      {
        PASTIX_INT  itervert;
        PASTIX_INT  baseval;
        PASTIX_INT  opt[8];

        baseval = 1;

        if (sizeof(PASTIX_INT) != sizeof(int))
          {
            errorPrint("Inconsistent integer type\n");
            retval = INTEGER_TYPE_ERR;
            break;
          }

        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
          print_onempi("%s", "calling metis...\n");

        /* call METIS and fill ordemesh (provide a partition) */
        opt[METIS_OPTION_PTYPE  ] = (iparm[IPARM_DEFAULT_ORDERING]==API_YES)?0:1;

        /* TODO: tester sans cette ligne... 0 if default */
        opt[METIS_OPTION_PTYPE  ] = 0;

        opt[METIS_OPTION_CTYPE  ] = iparm[IPARM_ORDERING_SWITCH_LEVEL];
        opt[METIS_OPTION_IPTYPE  ] = iparm[IPARM_ORDERING_CMIN];
        opt[METIS_OPTION_RTYPE  ] = iparm[IPARM_ORDERING_CMAX];
        opt[METIS_OPTION_DBGLVL ] = iparm[IPARM_ORDERING_FRAT];
        /*
         * NOTE Daniel Wirtz ABI Auckland: Disabled option so that Pastix would work with METIS 5.1.
         * As METIS will use defaults this doesnt break anything but may affect performance.
         */
        //opt[METIS_OPTION_OFLAGS ] = iparm[IPARM_STATIC_PIVOTING];
        opt[METIS_OPTION_PFACTOR] = iparm[IPARM_METIS_PFACTOR];
        opt[METIS_OPTION_NSEPS  ] = iparm[IPARM_NNZEROS];

        /*METIS_NodeND(&n,verttab,edgetab,&baseval,opt,
          ordemesh->permtab,ordemesh->peritab);*/
        METIS_NodeND(&n, *col2, *row2, &baseval,
                     opt, ordemesh->peritab, ordemesh->permtab);

        for (itervert=0; itervert<n+1; itervert++)
          ordemesh->rangtab[itervert] = itervert;
        ordemesh->cblknbr = n;
      }
#endif /* METIS */
      break;

      /*
       * Personal Ordering
       */
    case API_ORDER_PERSONAL:

      memcpy(ordemesh->permtab, perm, n*sizeof(PASTIX_INT));
      memcpy(ordemesh->peritab, invp, n*sizeof(PASTIX_INT));

      for (iter=0; iter<n+1; iter++)
        ordemesh->rangtab[iter] = iter;
      break;

      /*
       * Load ordering with Scotch Format
       */
    case API_ORDER_LOAD:
#ifdef WITH_SCOTCH
      pastix_order_load(ordemesh,
                        grafmesh,
                        procnum,
                        &((*pastix_data)->n2),
                        &((*pastix_data)->col2),
                        &((*pastix_data)->row2),
                        iparm[IPARM_IO_STRATEGY],
                        pastix_comm);
#endif
      break;

    default:
      errorPrint("Ordering not available");
      retval = BADPARAMETER_ERR;
      break;
    }

  MPI_Allreduce(&retval, &retval_rcv, 1, MPI_INT, MPI_MAX, pastix_comm);
  if (retval_rcv != NO_ERR)
    RETURN_ERROR(retval_rcv);

  orderBase(ordemesh, 0);

  clockStop(&timer1);
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    print_onempi(TIME_COMPUTE_ORDERING,clockVal(&timer1));

  /* Save i/o strategy */
  if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
    {

#ifdef WITH_SCOTCH
      retval = pastix_order_save(ordemesh,
                                 grafmesh,
                                 procnum,
                                 (*pastix_data)->n2,
                                 (*pastix_data)->col2,
                                 (*pastix_data)->row2,
                                 iparm[IPARM_IO_STRATEGY]);
      if (retval != NO_ERR)
        return retval;
#endif
    }

  /*
   * Return the ordering to user
   */
  if (iparm[IPARM_ORDERING] != API_ORDER_PERSONAL)
    {
      memcpy(perm, ordemesh->permtab, n*sizeof(PASTIX_INT));
      memcpy(invp, ordemesh->peritab, n*sizeof(PASTIX_INT));
    }

  iparm[IPARM_START_TASK]++;
  return NO_ERR;
}

#ifdef DISTRIBUTED
/*
  Function: dpastix_order_prepare_cscd

  Create a copy of user's CSCd and prepare it for ordering step.

  Symmetrize the graph and removes diagonal coefficients.

  Symetrize the graph, removes diagonal elements.

  Redistribute the CSCd to be abble to give it to scotch if needed.
  Indeed, PT-Scotch only allows consecutive columns.
  *PTS_permtab* is the permutation tabular from the user's distribution
  to PT-Scotch one, *PTS_peritab* is the reverse permutation.

  Parameters:
  pastix_data  - PaStiX internal data structure
  n            - Number of column in the CSC.
  colptr       - Start of each column in *rows* array.
  rows         - Row number of each non zeros.
  loc2glob     - local to global column number array.
  pastix_comm  - MPI communicator.
*/
int dpastix_order_prepare_cscd(pastix_data_t * pastix_data,
                               PASTIX_INT             n,
                               PASTIX_INT           * colptr,
                               PASTIX_INT           * rows,
                               PASTIX_INT           * loc2glob,
                               MPI_Comm        pastix_comm)
{
  PASTIX_INT   gN;
  PASTIX_INT * iparm;
  PASTIX_INT   i;
  int   OK, OKRecv;
  int * alln         = NULL;
  int   nlocal;
  int * displs       = NULL;

  iparm = pastix_data->iparm;
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX)
    fprintf(stdout, "> dpastix_order_prepare_cscd\n");

  if (pastix_data->PTS_permtab != NULL)
    {
      memFree_null(pastix_data->PTS_permtab);
      memFree_null(pastix_data->PTS_peritab);
    }
  /* Copy the cscd */
  pastix_data->n2 = n;
  MALLOC_INTERN(pastix_data->col2, n+1, PASTIX_INT);
  memcpy(pastix_data->col2, colptr, (n+1)*sizeof(PASTIX_INT));
  MALLOC_INTERN(pastix_data->loc2glob2, n, PASTIX_INT);
  memcpy((pastix_data->loc2glob2), loc2glob, n*sizeof(PASTIX_INT));
  MALLOC_INTERN(pastix_data->row2, colptr[n]-1, PASTIX_INT);
  memcpy(pastix_data->row2, rows, (colptr[n]-1)*sizeof(PASTIX_INT));

  gN = 0;
  MPI_Allreduce(&n, &gN, 1, COMM_INT, MPI_SUM, pastix_comm);

  /* Symmetrize the graph */
  if (iparm[IPARM_SYM]==API_SYM_YES || iparm[IPARM_SYM] == API_SYM_HER)
    {
      /* Symetric problem */
      /* Build non oriented graph */
      /* build non symmetric csc from symmetric csc */
      /*maillage global*/
      PASTIX_INT *tmpia;
      PASTIX_INT *tmpja;
      PASTIX_INT  tmpn;

      cscd_symgraph_int(pastix_data->n2,   pastix_data->col2,  pastix_data->row2 , NULL,
                        &tmpn, &tmpia, &tmpja, NULL,
                        pastix_data->loc2glob2, pastix_comm, API_YES);

      memFree_null(pastix_data->col2);
      pastix_data->col2 = tmpia;
      memFree_null(pastix_data->row2);
      pastix_data->row2 = tmpja;
      pastix_data->n2   = tmpn;
    }


#  ifdef DEBUG_DPASTIX
  cscd_save(pastix_data->n2,
            pastix_data->col2,
            pastix_data->row2,
            NULL, NULL, pastix_data->loc2glob2,
            iparm[IPARM_DOF_NBR], "cscd_after_sym", pastix_comm);
#  endif


  /* Remove diagonal coefficients */
  cscd_noDiag(pastix_data->n2, pastix_data->col2, pastix_data->row2, NULL, pastix_data->loc2glob2);
#  ifdef DEBUG_DPASTIX
  cscd_save(pastix_data->n2, pastix_data->col2, pastix_data->row2, NULL, NULL,
            pastix_data->loc2glob2, iparm[IPARM_DOF_NBR], "cscd_after_nodiag", pastix_comm);
#  endif
  if (pastix_data->procnbr > 1)
    {
      /* Check if matrix is not allready correct for scotch */
      /* PT-Scotch needs consecutives column blocs */
      OK = 0;
      OKRecv = 0;
      for (i = 0; i < n-1; i++)
        if (loc2glob[i] != loc2glob[i+1] -1)
          OK = 1;

      MPI_Allreduce(&OK, &OKRecv, 1, MPI_INT, MPI_SUM, pastix_comm);
      /* If it is not correct, permut it */
      if (OKRecv != 0)
        {

          /* Correct the cscd for scotch */
          MALLOC_INTERN(alln, pastix_data->procnbr, int);
          nlocal = n;
          MPI_Allgather(&nlocal, 1, MPI_INT,
                        alln, 1, MPI_INT,
                        pastix_comm);

          MALLOC_INTERN(displs, pastix_data->procnbr, int);
          displs[0] = 0;
          for (i = 1; i < pastix_data->procnbr; i++)
            displs[i] = displs[i-1] + alln[i-1];
          for (i = 0; i < n; i++)
            pastix_data->loc2glob2[i] = displs[pastix_data->procnum]+1+i;

          MALLOC_INTERN(pastix_data->PTS_peritab, gN, PASTIX_INT);
          MPI_Allgatherv(loc2glob, n, COMM_INT,
                         pastix_data->PTS_peritab, alln, displs, COMM_INT,
                         pastix_comm);

          memFree_null(displs);
          memFree_null(alln);
          MALLOC_INTERN(pastix_data->PTS_permtab, gN, PASTIX_INT);
          for (i = 0; i < gN; i++)
            pastix_data->PTS_permtab[pastix_data->PTS_peritab[i]-1] = i+1;

          /* Redistribue la cscd existante */
          for (i = 0; i < (pastix_data->col2)[n] - 1; i++)
            pastix_data->row2[i] = pastix_data->PTS_permtab[(pastix_data->row2)[i]-1];
        }
    }
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX)
    fprintf(stdout, "< dpastix_order_prepare_cscd\n");

  return NO_ERR;
}

/*
  Function: dpastix_task_scotch

  Execute ordering task, with a distributed graph.

  In LOAD mode, only load the graph from disk.

  Else, Clean col2, row2, loc2glob2 and grafmesh if ordering task
  has been called before.


  Build the graph and calls PT-Scotch ordering.

  Gather the graph to be abble to perform centralised version of symbolic
  factorisation.

  Save the graph if asked.

  Parameters:
  pastix_data - PaStiX data structure.
  pastix_comm - PaStiX MPI communicator.
  n           - size of the matrix/number of vertices.
  colptr      - starting index of each column in row
  row         - row of each element
  perm        - permutation tabular
  invp        - reverse permutation tabular
  loc2glob    - global index of local columns
*/
int dpastix_task_scotch(pastix_data_t ** pastix_data,
                        MPI_Comm         pastix_comm,
                        PASTIX_INT              n,
                        PASTIX_INT            * colptr,
                        PASTIX_INT            * row,
                        PASTIX_INT            * perm,
                        PASTIX_INT            * invp,
                        PASTIX_INT            * loc2glob)
{
#  ifndef WITH_SCOTCH
  errorPrint("Distributed PaStiX calls only works with -DWITH_SCOTCH");
  RETURN_ERROR(BADPARAMETER_ERR);
#  else
  PASTIX_INT              * iparm       = (*pastix_data)->iparm;
  PASTIX_INT              * n2;
  PASTIX_INT             ** col2;
  PASTIX_INT             ** row2;
  PASTIX_INT             ** loc2glob2;
  SCOTCH_Graph     * grafmesh;
  SCOTCH_Dordering * ordedat;
  SCOTCH_Ordering    ordering;
  SCOTCH_Dgraph    * dgraph;
  SCOTCH_Strat       stratdat;
  Order            * ordemesh;
  PASTIX_INT                gN;
  Clock              timer1;
  char               strat[550];
  PASTIX_INT                i;
  PASTIX_INT                procnum;
  int                retval     = NO_ERR;
  int                retval_rcv;

  print_debug(DBG_STEP,"-> pastix_task_scotch\n");

  if (sizeof(PASTIX_INT) != sizeof(SCOTCH_Num))
    {
      errorPrint("Inconsistent integer type\n");
      RETURN_ERROR(INTEGER_TYPE_ERR);
    }
  grafmesh  = &((*pastix_data)->grafmesh);
  ordedat   = &((*pastix_data)->ordedat);
  dgraph    = &((*pastix_data)->dgraph);
  ordemesh  = &((*pastix_data)->ordemesh);
  procnum   =   (*pastix_data)->procnum;
  n2        = &((*pastix_data)->n2);
  col2      = &((*pastix_data)->col2);
  row2      = &((*pastix_data)->row2);
  loc2glob2 = &((*pastix_data)->loc2glob2);

  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_onempi("%s", OUT_STEP_ORDER);

  /* On nettoie grafmesh et col2/row2 si ils sont d�ja allou�s */
#    ifdef WITH_SCOTCH
  if ((*pastix_data)->malgrf == 1)
    {
      SCOTCH_graphExit(grafmesh);
      (*pastix_data)->malgrf = 0;
    }
#    endif
  if ((*pastix_data)->bmalcolrow == 1)
    {
      if ((*col2)      != NULL) memFree_null(*col2);
      if ((*row2)      != NULL) memFree_null(*row2);
      if ((*loc2glob2) != NULL) memFree_null(*loc2glob2);
      (*pastix_data)->bmalcolrow = 0;
    }
  if (iparm[IPARM_ORDERING]  != API_ORDER_LOAD)
    {
      dpastix_order_prepare_cscd(*pastix_data,
                                 n,
                                 colptr,
                                 row,
                                 loc2glob,
                                 pastix_comm);
      /* Allocate a copy of col, row loc2glob */
      (*pastix_data)->bmalcolrow = 1;

    }

  switch (iparm[IPARM_ORDERING])
    {
      /*
       * Scotch Ordering
       */
    case API_ORDER_SCOTCH:
      /* construction du graphe */
      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphInit <\n");
      SCOTCH_dgraphInit(dgraph, pastix_comm);

      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphBuild <\n");
      if ( SCOTCH_dgraphBuild (dgraph,
                               1,              /* baseval */
                               *n2,            /* number of local vertices */
                               *n2,            /* Maximum number of local vertices     */
                               (*col2),
                               NULL,
                               NULL,           /* Local vertex load array (if any)     */
                               NULL,           /* Local vertex label array (if any)    */
                               (*col2)[(*n2)] - 1,
                               (*col2)[(*n2)] - 1,
                               (*row2),        /* Local edge array                     */
                               NULL,           /* Ghost edge array (if any); not const */
                               NULL))
        {
          errorPrint("SCOTCH_dgraphBuild");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }
      (*pastix_data)->malgrf = 1;

      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphCheck <\n");
      if (SCOTCH_dgraphCheck(dgraph))
        {
          errorPrint("pastix : SCOTCH_dgraphCheck");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

      print_debug(DBG_SCOTCH, "> SCOTCH_stratInit <\n");
      if (SCOTCH_stratInit(&stratdat))
        {
          errorPrint("pastix : SCOTCH_stratInit");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

      /* TODO : mettre des strategies par d�faut */
      if (iparm[IPARM_DEFAULT_ORDERING] == API_YES)
        {
          if (iparm[IPARM_INCOMPLETE] == API_NO)
            sprintf(strat, PTSCOTCH_STRAT_DIRECT);
          else
            sprintf(strat, PTSCOTCH_STRAT_INCOMP);
        }
      else /* personal ordering */
        {
          sprintf(strat, PTSCOTCH_STRAT_PERSO,
                  (long) iparm[IPARM_ORDERING_SWITCH_LEVEL],
                  (long) iparm[IPARM_ORDERING_CMIN],
                  (long) iparm[IPARM_ORDERING_CMAX],
                  ((float)iparm[IPARM_ORDERING_FRAT])/100.,
                  (long) iparm[IPARM_ORDERING_SWITCH_LEVEL],
                  (long) iparm[IPARM_ORDERING_CMIN],
                  (long) iparm[IPARM_ORDERING_CMAX],
                  ((float)iparm[IPARM_ORDERING_FRAT])/100.);

          if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
            print_onempi("Scotch Strategy |%s|\n", strat);
        }

      clockInit(&timer1);
      clockStart(&timer1);

      /*    print_debug(DBG_SCOTCH, "> SCOTCH_stratDgraphOrder <\n"); */
      /*    if (SCOTCH_stratDgraphOrder(&stratdat, strat)) */
      /*      { */
      /*        errorPrint("pastix : SCOTCH_stratDgraphOrder"); */
      /*        EXIT(MOD_SOPALIN,INTERNAL_ERR); */
      /*      } */
      if (procnum == 0 && iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        errorPrintW("PaStiX works only with PT-Scotch default strategy");

      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderInit <\n");
      if (0 != SCOTCH_dgraphOrderInit(dgraph, ordedat))
        {
          errorPrint("pastix : SCOTCH_dgraphOrderInit");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderCompute <\n");
      if (0 != SCOTCH_dgraphOrderCompute(dgraph, ordedat, &stratdat))
        {
          errorPrint("pastix : SCOTCH_dgraphOrderCompute");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

      print_debug(DBG_SCOTCH, "> SCOTCH_stratExit <\n");
      SCOTCH_stratExit(&stratdat);


      /* print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderPerm <\n"); */
      /*       if (0 != SCOTCH_dgraphOrderPerm(dgraph, ordedat, perm)) */
      /*  { */
      /*    errorPrint("pastix : SCOTCH_dgraphOrderPerm"); */
      /*    EXIT(MOD_SOPALIN,INTERNAL_ERR); */
      /*  } */

      clockStop(&timer1);
      if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        print_onempi(TIME_COMPUTE_ORDERING,clockVal(&timer1));


      /* Clean ordering if it exists */
      if ((*pastix_data)->malord)
        {
          orderExit(ordemesh);
          (*pastix_data)->malord=0;
        }

      if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        print_onempi("%s", OUT_ORDERINIT);
      orderInit(ordemesh);
      gN = 0;
      MPI_Allreduce(&n, &gN, 1, COMM_INT, MPI_SUM, pastix_comm);

      MALLOC_INTERN(ordemesh->rangtab, gN+1, PASTIX_INT);
      MALLOC_INTERN(ordemesh->permtab, gN,   PASTIX_INT);
      MALLOC_INTERN(ordemesh->peritab, gN,   PASTIX_INT);

      (*pastix_data)->malord=1;



      for (i=0;i<gN+1;i++)
        ordemesh->rangtab[i] = 0;

      SCOTCH_dgraphCorderInit (dgraph,
                               &ordering,
                               (SCOTCH_Num *)ordemesh->permtab,
                               (SCOTCH_Num *)ordemesh->peritab,
                               &ordemesh->cblknbr,
                               ordemesh->rangtab,
                               NULL);

      if (procnum == 0)
        {
          SCOTCH_dgraphOrderGather (dgraph, ordedat, &ordering);
        }
      else
        {
          SCOTCH_dgraphOrderGather (dgraph, ordedat, NULL);
        }
      MPI_Bcast(&ordemesh->cblknbr, 1                   , COMM_INT, 0, pastix_comm);
      MPI_Bcast(ordemesh->rangtab, (ordemesh->cblknbr+1), COMM_INT, 0, pastix_comm);
      MPI_Bcast(ordemesh->permtab, gN                   , COMM_INT, 0, pastix_comm);
      MPI_Bcast(ordemesh->peritab, gN                   , COMM_INT, 0, pastix_comm);

      global2localperm(n, perm, ((*pastix_data)->ordemesh).permtab, loc2glob);
      /* Gathering graph */
      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphGather <\n");
      SCOTCH_dgraphGather (dgraph, grafmesh);
      SCOTCH_dgraphCorderExit(dgraph, &ordering);
      SCOTCH_dgraphOrderExit(dgraph, ordedat);
      SCOTCH_dgraphExit(dgraph);

      SCOTCH_graphBase(grafmesh, 0);
      orderBase(ordemesh, 0);

      if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
        {
          PASTIX_INT nsave = 0;
          PASTIX_INT * colptrsave = NULL;
          PASTIX_INT * rowsave    = NULL;
          if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE_CSC))
            {
              cscd2csc_int((*pastix_data)->n2, (*pastix_data)->col2, (*pastix_data)->row2, NULL,
                           NULL, NULL, NULL,
                           &nsave, &colptrsave, &rowsave, NULL,
                           NULL, NULL, NULL,
                           (*pastix_data)->loc2glob2, pastix_comm, iparm[IPARM_DOF_NBR], API_YES);
            }

          retval = pastix_order_save(ordemesh,
                                     grafmesh,
                                     procnum,
                                     nsave,
                                     colptrsave,
                                     rowsave,
                                     iparm[IPARM_IO_STRATEGY]);
          if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE_CSC))
            {
              memFree_null(colptrsave);
              memFree_null(rowsave);
            }
          if (retval != NO_ERR)
            break;
        }


#    ifdef FORGET_PARTITION
      {
        PASTIX_INT iter;
        memFree_null(ordemesh->rangtab);
        MALLOC_INTERN(ordemesh->rangtab, n+1, PASTIX_INT);
        for (iter=0;iter<n+1;iter++)
          ordemesh->rangtab[iter]=iter;
        ordemesh->cblknbr=n;
      }
#    endif

      break;
    case API_ORDER_PERSONAL:
      {
        PASTIX_INT iter;
        PASTIX_INT * tmpperm = NULL;
        gN = 0;
        MPI_Allreduce(&n, &gN, 1, COMM_INT, MPI_SUM, pastix_comm);
        MALLOC_INTERN(tmpperm, gN, PASTIX_INT);
        /* Clean ordering if it exists */
        if ((*pastix_data)->malord)
          {
            orderExit(ordemesh);
            (*pastix_data)->malord=0;
          }

        if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
          print_onempi("%s", OUT_ORDERINIT);
        orderInit(ordemesh);
        gN = 0;
        MPI_Allreduce(&n, &gN, 1, COMM_INT, MPI_SUM, pastix_comm);

        MALLOC_INTERN(ordemesh->rangtab, gN+1, PASTIX_INT);
        memset(ordemesh->rangtab, 0,(gN+1)*sizeof(PASTIX_INT));
        MALLOC_INTERN(ordemesh->permtab, gN,   PASTIX_INT);
        MALLOC_INTERN(ordemesh->peritab, gN,   PASTIX_INT);

        (*pastix_data)->malord=1;

        memset(tmpperm, 0, gN*sizeof(PASTIX_INT));
        memset(ordemesh->permtab, 0, gN*sizeof(PASTIX_INT));
        for (iter = 0; iter < n; iter++)
          tmpperm[loc2glob[iter]-1] = perm[iter]-1;

        MPI_Allreduce(tmpperm, ordemesh->permtab, gN, COMM_INT, MPI_SUM, pastix_comm);
        for (iter = 0; iter < gN; iter++)
          ordemesh->peritab[ordemesh->permtab[iter]] = iter;

        for (iter=0; iter<gN+1; iter++)
          ordemesh->rangtab[iter] = iter;
        ordemesh->cblknbr = n;


        break;
      }
    case API_ORDER_LOAD:
      {
        PASTIX_INT   nload;
        PASTIX_INT * colptrload;
        PASTIX_INT * rowload;


        pastix_order_load(ordemesh,
                          grafmesh,
                          procnum,
                          &(nload),
                          &(colptrload),
                          &(rowload),
                          iparm[IPARM_IO_STRATEGY],
                          pastix_comm);
        break;
      }

    default:
      errorPrint("Ordering not available with distributed interface");
      retval = BADPARAMETER_ERR;
      break;
    }

  MPI_Allreduce(&retval, &retval_rcv, 1, MPI_INT, MPI_MAX, pastix_comm);
  if (retval_rcv != NO_ERR)
    RETURN_ERROR(retval_rcv);

  (*pastix_data)->malord=1;

  iparm[IPARM_START_TASK]++;
#  endif /* WITH_SCOTCH */
  return NO_ERR;
}
#endif /* DISTRIBUTED */

/*
  Function: pastix_task_fax

  Symbolic factorisation.

  Parameters:
  pastix_data - PaStiX data structure
  pastix_comm - PaStiX MPI communicator
  n           - Size of the matrix
  perm        - permutation tabular
  invp        - reverse permutation tabular
  flagWinvp   - flag to indicate if we have to print warning concerning perm and invp modification.

*/
void pastix_task_fax(pastix_data_t ** pastix_data, MPI_Comm pastix_comm,
                     PASTIX_INT * perm, PASTIX_INT * invp, int flagWinvp)
{
  PASTIX_INT           * iparm      =   (*pastix_data)->iparm;
  Order         * ordemesh   = &((*pastix_data)->ordemesh);
#ifdef WITH_SCOTCH
  SCOTCH_Graph  * grafmesh   = &((*pastix_data)->grafmesh);
#endif
  PASTIX_INT             procnum    =   (*pastix_data)->procnum;;
  FILE          * stream;
#ifdef DISTRIBUTED
  PASTIX_INT           * PTS_perm     = (*pastix_data)->PTS_permtab;
  PASTIX_INT           * PTS_rev_perm = (*pastix_data)->PTS_peritab;
  PASTIX_INT           * tmpperm      = NULL;
  PASTIX_INT           * tmpperi      = NULL;
  PASTIX_INT             gN;
  PASTIX_INT             i;
#endif
  PASTIX_INT             n;

  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_onempi("%s", OUT_STEP_FAX);

  print_debug(DBG_STEP,"->API_TASK_SYMBFACT\n");

  /* Force Load of symbmtx */
#ifdef ONLY_LOAD_SYMBMTX
  iparm[IPARM_IO_STRATEGY] = API_IO_LOAD;
#endif

  MALLOC_INTERN((*pastix_data)->symbmtx, 1, SymbolMatrix);
  /*
   * Load i/o strategy
   */
  if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD))
    {
      /* Load symbol & order */
      PASTIX_FOPEN(stream, "symbname","r");
      symbolLoad((*pastix_data)->symbmtx,stream);
      fclose(stream);
      symbolBase((*pastix_data)->symbmtx,0);
#ifdef WITH_SCOTCH
      if (ordemesh == NULL)
        {
          pastix_order_load(ordemesh,
                            grafmesh,
                            procnum,
                            &((*pastix_data)->n2),
                            &((*pastix_data)->col2),
                            &((*pastix_data)->row2),
                            iparm[IPARM_IO_STRATEGY],
                            pastix_comm);
        }
#endif
    }
  else /* not API_IO_LOAD */
    {
      n = ordemesh->rangtab[ordemesh->cblknbr];
#ifdef COMPACT_SMX
      if (iparm[IPARM_INCOMPLETE] == API_NO)
        {
          if (procnum == 0)
            errorPrintW("COMPACT_SMX only works with incomplete factorisation, forcing incomplete factorisation.");
          iparm[IPARM_INCOMPLETE] = API_YES;
        }
#endif /* COMPACT_SMX */

      if (iparm[IPARM_INCOMPLETE] == API_NO)
        {
          symbolInit((*pastix_data)->symbmtx);

          if ((iparm[IPARM_ORDERING] == API_ORDER_PERSONAL) ||
#ifdef METIS
              (iparm[IPARM_ORDERING] == API_ORDER_METIS)    ||
#endif
              0)
            {
              if ((procnum == 0) && (iparm[IPARM_LEVEL_OF_FILL] != -1))
                errorPrintW("metis or personal ordering can't be used without kass, forced use of kass.");
              iparm[IPARM_LEVEL_OF_FILL] = -1;
            }
#ifdef FORGET_PARTITION
          {
            if ((procnum == 0) && (iparm[IPARM_LEVEL_OF_FILL] != -1))
              errorPrintW("FORGET_PARTITION can't be used without kass, forced use of kass.");
            iparm[IPARM_LEVEL_OF_FILL] = -1;
          }
#endif

          if(iparm[IPARM_LEVEL_OF_FILL] == -1)
            {
              PASTIX_INT nkass;
              PASTIX_INT * colptrkass;
              PASTIX_INT * rowkass;

              if (iparm[IPARM_GRAPHDIST] == API_YES)
                {


                  CSC_sort((*pastix_data)->n2, (*pastix_data)->col2, (*pastix_data)->row2, NULL, 0);

                  cscd2csc_int((*pastix_data)->n2, (*pastix_data)->col2, (*pastix_data)->row2, NULL,
                               NULL, NULL, NULL,
                               &nkass, &colptrkass, &rowkass, NULL,
                               NULL, NULL, NULL,
#ifdef DISTRIBUTED
                               (*pastix_data)->loc2glob2,
#else
                               NULL,
#endif
                               (*pastix_data)->pastix_comm, iparm[IPARM_DOF_NBR], API_YES);
                  CSC_Fnum2Cnum(rowkass,
                                colptrkass,
                                nkass);
                }
              else
                {
                  nkass      = (*pastix_data)->n2;
                  colptrkass = (*pastix_data)->col2;
                  rowkass    = (*pastix_data)->row2;
                }

              kass(iparm[IPARM_LEVEL_OF_FILL],
                   iparm[IPARM_AMALGAMATION_LEVEL],
                   (*pastix_data)->symbmtx,
                   colptrkass[0],/* baseval*/
                   nkass,
                   colptrkass[nkass]-1,
                   colptrkass,
                   rowkass,
                   ordemesh,
                   pastix_comm);
              if (iparm[IPARM_GRAPHDIST] == API_YES)
                {
                  memFree_null(colptrkass);
                  memFree_null(rowkass);
                }
            }
          else /* iparm[IPARM_LEVEL_OF_FILL] != -1 */
            {
#ifdef WITH_SCOTCH
              symbolFaxGraph((*pastix_data)->symbmtx, grafmesh, ordemesh);
#endif
            }
          symbolBase((*pastix_data)->symbmtx,0);
        }
      else /* iparm[IPARM_INCOMPLETE] != API_NO */
        {
          PASTIX_INT nkass;
          PASTIX_INT * colptrkass;
          PASTIX_INT * rowkass;

          if (iparm[IPARM_GRAPHDIST] == API_YES)
            {

              CSC_sort((*pastix_data)->n2, (*pastix_data)->col2, (*pastix_data)->row2, NULL, 0);

              cscd2csc_int((*pastix_data)->n2, (*pastix_data)->col2, (*pastix_data)->row2, NULL,
                           NULL, NULL, NULL,
                           &nkass, &colptrkass, &rowkass, NULL,
                           NULL, NULL, NULL,
#ifdef DISTRIBUTED
                           (*pastix_data)->loc2glob2,
#else
                           NULL,
#endif
                           (*pastix_data)->pastix_comm, iparm[IPARM_DOF_NBR], API_YES);
              CSC_Fnum2Cnum(rowkass,
                            colptrkass,
                            nkass);

            }
          else
            {
              nkass = (*pastix_data)->n2;
              colptrkass = (*pastix_data)->col2;
              rowkass    = (*pastix_data)->row2;
            }

          symbolInit((*pastix_data)->symbmtx);
          /*bordi(iparm[IPARM_LEVEL_OF_FILL],&(solvmatr->symbmtx),grafmesh,ordemesh);*/
          kass(iparm[IPARM_LEVEL_OF_FILL],
               iparm[IPARM_AMALGAMATION_LEVEL],
               (*pastix_data)->symbmtx,
               colptrkass[0], /* baseval */
               nkass,
               colptrkass[nkass]-1,
               colptrkass,
               rowkass,
               ordemesh,
               pastix_comm);

          if (iparm[IPARM_GRAPHDIST] == API_YES)
            {
              memFree_null(colptrkass);
              memFree_null(rowkass);
            }
        }  /* iparm[IPARM_INCOMPLETE] != API_NO */
#ifdef DISTRIBUTED
      if (PTS_perm != NULL)
        {
          gN = n;

          MALLOC_INTERN(tmpperm, gN, PASTIX_INT);
          MALLOC_INTERN(tmpperi, gN, PASTIX_INT);
          for (i = 0; i < gN; i++)
            tmpperm[i] = ordemesh->permtab[PTS_perm[i]-1];

          memFree_null(ordemesh->permtab);
          ordemesh->permtab = tmpperm;

          for (i = 0; i < gN; i++)
            tmpperi[i] = PTS_rev_perm[ordemesh->peritab[i]]-1;
          memFree_null(ordemesh->peritab);
          ordemesh->peritab = tmpperi;

          memFree_null(PTS_perm);
          memFree_null(PTS_rev_perm);
        }
#endif /* DISTRIBUTED */
      /* WARNING : perm and invp can now be modified during symbolic factorization ??? */
      if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        if ( flagWinvp)
          if (procnum == 0)
            errorPrintW("perm and invp can be modified during symbolic factorization.");


      memcpy(perm, ordemesh->permtab, n*sizeof(PASTIX_INT));
      memcpy(invp, ordemesh->peritab, n*sizeof(PASTIX_INT));

      if ((*pastix_data)->bmalcolrow == 1)
        {
          if ((*pastix_data)->col2      != NULL) memFree_null((*pastix_data)->col2);
          if ((*pastix_data)->row2      != NULL) memFree_null((*pastix_data)->row2);
#ifdef DISTRIBUTED
          if ((*pastix_data)->loc2glob2 != NULL) memFree_null((*pastix_data)->loc2glob2);
#endif
          (*pastix_data)->bmalcolrow = 0;
        }
#ifdef WITH_SCOTCH
      if ((*pastix_data)->malgrf)
        {
          SCOTCH_graphExit(grafmesh);
          (*pastix_data)->malgrf=0;
        }
#endif
      if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE)) /* Save i/o strategy */
        {
          /* Save symbol */
          PASTIX_FOPEN(stream, "symbgen","w");
          symbolSave((*pastix_data)->symbmtx,stream);
          fclose(stream);
        }
    } /* not API_IO_LOAD */

#ifdef DUMP_SYMBOLMATRIX
  if ((*pastix_data)->procnum == 0)
    {
      PASTIX_FOPEN(stream, "symbol.eps", "w");
      symbolDraw((*pastix_data)->symbmtx,
                 stream);
      fclose(stream);
    }
#endif
  iparm[IPARM_START_TASK]++;

}

/*
  Function: dpastix_task_fax

  Symbolic factorisation.

  Parameters:
  pastix_data - PaStiX data structure
  pastix_comm - PaStiX MPI communicator
  n           - Size of the local matrix
  perm        - local permutation tabular
  loc2glob    - Global number of local columns (NULL if not ditributed)
  flagWinvp   - flag to indicate if we have to print warning concerning perm and invp modification.

*/
void dpastix_task_fax(pastix_data_t *pastix_data, MPI_Comm pastix_comm, PASTIX_INT n,
                      PASTIX_INT * perm, PASTIX_INT * loc2glob, int flagWinvp)
{
  PASTIX_INT   * gperm   = NULL;
  PASTIX_INT   * ginvp   = NULL;
  PASTIX_INT     gN;
  PASTIX_INT     my_n;
  PASTIX_INT   * my_perm = NULL;
  (void)pastix_comm;

  /* Note: for AUTOSPLIT_COMM
   *
   * perm is given by user, of size n,
   * we have to allocate it so that fax can write into it,
   * only written, data can be anything...
   *
   * loc2glob is also given by user and we need
   * a gathered one to build gperm from perm returned by fax.
   *
   * Anyway perm can't be used by user as it does not correspond to
   * the columns it gaves us. We just need to allocate data to write trash.
   */
  MPI_Allreduce(&n, &my_n, 1, COMM_INT, MPI_SUM,
                pastix_data->intra_node_comm);
  if (pastix_data->intra_node_procnum == 0)
    {
      if (my_n != n)
        {
          MALLOC_INTERN(my_perm, my_n, PASTIX_INT);
          perm = my_perm;
          if (pastix_data->procnum == 0)
            errorPrintW("User's perm array is invalid and will be unusable with IPARM_AUTOSPLIT_COMM");
        }

      if (!(PASTIX_MASK_ISTRUE(pastix_data->iparm[IPARM_IO_STRATEGY], API_IO_LOAD)))
        {
          gN = 0;
          MPI_Allreduce(&my_n, &gN, 1, COMM_INT, MPI_SUM, pastix_data->inter_node_comm);
          MALLOC_INTERN(gperm, gN, PASTIX_INT);
          MALLOC_INTERN(ginvp, gN, PASTIX_INT);
        }
      pastix_task_fax(&pastix_data, pastix_data->inter_node_comm,
                      gperm, ginvp, flagWinvp);

      if (my_n == n)
        {
          if (!(PASTIX_MASK_ISTRUE(pastix_data->iparm[IPARM_IO_STRATEGY], API_IO_LOAD)))
            {
              /* permtab may have been changed */
              global2localperm(my_n, perm, gperm, loc2glob);
            }
        }
      else
        {
          memFree_null(my_perm);
        }
      if (!(PASTIX_MASK_ISTRUE(pastix_data->iparm[IPARM_IO_STRATEGY], API_IO_LOAD))) {
        memFree_null(gperm);
        memFree_null(ginvp);
      }

    }
  else
    {
      if ( ( pastix_data->iparm[IPARM_INCOMPLETE] == API_NO &&
             ( pastix_data->iparm[IPARM_ORDERING] == API_ORDER_PERSONAL ||
               pastix_data->iparm[IPARM_ORDERING] == API_ORDER_METIS ||
               pastix_data->iparm[IPARM_LEVEL_OF_FILL] == -1) &&
             pastix_data->iparm[IPARM_GRAPHDIST] == API_YES ) ||
           ( pastix_data->iparm[IPARM_INCOMPLETE] == API_YES &&
             pastix_data->iparm[IPARM_GRAPHDIST] == API_YES ) )
        {
          PASTIX_INT nkass;
          PASTIX_INT * colptrkass;
          PASTIX_INT * rowkass;

          CSC_sort(pastix_data->n2, pastix_data->col2, pastix_data->row2, NULL, 0);

          cscd2csc_int(pastix_data->n2, pastix_data->col2, pastix_data->row2,
                       NULL,
                       NULL, NULL, NULL,
                       &nkass, &colptrkass, &rowkass, NULL,
                       NULL, NULL, NULL,
#ifdef DISTRIBUTED
                       pastix_data->loc2glob2,
#else
                       NULL,
#endif
                       pastix_data->pastix_comm, pastix_data->iparm[IPARM_DOF_NBR], API_YES);
          memFree_null(colptrkass);
          memFree_null(rowkass);
        }

      if (pastix_data->bmalcolrow == 1)
        {
          if (pastix_data->col2      != NULL) memFree_null(pastix_data->col2);
          if (pastix_data->row2      != NULL) memFree_null(pastix_data->row2);
#ifdef DISTRIBUTED
          if (pastix_data->loc2glob2 != NULL) memFree_null(pastix_data->loc2glob2);
#endif
          pastix_data->bmalcolrow = 0;
        }

    }
}


/*
  Function: pastix_task_blend

  Distribution task.

  Parameters:
  pastix_data - PaStiX data structure
  pastix_comm - PaStiX MPI communicator

*/
void pastix_task_blend(pastix_data_t **pastix_data,
                       MPI_Comm        pastix_comm)
{
  Dof            dofstr;
  BlendParam     blendpar;
  Assembly1D     assemb1D;
  Assembly2D     assemb2D;
  int            procnbr  = (*pastix_data)->inter_node_procnbr;
  PASTIX_INT           *iparm    = (*pastix_data)->iparm;
  double        *dparm    = (*pastix_data)->dparm;
  SolverMatrix  *solvmatr = &((*pastix_data)->solvmatr);
  PASTIX_INT            procnum  = (*pastix_data)->inter_node_procnum;;
  (void)pastix_comm;

  /* Si on refait blend on doit jeter nos ancien coefs */
  if ((*pastix_data)->malcof)
    {
      CoefMatrix_Free( &((*pastix_data)->sopar), solvmatr, iparm[IPARM_FACTORIZATION]);
      (*pastix_data)->malcof=0;
    }

  print_debug(DBG_STEP,"->API_TASK_ANALYSE\n");
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_onempi("%s", OUT_STEP_BLEND);

  dofInit(&dofstr);
  dofConstant(&dofstr, 0, (*pastix_data)->symbmtx->nodenbr,
              ((iparm[IPARM_DOF_COST] == 0)?iparm[IPARM_DOF_NBR]:iparm[IPARM_DOF_COST]));

  /*
    Warning: Never used
    nonzero=recursive_sum(0 ,solvmatr->symbmtx.cblknbr-1, nnz,
    &(solvmatr->symbmtx), &dofstr);
  */
  blendParamInit(&blendpar);

  blendpar.n         = (*pastix_data)->n2;
  blendpar.smpnbr    = iparm[IPARM_NB_SMP_NODE_USED];
  blendpar.ricar     = iparm[IPARM_INCOMPLETE];
  blendpar.abs       = iparm[IPARM_ABS];
  blendpar.blcolmin  = iparm[IPARM_MIN_BLOCKSIZE];
  blendpar.blcolmax  = iparm[IPARM_MAX_BLOCKSIZE];
  blendpar.blblokmin = iparm[IPARM_MIN_BLOCKSIZE];
  blendpar.blblokmax = iparm[IPARM_MAX_BLOCKSIZE];

  if(blendpar.blcolmin > blendpar.blcolmax)
    {
      errorPrint("Parameter error : blocksize max < blocksize min (cf. iparm.txt).");
      ASSERT(blendpar.blcolmin <=  blendpar.blcolmax, MOD_SOPALIN);
    }

  blendpar.level2D    = iparm[IPARM_DISTRIBUTION_LEVEL];
  blendpar.ratiolimit = (double)(iparm[IPARM_DISTRIBUTION_LEVEL]);

  if (blendpar.autolevel)
    printf("ratiolimit=%lf\n",blendpar.ratiolimit);
  else
    printf("level2D=%ld\n", (long) blendpar.level2D);
  printf("malt_limit=%ld\n", (long) blendpar.malt_limit);

#ifdef TRACE_SOPALIN
  if (procnum == 0)
    blendpar.tracegen = 1;
#endif

  blendpar.iparm = iparm;
  blendpar.dparm = dparm;

#ifdef FORCE_NOSMP
  iparm[IPARM_THREAD_NBR] = 1;
#endif
#ifdef STARPU_CONTEXT
  {
    compute_context_nbr(&(iparm[IPARM_STARPU_CTX_NBR]),
                        iparm[IPARM_THREAD_NBR],
                        iparm[IPARM_VERBOSE]);
    /* iparm[IPARM_THREAD_NBR] = iparm[IPARM_STARPU_CTX_NBR]-1; */
  }
#endif
  solverBlend(solvmatr, (*pastix_data)->symbmtx, &assemb1D,&assemb2D,
              procnbr, iparm[IPARM_THREAD_NBR], iparm[IPARM_CUDA_NBR],
              procnum,&blendpar,&dofstr);
  symbolExit((*pastix_data)->symbmtx);
  memFree_null((*pastix_data)->symbmtx);
  (*pastix_data)->malslv = 1;

  if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
    {
      iparm[IPARM_NNZEROS]       *=2;
      dparm[DPARM_PRED_FACT_TIME]*=2.;
    }
  dparm[DPARM_SOLV_FLOPS] = (double)iparm[IPARM_NNZEROS]; /* number of operations for solve */

  iparm[IPARM_NNZEROS_BLOCK_LOCAL] = solvmatr->coefnbr;

  /* Affichage */
  dparm[DPARM_FILL_IN]       = dparm[DPARM_FILL_IN]      *(double)(iparm[IPARM_NNZEROS]/(iparm[IPARM_DOF_NBR]*iparm[IPARM_DOF_NBR]));
  if ((procnum==0) && (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT))
    {
      fprintf(stdout,TIME_TO_ANALYSE,    (double)dparm[DPARM_ANALYZE_TIME]);
      fprintf(stdout, NNZERO_WITH_FILLIN_TH, (long)iparm[IPARM_NNZEROS]);
      fprintf(stdout, OUT_FILLIN_TH,         (double)dparm[DPARM_FILL_IN]);
      if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
        fprintf(stdout,NUMBER_OP_LU,     (double)dparm[DPARM_FACT_FLOPS]);
      else
        fprintf(stdout, NUMBER_OP_LLT,    (double)dparm[DPARM_FACT_FLOPS]);
      fprintf(stdout,TIME_FACT_PRED,PERF_MODEL,     (double)dparm[DPARM_PRED_FACT_TIME]);

    }
  if ((iparm[IPARM_VERBOSE] > API_VERBOSE_NO))
    {
      PASTIX_INT solversize = sizeofsolver(solvmatr, iparm);

      fprintf(stdout,SOLVMTX_WITHOUT_CO,  (int)procnum, (double)MEMORY_WRITE(solversize),
              MEMORY_UNIT_WRITE(solversize));
      fprintf(stdout, NNZERO_WITH_FILLIN, (int)procnum, (long)iparm[IPARM_NNZEROS_BLOCK_LOCAL]);
    }
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
    {
      PASTIX_INT sizeL = solvmatr->coefnbr;
      PASTIX_INT sizeG = 0;

      MPI_Reduce(&sizeL, &sizeG, 1, COMM_INT, MPI_MAX, 0, pastix_comm);

      if (procnum == 0)
        {
          sizeG *= sizeof(PASTIX_FLOAT);
          if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
            sizeG *= 2;

          fprintf(stdout, OUT_COEFSIZE, (double)MEMORY_WRITE(sizeG),
                  MEMORY_UNIT_WRITE(sizeG));
        }
    }


#if (defined X_ARCHpower_ibm_aix)
#  ifdef ENABLE_INTR
  mpc_enableintr();

  printf("mpc : enable=%ld\n",mpc_enableintr());
  printf("mpc : queryintrdelay=%ld query_intr=%ld\n",mpc_queryintrdelay(),
         mpc_queryintr());
#  endif
#endif /* (defined X_ARCHpower_ibm_aix) */

  iparm[IPARM_START_TASK]++;
}

/* Function: sopalin_check_param

   Check parameters consistency.

   Parameters:
   pastix_data - PaStiX data structure.

   Return:
   NO_ERR           - if no error occured
   BADPARAMETER_ERR - if Parameters are not correct on one proc.
*/
#define sopalin_check_param PASTIX_PREFIX_F(sopalin_check_param)
static inline
int sopalin_check_param(pastix_data_t *pastix_data)
{
  PASTIX_INT           * iparm    = pastix_data->iparm;
  int             ret      = NO_ERR;
  int             ret_recv = NO_ERR;

  if ((iparm[IPARM_SYM]           == API_SYM_NO) &&
      (iparm[IPARM_FACTORIZATION] != API_FACT_LU))
    {
      errorPrint("With unsymmetric patterns LU decomposition should be used");
      ret = BADPARAMETER_ERR;
    }
#ifndef TYPE_COMPLEX
  if (iparm[IPARM_FACTORIZATION] == API_FACT_LDLH)
    {
      errorPrint("LDLH only available with complex");
      ret = BADPARAMETER_ERR;
    }
#endif
  MPI_Allreduce(&ret, &ret_recv, 1, MPI_INT, MPI_MAX, pastix_data->inter_node_comm);
  RETURN_ERROR(ret_recv);
}


/* Function: pastix_check_param

   Check parameters consistency.

   Parameters:
   pastix_data - PaStiX data structure.

   Return:
   NO_ERR           - if no error occured
   BADPARAMETER_ERR - if Parameters are not correct on one proc.
*/
#define pastix_check_param PASTIX_PREFIX_F(pastix_check_param)
static inline
int pastix_check_param(pastix_data_t * pastix_data, int rhsnbr)
{
  int ret = NO_ERR, ret_recv;
  PASTIX_INT * iparm = pastix_data->iparm;

  if (iparm[IPARM_THREAD_NBR] < 1) {
    errorPrint("iparm[IPARM_THREAD_NBR] must be strictly positive.");
    ret = BADPARAMETER_ERR;
  }

  if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD_CSC) ||
      PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD_GRAPH))
    iparm[IPARM_IO_STRATEGY] |= API_IO_LOAD;

  if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE_CSC) ||
      PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE_GRAPH))
    iparm[IPARM_IO_STRATEGY] |= API_IO_SAVE;

  if ((PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD)) &&
      (iparm[IPARM_ORDERING]    != API_ORDER_LOAD))
    iparm[IPARM_ORDERING]  = API_ORDER_LOAD;

  if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD) &&
      PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
    {
      errorPrint("Save and load strategy are not compatible\n");
      ret = BADPARAMETER_ERR;
    }


#ifndef MULT_SMX
    if (iparm[IPARM_START_TASK] < API_TASK_CLEAN &&
        iparm[IPARM_END_TASK] > API_TASK_NUMFACT &&
        rhsnbr > 1)
      {
        errorPrint("Need to rebuild with -DMULT_SMX to use multiple RHS\n");
        ret = BADPARAMETER_ERR;
      }
#endif
    if (iparm[IPARM_START_TASK] < API_TASK_CLEAN &&
        iparm[IPARM_END_TASK] > API_TASK_NUMFACT &&
        rhsnbr < 1)
      {
        errorPrint("Number of RHS must be a positive number\n");
        ret = BADPARAMETER_ERR;
      }

    if (iparm[IPARM_STARPU] == API_YES)
      {
#ifndef WITH_STARPU
        errorPrintW("To use StarPU scheduler please build"
                    " PaStiX with -DWITH_STARPU");
        ret = BADPARAMETER_ERR;
#endif
        if (pastix_data->inter_node_procnbr > 1)
          {
            errorPrintW("StarPU does not works in MPI mode yet");
            ret = NOTIMPLEMENTED_ERR;
          }
      }

    if (iparm[IPARM_START_TASK] <= API_TASK_SOLVE  &&
        iparm[IPARM_END_TASK]   >= API_TASK_REFINE &&
        iparm[IPARM_ONLY_RAFF]  == API_YES) {
      errorPrintW("IPARM_ONLY_RAFF ignored, only possible if UPDO and RAFF are called in 2 steps");
      iparm[IPARM_ONLY_RAFF] = API_NO;
     }


    if (iparm[IPARM_FACTORIZATION] == API_FACT_LU &&
        iparm[IPARM_ONLY_RAFF]     == API_YES     &&
        iparm[IPARM_REFINEMENT]    == API_RAF_PIVOT)
      {
        errorPrint("Reffinement only is not compatible with the"
                   " simple iterative refinement.");
        ret = BADPARAMETER_ERR;
      }
  MPI_Allreduce(&ret, &ret_recv, 1, MPI_INT, MPI_MAX, pastix_data->pastix_comm);
  RETURN_ERROR(ret_recv);
}

#ifdef COMPUTE
/*
 * Function: pastix_fake_fillin_csc
 *
 * Fill in the internal csc based on the user csc and fill in the coeftab structure
 *
 * Parameters:
 * pastix_data - PaStiX data structure.
 * pastix_comm - PaStiX MPI communicator.
 * n           - Size of the matrix.
 * colptr      - starting index of each column in row and avals.
 * row         - row of each element of the matrix.
 * avals       - value of each element of the matrix.
 * b           - Right hand side.
 * nrhs        - Number of right-hand-sides.
 * loc2glob    - global number of local columns, NULL if not distributed.
 */
int pastix_fake_fillin_csc( pastix_data_t *pastix_data,
                            MPI_Comm       pastix_comm,
                            PASTIX_INT            n,
                            PASTIX_INT           *colptr,
                            PASTIX_INT           *row,
                            PASTIX_FLOAT         *avals,
                            PASTIX_FLOAT         *b,
                            PASTIX_INT            nrhs,
                            PASTIX_INT           *loc2glob)
{
  PASTIX_INT            *iparm    = pastix_data->iparm;
  PASTIX_INT            *l_colptr = NULL;
  PASTIX_INT            *l_row    = NULL;
  PASTIX_FLOAT          *l_val    = NULL;
  PASTIX_INT            *l_l2g    = NULL;
  PASTIX_FLOAT          *l_b      = NULL;
  PASTIX_INT             l_n      = 0;
  int             mal_l_l2g = API_NO;
  int             mal_l_b   = API_NO;
  (void)pastix_comm; (void)nrhs;
#  ifdef DISTRIBUTED
  int             OK       = 0;
  int             OK_RECV  = 0;
  int             retval = NO_ERR;
  int             retval_recv;
  PASTIX_INT             iter;
  PASTIX_INT             gN = -1;

  if (iparm[IPARM_GRAPHDIST] == API_YES)
    {
      /* If user has not specified that he is
         absolutly certain that is CSCd is
         correctly distributed */
      if (iparm[IPARM_CSCD_CORRECT] == API_NO)
        {
          /* Test que la cscd utilisateur correspond a la cscd pastix */
          l_n = 0;
          l_l2g = NULL;

          OK = 0;
          if (l_n != n)
            {
              OK = 1;
            }
          else
            {
              for (iter = 0; iter < l_n; iter++)
                {
                  if (l_l2g[iter] != loc2glob[iter])
                    {
                      OK = 1;
                      break;
                    }
                }
            }
          MPI_Allreduce(&OK, &OK_RECV, 1, MPI_INT, MPI_SUM, pastix_comm);
        }

      if (NULL != pastix_data->glob2loc)
        memFree_null(pastix_data->glob2loc);

      /* Building global to local column number correspondance */
      cscd_build_g2l(l_n,
                     l_l2g,
                     pastix_comm,
                     &gN,
                     &(pastix_data->glob2loc));

    }

  /* Si la cscd correspond pas, on la corrige */
  if (OK_RECV > 0)
    {
      if (pastix_data->procnum == 0 &&
          iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)
        fprintf(stdout,OUT_REDIS_CSC);

      /* redistributing cscd with correct local to global correspondance */
      retval = cscd_redispatch_int(  n,    colptr,    row,  avals,    b, nrhs, loc2glob,
                                     l_n, &l_colptr, &l_row, &l_val, &l_b, l_l2g,
                                     API_YES, pastix_comm, iparm[IPARM_DOF_NBR]);
      memFree_null(l_colptr); /* in fake fillin we do nothing with that */
      MPI_Allreduce(&retval, &retval_recv, 1, MPI_INT, MPI_MAX, pastix_comm);
      if (retval_recv != NO_ERR)
        RETURN_ERROR(retval_recv);

      mal_l_b = API_YES;
    }
  else
#  endif /* DISTRIBUTED */
    {
      /* the user CSCd is correctly distributed,
         only the pointers to the user CSCd arrays are used.
      */
      if ((iparm[IPARM_CSCD_CORRECT] == API_NO)
          && (l_l2g != NULL))
        memFree_null(l_l2g);

      l_n      = n;
      l_colptr = colptr;
      l_row    = row;
      l_val    = avals;
      l_l2g    = loc2glob;
      mal_l_l2g = API_NO;
      mal_l_b   = API_NO;
      l_b      = b;
    }

#  ifdef DISTRIBUTED
  if (pastix_data->mal_l2g_int == API_YES)
    memFree_null(pastix_data->l2g_int);
  if (pastix_data->malrhsd_int == API_YES)
    memFree_null(pastix_data->b_int);

  pastix_data->ncol_int = l_n;
  pastix_data->l2g_int  = l_l2g;
  pastix_data->b_int    = l_b;
  pastix_data->mal_l2g_int = mal_l_l2g;
  pastix_data->malrhsd_int = mal_l_b;
#  endif

  return NO_ERR;
}

/*
 *  Function: pastix_fillin_csc
 *
 *  Fill in the internal csc based on the user csc and fill in the coeftab structure
 *
 *  Parameters:
 *  pastix_data - PaStiX data structure.
 *  pastix_comm - PaStiX MPI communicator.
 *  n           - Size of the matrix.
 *  colptr      - starting index of each column in row and avals.
 *  row         - row of each element of the matrix.
 *  avals       - value of each element of the matrix.
 *  b           - Right hand side.
 *  nrhs        - Number of right-hand-sides.
 *  loc2glob    - global number of local columns, NULL if not distributed.
 */
int pastix_fillin_csc( pastix_data_t *pastix_data,
                       MPI_Comm       pastix_comm,
                       PASTIX_INT            n,
                       PASTIX_INT           *colptr,
                       PASTIX_INT           *row,
                       PASTIX_FLOAT         *avals,
                       PASTIX_FLOAT         *b,
                       PASTIX_INT            nrhs,
                       PASTIX_INT           *loc2glob)
{
  PASTIX_INT            *iparm    = pastix_data->iparm;
  SolverMatrix   *solvmatr = &(pastix_data->solvmatr);
  Order          *ordemesh = &(pastix_data->ordemesh);
  PASTIX_INT            *l_colptr = NULL;
  PASTIX_INT            *l_row    = NULL;
  PASTIX_FLOAT          *l_val    = NULL;
  PASTIX_INT            *l_l2g    = loc2glob;
  PASTIX_FLOAT          *l_b      = NULL;
  PASTIX_INT             l_n      = n;
  PASTIX_FLOAT         **transcsc = NULL;
  PASTIX_INT             procnum  = pastix_data->procnum;
  int             malcsc   = 0;
  int             forcetr  = 0;
  char            Type[4];
  Clock           clk;
  int             mal_l_l2g = API_NO;
  int             mal_l_b   = API_NO;
  (void)pastix_comm;
  (void)nrhs;
#  ifdef DISTRIBUTED
  int             OK       = 0;
  int             OK_RECV  = 0;
  int             retval = NO_ERR;
  int             retval_recv;
  PASTIX_INT             iter;
  PASTIX_INT             gN = -1;

  if (iparm[IPARM_GRAPHDIST] == API_YES)
    {
      /* If user has not specified that he is
         absolutly certain that is CSCd is
         correctly distributed */
      if (iparm[IPARM_CSCD_CORRECT] == API_NO)
        {
          /* Test que la cscd utilisateur correspond a la cscd pastix */
          l_n = pastix_getLocalNodeNbr(&pastix_data);
          MALLOC_INTERN(l_l2g, l_n, PASTIX_INT);
          if (l_l2g != NULL) mal_l_l2g = API_YES;
          pastix_getLocalNodeLst(&pastix_data, l_l2g);

          OK = 0;
          if (l_n != n)
            {
              OK = 1;
            }
          else
            {
              for (iter = 0; iter < l_n; iter++)
                {
                  if (l_l2g[iter] != loc2glob[iter])
                    {
                      OK = 1;
                      break;
                    }
                }
            }
          MPI_Allreduce(&OK, &OK_RECV, 1, MPI_INT, MPI_SUM, pastix_comm);
        }

      if (NULL != pastix_data->glob2loc)
        memFree_null(pastix_data->glob2loc);

      /* Building global to local column number correspondance */
      cscd_build_g2l(l_n,
                     l_l2g,
                     pastix_comm,
                     &gN,
                     &(pastix_data->glob2loc));

    }

  /* Si la cscd correspond pas, on la corrige */
  if (OK_RECV > 0)
    {
      if (procnum == 0 &&
          iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)
        fprintf(stdout,OUT_REDIS_CSC);

      /* redistributing cscd with correct local to global correspondance */
      clockInit(&clk);
      clockStart(&clk);
      retval = cscd_redispatch_int(  n,    colptr,    row,  avals,    b, nrhs, loc2glob,
                                     l_n, &l_colptr, &l_row, &l_val, &l_b, l_l2g,
                                     API_YES, pastix_comm, iparm[IPARM_DOF_NBR]);
      clockStop(&(clk));

      if (iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)
        print_onempi(OUT_REDISCSCDTIME, (double)clockVal(&clk));

      MPI_Allreduce(&retval, &retval_recv, 1, MPI_INT, MPI_MAX, pastix_comm);
      if (retval_recv != NO_ERR)
        RETURN_ERROR(retval_recv);

      malcsc = API_YES;
      mal_l_b = API_YES;
    }
  else
#  endif /* DISTRIBUTED */
    {
      /* the user CSCd is correctly distributed,
         only the pointers to the user CSCd arrays are used.
      */
      malcsc = API_NO;
      if ((iparm[IPARM_CSCD_CORRECT] == API_NO)
          && (l_l2g != NULL))
        memFree_null(l_l2g);

      l_n      = n;
      l_colptr = colptr;
      l_row    = row;
      l_val    = avals;
      l_l2g    = loc2glob;
      mal_l_l2g = API_NO;
      mal_l_b   = API_NO;
      l_b      = b;
    }

#  ifdef DISTRIBUTED
  if (pastix_data->mal_l2g_int == API_YES)
    memFree_null(pastix_data->l2g_int);
  if (pastix_data->malrhsd_int == API_YES)
    memFree_null(pastix_data->b_int);

  pastix_data->ncol_int = l_n;
  pastix_data->l2g_int  = l_l2g;
  pastix_data->b_int    = l_b;
  pastix_data->mal_l2g_int = mal_l_l2g;
  pastix_data->malrhsd_int = mal_l_b;
#  endif

  Type[0] = 'R';
  Type[1] = 'S';
  Type[2] = 'A';
  Type[3] = '\0';

  /* Fill in of the internal CSC */
  if (iparm[IPARM_FILL_MATRIX] == API_NO) /* not false factorisation */
    {
      if (pastix_data->malcsc)
        {
          CscExit(&(pastix_data->cscmtx));
          pastix_data->malcsc=0;
        }

      clockInit(&clk);
      clockStart(&clk);

      /* Choix des parametres pour CscOrdistrib */
      if (iparm[IPARM_SYM] == API_SYM_YES || iparm[IPARM_SYM] == API_SYM_HER) /* symmetric mtx */
        {
          if (iparm[IPARM_FACTORIZATION] == API_FACT_LU) /* LU */
            {
              /* Lu on RSA */
              printf("LU on RSA\n");
              forcetr = 1;
              Type[1] = 'S';
            }
          else
            {
              printf("LLt on RSA\n");
              forcetr = 0;
              Type[1] = 'S';
            }
          if (iparm[IPARM_SYM] == API_SYM_HER)
            Type[1] = 'H';
        }
      else
        {
          printf("LU on RUA\n");
          forcetr = 0;
          Type[1] = 'U';
        }
      transcsc = &(pastix_data->sopar.transcsc);

      if (iparm[IPARM_ONLY_RAFF] == API_YES)
        transcsc = NULL;

      /* Build internal CSCD from user CSC */
      if (iparm[IPARM_GRAPHDIST] == API_NO)
        {
          CscOrdistrib(&(pastix_data->cscmtx), Type,
                       transcsc, ordemesh,
                       l_n, l_n, l_colptr[l_n]-1, l_colptr,
                       l_row, l_val, forcetr,
                       solvmatr, procnum, iparm[IPARM_DOF_NBR]);
        }
#  ifdef DISTRIBUTED
      else
        {
          /* Build internal CSCD from user CSCD */
          CscdOrdistrib(&(pastix_data->cscmtx), Type,
                        transcsc, ordemesh,
                        l_n, l_colptr,
                        l_row, l_val,
                        l_l2g,
                        gN, pastix_data->glob2loc,
                        forcetr, solvmatr, procnum,
                        iparm[IPARM_DOF_NBR], pastix_data->inter_node_comm);
        }
#  endif /* DISTRIBUTED */
      pastix_data->malcsc = 1;

      /* Lib�ration de la csc interne temporaire */
      if (malcsc == API_YES)
        {
          /* l2g and rhs are not freed here beacause
             they will be used for solution filling */
          memFree_null(l_colptr);
          memFree_null(l_row);
          memFree_null(l_val);
        }

      clockStop(&(clk));

      if (iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)
        print_onempi(OUT_FILLCSCTIME, (double)clockVal(&clk));

      /* User Csc is useless after cscordistrib */
      if (iparm[IPARM_FREE_CSCUSER] == API_CSC_FREE)
        {
          free(colptr);
          free(row);
          free(avals);
          colptr = NULL;
          row    = NULL;
          avals  = NULL;
        }
    }

  if (pastix_data->malcof)
    {
      if (iparm[IPARM_SCHUR] == API_YES && pastix_data->schur_tab_set == API_YES)
        {
          SolverMatrix * datacode = &(pastix_data->solvmatr);
          PASTIX_INT            cblk;

          if (SOLV_TASKNBR > 0)
            {
              cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
              if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
                {
                  SOLV_COEFTAB(cblk) = NULL;
                }
            }
        }
      CoefMatrix_Free( &(pastix_data->sopar), solvmatr, iparm[IPARM_FACTORIZATION]);
      pastix_data->malcof=0;
    }

  /* On alloue par bloc colonne et pas tte la matrice */
  /* L'initialisation de la matrice est faite dans sopalin_init_smp              */
  pastix_data->sopar.iparm = iparm;
  if (iparm[IPARM_ONLY_RAFF] == API_NO)
    {
      if (iparm[IPARM_SCHUR] == API_YES && pastix_data->schur_tab_set == API_YES)
        {
          SolverMatrix * datacode = &(pastix_data->solvmatr);
          PASTIX_INT            cblk;

          if (SOLV_TASKNBR > 0)
            {
              cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
              if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
                {
                  SOLV_COEFTAB(cblk) = pastix_data->schur_tab;
                }
            }
        }

      pastix_data->malcof = 1;
    }

  return NO_ERR;
}
#else
#  define pastix_fillin_csc(ptx_data, comm, n, col, row, a, b, nrhs, l2g) NO_ERR
#endif /* COMPUTE */

/*
  Function: pastix_task_sopalin

  Factorisation, updown and raffinement tasks.

  Parameters:
  pastix_data - PaStiX data structure.
  pastix_comm - PaStiX MPI communicator.
  n           - Size of the matrix.
  colptr      - starting index of each column in row and avals.
  row         - row of each element of the matrix.
  avals       - value of each element of the matrix.
  b           - Right hand side.
  loc2glob    - global number of local columns, NULL if not distributed.
*/
int pastix_task_sopalin( pastix_data_t *pastix_data,
                         MPI_Comm       pastix_comm,
                         PASTIX_INT            n,
                         PASTIX_INT           *colptr,
                         PASTIX_INT           *row,
                         PASTIX_FLOAT         *avals,
                         PASTIX_FLOAT         *b,
                         PASTIX_INT            rhsnbr,
                         PASTIX_INT           *loc2glob)
{
  long            spivot, rpivot;
  double          sfacttime, rfacttime;
  double          ssolvtime,rsolvtime;
  double          srafftime,rrafftime;
  PASTIX_INT           * iparm    = pastix_data->iparm;
  double        * dparm    = pastix_data->dparm;
  SolverMatrix  * solvmatr = &(pastix_data->solvmatr);
  Order         * ordemesh = &(pastix_data->ordemesh);
  SopalinParam  * sopar    = &(pastix_data->sopar);
  PASTIX_INT             procnum  = pastix_data->inter_node_procnum;
  int             ret = 0;

  print_debug(DBG_STEP, "->API_TASK_NUMFACT\n");
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    {
      switch(iparm[IPARM_FACTORIZATION])
        {
        case API_FACT_LU:
          print_onempi("%s", OUT_STEP_NUMFACT_LU);
          break;
        case API_FACT_LLT:
          print_onempi("%s", OUT_STEP_NUMFACT_LLT);
          break;
        case API_FACT_LDLH:
          print_onempi("%s", OUT_STEP_NUMFACT_LDLH);
          break;
        case API_FACT_LDLT:
        default:
          print_onempi("%s", OUT_STEP_NUMFACT_LDLT);
        }
    }

  /* the error value has been reduced */
  if (NO_ERR != (ret = sopalin_check_param(pastix_data)))
    RETURN_ERROR(ret);

  /* Remplissage de la csc interne */
  if (pastix_data->cscInternFilled == API_NO) {
    pastix_fillin_csc(pastix_data, (pastix_data)->pastix_comm, n,
                      colptr, row, avals, b, rhsnbr, loc2glob);
  }
  /* sopalin */
  sopar->cscmtx      = &(pastix_data->cscmtx);
  sopar->itermax     = iparm[IPARM_ITERMAX];
  sopar->diagchange  = 0;
  sopar->epsilonraff = dparm[DPARM_EPSILON_REFINEMENT];
  sopar->rberror     = 0;
  sopar->espilondiag = dparm[DPARM_EPSILON_MAGN_CTRL];
  sopar->fakefact    = (iparm[IPARM_FILL_MATRIX] == API_YES) ? API_YES : API_NO;
  sopar->usenocsc    = 0;
  sopar->factotype   = iparm[IPARM_FACTORIZATION];
  sopar->symmetric   = iparm[IPARM_SYM];
  sopar->pastix_comm = pastix_comm;
  sopar->iparm       = iparm;
  sopar->dparm       = dparm;
  sopar->schur       = iparm[IPARM_SCHUR];
#ifdef DISTRIBUTED
  if (iparm[IPARM_GRAPHDIST] == API_YES) {
    sopar->n           = pastix_data->ncol_int;
    MPI_Allreduce(&(pastix_data->ncol_int), &sopar->gN, 1, COMM_INT, MPI_SUM, pastix_comm);
  } else
#endif
    {
    sopar->n           = n;
    sopar->gN          = n;
    }


  if (sopar->b != NULL)
    memFree_null(sopar->b);
  sopar->bindtab     = pastix_data->bindtab;

#ifdef FORCE_NOMPI
  iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
#elif (defined PASTIX_DYNSCHED)
  if ((iparm[IPARM_THREAD_COMM_MODE] != API_THREAD_COMM_ONE))
    {
      iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_COMM_ONE;
      if (procnum == 0 && iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        errorPrintW("Dynsched require API_THREAD_COMM_ONE, forced.");
    }
#endif

  switch (iparm[IPARM_THREAD_COMM_MODE])
    {
    case API_THREAD_COMM_ONE:
    case API_THREAD_FUNNELED:
      iparm[IPARM_NB_THREAD_COMM] = 1;

    case API_THREAD_COMM_DEFINED:
      iparm[IPARM_NB_THREAD_COMM] = MAX(iparm[IPARM_NB_THREAD_COMM],1);
      break;

    case API_THREAD_COMM_NBPROC:
      iparm[IPARM_NB_THREAD_COMM] = iparm[IPARM_THREAD_NBR];
      break;

    default:
      iparm[IPARM_NB_THREAD_COMM] = 0;
    }
  sopar->type_comm  = iparm[IPARM_THREAD_COMM_MODE];
  sopar->nbthrdcomm = iparm[IPARM_NB_THREAD_COMM];

  switch (iparm[IPARM_END_TASK])
    {
    case API_TASK_NUMFACT: /* Only sopalin */

      print_debug(DBG_STEP,"FACTO SEULE\n");

      /* no facto if only raff */
      if (iparm[IPARM_ONLY_RAFF] == API_NO)
        {
          switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
              ge_sopalin_thread(solvmatr, sopar);
              break;
            case API_FACT_LLT:
              po_sopalin_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLH:
              he_sopalin_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLT:
            default:
              sy_sopalin_thread(solvmatr, sopar);
            }
        }
      break;

    case API_TASK_SOLVE: /* Sopalin and updown */

      print_debug(DBG_STEP,"FACTO + UPDO\n");
 #ifndef PASTIX_FUNNELED
      if (THREAD_FUNNELED_ON)
        {
          if (procnum == 0)
            errorPrintW("API_THREAD_FUNNELED require -DPASTIX_FUNNELED,"
                        " force API_THREAD_MULTIPLE");
          sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
        }
#endif /* PASTIX_FUNNELED */
#ifndef PASTIX_THREAD_COMM
      if (THREAD_COMM_ON)
        {
          if (procnum == 0)
            errorPrintW("API_THREAD_COMM_* require -DPASTIX_THREAD_COMM,"
                        " force API_THREAD_MULTIPLE");
          sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
        }
#endif /* PASTIX_THREAD_COMM */

      /* Pour l'instant uniquement si on est en 1d */
      if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
        {
          /* attention MRHS_ALLOC */
          if (pastix_data->malsmx)
            {
              memFree_null(solvmatr->updovct.sm2xtab);
              pastix_data->malsmx=0;
            }
          solvmatr->updovct.sm2xnbr = rhsnbr;
          MALLOC_INTERN(solvmatr->updovct.sm2xtab,
                        solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                        PASTIX_FLOAT);

          pastix_data->malsmx=1;

          buildUpdoVect(pastix_data,
#ifdef DISTRIBUTED
                        pastix_data->l2g_int,
                        pastix_data->b_int,
#else
                        NULL,
                        b,
#endif
                        pastix_comm);
        }

      if (iparm[IPARM_ONLY_RAFF] == API_NO)
        {
          /* Pour l'instant uniquement si on est en 1d */
          if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
            {
              /* setting sopar->b for reffinement */
              /* Only 1 rhs is saved in sopar->b */
              if (sopar->b == NULL)
                {
                  MALLOC_INTERN(sopar->b, solvmatr->updovct.sm2xsze, PASTIX_FLOAT);
                }
              memcpy(sopar->b, solvmatr->updovct.sm2xtab,
                     solvmatr->updovct.sm2xsze*sizeof(PASTIX_FLOAT));
            }

          switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
              ge_sopalin_updo_thread(solvmatr, sopar);
              break;
            case API_FACT_LLT:
              po_sopalin_updo_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLH:
              he_sopalin_updo_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLT:
            default:
              sy_sopalin_updo_thread(solvmatr, sopar);
            }
        }

      iparm[IPARM_START_TASK]++;
      break;
    case API_TASK_REFINE: /* Sopalin, updown and raff */
    case API_TASK_CLEAN:

      print_debug(DBG_STEP,"FACTO + UPDO + RAFF (+ CLEAN)\n");
 #ifndef PASTIX_UPDO_ISEND
      if (THREAD_COMM_ON)
        {
          if (procnum == 0)
            errorPrintW("THREAD_COMM require -DPASTIX_UPDO_ISEND,"
                        " force API_THREAD_MULTIPLE");
          sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
        }
#endif /* PASTIX_UPDO_ISEND */
#ifndef STORAGE
      if (THREAD_COMM_ON)
        {
          if (procnum == 0)
            errorPrintW("THREAD_COMM require -DSTORAGE,"
                        " force API_THREAD_MULTIPLE");
          sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
        }
#endif /* STORAGE */

      /* Pour l'instant uniquement si on est en 1d */
      if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
        {
          /* attention MRHS_ALLOC */
          if (pastix_data->malsmx)
            {
              memFree_null(solvmatr->updovct.sm2xtab);
              pastix_data->malsmx=0;
            }
          if (rhsnbr > 1)
            errorPrintW("Reffinement works only with 1 rhs, please call them one after the other.");
          solvmatr->updovct.sm2xnbr = 1;
          MALLOC_INTERN(solvmatr->updovct.sm2xtab,
                        solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                        PASTIX_FLOAT);
          pastix_data->malsmx=1;

          buildUpdoVect(pastix_data,
#ifdef DISTRIBUTED
                        pastix_data->l2g_int,
                        pastix_data->b_int,
#else
                        NULL,
                        b,
#endif
                        pastix_comm);

          /* setting sopar->b for reffinement */
          if (sopar->b == NULL)
            {
              MALLOC_INTERN(sopar->b,
                            solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                            PASTIX_FLOAT);
            }
          memcpy(sopar->b, solvmatr->updovct.sm2xtab,
                 solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze*sizeof(PASTIX_FLOAT));
        }

      sopar->itermax     = iparm[IPARM_ITERMAX];
      sopar->epsilonraff = dparm[DPARM_EPSILON_REFINEMENT];
#ifdef OOC
      if (iparm[IPARM_GMRES_IM] != 1)
        {
          iparm[IPARM_GMRES_IM] = 1;
          if (procnum == 0)
            errorPrintW("IPARM_GMRES_IM force to 1 when using OOC");
        }
#endif
      sopar->gmresim = iparm[IPARM_GMRES_IM];

      switch (iparm[IPARM_REFINEMENT])
        {
        case API_RAF_GMRES:
          switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
              ge_sopalin_updo_gmres_thread(solvmatr, sopar);
              break;
            case API_FACT_LLT:
              po_sopalin_updo_gmres_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLH:
              he_sopalin_updo_gmres_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLT:
              sy_sopalin_updo_gmres_thread(solvmatr, sopar);
              break;
            default:
              errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
              return BADPARAMETER_ERR;
            }
          break;
        case API_RAF_PIVOT:
          switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
              ge_sopalin_updo_pivot_thread(solvmatr, sopar);
              break;
            case API_FACT_LLT:
            case API_FACT_LDLH:
            case API_FACT_LDLT:
              errorPrint("Refinement method and factorization type are incompatibles");
              return BADPARAMETER_ERR;
            default:
              errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
              return BADPARAMETER_ERR;
            }
          break;
        case API_RAF_GRAD:
          switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
              errorPrint("Refinement method and factorization type are incompatibles");
              return BADPARAMETER_ERR;
            case API_FACT_LLT:
              po_sopalin_updo_grad_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLH:
              he_sopalin_updo_grad_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLT:
              sy_sopalin_updo_grad_thread(solvmatr, sopar);
              break;
            default:
              errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
              return BADPARAMETER_ERR;
            }
          break;
        case API_RAF_BICGSTAB:
          switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
              ge_sopalin_updo_bicgstab_thread(solvmatr, sopar);
              break;
            case API_FACT_LLT:
            case API_FACT_LDLH:
            case API_FACT_LDLT:
              errorPrint("Refinement method and factorization type are incompatibles");
              return BADPARAMETER_ERR;
            default:
              errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
              return BADPARAMETER_ERR;
            }
          break;
        default:
          errorPrint("Undefined refinement method : %ld", (long)iparm[IPARM_REFINEMENT]);
          return BADPARAMETER_ERR;
        }

      /* sopar->b was only needed for raff */
      memFree_null(sopar->b);
      iparm[IPARM_START_TASK]++;
      iparm[IPARM_START_TASK]++;
      iparm[IPARM_NBITER]         = sopar->itermax;
      dparm[DPARM_RELATIVE_ERROR] = sopar->rberror;
    }

  if ((iparm[IPARM_END_TASK] > API_TASK_NUMFACT) /* Not only sopalin */
      && (iparm[IPARM_DISTRIBUTION_LEVEL] == 0))
    {
      /* b <- solution */
      if (iparm[IPARM_GRAPHDIST] == API_NO)
        {
          if (iparm[IPARM_ONLY_RAFF] == API_NO)
            {
              CscRhsUpdown(&(solvmatr->updovct),
                           solvmatr,
                           b, n, ordemesh->peritab,
                           iparm[IPARM_DOF_NBR],
                           iparm[IPARM_RHS_MAKING],
                           pastix_comm);
            }
        }
#ifdef DISTRIBUTED
      else
        {
          CscdRhsUpdown(&(solvmatr->updovct),
                        solvmatr,
                        pastix_data->b_int,
                        pastix_data->ncol_int,
                        pastix_data->glob2loc,
                        ordemesh->peritab,
                        (int)iparm[IPARM_DOF_NBR],
                        pastix_comm);
        }
#endif
    }

  iparm[IPARM_STATIC_PIVOTING] = sopar->diagchange;

  /*
   * Memory statistics
   */
#ifdef MEMORY_USAGE
  {
    unsigned long smem[2], rmem[2];

    smem[0] = memAllocGetMax();
    smem[1] = memAllocGetCurrent();

    dparm[DPARM_MEM_MAX] = (double)smem[0];
    MPI_Reduce(smem, rmem, 2, MPI_LONG, MPI_MAX, 0, pastix_comm);

    if (procnum == 0)
      {
        dparm[DPARM_MEM_MAX] = (double)rmem[0];

        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
          {
            fprintf(stdout, OUT_MAX_MEM_AF_SOP,  MEMORY_WRITE(rmem[0]), MEMORY_UNIT_WRITE(rmem[0]));
            fprintf(stdout, OUT_MEM_USED_AF_SOP, MEMORY_WRITE(rmem[1]), MEMORY_UNIT_WRITE(rmem[1]));
          }
      }
  }
#endif /* MEMORY_USAGE */


  spivot    = (long)  iparm[IPARM_STATIC_PIVOTING];
  sfacttime = (double)dparm[DPARM_FACT_TIME];
  MPI_Reduce(&spivot,   &rpivot,   1,MPI_LONG,  MPI_SUM,0,pastix_comm);
  MPI_Reduce(&sfacttime,&rfacttime,1,MPI_DOUBLE,MPI_MAX,0,pastix_comm);

  if (iparm[IPARM_ONLY_RAFF] == API_NO)
    {
      /*
       * Factorization Time
       */
      if ((procnum == 0) && (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT))
        {
          fprintf(stdout, OUT_STATIC_PIVOTING, rpivot);
          if (sopar->iparm[IPARM_INERTIA] != -1)
            {
              if (rpivot == 0)
                {
                  fprintf(stdout, OUT_INERTIA, (long)sopar->iparm[IPARM_INERTIA]);
                }
              else
                {
                  fprintf(stdout, OUT_INERTIA_PIVOT, (long)sopar->iparm[IPARM_INERTIA]);
                }
            }
          if (sopar->iparm[IPARM_ESP_NBTASKS] != -1)
            fprintf(stdout, OUT_ESP_NBTASKS,     (long)sopar->iparm[IPARM_ESP_NBTASKS]);
          fprintf(stdout, OUT_TIME_FACT,       rfacttime);
          fprintf(stdout, OUT_FLOPS_FACT,
                  PRINT_FLOPS((dparm[DPARM_FACT_FLOPS]/rfacttime)),
                  PRINT_FLOPS_UNIT((dparm[DPARM_FACT_FLOPS]/rfacttime)));
        }

      /*fprintf(stdout," Terms allocated during factorization %ld\n",(long)iparm[IPARM_ALLOCATED_TERMS]);*/

      /*
       * Solve Time
       */
      if (iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
        {
          ssolvtime = dparm[DPARM_SOLV_TIME];
          MPI_Reduce(&ssolvtime,&rsolvtime,1,MPI_DOUBLE,MPI_MAX,0,pastix_comm);

          if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
            print_onempi(OUT_TIME_SOLV, rsolvtime);
        }

      /*
       * Refinement Time
       */
      if (iparm[IPARM_END_TASK] > API_TASK_SOLVE)
        {
          srafftime = dparm[DPARM_RAFF_TIME];
          MPI_Reduce(&srafftime, &rrafftime, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);

          if ((procnum == 0) && (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT))
            {
              fprintf(stdout, OUT_RAFF_ITER_NORM,
                      (long)  iparm[IPARM_NBITER],
                      (double)dparm[DPARM_RELATIVE_ERROR]);
              fprintf(stdout, OUT_TIME_RAFF, rrafftime);
              if (iparm[IPARM_PRODUCE_STATS] == API_YES) {
                if (dparm[DPARM_RELATIVE_ERROR] > 0)
                  print_onempi(OUT_PREC1, dparm[DPARM_RELATIVE_ERROR]);
                if (dparm[DPARM_SCALED_RESIDUAL] > 0)
                print_onempi(OUT_PREC2, dparm[DPARM_SCALED_RESIDUAL]);
              }

            }
        }
    }

  printf("FACTO FIN\n");
  iparm[IPARM_START_TASK]++;
  return NO_ERR;
}



/*
  Function: pastix_task_updown

  Updown task.

  Parameters:
  pastix_data - PaStiX data structure.
  pastix_comm - PaStiX MPI communicator.
  n           - Matrix size.
  b           - Right hand side.
  loc2glob    - local to global column number.

*/
void pastix_task_updown(pastix_data_t *pastix_data,
                        MPI_Comm       pastix_comm,
                        PASTIX_INT            n,
                        PASTIX_FLOAT         *b,
                        PASTIX_INT            rhsnbr,
                        PASTIX_INT           *loc2glob)
{
  PASTIX_INT           * iparm    = pastix_data->iparm;
  double        * dparm    = pastix_data->dparm;
  SolverMatrix  * solvmatr = &(pastix_data->solvmatr);
  SopalinParam  * sopar    = &(pastix_data->sopar);
  Order         * ordemesh = &(pastix_data->ordemesh);
  PASTIX_INT             procnum  = pastix_data->procnum;
  double          ssolvtime,rsolvtime;

  print_debug(DBG_STEP,"-> pastix_task_updown\n");
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_onempi("%s", OUT_STEP_SOLVE);

  if (sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
      if (procnum == 0)
        errorPrintW("Updown step incompatible with 2D distribution");
      return;
    }

 #ifndef PASTIX_UPDO_ISEND
  if (THREAD_COMM_ON)
    {
      if (procnum == 0)
        errorPrintW("THREAD_COMM require -DPASTIX_UPDO_ISEND,"
                    " force API_THREAD_MULTIPLE");
      sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
    }
#endif /* PASTIX_UPDO_ISEND */
#ifndef STORAGE
  if (THREAD_COMM_ON)
    {
      if (procnum == 0)
        errorPrintW("THREAD_COMM require -DSTORAGE,"
                    " force API_THREAD_MULTIPLE");
      sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
    }
#endif /* STORAGE */

  printf("UPDO DEBUT\n");

  if ((iparm[IPARM_ONLY_RAFF] == API_YES) && (iparm[IPARM_END_TASK] > API_TASK_SOLVE))
    {
      errorPrintW("IPARM_ONLY_RAFF ignored, only possible if UPDO and RAFF are called in 2 steps");
      iparm[IPARM_ONLY_RAFF] = API_NO;
    }

  /* attention MRHS_ALLOC */
  if (pastix_data->malsmx)
    {
      memFree_null(solvmatr->updovct.sm2xtab);
      pastix_data->malsmx=0;
    }

  solvmatr->updovct.sm2xnbr = rhsnbr;
  MALLOC_INTERN(solvmatr->updovct.sm2xtab,
                solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                PASTIX_FLOAT);

  pastix_data->malsmx=1;

  buildUpdoVect(pastix_data,
                loc2glob,
                b,
                pastix_comm);

  if (iparm[IPARM_ONLY_RAFF] == API_NO)
    {
      /* setting sopar->b for reffinement */
      /* Only 1 rhs is saved in sopar->b */
      if (sopar->b == NULL)
        {
          MALLOC_INTERN(sopar->b, solvmatr->updovct.sm2xsze, PASTIX_FLOAT);
        }
      memcpy(sopar->b, solvmatr->updovct.sm2xtab,
             solvmatr->updovct.sm2xsze*sizeof(PASTIX_FLOAT));

      sopar->iparm = iparm;
      sopar->dparm = dparm;

      switch(iparm[IPARM_FACTORIZATION])
        {
        case API_FACT_LU:
          ge_updo_thread(solvmatr, sopar);
          break;
        case API_FACT_LLT:
          po_updo_thread(solvmatr, sopar);
          break;
        case API_FACT_LDLH:
          he_updo_thread(solvmatr, sopar);
          break;
        case API_FACT_LDLT:
        default:
          sy_updo_thread(solvmatr, sopar);
        }

      /*
        if ((procnum == 0) && (iparm[IPARM_END_TASK] < API_TASK_REFINE))
        errorPrintW("Need a call to step 6 (refinement) to put the solution in the user vector.");
      */
      if ((iparm[IPARM_END_TASK] < API_TASK_REFINE ||
           iparm[IPARM_TRANSPOSE_SOLVE] == API_YES)
          && (iparm[IPARM_DISTRIBUTION_LEVEL] == 0))
        {
          /* b <- solution */
          if (iparm[IPARM_GRAPHDIST] == API_NO)
            {
              CscRhsUpdown(&(solvmatr->updovct),
                           solvmatr,
                           b, n, ordemesh->peritab,
                           iparm[IPARM_DOF_NBR],
                           iparm[IPARM_RHS_MAKING],
                           pastix_comm);
            }
#ifdef DISTRIBUTED
          else
            {
              CscdRhsUpdown(&(solvmatr->updovct),
                            solvmatr,
                            b, n,
                            pastix_data->glob2loc,
                            ordemesh->peritab,
                            iparm[IPARM_DOF_NBR], pastix_comm);
            }
#endif
        }

      if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
        ssolvtime = (double)dparm[DPARM_SOLV_TIME];
        MPI_Reduce(&ssolvtime,&rsolvtime,1,MPI_DOUBLE,MPI_MAX,0,pastix_comm);
        print_onempi(OUT_TIME_SOLV,rsolvtime);
        if (iparm[IPARM_PRODUCE_STATS] == API_YES) {
          if (dparm[DPARM_RELATIVE_ERROR] > 0)
            print_onempi(OUT_PREC1, dparm[DPARM_RELATIVE_ERROR]);
          if (dparm[DPARM_SCALED_RESIDUAL] > 0)
            print_onempi(OUT_PREC2, dparm[DPARM_SCALED_RESIDUAL]);
        }
      }
    }
  printf("UPDO FIN\n");
  iparm[IPARM_START_TASK]++;
}

/*
  Function: pastix_task_raff

  Reffinement task

  Parameters:
  pastix_data - PaStiX data structure.
  pastix_comm - PaStiX MPI communicator.
  n           - Matrix size.
  b           - Right hand side.
  loc2glob    - local to global column number.
*/
void pastix_task_raff(pastix_data_t *pastix_data,
                      MPI_Comm       pastix_comm,
                      PASTIX_INT            n,
                      PASTIX_FLOAT         *b,
                      PASTIX_INT            rhsnbr,
                      PASTIX_INT           *loc2glob)
{
  PASTIX_INT           * iparm    = pastix_data->iparm;
  double        * dparm    = pastix_data->dparm;
  SopalinParam  * sopar    = &(pastix_data->sopar);
  SolverMatrix  * solvmatr = &(pastix_data->solvmatr);
  Order         * ordemesh = &(pastix_data->ordemesh);
  double          srafftime,rrafftime;
  PASTIX_INT             procnum  = pastix_data->procnum;;
  PASTIX_FLOAT         * tmp;

  print_debug(DBG_STEP, "->pastix_task_raff\n");
 #ifndef PASTIX_UPDO_ISEND
  if (THREAD_COMM_ON)
    {
      if (procnum == 0)
        errorPrintW("THREAD_COMM require -DPASTIX_UPDO_ISEND,"
                    " force API_THREAD_MULTIPLE");
      sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
    }
#endif /* PASTIX_UPDO_ISEND */
#ifndef STORAGE
  if (THREAD_COMM_ON)
    {
      if (procnum == 0)
        errorPrintW("THREAD_COMM require -DSTORAGE,"
                    " force API_THREAD_MULTIPLE");
      sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
    }
#endif /* STORAGE */

  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_onempi("%s", OUT_STEP_REFF);

  if (sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
      if (procnum == 0)
        errorPrintW("Refinment step incompatible with 2D distribution");
      return;
    }

  if (rhsnbr > 1)
    {
      errorPrintW("Reffinement works only with 1 rhs, please call them one after the other.");
      solvmatr->updovct.sm2xnbr = 1;
    }

  if (iparm[IPARM_ONLY_RAFF] == API_YES )
    {

      /* setting sopar->b for reffinement */
      if (sopar->b == NULL)
        {
          MALLOC_INTERN(sopar->b,
                        solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                        PASTIX_FLOAT);
        }

      tmp = solvmatr->updovct.sm2xtab;
      solvmatr->updovct.sm2xtab = sopar->b;

      buildUpdoVect(pastix_data,
                    loc2glob,
                    b,
                    pastix_comm);

      sopar->b = solvmatr->updovct.sm2xtab;
      solvmatr->updovct.sm2xtab = tmp;

    }

  sopar->itermax     = iparm[IPARM_ITERMAX];
  sopar->epsilonraff = dparm[DPARM_EPSILON_REFINEMENT];
#ifdef OOC
  if (iparm[IPARM_GMRES_IM] != 1)
    {
      iparm[IPARM_GMRES_IM] = 1;
      if (procnum == 0)
        errorPrintW("IPARM_GMRES_IM force to 1 when using OOC");
    }
#endif
  sopar->gmresim = iparm[IPARM_GMRES_IM];

      switch (iparm[IPARM_REFINEMENT])
        {
        case API_RAF_GMRES:
          switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
              ge_gmres_thread(solvmatr, sopar);
              break;
            case API_FACT_LLT:
              po_gmres_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLH:
              he_gmres_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLT:
              sy_gmres_thread(solvmatr, sopar);
              break;
            default:
              errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
              iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
              return;

            }
          break;
        case API_RAF_PIVOT:
          switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
              ge_pivot_thread(solvmatr, sopar);
              break;
            case API_FACT_LLT:
            case API_FACT_LDLH:
            case API_FACT_LDLT:
              errorPrint("Refinement method and factorization type are incompatibles");
              iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
              return;
            default:
              errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
              iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
              return;
            }
          break;
        case API_RAF_GRAD:
          switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
              errorPrint("Refinement method and factorization type are incompatibles");
              iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
              return;
            case API_FACT_LLT:
              po_grad_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLH:
              he_grad_thread(solvmatr, sopar);
              break;
            case API_FACT_LDLT:
              sy_grad_thread(solvmatr, sopar);
              break;
            default:
              errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
              iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
              return;
            }
          break;
        case API_RAF_BICGSTAB:
          switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
              ge_bicgstab_thread(solvmatr, sopar);
              break;
            case API_FACT_LLT:
            case API_FACT_LDLH:
            case API_FACT_LDLT:
              errorPrint("Refinement method and factorization type are incompatibles");
              iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
              return;
            default:
              errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
              iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
              return;
            }
          break;
        default:
          errorPrint("Undefined refinement method : %ld", (long)iparm[IPARM_REFINEMENT]);
          iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
          return;
        }

  dparm[DPARM_RELATIVE_ERROR] = sopar->rberror;
  iparm[IPARM_NBITER]         = sopar->itermax;

  /* sopar->b was only needed for raff */
  memFree_null(sopar->b);

  /* b <- solution */
  if (iparm[IPARM_GRAPHDIST] == API_NO)
    {
      CscRhsUpdown(&(solvmatr->updovct),
                   solvmatr,
                   b, n, ordemesh->peritab,
                   iparm[IPARM_DOF_NBR],
                   iparm[IPARM_RHS_MAKING],
                   pastix_comm);
    }
#ifdef DISTRIBUTED
  else
    {
      CscdRhsUpdown(&(solvmatr->updovct),
                    solvmatr,
                    b, n,
                    pastix_data->glob2loc,
                    ordemesh->peritab,
                    iparm[IPARM_DOF_NBR], pastix_comm);
    }
#endif

  /* Fin du roerdering */

  srafftime = (double)dparm[DPARM_RAFF_TIME];
  MPI_Reduce(&srafftime,&rrafftime,1,MPI_DOUBLE,MPI_MAX,0,pastix_comm);

  if ((procnum == 0) && (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT))
    {
      fprintf(stdout, OUT_RAFF_ITER_NORM, (long)iparm[IPARM_NBITER], (double)dparm[DPARM_RELATIVE_ERROR]);
      if (iparm[IPARM_PRODUCE_STATS] == API_YES) {
        if (dparm[DPARM_RELATIVE_ERROR] > 0)
          print_onempi(OUT_PREC1, dparm[DPARM_RELATIVE_ERROR]);
        if (dparm[DPARM_SCALED_RESIDUAL] > 0)
          print_onempi(OUT_PREC2, dparm[DPARM_SCALED_RESIDUAL]);
      }

      fprintf(stdout, OUT_TIME_RAFF, rrafftime);
    }
  printf("RAFF FIN\n");
  iparm[IPARM_START_TASK]++;

  return;
}

/*
 * Function: pastix_task_clean
 *
 * Cleaning task
 *
 * Parameters:
 *
 */
void pastix_task_clean(pastix_data_t **pastix_data,
                       MPI_Comm        pastix_comm)
{
  PASTIX_INT             i;
  PASTIX_INT           * iparm    = (*pastix_data)->iparm;
  Order         * ordemesh = &((*pastix_data)->ordemesh);
  SopalinParam  * sopar    = &((*pastix_data)->sopar);
  SolverMatrix  * solvmatr = &((*pastix_data)->solvmatr);
#ifdef PASTIX_DEBUG
  int             procnum  = (*pastix_data)->procnum;
  double        * dparm    = (*pastix_data)->dparm;
  FILE          * stream;
#endif
  (void)pastix_comm;

  print_debug(DBG_STEP, "->pastix_task_clean\n");
#ifdef DISTRIBUTED
  if ((*pastix_data)->mal_l2g_int == API_YES)
    memFree_null((*pastix_data)->l2g_int   );

  if ((*pastix_data)->malrhsd_int == API_YES)
    {
      if (iparm[IPARM_RHSD_CHECK] == API_YES)
        {
          memFree_null((*pastix_data)->b_int);
        }

      (*pastix_data)->malrhsd_int = API_NO;
    }
#endif

  if ((*pastix_data)->malcof)
    {
      if (iparm[IPARM_SCHUR] == API_YES && (*pastix_data)->schur_tab_set == API_YES)
        {
          SolverMatrix * datacode = &((*pastix_data)->solvmatr);
          PASTIX_INT            cblk;

          if (SOLV_TASKNBR > 0)
            {
              cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
              if (SYMB_LCOLNUM(cblk) == (*pastix_data)->n2*(*pastix_data)->iparm[IPARM_DOF_NBR]-1)
                {
                  SOLV_COEFTAB(cblk) = NULL;
                }
            }
        }
    }

#ifdef OOC
  {

    char str[STR_SIZE];
    struct stat stFileInfo;

#  ifndef OOC_DIR
#    define OOC_DIR                "/tmp/pastix"
#  endif

    for (i = 0; i <solvmatr->cblknbr; i++)
      {
        sprintf(str,"%s/pastix_coef_%d/%d",OOC_DIR,(int)iparm[IPARM_OOC_ID],(int)i);
        if (-1 != stat(str,&stFileInfo) && (-1 == remove(str)))
          {
            perror("remove");
            EXIT(MOD_SOPALIN,UNKNOWN_ERR);
          }
      }
    sprintf(str,"%s/pastix_coef_%d",OOC_DIR,(int)iparm[IPARM_OOC_ID]);
    if (-1 == remove(str))
      {
        perror("remove");
        EXIT(MOD_SOPALIN,UNKNOWN_ERR);
      }

    if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
      {
        for (i = 0; i <solvmatr->cblknbr; i++)
          {
            sprintf(str,"%s/pastix_ucoef_%d/%d",OOC_DIR, (int)iparm[IPARM_OOC_ID], (int)i);
            if (-1 != stat(str,&stFileInfo) && (-1 == remove(str)))
              {
                perror("remove");
                EXIT(MOD_SOPALIN,UNKNOWN_ERR);
              }
          }
        sprintf(str,"%s/pastix_ucoef_%d",OOC_DIR, (int)iparm[IPARM_OOC_ID]);
        if (-1 == remove(str))
          {
            perror("remove");
            EXIT(MOD_SOPALIN,UNKNOWN_ERR);
          }
      }

#  ifdef OOC_FTGT
    for (i = 0; i <solvmatr->ftgtnbr; i++)
      {
        sprintf(str,"%s/pastix_ftgt_%d/%d",OOC_DIR,(int)iparm[IPARM_OOC_ID],(int)i);
        if (-1 != stat(str,&stFileInfo) && -1 == remove(str))
          {
            perror("remove");
            EXIT(MOD_SOPALIN,UNKNOWN_ERR);
          }
      }

    sprintf(str,"%s/pastix_ftgt_%d",OOC_DIR,(int)iparm[IPARM_OOC_ID]);
    if (-1 == remove(str))
      {
        perror("remove");
        EXIT(MOD_SOPALIN,UNKNOWN_ERR);
      }
#  endif /* OOC_FTGT */
  }
#endif /* OOC */

#ifdef PASTIX_DEBUG
  if (iparm && dparm) {
    char filename[256];
    sprintf(filename, "parm%ld.dump", (long) procnum);
    PASTIX_FOPEN(stream, filename, "w");
    api_dumparm(stream,iparm,dparm);
    fclose(stream);
  }
#endif

  if ((*pastix_data)->malord)
    {
      orderExit(ordemesh);
      (*pastix_data)->malord=0;
    }

  if ((*pastix_data)->malcsc)
    {
      CscExit(&((*pastix_data)->cscmtx));
      (*pastix_data)->malcsc=0;
    }

  if ((*pastix_data)->malsmx)
    {
      memFree_null(solvmatr->updovct.sm2xtab);
      (*pastix_data)->malsmx=0;
    }

  /* Pour l'instant uniquement si on est en 1d */
  if (iparm && iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
    {
      if (solvmatr->updovct.cblktab)
        for (i=0; i<solvmatr->cblknbr; i++)
          {
            if (solvmatr->updovct.cblktab[i].browcblktab)
              memFree_null(solvmatr->updovct.cblktab[i].browcblktab);

            if (solvmatr->updovct.cblktab[i].browproctab)
              memFree_null(solvmatr->updovct.cblktab[i].browproctab);
          }

      memFree_null(solvmatr->updovct.lblk2gcblk);
      memFree_null(solvmatr->updovct.listblok);
      memFree_null(solvmatr->updovct.listcblk);
      memFree_null(solvmatr->updovct.gcblk2list);
      memFree_null(solvmatr->updovct.loc2glob);
      memFree_null(solvmatr->updovct.cblktab);
      memFree_null(solvmatr->updovct.listptr);
    }

  if (NULL != sopar->b)
    memFree_null(sopar->b);

  if ((*pastix_data)->malslv)
    {
      solverExit(solvmatr);
#ifdef PASTIX_DYNSCHED
      Bubble_Free(solvmatr->btree);
      memFree_null(solvmatr->btree);
#endif
      (*pastix_data)->malslv=0;
    }

  if ((*pastix_data)->sopar.bindtab != NULL)
    memFree_null((*pastix_data)->sopar.bindtab);

  if ((*pastix_data)->listschur != NULL)
    memFree_null((*pastix_data)->listschur);
#ifdef DISTRIBUTED
  if ((*pastix_data)->glob2loc != NULL)
    memFree_null((*pastix_data)->glob2loc);
#endif
#ifdef WITH_SEM_BARRIER
  if ((*pastix_data)->intra_node_procnbr > 1)
    {
      if (sem_close((*pastix_data)->sem_barrier) < 0)
        {
          perror("sem_close");
        }
      if ((*pastix_data)->intra_node_procnum == 0)
        {
          char sem_name[256];
          sprintf(sem_name, "/pastix_%d", (*pastix_data)->pastix_id);
          if (sem_unlink(sem_name) < 0)
            {
              perror("sem_unlink");
            }
        }
    }
#endif

  if (*pastix_data != NULL)
    memFree_null(*pastix_data);

  FreeMpiType();
  FreeMpiSum();
  pastix_print_memory_usage(iparm,pastix_comm);

  printf("CLEAN FIN\n");
}

void pastix_unscale(pastix_data_t *pastix_data, PASTIX_INT sym) {
#ifndef FORCE_MPI
  int size;
  MPI_Comm_size(pastix_data->pastix_comm, &size);
  if(size > 1) {
    errorPrint("pastix_task_unscale is not implemented in the distributed case (scaletabs must be distributed).");
    exit(1);
  }
#endif
  if(pastix_data->scaling == API_YES) {
    if(sym == API_YES)
      Matrix_Unscale_Sym(pastix_data, &pastix_data->solvmatr, pastix_data->scalerowtab, pastix_data->iscalerowtab);
    else
      Matrix_Unscale_Unsym(pastix_data, &pastix_data->solvmatr, pastix_data->scalerowtab, pastix_data->iscalerowtab, pastix_data->scalecoltab, pastix_data->iscalecoltab);
  }
}

#ifdef WITH_SEM_BARRIER
#  define SEM_BARRIER do {                                              \
    if ((*pastix_data)->intra_node_procnum == 0)                        \
      {                                                                 \
        int si_iter;                                                    \
        for (si_iter = 0;                                               \
             si_iter < (*pastix_data)->intra_node_procnbr-1;            \
             si_iter++)                                                 \
          {                                                             \
            sem_post((*pastix_data)->sem_barrier);                      \
          }                                                             \
      }                                                                 \
    else                                                                \
      {                                                                 \
        sem_wait((*pastix_data)->sem_barrier);                          \
      }                                                                 \
  } while(0)
#else
#  define SEM_BARRIER do {} while (0)
#endif
#define SYNC_IPARM do {                                                 \
  if ((*pastix_data)->intra_node_procnbr > 1)                           \
    {                                                                   \
      SEM_BARRIER;                                                      \
      MPI_Bcast(iparm,                                                  \
                IPARM_SIZE, COMM_INT,                                   \
                0, (*pastix_data)->intra_node_comm);                    \
    }                                                                   \
  } while (0)

#define WAIT_AND_RETURN do {                                    \
    if ( *pastix_data != NULL )                                 \
      {                                                         \
        SYNC_IPARM;                                             \
        if (iparm[IPARM_START_TASK] > API_TASK_SOLVE)           \
          {                                                     \
            MPI_Bcast(b, n, COMM_FLOAT, 0,                      \
                      (*pastix_data)->intra_node_comm );        \
          }                                                     \
      }                                                         \
    else {                                                      \
      MPI_Barrier((*pastix_data)->intra_node_comm);             \
    }                                                           \
    return;                                                     \
  } while (0)

/** ****************************************************************************************
 *
 *  Function: pastix
 *
 *  Computes one to all steps of the resolution of Ax=b linear system, using direct methods.
 *
 *  The matrix is given in CSC format.
 *
 *  Parameters:
 *  pastix_data - Data used for a step by step execution.
 *  pastix_comm - MPI communicator which compute the resolution.
 *  n           - Size of the system.
 *  colptr      - Tabular containing the start of each column in row and avals tabulars.
 *  row         - Tabular containing the row number for each element sorted by column.
 *  avals       - Tabular containing the values of each elements sorted by column.
 *  perm        - Permutation tabular for the renumerotation of the unknowns.
 *  invp        - Reverse permutation tabular for the renumerotation of the unknowns.
 *  b           - Right hand side vector(s).
 *  rhs         - Number of right hand side vector(s).
 *  iparm       - Integer parameters given to pastix.
 *  dparm       - Double parameters given to p�stix.
 *
 *  About: Example
 *
 *  from file <simple.c> :
 *
 *  > /\*******************************************\/
 *  > /\*    Check Matrix format                  *\/
 *  > /\*******************************************\/
 *  > /\*
 *  >  * Matrix needs :
 *  >  *    - to be in fortran numbering
 *  >  *    - to have only the lower triangular part in symmetric case
 *  >  *    - to have a graph with a symmetric structure in unsymmetric case
 *  >  *\/
 *  > mat_type = API_SYM_NO;
 *  > if (MTX_ISSYM(type)) mat_type = API_SYM_YES;
 *  > if (MTX_ISHER(type)) mat_type = API_SYM_HER;
 *  > pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
 *  >                    mat_type,  API_YES,
 *  >                    ncol, &colptr, &rows, &values, NULL);
 *  >
 *  > /\*******************************************\/
 *  > /\* Initialize parameters to default values *\/
 *  > /\*******************************************\/
 *  > iparm[IPARM_MODIFY_PARAMETER] = API_NO;
 *  > pastix(&pastix_data, MPI_COMM_WORLD,
 *  >        ncol, colptr, rows, values,
 *  >        perm, invp, rhs, 1, iparm, dparm);
 *  >
 *  > /\*******************************************\/
 *  > /\*       Customize some parameters         *\/
 *  > /\*******************************************\/
 *  > iparm[IPARM_THREAD_NBR] = nbthread;
 *  > iparm[IPARM_SYM] = mat_type;
 *  > switch (mat_type)
 *  >   {
 *  >     case API_SYM_YES:
 *  >       iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
 *  >       break;
 *  >     case API_SYM_HER:
 *  >       iparm[IPARM_FACTORIZATION] = API_FACT_LDLH;
 *  >       break;
 *  >     default:
 *  >       iparm[IPARM_FACTORIZATION] = API_FACT_LU;
 *  >    }
 *  > iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
 *  > iparm[IPARM_END_TASK]            = API_TASK_CLEAN;
 *  >
 *  > /\*******************************************\/
 *  > /\*           Save the rhs                  *\/
 *  > /\*    (it will be replaced by solution)    *\/
 *  > /\*******************************************\/
 *  > rhssaved = malloc(ncol*sizeof(pastix_float_t));
 *  > memcpy(rhssaved, rhs, ncol*sizeof(pastix_float_t));
 *  >
 *  > /\*******************************************\/
 *  > /\*           Call pastix                   *\/
 *  > /\*******************************************\/
 *  > perm = malloc(ncol*sizeof(pastix_int_t));
 *  > invp = malloc(ncol*sizeof(pastix_int_t));
 *  >
 *  > pastix(&pastix_data, MPI_COMM_WORLD,
 *  >  ncol, colptr, rows, values,
 *  >  perm, invp, rhs, 1, iparm, dparm);
 */
void pastix(pastix_data_t **pastix_data,
            MPI_Comm        pastix_comm,
            PASTIX_INT             n,
            PASTIX_INT            *colptr,
            PASTIX_INT            *row,
            PASTIX_FLOAT          *avals,
            PASTIX_INT            *perm,
            PASTIX_INT            *invp,
            PASTIX_FLOAT          *b,
            PASTIX_INT             rhs,
            PASTIX_INT            *iparm,
            double         *dparm)
{
  int flagWinvp = 1;
  int ret = NO_ERR;
#ifdef FIX_SCOTCH /* Pour le debug au cines */
  _SCOTCHintRandInit();
#endif

  iparm[IPARM_GRAPHDIST] = API_NO;

  if (iparm[IPARM_MODIFY_PARAMETER] == API_NO) /* init task */
    {
      /* init with default for iparm & dparm */
      pastix_initParam(iparm, dparm);
      iparm[IPARM_GRAPHDIST] = API_NO;
      return;
    }
  /*
   * Init : create pastix_data structure if it's first time
   */
  if (*pastix_data == NULL)
    {
      /* Need to be set to -1 in every cases */
      iparm[IPARM_OOC_ID] = -1;

      /* Allocation de la structure pastix_data qd on rentre dans
         pastix pour la premi�re fois */
      pastix_task_init(pastix_data, pastix_comm, iparm, dparm);
      if ((*pastix_data)->intra_node_procnum == 0) {
        /* Affichage des options */
        pastix_welcome_print(*pastix_data, colptr, n);

        /* Matrix verification */
        if (iparm[IPARM_MATRIX_VERIFICATION] == API_YES)
          if ( NO_ERR != (ret = pastix_checkMatrix((*pastix_data)->inter_node_comm,
                                                   iparm[IPARM_VERBOSE], iparm[IPARM_SYM],
                                                   API_NO, n, &colptr, &row,
                                                   (avals == NULL)?NULL:(&avals),
                                                   NULL, iparm[IPARM_DOF_NBR])))
            {
              errorPrint("The matrix is not in the correct format");
              iparm[IPARM_ERROR_NUMBER] = ret;
              return;
            }
      }
    }

  (*pastix_data)->n  = n;

  if (NO_ERR != (ret = pastix_check_param(*pastix_data, rhs)))
    {
      iparm[IPARM_ERROR_NUMBER] = ret;
      return;
    }


  if ((*pastix_data)->intra_node_procnum == 0) {
    /* only master node do computations */
    if (iparm[IPARM_END_TASK]<API_TASK_ORDERING) {
      WAIT_AND_RETURN;
    }

    /* On n'affiche pas le warning si on enchaine scotch et fax */
    if ( (iparm[IPARM_START_TASK] <= API_TASK_ORDERING)
         && (iparm[IPARM_END_TASK] >= API_TASK_SYMBFACT)) {
      flagWinvp = 0;
    }

#ifdef SEPARATE_ZEROS
    {
      PASTIX_INT  n_nz       = 0;
      PASTIX_INT *colptr_nz  = NULL;
      PASTIX_INT *rows_nz    = NULL;
      PASTIX_INT  n_z        = 0;
      PASTIX_INT *colptr_z   = NULL;
      PASTIX_INT *rows_z     = NULL;
      PASTIX_INT *permz      = NULL;
      PASTIX_INT *revpermz      = NULL;
      int  ret;
      Order *ordemesh = &((*pastix_data)->ordemesh);
      Order *order    = NULL;
      PASTIX_INT  iter;

      MALLOC_INTERN(order, 1, Order);
      orderInit(order);
      MALLOC_INTERN(order->rangtab, n+1, PASTIX_INT);
      MALLOC_INTERN(order->permtab, n,   PASTIX_INT);

      MALLOC_INTERN(permz, n, PASTIX_INT);
      MALLOC_INTERN(revpermz, n, PASTIX_INT);

      CSC_buildZerosAndNonZerosGraphs(n,
                                      colptr,
                                      row,
                                      avals,
                                      &n_nz,
                                      &colptr_nz,
                                      &rows_nz,
                                      &n_z,
                                      &colptr_z,
                                      &rows_z,
                                      permz,
                                      revpermz,
                                      dparm[DPARM_EPSILON_MAGN_CTRL]);

      if (n_z != 0 && n_nz != 0)
        {
          fprintf(stdout, ":: Scotch on non zeros\n");
          /*
           * Scotch : Ordering
           */
          if (iparm[IPARM_START_TASK] == API_TASK_ORDERING) /* scotch task */
            if (NO_ERR != (ret = pastix_task_scotch(pastix_data,
                                                    (*pastix_data)->inter_node_comm,
                                                    n_nz, colptr_nz, rows_nz,
                                                    perm, invp)))
              {
                iparm[IPARM_ERROR_NUMBER] = ret;
                WAIT_AND_RETURN;
              }
          for (iter = 0; iter < n_nz; iter++)
            order->permtab[revpermz[iter]-1] = ordemesh->permtab[iter];
          memcpy(order->rangtab,
                 ordemesh->rangtab,
                 (ordemesh->cblknbr +1)*sizeof(PASTIX_INT));
          order->cblknbr = ordemesh->cblknbr;
          iparm[IPARM_START_TASK]--;
          fprintf(stdout, ":: Scotch on zeros\n");
          if (iparm[IPARM_START_TASK] == API_TASK_ORDERING) /* scotch task */
            if (NO_ERR != (ret = pastix_task_scotch(pastix_data,
                                                    (*pastix_data)->inter_node_comm,
                                                    n_z, colptr_z, rows_z,
                                                    perm, invp)))
              {
                iparm[IPARM_ERROR_NUMBER] = ret;
                WAIT_AND_RETURN;
              }

          for (iter = 0; iter < n_z; iter++)
            order->permtab[revpermz[n_nz+iter]-1] = n_nz+ordemesh->permtab[iter];

          for (iter = 0; iter < ordemesh->cblknbr+1; iter++)
            order->rangtab[order->cblknbr+iter] = order->rangtab[order->cblknbr] +
              ordemesh->rangtab[iter];

          iparm[IPARM_START_TASK]--;

          if (iparm[IPARM_START_TASK] == API_TASK_ORDERING) /* scotch task */
            if (NO_ERR != (ret = pastix_task_scotch(pastix_data,
                                                    (*pastix_data)->inter_node_comm,
                                                    n, colptr, row, perm, invp)))
              {
                iparm[IPARM_ERROR_NUMBER] = ret;
                WAIT_AND_RETURN;
              }
          MALLOC_INTERN(order->peritab, n,   PASTIX_INT);
          for (iter = 0; iter < n; iter++)
            {
              ASSERT(order->permtab[iter] < n, MOD_SOPALIN);
              ASSERT(order->permtab[iter]  >= 0, MOD_SOPALIN);
            }
          for (iter = 0; iter < n; iter++)
            order->peritab[order->permtab[iter]]=iter;
          for (iter = 0; iter < n; iter++)
            {
              ASSERT(order->peritab[iter] < n, MOD_SOPALIN);
              ASSERT(order->peritab[iter]  >= 0, MOD_SOPALIN);
            }
          memFree_null(ordemesh->rangtab);
          memFree_null(ordemesh->permtab);
          memFree_null(ordemesh->peritab);
          memcpy(ordemesh, order, sizeof(Order));
          memFree_null(order);
          memFree_null(colptr_nz);
          memFree_null(colptr_z);
          memFree_null(rows_nz);
          memFree_null(rows_z);
        }
    }
#endif /* SEPARATE_ZEROS */

    if (iparm[IPARM_ISOLATE_ZEROS] == API_YES &&
        iparm[IPARM_SCHUR] == API_YES)
      {
        errorPrint("Schur complement is incompatible with diagonal zeros isolation.");
        iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
        WAIT_AND_RETURN;
      }

    /*
     * Scotch : Ordering
     */
    if (iparm[IPARM_START_TASK] == API_TASK_ORDERING) /* scotch task */
      {
        if (iparm[IPARM_ISOLATE_ZEROS] == API_YES)
          {
            PASTIX_INT itercol;
            PASTIX_INT iterrow;
            PASTIX_INT iterschur = 0;
            int found;

            (*pastix_data)->nschur=0;
            for (itercol = 0; itercol < n; itercol++)
              {
                found = API_NO;
                for (iterrow = colptr[itercol]-1; iterrow <  colptr[itercol+1]-1; iterrow++)
                  {
                    if (row[iterrow]-1 == itercol)
                      {
                        if (ABS_FLOAT(avals[iterrow]) < dparm[DPARM_EPSILON_REFINEMENT])
                          {
                            (*pastix_data)->nschur++;
                          }
                        found = API_YES;
                        break;
                      }
                  }
                if (found == API_NO)
                  {
                    (*pastix_data)->nschur++;
                  }
              }

            if ((*pastix_data)->nschur > 0) {
            MALLOC_INTERN((*pastix_data)->listschur,
                          (*pastix_data)->nschur,
                          PASTIX_INT);
	      for (itercol = 0; itercol < n; itercol++) {
                found = API_NO;
                for (iterrow = colptr[itercol]-1; iterrow <  colptr[itercol+1]-1; iterrow++) {
		  if (row[iterrow]-1 == itercol) {
		    if (ABS_FLOAT(avals[iterrow]) < dparm[DPARM_EPSILON_REFINEMENT]) {
                            (*pastix_data)->listschur[iterschur] = itercol+1;
                            iterschur++;
                          }
                        found = API_YES;
                        break;
                      }
                  }
                if (found == API_NO) {
                    (*pastix_data)->listschur[iterschur] = itercol+1;
                    iterschur++;
                  }
              }
	    } else {
	      iparm[IPARM_ISOLATE_ZEROS] = API_NO;
	    }
          }
        if (NO_ERR != (ret = pastix_task_scotch(pastix_data,
                                                (*pastix_data)->inter_node_comm,
                                                n, colptr, row, perm, invp)))
          {
            iparm[IPARM_ERROR_NUMBER] = ret;
            WAIT_AND_RETURN;
          }
        if (iparm[IPARM_ISOLATE_ZEROS] == API_YES)
          {
            memFree_null((*pastix_data)->listschur);
          }
      }
    if (iparm[IPARM_END_TASK]<API_TASK_SYMBFACT) {
      WAIT_AND_RETURN;
    }

    /*
     * Fax : Facto symbolic
     */
    if (iparm[IPARM_START_TASK] == API_TASK_SYMBFACT) /* Fax task */
      pastix_task_fax(pastix_data, pastix_comm, perm, invp, flagWinvp);

    if (iparm[IPARM_END_TASK]<API_TASK_ANALYSE) {
      WAIT_AND_RETURN;
    }

    /*
     * Blend : Scheduling
     */
    if (iparm[IPARM_START_TASK] == API_TASK_ANALYSE) /* Blend task */
      pastix_task_blend(pastix_data, (*pastix_data)->inter_node_comm);

    if (iparm[IPARM_END_TASK]<API_TASK_NUMFACT) {
      WAIT_AND_RETURN;
    }

#if defined(PROFILE) && defined(MARCEL)
    profile_activate(FUT_ENABLE, MARCEL_PROF_MASK, 0);
    marcel_printf("DEBUT profil marcel\n");
#endif

    /*
     * Sopalin : Factorisation
     */
    if (iparm[IPARM_START_TASK] == API_TASK_NUMFACT) /* Sopalin task */
      {
        ret = pastix_task_sopalin(*pastix_data,
                                  (*pastix_data)->inter_node_comm, n,
                                  colptr, row, avals, b, rhs, NULL);


        MPI_Bcast(&ret, 1, MPI_INT, 0, (*pastix_data)->inter_node_comm);
        if (NO_ERR != ret) {
          iparm[IPARM_ERROR_NUMBER] = ret;
          WAIT_AND_RETURN;
        }
      }
    if (iparm[IPARM_END_TASK]<iparm[IPARM_START_TASK]) {
      WAIT_AND_RETURN;
    }

    /*
     * Updo : solve
     */
    if (iparm[IPARM_START_TASK] == API_TASK_SOLVE) /* Updown task */
      {
        /* For thread comm */
        (*pastix_data)->sopar.stopthrd = API_YES;
        pastix_task_updown(*pastix_data, (*pastix_data)->inter_node_comm,
                           n, b, rhs, NULL);
        /* For thread comm */
        (*pastix_data)->sopar.stopthrd = API_NO;
      }
    if (iparm[IPARM_END_TASK]<API_TASK_REFINE) {
      WAIT_AND_RETURN;
    }

    /*
     * Raff
     */
    if (iparm[IPARM_START_TASK] == API_TASK_REFINE) /* Refinement task */
      {
        pastix_task_raff(*pastix_data, (*pastix_data)->inter_node_comm,
                         n, b, rhs, NULL);
      }
    if (iparm[IPARM_END_TASK]<API_TASK_CLEAN) {
      WAIT_AND_RETURN;
    }

    /*
     * Clean
     */
  } /* (*pastix_data)->intra_node_procnum == 0 */

  SYNC_IPARM;
  if (iparm[IPARM_END_TASK]<API_TASK_CLEAN)
    return;

  if (iparm[IPARM_START_TASK] == API_TASK_CLEAN)
    pastix_task_clean(pastix_data, pastix_comm);

#if defined(PROFILE) && defined(MARCEL)
  profile_stop();
  marcel_printf("FIN profil marcel\n");
#endif
}

#define REDISTRIBUTE_RHS                                        \
  {                                                             \
                                                                \
    if (b != NULL && rhsHasBeenRedistributed == API_NO)         \
      {                                                         \
        rhsHasBeenRedistributed = API_YES;                      \
        if (rhs_need_redispatch == API_YES)                     \
          {                                                     \
            if ((*pastix_data)->procnum == 0 &&                 \
                iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)        \
              fprintf(stdout,OUT_REDIS_RHS);                    \
                                                                \
            /* Distribute the user RHS into                     \
             intern distribution */                             \
            if (b_int == NULL)                                  \
              {                                                 \
                MALLOC_INTERN(b_int, ncol_int*rhs, PASTIX_FLOAT);      \
                (*pastix_data)->b_int = b_int;                  \
              }                                                 \
            redispatch_rhs(n,                                   \
                           b,                                   \
                           rhs,                                 \
                           loc2glob,                            \
                           ncol_int,                            \
                           b_int,                               \
                           l2g_int,                             \
                           (*pastix_data)->procnbr,             \
                           (*pastix_data)->procnum,             \
                           pastix_comm,                         \
                           iparm[IPARM_DOF_NBR]);               \
                                                                \
          }                                                     \
        else                                                    \
          {                                                     \
            b_int = b;                                          \
          }                                                     \
      }                                                         \
  }
#define REDISTRIBUTE_SOL                                \
  {                                                     \
   if (rhs_need_redispatch == API_YES)                  \
     {                                                  \
      if ((*pastix_data)->procnum == 0 &&               \
            iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)    \
        fprintf(stdout,OUT_REDIS_SOL);                  \
                                                        \
      redispatch_rhs(ncol_int,                          \
                     b_int,                             \
                     rhs,                               \
                     l2g_int,                           \
                     n,                                 \
                     b,                                 \
                     loc2glob,                          \
                     (*pastix_data)->procnbr,           \
                     (*pastix_data)->procnum,           \
                     pastix_comm,                       \
                     iparm[IPARM_DOF_NBR]);             \
     }                                                  \
  }

/** Function: dpastix

    Computes one to all steps of the resolution of
    Ax=b linear system, using direct methods.
    Here the matrix is given distributed.

    The matrix is given in CSCD format.

    Parameters:
    pastix_data - Data used for a step by step execution.
    pastix_comm - MPI communicator which compute the resolution.
    n           - Size of the system.
    colptr      - Tabular containing the start of each column in row and avals tabulars.
    row         - Tabular containing the row number for each element sorted by column.
    avals       - Tabular containing the values of each elements sorted by column.
    loc2glob    - Global column number of the local columns.
    perm        - Permutation tabular for the renumerotation of the unknowns.
    invp        - Reverse permutation tabular for the renumerotation of the unknowns.
    b           - Right hand side vector(s).
    rhs         - Number of right hand side vector(s).
    iparm       - Integer parameters given to pastix.
    dparm       - Double parameters given to p�stix.
*/
#ifndef dpastix
#  error "REDEFINIR dpastix dans redefine_functions.h"
#endif
void dpastix(pastix_data_t **pastix_data,
             MPI_Comm        pastix_comm,
             PASTIX_INT             n,
             PASTIX_INT            *colptr,
             PASTIX_INT            *row,
             PASTIX_FLOAT          *avals,
             PASTIX_INT            *loc2glob,
             PASTIX_INT            *perm,
             PASTIX_INT            *invp,
             PASTIX_FLOAT          *b,
             PASTIX_INT             rhs,
             PASTIX_INT            *iparm,
             double         *dparm)
{
#ifdef DISTRIBUTED
  int    flagWinvp               = 1;
  PASTIX_INT    ncol_int                = 0;
  PASTIX_INT   *l2g_int                 = NULL;
  PASTIX_FLOAT *b_int                   = NULL;
  int    ret                     = NO_ERR;
  int    ret_rcv                 = NO_ERR;
  PASTIX_INT    gN                      = -1;
  int    rhs_need_redispatch     = API_NO;
  int    rhsHasBeenRedistributed = API_NO;
  int    mayNeedReturnSol        = API_NO;
#  ifdef FIX_SCOTCH /* Pour le debug au cines */
  _SCOTCHintRandInit();
#  endif


  if (iparm[IPARM_MODIFY_PARAMETER] == API_NO) /* init task */
    {
      /* init with default for iparm & dparm */
      pastix_initParam(iparm, dparm);
      return;
    }

  /* Si pastix_data est nul, c'est qu'on rentre dans
     la fonction pour la premi�re fois */
  if (*pastix_data == NULL)
    {
      iparm[IPARM_OOC_ID]         = -1;
      /* initialisation et allocation de pastix_data */
      pastix_task_init(pastix_data,pastix_comm,iparm,dparm);

      /* Affichage des options */
      pastix_welcome_print(*pastix_data, colptr, n);

      /* multiple RHS see MRHS_ALLOC */
      if ( ((*pastix_data)->procnum == 0)  && (rhs!=1) )
        errorPrintW("multiple right-hand-side not tested...");

      /* Matrix verification */
      if (iparm[IPARM_MATRIX_VERIFICATION] == API_YES)
        if (NO_ERR != (ret = pastix_checkMatrix(pastix_comm, iparm[IPARM_VERBOSE],
                                                iparm[IPARM_SYM], API_NO,
                                                n, &colptr, &row , (avals == NULL)?NULL:(&avals),
                                                ( (iparm[IPARM_GRAPHDIST] == API_NO)?
                                                  NULL:(&loc2glob) ),
                                                iparm[IPARM_DOF_NBR])))
          {
            errorPrint("The matrix is not in the correct format");
            iparm[IPARM_ERROR_NUMBER] = ret;
            return;
          }

    }
  (*pastix_data)->n  = n;


  if (NO_ERR != (ret = pastix_check_param(*pastix_data, rhs)))
    {
      iparm[IPARM_ERROR_NUMBER] = ret;
      return;
    }


  /* On n'affiche pas le warning si on enchaine scotch et fax */
  if ( (iparm[IPARM_START_TASK] <= API_TASK_ORDERING)
       && (iparm[IPARM_END_TASK] >= API_TASK_SYMBFACT))
    flagWinvp = 0;


  /*
   * Scotch : Ordering
   */
  if (iparm[IPARM_START_TASK] == API_TASK_ORDERING) /* scotch task */
    {

      if (iparm[IPARM_GRAPHDIST] == API_YES)
        {
          if (NO_ERR != (ret = dpastix_task_scotch(pastix_data, pastix_comm,
                                                   n, colptr, row,
                                                   perm, invp, loc2glob)))
            {
              errorPrint("Error in ordering task\n");
              iparm[IPARM_ERROR_NUMBER] = ret;
              return;
            }
        }
      else
        {
          if ((*pastix_data)->intra_node_procnum == 0) {
            if (NO_ERR != (ret = pastix_task_scotch(pastix_data,
                                                    (*pastix_data)->inter_node_comm,
                                                    n, colptr, row,
                                                    perm, invp)))
              {
                errorPrint("Error in ordering task\n");
                iparm[IPARM_ERROR_NUMBER] = ret;
              }
          }
          SYNC_IPARM;
          if (iparm[IPARM_ERROR_NUMBER] != NO_ERR)
            return;
        }

    }
  if (iparm[IPARM_END_TASK]<API_TASK_SYMBFACT)
    return;

  /*
   * Fax : Facto symbolic
   */
  if (iparm[IPARM_START_TASK] == API_TASK_SYMBFACT) /* Fax task */
    {
      if (iparm[IPARM_GRAPHDIST] == API_YES)
        {
          dpastix_task_fax(*pastix_data,
                           (*pastix_data)->inter_node_comm,
                           n, perm, loc2glob, flagWinvp);
        }
      else
        {
          if ((*pastix_data)->intra_node_procnum == 0)
            {
              pastix_task_fax(pastix_data,
                              (*pastix_data)->inter_node_comm,
                              perm, invp, flagWinvp);
            }
        }
      SYNC_IPARM;
    }

  if (iparm[IPARM_END_TASK]<API_TASK_ANALYSE)
    return;

  if ((*pastix_data)->intra_node_procnum == 0)
    {
      if (iparm[IPARM_START_TASK] == API_TASK_ANALYSE) /* Blend task */
        {
          pastix_task_blend(pastix_data, (*pastix_data)->inter_node_comm);
        }
    }
  SYNC_IPARM;
  if (iparm[IPARM_END_TASK]<API_TASK_NUMFACT)
    return;

#  if defined(PROFILE) && defined(MARCEL)
  profile_activate(FUT_ENABLE, MARCEL_PROF_MASK, 0);
  marcel_printf("DEBUT profil marcel\n");
#  endif

  if (iparm[IPARM_START_TASK] == API_TASK_NUMFACT) /* Sopalin task */
    {
      mayNeedReturnSol = API_YES;
      if ((*pastix_data)->intra_node_procnum == 0)
        {
          ret = pastix_task_sopalin(*pastix_data,
                                    (*pastix_data)->inter_node_comm,
                                    n,
                                    colptr,
                                    row,
                                    avals,
                                    b,
                                    rhs,
                                    loc2glob);

          MPI_Allreduce(&ret, &ret_rcv, 1, MPI_INT, MPI_MAX,
                        (*pastix_data)->inter_node_comm);
          if (NO_ERR != ret_rcv)
            {
              errorPrint("Error in numeric factorisation task\n");
              iparm[IPARM_ERROR_NUMBER] = ret_rcv;
            }
        }
      else
        {
          /* Remplissage de la csc interne */
          if ((*pastix_data)->cscInternFilled == API_NO) {
            pastix_fake_fillin_csc(*pastix_data, pastix_comm, n,
                                   colptr, row, avals, b, rhs, loc2glob);
          }
        }
      SYNC_IPARM;

      if (iparm[IPARM_ERROR_NUMBER] != NO_ERR)
        return;

      if (iparm[IPARM_START_TASK] > API_TASK_SOLVE)
        {
          ncol_int = (*pastix_data)->ncol_int;
          b_int    = (*pastix_data)->b_int;
          l2g_int  = (*pastix_data)->l2g_int;
          if ((*pastix_data)->malrhsd_int)
            rhs_need_redispatch = API_YES;
          REDISTRIBUTE_SOL;
        }

      if (iparm[IPARM_END_TASK] < API_TASK_CLEAN)
        return;

    }



  if ((*pastix_data)->intra_node_procnum == 0)
    {
      ncol_int = (*pastix_data)->ncol_int;
      b_int    = (*pastix_data)->b_int;
      l2g_int  = (*pastix_data)->l2g_int;

      /* User can change CSCD after blend */
      if ((iparm[IPARM_GRAPHDIST] == API_YES) &&
          ((*pastix_data)->glob2loc == NULL))
        {
          cscd_build_g2l(ncol_int,
                         l2g_int,
                         (*pastix_data)->inter_node_comm,
                         &gN,
                         &((*pastix_data)->glob2loc));
        }


      /* Updown task */

      /* If user has not specified that he is
         absolutly certain that is CSCd is
         correctly distributed */
      if ( ( iparm[IPARM_START_TASK] == API_TASK_SOLVE ||
             iparm[IPARM_START_TASK] == API_TASK_REFINE ) &&
           iparm[IPARM_CSCD_CORRECT] == API_NO )
        {
          PASTIX_INT my_n;
          PASTIX_INT * my_l2g = NULL;
          int OK, OK_RECV;
          PASTIX_INT iter;
          /* Test que la cscd utilisateur correspond a la cscd pastix */
          my_n = pastix_getLocalNodeNbr(pastix_data);

          OK = 0;
          if (my_n != n)
            {
              OK = 1;
            }
          else
            {
              if ((*pastix_data)->l2g_int) {
                my_l2g = (*pastix_data)->l2g_int;
              } else {
                MALLOC_INTERN(my_l2g, my_n, PASTIX_INT);
                pastix_getLocalNodeLst(pastix_data, my_l2g);
              }
              for (iter = 0; iter < my_n; iter++)
                {
                  if (my_l2g[iter] != loc2glob[iter])
                    {
                      OK = 1;
                      break;
                    }
                }
              if (!(*pastix_data)->l2g_int) {
                memFree_null(my_l2g);
              }
            }
          MPI_Allreduce(&OK, &OK_RECV, 1, MPI_INT, MPI_SUM, pastix_comm);
          if (OK_RECV != 0)
            rhs_need_redispatch = API_YES;
        }

      if (iparm[IPARM_START_TASK] == API_TASK_SOLVE)
        {
          mayNeedReturnSol = API_YES;
          REDISTRIBUTE_RHS;
          /* For thread comm */
          (*pastix_data)->sopar.stopthrd = API_YES;
          pastix_task_updown(*pastix_data,
                             (*pastix_data)->inter_node_comm,
                             ncol_int, b_int, rhs, l2g_int);
          /* For thread comm */
          (*pastix_data)->sopar.stopthrd = API_NO;
        }
    }
  else
    {
      /* If user has not specified that he is
         absolutly certain that is CSCd is
         correctly distributed */
      if ( ( iparm[IPARM_START_TASK] == API_TASK_SOLVE ||
             iparm[IPARM_START_TASK] == API_TASK_REFINE ) &&
           iparm[IPARM_CSCD_CORRECT] == API_NO )
        {
          int OK, OK_RECV;
          OK = 0;

          MPI_Allreduce(&OK, &OK_RECV, 1, MPI_INT, MPI_SUM, pastix_comm);
          if (OK_RECV != 0)
            rhs_need_redispatch = API_YES;
        }

      if (iparm[IPARM_START_TASK] == API_TASK_SOLVE)
        {
          mayNeedReturnSol = API_YES;
          REDISTRIBUTE_RHS;
        }
    }

  if ( mayNeedReturnSol &&
       iparm[IPARM_END_TASK]<API_TASK_REFINE)
    {
      REDISTRIBUTE_SOL;
      return;
    }

  if (iparm[IPARM_START_TASK] == API_TASK_REFINE) /* Refinement task */
    {
      /* If it wasn't done just after solve */
      REDISTRIBUTE_RHS;
      if ((*pastix_data)->intra_node_procnum == 0)
        {
          pastix_task_raff(*pastix_data, (*pastix_data)->inter_node_comm,
                           ncol_int, b_int, rhs, l2g_int);
        }
      REDISTRIBUTE_SOL;
    }

  SYNC_IPARM;
  if (iparm[IPARM_END_TASK]<API_TASK_CLEAN)
    return;

  if (iparm[IPARM_START_TASK] == API_TASK_CLEAN)
    pastix_task_clean(pastix_data, pastix_comm);

#  if defined(PROFILE) && defined(MARCEL)
  profile_stop();
  marcel_printf("FIN profil marcel\n");
#  endif
#else
  (void)pastix_data; (void)pastix_comm; (void)n; (void)colptr; (void)row;
  (void)avals; (void)loc2glob; (void)perm; (void)invp; (void)b; (void)rhs;
  (void)iparm; (void)dparm;
  errorPrint("To use dpastix please compile with -DDISTRIBUTED");
  iparm[IPARM_ERROR_NUMBER] = BAD_DEFINE_ERR;
#endif /* DISTRIBUTED */
}



/*
  Function: pastix_bindThreads

  Set bindtab in pastix_data, it gives for each thread the CPU to bind in to.
  bindtab follows this organisation :

  bindtab[threadnum] = cpu to set thread threadnum.

  Parameters:
  pastix_data - Structure de donn�e pour l'utilisation step by step
  thrdnbr     - Nombre de threads / Taille du tableau
  bindtab     - Tableau de correspondance entre chaque thread et coeur de la machine
*/

void pastix_bindThreads ( pastix_data_t *pastix_data, PASTIX_INT thrdnbr, PASTIX_INT *bindtab)
{
  int i;

  if ( pastix_data == NULL )
    {
      errorPrint("Pastix_data need to be initialized before to try to set bindtab.");
      EXIT(MOD_SOPALIN, BADPARAMETER_ERR);
    }

  /* Copy association tab between threads and cores */
  MALLOC_INTERN(pastix_data->sopar.bindtab, thrdnbr, int);
  for (i = 0; i < thrdnbr; i++)
    {
      pastix_data->sopar.bindtab[i] = bindtab[i];
    }
  /* Check values in bindtab */
  {
    int nbproc;
#ifdef MARCEL
    nbproc = marcel_nbvps();
#else
    nbproc = sysconf(_SC_NPROCESSORS_ONLN);
#endif

    for (i=0; i< thrdnbr; i++)
      if (!(pastix_data->sopar.bindtab[i] < nbproc))
        {
          errorPrint("Try to bind thread on an unavailable core.");
          EXIT(MOD_SOPALIN, BADPARAMETER_ERR);
        }

  }

  pastix_data->bindtab = pastix_data->sopar.bindtab;
  return;
}
/*
 * Function: pastix_checkMatrix_int
 *
 * Check the matrix :
 * - Renumbers in Fortran numerotation (base 1) if needed (base 0)
 * - Check that the matrix contains no doubles,  with flagcor == API_YES,
 *   correct it.
 * - Can scale the matrix if compiled with -DMC64 -DSCALING (untested)
 * - Checks the symetry of the graph in non symmetric mode.
 *   With non distributed matrices, with flagcor == API_YES,
 *   correct the matrix.
 * - sort the CSC.
 *
 * Parameters:
 *   pastix_comm - PaStiX MPI communicator
 *   verb        - Level of prints (API_VERBOSE_[NOT|NO|YES])
 *   flagsym     - Indicate if the given matrix is symetric
 *                 (API_SYM_YES or API_SYM_NO)
 *   flagcor     - Indicate if we permit the function to reallocate the matrix.
 *   n           - Number of local columns.
 *   colptr      - First element of each row in *row* and *avals*.
 *   row         - Row of each element of the matrix.
 *   avals       - Value of each element of the matrix.
 *   loc2glob    - Global column number of local columns
 *                 (NULL if not distributed).
 *   dof         - Number of degrees of freedom.
 *   flagalloc   - indicate if allocation on CSC uses internal malloc.
 */
#pragma intel optimization_level 0
PASTIX_INT pastix_checkMatrix_int(MPI_Comm pastix_comm,
                                  PASTIX_INT      verb,
                                  PASTIX_INT      flagsym,
                                  PASTIX_INT      flagcor,
                                  PASTIX_INT      n,
                                  PASTIX_INT    **colptr,
                                  PASTIX_INT    **row,
                                  PASTIX_FLOAT  **avals,
                                  PASTIX_INT    **loc2glob,
                                  PASTIX_INT      dof,
                                  PASTIX_INT      flagalloc)
{
  int  procnum;
  int  ret;
  int  OK;
  int  OK_RECV;
  PASTIX_INT  old;
  PASTIX_INT  i;
  PASTIX_INT  l2g_sum_n[2];
  PASTIX_INT  l2g_sum_n_reduced[2];

  MPI_Comm_rank(pastix_comm, &procnum);

  if (verb > API_VERBOSE_NOT)
    print_onempi("%s","Check : Numbering");

  if (!(*colptr)[0])
    {
      /* fortran-style numbering */
      if (verb > API_VERBOSE_NOT)
        print_onempi("%s", "\n\tC numbering to Fortran Numbering\tOK\n");
      CSC_Cnum2Fnum(*row,*colptr,n);
      if (loc2glob != NULL)
        for (i = 0; i <  n; i++)
          (*loc2glob)[i]++;
    }
  else
    {
      if (verb > API_VERBOSE_NOT)
        print_onempi("%s", "\t\tOK\n");
    }

  if (loc2glob != NULL)
    {
      PASTIX_INT l2g_OK;
      PASTIX_INT l2g_OK_rcv;
      l2g_sum_n[0] = 0;
      l2g_sum_n[1] = n;
      l2g_OK       = 0;

      for (i = 0; i < n ; i++)
        {
          l2g_sum_n[0] += (*loc2glob)[i];
          if ( i > 0 && (*loc2glob)[i] <= (*loc2glob)[i-1] )
            {
              l2g_OK = 1;
            }
        }
      MPI_Allreduce(&l2g_OK,  &l2g_OK_rcv,        1,
                    COMM_INT, MPI_SUM, pastix_comm);
      if (l2g_OK_rcv > 0)
        {
          print_onempi("%s", "Local column must be ordered increasingly\n");
          RETURN_ERROR(BADPARAMETER_ERR);
        }
      MPI_Allreduce(l2g_sum_n, l2g_sum_n_reduced, 2, COMM_INT, MPI_SUM, pastix_comm);
      /* All column have been distributed */
      if (2*l2g_sum_n_reduced[0] != (l2g_sum_n_reduced[1]*(l2g_sum_n_reduced[1]+1)))
        {
          print_onempi("%s", "All column must be destributed once and only once\n");
          RETURN_ERROR(BADPARAMETER_ERR);
        }
    }


  /* sorting */
  if (verb > API_VERBOSE_NOT)
    print_onempi("%s", "Check : Sort CSC");

  if (avals != NULL)
    CSC_sort(n,*colptr,*row,*avals, dof);
  else
    CSC_sort(n,*colptr,*row,NULL, 0);
  if (verb > API_VERBOSE_NOT)
    print_onempi("%s","\t\tOK\n");


  if (verb > API_VERBOSE_NOT)
    print_onempi("%s", "Check : Duplicates");

  old = (*colptr)[n]-1;
  /* Preserve sorting */
  ret = csc_check_doubles(n,
                          *colptr,
                          row,
                          avals,
                          dof,
                          flagcor,
                          flagalloc);

  if (loc2glob != NULL)
    {
      if (ret == API_YES) {OK = 0;}
      MPI_Allreduce(&OK, &OK_RECV, 1, MPI_INT, MPI_SUM, pastix_comm);
      if (OK_RECV>0) {ret = API_NO;}
    }

  if (ret == API_YES)
    {
      if (verb > API_VERBOSE_NOT)
        print_onempi("%s", "\t\tOK\n");
    }
  else
    {
      if (verb > API_VERBOSE_NOT)
        print_onempi("%s", "\t\t\tKO\n");
      RETURN_ERROR(MATRIX_ERR);
    }

  if (verb > API_VERBOSE_NOT)
    {
      if (old != ((*colptr)[n] - 1))
        {
          if (loc2glob != NULL)
            {
              fprintf(stdout,
                      "\n\t%ld double terms merged on proc %ld \n",
                      (long)(old - (*colptr)[n]-1),(long)procnum);
            }
          else
            {
              print_onempi("\t%ld double terms merged\n",
                           (long)(old - ((*colptr)[n]-1)));
            }
        }
    }
  {
    PASTIX_INT cnt_lower     = 0;
    PASTIX_INT cnt_upper     = 0;
    PASTIX_INT cnt_diag      = 0;
    PASTIX_INT cnt_num_zeros = 0;
    PASTIX_INT globn         = n;
    PASTIX_INT itercol;
    PASTIX_INT iterrow;

    for (itercol = 0; itercol < n; itercol ++) {
      for ( iterrow = (*colptr)[itercol] - 1;
            iterrow < (*colptr)[itercol+1] - 1;
            iterrow++) {
        PASTIX_INT column = itercol;
        if (loc2glob != NULL)
          column = (*loc2glob)[itercol]-1;
        if ((*row)[iterrow]-1 > column) {
          cnt_lower++;
        } else {
          if ((*row)[iterrow]-1 < column) {
            cnt_upper++;
          } else {
            cnt_diag++;
            if (avals != NULL && (*avals)[iterrow*dof*dof] == 0.)
              cnt_num_zeros++;
          }
        }
      }
    }


    if (loc2glob != NULL)
      {
        PASTIX_INT send_data[5];
        PASTIX_INT recv_data[5];

        send_data[0] = cnt_lower;
        send_data[1] = cnt_upper;
        send_data[2] = cnt_diag;
        send_data[3] = cnt_num_zeros;
        send_data[4] = globn;
        MPI_Allreduce(send_data, recv_data, 5, COMM_INT, MPI_SUM, pastix_comm);
        cnt_lower      = recv_data[0];
        cnt_upper      = recv_data[1];
        cnt_diag       = recv_data[2];
        cnt_num_zeros  = recv_data[3];
        globn          = recv_data[4];
      }

    if (cnt_diag != globn)
      if (verb > API_VERBOSE_NOT)
        errorPrintW("%d/%d structural zeros found on the diagonal.",
                    globn-cnt_diag, globn);

    if (cnt_num_zeros != 0)
      if (verb > API_VERBOSE_NOT)
        errorPrintW("%d numerical zeros found on the diagonal.", cnt_num_zeros);

    if (cnt_upper == cnt_lower &&
        cnt_lower != 0 &&
        ( flagsym == API_SYM_YES  ||
          flagsym == API_SYM_HER ) &&
        flagcor == API_YES)
      {
        PASTIX_INT index = 0;
        PASTIX_INT lastindex = 0;
        PASTIX_INT   * tmprows;
        PASTIX_FLOAT * tmpvals;
        errorPrintW("Upper and lower part given on a symmetric matrix, dropping upper");
        for (itercol = 0; itercol < n; itercol++) {
          PASTIX_INT column = itercol;
          if (loc2glob != NULL) column = (*loc2glob)[itercol]-1;
            for (iterrow = (*colptr)[itercol]-1;
                 iterrow < (*colptr)[itercol+1]-1;
                 iterrow++) {
                if ((*row)[iterrow]-1 >= column) {
                  (*row)[index] = (*row)[iterrow];
                  if (avals != NULL)
                    (*avals)[index] = (*avals)[iterrow];
                  index++;
                }
            }
            (*colptr)[itercol] = lastindex+1;
            lastindex = index;
        }
        (*colptr)[n] = lastindex+1;
        MALLOC_EXTERN(tmprows, lastindex, PASTIX_INT);
        memcpy(tmprows, (*row),   lastindex*sizeof(PASTIX_INT));
        free((*row));
        (*row) = tmprows;
        if (avals != NULL)
          {
            MALLOC_EXTERN(tmpvals, lastindex, PASTIX_FLOAT);
            memcpy(tmpvals, (*avals), lastindex*sizeof(PASTIX_FLOAT));
            free((*avals));
            (*avals) = tmpvals;
          }
      }
    else
      {
        if ( ( flagsym == API_SYM_YES ||
               flagsym == API_SYM_HER ) &&
             cnt_lower != 0 && cnt_upper != 0 )
          {
            errorPrint("Only lower or upper part should be given (lower %d upper %d diag %d)",
                       cnt_lower, cnt_upper, cnt_diag);
            RETURN_ERROR(MATRIX_ERR);
          }
      }
  }


  /* Pre-conditionnement mc64 */
#ifdef SCALING
#  ifdef MC64
  if (sizeof(int) != sizeof(PASTIX_INT))
    {
      errorPrint("MC64 only works with classical integers\n");
      RETURN_ERROR(INTEGER_TYPE_ERR);
    }

  errorPrintW("NOT TESTED");
  for (i = 0; i > (*colptr)[n]-1; i++)
    if ((*row)[i] == 0)
      errorPrint("Et Merde\n");

  if ((flagcor == API_YES) && (iparm[IPARM_MC64] == 1))
    {
      PASTIX_INT    job;
      PASTIX_INT    m      = n;
      PASTIX_INT    ne     = (*colptr)[n]-1;
      PASTIX_INT    num;
      PASTIX_INT   *p, *ip;
      PASTIX_INT    liw    = 3*m+2*n+ne;
      PASTIX_INT   *iw;
      PASTIX_INT    ldw    = n+3*m+ne;
      PASTIX_INT    nicntl = 10;
      PASTIX_INT    ncntl  = 10;
      PASTIX_INT   *icntl;
      PASTIX_INT    info;
      PASTIX_FLOAT *dw;
      PASTIX_FLOAT *cntl;

      print_onempi("%s", "Preconditioning...\n");

      MALLOC_INTERN(p,     m,      PASTIX_INT);
      MALLOC_INTERN(ip,    m,      PASTIX_INT);
      MALLOC_INTERN(iw,    liw,    PASTIX_INT);
      MALLOC_INTERN(dw,    ldw,    PASTIX_FLOAT);
      MALLOC_INTERN(icntl, nicntl, PASTIX_INT);
      MALLOC_INTERN(cntl,  ncntl,  PASTIX_FLOAT);

      for (i=0;i<m;i++)
        {
          p[i]  = i+1;
          ip[i] = i+1;
        }
      for (i=0;i<ldw;i++)
        dw[i] = 0.0;
      for (i=0;i<liw;i++)
        iw[i] = 0;
      for (i=0;i<nicntl;i++)
        icntl[i] = 0;
      for (i=0;i<ncntl;i++)
        cntl[i] = 0.0;

      /* compute scaling and unsymmetric column permutation */
      print_onempi("%s", "compute scaling and unsymmetric column permutation...\n");

      FORTRAN_CALL(mc64id)(icntl,cntl);

      cntl[1] = DBL_MAX;
      job     = 6;
      num     = n;

      printf("adresse1 job=%p m=%p n=%p ne=%p ia=%p\n",&job,&m,&n,&ne,*colptr);
      FORTRAN_CALL(mc64ad)(&job,&m,&n,&ne,*colptr,*row,*avals,&num,p,&liw,iw,
                           &ldw,dw,icntl,cntl,&info);
      fprintf(stdout,"return info=%ld (num=%ld n=%ld)\n",info,num,n);
      printf("adresse1 job=%p m=%p n=%p ne=%p ia=%p\n",&job,&m,&n,&ne,*colptr);

      if (num<0 || info<0)
        {
          errorPrint("Error in MC64AD !!!");
          memFree_null(p);
          memFree_null(ip);
          memFree_null(iw);
          memFree_null(dw);
          memFree_null(icntl);
          memFree_null(cntl);
        }
      else
        {
          /* scaling */
          for (i=0;i<m+n;i++)
            dw[i]=exp(dw[i]); /* a_ij := aij * exp(u_i + u_j) */

          print_onempi("%s", "scaling rows...\n");
          CSC_rowScale(n,*colptr,*row,*avals,dw);

          print_onempi("%s", "scaling columns...\n");
          CSC_colScale(n,*colptr,*row,*avals,dw+m);

          /* apply unsymmetric column permutation */
          for (i=0;i<m;i++)
            ip[p[i]-1]=i+1; /* inverse permutation */

          print_onempi("s%", "apply unsymmetric column permutation...\n");
          CSC_colPerm(n,*colptr,*row,*avals,ip);

          memFree_null(p);
          memFree_null(ip);
          memFree_null(iw);
          memFree_null(dw);
          memFree_null(icntl);
          memFree_null(cntl);
        }
    }
#  endif /* MC64 */
#endif /* SCALING */

  /* Symmetrisation du graphe des matrices non-symm�triques */
  if (flagsym == API_SYM_NO)
    {
      if (verb > API_VERBOSE_NOT)
        print_onempi("%s", "Check : Graph symmetry");

      old = (*colptr)[n]-1;

      /* Version distribu�e */
      if ((loc2glob != NULL))
        {
          /* Preserve sorting */
          if (EXIT_SUCCESS == cscd_checksym(n, *colptr, row, avals, *loc2glob, flagcor,
                                            flagalloc, dof,  pastix_comm))
            {
              if (verb > API_VERBOSE_NOT)
                {
                  print_onempi("%s", "\t\tOK\n");
                }
              if (verb > API_VERBOSE_NOT)
                {
                  if (old != ((*colptr)[n] - 1))
                    {
                      fprintf(stdout, "\tAdd %ld null terms on proc %ld\n",
                              (long)((*colptr)[n]-1-old), (long)procnum);
                    }
                }

            }
          else
            {
              if (verb > API_VERBOSE_NOT)
                {
                  print_onempi("%s", "\t\tKO\n");
                }
              RETURN_ERROR(MATRIX_ERR);
            }
        }
      /* Version non distribu�e */
      else
        {
          /* Preserve sorting */
          if (EXIT_SUCCESS == csc_checksym(n, *colptr, row, avals, flagcor, flagalloc, dof))
            {
              if (verb > API_VERBOSE_NOT)
                {
                  print_onempi("%s", "\t\tOK\n");
                }
            }
          else
            {
              if (verb > API_VERBOSE_NOT)
                {
                  print_onempi("%s", "\t\tKO\n");
                }
              RETURN_ERROR(MATRIX_ERR);
            }
          if (verb > API_VERBOSE_NOT)
            {
              if (old != ((*colptr)[n] - 1))
                {
                  print_onempi("\tAdd %ld null terms\t\n",(long)((*colptr)[n]-1-old));
                }
            }


        }
    }

  return NO_ERR;
}

/*
  Function: pastix_getLocalUnknownNbr

  Return the node number in the new distribution computed by blend.
  Needs blend to be runned with pastix_data before.

  Parameters:
  pastix_data - Data used for a step by step execution.

  Returns:
  Number of local nodes/columns in new distribution.
*/
PASTIX_INT pastix_getLocalUnknownNbr(pastix_data_t ** pastix_data)
{
  SolverMatrix  * solvmatr = &((*pastix_data)->solvmatr);
  PASTIX_INT index;
  PASTIX_INT nodenbr;

  nodenbr = 0;
  if ((*pastix_data)->intra_node_procnum == 0)
    for (index=0; index<solvmatr->cblknbr; index++)
      {
        nodenbr += solvmatr->cblktab[index].lcolnum-solvmatr->cblktab[index].fcolnum+1;
      }
  return nodenbr;
}

/*
  Function: pastix_getLocalNodeNbr

  Return the node number in the new distribution computed by blend.
  Needs blend to be runned with pastix_data before.

  Parameters:
  pastix_data - Data used for a step by step execution.

  Returns:
  Number of local nodes/columns in new distribution.
*/
PASTIX_INT pastix_getLocalNodeNbr(pastix_data_t ** pastix_data)
{
  SolverMatrix  * solvmatr = &((*pastix_data)->solvmatr);
  PASTIX_INT index;
  PASTIX_INT nodenbr;

  nodenbr = 0;
  if ((*pastix_data)->intra_node_procnum == 0)
    for (index=0; index<solvmatr->cblknbr; index++)
      {
        nodenbr += solvmatr->cblktab[index].lcolnum-solvmatr->cblktab[index].fcolnum+1;
      }
  nodenbr = nodenbr/(*pastix_data)->iparm[IPARM_DOF_NBR];
  return nodenbr;
}
/* qsort int comparison function */
int
cmpint(const void *p1, const void *p2)
{
  const PASTIX_INT *a = (const PASTIX_INT *)p1;
  const PASTIX_INT *b = (const PASTIX_INT *)p2;

  return (int) *a - *b;
}
/*
  Function: pastix_getLocalUnknownLst

  Fill in unknowns with the list of local nodes/clumns.
  Needs nodelst to be allocated with nodenbr*sizeof(pastix_int_t),
  where nodenbr has been computed by <pastix_getLocalUnknownNbr>.

  Parameters:
  pastix_data - Data used for a step by step execution.
  nodelst     - An array where to write the list of local nodes/columns.
*/
PASTIX_INT pastix_getLocalUnknownLst(pastix_data_t **pastix_data,
                              PASTIX_INT            *nodelst)
{

  SolverMatrix  * solvmatr = &((*pastix_data)->solvmatr);
  Order         * ordemesh = &((*pastix_data)->ordemesh);
  PASTIX_INT index;
  PASTIX_INT index2;
  PASTIX_INT index3;
  int dof = (*pastix_data)->iparm[IPARM_DOF_NBR];

  index3 = 0;
  if ((*pastix_data)->intra_node_procnum == 0)
    for (index=0; index<solvmatr->cblknbr; index++)
      {
        for (index2 = solvmatr->cblktab[index].fcolnum;
             index2 < solvmatr->cblktab[index].lcolnum + 1;
             index2++)
          nodelst[index3++] = dof*ordemesh->peritab[(index2-index2%dof)/dof]+1+index2%dof;
      }
  qsort(nodelst, index3, sizeof(PASTIX_INT), cmpint);

  return NO_ERR;
}


/*
  Function: pastix_getLocalNodeLst

  Fill in nodelst with the list of local nodes/clumns.
  Needs nodelst to be allocated with nodenbr*sizeof(pastix_int_t),
  where nodenbr has been computed by <pastix_getLocalNodeNbr>.

  Parameters:
  pastix_data - Data used for a step by step execution.
  nodelst     - An array where to write the list of local nodes/columns.
*/
PASTIX_INT pastix_getLocalNodeLst(pastix_data_t **pastix_data,
                           PASTIX_INT            *nodelst)
{

  SolverMatrix  * solvmatr = &((*pastix_data)->solvmatr);
  Order         * ordemesh = &((*pastix_data)->ordemesh);
  PASTIX_INT index;
  PASTIX_INT index2;
  PASTIX_INT index3;
  int dof = (*pastix_data)->iparm[IPARM_DOF_NBR];

  index3 = 0;
  if ((*pastix_data)->intra_node_procnum == 0)
    for (index=0; index<solvmatr->cblknbr; index++)
      {
        for (index2 = solvmatr->cblktab[index].fcolnum;
             index2 < solvmatr->cblktab[index].lcolnum + 1;
             index2+=dof)
          nodelst[index3++] = ordemesh->peritab[index2/dof]+1;
      }
  qsort(nodelst, index3, sizeof(PASTIX_INT), cmpint);

  return NO_ERR;
}

/*
  Function: pastix_setSchurUnknownList

  Set the list of unknowns to isolate at the end
  of the matrix via permutations.

  Parameters:
  pastix_data - Data used for a step by step execution.
  n           - Number of unknowns.
  list        - List of unknowns.
*/
PASTIX_INT pastix_setSchurUnknownList(pastix_data_t * pastix_data,
                               PASTIX_INT  n,
                               PASTIX_INT *list)
{
  if (pastix_data == NULL)
    return STEP_ORDER_ERR;

  if (n == 0 || list == NULL)
    return BADPARAMETER_ERR;

  pastix_data->nschur = n;
  MALLOC_INTERN(pastix_data->listschur, n, PASTIX_INT);
  memcpy(pastix_data->listschur, list, n*sizeof(PASTIX_INT));
  return NO_ERR;
}
/*
  Function: pastix_getSchurLocalNodeNbr

  Compute the number of nodes in the local part of the Schur.

  Parameters:
  pastix_data - Common data structure for PaStiX calls.
  nodeNbr     - (out) Number of nodes in schur (local).

  Returns:
  NO_ERR      - For the moment

  TODO: Error management.
*/
PASTIX_INT pastix_getSchurLocalNodeNbr(pastix_data_t * pastix_data, PASTIX_INT * nodeNbr)
{
  SolverMatrix * datacode = &(pastix_data->solvmatr);
  int            owner = API_NO;
  PASTIX_INT            cblk;

  if (SOLV_TASKNBR > 0)
    {
      cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
      if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
        {
          owner = API_YES;
        }

    }

  if (owner == API_YES)
    {
      *nodeNbr = SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1;
    }
  else
    {
      *nodeNbr = 0;
    }
  return NO_ERR;
}

/*
  Function: pastix_getSchurLocalUnkownNbr

  Compute the number of unknowns in the local part of the Schur.

  Parameters:
  pastix_data - Common data structure for PaStiX calls.
  unknownNbr  - (out) Number of unknowns in schur (local).

  Returns:
  NO_ERR      - For the moment

  TODO: Error management.
*/
PASTIX_INT pastix_getSchurLocalUnkownNbr(pastix_data_t * pastix_data, PASTIX_INT * unknownNbr)
{
  SolverMatrix * datacode = &(pastix_data->solvmatr);
  int            owner = API_NO;
  PASTIX_INT            cblk;

  if (SOLV_TASKNBR > 0)
    {
      cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
      if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
        {
          owner = API_YES;
        }

    }

  if (owner == API_YES)
    {
      *unknownNbr = (SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1)*pastix_data->iparm[IPARM_DOF_NBR];
    }
  else
    {
      *unknownNbr = 0;
    }
  return NO_ERR;
}

/*
  Function: pastix_getSchurLocalNodeList

  Compute the list of nodes in the local part of the Schur.

  Parameters:
  pastix_data - Common data structure for PaStiX calls.
  nodes     - (out) Nodes in schur (local).

  Returns:
  NO_ERR      - For the moment

  TODO: Error management.
*/
PASTIX_INT pastix_getSchurLocalNodeList(pastix_data_t * pastix_data, PASTIX_INT * nodes)
{
  SolverMatrix * datacode = NULL;
  Order        * ordemesh = NULL;
  int            owner = API_NO;
  PASTIX_INT            cblk;
  PASTIX_INT            intern_index;
  PASTIX_INT            dof;
  PASTIX_INT            intern_index_dof;
  datacode = &(pastix_data->solvmatr);
  ordemesh = &(pastix_data->ordemesh);

  if (SOLV_TASKNBR > 0)
    {
      cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
      if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
        {
          owner = API_YES;
        }

    }

  if (owner == API_YES)
    {
      PASTIX_INT iter;
      for (iter = 0; iter < SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1; iter+=pastix_data->iparm[IPARM_DOF_NBR])
        {
          intern_index = iter + SYMB_FCOLNUM(cblk);
          dof = intern_index % pastix_data->iparm[IPARM_DOF_NBR];
          intern_index_dof = (intern_index - dof) / pastix_data->iparm[IPARM_DOF_NBR];
          nodes[iter/pastix_data->iparm[IPARM_DOF_NBR]] = ordemesh->peritab[intern_index_dof];
        }
    }

  return NO_ERR;
}


/*
  Function: pastix_getSchurLocalUnkownList

  Compute the list of unknowns in the local part of the Schur.

  Parameters:
  pastix_data - Common data structure for PaStiX calls.
  unknowns    - (out) Unknowns in schur (local).

  Returns:
  NO_ERR      - For the moment

  TODO: Error management.
*/
PASTIX_INT pastix_getSchurLocalUnknownList(pastix_data_t * pastix_data, PASTIX_INT * unknowns)
{
  SolverMatrix * datacode = NULL;
  Order        * ordemesh = NULL;
  int            owner = API_NO;
  PASTIX_INT            cblk;
  PASTIX_INT            intern_index;
  PASTIX_INT            dof;
  PASTIX_INT            intern_index_dof;
  datacode = &(pastix_data->solvmatr);
  ordemesh = &(pastix_data->ordemesh);

  if (SOLV_TASKNBR > 0)
    {
      cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
      if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
        {
          owner = API_YES;
        }

    }

  if (owner == API_YES)
    {
      PASTIX_INT iter;
      for (iter = 0; iter < SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1; iter++)
        {
          intern_index = iter + SYMB_FCOLNUM(cblk);
          dof = intern_index % pastix_data->iparm[IPARM_DOF_NBR];
          intern_index_dof = (intern_index - dof) / pastix_data->iparm[IPARM_DOF_NBR];
          unknowns[iter] = (ordemesh->peritab[intern_index_dof]*pastix_data->iparm[IPARM_DOF_NBR])+dof;
        }
    }

  return NO_ERR;
}


/*
  Function: pastix_getSchurLocalUnkownList

  Give user memory area to store schur in PaStiX.

  Parameters:
  pastix_data - Common data structure for PaStiX calls.
  array       - Memory area to store the schur.

  Returns:
  NO_ERR      - For the moment

  TODO: Error management.
*/
PASTIX_INT pastix_setSchurArray(pastix_data_t * pastix_data, PASTIX_FLOAT * array)
{
  pastix_data->schur_tab = array;
  pastix_data->schur_tab_set = API_YES;
  return NO_ERR;
}
/*
  Function: pastix_getSchur

  Get the Schur complement from PaStiX.

  Schur complement is a dense block in a
  column scheme.

  Parameters:
  pastix_data - Data used for a step by step execution.
  schur - Array to fill-in with Schur complement.

*/
PASTIX_INT pastix_getSchur(pastix_data_t * pastix_data,
                    PASTIX_FLOAT * schur)
{
  SolverMatrix * datacode = &(pastix_data->solvmatr);
  int            owner = API_NO;
  PASTIX_INT            send[2];
  PASTIX_INT            recv[2];
  PASTIX_INT            cblk;

  if (SOLV_TASKNBR > 0)
    {
      cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
      if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
        {
          owner = API_YES;
        }

    }

  if (owner == API_YES)
    {
      PASTIX_INT coefnbr  = SOLV_STRIDE(cblk) * (SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1);
      memcpy(schur, SOLV_COEFTAB(cblk), coefnbr*sizeof(PASTIX_FLOAT));
      send[0] = coefnbr;
      send[1] = SOLV_PROCNUM;

      MPI_Allreduce(&send, &recv, 2, COMM_INT, MPI_SUM, pastix_data->pastix_comm);
    }
  else
    {
      send[0] = 0;
      send[1] = 0;

      MPI_Allreduce(&send, &recv, 2, COMM_INT, MPI_SUM, pastix_data->pastix_comm);
    }
  MPI_Bcast(schur, recv[0], COMM_FLOAT, recv[1], pastix_data->pastix_comm);
  return NO_ERR;
}
/*
 * Function: pastix_checkMatrix
 *
 * Check the matrix :
 * - Renumbers in Fortran numerotation (base 1) if needed (base 0)
 * - Check that the matrix contains no doubles,  with flagcor == API_YES,
 *   correct it.
 * - Can scale the matrix if compiled with -DMC64 -DSCALING (untested)
 * - Checks the symetry of the graph in non symmetric mode.
 *   With non distributed matrices, with flagcor == API_YES,
 *   correct the matrix.
 * - sort the CSC.
 *
 * Parameters:
 *   pastix_comm - PaStiX MPI communicator
 *   verb        - Level of prints (API_VERBOSE_[NOT|NO|YES])
 *   flagsym     - Indicate if the given matrix is symetric
 *                 (API_SYM_YES or API_SYM_NO)
 *   flagcor     - Indicate if we permit the function to reallocate the matrix.
 *   n           - Number of local columns.
 *   colptr      - First element of each row in *row* and *avals*.
 *   row         - Row of each element of the matrix.
 *   avals       - Value of each element of the matrix.
 *   loc2glob    - Global column number of local columns
 *                 (NULL if not distributed).
 *   dof         - Number of degrees of freedom.
 */
PASTIX_INT pastix_checkMatrix(MPI_Comm pastix_comm,
                       PASTIX_INT      verb,
                       PASTIX_INT      flagsym,
                       PASTIX_INT      flagcor,
                       PASTIX_INT      n,
                       PASTIX_INT    **colptr,
                       PASTIX_INT    **row,
                       PASTIX_FLOAT  **avals,
                       PASTIX_INT    **loc2glob,
                       PASTIX_INT      dof)
{
  return pastix_checkMatrix_int(pastix_comm,
                                verb,
                                flagsym,
                                flagcor,
                                n,
                                colptr,
                                row,
                                avals,
                                loc2glob,
                                dof,
                                API_NO);
}

#define pastix_getMemoryUsage PASTIX_EXTERN_F(pastix_getMemoryUsage)
unsigned long pastix_getMemoryUsage() {
#ifdef MEMORY_USAGE
  return memAllocGetCurrent();
#else
  return -1;
#endif
}

#define pastix_getMaxMemoryUsage PASTIX_EXTERN_F(pastix_getMaxMemoryUsage)
unsigned long pastix_getMaxMemoryUsage() {
#ifdef MEMORY_USAGE
  return memAllocGetMax();
#else
  return -1;
#endif
}
