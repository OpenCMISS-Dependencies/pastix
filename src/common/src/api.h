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
  header: api.h

  Header file containing constants used in PaStiX and provided to users.

  Authors:
    Mathieu Faverge - faverge@labri.fr
    Xavier   Lacoste - lacoste@labri.fr
    Pierre Ramet     - ramet@labri.fr

 */

#ifndef API_H
#define API_H
/* Acces au tableau iparm*/
/*
   enum: IPARM_ACCESS

   Integer parameters tabular accessors

   IPARM_MODIFY_PARAMETER      - Indicate if parameters have been set by user             Default: API_YES             IN
   IPARM_START_TASK            - Indicate the first step to execute (see PaStiX steps)    Default: API_TASK_ORDERING   IN
   IPARM_END_TASK              - Indicate the last step to execute (see PaStiX steps)     Default: API_TASK_CLEAN      IN
   IPARM_VERBOSE               - Verbose mode (see Verbose modes)                         Default: API_VERBOSE_NO      IN
   IPARM_DOF_NBR               - Degree of freedom per node                               Default: 1                   IN
   IPARM_ITERMAX               - Maximum iteration number for refinement                  Default: 250                 IN
   IPARM_MATRIX_VERIFICATION   - Check the input matrix                                   Default: API_NO              IN
   IPARM_MC64                  - MC64 operation <pastix.h> IGNORE                         Default: 0                   IN
   IPARM_ONLY_RAFF             - Refinement only                                          Default: API_NO              IN
   IPARM_CSCD_CORRECT          - Indicate if the cscd has been redistributed after blend  Default: API_NO              IN
   IPARM_NBITER                - Number of iterations performed in refinement       Default: -                   OUT
   IPARM_TRACEFMT              - Trace format (see Trace modes)                           Default: API_TRACE_PICL      IN
   IPARM_GRAPHDIST             - Specify if the given graph is distributed or not         Default: API_YES             IN
   IPARM_AMALGAMATION_LEVEL    - Amalgamation level                                       Default: 5                   IN
   IPARM_ORDERING              - Choose ordering                                          Default: API_ORDER_SCOTCH    IN
   IPARM_DEFAULT_ORDERING      - Use default ordering parameters with \scotch{} or \metis{} Default: API_YES             IN
   IPARM_ORDERING_SWITCH_LEVEL - Ordering switch level    (see \scotch{} User's Guide)    Default: 120                 IN
   IPARM_ORDERING_CMIN         - Ordering cmin parameter  (see \scotch{} User's Guide)    Default: 0                   IN
   IPARM_ORDERING_CMAX         - Ordering cmax parameter  (see \scotch{} User's Guide)    Default: 100000              IN
   IPARM_ORDERING_FRAT         - Ordering frat parameter  (see \scotch{} User's Guide)    Default: 8                   IN
   IPARM_STATIC_PIVOTING       - Static pivoting                                          Default: -                   OUT
   IPARM_METIS_PFACTOR         - \metis{} pfactor                                         Default: 0                   IN
   IPARM_NNZEROS               - Number of nonzero entries in the factorized matrix       Default: -                   OUT
   IPARM_ALLOCATED_TERMS       - Maximum memory allocated for matrix terms                Default: -                   OUT
   IPARM_BASEVAL               - Baseval used for the matrix                              Default: 0                   IN
   IPARM_MIN_BLOCKSIZE         - Minimum block size                                       Default: 60                  IN
   IPARM_MAX_BLOCKSIZE         - Maximum block size                                       Default: 120                 IN
   IPARM_SCHUR                 - Schur mode                                               Default: API_NO              IN
   IPARM_ISOLATE_ZEROS         - Isolate null diagonal terms at the end of the matrix     Default: API_NO              IN
   IPARM_RHSD_CHECK            - Set to API_NO to avoid RHS redistribution                Default: API_YES             IN
   IPARM_FACTORIZATION         - Factorization mode (see Factorization modes)             Default: API_FACT_LDLT       IN
   IPARM_NNZEROS_BLOCK_LOCAL   - Number of nonzero entries in the local block factorized matrix Default: -                   OUT
   IPARM_CPU_BY_NODE           - Number of CPUs per SMP node                              Default: 0                   IN
   IPARM_BINDTHRD              - Thread binding mode (see Thread binding modes)           Default: API_BIND_AUTO       IN
   IPARM_THREAD_NBR            - Number of threads per MPI process                        Default: 1                   IN
   IPARM_DISTRIBUTION_LEVEL    - Distribution level IGNORE                                Default:                     IN
   IPARM_LEVEL_OF_FILL         - Level of fill for incomplete factorization               Default: 1                   IN
   IPARM_IO_STRATEGY           - IO strategy (see Checkpoints modes)                      Default: API_IO_NO           IN
   IPARM_RHS_MAKING            - Right-hand-side making (see Right-hand-side modes)      Default: API_RHS_B           IN
   IPARM_REFINEMENT            - Refinement type (see Refinement modes)                   Default: API_RAF_GMRES       IN
   IPARM_SYM                   - Symmetric matrix mode (see Symmetric modes)              Default: API_SYM_YES         IN
   IPARM_INCOMPLETE            - Incomplete factorization                                 Default: API_NO              IN
   IPARM_ABS                   - ABS level (Automatic Blocksize Splitting)                Default: 1                   IN
   IPARM_ESP                   - ESP (Enhanced Sparse Parallelism)                        Default: API_NO              IN
   IPARM_GMRES_IM              - GMRES restart parameter                                  Default: 25                  IN
   IPARM_FREE_CSCUSER          - Free user CSC                                            Default: API_CSC_PRESERVE    IN
   IPARM_FREE_CSCPASTIX        - Free internal CSC (Use only without call to Refin. step) Default: API_CSC_PRESERVE    IN
   IPARM_OOC_LIMIT             - Out of core memory limit (Mo)                            Default: 2000                IN
   IPARM_OOC_THREAD            - Out of core thread number IGNORE                         Default: 1                   IN
   IPARM_OOC_ID                - Out of core run ID        IGNORE                         Default: -                   OUT
   IPARM_NB_SMP_NODE_USED      - Number of SMP node used   IGNORE                         Default:                     IN
   IPARM_THREAD_COMM_MODE      - Threaded communication mode (see Communication modes)    Default: API_THREAD_MULT     IN
   IPARM_NB_THREAD_COMM        - Number of thread(s) for communication                    Default: 1                   IN
   IPARM_FILL_MATRIX           - Initialize matrix coefficients (for test only)  IGNORE   Default:                     IN
   IPARM_INERTIA               - Return the inertia (symmetric matrix without pivoting)   Default: -                   OUT
   IPARM_ESP_NBTASKS           - Return the number of tasks generated by ESP              Default: -                   OUT
   IPARM_ESP_THRESHOLD         - Minimal block sizee to switch in ESP mode (128 * 128)    Default: 16384               IN
   IPARM_DOF_COST              - Degree of freedom for cost computation (If different from IPARM_DOF_NBR) Default: 0                    IN
   IPARM_MURGE_REFINEMENT      - Enable refinement in MURGE                               Default: API_YES             IN
   IPARM_STARPU                - Use StarPU runtime                                       Default: API_NO              IN
   IPARM_AUTOSPLIT_COMM        - Automaticaly split communicator to have one MPI task by node             Default: API_NO               IN
   IPARM_FLOAT                 - Indicate the floating point type  IGNORE                 Default: -                   INOUT
   IPARM_PID                   - Pid of the first process (used for naming the log directory) Default: -1                  OUT
   IPARM_ERROR_NUMBER          - Return value                                             Default: -                   OUT
   IPARM_CUDA_NBR              - Number of cuda devices                                   Default: 0                   IN
   IPARM_TRANSPOSE_SOLVE       - Use transposed matrix during solve                       Default: API_NO              IN
   IPARM_STARPU_CTX_DEPTH      - Tree depth of the contexts given to StarPU               Default:3                    IN
   IPARM_STARPU_CTX_NBR        - Number of contexts created                               Default:-1                   INOUT
   IPARM_PRODUCE_STATS         - Compute some statistiques (such as precision error)      Default:API_NO               IN
   IPARM_GPU_CRITERIUM         - Criterium for sorting GPU                                Default:0                    IN
   IPARM_MURGE_MAY_REFINE      - Enable refinement in MURGE                               Default: API_NO             IN
   IPARM_SIZE                  - Iparm Size                IGNORE                         Default:                     IN
*/
enum IPARM_ACCESS {
  IPARM_MODIFY_PARAMETER        = 0,
  IPARM_START_TASK              = 1,
  IPARM_END_TASK                = 2,
  IPARM_VERBOSE                 = 3,
  IPARM_DOF_NBR                 = 4,
  IPARM_ITERMAX                 = 5,
  IPARM_MATRIX_VERIFICATION     = 6,
  IPARM_MC64                    = 7,
  IPARM_ONLY_RAFF               = 8,
  IPARM_CSCD_CORRECT            = 9,
  IPARM_NBITER                  = 10,
  IPARM_TRACEFMT                = 11,
  IPARM_GRAPHDIST               = 12,
  IPARM_AMALGAMATION_LEVEL      = 13,
  IPARM_ORDERING                = 14,
  IPARM_DEFAULT_ORDERING        = 15,
  IPARM_ORDERING_SWITCH_LEVEL   = 16,
  IPARM_ORDERING_CMIN           = 17,
  IPARM_ORDERING_CMAX           = 18,
  IPARM_ORDERING_FRAT           = 19,
  IPARM_STATIC_PIVOTING         = 20,
  IPARM_METIS_PFACTOR           = 21,
  IPARM_NNZEROS                 = 22,
  IPARM_ALLOCATED_TERMS         = 23,
  IPARM_BASEVAL                 = 24,
  IPARM_MIN_BLOCKSIZE           = 25,
  IPARM_MAX_BLOCKSIZE           = 26,
  IPARM_SCHUR                   = 27,
  IPARM_ISOLATE_ZEROS           = 28,
  IPARM_RHSD_CHECK              = 29,
  IPARM_FACTORIZATION           = 30,
  IPARM_NNZEROS_BLOCK_LOCAL     = 31,
  IPARM_CPU_BY_NODE             = 32,
  IPARM_BINDTHRD                = 33,
  IPARM_THREAD_NBR              = 34,
  IPARM_DISTRIBUTION_LEVEL      = 35,
  IPARM_LEVEL_OF_FILL           = 36,
  IPARM_IO_STRATEGY             = 37,
  IPARM_RHS_MAKING              = 38,
  IPARM_REFINEMENT              = 39,
  IPARM_SYM                     = 40,
  IPARM_INCOMPLETE              = 41,
  IPARM_ABS                     = 42,
  IPARM_ESP                     = 43,
  IPARM_GMRES_IM                = 44,
  IPARM_FREE_CSCUSER            = 45,
  IPARM_FREE_CSCPASTIX          = 46,
  IPARM_OOC_LIMIT               = 47,
  IPARM_OOC_THREAD              = 48,
  IPARM_OOC_ID                  = 49,
  IPARM_NB_SMP_NODE_USED        = 50,
  IPARM_THREAD_COMM_MODE        = 51,
  IPARM_NB_THREAD_COMM          = 52,
  IPARM_FILL_MATRIX             = 53,
  IPARM_INERTIA                 = 54,
  IPARM_ESP_NBTASKS             = 55,
  IPARM_ESP_THRESHOLD           = 56,
  IPARM_DOF_COST                = 57,
  IPARM_MURGE_REFINEMENT        = 58,
  IPARM_STARPU                  = 59,
  IPARM_AUTOSPLIT_COMM          = 60,
  IPARM_FLOAT                   = 61,
  IPARM_PID                     = 62,
  IPARM_ERROR_NUMBER            = 63,
  IPARM_CUDA_NBR                = 64,
  IPARM_TRANSPOSE_SOLVE         = 65,
  IPARM_STARPU_CTX_DEPTH        = 66,
  IPARM_STARPU_CTX_NBR          = 67,
  IPARM_PRODUCE_STATS           = 68,
  IPARM_GPU_CRITERIUM           = 69,
  IPARM_MURGE_MAY_REFINE        = 70,
  IPARM_SIZE                    = 128  /* Need to be greater or equal to 64 for backward compatibility */
};

/* Acces au tableau dparm */
/*
   Enum: DPARM_ACCESS

   Floating point parameters tabular accossors

   DPARM_FILL_IN            - Fill-in                                           Default: -                OUT
   DPARM_MEM_MAX            - Maximum memory (-DMEMORY_USAGE)                   Default: -                OUT
   DPARM_EPSILON_REFINEMENT - Epsilon for refinement                            Default: 1e^{-12}         IN
   DPARM_RELATIVE_ERROR     - Relative backward error                           Default: -                OUT
   DPARM_EPSILON_MAGN_CTRL  - Epsilon for magnitude control                     Default: 1e^{-31}         IN
   DPARM_ANALYZE_TIME       - Time for Analyse step (wallclock)                 Default: -                OUT
   DPARM_PRED_FACT_TIME     - Predicted factorization time                      Default: -                OUT
   DPARM_FACT_TIME          - Time for Numerical Factorization step (wallclock) Default: -                OUT
   DPARM_SOLV_TIME          - Time for Solve step (wallclock)                   Default: -                OUT
   DPARM_FACT_FLOPS         - Numerical Factorization flops (rate!)             Default: -                OUT
   DPARM_SOLV_FLOPS         - Solve flops (rate!)                               Default: -                OUT
   DPARM_RAFF_TIME          - Time for Refinement step (wallclock)              Default: -                OUT
   DPARM_SIZE               - Dparm Size         IGNORE                         Default: -                IN
 */
enum DPARM_ACCESS {
  DPARM_FILL_IN                 = 1,
  DPARM_MEM_MAX                 = 2,
  DPARM_EPSILON_REFINEMENT      = 5,
  DPARM_RELATIVE_ERROR          = 6,
  DPARM_SCALED_RESIDUAL         = 7,
  DPARM_EPSILON_MAGN_CTRL       = 10,
  DPARM_ANALYZE_TIME            = 18,
  DPARM_PRED_FACT_TIME          = 19,
  DPARM_FACT_TIME               = 20,
  DPARM_SOLV_TIME               = 21,
  DPARM_FACT_FLOPS              = 22,
  DPARM_SOLV_FLOPS              = 23,
  DPARM_RAFF_TIME               = 24,
  DPARM_SIZE                    = 64 /* Need to be greater or equal to 64 for backward compatibility */
};

/** Etapes de résolution de PaStiX */
/*
  Enum: API_TASK

  PaStiX step modes (index IPARM_START_TASK and IPARM_END_TASK)

  API_TASK_INIT       - Set default parameters
  API_TASK_ORDERING   - Ordering
  API_TASK_SYMBFACT   - Symbolic factorization
  API_TASK_ANALYSE    - Tasks mapping and scheduling
  API_TASK_NUMFACT    - Numerical factorization
  API_TASK_SOLVE      - Numerical solve
  API_TASK_REFINE     - Numerical refinement
  API_TASK_CLEAN      - Clean
 */
/* _POS_ 1 */
enum API_TASK {
  API_TASK_INIT       = 0,
  API_TASK_ORDERING   = 1,
  API_TASK_SYMBFACT   = 2,
  API_TASK_ANALYSE    = 3,
  API_TASK_NUMFACT    = 4,
  API_TASK_SOLVE      = 5,
  API_TASK_REFINE     = 6,
  API_TASK_CLEAN      = 7
};

/** Etapes de résolution de PaStiX pour compatibilte avec les anciennes version */
/*
  Enum: API_TASK_OLD

  API_TASK_SCOTCH     - Ordering
  API_TASK_FAX        - Symbolic factorization
  API_TASK_BLEND      - Tasks mapping and scheduling
  API_TASK_SOPALIN    - Numerical factorization
  API_TASK_UPDOWN     - Numerical solve
  API_TASK_REFINEMENT - Numerical refinement
 */
/* _POS_ -1 */
enum API_TASK_OLD {
  API_TASK_SCOTCH     = 1,
  API_TASK_FAX        = 2,
  API_TASK_BLEND      = 3,
  API_TASK_SOPALIN    = 4,
  API_TASK_UPDOWN     = 5,
  API_TASK_REFINEMENT = 6
};

/** Affichage de PaStiX */
/*
   Enum: API_VERBOSE

   Verbose modes (index IPARM_VERBOSE)

   API_VERBOSE_NOT        - Silent mode, no messages
   API_VERBOSE_NO         - Some messages
   API_VERBOSE_YES        - Many messages
   API_VERBOSE_CHATTERBOX - Like a gossip
   API_VERBOSE_UNBEARABLE - Really talking too much...
*/
/* _POS_ 5 */
enum API_VERBOSE {
  API_VERBOSE_NOT        = 0, /* Nothing  */
  API_VERBOSE_NO         = 1, /* Default  */
  API_VERBOSE_YES        = 2, /* Extended */
  API_VERBOSE_CHATTERBOX = 3,
  API_VERBOSE_UNBEARABLE = 4
};

/** Load strategy for graph and ordering */
/*
  Enum: API_IO

  Check-points modes (index IPARM_IO)

  API_IO_NO         - No output or input
  API_IO_LOAD       - Load ordering during ordering step and symbol matrix instead of symbolic factorisation.
  API_IO_SAVE       - Save ordering during ordering step and symbol matrix instead of symbolic factorisation.
  API_IO_LOAD_GRAPH - Load graph during ordering step.
  API_IO_SAVE_GRAPH - Save graph during ordering step.
  API_IO_LOAD_CSC   - Load CSC(d) during ordering step.
  API_IO_SAVE_CSC   - Save CSC(d) during ordering step.
 */
/* _POS_ 6 */
enum API_IO {
  API_IO_NO         = 0,
  API_IO_LOAD       = 1,
  API_IO_SAVE       = 2,
  API_IO_LOAD_GRAPH = 4,
  API_IO_SAVE_GRAPH = 8,
  API_IO_LOAD_CSC   = 16,
  API_IO_SAVE_CSC   = 32
};

/** Genération du second membre */
/*
  Enum: API_RHS

  Right-hand-side modes (index IPARM_RHS)

  API_RHS_B - User's right hand side
  API_RHS_1 - $ \forall i, X_i = 1 $
  API_RHS_I - $ \forall i, X_i = i $

 */
/* _POS_ 7 */
enum API_RHS {
  API_RHS_B = 0, /* Utilisation du second membre fournit */
  API_RHS_1 = 1, /* Utilisation d'un second membre dont tous les coefficients valent 1 */
  API_RHS_I = 2, /* Utilisation d'un second membre tel que RHS(i) = i */
  API_RHS_0 = 3  /* Initialisation en mode ONLY_RAFF d'une solution X0(i) = 0 */
};

/** Type de raffinement utilisé */
/*
  Enum: API_RAF

  Refinement modes (index IPARM_REFINEMENT)

  API_RAF_GMRES   - GMRES
  API_RAF_GRAD    - Conjugate Gradient ($LL^t$ or $LDL^t$ factorization)
  API_RAF_PIVOT   - Iterative Refinement (only for $LU$ factorization)
  API_RAF_BICGSTAB - BICGSTAB
 */
/* _POS_ 8 */
enum API_RAF {
  API_RAF_GMRES   = 0, /* Utilisation de GMRES */
  API_RAF_GRAD    = 1, /* Utilisation du gradient conjugue */
  API_RAF_PIVOT   = 2, /* Utilisation de la methode du pivot */
  API_RAF_BICGSTAB = 3
};

/** Type de facto utilisée (LLT,LDLT,LU)*/
/*
  Enum: API_FACT

  Factorization modes (index IPARM_FACTORISATION)

  API_FACT_LLT  - $LL^t$ Factorization
  API_FACT_LDLT - $LDL^t$ Factorization
  API_FACT_LU   - $LU$ Factorization
  API_FACT_LDLH - $LDL^h$ hermitian factorization
*/
/* _POS_ 4 */
enum API_FACT {
  API_FACT_LLT  = 0, /* Factorisation de Cholesky */
  API_FACT_LDLT = 1, /* Factorisation de Crout */
  API_FACT_LU   = 2, /* Factorisation LU */
  API_FACT_LDLH  = 3
};

/** Matrice symétrique ou non (0 : symétrique, 1 : non) */
/*
  Enum: API_SYM

  Symmetric modes (index IPARM_SYM)

  API_SYM_YES - Symmetric matrix
  API_SYM_NO  - Nonsymmetric matrix
  API_SYM_HER - Hermitian

 */
/* _POS_ 3 */
enum API_SYM {
  API_SYM_YES = 0, /* Matrice symetrique     */
  API_SYM_NO  = 1,  /* Matrice non symetrique */
  API_SYM_HER = 2
};

/** Supressing user CSC(D) when not usefull anymore */
/*
  Enum: API_ERASE_CSC

  CSC Management modes (index IPARM_FREE_CSCUSER and IPARM_FREE_CSCPASTIX)

  API_CSC_PRESERVE - Do not free the CSC
  API_CSC_FREE     - Free the CSC when no longer needed
 */
/* _POS_ 11 */
enum API_ERASE_CSC{
  API_CSC_PRESERVE = 0,
  API_CSC_FREE     = 1
};

/** DMP communication mode */
/*
  Enum: API_THREAD_MODE

  Comunication modes (index IPARM_THREAD_COMM_MODE)

  API_THREAD_MULTIPLE      - All threads communicate.
  API_THREAD_FUNNELED      - One thread perform all the MPI Calls.
  API_THREAD_COMM_ONE      - One dedicated communication thread will receive messages.
  API_THREAD_COMM_DEFINED  - Then number of threads receiving the messages is given by IPARM_NB_THREAD_COMM.
  API_THREAD_COMM_NBPROC   - One communication thread per computation thread will receive messages.
 */
/* _POS_ 9 */
enum API_THREAD_MODE {
  API_THREAD_MULTIPLE      = 1,
  API_THREAD_FUNNELED      = 2,
  API_THREAD_COMM_ONE      = 4,
  API_THREAD_COMM_DEFINED  = 8,
  API_THREAD_COMM_NBPROC   = 16
};

/** Thread binding */
/*
  Enum: API_BIND_MODE

  Thread-binding modes (index IPARM_BINTHRD)

  API_BIND_NO   - Do not bind thread
  API_BIND_AUTO - Default binding
  API_BIND_TAB  - Use vector given by pastix_setBind
*/
/* _POS_ 12 */
enum API_BIND_MODE {
  API_BIND_NO      = 0, /* Do not bind threads                            */
  API_BIND_AUTO    = 1, /* Default thread binding                         */
  API_BIND_TAB     = 2  /* Use tabular given by pastix_setBind to bind */
};

/** Boolean */
/*
  Enum: API_BOOLEAN

  Boolean modes (All boolean except IPARM_SYM)

  API_NO  - No
  API_YES - Yes
 */
/* _POS_ 2 */
enum API_BOOLEAN {
  API_NO  = 0,
  API_YES = 1
};

/** Trace format */
/*
  Enum: API_TRACEFMT

  Trace modes (index IPARM_TRACEFMT)

  API_TRACE_PICL       - Use PICL trace format
  API_TRACE_PAJE       - Use Paje trace format
  API_TRACE_HUMREAD    - Use human-readable text trace format
  API_TRACE_UNFORMATED - Unformated trace format
 */
/* _POS_ 10 */
enum API_TRACEFMT {
  API_TRACE_PICL       = 0, /* Use PICL trace format       */
  API_TRACE_PAJE       = 1, /* Use Paje trace format       */
  API_TRACE_HUMREAD    = 2, /* Use text trace format       */
  API_TRACE_UNFORMATED = 3  /* Use unformated trace format */
};

/*
  Enum: API_ORDER

  Ordering modes (index IPARM_ORDERING)

  API_ORDER_SCOTCH   - Use \scotch{} ordering
  API_ORDER_METIS    - Use \metis{} ordering
  API_ORDER_PERSONAL - Apply user's permutation
  API_ORDER_LOAD     - Load ordering from disk
 */
/* _POS_ 11 */
enum API_ORDER {
  API_ORDER_SCOTCH    = 0,
  API_ORDER_METIS     = 1,
  API_ORDER_PERSONAL  = 2,
  API_ORDER_LOAD      = 3
};

/*
  Enum: API_FLOAT

  Ordering modes (index IPARM_ORDERING)

  API_REALSINGLE    - Use \scotch{} ordering
  API_REALDOUBLE    - Use \metis{} ordering
  API_COMPLEXSINGLE - Apply user's permutation
  API_COMPLEXDOUBLE - Load ordering from disk
 */
/* _POS_ 61 */
enum API_FLOAT {
  API_REALSINGLE    = 0,
  API_REALDOUBLE    = 1,
  API_COMPLEXSINGLE = 2,
  API_COMPLEXDOUBLE = 3
};

/*
 * Enum: API_GPU_CRITERIUM
 *
 * Criterium used to decide to put tasks on GPUs.
 *
 * API_GPU_CRITERION_UPDATES  - Number of updates on the panel.
 * API_GPU_CRITERION_CBLKSIZE - Size of the target panel.
 * API_GPU_CRITERION_FLOPS    - Number of FLOP involved in updates.
 * API_GPU_CRITERION_PRIORITY - Priority computed in static scheduler.
 */
enum API_GPU_CRITERIUM {
  API_GPU_CRITERION_UPDATES  = 0,
  API_GPU_CRITERION_CBLKSIZE = 1,
  API_GPU_CRITERION_FLOPS    = 2,
  API_GPU_CRITERION_PRIORITY = 3
};
/*
  Enum: MODULES

  Module Identification number.

  If an error occurs, error value is set to
  MODULE + EER_NUMBER.

  User can catch error by computing iparm[IPARM_ERROR_NUMBER]%100.

  MODULE can be catch by computing iparm[IPARM_ERROR_NUMBER] - iparm[IPARM_ERROR_NUMBER]%100.

  MOD_UNKNOWN - Unknown module
  MOD_SOPALIN - Numerical factorisation module
  MOD_BLEND   - Analysing module
  MOD_SCOTCH  - Scotch module
  MOD_FAX     - Symbolic factorisation module
  MOD_ORDER   - Order module
  MOD_COMMON  - Common module
  MOD_SI      -
  MOD_GRAPH   - Graph module
  MOD_SYMBOL  - Symbol structure module
  MOD_KASS    - Kass module
  MOD_BUBBLE  - Bubble
  MOD_MURGE   - Murge

*/
enum MODULES {
  MOD_UNKNOWN =    0,
  MOD_SOPALIN =  100,
  MOD_BLEND   =  200,
  MOD_SCOTCH  =  300,
  MOD_FAX     =  400,
  MOD_ORDER   =  500,
  MOD_COMMON  =  600,
  MOD_SI      =  700,
  MOD_GRAPH   =  800,
  MOD_SYMBOL  =  900,
  MOD_KASS    = 1000,
  MOD_BUBBLE  = 1100,
  MOD_MURGE   = 1200
};

/* Enum: ERR_NUMBERS

   Error Numbers

   NO_ERR             - No error
   UNKNOWN_ERR        - Unknown error
   ALLOC_ERR          - Allocation error
   ASSERT_ERR         - Error in one assertion
   NOTIMPLEMENTED_ERR - Not implemented feature
   OUTOFMEMORY_ERR    - Not enough memory (OOC)
   THREAD_ERR         - Error with threads
   INTERNAL_ERR       - Internal error
   BADPARAMETER_ERR   - Bad parameters given
   FILE_ERR           - Error in In/Out operations
   BAD_DEFINE_ERROR   - Error with defines during compilation
   INTEGER_TYPE_ERR   - Error with integer types
   IO_ERR             - Error with input/output
   MATRIX_ERR         - Wrongly defined matrix
   FLOAT_TYPE_ERR     - Wrong type of floating point values
   STEP_ORDER_ERR     - Error in ordering
   MPI_ERR            - Error with MPI calls
*/
/* Need to conserve it MURGE compliant */
enum ERR_NUMBERS {
  NO_ERR             = 0,
  UNKNOWN_ERR        = 1,
  ALLOC_ERR          = 2,
  ASSERT_ERR         = 3,
  NOTIMPLEMENTED_ERR = 4,
  OUTOFMEMORY_ERR    = 5,
  THREAD_ERR         = 6,
  INTERNAL_ERR       = 7,
  BADPARAMETER_ERR   = 8,
  FILE_ERR           = 9,
  BAD_DEFINE_ERR     = 10,
  INTEGER_TYPE_ERR   = 11,
  IO_ERR             = 12,
  MATRIX_ERR         = 13,
  FLOAT_TYPE_ERR     = 14,
  STEP_ORDER_ERR     = 15,
  MPI_ERR            = 16
};

/** Matrix type */
#define MTX_ISSYM(a) ((a)[1]=='S')
#define MTX_ISHER(a) ((a)[1]=='H')
#define MTX_ISCOM(a) ((a)[0]=='C')
#define MTX_ISRHX(a) ((a)[2]=='X')
#define MTX_ISRHS(a) ((a)[0]!='\0')

/* **************************************** */
#endif /* not API_H */
