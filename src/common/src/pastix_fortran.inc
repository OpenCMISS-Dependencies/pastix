!  Copyright 2008 BORDEAUX I UNIVERSITY & INRIA ! **
! ** This file is part of the PaStiX parallel sparse matrix package.
! **
! ** This software is governed by the CeCILL-C license under French law
! ** and abiding by the rules of distribution of free software. You can
! ** use, modify and/or redistribute the software under the terms of the
! ** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
! ** URL: "http://www.cecill.info".
! ** 
! ** As a counterpart to the access to the source code and rights to copy,
! ** modify and redistribute granted by the license, users are provided
! ** only with a limited warranty and the software's author, the holder of
! ** the economic rights, and the successive licensors have only limited
! ** liability.
! ** 
! ** In this respect, the user's attention is drawn to the risks associated
! ** with loading, using, modifying and/or developing or reproducing the
! ** software by the user in light of its specific status of free software,
! ** that may mean that it is complicated to manipulate, and that also
! ** therefore means that it is reserved for developers and experienced
! ** professionals having in-depth computer knowledge. Users are therefore
! ** encouraged to load and test the software's suitability as regards
! ** their requirements in conditions enabling the security of their
! ** systems and/or data to be ensured and, more generally, to use and
! ** operate it in the same conditions as regards security.
! ** 
! ** The fact that you are presently reading this means that you have had
! ** knowledge of the CeCILL-C license and that you accept its terms.
! !   header: api.h
! 
!   Header file containing constants used in PaStiX and provided to users.
! 
!   Authors:
!     Mathieu Faverge - faverge@labri.fr
!     Xavier   Lacoste - lacoste@labri.fr
!     Pierre Ramet     - ramet@labri.fr
! 
! !    enum: IPARM_ACCESS
! 
!    Integer parameters tabular accessors
! 
!    IPARM_MODIFY_PARAMETER      - Indicate if parameters have been set by user             Default: API_YES             IN
!    IPARM_START_TASK            - Indicate the first step to execute (see PaStiX steps)    Default: API_TASK_ORDERING   IN
!    IPARM_END_TASK              - Indicate the last step to execute (see PaStiX steps)     Default: API_TASK_CLEAN      IN
!    IPARM_VERBOSE               - Verbose mode (see Verbose modes)                         Default: API_VERBOSE_NO      IN
!    IPARM_DOF_NBR               - Degree of freedom per node                               Default: 1                   IN
!    IPARM_ITERMAX               - Maximum iteration number for refinement                  Default: 250                 IN
!    IPARM_MATRIX_VERIFICATION   - Check the input matrix                                   Default: API_NO              IN
!    IPARM_MC64                  - MC64 operation <pastix.h> IGNORE                         Default: 0                   IN
!    IPARM_ONLY_RAFF             - Refinement only                                          Default: API_NO              IN
!    IPARM_CSCD_CORRECT          - Indicate if the cscd has been redistributed after blend  Default: API_NO              IN
!    IPARM_NBITER                - Number of iterations performed in refinement       Default: -                   OUT
!    IPARM_TRACEFMT              - Trace format (see Trace modes)                           Default: API_TRACE_PICL      IN
!    IPARM_GRAPHDIST             - Specify if the given graph is distributed or not         Default: API_YES             IN
!    IPARM_AMALGAMATION_LEVEL    - Amalgamation level                                       Default: 5                   IN
!    IPARM_ORDERING              - Choose ordering                                          Default: API_ORDER_SCOTCH    IN
!    IPARM_DEFAULT_ORDERING      - Use default ordering parameters with \scotch{} or \metis{} Default: API_YES             IN
!    IPARM_ORDERING_SWITCH_LEVEL - Ordering switch level    (see \scotch{} User's Guide)    Default: 120                 IN
!    IPARM_ORDERING_CMIN         - Ordering cmin parameter  (see \scotch{} User's Guide)    Default: 0                   IN
!    IPARM_ORDERING_CMAX         - Ordering cmax parameter  (see \scotch{} User's Guide)    Default: 100000              IN
!    IPARM_ORDERING_FRAT         - Ordering frat parameter  (see \scotch{} User's Guide)    Default: 8                   IN
!    IPARM_STATIC_PIVOTING       - Static pivoting                                          Default: -                   OUT
!    IPARM_METIS_PFACTOR         - \metis{} pfactor                                         Default: 0                   IN
!    IPARM_NNZEROS               - Number of nonzero entries in the factorized matrix       Default: -                   OUT
!    IPARM_ALLOCATED_TERMS       - Maximum memory allocated for matrix terms                Default: -                   OUT
!    IPARM_BASEVAL               - Baseval used for the matrix                              Default: 0                   IN
!    IPARM_MIN_BLOCKSIZE         - Minimum block size                                       Default: 60                  IN
!    IPARM_MAX_BLOCKSIZE         - Maximum block size                                       Default: 120                 IN
!    IPARM_SCHUR                 - Schur mode                                               Default: API_NO              IN
!    IPARM_ISOLATE_ZEROS         - Isolate null diagonal terms at the end of the matrix     Default: API_NO              IN
!    IPARM_RHSD_CHECK            - Set to API_NO to avoid RHS redistribution                Default: API_YES             IN
!    IPARM_FACTORIZATION         - Factorization mode (see Factorization modes)             Default: API_FACT_LDLT       IN
!    IPARM_NNZEROS_BLOCK_LOCAL   - Number of nonzero entries in the local block factorized matrix Default: -                   OUT
!    IPARM_CPU_BY_NODE           - Number of CPUs per SMP node                              Default: 0                   IN
!    IPARM_BINDTHRD              - Thread binding mode (see Thread binding modes)           Default: API_BIND_AUTO       IN
!    IPARM_THREAD_NBR            - Number of threads per MPI process                        Default: 1                   IN
!    IPARM_DISTRIBUTION_LEVEL    - Distribution level IGNORE                                Default:                     IN
!    IPARM_LEVEL_OF_FILL         - Level of fill for incomplete factorization               Default: 1                   IN
!    IPARM_IO_STRATEGY           - IO strategy (see Checkpoints modes)                      Default: API_IO_NO           IN
!    IPARM_RHS_MAKING            - Right-hand-side making (see Right-hand-side modes)      Default: API_RHS_B           IN
!    IPARM_REFINEMENT            - Refinement type (see Refinement modes)                   Default: API_RAF_GMRES       IN
!    IPARM_SYM                   - Symmetric matrix mode (see Symmetric modes)              Default: API_SYM_YES         IN
!    IPARM_INCOMPLETE            - Incomplete factorization                                 Default: API_NO              IN
!    IPARM_ABS                   - ABS level (Automatic Blocksize Splitting)                Default: 1                   IN
!    IPARM_ESP                   - ESP (Enhanced Sparse Parallelism)                        Default: API_NO              IN
!    IPARM_GMRES_IM              - GMRES restart parameter                                  Default: 25                  IN
!    IPARM_FREE_CSCUSER          - Free user CSC                                            Default: API_CSC_PRESERVE    IN
!    IPARM_FREE_CSCPASTIX        - Free internal CSC (Use only without call to Refin. step) Default: API_CSC_PRESERVE    IN
!    IPARM_OOC_LIMIT             - Out of core memory limit (Mo)                            Default: 2000                IN
!    IPARM_OOC_THREAD            - Out of core thread number IGNORE                         Default: 1                   IN
!    IPARM_OOC_ID                - Out of core run ID        IGNORE                         Default: -                   OUT
!    IPARM_NB_SMP_NODE_USED      - Number of SMP node used   IGNORE                         Default:                     IN
!    IPARM_THREAD_COMM_MODE      - Threaded communication mode (see Communication modes)    Default: API_THREAD_MULT     IN
!    IPARM_NB_THREAD_COMM        - Number of thread(s) for communication                    Default: 1                   IN
!    IPARM_FILL_MATRIX           - Initialize matrix coefficients (for test only)  IGNORE   Default:                     IN
!    IPARM_INERTIA               - Return the inertia (symmetric matrix without pivoting)   Default: -                   OUT
!    IPARM_ESP_NBTASKS           - Return the number of tasks generated by ESP              Default: -                   OUT
!    IPARM_ESP_THRESHOLD         - Minimal block sizee to switch in ESP mode (128 * 128)    Default: 16384               IN
!    IPARM_DOF_COST              - Degree of freedom for cost computation (If different from IPARM_DOF_NBR) Default: 0                    IN
!    IPARM_MURGE_REFINEMENT      - Enable refinement in MURGE                               Default: API_YES             IN
!    IPARM_STARPU                - Use StarPU runtime                                       Default: API_NO              IN
!    IPARM_AUTOSPLIT_COMM        - Automaticaly split communicator to have one MPI task by node             Default: API_NO               IN
!    IPARM_FLOAT                 - Indicate the floating point type  IGNORE                 Default: -                   INOUT
!    IPARM_PID                   - Pid of the first process (used for naming the log directory) Default: -1                  OUT
!    IPARM_ERROR_NUMBER          - Return value                                             Default: -                   OUT
!    IPARM_CUDA_NBR              - Number of cuda devices                                   Default: 0                   IN
!    IPARM_TRANSPOSE_SOLVE       - Use transposed matrix during solve                       Default: API_NO              IN
!    IPARM_STARPU_CTX_DEPTH      - Tree depth of the contexts given to StarPU               Default:3                    IN
!    IPARM_STARPU_CTX_NBR        - Number of contexts created                               Default:-1                   INOUT
!    IPARM_PRODUCE_STATS         - Compute some statistiques (such as precision error)      Default:API_NO               IN
!    IPARM_GPU_CRITERIUM         - Criterium for sorting GPU                                Default:0                    IN
!    IPARM_MURGE_MAY_REFINE      - Enable refinement in MURGE                               Default: API_NO             IN
!    IPARM_SIZE                  - Iparm Size                IGNORE                         Default:                     IN
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MODIFY_PARAMETER         = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_START_TASK               = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_END_TASK                 = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_VERBOSE                  = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_DOF_NBR                  = 5
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ITERMAX                  = 6
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MATRIX_VERIFICATION      = 7
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MC64                     = 8
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ONLY_RAFF                = 9
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_CSCD_CORRECT             = 10
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_NBITER                   = 11
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_TRACEFMT                 = 12
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_GRAPHDIST                = 13
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_AMALGAMATION_LEVEL       = 14
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ORDERING                 = 15
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_DEFAULT_ORDERING         = 16
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ORDERING_SWITCH_LEVEL    = 17
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ORDERING_CMIN            = 18
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ORDERING_CMAX            = 19
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ORDERING_FRAT            = 20
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_STATIC_PIVOTING          = 21
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_METIS_PFACTOR            = 22
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_NNZEROS                  = 23
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ALLOCATED_TERMS          = 24
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_BASEVAL                  = 25
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MIN_BLOCKSIZE            = 26
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MAX_BLOCKSIZE            = 27
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_SCHUR                    = 28
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ISOLATE_ZEROS            = 29
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_RHSD_CHECK               = 30
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_FACTORIZATION            = 31
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_NNZEROS_BLOCK_LOCAL      = 32
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_CPU_BY_NODE              = 33
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_BINDTHRD                 = 34
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_THREAD_NBR               = 35
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_DISTRIBUTION_LEVEL       = 36
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_LEVEL_OF_FILL            = 37
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_IO_STRATEGY              = 38
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_RHS_MAKING               = 39
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_REFINEMENT               = 40
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_SYM                      = 41
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_INCOMPLETE               = 42
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ABS                      = 43
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ESP                      = 44
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_GMRES_IM                 = 45
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_FREE_CSCUSER             = 46
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_FREE_CSCPASTIX           = 47
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_OOC_LIMIT                = 48
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_OOC_THREAD               = 49
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_OOC_ID                   = 50
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_NB_SMP_NODE_USED         = 51
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_THREAD_COMM_MODE         = 52
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_NB_THREAD_COMM           = 53
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_FILL_MATRIX              = 54
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_INERTIA                  = 55
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ESP_NBTASKS              = 56
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ESP_THRESHOLD            = 57
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_DOF_COST                 = 58
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MURGE_REFINEMENT         = 59
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_STARPU                   = 60
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_AUTOSPLIT_COMM           = 61
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_FLOAT                    = 62
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_PID                      = 63
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ERROR_NUMBER             = 64
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_CUDA_NBR                 = 65
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_TRANSPOSE_SOLVE          = 66
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_STARPU_CTX_DEPTH         = 67
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_STARPU_CTX_NBR           = 68
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_PRODUCE_STATS            = 69
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_GPU_CRITERIUM            = 70
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MURGE_MAY_REFINE         = 71
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_SIZE                     = 128
! !    Enum: DPARM_ACCESS
! 
!    Floating point parameters tabular accossors
! 
!    DPARM_FILL_IN            - Fill-in                                           Default: -                OUT
!    DPARM_MEM_MAX            - Maximum memory (-DMEMORY_USAGE)                   Default: -                OUT
!    DPARM_EPSILON_REFINEMENT - Epsilon for refinement                            Default: 1e^{-12}         IN
!    DPARM_RELATIVE_ERROR     - Relative backward error                           Default: -                OUT
!    DPARM_EPSILON_MAGN_CTRL  - Epsilon for magnitude control                     Default: 1e^{-31}         IN
!    DPARM_ANALYZE_TIME       - Time for Analyse step (wallclock)                 Default: -                OUT
!    DPARM_PRED_FACT_TIME     - Predicted factorization time                      Default: -                OUT
!    DPARM_FACT_TIME          - Time for Numerical Factorization step (wallclock) Default: -                OUT
!    DPARM_SOLV_TIME          - Time for Solve step (wallclock)                   Default: -                OUT
!    DPARM_FACT_FLOPS         - Numerical Factorization flops (rate!)             Default: -                OUT
!    DPARM_SOLV_FLOPS         - Solve flops (rate!)                               Default: -                OUT
!    DPARM_RAFF_TIME          - Time for Refinement step (wallclock)              Default: -                OUT
!    DPARM_SIZE               - Dparm Size         IGNORE                         Default: -                IN
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_FILL_IN                  = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_MEM_MAX                  = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_EPSILON_REFINEMENT       = 6
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_RELATIVE_ERROR           = 7
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_SCALED_RESIDUAL          = 8
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_EPSILON_MAGN_CTRL        = 11
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_ANALYZE_TIME             = 19
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_PRED_FACT_TIME           = 20
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_FACT_TIME                = 21
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_SOLV_TIME                = 22
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_FACT_FLOPS               = 23
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_SOLV_FLOPS               = 24
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_RAFF_TIME                = 25
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_SIZE                     = 64
! !   Enum: API_TASK
! 
!   PaStiX step modes (index IPARM_START_TASK and IPARM_END_TASK)
! 
!   API_TASK_INIT       - Set default parameters
!   API_TASK_ORDERING   - Ordering
!   API_TASK_SYMBFACT   - Symbolic factorization
!   API_TASK_ANALYSE    - Tasks mapping and scheduling
!   API_TASK_NUMFACT    - Numerical factorization
!   API_TASK_SOLVE      - Numerical solve
!   API_TASK_REFINE     - Numerical refinement
!   API_TASK_CLEAN      - Clean

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_INIT                  = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_ORDERING              = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_SYMBFACT              = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_ANALYSE               = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_NUMFACT               = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_SOLVE                 = 5
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_REFINE                = 6
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_CLEAN                 = 7
! !   Enum: API_TASK_OLD
! 
!   API_TASK_SCOTCH     - Ordering
!   API_TASK_FAX        - Symbolic factorization
!   API_TASK_BLEND      - Tasks mapping and scheduling
!   API_TASK_SOPALIN    - Numerical factorization
!   API_TASK_UPDOWN     - Numerical solve
!   API_TASK_REFINEMENT - Numerical refinement

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_SCOTCH                = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_FAX                   = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_BLEND                 = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_SOPALIN               = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_UPDOWN                = 5
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_REFINEMENT            = 6
! !    Enum: API_VERBOSE
! 
!    Verbose modes (index IPARM_VERBOSE)
! 
!    API_VERBOSE_NOT        - Silent mode, no messages
!    API_VERBOSE_NO         - Some messages
!    API_VERBOSE_YES        - Many messages
!    API_VERBOSE_CHATTERBOX - Like a gossip
!    API_VERBOSE_UNBEARABLE - Really talking too much...

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_VERBOSE_NOT                = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_VERBOSE_NO                 = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_VERBOSE_YES                = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_VERBOSE_CHATTERBOX         = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_VERBOSE_UNBEARABLE         = 4
! !   Enum: API_IO
! 
!   Check-points modes (index IPARM_IO)
! 
!   API_IO_NO         - No output or input
!   API_IO_LOAD       - Load ordering during ordering step and symbol matrix instead of symbolic factorisation.
!   API_IO_SAVE       - Save ordering during ordering step and symbol matrix instead of symbolic factorisation.
!   API_IO_LOAD_GRAPH - Load graph during ordering step.
!   API_IO_SAVE_GRAPH - Save graph during ordering step.
!   API_IO_LOAD_CSC   - Load CSC(d) during ordering step.
!   API_IO_SAVE_CSC   - Save CSC(d) during ordering step.

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_NO                      = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_LOAD                    = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_SAVE                    = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_LOAD_GRAPH              = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_SAVE_GRAPH              = 8
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_LOAD_CSC                = 16
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_SAVE_CSC                = 32
! !   Enum: API_RHS
! 
!   Right-hand-side modes (index IPARM_RHS)
! 
!   API_RHS_B - User's right hand side
!   API_RHS_1 - $ \forall i, X_i = 1 $
!   API_RHS_I - $ \forall i, X_i = i $
! 

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RHS_B                      = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RHS_1                      = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RHS_I                      = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RHS_0                      = 3
! !   Enum: API_RAF
! 
!   Refinement modes (index IPARM_REFINEMENT)
! 
!   API_RAF_GMRES   - GMRES
!   API_RAF_GRAD    - Conjugate Gradient ($LL^t$ or $LDL^t$ factorization)
!   API_RAF_PIVOT   - Iterative Refinement (only for $LU$ factorization)
!   API_RAF_BICGSTAB - BICGSTAB

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RAF_GMRES                  = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RAF_GRAD                   = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RAF_PIVOT                  = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RAF_BICGSTAB               = 3
! !   Enum: API_FACT
! 
!   Factorization modes (index IPARM_FACTORISATION)
! 
!   API_FACT_LLT  - $LL^t$ Factorization
!   API_FACT_LDLT - $LDL^t$ Factorization
!   API_FACT_LU   - $LU$ Factorization
!   API_FACT_LDLH - $LDL^h$ hermitian factorization

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_FACT_LLT                   = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_FACT_LDLT                  = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_FACT_LU                    = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_FACT_LDLH                  = 3
! !   Enum: API_SYM
! 
!   Symmetric modes (index IPARM_SYM)
! 
!   API_SYM_YES - Symmetric matrix
!   API_SYM_NO  - Nonsymmetric matrix
!   API_SYM_HER - Hermitian
! 

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_SYM_YES                    = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_SYM_NO                     = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_SYM_HER                    = 2
! !   Enum: API_ERASE_CSC
! 
!   CSC Management modes (index IPARM_FREE_CSCUSER and IPARM_FREE_CSCPASTIX)
! 
!   API_CSC_PRESERVE - Do not free the CSC
!   API_CSC_FREE     - Free the CSC when no longer needed

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_CSC_PRESERVE               = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_CSC_FREE                   = 1
! !   Enum: API_THREAD_MODE
! 
!   Comunication modes (index IPARM_THREAD_COMM_MODE)
! 
!   API_THREAD_MULTIPLE      - All threads communicate.
!   API_THREAD_FUNNELED      - One thread perform all the MPI Calls.
!   API_THREAD_COMM_ONE      - One dedicated communication thread will receive messages.
!   API_THREAD_COMM_DEFINED  - Then number of threads receiving the messages is given by IPARM_NB_THREAD_COMM.
!   API_THREAD_COMM_NBPROC   - One communication thread per computation thread will receive messages.

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_THREAD_MULTIPLE            = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_THREAD_FUNNELED            = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_THREAD_COMM_ONE            = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_THREAD_COMM_DEFINED        = 8
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_THREAD_COMM_NBPROC         = 16
! !   Enum: API_BIND_MODE
! 
!   Thread-binding modes (index IPARM_BINTHRD)
! 
!   API_BIND_NO   - Do not bind thread
!   API_BIND_AUTO - Default binding
!   API_BIND_TAB  - Use vector given by pastix_setBind

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_BIND_NO                    = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_BIND_AUTO                  = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_BIND_TAB                   = 2
! !   Enum: API_BOOLEAN
! 
!   Boolean modes (All boolean except IPARM_SYM)
! 
!   API_NO  - No
!   API_YES - Yes

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_NO                         = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_YES                        = 1
! !   Enum: API_TRACEFMT
! 
!   Trace modes (index IPARM_TRACEFMT)
! 
!   API_TRACE_PICL       - Use PICL trace format
!   API_TRACE_PAJE       - Use Paje trace format
!   API_TRACE_HUMREAD    - Use human-readable text trace format
!   API_TRACE_UNFORMATED - Unformated trace format

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TRACE_PICL                 = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TRACE_PAJE                 = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TRACE_HUMREAD              = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TRACE_UNFORMATED           = 3
! !   Enum: API_ORDER
! 
!   Ordering modes (index IPARM_ORDERING)
! 
!   API_ORDER_SCOTCH   - Use \scotch{} ordering
!   API_ORDER_METIS    - Use \metis{} ordering
!   API_ORDER_PERSONAL - Apply user's permutation
!   API_ORDER_LOAD     - Load ordering from disk

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_ORDER_SCOTCH               = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_ORDER_METIS                = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_ORDER_PERSONAL             = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_ORDER_LOAD                 = 3
! !   Enum: API_FLOAT
! 
!   Ordering modes (index IPARM_ORDERING)
! 
!   API_REALSINGLE    - Use \scotch{} ordering
!   API_REALDOUBLE    - Use \metis{} ordering
!   API_COMPLEXSINGLE - Apply user's permutation
!   API_COMPLEXDOUBLE - Load ordering from disk

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_REALSINGLE                 = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_REALDOUBLE                 = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_COMPLEXSINGLE              = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_COMPLEXDOUBLE              = 3
! !  * Enum: API_GPU_CRITERIUM
!  *
!  * Criterium used to decide to put tasks on GPUs.
!  *
!  * API_GPU_CRITERION_UPDATES  - Number of updates on the panel.
!  * API_GPU_CRITERION_CBLKSIZE - Size of the target panel.
!  * API_GPU_CRITERION_FLOPS    - Number of FLOP involved in updates.
!  * API_GPU_CRITERION_PRIORITY - Priority computed in static scheduler.
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_GPU_CRITERION_UPDATES      = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_GPU_CRITERION_CBLKSIZE     = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_GPU_CRITERION_FLOPS        = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_GPU_CRITERION_PRIORITY     = 3
! !   Enum: MODULES
! 
!   Module Identification number.
! 
!   If an error occurs, error value is set to
!   MODULE + EER_NUMBER.
! 
!   User can catch error by computing iparm[IPARM_ERROR_NUMBER]%100.
! 
!   MODULE can be catch by computing iparm[IPARM_ERROR_NUMBER] - iparm[IPARM_ERROR_NUMBER]%100.
! 
!   MOD_UNKNOWN - Unknown module
!   MOD_SOPALIN - Numerical factorisation module
!   MOD_BLEND   - Analysing module
!   MOD_SCOTCH  - Scotch module
!   MOD_FAX     - Symbolic factorisation module
!   MOD_ORDER   - Order module
!   MOD_COMMON  - Common module
!   MOD_SI      -
!   MOD_GRAPH   - Graph module
!   MOD_SYMBOL  - Symbol structure module
!   MOD_KASS    - Kass module
!   MOD_BUBBLE  - Bubble
!   MOD_MURGE   - Murge
! 
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_UNKNOWN                    = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_SOPALIN                    = 100
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_BLEND                      = 200
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_SCOTCH                     = 300
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_FAX                        = 400
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_ORDER                      = 500
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_COMMON                     = 600
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_SI                         = 700
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_GRAPH                      = 800
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_SYMBOL                     = 900
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_KASS                       = 1000
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_BUBBLE                     = 1100
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_MURGE                      = 1200
!  Enum: ERR_NUMBERS! 
!    Error Numbers
! 
!    NO_ERR             - No error
!    UNKNOWN_ERR        - Unknown error
!    ALLOC_ERR          - Allocation error
!    ASSERT_ERR         - Error in one assertion
!    NOTIMPLEMENTED_ERR - Not implemented feature
!    OUTOFMEMORY_ERR    - Not enough memory (OOC)
!    THREAD_ERR         - Error with threads
!    INTERNAL_ERR       - Internal error
!    BADPARAMETER_ERR   - Bad parameters given
!    FILE_ERR           - Error in In/Out operations
!    BAD_DEFINE_ERROR   - Error with defines during compilation
!    INTEGER_TYPE_ERR   - Error with integer types
!    IO_ERR             - Error with input/output
!    MATRIX_ERR         - Wrongly defined matrix
!    FLOAT_TYPE_ERR     - Wrong type of floating point values
!    STEP_ORDER_ERR     - Error in ordering
!    MPI_ERR            - Error with MPI calls
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: NO_ERR                         = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: UNKNOWN_ERR                    = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: ALLOC_ERR                      = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: ASSERT_ERR                     = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: NOTIMPLEMENTED_ERR             = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: OUTOFMEMORY_ERR                = 5
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: THREAD_ERR                     = 6
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: INTERNAL_ERR                   = 7
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: BADPARAMETER_ERR               = 8
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: FILE_ERR                       = 9
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: BAD_DEFINE_ERR                 = 10
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: INTEGER_TYPE_ERR               = 11
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IO_ERR                         = 12
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MATRIX_ERR                     = 13
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: FLOAT_TYPE_ERR                 = 14
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: STEP_ORDER_ERR                 = 15
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MPI_ERR                        = 16

