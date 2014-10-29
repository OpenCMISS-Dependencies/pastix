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
#define IPARM_MODIFY_PARAMETER         INT(1 , KIND=PASTIX_INT_KIND)
#define IPARM_START_TASK               INT(2 , KIND=PASTIX_INT_KIND)
#define IPARM_END_TASK                 INT(3 , KIND=PASTIX_INT_KIND)
#define IPARM_VERBOSE                  INT(4 , KIND=PASTIX_INT_KIND)
#define IPARM_DOF_NBR                  INT(5 , KIND=PASTIX_INT_KIND)
#define IPARM_ITERMAX                  INT(6 , KIND=PASTIX_INT_KIND)
#define IPARM_MATRIX_VERIFICATION      INT(7 , KIND=PASTIX_INT_KIND)
#define IPARM_MC64                     INT(8 , KIND=PASTIX_INT_KIND)
#define IPARM_ONLY_RAFF                INT(9 , KIND=PASTIX_INT_KIND)
#define IPARM_CSCD_CORRECT             INT(10 , KIND=PASTIX_INT_KIND)
#define IPARM_NBITER                   INT(11 , KIND=PASTIX_INT_KIND)
#define IPARM_TRACEFMT                 INT(12 , KIND=PASTIX_INT_KIND)
#define IPARM_GRAPHDIST                INT(13 , KIND=PASTIX_INT_KIND)
#define IPARM_AMALGAMATION_LEVEL       INT(14 , KIND=PASTIX_INT_KIND)
#define IPARM_ORDERING                 INT(15 , KIND=PASTIX_INT_KIND)
#define IPARM_DEFAULT_ORDERING         INT(16 , KIND=PASTIX_INT_KIND)
#define IPARM_ORDERING_SWITCH_LEVEL    INT(17 , KIND=PASTIX_INT_KIND)
#define IPARM_ORDERING_CMIN            INT(18 , KIND=PASTIX_INT_KIND)
#define IPARM_ORDERING_CMAX            INT(19 , KIND=PASTIX_INT_KIND)
#define IPARM_ORDERING_FRAT            INT(20 , KIND=PASTIX_INT_KIND)
#define IPARM_STATIC_PIVOTING          INT(21 , KIND=PASTIX_INT_KIND)
#define IPARM_METIS_PFACTOR            INT(22 , KIND=PASTIX_INT_KIND)
#define IPARM_NNZEROS                  INT(23 , KIND=PASTIX_INT_KIND)
#define IPARM_ALLOCATED_TERMS          INT(24 , KIND=PASTIX_INT_KIND)
#define IPARM_BASEVAL                  INT(25 , KIND=PASTIX_INT_KIND)
#define IPARM_MIN_BLOCKSIZE            INT(26 , KIND=PASTIX_INT_KIND)
#define IPARM_MAX_BLOCKSIZE            INT(27 , KIND=PASTIX_INT_KIND)
#define IPARM_SCHUR                    INT(28 , KIND=PASTIX_INT_KIND)
#define IPARM_ISOLATE_ZEROS            INT(29 , KIND=PASTIX_INT_KIND)
#define IPARM_RHSD_CHECK               INT(30 , KIND=PASTIX_INT_KIND)
#define IPARM_FACTORIZATION            INT(31 , KIND=PASTIX_INT_KIND)
#define IPARM_NNZEROS_BLOCK_LOCAL      INT(32 , KIND=PASTIX_INT_KIND)
#define IPARM_CPU_BY_NODE              INT(33 , KIND=PASTIX_INT_KIND)
#define IPARM_BINDTHRD                 INT(34 , KIND=PASTIX_INT_KIND)
#define IPARM_THREAD_NBR               INT(35 , KIND=PASTIX_INT_KIND)
#define IPARM_DISTRIBUTION_LEVEL       INT(36 , KIND=PASTIX_INT_KIND)
#define IPARM_LEVEL_OF_FILL            INT(37 , KIND=PASTIX_INT_KIND)
#define IPARM_IO_STRATEGY              INT(38 , KIND=PASTIX_INT_KIND)
#define IPARM_RHS_MAKING               INT(39 , KIND=PASTIX_INT_KIND)
#define IPARM_REFINEMENT               INT(40 , KIND=PASTIX_INT_KIND)
#define IPARM_SYM                      INT(41 , KIND=PASTIX_INT_KIND)
#define IPARM_INCOMPLETE               INT(42 , KIND=PASTIX_INT_KIND)
#define IPARM_ABS                      INT(43 , KIND=PASTIX_INT_KIND)
#define IPARM_ESP                      INT(44 , KIND=PASTIX_INT_KIND)
#define IPARM_GMRES_IM                 INT(45 , KIND=PASTIX_INT_KIND)
#define IPARM_FREE_CSCUSER             INT(46 , KIND=PASTIX_INT_KIND)
#define IPARM_FREE_CSCPASTIX           INT(47 , KIND=PASTIX_INT_KIND)
#define IPARM_OOC_LIMIT                INT(48 , KIND=PASTIX_INT_KIND)
#define IPARM_OOC_THREAD               INT(49 , KIND=PASTIX_INT_KIND)
#define IPARM_OOC_ID                   INT(50 , KIND=PASTIX_INT_KIND)
#define IPARM_NB_SMP_NODE_USED         INT(51 , KIND=PASTIX_INT_KIND)
#define IPARM_THREAD_COMM_MODE         INT(52 , KIND=PASTIX_INT_KIND)
#define IPARM_NB_THREAD_COMM           INT(53 , KIND=PASTIX_INT_KIND)
#define IPARM_FILL_MATRIX              INT(54 , KIND=PASTIX_INT_KIND)
#define IPARM_INERTIA                  INT(55 , KIND=PASTIX_INT_KIND)
#define IPARM_ESP_NBTASKS              INT(56 , KIND=PASTIX_INT_KIND)
#define IPARM_ESP_THRESHOLD            INT(57 , KIND=PASTIX_INT_KIND)
#define IPARM_DOF_COST                 INT(58 , KIND=PASTIX_INT_KIND)
#define IPARM_MURGE_REFINEMENT         INT(59 , KIND=PASTIX_INT_KIND)
#define IPARM_STARPU                   INT(60 , KIND=PASTIX_INT_KIND)
#define IPARM_AUTOSPLIT_COMM           INT(61 , KIND=PASTIX_INT_KIND)
#define IPARM_FLOAT                    INT(62 , KIND=PASTIX_INT_KIND)
#define IPARM_PID                      INT(63 , KIND=PASTIX_INT_KIND)
#define IPARM_ERROR_NUMBER             INT(64 , KIND=PASTIX_INT_KIND)
#define IPARM_CUDA_NBR                 INT(65 , KIND=PASTIX_INT_KIND)
#define IPARM_TRANSPOSE_SOLVE          INT(66 , KIND=PASTIX_INT_KIND)
#define IPARM_STARPU_CTX_DEPTH         INT(67 , KIND=PASTIX_INT_KIND)
#define IPARM_STARPU_CTX_NBR           INT(68 , KIND=PASTIX_INT_KIND)
#define IPARM_PRODUCE_STATS            INT(69 , KIND=PASTIX_INT_KIND)
#define IPARM_GPU_CRITERIUM            INT(70 , KIND=PASTIX_INT_KIND)
#define IPARM_MURGE_MAY_REFINE         INT(71 , KIND=PASTIX_INT_KIND)
#define IPARM_SIZE                     INT(128 , KIND=PASTIX_INT_KIND)
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
#define DPARM_FILL_IN                  INT(2 , KIND=PASTIX_INT_KIND)
#define DPARM_MEM_MAX                  INT(3 , KIND=PASTIX_INT_KIND)
#define DPARM_EPSILON_REFINEMENT       INT(6 , KIND=PASTIX_INT_KIND)
#define DPARM_RELATIVE_ERROR           INT(7 , KIND=PASTIX_INT_KIND)
#define DPARM_SCALED_RESIDUAL          INT(8 , KIND=PASTIX_INT_KIND)
#define DPARM_EPSILON_MAGN_CTRL        INT(11 , KIND=PASTIX_INT_KIND)
#define DPARM_ANALYZE_TIME             INT(19 , KIND=PASTIX_INT_KIND)
#define DPARM_PRED_FACT_TIME           INT(20 , KIND=PASTIX_INT_KIND)
#define DPARM_FACT_TIME                INT(21 , KIND=PASTIX_INT_KIND)
#define DPARM_SOLV_TIME                INT(22 , KIND=PASTIX_INT_KIND)
#define DPARM_FACT_FLOPS               INT(23 , KIND=PASTIX_INT_KIND)
#define DPARM_SOLV_FLOPS               INT(24 , KIND=PASTIX_INT_KIND)
#define DPARM_RAFF_TIME                INT(25 , KIND=PASTIX_INT_KIND)
#define DPARM_SIZE                     INT(64 , KIND=PASTIX_INT_KIND)
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

#define API_TASK_INIT                  INT(0 , KIND=PASTIX_INT_KIND)
#define API_TASK_ORDERING              INT(1 , KIND=PASTIX_INT_KIND)
#define API_TASK_SYMBFACT              INT(2 , KIND=PASTIX_INT_KIND)
#define API_TASK_ANALYSE               INT(3 , KIND=PASTIX_INT_KIND)
#define API_TASK_NUMFACT               INT(4 , KIND=PASTIX_INT_KIND)
#define API_TASK_SOLVE                 INT(5 , KIND=PASTIX_INT_KIND)
#define API_TASK_REFINE                INT(6 , KIND=PASTIX_INT_KIND)
#define API_TASK_CLEAN                 INT(7 , KIND=PASTIX_INT_KIND)
! !   Enum: API_TASK_OLD
! 
!   API_TASK_SCOTCH     - Ordering
!   API_TASK_FAX        - Symbolic factorization
!   API_TASK_BLEND      - Tasks mapping and scheduling
!   API_TASK_SOPALIN    - Numerical factorization
!   API_TASK_UPDOWN     - Numerical solve
!   API_TASK_REFINEMENT - Numerical refinement

#define API_TASK_SCOTCH                INT(1 , KIND=PASTIX_INT_KIND)
#define API_TASK_FAX                   INT(2 , KIND=PASTIX_INT_KIND)
#define API_TASK_BLEND                 INT(3 , KIND=PASTIX_INT_KIND)
#define API_TASK_SOPALIN               INT(4 , KIND=PASTIX_INT_KIND)
#define API_TASK_UPDOWN                INT(5 , KIND=PASTIX_INT_KIND)
#define API_TASK_REFINEMENT            INT(6 , KIND=PASTIX_INT_KIND)
! !    Enum: API_VERBOSE
! 
!    Verbose modes (index IPARM_VERBOSE)
! 
!    API_VERBOSE_NOT        - Silent mode, no messages
!    API_VERBOSE_NO         - Some messages
!    API_VERBOSE_YES        - Many messages
!    API_VERBOSE_CHATTERBOX - Like a gossip
!    API_VERBOSE_UNBEARABLE - Really talking too much...

#define API_VERBOSE_NOT                INT(0 , KIND=PASTIX_INT_KIND)
#define API_VERBOSE_NO                 INT(1 , KIND=PASTIX_INT_KIND)
#define API_VERBOSE_YES                INT(2 , KIND=PASTIX_INT_KIND)
#define API_VERBOSE_CHATTERBOX         INT(3 , KIND=PASTIX_INT_KIND)
#define API_VERBOSE_UNBEARABLE         INT(4 , KIND=PASTIX_INT_KIND)
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

#define API_IO_NO                      INT(0 , KIND=PASTIX_INT_KIND)
#define API_IO_LOAD                    INT(1 , KIND=PASTIX_INT_KIND)
#define API_IO_SAVE                    INT(2 , KIND=PASTIX_INT_KIND)
#define API_IO_LOAD_GRAPH              INT(4 , KIND=PASTIX_INT_KIND)
#define API_IO_SAVE_GRAPH              INT(8 , KIND=PASTIX_INT_KIND)
#define API_IO_LOAD_CSC                INT(16 , KIND=PASTIX_INT_KIND)
#define API_IO_SAVE_CSC                INT(32 , KIND=PASTIX_INT_KIND)
! !   Enum: API_RHS
! 
!   Right-hand-side modes (index IPARM_RHS)
! 
!   API_RHS_B - User's right hand side
!   API_RHS_1 - $ \forall i, X_i = 1 $
!   API_RHS_I - $ \forall i, X_i = i $
! 

#define API_RHS_B                      INT(0 , KIND=PASTIX_INT_KIND)
#define API_RHS_1                      INT(1 , KIND=PASTIX_INT_KIND)
#define API_RHS_I                      INT(2 , KIND=PASTIX_INT_KIND)
#define API_RHS_0                      INT(3 , KIND=PASTIX_INT_KIND)
! !   Enum: API_RAF
! 
!   Refinement modes (index IPARM_REFINEMENT)
! 
!   API_RAF_GMRES   - GMRES
!   API_RAF_GRAD    - Conjugate Gradient ($LL^t$ or $LDL^t$ factorization)
!   API_RAF_PIVOT   - Iterative Refinement (only for $LU$ factorization)
!   API_RAF_BICGSTAB - BICGSTAB

#define API_RAF_GMRES                  INT(0 , KIND=PASTIX_INT_KIND)
#define API_RAF_GRAD                   INT(1 , KIND=PASTIX_INT_KIND)
#define API_RAF_PIVOT                  INT(2 , KIND=PASTIX_INT_KIND)
#define API_RAF_BICGSTAB               INT(3 , KIND=PASTIX_INT_KIND)
! !   Enum: API_FACT
! 
!   Factorization modes (index IPARM_FACTORISATION)
! 
!   API_FACT_LLT  - $LL^t$ Factorization
!   API_FACT_LDLT - $LDL^t$ Factorization
!   API_FACT_LU   - $LU$ Factorization
!   API_FACT_LDLH - $LDL^h$ hermitian factorization

#define API_FACT_LLT                   INT(0 , KIND=PASTIX_INT_KIND)
#define API_FACT_LDLT                  INT(1 , KIND=PASTIX_INT_KIND)
#define API_FACT_LU                    INT(2 , KIND=PASTIX_INT_KIND)
#define API_FACT_LDLH                  INT(3 , KIND=PASTIX_INT_KIND)
! !   Enum: API_SYM
! 
!   Symmetric modes (index IPARM_SYM)
! 
!   API_SYM_YES - Symmetric matrix
!   API_SYM_NO  - Nonsymmetric matrix
!   API_SYM_HER - Hermitian
! 

#define API_SYM_YES                    INT(0 , KIND=PASTIX_INT_KIND)
#define API_SYM_NO                     INT(1 , KIND=PASTIX_INT_KIND)
#define API_SYM_HER                    INT(2 , KIND=PASTIX_INT_KIND)
! !   Enum: API_ERASE_CSC
! 
!   CSC Management modes (index IPARM_FREE_CSCUSER and IPARM_FREE_CSCPASTIX)
! 
!   API_CSC_PRESERVE - Do not free the CSC
!   API_CSC_FREE     - Free the CSC when no longer needed

#define API_CSC_PRESERVE               INT(0 , KIND=PASTIX_INT_KIND)
#define API_CSC_FREE                   INT(1 , KIND=PASTIX_INT_KIND)
! !   Enum: API_THREAD_MODE
! 
!   Comunication modes (index IPARM_THREAD_COMM_MODE)
! 
!   API_THREAD_MULTIPLE      - All threads communicate.
!   API_THREAD_FUNNELED      - One thread perform all the MPI Calls.
!   API_THREAD_COMM_ONE      - One dedicated communication thread will receive messages.
!   API_THREAD_COMM_DEFINED  - Then number of threads receiving the messages is given by IPARM_NB_THREAD_COMM.
!   API_THREAD_COMM_NBPROC   - One communication thread per computation thread will receive messages.

#define API_THREAD_MULTIPLE            INT(1 , KIND=PASTIX_INT_KIND)
#define API_THREAD_FUNNELED            INT(2 , KIND=PASTIX_INT_KIND)
#define API_THREAD_COMM_ONE            INT(4 , KIND=PASTIX_INT_KIND)
#define API_THREAD_COMM_DEFINED        INT(8 , KIND=PASTIX_INT_KIND)
#define API_THREAD_COMM_NBPROC         INT(16 , KIND=PASTIX_INT_KIND)
! !   Enum: API_BIND_MODE
! 
!   Thread-binding modes (index IPARM_BINTHRD)
! 
!   API_BIND_NO   - Do not bind thread
!   API_BIND_AUTO - Default binding
!   API_BIND_TAB  - Use vector given by pastix_setBind

#define API_BIND_NO                    INT(0 , KIND=PASTIX_INT_KIND)
#define API_BIND_AUTO                  INT(1 , KIND=PASTIX_INT_KIND)
#define API_BIND_TAB                   INT(2 , KIND=PASTIX_INT_KIND)
! !   Enum: API_BOOLEAN
! 
!   Boolean modes (All boolean except IPARM_SYM)
! 
!   API_NO  - No
!   API_YES - Yes

#define API_NO                         INT(0 , KIND=PASTIX_INT_KIND)
#define API_YES                        INT(1 , KIND=PASTIX_INT_KIND)
! !   Enum: API_TRACEFMT
! 
!   Trace modes (index IPARM_TRACEFMT)
! 
!   API_TRACE_PICL       - Use PICL trace format
!   API_TRACE_PAJE       - Use Paje trace format
!   API_TRACE_HUMREAD    - Use human-readable text trace format
!   API_TRACE_UNFORMATED - Unformated trace format

#define API_TRACE_PICL                 INT(0 , KIND=PASTIX_INT_KIND)
#define API_TRACE_PAJE                 INT(1 , KIND=PASTIX_INT_KIND)
#define API_TRACE_HUMREAD              INT(2 , KIND=PASTIX_INT_KIND)
#define API_TRACE_UNFORMATED           INT(3 , KIND=PASTIX_INT_KIND)
! !   Enum: API_ORDER
! 
!   Ordering modes (index IPARM_ORDERING)
! 
!   API_ORDER_SCOTCH   - Use \scotch{} ordering
!   API_ORDER_METIS    - Use \metis{} ordering
!   API_ORDER_PERSONAL - Apply user's permutation
!   API_ORDER_LOAD     - Load ordering from disk

#define API_ORDER_SCOTCH               INT(0 , KIND=PASTIX_INT_KIND)
#define API_ORDER_METIS                INT(1 , KIND=PASTIX_INT_KIND)
#define API_ORDER_PERSONAL             INT(2 , KIND=PASTIX_INT_KIND)
#define API_ORDER_LOAD                 INT(3 , KIND=PASTIX_INT_KIND)
! !   Enum: API_FLOAT
! 
!   Ordering modes (index IPARM_ORDERING)
! 
!   API_REALSINGLE    - Use \scotch{} ordering
!   API_REALDOUBLE    - Use \metis{} ordering
!   API_COMPLEXSINGLE - Apply user's permutation
!   API_COMPLEXDOUBLE - Load ordering from disk

#define API_REALSINGLE                 INT(0 , KIND=PASTIX_INT_KIND)
#define API_REALDOUBLE                 INT(1 , KIND=PASTIX_INT_KIND)
#define API_COMPLEXSINGLE              INT(2 , KIND=PASTIX_INT_KIND)
#define API_COMPLEXDOUBLE              INT(3 , KIND=PASTIX_INT_KIND)
! !  * Enum: API_GPU_CRITERIUM
!  *
!  * Criterium used to decide to put tasks on GPUs.
!  *
!  * API_GPU_CRITERION_UPDATES  - Number of updates on the panel.
!  * API_GPU_CRITERION_CBLKSIZE - Size of the target panel.
!  * API_GPU_CRITERION_FLOPS    - Number of FLOP involved in updates.
!  * API_GPU_CRITERION_PRIORITY - Priority computed in static scheduler.
#define API_GPU_CRITERION_UPDATES      INT(0 , KIND=PASTIX_INT_KIND)
#define API_GPU_CRITERION_CBLKSIZE     INT(1 , KIND=PASTIX_INT_KIND)
#define API_GPU_CRITERION_FLOPS        INT(2 , KIND=PASTIX_INT_KIND)
#define API_GPU_CRITERION_PRIORITY     INT(3 , KIND=PASTIX_INT_KIND)
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
#define MOD_UNKNOWN                    INT(0 , KIND=PASTIX_INT_KIND)
#define MOD_SOPALIN                    INT(100 , KIND=PASTIX_INT_KIND)
#define MOD_BLEND                      INT(200 , KIND=PASTIX_INT_KIND)
#define MOD_SCOTCH                     INT(300 , KIND=PASTIX_INT_KIND)
#define MOD_FAX                        INT(400 , KIND=PASTIX_INT_KIND)
#define MOD_ORDER                      INT(500 , KIND=PASTIX_INT_KIND)
#define MOD_COMMON                     INT(600 , KIND=PASTIX_INT_KIND)
#define MOD_SI                         INT(700 , KIND=PASTIX_INT_KIND)
#define MOD_GRAPH                      INT(800 , KIND=PASTIX_INT_KIND)
#define MOD_SYMBOL                     INT(900 , KIND=PASTIX_INT_KIND)
#define MOD_KASS                       INT(1000 , KIND=PASTIX_INT_KIND)
#define MOD_BUBBLE                     INT(1100 , KIND=PASTIX_INT_KIND)
#define MOD_MURGE                      INT(1200 , KIND=PASTIX_INT_KIND)
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
#define NO_ERR                         INT(0 , KIND=PASTIX_INT_KIND)
#define UNKNOWN_ERR                    INT(1 , KIND=PASTIX_INT_KIND)
#define ALLOC_ERR                      INT(2 , KIND=PASTIX_INT_KIND)
#define ASSERT_ERR                     INT(3 , KIND=PASTIX_INT_KIND)
#define NOTIMPLEMENTED_ERR             INT(4 , KIND=PASTIX_INT_KIND)
#define OUTOFMEMORY_ERR                INT(5 , KIND=PASTIX_INT_KIND)
#define THREAD_ERR                     INT(6 , KIND=PASTIX_INT_KIND)
#define INTERNAL_ERR                   INT(7 , KIND=PASTIX_INT_KIND)
#define BADPARAMETER_ERR               INT(8 , KIND=PASTIX_INT_KIND)
#define FILE_ERR                       INT(9 , KIND=PASTIX_INT_KIND)
#define BAD_DEFINE_ERR                 INT(10 , KIND=PASTIX_INT_KIND)
#define INTEGER_TYPE_ERR               INT(11 , KIND=PASTIX_INT_KIND)
#define IO_ERR                         INT(12 , KIND=PASTIX_INT_KIND)
#define MATRIX_ERR                     INT(13 , KIND=PASTIX_INT_KIND)
#define FLOAT_TYPE_ERR                 INT(14 , KIND=PASTIX_INT_KIND)
#define STEP_ORDER_ERR                 INT(15 , KIND=PASTIX_INT_KIND)
#define MPI_ERR                        INT(16 , KIND=PASTIX_INT_KIND)

