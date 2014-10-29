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
  File: sopalin3d.h

  Contains structures used in sopalin and declarations of functions of <sopalin3d.h>
 */

#ifndef SOPALIN_3D_H
#define SOPALIN_3D_H

#ifndef SOPALIN_THREAD_H
#error "sopalin_thread.h must be included before sopalin3d.h"
#endif

#ifndef SOPALIN_DEFINE_H
#error "sopalin_define.h must be included before sopalin3d.h"
#endif

typedef struct Sopalin_Data_ Sopalin_Data_t;

#ifdef WITH_STARPU
#include "starpu_submit_tasks.h"
#endif

#ifdef OOC
typedef struct ooc_thread_ ooc_thread_t;
typedef struct ooc_ ooc_t;
#endif

/*
  Enum: SOPALIN_TASK

  tasks wich will be computed in current PaStiX call

  constants:

  SOPALIN_ONLY       - Only runs factorisation.
  SOPALIN_UPDO       - Runs factorisation and up down.
  SOPALIN_UPDO_GMRES - Runs factorisation, up down and GMRES.
  SOPALIN_UPDO_GRAD  - Runs factorisation, up down and conjugate gradient.
  SOPALIN_UPDO_PIVOT - Runs factorisation, up down and simple iterative raffinement.
  UPDO_ONLY          - Only up down.
  RAFF_GMRES         - Only GMRES.
  RAFF_GRAD          - Only conjugate gradient.
  RAFF_PIVOT         - Only simple iterative raffinement.
  SOPALIN_NBTASKS    - Number of existing tasks.
*/
enum SOPALIN_TASK {
  SOPALIN_ONLY = 0,
  SOPALIN_UPDO,
  SOPALIN_UPDO_GMRES,
  SOPALIN_UPDO_GRAD,
  SOPALIN_UPDO_PIVOT,
  UPDO_ONLY,
  RAFF_GMRES,
  RAFF_GRAD,
  RAFF_PIVOT,
  SOPALIN_NBTASKS
};


/************************************************/
/*         Parametres de sopalin                */
/************************************************/
/*
  struct: SopalinParam_

  Parameters for factorisation, updown and reffinement.
 */
typedef struct SopalinParam_ {
  CscMatrix *cscmtx;          /*+ Compress Sparse Column matrix                    +*/
  double     epsilonraff;     /*+ epsilon to stop reffinement                      +*/
  double     rberror;         /*+ ||r||/||b||                                      +*/
  double     espilondiag;     /*+ epsilon critere for diag control                 +*/
  PASTIX_FLOAT     *b;               /*+ b vector (RHS and solution)                      +*/
  PASTIX_FLOAT     *transcsc;        /*+ transpose csc                                    +*/
  PASTIX_INT        itermax;         /*+ max number of iteration                          +*/
  PASTIX_INT        diagchange;      /*+ number of change of diag                         +*/
  PASTIX_INT        gmresim;         /*+ Krylov subspace size for GMRES                   +*/
  PASTIX_INT        fakefact;        /*+ Flag indicating if we want fake factorisation    +*/
  PASTIX_INT        usenocsc;        /*+ Flag indicating if we want to use the intern CSC +*/
  int        factotype;       /*+ Type of factorization                            +*/
  int        symmetric;       /*+ Symmetric                                        +*/
  MPI_Comm   pastix_comm;     /*+ MPI communicator                                 +*/
  int        type_comm;       /*+ Communication mode                               +*/
  int        nbthrdcomm;      /*+ Communication's thread number                    +*/
  PASTIX_INT       *iparm;           /*+ In/Out integer parameters                        +*/
  double    *dparm;           /*+ In/Out float parameters                          +*/
  int       *bindtab;         /*+ Define where to bin threads                      +*/
  int        stopthrd;        /*+ Boolean for communication thread controlling     +*/
  int        schur;           /*+ If API_YES won't compute last diag               +*/
  PASTIX_INT        n;               /*+ size of the matrix                               +*/
  PASTIX_INT        gN;
} SopalinParam;

/************************************************/
/*  Structure pour le AllReduce Funneled        */
/************************************************/
/*
  Struct: Pastix_Allreduce_

  Structure used for MPI_Allreduce operations in Thread Funneled mode.
*/
typedef struct Pastix_Allreduce_ {
  void *           sendbuf;        /*+ Sending buffer                  +*/
  void *           recvbuf;        /*+ receiving buffer                +*/
  int              count;          /*+ Number of elements to reduce    +*/
  MPI_Datatype     datatype;       /*+ MPI Datatype to reduce          +*/
  MPI_Op           op;             /*+ MPI operation                   +*/
} Pastix_Allreduce_t;


/************************************************/
/*         Thread Data                          */
/*   Données locales à chaque threads           */
/*   nécessitant pas de mutex pour modification */
/************************************************/
/*
   Struct: Thread_Data_

   Structure used to contain data local to each thread.
   These datas do not need to be protected by mutexes.
 */
typedef struct Thread_Data_ {
  Clock            sop_clk;                  /*+ Clock                                               +*/
  Clock            sop_clk_comm;             /*+ Communication clock                                 +*/
  PASTIX_INT              nbpivot;                  /*+ Number of pivoting performed                        +*/
  PASTIX_INT              flag_bind;                /*+ Indicate if threads are binded on processors        +*/
#ifdef TRYLOCK
  PASTIX_INT              ptbusy;                   /*+ Number of mutexes in use                            +*/
  PASTIX_INT              ptfree;                   /*+ Nomber of free mutexes                              +*/
  PASTIX_INT              ptwait;                   /*+ Nombre de cond_wait appele                          +*/
#endif
  PASTIX_INT              firstbloktab;             /*+ First block index to compute stride in maxbloktab   +*/
  PASTIX_INT              stridebloktab;            /*+ Temporary tabulars maxblokttab's stride copy        +*/
  PASTIX_FLOAT           *maxbloktab1;              /*+ Temporary tabular to add contributions              +*/
  PASTIX_FLOAT           *maxbloktab2;              /*+ Temporary tabular to add contributions (for LU)     +*/
  MPI_Request     *send_block_requests;      /*+ sent blocks requests                                +*/
  PASTIX_INT             *send_block_target;        /*+ sent blocks targets                                 +*/
  MPI_Request     *send_fanin_requests;      /*+ sent fanins requests                                +*/
#if (!defined NO_MPI_TYPE)
  MPI_Datatype    *send_fanin_mpitypes;      /*+ sent fanins mpi types                               +*/
  PASTIX_INT            **send_fanin_infotab;       /*+ sent fanins mpi types                               +*/
#endif /* not NO_MPI_TYPE */
  PASTIX_INT             *send_fanin_target;        /*+ sent fanins targets                                 +*/
  PASTIX_INT            **send_fanin_target_extra;  /*+ extra sent fanin targets                            +*/
#ifdef TEST_IRECV
  MPI_Request     *recv_fanin_request;       /*+ receiving fanin requests                            +*/
  MPI_Request     *recv_block_request;       /*+ receiving blocks requests                           +*/
  void           **recv_fanin_buffer;        /*+ fanin reception buffers                             +*/
  void           **recv_block_buffer;        /*+ blocks reception buffers                            +*/
#endif
  PASTIX_INT              maxsrequest_fanin;        /*+ Maximum number of send requests used                +*/
  PASTIX_INT              maxsrequest_block;        /*+ Maximum number of send requests used                +*/
#ifndef FORCE_NOMPI
  MPI_Status      *srteststatus;
  int             *srtestindices;
#endif
  void            *recv_buffer;              /*+ reception buffer                                    +*/
  int             *gtabsize;                 /*+ size of the data type to send                       +*/
#ifndef NO_MPI_TYPE
  MPI_Aint        *gtaboffs;                 /*+ offsize of the data type to send                    +*/
  MPI_Datatype    *gtabtype;                 /*+ Type of data to send                                +*/
#else
  void           **gtaboffs;                 /*+ pointer to the data to send                         +*/
  int             *gtabtype;                 /*+ Size of the data to send                            +*/
  void           **send_fanin_buffer;        /*+ Fanins sending buffers                              +*/
  void           **send_block_buffer;        /*+ Blocks sending buffers                              +*/
  PASTIX_INT             *send_fanin_buffer_size;   /*+ Fanins sending buffers size                         +*/
  PASTIX_INT             *send_block_buffer_size;   /*+ Blocks sending buffers size                         +*/
#endif /* NO_MPI_TYPE */
#ifdef TRACE_SOPALIN
  FILE            *tracefile;                /*+ Tracing file for the solver                         +*/
  int              traceactive;              /*+ Flag indicating if trace is active                  +*/
  int              traceid;                  /*+ Flag indicating if trace is active                  +*/
#endif
#ifdef COMPACT_SMX
  PASTIX_INT             stridebloktab2;
#endif /* COMPACT_SMX */
#ifdef PASTIX_DYNSCHED
  PASTIX_INT             esp;
#ifndef PASTIX_DYNSCHED_WITH_TREE
  PASTIX_INT            *tabtravel;
#endif
#endif
} Thread_Data_t;

/************************************************/
/*          Sopalin Data                        */
/*   Données Communes à tous les threads        */
/************************************************/
/*
   Struct: Sopalin_Data_

   Data common to all threads.
*/
struct Sopalin_Data_ {
  SolverMatrix    *datacode;                 /*+ Solver matrix                      +*/
  SopalinParam    *sopar;                    /*+ Sopalin parameters                 +*/
  Thread_Data_t  **thread_data;              /*+ Threads data                       +*/
  Queue           *fanintgtsendqueue;        /*+ Fanins to send queue               +*/
  Queue           *blocktgtsendqueue;        /*+ Blocks to send queue               +*/
  PASTIX_INT             *taskmark;                 /*+ Task marking for 2D                +*/
#ifdef TRACE_SOPALIN
  FILE            *tracefile;                /*+ Tracing file for the solver        +*/
  double           timestamp;                /*+ Original time for trace            +*/
#endif
#if (defined COMPUTE_ALLOC) || (defined STATS_SOPALIN)
  PASTIX_INT              current_alloc;            /*+ Current allocated memory           +*/
#endif
#ifdef ALLOC_FTGT
  PASTIX_INT              max_alloc;                /*+ Maximum allocated memory           +*/
  PASTIX_INT              alloc_init;               /*+ Initial allocated memory           +*/
#ifdef STATS_SOPALIN
  pthread_mutex_t  mutex_alloc;              /*+ Mutex on allocated memory variable +*/
#endif
#endif
#ifdef USE_CSC
  double volatile  critere;                  /*+ Stopping threshold for refinement   +*/
  double volatile  stop;                     /*+ Flag to stop threads on refinement  +*/
  double volatile  berr;                     /*+ Error in refinement                 +*/
  double volatile  lberr;                    /*+ Last error in refinement            +*/
  PASTIX_INT    volatile  raffnbr;                  /*+ Number of iterations in refinement  +*/
  PASTIX_INT    volatile  count_iter;               /*+ Number of iterations in refinement  +*/
  PASTIX_INT    volatile  flag_gmres;               /*+ Flag to continue in static pivoting +*/
#endif
#ifdef SMP_SOPALIN
  pthread_mutex_t *mutex_task;               /*+ Mutex on each task                               +*/
  pthread_cond_t  *cond_task;                /*+ Cond for each task                               +*/
  pthread_mutex_t *mutex_fanin;              /*+ Mutex on each fanin                              +*/
  pthread_cond_t  *cond_fanin;               /*+ Cond for each fanin                              +*/
  pthread_mutex_t *mutex_blok;               /*+ Mutex on each block                              +*/
  pthread_mutex_t *mutex_queue_fanin;        /*+ Mutex on the fanins queue                        +*/
  pthread_mutex_t *mutex_queue_block;        /*+ Mutex on the blocks queue                        +*/
#else /* SMP_SOPALIN */
  Queue            taskqueue;                /*+ Task queue for NOSMP version                     +*/
#endif
  sopthread_barrier_t barrier;               /*+ Threads synchronisation barrier                  +*/
  pthread_mutex_t  mutex_comm;               /*+ Mutex on communication variables                 +*/
  pthread_cond_t   cond_comm;                /*+ Condition on step_comm                           +*/
  int              step_comm;                /*+ Current step indicator                           +*/
  Pastix_Allreduce_t allreduce;              /*+ Data structure for MPi_Allreduce                 +*/
  Queue             *sendqueue;              /*+ Ready to send data queue                         +*/
#ifdef STORAGE
  PASTIX_FLOAT           *grhs;                     /*+ Data storage tabular                             +*/
  volatile PASTIX_INT    *flagtab;                  /*+ Indicate received cblk in up step                +*/
  pthread_mutex_t *mutex_flagtab;            /*+ Mutex on flagtab                                 +*/
  pthread_cond_t  *cond_flagtab;             /*+ cond on flagtab                                  +*/
#endif
  volatile void   *ptr_raff[15];             /*+ pointers used in refinement                      +*/
  void           **ptr_csc;                  /*+ pointer to data used by each threads in csc_code +*/
  volatile PASTIX_FLOAT    *common_flt;      /*+ Common pointer to share a float between threads  +*/
  volatile double          *common_dbl;      /*+ Common pointer to share a double between threads  +*/
  pthread_mutex_t  mutex_raff;               /*+ mutex on common tabulars during csc operations   +*/
  pthread_cond_t   cond_raff;                /*+ cond corresponding to mutex_raff                 +*/
#ifdef PASTIX_DYNSCHED
  pthread_mutex_t *tasktab_mutex;            /*+ +*/
  pthread_cond_t  *tasktab_cond;             /*+ +*/
  volatile PASTIX_INT    *tasktab_indice;           /*+ +*/
  volatile PASTIX_INT    *tasktab_nbthrd;           /*+ +*/
  Queue           *taskqueue;                /*+ +*/
#endif
#ifdef OOC
  ooc_t           *ooc;                      /*+ Data structure needed for Out-of-core            +*/
#endif
#ifndef WITH_HWLOC
#  ifdef PASTIX_GET_SCHED_AFFINITY
  int             *allowed_cpus;             /*+ List of authorized CPUs for binding +*/
  int             ncore_avail;               /*+ number of cores available +*/
#  endif
#endif
#ifdef WITH_STARPU
  starpu_loop_data_t  *starpu_loop_data;
#endif
};

#ifdef WITH_STARPU
struct starpu_trf_data_ {
  PASTIX_INT              cblknum;
#ifdef STARPU_SUBMIT_READY
  PASTIX_INT              tasknum;
#endif
  Sopalin_Data_t * sopalin_data;
};
typedef struct starpu_trf_data_ starpu_trf_data_t;

struct starpu_gemm_data_ {
  PASTIX_INT              cblknum;
#ifdef STARPU_SUBMIT_READY
  PASTIX_INT              tasknum;
#endif
  PASTIX_INT              bloknum;
  PASTIX_INT              fcblknum;
  PASTIX_INT              nblocs;
  Sopalin_Data_t * sopalin_data;
  int           ** d_blocktab;
};
typedef struct starpu_gemm_data_ starpu_gemm_data_t;
#endif
/************************************************/
/*     Fonctions publiques de sopalin3d         */
/************************************************/
#define he_sopalin_thread PASTIX_PREFIX_F(he_sopalin_thread)
#define sy_sopalin_thread PASTIX_PREFIX_F(sy_sopalin_thread)
#define ge_sopalin_thread PASTIX_PREFIX_F(ge_sopalin_thread)

#define po_sopalin_updo_thread PASTIX_PREFIX_F(po_sopalin_updo_thread)
#define he_sopalin_updo_thread PASTIX_PREFIX_F(he_sopalin_updo_thread)
#define sy_sopalin_updo_thread PASTIX_PREFIX_F(sy_sopalin_updo_thread)
#define ge_sopalin_updo_thread PASTIX_PREFIX_F(ge_sopalin_updo_thread)

#define po_sopalin_updo_gmres_thread PASTIX_PREFIX_F(po_sopalin_updo_gmres_thread)
#define he_sopalin_updo_gmres_thread PASTIX_PREFIX_F(he_sopalin_updo_gmres_thread)
#define sy_sopalin_updo_gmres_thread PASTIX_PREFIX_F(sy_sopalin_updo_gmres_thread)
#define ge_sopalin_updo_gmres_thread PASTIX_PREFIX_F(ge_sopalin_updo_gmres_thread)

#define po_sopalin_updo_grad_thread	PASTIX_PREFIX_F(po_sopalin_updo_grad_thread)
#define he_sopalin_updo_grad_thread	PASTIX_PREFIX_F(he_sopalin_updo_grad_thread)
#define sy_sopalin_updo_grad_thread	PASTIX_PREFIX_F(sy_sopalin_updo_grad_thread)
#define ge_sopalin_updo_pivot_thread	PASTIX_PREFIX_F(ge_sopalin_updo_pivot_thread)

#define po_updo_thread PASTIX_PREFIX_F(po_updo_thread)
#define he_updo_thread PASTIX_PREFIX_F(he_updo_thread)
#define sy_updo_thread PASTIX_PREFIX_F(sy_updo_thread)
#define ge_updo_thread PASTIX_PREFIX_F(ge_updo_thread)

#define po_gmres_thread PASTIX_PREFIX_F(po_gmres_thread)
#define he_gmres_thread PASTIX_PREFIX_F(he_gmres_thread)
#define sy_gmres_thread PASTIX_PREFIX_F(sy_gmres_thread)
#define ge_gmres_thread PASTIX_PREFIX_F(ge_gmres_thread)

#define po_grad_thread	PASTIX_PREFIX_F(po_grad_thread)
#define he_grad_thread	PASTIX_PREFIX_F(he_grad_thread)
#define sy_grad_thread	PASTIX_PREFIX_F(sy_grad_thread)
#define ge_pivot_thread PASTIX_PREFIX_F(ge_pivot_thread)

/*
  Functions: <Sopalin3d.c> functions declarations.
 */
#define ge_sopalin_launch                 PASTIX_PREFIX_F(ge_sopalin_launch)
#define ge_sopalin_thread                 PASTIX_PREFIX_F(ge_sopalin_thread)
#define ge_sopalin_updo_thread            PASTIX_PREFIX_F(ge_sopalin_updo_thread)
#define ge_sopalin_updo_gmres_thread      PASTIX_PREFIX_F(ge_sopalin_updo_gmres_thread)
#define ge_sopalin_updo_grad_thread       PASTIX_PREFIX_F(ge_sopalin_updo_grad_thread)
#define ge_sopalin_updo_pivot_thread      PASTIX_PREFIX_F(ge_sopalin_updo_pivot_thread)
#define ge_sopalin_updo_bicgstab_thread    PASTIX_PREFIX_F(ge_sopalin_updo_bicgstab_thread)
#define ge_updo_thread                    PASTIX_PREFIX_F(ge_updo_thread)
#define ge_pivot_thread                   PASTIX_PREFIX_F(ge_pivot_thread)
#define ge_gmres_thread                   PASTIX_PREFIX_F(ge_gmres_thread)
#define ge_grad_thread                    PASTIX_PREFIX_F(ge_grad_thread)
#define ge_bicgstab_thread                 PASTIX_PREFIX_F(ge_bicgstab_thread)
void    ge_sopalin_launch                (SolverMatrix *, SopalinParam *, PASTIX_INT cas);
void    ge_sopalin_thread                (SolverMatrix *, SopalinParam *);
void    ge_sopalin_updo_thread           (SolverMatrix *, SopalinParam *);
void    ge_sopalin_updo_gmres_thread     (SolverMatrix *, SopalinParam *);
void    ge_sopalin_updo_grad_thread      (SolverMatrix *, SopalinParam *);
void    ge_sopalin_updo_pivot_thread     (SolverMatrix *, SopalinParam *);
void    ge_sopalin_updo_bicgstab_thread   (SolverMatrix *, SopalinParam *);
void    ge_updo_thread                   (SolverMatrix *, SopalinParam *);
void    ge_pivot_thread                  (SolverMatrix *, SopalinParam *);
void    ge_gmres_thread                  (SolverMatrix *, SopalinParam *);
void    ge_grad_thread                   (SolverMatrix *, SopalinParam *);
void    ge_bicgstab_thread                (SolverMatrix *, SopalinParam *);

#define sy_sopalin_launch                 PASTIX_PREFIX_F(sy_sopalin_launch)
#define sy_sopalin_thread                 PASTIX_PREFIX_F(sy_sopalin_thread)
#define sy_sopalin_updo_thread            PASTIX_PREFIX_F(sy_sopalin_updo_thread)
#define sy_sopalin_updo_gmres_thread      PASTIX_PREFIX_F(sy_sopalin_updo_gmres_thread)
#define sy_sopalin_updo_grad_thread       PASTIX_PREFIX_F(sy_sopalin_updo_grad_thread)
#define sy_sopalin_updo_pivot_thread      PASTIX_PREFIX_F(sy_sopalin_updo_pivot_thread)
#define sy_sopalin_updo_bicgstab_thread    PASTIX_PREFIX_F(sy_sopalin_updo_bicgstab_thread)
#define sy_updo_thread                    PASTIX_PREFIX_F(sy_updo_thread)
#define sy_pivot_thread                   PASTIX_PREFIX_F(sy_pivot_thread)
#define sy_gmres_thread                   PASTIX_PREFIX_F(sy_gmres_thread)
#define sy_grad_thread                    PASTIX_PREFIX_F(sy_grad_thread)
#define sy_bicgstab_thread                 PASTIX_PREFIX_F(sy_bicgstab_thread)
void    sy_sopalin_launch                (SolverMatrix *, SopalinParam *, PASTIX_INT cas);
void    sy_sopalin_thread                (SolverMatrix *, SopalinParam *);
void    sy_sopalin_updo_thread           (SolverMatrix *, SopalinParam *);
void    sy_sopalin_updo_gmres_thread     (SolverMatrix *, SopalinParam *);
void    sy_sopalin_updo_grad_thread      (SolverMatrix *, SopalinParam *);
void    sy_sopalin_updo_pivot_thread     (SolverMatrix *, SopalinParam *);
void    sy_sopalin_updo_bicgstab_thread   (SolverMatrix *, SopalinParam *);
void    sy_updo_thread                   (SolverMatrix *, SopalinParam *);
void    sy_pivot_thread                  (SolverMatrix *, SopalinParam *);
void    sy_gmres_thread                  (SolverMatrix *, SopalinParam *);
void    sy_grad_thread                   (SolverMatrix *, SopalinParam *);
void    sy_bicgstab_thread                (SolverMatrix *, SopalinParam *);

#define he_sopalin_launch                 PASTIX_PREFIX_F(he_sopalin_launch)
#define he_sopalin_thread                 PASTIX_PREFIX_F(he_sopalin_thread)
#define he_sopalin_updo_thread            PASTIX_PREFIX_F(he_sopalin_updo_thread)
#define he_sopalin_updo_gmres_thread      PASTIX_PREFIX_F(he_sopalin_updo_gmres_thread)
#define he_sopalin_updo_grad_thread       PASTIX_PREFIX_F(he_sopalin_updo_grad_thread)
#define he_sopalin_updo_pivot_thread      PASTIX_PREFIX_F(he_sopalin_updo_pivot_thread)
#define he_sopalin_updo_bicgstab_thread    PASTIX_PREFIX_F(he_sopalin_updo_bicgstab_thread)
#define he_updo_thread                    PASTIX_PREFIX_F(he_updo_thread)
#define he_pivot_thread                   PASTIX_PREFIX_F(he_pivot_thread)
#define he_gmres_thread                   PASTIX_PREFIX_F(he_gmres_thread)
#define he_grad_thread                    PASTIX_PREFIX_F(he_grad_thread)
#define he_bicgstab_thread                 PASTIX_PREFIX_F(he_bicgstab_thread)
void    he_sopalin_launch                (SolverMatrix *, SopalinParam *, PASTIX_INT cas);
void    he_sopalin_thread                (SolverMatrix *, SopalinParam *);
void    he_sopalin_updo_thread           (SolverMatrix *, SopalinParam *);
void    he_sopalin_updo_gmres_thread     (SolverMatrix *, SopalinParam *);
void    he_sopalin_updo_grad_thread      (SolverMatrix *, SopalinParam *);
void    he_sopalin_updo_pivot_thread     (SolverMatrix *, SopalinParam *);
void    he_sopalin_updo_bicgstab_thread   (SolverMatrix *, SopalinParam *);
void    he_updo_thread                   (SolverMatrix *, SopalinParam *);
void    he_pivot_thread                  (SolverMatrix *, SopalinParam *);
void    he_gmres_thread                  (SolverMatrix *, SopalinParam *);
void    he_grad_thread                   (SolverMatrix *, SopalinParam *);
void    he_bicgstab_thread                (SolverMatrix *, SopalinParam *);

#define po_sopalin_launch                 PASTIX_PREFIX_F(po_sopalin_launch)
#define po_sopalin_thread                 PASTIX_PREFIX_F(po_sopalin_thread)
#define po_sopalin_updo_thread            PASTIX_PREFIX_F(po_sopalin_updo_thread)
#define po_sopalin_updo_gmres_thread      PASTIX_PREFIX_F(po_sopalin_updo_gmres_thread)
#define po_sopalin_updo_grad_thread       PASTIX_PREFIX_F(po_sopalin_updo_grad_thread)
#define po_sopalin_updo_pivot_thread      PASTIX_PREFIX_F(po_sopalin_updo_pivot_thread)
#define po_sopalin_updo_bicgstab_thread    PASTIX_PREFIX_F(po_sopalin_updo_bicgstab_thread)
#define po_updo_thread                    PASTIX_PREFIX_F(po_updo_thread)
#define po_pivot_thread                   PASTIX_PREFIX_F(po_pivot_thread)
#define po_gmres_thread                   PASTIX_PREFIX_F(po_gmres_thread)
#define po_grad_thread                    PASTIX_PREFIX_F(po_grad_thread)
#define po_bicgstab_thread                 PASTIX_PREFIX_F(po_bicgstab_thread)
void    po_sopalin_launch                (SolverMatrix *, SopalinParam *, PASTIX_INT cas);
void    po_sopalin_thread                (SolverMatrix *, SopalinParam *);
void    po_sopalin_updo_thread           (SolverMatrix *, SopalinParam *);
void    po_sopalin_updo_gmres_thread     (SolverMatrix *, SopalinParam *);
void    po_sopalin_updo_grad_thread      (SolverMatrix *, SopalinParam *);
void    po_sopalin_updo_pivot_thread     (SolverMatrix *, SopalinParam *);
void    po_sopalin_updo_bicgstab_thread   (SolverMatrix *, SopalinParam *);
void    po_updo_thread                   (SolverMatrix *, SopalinParam *);
void    po_pivot_thread                  (SolverMatrix *, SopalinParam *);
void    po_gmres_thread                  (SolverMatrix *, SopalinParam *);
void    po_grad_thread                   (SolverMatrix *, SopalinParam *);
void    po_bicgstab_thread                (SolverMatrix *, SopalinParam *);

#endif
