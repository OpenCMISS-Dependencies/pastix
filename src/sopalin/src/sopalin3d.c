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
  File: sopalin3d.c

  sopalin 3d main program.

  Authors:
        Mathieu Faverge - faverge@labri.fr
        Xavier  Lacoste - lacoste@labri.fr
        Pierre  Ramet   - ramet@labri.fr

  Date:
        Version 0.0 - february 2003
*/

/***********************************************/
/*                HEADERS                      */
/***********************************************/

/*
 * Redefinition du test pour madeleine, car les
 * requetes madeleine sont des structures
 */

#ifndef MAD_MPI
#define MPI_Request_is_equal(r1, r2) ((r1) == (r2))
#endif

/* Utilisation des blas IBM sur power */

#if (defined X_ARCHpower_ibm_aix)
#if (!(defined X_INCLUDE_ESSL))
#define X_INCLUDE_ESSL
#endif
#endif

/* Include pour compaq */
#if (defined X_ARCHalpha_compaq_osf1)
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/processor.h>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#define X_INCLUDE_CXML
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>

#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#ifdef X_INCLUDE_ESSL
#include <essl.h>
/*#include <pessl.h>*/
#endif /* X_INCLUDE_ESSL */

#include <signal.h>
#include "common_pastix.h"
#include "tools.h"
#ifdef PASTIX_EZTRACE
#include "pastix_eztrace.h"
#else
#include "trace.h"
#endif
#include "sopalin_define.h"
#include "symbol.h"
#include "ftgt.h"
#include "csc.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
#include "sopalin_thread.h"
#include "stack.h"
#include "sopalin3d.h"
#include "sopalin_init.h"
#include "perf.h"
#include "out.h"
#include "coefinit.h"
#include "ooc.h"
#include "order.h"
#include "debug_dump.h"
#ifdef WITH_STARPU
#include "starpu_submit_tasks.h"
#endif

/**************************************/
/*          A Lire                    */
/**************************************/

/* ATTENTION cast alignement des buffers de reception (facto+updo) complexe */
/* SMP : mutex_queue necessaire ??? (si les files sont par thread) */
/* SMP : PB_SMP commenter SYMB_CBLKNUM(-t)=... 2D sans SOLVER_UPDOWN ??? */
/* SMP : PB_SMP free blocktarget ??? */
/* EXACT_TAG risque debordement SEPFB -> EXACT_THREAD pb dans updown ??? */


/**************************************/
/*           Options de compil        */
/**************************************/


#ifdef HPM_SOPALIN
#include <libhpm.h>
#endif


/* Définition d'un décalage pour les tags fanin et block */
#define SEPFB MAXTHRDS
#if (defined EXACT_TAG)
#undef  SEPFB
#define SEPFB 40000
#endif


/*static PASTIX_INT iterator;*/
#define DO_ITER_MAX 10
/*#define DO_ITER(x) {for(iterator=0;iterator<DO_ITER_MAX;iterator++){(x);}}*/
#define DO_ITER(x) {if (SOLV_PROCNBR > 1) {(x);};}


#include "sopalin_time.h"

#ifdef COMPUTE

#include "sopalin_compute.h"

/*
   Section: Global variables
*/
/*
   int: iun
   Integer 1
*/
static PASTIX_INT   iun   = 1;
/* static PASTIX_INT izero=0; */
/*
  float: fun
  Floating point   1.0
*/
#ifdef CPLX
static PASTIX_FLOAT fun   = 1.0+0.0*I;
#else
static PASTIX_FLOAT fun   = 1.0;
#endif
/*
  Float: fzero
  Floating point   0.0
*/
static PASTIX_FLOAT fzero = 0.0;


#else /* COMPUTE */

/* BLAS computations are disabled */
#define SOPALIN_GEAM
#define SOPALIN_GESM
#define SOPALIN_COPY
#define SOPALIN_SCAL
#define SOPALIN_PPF
#define SOPALIN_POF
#define SOPALIN_TRSM
#define SOPALIN_TRSV
#define SOPALIN_GEMM
#define SOPALIN_GEMV
#define SOPALIN_AXPY
#define SOPALIN_GER
#define SOPALIN_SYR

#endif /* COMPUTE */

/* Definition des debug */
#include "sopalin_debug.h"

/* Acces aux champs de datacode */
#include "sopalin_acces.h"

/***********************************/
/*      Affichage                  */
/***********************************/
/*
  Section: Macros

  Macros: Printing maccros.

  print_onempi(fmt, ...) - Print message by one MPI task.
  print_one(fmt, ...)    - Print message by one thread of one MPI task
  print_error(...)       - Print error messages (ignored).

 */
#define print_onempi(fmt, ...) if( SOLV_PROCNUM == 0 )           fprintf(stdout, fmt, __VA_ARGS__)
#define print_one(fmt, ...)    if( me == 0 && SOLV_PROCNUM == 0) fprintf(stdout, fmt, __VA_ARGS__)
#define print_error(...)


#ifdef VERIF_MPI
int err_mpi;
#endif

#define LOCAL_ALLOC_BTAG -100

/* ??? extra-diag blocks in 1D column-block (computed by blend) */
/*
#undef PACKMAX
#define PACKMAX 32
#undef PACKAREA
#define PACKAREA 200000
*/

/************************************************/
/*       Déclaration des fonctions              */
/************************************************/
/* Section: Prototypes */
#define dump_all                  API_CALL(dump_all)
#define init_struct_sopalin       API_CALL(init_struct_sopalin)
#define sopalin_launch            API_CALL(sopalin_launch)
#define sopalin_updo_comm         API_CALL(sopalin_updo_comm)
#define sopalin_thread            API_CALL(sopalin_thread)
#define sopalin_smp               API_CALL(sopalin_smp)
#define sopalin_updo_thread       API_CALL(sopalin_updo_thread)
#define sopalin_updo_smp          API_CALL(sopalin_updo_smp)
#define sopalin_updo_gmres_thread API_CALL(sopalin_updo_gmres_thread)
#define sopalin_updo_gmres_smp    API_CALL(sopalin_updo_gmres_smp)
#define sopalin_updo_grad_thread  API_CALL(sopalin_updo_grad_thread)
#define sopalin_updo_grad_smp     API_CALL(sopalin_updo_grad_smp)
#define sopalin_updo_pivot_thread API_CALL(sopalin_updo_pivot_thread)
#define sopalin_updo_pivot_smp    API_CALL(sopalin_updo_pivot_smp)
#define up_down                   API_CALL(up_down)
#define up_down_smp               API_CALL(up_down_smp)



void  dump_all                 (SolverMatrix *, CscMatrix * cscmtx, int);
void  init_struct_sopalin      (Sopalin_Data_t *sopalin_data, SolverMatrix *m,
                                SopalinParam *sopar);
void  sopalin_launch           (SolverMatrix *m, SopalinParam *sopaparam, PASTIX_INT cas);
void* sopalin_updo_comm        (void *arg);
void  sopalin_thread           (SolverMatrix *m, SopalinParam *sopaparam);
void* sopalin_smp              (void *arg);
void  sopalin_updo_thread      (SolverMatrix *m, SopalinParam *sopaparam);
void* sopalin_updo_smp         (void *arg);
void  sopalin_updo_gmres_thread(SolverMatrix *m, SopalinParam *sopaparam);
void* sopalin_updo_gmres_smp   (void *arg);
void  sopalin_updo_grad_thread (SolverMatrix *m, SopalinParam *sopaparam);
void* sopalin_updo_grad_smp    (void *arg);
void  sopalin_updo_pivot_thread(SolverMatrix *m, SopalinParam *sopaparam);
void* sopalin_updo_pivot_smp   (void *arg);
void  up_down                  (void);
void* up_down_smp              (void * arg);

/************************************************/
/*              Dump des matrices               */
/************************************************/
/* Section: Debug functions */
/*
   Function: API_CALL(dump_all)

   Dumps the matrix and right-hand-side on disk.

   This function can dump the internal distributed CSC matrix,
   or the solvermatrix, or the Up-down vector.

   This function must be called by only one thread for
   each MPI process.

   The *x* value can be defined using *DUMP_CSC*, *DUMP_SOLV* and *DUMP_SMB*,
   combined like *DUMP_CSC | DUMP_SOLV | DUMP_SMB*

   Parameters:
         datacode - SolverMatrix
         x        - value indicating what to dump.

   Returns:
         void
*/
void dump_all(SolverMatrix *datacode,
              CscMatrix    *cscmtx,
              int           x)
{
  /*
  ** Fonction appellée par un seul thread par processus MPI
  **   instance : étape où sont dumpées les matrices
  **   x        : vecteurs a dumper.
  */
        static int instance = 0;
        FILE *file;
        char  filename[250];
#ifdef SOPALIN_LU
        FILE *fileu;
        char  filenameu[250];
#endif

        printf("Dump CSC SOLV SMB (%ld)...\n",(long) instance);

        /* CSC */
#ifdef USE_CSC
        if (x & DUMP_CSC)
          {
        sprintf(filename, "csc%ld.%ld",
                (long) instance,
                (long) SOLV_PROCNUM);
        file = fopen(filename, "w");
        dump2(datacode, cscmtx, NULL, file);
        fclose(file);
          }
#endif

        /* SOLV */
        if (x & DUMP_SOLV)
          {
        sprintf(filename, "solv%ld.%ld",(long) instance,(long) SOLV_PROCNUM);
        file = fopen(filename, "w");
#ifdef SOPALIN_LU
        sprintf(filenameu, "solvu%ld.%ld",(long) instance,(long) SOLV_PROCNUM);
        fileu = fopen(filenameu, "w");
        dump3_LU(datacode, file, fileu);
        fclose(fileu);
#else
        dump3(datacode, file);
#endif
        fclose(file);
          }

        /* SMB */
        if (x & DUMP_SMB)
          {
        sprintf(filename, "smb%ld.%ld",(long) instance,(long) SOLV_PROCNUM);
        file = fopen( filename, "w");
        dump5(datacode, file);
        fclose(file);
          }
        instance++;
}

/****************************************************************************/
/*                     COMMUNICATION ROUTINES                               */
/****************************************************************************/

#include "./sopalin_sendrecv.c"

/****************************************************************************/
/*                       COMPUTE ROUTINES                                   */
/****************************************************************************/

/* Section : Computation routines prototypes */
void API_CALL(compute_diag)  (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task);
void API_CALL(compute_1d)    (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task);
void API_CALL(compute_1dgemm)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task, PASTIX_INT i, PASTIX_INT b2);
void API_CALL(compute_e1)    (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task);
void API_CALL(compute_e2)    (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task);

#include "./sopalin_compute.c"

#if defined(PASTIX_DYNSCHED)
#include "./dynsched.h"
#endif

#ifdef SOLVER_UPDOWN

/****************************************************************************/
/*                       UP AND DOWN STEP                                   */
/****************************************************************************/
/* Section : Up-down step routines prototypes */
/* Initialisation / Nettoyage */
void  API_CALL(updo_init)(Sopalin_Data_t *sopalin_data, SolverMatrix *datacode, SopalinParam *sopaparam);

/* Thread de communication */
void* API_CALL(updo_thread_comm)(void *);

/* Lancement de updo seul */
void  API_CALL(updo_thread)(SolverMatrix *datacode, SopalinParam *sopaparam);

#include "updo.c"
/* Section : Reffinement step routines prototypes */
/* Raffinement du second membre */
void* API_CALL(pivotstatique_smp)(void *arg);
void* API_CALL(gmres_smp)        (void *arg);
void* API_CALL(grad_smp)         (void *arg);
void* API_CALL(bicgstab_smp)      (void *arg);

/* Lancement d'une des fonctions seules */
void  API_CALL(pivot_thread)   (SolverMatrix *datacode, SopalinParam *sopaparam);
void  API_CALL(gmres_thread)   (SolverMatrix *datacode, SopalinParam *sopaparam);
void  API_CALL(grad_thread)    (SolverMatrix *datacode, SopalinParam *sopaparam);
void  API_CALL(bicgstab_thread) (SolverMatrix *datacode, SopalinParam *sopaparam);

#include "csc_intern_compute.h"

#define RAFF_CLOCK_INIT {clockInit(&raff_clk);clockStart(&raff_clk);}
#define RAFF_CLOCK_STOP {clockStop(&(raff_clk));}
#define RAFF_CLOCK_GET  clockVal(&(raff_clk))

/* #define DEBUG_RAFF */
#include "raff_functions.h"
#include "raff_grad.c"
#include "raff_gmres.c"
#include "raff_pivot.c"
#include "raff_bicgstab.c"
#endif

/****************************************************************************/
/*                  INITIALIZATION ROUTINES                                 */
/****************************************************************************/

/* Initialisation de la CSC */
/* Section: Sopalin functions */

/*
 * Function: init_struct_sopalin
 *
 * Initialization routine.
 *
 * Set the sopalin_data->critere, allocate FANIN_COEFTAB and TASK_BTAGPTR.
 *
 * Parameters:
 *       sopalin_data - Structure used during factorisation and resolution.
 *       datacode     - SolverMatrix structure.
 *	sopar        - Factorisation parameters.
 */
#define init_struct_sopalin API_CALL(init_struct_sopalin)
void init_struct_sopalin (Sopalin_Data_t *sopalin_data,
                          SolverMatrix   *datacode,
                          SopalinParam   *sopar)
{
  MPI_Comm pastix_comm = PASTIX_COMM;
  PASTIX_INT      i;
#if (defined DEADCODE) && !(defined ALLOC_FTGT)
  PASTIX_INT j;
#endif

  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_onempi("%s", OUT2_SOP_TABG);

  /* A la louche ??? */
#ifdef SOPALIN_LU
  SOLV_CPFTMAX *= 2;
  PACKAREA     *= 2;
#endif

#if (defined COMM_REORDER) && (defined PASTIX_NMAD_AGGREG)
  PACKMAX = 1;
#endif

  /*PACKMAX=(PACKMAX<32)?PACKMAX:32;*/
  /*PACKAREA=PACKMAX*SOLV_CPFTMAX;*/

  print_debug(DBG_SOPALIN_MAIN,"--- PACKMAX=%ld PACKAREA=%ld\n", (long) PACKMAX, (long) PACKAREA);
  print_debug(DBG_SOPALIN_MAIN,"--- ** sopalin coefnbr=%ld ftgtnbr=%ld bpftmax=%ld cpftmax=%ld coefmax=%ld cblknbr=%ld bloknbr=%ld procnum=%ld procnbr=%ld\n",
                  (long) SOLV_COEFNBR, (long) SOLV_FTGTNBR, (long) SOLV_BPFTMAX,
                  (long) SOLV_CPFTMAX, (long) SOLV_COEFMAX, (long) SYMB_CBLKNBR,
                  (long) SYMB_BLOKNBR, (long) SOLV_PROCNUM, (long) SOLV_PROCNBR);

  /* Redéfinition des erreurs pour Compaq */
#ifdef X_INCLUDE_ESSL
  {
        int ierno,inoal,inomes,itrace,iusadr,irange,dummy;

        /* Special pour le pere Goudin ... */

#define NUM_ERROR_PPF 2104
#define NUM_ERROR_POF 2115
#define NOT_ALTERED      0
#define LIMIT_ERROR    255
#define UNLIMITED_ERROR LIMIT_ERROR+1
#define IGNORE_ERROR    -1
#define IGNORE_MESSAGE  -1

        /*
        printf("les modifs du pere goudin\n");

        einfo(0,&dummy,&dummy);

        ierno=NUM_ERROR_PPF;
        inoal=UNLIMITED_ERROR;
        inomes=NOT_ALTERED;
        itrace=NOT_ALTERED;
        iusadr=NOT_ALTERED;
        irange=NUM_ERROR_PPF;

        errset(&ierno,&inoal,&inomes,&itrace,&iusadr,&irange);

        ierno=NUM_ERROR_POF;
        inoal=UNLIMITED_ERROR;
        inomes=NOT_ALTERED;
        itrace=NOT_ALTERED;
        iusadr=NOT_ALTERED;
        irange=NUM_ERROR_POF;

        errset(&ierno,&inoal,&inomes,&itrace,&iusadr,&irange);
        */
  }
#endif /* X_INCLUDE_ESSL */

  /* Statistiques d'allocation */
#ifdef ALLOC_FTGT
  {
        double factor = 0.0;
        PASTIX_INT    alloc_init =
          SYMB_CBLKNBR*3*sizeof(PASTIX_INT)+
          SYMB_BLOKNBR*3*sizeof(PASTIX_INT)+
          SYMB_CBLKNBR*1*sizeof(PASTIX_INT)+
          SYMB_BLOKNBR*1*sizeof(PASTIX_INT)+
          SOLV_TASKNBR  *sizeof(Task)+
          SOLV_FTGTNBR  *sizeof(FanInTarget)+
          SOLV_COEFNBR  *sizeof(PASTIX_FLOAT)+
          SOLV_BTAGNBR  *sizeof(BlockTarget)+
          SOLV_BCOFNBR  *sizeof(BlockCoeff)+
          SOLV_INDNBR   *sizeof(PASTIX_INT);

        factor = 100.0 / (double)alloc_init;
        (void)factor;
        printf("symbol.cblk %12ld %2.2lf %%\n",
           (long)  (SYMB_CBLKNBR*3*sizeof(PASTIX_INT)),
           (double)(SYMB_CBLKNBR*3*sizeof(PASTIX_INT))*factor);
        printf("symbol.blok %12ld %2.2lf %%\n",
           (long)  (SYMB_BLOKNBR*3*sizeof(PASTIX_INT)),
           (double)(SYMB_BLOKNBR*3*sizeof(PASTIX_INT))*factor);
        printf("solver.cblk %12ld %2.2lf %%\n",
           (long)  (SYMB_CBLKNBR*1*sizeof(PASTIX_INT)),
           (double)(SYMB_CBLKNBR*1*sizeof(PASTIX_INT))*factor);
        printf("solver.blok %12ld %2.2lf %%\n",
           (long)  (SYMB_BLOKNBR*1*sizeof(PASTIX_INT)),
           (double)(SYMB_BLOKNBR*1*sizeof(PASTIX_INT))*factor);
        printf("solver.task %12ld %2.2lf %%\n",
           (long)  (SOLV_TASKNBR*1*sizeof(Task)),
           (double)(SOLV_TASKNBR*1*sizeof(Task))*factor);
        printf("solver.ftgt %12ld %2.2lf %%\n",
           (long)  (SOLV_FTGTNBR*1*sizeof(FanInTarget)),
           (double)(SOLV_FTGTNBR*1*sizeof(FanInTarget))*factor);
        printf("solver.coef %12ld %2.2lf %%\n",
           (long)  (SOLV_COEFNBR*1*sizeof(PASTIX_FLOAT)),
           (double)(SOLV_COEFNBR*1*sizeof(PASTIX_FLOAT))*factor);
        printf("solver.btag %12ld %2.2lf %%\n",
           (long)  (SOLV_BTAGNBR*1*sizeof(BlockTarget)),
           (double)(SOLV_BTAGNBR*1*sizeof(BlockTarget))*factor);
        printf("solver.bcof %12ld %2.2lf %%\n",
           (long)  (SOLV_BCOFNBR*1*sizeof(BlockCoeff)),
           (double)(SOLV_BCOFNBR*1*sizeof(BlockCoeff))*factor);
        printf("solver.ind  %12ld %2.2lf %%\n",
           (long)  (SOLV_INDNBR *1*sizeof(PASTIX_INT)),
           (double)(SOLV_INDNBR *1*sizeof(PASTIX_INT))*factor);
  }
#endif

  /* Computing sopalin_data->critere */
  /* Computing sopalin_data->critere */
#ifdef COMPUTE
#ifdef USE_CSC

  sopalin_data->critere = sopar->espilondiag;
  if (sopalin_data->critere<0.0) {
    /* sopalin_data->critere absolu */
    sopalin_data->critere=-sopalin_data->critere;
    if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_YES) {
      fprintf(stdout, "Pivoting criterium (epsilon) = %g\n", sopalin_data->critere);
    }
  } else {
    if (sopar->usenocsc == 1) {
      sopalin_data->critere = sopar->espilondiag;
    } else {
      if (sopar->fakefact == 1) {
        sopalin_data->critere = (UPDOWN_GNODENBR*UPDOWN_GNODENBR+UPDOWN_GNODENBR)*sqrt(sopar->espilondiag);
      } else {
        sopalin_data->critere = CscNorm1(sopalin_data->sopar->cscmtx, pastix_comm)*sqrt(sopar->espilondiag);
      }
    }
    if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_YES) {
      fprintf(stdout,"Pivoting criterium (||A||*sqrt(epsilon)) = %g\n", sopalin_data->critere);
    }
  }

  /* Allocating FANIN_COEFTAB */
#endif /* USE_CSC */
#endif

  for (i=0;i<SOLV_FTGTNBR;i++)
    {
      PASTIX_INT ftgtsize;
      (void)ftgtsize;

#ifdef DEBUG_SOPALIN_INIT
      printf("ftgt %ld : %ld %ld %ld %ld\n",i,FANIN_LROWNUM(i),FANIN_FROWNUM(i),FANIN_LCOLNUM(i),FANIN_FCOLNUM(i));
#endif

      ftgtsize = (FANIN_LROWNUM(i)-FANIN_FROWNUM(i)+1)
        *(FANIN_LCOLNUM(i)-FANIN_FCOLNUM(i)+1);

#ifdef SOPALIN_LU
      ftgtsize *= 2;
#endif
#if (defined ALLOC_FTGT )
#if !(defined OOC_FTGT)
      FANIN_COEFTAB(i) = NULL;
#endif
#else
      MALLOC_INTERN(FANIN_COEFTAB(i), ftgtsize, PASTIX_FLOAT);
      for (j=0;j<ftgtsize;j++)
        FANIN_COEFTAB(i)[j] = 0.0;
#endif
    }
#ifdef DEBUG_SOPALIN_INIT
  printf("end init ftgttab\n");
#endif

  /* ???
         SYMB_NODENBR = UPDOWN_GNODENBR;
         fprintf(stdout, "Node NBR : %ld\n",SYMB_NODENBR);
  */

}

#include "./contrib.c"

/****************************************************************************/
/*                    FACTORIZATION ROUTINE                                 */
/****************************************************************************/
/*
 * Function: sopalin_smp
 *
 * Factorization routine.
 *
 * This function is meant to be used when launching a thread.
 *
 * Parameters:
 *       arg - Pointer to a data structure <sopthread_data_t> with a
 *                <Sopalin_Data_t> pointer as *data*.
 *
 */
#define sopalin_smp API_CALL(sopalin_smp)
void* sopalin_smp(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  SolverMatrix     *datacode     = sopalin_data->datacode;
  SopalinParam     *sopar        = sopalin_data->sopar;
  Thread_Data_t    *thread_data;
  PASTIX_INT               me           = argument->me;
  PASTIX_INT               i            = 0;
#ifndef PASTIX_DYNSCHED
  PASTIX_INT               ii           = 0;
#endif
  PASTIX_INT               nbpivotT     = 0;
  int               init;
  double            mintime, maxtime;
/*   PASTIX_INT               cptinv = 0; */
#if (!(defined FORCE_NOMPI))
  MPI_Comm          pastix_comm = PASTIX_COMM;
#ifdef TEST_IRECV
  PASTIX_INT               size;
#endif
#endif

#if (defined PASTIX_DYNSCHED)
  PASTIX_INT bloknum;
  PASTIX_INT itasktab  = me;
  PASTIX_INT itasktab2 = me;
  int stolen = 0;
#endif

  MONOTHREAD_BEGIN;
  trace_begin_task(sopalin_data->tracefile,
                                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                                   STATE_L0_FACTOINIT, 0);
  MONOTHREAD_END;

  /* Initialisation des données propres à chaque thread */
  print_debug(DBG_SOPALIN_MAIN, "----- %2d : START init sopalin smp\n",
              (int)me);
  init = INIT_COMPUTE;
  if (THREAD_FUNNELED_OFF)
    {
      init = init | INIT_SEND;
    }
  if (THREAD_COMM_OFF)
    {
      init = init | INIT_RECV;
    }

  sopalin_init_smp(sopalin_data, me, 1, init);
  thread_data = sopalin_data->thread_data[me];
  print_debug(DBG_SOPALIN_MAIN, "----- %2d : END init sopalin smp\n", (int)me);

#ifdef TEST_IRECV
  size = PACKMAX*(sizeof(PASTIX_INT)*MAXINFO)+PACKAREA*sizeof(PASTIX_FLOAT);
  for (i=0;i<MAX_R_REQUESTS;i++)
    {
      CALL_MPI MPI_Irecv(thread_data->recv_fanin_buffer[i],size,MPI_BYTE,
                         MPI_ANY_SOURCE,me,pastix_comm,
                         &(thread_data->recv_fanin_request[i]));
      TEST_MPI("MPI_Irecv");
    }
  size=sizeof(PASTIX_INT)*(BTAGINFO+BCOFINFO)+sizeof(PASTIX_FLOAT)*SOLV_BPFTMAX;
  for (i=0;i<MAX_R_REQUESTS;i++)
    {
      CALL_MPI MPI_Irecv(thread_data->recv_block_buffer[i],size,MPI_BYTE,MPI_ANY_SOURCE,
                         SEPFB+me,pastix_comm,&(thread_data->recv_block_request[i]));
      TEST_MPI("MPI_Irecv");
    }
#endif /* TEST_IRECV */

#ifdef HPM_SOPALIN
  hpmInit(SOLV_PROCNUM,"sopalin");
#endif

  /* Synchro de fin d'initialisation */
  SYNCHRO_THREAD;
  MONOTHREAD_BEGIN;

#ifdef PASTIX_DUMP_FACTO
  API_CALL(dump_all)(datacode, sopar->cscmtx,
                     ((datacode->updovct.sm2xtab!=NULL)?
                      (DUMP_CSC | DUMP_SOLV | DUMP_SMB):
                      (DUMP_CSC | DUMP_SOLV)));
#endif

  if (THREAD_FUNNELED_OFF)
    {
      CALL_MPI MPI_Barrier(pastix_comm);
      TEST_MPI("MPI_Barrier");
    }

  if (THREAD_COMM_ON)
    {
      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      sopalin_data->step_comm = COMMSTEP_FACTO;
      print_debug(DBG_THCOMM, "%s:%d FACTO\n", __FILE__, __LINE__);
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
      pthread_cond_broadcast(&(sopalin_data->cond_comm));
    }

  ooc_set_step(sopalin_data, OOCSTEP_SOPALIN);

  trace_begin_task(sopalin_data->tracefile,
        SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
        STATE_L0_FACTO, 0);

  MONOTHREAD_END;
  SYNCHRO_THREAD;
  SOPALIN_CLOCK_INIT; /* Debut du compteur pour le temps de facto */
  COMM_CLOCK_INIT;
  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_one("%s", OUT2_SOP_BSOP);
  print_debug(DBG_SOPALIN_MAIN,"----- [%d]%2d: sopalin starting...\n",
        (int) SOLV_PROCNUM, (int)me);

#ifdef COMPUTE_ALLOC
  ALLOC_CLOCK_INIT;
#endif

  /*****************************************************/
  /*            Main Loop                              */
  /*****************************************************/

#if defined PASTIX_DYNSCHED
  while(1){
    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                     STATE_WAITTSK, i);

    i = API_CALL(sopalin_dynsched_getNexTask)( sopalin_data, datacode, thread_data,
                                               &itasktab, &itasktab2, &bloknum, me );

    stolen = itasktab != itasktab2;
    if ( i == -1 )
      break;

    /* GEMM tasks from ESP option */
    else if (i < -1)
      {
        i = TASK_ESP2TASK( i );
#ifdef ESP_WRITE
        trace_begin_task1(thread_data->tracefile,
                          SOPALIN_CLOCK_TRACE,
                          SOLV_PROCNUM, me,
                          STATE_COMP1DGEMM,
                          SOLV_PROCDIAG( SYMB_CBLKNUM( bloknum ) ),
                          SOLV_TASKTAB[i],
                          stolen );
#else
        trace_begin_task1(thread_data->tracefile,
                          SOPALIN_CLOCK_TRACE,
                          SOLV_PROCNUM, me,
                          STATE_COMP1DGEMM,
                          TASK_PROC( i ),
                          SOLV_TASKTAB[i],
                          stolen );
#endif

        API_CALL(compute_1dgemm)(sopalin_data, me, i, bloknum, -1);
        trace_end_task1();

        MUTEX_LOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
        sopalin_data->tasktab_indice[itasktab2]++;
        MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
        continue;
      }

#elif (defined SMP_SOPALIN)
    for (ii=0;ii<SOLV_TTSKNBR;ii++){
        i = SOLV_TTSKTAB(ii);

#else /* DYNSCHED, SMP_SOPALIN */
        for (ii=0;ii<SOLV_TASKNBR;ii++){
        i = queueGet(&(sopalin_data->taskqueue));
#endif /* DYNSCHED, SMP_SOPALIN */

#ifdef COMPUTE_ALLOC
        ALLOC_CLOCK_STOP;
        printf("Step %lf memsize %lf\n",ALLOC_CLOCK_GET,
                   100.0*((double)(sopalin_data->current_alloc))/((double)(SOLV_COEFNBR)));
#endif /* COMPUTE_ALLOC */

        print_debug(DBG_SOPALIN_MAIN,
                                "[%ld]%ld: Task %ld\n"
                                "[%ld]%ld: taskid prionum cblknum bloknum ctrcnt btagptr"
                                " indnum tasknext\n"
                                "[%ld]%ld: %ld %ld %ld %ld %ld %7x %ld %ld (%ld %ld %ld %ld)\n",
                                (long)SOLV_PROCNUM,(long)me,
                                (long)i,(long)SOLV_PROCNUM,(long)me,
                                (long)SOLV_PROCNUM,(long)me,
                                (long)TASK_TASKID(i),(long)TASK_PRIONUM(i),
                                (long)TASK_CBLKNUM(i),(long)TASK_BLOKNUM(i),
                                (long)TASK_CTRBCNT(i),(unsigned int)(intptr_t)TASK_BTAGPTR(i),
                                (long)TASK_INDNUM(i),(long)TASK_TASKNEXT(i),
                                (long)SYMB_FCOLNUM(TASK_CBLKNUM(i)),
                                (long)SYMB_LCOLNUM(TASK_CBLKNUM(i)),
                                (long)SYMB_FROWNUM(TASK_BLOKNUM(i)),
                                (long)SYMB_LROWNUM(TASK_BLOKNUM(i)));

        trace_begin_task(thread_data->tracefile,
                                         SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                         STATE_WAITREM, i);

        /* Compute task */
        switch(TASK_TASKID(i))
          {
          case COMP_1D:
#ifdef HPM_SOPALIN
        hpmStart(COMP_1D+1,"COMP_1D");
        hpmStart(6,"RECV");
#endif /* HPM_SOPALIN */

        /* Wait for contributions */
        API_CALL(wait_contrib_comp_1d)(sopalin_data, me, i);

#ifdef HPM_SOPALIN
        hpmStop(6);
#endif /* HPM_SOPALIN */
        ooc_wait_task(sopalin_data,i, me);
        ooc_wait_for_cblk(sopalin_data, TASK_CBLKNUM(i),me);


        trace_begin_task1(thread_data->tracefile,
                          SOPALIN_CLOCK_TRACE,
                          SOLV_PROCNUM, me,
                          STATE_COMP1D,
                          TASK_PROC( i ),
                          SOLV_TASKTAB[i],
                          stolen );

        /* Compute */
        API_CALL(compute_1d)(sopalin_data, me, i);

        trace_end_task1();

        ooc_save_coef(sopalin_data, i, TASK_CBLKNUM(i), me);

#ifdef HPM_SOPALIN
        hpmStop(COMP_1D+1);
#endif /* HPM_SOPALIN */
        break;

          case DIAG:
#ifdef HPM_SOPALIN
        hpmStart(DIAG+1,"DIAG");
        hpmStart(6,"RECV");
#endif /* HPM_SOPALIN */
        API_CALL(wait_contrib_comp_1d)(sopalin_data, me,i);

#ifdef HPM_SOPALIN
        hpmStop(6);
#endif /* HPM_SOPALIN */

        ooc_wait_for_cblk(sopalin_data, TASK_CBLKNUM(i),me);

        trace_begin_task1(thread_data->tracefile,
                          SOPALIN_CLOCK_TRACE,
                          SOLV_PROCNUM, me,
                          STATE_DIAG,
                          TASK_PROC( i ),
                          SOLV_TASKTAB[i],
                          stolen );

        API_CALL(compute_diag)(sopalin_data, me, i);

        trace_end_task1();

        ooc_save_coef(sopalin_data, i, TASK_CBLKNUM(i),me);

#ifdef HPM_SOPALIN
        hpmStop(DIAG+1);
#endif /* HPM_SOPALIN */
        break;
          case E1:
#ifdef HPM_SOPALIN
        hpmStart(E1+1,"E1");
        hpmStart(6,"RECV");
#endif

        API_CALL(wait_contrib_comp_1d)(sopalin_data, me, i);

        API_CALL(wait_contrib_comp_2d)(sopalin_data, me, i);

#ifdef HPM_SOPALIN
        hpmStop(6);
#endif

        ooc_wait_for_cblk(sopalin_data, TASK_CBLKNUM(i),me);

        trace_begin_task1(thread_data->tracefile,
                          SOPALIN_CLOCK_TRACE,
                          SOLV_PROCNUM, me,
                          STATE_E1,
                          TASK_PROC( i ),
                          SOLV_TASKTAB[i],
                          stolen );

        API_CALL(compute_e1)(sopalin_data, me, i);

        trace_end_task1();

        ooc_save_coef(sopalin_data, i, TASK_CBLKNUM(i),me);

#ifdef HPM_SOPALIN
        hpmStop(E1+1);
#endif
        break;

          case E2:
#ifdef HPM_SOPALIN
        hpmStart(E2+1,"E2");
        hpmStart(6,"RECV");
#endif

        API_CALL(wait_contrib_comp_2d)(sopalin_data, me, i);

#ifdef HPM_SOPALIN
        hpmStop(6);
#endif

        ooc_wait_for_cblk(sopalin_data, TASK_CBLKNUM(i), me);

        trace_begin_task1(thread_data->tracefile,
                          SOPALIN_CLOCK_TRACE,
                          SOLV_PROCNUM, me,
                          STATE_E2,
                          TASK_PROC( i ),
                          SOLV_TASKTAB[i],
                          stolen );

        API_CALL(compute_e2)(sopalin_data, me, i);

        trace_end_task1();

        ooc_save_coef(sopalin_data, i, TASK_CBLKNUM(i), me);

#ifdef HPM_SOPALIN
        hpmStop(E2+1);
#endif
        break;

          default:
        errorPrint("Taskid unknown for task %ld\n", (long)i);
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
          }

#ifdef FORCE_CONSO
        API_CALL(rcsd_testall_fab)(sopalin_data, me);
#endif

#ifdef PASTIX_DYNSCHED
        MUTEX_LOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
        sopalin_data->tasktab_indice[itasktab2]++;
        MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
#endif
        trace_begin_task(thread_data->tracefile,
                         SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                         STATE_IDLE, i);
  } /* FIN boucle Principale */

#ifdef _UNUSED_
  }}}
#endif
  /* Sauvegarde du temps de facto */
  SOPALIN_CLOCK_STOP;

  trace_begin_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                   STATE_IDLE, i);

  print_debug(DBG_SOPALIN_MAIN,"----- %2d-%2d : sopalin time %lf\n",
                  (int)SOLV_PROCNUM, (int)me, SOPALIN_CLOCK_GET);

/*   fprintf(stdout, "%d - %d : Nombre d'inversion %d\n", (int)SOLV_PROCNUM, (int)me, (int)cptinv); */

  /* Suppression des Comms lancées inutilement */
#if (defined TEST_IRECV) && !(defined FORCE_NOMPI)
  for (i=0; i<MAX_R_REQUESTS; i++)
        {
          int flag;
          MPI_Status status;

          CALL_MPI MPI_Cancel(&thread_data->recv_fanin_request[i]);
          TEST_MPI("MPI_Cancel");
          CALL_MPI MPI_Test(&thread_data->recv_fanin_request[i], &flag, &status);
          TEST_MPI("MPI_Test");

          CALL_MPI MPI_Cancel(&thread_data->recv_block_request[i]);
          TEST_MPI("MPI_Cancel");
          CALL_MPI MPI_Test(&thread_data->recv_block_request[i], &flag, &status);
          TEST_MPI("MPI_Test");
        }
#endif

  /* Fin HPM */
#ifdef HPM_SOPALIN
  hpmTerminate(SOLV_PROCNUM);
#endif

#if (defined TEST_ISEND) && !(defined FORCE_NOMPI)
if (THREAD_FUNNELED_OFF)
  {
    /* Attente de la fin des communications en envoi */
    if (SOLV_PROCNBR > 1)
      API_CALL(send_waitall_fab)(sopalin_data, me);
  }
#endif /* Fin attente comm */

#ifdef TRYLOCK
  print_debug(DBG_SOPALIN_MAIN,"----- %2d : TRYLOCK free = %4ld / busy = %4ld / wait = %4ld\n",
                  (int)me, (int)thread_data->ptfree, (int)thread_data->ptbusy, (int)thread_data->ptwait);
#endif

  trace_finish(thread_data->tracefile, SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me);
  /* Nettoyage des variables associées aux threads */
  sopalin_clean_smp(sopalin_data, me);

  /* Synchro de tous les threads avant calcul globaux par le thread principal */
  SYNCHRO_THREAD;

  print_debug(DBG_SOPALIN_MAIN,"%d-%d : Synchro Avant reduction resultats\n",
                  (int)SOLV_PROCNUM, (int)me);

  /*******************************************/
  /*           Partie Monothread             */
  /*******************************************/
  MONOTHREAD_BEGIN;

  trace_begin_task(sopalin_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                   STATE_L0_FACTOCLEAN, 0);

  /* Envoi du message de fin aux threads de comm */
#ifndef FORCE_NOMPI
  if (THREAD_COMM_ON && THREAD_FUNNELED_OFF)
    {
      for (i=0; i < SOLV_PROCNBR; i++)
        {
          int tag, iterator;
          int miun = -1;
          if (i == SOLV_PROCNUM) continue;
          for (iterator=0; iterator < sopar->nbthrdcomm; iterator++)
            {
              tag = (sopar->type_comm == 3) ? iterator : TAG_FANIN;
              CALL_MPI MPI_Send(&miun, 1, MPI_INT, i, tag, PASTIX_COMM);
              TEST_MPI("MPI_Send");
              /*fprintf(stderr," %d : Envoi %d à %d\n", SOLV_PROCNUM, iterator, i);*/
            }
        }
    }
#endif

  /* Calcul du nombre de pivotage réalisé */
  for(i= 0; i< SOLV_THRDNBR; i++)
        nbpivotT += sopalin_data->thread_data[i]->nbpivot;
  sopar->diagchange = nbpivotT;

  /* Calcul du temps de facto */
  mintime = thread_data->sop_clk.time[0];
  maxtime = thread_data->sop_clk.time[1];
  for(i=1; i<SOLV_THRDNBR; i++)
        {
          mintime = MIN(sopalin_data->thread_data[i]->sop_clk.time[0], mintime);
          maxtime = MAX(sopalin_data->thread_data[i]->sop_clk.time[1], maxtime);
        }
  sopar->dparm[DPARM_FACT_TIME] = (maxtime - mintime);

  /* WARNING : Don't put one (All)Reduce before thread synchronization FACTOEND */
  if (THREAD_COMM_ON)
    {
      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      while(sopalin_data->step_comm != COMMSTEP_FACTOEND)
        COND_WAIT(&(sopalin_data->cond_comm), &(sopalin_data->mutex_comm));
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
    }

  /* Calcul de l'inertie de la matrice (pour CROUT seulement, sans pivotage) */
  sopar->iparm[IPARM_INERTIA] = -1;
#if (!defined TYPE_COMPLEX) && (!defined CHOL_SOPALIN) && (!defined OOC)
  {
        PASTIX_FLOAT *ga;
        PASTIX_INT c, k, stride, size, inertia;
        inertia=0;
        for (c=0;c<SYMB_CBLKNBR;c++)
          {
        ga     =&(SOLV_COEFTAB(c)[SOLV_COEFIND(SYMB_BLOKNUM(c))]);
        stride = SOLV_STRIDE(c);
        size   = SYMB_LCOLNUM(c)-SYMB_FCOLNUM(c)+1;
        for (k=0;k<size;k++)
          if (ga[k+k*stride]>fzero) inertia++;
          }
        MyMPI_Allreduce(&inertia,&(sopar->iparm[IPARM_INERTIA]),1,
                        COMM_INT,MPI_SUM,pastix_comm);
  }
#endif

#ifdef PASTIX_DYNSCHED
  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        {
          for(i=1; i<SOLV_THRDNBR; i++)
        {
          thread_data->esp += sopalin_data->thread_data[i]->esp;
        }
          MyMPI_Reduce(&(thread_data->esp), &(sopar->iparm[IPARM_ESP_NBTASKS]), 1,
                   COMM_INT, MPI_SUM, 0, pastix_comm);
        }
#endif

if (THREAD_FUNNELED_OFF)
  {
    CALL_MPI MPI_Barrier(pastix_comm);
    TEST_MPI("MPI_Barrier");
  }

#ifdef STATS_SOPALIN
  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        {
          double overhead     = (double)(sopalin_data->max_alloc + SOLV_COEFNBR)/(double)(SOLV_COEFNBR);

          fprintf(stdout, "%2d - Local number of terms allocated\t Cblk+Ftgt : %10ld,\t Cblk : %10ld,\t Overhead : %.2lf (%.2lf%%)\n",
                  (int)SOLV_PROCNUM,
                  (long)(sopalin_data->max_alloc + SOLV_COEFNBR),
                  (long)(SOLV_COEFNBR),
                  overhead,
                  (overhead - 1.0) * 100.0);
          {
        PASTIX_INT    tmp_max_alloc = sopalin_data->max_alloc;
        PASTIX_INT    tmp_coefnbr   = SOLV_COEFNBR;
        PASTIX_INT    max_max_alloc = 0;
        PASTIX_INT    max_coefnbr   = 0;
        PASTIX_INT    sum_max_alloc = 0;
        PASTIX_INT    sum_coefnbr   = 0;
        double overhead2;

        MyMPI_Allreduce(&tmp_max_alloc,&max_max_alloc,1,COMM_INT,MPI_MAX,pastix_comm);
        MyMPI_Allreduce(&tmp_coefnbr,  &max_coefnbr,  1,COMM_INT,MPI_MAX,pastix_comm);
        MyMPI_Allreduce(&tmp_max_alloc,&sum_max_alloc,1,COMM_INT,MPI_SUM,pastix_comm);
        MyMPI_Allreduce(&tmp_coefnbr,  &sum_coefnbr,  1,COMM_INT,MPI_SUM,pastix_comm);

        overhead2 = (double)(max_max_alloc+max_coefnbr)/(double)(max_coefnbr);
        print_one("Maximum number of terms allocated\t Cblk+Ftgt : %10ld,\t Cblk : %10ld,\t Overhead : %.2lf (%.2lf%%)\n",
                  (long)(max_max_alloc+max_coefnbr),
                  (long) max_coefnbr,
                  overhead2,
                  (overhead2 - 1.0) * 100.0);

        overhead2 = (double)(sum_max_alloc+sum_coefnbr)/(double)(sum_coefnbr);
        print_one("Total number of terms allocated\t\t Cblk+Ftgt : %10ld,\t Cblk : %10ld,\t Overhead : %.2lf (%.2lf%%)\n",
                  (long)(sum_max_alloc+sum_coefnbr),
                  (long) sum_coefnbr,
                  overhead2,
                  (overhead2 - 1.0) * 100.0);
        sopar->iparm[IPARM_ALLOCATED_TERMS]=sum_max_alloc+sum_coefnbr;
          }
        }
#endif

  /* free file structures */
  sopalin_clean(sopalin_data, 1);

#ifdef PASTIX_DUMP_FACTO
API_CALL(dump_all)(datacode, sopar->cscmtx, DUMP_SOLV);
#endif

  /* Fin des threads de comms et d'OOC */
  if (THREAD_COMM_ON)
    {
      if ((sopar->iparm[IPARM_END_TASK] == API_TASK_NUMFACT)
          || (sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0))
        {
          MUTEX_LOCK(&(sopalin_data->mutex_comm));
          sopalin_data->step_comm = COMMSTEP_END;
          print_debug(DBG_THCOMM, "%s:%d END\n", __FILE__, __LINE__);
          MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
          pthread_cond_broadcast(&(sopalin_data->cond_comm));
        }
    }

#ifdef OOC
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] < API_TASK_SOLVE)
        ooc_stop_thread(sopalin_data);
#endif

  trace_begin_task(sopalin_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                   STATE_L0_IDLE, 0);

  MONOTHREAD_END;
  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_one("%s", OUT2_SOP_ESOP);

  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX)
        fprintf(stdout, OUT4_FACT_COMM_TIME,
                (int)SOLV_PROCNUM, (int)me, COMM_CLOCK_GET);
  return 0;
}

/******************************************************************************/
/*         Fonction pour les threads de communication                         */
/******************************************************************************/
/*
  Function: API_CALL(sopalin_updo_comm)

  Updown function used for the communication thread.

  Parameters:
        arg - Pointer to a data structure <sopthread_data_t> with a
                  <Sopalin_Data_t> pointer as *data*.
*/
#define sopalin_updo_comm API_CALL(sopalin_updo_comm)
void *sopalin_updo_comm ( void *arg )
{
#ifndef FORCE_NOMPI
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  SolverMatrix     *datacode     = sopalin_data->datacode;
  PASTIX_INT               me           = argument->me;
  if (THREAD_COMM_ON)
    {
      /* Thread_Data_t    *thread_data  = sopalin_data->thread_data[me]; */

      /* On se met en attente du debut de la descente
       * ou on quitte si on ne reitere pas */
      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      while(1)
        {
          switch(sopalin_data->step_comm)
            {
              /* Il n'y a plus de comm */
            case COMMSTEP_END:
              sopalin_data->step_comm = COMMSTEP_INIT;
              print_debug(DBG_THCOMM, "%s:%d INIT\n", __FILE__, __LINE__);
              MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
              print_debug(DBG_SOPALIN_THREADCOMM,
                          "%d : je quitte\n", (int)SOLV_PROCNUM);
              return (void *)1;
              break;

              /* Factorisation */
            case COMMSTEP_FACTO:
              MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
              API_CALL(sendrecv_smp)(arg);
              /* On ne lance qu'une facto et
               * le reste ne necessite qu'un thread pour l'instant */
              if (me > SOLV_THRDNBR) return (void *)1;
              MUTEX_LOCK(&(sopalin_data->mutex_comm));
              break;

              /* Udpo */
            case COMMSTEP_DOWN:
              MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
              API_CALL(updo_thread_comm)(arg);
              MUTEX_LOCK(&(sopalin_data->mutex_comm));
              break;

              /* Un AllReduce est a faire en funneled */
            case COMMSTEP_ALLREDUCE:
              if (THREAD_FUNNELED_ON)
                {
                  {
                    Pastix_Allreduce_t *allreduce = &(sopalin_data->allreduce);
                    MPI_Allreduce(allreduce->sendbuf,
                                  allreduce->recvbuf,
                                  allreduce->count,
                                  allreduce->datatype,
                                  allreduce->op,
                                  PASTIX_COMM);
                    sopalin_data->step_comm = COMMSTEP_INIT;
                    print_debug(DBG_THCOMM, "%s:%d INIT\n", __FILE__, __LINE__);
                    pthread_cond_broadcast(&(sopalin_data->cond_comm));
                  }
                  break;
                }
            case COMMSTEP_REDUCE:
              if (THREAD_FUNNELED_ON)
                {
                  {
                    Pastix_Allreduce_t *reduce = &(sopalin_data->allreduce);
                    MPI_Reduce(reduce->sendbuf,
                               reduce->recvbuf,
                               reduce->count,
                               reduce->datatype,
                               reduce->op,
                               0,
                               PASTIX_COMM);
                    sopalin_data->step_comm = COMMSTEP_INIT;
                    print_debug(DBG_THCOMM, "%s:%d INIT\n", __FILE__, __LINE__);
                    pthread_cond_broadcast(&(sopalin_data->cond_comm));
                  }
                  break;
                }
              /* Sortie spontannée du wait */
            default:
              COND_WAIT(&(sopalin_data->cond_comm),
                        &(sopalin_data->mutex_comm));
            }
        }
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
    }
#endif
  (void)arg;
  return (void*)NULL;
}

/******************************************************************************/
/*         Fonction pour les threads de calculs                               */
/******************************************************************************/

/*
 * Function: sopalin_thread
 *
 * Function launching computing, communictating and out of core threads on
 * the factorization step only.
 *
 * Initiate the <Sopalin_Data_t> structure, launch threads, clean and restore.
 *
 * Parameters:
 *       m         - The <SolverMatrix> structure.
 *	sopaparam - Sopalin parameters in the <SopalinParam> stucture.
 */
#define sopalin_thread API_CALL(sopalin_thread)
void sopalin_thread(SolverMatrix *m,
                    SopalinParam *sopaparam)
{
  Backup b;
  Sopalin_Data_t *sopalin_data = NULL;
  SolverMatrix  *datacode = NULL;

  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  sopalin_backup(m,&b);
  sopalin_init(sopalin_data, m, sopaparam, 1);
  API_CALL(init_struct_sopalin)(sopalin_data, m, sopaparam);

  datacode = sopalin_data->datacode;

#ifdef WITH_STARPU
  if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {
      starpu_submit_tasks(sopalin_data);
    }
  else
#endif
    {
      sopalin_launch_thread(sopalin_data,
                            SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree,
                            sopalin_data->sopar->iparm[IPARM_VERBOSE],
                            SOLV_THRDNBR,          API_CALL(sopalin_smp),       sopalin_data,
                            sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                            OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);
    }
  sopalin_clean(sopalin_data, 2);
  sopalin_restore(m,&b);

  memFree_null(sopalin_data);
}

/*
  Function: sopalin_updo_smp

  Function used for computing thread creation to compute factorisation and
  resolution.

  Parameters:
        arg - Pointer to a data structure <sopthread_data_t> with a
                  <Sopalin_Data_t> pointer as *data*.
*/
void* API_CALL(sopalin_updo_smp)(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  PASTIX_INT               me           = argument->me;

  API_CALL(sopalin_smp)(argument);
  if (sopalin_data->sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
        {
          if ((sopalin_data->datacode->clustnum == 0) && (me == 0))
      errorPrintW("Updown incompatible with 2D distribution");
          return 0;
        }

  MONOTHREAD_BEGIN;
  sopalin_init(sopalin_data, NULL, NULL, 0);
  MONOTHREAD_END;
  API_CALL(up_down_smp)(argument);

  return 0;
}
/*
  Function: API_CALL(sopalin_updo_thread)

  Function launching computing, communicating and out of core threads on
  the factorization and solve steps.

  Initiate the <Sopalin_Data_t> structure, launch threads, clean and restore.

  Parameters:
        m         - The <SolverMatrix> structure.
        sopaparam - Sopalin parameters in the <SopalinParam> stucture.
*/
void API_CALL(sopalin_updo_thread)(SolverMatrix *m,
                                   SopalinParam *sopaparam)
{
  Backup b;
  Sopalin_Data_t *sopalin_data = NULL;
  SolverMatrix   *datacode = NULL;

  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  sopalin_backup(m,&b);
  sopalin_init(sopalin_data, m, sopaparam, 1);
  API_CALL(init_struct_sopalin)(sopalin_data, m, sopaparam);
  datacode = sopalin_data->datacode;
#ifdef WITH_STARPU
  if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {
      starpu_submit_tasks(sopalin_data);

    }
  else
#endif
    {
      sopalin_launch_thread(sopalin_data,
                            SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                            SOLV_THRDNBR,          API_CALL(sopalin_updo_smp),  sopalin_data,
                            sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                            OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);
    }

  sopalin_clean(sopalin_data, 2);
  sopalin_restore(m,&b);

  memFree_null(sopalin_data);
}

/*
  Function: API_CALL(sopalin_updo_gmres_smp)

  Function used for computing thread creation to compute factorisation,
  resolution and gmres.

  Parameters:
        arg - Pointer to a data structure <sopthread_data_t> with a
                  <Sopalin_Data_t> pointer as *data*.
*/
void* API_CALL(sopalin_updo_gmres_smp)(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  PASTIX_INT               me           = argument->me;

  API_CALL(sopalin_smp)(argument);
  if (sopalin_data->sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
        {
          if ((sopalin_data->datacode->clustnum == 0) && (me == 0))
            errorPrintW("Updown incompatible with 2D distribution");
          return 0;
        }

  MONOTHREAD_BEGIN;
  sopalin_init(sopalin_data, NULL, NULL, 0);
  MONOTHREAD_END;
  API_CALL(gmres_smp)(argument);

  return 0;
}
/*
  Function: API_CALL(sopalin_updo_gmres_thread)

  Function launching computing, communicating and out of core threads on
  the factorization, solve and reffinement (using GMRES) steps.

  Initiate the <Sopalin_Data_t> structure, launch threads, clean and restore.

  Parameters:
        m         - The <SolverMatrix> structure.
        sopaparam - Sopalin parameters in the <SopalinParam> stucture.
*/
void API_CALL(sopalin_updo_gmres_thread)(SolverMatrix *m, SopalinParam *sopaparam)
{
  Backup b;
  Sopalin_Data_t *sopalin_data;
  SolverMatrix   *datacode = m;



  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  sopalin_backup(m,&b);
  sopalin_init(sopalin_data, m, sopaparam, 1);
  API_CALL(init_struct_sopalin)(sopalin_data, m, sopaparam);
  datacode = sopalin_data->datacode;
#ifdef WITH_STARPU
  if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {

      starpu_submit_tasks(sopalin_data);

    }
  else
#endif
    {
      sopalin_launch_thread(sopalin_data,
                            SOLV_PROCNUM,          SOLV_PROCNBR,                     datacode->btree,
                            sopalin_data->sopar->iparm[IPARM_VERBOSE],
                            SOLV_THRDNBR,          API_CALL(sopalin_updo_gmres_smp), sopalin_data,
                            sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm),      sopalin_data,
                            OOC_THREAD_NBR,        ooc_thread,                       sopalin_data);
    }
  sopalin_clean(sopalin_data, 2);
  sopalin_restore(m,&b);

  memFree_null(sopalin_data);
}

/*
  Function: API_CALL(sopalin_updo_grad_smp)

  Function used for computing thread creation to compute factorisation,
  resolution and conjugate gradient.

  Parameters:
        arg - Pointer to a data structure <sopthread_data_t> with a
                  <Sopalin_Data_t> pointer as *data*.
*/
void* API_CALL(sopalin_updo_grad_smp)(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  PASTIX_INT               me           = argument->me;

  API_CALL(sopalin_smp)(argument);
  if (sopalin_data->sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
        {
          if ((sopalin_data->datacode->clustnum == 0) && (me == 0))
        errorPrintW("Updown incompatible with 2D distribution");
          return 0;
        }

  MONOTHREAD_BEGIN;
  sopalin_init(sopalin_data, NULL, NULL, 0);
  MONOTHREAD_END;
  API_CALL(up_down_smp)(argument);
  API_CALL(grad_smp)   (argument);

  return 0;
}

/*
  Function: API_CALL(sopalin_updo_grad_thread)

  Function launching computing, communicating and out of core threads on
  the factorization, solve and reffinement (using conjugate grandient) steps.

  Initiate the <Sopalin_Data_t> structure, launch threads, clean and restore.

  Parameters:
        m         - The <SolverMatrix> structure.
        sopaparam - Sopalin parameters in the <SopalinParam> stucture.
*/
void API_CALL(sopalin_updo_grad_thread)(SolverMatrix *m, SopalinParam *sopaparam)
{
  Backup b;
  Sopalin_Data_t *sopalin_data = NULL;
  SolverMatrix   *datacode = NULL;

  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  sopalin_backup(m,&b);
  sopalin_init(sopalin_data, m, sopaparam, 1);
  API_CALL(init_struct_sopalin)(sopalin_data, m, sopaparam);
  datacode = sopalin_data->datacode;
#ifdef WITH_STARPU
  if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {

      starpu_submit_tasks(sopalin_data);

    }
  else
#endif
    {
      sopalin_launch_thread(sopalin_data,
                            SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                            SOLV_THRDNBR,          API_CALL(sopalin_updo_grad_smp), sopalin_data,
                            sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm),     sopalin_data,
                            OOC_THREAD_NBR,        ooc_thread,                      sopalin_data);
    }
  sopalin_clean(sopalin_data, 2);
  sopalin_restore(m,&b);

  memFree_null(sopalin_data);
}

/*
  Function: API_CALL(sopalin_updo_pivot_smp)

  Function used for computing thread creation to compute factorisation,
  resolution and pivoting refinement.

  Parameters:
        arg - Pointer to a data structure <sopthread_data_t> with a
                  <Sopalin_Data_t> pointer as *data*.
*/
void* API_CALL(sopalin_updo_pivot_smp)(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  PASTIX_INT               me           = argument->me;

  API_CALL(sopalin_smp)(argument);
  if (sopalin_data->sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
        {
          if ((sopalin_data->datacode->clustnum == 0) && (me == 0))
        errorPrintW("Updown incompatible with 2D distribution");
          return 0;
        }

  MONOTHREAD_BEGIN;
  sopalin_init(sopalin_data, NULL, NULL, 0);
  MONOTHREAD_END;
  API_CALL(up_down_smp)(argument);
  API_CALL(pivotstatique_smp)(argument);

  return 0;
}

/*
  Function: API_CALL(sopalin_updo_pivot_thread)

  Function launching computing, communicating and out of core threads on
  the factorization, solve and reffinement (using pivoting refinement) steps.

  Initiate the <Sopalin_Data_t> structure, launch threads, clean and restore.

  Parameters:
        m         - The <SolverMatrix> structure.
        sopaparam - Sopalin parameters in the <SopalinParam> stucture.
*/
void API_CALL(sopalin_updo_pivot_thread)(SolverMatrix *m, SopalinParam *sopaparam)
{
  Backup b;
  Sopalin_Data_t *sopalin_data = NULL;
  SolverMatrix   *datacode = NULL;

  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  sopalin_backup(m,&b);
  sopalin_init(sopalin_data, m, sopaparam, 1);
  datacode = sopalin_data->datacode;
#ifdef WITH_STARPU
  if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {

      starpu_submit_tasks(sopalin_data);

    }
  else
#endif
    {
      sopalin_launch_thread(sopalin_data,
                            SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                            SOLV_THRDNBR,          API_CALL(sopalin_updo_pivot_smp), sopalin_data,
                            sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm),      sopalin_data,
                            OOC_THREAD_NBR,        ooc_thread,                       sopalin_data);
    }
  sopalin_clean(sopalin_data, 2);
  sopalin_restore(m,&b);
}


/*
  Function: API_CALL(sopalin_updo_bicgstab_smp)

  Function used for computing thread creation to compute factorisation,
  resolution and bicgstab refinement.

  Parameters:
        arg - Pointer to a data structure <sopthread_data_t> with a
                  <Sopalin_Data_t> pointer as *data*.
*/
void* API_CALL(sopalin_updo_bicgstab_smp)(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  PASTIX_INT               me           = argument->me;

  API_CALL(sopalin_smp)(argument);
  if (sopalin_data->sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
        {
          if ((sopalin_data->datacode->clustnum == 0) && (me == 0))
        errorPrintW("Updown incompatible with 2D distribution");
          return 0;
        }

  MONOTHREAD_BEGIN;
  sopalin_init(sopalin_data, NULL, NULL, 0);
  MONOTHREAD_END;
  API_CALL(up_down_smp)(argument);
  API_CALL(bicgstab_smp)(argument);

  return 0;
}
/*
  Function: API_CALL(sopalin_updo_bicgstab_thread)

  Function launching computing, communicating and out of core threads on
  the factorization, solve and reffinement (using bicgstab refinement) steps.

  Initiate the <Sopalin_Data_t> structure, launch threads, clean and restore.

  Parameters:
        m         - The <SolverMatrix> structure.
        sopaparam - Sopalin parameters in the <SopalinParam> stucture.
*/
void API_CALL(sopalin_updo_bicgstab_thread)(SolverMatrix *m, SopalinParam *sopaparam)
{
  Backup b;
  Sopalin_Data_t *sopalin_data = NULL;
  SolverMatrix   *datacode = NULL;

  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  sopalin_backup(m,&b);
  sopalin_init(sopalin_data, m, sopaparam, 1);
  datacode = sopalin_data->datacode;
#ifdef WITH_STARPU
  if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {

      starpu_submit_tasks(sopalin_data);

    }
  else
#endif
    {
      sopalin_launch_thread(sopalin_data,
                            SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                            SOLV_THRDNBR,          API_CALL(sopalin_updo_bicgstab_smp), sopalin_data,
                            sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm),      sopalin_data,
                            OOC_THREAD_NBR,        ooc_thread,                       sopalin_data);
    }
  sopalin_clean(sopalin_data, 2);
  sopalin_restore(m,&b);
}

/*
  Function: API_CALL(sopalin_launch)

  TODO: Comment (unused ?)
 */
void API_CALL(sopalin_launch)(SolverMatrix *m,
                                  SopalinParam *sopaparam,
                                  PASTIX_INT cas)
{
  Backup b;
  Sopalin_Data_t *sopalin_data = NULL;
  SolverMatrix   *datacode     = NULL;

  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  if (cas < UPDO_ONLY)
        {
          sopalin_backup(m,&b);
          sopalin_init(sopalin_data, m, sopaparam, 1);
          API_CALL(init_struct_sopalin)(sopalin_data, m, sopaparam);
        }
  else
        {
          sopalin_init(sopalin_data, m, sopaparam, 0);
        }

  datacode = sopalin_data->datacode;
  switch(cas){
  case SOPALIN_ONLY:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(sopalin_smp),       sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);
        break;
  case SOPALIN_UPDO:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(sopalin_updo_smp),  sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);
        break;
  case SOPALIN_UPDO_GMRES:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(sopalin_updo_gmres_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm),      sopalin_data,
                              OOC_THREAD_NBR,        ooc_thread,                       sopalin_data);
        break;
  case SOPALIN_UPDO_GRAD:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(sopalin_updo_grad_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm),     sopalin_data,
                              OOC_THREAD_NBR,        ooc_thread,                      sopalin_data);
        break;
  case SOPALIN_UPDO_PIVOT:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(sopalin_updo_pivot_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm),      sopalin_data,
                              OOC_THREAD_NBR,        ooc_thread,                       sopalin_data);
        break;
  case UPDO_ONLY:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(up_down_smp),       sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);
        break;
  case RAFF_GMRES:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(gmres_smp),         sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);
        break;
  case RAFF_GRAD:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(grad_smp),          sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);
        break;
  case RAFF_PIVOT:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(pivotstatique_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);
        break;
  default:
        if( SOLV_PROCNUM == 0 )
          {
        errorPrint("undefined case.");
        EXIT(MOD_SOPALIN,BADPARAMETER_ERR);
          }
  }

  sopalin_clean(sopalin_data, 2);
  if (cas < UPDO_ONLY)
          sopalin_restore(m,&b);

  memFree_null(sopalin_data);
}
