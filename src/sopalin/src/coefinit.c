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
 * File: coefinit.c
 *
 * Allocation and initialisation of the coeficient of the solver matrix.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include "common_pastix.h"
#ifndef FORCE_NOSMP
#  include <pthread.h>
#endif
#ifdef FORCE_NOMPI
#  include "nompi.h"
#else
#  include <mpi.h>
#endif
#include "sopalin_define.h"
#include "dof.h"
#include "ftgt.h"
#include "symbol.h"
#include "csc.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"

/*
 * Section: Defines
 *
 * Define: COEFTAB_TYPE
 *
 * Type allocated in coeftab: PASTIX_FLOAT *
 *
 * Define: COEFTAB_SIZE
 *
 * Size of the array coeftab: number of column blocks
 */
#define COEFTAB_TYPE PASTIX_FLOAT *
#define COEFTAB_SIZE SYMB_CBLKNBR

#include "sopalin_thread.h"
#include "stack.h"
#include "sopalin3d.h"
#include "order.h"
#include "debug_dump.h"
#include "csc_intern_solve.h"
#include "csc_intern_build.h"
#include "sopalin_time.h"
#include "ooc.h"
#include "sopalin_acces.h" /* ATTENTION : inclure apres define SMP_SOPALIN */
#include "coefinit.h"

/* Section: Functions */

/*
 * Function: CoefMatrix_Allocate
 *
 * Allocate matrix coefficients in coeftab and ucoeftab.
 *
 * Should be first called with me = -1 to allocated coeftab.
 * Then, should be called with me set to thread ID
 * to allocate column blocks coefficients arrays.
 *
 * Parameters
 *
 *    datacode  - solverMatrix
 *    factotype - factorization type (LU, LLT ou LDLT)
 *    me        - thread number. (-1 for first call,
 *                from main thread. >=0 to allocate column blocks
 *     assigned to each thread.)
 */
void CoefMatrix_Allocate(SopalinParam    *sopar,
                         SolverMatrix    *datacode,
                         pthread_mutex_t *mutex,
                         PASTIX_INT              factotype,
                         PASTIX_INT              me)
{
  PASTIX_INT i;
#ifndef OOC
  PASTIX_INT itercblk, coefnbr;
#endif
#ifdef WITH_STARPU
  /* For CUDA devices we have no allocation (yet?) */
  if ( sopar->iparm[IPARM_STARPU] == API_YES && me >= SOLV_THRDNBR)
    return;
#endif /* WITH_STARPU */
#ifndef OOC
  {
    /* On ne passe pas ici en OOC */
    PASTIX_INT bubnum  = me;
    PASTIX_INT task;

#  ifdef PASTIX_DYNSCHED
    while (bubnum != -1)
      {
        PASTIX_INT fcandnum = datacode->btree->nodetab[bubnum].fcandnum;
        PASTIX_INT lcandnum = datacode->btree->nodetab[bubnum].lcandnum;
        for (i=(me-fcandnum);i < datacode->ttsknbr[bubnum]; i+=(lcandnum-fcandnum+1))
#  else
    for (i=0; i < datacode->ttsknbr[bubnum]; i++)
#  endif /* PASTIX_DYNSCHED */

      {
        task = datacode->ttsktab[bubnum][i];
        itercblk = TASK_CBLKNUM(task);
        coefnbr  = SOLV_STRIDE(itercblk) * (SYMB_LCOLNUM(itercblk) -
                                            SYMB_FCOLNUM(itercblk) + 1);
        if ((TASK_TASKID(task) == COMP_1D)
            || (TASK_TASKID(task) == DIAG))
          {
            datacode->cblktab[itercblk].procdiag = me;
            if (SOLV_COEFTAB(itercblk) == NULL)
              { /* If not NULL it should be the schur */
                MALLOC_INTERN(SOLV_COEFTAB(itercblk), coefnbr, PASTIX_FLOAT);
              }
          }
        else if ( SOLV_COEFIND(TASK_BLOKNUM(task)) == 0 )
          {
            MUTEX_LOCK(mutex);
            datacode->cblktab[itercblk].procdiag = me;
            if (SOLV_COEFTAB(itercblk) == NULL)
              {
                MALLOC_INTERN(SOLV_COEFTAB(itercblk), coefnbr, PASTIX_FLOAT);
              }
            MUTEX_UNLOCK(mutex);
          }
      }

#  ifdef PASTIX_DYNSCHED
        bubnum = BFATHER(datacode->btree, bubnum);
      }
#  endif /* PASTIX_DYNSCHED */
  }
#endif /* OOC */

  /*
   * Allocate LU coefficient arrays
   * We also use it to store the diagonal in LDLt factorization using esp
   */
  if ( (factotype == API_FACT_LU) /* LU */
       || ( (factotype == API_FACT_LDLT) && sopar->iparm[IPARM_ESP] ) )
    {
#ifndef OOC
      {
        /* On ne passe pas ici en OOC */
        PASTIX_INT bubnum  = me;
        PASTIX_INT task;

#  ifdef PASTIX_DYNSCHED
        while (bubnum != -1)
          {
            PASTIX_INT fcandnum = datacode->btree->nodetab[bubnum].fcandnum;
            PASTIX_INT lcandnum = datacode->btree->nodetab[bubnum].lcandnum;
            for (i=(me-fcandnum); i < datacode->ttsknbr[bubnum]; i+=(lcandnum-fcandnum+1))
#  else
        for (i=0; i < datacode->ttsknbr[bubnum]; i++)
#  endif /* PASTIX_DYNSCHED */

          {
            task = datacode->ttsktab[bubnum][i];
            itercblk = TASK_CBLKNUM(task);
            if ( ((TASK_MASTER(task) != -1) && (task != TASK_MASTER(task)))
                 || (me != datacode->cblktab[itercblk].procdiag)
                 || (SOLV_UCOEFTAB(itercblk) != NULL) )
              {
                continue;
              }

            if ( (factotype == API_FACT_LDLT) && sopar->iparm[IPARM_ESP] )
              {
                coefnbr  = SYMB_LCOLNUM(itercblk) - SYMB_FCOLNUM(itercblk) + 1;
              }
            else
              {
                coefnbr  = SOLV_STRIDE(itercblk) * (SYMB_LCOLNUM(itercblk) -
                                                    SYMB_FCOLNUM(itercblk) + 1);
              }

            MALLOC_INTERN(SOLV_UCOEFTAB(itercblk), coefnbr, PASTIX_FLOAT);
          }

#  ifdef PASTIX_DYNSCHED
            bubnum = BFATHER(datacode->btree, bubnum);
          }
#  endif /* PASTIX_DYNSCHED */
      }
#endif /* OOC */
    }
  printf("fin alloc (u)coeff\n");
}

/*
 * Function: CoefMatrix_Init
 *
 * Init coeftab and ucoeftab coefficients.
 *
 * Parameters:
 *    datacode     - solverMatrix
 *    barrier      - Barrier used for thread synchronisation.
 *    me           - Thread ID
 *    iparm        - Integer parameters array.
 *    transcsc     - vecteur transcsc
 *    sopalin_data - <Sopalin_Data_t> structure.
 */
void CoefMatrix_Init(SolverMatrix         *datacode,
                     sopthread_barrier_t  *barrier,
                     PASTIX_INT                   me,
                     PASTIX_INT                  *iparm,
                     PASTIX_FLOAT               **transcsc,
                     Sopalin_Data_t       *sopalin_data)
{

  PASTIX_INT j, itercblk;
  PASTIX_INT i, coefnbr;

#ifdef WITH_STARPU
  /* For CUDA devices we have no allocation (yet?) */
  if ( iparm[IPARM_STARPU] == API_YES && me >= SOLV_THRDNBR)
    return;
#endif /* WITH_STARPU */

  /* Remplissage de la matrice */
  if (iparm[IPARM_FILL_MATRIX] == API_NO)
  {
    /* Remplissage par bloc */
    PASTIX_INT bubnum  = me;
#ifdef PASTIX_DYNSCHED
    while (bubnum != -1)
    {
      PASTIX_INT fcandnum = datacode->btree->nodetab[bubnum].fcandnum;
      PASTIX_INT lcandnum = datacode->btree->nodetab[bubnum].lcandnum;
      for (i=(me-fcandnum);i < datacode->ttsknbr[bubnum]; i+=(lcandnum-fcandnum+1))
#else
        for (i=0; i < datacode->ttsknbr[bubnum]; i++)
#endif /* PASTIX_DYNSCHED */

        {
          PASTIX_INT task;
          PASTIX_INT k = i;
#ifdef OOC
          /* En OOC, on inverse la boucle pour conserver les premiers blocs en mémoire */
          k = datacode->ttsknbr[bubnum]-i-1;
#endif
          task = datacode->ttsktab[bubnum][k];
          itercblk = TASK_CBLKNUM(task);

          if ( ((TASK_MASTER(task) != -1) && (task != TASK_MASTER(task)))
               || (me != datacode->cblktab[itercblk].procdiag))
            continue;

          coefnbr  = SOLV_STRIDE(itercblk) * (SYMB_LCOLNUM(itercblk) - SYMB_FCOLNUM(itercblk) + 1);

          ooc_wait_for_cblk(sopalin_data, itercblk, me);

          /* initialisation du bloc colonne */
          for (j=0 ; j < coefnbr ; j++)
          {
            SOLV_COEFTAB(itercblk)[j] = ZERO;
            if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
              SOLV_UCOEFTAB(itercblk)[j] = ZERO;
          }

          /* remplissage */
          Csc2solv_cblk(sopalin_data->sopar->cscmtx, datacode, *transcsc, itercblk);

          ooc_save_coef(sopalin_data, task, itercblk, me);

#if defined(PASTIX_DYNSCHED) || defined(TRACE_SOPALIN)
          TASK_CAND(task) = me;
#endif
        }

#ifdef PASTIX_DYNSCHED
      bubnum = BFATHER(datacode->btree, bubnum);
    }
#endif /* PASTIX_DYNSCHED */

    printf("fin fill-in\n");

#ifdef DEBUG_COEFINIT
    if (me == 0)
    {
      FILE *transfile;
      char transfilename[10];
      sprintf(transfilename, "trans%ld.%ld",(long) me,(long) SOLV_PROCNUM);
      transfile = fopen(transfilename, "w");
      dump7(*transcsc, transfile);
    }
#endif
    /* Libération de mémoire */
#ifdef STARPU_INIT_SMP
    if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_NO)
#endif /* STARPU_INIT_SMP */
      SYNCHRO_X_THREAD(SOLV_THRDNBR, *barrier);
    if (me == 0)
    {
      if (iparm[IPARM_FACTORIZATION] != API_FACT_LU)
      {
        if (*transcsc != NULL)
          memFree_null(*transcsc);
      }
      else
      {
        if (iparm[IPARM_SYM] == API_SYM_YES || iparm[IPARM_SYM] == API_SYM_HER) /* Symmetric */
          *transcsc = NULL;
        else /* Unsymmetric */
          memFree_null(*transcsc);
      }
    }
  }
  else  /* fake factorisation */
  {

    /* deadcode */
    PASTIX_INT itercol;

    /* Initialisation de la matrice à 0 et 1 ou 2 */
    PASTIX_INT task;
    PASTIX_INT bubnum  = me;

#ifdef PASTIX_DYNSCHED
    while (bubnum != -1)
    {
      PASTIX_INT fcandnum = datacode->btree->nodetab[bubnum].fcandnum;
      PASTIX_INT lcandnum = datacode->btree->nodetab[bubnum].lcandnum;
      for (i=(me-fcandnum);i < datacode->ttsknbr[bubnum]; i+=(lcandnum-fcandnum+1))
#else
        for (i=0; i < datacode->ttsknbr[bubnum]; i++)
#endif /* PASTIX_DYNSCHED */

        {
          task = datacode->ttsktab[bubnum][i];
          itercblk = TASK_CBLKNUM(task);
          if ( ((TASK_MASTER(task) != -1) && (task != TASK_MASTER(task)))
               || (me != datacode->cblktab[itercblk].procdiag))
            continue;
          coefnbr  = SOLV_STRIDE(itercblk) * (SYMB_LCOLNUM(itercblk) - SYMB_FCOLNUM(itercblk) + 1);

          for (j=0 ; j < coefnbr ; j++)
          {
            if (iparm[IPARM_FILL_MATRIX] == API_NO)
              SOLV_COEFTAB(itercblk)[j] = ZERO;
            else
              SOLV_COEFTAB(itercblk)[j] = UN;
            if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
            {
              if (iparm[IPARM_FILL_MATRIX] == API_NO)
                SOLV_UCOEFTAB(itercblk)[j] = ZERO;
              else
                SOLV_UCOEFTAB(itercblk)[j] = DEUX;
            }
          }
#ifdef _UNUSED_
        }
#endif
    }

#ifdef PASTIX_DYNSCHED
    bubnum = BFATHER(datacode->btree, bubnum);
  }
#endif /* PASTIX_DYNSCHED */

  /* 2 eme phase de l'initialisation de la matrice */
  for (i=0; i < SOLV_TTSKNBR; i++)
  {
    itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(i));
    coefnbr  = SOLV_STRIDE(itercblk) * (SYMB_LCOLNUM(itercblk) - SYMB_FCOLNUM(itercblk) + 1);

    ooc_wait_for_cblk(sopalin_data, itercblk, me);

    /* initialisation du bloc colonne */
    for (j=0 ; j < coefnbr ; j++)
    {
      SOLV_COEFTAB(itercblk)[j] = UN;
      if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
        SOLV_UCOEFTAB(itercblk)[j] = DEUX;
    }

    /* if we are on a diagonal bloc */
    if (SYMB_FCOLNUM(itercblk) == SYMB_FROWNUM(SYMB_BLOKNUM(itercblk)))
    {
      PASTIX_INT index  = SOLV_COEFIND(SYMB_BLOKNUM(itercblk));
      PASTIX_INT size   = SYMB_LCOLNUM(itercblk) - SYMB_FCOLNUM(itercblk) + 1;
      PASTIX_INT stride = SOLV_STRIDE(itercblk);

      for (itercol=0; itercol<size; itercol++)
      {
        /* On s'assure que la matrice est diagonale dominante */
        SOLV_COEFTAB(itercblk)[index+itercol*stride+itercol] = (PASTIX_FLOAT) (UPDOWN_GNODENBR*UPDOWN_GNODENBR);
      }
      /* copie de la partie de block diag de U dans L */
      if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
      {
        PASTIX_INT iterrow;
        for (itercol=0; itercol<size; itercol++)
        {
          for (iterrow=itercol+1; iterrow<size; iterrow++)
          {
            SOLV_COEFTAB(itercblk)[index+iterrow*stride+itercol] = SOLV_UCOEFTAB(itercblk)[index+itercol*stride+iterrow];
          }
        }
      }
    }

    ooc_save_coef(sopalin_data, SOLV_TTSKTAB(i), itercblk, me);
  }
#ifdef _UNUSED_
}
#endif
printf("fin false fill-in\n");
}


if (iparm[IPARM_FREE_CSCPASTIX] == API_CSC_FREE)
 {
   SYNCHRO_X_THREAD(SOLV_THRDNBR, *barrier);

   /* Internal csc is useless if we don't want to do refinement step */
   if (me == 0)
   {
     if ((iparm[IPARM_END_TASK] < API_TASK_SOLVE) ||
         ((iparm[IPARM_END_TASK] < API_TASK_REFINE) &&
          (iparm[IPARM_RHS_MAKING] == API_RHS_B)))
     {
       CscExit(sopalin_data->sopar->cscmtx);
     }
     else
     {
       errorPrintW("The internal CSC can't be freed if you want to use refinement or if you don't give one RHS.\n");
     }
   }
 }
}

/*
 Function: CoefMatrix_Free

 Free the solver matrix coefficient tabular : coeftab and ucoeftab.

 WARNING: Call it with one unnique thread.

 Parameters:
 datacode   - solverMatrix
 factotype  - factorisation type (<API_FACT>)

 */
void CoefMatrix_Free(SopalinParam *sopar,
                     SolverMatrix *datacode,
                     PASTIX_INT           factotype)
{
  PASTIX_INT i;

  if ( (factotype == API_FACT_LU)
       || ( (factotype == API_FACT_LDLT) && sopar->iparm[IPARM_ESP]) )
  {
    for (i=0 ; i < SYMB_CBLKNBR; i++)
      if (SOLV_UCOEFTAB(i) != NULL)
        memFree_null(SOLV_UCOEFTAB(i));
  }
  for (i=0 ; i < SYMB_CBLKNBR; i++)
    if (SOLV_COEFTAB(i) != NULL)
      memFree_null(SOLV_COEFTAB(i));
}
