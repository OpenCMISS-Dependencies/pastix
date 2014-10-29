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
  File: sopalin_compute.c

  Computation functions.

  Pierre Ramet : fev 2003

*/

void API_CALL(factor_diag)  (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT c);
void API_CALL(factor_trsm1d)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT c);
void API_CALL(compute_contrib_compact)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT c, PASTIX_INT b1, PASTIX_INT b2, PASTIX_INT usediag);
void API_CALL(add_contrib_local)      (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT b1,PASTIX_INT b2,PASTIX_INT c,PASTIX_INT b3,PASTIX_INT cbl);
void API_CALL(add_contrib_target)     (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT b1,PASTIX_INT b2,PASTIX_INT task,PASTIX_INT t);

/*
 * Compute tasks
 */
void API_CALL(compute_diag)  (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task);
void API_CALL(compute_1d)    (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task);
void API_CALL(compute_1dgemm)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task, PASTIX_INT i, PASTIX_INT b2);
void API_CALL(compute_e1)    (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task);
void API_CALL(compute_e2)    (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task);
void API_CALL(compute_unlock_after_DiagE1)(Sopalin_Data_t * sopalin_data, PASTIX_INT task);

#include "./compute_gemdm.c"
#include "./compute_diag.c"
#include "./compute_trsm.c"

#if (defined COMM_REORDER) || (defined PASTIX_DYNSCHED)
void API_CALL(compute_unlock_after_DiagE1)(Sopalin_Data_t * sopalin_data, PASTIX_INT task)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  PASTIX_INT            btagnum, sendcnt;
  PASTIX_INT            dest, i;

  /* Add communication in adapted tasks list */
  /* And add local tasks in lists */

  btagnum = SOLV_INDTAB[TASK_INDNUM(task)];
  sendcnt = STASK_SENDCNT(task);

  for(i=0; i<sendcnt; i++, btagnum++)
  {
    dest = BTAG_PROCDST(btagnum);

    if (dest != SOLV_PROCNUM)
  {
#ifdef COMM_REORDER
  if (THREAD_FUNNELED_ON)
    {
      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      queueAdd2(sopalin_data->sendqueue, btagnum,
                (double)(-(dest+1)), (PASTIX_INT)BTAG_PRIONUM(btagnum));
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
    }
  else
    {
      MUTEX_LOCK(&(sopalin_data->mutex_queue_block[dest]));
      queueAdd2(&(sopalin_data->blocktgtsendqueue[dest]), btagnum,
                ((double)BTAG_PRIONUM(btagnum)), btagnum);
      MUTEX_UNLOCK(&(sopalin_data->mutex_queue_block[dest]));
    }
#endif /* COMM_REORDER */
  }
#ifdef PASTIX_DYNSCHED
    else
  {
    PASTIX_INT j;
    PASTIX_INT firsttask = BTAG_TASKDST(btagnum);
    PASTIX_INT localtask = firsttask;

    MUTEX_LOCK(&(sopalin_data->mutex_task[localtask]));
    if ((!TASK_CTRBCNT(localtask))
      && (sopalin_data->taskmark[localtask] == -1))
    {
      j = TASK_THREADID(localtask);

#if (DBG_PASTIX_DYNSCHED > 0)
      ASSERTDBG(sopalin_data->taskmark[localtask] == -1, MOD_SOPALIN);
      ASSERTDBG(TASK_BTAGPTR(localtask) != NULL, MOD_SOPALIN);
      ASSERTDBG(RTASK_COEFTAB(localtask) != NULL, MOD_SOPALIN);
#endif
      sopalin_data->taskmark[localtask]++;

      MUTEX_LOCK(&(sopalin_data->tasktab_mutex[j]));
      queueAdd(&(sopalin_data->taskqueue[j]),
         localtask,
         (double)TASK_PRIONUM(localtask));
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[j]));
      pthread_cond_broadcast(&(sopalin_data->tasktab_cond[j]));
    }
    MUTEX_UNLOCK(&(sopalin_data->mutex_task[localtask]));

    while (TASK_TASKNEXT(localtask) != firsttask)
    {
      localtask = TASK_TASKNEXT(localtask);

      MUTEX_LOCK(&(sopalin_data->mutex_task[localtask]));
      if ((!TASK_CTRBCNT(localtask))
      && (sopalin_data->taskmark[localtask] == -1))
    {
      j = TASK_THREADID(localtask);

#if (DBG_PASTIX_DYNSCHED > 0)
      ASSERTDBG(sopalin_data->taskmark[localtask] == -1, MOD_SOPALIN);
      ASSERTDBG(TASK_BTAGPTR(localtask) != NULL, MOD_SOPALIN);
      ASSERTDBG(RTASK_COEFTAB(localtask) != NULL, MOD_SOPALIN);
#endif
      sopalin_data->taskmark[localtask]++;

      MUTEX_LOCK(&(sopalin_data->tasktab_mutex[j]));
      queueAdd(&(sopalin_data->taskqueue[j]),
         localtask,
         (double)TASK_PRIONUM(localtask));
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[j]));
      pthread_cond_broadcast(&(sopalin_data->tasktab_cond[j]));
    }
      MUTEX_UNLOCK(&(sopalin_data->mutex_task[localtask]));
    }
  }
#endif
  }
}
#endif /* REORDER/BUBBLE */

/****************************************************************************/
/* COMPUTE TASK DIAG                                                        */
/****************************************************************************/
/*
 * Task to factorize the diagonal block
 */
void API_CALL(compute_diag)(Sopalin_Data_t *sopalin_data,
        PASTIX_INT             me,
        PASTIX_INT             task)
{
  PASTIX_INT            c;
  SolverMatrix  *datacode    = sopalin_data->datacode;
#ifdef TRACE_SOPALIN
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif

  c = TASK_CBLKNUM(task);

  /*  For Schur complement, last diagonal block is not factorized */
  if ( (sopalin_data->sopar->schur == API_YES) &&
     (SOLV_INDTAB[TASK_INDNUM(task)] == -1) )
  {
    return;
  }

  trace_begin_task(thread_data->tracefile,
       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
       STATE_L2_DIAG, task);

  /*
   * Compute
   */
  API_CALL(factor_diag)(sopalin_data, me, c);

  /* if not last diag, we copy the factorized diagonal block for local and remote task */
  if (SOLV_INDTAB[TASK_INDNUM(task)] != -1)
  {
    PASTIX_INT N      = SYMB_LCOLNUM(c) - SYMB_FCOLNUM(c) + 1;
    PASTIX_INT offset = SOLV_COEFIND(SYMB_BLOKNUM(c));
    PASTIX_INT size   = N * N;
    PASTIX_INT stride = SOLV_STRIDE(c);
    PASTIX_FLOAT *ga = NULL;
    PASTIX_FLOAT *gc = NULL;
#ifdef SOPALIN_LU
    PASTIX_FLOAT *gc2;

    MALLOC_INTERN(gc, 2*size, PASTIX_FLOAT);
    STATS_ADD( 2*size );
    gc2 = gc + size;
#else
    MALLOC_INTERN(gc,   size, PASTIX_FLOAT);
    STATS_ADD( size );
#endif
    print_debug(DBG_SOPALIN_ALLOC,
          "alloc block coeff %x\n",
          (unsigned int)(intptr_t)gc);

    /* OOC coeftab loaded and locked above */
    ga = &(SOLV_COEFTAB(c)[offset]);
    SOPALIN_LACPY(N, N, ga, stride, gc, N);

#ifdef SOPALIN_LU
    ga = &(SOLV_UCOEFTAB(c)[offset]);
    SOPALIN_LACPY(N, N, ga, stride, gc2, N);
#endif
    print_debug(DBG_SOPALIN_DIAG, "STASK_BCOFPTR : %x\n",
      (unsigned int)(intptr_t)STASK_BCOFPTR(task));

    MUTEX_LOCK(&(sopalin_data->mutex_task[task]));
    STASK_COEFTAB(task) = gc;
    MUTEX_UNLOCK(&(sopalin_data->mutex_task[task]));

    /* Signal for E1 tasks */
    pthread_cond_broadcast(&(sopalin_data->cond_task[task]));

    print_debug(DBG_SOPALIN_DIAG, "STASK_COEFTAB : %x\n",
      (unsigned int)(intptr_t)STASK_COEFTAB(task));
    print_debug(DBG_SOPALIN_DIAG, "btagtab[%ld]\n",
      (long)SOLV_INDTAB[TASK_INDNUM(task)]);

#if (defined COMM_REORDER) || (defined PASTIX_DYNSCHED)
    /* Add communication in adapted tasks list */
    /* And add local tasks in lists */
    API_CALL(compute_unlock_after_DiagE1)(sopalin_data, task);
#endif

  }

  /* Send what needs to be sent */
  if (THREAD_FUNNELED_OFF)
    {
      API_CALL(send_all_block)(sopalin_data, me);
    }
  else
    {
      pthread_cond_broadcast(&(sopalin_data->cond_comm));
    }

  trace_end_task(thread_data->tracefile,
     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
     STATE_L2_DIAG, task);
}

/****************************************************************************/
/* COMPUTE TASK 1D                                                          */
/****************************************************************************/

/*
 * Compute the update in one big block for all the left block column
 * multiply by the one in regards with the diagonal block
 * The result is stored in a temporary buffer to be added part
 * by part to the target cblk
 */
void API_CALL(compute_contrib_compact)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT c, PASTIX_INT b1, PASTIX_INT b2, PASTIX_INT usediag)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  PASTIX_FLOAT         *gaik, *gb, *gc;
#ifdef CHOL_SOPALIN
  PASTIX_FLOAT         *gajk;
#endif
  PASTIX_INT            dima, dimi, dimj, stride, k;
  (void)b2; (void)gb; (void)usediag;

  gb = thread_data->maxbloktab1; /* C in U for LU, B for LDLt */
  gc = thread_data->maxbloktab2; /* C in L */
  stride = SOLV_STRIDE(c);

  /* Matrix A = Aik */
  gaik = &(SOLV_COEFTAB(c)[SOLV_COEFIND(b1)]);

  /* Compute M */
  dimi = 0;
  for (k=b1; k<SYMB_BLOKNUM(c+1); k++)
  dimi += SYMB_LROWNUM(k) - SYMB_FROWNUM(k) + 1;

  /* Compute N */
  dimj = SYMB_LROWNUM(b1) - SYMB_FROWNUM(b1) + 1;/* ATTENTION stride-dima;*/
#ifdef PASTIX_ESC
  for (k=b2; k>b1; k--)
    dimj += SYMB_LROWNUM(k) - SYMB_FROWNUM(k) + 1;
#endif

  /* Compute K */
  dima = SYMB_LCOLNUM(c) - SYMB_FCOLNUM(c) + 1;

  thread_data->firstbloktab  = b1;
  thread_data->stridebloktab = dimi;

  ASSERTDBG(dimj*thread_data->stridebloktab <= SOLV_COEFMAX, MOD_SOPALIN);

  /* multAikB(gc,gaik,gb,stride,dimj,dima,dimi) */
#ifdef CHOL_SOPALIN
#ifdef SOPALIN_LU
  /* L update */
  gajk = &(SOLV_UCOEFTAB(c)[SOLV_COEFIND(b1)]);
  SOPALIN_GEMM("N", "T", dimi, dimj, dima,
         fun,   gaik, stride,
            gajk, stride,
       fzero, gc,   dimi);

  /* U update */
  gaik = &(SOLV_UCOEFTAB(c)[SOLV_COEFIND(b1)]);
  gajk = &(SOLV_COEFTAB(c)[ SOLV_COEFIND(b1)]);
  SOPALIN_GEMM("N", "T", dimi, dimj, dima,
         fun,   gaik, stride,
            gajk, stride,
       fzero, gb,   dimi);
#else /* SOPALIN_LU */
  gajk = &(SOLV_COEFTAB(c)[ SOLV_COEFIND(b1)]);
  SOPALIN_GEMM("N",
               "C",
               dimi, dimj, dima,
               fun,   gaik, stride,
               gajk, stride,
               fzero, gc,   dimi);
#endif /* SOPALIN_LU */
#else /* CHOL_SOPALIN */
  if ( usediag == 1 )
  {
  int ldw = SOLV_COEFMAX;
  PASTIX_FLOAT *D = &(SOLV_UCOEFTAB(c)[SOLV_COEFIND(SYMB_BLOKNUM(c))]);
  gb = &(SOLV_COEFTAB(c)[SOLV_COEFIND(b1)]);

  API_CALL(CORE_gemdm)(PastixNoTrans,
#ifdef HERMITIAN
                       PastixConjTrans,
#else
                       PastixTrans,
#endif
                       dimi, dimj, dima,
                       fun,   gaik, stride,
                       gb,   stride,
                       fzero, gc,   dimi,
                       D, iun,
                       thread_data->maxbloktab1, ldw );
  }
  else
  {
  gb = thread_data->maxbloktab1+(stride-dima-dimi); /* Correction pour le pere Goudin */
#ifdef HERMITIAN
  SOPALIN_GEMM("N",
               "C",
               dimi, dimj, dima,
               fun,   gaik, stride,
               gb,   stride-dima,
               fzero, gc,   dimi);
#else
  SOPALIN_GEMM("N",
               "T",
               dimi, dimj, dima,
               fun,   gaik, stride,
               gb,   stride-dima,
               fzero, gc,   dimi);
#endif
  }
#endif /* CHOL_SOPALIN */
}

/*
   Function: API_CALL(add_contrib_local)

   Add contribution from b2 on b3.

   Parameters:
   sopalin_data - Global sopalin structure
   me           - Computing thread number
   b1           - extra digonal block of c
   b2           - extra diagonal block after b1 on c
   c            - Column block
   b3           - Extra digonal block on cbl facing b2
   cbl          - Column block facing b1

 */
void API_CALL(add_contrib_local)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT b1, PASTIX_INT b2, PASTIX_INT c, PASTIX_INT b3, PASTIX_INT cbl)
{
  PASTIX_INT frownum, lrownum, ofrownum, olrownum;
  PASTIX_INT dimi, dimj, stridea, strideb, step, k;
  PASTIX_FLOAT *ga,*gb;
#ifdef SOPALIN_LU
  PASTIX_FLOAT *ga2, *gb2;
#endif /* SOPALIN_LU */
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#ifdef DEBUG_SOPALIN_NAPA
  int flag;
#endif
  (void)c;

  /* column block in front of b1 */
  /* cbl=SYMB_CBLKNUM(b1); ??? */
  stridea = SOLV_STRIDE(cbl);

  print_debug(DBG_SOPALIN_COMP1D,
              "FR(b2) %ld FR(b3) %ld\nLR(b3) %ld LR(b2) %ld\nFR(b1) %ld FC(cb) %ld\nLC(cb) %ld LR(b1) %ld\n",
              (long)SYMB_FROWNUM(b2), (long)SYMB_FROWNUM(b3),
              (long)SYMB_LROWNUM(b3), (long)SYMB_LROWNUM(b2),
              (long)SYMB_FROWNUM(b1), (long)SYMB_FCOLNUM(cbl),
              (long)SYMB_LCOLNUM(cbl),(long)SYMB_LROWNUM(b1));

#ifdef NAPA_SOPALIN
  ASSERTDBG((SYMB_FROWNUM(b1)  >= SYMB_FCOLNUM(cbl)) &&
            (SYMB_LCOLNUM(cbl) >= SYMB_LROWNUM(b1) ), MOD_SOPALIN);
#else
  ASSERTDBG((SYMB_FROWNUM(b2)  >= SYMB_FROWNUM(b3) ) &&
            (SYMB_LROWNUM(b3)  >= SYMB_LROWNUM(b2) ) &&
            (SYMB_FROWNUM(b1)  >= SYMB_FCOLNUM(cbl)) &&
            (SYMB_LCOLNUM(cbl) >= SYMB_LROWNUM(b1) ), MOD_SOPALIN);
#endif

  ga = &(SOLV_COEFTAB(cbl)[ SOLV_COEFIND(b3)+
                            (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                            (SYMB_FROWNUM(b2)-SYMB_FROWNUM(b3))]);
#ifdef SOPALIN_LU
  if ((b3 == SYMB_BLOKNUM(cbl)) && (b1 != b2)) {
    ga2 = &(SOLV_COEFTAB(cbl)[SOLV_COEFIND(b3)+
                              (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))+
                              (SYMB_FROWNUM(b2)-SYMB_FROWNUM(b3))*stridea]);
  }
  else {
    ga2 = &(SOLV_UCOEFTAB(cbl)[SOLV_COEFIND(b3)+
                               (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                               (SYMB_FROWNUM(b2)-SYMB_FROWNUM(b3))]);
  }
#endif /* SOPALIN_LU */

  /* vertical dimension */
  dimj = SYMB_LROWNUM(b1) - SYMB_FROWNUM(b1) + 1;
  /* vertical dimension */
  dimi = SYMB_LROWNUM(b2) - SYMB_FROWNUM(b2) + 1;

  strideb = thread_data->stridebloktab;

  step = b2 - b1;
  for (k=b1; k<b2; k++)
    step += SYMB_LROWNUM(k) - SYMB_FROWNUM(k);

  ASSERTDBG(step+dimi    <= strideb,     MOD_SOPALIN);
  ASSERTDBG(strideb*dimj <= SOLV_COEFMAX,MOD_SOPALIN);

#ifdef PASTIX_ESC
  /* Décalage du bloc colonne de k + decalage du bloc k inutile */
  for (k=thread_data->firstbloktab; k<b1; k++)
    step += (SYMB_LROWNUM(k) - SYMB_FROWNUM(k) + 1) * (strideb + 1);
#endif

  gb  = &(thread_data->maxbloktab2[step]);
#ifdef SOPALIN_LU
  gb2 = &(thread_data->maxbloktab1[step]);
#endif /* SOPALIN_LU */

#ifdef NAPA_SOPALIN
  ofrownum=SYMB_FROWNUM(b2);
  olrownum=SYMB_LROWNUM(b2);
  b3--;
#ifdef DEBUG_SOPALIN_NAPA
  flag = 1;
#endif
  do {
#ifdef DEBUG_SOPALIN_NAPA
    PASTIX_INT trace = 0;
    /* il peut y avoir plusieurs cibles partielles */
    if (!flag)
      {
        print_debug(DBG_SOPALIN_NAPA, "ILU: plusieurs cibles locales\n");
      }
#endif
    frownum=ofrownum;
    lrownum=olrownum;
    b3++;
#ifdef DEBUG_SOPALIN_NAPA
    if ((!flag) || (SYMB_FROWNUM(b3)>frownum) ||
        (SYMB_LROWNUM(b3)<lrownum))
      {
        trace = 1;
        if (flag)
          {
            print_debug(DBG_SOPALIN_NAPA, "\nILU: debug local SF=%ld F=%ld SL=%ld L=%ld gb=%ld ga=%ld\n",
                        (long)SYMB_FROWNUM(b3), (long)frownum,
                        (long)SYMB_LROWNUM(b3), (long)lrownum,
                        (long)0, (long)(SOLV_COEFIND(b3)+
                                        (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                                        (SYMB_FROWNUM(b2)-SYMB_FROWNUM(b3))));
          }
      }
#endif
    if (SYMB_FROWNUM(b3)>frownum)
      {
        frownum=SYMB_FROWNUM(b3);
        print_debug(DBG_SOPALIN_NAPA, "ILU: tronque frownum\n");
      }
    if (SYMB_LROWNUM(b3)<lrownum)
      {
        lrownum=SYMB_LROWNUM(b3);
        print_debug(DBG_SOPALIN_NAPA, "ILU: tronque lrownum\n");
      }
    dimi=lrownum-frownum+1;
    /*
     gb=&(maxbloktab2[me][step])+frownum-ofrownum;
     */
    gb=&(thread_data->maxbloktab2[step+frownum-ofrownum]);
    ga=&(SOLV_COEFTAB(cbl)[SOLV_COEFIND(b3)+
                           (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                           (frownum-SYMB_FROWNUM(b3))]);
#ifdef SOPALIN_LU
    if ( b3!=SYMB_BLOKNUM(cbl) )
      {
        ga2=&(SOLV_UCOEFTAB(cbl)[SOLV_COEFIND(b3)+
                                 (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                                 (frownum-SYMB_FROWNUM(b3))]);
      }
    else if ( b1 != b2 ) /* Store U and L directly in L for factorization */
      {
        ga2=&(SOLV_COEFTAB(cbl)[SOLV_COEFIND(b3)+
                                (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))+
                                (frownum-SYMB_FROWNUM(b3))*stridea]);
      }
    else {
      /* We are on diagonal block => we don't do anything */
    }
    gb2=&(thread_data->maxbloktab1[step+frownum-ofrownum]);
#endif /* SOPALIN_LU */
#ifdef DEBUG_SOPALIN_NAPA
    if (trace)
      {
        print_debug(DBG_SOPALIN_NAPA,
                    "ILU: debug local SF=%ld F=%ld SL=%ld L=%ld gb=%ld ga=%ld\n",
                    (long)SYMB_FROWNUM(b3), (long)frownum,
                    (long)SYMB_LROWNUM(b3), (long)lrownum,
                    (long)(frownum-ofrownum),
                    (long)(SOLV_COEFIND(b3) +
                           (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                           (frownum-SYMB_FROWNUM(b3))));
      }
#endif

    ASSERTDBG((SYMB_FROWNUM(b3)<=frownum) &&
              (SYMB_LROWNUM(b3)>=lrownum),MOD_SOPALIN);

#else /* NAPA_SOPALIN */

    ASSERTDBG((SYMB_FROWNUM(b3)<=SYMB_FROWNUM(b2)) &&
              (SYMB_LROWNUM(b3)>=SYMB_LROWNUM(b2)),MOD_SOPALIN);

#endif /* NAPA_SOPALIN */

    MUTEX_LOCK(&(sopalin_data->mutex_blok[b3]));
    SOPALIN_GESM( "N", "N", dimi, dimj,
                  fun, gb, strideb,
                  ga, stridea);
#ifdef SOPALIN_LU
    if ( b3!=SYMB_BLOKNUM(cbl) )
      {
        SOPALIN_GESM("N","N",dimi,dimj,fun,gb2,strideb,ga2,stridea);
      }
    else if ( b1 != b2 ) /* Store U and L directly in L for factorization */
      {
        SOPALIN_GESM("T","N",dimj,dimi,fun,gb2,strideb,ga2,stridea);
      }
    else {
      /* We are on diagonal block => we don't do anything */
    }
#endif
    MUTEX_UNLOCK(&(sopalin_data->mutex_blok[b3]));

#ifdef NAPA_SOPALIN
#ifdef DEBUG_SOPALIN_NAPA
    flag = 0;
#endif
  } while  ((b3+1<SYMB_BLOKNUM(cbl+1)) &&
            (
              ((SYMB_FROWNUM(b3+1)<=SYMB_FROWNUM(b2)) &&
               (SYMB_LROWNUM(b3+1)>=SYMB_FROWNUM(b2))) ||
              ((SYMB_LROWNUM(b3+1)>=SYMB_LROWNUM(b2)) &&
               (SYMB_FROWNUM(b3+1)<=SYMB_LROWNUM(b2))) ||
              ((SYMB_FROWNUM(b3+1)>=SYMB_FROWNUM(b2)) &&
               (SYMB_LROWNUM(b3+1)<=SYMB_LROWNUM(b2))) ||
              ((SYMB_FROWNUM(b3+1)<=SYMB_FROWNUM(b2)) &&
               (SYMB_LROWNUM(b3+1)>=SYMB_LROWNUM(b2)))
             ));
#endif
}

void API_CALL(add_contrib_target)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT b1,PASTIX_INT b2, PASTIX_INT task, PASTIX_INT t)
{
  PASTIX_INT dimi,dimj,stridea,strideb,step,k;
  PASTIX_FLOAT *ga,*gb;
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  (void)task;

  stridea = FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1;

  /* Allocation du buffer pour la contribution a envoyer */
#ifdef ALLOC_FTGT
  MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
#ifdef OOC_FTGT
  print_debug(DBG_OOC_FTGT, "WAIT %4d %4d\n", (int)t, (int) task);
  ooc_wait_for_ftgt(sopalin_data, t, me);

  ASSERTDBG(((unsigned long)(*(((double*)FANIN_COEFTAB(t))-1))) ==
            sizeof(PASTIX_FLOAT) * ((FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1) * (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1)   *
                                    ((sopalin_data->sopar->factotype == API_FACT_LU)?2:1))
            , MOD_SOPALIN);
#else
  if (FANIN_COEFTAB(t)==NULL)
    {
      PASTIX_INT j;
      PASTIX_INT ftgtsize = (FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1)
        *(FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);
#ifdef SOPALIN_LU
      ftgtsize *= 2;
#endif

      MALLOC_INTERN(FANIN_COEFTAB(t), ftgtsize, PASTIX_FLOAT);
      for (j=0;j<ftgtsize;j++)
        FANIN_COEFTAB(t)[j] = fzero;

      print_debug(DBG_SOPALIN_ALLOC, "alloc fanin coeff %x\n",(unsigned int)(intptr_t)FANIN_COEFTAB(t));

      STATS_ADD( ftgtsize );
    }
#endif /* OOC_FTGT */
  MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
#endif /* ALLOC_FTGT */

  print_debug(DBG_SOPALIN_COMP1D,
              "FR(b2) %ld FR(tg) %ld\nLR(tg) %ld LR(b2) %ld\nFR(b1) %ld FC(tg) %ld\nLC(tg) %ld LR(b1) %ld\n",
              (long)SYMB_FROWNUM(b2), (long)FANIN_FROWNUM(t),
              (long)FANIN_LROWNUM(t), (long)SYMB_LROWNUM(b2),
              (long)SYMB_FROWNUM(b1), (long)FANIN_FCOLNUM(t),
              (long)FANIN_LCOLNUM(t), (long)SYMB_LROWNUM(b1));

  ASSERTDBG((SYMB_FROWNUM(b2)>=FANIN_FROWNUM(t))&&
            (FANIN_LROWNUM(t)>=SYMB_LROWNUM(b2))&&
            (SYMB_FROWNUM(b1)>=FANIN_FCOLNUM(t))&&
            (FANIN_LCOLNUM(t)>=SYMB_LROWNUM(b1)),MOD_SOPALIN);

  ga=&(FANIN_COEFTAB(t)[(SYMB_FROWNUM(b1)-FANIN_FCOLNUM(t))*stridea+
                        (SYMB_FROWNUM(b2)-FANIN_FROWNUM(t))]);

  /* vertical dimension */
  dimj=SYMB_LROWNUM(b1)-SYMB_FROWNUM(b1)+1;
  /* vertical dimension */
  dimi=SYMB_LROWNUM(b2)-SYMB_FROWNUM(b2)+1;

  strideb=thread_data->stridebloktab;

  /* Indice de debut de bloc dans maxbloktab */
  step = b2 - b1;
  for (k=b1; k<b2; k++)
    step += SYMB_LROWNUM(k) - SYMB_FROWNUM(k);

  ASSERTDBG(step+dimi<=strideb,MOD_SOPALIN);
  ASSERTDBG(strideb*dimj<=SOLV_COEFMAX,MOD_SOPALIN);

#ifdef PASTIX_ESC
  /* Décalage du bloc colonne de k + decalage du bloc k inutile */
  for (k=thread_data->firstbloktab; k<b1; k++)
    step += (SYMB_LROWNUM(k) - SYMB_FROWNUM(k) + 1) * (strideb + 1);
#endif

  gb=&(thread_data->maxbloktab2[step]);

  MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));

  SOPALIN_GESM("N","N",dimi,dimj,fun,gb,strideb,ga,stridea);

#ifdef SOPALIN_LU
  if ( (b1!=b2) &&
       (SYMB_FROWNUM(b2)>= FANIN_FCOLNUM(t)) &&
       (SYMB_LROWNUM(b2)<= FANIN_LCOLNUM(t)))
    {
      ga=&(FANIN_COEFTAB(t)[(SYMB_FROWNUM(b1)-FANIN_FCOLNUM(t))+
                            (SYMB_FROWNUM(b2)-FANIN_FROWNUM(t))*stridea]);
      gb=&(thread_data->maxbloktab1[step]);
      SOPALIN_GESM("T","N",dimj,dimi,fun,gb,strideb,ga,stridea);
    }
  else
    {
      /* ga += ftgtsize/2 */
      ga+=(FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1)
        *(FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);
      gb=&(thread_data->maxbloktab1[step]);
      SOPALIN_GESM("N","N",dimi,dimj,fun,gb,strideb,ga,stridea);
    }
#endif

  /* WARNING : sauvegarder que si necessaire */
  ooc_save_ftgt(sopalin_data, task, t, me);

  FANIN_CTRBCNT(t)--;
  /* Add fanin to the correct ready heap */
#ifdef COMM_REORDER
  if (FANIN_CTRBCNT(t) == 0)
    {
      PASTIX_INT dest = FANIN_PROCDST(t);
      MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
      if (THREAD_FUNNELED_ON)
        {
          MUTEX_LOCK(&(sopalin_data->mutex_comm));
          queueAdd2(sopalin_data->sendqueue, t,
                    (double)(dest+1), (PASTIX_INT)FANIN_PRIONUM(t));
          MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
        }
      else
        {
          MUTEX_LOCK(&(sopalin_data->mutex_queue_fanin[dest]));
          queueAdd2(&(sopalin_data->fanintgtsendqueue[dest]), t,
                    ((double)FANIN_PRIONUM(t)), t);
          MUTEX_UNLOCK(&(sopalin_data->mutex_queue_fanin[dest]));
        }
    }
  else
#endif
    MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
}


/*
   Function: API_CALL(compute_1d)

   Factorisation of one block column

  Parameters:
   sopalin_data -
   me           - Thread number
   task         - Task number

*/
void API_CALL(compute_1d)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
#if defined(TRACE_SOPALIN) || defined(PASTIX_DYNSCHED)
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif
#ifdef PASTIX_DYNSCHED
  int            esp         = sopalin_data->sopar->iparm[IPARM_ESP];
  int            esparam     = sopalin_data->sopar->iparm[IPARM_ESP_THRESHOLD];
  PASTIX_INT     t;
#endif
  PASTIX_INT            c, fblknum, lblknum;
  PASTIX_INT            i, ii, jj, n;
  PASTIX_INT            dimb, dimb2;
#ifdef OOC
  PASTIX_INT            tooc;
#endif

  c = TASK_CBLKNUM(task);

  if (sopalin_data->sopar->schur == API_YES)
  {
    PASTIX_INT lN = sopalin_data->sopar->gN * sopalin_data->sopar->iparm[IPARM_DOF_NBR] - 1;
    if ( SYMB_LCOLNUM(c) == lN )
      return;
  }

  trace_begin_task(thread_data->tracefile,
       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
       STATE_L2_COMP1D, task);

  API_CALL(factor_diag)(sopalin_data, me, c);

  fblknum = SYMB_BLOKNUM(c);
  lblknum = SYMB_BLOKNUM(c + 1);

  /* Maximal M dimension for the GEMMs */
  dimb = SOLV_STRIDE(c) - ( SYMB_LROWNUM(fblknum) - SYMB_FROWNUM(fblknum) + 1 );
  fblknum++;

  /* if there is an extra-diagonal bloc in column block */
  if ( fblknum < lblknum )
  {
    API_CALL(factor_trsm1d)(sopalin_data, me, c);
  }

  /* for all extra-diagonal column blocks */
  n = 0;
  for (i=fblknum; i<lblknum; )
  {
    ii = 1;

    /* N dimension of the GEMM */
    dimb2 = SYMB_LROWNUM(i) - SYMB_FROWNUM(i) + 1;
#ifdef PASTIX_ESC
    while( ( (i+ii) < lblknum )
           && ( SYMB_CBLKNUM(i) == SYMB_CBLKNUM(i+ii) ) ) {
      dimb2 += SYMB_LROWNUM(i+ii) - SYMB_FROWNUM(i+ii) + 1;
      ii++;
    }
#endif

#ifdef PASTIX_DYNSCHED
    t  = SOLV_INDTAB[TASK_INDNUM(task)+n];
#ifdef ESP_A
    if (esp && (dimb2*dimb2 > esparam) )
#else
    if (esp && (dimb2*dimb  > esparam) )
#endif
    {
      PASTIX_INT prionum;
#ifdef ESP_WRITE
      PASTIX_INT tid = TASK_THREADID(-t);
#else
      PASTIX_INT tid = TASK_THREADID(task);
#endif
      if (t > 0) {
        prionum = (double)(FANIN_PRIONUM(t));
      }
      else {
        prionum = (double)(TASK_PRIONUM(-t));
      }

      /*
       * We take the average priority of the current task and the
       * target task to avoid to delay all the tasks at the same
       * moment, creating conflicts to acces the mutex
       */
      prionum = (PASTIX_INT)floor( (( 2.* prionum + TASK_PRIONUM(task)) / 3.) - 1. );

      MUTEX_LOCK(&(sopalin_data->tasktab_mutex[tid]));
      queueAdd2(&(sopalin_data->taskqueue[tid]), TASK_TASK2ESP(task), prionum, i);
      sopalin_data->tasktab_indice[tid]--;
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[tid]));
      pthread_cond_broadcast(&(sopalin_data->tasktab_cond[tid]));

      thread_data->esp++;

      print_debug(DBG_SOPALIN_COMP1D,
                  "COMP1D %05d Ajout bloc %05d / taskdst %05d / priorite %05d / me %d / n %03d\n",
                  (int)task, (int)i, (int)-t, (int)prionum, (int)tid, (int)n);
    }
    else
#endif /* PASTIX_DYNSCHED */
    {
      API_CALL(compute_1dgemm)(sopalin_data, me, task, i, i+ii);
    }

    dimb -= dimb2;
    for (jj=0; jj<ii; jj++,i++)
      n += (lblknum - i);
  }

  trace_end_task(thread_data->tracefile,
       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
       STATE_L2_COMP1D, task);
}

void API_CALL(compute_1dgemm)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task, PASTIX_INT i, PASTIX_INT b2)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
#ifdef TRACE_SOPALIN
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif
  PASTIX_INT            j, t, fblknum, lblknum, n, usediag = 0;
  PASTIX_INT            c           = TASK_CBLKNUM(task);
#ifdef OOC
  PASTIX_INT            tooc;
#endif

  fblknum = SYMB_BLOKNUM(c); /* Be careful, not the same fblknum than in comp1d[plus] */
  lblknum = SYMB_BLOKNUM(c + 1);

  n = i - fblknum;
  n = (n * (n - 1)) / 2;
  n = (i - fblknum - 1) * (lblknum - fblknum) - n;

  trace_begin_task(thread_data->tracefile,
       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
       STATE_L2_COMP1DGEMM, task);

  /* Si on arrive par une tache qui n'a pas l'info on la recalcule */
  if (b2 == -1)
  {
    b2 = i+1;
#ifdef PASTIX_ESC
    while( (b2 < lblknum)
           && (SYMB_CBLKNUM(i) == SYMB_CBLKNUM(b2)))
      b2++;
#endif
    usediag = 1;
  }

  /* Compute contributions (GEMM) */
  API_CALL(compute_contrib_compact)(sopalin_data, me, c, i, b2-1, usediag);

#ifdef OOC
  if ((tooc = SOLV_INDTAB[TASK_INDNUM(task)+(n)]) < 0)
  {
    ooc_wait_for_cblk(sopalin_data, TASK_CBLKNUM(-tooc), me);
  }
#endif /* OOC */

  /* for all following blocks in block column */
  for (; i<b2; i++)
  for (j=i; j<lblknum; j++)
  {
    t = SOLV_INDTAB[TASK_INDNUM(task)+(n++)];
    if (t < 0)
  {
    PASTIX_INT b;
    /* if the contribution is local */
#ifdef NAPA_SOPALIN
    /* if the contribution really exist */
    ASSERT(-t<SOLV_TASKNBR,MOD_SOPALIN); /* pas possible dans blend */
#endif
    b = TASK_BLOKNUM(-t);
    if (TASK_TASKID(-t)==COMP_1D)
    {
      b = SYMB_BLOKNUM(TASK_CBLKNUM(-t));
      /* if TASK 1D -> look for bloknum */
#ifdef NAPA_SOPALIN
      while (!(((SYMB_FROWNUM(j)>=SYMB_FROWNUM(b)) &&
      (SYMB_LROWNUM(j)<=SYMB_LROWNUM(b))) ||
         ((SYMB_FROWNUM(j)<=SYMB_FROWNUM(b)) &&
      (SYMB_LROWNUM(j)>=SYMB_LROWNUM(b))) ||
         ((SYMB_FROWNUM(j)<=SYMB_FROWNUM(b)) &&
      (SYMB_LROWNUM(j)>=SYMB_FROWNUM(b))) ||
         ((SYMB_FROWNUM(j)<=SYMB_LROWNUM(b)) &&
      (SYMB_LROWNUM(j)>=SYMB_LROWNUM(b)))))
#else
    while (!((SYMB_FROWNUM(j)>=SYMB_FROWNUM(b)) &&
       (SYMB_LROWNUM(j)<=SYMB_LROWNUM(b))))
#endif
      {
        b++;
        ASSERTDBG(b<SYMB_BLOKNUM(TASK_CBLKNUM(-t)+1),MOD_SOPALIN);
      }
    }

    print_debug(DBG_SOPALIN_COMPUTE,
        "%ld add local contrib %ld %ld %ld %ld (%ld %ld) %ld\n",
        (long)me, (long)i, (long)j, (long)c, (long)b, (long)-t,
        (long)TASK_TASKID(-t), (long)TASK_CTRBCNT(-t));

    ooc_wait_for_cblk(sopalin_data, TASK_CBLKNUM(-t), me);

    API_CALL(add_contrib_local)(sopalin_data, me, i, j, c, b, TASK_CBLKNUM(-t));

    ooc_save_coef(sopalin_data, task, TASK_CBLKNUM(-t), me);

    MUTEX_LOCK(&(sopalin_data->mutex_task[-t]));
    TASK_CTRBCNT(-t)--;
    ASSERTDBG((TASK_CTRBCNT(-t) >= 0), MOD_SOPALIN);
#ifdef PASTIX_DYNSCHED
    if (((!TASK_CTRBCNT(-t)) && (sopalin_data->taskmark[-t] == -1)) &&
      (!((TASK_TASKID(-t) == E1) && ((TASK_BTAGPTR(-t) == NULL) || (RTASK_COEFTAB(-t) == NULL)))))
    {
      PASTIX_INT iter;

#if (DBG_PASTIX_DYNSCHED > 0)
      ASSERTDBG(sopalin_data->taskmark[-t] == -1, MOD_SOPALIN);
#endif
      sopalin_data->taskmark[-t]++;
      MUTEX_UNLOCK(&(sopalin_data->mutex_task[-t]));

      iter = TASK_THREADID(-t);
      MUTEX_LOCK(&(sopalin_data->tasktab_mutex[iter]));
      queueAdd(&(sopalin_data->taskqueue[iter]),
         -t, (double)(TASK_PRIONUM(-t)));
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[iter]));
      pthread_cond_broadcast(&(sopalin_data->tasktab_cond[iter]));
    }
    else
    MUTEX_UNLOCK(&(sopalin_data->mutex_task[-t]));
#else
    MUTEX_UNLOCK(&(sopalin_data->mutex_task[-t]));
    pthread_cond_broadcast(&(sopalin_data->cond_task[-t]));
#endif
  }
    else
  {
    /* if the contribution is not local */
#ifdef NAPA_SOPALIN
    /* if the contribution really exist */
    if (t<SOLV_FTGTNBR) {
#endif

    print_debug(DBG_SOPALIN_COMPUTE,
      "%ld add fanin contrib %ld %ld %ld %ld %ld\n",
      (long)me, (long)i, (long)j, (long)c, (long)t,
      (long)FANIN_CTRBCNT(t));

    API_CALL(add_contrib_target)(sopalin_data, me, i, j, task, t);

#ifdef NAPA_SOPALIN
    }
    else
    {
      print_debug(DBG_SOPALIN_NAPA, "ILU: drop (c=%ld,b=%ld)\n", (long)i, (long)j);
    }
#endif
  }
  }

#ifdef OOC
  if (tooc < 0)
  {
    ooc_save_coef(sopalin_data, task, TASK_CBLKNUM(-tooc), me);
  }
#endif

  trace_end_task(thread_data->tracefile,
     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
     STATE_L2_COMP1DGEMM, task);

  if (THREAD_FUNNELED_OFF)
    {
      PASTIX_INT dest;
      for (dest=0;dest<SOLV_PROCNBR;dest++)
        {
          if (dest == SOLV_PROCNUM) continue;
          API_CALL(send_all_fanin)(sopalin_data, me, dest);
        }
    }
}

/****************************************************************************/
/* COMPUTE CLEAN buffer for E1/E2                                           */
/****************************************************************************/

static inline
void API_CALL(compute_cleane1e2)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task)
{
  SolverMatrix *datacode = sopalin_data->datacode;
  (void)me;
#ifdef STATS_SOPALIN
  PASTIX_INT size = (RTASK_LROWNUM(task) - RTASK_FROWNUM(task) + 1)
       * (RTASK_LCOLNUM(task) - RTASK_FCOLNUM(task) + 1);
#ifdef SOPALIN_LU
  size *= 2;
#endif
#endif

  MUTEX_LOCK(&(sopalin_data->mutex_task[TASK_MASTER(task)]));
  RTASK_TASKCNT(task)--;
  ASSERTDBG(RTASK_TASKCNT(task)>=0,MOD_SOPALIN);

  /*
   * If it is the last E1 task
   *    -> Need to free the copy of the diagonal block from DIAG task
   * If it is the last E2 task
   *    -> Need to free the copy of the extra-diagonal block from E1 task
   */
  if (RTASK_TASKCNT(task)==0)
  {
    if (RTASK_SENDCNT(task) != LOCAL_ALLOC_BTAG) /* not marked btag */
  {
    /* MUTEX_LOCK(&(sopalin_data->mutex_task[task])); */
    RTASK_SENDCNT(task)--;
    /* MUTEX_UNLOCK(&(sopalin_data->mutex_task[task])); */
    ASSERTDBG(RTASK_SENDCNT(task)>=0,MOD_SOPALIN);

    if (RTASK_SENDCNT(task)==0)
    {
      memFree_null(RTASK_COEFTAB(task));
      print_debug(DBG_SOPALIN_ALLOC, "free block coeff %x\n",
              (unsigned int)(intptr_t)RTASK_COEFTAB(task));
        STATS_SUB( size );
    }
  }
    else
  {
    memFree_null(RTASK_COEFTAB(task));
    print_debug(DBG_SOPALIN_ALLOC, "free block coeff %x\n",
            (unsigned int)(intptr_t)RTASK_COEFTAB(task));

      STATS_SUB( size );

    memFree_null(RTASK_BCOFPTR(task));
    print_debug(DBG_SOPALIN_ALLOC, "free block bcof %x\n",
        (unsigned int)(intptr_t)RTASK_BCOFPTR(task));

    memFree_null(TASK_BTAGPTR(task));
    print_debug(DBG_SOPALIN_ALLOC, "free block btag %x\n",
        (unsigned int)(intptr_t)TASK_BTAGPTR(task));

  }
  }
  MUTEX_UNLOCK(&(sopalin_data->mutex_task[TASK_MASTER(task)]));
}

/****************************************************************************/
/* COMPUTE TASK E1                                                          */
/****************************************************************************/

void API_CALL(compute_e1)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task)
{
#if defined(FLAG_ASSERT) || defined(PASTIX_DEBUG) || defined(TRACE_SOPALIN)
  SolverMatrix  *datacode    = sopalin_data->datacode;
#endif
#ifdef TRACE_SOPALIN
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif

  /*
   * Check parameters
   */
  print_debug(DBG_SOPALIN_E1, "RTASK_BCOFPTR :%x\n", (unsigned int)(intptr_t)RTASK_BCOFPTR(task));
  ASSERTDBG(TASK_BTAGPTR(task),  MOD_SOPALIN);
  ASSERTDBG(RTASK_COEFTAB(task), MOD_SOPALIN);

  trace_begin_task(thread_data->tracefile,
       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
       STATE_L2_E1, task);

  /* Compute the TRSM */
  API_CALL(factor_trsm2d)(sopalin_data, me, task);

  trace_end_task(thread_data->tracefile,
     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
     STATE_L2_E1, task);

  print_debug(DBG_SOPALIN_E1,
      "STASK_TASKDST %ld STASK_PROCDST %ld STASK_BCOFPTR %x (btagtab[%ld])\n",
      (long)STASK_TASKDST(task), (long)STASK_PROCDST(task),
      (unsigned int)(intptr_t)STASK_BCOFPTR(task),
        (long)SOLV_INDTAB[TASK_INDNUM(task)]);

  /*
   * Unlock communication and local tasks linked to this
   * computation and add it in the adapted tasks list
   */
#if (defined COMM_REORDER) || (defined PASTIX_DYNSCHED)
  API_CALL(compute_unlock_after_DiagE1)(sopalin_data, task);
#endif

  /*
   * Try to send available block or signal to
   * communication thread that he has some work to do
   */
  if (THREAD_FUNNELED_OFF)
    {
      API_CALL(send_all_block)(sopalin_data, me);
    }
  else
    {
      pthread_cond_broadcast(&(sopalin_data->cond_comm));
    }

  trace_begin_task(thread_data->tracefile,
       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
       STATE_L2_E1, task);

  /* Clean Buffers if last task using it */
  API_CALL(compute_cleane1e2)(sopalin_data, me, task);

  trace_end_task(thread_data->tracefile,
     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
     STATE_L2_E1, task);
}

/****************************************************************************/
/* COMPUTE TASK E2                                                          */
/****************************************************************************/

void API_CALL(compute_e2)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task)
{
  PASTIX_INT c, b2, dima, dimi, dimj, stride, stridec, t;
  PASTIX_INT offsetC, offsetA;
  PASTIX_FLOAT         *gaik, *gajk, *gc;
#ifdef SOPALIN_LU
  PASTIX_FLOAT         *gaik2, *gajk2, *gc2;
#endif
  SolverMatrix  *datacode    = sopalin_data->datacode;
#ifdef TRACE_SOPALIN
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif

  c  = TASK_CBLKNUM(task);
  b2 = TASK_BLOKNUM(task);

  /* even if local gajk=RTASK_COEFTAB */

  trace_begin_task(thread_data->tracefile,
       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
       STATE_L2_E2, task);

  print_debug(DBG_SOPALIN_E2, "RTASK_BCOFPTR :%x\n",
      (unsigned int)(intptr_t)RTASK_BCOFPTR(task));
  print_debug(DBG_SOPALIN_E2, "RTASK_COEFTAB :%x\n",
      (unsigned int)(intptr_t)RTASK_COEFTAB(task));

  ASSERTDBG(TASK_BTAGPTR(task),MOD_SOPALIN);
  ASSERTDBG(RTASK_COEFTAB(task),MOD_SOPALIN);

  offsetA = SOLV_COEFIND(b2);
  gajk = (PASTIX_FLOAT *)RTASK_COEFTAB(task);
  gaik = &(SOLV_COEFTAB(c)[offsetA]);

  /* vertical dimension */
  dimj = RTASK_LROWNUM(task) - RTASK_FROWNUM(task) + 1;

  /* horizontal dimension */
  dima = SYMB_LCOLNUM(c) - SYMB_FCOLNUM(c) + 1;

  /* vertical dimension */
  dimi = SYMB_LROWNUM(b2) - SYMB_FROWNUM(b2) + 1;

  stride = SOLV_STRIDE(c);

#ifdef SOPALIN_LU
  /* Swap L and U */
  gajk2 = gajk;
  gajk  = gajk2 + dima * dimj;
  gaik2 = &(SOLV_UCOEFTAB(c)[offsetA]);
#endif

  t = SOLV_INDTAB[TASK_INDNUM(task)];
  if ( t < 0 )
  {
    PASTIX_INT cbl = TASK_CBLKNUM(-t);
    PASTIX_INT b3  = TASK_BLOKNUM(-t);

    /* gaij is local */
    print_debug(DBG_SOPALIN_E2, "add local\n");

    stridec = SOLV_STRIDE(cbl);

    ooc_wait_for_cblk(sopalin_data, cbl, me);

    offsetC  = SOLV_COEFIND(b3);
    offsetC += (RTASK_FROWNUM(task) - SYMB_FCOLNUM(cbl)) * stridec;
    offsetC += SYMB_FROWNUM(b2) - SYMB_FROWNUM(b3);

    gc = &(SOLV_COEFTAB(cbl)[offsetC]);

    ASSERTDBG(( SYMB_FROWNUM(b2)    >= SYMB_FROWNUM(b3) ) &&
        ( SYMB_LROWNUM(b2)    <= SYMB_LROWNUM(b3) ) &&
        ( RTASK_FROWNUM(task) >= SYMB_FCOLNUM(b3) ) &&
        ( RTASK_FROWNUM(task) <= SYMB_LCOLNUM(b3) ), MOD_SOPALIN);

    MUTEX_LOCK(&(sopalin_data->mutex_blok[b3]));
    /* multAikB(gc,gaik,gb,stride,dimj,dima,dimi) */
#ifdef HERMITIAN
    SOPALIN_GEMM("N",
                 "C",
                 dimi, dimj, dima,
                 -fun, gaik, stride,
                 gajk, dimj,
                 fun,  gc,   stridec);
#else
    SOPALIN_GEMM("N",
                 "T",
                 dimi, dimj, dima,
                 -fun, gaik, stride,
                 gajk, dimj,
                 fun,  gc,   stridec);
#endif
#ifdef SOPALIN_LU
    if ( b3 != SYMB_BLOKNUM(cbl) )
    {
      gc2 = &(SOLV_UCOEFTAB(cbl)[offsetC]);
      SOPALIN_GEMM("N", "T", dimi, dimj, dima,
             -fun, gaik2, stride,
               gajk2, dimj,
             fun,  gc2,   stridec);
    }
    /* Replace b1 != b2 since b1 is not available */
    /* no need to test LROWNUM */
    else if ( RTASK_FROWNUM(task) != SYMB_FROWNUM(b2) )
    {
      gc2 = &(SOLV_COEFTAB(cbl)[offsetC]);
      /* Store U and L directly in L for factorization */
      SOPALIN_GEMM("N", "T", dimj, dimi, dima,
             -fun, gajk2, dimj,
               gaik2, stride,
             fun,  gc2,   stridec);
    }
#endif
    MUTEX_UNLOCK(&(sopalin_data->mutex_blok[b3]));

    /* update task E1 or DIAG cnt */
    MUTEX_LOCK(&(sopalin_data->mutex_task[-t]));
    TASK_CTRBCNT(-t)--;
    ASSERTDBG((TASK_CTRBCNT(-t) >= 0), MOD_SOPALIN);

#ifdef PASTIX_DYNSCHED
    if ( ( (!TASK_CTRBCNT(-t)) && (sopalin_data->taskmark[-t] == -1) ) &&
    (!((TASK_TASKID(-t) == E1) && ((TASK_BTAGPTR(-t) == NULL) || (RTASK_COEFTAB(-t) == NULL)))))
  {
    PASTIX_INT i;

#if (DBG_PASTIX_DYNSCHED > 0)
    ASSERTDBG(sopalin_data->taskmark[-t] == -1, MOD_SOPALIN);
#endif
    sopalin_data->taskmark[-t]++;
    MUTEX_UNLOCK(&(sopalin_data->mutex_task[-t]));

    i = TASK_THREADID(-t);
    MUTEX_LOCK(&(sopalin_data->tasktab_mutex[i]));
    queueAdd(&(sopalin_data->taskqueue[i]),
       -t,
       (double)TASK_PRIONUM(-t));
    MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[i]));
    pthread_cond_broadcast(&(sopalin_data->tasktab_cond[i]));
  }
    else
  MUTEX_UNLOCK(&(sopalin_data->mutex_task[-t]));
#else
    MUTEX_UNLOCK(&(sopalin_data->mutex_task[-t]));
    pthread_cond_broadcast(&(sopalin_data->cond_task[-t]));
#endif

    ooc_save_coef(sopalin_data, task, TASK_CBLKNUM(-t), me);
  }
  else
  {
    /* gaij is not local */
    print_debug(DBG_SOPALIN_E2, "add fanin\n");

    stridec = FANIN_LROWNUM(t) - FANIN_FROWNUM(t) + 1;
#ifdef ALLOC_FTGT
    MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
#ifdef OOC_FTGT
    print_debug(DBG_OOC_FTGT, "WAIT %4d %4d\n", (int)t, (int) task);
    ooc_wait_for_ftgt(sopalin_data, t, me);
    ASSERTDBG(((unsigned long)(*(((double*)FANIN_COEFTAB(t))-1))) ==
        sizeof(PASTIX_FLOAT)*((FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1) * (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1)   *
                 ((sopalin_data->sopar->factotype == API_FACT_LU)?2:1))
        , MOD_SOPALIN);
#else
    if (FANIN_COEFTAB(t) == NULL)
  {
    PASTIX_INT j;
    PASTIX_INT ftgtsize = stridec*(FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);

#ifdef SOPALIN_LU
      ftgtsize *= 2;
#endif
    MALLOC_INTERN(FANIN_COEFTAB(t), ftgtsize, PASTIX_FLOAT);

    for (j=0;j<ftgtsize;j++)
    FANIN_COEFTAB(t)[j] = fzero;

    print_debug(DBG_SOPALIN_ALLOC, "alloc fanin coeff %x\n",(unsigned int)(intptr_t)FANIN_COEFTAB(t));

      STATS_ADD( ftgtsize );
  }
#endif /* OOC_FTGT */
#endif /* ALLOC_FTGT */

    offsetC = ( RTASK_FROWNUM(task) - FANIN_FCOLNUM(t) ) * stridec
    +         SYMB_FROWNUM(b2) - FANIN_FROWNUM(t);

    gc = &(FANIN_COEFTAB(t)[offsetC]);

    ASSERTDBG((SYMB_FROWNUM(b2)>=FANIN_FROWNUM(t))&&
        (FANIN_LROWNUM(t)>=SYMB_LROWNUM(b2))&&
        (RTASK_FROWNUM(task)>=FANIN_FCOLNUM(t))&&
        (FANIN_LCOLNUM(t)>=RTASK_LROWNUM(task)),MOD_SOPALIN);

    /* multAikB(gc,gaik,gb,stride,dimj,dima,dimi) */
#ifdef HERMITIAN
    SOPALIN_GEMM("N",
                 "C",
                 dimi, dimj, dima,
                 -fun,  gaik, stride,
                 gajk, dimj,
                 fzero, gc,   stridec);
#else
    SOPALIN_GEMM("N",
                 "T",
                 dimi, dimj, dima,
                 -fun,  gaik, stride,
                 gajk, dimj,
                 fzero, gc,   stridec);
#endif

#ifdef SOPALIN_LU
    /* Block extra-diagonal => contribution on U */
    if ( FANIN_FROWNUM(t) != FANIN_FCOLNUM(t) )
    {
      gc2 = gc + stridec*(FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);
      SOPALIN_GEMM("N", "T", dimi, dimj, dima,
                   -fun, gaik2, stride,
                   gajk2, dimj,
                   fun,  gc2,   stridec);
    }
    /* Fanin target diagonal and contribution not diagonal */
    /* => Contribution on L */
    else if ( RTASK_FROWNUM(task) != SYMB_FROWNUM(b2) )
    {
      /* Store U and L directly in L for factorization */
      offsetC = ( RTASK_FROWNUM(task) - FANIN_FCOLNUM(t) )
        +       ( SYMB_FROWNUM(b2)    - FANIN_FROWNUM(t) ) * stridec;
      gc2 = &(FANIN_COEFTAB(t)[offsetC]);
      SOPALIN_GEMM("N", "T", dimj, dimi, dima,
                   -fun, gajk2, dimj,
                   gaik2, stride,
                   fun,  gc2,   stridec);
    }
    else {
    /* We are on diagonal block => we don't do anything */
    }
#endif

#ifdef ALLOC_FTGT
    ooc_save_ftgt(sopalin_data, task, t, me);
    MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
#endif

    /* update fanin target cnt */
    MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
    FANIN_CTRBCNT(t)--;
    /* Put fanin in ready heap */
#ifdef COMM_REORDER
    if (FANIN_CTRBCNT(t) == 0)
  {
    PASTIX_INT dest = FANIN_PROCDST(t);
    MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
  if (THREAD_FUNNELED_ON)
    {
      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      queueAdd2(sopalin_data->sendqueue, t,
                (double)(dest+1), (PASTIX_INT)FANIN_PRIONUM(t));
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
    }
  else
    {
      MUTEX_LOCK(&(sopalin_data->mutex_queue_fanin[dest]));
      queueAdd2(&(sopalin_data->fanintgtsendqueue[dest]), t,
                ((double)FANIN_PRIONUM(t)), t);
      MUTEX_UNLOCK(&(sopalin_data->mutex_queue_fanin[dest]));
    }
  }
    else
#endif
  MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));

    trace_end_task(thread_data->tracefile,
       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
       STATE_L2_E2, task);

    if (THREAD_FUNNELED_OFF)
      {
        API_CALL(send_all_fanin)(sopalin_data, me, FANIN_PROCDST(t));
      }
    trace_begin_task(thread_data->tracefile,
         SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
         STATE_L2_E2, task);
  }

  /* Clean Buffers if last task using it */
  API_CALL(compute_cleane1e2)(sopalin_data, me, task);

  trace_end_task(thread_data->tracefile,
                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                 STATE_L2_E2, task);
}

#ifdef WITH_STARPU
#include "./starpu_kernels.c"
#include "./starpu_updo_kernels.c"
#endif
