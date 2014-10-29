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
#undef EXACT_TAG
#undef UPDOWN_SM2XMAX
#define UPDOWN_SM2XMAX 4096
*/
#include "csc_intern_compute.h"

#ifdef PASTIX_EZTRACE
#  include "pastix_eztrace.h"
#else /* not PASTIX_EZTRACE */
#  include "trace.h"
#endif /* not PASTIX_EZTRACE */

/* ??? Attention a TASK_PRIONUM en 2D                 */
/* ??? Attention free apres Isend                     */
/* ??? Attention ne fonctionne pas sans l'aggregation */

#include "updo_sendrecv.h"

#define updo_thread API_CALL(updo_thread)
/* Lancement de updo seul */
void updo_thread ( SolverMatrix * datacode, SopalinParam * sopar);

/*********************************/
/*
 * Function: API_CALL(updo_thread)
 *
 * Launch threads for solving step.
 *
 * Parameters:
 *   datacode  - SolverMatrix structure (common data)
 *   sopaparam - sopalin parameters.
 *
 * Returns:
 *   void
 */
void updo_thread ( SolverMatrix *datacode, SopalinParam *sopaparam )
{
  Sopalin_Data_t *sopalin_data = NULL;
  BackupSolve_t b;

  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  solve_backup(datacode,&b);
  /*   if (sopaparam->iparm[IPARM_ONLY_RAFF] == 1) */
  /*     sopalin_init(sopalin_data, datacode, sopaparam); */
  sopalin_init(sopalin_data, datacode, sopaparam, 0);
#ifdef WITH_STARPU
  if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {
      starpu_submit_tasks(sopalin_data);
    }
  else
#endif
    {
      sopalin_launch_thread(sopalin_data,
                            SOLV_PROCNUM,          SOLV_PROCNBR,                datacode->btree,
                            sopalin_data->sopar->iparm[IPARM_VERBOSE],
                            SOLV_THRDNBR,          API_CALL(up_down_smp),       sopalin_data,
                            sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                            OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);
    }
    sopalin_clean(sopalin_data, 2);
    solve_restore(datacode,&b);

    memFree_null(sopalin_data);
}

/*
 * Function: API_CALL(updo_down_smp)
 *
 * Compute solve step in shared memory
 * processus mode.
 *
 * Parameters:
 * arg - address of a *sopthread_data_t* structure
 * containing the thread number and the *sopalin_data*.
 *
 * Returns:
 * EXIT_SUCCESS
 *
 */
#define updo_down_smp API_CALL(updo_down_smp)
void* up_down_smp ( void *arg )
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  PASTIX_INT               me           = argument->me;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  SolverMatrix     *datacode     = sopalin_data->datacode;
  SopalinParam     *sopar        = sopalin_data->sopar;
  Thread_Data_t    *thread_data  = sopalin_data->thread_data[me];
  /* Taille et buffer de communication (inutilisé si thread comm) */
  void *            updo_buffer      = NULL;
  PASTIX_INT               updo_buffer_size = 0;
  PASTIX_FLOAT            *ga,*gb,*gc;
  PASTIX_INT               i,ii,j,size,sizea,stride;
  int               init;
#ifdef DEP_SMX
  PASTIX_INT               count_cblk;
#endif
#ifndef STORAGE
  PASTIX_INT               count;
#else
  int               itersmx;
  PASTIX_FLOAT            *lgb;
  PASTIX_FLOAT            *lgrhs;
#endif
#ifdef PASTIX_DYNSCHED
  PASTIX_INT               itasktab, itasktab2, nblevel;
  PASTIX_INT               l;
  int              *indbubble = NULL;
#else
  Queue             cblreadyqueue;
#endif
  PASTIX_FLOAT *b = NULL;
  MONOTHREAD_BEGIN;

  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX)
    fprintf(stdout,  "%d:%d up_down_smp\n", (int)SOLV_PROCNUM, (int)me);

  ooc_set_step(sopalin_data, OOCSTEP_DOWN);

  trace_begin_task(sopalin_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                   STATE_L0_UPDOINIT, 0);

  if (sopalin_data->sopar->iparm[IPARM_END_TASK] == API_TASK_SOLVE &&
      sopalin_data->sopar->iparm[IPARM_PRODUCE_STATS] == API_YES) {
    MALLOC_INTERN(b, UPDOWN_SM2XSZE, PASTIX_FLOAT);
    SOPALIN_COPY(UPDOWN_SM2XSZE,sopar->b,iun,b,iun);
  }

  MONOTHREAD_END;

#ifdef SOPALIN_LU
#  ifdef PASTIX_DUMP_SOLV
  API_CALL(dump_all)(datacode, sopar->cscmtx, DUMP_SOLV);
#  endif
  MONOTHREAD_BEGIN
  if (sopalin_data->sopar->iparm[IPARM_TRANSPOSE_SOLVE] == API_YES)
    {
      PASTIX_INT itercblk;
      PASTIX_INT iun = 1;
      PASTIX_FLOAT * tmp;

      /* U <- 1/diag(U) * U */
      for (itercblk=0; itercblk<SYMB_CBLKNBR; itercblk++)
        {
          PASTIX_INT col;

          stride = SOLV_STRIDE(itercblk);
          for (col = 0;
               col < SYMB_LCOLNUM(itercblk)-SYMB_FCOLNUM(itercblk)+1;
               col++)
            {
              if (stride-(col+1) > 0)
                {
                  PASTIX_FLOAT d  = ((PASTIX_FLOAT)1.)/SOLV_UCOEFTAB(itercblk)[
                    SOLV_COEFIND(SYMB_BLOKNUM(itercblk)) +
                    col*(stride+1)];
                  PASTIX_FLOAT *v = &(SOLV_UCOEFTAB(itercblk)[
                                 SOLV_COEFIND(SYMB_BLOKNUM(itercblk)) +
                                 col*(stride+1)+1]);
                  SOPALIN_SCAL(stride-(col+1), d, v, iun);
                }
            }
        }

      /* L <-  L * diag(U) */
      for (itercblk=0; itercblk<SYMB_CBLKNBR; itercblk++)
        {
          PASTIX_INT col;
          stride = SOLV_STRIDE(itercblk);
          for (col = 0;
               col < SYMB_LCOLNUM(itercblk)-SYMB_FCOLNUM(itercblk)+1;
               col++)
            {
              if (stride-(col+1) > 0)
                {
                  PASTIX_FLOAT d  = SOLV_COEFTAB(itercblk)[
                    SOLV_COEFIND(SYMB_BLOKNUM(itercblk))+
                    col*(stride+1)];
                  PASTIX_FLOAT *v = &(SOLV_COEFTAB(itercblk)[
                                 SOLV_COEFIND(SYMB_BLOKNUM(itercblk))+
                                 col*(stride+1)+1]);
                  SOPALIN_SCAL(stride-(col+1), d, v, iun);
                }
            }
        }

      for (itercblk=0; itercblk<SYMB_CBLKNBR; itercblk++)
        {
          tmp = SOLV_COEFTAB(itercblk);
          SOLV_COEFTAB(itercblk) = SOLV_UCOEFTAB(itercblk);
          SOLV_UCOEFTAB(itercblk) = tmp;
        }
    }
  MONOTHREAD_END

#  ifdef PASTIX_DUMP_SOLV
  API_CALL(dump_all)(datacode, sopar->cscmtx, DUMP_SOLV);
#  endif /* PASTIX_DUMP_SOLV */

#endif
  /* Initialisation des données par thread */
  init = INIT_COMPUTE;
  if (THREAD_FUNNELED_OFF)
    {
      init = init | INIT_SEND;
      if (THREAD_COMM_OFF)
        {
          init = init | INIT_RECV;
        }
    }
  sopalin_init_smp(sopalin_data, me, 0, init);
  thread_data = sopalin_data->thread_data[me];

  if (THREAD_COMM_OFF)
    {
      /* Initialisation buffer communication */
      updo_buffer_size = UPDOWN_SIZETAB*sizeof(PASTIX_INT)+
        UPDOWN_SM2XMAX*sizeof(PASTIX_FLOAT)*UPDOWN_SM2XNBR;

      MALLOC_INTERN(updo_buffer, updo_buffer_size, char);

      print_debug(DBG_SOPALIN_UPDO, "  buffer size = %d\n", (int)updo_buffer_size);
    }

  /******************************/
  /*    INITIALISATION DOWN     */
  /******************************/

  /* Init Down monothread */
  SYNCHRO_THREAD;
  SOPALIN_CLOCK_INIT;
  MONOTHREAD_BEGIN;
  /* Init SM2X */
#ifdef UPDO_DEADCODE
  for (i=0;i<SYMB_CBLKNBR;i++)
    for (j=0;j<SYMB_LCOLNUM(i)-SYMB_FCOLNUM(i)+1;j++)
    {
      UPDOWN_SM2XTAB[UPDOWN_SM2XIND(i)+j] += 3.0*(PASTIX_FLOAT)SYMB_NODENBR;
    }
#endif /* UPDO_DEADCODE */

  /* Waiting for CTRBNBR contributions on Xi */
  for (i=0;i<SYMB_CBLKNBR;i++)
  {
    UPDOWN_CTRBCNT(i) = UPDOWN_CTRBNBR(i); /* SOLV_CTRBNBR(i); ??? */
    UPDOWN_MSGCNT(i)  = UPDOWN_MSGNBR(i);
  }

  for (i=0;i<SOLV_FTGTNBR;i++)
  {
    FANIN_CTRBCNT(i) = 0;
#ifndef OOC_FTGT
    FANIN_COEFTAB(i) = NULL; /* Should be ok ... ??? */
#endif
  }

  for (i=0;i<SYMB_CBLKNBR;i++)
    for (j=SYMB_BLOKNUM(i)+1;j<SYMB_BLOKNUM(i+1);j++)
      if (SYMB_CBLKNUM(j)<=0)
      {
        FANIN_CTRBCNT(-SYMB_CBLKNUM(j))++;
      }

  for (i=0;i<SOLV_FTGTNBR;i++)
    FANIN_CTRBNBR(i) = FANIN_CTRBCNT(i);

#ifdef STORAGE
  for (i=0;i<UPDOWN_GNODENBR*UPDOWN_SM2XNBR;i++) sopalin_data->grhs[i]    = 0.0;
  for (i=0;i<UPDOWN_GCBLKNBR;i++) sopalin_data->flagtab[i] = 0;
#endif

#ifdef PASTIX_DUMP_SOLV
  API_CALL(dump_all)(datacode, sopar->cscmtx, DUMP_SMB);
#endif

  MONOTHREAD_END;
  SYNCHRO_THREAD;

#ifndef SMP_SOPALIN
  /* allocate file structures */
  queueInit(&cblreadyqueue,SYMB_CBLKNBR);

  for (i=0;i<SYMB_CBLKNBR;i++)
#  ifdef DEP_SMX
    /* Dans le cas de dep_smx, on ne met que les taches pretes */
    if (!UPDOWN_CTRBCNT(i))
#  endif
      queueAdd(&cblreadyqueue,i,(double)(TASK_PRIONUM(i)));

#else /* SMP_SOPALIN */
#  ifdef PASTIX_DYNSCHED
  {
    int bubnum = me;
    int bubnum2 = me;

    while (bubnum != -1)
    {
      /* Allocation de la file de tâches */
      if (sopalin_data->taskqueue[bubnum].size == 0)
        queueInit(&(sopalin_data->taskqueue[bubnum]), datacode->ttsknbr[bubnum]);

      for (i=0;i<datacode->ttsknbr[bubnum];i++)
      {
        PASTIX_INT task = datacode->ttsktab[bubnum][i];
        if (!UPDOWN_CTRBCNT(task))
        {
          queueAdd(&(sopalin_data->taskqueue[bubnum]), task, ((double)TASK_PRIONUM(task)));
        }
        if (task > (SOLV_TASKNBR-1))
          errorPrint("Pb task trop grand\n");
      }
      bubnum = BFATHER(datacode->btree, bubnum2);
      if ((bubnum != -1) &&
          (datacode->btree->nodetab[bubnum].fcandnum != datacode->btree->nodetab[bubnum2].fcandnum))
        bubnum = -1;
      bubnum2 = bubnum;
    }
  }
  MONOTHREAD_BEGIN;
  for (i=0; i<SOLV_BUBLNBR; i++)
    sopalin_data->tasktab_indice[i] = 0;
  MONOTHREAD_END;
#  else /* PASTIX_DYNSCHED */

  queueInit(&cblreadyqueue,SOLV_TTSKNBR);
  for (ii=0;ii<SOLV_TTSKNBR;ii++)
  {
    if ((TASK_TASKID(SOLV_TTSKTAB(ii))!=COMP_1D) &&
        (TASK_TASKID(SOLV_TTSKTAB(ii))!=DIAG)) continue;
    i = TASK_CBLKNUM(SOLV_TTSKTAB(ii));
    queueAdd(&cblreadyqueue,i,(double)(TASK_PRIONUM(i)));
  }
#  endif
#endif /* SMP_SOPALIN */

  /* Synchro Fin Initialisation */
  SYNCHRO_THREAD;

  MONOTHREAD_BEGIN;
  trace_begin_task(sopalin_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me,
                   0, STATE_L0_DOWN, 0);
  if (THREAD_COMM_ON)
    {
      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      sopalin_data->step_comm = COMMSTEP_DOWN;
      print_debug(DBG_THCOMM, "%s:%d DOWN\n", __FILE__, __LINE__);
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
      pthread_cond_broadcast(&(sopalin_data->cond_comm));
    }
  MONOTHREAD_END;

  /******************************/
  /*          DOWN              */
  /******************************/
  SOPALIN_CLOCK_STOP;
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX)
    fprintf(stdout, OUT4_UPDO_TIME_INIT,
            (int)SOLV_PROCNUM, (int)me, SOPALIN_CLOCK_GET);

  /* Initialisation du compteur solve */
  SOPALIN_CLOCK_INIT;
  COMM_CLOCK_INIT;

  /* Down Step */
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_one("%s", OUT2_SOP_DOWN);

  /*
   * Différentes boucles suivant les options
   * la premiere partie ne contient le type de boucle et l'attente
   * des contributions sur la tâche à traiter
   *
   */
#ifdef PASTIX_DYNSCHED
  itasktab  = me;
  itasktab2 = me;
  while(1){
    PASTIX_INT stolen, bloknum;
    ii = API_CALL(sopalin_dynsched_getNexTask)( sopalin_data, datacode, thread_data,
                                                &itasktab, &itasktab2, &bloknum, me );

    stolen = itasktab != itasktab2;
    if ( ii == -1 )
      break;

    if ((TASK_TASKID(ii)!=COMP_1D) &&
        (TASK_TASKID(ii)!=DIAG)) continue;

    i = TASK_CBLKNUM(ii);
    /* Ignore schur */
    if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
        UPDOWN_LOC2GLOB(i) == UPDOWN_GCBLKNBR-1)
      continue;

#elif (defined DEP_SMX)

  count_cblk = SYMB_CBLKNBR;
  while (count_cblk)
  {
    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                     STATE_WAITREM, 0);

    while (!queueSize(&cblreadyqueue))
    {
      MPI_Status status;
      COMM_CLOCK_START;

      /* receive and add contributions */
      CALL_MPI MPI_Recv(updo_buffer,updo_buffer_size,MPI_BYTE,
                        MPI_ANY_SOURCE,0,PASTIX_COMM,&status);
      TEST_MPI("MPI_Recv");
      COMM_CLOCK_STOP;
      API_CALL(updo_down_recv)(sopalin_data, updo_buffer, status, me);
    }

    i = queueGet(&cblreadyqueue);
#  error "Not implemented"

#elif (defined SMP_SOPALIN)

  for (ii=0;ii<SOLV_TTSKNBR;ii++)
  {
    i = queueGet(&cblreadyqueue);

    /* Ignore Schur */
    if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
        UPDOWN_LOC2GLOB(i) == UPDOWN_GCBLKNBR -1)
      continue;
    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                     STATE_WAITREM, i);

#  if (!defined FORCE_NOMPI)
    if (THREAD_COMM_OFF)
      {
        /* On attend les messages */
        while (UPDOWN_MSGCNT(i))
          {
            MPI_Status status;
            int        tag;

#    ifdef EXACT_TAG
            tag = i;
#    elif (defined EXACT_THREAD)
            tag = me;
#    else /* not EXACT_TAG nor EXACT_THREAD */
            tag = TAG_DOWN;
#    endif /* not EXACT_TAG nor EXACT_THREAD */
#    if (defined PASTIX_UPDO_ISEND) && (defined FORCE_CONSO)
            API_CALL(send_testall_down)(sopalin_data, me);
#    endif /* PASTIX_UPDO_ISEND && FORCE_CONSO */
            COMM_CLOCK_START;
            CALL_MPI MPI_Recv(updo_buffer, updo_buffer_size, MPI_BYTE,
                              MPI_ANY_SOURCE, tag, PASTIX_COMM, &status);
            TEST_MPI("MPI_Recv");
            COMM_CLOCK_STOP;
            API_CALL(updo_down_recv)(sopalin_data, updo_buffer, status, me);
          }
      }
#  endif /* !FORCE_NOMPI */
    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                     STATE_WAITLOC, i);

    /* On attend toutes les contributions */
    MUTEX_LOCK(&(sopalin_data->mutex_task[i]));
    while(UPDOWN_CTRBCNT(i))
      COND_WAIT(&(sopalin_data->cond_task[i]),
                &(sopalin_data->mutex_task[i]));
    MUTEX_UNLOCK(&(sopalin_data->mutex_task[i]));

#else /* not SMP_SOPALIN */

  for (ii=0;ii<SYMB_CBLKNBR;ii++)
  {
    i = queueGet(&cblreadyqueue);

    /* ignore schur */
    if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
        UPDOWN_LOC2GLOB(i) == UPDOWN_GCBLKNBR -1) {
      continue;
    }
    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                     STATE_WAITREM, i);

#  ifndef FORCE_NOMPI
    while(UPDOWN_CTRBCNT(i))
    {
      MPI_Status status;
      int        tag;

      tag = TAG_DOWN;
      if (THREAD_COMM_OFF)
        {
#    ifdef EXACT_TAG
          tag = i;
#    elif (defined EXACT_THREAD)
          tag = me;
#    endif /* EXACT_THREAD */
        }

#    if (defined PASTIX_UPDO_ISEND) && (defined FORCE_CONSO)
      API_CALL(send_testall_down)(sopalin_data, me);
#    endif
      COMM_CLOCK_START;
      CALL_MPI MPI_Recv(updo_buffer, updo_buffer_size, MPI_BYTE,
                        MPI_ANY_SOURCE, tag, PASTIX_COMM, &status);
      TEST_MPI("MPI_Recv");
      COMM_CLOCK_STOP;
      API_CALL(updo_down_recv)(sopalin_data, updo_buffer, status, me);
    }
#  endif /* FORCE_NOMPI */

#endif /* Fin des differentes boucles Down */

    /*********************************/
    /* Block de calcul de la tache i */
    /*********************************/
    {

      /* ignore schur */
      if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
          UPDOWN_LOC2GLOB(i) == UPDOWN_GCBLKNBR-1)
        continue;

      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                       STATE_DOWN, TASK_COLOR(i));
      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                       STATE_L2_DOWN, i);

      print_debug(DBG_SOPALIN_DOWN, "Calcul de X%d\n",(int)i);

      ooc_wait_for_cblk(sopalin_data, i, me);

      ga     =&(SOLV_COEFTAB(i)[SOLV_COEFIND(SYMB_BLOKNUM(i))]);
      gb     =&(UPDOWN_SM2XTAB[UPDOWN_SM2XIND(i)]);
      stride = SOLV_STRIDE(i);
      size   = SYMB_LCOLNUM(i) - SYMB_FCOLNUM(i) + 1;

      /* Xi=(Lii)-1Xi; */
      MUTEX_LOCK(&(sopalin_data->mutex_task[i]));

#if (defined CHOL_SOPALIN) && (!defined SOPALIN_LU)
#  ifdef MULT_SMX
      SOPALIN_TRSM("L","L","N","N",size,UPDOWN_SM2XNBR,fun,ga,stride,gb,UPDOWN_SM2XSZE);
#  else /* MULT_SMX */
      SOPALIN_TRSV("L","N","N",size,ga,stride,gb,iun);
#  endif
#else
#  ifdef MULT_SMX
          SOPALIN_TRSM("L","L","N","U",size,UPDOWN_SM2XNBR,fun,ga,stride,gb,UPDOWN_SM2XSZE);
#  else /* MULT_SMX */
          SOPALIN_TRSV("L","N","U",size,ga,stride,gb,iun);
#  endif
#endif

      MUTEX_UNLOCK(&(sopalin_data->mutex_task[i]));

#ifdef DEP_SMX
      count_cblk--;
#endif

      trace_end_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                     STATE_L2_DOWN, i);
      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                       STATE_L2_ADD, i);

#ifdef COMPACT_SMX

      ooc_wait_for_cblk(sopalin_data, i, me);

      ga    = &(SOLV_COEFTAB(i)[SOLV_COEFIND(SYMB_BLOKNUM(i)+1)]);
      sizea = stride - size;
      gc    = thread_data->maxbloktab2;
      thread_data->stridebloktab=0;
#  ifdef MULT_SMX
      ASSERTDBG(sizea*UPDOWN_SM2XNBR<=SOLV_COEFMAX,MOD_SOPALIN);
      ASSERTDBG(sizea*UPDOWN_SM2XSZE<=SOLV_COEFMAX,MOD_SOPALIN);

      SOPALIN_GEMM("N","N",sizea,UPDOWN_SM2XNBR,size,fun,ga,stride,gb,UPDOWN_SM2XSZE,fzero,gc,UPDOWN_SM2XSZE);
#  else /* MULT_SMX */
      SOPALIN_GEMV("N",sizea,size,fun,ga,stride,gb,iun,fzero,gc,iun);
#  endif /* MULT_SMX */
#endif /* COMPACT_SMX */

      /******************************/
      /*  AJOUT/ENVOI CONTRIBUTIONS */
      /******************************/

      for (j=SYMB_BLOKNUM(i)+1;j<SYMB_BLOKNUM(i+1);j++)
      {
        /********************************/
        /* if the contribution is LOCAL */
        /********************************/

        if (SYMB_CBLKNUM(j)>0) {
          PASTIX_INT cblknum = SYMB_CBLKNUM(j);
          /* ignore schur */
          if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
              UPDOWN_LOC2GLOB(cblknum) == UPDOWN_GCBLKNBR-1) {
            UPDOWN_CTRBCNT(cblknum)--;
            continue;
          }


          sizea = SYMB_LROWNUM(j) - SYMB_FROWNUM(j) + 1;

#ifdef COMPACT_SMX
          ga = thread_data->maxbloktab2 + thread_data->stridebloktab;
          thread_data->stridebloktab += sizea;
#else /* COMPACT_SMX */
          ga    =&(SOLV_COEFTAB(i)[SOLV_COEFIND(j)]);
#endif /* COMPACT_SMX */

          gc    =&(UPDOWN_SM2XTAB[UPDOWN_SM2XIND(cblknum)+
                                  SYMB_FROWNUM(j)-
                                  SYMB_FCOLNUM(cblknum)]);

          /* Xk=Xk-LkiXi */
          print_debug(DBG_SOPALIN_DOWN,"Mise A Jour locale de X%d\n",(int)cblknum);

          MUTEX_LOCK(&(sopalin_data->mutex_task[cblknum]));

#ifdef COMPACT_SMX
#  ifdef MULT_SMX
          SOPALIN_GESM("N","N",sizea,UPDOWN_SM2XNBR,-fun,ga,UPDOWN_SM2XSZE,gc,UPDOWN_SM2XSZE);
#  else /* MULT_SMX */
          SOPALIN_AXPY(sizea,-fun,ga,iun,gc,iun);
#  endif /* MULT_SMX */
#else /* COMPACT_SMX */
#  ifdef MULT_SMX
          SOPALIN_GEMM("N","N",sizea,UPDOWN_SM2XNBR,size,-fun,ga,stride,gb,UPDOWN_SM2XSZE,fun,gc,UPDOWN_SM2XSZE);
#  else /* MULT_SMX */
          SOPALIN_GEMV("N",sizea,size,-fun,ga,stride,gb,iun,fun,gc,iun);
#  endif /* MULT_SMX */
#endif /* COMPACT_SMX */

          UPDOWN_CTRBCNT(cblknum)--;
          if (!UPDOWN_CTRBCNT(cblknum))
          {
            MUTEX_UNLOCK(&(sopalin_data->mutex_task[cblknum]));
#ifdef DEP_SMX
            queueAdd(&cblreadyqueue, cblknum, (double)(TASK_PRIONUM(cblknum)));
#elif (defined PASTIX_DYNSCHED)
            {
              int threadid = TASK_THREADID(cblknum);
              MUTEX_LOCK(&(sopalin_data->tasktab_mutex[threadid]));
              queueAdd(&(sopalin_data->taskqueue[threadid]),
                       cblknum, (double)(TASK_PRIONUM(cblknum)));
              MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[threadid]));
              pthread_cond_broadcast(&(sopalin_data->tasktab_cond[threadid]));
            }
#else
            pthread_cond_broadcast(&(sopalin_data->cond_task[cblknum]));
#endif
          }
          else
            MUTEX_UNLOCK(&(sopalin_data->mutex_task[cblknum]));

        }
        /************************************/
        /* if the contribution is NOT LOCAL */
        /************************************/
#ifndef FORCE_NOMPI
        else
        {
          PASTIX_INT infotab[UPDOWN_SIZETAB];
          int tabsize[2];

          /* ignore schur */
          if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
              FANIN_CBLKDST(SOLV_FTGTIND(j)) ==  UPDOWN_GCBLKNBR-1)
            continue;

          /* En-tête du message */
          infotab[0] = FANIN_FCOLNUM(SOLV_FTGTIND(j));
          infotab[1] = FANIN_LCOLNUM(SOLV_FTGTIND(j));
          infotab[2] = FANIN_CBLKDST(SOLV_FTGTIND(j));
          infotab[3] = FANIN_CTRBNBR(SOLV_FTGTIND(j));

          tabsize[0] = UPDOWN_SIZETAB;
          tabsize[1] = (infotab[1]-infotab[0]+1)*UPDOWN_SM2XNBR;

          sizea = SYMB_LROWNUM(j)-SYMB_FROWNUM(j)+1;

#  ifdef COMPACT_SMX
          ga    = thread_data->maxbloktab2 + thread_data->stridebloktab;
          thread_data->stridebloktab += sizea;
#  else /* COMPACT_SMX */
          ga    = &(SOLV_COEFTAB(i)[SOLV_COEFIND(j)]);
#  endif /* COMPACT_SMX */


          MUTEX_LOCK(&(sopalin_data->mutex_fanin[SOLV_FTGTIND(j)]));
          ooc_wait_for_ftgt(sopalin_data, SOLV_FTGTIND(j), me);

          if (FANIN_COEFTAB(SOLV_FTGTIND(j)) == NULL)
          {
            PASTIX_INT k;
            MALLOC_INTERN(FANIN_COEFTAB(SOLV_FTGTIND(j)),
                          tabsize[1], PASTIX_FLOAT);
            for (k=0; k<tabsize[1]; k++)
              FANIN_COEFTAB(SOLV_FTGTIND(j))[k] = fzero;
          }
#  ifdef MEMORY_USAGE
          ASSERTDBG(((unsigned long)(*(((double*)FANIN_COEFTAB(SOLV_FTGTIND(j)))-1))) ==
                    sizeof(PASTIX_FLOAT) * UPDOWN_SM2XNBR *
                    (FANIN_LCOLNUM(SOLV_FTGTIND(j)) - FANIN_FCOLNUM(SOLV_FTGTIND(j)) + 1)
                    , MOD_SOPALIN);
#  endif
          gc = FANIN_COEFTAB(SOLV_FTGTIND(j));

          /* Y=LkiXi */
          print_debug(DBG_SOPALIN_DOWN, "Calcul Contribution pour %d\n",(int)j);

#  ifdef COMPACT_SMX
#    ifdef MULT_SMX
          SOPALIN_GEAM("N","N",sizea,UPDOWN_SM2XNBR,-fun,ga,UPDOWN_SM2XSZE,&(gc[SYMB_FROWNUM(j)-FANIN_FCOLNUM(SOLV_FTGTIND(j))]),infotab[1]-infotab[0]+1);
#    else /* MULT_SMX */
          SOPALIN_AXPY(sizea,-fun,ga,iun,&(gc[SYMB_FROWNUM(j)-FANIN_FCOLNUM(SOLV_FTGTIND(j))]),iun);
#    endif
#  else /* COMPACT_SMX */
#    ifdef MULT_SMX
          SOPALIN_GEMM("N","N",sizea,UPDOWN_SM2XNBR,size,fun,ga,stride,gb,UPDOWN_SM2XSZE,fun,
                       &(gc[SYMB_FROWNUM(j)-FANIN_FCOLNUM(SOLV_FTGTIND(j))]),infotab[1]-infotab[0]+1);
#    else /* MULT_SMX */
          SOPALIN_GEMV("N",sizea,size,fun,ga,stride,gb,iun,fun,
                       &(gc[SYMB_FROWNUM(j)-FANIN_FCOLNUM(SOLV_FTGTIND(j))]),iun);
#    endif /* MULT_SMX */
#  endif /* COMPACT_SMX */

          FANIN_CTRBCNT(SOLV_FTGTIND(j))--;
          if (!(FANIN_CTRBCNT(SOLV_FTGTIND(j))))
          {
            MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[SOLV_FTGTIND(j)]));
            if (THREAD_FUNNELED_OFF)
              {
                API_CALL(updo_down_send)(sopalin_data, me, i, j);
              }
            else
              {
                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                queueAdd2(sopalin_data->sendqueue, i, 0., j);
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
              }
          }
          else
            MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[SOLV_FTGTIND(j)]));

        }
#endif /* FORCE_NOMPI */
      }/* Fin Boucle Ajout contributions */

      trace_end_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                     STATE_L2_ADD, i);
      trace_end_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                     STATE_DOWN, TASK_COLOR(i));

      /*******************************/
      /*            DIAG OOC         */
      /*******************************/

      /* TODO: Pas bon au niveau des comms : A revoir */
#if (!defined CHOL_SOPALIN) && (defined OOC)
      if ((TASK_TASKID(SOLV_TTSKTAB(ii))==COMP_1D) ||
          (TASK_TASKID(SOLV_TTSKTAB(ii))==DIAG))
      {

        trace_begin_task(thread_data->tracefile,
                         SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                         STATE_L2_DIAG, i);

        print_debug(DBG_SOPALIN_DOWN, "Calcul des X%d=D-1X%d\n",(int)i,(int)i);

        ga     =&(SOLV_COEFTAB(i)[SOLV_COEFIND(SYMB_BLOKNUM(i))]);
        gb     =&(UPDOWN_SM2XTAB[UPDOWN_SM2XIND(i)]);
        stride = SOLV_STRIDE(i);
        size   = SYMB_LCOLNUM(i)-SYMB_FCOLNUM(i)+1;
        PASTIX_INT jdiag;
        PASTIX_INT kdiag;
        for (jdiag=0;jdiag<UPDOWN_SM2XNBR;jdiag++)
        {
          PASTIX_INT dec = jdiag*UPDOWN_SM2XSZE;
          for (kdiag=0;kdiag<size;kdiag++)
            gb[dec+kdiag] /= ga[kdiag+kdiag*stride];
        }

        trace_end_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                       STATE_L2_DIAG, i);
      }
#endif /* CHOL_SOPALIN */

      ooc_save_coef(sopalin_data, i, i, me);


    } /* Fin bloc de calcul de la tache i */

#if defined(PASTIX_DYNSCHED)
    MUTEX_LOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
    sopalin_data->tasktab_indice[itasktab2]++;
    MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
#endif
  } /* Fin Boucle principale */
#ifdef _UNUSED_
  }}}
#endif

  /* Fin des communications avant diagonales */
#if (defined PASTIX_UPDO_ISEND && !(defined FORCE_NOMPI))
  if (THREAD_FUNNELED_OFF)
    {
      /* Attente de la fin des communications en envoi */
      if (SOLV_PROCNBR > 1)
        {
          API_CALL(send_waitall_down)(sopalin_data, me);
        }
    }
#endif /* PASTIX_UPDO_ISEND */


  SOPALIN_CLOCK_STOP;
  print_debug(DBG_SOPALIN_UPDO, "%d : down time %lf\n", (int)me, SOPALIN_CLOCK_GET);

  trace_begin_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                   STATE_IDLE, 0);

  /* Synchro Avant diagonale */
  if (THREAD_COMM_ON)
    {
      MONOTHREAD_BEGIN;
      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      while(sopalin_data->step_comm != COMMSTEP_UPDOEND)
        COND_WAIT(&(sopalin_data->cond_comm), &(sopalin_data->mutex_comm));
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
      MONOTHREAD_END;
    }

#ifndef PASTIX_DYNSCHED
  queueExit(&cblreadyqueue);
#endif

  SYNCHRO_THREAD;
#ifdef PASTIX_DUMP_SOLV
  MONOTHREAD_BEGIN;
  API_CALL(dump_all)(datacode, sopar->cscmtx, DUMP_SMB);
  MONOTHREAD_END;
  SYNCHRO_THREAD;
#endif

  /******************************/
  /*        DIAGONALE           */
  /******************************/

  MONOTHREAD_BEGIN;
  trace_begin_task(sopalin_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                   STATE_L0_DIAG, 0);
  MONOTHREAD_END;

  /* On ne fait pas la diagonale dans la facto de cholesky et en OOC */
#if (!defined CHOL_SOPALIN) && (!defined OOC)
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_one("%s", OUT2_SOP_DIAG);


#  define DYNSCHED_START_LOOP                                           \
  while (bubnum != -1)                                                  \
    {                                                                   \
    PASTIX_INT fcandnum = datacode->btree->nodetab[bubnum].fcandnum;    \
    PASTIX_INT lcandnum = datacode->btree->nodetab[bubnum].lcandnum;    \
    for (ii=(me-fcandnum);                                              \
         ii < datacode->ttsknbr[bubnum];                                \
         ii+=(lcandnum-fcandnum+1))
#  define DYNSCHED_END_LOOP }

#  define STATIC_START_LOOP                                \
  for (ii=0; ii < datacode->ttsknbr[bubnum]; ii++)
#  define STATIC_END_LOOP

#  ifdef PASTIX_DYNSCHED
#    define SMP_START_LOOP DYNSCHED_START_LOOP
#    define SMP_END_LOOP DYNSCHED_END_LOOP
#    undef DYNSCHED_END_LOOP
#  else
#    define SMP_START_LOOP STATIC_START_LOOP
#    define SMP_END_LOOP STATIC_END_LOOP
#  endif

#  ifdef SMP_SOPALIN
#    define START_LOOP                          \
  int bubnum  = me;                             \
  SMP_START_LOOP {                              \
    i = datacode->ttsktab[bubnum][ii];          \
    if ((TASK_TASKID(i)!=COMP_1D) &&            \
        (TASK_TASKID(i)!=DIAG)) continue;       \
    i = TASK_CBLKNUM(i);
#    define END_LOOP SMP_END_LOOP }
#  else /* SMP_SOPALIN */
#    define START_LOOP for (i=0;i<SYMB_CBLKNBR;i++)
#    define END_LOOP
#  endif

  /* Xi=(Dii)-1Xi; */
  {
    START_LOOP {
      /* ignore schur */
      if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
          UPDOWN_LOC2GLOB(i) == UPDOWN_GCBLKNBR-1)
        continue;
      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                       STATE_L2_DIAG, i);

      print_debug(DBG_SOPALIN_DOWN, "Calcul des X%d=D-1X%d\n",(int)i,(int)i);

      ooc_wait_for_cblk(sopalin_data, i, me);

      ga     =&(SOLV_COEFTAB(i)[SOLV_COEFIND(SYMB_BLOKNUM(i))]);
      gb     =&(UPDOWN_SM2XTAB[UPDOWN_SM2XIND(i)]);
      stride = SOLV_STRIDE(i);
      size   = SYMB_LCOLNUM(i)-SYMB_FCOLNUM(i)+1;

      for (j=0;j<UPDOWN_SM2XNBR;j++) {
        PASTIX_INT k, dec = j*UPDOWN_SM2XSZE;
        for (k=0;k<size;k++)
          gb[dec+k] /= ga[k+k*stride];
      }

      ooc_save_coef(sopalin_data,i, i, me);

      trace_end_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                     STATE_L2_DIAG, i);
    }
#  ifdef PASTIX_DYNSCHED
    bubnum = BFATHER(datacode->btree, bubnum);
#  endif
    END_LOOP;
  }

#undef START_LOOP
#undef END_LOOP
#undef DYNSCHED_START_LOOP
#undef DYNSCHED_END_LOOP
#undef STATIC_START_LOOP
#undef STATIC_END_LOOP
#undef SMP_START_LOOP
#undef SMP_END_LOOP


  SOPALIN_CLOCK_STOP;
  print_debug(DBG_SOPALIN_UPDO, "%d : diag time %lf\n", (int)me, SOPALIN_CLOCK_GET);

#endif /* CHOL_SOPALIN */

  /******************************/
  /*    INITIALISATION UP       */
  /******************************/

  /* Init Up monothread */
  SYNCHRO_THREAD; /* INUTILE sauf peut-etre pour ooc_set_step */

  MONOTHREAD_BEGIN;
  ooc_set_step(sopalin_data, OOCSTEP_UP);
  trace_begin_task(sopalin_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                   STATE_L0_UPDOINIT, 0);

  for (i=0;i<SYMB_CBLKNBR;i++)
    ASSERTDBG(!UPDOWN_CTRBCNT(i),MOD_SOPALIN);

  /* Waiting for BCOLNBR contributions on Xi */
  for (i=0;i<SYMB_CBLKNBR;i++)
  {
    UPDOWN_CTRBCNT(i) = SYMB_BLOKNUM(i+1)-SYMB_BLOKNUM(i)-1;
  }

#ifdef PASTIX_DUMP_SOLV
  API_CALL(dump_all)(datacode, sopar->cscmtx, DUMP_SMB);
#endif
  MONOTHREAD_END;

  SYNCHRO_THREAD;

#ifndef SMP_SOPALIN
  /* allocate file structures */
  queueInit(&cblreadyqueue,SYMB_CBLKNBR);

  for (i=0;i<SYMB_CBLKNBR;i++)
#  ifdef DEP_SMX
    if (!UPDOWN_CTRBCNT(i))
#  endif
      queueAdd(&cblreadyqueue,i,-(double)(TASK_PRIONUM(i)));

#else /* SMP_SOPALIN */

#  ifdef PASTIX_DYNSCHED
  {
    int bubnum;
    int bubnum2;

    nblevel = 1; bubnum = BFATHER(datacode->btree, me);
    while (bubnum != -1)
    {
      nblevel++;
      bubnum = BFATHER(datacode->btree, bubnum);
    }
    MALLOC_INTERN(indbubble, nblevel, int);
    bubnum = me;
    for (i=0; i<nblevel; i++)
    {
      indbubble[nblevel-i-1] = bubnum;
      bubnum = BFATHER(datacode->btree, bubnum);
    }

    /* Initialisation des files de taches pretes */
    {
      bubnum  = me;
      bubnum2 = me;

      while (bubnum != -1)
      {
        for (i=0;i<datacode->ttsknbr[bubnum];i++)
        {
          PASTIX_INT task = datacode->ttsktab[bubnum][i];
          queueAdd(&(sopalin_data->taskqueue[bubnum]), task, -((double)TASK_PRIONUM(task)));
          if (task > (SOLV_TASKNBR-1))
            errorPrint("Pb task trop grand\n");
        }

        bubnum = BFATHER(datacode->btree, bubnum2);
        if ((bubnum != -1) &&
            (datacode->btree->nodetab[bubnum].fcandnum != datacode->btree->nodetab[bubnum2].fcandnum))
          bubnum = -1;
        bubnum2 = bubnum;
      }
    }
  }
  MONOTHREAD_BEGIN;
  for (i=0; i<SOLV_BUBLNBR; i++)
    sopalin_data->tasktab_indice[i] = 0;
  MONOTHREAD_END;
#  else /* not PASTIX_DYNSCHED */
  queueInit(&cblreadyqueue,SOLV_TTSKNBR);
  for (ii=0;ii<SOLV_TTSKNBR;ii++)
  {
    if ((TASK_TASKID(SOLV_TTSKTAB(ii))!=COMP_1D) &&
        (TASK_TASKID(SOLV_TTSKTAB(ii))!=DIAG)) continue;
    i = TASK_CBLKNUM(SOLV_TTSKTAB(ii));
    queueAdd(&cblreadyqueue,i,-(double)(TASK_PRIONUM(i)));
  }
#  endif
#endif /* SMP_SOPALIN */

  SYNCHRO_THREAD;

  MONOTHREAD_BEGIN;
  trace_begin_task(sopalin_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                   STATE_L0_UP, 0);
  if (THREAD_COMM_ON)
    {
      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      sopalin_data->step_comm = COMMSTEP_UP;
      print_debug(DBG_THCOMM, "%s:%d UP\n", __FILE__, __LINE__);
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
      pthread_cond_broadcast(&(sopalin_data->cond_comm));
    }
  MONOTHREAD_END;

  /******************************/
  /*            UP              */
  /******************************/
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_one("%s", OUT2_SOP_UP);

  /*
   * Différente boucles suivant les options
   */
#ifdef PASTIX_DYNSCHED
  itasktab = 0;
  while(1){
    deb2:

    l = indbubble[itasktab];
    MUTEX_LOCK(&(sopalin_data->tasktab_mutex[l]));

    /* On descend dans l'arbre si on a rien a faire a ce niveau */
    while (queueSize(&(sopalin_data->taskqueue[l])) == 0) {
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[l]));
      itasktab++;
      if (itasktab == nblevel) break;
      l = indbubble[itasktab];
      MUTEX_LOCK(&(sopalin_data->tasktab_mutex[l]));
    }

    /* Il n'y a plus rien a faire, on quitte */
    if (itasktab == nblevel) break;

    i = queueGet(&(sopalin_data->taskqueue[l]));
    if (i == -1) {
      errorPrint("Probleme de file vide");
      goto deb2;
    }
    else
      sopalin_data->tasktab_indice[l]++;
    MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[l]));

    /* ignore schur */
    if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
        UPDOWN_LOC2GLOB(i) == UPDOWN_GCBLKNBR-1) {
      for (j=UPDOWN_BROWPROCNBR(i)-1;j>=0;j--) {
        /* if the contribution is not local */
        if (UPDOWN_BROWPROCTAB(i)[j] != SOLV_PROCNUM) {
#  ifndef FORCE_NOMPI
          if (THREAD_FUNNELED_OFF) {
            API_CALL(updo_up_send)(sopalin_data, me, i, j);
          } else {
            MUTEX_LOCK(&(sopalin_data->mutex_comm));
            queueAdd2(sopalin_data->sendqueue, i, 1., j);
            MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
          }
#  endif /* FORCE_NOMPI */
        }
#  ifndef STORAGE /* On calcul au fur et a mesure les contributions */
        else {
          for (count=UPDOWN_LISTPTR(UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(i)));
               count<UPDOWN_LISTPTR(UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(i))+1);
               count++) {
            PASTIX_INT kk = UPDOWN_LISTCBLK(count);
            UPDOWN_CTRBCNT(kk)--;
          }
        }
#  endif /* STORAGE */
      }
      continue;
    }
    ooc_wait_for_cblk(sopalin_data, i, me);

#elif (defined DEP_SMX)

  count_cblk = SYMB_CBLKNBR;
  while (count_cblk) {
    i = 0;
    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                     STATE_WAITREM, i);
    while (!queueSize(&cblreadyqueue)) {
#  error "Not implemented"
    }}
#elif (defined SMP_SOPALIN)
  for (ii=0;ii<SOLV_TTSKNBR;ii++) {
    i = queueGet(&cblreadyqueue);
    print_debug(DBG_SOPALIN_UP, "%d : Task %4d\n", (int)me, (int)i);

    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                     STATE_WAITREM, i);

    /* ignore schur */
    if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
        UPDOWN_LOC2GLOB(i) == UPDOWN_GCBLKNBR-1) {
      for (j=UPDOWN_BROWPROCNBR(i)-1;j>=0;j--) {
        /* if the contribution is not local */
        if (UPDOWN_BROWPROCTAB(i)[j] != SOLV_PROCNUM) {
#  ifndef FORCE_NOMPI
          if (THREAD_FUNNELED_OFF) {
            API_CALL(updo_up_send)(sopalin_data, me, i, j);
          } else {
            MUTEX_LOCK(&(sopalin_data->mutex_comm));
            queueAdd2(sopalin_data->sendqueue, i, 1., j);
            MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
          }
#  endif /* FORCE_NOMPI */
        }
#  ifndef STORAGE /* On calcul au fur et a mesure les contributions */
        else {
          for (count=UPDOWN_LISTPTR(UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(i)));
               count<UPDOWN_LISTPTR(UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(i))+1);
               count++) {
            PASTIX_INT kk = UPDOWN_LISTCBLK(count);
            UPDOWN_CTRBCNT(kk)--;
          }
        }
#  endif /* STORAGE */
      }
      continue;
    }
    ooc_wait_for_cblk(sopalin_data, i, me);

#  ifndef STORAGE
    while (UPDOWN_CTRBCNT(i))
      while (API_CALL(probe_updown)(PASTIX_COMM, UPDOWN_LOC2GLOB(i)))
#  endif /* STORAGE */

#else /* not SMP_SOPALIN */

  for (ii=0;ii<SYMB_CBLKNBR;ii++) {
    i = queueGet(&cblreadyqueue);

    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                     STATE_WAITREM, i);
    /* ignore schur */
    if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
        SYMB_LCOLNUM(i) ==
        sopalin_data->sopar->gN*sopalin_data->sopar->iparm[IPARM_DOF_NBR]-1) {
      for (j=UPDOWN_BROWPROCNBR(i)-1;j>=0;j--) {
        /* if the contribution is not local */
        if (UPDOWN_BROWPROCTAB(i)[j] != SOLV_PROCNUM) {
#  ifndef FORCE_NOMPI
          if (THREAD_FUNNELED_OFF) {
            API_CALL(updo_up_send)(sopalin_data, me, i, j);
          } else {
            MUTEX_LOCK(&(sopalin_data->mutex_comm));
            queueAdd2(sopalin_data->sendqueue, i, 1., j);
            MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
          }
#  endif /* FORCE_NOMPI */
        }
#  ifndef STORAGE /* On calcul au fur et a mesure les contributions */
        else {
          for (count=UPDOWN_LISTPTR(UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(i)));
               count<UPDOWN_LISTPTR(UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(i))+1);
               count++) {
            PASTIX_INT kk = UPDOWN_LISTCBLK(count);
            UPDOWN_CTRBCNT(kk)--;
          }
        }
#  endif /* STORAGE */ 
      }
     continue;
    }
    ooc_wait_for_cblk(sopalin_data, i, me);

    while (UPDOWN_CTRBCNT(i))

#endif /* Fin des diffrentes boucles Up */

    {
      /* Attente contributions */
#ifdef STORAGE
      API_CALL(updo_up_WaitCtrb_storage)  (sopalin_data, updo_buffer_size, updo_buffer, me, i);
#else
      API_CALL(updo_up_WaitCtrb_nostorage)(sopalin_data, updo_buffer_size, updo_buffer, me, i);
#endif
    }

    /* CALCUL de la tache i */
    {
#ifdef DEP_SMX
      i = queueGet(&cblreadyqueue);
#endif

      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                       STATE_UP, TASK_COLOR(i));
      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                       STATE_L2_UP, i);

      print_debug(DBG_SOPALIN_UP,"Calcul de X%d\n", (int)i);

#ifdef SOPALIN_LU
      ga=&(SOLV_UCOEFTAB(i)[SOLV_COEFIND(SYMB_BLOKNUM(i))]);
#else
      ga=&(SOLV_COEFTAB(i)[SOLV_COEFIND(SYMB_BLOKNUM(i))]);
#endif

      gb     = &(UPDOWN_SM2XTAB[UPDOWN_SM2XIND(i)]);
      stride = SOLV_STRIDE(i);
      size   = SYMB_LCOLNUM(i)-SYMB_FCOLNUM(i)+1;

      /* Xi=(Liit)-1Xi; */
      MUTEX_LOCK(&(sopalin_data->mutex_task[i]));
#ifdef CHOL_SOPALIN
#  ifdef MULT_SMX
      SOPALIN_TRSM("L","L","T","N",size,UPDOWN_SM2XNBR,fun,ga,stride,gb,UPDOWN_SM2XSZE);
#  else /* MULT_SMX */
      SOPALIN_TRSV("L","T","N",size,ga,stride,gb,iun);
#  endif
#else /* CHOL_SOPALIN */
#  ifdef HERMITIAN
#    ifdef MULT_SMX
      SOPALIN_TRSM("L","L","C","U",size,UPDOWN_SM2XNBR,fun,ga,stride,gb,UPDOWN_SM2XSZE);
#    else /* MULT_SMX */
      SOPALIN_TRSV("L","C","U",size,ga,stride,gb,iun);
#    endif
#  else /* not HERMITIAN */
#    ifdef MULT_SMX
      SOPALIN_TRSM("L","L","T","U",size,UPDOWN_SM2XNBR,fun,ga,stride,gb,UPDOWN_SM2XSZE);
#    else /* MULT_SMX */
      SOPALIN_TRSV("L","T","U",size,ga,stride,gb,iun);
#    endif
#  endif /* not HERMITIAN */
#endif
      MUTEX_UNLOCK(&(sopalin_data->mutex_task[i]));

#ifdef DEP_SMX
      count_cblk--;
#endif

      trace_end_task  (thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                       STATE_L2_UP, i);
      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                       STATE_L2_ADD, i);

#ifdef STORAGE
      /* deposer dans grhs; */
      print_debug(DBG_SOPALIN_UP, "me=%ld depose i=%ld gi=%ld\n",
                  (long)me,(long)i,(long)UPDOWN_LOC2GLOB(i));

      ASSERTDBG(SYMB_FCOLNUM(i)<UPDOWN_GNODENBR,MOD_SOPALIN);

      /* WARNING : Est-ce que c'est normal que cette copie ne soit pas protégée ? */
      lgb   = gb;
      lgrhs = &(sopalin_data->grhs[SYMB_FCOLNUM(i)]);
      for (itersmx=0; itersmx<UPDOWN_SM2XNBR;
           itersmx++, lgb = lgb + UPDOWN_SM2XSZE, lgrhs = lgrhs + UPDOWN_GNODENBR)
        SOPALIN_COPY(size, lgb, iun, lgrhs, iun);
      MUTEX_LOCK(&(sopalin_data->mutex_flagtab[UPDOWN_LOC2GLOB(i)]));
      sopalin_data->flagtab[UPDOWN_LOC2GLOB(i)] = 1;
      MUTEX_UNLOCK(&(sopalin_data->mutex_flagtab[UPDOWN_LOC2GLOB(i)]));
      pthread_cond_broadcast(&(sopalin_data->cond_flagtab[UPDOWN_LOC2GLOB(i)]));
#endif

      /******************************/
      /*  AJOUT/ENVOI CONTRIBUTIONS */
      /******************************/

      for (j=UPDOWN_BROWPROCNBR(i)-1;j>=0;j--)
      {

        print_debug(DBG_SOPALIN_UP, "Contribution %d de X%d vers (proc %d)\n",
                    (int)j, (int)i, (int)UPDOWN_BROWPROCTAB(i)[j]);

        /* if the contribution is not local */
        if (UPDOWN_BROWPROCTAB(i)[j] != SOLV_PROCNUM)
        {
#ifndef FORCE_NOMPI
          if (THREAD_FUNNELED_OFF)
            {
             API_CALL(updo_up_send)(sopalin_data, me, i, j);
            }
          else
            {
              MUTEX_LOCK(&(sopalin_data->mutex_comm));
              queueAdd2(sopalin_data->sendqueue, i, 1., j);
              MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
            }
#endif /* FORCE_NOMPI */
        }

#ifndef STORAGE /* On calcul au fur et a mesure les contributions */
        /* if the contribution is local */
        else
        {
          PASTIX_INT infotab[UPDOWN_SIZETAB];

          infotab[0] = SYMB_FCOLNUM(i);
          infotab[1] = SYMB_LCOLNUM(i);
          infotab[2] = UPDOWN_LOC2GLOB(i);
          infotab[3] = 0;

          for (count=UPDOWN_LISTPTR(UPDOWN_GCBLK2LIST(infotab[2]));
               count<UPDOWN_LISTPTR(UPDOWN_GCBLK2LIST(infotab[2])+1);
               count++)
          {
            PASTIX_INT kk = UPDOWN_LISTCBLK(count);
            PASTIX_INT k  = UPDOWN_LISTBLOK(count);

            ASSERTDBG((SYMB_FROWNUM(k)>=infotab[0]) &&
                      (SYMB_LROWNUM(k)<=infotab[1]),MOD_SOPALIN);

            print_debug(DBG_SOPALIN_UP, "trouve blok %d -> MAJ\n", (int)k);

            ooc_wait_for_cblk(sopalin_data, kk, me);

#  ifdef SOPALIN_LU
            ga =&(SOLV_UCOEFTAB(kk)[SOLV_COEFIND(k)]);
#  else
            ga =&(SOLV_COEFTAB(kk)[SOLV_COEFIND(k)]);
#  endif
            stride = SOLV_STRIDE(kk);
            size   = SYMB_LCOLNUM(kk) - SYMB_FCOLNUM(kk)+1;
            sizea  = SYMB_LROWNUM(k)  - SYMB_FROWNUM(k) +1;
            gc     =&(UPDOWN_SM2XTAB[UPDOWN_SM2XIND(kk)]);

            MUTEX_LOCK(&(sopalin_data->mutex_task[kk]));
#  ifdef HERMITIAN
#    ifdef MULT_SMX
            SOPALIN_GEMM("C","N",size,UPDOWN_SM2XNBR,sizea,-fun,ga,stride,
                         gb+SYMB_FROWNUM(k)-infotab[0],UPDOWN_SM2XSZE,fun,gc,UPDOWN_SM2XSZE);
#    else /* MULT_SMX */
            SOPALIN_GEMV("C",sizea,size,-fun,ga,stride,
                         gb+SYMB_FROWNUM(k)-infotab[0],iun,fun,gc,iun);
#    endif
#  else /* not HERMITIAN */
#    ifdef MULT_SMX
            SOPALIN_GEMM("T","N",size,UPDOWN_SM2XNBR,sizea,-fun,ga,stride,
                         gb+SYMB_FROWNUM(k)-infotab[0],UPDOWN_SM2XSZE,fun,gc,UPDOWN_SM2XSZE);
#    else /* MULT_SMX */
            SOPALIN_GEMV("T",sizea,size,-fun,ga,stride,
                         gb+SYMB_FROWNUM(k)-infotab[0],iun,fun,gc,iun);
#    endif
#  endif /* not HERMITIAN */
            UPDOWN_CTRBCNT(kk)--;
            if (!UPDOWN_CTRBCNT(kk))
            {
#  ifdef DEP_SMX
              queueAdd(&cblreadyqueue,kk,-(double)(TASK_PRIONUM(kk)));
#  endif
              pthread_cond_broadcast(&(sopalin_data->cond_task[kk]));
            }
            MUTEX_UNLOCK(&(sopalin_data->mutex_task[kk]));
          }
        }
#endif /* STORAGE */
      } /* Fin Boucle Contributions */
    } /* Fin calcul tache i */

    ooc_save_coef(sopalin_data, i, i, me);

    trace_end_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                   STATE_L2_ADD, i);
    trace_end_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                   STATE_UP, TASK_COLOR(i));
  }/* Fin boucle principale */
#ifdef _UNUSED_
  }}}
#endif

#if (defined PASTIX_UPDO_ISEND  && !(defined FORCE_NOMPI))
  if (THREAD_FUNNELED_OFF)
    {
      /* Attente de la fin des communications en envoi */
      if (SOLV_PROCNBR > 1)
        {
          API_CALL(send_waitall_up)(sopalin_data, me);
        }
    }
#endif /* PASTIX_UPDO_ISEND */

  SOPALIN_CLOCK_STOP;
  print_debug(DBG_SOPALIN_UPDO, "%d : updown time %lf\n", (int)me, SOPALIN_CLOCK_GET);
  set_dparm(sopar->dparm, DPARM_SOLV_TIME, SOPALIN_CLOCK_GET);

  trace_begin_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                   STATE_IDLE, 0);

#ifdef OOC
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] < API_TASK_REFINE ||
      sopalin_data->sopar->iparm[IPARM_START_TASK] > API_TASK_NUMFACT)
    ooc_stop_thread(sopalin_data);
#endif

  SYNCHRO_THREAD;

  MONOTHREAD_BEGIN;
  trace_begin_task(sopalin_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                   STATE_L0_UPDOCLEAN, 0);


  for (i=0;i<SYMB_CBLKNBR;i++)
    ASSERTDBG(!UPDOWN_CTRBCNT(i),MOD_SOPALIN);


  if (THREAD_COMM_ON)
    {
      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      while(sopalin_data->step_comm != COMMSTEP_UPDOEND)
        COND_WAIT(&(sopalin_data->cond_comm), &(sopalin_data->mutex_comm));
      if ((sopar->stopthrd == API_YES)
          || (sopar->iparm[IPARM_END_TASK] == API_TASK_SOLVE))
        {
          sopalin_data->step_comm = COMMSTEP_END;
          print_debug(DBG_THCOMM, "%s:%d END\n", __FILE__, __LINE__);
        }
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
      pthread_cond_broadcast(&(sopalin_data->cond_comm));
    }
  MONOTHREAD_END;
  SYNCHRO_THREAD;

#ifndef PASTIX_DYNSCHED
  queueExit(&cblreadyqueue);
#else
  memFree_null(indbubble);
#endif

  if (THREAD_COMM_OFF)
    memFree_null(updo_buffer);


  sopalin_clean_smp(sopalin_data, me);

  MONOTHREAD_BEGIN;
  trace_begin_task(sopalin_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                   STATE_L0_IDLE, 0);
#ifdef SOPALIN_LU
  if (sopalin_data->sopar->iparm[IPARM_TRANSPOSE_SOLVE] == API_YES)
    {
      PASTIX_INT itercblk;
      PASTIX_INT iun = 1;
      PASTIX_FLOAT * tmp;
      for (itercblk=0; itercblk<SYMB_CBLKNBR; itercblk++)
        {
          tmp = SOLV_COEFTAB(itercblk);
          SOLV_COEFTAB(itercblk) = SOLV_UCOEFTAB(itercblk);
          SOLV_UCOEFTAB(itercblk) = tmp;
        }

      /* U <- diag(u) * U */
      for (itercblk=0; itercblk<SYMB_CBLKNBR; itercblk++)
        {
          PASTIX_INT col;
          stride = SOLV_STRIDE(itercblk);
          for (col = 0;
               col < SYMB_LCOLNUM(itercblk)-SYMB_FCOLNUM(itercblk)+1;
               col++)
            {
              if (stride-(col+1) > 0)
                {
                  PASTIX_FLOAT d = SOLV_UCOEFTAB(itercblk)[SOLV_COEFIND(SYMB_BLOKNUM(itercblk)) +
                                                    col*(stride+1)];

                  PASTIX_FLOAT *v = &(SOLV_UCOEFTAB(itercblk)[SOLV_COEFIND(SYMB_BLOKNUM(itercblk)) +
                                                       col*(stride+1)+1]);
                  SOPALIN_SCAL(stride-(col+1), d, v, iun);
                }
            }
        }

      /* L <- L * 1/diag(U) */
      for (itercblk=0; itercblk<SYMB_CBLKNBR; itercblk++)
        {
          PASTIX_INT col;
          stride = SOLV_STRIDE(itercblk);
          for (col = 0;
               col < SYMB_LCOLNUM(itercblk)-SYMB_FCOLNUM(itercblk)+1;
               col++)
            {
              if (stride-(col+1) > 0)
                {
                  PASTIX_FLOAT d = ((PASTIX_FLOAT)1.)/SOLV_COEFTAB(itercblk)[SOLV_COEFIND(SYMB_BLOKNUM(itercblk)) +
                                                               col*(stride+1)];
                  PASTIX_FLOAT *v = &(SOLV_COEFTAB(itercblk)[SOLV_COEFIND(SYMB_BLOKNUM(itercblk)) +
                                                      col*(stride+1)+1]);
                  SOPALIN_SCAL(stride-(col+1), d, v, iun);
                }
            }
        }
    }
#endif
#ifdef PASTIX_DUMP_SOLV
  API_CALL(dump_all)(datacode, sopar->cscmtx, DUMP_SMB);
#endif
  MONOTHREAD_END;

  if (sopalin_data->sopar->iparm[IPARM_END_TASK] == API_TASK_SOLVE &&
      sopalin_data->sopar->iparm[IPARM_PRODUCE_STATS] == API_YES) {
    double prec;
    PASTIX_FLOAT *r, *s;
    MONOTHREAD_BEGIN;
    MALLOC_INTERN(r, UPDOWN_SM2XSZE, PASTIX_FLOAT);
    MALLOC_INTERN(s, UPDOWN_SM2XSZE, PASTIX_FLOAT);
    sopalin_data->ptr_raff[0] = (void *)r;
    sopalin_data->ptr_raff[1] = (void *)s;
    sopalin_data->ptr_raff[2] = (void *)b;
    MONOTHREAD_END;
    SYNCHRO_THREAD;

    r = (PASTIX_FLOAT *)sopalin_data->ptr_raff[0];
    s = (PASTIX_FLOAT *)sopalin_data->ptr_raff[1];
    b = (PASTIX_FLOAT *)sopalin_data->ptr_raff[2];

    /* compute r = b - Ax */
    CscbMAx(sopalin_data, me, r, b, sopalin_data->sopar->cscmtx,
            &(datacode->updovct), datacode, PASTIX_COMM,
            sopar->iparm[IPARM_TRANSPOSE_SOLVE]);
    /* |A||x| + |b| */
    CscAxPb( sopalin_data, me, s, b, sopalin_data->sopar->cscmtx,
             &(datacode->updovct), datacode, PASTIX_COMM,
             sopar->iparm[IPARM_TRANSPOSE_SOLVE]);

    CscBerr(sopalin_data, me, r, s, UPDOWN_SM2XSZE,
            1, &prec , PASTIX_COMM);
    sopalin_data->sopar->dparm[DPARM_SCALED_RESIDUAL] = prec;

    prec = CscNormErr(sopalin_data,
                      me,
                      r,
                      b,
                      UPDOWN_SM2XSZE,
                      1,
                      PASTIX_COMM);
    sopalin_data->sopar->dparm[DPARM_RELATIVE_ERROR] = prec;
    if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO) {
      print_onempi(OUT_PREC1, sopalin_data->sopar->dparm[DPARM_RELATIVE_ERROR]);
      print_onempi(OUT_PREC2, sopalin_data->sopar->dparm[DPARM_SCALED_RESIDUAL]);
    }
    memFree_null(r);
    memFree_null(s);
    memFree_null(b);
  }

  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX)
    fprintf(stdout, OUT4_UPDO_COMM_TIME,
            (int)SOLV_PROCNUM, (int)me, COMM_CLOCK_GET);

  return 0;
}

#include "./updo_sendrecv.c"
