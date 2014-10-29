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
/*************************************/
/*          RECEIVE FUNCTIONS        */
/*************************************/
/* Handle received data */
#define recv_handle_fanin API_CALL(recv_handle_fanin)
#define recv_handle_block API_CALL(recv_handle_block)
void  recv_handle_fanin (Sopalin_Data_t *, PASTIX_INT, void *buffer,
                         MPI_Status status, int elected);
void  recv_handle_block (Sopalin_Data_t *, PASTIX_INT, void *buffer,
                         MPI_Status status, int elected);

/* Wait for one reception */
#define recv_waitone_fanin API_CALL(recv_waitone_fanin)
#define recv_waitone_block API_CALL(recv_waitone_block)
#define recv_waitone_fob   API_CALL(recv_waitone_fob)
void  recv_waitone_fanin(Sopalin_Data_t *, PASTIX_INT, PASTIX_INT tag);
void  recv_waitone_block(Sopalin_Data_t *, PASTIX_INT, PASTIX_INT tag);
void  recv_waitone_fob  (Sopalin_Data_t *, PASTIX_INT);

/* Test reception */
#define recv_testone_fob API_CALL(recv_testone_fob)
#define recv_testall_fab API_CALL(recv_testall_fab)
void  recv_testone_fob  (Sopalin_Data_t *, PASTIX_INT);
void  recv_testall_fab  (Sopalin_Data_t *, PASTIX_INT);

/*************************************/
/*          SEND FUNCTIONS           */
/*************************************/
/* Send one communication */
#define send_one_fanin API_CALL(send_one_fanin)
#define send_one_block API_CALL(send_one_block)
int   send_one_fanin    (Sopalin_Data_t *, PASTIX_INT, PASTIX_INT t);
int   send_one_block    (Sopalin_Data_t *, PASTIX_INT, PASTIX_INT t);

/* Send all available communications */
#define send_all_fanin API_CALL(send_all_fanin)
#define send_all_block API_CALL(send_all_block)
void  send_all_fanin    (Sopalin_Data_t *, PASTIX_INT, PASTIX_INT dest);
void  send_all_block    (Sopalin_Data_t *, PASTIX_INT);

/* Free data structure associate with send request */
#define send_free_fanin API_CALL(send_free_fanin)
#define send_free_block API_CALL(send_free_block)
void  send_free_fanin   (Sopalin_Data_t *, PASTIX_INT, PASTIX_INT s_index);
void  send_free_block   (Sopalin_Data_t *, PASTIX_INT, PASTIX_INT s_index);

/* Test and wait send requests */
#define send_testall       API_CALL(send_testall)
#define send_testall_fanin API_CALL(send_testall_fanin)
#define send_testall_block API_CALL(send_testall_block)
#define send_testall_fab   API_CALL(send_testall_fab)
#define send_waitone       API_CALL(send_waitone)
#define send_waitone_fanin API_CALL(send_waitone_fanin)
#define send_waitone_block API_CALL(send_waitone_block)
#define send_waitall_fab   API_CALL(send_waitall_fab)
void  send_testall      (Sopalin_Data_t *,
                         PASTIX_INT, void (*funcfree)(Sopalin_Data_t*,
                                               PASTIX_INT, PASTIX_INT));
void  send_testall_fanin(Sopalin_Data_t *, PASTIX_INT);
void  send_testall_block(Sopalin_Data_t *, PASTIX_INT);
void  send_testall_fab  (Sopalin_Data_t *, PASTIX_INT);
int   send_waitone      (Sopalin_Data_t *, PASTIX_INT,
                         void (*funcfree)(Sopalin_Data_t*, PASTIX_INT, PASTIX_INT));
int   send_waitone_fanin(Sopalin_Data_t *, PASTIX_INT);
int   send_waitone_block(Sopalin_Data_t *, PASTIX_INT);
void  send_waitall_fab  (Sopalin_Data_t *, PASTIX_INT);

/* Test all receive and send requests */
#define rcsd_testall_fab API_CALL(rcsd_testall_fab)
void  rcsd_testall_fab  (Sopalin_Data_t *, PASTIX_INT);

/* Fonction pour thread de comm */
#define sendrecv_smp API_CALL(sendrecv_smp)
void* sendrecv_smp(void *arg);

#ifdef FORCE_NOMPI
#define recv_waitone_fanin API_CALL(recv_waitone_fanin)
#define recv_waitone_block API_CALL(recv_waitone_block)
#define recv_waitone_fob   API_CALL(recv_waitone_fob)
#define recv_testone_fob   API_CALL(recv_testone_fob)
#define recv_testall_fab   API_CALL(recv_testall_fab)
#define send_all_fanin     API_CALL(send_all_fanin)
#define send_all_block     API_CALL(send_all_block)
#define send_free_fanin    API_CALL(send_free_fanin)
#define send_free_block    API_CALL(send_free_block)
#define rcsd_testall_fab   API_CALL(rcsd_testall_fab)
#define sendrecv_smp       API_CALL(sendrecv_smp)
void  recv_waitone_fanin(Sopalin_Data_t *sopalin_data,
                         PASTIX_INT me, PASTIX_INT tag)
{
  (void)sopalin_data; (void)me; (void)tag;
}
void  recv_waitone_block(Sopalin_Data_t *sopalin_data,
                         PASTIX_INT me, PASTIX_INT tag)
{
  (void)sopalin_data; (void)me; (void)tag;
}
void  recv_waitone_fob  (Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  (void)sopalin_data; (void)me;
}
void  recv_testone_fob  (Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  (void)sopalin_data; (void)me;
}
void  recv_testall_fab  (Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  (void)sopalin_data; (void)me;
}
void  send_all_fanin    (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT dest)
{
  (void)sopalin_data; (void)me; (void)dest;
}
void  send_all_block    (Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  (void)sopalin_data; (void)me;
}
void  send_free_fanin   (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT s_index)
{
  (void)sopalin_data; (void)me; (void)s_index;
}

void  send_free_block   (Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT s_index)
{
  (void)sopalin_data; (void)me; (void)s_index;
}

void  rcsd_testall_fab  (Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  (void)sopalin_data; (void)me;
}
void* sendrecv_smp(void *arg){ (void)arg; return NULL;}
#else

/****************************************************************************/
/* RECEIVE ROUTINES                                                         */
/****************************************************************************/

/*
 * Function: recv_handle_fanin
 *
 * Add fanin contribution received in recv_buffer.
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    recv_buffer  - Received data
 *    status       - MPI communication status
 *    elected      - Index of communication used
 *
 */
void
recv_handle_fanin(Sopalin_Data_t *sopalin_data,
                  PASTIX_INT             me,
                  void           *recv_buffer,
                  MPI_Status      status,
                  int             elected)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
#  if (defined TRACE_SOPALIN) || (defined TEST_IRECV)
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#  endif
  void  *buffer      = recv_buffer;
  PASTIX_FLOAT *ga,  *gb;
#  ifdef SOPALIN_LU
  PASTIX_FLOAT *ga2, *gb2;
#  endif
  PASTIX_INT    packnbr, pack, ind;
  PASTIX_INT    ctrbnbr, taskdst, cblkdst, blokdst, stride, str, dimi, dimj;
  PASTIX_INT    fcolnum, lcolnum, frownum, lrownum, ofrownum, olrownum;
#ifdef TEST_IRECV
  PASTIX_INT    size;
#endif
  int    flag = 0;

  print_debug(DBG_SOPALIN_RECV,"%ld: receive fanin target\n",(long)me);

#ifdef TEST_IRECV
  size = PACKMAX*(sizeof(PASTIX_INT)*MAXINFO)+PACKAREA*sizeof(PASTIX_FLOAT);
#endif
  taskdst = ((PASTIX_INT*)buffer)[FTGT_TASKDST];
  packnbr = ((PASTIX_INT*)buffer)[FTGT_PRIONUM];

  trace_begin_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                   STATE_L2_RECVF, taskdst);
  trace_recv(thread_data->tracefile,
             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, status.MPI_SOURCE,
             COMM_FANIN, taskdst, 0, ((PASTIX_INT*)buffer)[FTGT_IDTRACE]);

  /* Loop on all ftgt */
  for (pack=0;pack<packnbr;pack++)
    {
      /* add fanintarget */
      ctrbnbr = ((PASTIX_INT*)buffer)[FTGT_CTRBNBR];
      taskdst = ((PASTIX_INT*)buffer)[FTGT_TASKDST];
      blokdst = ((PASTIX_INT*)buffer)[FTGT_BLOKDST];
      fcolnum = ((PASTIX_INT*)buffer)[FTGT_FCOLNUM];
      lcolnum = ((PASTIX_INT*)buffer)[FTGT_LCOLNUM];
      frownum = ((PASTIX_INT*)buffer)[FTGT_FROWNUM];
      lrownum = ((PASTIX_INT*)buffer)[FTGT_LROWNUM];
      dimj    = lcolnum - fcolnum + 1;
      dimi    = lrownum - frownum + 1;

#  ifdef DRUNK_SOPALIN
      if (taskdst == -DRUNK)
        {
          PASTIX_INT c;

          /* find cblkdst without taskdst */
          taskdst = SOLV_TASKNBR-1;
          for (c=TASK_CBLKNUM(taskdst); c<SYMB_CBLKNBR; c++)
            if ((blokdst>=SYMB_BLOKNUM(c)) && (blokdst<SYMB_BLOKNUM(c+1)))
              break;

          ASSERTDBG(c!=SYMB_CBLKNBR,MOD_SOPALIN);

          cblkdst = c;

          print_debug(DBG_SOPALIN_DRUNK, "add on DRUNK ctrncnt=%ld\n",
                      (long)TASK_CTRBCNT(SOLV_TASKNBR-1));
          print_debug(DBG_SOPALIN_DRUNK, "cblkdst blokdst DRUNK = %ld %ld \n",
                      (long)c, (long)blokdst);
          print_debug(DBG_SOPALIN_DRUNK, "stride dimi dimj"
                      " DRUNK = %ld %ld %ld \n",
                      (long)stride, (long)dimi, (long)dimj);
        }
      else
#  endif /* DRUNK_SOPALIN */
        {
          cblkdst = TASK_CBLKNUM(taskdst);
        }

      stride  = SOLV_STRIDE(cblkdst);

      /* Load unpredicted cblk */
      ooc_hack_load(sopalin_data, cblkdst, me);
      ooc_wait_for_cblk(sopalin_data, cblkdst, me);

      ind = SOLV_COEFIND(blokdst)
        + (fcolnum - SYMB_FCOLNUM(cblkdst)) * stride
        +  frownum - SYMB_FROWNUM(blokdst);

      ga  = &(SOLV_COEFTAB(cblkdst)[ind]);
#  ifdef SOPALIN_LU
      ga2 = &(SOLV_UCOEFTAB(cblkdst)[ind]);
#  endif

      gb  =(PASTIX_FLOAT *)(((char *) buffer) + MAXINFO*sizeof(PASTIX_INT));

      print_debug(DBG_SOPALIN_RECV,
                  "%ld: Recv fanintarget\n"
                  "%ld: ctrbnbr ctrbcnt procdst taskdst blokdst prionum"
                  " fcolnum lcolnum frownum lrownum (cblkdst)\n"
                  "%ld: %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld (%ld)\n",
                  (long)me, (long)me, (long)me, (long)ctrbnbr,
                  (long)((PASTIX_INT*)buffer)[FTGT_CTRBCNT],
                  (long)((PASTIX_INT*)buffer)[FTGT_PROCDST],
                  (long)taskdst, (long)blokdst,
                  (long)((PASTIX_INT*)buffer)[FTGT_PRIONUM],
                  (long)fcolnum, (long)lcolnum,
                  (long)frownum, (long)lrownum, (long)cblkdst);

      /* pas de test possible sur prionum (=packnbr) */
#  ifdef EXACT_THREAD
      if (THREAD_COMM_OFF)
        ASSERTDBG(me==((PASTIX_INT*)buffer)[FTGT_PROCDST]%SOLV_THRDNBR,MOD_SOPALIN);
#  endif
      ASSERTDBG((SYMB_FCOLNUM(cblkdst)<=fcolnum) &&
                (SYMB_LCOLNUM(cblkdst)>=lcolnum),MOD_SOPALIN);

      str = dimi;

#  ifdef NAPA_SOPALIN
      ofrownum = frownum;
      olrownum = lrownum;
      blokdst--;
      flag = 1;
      do {

        PASTIX_INT trace = 0;

        /* il peut y avoir plusieurs cibles partielles */
        if (!flag)
          print_debug(DBG_SOPALIN_NAPA, "ILU: plusieurs cibles distantes\n");

        frownum = ofrownum;
        lrownum = olrownum;
        blokdst++;

        if ((!flag)
            || (SYMB_FROWNUM(blokdst) > frownum)
            || (SYMB_LROWNUM(blokdst) < lrownum))
          {
            trace = 1;
            if (flag)
              print_debug(DBG_SOPALIN_NAPA,
                          "\nILU: debug fanin"
                          " SF=%ld F=%ld SL=%ld L=%ld gb=%ld ga=%ld"
                          " (stride=%ld SFC=%ld FC=%ld)\n",
                          (long)SYMB_FROWNUM(blokdst), (long)frownum,
                          (long)SYMB_LROWNUM(blokdst), (long)lrownum, (long)0,
                          (long)(SOLV_COEFIND(blokdst)
                                 + (fcolnum-SYMB_FCOLNUM(cblkdst))*stride
                                 + frownum-SYMB_FROWNUM(blokdst)),
                          (long)stride, (long)SYMB_FCOLNUM(cblkdst),
                          (long)fcolnum);
          }

        if (SYMB_FROWNUM(blokdst)>frownum)
          {
            frownum = SYMB_FROWNUM(blokdst);
            print_debug(DBG_SOPALIN_NAPA, "ILU: tronque frownum\n");
          }
        if (SYMB_LROWNUM(blokdst)<lrownum)
          {
            lrownum=SYMB_LROWNUM(blokdst);
            print_debug(DBG_SOPALIN_NAPA, "ILU: tronque lrownum\n");
          }
        dimi = lrownum - frownum + 1;

        ind = SOLV_COEFIND(blokdst)
          + (fcolnum - SYMB_FCOLNUM(cblkdst)) * stride
          +  frownum - SYMB_FROWNUM(blokdst);

        ga  = &(SOLV_COEFTAB(cblkdst)[ind]);
#    ifdef SOPALIN_LU
        ga2 = &(SOLV_UCOEFTAB(cblkdst)[ind]);
#    endif
        gb  = (PASTIX_FLOAT *) ( ((char *) buffer)
                          + MAXINFO*sizeof(PASTIX_INT)
                          +(frownum-ofrownum)*sizeof(PASTIX_FLOAT));

        if (trace)
          {
            print_debug(DBG_SOPALIN_NAPA,
                        "ILU: debug fanin"
                        " SF=%ld F=%ld SL=%ld L=%ld gb=%ld ga=%ld"
                        " (stride=%ld SFC=%ld FC=%ld)\n",
                        (long)SYMB_FROWNUM(blokdst), (long)frownum,
                        (long)SYMB_LROWNUM(blokdst), (long)lrownum,
                        (long)(frownum - ofrownum),
                        (long)(SOLV_COEFIND(blokdst)
                               +(fcolnum-SYMB_FCOLNUM(cblkdst))*stride
                               +frownum-SYMB_FROWNUM(blokdst)),
                        (long)stride, (long)SYMB_FCOLNUM(cblkdst),
                        (long)fcolnum);
          }
#  endif /* NAPA_SOPALIN */

        ASSERTDBG((SYMB_FROWNUM(blokdst) <= frownum) &&
                  (SYMB_LROWNUM(blokdst) >= lrownum),MOD_SOPALIN);

        MUTEX_LOCK(&(sopalin_data->mutex_blok[blokdst]));
        SOPALIN_GEAM("N","N",dimi,dimj,fun,gb,str,ga,stride);
#  ifdef SOPALIN_LU
        gb2 = gb + str*dimj;
        SOPALIN_GEAM("N","N",dimi,dimj,fun,gb2,str,ga2,stride);
#  endif

        ooc_save_coef(sopalin_data, -1, cblkdst, me);
        MUTEX_UNLOCK(&(sopalin_data->mutex_blok[blokdst]));

        /* updatecontrib cnt */
        ASSERTDBG(ctrbnbr > 0, MOD_SOPALIN);

#  ifdef NAPA_SOPALIN
        if (flag)
          {
#  endif

            MUTEX_LOCK(&(sopalin_data->mutex_task[taskdst]));
            TASK_CTRBCNT(taskdst)-=ctrbnbr;
            TASK_FTGTCNT(taskdst)--; /*-=ctrbnbr*/

            /* Unlock taskdst if counter is null */
#  ifdef PASTIX_DYNSCHED
            if ( ( (!TASK_CTRBCNT(taskdst)) &&
                   (sopalin_data->taskmark[taskdst] == -1)) &&
                 (!( (TASK_TASKID(taskdst) == E1) &&
                     ( (TASK_BTAGPTR(taskdst) == NULL) ||
                       (RTASK_COEFTAB(taskdst) == NULL)))))
              {
                PASTIX_INT i;

                ASSERTDBG(TASK_TASKID(taskdst) != E2, MOD_SOPALIN);

#    if (DBG_PASTIX_DYNSCHED > 0)
                ASSERTDBG(sopalin_data->taskmark[taskdst] == -1, MOD_SOPALIN);
#    endif
                sopalin_data->taskmark[taskdst]++;
                MUTEX_UNLOCK(&(sopalin_data->mutex_task[taskdst]));

                i = TASK_THREADID(taskdst);
                MUTEX_LOCK(&(sopalin_data->tasktab_mutex[i]));
                queueAdd(&(sopalin_data->taskqueue[i]),
                         taskdst, (double)TASK_PRIONUM(taskdst));
                MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[i]));
                pthread_cond_broadcast(&(sopalin_data->tasktab_cond[i]));
              }
            else
              MUTEX_UNLOCK(&(sopalin_data->mutex_task[taskdst]));
#  else
            MUTEX_UNLOCK(&(sopalin_data->mutex_task[taskdst]));
            pthread_cond_broadcast(&(sopalin_data->cond_task[taskdst]));
#  endif

            if (THREAD_FUNNELED_ON)
              {
                /* MUTEX_LOCK(&(sopalin_data->mutex_comm));   */
                SOLV_FTGTCNT--;
                /* MUTEX_UNLOCK(&(sopalin_data->mutex_comm)); */
              }

#ifdef NAPA_SOPALIN
          }
        flag = 0;
      } while ((blokdst+1<SYMB_BLOKNUM(cblkdst+1)) &&
               (
                 ((SYMB_FROWNUM(blokdst+1)<=ofrownum) &&
                  (SYMB_LROWNUM(blokdst+1)>=ofrownum)) ||
                 ((SYMB_LROWNUM(blokdst+1)>=olrownum) &&
                  (SYMB_FROWNUM(blokdst+1)<=olrownum)) ||
                 ((SYMB_FROWNUM(blokdst+1)>=ofrownum) &&
                  (SYMB_LROWNUM(blokdst+1)<=olrownum)) ||
                 ((SYMB_FROWNUM(blokdst+1)<=ofrownum) &&
                  (SYMB_LROWNUM(blokdst+1)>=olrownum))
                ));
#endif

      /* Next contribution */
#ifdef SOPALIN_LU
      buffer = ((char *) buffer)+MAXINFO*sizeof(PASTIX_INT)+str*dimj*sizeof(PASTIX_FLOAT)*2;
#else
      buffer = ((char *) buffer)+MAXINFO*sizeof(PASTIX_INT)+str*dimj*sizeof(PASTIX_FLOAT);
#endif

      print_debug(DBG_SOPALIN_RECV, "%ld: fin ajout fanin target\n",(long)me);
    }


  /* If we use Irecv, we launch a new communication on the elected buffer */
#ifdef TEST_IRECV
  {
    CALL_MPI MPI_Irecv(thread_data->recv_fanin_buffer[elected], size, MPI_BYTE,
                       MPI_ANY_SOURCE, me, PASTIX_COMM,
                       &(thread_data->recv_fanin_request[elected]));
    TEST_MPI("MPI_Irecv");
  }
#endif /* TEST_IRECV */

  trace_end_task(thread_data->tracefile,
                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                 STATE_L2_RECVF, taskdst);
}

/*
 * Function: recv_handle_block
 *
 * Add block contribution received in recv_buffer.
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    recv_buffer  - Received data
 *    status       - MPI communication status
 *    elected      - Index of communication used
 *
 */
void
recv_handle_block(Sopalin_Data_t *sopalin_data,
                  PASTIX_INT             me,
                  void           *buffer,
                  MPI_Status      status,
                  int             elected)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
#if (defined TRACE_SOPALIN) || (defined TEST_IRECV)
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif
  BlockTarget   *btag;
  BlockCoeff    *bcof;
  PASTIX_INT            task, bloksize, taskcnt;
#ifdef TEST_IRECV
  PASTIX_INT            size;
#endif
#ifdef FLAG_ASSERT
  PASTIX_INT            taskcntT;
#endif
  PASTIX_INT            bttaskdst;
#ifdef PASTIX_DEBUG
  PASTIX_INT            btprocdst, bttaskcnt, btprionum;
#endif
  PASTIX_INT            bcfrownum, bclrownum, bcfcolnum, bclcolnum;

  print_debug(DBG_SOPALIN_RECV,"%ld: receive block target\n", (long)me);
#ifdef TEST_IRECV
  size = sizeof(PASTIX_INT)*(BTAGINFO+BCOFINFO)+sizeof(PASTIX_FLOAT)*SOLV_BPFTMAX;
#endif
  /* build blocktarget */
  MALLOC_INTERN(btag, 1, BlockTarget);
  MALLOC_INTERN(bcof, 1, BlockCoeff);

  bttaskdst = btag->infotab[BTAG_TASKDST] = ((PASTIX_INT*)buffer)[BTAG_TASKDST];
#ifdef PASTIX_DEBUG
  btprocdst = btag->infotab[BTAG_PROCDST] = ((PASTIX_INT*)buffer)[BTAG_PROCDST];
  btprionum = btag->infotab[BTAG_PRIONUM] = ((PASTIX_INT*)buffer)[BTAG_PRIONUM];
  bttaskcnt = btag->infotab[BTAG_TASKCNT] = ((PASTIX_INT*)buffer)[BTAG_TASKCNT];
#endif
  btag->bcofptr = bcof;

  btag->bcofptr->sendcnt = LOCAL_ALLOC_BTAG; /* need a mark to be free */
  bcfrownum = btag->bcofptr->infotab[BCOF_FROWNUM] = ((PASTIX_INT*)buffer)[BTAGINFO+BCOF_FROWNUM];
  bclrownum = btag->bcofptr->infotab[BCOF_LROWNUM] = ((PASTIX_INT*)buffer)[BTAGINFO+BCOF_LROWNUM];
  bcfcolnum = btag->bcofptr->infotab[BCOF_FCOLNUM] = ((PASTIX_INT*)buffer)[BTAGINFO+BCOF_FCOLNUM];
  bclcolnum = btag->bcofptr->infotab[BCOF_LCOLNUM] = ((PASTIX_INT*)buffer)[BTAGINFO+BCOF_LCOLNUM];

  bloksize = (bclrownum - bcfrownum + 1)*
    (bclcolnum - bcfcolnum + 1);

  trace_begin_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                   STATE_L2_RECVB, bttaskdst);

  trace_recv(thread_data->tracefile,
             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, status.MPI_SOURCE,
             COMM_BLOCK, bttaskdst, bloksize, ((PASTIX_INT*)buffer)[BTAG_IDTRACE]);

  print_debug(DBG_SOPALIN_RECV,
              "%ld: Recv blocktarget\n"
              "%ld: prionum taskdst procdst taskcnt frownum"
              " lrownum fcolnum lcolnum\n"
              "%ld: %ld %ld %ld %ld %ld %ld %ld %ld\n",
              (long)me, (long)me, (long)me,
              (long)btprionum, (long)bttaskdst,
              (long)btprocdst, (long)bttaskcnt,
              (long)bcfrownum, (long)bclrownum,
              (long)bcfcolnum, (long)bclcolnum);

#ifdef EXACT_THREAD
  ASSERTDBG((SEPFB+me) == (SEPFB+btprocdst%SOLV_THRDNBR), MOD_SOPALIN);
#endif
#ifndef PASTIX_DYNSCHED
  ASSERTDBG((btprocdst/SOLV_THRDNBR) == SOLV_PROCNUM,     MOD_SOPALIN);
#endif
  MALLOC_INTERN(btag->bcofptr->coeftab, bloksize, PASTIX_FLOAT);

  STATS_ADD(bloksize);

  memcpy((void *) btag->bcofptr->coeftab,
         (void *) (((char *) buffer)+(BTAGINFO+BCOFINFO)*sizeof(PASTIX_FLOAT)),
         bloksize*sizeof(PASTIX_FLOAT));

  /* link  blocktarget -> tasks */
  task = bttaskdst;

  ASSERTDBG(TASK_BTAGPTR(task)==NULL,MOD_SOPALIN);
  ASSERTDBG(((bclrownum == bclcolnum) && (bcfrownum == bcfcolnum)) ?
            (TASK_TASKID(task)==E1):(TASK_TASKID(task)==E2),MOD_SOPALIN);

  /* Lock on the first task from the cycle */
  MUTEX_LOCK(&(sopalin_data->mutex_task[task]));
  TASK_BTAGPTR(task) = btag;

  taskcnt  = 1;
#ifdef FLAG_ASSERT
  taskcntT = RTASK_TASKCNT(task);
#endif

#ifdef PASTIX_DYNSCHED
  if ( ( (!TASK_CTRBCNT(task)) &&
         (sopalin_data->taskmark[task] == -1)) &&
       (!( (TASK_TASKID(task) == E1) &&
           ( (TASK_BTAGPTR(task) == NULL) ||
             (RTASK_COEFTAB(task) == NULL)))) &&
       (!( (TASK_TASKID(task) == E2) &&
           ( (TASK_BTAGPTR(task) == NULL) ||
             (RTASK_COEFTAB(task) == NULL)))))
    {
      PASTIX_INT i = TASK_THREADID(task);

#  ifdef PASTIX_DYNSCHED
      ASSERTDBG(sopalin_data->taskmark[task] == -1, MOD_SOPALIN);
#  endif
      sopalin_data->taskmark[task]++;

      MUTEX_LOCK(&(sopalin_data->tasktab_mutex[i]));
      queueAdd(&(sopalin_data->taskqueue[i]),
               task, (double)TASK_PRIONUM(task));
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[i]));
      pthread_cond_broadcast(&(sopalin_data->tasktab_cond[i]));
    }
#endif
  MUTEX_UNLOCK(&(sopalin_data->mutex_task[task]));

  print_debug(DBG_SOPALIN_RECV, "%ld: link (%ld[%ld]) ",
              (long)me, (long)task, (long)TASK_TASKID(task));

  /* Add pointer to receive buffer for tasks in the cycle */
  while (TASK_TASKNEXT(task) != bttaskdst)
    {
      /* same id for all links E1 or E2 */
      ASSERTDBG(TASK_TASKID(task) == TASK_TASKID(TASK_TASKNEXT(task)),
                MOD_SOPALIN);

      task = TASK_TASKNEXT(task);
      MUTEX_LOCK(&(sopalin_data->mutex_task[task]));
      TASK_BTAGPTR(task) = btag;

      print_debug(DBG_SOPALIN_RECV, "%ld[%ld] ",
                  (long)task,(long)TASK_TASKID(task));

#ifdef PASTIX_DYNSCHED
      if (((!TASK_CTRBCNT(task)) && (sopalin_data->taskmark[task] == -1)) &&
          (!((TASK_TASKID(task) == E1) && ((TASK_BTAGPTR(task) == NULL) ||
                                           (RTASK_COEFTAB(task) == NULL)))) &&
          (!((TASK_TASKID(task) == E2) && ((TASK_BTAGPTR(task) == NULL) ||
                                           (RTASK_COEFTAB(task) == NULL)))))
        {
          PASTIX_INT i = TASK_THREADID(task);

#  if (DBG_PASTIX_DYNSCHED > 0)
          ASSERTDBG(sopalin_data->taskmark[task] == -1, MOD_SOPALIN);
#  endif
          sopalin_data->taskmark[task]++;

          MUTEX_LOCK(&(sopalin_data->tasktab_mutex[i]));
          queueAdd(&(sopalin_data->taskqueue[i]),
                   task, (double)TASK_PRIONUM(task));
          MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[i]));
          pthread_cond_broadcast(&(sopalin_data->tasktab_cond[i]));
        }
#endif
      MUTEX_UNLOCK(&(sopalin_data->mutex_task[task]));
      taskcnt++;
    }

  ASSERTDBG(taskcnt == taskcntT, MOD_SOPALIN);
  pthread_cond_broadcast(&(sopalin_data->cond_task[TASK_MASTER(bttaskdst)]));

  print_debug(DBG_SOPALIN_RECV, "\n");
  print_debug(DBG_SOPALIN_RECV, "%ld: fin ajout block target\n",(long)me);

#ifdef TEST_IRECV
  /* Restart block reception */
  {
    CALL_MPI MPI_Irecv(thread_data->recv_block_buffer[elected],size,MPI_BYTE,
                       MPI_ANY_SOURCE,SEPFB+me,PASTIX_COMM,
                       &(thread_data->recv_block_request[elected]));
    TEST_MPI("MPI_Irecv");
  }
#endif /* TEST_IRECV */

  trace_end_task(thread_data->tracefile,
                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                 STATE_L2_RECVB, bttaskdst);
}

/*
 * Function: recv_waitone_fanin
 *
 * Wait one fanin communication and call recv_handle_fanin
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    tag          - tag used for communication
 *
 */
void
recv_waitone_fanin(Sopalin_Data_t *sopalin_data,
                   PASTIX_INT             me,
                   PASTIX_INT             tag)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  MPI_Status     status;
  int            elected     = 0;

  print_debug(DBG_SOPALIN_RECV, "%ld: recv_waitone_fanin\n", (long)me);

#ifdef TEST_IRECV

  /* Si on est en irecv, on attend la reception */
  CALL_MPI MPI_Waitany(MAX_R_REQUESTS, thread_data->recv_fanin_request,
                       &elected, &status);
  TEST_MPI("MPI_Waitany");
  thread_data->recv_buffer = thread_data->recv_fanin_buffer[elected];

#else /* TEST_IRECV */

  {
    PASTIX_INT size;

    /* sinon on prepare le mpi_recv */
    tag = TAG_FANIN;
    if (THREAD_COMM_OFF)
      {
#  if (defined EXACT_THREAD)
        tag = me;
#  elif (defined EXACT_TAG)
        tag = tag;
#  endif
      }

    size = PACKMAX*(sizeof(PASTIX_INT)*MAXINFO)+PACKAREA*sizeof(PASTIX_FLOAT);
    OOC_RECEIVING;
    CALL_MPI MPI_Recv(thread_data->recv_buffer, size, MPI_BYTE,
                      MPI_ANY_SOURCE, tag, PASTIX_COMM, &status);
    OOC_RECEIVED;
    TEST_MPI("MPI_Recv");
  }
#endif /* TEST_IRECV */

  recv_handle_fanin(sopalin_data, me, thread_data->recv_buffer,
                    status, elected);
}


/*
 * Function: recv_waitone_block
 *
 * Wait one block communication and call recv_handle_block
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    tag          - tag used for communication
 *
 */
void
recv_waitone_block(Sopalin_Data_t *sopalin_data,
                   PASTIX_INT             me,
                   PASTIX_INT             tag)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  MPI_Status     status;
  int            elected     = 0;

  print_debug(DBG_SOPALIN_RECV, "%ld: recv_waitone_block\n", (long)me);

  /* receive blocktarget */
#ifdef TEST_IRECV

  CALL_MPI MPI_Waitany(MAX_R_REQUESTS, thread_data->recv_block_request,
                       &elected, &status);
  TEST_MPI("MPI_Waitany");
  thread_data->recv_buffer = thread_data->recv_block_buffer[elected];

#else /*TEST_IRECV */

  {
    PASTIX_INT size;

    tag = TAG_BLOCK;
    if (THREAD_COMM_OFF)
      {
#  if (defined EXACT_THREAD)
        tag = SEPFB + me;
#  elif (defined EXACT_TAG)
        tag = SEPFB + tag;
#  endif
      }
    size = sizeof(PASTIX_INT)*(BTAGINFO+BCOFINFO) + sizeof(PASTIX_FLOAT)*SOLV_BPFTMAX;
    CALL_MPI MPI_Recv(thread_data->recv_buffer, size, MPI_BYTE,
                      MPI_ANY_SOURCE, tag, PASTIX_COMM, &status);
    TEST_MPI("MPI_Recv");
  }

#endif /* TEST_IRECV */

  recv_handle_block(sopalin_data, me, thread_data->recv_buffer,
                    status, elected);
}

/*
 * Function: recv_waitone_fob
 *
 * Wait one fanin or one block communication and call
 * associate recv_handle_(fanin|block)
 * Works only without EXACT_TAG or EXACT_THREAD
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *
 */
void
recv_waitone_fob(Sopalin_Data_t *sopalin_data,
                 PASTIX_INT             me)
{
  print_debug(DBG_SOPALIN_RECV, "%ld: recv_waitone_fob\n",(long)me);

#ifdef TEST_IRECV
  {
    SolverMatrix  *datacode    = sopalin_data->datacode;
    Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    MPI_Status    *statuses    = thread_data->srteststatus;
    int           *indices     = thread_data->srtestindices;
    int            i, elected, outcount1, outcount2;

    if (thread_data->srteststatus == NULL)
      {
        MALLOC_INTERN(thread_data->srteststatus,  MAX(MAX_R_REQUESTS,
                                                      MAX_S_REQUESTS),
                      MPI_Status);
        MALLOC_INTERN(thread_data->srtestindices, MAX(MAX_R_REQUESTS,
                                                      MAX_S_REQUESTS), int);
        statuses = thread_data->srteststatus;
        indices  = thread_data->srtestindices;
      }

    outcount1 = 0;
    outcount2 = 0;
    while ((!outcount1) && (!outcount2))
      {
        CALL_MPI MPI_Testsome(MAX_R_REQUESTS, thread_data->recv_fanin_request,
                              &outcount1, indices, statuses);
        TEST_MPI("MPI_Testsome");

        for (i=0;i<outcount1;i++)
          {
            elected = indices[i];
            recv_handle_fanin(sopalin_data, me,
                              thread_data->recv_fanin_buffer[elected],
                              statuses[elected], elected);
          }

        CALL_MPI MPI_Testsome(MAX_R_REQUESTS, thread_data->recv_block_request,
                              &outcount2, indices, statuses);
        TEST_MPI("MPI_Testsome");

        for (i=0;i<outcount2;i++)
          {
            elected = indices[i];
            recv_handle_block(sopalin_data, me,
                              thread_data->recv_block_buffer[elected],
                              statuses[elected],elected);
          }
      }
  }
#else

#  if (defined EXACT_TAG) || (defined EXACT_THREAD)
  /* Can add some deadlock if we wait on any_tag in smp */
#    ifdef SMP_SOPALIN
  errorPrintW("Tag EXACT or THREAD are incompatible with the function"
              " recv_waitone_fob");
#    else
  {
    SolverMatrix  *datacode    = sopalin_data->datacode;
    Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    MPI_Status     status;
    PASTIX_INT            size;
    size = MAX(PACKMAX*(sizeof(PASTIX_INT)*MAXINFO)+
               PACKAREA*sizeof(PASTIX_FLOAT),
               sizeof(PASTIX_INT)*(BTAGINFO+BCOFINFO)+
               sizeof(PASTIX_FLOAT)*SOLV_BPFTMAX);


    /* Test one fanin */
    CALL_MPI MPI_Recv(thread_data->recv_buffer, size, MPI_BYTE,
                      MPI_ANY_SOURCE, MPI_ANY_TAG, PASTIX_COMM, &status);
    TEST_MPI("MPI_Recv");

    if (status.MPI_TAG < SEPFB)
      recv_handle_fanin(sopalin_data, me,
                        thread_data->recv_buffer,
                        status, 0);
    else
      recv_handle_block(sopalin_data, me,
                        thread_data->recv_buffer,
                        status, 0);
  }
#    endif
#  else
  {
    SolverMatrix  *datacode    = sopalin_data->datacode;
    Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    MPI_Status     status;
    int            elected     = 0;
    PASTIX_INT            size;
    size = MAX(PACKMAX*(sizeof(PASTIX_INT)*MAXINFO)+
               PACKAREA*sizeof(PASTIX_FLOAT),
               sizeof(PASTIX_INT)*(BTAGINFO+BCOFINFO)+
               sizeof(PASTIX_FLOAT)*SOLV_BPFTMAX);


    CALL_MPI MPI_Recv(thread_data->recv_buffer, size, MPI_BYTE,
                      MPI_ANY_SOURCE, MPI_ANY_TAG, PASTIX_COMM, &status);
    TEST_MPI("MPI_Recv");
    switch(status.MPI_TAG)
      {
      case TAG_FANIN:
        recv_handle_fanin(sopalin_data, me,
                          thread_data->recv_buffer,
                          status, elected);
        break;
      case TAG_BLOCK:
        recv_handle_block(sopalin_data, me,
                          thread_data->recv_buffer,
                          status, elected);
        break;
      default:
        print_debug(DBG_SOPALIN_COMM, "tag unknown\n");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
      }
  }
#  endif
#endif
}

/*
 * Function: recv_testone_fob
 *
 * Test one fanin or one block communication and call associate
 * recv_handle_(fanin|block)
 * Works only without EXACT_TAG or EXACT_THREAD
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *
 */
void recv_testone_fob(Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  print_debug(DBG_SOPALIN_RECV, "%ld: recv_testone_fob\n",(long)me);

#ifdef TEST_IRECV
  {
    SolverMatrix  *datacode    = sopalin_data->datacode;
    Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    MPI_Status     status;
    int            flag;
    int            elected;

    CALL_MPI MPI_Testany(MAX_R_REQUESTS, thread_data->recv_fanin_request,
                         &elected, &flag, &status);
    TEST_MPI("MPI_Testany");
    if (flag)
      {
        recv_handle_fanin(sopalin_data, me,
                          thread_data->recv_fanin_buffer[elected],
                          status, elected);
      }

    CALL_MPI MPI_Testany(MAX_R_REQUESTS, thread_data->recv_block_request,
                         &elected, &flag, &status);
    TEST_MPI("MPI_Testany");
    if (flag)
      {
        recv_handle_block(sopalin_data, me,
                          thread_data->recv_block_buffer[elected],
                          status, elected);
      }
  }

#else /* Test_Irecv */

#  if (defined EXACT_TAG)
  /* Can add some deadlock if we wait on any_tag in smp */
#    ifdef SMP_SOPALIN
  errorPrintW("Tag EXACT is incompatible with the function recv_testall_fab");
#    else
  {
    MPI_Status    status;
    int           flag;

    flag = 0;
    while ((!flag))
      {
        /* Test one fanin */
        CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, PASTIX_COMM,
                            &flag, &status);
        TEST_MPI("MPI_Iprobe");
        if (flag)
          {
            if (flag < SEPFB)
              recv_waitone_fanin(sopalin_data, me, status.MPI_TAG);
            else
              recv_waitone_block(sopalin_data, me, status.MPI_TAG);
          }
      }
  }
#    endif
#  elif (defined EXACT_THREAD)
  {
#    ifdef SMP_SOPALIN
    SolverMatrix *datacode    = sopalin_data->datacode;
#    endif
    MPI_Status    status;
    int           flag1, flag2;

    flag1 = 0; flag2 = 0;
    while ((!flag1) && (!flag2))
      {
        /* Test one fanin */
        CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, me, PASTIX_COMM,
                            &flag1, &status);
        TEST_MPI("MPI_Iprobe");
        if (flag1)
          recv_waitone_fanin(sopalin_data, me, status.MPI_TAG);

        /* Test one block */
        CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, SEPFB+me, PASTIX_COMM,
                            &flag2, &status);
        TEST_MPI("MPI_Iprobe");

        if (flag2)
          recv_waitone_block(sopalin_data, me, status.MPI_TAG);
      }
  }
#  else
  {
    MPI_Status status;
    int        flag;

    CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, PASTIX_COMM,
                        &flag, &status);
    TEST_MPI("MPI_Iprobe");
    if (flag)
      {
        switch (status.MPI_TAG)
          {
          case TAG_FANIN:
            recv_waitone_fanin(sopalin_data, me, status.MPI_TAG);
            break;
          case TAG_BLOCK:
            recv_waitone_block(sopalin_data, me, status.MPI_TAG);
            break;
          default:
            print_debug(DBG_SOPALIN_COMM, "tag unknown\n");
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
          }
      }
  }
#  endif
#endif /* TEST_IRECV */
}

/*
 * Function: recv_testall_fab
 *
 * Test all active receive communication and call associate
 * recv_handle_(fanin|block)
 * Works only without EXACT_TAG or EXACT_THREAD
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *
 */
void
recv_testall_fab(Sopalin_Data_t *sopalin_data,
                 PASTIX_INT             me)
{

  print_debug(DBG_SOPALIN_RECV, "%ld: recv_testall_fab\n",(long)me);

  /* We don't care about the tag since there is launched requests */
#ifdef TEST_IRECV
  {
    SolverMatrix  *datacode    = sopalin_data->datacode;
    Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    MPI_Status    *statuses;
    int           *indices;
    int            i, elected, outcount;

    if (thread_data->srteststatus == NULL)
      {
        MALLOC_INTERN(thread_data->srteststatus,  MAX(MAX_R_REQUESTS,
                                                      MAX_S_REQUESTS),
                      MPI_Status);
        MALLOC_INTERN(thread_data->srtestindices, MAX(MAX_R_REQUESTS,
                                                      MAX_S_REQUESTS), int);
      }
    statuses = thread_data->srteststatus;
    indices  = thread_data->srtestindices;

    CALL_MPI MPI_Testsome(MAX_R_REQUESTS, thread_data->recv_fanin_request,
                          &outcount, indices, statuses);
    TEST_MPI("MPI_Testsome");

    for (i=0;i<outcount;i++)
      {
        elected = indices[i];
        recv_handle_fanin(sopalin_data, me,
                          thread_data->recv_fanin_buffer[elected],
                          statuses[elected], elected);
      }

    CALL_MPI MPI_Testsome(MAX_R_REQUESTS, thread_data->recv_block_request,
                          &outcount, indices, statuses);
    TEST_MPI("MPI_Testsome");

    for (i=0;i<outcount;i++)
      {
        elected = indices[i];
        recv_handle_block(sopalin_data, me,
                          thread_data->recv_block_buffer[elected],
                          statuses[elected],elected);
      }
  }
#else
#  if (defined EXACT_TAG)
  /* Can add some deadlock if we wait on any_tag */
  errorPrintW("Tag EXACT is incompatible with the function recv_testall_fab");
#  elif (defined EXACT_THREAD)
  {
#    ifdef SMP_SOPALIN
    SolverMatrix *datacode    = sopalin_data->datacode;
#    endif
    MPI_Status    status;
    int           flag1, flag2;

    flag1 = 1;
    flag2 = 1;
    while (flag1 || flag2)
      {
        /* Test one fanin */
        CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, me, PASTIX_COMM,
                            &flag1, &status);
        TEST_MPI("MPI_Iprobe");
        if (flag1)
          recv_waitone_fanin(sopalin_data, me, status.MPI_TAG);

        /* Test one block */
        CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, SEPFB+me, PASTIX_COMM,
                            &flag2, &status);
        TEST_MPI("MPI_Iprobe");

        if (flag2)
          recv_waitone_block(sopalin_data, me, status.MPI_TAG);
      }
  }
#  else
  /* Normally don't use in this section exept if I modifiy some code
   in thread_comm */
  {
    MPI_Status status;
    int        flag;

    flag = 1;
    while (flag)
      {
        CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, PASTIX_COMM,
                            &flag, &status);
        TEST_MPI("MPI_Iprobe");
        if (flag)
          {
            switch (status.MPI_TAG)
              {
              case TAG_FANIN:
                recv_waitone_fanin(sopalin_data, me, status.MPI_TAG);
                break;
              case TAG_BLOCK:
                recv_waitone_block(sopalin_data, me, status.MPI_TAG);
                break;
              default:
                print_debug(DBG_SOPALIN_COMM, "tag unknown\n");
                EXIT(MOD_SOPALIN,INTERNAL_ERR);
              }
          }
      }
  }
#  endif /* tag */
#endif /* TEST_IRECV */
}


/****************************************************************************/
/* SEND ROUTINES MPI                                                        */
/****************************************************************************/

/*
 * Function: send_one_fanin
 *
 * Send all contribution for a same task on a same destination.
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    t            - First fanin number
 *
 * Returns:
 *    Number of fanin sent
 */
int
send_one_fanin ( Sopalin_Data_t *sopalin_data,
                 PASTIX_INT             me,
                 PASTIX_INT             t)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  Queue         *sendqueue;
  PASTIX_INT            id_req  = 0;
  PASTIX_INT            packnbr = 0;
  PASTIX_INT            tdeb, tag;
  PASTIX_INT           *extra; /* contient la liste supplementaire
                         * des fanintgt a liberer */
  PASTIX_INT            sum_pack;
#ifdef NO_MPI_TYPE
  PASTIX_INT            iter;
  PASTIX_INT            copied;
#else  /* NO_MPI_TYPE */
  MPI_Datatype   newtype;
#endif /* NO_MPI_TYPE  */

  print_debug(DBG_SOPALIN_BLEND,
              "%ld: Send fanintarget\n %ld:"
              " ctrnbr ctrbcnt procdst taskdst blokdst"
              " prionum fcolnum lcolnum frownum lrownum ->"
              " tag\n %ld: %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",
              (long)me,(long)me,(long)me,(long)FANIN_CTRBNBR(t),
              (long)FANIN_CTRBCNT(t),(long)FANIN_PROCDST(t),
              (long)FANIN_TASKDST(t),(long)FANIN_BLOKDST(t),
              (long)FANIN_PRIONUM(t),(long)FANIN_FCOLNUM(t),
              (long)FANIN_LCOLNUM(t),(long)FANIN_FROWNUM(t),
              (long)FANIN_LROWNUM(t));

  tdeb  = t;
  extra = NULL;

  /*************************/
  /* add first contibution */
  /*************************/

  /* Number of elements for each type */
  thread_data->gtabsize[packnbr]   = MAXINFO;
  thread_data->gtabsize[packnbr+1] = (FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1) *
    (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);
#ifdef SOPALIN_LU
  thread_data->gtabsize[packnbr+1]*= 2;
#endif

  sum_pack = thread_data->gtabsize[packnbr+1];

  /* Type of each vector */
#ifdef NO_MPI_TYPE
  thread_data->gtabtype[packnbr]   = sizeof(PASTIX_INT);
  thread_data->gtabtype[packnbr+1] = sizeof(PASTIX_FLOAT);
#else /* NO_MPI_TYPE */
  thread_data->gtabtype[packnbr]   = COMM_INT;
  thread_data->gtabtype[packnbr+1] = COMM_FLOAT;
#endif /* NO_MPI_TYPE */

#ifdef OOC_FTGT
  print_debug(DBG_OOC_FTGT, "WAIT %4d %4d\n", (int)t, (int) -1);
  MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
  ooc_wait_for_ftgt(sopalin_data, t, me);
  ASSERTDBG(((unsigned long)(*(((double*)FANIN_COEFTAB(t))-1))) ==
            sizeof(PASTIX_FLOAT)*((FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1) *
                           (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1) *
                           ((sopalin_data->sopar->factotype ==
                             API_FACT_LU)?2:1))
            , MOD_SOPALIN);
  MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
#endif

  print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_INFOTAB : %x\n", (long)me,
              (unsigned int)(intptr_t)FANIN_INFOTAB(t));
  print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_COEFTAB : %x\n", (long)me,
              (unsigned int)(intptr_t)FANIN_COEFTAB(t));
  print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_CTRBCNT : %x\n", (long)me,
              (unsigned int)(intptr_t)FANIN_CTRBCNT(t));

  /* Adress of each vector */
#ifdef NO_MPI_TYPE
  thread_data->gtaboffs[packnbr]   = FANIN_INFOTAB(t);
  thread_data->gtaboffs[packnbr+1] = FANIN_COEFTAB(t);
#else /* NO_MPI_TYPE */
  CALL_MPI MPI_Address(FANIN_INFOTAB(t),&(thread_data->gtaboffs[packnbr]));
  TEST_MPI("MPI_Address");
  CALL_MPI MPI_Address(FANIN_COEFTAB(t),&(thread_data->gtaboffs[packnbr+1]));
  TEST_MPI("MPI_Address");
#endif /* NO_MPI_TYPE */

  /* Add other contribution for the same task */
  if (THREAD_FUNNELED_ON)
    {
      sendqueue = sopalin_data->sendqueue;
    }
  else
    {
      sendqueue = &(sopalin_data->fanintgtsendqueue[FANIN_PROCDST(tdeb)]);
    }

  if (queueSize(sendqueue))
    {
      t = queueRead(sendqueue);

      print_debug(DBG_SOPALIN_COMM,
                  "%ld-%ld C: ftgt %ld / dest %ld / key %ld / task %ld\n",
                  (long)SOLV_PROCNUM, (long)me, (long)t, (long)FANIN_PROCDST(t),
                  (long)FANIN_PRIONUM(t), (long)FANIN_TASKDST(t));

      /* allocation du tableau listant les contributions qui seront envoyees */
      /* Attention : FANIN_CTRBCNT est testé en 2 fois et doit donc
       *             etre protégé */
      MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
      if ((FANIN_PROCDST(tdeb) == FANIN_PROCDST(t)) &&
          (FANIN_TASKDST(tdeb) == FANIN_TASKDST(t)) &&
          (!(FANIN_CTRBCNT(t))) && ((packnbr/2)<(PACKMAX-1)))
        MALLOC_INTERN(extra, PACKMAX, PASTIX_INT);

      while ((FANIN_PROCDST(tdeb) == FANIN_PROCDST(t)) &&
             (FANIN_TASKDST(tdeb) == FANIN_TASKDST(t)) &&
             (!(FANIN_CTRBCNT(t))) && ((packnbr/2)<(PACKMAX-1)))
        {
          PASTIX_INT next;

          MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
          next = (FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1)*
            (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);

#ifdef SOPALIN_LU
          next*=2;
#endif

          /* Si cela fait trop de donnees a envoyer, on sort de la boucle */
          if (sum_pack+next>PACKAREA) break;

          t = queueGet(sendqueue);

          print_debug(DBG_SOPALIN_SEND, "%ld: Extra fanintarget\n"
                      "%ld: ctrnbr ctrbcnt procdst taskdst blokdst prionum"
                      " fcolnum lcolnum frownum lrownum\n"
                      "%ld: %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",
                      (long)me,(long)me,(long)me,(long)FANIN_CTRBNBR(t),
                      (long)FANIN_CTRBCNT(t),(long)FANIN_PROCDST(t),
                      (long)FANIN_TASKDST(t),
                      (long)FANIN_BLOKDST(t),(long)FANIN_PRIONUM(t),
                      (long)FANIN_FCOLNUM(t),
                      (long)FANIN_LCOLNUM(t),(long)FANIN_FROWNUM(t),
                      (long)FANIN_LROWNUM(t));

          extra[packnbr/2] = t;
          packnbr         += 2;

          /* create MPI type */
          thread_data->gtabsize[packnbr]   = MAXINFO;
          thread_data->gtabsize[packnbr+1] =
            (FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1)*
            (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);
#ifdef SOPALIN_LU
          thread_data->gtabsize[packnbr+1]*= 2;
#endif

          sum_pack += thread_data->gtabsize[packnbr+1];

#ifdef NO_MPI_TYPE
          thread_data->gtabtype[packnbr]   = sizeof(PASTIX_INT);
          thread_data->gtabtype[packnbr+1] = sizeof(PASTIX_FLOAT);
#else /* NO_MPI_TYPE */
          thread_data->gtabtype[packnbr]   = COMM_INT;
          thread_data->gtabtype[packnbr+1] = COMM_FLOAT;
#endif /* NO_MPI_TYPE */

#ifdef OOC_FTGT
          print_debug(DBG_OOC_FTGT, "WAIT %4d %4d\n", (int)t, (int) -1);
          MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
          ooc_wait_for_ftgt(sopalin_data, t, me);
          ASSERTDBG(((unsigned long)(*(((double*)FANIN_COEFTAB(t))-1))) ==
                    sizeof(PASTIX_FLOAT)*((FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1) *
                                   (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1)   *
                                   ((sopalin_data->sopar->factotype == API_FACT_LU)?2:1))
                    , MOD_SOPALIN);
          MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
#endif

          print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_INFOTAB : %x\n",
                      (long)me, (unsigned int)(intptr_t)FANIN_INFOTAB(t));
          print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_COEFTAB : %x\n",
                      (long)me, (unsigned int)(intptr_t)FANIN_COEFTAB(t));
          print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_CTRBCNT : %x\n",
                      (long)me, (unsigned int)(intptr_t)FANIN_CTRBCNT(t));

#ifdef NO_MPI_TYPE
          thread_data->gtaboffs[packnbr]   = FANIN_INFOTAB(t);
          thread_data->gtaboffs[packnbr+1] = FANIN_COEFTAB(t);
#else /* NO_MPI_TYPE */
          CALL_MPI MPI_Address(FANIN_INFOTAB(t),
                               &(thread_data->gtaboffs[packnbr]));
          TEST_MPI("MPI_Address");
          CALL_MPI MPI_Address(FANIN_COEFTAB(t),
                               &(thread_data->gtaboffs[packnbr+1]));
          TEST_MPI("MPI_Address");
#endif /* NO_MPI_TYPE */

          if (queueSize(sendqueue))
            {
              t = queueRead(sendqueue);
              MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
            }
          else
            {
              MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
              break;
            }
        }

      MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
      if (extra)
        extra[packnbr/2] = -1;
    }

  if (packnbr/2 > PACKMAX)
    print_debug(DBG_SOPALIN_SEND, "packnbr/2=%ld<?(%ld)\n",
                (long)(packnbr/2),(long)PACKMAX);
  if (sum_pack > PACKAREA)
    print_debug(DBG_SOPALIN_SEND, "sum_pack=%ld<?(%ld)\n",
                (long)(sum_pack),(long)PACKAREA);
  ASSERTDBG(packnbr/2 <= PACKMAX, MOD_SOPALIN);
  ASSERTDBG(sum_pack  <= PACKAREA,MOD_SOPALIN);

  /*   fprintf(thread_data->tracefile,"\n"); */
  print_debug(DBG_SOPALIN_SEND, "%ld: packnbr/2=%ld sum_pack=%ld\n",
              (long)me,(long)(packnbr/2),(long)sum_pack);

  /* Choix du tag suivant la version */
  if (sopalin_data->sopar->type_comm != 3)
    tag = TAG_FANIN;
  else
    tag = FANIN_INFOTAB(tdeb)[FTGT_PROCDST]%SOLV_THRDNBR;
  if (THREAD_COMM_OFF)
    {
#if (defined EXACT_THREAD)
      tag = FANIN_INFOTAB(tdeb)[FTGT_PROCDST]%SOLV_THRDNBR;
#elif (defined EXACT_TAG)
      tag = FANIN_PRIONUM(tdeb);
#endif
    }

  /* le nombre de paquets est code dans le champ prionum de la 1ere fanin */
  FANIN_PRIONUM(tdeb) = packnbr/2+1;

  /* On recherche la premiere requete disponible */
#ifdef TEST_ISEND
  id_req = send_waitone_fanin(sopalin_data, me);
  ASSERTDBG(id_req<MAX_S_REQUESTS,MOD_SOPALIN);
#endif /* TEST_ISEND */

  /* Envoi des donnees */
  trace_send(thread_data->tracefile,
             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, FANIN_PROCDST(tdeb),
             COMM_FANIN, FANIN_TASKDST(tdeb), thread_data->gtabsize[1],
             &(FANIN_IDTRACE(tdeb)));

#ifdef NO_MPI_TYPE
  thread_data->send_fanin_buffer_size[id_req] = 0;
  for (iter = 0; iter < 2*(packnbr/2+1); iter++) {
    thread_data->send_fanin_buffer_size[id_req] +=
      thread_data->gtabsize[iter]*thread_data->gtabtype[iter];
  }
  MALLOC_INTERN(thread_data->send_fanin_buffer[id_req],
                thread_data->send_fanin_buffer_size[id_req],
                char);
  copied = 0;
  for (iter = 0; iter < 2*(packnbr/2+1); iter++) {
    memcpy(thread_data->send_fanin_buffer[id_req]+copied,
           thread_data->gtaboffs[iter],
           thread_data->gtabtype[iter]*thread_data->gtabsize[iter] );
    copied += thread_data->gtabsize[iter]*thread_data->gtabtype[iter];
  }

#  ifdef TEST_ISEND
  CALL_MPI MPI_Isend(thread_data->send_fanin_buffer[id_req],
                     thread_data->send_fanin_buffer_size[id_req],
                     MPI_CHAR,FANIN_PROCDST(tdeb),tag,
                     PASTIX_COMM,&((thread_data->send_fanin_requests[id_req])));
  TEST_MPI("MPI_Isend");
#  else
  CALL_MPI MPI_Rsend(thread_data->send_fanin_buffer[id_req],
                     thread_data->send_fanin_buffer_size[id_req],
                     MPI_CHAR,FANIN_PROCDST(tdeb),tag,
                     PASTIX_COMM);
  TEST_MPI("MPI_Rsend");
#  endif
#else /* NO_MPI_TYPE */
  CALL_MPI MPI_Type_struct(2*(packnbr/2+1), thread_data->gtabsize,
                           thread_data->gtaboffs,
                           thread_data->gtabtype, &newtype);
  TEST_MPI("MPI_Type_struct");
  CALL_MPI MPI_Type_commit(&newtype);
  TEST_MPI("MPI_Type_commit");
#  ifdef TEST_ISEND
  thread_data->send_fanin_mpitypes[id_req] = newtype;
  CALL_MPI MPI_Isend(MPI_BOTTOM,1,thread_data->send_fanin_mpitypes[id_req],
                     FANIN_PROCDST(tdeb),tag,
                     PASTIX_COMM,&((thread_data->send_fanin_requests[id_req])));
  TEST_MPI("MPI_Isend");
#  else
  CALL_MPI MPI_Rsend(MPI_BOTTOM,1,newtype,FANIN_PROCDST(tdeb),tag,
                     PASTIX_COMM);
  TEST_MPI("MPI_Rsend");
  CALL_MPI MPI_Type_free(&newtype);
  TEST_MPI("MPI_Type_free");
#  endif
#endif /* NO_MPI_TYPE */

  thread_data->send_fanin_target[id_req]       = tdeb;
  thread_data->send_fanin_target_extra[id_req] = extra;

  /* Try to free some buffer */
#ifdef TEST_ISEND
  send_testall_fanin(sopalin_data, me);
#else
  send_free_fanin(sopalin_data, me, id_req);
#endif

  return (packnbr/2+1);
}

/*
 * Function: send_one_block
 *
 * Send the block t to the correct destination.
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    t            - Block Target number
 *
 * Returns:
 *    Number of blocks sent (always 1)
 */
int
send_one_block(Sopalin_Data_t *sopalin_data,
               PASTIX_INT             me,
               PASTIX_INT             t)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  int            tabsize[3];
#ifdef NO_MPI_TYPE
  void *         taboffs[3];
  int            tabtype[3]  = {sizeof(PASTIX_INT),sizeof(PASTIX_INT),sizeof(PASTIX_FLOAT)};
  PASTIX_INT            copied;
  PASTIX_INT            iter;
#else /* NO_MPI_TYPE */
  MPI_Aint       taboffs[3];
  MPI_Datatype   tabtype[3]  = {COMM_INT,COMM_INT,COMM_FLOAT};
  MPI_Datatype   newtype;
#endif /* NO_MPI_TYPE */
  PASTIX_INT tag;
  PASTIX_INT id_req;

  print_debug(DBG_SOPALIN_SEND, "%ld: Send blocktarget\n"
              "%ld: prionum taskdst procdst taskcnt frownum"
              " lrownum fcolnum lcolnum\n"
              "%ld: %ld %ld %ld %ld %ld %ld %ld %ld\n",
              (long)me, (long)me, (long)me,
              (long)BTAG_PRIONUM(t), (long)BTAG_TASKDST(t),
              (long)BTAG_PROCDST(t), (long)BTAG_TASKCNT(t),
              (long)BTAG_FROWNUM(t), (long)BTAG_LROWNUM(t),
              (long)BTAG_FCOLNUM(t), (long)BTAG_LCOLNUM(t));

  /* create MPI type */
  tabsize[0] = BTAGINFO;
  tabsize[1] = BCOFINFO;
  tabsize[2] = (BTAG_LROWNUM(t)-BTAG_FROWNUM(t)+1)*
    (BTAG_LCOLNUM(t)-BTAG_FCOLNUM(t)+1);

  print_debug(DBG_SOPALIN_SEND, "%ld: BTAG_BTAGTAB : %x\n",
              (long)me, (unsigned int)(intptr_t)BTAG_BTAGTAB(t));
  print_debug(DBG_SOPALIN_SEND, "%ld: BTAG_BCOFTAB : %x\n",
              (long)me, (unsigned int)(intptr_t)BTAG_BCOFTAB(t));
  print_debug(DBG_SOPALIN_SEND, "%ld: BTAG_BCOFPTR : %x\n",
              (long)me, (unsigned int)(intptr_t)BTAG_BCOFPTR(t));
  print_debug(DBG_SOPALIN_SEND, "%ld: BTAG_COEFTAB : %x\n",
              (long)me, (unsigned int)(intptr_t)BTAG_COEFTAB(t));

#ifndef NO_MPI_TYPE
  CALL_MPI MPI_Address(BTAG_BTAGTAB(t),&(taboffs[0]));
  TEST_MPI("MPI_Address");
  CALL_MPI MPI_Address(BTAG_BCOFTAB(t),&(taboffs[1]));
  TEST_MPI("MPI_Address");
  CALL_MPI MPI_Address((void *)BTAG_COEFTAB(t),&(taboffs[2]));
  TEST_MPI("MPI_Address");

  CALL_MPI MPI_Type_struct(3,tabsize,taboffs,tabtype,&newtype);
  TEST_MPI("MPI_Type_struct");

  CALL_MPI MPI_Type_commit(&newtype);
  TEST_MPI("MPI_Type_commit");
#endif /* NO_MPI_TYPE */

  /* On recherche la premiere requete disponible */
#ifdef TEST_ISEND
  id_req = send_waitone_block(sopalin_data, me);
  ASSERTDBG(id_req<MAX_S_REQUESTS,MOD_SOPALIN);
#endif /* TEST_ISEND */

  trace_send(thread_data->tracefile,
             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me,
             BTAG_BTAGTAB(t)[BTAG_PROCDST],
             COMM_BLOCK, (SEPFB+BTAG_BTAGTAB(t)[BTAG_TASKDST]), tabsize[2],
             &(BTAG_IDTRACE(t)));

  print_debug(DBG_SOPALIN_SEND, "%ld: send block cpu=%ld thrd=%ld(%ld)\n",
              (long)me, (long)BTAG_PROCDST(t),
              (long)(BTAG_BTAGTAB(t)[BTAG_PROCDST]%SOLV_THRDNBR),
              (long)(BTAG_BTAGTAB(t)[BTAG_PROCDST]));

  if (sopalin_data->sopar->type_comm != API_THREAD_COMM_NBPROC)
    tag = TAG_BLOCK;
  else
    tag = SOLV_THRDNBR + BTAG_BTAGTAB(t)[BTAG_PROCDST]%SOLV_THRDNBR;
  if (THREAD_COMM_OFF)
    {
#if (defined EXACT_THREAD)
      tag = SEPFB+BTAG_BTAGTAB(t)[BTAG_PROCDST]%SOLV_THRDNBR;
#elif (defined EXACT_TAG)
      tag = SEPFB+BTAG_PRIONUM(t);
#endif
    }
#ifdef NO_MPI_TYPE
  taboffs[0] = (void*)BTAG_BTAGTAB(t);
  taboffs[1] = (void*)BTAG_BCOFTAB(t);
  taboffs[2] = (void*)BTAG_COEFTAB(t);

  thread_data->send_block_buffer_size[id_req] = 0;
  for(iter = 0; iter < 3; iter ++){
    thread_data->send_block_buffer_size[id_req] += tabsize[iter]*tabtype[iter];
  }

  MALLOC_INTERN(thread_data->send_block_buffer[id_req],
                thread_data->send_block_buffer_size[id_req],
                char);

  copied = 0;
  for(iter = 0; iter < 3; iter ++)
    {
      memcpy(thread_data->send_block_buffer+copied,
             taboffs[iter], tabsize[iter]*tabtype[iter]);
      copied +=tabsize[iter]*tabtype[iter];
    }
#  ifdef TEST_ISEND
  CALL_MPI MPI_Isend(thread_data->send_block_buffer[id_req],
                     thread_data->send_block_buffer_size[id_req],
                     MPI_CHAR,BTAG_PROCDST(t),tag,
                     PASTIX_COMM,&(thread_data->send_block_requests[id_req]));
  TEST_MPI("MPI_Isend");
#  else
  CALL_MPI MPI_Rsend(thread_data->send_block_buffer[id_req],
                     thread_data->send_block_buffer_size[id_req],
                     MPI_CHAR,BTAG_PROCDST(t),tag,
                     PASTIX_COMM);
  TEST_MPI("MPI_Rsend");
#  endif
#else /* NO_MPI_TYPE */
#  ifdef TEST_ISEND
  CALL_MPI MPI_Isend(MPI_BOTTOM,1,newtype,BTAG_PROCDST(t),tag,
                     PASTIX_COMM,&(thread_data->send_block_requests[id_req]));
  TEST_MPI("MPI_Isend");
#  else
  CALL_MPI MPI_Rsend(MPI_BOTTOM,1,newtype,BTAG_PROCDST(t),tag,
                     PASTIX_COMM);
  TEST_MPI("MPI_Send");
#  endif
#endif /* NO_MPI_TYPE */

  thread_data->send_block_target[id_req] = t;

#ifndef NO_MPI_TYPE
  /* free MPI type */
  CALL_MPI MPI_Type_free(&newtype);
  TEST_MPI("MPI_Type_free");
#endif /* NO_MPI_TYPE */

  /* Try to free some buffer */
#ifdef TEST_ISEND
  send_testall_block(sopalin_data, me);
#else
  send_free_block(sopalin_data, me, id_req);
#endif

  return 1;
}

/*
 * Function: send_all_fanin
 *
 * Send all contribution for a different task on a same destination.
 * Can't be called with API_THREAD_FUNNELED
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    dest         - First fanin number
 *
 */
void
send_all_fanin(Sopalin_Data_t *sopalin_data,
               PASTIX_INT             me,
               PASTIX_INT             dest)
{
  if (THREAD_FUNNELED_ON)
    {
      errorPrintW("API_THREAD_FUNNELED is incompatible with the function"
                  " send_all_fanin");
    }
  else
    {
      SolverMatrix  *datacode    = sopalin_data->datacode;
#ifdef TRACE_SOPALIN
      Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif
      int            flag;

      print_debug(DBG_SOPALIN_SEND, "%ld : --> send_all_fanin\n", (long)me);

      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                       STATE_L2_SENDF, 0);

      flag = 1;
      if (dest != SOLV_PROCNUM)
        {
          MUTEX_LOCK(&(sopalin_data->mutex_queue_fanin[dest]));
          while ( (flag) &&
                  (queueSize(&(sopalin_data->fanintgtsendqueue[dest]))))
            {
              PASTIX_INT t = queueRead(&(sopalin_data->fanintgtsendqueue[dest]));

              print_debug(DBG_SOPALIN_SEND,
                          "send dest %ld fanintarget %ld\n",(long)dest,(long)t);

              /* If the target is ready, we send it */
              /* With -DCOMM_REORDER : Always true  */
              if (!(FANIN_CTRBCNT(t)))
                {
                  double key;

                  t = queueGet2(&(sopalin_data->fanintgtsendqueue[dest]),
                                &key, NULL);

                  ASSERTDBG(FANIN_PROCDST(t) != SOLV_PROCNUM, MOD_SOPALIN);

                  /* send fanin target t */
                  send_one_fanin(sopalin_data, me, t);

                }
              else
                flag = 0;
            }
          MUTEX_UNLOCK(&(sopalin_data->mutex_queue_fanin[dest]));
        }

      trace_end_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                     STATE_L2_SENDF, 0);

      print_debug(DBG_SOPALIN_SEND, "%ld : <-- send_all_fanin\n", (long)me);
    }
  return;
}

/*
 * Function: send_all_block
 *
 * Send all contribution for a same task on a same destination.
 * Can't be called in API_THREAD_FUNNELED
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *
 */
void
send_all_block ( Sopalin_Data_t *sopalin_data,
                 PASTIX_INT             me)
{
  if (THREAD_FUNNELED_ON)
    {
      errorPrintW("API_THREAD_FUNNELED is incompatible with the function"
                  " send_all_block");
    }
  else
    {
      SolverMatrix  *datacode    = sopalin_data->datacode;
#ifdef TRACE_SOPALIN
      Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif
      int            flag,dest;

      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                       STATE_L2_SENDB, 0);

      print_debug(DBG_SOPALIN_SEND, "%ld: send_all_block\n", (long)me);

      for (dest=0;dest<SOLV_PROCNBR;dest++)
        {
          flag = 1;
          if (dest == SOLV_PROCNUM)
            continue;

          MUTEX_LOCK(&(sopalin_data->mutex_queue_block[dest]));
          while ( (flag) &&
                  (queueSize(&(sopalin_data->blocktgtsendqueue[dest]))) )
            {
              PASTIX_INT t = queueRead(&(sopalin_data->blocktgtsendqueue[dest]));

              print_debug(DBG_SOPALIN_SEND,
                          "send dest %ld blocktarget %ld bcofptr %x (ready %x)"
                          " prio %ld task %ld\n",
                          (long)dest, (long)t,
                          (unsigned int)(intptr_t)BTAG_BCOFPTR(t),
                          (unsigned int)(intptr_t)BTAG_COEFTAB(t),
                          (long)BTAG_PRIONUM(t), (long)BTAG_TASKDST(t));

              /* if the target is ready */
              if (BTAG_COEFTAB(t))
                {
                  double key;

                  t = queueGet2(&(sopalin_data->blocktgtsendqueue[dest]),
                                &key, NULL);

                  ASSERTDBG(BTAG_PROCDST(t)!=SOLV_PROCNUM,MOD_SOPALIN);

                  /* send block target t */
                  send_one_block(sopalin_data, me, t);
                }
              else
                flag = 0;
            }
          MUTEX_UNLOCK(&(sopalin_data->mutex_queue_block[dest]));
        }

      trace_end_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                     STATE_L2_SENDB, 0);

    }
  return;
}

/*
 * Function: send_free_fanin
 *
 * Free associated structure to fanin sent.
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
void
send_free_fanin ( Sopalin_Data_t *sopalin_data,
                  PASTIX_INT             me,
                  PASTIX_INT             s_index)
{
#if ((!defined OOC_FTGT) || defined PASTIX_DEBUG )
  SolverMatrix  *datacode    = sopalin_data->datacode;
#endif
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  PASTIX_INT f = 0;
  PASTIX_INT i;

  i = thread_data->send_fanin_target[s_index];
  print_debug(DBG_SOPALIN_SEND, "->send_free_fanin %ld \n", (long)i);
#ifdef OOC_FTGT
  ooc_reset_ftgt(sopalin_data, i, me);
#else
  memFree_null(FANIN_COEFTAB(i));
#endif
#ifndef NO_MPI_TYPE
  /* free MPI type */
  CALL_MPI MPI_Type_free(&(thread_data->send_fanin_mpitypes[s_index]));
  TEST_MPI("MPI_Type_free");
#endif /* NO_MPI_TYPE */

  print_debug(DBG_SOPALIN_ALLOC, "free fanin coeff %x\n",
              (unsigned int)(intptr_t)FANIN_COEFTAB(i));

  STATS_SUB((FANIN_LROWNUM(i)-FANIN_FROWNUM(i)+1)*
            (FANIN_LCOLNUM(i)-FANIN_FCOLNUM(i)+1));

  if (thread_data->send_fanin_target_extra[s_index])
    {
      while((i = thread_data->send_fanin_target_extra[s_index][f]) != -1)
        {
#ifdef OOC_FTGT
          ooc_reset_ftgt(sopalin_data,i,me);
#else
          memFree_null(FANIN_COEFTAB(i));
#endif

          print_debug(DBG_SOPALIN_ALLOC, "free fanin coeff %x\n",
                      (unsigned int)(intptr_t)FANIN_COEFTAB(i));

          STATS_SUB((FANIN_LROWNUM(i)-FANIN_FROWNUM(i)+1)*
                    (FANIN_LCOLNUM(i)-FANIN_FCOLNUM(i)+1));

          f++;
        }

      memFree_null(thread_data->send_fanin_target_extra[s_index]);
      thread_data->send_fanin_target_extra[s_index] = NULL;
    }

#ifdef NO_MPI_TYPE
  memFree_null(thread_data->send_fanin_buffer[s_index]);
  thread_data->send_fanin_buffer[s_index] = NULL;
#endif /* NO_MPI_TYPE */
  print_debug(DBG_SOPALIN_SEND, "<-send_free_fanin %ld \n", (long)i);
}

/*
 * Function: send_free_block
 *
 * Free associated structure to block sent.
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
void
send_free_block(Sopalin_Data_t *sopalin_data,
                PASTIX_INT             me,
                PASTIX_INT             s_index)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  PASTIX_INT            i           = thread_data->send_block_target[s_index];

  BTAG_SENDCNT(i)--;
  ASSERTDBG(BTAG_SENDCNT(i)>=0,MOD_SOPALIN);

  if (BTAG_SENDCNT(i)==0)
    {
      memFree_null(BTAG_COEFTAB(i));
      print_debug(DBG_SOPALIN_ALLOC, "free block coeff %x\n",
                  (unsigned int)(intptr_t)BTAG_COEFTAB(i));

      STATS_SUB((BTAG_LROWNUM(i)-BTAG_FROWNUM(i)+1)*
                (BTAG_LCOLNUM(i)-BTAG_FCOLNUM(i)+1));
    }

#ifdef NO_MPI_TYPE
  memFree_null(thread_data->send_block_buffer[s_index]);
#endif /* NO_MPI_TYPE */
}

/*********************************/
/*
 * Function: send_testall_fanin
 *
 * Test all fanin sent to progress communications
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
void send_testall ( Sopalin_Data_t *sopalin_data, PASTIX_INT me,
                    void (*funcfree)(Sopalin_Data_t*, PASTIX_INT, PASTIX_INT))
{
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  int            i;

  print_debug(DBG_SOPALIN_UPDO, "%ld: test_all_downsend\n", (long)me);

#ifndef PASTIX_TERA
  {
#  if (defined TEST_IRECV) || ((defined TEST_ISEND) && (defined SMP_SOPALIN))
    SolverMatrix  *datacode = sopalin_data->datacode;
#  endif
    MPI_Status    *statuses;
    int           *indices;
    int            outcount = 0;

    if (thread_data->srteststatus == NULL)
      {
        MALLOC_INTERN(thread_data->srteststatus,
                      MAX(MAX_R_REQUESTS, MAX_S_REQUESTS), MPI_Status);
        MALLOC_INTERN(thread_data->srtestindices,
                      MAX(MAX_R_REQUESTS, MAX_S_REQUESTS), int);
      }
    statuses = thread_data->srteststatus;
    indices  = thread_data->srtestindices;

    if (thread_data->maxsrequest_fanin > 0)
      {
        CALL_MPI MPI_Testsome(thread_data->maxsrequest_fanin,
                              thread_data->send_fanin_requests,
                              &outcount, indices, statuses);
        TEST_MPI("MPI_Testsome");

        for(i=0; i<outcount; i++)
          {
            funcfree(sopalin_data, me, indices[i]);
          }
      }
  }
#else /* PASTIX_TERA */
  /* Can be removed in next release, if the first version is ok on TERA10 */
  {
    MPI_Status s_status;
    int        s_flag;

    for(i=0; i<thread_data->maxsrequest_fanin; i++)
      {
        if (!MPI_Request_is_equal(thread_data->send_fanin_requests[i],
                                  MPI_REQUEST_NULL))
          {
            CALL_MPI MPI_Test(&(thread_data->send_fanin_requests[i]),
                              &s_flag, &s_status);
            TEST_MPI("MPI_Test");
            if (s_flag)
              {
                funcfree(sopalin_data, me, i);
              }
          }
      }
  }
#endif /* PASTIX_TERA */
}

void send_testall_fanin(Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  send_testall(sopalin_data, me, send_free_fanin);
  return;
}

/*********************************/
/*
 * Function: send_testall_block
 *
 * Test all block sent to progress communications
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
void send_testall_block ( Sopalin_Data_t *sopalin_data, PASTIX_INT me )
{
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];

  print_debug(DBG_SOPALIN_RECV, "%ld: send_testall_block\n", (long)me);

#ifndef PASTIX_TERA
  {
#  if (defined TEST_IRECV) || ((defined TEST_ISEND) && (defined SMP_SOPALIN))
    SolverMatrix  *datacode    = sopalin_data->datacode;
#  endif
    MPI_Status    *statuses    = thread_data->srteststatus;
    int           *indices     = thread_data->srtestindices;
    int            outcount = 0;
    int            i;

    if (thread_data->srteststatus == NULL)
      {
        MALLOC_INTERN(thread_data->srteststatus,
                      MAX(MAX_R_REQUESTS, MAX_S_REQUESTS), MPI_Status);
        MALLOC_INTERN(thread_data->srtestindices,
                      MAX(MAX_R_REQUESTS, MAX_S_REQUESTS), int);
        statuses = thread_data->srteststatus;
        indices  = thread_data->srtestindices;
      }

    if (thread_data->maxsrequest_block > 0)
      {
        CALL_MPI MPI_Testsome(thread_data->maxsrequest_block,
                              thread_data->send_block_requests,
                              &outcount, indices, statuses);
        TEST_MPI("MPI_Testsome");

        for(i=0; i<outcount; i++)
          {
            send_free_block(sopalin_data, me, indices[i]);
          }
      }
  }
#else
  /* Can be removed in next release, if the first version is ok on TERA10 */
  {
    MPI_Status s_status;
    int        s_flag;
    PASTIX_INT        j;

    for(j=0; j<thread_data->maxsrequest_block; j++)
      {
        if (!MPI_Request_is_equal(thread_data->send_block_requests[j],
                                  MPI_REQUEST_NULL))
          {
            CALL_MPI MPI_Test(&(thread_data->send_block_requests[j]),
                              &s_flag, &s_status);
            TEST_MPI("MPI_Test");
            if (s_flag)
              {
                send_free_block(sopalin_data, me, j);
              }
          }
      }
  }
#endif
  return;
}

/*********************************/
/*
 * Function: send_testall_fab
 *
 * Test all block sent to progress communications
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
void send_testall_fab(Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  send_testall_fanin(sopalin_data, me);
  send_testall_block(sopalin_data, me);
}

/*********************************/
/*
 * Function: send_waitone_fanin
 *
 * Test fanin sent to return an id or wait until one fanin finished.
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
int send_waitone ( Sopalin_Data_t *sopalin_data, PASTIX_INT me,
                   void (*funcfree)(Sopalin_Data_t*, PASTIX_INT, PASTIX_INT) )
{
#if (defined TEST_ISEND) && (defined SMP_SOPALIN)
  SolverMatrix  *datacode    = sopalin_data->datacode;
#endif
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  MPI_Status     status;
  int            flag = 0;
  int            req  = 0;

  print_debug(DBG_SOPALIN_RECV, "%ld: send_waitone_fanin\n",(long)me);

  while((!flag) && (thread_data->maxsrequest_fanin > 0))
    {
      CALL_MPI MPI_Testany(thread_data->maxsrequest_fanin,
                           thread_data->send_fanin_requests,
                           &req, &flag, &status);
      TEST_MPI("MPI_Testany");

      if ( flag )
        {
          if (req != MPI_UNDEFINED)
            {
              funcfree(sopalin_data, me, req);
              return req;
            }
          else
            /* Case where all requests are finished */
            return 0;
        }
      else
        {
          if (thread_data->maxsrequest_fanin < MAX_S_REQUESTS)
            {
              req = thread_data->maxsrequest_fanin;
              thread_data->maxsrequest_fanin++;
              return req;
            }
        }
    }

  /* On n'est pas rentré dans la boucle */
  thread_data->maxsrequest_fanin++;
  return req;
}

int send_waitone_fanin(Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  return send_waitone(sopalin_data, me, send_free_fanin);
}

/*********************************/
/*
 * Function: send_waitone_block
 *
 * Tests send block. Must be called only when all requests are used.
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
int send_waitone_block ( Sopalin_Data_t *sopalin_data, PASTIX_INT me )
{
#if (defined TEST_ISEND) && (defined SMP_SOPALIN)
  SolverMatrix  *datacode    = sopalin_data->datacode;
#endif
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  MPI_Status     status;
  int            flag = 0;
  int            req  = 0;

  print_debug(DBG_SOPALIN_RECV, "%ld: send_waitone_block\n",(long)me);

  while((!flag) && (thread_data->maxsrequest_block > 0))
    {
      CALL_MPI MPI_Testany(thread_data->maxsrequest_block,
                           thread_data->send_block_requests,
                           &req, &flag, &status);
      TEST_MPI("MPI_Testany");

      if (flag)
        {
          if (req != MPI_UNDEFINED)
            {
              send_free_block(sopalin_data, me, req);
              return req;
            }
          else
            /* Case where all requests are finished */
            return 0;
        }
      else
        {
          if (thread_data->maxsrequest_block < MAX_S_REQUESTS)
            {
              req = thread_data->maxsrequest_block;
              thread_data->maxsrequest_block++;
              return req;
            }
        }
    }

  thread_data->maxsrequest_block++;
  return req;
}

/*********************************/
/*
 * Function: send_waitall_fab
 *
 * Wait for all pending communications (fanin and block).
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
void send_waitall_fab(Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  MPI_Status     s_status;
  int            i;

#ifdef FORCE_CONSO
  PASTIX_INT        nb_envois_fanin = 0;
  PASTIX_INT        nb_envois_block = 0;
  int        s_flag          = 0;

  while((nb_envois_fanin != thread_data->maxsrequest_fanin) ||
        (nb_envois_block != thread_data->maxsrequest_block))
    {
      nb_envois_fanin = 0;
      for(i=0; i<thread_data->maxsrequest_fanin; i++)
        {
          if (!MPI_Request_is_equal(thread_data->send_fanin_requests[i], MPI_REQUEST_NULL))
            {
              CALL_MPI MPI_Test(&(thread_data->send_fanin_requests[i]),
                                &s_flag, &s_status);
              TEST_MPI("MPI_Test");
              if (s_flag)
                {
                  send_free_fanin(sopalin_data, me, i);
                  nb_envois_fanin++;
                }
            }
          else
            nb_envois_fanin++;
        }

      nb_envois_block = 0;
      for(i=0; i<thread_data->maxsrequest_block; i++)
        {
          if (!MPI_Request_is_equal(thread_data->send_block_requests[i],
                                    MPI_REQUEST_NULL))
            {
              CALL_MPI MPI_Test(&(thread_data->send_block_requests[i]),
                                &s_flag, &s_status);
              TEST_MPI("MPI_Test");
              if (s_flag)
                {
                  send_free_block(sopalin_data, me, i);
                  nb_envois_block++;
                }
            }
          else
            nb_envois_block++;
        }
    }
#endif

  for (i=0;i<thread_data->maxsrequest_fanin;i++)
    if (!MPI_Request_is_equal(thread_data->send_fanin_requests[i],
                              MPI_REQUEST_NULL))
      {
        CALL_MPI MPI_Wait(&thread_data->send_fanin_requests[i], &s_status);
        TEST_MPI("MPI_Wait");

        send_free_fanin(sopalin_data, me, i);
      }

  for (i=0;i<thread_data->maxsrequest_block;i++)
    if (!MPI_Request_is_equal(thread_data->send_block_requests[i],
                              MPI_REQUEST_NULL))
      {
        CALL_MPI MPI_Wait(&thread_data->send_block_requests[i], &s_status);
        TEST_MPI("MPI_Wait");

        send_free_block(sopalin_data, me, i);
      }
}

/*********************************/
/*
 * Function: rcsd_testall_fab
 *
 * Launch threads for solving step.
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
void rcsd_testall_fab(Sopalin_Data_t *sopalin_data, PASTIX_INT me)
{
#ifdef TEST_ISEND
  /* Test des envois */
  send_testall_fab(sopalin_data, me);
#endif

  /* Test des receptions */
  recv_testall_fab(sopalin_data, me);
}


/*
 * Function: sendrecv_smp
 *
 * fonction de réception des comms dans la version threads séparés
 *
 */
void* sendrecv_smp ( void *arg )
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  if (THREAD_COMM_ON)
    {
      SolverMatrix     *datacode     = sopalin_data->datacode;
      MPI_Comm          pastix_comm  = PASTIX_COMM;
      int               type_thcomm  = sopalin_data->sopar->type_comm;
      int               me           = argument->me;
#ifdef TRACE_SOPALIN
      Thread_Data_t    *thread_data;
#endif
      void            **receive_buffer;
      MPI_Request      *request;
      MPI_Status        status;
      int               type_comm;
      int               tag_fanin, tag_block;
      int               nbrequest   = 1;
      int               nbrequesttot= 1;
      PASTIX_INT               size[2];
      int               i;
      int               init;
      int               flag, wait;
      int               nbsend, nbsend_fanin, nbsend_block;
      int               nbrecv, nbrecv_block;
      PASTIX_INT               save_ftgtcnt = -1;
      PASTIX_INT               ftgt, key;
      double            dest;
      int               nb_proc_end = 1;


      print_debug(DBG_SOPALIN_THREADCOMM, "%d - %d : --> SendRecv\n",
                  (int)SOLV_PROCNUM, (int)me);

      if (SOLV_PROCNBR == 1) goto end;

      /* Allocation de la structure de données spécifique a ce thread */
      init = 0;
      if (THREAD_FUNNELED_ON)
        {
          init = init | INIT_SEND;
        }

      sopalin_init_smp(sopalin_data, me, 1, init);
#ifdef TRACE_SOPALIN
      thread_data = sopalin_data->thread_data[me];
#endif
      /***********************************/
      /*         Reception               */
      /***********************************/

      /* Taille des buffers de réception */
      size[0] = PACKMAX*(sizeof(PASTIX_INT)*MAXINFO)+PACKAREA*sizeof(PASTIX_FLOAT);
      size[1] = sizeof(PASTIX_INT)*(BTAGINFO+BCOFINFO)+sizeof(PASTIX_FLOAT)*SOLV_BPFTMAX;

#ifdef THREAD_COMM_MULTIPLE
      nbrequest = MAX(SOLV_PROCNBR-1, 1);
#endif
      nbrequesttot = 2*nbrequest;

      /* Allocation des buffers de reception et des requetes */
      MALLOC_INTERN(receive_buffer, nbrequesttot, void*);
      MALLOC_INTERN(request       , nbrequesttot, MPI_Request);

      for(i=0; i< nbrequest; i++)
        {
          MALLOC_INTERN(receive_buffer[2*i],   size[0], char);
          MALLOC_INTERN(receive_buffer[2*i+1], size[1], char);
        }

      /* Initialisation des requêtes */
      for (i=0; i<nbrequesttot; i++)
        request[i] = MPI_REQUEST_NULL;

      /* Exact thread ou Tag classique
       ( cas 1 thread de comm par thread de calcul )*/
      if (type_thcomm == 3)
        {
          tag_fanin = me-SOLV_THRDNBR;
          tag_block = me;
        }
      else
        {
          tag_fanin = TAG_FANIN;
          tag_block = TAG_BLOCK;
        }

      /* Initialisation des comms persistantes en réception */
      /* Pas performant, pourrait etre vire */
#ifdef THREAD_COMM_MULTIPLE
      /* Proc de rang inferieur au proc local */
      for(i=0; i < SOLV_PROCNUM; i++){
        CALL_MPI MPI_Recv_init(receive_buffer[2*i],size[0],MPI_BYTE,
                               i,tag_fanin,pastix_comm,&request[2*i]);
        TEST_MPI("MPI_Recv_init");

        CALL_MPI MPI_Recv_init(receive_buffer[2*i+1],size[1],MPI_BYTE,
                               i,tag_block,pastix_comm,&request[2*i+1]);
        TEST_MPI("MPI_Recv_init");
      }
      /* Proc de rang supérieur au proc local */
      for(i=SOLV_PROCNUM+1; i<SOLV_PROCNBR; i++){
        CALL_MPI MPI_Recv_init(receive_buffer[2*(i-1)],size[0],MPI_BYTE,
                               i,tag_fanin,pastix_comm,&request[2*(i-1)]);
        TEST_MPI("MPI_Recv_init");

        CALL_MPI MPI_Recv_init(receive_buffer[2*(i-1)+1],size[1],MPI_BYTE,
                               i,tag_block,pastix_comm,&request[2*(i-1)+1]);
        TEST_MPI("MPI_Recv_init");
      }
#else
      CALL_MPI MPI_Recv_init(receive_buffer[0], size[0], MPI_BYTE,
                             MPI_ANY_SOURCE, tag_fanin, pastix_comm,
                             &request[0]);
      TEST_MPI("MPI_Recv_init");

      CALL_MPI MPI_Recv_init(receive_buffer[1], size[1], MPI_BYTE,
                             MPI_ANY_SOURCE, tag_block, pastix_comm,
                             &request[1]);
      TEST_MPI("MPI_Recv_init");
#endif

      /***********************************/
      /*       Boucle Principale         */
      /***********************************/

      print_debug(DBG_SOPALIN_THREADCOMM,
                  "%d - %d : Boucle Emission-Reception\n",
                  (int)SOLV_PROCNUM, (int)me);

      /* Lancement des communications */
      CALL_MPI MPI_Startall(nbrequesttot, request);
      TEST_MPI("MPI_Startall");

      if (THREAD_FUNNELED_OFF)
        {
          /*
           * Pas de test sur le nombre de comm restant a recevoir
           * car ca pourrait poser pb en multi-thread de comm
           */
          while(1){

            trace_begin_task(thread_data->tracefile,
                             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                             STATE_WAITREM, 1);

            /* On attend une comm */
            CALL_MPI MPI_Waitany(nbrequesttot, request, &type_comm, &status);
            TEST_MPI("MPI_Waitany");

            print_debug(DBG_SOPALIN_THREADCOMM,
                        "%d - %d : 1 reception / nb_proc_end %d\n",
                        (int)SOLV_PROCNUM, (int)me, nb_proc_end);

            /* Recuperation des comms de fin */

            trace_begin_task(thread_data->tracefile,
                             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                             STATE_COMPUTE, 1);

            if ((status.MPI_TAG == tag_fanin)
                && ( ((int *)(receive_buffer[type_comm]))[0] == -1 ))
              {
                nb_proc_end++;
                /* fprintf(stderr, "%d - %d : reçu msg end de %d (%d, %d)\n",
                 SOLV_PROCNUM, me, status.MPI_SOURCE,
                 status.MPI_TAG, type_comm);*/
                if (nb_proc_end < SOLV_PROCNBR)
                  {
                    CALL_MPI MPI_Start(&request[type_comm]);
                    TEST_MPI("MPI_Start");
                    continue;
                  }
                else
                  break;
              }

            /* On fait le calcul associé */
            if (type_comm%2 == 1)
              recv_handle_block(sopalin_data, me,
                                receive_buffer[type_comm],
                                status, 0);
            else
              recv_handle_fanin(sopalin_data, me,
                                receive_buffer[type_comm],
                                status, 0);

            /* On relance l'attente sur une comm */
            CALL_MPI MPI_Start(&request[type_comm]);
            TEST_MPI("MPI_Start");
          }
        }
      else
        {
          /* THREAD_FUNNELED_ON */
          wait   = 0;
          nbsend = 2;
          nbsend_fanin = SOLV_FTGTNBR;
          nbsend_block = SOLV_BTGSNBR;
          nbrecv_block = SOLV_BTGRNBR;
          save_ftgtcnt = SOLV_FTGTCNT;
          nbrecv = SOLV_FTGTCNT + nbrecv_block;

          if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
            {
              fprintf(stdout, OUT2_FUN_STATS,
                      (long)SOLV_PROCNUM, (long)nbrecv,
                      (long)(nbsend_block+nbsend_fanin));
            }

          trace_begin_task(thread_data->tracefile,
                           SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                           STATE_WAITREM, 1);

          while( nbsend || nbrecv ) {

            /*
             * Attente si rien à faire
             */
            if (wait)
              {
                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                COND_TIMEWAIT(&(sopalin_data->cond_comm),
                              &(sopalin_data->mutex_comm));
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
              }
            wait = 1;

            /*
             * Réception des données
             */
            CALL_MPI MPI_Testany(nbrequesttot, request, &type_comm,
                                 &flag, &status);
            TEST_MPI("MPI_Testany");

            if(flag && (type_comm != MPI_UNDEFINED))
              {
                wait = 0;

                trace_begin_task(thread_data->tracefile,
                                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                 STATE_COMPUTE, 1);

                /* On fait le calcul associé */
                if (type_comm%2)
                  {
                    recv_handle_block(sopalin_data, me,
                                      receive_buffer[type_comm],
                                      status, 0);
                    nbrecv_block--;
                  }
                else
                  recv_handle_fanin(sopalin_data, me,
                                    receive_buffer[type_comm],
                                    status, 0);

                nbrecv = nbrecv_block + SOLV_FTGTCNT;

                /* On relance l'attente sur une comm */
                CALL_MPI MPI_Start(&request[type_comm]);
                TEST_MPI("MPI_Start");

                print_debug(DBG_SOPALIN_THREADCOMM,
                            "%d : %d Reception solv_ftgtcnt %d\n",
                            (int)SOLV_PROCNUM, (int)me, (int)SOLV_FTGTCNT);

                trace_begin_task(thread_data->tracefile,
                                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                 STATE_WAITREM, 1);
              } /* Fin Boucle réception */

            /*
             * Progression des envois
             */
            send_testall_fab(sopalin_data, me);

            /*
             * Envoi des données
             */
            if (nbsend > 0)
              {
                trace_begin_task(thread_data->tracefile,
                                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                 STATE_COMPUTE, 2);

                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                if(queueSize(sopalin_data->sendqueue))
                  {
                    ftgt = queueGet2(sopalin_data->sendqueue, &dest, &key);

                    print_debug(DBG_FUNNELED,
                                "%d-%d C : ftgt %ld / dest %ld / key %ld\n",
                                (int)SOLV_PROCNUM, (int)me,
                                (long)ftgt, (long)dest,
                                (long)key);
                    print_debug(DBG_FUNNELED,
                                "%d-%d C :"
                                " ftgt %ld / dest %ld / key %ld / task %ld\n",
                                (int)SOLV_PROCNUM, (int)me, (long)ftgt,
                                (long)FANIN_PROCDST(ftgt),
                                (long)FANIN_PRIONUM(ftgt),
                                (long)FANIN_TASKDST(ftgt));
                    print_debug(DBG_FUNNELED,
                                "%d-%d C : fanin %ld / block %ld\n",
                                (int)SOLV_PROCNUM, (int)me, (long)nbsend_fanin,
                                (long)nbsend_block);
                    if (dest > 0)
                      nbsend_fanin -= send_one_fanin(sopalin_data,
                                                     me, ftgt);
                    else
                      nbsend_block -= send_one_block(sopalin_data,
                                                     me, ftgt);

                    wait = 0;
                  }
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));

                nbsend = nbsend_fanin + nbsend_block;

                trace_begin_task(thread_data->tracefile,
                                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                 STATE_WAITREM, 2);
              }
          }
        }

      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                       STATE_IDLE, 1);

      print_debug(DBG_SOPALIN_THREADCOMM,
                  "%d - %d : FIN Boucle Emission-Reception\n",
                  (int)SOLV_PROCNUM, (int)me);


      /***********************************/
      /*         Reception               */
      /***********************************/

      /* Annulation des comms inutiles */
      /* et Liberation de la memoire */
      for(i=0; i<2*nbrequest; i++)
        {
          int fflag;
          if (THREAD_FUNNELED_OFF)
            {
              if (i == type_comm)
                {
                  memFree_null(receive_buffer[i]);
                  continue;
                }
            }
          /* Annulation des comms */
          CALL_MPI MPI_Cancel(&request[i]);
          TEST_MPI("MPI_Cancel");
          /* Test pour rendre inactive la requete */
          CALL_MPI MPI_Test(&request[i], &fflag, &status);
          TEST_MPI("MPI_Test");
          /* Liberation de la requete persistante */
          CALL_MPI MPI_Request_free(&request[i]);
          TEST_MPI("MPI_Request_free");
          /* Liberation du buffer */
          memFree_null(receive_buffer[i]);
        }

      memFree_null(receive_buffer);
      memFree_null(request);

      /***********************************/
      /*         Emission                */
      /***********************************/
      if (THREAD_FUNNELED_ON)
        {
          /* Attente de la fin des communications en envoi */
          send_waitall_fab(sopalin_data, me);

          CALL_MPI MPI_Barrier(PASTIX_COMM);
          TEST_MPI("MPI_Barrier");

          /* Restoration de SOLV_FTGTCNT */
          SOLV_FTGTCNT = save_ftgtcnt;
        }

      sopalin_clean_smp(sopalin_data, me);

      end:
      if (me == SOLV_THRDNBR)
        {
          MUTEX_LOCK(&(sopalin_data->mutex_comm));
          sopalin_data->step_comm = COMMSTEP_FACTOEND;
          print_debug(DBG_THCOMM, "%s:%d FACTOEND\n", __FILE__, __LINE__);
          MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
          pthread_cond_broadcast(&(sopalin_data->cond_comm));
        }
      print_debug(DBG_SOPALIN_THREADCOMM,
                  "%d - %d : <-- SendRecv\n", (int)SOLV_PROCNUM, (int)me);

      return 0;
    }
  else
    {
      errorPrint("sendrecv_smp sould not be called in THREAD_COMM mode.\n");
      return NULL;
    }
}
#endif /* FORCE_NOMPI */
