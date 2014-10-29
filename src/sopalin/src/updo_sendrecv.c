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
/******************************************************************************
 * File: updo_sendrecv.c                                                      *
 *                                                                            *
 * Implement the communications routines for updown step                      *
 *                                                                            *
 ******************************************************************************/

#include "updo_sendrecv.h"

/******************************************/
/*                                        */
/*              DOWN                      */
/*                                        */
/******************************************/
/*
 * Backward substitutions routines are only available in MPI version.
 * They don't act on the code contrary to forward substitutions
 * routines in sequential mode.
 */

#ifndef FORCE_NOMPI

void send_free_down ( Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT s_index )
{
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#  ifndef NO_MPI_TYPE
  SolverMatrix  *datacode    = sopalin_data->datacode;

  memFree_null(FANIN_COEFTAB(thread_data->send_fanin_target[s_index]));
  memFree_null(thread_data->send_fanin_infotab[s_index]);
  CALL_MPI MPI_Type_free(&(thread_data->send_fanin_mpitypes[s_index]));
  TEST_MPI("MPI_Type_free");
#  else /* NO_MPI_TYPE */
  memFree_null(thread_data->send_fanin_buffer[s_index]);
#  endif /* NO_MPI_TYPE */
  thread_data->send_fanin_requests[s_index] = MPI_REQUEST_NULL;
}

void send_free_up ( Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT s_index )
{
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#  ifndef NO_MPI_TYPE
  memFree_null(thread_data->send_fanin_infotab[s_index]);
  CALL_MPI MPI_Type_free(&(thread_data->send_fanin_mpitypes[s_index]));
  TEST_MPI("MPI_Type_free");
#  else /* NO_MPI_TYPE */
  memFree_null(thread_data->send_fanin_buffer[s_index]);
#  endif /* NO_MPI_TYPE */
  thread_data->send_fanin_requests[s_index] = MPI_REQUEST_NULL;
}

void send_testall_down ( Sopalin_Data_t *sopalin_data, PASTIX_INT me )
{
  API_CALL(send_testall)(sopalin_data, me, API_CALL(send_free_down));
  return;
}

void send_testall_up ( Sopalin_Data_t *sopalin_data, PASTIX_INT me )
{
  API_CALL(send_testall)(sopalin_data, me, API_CALL(send_free_up));
  return;
}

void send_waitall ( Sopalin_Data_t *sopalin_data, PASTIX_INT me,
                    void (*funcfree)(Sopalin_Data_t*, PASTIX_INT, PASTIX_INT) )
{
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  MPI_Status     s_status;
  int            i;

#  ifdef FORCE_CONSO
  int            s_flag;
  int            nb_envois_fanin = 0;
#  endif /* FORCE_CONSO */

  print_debug(DBG_SOPALIN_UPDO, "%ld: %s", (long)me, __func__);

#  ifdef FORCE_CONSO
  while(nb_envois_fanin != thread_data->maxsrequest_fanin)
    {
      nb_envois_fanin = 0;
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
                  nb_envois_fanin++;
                }
            }
          else
            nb_envois_fanin++;
        }
    }
#  else /* not FORCE_CONSO */
  for (i=0; i<thread_data->maxsrequest_fanin; i++)
    if (!MPI_Request_is_equal(thread_data->send_fanin_requests[i],
                              MPI_REQUEST_NULL))
      {
        CALL_MPI MPI_Wait(&thread_data->send_fanin_requests[i],
                          &s_status);
        TEST_MPI("MPI_Wait");

        funcfree(sopalin_data, me, i);
      }
#  endif /* not FORCE_CONSO */
}

void send_waitall_down ( Sopalin_Data_t *sopalin_data, PASTIX_INT me )
{
  API_CALL(send_waitall)(sopalin_data, me, API_CALL(send_free_down));
  return;
}

void send_waitall_up ( Sopalin_Data_t *sopalin_data, PASTIX_INT me )
{
  API_CALL(send_waitall)(sopalin_data, me, API_CALL(send_free_up));
  return;
}

int send_waitone_down ( Sopalin_Data_t *sopalin_data, PASTIX_INT me )
{
  return API_CALL(send_waitone)(sopalin_data, me, API_CALL(send_free_down));
}

int send_waitone_up ( Sopalin_Data_t *sopalin_data, PASTIX_INT me )
{
  return API_CALL(send_waitone)(sopalin_data, me, API_CALL(send_free_up));
}

/*
 * Function: updo_down_recv
 *
 * Compute the contribution receipts in the updo_buffer for the
 * backward substitution and unlock tasks which waited for it.
 *
 * Parameters:
 *   sopalin_data - Sopalin_data structure to clean
 *   updo_buffer  - Receipt buffer
 *   status       - MPI communication status
 *   me           - Thread Id
 *
 * Returns:
 *   void
 */
void updo_down_recv ( Sopalin_Data_t *sopalin_data,
                      void           *updo_buffer,
                      MPI_Status      status,
                      PASTIX_INT             me )
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  PASTIX_FLOAT         *gb, *gc;
  PASTIX_INT           *infotab;
  PASTIX_INT            size;
#  ifdef TRACE_SOPALIN
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#  endif /* TRACE_SOPALIN */

  infotab = ((PASTIX_INT*)(updo_buffer));

  trace_recv ( thread_data->tracefile,
               SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me,
               status.MPI_SOURCE, COMM_DOWN, infotab[2], 0, infotab[5]);
  trace_begin_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                   STATE_L2_RECVDOWN, infotab[2]);

  print_debug(DBG_SOPALIN_DOWN, "Recoit infotab %d %d %d\n",
              (int)infotab[0], (int)infotab[1], (int)infotab[2]);


  gb   = ((PASTIX_FLOAT*)(updo_buffer)) + UPDOWN_SIZETAB*sizeof(PASTIX_INT)/sizeof(PASTIX_FLOAT);
  size = infotab[1] - infotab[0] + 1;
  gc   =&(UPDOWN_SM2XTAB[UPDOWN_SM2XIND(infotab[2])+
                         infotab[0] - SYMB_FCOLNUM(infotab[2])]);

#  ifdef PASTIX_DUMP_SOLV_COMM
  {
    FILE *file;
    char  filename[30];
    int   i;

    sprintf(filename, "recv/down_%02ld_%02ld_%05ld",
            (long)status.MPI_SOURCE,
            (long)SOLV_PROCNUM,
            (long)infotab[5]);

    file = fopen(filename, "a+");
    for (i=0; i<UPDOWN_SIZETAB; i++)
      fprintf(file, "%ld ", (long)infotab[i]);
    fprintf(file, "\n");
    for (i=0; i<size; i++)
      fprintf(file, "%lf ", gb[i]);
    fprintf(file, "\n");
    fclose(file);
  }
#  endif /* PASTIX_DUMP_SOLV_COMM */

  /* Xi=Xi-Y */
  MUTEX_LOCK(&(sopalin_data->mutex_task[infotab[2]]));
#  ifdef MULT_SMX
  SOPALIN_GESM("N", "N", size, UPDOWN_SM2XNBR, fun,
               gb, size, gc, UPDOWN_SM2XSZE);
#  else /* MULT_SMX */
  SOPALIN_AXPY(size, -fun, gb, iun, gc, iun);
#  endif /* MULT_SMX */

  trace_end_task(thread_data->tracefile,
                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                 STATE_L2_RECVDOWN, infotab[2]);

  UPDOWN_CTRBCNT(infotab[2])-=infotab[3];
  UPDOWN_MSGCNT(infotab[2])--;

  if (!UPDOWN_CTRBCNT(infotab[2]))
    {
      MUTEX_UNLOCK(&(sopalin_data->mutex_task[infotab[2]]));

#  ifdef DEP_SMX
      queueAdd(&cblreadyqueue,infotab[2], (double)(TASK_PRIONUM(infotab[2])));
#  elif (defined PASTIX_DYNSCHED)
      {
        int threadid = TASK_THREADID(infotab[2]);

        MUTEX_LOCK(&(sopalin_data->tasktab_mutex[threadid]));
        queueAdd(&(sopalin_data->taskqueue[threadid]),
                 infotab[2], (double)(TASK_PRIONUM(infotab[2])));
        MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[threadid]));
        pthread_cond_broadcast(&(sopalin_data->tasktab_cond[threadid]));
      }
#  else /* not DEP_SMX nor PASTIX_DYNSCHED */
      pthread_cond_broadcast(&(sopalin_data->cond_task[infotab[2]]));
#  endif /* not DEP_SMX nor PASTIX_DYNSCHED */
    }
  else
    MUTEX_UNLOCK(&(sopalin_data->mutex_task[infotab[2]]));
}


/*
 * Function: updo_down_send
 *
 * Send the contribution to the other process.
 *
 * Parameters:
 *   sopalin_data - Sopalin_data structure to clean
 *   me           - Thread Id
 *   i            -
 *   j            -
 *
 * Returns:
 *   void
 */
void updo_down_send ( Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT i, PASTIX_INT j )
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#  if (defined PASTIX_UPDO_ISEND) || (defined NO_MPI_TYPE)
  PASTIX_INT            id_req = 0;
#  endif /* PASTIX_UPDO_ISEND || NO_MPI_TYPE */
  int           *tabsize;
#  ifdef NO_MPI_TYPE
  void         **taboffs;
  int           *tabtype;
  char          *send_buffer;
  int            send_buffer_size;
#  else /* NO_MPI_TYPE */
  MPI_Aint      *taboffs;
  MPI_Datatype  *tabtype;
  MPI_Datatype   newtype;
#  endif /* NO_MPI_TYPE */

  PASTIX_INT            sizeb, tag;
  PASTIX_INT            iter;
#  if (!defined PASTIX_UPDO_ISEND) || (defined NO_MPI_TYPE)
  PASTIX_INT            infotab[UPDOWN_SIZETAB];
#  else  /* !(!PASTIX_UPDO_ISEND || NO_MPI_TYPE) */
  PASTIX_INT           *infotab;

  MALLOC_INTERN(infotab, UPDOWN_SIZETAB, PASTIX_INT);
#  endif  /* !(!PASTIX_UPDO_ISEND || NO_MPI_TYPE) */
  tabsize = thread_data->gtabsize;
  taboffs = thread_data->gtaboffs;
  tabtype = thread_data->gtabtype;

  /*********************************/
  /*   Data to send preparation    */
  /*********************************/

  /* En-tête du message */
  for (iter = 0; iter < UPDOWN_SIZETAB; iter++)
    infotab[iter] = 0;
  infotab[0] = FANIN_FCOLNUM(SOLV_FTGTIND(j));
  infotab[1] = FANIN_LCOLNUM(SOLV_FTGTIND(j));
  infotab[2] = FANIN_CBLKDST(SOLV_FTGTIND(j));
  infotab[3] = FANIN_CTRBNBR(SOLV_FTGTIND(j));

  sizeb = infotab[1]-infotab[0]+1;

  tabsize[0] = UPDOWN_SIZETAB;
  tabsize[1] = UPDOWN_SM2XNBR * sizeb;

  /* Envoi de la contribution */
#  ifdef NO_MPI_TYPE
  tabtype[0] = sizeof(PASTIX_INT);
  tabtype[1] = sizeof(PASTIX_FLOAT);
  taboffs[0] = infotab;
  taboffs[1] = FANIN_COEFTAB(SOLV_FTGTIND(j));
  send_buffer_size = tabsize[0] * tabtype[0] +
    tabsize[1] * tabtype[1];
#  else /* NO_MPI_TYPE */
  tabtype[0] = COMM_INT;
  tabtype[1] = COMM_FLOAT;

  CALL_MPI MPI_Address(infotab,&(taboffs[0]));
  TEST_MPI("MPI_Address");
  CALL_MPI MPI_Address(FANIN_COEFTAB(SOLV_FTGTIND(j)),&(taboffs[1]));
  TEST_MPI("MPI_Address");
  CALL_MPI MPI_Type_struct(2,tabsize,taboffs,tabtype,&newtype);
  TEST_MPI("MPI_Type_struct");
  CALL_MPI MPI_Type_commit(&newtype);
  TEST_MPI("MPI_Type_commit");
#  endif /* NO_MPI_TYPE */

#  ifdef PASTIX_DUMP_SOLV_COMM
  {
    static int send;
    FILE *file;
    char  filename[30];
    int   i;

    infotab[5] = send++;
    sprintf(filename, "send/down_%02ld_%02ld_%05ld",
            (long)SOLV_PROCNUM,
            (long)FANIN_PROCDST(SOLV_FTGTIND(j)),
            (long)infotab[5]);

    file = fopen(filename, "a+");
    for (i=0; i<UPDOWN_SIZETAB; i++)
      fprintf(file, "%ld ", (long)infotab[i]);
    fprintf(file, "\n");
    for (i=0; i<sizeb; i++)
      fprintf(file, "%lf ", FANIN_COEFTAB(SOLV_FTGTIND(j))[i]);
    fprintf(file, "\n");
    fclose(file);
  }
#  endif /* PASTIX_DUMP_SOLV_COMM */

  tag = TAG_DOWN;
  if (THREAD_COMM_OFF)
    {
#  if (defined EXACT_TAG)
      tag = infotab[2];
#  elif (defined EXACT_THREAD)
      tag = FANIN_INFOTAB(SOLV_FTGTIND(j))[FTGT_PROCDST]%SOLV_THRDNBR;
#  endif /* EXACT_THREAD */
    }
  /*********************************/
  /*   Fin de précedents envois    */
  /*********************************/
#  ifdef PASTIX_UPDO_ISEND
  id_req = API_CALL(send_waitone_down)(sopalin_data, me);
  ASSERTDBG(id_req<MAX_S_REQUESTS, MOD_SOPALIN);
#  endif /* PASTIX_UPDO_ISEND */

  /*********************************/
  /*    Send                       */
  /*********************************/
  print_debug(DBG_SOPALIN_DOWN, "Envoi Contribution %d vers %d (%d,%d,%d)\n",
              (int)j, (int)FANIN_PROCDST(SOLV_FTGTIND(j)),
              (int)infotab[0], (int)infotab[1], (int)infotab[2]);

  trace_send(thread_data->tracefile,
             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me,
             FANIN_PROCDST(SOLV_FTGTIND(j)),
             COMM_DOWN, infotab[2], tabsize[0]+tabsize[1], &(infotab[5]));

#  ifdef NO_MPI_TYPE
  MALLOC_INTERN(thread_data->send_fanin_buffer[id_req], send_buffer_size, char);
  send_buffer = thread_data->send_fanin_buffer[id_req];
  memcpy(send_buffer,                         taboffs[0],
         tabsize[0]*tabtype[0]);
  memcpy(send_buffer+(tabsize[0]*tabtype[0]), taboffs[1],
         tabsize[1]*tabtype[1]);
  COMM_CLOCK_START;
#    ifndef PASTIX_UPDO_ISEND
  CALL_MPI MPI_Send(send_buffer, send_buffer_size, MPI_CHAR,
                    FANIN_PROCDST(SOLV_FTGTIND(j)), tag,
                    PASTIX_COMM);
  TEST_MPI("MPI_Send");

  memFree_null(thread_data->send_fanin_buffer[id_req]);
#    else /* not PASTIX_UPDO_ISEND */
  CALL_MPI MPI_Isend(send_buffer, send_buffer_size, MPI_CHAR,
                     FANIN_PROCDST(SOLV_FTGTIND(j)), tag,
                     PASTIX_COMM, &(thread_data->send_fanin_requests[id_req]));
  TEST_MPI("MPI_Isend");
#    endif /* PASTIX_UPDO_ISEND */
  COMM_CLOCK_STOP;
#    ifdef OOC_FTGT
  ooc_reset_ftgt(sopalin_data, SOLV_FTGTIND(j), me);
#    else /* not OOC_FTGT */
  memFree_null(FANIN_COEFTAB(SOLV_FTGTIND(j)));
#    endif /* not OOC_FTGT */

#  else /* NO_MPI_TYPE */

  COMM_CLOCK_START;
#    ifndef PASTIX_UPDO_ISEND
  CALL_MPI MPI_Send(MPI_BOTTOM, 1, newtype,
                    FANIN_PROCDST(SOLV_FTGTIND(j)), tag,
                    PASTIX_COMM);
  TEST_MPI("MPI_Send");
  CALL_MPI MPI_Type_free(&newtype);
  TEST_MPI("MPI_Type_free");
#      ifdef OOC_FTGT
  ooc_reset_ftgt(sopalin_data,SOLV_FTGTIND(j),me);
#      else /* not OOC_FTGT */
  memFree_null(FANIN_COEFTAB(SOLV_FTGTIND(j)));
#      endif /* not OOC_FTGT */
#    else /* not PASTIX_UPDO_ISEND */
  thread_data->send_fanin_mpitypes[id_req] = newtype;
  thread_data->send_fanin_infotab[id_req]  = infotab;

  CALL_MPI MPI_Isend(MPI_BOTTOM, 1, thread_data->send_fanin_mpitypes[id_req],
                     FANIN_PROCDST(SOLV_FTGTIND(j)), tag,
                     PASTIX_COMM, &(thread_data->send_fanin_requests[id_req]));
  TEST_MPI("MPI_Isend");
  thread_data->send_fanin_target[id_req] = SOLV_FTGTIND(j);
#    endif /* PASTIX_UPDO_ISEND */
  COMM_CLOCK_STOP;
#  endif /* NO_MPI_TYPE */

}
#endif /* FORCE_NOMPI */


/******************************************/
/*                                        */
/*                UP                      */
/*                                        */
/******************************************/


/*
 * Function: updo_up_WaitCtrb_storage
 *
 * Wait contribution in case where a global rhs is stored.
 *
 * Parameters:
 *   sopalin_data     - Sopalin_data structure to clean
 *   updo_buffer_size - Size of updo_buffer
 *   updo_buffer      - Receipt buffer
 *   me               - Thread Id
 *   i                - Wait contributions fot task i
 *
 * Returns:
 *   void
 */
#ifdef STORAGE
void updo_up_WaitCtrb_storage ( Sopalin_Data_t *sopalin_data,
                                PASTIX_INT      updo_buffer_size,
                                void           *updo_buffer,
                                PASTIX_INT      me,
                                PASTIX_INT      i )
{
  SolverMatrix  *datacode = sopalin_data->datacode;
  PASTIX_FLOAT         *ga, *gb, *gc;
  PASTIX_INT            j, stride;
  PASTIX_INT            size, sizea;
  (void)updo_buffer; (void)updo_buffer_size; (void)me;

  print_debug(DBG_SOPALIN_UP, "%d : Task %4d Wait %4d\n",
              (int)me, (int)i, (int)UPDOWN_CTRBCNT(i));

  /* Boucle sur les blocs du bloc colonne courant à partir du dernier */
  for (j=SYMB_BLOKNUM(i+1)-1; j>=SYMB_BLOKNUM(i)+1; j--)
    {

      print_debug(DBG_SOPALIN_UP, "me=%ld j=%ld gj=%ld\n",
                  (long)me,(long)j,(long)UPDOWN_LBLK2GCBLK(j));

      /* Ignore data comming from schur */
      if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
          UPDOWN_LBLK2GCBLK(j) == UPDOWN_GCBLKNBR-1) {
        UPDOWN_CTRBCNT(i)--;
        continue;
      }

      if (THREAD_COMM_ON)
        {
          /* On attend toutes les contributions locales ou non */
          MUTEX_LOCK(&(sopalin_data->mutex_flagtab[UPDOWN_LBLK2GCBLK(j)]));
          while (!sopalin_data->flagtab[UPDOWN_LBLK2GCBLK(j)])
            COND_WAIT(&(sopalin_data->cond_flagtab[UPDOWN_LBLK2GCBLK(j)]),
                      &(sopalin_data->mutex_flagtab[UPDOWN_LBLK2GCBLK(j)]));
          MUTEX_UNLOCK(&(sopalin_data->mutex_flagtab[UPDOWN_LBLK2GCBLK(j)]));
        }
      else
        {
          /*
           * if the contribution is local
           */
          if (SYMB_CBLKNUM(j)>0)
            {
              MUTEX_LOCK(&(sopalin_data->mutex_flagtab[UPDOWN_LBLK2GCBLK(j)]));
              while (!sopalin_data->flagtab[UPDOWN_LBLK2GCBLK(j)])
                COND_WAIT(&(sopalin_data->cond_flagtab[UPDOWN_LBLK2GCBLK(j)]),
                          &(sopalin_data->mutex_flagtab[UPDOWN_LBLK2GCBLK(j)]));
              MUTEX_UNLOCK(&(sopalin_data->mutex_flagtab[UPDOWN_LBLK2GCBLK(j)]));
            }
          /*
           * if the contribution is not local
           */
#  ifndef FORCE_NOMPI
          else
            {
              MUTEX_LOCK(&(sopalin_data->mutex_flagtab[UPDOWN_LBLK2GCBLK(j)]));

#    if ((defined EXACT_TAG) || (defined EXACT_THREAD)) && !(defined FORCE_NOSMP)
              if(THREAD_COMM_OFF &&
                 !sopalin_data->flagtab[UPDOWN_LBLK2GCBLK(j)])
                {
                  MPI_Status status;
                  Thread_Data_t *thread_data = sopalin_data->thread_data[me];

#      if (defined PASTIX_UPDO_ISEND) && (defined FORCE_CONSO)
                  API_CALL(send_testall_up)(sopalin_data, me);
#      endif /* PASTIX_UPDO_ISEND && FORCE_CONSO */
                  COMM_CLOCK_START;
                  CALL_MPI MPI_Recv(updo_buffer, updo_buffer_size, MPI_BYTE,
                                    MPI_ANY_SOURCE, UPDOWN_LBLK2GCBLK(j),
                                    PASTIX_COMM, &status);
                  TEST_MPI("MPI_Recv");
                  COMM_CLOCK_STOP;
                  API_CALL(updo_up_recv)(sopalin_data, updo_buffer, status, me);
                }
#    else /* EXACT_TAG */
              while(!sopalin_data->flagtab[UPDOWN_LBLK2GCBLK(j)])
                {
                  MUTEX_UNLOCK(&(sopalin_data->mutex_flagtab[UPDOWN_LBLK2GCBLK(j)]));
                  MPI_Status status;
                  Thread_Data_t *thread_data = sopalin_data->thread_data[me];

#      if (defined PASTIX_UPDO_ISEND) && (defined FORCE_CONSO)
                  API_CALL(send_testall_up)(sopalin_data, me);
#      endif /* PASTIX_UPDO_ISEND && FORCE_CONSO */
                  COMM_CLOCK_START;
                  CALL_MPI MPI_Recv(updo_buffer, updo_buffer_size, MPI_BYTE,
                                    MPI_ANY_SOURCE, TAG_UP, PASTIX_COMM,
                                    &status);
                  TEST_MPI("MPI_Recv");
                  COMM_CLOCK_STOP;
                  MUTEX_LOCK(&(sopalin_data->mutex_flagtab[UPDOWN_LBLK2GCBLK(j)]));
                  API_CALL(updo_up_recv)(sopalin_data, updo_buffer, status, me);
                }
#    endif /* EXACT_TAG */
              MUTEX_UNLOCK(&(sopalin_data->mutex_flagtab[UPDOWN_LBLK2GCBLK(j)]));
            }
#  endif /* FORCE_NOMPI */
        }

      ASSERTDBG(SYMB_FROWNUM(j)<UPDOWN_GNODENBR,MOD_SOPALIN);

      gb =&(sopalin_data->grhs[SYMB_FROWNUM(j)]);

      /* OOC : SOLV_COEFTAB loaded and locked above */
#  ifdef SOPALIN_LU
      ga =&(SOLV_UCOEFTAB(i)[SOLV_COEFIND(j)]);
#  else /* not SOPALIN_LU */
      ga =&(SOLV_COEFTAB(i)[SOLV_COEFIND(j)]);
#  endif /* not SOPALIN_LU */

      stride = SOLV_STRIDE(i);
      size   = SYMB_LCOLNUM(i)-SYMB_FCOLNUM(i)+1;
      sizea  = SYMB_LROWNUM(j)-SYMB_FROWNUM(j)+1;
      gc     =&(UPDOWN_SM2XTAB[UPDOWN_SM2XIND(i)]);

      MUTEX_LOCK(&(sopalin_data->mutex_task[i]));
#  ifdef HERMITIAN
#    ifdef MULT_SMX
      SOPALIN_GEMM("C", "N", size, UPDOWN_SM2XNBR, sizea, -fun, ga, stride,
                   gb, UPDOWN_GNODENBR,
                   fun, gc, UPDOWN_SM2XSZE);
#    else /* not MULT_SMX */
      SOPALIN_GEMV("C",sizea,size,-fun,ga,stride,gb,iun,fun,gc,iun);
#    endif /* not MULT_SMX */
#  else /* not HERMITIAN */
#    ifdef MULT_SMX
      SOPALIN_GEMM("T", "N", size, UPDOWN_SM2XNBR, sizea, -fun, ga, stride,
                   gb, UPDOWN_GNODENBR,
                   fun, gc, UPDOWN_SM2XSZE);
#    else /* not MULT_SMX */
      SOPALIN_GEMV("T",sizea,size,-fun,ga,stride,gb,iun,fun,gc,iun);
#    endif /* not MULT_SMX */
#  endif /* not HERMITIAN */
      UPDOWN_CTRBCNT(i)--;

      MUTEX_UNLOCK(&(sopalin_data->mutex_task[i]));
    } /* Fin boucle sur les blocs */

  pthread_cond_broadcast(&(sopalin_data->cond_task[i]));
  ASSERTDBG(UPDOWN_CTRBCNT(i)==0,MOD_SOPALIN);
}

/*
 * Function: updo_up_recv
 *
 * Copy the receipts data in place in the global rhs
 *
 * Parameters:
 *   sopalin_data     - Sopalin_data structure to clean
 *   updo_buffer      - Receipt buffer
 *   status           - Communication status
 *   me               - Thread Id
 *
 * Returns:
 *   void
 */
#  ifndef FORCE_NOMPI
void updo_up_recv ( Sopalin_Data_t *sopalin_data,
                    void           *updo_buffer,
                    MPI_Status      status,
                    PASTIX_INT             me)
{
  PASTIX_INT            size, itersmx;
  PASTIX_INT           *infotab;
  PASTIX_FLOAT         *gb, *lgb, *lgrhs;
  /* #if (defined FLAG_ASSERT) || (defined TRACE_SOPALIN) || (defined PASTIX_DUMP_SOLV_COMM) */
  SolverMatrix  *datacode     = sopalin_data->datacode;
  /* #endif */
#    ifdef TRACE_SOPALIN
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#    endif /* TRACE_SOPALIN */

  infotab = ((PASTIX_INT*)(updo_buffer));

  trace_recv(thread_data->tracefile,
             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me,
             status.MPI_SOURCE, COMM_UP, infotab[2], 0, infotab[5]);

  gb      = ((PASTIX_FLOAT*)(updo_buffer)) +
    UPDOWN_SIZETAB*sizeof(PASTIX_INT)/sizeof(PASTIX_FLOAT);

  /* deposer dans grhs; */
  size = infotab[1] - infotab[0] + 1;

  ASSERTDBG(infotab[0]<UPDOWN_GNODENBR,MOD_SOPALIN);

#    ifdef PASTIX_DUMP_SOLV_COMM
  {
    FILE *file;
    char  filename[30];
    int   i;

    sprintf(filename, "recv/up_%02ld_%02ld_%05ld",
            (long)status.MPI_SOURCE,
            (long)SOLV_PROCNUM,
            (long)infotab[5]);

    file = fopen(filename, "a+");
    for (i=0; i<UPDOWN_SIZETAB; i++)
      fprintf(file, "%ld ", (long)infotab[i]);
    fprintf(file, "\n");
    for (i=0; i<size; i++)
      fprintf(file, "%lf ", gb[i]);
    fprintf(file, "\n");
    fclose(file);
  }
#    endif /* PASTIX_DUMP_SOLV_COMM */

  if (THREAD_COMM_ON)
    MUTEX_LOCK(&(sopalin_data->mutex_flagtab[infotab[2]]));


  lgb   = gb;
  lgrhs = &(sopalin_data->grhs[infotab[0]]);
  for (itersmx=0; itersmx<UPDOWN_SM2XNBR;
       itersmx++, lgb = lgb + size, lgrhs = lgrhs + UPDOWN_GNODENBR)
    SOPALIN_COPY(size, lgb, iun, lgrhs, iun);
  sopalin_data->flagtab[infotab[2]] = 1;

  if (THREAD_COMM_ON)
    MUTEX_UNLOCK(&(sopalin_data->mutex_flagtab[infotab[2]]));

  pthread_cond_broadcast(&(sopalin_data->cond_flagtab[infotab[2]]));
}
#  endif /* FORCE_NOMPI */


#else /* STORAGE */
void updo_up_recv ( Sopalin_Data_t *sopalin_data,
                    void           *updo_buffer,
                    MPI_Status      status,
                    PASTIX_INT             me)
{}
#  ifndef FORCE_NOMPI
/*
 * Function: probe_updown
 *
 * Test if some communications is arrived
 *
 * Parameters:
 *   comm - MPI comunicator
 *   tag  - communication tag
 *
 * Returns:
 *   flag - result of MPI_Iprobe
 */
int probe_updown ( MPI_Comm comm, PASTIX_INT tag )
{
  MPI_Status status;
  int flag;
#    if ((defined EXACT_TAG) || (defined EXACT_THREAD)) && (!defined FORCE_NOSMP)
  PASTIX_INT mytag = tag;
#    else /* not ((EXACT_TAG||EXACT_THREAD) && ! FORCE_NOSMP) */
  PASTIX_INT mytag = TAG_UP;
#    endif /* not ((EXACT_TAG||EXACT_THREAD) && ! FORCE_NOSMP) */

  CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, mytag, comm, &flag, &status);
  TEST_MPI("MPI_Iprobe");
  return flag;
}

/*
 * Function: updo_up_WaitCtrb_nostorage
 *
 * Wait contribution in case where a global rhs is stored.
 *
 * Parameters:
 *   sopalin_data     - Sopalin_data structure to clean
 *   updo_buffer_size - Size of updo_buffer
 *   updo_buffer      - Receipt buffer
 *   me               - Thread Id
 *   i                - Wait contributions fot task i
 *
 * Returns:
 *   void
 */
void updo_up_WaitCtrb_nostorage ( Sopalin_Data_t *sopalin_data,
                                  PASTIX_INT             updo_buffer_size,
                                  void           *updo_buffer,
                                  PASTIX_INT             me,
                                  PASTIX_INT             i )
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  MPI_Status     status;
  PASTIX_FLOAT         *ga, *gb, *gc;
  PASTIX_INT           *infotab;
  PASTIX_INT            k, kk, count, tag;
  PASTIX_INT            size, sizea, stride;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];

  print_debug(DBG_SOPALIN_UP, "attendre %d contrib pour %d\n",
              (int)UPDOWN_CTRBCNT(i), (int)i);

  tag = TAG_UP;
#    if ((defined EXACT_TAG) || (defined EXACT_THREAD)) && (!defined FORCE_NOSMP)
  if (THREAD_COMM_OFF)
    tag = UPDOWN_LOC2GLOB(i);
#    endif /* ((EXACT_TAG||EXACT_THREAD) && ! FORCE_NOSMP) */

  /* receive and add contributions */
  COMM_CLOCK_START;
  CALL_MPI MPI_Recv(updo_buffer,updo_buffer_size,MPI_BYTE,
                    MPI_ANY_SOURCE,tag,PASTIX_COMM,&status);
  TEST_MPI("MPI_Recv");
  COMM_CLOCK_STOP;
  infotab = ((PASTIX_INT*)(updo_buffer));

  trace_begin_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                   STATE_L2_RECVUP, infotab[2]);
  trace_recv      (thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, ,me, status.MPI_SOURCE,
                   COMM_UP, infotab[2], updo_buffer_size, infotab[5]);

  print_debug(DBG_SOPALIN_UP, "Recoit infotab %d %d %d\n",
              (int)infotab[0], (int)infotab[1], (int)infotab[2]);

  gb = ((PASTIX_FLOAT*)(updo_buffer))+
    UPDOWN_SIZETAB*sizeof(PASTIX_INT)/sizeof(PASTIX_FLOAT);

  /* Xi=Xi-LkitXk */
  for (count=UPDOWN_LISTPTR(UPDOWN_GCBLK2LIST(infotab[2]));
       count<UPDOWN_LISTPTR(UPDOWN_GCBLK2LIST(infotab[2])+1);
       count++)
    {
      kk = UPDOWN_LISTCBLK(count);
      k  = UPDOWN_LISTBLOK(count);

      /* Ignore data comming from Schur */
      if (sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
          infotab[2] == UPDOWN_GCBLKNBR-1) {
        /* We discard contribution from shcur */
        UPDOWN_CTRBCNT(kk)--;
        continue;
      }
      ASSERTDBG((SYMB_FROWNUM(k)>=infotab[0]) &&
                (SYMB_LROWNUM(k)<=infotab[1]),MOD_SOPALIN);

      print_debug(DBG_SOPALIN_UP, "trouve blok %d -> MAJ\n", (int)k);

#    ifdef SOPALIN_LU
      ga =&(SOLV_UCOEFTAB(kk)[SOLV_COEFIND(k)]);
#    else /* not SOPALIN_LU */
      ga =&(SOLV_COEFTAB(kk)[SOLV_COEFIND(k)]);
#    endif /* not SOPALIN_LU */
      stride = SOLV_STRIDE(kk);
      size   = SYMB_LCOLNUM(kk) - SYMB_FCOLNUM(kk)+1;
      sizea  = SYMB_LROWNUM(k)  - SYMB_FROWNUM(k) +1;
      gc     =&(UPDOWN_SM2XTAB[UPDOWN_SM2XIND(kk)]);

      MUTEX_LOCK(&(sopalin_data->mutex_task[kk]));
#    ifdef HERMITIAN
#      ifdef MULT_SMX
      SOPALIN_GEMM("C","N",size,UPDOWN_SM2XNBR,sizea,-fun,ga,stride,
                   gb+SYMB_FROWNUM(k)-infotab[0],infotab[1]-infotab[0]+1,
                   fun,gc,UPDOWN_SM2XSZE);
#      else /* not MULT_SMX */
      SOPALIN_GEMV("C",sizea,size,-fun,ga,stride,
                   gb+SYMB_FROWNUM(k)-infotab[0],iun,fun,gc,iun);
#      endif /* not MULT_SMX */
#    else /* not HERMITIAN */
#      ifdef MULT_SMX
      SOPALIN_GEMM("T","N",size,UPDOWN_SM2XNBR,sizea,-fun,ga,stride,
                   gb+SYMB_FROWNUM(k)-infotab[0],infotab[1]-infotab[0]+1,
                   fun,gc,UPDOWN_SM2XSZE);
#      else /* not MULT_SMX */
      SOPALIN_GEMV("T",sizea,size,-fun,ga,stride,
                   gb+SYMB_FROWNUM(k)-infotab[0],iun,fun,gc,iun);
#      endif /* not MULT_SMX */
#    endif /* not HERMITIAN */
      UPDOWN_CTRBCNT(kk)--;
      if (!UPDOWN_CTRBCNT(kk))
        {
          MUTEX_UNLOCK(&(sopalin_data->mutex_task[kk]));
#    ifdef DEP_SMX
          queueAdd(&cblreadyqueue,kk,-(double)(TASK_PRIONUM(kk)));
#    endif /* DEP_SMX */
          pthread_cond_broadcast(&(sopalin_data->cond_task[kk]));
        }
      else
        MUTEX_UNLOCK(&(sopalin_data->mutex_task[kk]));

    }

  trace_end_task(thread_data->tracefile,
                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                 STATE_L2_RECVUP, infotab[2]);
}

#  endif /* FORCE_NOMPI */
#endif /* STORAGE */



/*
 * Function: updo_up_send
 *
 * Send contribution associated to (i,j)
 *
 * Parameters:
 *   sopalin_data     - Sopalin_data structure to clean
 *   me               - Thread Id
 *   i                -
 *   j                -
 *
 * Returns:
 *   void
 */
#ifndef FORCE_NOMPI
void updo_up_send ( Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT i, PASTIX_INT j)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#  if (defined PASTIX_UPDO_ISEND) || (defined NO_MPI_TYPE)
  PASTIX_INT            id_req = 0;
#  endif /* PASTIX_UPDO_ISEND || NO_MPI_TYPE */
  PASTIX_FLOAT         *gb;
  int           *tabsize;
#  ifdef NO_MPI_TYPE
  void         **taboffs;
  int           *tabtype;
  void          *send_buffer;
  PASTIX_INT            send_buffer_size;
  PASTIX_INT            copied;
#  else /* NO_MPI_TYPE */
  MPI_Aint      *taboffs;
  MPI_Datatype  *tabtype;
  MPI_Datatype   newtype;
#  endif /* NO_MPI_TYPE */
  PASTIX_INT            size, iter, tag;
#  if (!defined PASTIX_UPDO_ISEND) || (defined NO_MPI_TYPE)
  PASTIX_INT            infotab[UPDOWN_SIZETAB];
#  else /* !PASTIX_UPDO_ISEND || NO_MPI_TYPE */
  PASTIX_INT           *infotab;
  MALLOC_INTERN(infotab, UPDOWN_SIZETAB, PASTIX_INT);
#  endif /* !PASTIX_UPDO_ISEND || NO_MPI_TYPE */

  gb = &(UPDOWN_SM2XTAB[UPDOWN_SM2XIND(i)]);

  tabsize = thread_data->gtabsize;
  taboffs = thread_data->gtaboffs;
  tabtype = thread_data->gtabtype;

  for (iter = 0; iter < UPDOWN_SIZETAB; iter++)
    infotab[iter] = 0;
  infotab[0] = SYMB_FCOLNUM(i);
  infotab[1] = SYMB_LCOLNUM(i);
  infotab[2] = UPDOWN_LOC2GLOB(i);
  infotab[3] = 0;

  tabsize[0] = UPDOWN_SIZETAB;
  size = SYMB_LCOLNUM(i) - SYMB_FCOLNUM(i) + 1;

#  ifdef   NO_MPI_TYPE
  tabtype[0] = sizeof(PASTIX_INT);
  taboffs[0] = infotab;

  send_buffer_size = tabsize[0]*tabtype[0];

  /* If schur, send empty data */
  if ((sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
       UPDOWN_LOC2GLOB(i) == UPDOWN_GCBLKNBR-1)) {
    for (iter=1; iter<(UPDOWN_SM2XNBR+1); iter++) {
      tabsize[iter] = 0;
      tabtype[iter] = sizeof(PASTIX_FLOAT);
    }
  } else {
    for (iter=1; iter<(UPDOWN_SM2XNBR+1); iter++) {
      tabsize[iter] = size;
      tabtype[iter] = sizeof(PASTIX_FLOAT);
      taboffs[iter] = gb + ((iter-1)*UPDOWN_SM2XSZE);
      send_buffer_size += tabtype[iter]*tabsize[iter];
    }
  }
#  else  /* NO_MPI_TYPE */

  tabtype[0] = COMM_INT;

  CALL_MPI MPI_Address(infotab,&(taboffs[0]));
  TEST_MPI("MPI_Address");

  /* If schur, send empty data */
  if ((sopalin_data->sopar->iparm[IPARM_SCHUR] == API_YES &&
       UPDOWN_LOC2GLOB(i) == UPDOWN_GCBLKNBR-1)) {
    for (iter=1; iter<UPDOWN_SM2XNBR+1; iter++) {
      tabsize[iter] = 0;
      tabtype[iter] = COMM_FLOAT;
    }
  } else {
    for (iter=1; iter<UPDOWN_SM2XNBR+1; iter++) {
      tabsize[iter] = size;
      tabtype[iter] = COMM_FLOAT;
      CALL_MPI MPI_Address(gb+((iter-1)*UPDOWN_SM2XSZE),&(taboffs[iter]));
      TEST_MPI("MPI_Address");
    }
  }
  CALL_MPI MPI_Type_struct(UPDOWN_SM2XNBR+1,tabsize,taboffs,tabtype,&newtype);
  TEST_MPI("MPI_Type_struct");
  CALL_MPI MPI_Type_commit(&newtype);
  TEST_MPI("MPI_Type_commit");

#  endif /* NO_MPI_TYPE */

#  ifdef PASTIX_DUMP_SOLV_COMM
  {
    static int sendup;
    FILE *file;
    char  filename[30];
    int   k;

    infotab[5] = sendup++;
    sprintf(filename, "send/up_%02ld_%02ld_%05ld",
            (long)SOLV_PROCNUM,
            (long) UPDOWN_BROWPROCTAB(i)[j],
            (long)infotab[5]);

    file = fopen(filename, "a+");
    for (k=0; k<UPDOWN_SIZETAB; k++)
      fprintf(file, "%ld ", (long)infotab[k]);
    fprintf(file, "\n");
    for (k=0; k<size; k++)
      fprintf(file, "%lf ", gb[k]);
    fprintf(file, "\n");
    fclose(file);
  }
#  endif /* PASTIX_DUMP_SOLV_COMM */

  /* Choix du tag a utiliser */
#  if (defined STORAGE) && (!defined FORCE_NOSMP)
  if (THREAD_COMM_OFF)
    {
      tag = UPDOWN_LOC2GLOB(i);
    }
  else
#  endif
    {
#  if ((defined EXACT_TAG) || (defined EXACT_THREAD)) && (!defined FORCE_NOSMP)
      if (THREAD_COMM_OFF)
        {
          tag = UPDOWN_BROWCBLKTAB(i)[j];
        }
      else
#  endif /* (EXACT_TAG || EXACT_THREAD) && !FORCE_NOSMP */
        {
          tag = TAG_UP;
        }
    }

  /*********************************/
  /*   Fin de précedents envois    */
  /*********************************/
#  ifdef PASTIX_UPDO_ISEND
  id_req = API_CALL(send_waitone_up)(sopalin_data, me);
  ASSERTDBG(id_req<MAX_S_REQUESTS, MOD_SOPALIN);
#  endif /* PASTIX_UPDO_ISEND */

  /*********************************/
  /*    Send                       */
  /*********************************/

  trace_send(thread_data->tracefile,
             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, UPDOWN_BROWPROCTAB(i)[j],
             COMM_UP, infotab[2], tabsize[0]+tabsize[1], &(infotab[5]));

  print_debug(DBG_SOPALIN_UP, "me=%ld send tag=%ld pour %ld\n",
              (long)me,(long)tag,(long)UPDOWN_BROWPROCTAB(i)[j]);
  print_debug(DBG_SOPALIN_UP, "Envoi X%d vers %d (%d,%d,%d)\n",
              (int)i, (int)UPDOWN_BROWPROCTAB(i)[j],
              (int)infotab[0], (int)infotab[1], (int)infotab[2]);


#  ifdef NO_MPI_TYPE
  MALLOC_INTERN(thread_data->send_fanin_buffer[id_req], send_buffer_size, char);

  send_buffer = thread_data->send_fanin_buffer[id_req];

  copied = 0;
  for (iter = 0; iter <(UPDOWN_SM2XNBR+1); iter ++)
    {
      PASTIX_INT mySize = tabsize[iter]*tabtype[iter];
      memcpy(send_buffer + copied, taboffs[iter], mySize);
      copied += mySize;
    }
  COMM_CLOCK_START;
#    ifndef PASTIX_UPDO_ISEND
  CALL_MPI MPI_Send(send_buffer, send_buffer_size, MPI_CHAR,
                    UPDOWN_BROWPROCTAB(i)[j], tag,
                    PASTIX_COMM);
  TEST_MPI("MPI_Send");
  memFree_null(thread_data->send_fanin_buffer[id_req]);
#    else
  CALL_MPI MPI_Isend(send_buffer, send_buffer_size, MPI_CHAR,
                     UPDOWN_BROWPROCTAB(i)[j], tag,
                     PASTIX_COMM, &(thread_data->send_fanin_requests[id_req]));
  TEST_MPI("MPI_Isend");
#    endif
  COMM_CLOCK_STOP;
#  else /* NO_MPI_TYPE */

  COMM_CLOCK_START;
#    ifndef PASTIX_UPDO_ISEND
  CALL_MPI MPI_Send(MPI_BOTTOM, 1,newtype,
                    UPDOWN_BROWPROCTAB(i)[j], tag,
                    PASTIX_COMM);
  TEST_MPI("MPI_Send");
  CALL_MPI MPI_Type_free(&newtype);
  TEST_MPI("MPI_Type_free");
#    else
  thread_data->send_fanin_mpitypes[id_req] = newtype;
  thread_data->send_fanin_infotab[id_req]  = infotab;
  CALL_MPI MPI_Isend(MPI_BOTTOM, 1, thread_data->send_fanin_mpitypes[id_req],
                     UPDOWN_BROWPROCTAB(i)[j], tag,
                     PASTIX_COMM, &(thread_data->send_fanin_requests[id_req]));
  TEST_MPI("MPI_Isend");
#    endif
  COMM_CLOCK_STOP;

#  endif /* NO_MPI_TYPE */
}
#endif /* FORCE_NOMPI */

/******************************************************************************/
/*                                                                            */
/*                    THREAD DE COMM POUR UP/DOWN                             */
/*                                                                            */
/******************************************************************************/
#ifndef FORCE_NOMPI
void* updo_thread_comm ( void *arg )
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  if (THREAD_COMM_ON)
    {
      SolverMatrix     *datacode     = sopalin_data->datacode;
      MPI_Request      *request;
      PASTIX_INT               me           = argument->me;
      Thread_Data_t    *thread_data;
      MPI_Status        status;
      PASTIX_INT               nbmsg;
      PASTIX_INT               updo_buffer_size;
      void             *updo_buffer  = NULL;
      int               init;
      PASTIX_INT               i, j, nbmsgsend = 0;
      int               flag, wait = 0;

      /* Allocation de la structure de données spécifique a ce thread */
      init = 0;
      if (THREAD_FUNNELED_ON)
        {
          init = init | INIT_SEND;
        }
      sopalin_init_smp(sopalin_data, me, 0, init);
      thread_data = sopalin_data->thread_data[me];

      /* Initialisation buffer communication */
      updo_buffer_size = UPDOWN_SIZETAB*sizeof(PASTIX_INT) +
        UPDOWN_SM2XMAX*sizeof(PASTIX_FLOAT)*UPDOWN_SM2XNBR;

      MALLOC_INTERN(updo_buffer, updo_buffer_size, char);

      MALLOC_INTERN(request, 1, MPI_Request);
      request[0] = MPI_REQUEST_NULL;

      print_debug(DBG_SOPALIN_UPDO, "buffer size = %d\n", (int)updo_buffer_size);

      /********************************/
      /*      Boucle Down             */
      /********************************/
      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      while(sopalin_data->step_comm != COMMSTEP_DOWN)
        COND_WAIT(&(sopalin_data->cond_comm),
                  &(sopalin_data->mutex_comm));
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));

      print_debug(DBG_SOPALIN_UPDO, "--- %d-%d : Debut Boucle Down ---\n", (int)SOLV_PROCNUM, (int)me);
      nbmsg = UPDOWN_DOWNMSGNBR;

      if (UPDOWN_DOWNMSGNBR > 0)
        {
          CALL_MPI MPI_Recv_init(updo_buffer, updo_buffer_size, MPI_BYTE,
                                 MPI_ANY_SOURCE, TAG_DOWN, PASTIX_COMM, request);
          TEST_MPI("MPI_Recv_init");
        }

      if (THREAD_FUNNELED_OFF)
        {
          while(nbmsg){
            COMM_CLOCK_START;
            CALL_MPI MPI_Start(request);
            TEST_MPI("MPI_Start");

            CALL_MPI MPI_Wait(request, &status);
            TEST_MPI("MPI_Wait");
            COMM_CLOCK_STOP;
            API_CALL(updo_down_recv)(sopalin_data, updo_buffer, status, me);
            nbmsg--;
          }
        }
      else
        {
          nbmsgsend = UPDOWN_UPMSGNBR;

          /* Lancement de la premiere reception */
          if (UPDOWN_DOWNMSGNBR > 0)
            {
              CALL_MPI MPI_Start(request);
              TEST_MPI("MPI_Start");
            }

          while(nbmsg || nbmsgsend){

            /* print_debug(DBG_SOPALIN_UPDO, "%d-%d : Recv %d - Send %d\n",  */
            /*             (int)SOLV_PROCNUM, (int)me, nbmsg, nbmsgsend); */

            /*
             * Attente si rien à faire
             */
            if (wait)
              {
                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                COND_TIMEWAIT(&(sopalin_data->cond_comm), &(sopalin_data->mutex_comm));
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
              }
            wait = 1;

            /* Reception */
            if (nbmsg)
              {
                CALL_MPI MPI_Test(request, &flag, &status);
                TEST_MPI("MPI_Wait");

                if (flag)
                  {
                    API_CALL(updo_down_recv)(sopalin_data, updo_buffer, status, me);
                    nbmsg--;
                    if (nbmsg)
                      {
                        CALL_MPI MPI_Start(request);
                        TEST_MPI("MPI_Start");
                      }
                    wait = 0;
                  }
              }

            /* Envoi */
            if (nbmsgsend)
              {
                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                if (queueSize(sopalin_data->sendqueue))
                  {
                    i = queueGet2(sopalin_data->sendqueue, NULL, &j);
                    API_CALL(updo_down_send)(sopalin_data, me, i, j);
                    MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
                    nbmsgsend--;
                    wait = 0;
                  }
                else
                  MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
              }
          }
        }

      if (THREAD_FUNNELED_ON)
        {
          /* Attente de la fin des communications en envoi */
          if (SOLV_PROCNBR > 1)
            {
              API_CALL(send_waitall_down)(sopalin_data, me);
            }
        }

      if (UPDOWN_DOWNMSGNBR > 0)
        {
          CALL_MPI MPI_Request_free(request);
          TEST_MPI("MPI_Request_free");
        }

      if (me == SOLV_THRDNBR)
        {
          MUTEX_LOCK(&(sopalin_data->mutex_comm));
          sopalin_data->step_comm = COMMSTEP_UPDOEND;
          print_debug(DBG_THCOMM, "%s:%d UPDOEND\n", __FILE__, __LINE__);
          MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
          pthread_cond_broadcast(&(sopalin_data->cond_comm));
        }

      /********************************/
      /*      Boucle Up               */
      /********************************/

      MUTEX_LOCK(&(sopalin_data->mutex_comm));
      while(sopalin_data->step_comm != COMMSTEP_UP)
        COND_WAIT(&(sopalin_data->cond_comm),
                  &(sopalin_data->mutex_comm));
      MUTEX_UNLOCK(&(sopalin_data->mutex_comm));

      print_debug(DBG_SOPALIN_UPDO, "--- %d-%d : Debut Boucle Up ---\n", (int)SOLV_PROCNUM, (int)me);
      nbmsg = UPDOWN_UPMSGNBR;

      if (UPDOWN_UPMSGNBR > 0)
        {
          CALL_MPI MPI_Recv_init(updo_buffer, updo_buffer_size, MPI_BYTE,
                                 MPI_ANY_SOURCE, TAG_UP, PASTIX_COMM, request);
          TEST_MPI("MPI_Recv_init");
        }
      if (THREAD_FUNNELED_OFF)
        {
          while(nbmsg){
            print_debug(DBG_SOPALIN_UPDO, "--- NB UP to receive %d ---\n", (int)nbmsg);
            COMM_CLOCK_START;
            CALL_MPI MPI_Start(request);
            TEST_MPI("MPI_Start");

            CALL_MPI MPI_Wait(request, &status);
            TEST_MPI("MPI_Wait");
            COMM_CLOCK_STOP;
            API_CALL(updo_up_recv)(sopalin_data, updo_buffer, status, me);
            nbmsg--;
          }
        }
      else
        {
          nbmsgsend = UPDOWN_DOWNMSGNBR;
          wait = 0;

          /* Lancement de la premiere reception */
          if (UPDOWN_UPMSGNBR > 0)
            {
              CALL_MPI MPI_Start(request);
              TEST_MPI("MPI_Start");
            }
          while(nbmsg || nbmsgsend){

            /* print_debug(DBG_SOPALIN_UPDO, "%d-%d : Recv %d - Send %d\n", */
            /*                  (int)SOLV_PROCNUM, (int)me, nbmsg, nbmsgsend);      */

            /*
             * Attente si rien à faire
             */
            if (wait)
              {
                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                COND_TIMEWAIT(&(sopalin_data->cond_comm), &(sopalin_data->mutex_comm));
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
              }
            wait = 1;

            /* Reception */
            if (nbmsg)
              {
                CALL_MPI MPI_Test(request, &flag, &status);
                TEST_MPI("MPI_Wait");

                if (flag)
                  {
                    API_CALL(updo_up_recv)(sopalin_data, updo_buffer, status, me);
                    nbmsg--;

                    if (nbmsg)
                      {
                        CALL_MPI MPI_Start(request);
                        TEST_MPI("MPI_Start");
                      }
                    wait = 0;
                  }
              }

            /* Envoi */
            if (nbmsgsend)
              {
                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                if (queueSize(sopalin_data->sendqueue))
                  {
                    i = queueGet2(sopalin_data->sendqueue, NULL, &j);
                    API_CALL(updo_up_send)(sopalin_data, me, i, j);
                    nbmsgsend--;
                    wait = 0;
                  }
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
              }
          }
        }

      if (UPDOWN_UPMSGNBR > 0)
        {
          CALL_MPI MPI_Request_free(request);
          TEST_MPI("MPI_Request_free");
        }
      print_debug(DBG_SOPALIN_UPDO, "--- Fin Comm ---\n");

      if (THREAD_FUNNELED_ON)
        {
          /* Attente de la fin des communications en envoi */
          if (SOLV_PROCNBR > 1)
            {
              API_CALL(send_waitall_up)(sopalin_data, me);
            }
        }

      /* Liberation de la memoire */
      memFree_null(updo_buffer);
      memFree_null(request);
      sopalin_clean_smp(sopalin_data, me);

      if (me == SOLV_THRDNBR)
        {
          MUTEX_LOCK(&(sopalin_data->mutex_comm));
          sopalin_data->step_comm = COMMSTEP_UPDOEND;
          print_debug(DBG_THCOMM, "%s:%d UPDOEND\n", __FILE__, __LINE__);
          MUTEX_UNLOCK(&(sopalin_data->mutex_comm));



          pthread_cond_broadcast(&(sopalin_data->cond_comm));
        }
      return (void *)0;
    }
  else
    {
      errorPrint("updo_thread_comm should not be called in THREAD_COMM mode.");
      return (void *)1;
    }
}
#endif /* FORCE_NOMPI */
