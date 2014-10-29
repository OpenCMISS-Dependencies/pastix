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
  File: trace.h                                 

  trace instruction for different format

  Authors: 
    Mathieu Faverge - faverge@labri.fr

  Dates: 
    Varsion 0.0 - from 12 apr 2008
                  to   12 apr 2008
*/
#ifndef TRACE_H
#define TRACE_H

#ifdef TRACE_SOPALIN

typedef enum API_TRACEFMT TraceFmt_t;

typedef enum Trace_State {
  STATE_IDLE         =  0,
  STATE_WAITLOC      =  1,
  STATE_WAITREM      =  2,
  STATE_WAITTSK	     =  3,
  STATE_L0_IDLE      =  0,
  STATE_L0_FACTOINIT =  4,
  STATE_L0_FACTOCLEAN=  5,
  STATE_L0_FACTO     =  6,
  STATE_L0_UPDOINIT  =  7,
  STATE_L0_UPDOCLEAN =  8,
  STATE_L0_DOWN	     =  9,
  STATE_L0_DIAG	     = 10,
  STATE_L0_UP	     = 11,
  STATE_L0_RAFF	     = 12,
  STATE_COMPUTE	     = 13,
  STATE_BLOCK	     = 14,
  STATE_DIAG	     = 15,
  STATE_COMP1D	     = 16,
  STATE_COMP1DPLUS   = 17,
  STATE_COMP1DGEMM   = 18,
  STATE_E1	     = 19,
  STATE_E2	     = 20,
  STATE_DOWN	     = 21,
  STATE_UP	     = 22,
  STATE_L2_ADD	     = 23,
  STATE_L2_SENDF     = 24,
  STATE_L2_SENDB     = 25,
  STATE_L2_RECVF     = 26,
  STATE_L2_RECVB     = 27,
  STATE_L2_RECVDOWN  = 28,
  STATE_L2_RECVUP    = 29,
  STATE_COMMG	     = 30,
  STATE_ALLOC	     = 31,
  STATE_FREE         = 32,
  STATE_L2_COMP1D    = 33,
  STATE_L2_COMP1DGEMM= 34,
  STATE_L2_E1	     = 35,
  STATE_L2_E2	     = 36,
  STATE_L2_DIAG	     = 37,
  STATE_L2_DOWN	     = 38,
  STATE_L2_UP	     = 39,
  STATE_NBSTATES     = 40,
} Trace_State_t;       

typedef enum Trace_Comm {
  COMM_FANIN = 0,
  COMM_BLOCK,
  COMM_DOWN,
  COMM_UP,
  COMM_NBTYPECOMM
} Trace_Comm_t;

/*
 * Common part :
 *    - TraceFmt_t  fmt     : Format de trace
 *    - FILE       *file    : fichier de trace
 *    - double      time    : timestamp
 *    - PASTIX_INT         procnum
 *    - PASTIX_INT         thrdnum
 * State :
 *    - PASTIX_INT           level : niveau de trace (+il est grand plus c'est détaillé)
 *    - Trace_State_T state : Etat dans lequel se retrouve le proc
 *    - PASTIX_INT           id    : identifiant de la tâche
 * Comm :
 *    - PASTIX_INT          dest/src : destinataire ou source du message
 *    - Trace_Comm_t Type     : type de communication
 *    - PASTIX_INT          id       : identifiant du message
 *    - PASTIX_INT          size     : taille du message
 */   

extern volatile PASTIX_INT idg;

#define trace_start(file, time, procnum, thrdnum) \
  fprintf(file, "%9.9lf %ld %ld 0 -1 0\n", (double)time, (long)(procnum), (long)(thrdnum));
#define trace_finish(file, time, procnum, thrdnum) \
  fprintf(file, "%9.9lf %ld %ld 0 -1 1\n", (double)time, (long)(procnum), (long)(thrdnum));
#define trace_send(file, time, procnum, thrdnum, dest, type, id, size, idreq) \
  {									\
  PASTIX_INT idl = idg++;							\
  fprintf(file, "%9.9lf %ld %ld %ld 2 %ld %ld %ld\n",			\
	  (double)time, (long)procnum, (long)thrdnum, (long)dest, (long)type, (long)id, (long)idl);	\
  *idreq = idl;								\
  }
#define trace_recv(file, time, procnum, thrdnum, src, type, id, size, idreq) \
  fprintf(file, "%9.9lf %ld %ld %ld 3 %ld %ld %ld\n",			\
	  (double)time, (long)procnum, (long)thrdnum, (long)src, (long)type, (long)id, (long)idreq);

#define trace_begin_task(file, time, procnum, thrdnum, level, state, id) \
  {									\
  if (level != 2)							\
    fprintf(file, "%9.9lf %ld %ld %ld 0 %ld %ld\n",			\
	    (double)time, (long)(procnum), (long)(thrdnum),		\
	    (long)(level), (long)(state), (long)(id));			\
  }

#define trace_end_task(file, time, procnum, thrdnum, level, state, id)\
  {									\
    if (level != 2)							\
      fprintf(file, "%9.9lf %ld %ld %ld 1 %ld %ld\n",			\
	      (double)time, (long)(procnum), (long)(thrdnum),		\
	      (long)(level), (long)(state), (long)(id));		\
  }

#define trace_begin_task1(file, time, procnum, thrdnum, state, id2, task, stolen) \
  {									\
    fprintf(file, "%9.9lf %ld %ld 1 0 %ld %ld %ld %ld %ld %ld %ld %ld\n", \
	    (double)time, (long)(procnum), (long)(thrdnum), (long)(state), \
	    (long)(task.id), (long)(task.threadid), (long)(task.fcandnum), \
	    (long)(task.lcandnum), (long)(task.cand), (long)(task.taskid), \
	    (long)(id2));						\
  }

#define trace_begin_task2(file, time, procnum, thrdnum, state, id, fcand, lcand, cand) \
  {									\
    fprintf(file, "%9.9lf %ld %ld 1 0 %ld %ld %ld %ld %ld %ld\n",	\
	    (double)time,  (long)(procnum), (long)(thrdnum), (long)(state), \
	    (long)(id),    (long)0, (long)(fcand),			\
	    (long)(lcand), (long)(cand)					\
	    );								\
  }

#define trace_malloc(file, time, procnum, state, memory)	\
  fprintf(file, "%9.9lf %ld 0 0 4 %ld %f\n",				\
	  (double)time, (long)procnum, (long)state, (double)memory/(1<<20));

#else /* TRACE_SOPALIN */

#define trace_start(...)      do {} while(0)
#define trace_finish(...)     do {} while(0)
#define trace_begin_task(...) do {} while(0)
#define trace_end_task(...)   do {} while(0)
#define trace_begin_task1(...)do {} while(0)
#define trace_begin_task2(...)do {} while(0)
#define trace_send(...)       do {} while(0)
#define trace_recv(...)       do {} while(0)
#define trace_malloc(...)     do {} while(0)

#endif /* TRACE_SOPALIN */

#define trace_settask(...) do {} while(0)
#define trace_addtask(...) do {} while(0)
#define trace_deltask(...) do {} while(0)
#define trace_end_task1()  do {} while(0)

#endif /* TRACE_H */
