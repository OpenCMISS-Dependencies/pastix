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
#ifdef TRACE_SOPALIN

#include <stdio.h>
#include "common_pastix.h"
#include "trace.h"

char *textstate[STATE_NBSTATES] =
  {
    "Inactif",
    "WaitLocal",
    "WaitDist",
    "WaitTask",
    "FactoInit",
    "FactoClean",
    "Facto",
    "UpdoInit",
    "UpdoClean",
    "Down",
    "Up",
    "Raff",
    "Compute",
    "Block",
    "Diag",
    "Comp1d",
    "E1",
    "E2",
    "Add",
    "SendF",
    "SendB",
    "RecvF",
    "RecvB",
    "RecvDown",
    "RecvUp",
    "Commg"
  };

char *textcomm[COMM_NBTYPECOMM] =
  {
    "Fanin",
    "Block",
    "Down",
    "Up"
  };

volatile PASTIX_INT idg = 0;

/* /\* Write time when processor begin computation *\/ */
/* void trace_start(TraceFmt_t fmt, FILE *file, double time,  */
/* 		 PASTIX_INT procnum, PASTIX_INT thrdnum) */
/* { */
/*   switch(fmt) */
/*     { */
/*     case API_TRACE_PICL: */
/*       fprintf(file, "-3 -901 %10.10g %ld -1 0 \n",         (double)time, (long)procnum); */
/*       fprintf(file, "-2 -904 %10.10g %ld -1 3 2 4 4 0 \n", (double)time, (long)procnum); */
/*       fprintf(file, "-3 -402 %10.10g %ld -1 0 \n",         (double)time, (long)procnum); */
/*       fprintf(file, "-4 -402 %10.10g %ld -1 0 \n",         (double)time, (long)procnum); */
/*       break; */
/*     case API_TRACE_PAJE: */
/*     case API_TRACE_HUMREAD: */
/*     default: */
/*       fprintf(file, "%9.9lf %d %d 0 -1 0\n", */
/* 	      (double)time, (int)procnum, (int)thrdnum); */
/*       break; */
/*     } */
/* } */

/* /\* Write time when processor exit computation *\/ */
/* void trace_finish(TraceFmt_t fmt, FILE *file, double time,  */
/* 		  PASTIX_INT procnum, PASTIX_INT thrdnum) */
/* { */
/*   switch(fmt) */
/*     { */
/*     case API_TRACE_PICL: */
/*       fprintf(file, "-4 -901 %10.10g %ld -1 0 \n", (double)time, (long)procnum); */
/*       break; */
/*     case API_TRACE_PAJE: */
/*     case API_TRACE_HUMREAD: */
/*     default: */
/*       fprintf(file, "%9.9lf %d %d 0 -1 1\n", */
/* 	      (double)time, (int)procnum, (int)thrdnum); */
/*       break; */
/*     } */
/* } */


/* /\* Write time when the message (No "id") of size "size" in octet is sent by the processor "procnum"  */
/*    to the processor "dest" *\/ */
/* void trace_send(TraceFmt_t fmt, FILE *file, double time,  */
/* 		PASTIX_INT procnum, PASTIX_INT thrdnum,  */
/* 		PASTIX_INT dest, Trace_Comm_t type, PASTIX_INT id, PASTIX_INT size, PASTIX_INT *idreq) */
/* { */
/*   static volatile PASTIX_INT idg = 0; */
/*   PASTIX_INT idl = idg++; */

/*   switch(fmt) */
/*       { */
/*       case API_TRACE_PICL: */
/* 	fprintf(file, "-3 -27 %10.10g %ld -1 3 2 %ld 88 %ld \n", (double)time, (long)procnum, (long)size, (long)dest); */
/* 	fprintf(file, "-4 -27 %10.10g %ld -1 1 2 %ld \n",        (double)time, (long)procnum, (long)id); */
/* 	fprintf(file, "-3 -30 %10.10g %ld -1 1 2 %ld \n",        (double)time, (long)procnum, (long)id); */
/*         fprintf(file, "-4 -30 %10.10g %ld -1 0 \n",              (double)time, (long)procnum); */
/* 	break; */
/*       case API_TRACE_PAJE: */
/* 	fprintf(file,"42 %f L_1 C_Net0 C_P%d %ld %d_%d_%ld\n", */
/* 		(double)time, (int)procnum, (long)id, (int)procnum, (int)dest, (long)id); */
/* 	break; */
/*       case API_TRACE_HUMREAD: */
/* 	fprintf(file, "%9.9lf ( %02d - %02d ) SEND%s %02d %ld %ld\n", */
/* 		(double)time, (int)procnum, (int)thrdnum, textcomm[type], (int)dest, (long)id, (long)size); */
/* 	break; */
/*       default: */
/* 	fprintf(file, "%9.9lf %d %d %d 2 %d %ld %ld\n", */
/* 		(double)time, (int)procnum, (int)thrdnum, (int)dest, (int)type, (long)id, (long)idl); */
/*       } */
/*   *idreq = idl; */
/* } */

/* /\* Write time when the message (No "id") of size "size" in octet is received by the processor "procnum"  */
/*   from the processor "src" *\/ */
/* void trace_recv(TraceFmt_t fmt, FILE *file, double time,  */
/* 		PASTIX_INT procnum, PASTIX_INT thrdnum,  */
/* 		PASTIX_INT src, Trace_Comm_t type, PASTIX_INT id, PASTIX_INT size, PASTIX_INT idreq) */
/* { */
/*   switch(fmt) */
/*     { */
/*     case API_TRACE_PICL: */
/*       fprintf(file, "-3 -57 %10.10g %ld -1 3 2 %ld 88 %ld \n", (double)time, (long)procnum, (long)size, (long)procnum); */
/*       fprintf(file, "-4 -57 %10.10g %ld -1 1 2 %ld \n",        (double)time, (long)procnum, (long)id); */
/*       fprintf(file, "-3 -60 %10.10g %ld -1 1 2 %ld \n",        (double)time, (long)procnum, (long)id); */
/*       fprintf(file, "-4 -60 %10.10g %ld -1 3 2 %ld 88 %ld \n", (double)time, (long)procnum, (long)size, (long)src); */
/*       break; */
/*     case API_TRACE_PAJE: */
/*       fprintf(file,"43 %f L_1 C_Net0 C_P%d %ld %d_%d_%ld\n",  */
/* 	      (double)time, (int)procnum, (long)id, (int)src, (int)procnum, (long)id); */
/*       break; */
/*     case API_TRACE_HUMREAD: */
/*       fprintf(file, "%9.9lf ( %02d - %02d ) RECV%s %02d %ld %ld\n", */
/* 	      (double)time, (int)procnum, (int)thrdnum, textcomm[type], (int)src, (long)id, (long)size); */
/*       break; */
/*     default: */
/*       fprintf(file, "%9.9lf %d %d %d 3 %d %ld %ld\n", */
/* 	      (double)time, (int)procnum, (int)thrdnum, (int)src, (int)type, (long)id, (long)idreq); */
/*     } */
/* } */

/* /\* Write time when the computational task of type "no" on processor "procnum" is beginning *\/ */
/* void trace_begin_task(TraceFmt_t fmt, FILE *file, double time,  */
/* 		      PASTIX_INT procnum, PASTIX_INT thrdnum, PASTIX_INT level,  */
/* 		      Trace_State_t state, PASTIX_INT id) */
/* { */
/*   switch(fmt) */
/*     { */
/*     case API_TRACE_PICL: */
/*       if (level != 2) return; */
/*       fprintf(file, "-3 %ld %10.10g %ld -1 0 \n", (long)state, (double)time, (long)procnum); */
/*       break; */
/*     case API_TRACE_PAJE: */
/*       fprintf(file, "10 %f ST_ProcState C_P%d S_%ld\n", */
/* 	      (double)time, (int)procnum, (long)state); */
/*       break; */
/*     case API_TRACE_HUMREAD: */
/*       if (level != 1) return; */
/*       fprintf(file, "%9.9lf ( %02d - %02d ) %s %ld\n", */
/* 	      (double)time, (int)procnum, (int)thrdnum, textstate[state], (long)id); */
/*       break; */
/*     default: */
/*       if (level == 2) return; */
/*       fprintf(file, "%9.9lf %d %d %d 0 %d %ld\n", */
/* 	      (double)time, (int)procnum, (int)thrdnum, (int)level, (int)state, (long)id); */
/*     } */
/* } */

/* /\* Write time when the computational task of type "no" on processor "procnum" ends *\/ */
/* void trace_end_task(TraceFmt_t fmt, FILE *file, double time,  */
/* 		    PASTIX_INT procnum, PASTIX_INT thrdnum, PASTIX_INT level,  */
/* 		    Trace_State_t state, PASTIX_INT id) */
/* { */
/*   switch(fmt) */
/*     { */
/*     case API_TRACE_PICL: */
/*       if (level != 2) return; */
/*       fprintf(file, "-4 %ld %10.10g %ld -1 0 \n", (long)state, (double)time, (long)procnum); */
/*       break; */
/*     case API_TRACE_PAJE: */
/*       break; */
/*     case API_TRACE_HUMREAD: */
/*       if (level != 1) return; */
/*       fprintf(file, "%9.9lf ( %02d - %02d ) End%s %ld\n", */
/* 	      (double)time, (int)procnum, (int)thrdnum, textstate[state], (long)id); */
/*       break; */
/*     default: */
/*       if (level == 2) return; */
/*       fprintf(file, "%9.9lf %d %d %d 1 %d %ld\n", */
/* 	      (double)time, (int)procnum, (int)thrdnum, (int)level, (int)state, (long)id); */
/*     } */
/* } */


/* /\* Write time and memory after a malloc *\/ */
/* void trace_malloc(TraceFmt_t fmt, FILE *file, double time,  */
/* 		  PASTIX_INT procnum, */
/* 		  Trace_State_t state, PASTIX_INT memory) */
/* { */
/*   switch(fmt) */
/*     { */
/*     case API_TRACE_PICL: */
/*       fprintf(file, "-3 %ld %10.10g %ld -1 0 %ld\n", (long)state, (double)time, (long)procnum, (long)memory); */
/*       break; */
/*     case API_TRACE_PAJE: */
/*       fprintf(file, "51 %f V_Mem C_N%d %f\n", */
/* 	      (double)time, (int)procnum, (double)memory/(1<<20)); */
/*       break; */
/*     case API_TRACE_HUMREAD: */
/*       fprintf(file, "%9.9lf ( %02d ) %s %ld\n", */
/* 	      (double)time, (int)procnum,  textstate[state], (long)memory); */
/*       break; */
/*     default: */
/*       /\* On force le thrdnum a 0 et le level a 0 pour avoir la meme  */
/* 	 structure de base sur chaque ligne de la trace *\/ */
/*       fprintf(file, "%9.9lf %d 0 0 4 %d %f\n", */
/* 	      (double)time, (int)procnum, (int)state, (double)memory/(1<<20)); */
/*     } */
/* } */
#else
/* ISO C forbids an empty source file */
#include "not_empty.h"
NOT_EMPTY(trace)
#endif /* TRACE_SOPALIN */
