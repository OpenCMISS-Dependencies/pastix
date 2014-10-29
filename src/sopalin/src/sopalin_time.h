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
   File: sopalin_time.h
   
   Few macros to manage clocks in sopalin.

   macros: 
     SOPALIN_CLOCK_INIT  - Initiate clock and start it.
     SOPALIN_CLOCK_STOP  - Stop clock.
     SOPALIN_CLOCK_GET   - Get value from clock.
     SOPALIN_CLOCK_TRACE - Get clock relatively to tracing start.
 */
#define SOPALIN_CLOCK_INIT  {clockInit(&(thread_data->sop_clk));clockStart(&(thread_data->sop_clk));}
#define SOPALIN_CLOCK_STOP  {clockStop(&(thread_data->sop_clk));}
#define SOPALIN_CLOCK_GET   clockVal(&(thread_data->sop_clk))
#define COMM_CLOCK_INIT  clockInit(&(thread_data->sop_clk_comm))
#define COMM_CLOCK_START clockStart(&(thread_data->sop_clk_comm))
#define COMM_CLOCK_STOP  clockStop(&(thread_data->sop_clk_comm))
#define COMM_CLOCK_GET   clockVal(&(thread_data->sop_clk_comm))
#define SOPALIN_CLOCK_TRACE (clockGet() - (sopalin_data->timestamp))
/* #define SOPALIN_CLOCK_TRACE ((clockGet() - (thread_data->sop_clk).time[0])) */
/* #define SOPALIN_CLOCK_TRACE (clockGet()) */
/* #define SOPALIN_CLOCK_TRACE (MPI_Wtime()) */

#ifdef COMPUTE_ALLOC 

static Clock alloc_clk;

#define ALLOC_CLOCK_INIT {clockInit(&alloc_clk);clockStart(&alloc_clk);}
#define ALLOC_CLOCK_STOP {clockStop(&alloc_clk);}
#define ALLOC_CLOCK_GET  clockVal(&alloc_clk)

#endif /* COMPUTE_ALLOC */
