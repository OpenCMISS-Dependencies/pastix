/* Copyright INRIA 2004
**
** This file is part of the Scotch distribution.
**
** The Scotch distribution is libre/free software; you can
** redistribute it and/or modify it under the terms of the
** GNU Lesser General Public License as published by the
** Free Software Foundation; either version 2.1 of the
** License, or (at your option) any later version.
**
** The Scotch distribution is distributed in the hope that
** it will be useful, but WITHOUT ANY WARRANTY; without even
** the implied warranty of MERCHANTABILITY or FITNESS FOR A
** PARTICULAR PURPOSE. See the GNU Lesser General Public
** License for more details.
**
** You should have received a copy of the GNU Lesser General
** Public License along with the Scotch distribution; if not,
** write to the Free Software Foundation, Inc.,
** 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
**
** $Id: common_memory.c 176 2004-10-12 13:53:26Z goureman $
*/
/*
 * File: common_memory.c
 *
 * Part of a parallel direct block solver.
 *
 * This module handles errors.
 *
 * Authors:
 *   Mathieu Faverge    - faverge@labri.fr
 *   Xavier   LACOSTE    - lacoste@labri.fr
 *   Francois PELLEGRINI - .
 *
 * Dates:
 *   Version 0.0  - from 07 sep 2001
 *                  to   07 sep 2001
 *   Version 0.1  - from 14 apr 2001
 *                  to   24 mar 2003
 *   Version 1.3  - from 25 feb 2004
 *                  to   25 feb 2004
 */

/*
 **  The defines and includes.
 */

#define COMMON_MEMORY

#include "common_pastix.h"
#ifdef PASTIX_EZTRACE
#  include "pastix_eztrace.h"
#else
#  include "trace.h"
#endif
#ifdef MEMORY_USAGE
#  include <pthread.h>
#endif /* MEMORY_USAGE */

#ifndef TRACE_SOPALIN
#  undef  MEMORY_TRACE
#endif

/*
 * Group: Variables
 *
 * The static variables.
 *
 * int: memallocmutexflag
 *   Boolean indicating if <memallocmutexdat> mutex has been initialized.
 *
 * pthread_mutex_t: memallocmutexdat
 *   mutex protecting <memalloccurrent>, <memallocmax>, <memalloctraceflag>,
 *   <trace_file>, <trace_timestamp> and <trace_procnum>
 *
 * ulong: memalloccurrent
 *   Current memory allocated using <memAlloc_func>.
 *
 * ulong: memallocmax
 *   Maximum value of memalloccurrent since the program started.
 *
 * int: memalloctraceflag
 *   Boolean indicating if we want to trace allocation.
 *
 * stream: trace_file
 *   File into which to write traces.
 *
 * double: time_stamp
 *   Origin of traces.
 *
 * int: trace_procnum
 *   Processor tracing allocations.
 */

#ifdef MEMORY_USAGE
static int                  memallocmutexflag = 0;
static pthread_mutex_t      memallocmutexdat;     /*+ Local mutex +*/

unsigned long               memalloccurrent   = 0;
unsigned long               memallocmax       = 0;
int                         memalloctraceflag = 0;

#  ifdef MEMORY_TRACE
static FILE                *trace_file;
static double               trace_timestamp;
static int                  trace_procnum;
#  endif
#endif /* MEMORY_USAGE */

/*
 * Group: Functions
 *
 * The memory handling routines.
 */

void * mymalloc ( size_t size,
                  char * filename,
                  int    line)
{
  if (size > 0) {
    return malloc(size);
  }
  else {
    fprintf(stderr, "Pb Alloc 0 %s:%d\n", filename, line);
    return (void *)NULL;
  }
}

#ifdef MEMORY_USAGE

#  ifdef MEMORY_TRACE
/*
 * Function: memAllocaTrace
 *
 * Start tracing memory.
 *
 * Initialize <memallocmutexdat> if not done.
 *
 * Defines all tracing variables.
 *
 * Parameters:
 *   file      - Stream where to write traces, opened in write mode.
 *   timestamp - Traces origin.
 *   procnum   - Processor writting traces.
 *
 * Returns:
 *   void - In all cases.
 */
void memAllocTrace (FILE  *file,
                    double timestamp,
                    int    procnum)
{
  if (memallocmutexflag == 0) {
    /* If memory mutex not yet initialized */
    memallocmutexflag = 1;
    /* Initialize local mutex */
    pthread_mutex_init (&memallocmutexdat, NULL);
  }

  /* Lock local mutex */
  pthread_mutex_lock (&memallocmutexdat);
  memalloctraceflag = 1;
  trace_file        = file;
  trace_timestamp   = timestamp;
  trace_procnum     = procnum;
  /* UnLock local mutex */
  pthread_mutex_unlock (&memallocmutexdat);

}
/*
 * Function: memAllocUntrace
 *
 * Stop tracing allocations.
 *
 * Returns:
 *   void - in all cases.
 */
void memAllocUntrace ()
{
  /* Lock local mutex */
  pthread_mutex_lock (&memallocmutexdat);
  memalloctraceflag = 0;
  /* UnLock local mutex */
  pthread_mutex_unlock (&memallocmutexdat);
}
#  endif
/*
 * Function: memAllocGetCurrent
 *
 * Get the current memory allocated.
 *
 * Returns:
 *   <memalloccurrent> value.
 */
unsigned long memAllocGetCurrent ()
{
  return (memalloccurrent);
}

/*
 * Function: memAllocGetMax
 *
 * Get the maximu memory allocated.
 *
 * Returns:
 *   <memallocmax> value.
 */
unsigned long memAllocGetMax () {
  return (memallocmax);
}

/*
 * Function: memAllocTraceReset
 *
 * Restarts tracing allocation with reseting <memallocmax>.
 *
 * Returns:
 *   void - in all cases.
 */
void memAllocTraceReset () {
  if (memallocmutexflag == 0) {
    /* If memory mutex not yet initialized */
    memallocmutexflag = 1;
    pthread_mutex_init (&memallocmutexdat, NULL); /* Initialize local mutex */
  }
  pthread_mutex_lock (&memallocmutexdat); /* Lock local mutex */
  memalloctraceflag = 1;
  memallocmax = 0;
  pthread_mutex_unlock (&memallocmutexdat); /* UnLock local mutex */
}

/*
 *  Function: memAlloc_func
 *
 *  This is a thread-safe memory allocation routine.
 *
 *  Parameters:
 *    size     - Memory size wanted.
 *    filename - Used for error message, file where the function is called.
 *    line     - Used for erro message, line where the function is called.
 *
 *  Returns:
 *    !NULL - pointer to memory block.
 *    NULL  - no array allocated.
 */

void * memAlloc_func ( size_t size,
                       char * filename,
                       int    line)
{
  double *              memptr;
#  ifdef DEBUG_ALLOC
  if (size == 0)
    errorPrint("%s:%d allocating 0\n",filename,line);
#  endif

  if (memallocmutexflag == 0) {
    /* If memory mutex not yet initialized */
    memallocmutexflag = 1;
    pthread_mutex_init (&memallocmutexdat, NULL); /* Initialize local mutex */
  }

  pthread_mutex_lock (&memallocmutexdat); /* Lock local mutex */

  /* Add a double containing the size of the array */
  memptr = (double*)malloc (size + sizeof(double));
  if (memptr != NULL) {
    memptr[0] = (double)size;
    memalloccurrent += (unsigned long) size;
    memallocmax = MAX (memallocmax, memalloccurrent);
    memptr ++;
#  ifdef MEMORY_TRACE
    if (memalloctraceflag)
      trace_malloc(trace_file, (clockGet()-trace_timestamp),
                   trace_procnum, STATE_ALLOC, memalloccurrent);
#  endif
  }
  else {
    perror("malloc");
    fprintf (stderr, "File : %s:%d\n"
	     "ALLOCATION FAILURE, used : %lu (%.3g %s), asked : %lu (%.3g %s)\n",
             filename, line,
	     memalloccurrent, MEMORY_WRITE(memalloccurrent), MEMORY_UNIT_WRITE(memalloccurrent),
	     (unsigned long)size, MEMORY_WRITE(size), MEMORY_UNIT_WRITE(size));
  }

  /* Unlock local mutex */
  pthread_mutex_unlock (&memallocmutexdat);

  return ((void*)memptr);
}

/*
 *  Function: memRealloc_func
 *
 *  This is a thread-safe memory reallocation routine.
 *
 *  Parameters:
 *    memptr - address of the array to realloc.
 *    size   - New size wanted.
 *
 *  Returns:
 *    !NULL - pointer to memory block.
 *    NULL  - no array allocated.
 */

void * memRealloc_func ( void * memptr,
                         size_t size, char * filename, int line)
{
  double *              newmemptr;


  if (memptr == NULL) {
    newmemptr = memAlloc(size);
    PRINT_ALLOC(newmemptr, size, filename, line);
    return newmemptr;
  }

  PRINT_DEALLOC(memptr, filename, line);
  pthread_mutex_lock (&memallocmutexdat);         /* Lock local mutex */

  newmemptr = (double*) memptr;
  newmemptr --;
  memalloccurrent -= (unsigned long) newmemptr [0] ;
  newmemptr = realloc (newmemptr, size + sizeof (double));
  if (newmemptr != NULL) {
    newmemptr[0] = (double)size;
    memalloccurrent += (unsigned long) size;
    memallocmax = MAX (memallocmax, memalloccurrent);
    newmemptr ++;
#  ifdef MEMORY_TRACE
    if (memalloctraceflag)
      trace_malloc(trace_file, (clockGet()-trace_timestamp),
                   trace_procnum, STATE_ALLOC, memalloccurrent);
#  endif
  }
  else {
    perror("realloc");
    fprintf (stderr, "ALLOCATION FAILURE, used : %lu (%.3g %s), asked : %lu (%.3g %s)\n",
             memalloccurrent, MEMORY_WRITE(memalloccurrent), MEMORY_UNIT_WRITE(memalloccurrent),
	     (unsigned long)size, MEMORY_WRITE(size), MEMORY_UNIT_WRITE(size));
  }

  pthread_mutex_unlock (&memallocmutexdat);       /* Unlock local mutex */
  PRINT_ALLOC(newmemptr, size, filename, line);

  return ((void*)newmemptr);
}

/*
 *  Function: memFree
 *
 *  This is a thread-safe memory deallocation routine.
 *
 *  It returns:
 *    void - in all cases
 */
void memFree (void * memptr)
{
  double *              newmemptr;

  pthread_mutex_lock (&memallocmutexdat);         /* Lock local mutex */

  if (!(memptr == NULL)) {
    newmemptr = (double*)memptr;
    newmemptr --;
    memalloccurrent -= (unsigned long) newmemptr [0] ;
    free (newmemptr);
#  ifdef MEMORY_TRACE
    if (memalloctraceflag)
      trace_malloc(trace_file, (clockGet()-trace_timestamp),
                   trace_procnum, STATE_FREE, memalloccurrent);
#  endif
  }

  pthread_mutex_unlock (&memallocmutexdat);       /* Unlock local mutex */
}

#endif /* MEMORY_USAGE */

/*
 *  Function: memAllocGroup
 *
 *  This routine allocates a set of arrays in
 *  a single memAlloc()'ed array, the address
 *  of which is placed in the first argument.
 *  Arrays to be allocated are described as
 *  a duplet of ..., &ptr, size, ...,
 *  terminated by a NULL pointer.
 *
 *  Parameters:
 *    memptr - Pointer to first argument to allocate
 *    ...    - list of duplets &ptr, size, starting by size and ending by size.
 *
 *  Returns:
 *    !NULL - pointer to block, all arrays allocated.
 *    NULL  - no array allocated.
 */

void * memAllocGroup (void **  memptr,
                      ...)
{
  va_list     memlist;  /* Argument list of the call              */
  unsigned char    **  memloc;   /* Pointer to pointer of current argument */
  size_t      memoff;   /* Offset value of argument               */
  unsigned char    *   blkptr;   /* Pointer to memory chunk                */

  memoff = 0;
  memloc = (unsigned char **) memptr;             /* Point to first memory argument */
  va_start (memlist, memptr);            /* Start argument parsing         */
  while (memloc != NULL) {               /* As long as not NULL pointer    */
    memoff  = (memoff + (sizeof (double) - 1)) & (~ (sizeof (double) - 1));
    memoff += va_arg (memlist, size_t);
    memloc  = va_arg (memlist, unsigned char **);
  }

  if ((blkptr = (unsigned char *) memAlloc (memoff)) == NULL) {
    /* If cannot allocate   */
    *memptr = NULL;                      /* Set first pointer to NULL */
    return (NULL);
  }

  memoff = 0;
  memloc = (unsigned char **) memptr;             /* Point to first memory argument */
  va_start (memlist, memptr);            /* Restart argument parsing       */
  while (memloc != NULL) {               /* As long as not NULL pointer    */
    memoff  = (memoff + (sizeof (double) - 1)) &
      (~ (sizeof (double) - 1)); /* Pad  */
    *memloc = blkptr + memoff;           /* Set argument address           */
    memoff += va_arg (memlist, size_t);  /* Accumulate padded sizes        */
    /* Get next argument pointer      */
    memloc  = (unsigned char **) va_arg (memlist, void *);
  }

  return ((void *) blkptr);
}

/*
 *  Function: memAllocGroup
 *
 *  This routine reallocates a set of arrays in
 *  a single memRealloc()'ed array passed as
 *  first argument, and the address of which
 *  is placed in the second argument.
 *  Arrays to be allocated are described as
 *  a duplet of ..., &ptr, size, ...,
 *  terminated by a NULL pointer.
 *  WARNING: Because of memory alignment issues
 *  between int and double values, when arrays
 *  are not reallocated in place, offsets of
 *  arrays may vary, so that one should rather
 *  compute differences with respect to original
 *  offsets than rely on offsets returned by the
 *  routine. This routine should be used with
 *  extreme caution!
 *
 *  Parameters:
 *    oldptr - Pointer to first block to reallocate.
 *    ...    - list of duplets &ptr, size, starting by size and ending by size.
 *
 *  Returns:
 *    !NULL - pointer to block, all arrays allocated.
 *    NULL  - no array allocated.
*/
void * memReallocGroup (void * oldptr,
                        ...)
{
  va_list     memlist;    /* Argument list of the call              */
  unsigned char    **  memloc;     /* Pointer to pointer of current argument */
  size_t      memoff;     /* Offset value of argument               */
  unsigned char    *   blkptr;     /* Pointer to memory chunk                */

  memoff = 0;
  va_start (memlist, oldptr);      /* Start argument parsing */

  while ((memloc = va_arg (memlist, unsigned char **)) != NULL) {
    /* As long as not NULL pointer */
    /* Pad */
    memoff  = (memoff + (sizeof (double) - 1)) & (~ (sizeof (double) - 1));
    /* Accumulate padded sizes */
    memoff += va_arg (memlist, size_t);
  }
  /* Ensure reallocation of a padded area into a non-padded area will
   * never result in data loss */
  memoff += sizeof (double);

  /* If cannot allocate block */
  if ((blkptr = (unsigned char *) memRealloc (oldptr, memoff)) == NULL)
    return (NULL);

  memoff = 0;
  /* Restart argument parsing           */
  va_start (memlist, oldptr);
  while ((memloc = va_arg (memlist, unsigned char **)) != NULL) {
    /* As long as not NULL pointer */
    /* Pad */
    memoff  = (memoff + (sizeof (double) - 1)) & (~ (sizeof (double) - 1));
    /* Set argument address */
    *memloc = blkptr + memoff;
    /* Accumulate padded sizes */
    memoff += va_arg (memlist, size_t);
  }

  return ((void *) blkptr);
}

/*
 *  Function: memOffset
 *
 *  This routine computes the offsets of arrays
 *  of given sizes and types with respect to a
 *  given base address passed as first argument.
 *  Arrays the offsets of which are to be computed
 *  are described as a duplet of ..., &ptr, size, ...,
 *  terminated by a NULL pointer.
 *
 *  Parameters:
 *    memptr - Pointer to base address of memory area.
 *    ...    - list of duplets &ptr, size, starting by size and ending by size.
 *
 *  Returns:
 *    !NULL - in all cases, pointer to the end of the memory area.
 */
void * memOffset (void * memptr,
                  ...)
{
  va_list    memlist; /* Argument list of the call              */
  unsigned char    ** memloc;  /* Pointer to pointer of current argument */
  size_t     memoff;  /* Offset value of argument               */

  memoff = 0;
  va_start (memlist, memptr);                     /* Start argument parsing */

  while ((memloc = va_arg (memlist, unsigned char **)) != NULL) {
    /* As long as not NULL pointer */
    memoff  = (memoff + (sizeof (double) - 1)) & (~ (sizeof (double) - 1));
    *memloc = (unsigned char *) memptr + memoff;           /* Set argument address    */
    memoff += va_arg (memlist, size_t);           /* Accumulate padded sizes */
  }

  return ((void *) ((unsigned char *) memptr + memoff));
}
