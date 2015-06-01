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
** $Id: common.c 353 2005-11-04 15:57:47Z papin $
*/
/*
  File: common.c

  Part of a parallel direct block solver.

  These lines are common routines used   
  by all modules.                        

  Authors:
    Mathieu Faverge    - faverge@labri.fr
    David    GOUDIN     - .       
    Pascal   HENON      - henon@labri.fr
    Xavier   LACOSTE    - lacoste@labri.fr
    Francois PELLEGRINI - .
    Pierre   RAMET      - ramet@labri.fr 

  Dates:
    Version 0.0 - from : 08 may 1998
                  to     14 sep 1998
    Version 2.0 - from : 27 sep 2004
                  to     27 sep 2004
*/

/*
**  The defines and includes.
*/

#include "common_pastix.h"
#include <time.h>

#if __MINGW32__
#include "windows.h"
#define CLOCK_REALTIME 1
    LARGE_INTEGER
    getFILETIMEoffset()
    {
        SYSTEMTIME s;
        FILETIME f;
        LARGE_INTEGER t;

        s.wYear = 1970;
        s.wMonth = 1;
        s.wDay = 1;
        s.wHour = 0;
        s.wMinute = 0;
        s.wSecond = 0;
        s.wMilliseconds = 0;
        SystemTimeToFileTime(&s, &f);
        t.QuadPart = f.dwHighDateTime;
        t.QuadPart <<= 32;
        t.QuadPart |= f.dwLowDateTime;
        return (t);
    }

    int
    clock_gettime(int X, struct timespec *tv)
    {
        LARGE_INTEGER           t;
        FILETIME            f;
        double                  microseconds;
        static LARGE_INTEGER    offset;
        static double           frequencyToMicroseconds;
        static int              initialized = 0;
        static BOOL             usePerformanceCounter = 0;

        if (!initialized) {
            LARGE_INTEGER performanceFrequency;
            initialized = 1;
            usePerformanceCounter = QueryPerformanceFrequency(&performanceFrequency);
            if (usePerformanceCounter) {
                QueryPerformanceCounter(&offset);
                frequencyToMicroseconds = (double)performanceFrequency.QuadPart / 1000000.;
            } else {
                offset = getFILETIMEoffset();
                frequencyToMicroseconds = 10.;
            }
        }
        if (usePerformanceCounter) QueryPerformanceCounter(&t);
        else {
            GetSystemTimeAsFileTime(&f);
            t.QuadPart = f.dwHighDateTime;
            t.QuadPart <<= 32;
            t.QuadPart |= f.dwLowDateTime;
        }

        t.QuadPart -= offset.QuadPart;
        microseconds = (double)t.QuadPart / frequencyToMicroseconds;
        t.QuadPart = microseconds;
        tv->tv_sec = t.QuadPart / 1000000;
        tv->tv_nsec = t.QuadPart % 1000000;
        return (0);
    }
#endif

/*
  Function: clockGet

  Timing routine.
 
  Uses different timing routines depending on the machine architecture.

  Returns:
    Returns the time ellapsed since <clockStart>.
 */
double clockGet (void)
{
#if (defined X_ARCHalpha_compaq_osf1)
  struct rusage       data;
  getrusage (RUSAGE_SELF, &data);
  return (((double) data.ru_utime.tv_sec  + (double) data.ru_stime.tv_sec) +
          ((double) data.ru_utime.tv_usec + (double) data.ru_stime.tv_usec) * 
	  1.0e-6L);
#else
#if defined X_ARCHi686_mac
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double)tp.tv_sec + (double) tp.tv_usec * (double)1.0e-6L);
#else
  struct timespec tp;
  
  clock_gettime (CLOCK_REALTIME, &tp);            /* Elapsed time */

  return ((double) tp.tv_sec + (double) tp.tv_nsec * (double)1.0e-9L);

#endif
#endif
}

/*
  Function: usagePrint

  Usage printing routine.

  Prints usage into *stream* using strings inside *data* array.
  
  *data* array must be NULL terminated.

  Parameters:
    stream - File opened in write mode
    data   - tabular containing strings and NULL terminated.
 */
void
usagePrint (
FILE * const                stream,
const char ** const         data)
{
  const char **       cptr;

  fprintf (stream, "Usage is:\n");
  for (cptr = data; *cptr != NULL; cptr ++)
    fprintf (stream, "  %s\n", *cptr);
}

/*
  Function: api_iparmreader

  Reads integer parameters file from disk.
  
  The file must contain IPARM_SIZE lines starting with 
  the integer value corresponding.

  TODO: return values instead of exit(-1)...

  Parameters:
    filename - name of the file to read from.
    iparmtab - Array where to store parameters.

  Returns:
    1        - if file couldn't be read.
*/
int api_iparmreader (char * filename, 
		     PASTIX_INT *iparmtab)
{
  FILE*   m_File;
  int     i = 0;
  char    szbuff[MAX_CHAR_PER_LINE];
  char*   ret;
  char*   token;

#ifdef PASTIX_LOG
  fprintf(stderr, "-> api_iparmreader\n");
#endif
  m_File = fopen(filename,"rt");

  if(!m_File)
    {
#ifdef PASTIX_LOG
  fprintf(stderr, "<- api_iparmreader\n");
#endif
      return 1;
    }

  while(!feof(m_File) && i < IPARM_SIZE)
    { 
      ret = fgets(szbuff, MAX_CHAR_PER_LINE, m_File);
      if (ret == NULL)
	{
	  EXIT(MOD_UNKNOWN, UNKNOWN_ERR);
	}
      token = strtok(szbuff," ");
      iparmtab[i] = (PASTIX_INT)atol(token);
      i++;
    }

  fclose(m_File);

#ifdef PASTIX_LOG
  fprintf(stderr, "<- api_iparmreader\n");
#endif

#ifdef OOC
/*   if (iparmtab[IPARM_OOC_THREAD] > 1) */
    iparmtab[IPARM_OOC_THREAD] = 1;
#endif
  return 0;
}

/*
  Function: api_dparmreader

  Reads double parameters file from disk.
  
  The file must contain IPARM_SIZE lines starting with 
  the double value corresponding.

  See *atof* manual for the format required.

  TODO: return values instead of exit(-1)...

  Parameters:
    filename - name of the file to read from.
    dparmtab - Array where to store parameters.

  Returns:
    1        - if file couldn't be read.
*/
int api_dparmreader(char * filename, 
		    double *dparmtab)
{
  FILE*   m_File;
  int     i = 0;
  char    szbuff[MAX_CHAR_PER_LINE];
  char*   ret;
  char*   token;

#ifdef PASTIX_LOG
  fprintf(stderr, "-> api_dparmreader\n");
#endif
  m_File = fopen(filename,"rt");

  if(!m_File) return 1;

  while(!feof(m_File) && i < DPARM_SIZE)
    {   
      ret = fgets(szbuff, MAX_CHAR_PER_LINE, m_File);
      if (ret == NULL)
	{
	  EXIT(MOD_UNKNOWN, UNKNOWN_ERR);
	}
      token = strtok(szbuff," ");
      dparmtab[i] = atof(token);
      i++;
    }
  fclose(m_File);

#ifdef PASTIX_LOG
  fprintf(stderr, "<- api_dparmreader\n");
#endif
  return 0;
}

/*
  Function: api_dumparm

  Dump PaStiX parameters arrays to disk.

  Parameters:
    stream - File opened in write mode
    iparm  - integer parameters array
    dparm  - floating parameters array
  
 */
void api_dumparm(FILE *stream, 
		 PASTIX_INT *iparm, 
		 double *dparm)
{
  PASTIX_INT i;

  for (i=0; i<IPARM_SIZE; i++)
    {
      fprintf(stream, "iparm[%ld] = %ld\n", (long) i, (long) iparm[i]);
    }
  fprintf(stream, "----\n");
  for (i=0; i<DPARM_SIZE; i++)
    {
      fprintf(stream, "dparm[%ld] = %e\n", (long) i, dparm[i]);
    }
}

void set_iparm(PASTIX_INT *iparm, enum IPARM_ACCESS offset, PASTIX_INT value)
{
  if (iparm != NULL) iparm[offset] = (PASTIX_INT)value;
}

void set_dparm(double *dparm, enum DPARM_ACCESS offset, double value)
{
  if (dparm != NULL) dparm[offset] = (double)value;
}
/*******************************************************************************
 * Section: Sort functions
 */



/*
   Function: qsortIntFloatAsc

   Sort 2 arrays simultaneously, the first array is an
   array of PASTIX_INT and used as key for sorting.
   The second array is an array of PASTIX_FLOAT.

   Parameters:
     pbase       - Array of pointers to the first element of each array to sort.
     total_elems - Number of element in each array.

   Returns:
     Nothing

*/
#define INTSORTNAME            qsortIntFloatAsc
#define INTSORTSIZE(x)         ((x==0)?(sizeof (PASTIX_INT)):(          \
                                  (PASTIX_INT)(*(PASTIX_INT*)(*(pbase+2))) \
                                  *sizeof(PASTIX_FLOAT)))
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {                                     \
    PASTIX_INT     t;                                                   \
    long    disp_p   = (((PASTIX_INT*)p)-((PASTIX_INT*)base_ptr));      \
    long    disp_q   = (((PASTIX_INT*)q)-((PASTIX_INT*)base_ptr));      \
    PASTIX_FLOAT * floatptr = *(pbase+1);                               \
    PASTIX_FLOAT   f;                                                   \
    PASTIX_INT     dof, ndof2= (PASTIX_INT)(*(PASTIX_INT*)(*(pbase+2))); \
    /* swap integers */                                                 \
    t = *((PASTIX_INT *) (p));                                          \
    *((PASTIX_INT *) (p)) = *((PASTIX_INT *) (q));                      \
    *((PASTIX_INT *) (q)) = t;                                          \
    /* swap corresponding values */                                     \
    for (dof = 0; dof < ndof2; dof++) {                                 \
      f = floatptr[disp_p*ndof2+dof];                                   \
      floatptr[disp_p*ndof2+dof] = floatptr[disp_q*ndof2+dof];          \
      floatptr[disp_q*ndof2+dof] = f;                                   \
    }                                                                   \
  } while (0)
#define INTSORTCMP(p,q)             (*((PASTIX_INT *) (p)) < *((PASTIX_INT *) (q)))
#include "common_sort2.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB

/*
   Function: qsort2IntFloatAsc

   Sort 3 arrays simultaneously, the first array is an
   array of PASTIX_INT and used as primary key for sorting.
   The second array is an other array of PASTIX_INT used
   as secondary key.
   The third array is an array of PASTIX_FLOAT.

   Parameters:
     pbase       - Array of pointers to the first element of each array to sort.
     total_elems - Number of element in each array.

   Returns:
     Nothing

*/
#define INTSORTNAME            qsort2IntFloatAsc
#define INTSORTSIZE(x)         ((x<2)?(sizeof (PASTIX_INT)):(sizeof (PASTIX_FLOAT)))
#define INTSORTNTAB            3
#define INTSORTSWAP(p,q)       do {					\
    PASTIX_INT     t;								\
    long    disp_p   = (((PASTIX_INT*)p)-((PASTIX_INT*)base_ptr));			\
    long    disp_q   = (((PASTIX_INT*)q)-((PASTIX_INT*)base_ptr));			\
    PASTIX_INT   * int2ptr  = *(pbase+1);					\
    PASTIX_FLOAT * floatptr = *(pbase+2);					\
    PASTIX_FLOAT   f;								\
    /* swap integers */							\
    t = *((PASTIX_INT *) (p));							\
    *((PASTIX_INT *) (p)) = *((PASTIX_INT *) (q));					\
    *((PASTIX_INT *) (q)) = t;							\
    /* swap on secont integer array */					\
    t = int2ptr[disp_p];						\
    int2ptr[disp_p] = int2ptr[disp_q];					\
    int2ptr[disp_q] = t;			 			\
    /* swap corresponding values */					\
    f = floatptr[disp_p];						\
    floatptr[disp_p] = floatptr[disp_q];				\
    floatptr[disp_q] = f;						\
  } while (0)
#define INTSORTCMP(p,q)  ((*((PASTIX_INT *) (p)) < *((PASTIX_INT *) (q))) ||		\
			  ((*((PASTIX_INT *) (p)) == *((PASTIX_INT *) (q))) &&	\
			   ((( PASTIX_INT *)(*(pbase+1)))[(((PASTIX_INT*)p)-((PASTIX_INT*)base_ptr))] < \
			    (( PASTIX_INT *)(*(pbase+1)))[(((PASTIX_INT*)q)-((PASTIX_INT*)base_ptr))])))
#include "common_sort2.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB

/*
   Function: qsort2IntAsc

   Sort 2 arrays simultaneously, the first array is an
   array of PASTIX_INT and used as primary key for sorting.
   The second array is an other array of PASTIX_INT used
   as secondary key.

   Parameters:
     pbase       - Array of pointers to the first element of each array to sort.
     total_elems - Number of element in each array.

   Returns:
     Nothing

*/
#define INTSORTNAME            qsort2IntAsc
#define INTSORTSIZE(x)         (sizeof (PASTIX_INT))
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {                                     \
    PASTIX_INT     t;                                                   \
    long    disp_p   = (((PASTIX_INT*)p)-((PASTIX_INT*)base_ptr));			\
    long    disp_q   = (((PASTIX_INT*)q)-((PASTIX_INT*)base_ptr));			\
    PASTIX_INT   * int2ptr  = *(pbase+1);                               \
    /* swap integers */                                                 \
    t = *((PASTIX_INT *) (p));                                          \
    *((PASTIX_INT *) (p)) = *((PASTIX_INT *) (q));                      \
    *((PASTIX_INT *) (q)) = t;                                          \
    /* swap on secont integer array */                                  \
    t = int2ptr[disp_p];                                                \
    int2ptr[disp_p] = int2ptr[disp_q];                                  \
    int2ptr[disp_q] = t;                                                \
  } while (0)
#define INTSORTCMP(p,q)  ((*((PASTIX_INT *) (p)) < *((PASTIX_INT *) (q))) || \
                          ((*((PASTIX_INT *) (p)) == *((PASTIX_INT *) (q))) && \
                           ((( PASTIX_INT *)(*(pbase+1)))[(((PASTIX_INT*)p)-((PASTIX_INT*)base_ptr))] < \
                            (( PASTIX_INT *)(*(pbase+1)))[(((PASTIX_INT*)q)-((PASTIX_INT*)base_ptr))])))
#include "common_sort2.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB

/*
   Function: qsort2SmallIntAsc

   Sort 2 arrays simultaneously, the first array is an
   array of integers (int) and used as primary key for sorting.
   The second array is an other array of int used
   as secondary key.

   Parameters:
     pbase       - Array of pointers to the first element of each array to sort.
     total_elems - Number of element in each array.

   Returns:
     Nothing

*/
#define INTSORTNAME            qsort2SmallIntAsc
#define INTSORTSIZE(x)         (sizeof (int))
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {					\
    int     t;								\
    long    disp_p   = (((int*)p)-((int*)base_ptr));			\
    long    disp_q   = (((int*)q)-((int*)base_ptr));			\
    int   * int2ptr  = *(pbase+1);						\
    /* swap integers */							\
    t = *((int *) (p));							\
    *((int *) (p)) = *((int *) (q));					\
    *((int *) (q)) = t;							\
    /* swap on secont integer array */					\
    t = int2ptr[disp_p];						\
    int2ptr[disp_p] = int2ptr[disp_q];					\
    int2ptr[disp_q] = t;						\
  } while (0)
#define INTSORTCMP(p,q)  ((*((int *) (p)) < *((int *) (q))) ||		\
			  ((*((int *) (p)) == *((int *) (q))) &&	\
			   ((( int *)(*(pbase+1)))[(((int*)p)-((int*)base_ptr))] < \
			    (( int *)(*(pbase+1)))[(((int*)q)-((int*)base_ptr))])))
#include "common_sort2.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB
