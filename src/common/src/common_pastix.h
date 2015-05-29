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
 *  File: common_pastix.h
 *
 *  Part of a parallel direct block solver.
 *
 *  These lines are the common data
 *  declarations for all modules.
 *
 *  Authors:
 *    Mathieu Faverge    - faverge@labri.fr
 *    David    GOUDIN     - .
 *    Pascal   HENON      - henon@labri.fr
 *    Xavier   LACOSTE    - lacoste@labri.fr
 *    Francois PELLEGRINI - .
 *    Pierre   RAMET      - ramet@labri.fr
 *
 *  Dates:
 *    Version 0.0 - from 08 may 1998
 *                  to   08 jan 2001
 *    Version 1.0 - from 06 jun 2002
 *                  to   06 jun 2002
 */
#ifndef COMMON_PASTIX_H
#define COMMON_PASTIX_H
#include "api.h"
#include "debug.h"
#include "errors.h"

#if (defined DYNSCHED && defined WITH_STARPU)
#  error "STARPU and Dynsched are not compatible"
#endif

#ifdef __INTEL_COMPILER
/* Ignore icc remark : "operands are evaluated in unspecified order"*/
#  pragma warning(disable:981)
/* Ignore icc remark : "external function definition with no prior declaration" */
#  pragma warning(disable:1418)
/* Ignore icc remark : "external declaration in primary source file" */
#  pragma warning(disable:1419)
/* Ignore icc remark : " parameter "arg" was never referenced" */
#  pragma warning(disable:869)
/* Ignore icc remark : "variable "size" was set but never used" */
#  pragma warning(disable:593)
/* Ignore icc remark : "floating-point equality and inequality comparisons are unreliable" */
#  pragma warning(disable:1572)
/* Ignore icc remark : "statement is unreachable" */
#  pragma warning(disable:111)
#endif

#ifdef OOC_FTGT
#  ifndef OOC_FTGT_RESET
#    define OOC_FTGT_RESET
#  endif
#  ifndef OOC
#    define OOC
#  endif
#endif

#ifdef OOC
#  ifndef MEMORY_USAGE
#    define MEMORY_USAGE
#  endif
#endif


/*
** Machine configuration values.
** The end of the X_ARCH variable is built with parts of the
** `uname -m`, `uname -r`, and `uname -s` commands.
*/

#define X_C_NORESTRICT
#ifndef X_C_NORESTRICT
#  define X_C_RESTRICT
#endif /* X_C_NORESTRICT */

#if (defined X_ARCHi686_mac)
#  define X_ARCHi686_pc_linux
#endif

#if (defined X_ARCHpower_ibm_aix)
#  define X_INCLUDE_ESSL
#  undef  X_C_RESTRICT
#endif /* (defined X_ARCHpower_ibm_aix) */

#if (defined X_ARCHalpha_compaq_osf1)
#  define restrict
/*#define volatile*/
#endif /* (defined X_ARCHalpha_compaq_osf1) */


/*
** Compiler optimizations.
*/

#ifdef X_C_RESTRICT
#  ifdef __GNUC__
#    define restrict                    __restrict
#  endif /* __GNUC__ */
#else /* X_C_RESTRICT */
#  define restrict
#endif /* X_C_RESTRICT */

/*
** The includes.
*/

/* Redefinition de malloc,free,printf,fprintf */
#if defined(MARCEL) || defined(_MINGW_)
#  include <pthread.h>
#endif

#include            <ctype.h>
#include            <math.h>
#ifdef X_ARCHi686_mac
#  include            <malloc/malloc.h>
#else /* X_ARCHi686_mac */
#  ifdef __FreeBSD__
#    include            <stdlib.h>
#  else /* not __FreeBSD__ */
#    include            <malloc.h>
#  endif /* not __FreeBSD__ */
#endif /* X_ARCHi686_mac */
#include            <memory.h>
#include            <stdio.h>
#include            <stdarg.h>
#include            <stdlib.h>
#include            <string.h>
#include            <time.h>         /* For the effective calls to clock () */
#include            <limits.h>
#include            <sys/types.h>
#include            <sys/time.h>
#ifdef WIN32
	#include            <windows_fix.h>
#else
	#include            <sys/resource.h>	
#endif
#include            <unistd.h>
#include            <float.h>
#include            <stdint.h>
#ifdef X_INCLUDE_ESSL
#  include            <essl.h>
#endif /* X_INCLUDE_ESSL */

#ifdef X_ASSERT
#  include <assert.h>
#endif /* X_ASSERT */

#ifndef MIN
#  define MIN(x,y) (((x)<(y))?(x):(y))
#endif

#ifndef MAX
#  define MAX(x,y) (((x)<(y))?(y):(x))
#endif

/*
 * Checking incompatible options.
 */
#ifdef FORCE_NOMPI
#  ifdef DISTRIBUTED
#    error "-DFORCE_NOMPI is not compatible with -DDISTRIBUTED"
#  endif
#endif

/*
**  Handling of generic types.
*/


#ifdef INTSIZE32
#  ifndef FORCE_INT32
#    define FORCE_INT32
#  endif

#endif

#if (defined INTSIZE64 || defined INTSSIZE64)
#  ifndef FORCE_INT64
#    define FORCE_INT64
#  endif
#endif

#ifdef FORCE_INT32
#  ifndef INTSIZE32
#    define INTSIZE32
#  endif
#endif


#ifdef FORCE_INT64
#  if !(defined INTSIZE64) && !(defined INTSSIZE64)
#    define INTSIZE64
#  endif
#endif


#ifdef PREC_DOUBLE
#  define BLAS_DOUBLE
#  define BASE_FLOAT double
#else
#  define BASE_FLOAT float
#endif

#ifdef TYPE_COMPLEX
#  define CPLX
#else
#  define creal(x) x
#  define cimag(x) 0.0
#endif

#ifdef CPLX
#  if (defined X_ARCHalpha_compaq_osf1)

#    ifndef USE_CXX

#      ifndef   _RWSTD_HEADER_REQUIRES_HPP
#        include <complex>
#      else  /* _RWSTD_HEADER_REQUIRES_HPP */
#        include <complex.hpp>
#      endif /* _RWSTD_HEADER_REQUIRES_HPP */

#      define PASTIX_FLOAT complex<BASE_FLOAT>

#      ifdef    PREC_DOUBLE
#        define COMM_FLOAT MPI_DOUBLE_COMPLEX
#        define FLOAT_MAX DBL_MAX
#      else  /* PREC_DOUBLE */
#        define COMM_FLOAT MPI_COMPLEX
#        define FLOAT_MAX MAXFLOAT
#      endif /* PREC_DOUBLE */

#      define COMM_SUM GetMpiSum()
#      define ABS_FLOAT(x) abs(x)
#      define fabs(x) abs(x)
#      define cabs(x) abs(x)
#      define csqrt(x) sqrt(x)
#      define CONJ_FLOAT(x) conj(x)
#      define creal(x) real(x)
#      define cimag(x) imag(x)
#    endif /*USE_CXX*/

#  else /*X_ARCHalpha_compaq_osf1*/
#    include <complex.h>
#    define COMM_FLOAT GetMpiType()
#    define COMM_SUM GetMpiSum()

#    ifdef    PREC_DOUBLE
#      define PASTIX_FLOAT double complex
#      define ABS_FLOAT(x) cabs(x)
#      ifdef    _DCMPLX
#        define BLAS_FLOAT dcmplx
#      else  /* _DCMPLX */
#        define BLAS_FLOAT double complex
#      endif /* _DCMPLX */
#      define CONJ_FLOAT(x) conj(x)
#      define FLOAT_MAX DBL_MAX
#    else  /* PREC_DOUBLE */
#      define PASTIX_FLOAT float complex
#      define ABS_FLOAT(x) cabsf(x)
#      ifdef    _CMPLX
#        define BLAS_FLOAT cmplx
#      else  /* _CMPLX */
#        define BLAS_FLOAT float complex
#      endif /* _CMPLX */
#      define CONJ_FLOAT(x) conjf(x)
#      define FLOAT_MAX MAXFLOAT
#    endif /* PREC_DOUBLE */

#  endif /* X_ARCHalpha_compaq_osf1 */
#else /* CPLX */
#  define COMM_SUM MPI_SUM

#  define CONJ_FLOAT(x) x
#  ifdef PREC_DOUBLE
#    define PASTIX_FLOAT double
#    define ABS_FLOAT(x) fabs(x)
#    define COMM_FLOAT MPI_DOUBLE
#    define FLOAT_MAX  DBL_MAX
#  else /* PREC_DOUBLE */
#    define PASTIX_FLOAT float
#    define FLOAT_MAX MAXFLOAT
#    define COMM_FLOAT MPI_FLOAT
#    define ABS_FLOAT(x) fabsf(x)
#  endif /* PREC_DOUBLE */
#endif /* CPLX */

#ifndef BLAS_FLOAT
#  define BLAS_FLOAT PASTIX_FLOAT
#endif

#ifdef PREC_DOUBLE
#  define BLAS_REAL double
#else
#  define BLAS_REAL float
#endif

/*
 *  D�finition de la taille des entiers utilis�s
 */
#ifdef FORCE_LONG
#  define PASTIX_INT           long          /* Long integer type */
#  define PASTIX_UINT          unsigned long
#  define COMM_INT      MPI_LONG
#elif (defined FORCE_INT32)
#  define PASTIX_INT           int32_t
#  define PASTIX_UINT          uint32_t
#  define COMM_INT      MPI_INTEGER4
#elif (defined FORCE_INT64)
#  define PASTIX_INT           int64_t
#  define PASTIX_UINT          uint64_t
#  define COMM_INT      MPI_INTEGER8
#else
#  define PASTIX_INT           int           /* Default integer type     */
#  define PASTIX_UINT          unsigned int
#  define COMM_INT      MPI_INT       /* Generic MPI integer type */
#endif

#ifndef INTSIZEBITS
#  define INTSIZEBITS   (sizeof (PASTIX_INT) << 3)
#endif /* INTSIZEBITS */

#define INTVALMAX     ((PASTIX_INT) (((PASTIX_UINT) 1 << (INTSIZEBITS - 1)) - 1))

#include "redefine_functions.h"

#define MEMORY_WRITE(mem) ( ((mem) < 1<<10) ?                           \
                            ( (double)(mem) ) :                         \
                            ( ( (mem) < 1<<20 ) ?                       \
                              ( (double)(mem)/(double)(1<<10) ) :       \
                              ( ((mem) < 1<<30 ) ?                      \
                                ( (double)(mem)/(double)(1<<20) ) :     \
                                ( (double)(mem)/(double)(1<<30) ))))
#define MEMORY_UNIT_WRITE(mem) (((mem) < 1<<10) ?                       \
                                "o" :                                   \
                                ( ( (mem) < 1<<20 ) ?                   \
                                  "Ko" :                                \
                                  ( ( (mem) < 1<<30 ) ?                 \
                                    "Mo" :                              \
                                    "Go" )))

#define PRINT_FLOPS(flops) ( ((flops) < 1<<10) ?                        \
                             ( (double)(flops) ) :                      \
                             ( ( (flops) < 1<<20 ) ?                    \
                               ( (double)(flops)/(double)(1<<10) ) :    \
                               ( ((flops) < 1<<30 ) ?                   \
                                 ( (double)(flops)/(double)(1<<20) ) :  \
                                 ( (double)(flops)/(double)(1<<30) ))))
#define PRINT_FLOPS_UNIT(flops) ( ((flops) < 1<<10) ?                   \
                                  ( "FLOPS" ) :                         \
                                  ( ( (flops) < 1<<20 ) ?               \
                                    ( "KFLOPS" ) :                      \
                                    ( ((flops) < 1<<30 ) ?              \
                                      ( "MFLOPS" ) :                    \
                                      ( "GFLOPS" ))))
#ifdef PRINT_ALL_MALLOC
#  define PRINT_ALLOC(ptr, size, file, line) do                         \
    {                                                                 \
      fprintf(stdout, "%s:%d allocation size %.3g %s : %p (%s)\n",    \
              file, line, MEMORY_WRITE((size)),                       \
              MEMORY_UNIT_WRITE((size)), ptr, #ptr);                  \
    } while(0)

#  define PRINT_DEALLOC(ptr, file, line) do                             \
    {                                                                   \
      if (ptr != NULL) {                                                \
        double * pda_memptr = (double*)(ptr);                           \
        unsigned long pda_size;                                         \
        pda_memptr --;                                                  \
        pda_size = (unsigned long) pda_memptr[0];                       \
        fprintf(stdout, "%s:%d"                                         \
                " deallocation size %.3g %s : %p (%s)\n",               \
                __FILE__, __LINE__,                                     \
                MEMORY_WRITE(pda_size),                                 \
                MEMORY_UNIT_WRITE(pda_size), ptr, #ptr);                \
      }                                                                 \
    } while(0)
#  define PRINT_DEALLOC_EXT(ptr, size, file, line) do                   \
    {                                                                   \
      fprintf(stdout, "%s:%d"                                           \
              " deallocation size %.3g %s : %p (%s)\n",                 \
              __FILE__, __LINE__,                                       \
              MEMORY_WRITE(size),                                       \
              MEMORY_UNIT_WRITE(size), ptr, #ptr);                      \
    } while(0)
#else /* not PRINT_ALL_MALLOC */
#  define PRINT_ALLOC(ptr, size, file, line) do   \
    {                                           \
    } while(0)

#  define PRINT_DEALLOC(ptr, file, line) do       \
    {                                           \
    } while(0)
#  define PRINT_DEALLOC_EXT(ptr, size, file, line) do   \
    {                                                   \
    } while(0)
#endif /* not PRINT_ALL_MALLOC */

/* Working definitions. */
#define memAlloca(size)                alloca(size)
#ifndef MEMORY_USAGE   /* TODO : remove mutex protection in multi-thread mode */
#  ifdef X_ARCHpower_ibm_aix
#    define memAlloc(size)                 mymalloc(size, __FILE__,__LINE__)
#  else
#    define memAlloc(size)                 malloc(size)
#  endif
#  define memAlloca(size)                alloca(size)
#  define memRealloc(ptr,size)           realloc((ptr),(size))
#  define memFree(ptr)                   ( free ((char *) (ptr)) , 0)
#  define PTR_MEMSIZE(ptr) (size_t)(0)
#else
#  define memAlloc(size)                 (memAlloc_func(size,__FILE__,__LINE__))
#  define PTR_MEMSIZE(ptr) ((ptr==NULL)?0:((size_t)(((double*)(ptr))[-1])))
#endif
#define memFree_null(ptr)                  do   \
    {                                           \
      PRINT_DEALLOC(ptr, __FILE__, __LINE__);   \
      memFree ((char *) (ptr));                 \
      (ptr) = NULL;                             \
    } while(0)

/* Freeing function if memAlloca implemented as malloc */
#define memFreea(ptr,module)           (0)
#define memSet(ptr,val,siz)                memset((ptr),(val),(siz));
#define memCpy(dst,src,siz)                memcpy((dst),(src),(siz));
#define memMov(dst,src,siz)                memmove((dst),(src),(siz));

#ifdef WARNINGS_MALLOC
#  define CHECK_ALLOC(ptr, size, type)                                  \
  do {                                                                \
    if (ptr != NULL)                                                  \
      errorPrintW("non NULL pointer in allocation (line=%d,file=%s)", \
                  __LINE__,__FILE__);                                 \
    if ((size) * sizeof(type))                                        \
      errorPrintW("Allocation of size 0 (line=%d,file=%s)",           \
                  __LINE__,__FILE__);                                 \
  } while (0)
#else  /* not WARNINGS_MALLOC */
#  define CHECK_ALLOC(ptr, size, type)                                  \
  do {                                                                \
  } while (0)

#endif /* not WARNINGS_MALLOC */
/*
 * Macro: MALLOC_INTOREXTERN
 *
 * Choose between <MALLOC_INTERN> and <MALLOC_EXTERN>
 * following flag_int.
 *
 * Parameters:
 *   ptr      - address where to allocate.
 *   size     - Number of elements to allocate.
 *   types    - Type of the elements to allocate.
 *   flag_int - API_YES for internal allocation, API_NO for external.
 */
#define MALLOC_INTOREXTERN(ptr, size, type, flag_int) \
  do {                                                \
    if (flag_int == API_YES)                          \
      {                                               \
        MALLOC_INTERN(ptr, size, type);               \
      }                                               \
    else                                              \
      {                                               \
        MALLOC_EXTERN(ptr, size, type);               \
      }                                               \
  } while (0)

#define FREE_NULL_INTOREXT(ptr, flag_int)         \
  do {                                            \
    if (flag_int == API_YES)                      \
      {                                           \
        memFree_null(ptr);                        \
      }                                           \
    else                                          \
      {                                           \
        free(ptr);                                \
        ptr = NULL;                               \
      }                                           \
  } while (0)
/*
 * Macro: MALLOC_EXTERN
 *
 * Allocate a space of size *size* x sizeof(*type*)
 * at the adress indicated by ptr, using external *malloc*.
 *
 * Parameters:
 *   ptr   - address where to allocate.
 *   size  - Number of elements to allocate.
 *   types - Type of the elements to allocate.
 */
#define MALLOC_EXTERN(ptr, size, type)                                  \
  do {                                                                  \
    CHECK_ALLOC(ptr, size, type);                                       \
    if (((long long)(size) * sizeof(type)) == 0)                        \
      {                                                                 \
        ptr = NULL;                                                     \
      }                                                                 \
    else                                                                \
      {                                                                 \
        if (((long long)(((size)*sizeof(type))) < 0) ||                 \
            NULL == (ptr = (type *) malloc((size) * sizeof(type))))     \
          {                                                             \
            MALLOC_ERROR(#ptr);                                         \
          }                                                             \
      }                                                                 \
  } while(0)


/*
 * Macro: MALLOC_INTERN
 *
 * Allocate a space of size *size* x sizeof(*type*)
 * at the adress indicated by ptr, using internal *memAlloc*.
 *
 * Parameters:
 *   ptr   - address where to allocate.
 *   size  - Number of elements to allocate.
 *   types - Type of the elements to allocate.
 */
#define MALLOC_INTERN(ptr, size, type)                                  \
  {                                                                     \
    CHECK_ALLOC(ptr, size, type);                                       \
    if (((long long)(size) * sizeof(type)) == 0)                        \
      {                                                                 \
        ptr = NULL;                                                     \
      }                                                                 \
    else                                                                \
      {                                                                 \
        if ((((long long) ((size) * sizeof(type))) < 0) ||              \
            (NULL == (ptr = (type *) memAlloc((size) * sizeof(type))))) \
          {                                                             \
            MALLOC_ERROR(#ptr);                                         \
          }                                                             \
      }                                                                 \
    PRINT_ALLOC(ptr, (size*sizeof(type)), __FILE__,__LINE__);           \
  }

/*
  Macro: PASTIX_FOPEN

  Open a file and handle errors.

  Parameters:
  FILE      - Stream (FILE*) to link to the file.
  filenamne - String containing the path to the file.
  mode      - String containing the opening mode.

*/
#define PASTIX_FOPEN(FILE, filenamne, mode)                                \
  {                                                                 \
    FILE = NULL;                                                    \
    if (NULL == (FILE = fopen(filenamne, mode)))                    \
      {                                                             \
        errorPrint("%s:%d Couldn't open file : %s with mode %s\n",  \
                   __FILE__, __LINE__, filenamne, mode);            \
        EXIT(MOD_UNKNOWN,FILE_ERR);                                 \
      }                                                             \
  }
/*
  Macro: PASTIX_FREAD

  Calls fread function and test his return value

  Parameters:
  buff   - Memory area where to copy read data.
  size   - Size of an element to read.
  count  - Number of elements to read
  stream - Stream to read from
*/
#define PASTIX_FREAD(buff, size, count, stream)        \
  {                                             \
    if ( 0 == fread(buff, size, count, stream)) \
      {                                         \
        errorPrint("%s:%d fread error\n",       \
                   __FILE__, __LINE__);         \
        EXIT(MOD_UNKNOWN,FILE_ERR);             \
      }                                         \
  }
/*
 * Other working definitions
 */

#define MAX_CHAR_PER_LINE 1000


/*
**  Handling of timers.
*/

/** The clock type. **/

typedef struct Clock_ {
  double                    time[2];    /*+ The start and accumulated times +*/
} Clock;

/*
**  Handling of files.
*/

/** The file structure. **/

typedef struct File_ {
  char *                    name;                 /*+ File name    +*/
  FILE *                    pntr;                 /*+ File pointer +*/
  char *                    mode;                 /*+ Opening mode +*/
} File;

/*
**  The function prototypes.
*/

#ifdef X_ARCHalpha_compaq_osf1
#  ifndef USE_CXX
extern "C" {
#  endif
#endif

#ifdef MEMORY_USAGE
  void *         memAlloc_func       (size_t,char*,int);
  void *         memRealloc_func     (void *, size_t, char*, int);
#  define memRealloc(ptr,size)                      \
  memRealloc_func(ptr, size, __FILE__, __LINE__)

  void           memFree             (void *);
  unsigned long  memAllocGetCurrent  (void);
  unsigned long  memAllocGetMax      (void);
  void           memAllocTraceReset  (void);
#else
  void *         mymalloc            (size_t,char*,int);
#endif /* MEMORY_USAGE */
#ifdef MEMORY_TRACE
  void           memAllocTrace       (FILE *, double, int);
  void           memAllocUntrace     ();
#else
#  define          memAllocTrace(a, b, c) {}
#  define          memAllocUntrace()      {}
#endif
  void *         memAllocGroup       (void **, ...);
  void *         memReallocGroup     (void *, ...);
  void *         memOffset           (void *, ...);

  void           usagePrint          (FILE * const, const char ** const);

  void           errorProg           (const char * const);
  void           errorPrint          (const char * const, ...);
  void           errorPrintW         (const char * const, ...);

  int            intLoad             (FILE * const, PASTIX_INT * const);
  int            intSave             (FILE * const, const PASTIX_INT);
  void           intAscn             (PASTIX_INT * restrict const,
                                      const PASTIX_INT, const PASTIX_INT);
  void           intPerm             (PASTIX_INT * restrict const, const PASTIX_INT);
  void           intRandInit         (void);
  PASTIX_INT            intRandVal          (PASTIX_INT);
  void           intSort1asc1        (void * const, const PASTIX_INT);
  void           intSort2asc1        (void * const, const PASTIX_INT);
  void           intSort2asc2        (void * const, const PASTIX_INT);
  void           intSort3asc1dsc2    (void * const, const PASTIX_INT);

  void           clockInit           (Clock * const);
  void           clockStart          (Clock * const);
  void           clockStop           (Clock * const);
  double         clockVal            (Clock * const);
  double         clockGet            (void);

#ifdef X_ARCHalpha_compaq_osf1
#  ifndef USE_CXX
}
#  endif
#endif

/*
**  The macro definitions.
*/

#define clockInit(clk)        ((clk)->time[0]  = (clk)->time[1] = 0)
#define clockStart(clk)       ((clk)->time[0]  = clockGet ())
#define clockStop(clk)        ((clk)->time[1]  = clockGet ())
#define clockVal(clk)         ((clk)->time[1] - (clk)->time[0])

#define intRandVal(ival)      ((PASTIX_INT) (((PASTIX_UINT) random ()) % ((PASTIX_UINT) (ival))))

#define FORTRAN(nu,nl,pl,pc)                    \
  void nu pl;                                   \
  void nl pl                                    \
  { nu pc; }                                    \
  void nl##_ pl                                 \
  { nu pc; }                                    \
  void nl##__ pl                                \
  { nu pc; }                                    \
  void nu pl

#ifdef MARCEL
#  define marcel_printf(...) do {  } while(0)
/* #define marcel_printf(...) marcel_fprintf(stderr, __VA_ARGS__) */
#else
#  define printf(...) do {  } while(0)
/* #define printf(...) fprintf(stderr, __VA_ARGS__) */
#endif

void api_dumparm(FILE *stream, PASTIX_INT *iparm, double *dparm);
int  api_dparmreader(char * filename, double *dparmtab);
int  api_iparmreader(char * filename, PASTIX_INT    *iparmtab);
void set_iparm(PASTIX_INT    *iparm, enum IPARM_ACCESS offset, PASTIX_INT    value);
void set_dparm(double *dparm, enum DPARM_ACCESS offset, double value);

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
void qsortIntFloatAsc(void ** const pbase,
                      const PASTIX_INT     total_elems);

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
void qsort2IntFloatAsc(void ** const pbase,
                       const PASTIX_INT     total_elems);


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
void qsort2IntAsc(void ** const pbase,
                  const PASTIX_INT     total_elems);

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
void qsort2SmallIntAsc(void ** const pbase,
                       const PASTIX_INT     total_elems);



/*
 * Macro to write in file pastix.pid/fname
 * Don't forget to close the FILE out after
 */
#ifdef WIN32
#define OUT_OPENFILEINDIR(iparm, file, fname, mode)               \
  {                                                               \
    char  outdirname[255];                                        \
    char  outfilename[255];                                       \
    sprintf(outdirname, "./pastix.%d", (int)(iparm)[IPARM_PID]);  \
    mkdir(outdirname);                                      	  \
    sprintf(outfilename, "%s/%s", outdirname, fname);             \
    file = fopen(outfilename, mode);                              \
  }
#else  
#define OUT_OPENFILEINDIR(iparm, file, fname, mode)               \
  {                                                               \
    char  outdirname[255];                                        \
    char  outfilename[255];                                       \
    sprintf(outdirname, "./pastix.%d", (int)(iparm)[IPARM_PID]);  \
	mkdir(outdirname, 0755);                                      \
    sprintf(outfilename, "%s/%s", outdirname, fname);             \
    file = fopen(outfilename, mode);                              \
  }
#endif
#define OUT_CLOSEFILEINDIR(file) fclose(file);

#define PASTIX_MASK_ISTRUE(var, mask) (var == (var | mask))
/*
 * macro: CHECK_MPI
 *
 * Perform MPI call, check for error, print an error message
 * and call MPI_Abort.
 *
 */
#ifdef FORCE_NOMPI
#  define CHECK_MPI(call) call
#  define CHECK_THREAD_LEVEL(THREAD_MODE) do {  \
    THREAD_MODE = THREAD_MODE;                  \
  } while (0)
#else /* not FORCE_NOMPI */
#  define CHECK_MPI(call) do {                                          \
    int error_code;                                                     \
    error_code = call;                                                  \
    if (error_code != MPI_SUCCESS) {                                    \
                                                                        \
      char error_string[MPI_MAX_ERROR_STRING];                          \
      int length_of_error_string, error_class;                          \
      int my_rank = -1;                                                 \
                                                                        \
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);                          \
      MPI_Error_class(error_code, &error_class);                        \
      MPI_Error_string(error_class, error_string,                       \
                       &length_of_error_string);                        \
      fprintf(stderr, "%3d: %s\n", my_rank, error_string);              \
      MPI_Error_string(error_code, error_string,                        \
                       &length_of_error_string);                        \
      fprintf(stderr, "%3d: %s\n", my_rank, error_string);              \
      MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                               \
    }                                                                   \
  } while(0)

#  ifdef FORCE_NOSMP
#  define CHECK_THREAD_LEVEL(THREAD_MODE) do {  \
    THREAD_MODE = THREAD_MODE;                  \
  } while (0)
#  else  /* not FORCE_NOSMP */
#    define CHECK_THREAD_LEVEL(THREAD_MODE)                             \
  do {                                                                  \
    int      provided;                                                  \
                                                                        \
    CHECK_MPI(MPI_Query_thread(&provided));                             \
    if (THREAD_MODE == API_THREAD_FUNNELED)                              \
      {                                                                 \
        switch(provided) {                                              \
        case MPI_THREAD_SINGLE:                                         \
          errorPrint("This run only supports MPI_THREAD_SINGLE\n"       \
                     "  either use -DFORCE_NOSMP,\n"                    \
                     "  change your MPI Library\n"                      \
                     "  or check that MPI_Init_thread"                  \
                     " is correctly called\n");                         \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        case MPI_THREAD_FUNNELED:                                       \
        case MPI_THREAD_SERIALIZED:                                     \
        case MPI_THREAD_MULTIPLE:                                       \
          break;                                                        \
        default:                                                        \
          errorPrint("provided thread level support is unknown");       \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        }                                                               \
      }                                                                 \
    else                                                                \
      {                                                                 \
        switch(provided) {                                              \
        case MPI_THREAD_SINGLE:                                         \
          errorPrint("This run only supports MPI_THREAD_SINGLE\n"       \
                     "  either use -DFORCE_NOSMP,\n"                    \
                     "  change your MPI Library\n"                      \
                     "  or check that MPI_Init_thread"                  \
                     " is correctly called\n");                         \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        case MPI_THREAD_FUNNELED:                                       \
          errorPrint("This run only supports MPI_THREAD_FUNNELED\n"     \
                     "  either use API_THREAD_FUNNELED,\n"              \
                     "  change your MPI Library\n"                      \
                     "  or check that MPI_Init_thread"                  \
                     " is correctly called\n");                         \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        case MPI_THREAD_SERIALIZED:                                     \
          errorPrint("This run only supports MPI_THREAD_SERIALIZED\n"   \
                     "  either use API_THREAD_FUNNELED,\n"                \
                     "  change your MPI Library\n"                      \
                     "  or check that MPI_Init_thread"                  \
                     " is correctly called\n");                         \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        case MPI_THREAD_MULTIPLE:                                       \
          break;                                                        \
        default:                                                        \
          errorPrint("provided thread level support is unknown");       \
          MPI_Abort(MPI_COMM_WORLD, MPI_ERR);                           \
          break;                                                        \
        }                                                               \
      }                                                                 \
  } while(0)
#  endif /* not FORCE_NOSMP */
#endif /* not FORCE_NOMPI */

#endif /* COMMON_PASTIX_H */
