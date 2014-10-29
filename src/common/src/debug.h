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
  File: debug.h

  Defines debugs flags and <print_debug> macro which use them.

  Authors:
    Mathieu Faverge - faverge@labri.fr
    Xavier  LACOSTE - lacoste@labri.fr

 */
/* DEBUG FLAGS for print_debug */
/* Scotch */
#define DBG_SCOTCH              0
#define DBG_STEP                0
#define DBG_CSCD                0
/* blend */
#define DBG_BUBBLESPLIT         0
#define DBG_BUBBLES             0

/* Sopalin */
#define DBG_SOPALIN_NAN         0
#define DBG_SOPALIN_INF         0
#define DBG_SOPALIN_RAFF        0
#define DBG_SOPALIN_DEBUG       0
#define DBG_SOPALIN_MAIN        0
#define DBG_SOPALIN_THREADCOMM  0
#define DBG_SOPALIN_ALLOC       0
#define DBG_SOPALIN_BLEND       0
#define DBG_SOPALIN_DRUNK       0
#define DBG_SOPALIN_COMM        0
#define DBG_SOPALIN_COMPUTE     0
#define DBG_SOPALIN_COMP1D      0
#define DBG_SOPALIN_DIAG        0
#define DBG_SOPALIN_E1          0
#define DBG_SOPALIN_E2          0
#define DBG_SOPALIN_NAPA        0
#define DBG_FUNNELED            0
#define DBG_THCOMM              0
#define DBG_UPDO                0
#define DBG_PASTIX_DYNSCHED     0
#define DBG_SOPALIN_TIME        0
#define DBG_CSC_LOG             0
#define DBG_PASTIX_REVERTSTEAL  0

/* Sopalin SendRecv */
#define DBG_SOPALIN_SEND        0
#define DBG_SOPALIN_RECV        0


/* updown */		        
#define DBG_SOPALIN_UPDO        0
#define DBG_SOPALIN_UP          0
#define DBG_SOPALIN_DOWN        0

/* OOC */		        
#define DBG_OOC_TRACE_V1        0
#define DBG_OOC_TRACE_V2        0
#define DBG_OOC_DEBUG           0
#define DBG_OOC_PREDICTION      0
#define DBG_OOC_SAVE            0
#define DBG_OOC_MUTEX_ALLOCATED 0
#define DBG_OOC_WAIT_FOR_FTGT   0
#define DBG_OOC_FTGT            0

/* Raff */
#define DBG_RAFF_PIVOT          0
#define DBG_RAFF_GMRES          0
#define DBG_RAFF_GRAD           0

/* MURGE */
#define DBG_MURGE               0

/*
 * Macro: print_debug
 *
 * Prints debugging message if PASTIX_DEBUG is defined and if the
 * debug flag is the to 1.
 */
#ifdef PASTIX_DEBUG
#define print_debug(mod,...) {if (mod) fprintf(stderr, __VA_ARGS__);}
#else
#define print_debug(...)     {}
#endif

/*
 * Macros:
 *   FLOAT_FMT   - Format for printing a float.
 *   FLOAT_PRINT - Cast *val* to be usable with <FLOAT_FMT>.
 *   PRINT_COEF  - Print "FILE:LINE <variable name> <value>"
 *   PRINT_REAL  - Print "FILE:LINE <variable name> <value>", with a real variable
 */
#ifdef TYPE_COMPLEX
#  define FLOAT_FMT "%.20g %.20g"
#  ifdef PREC_DOUBLE
#    define FLOAT_PRINT(val) creal(val), cimag(val)
#  else
#    define FLOAT_PRINT(val) (double)crealf(val), (double)cimagf(val)
#  endif
#else
#  define FLOAT_FMT "%.20g"
#  ifdef PREC_DOUBLE
#    define FLOAT_PRINT(val) (val)
#  else
#    define FLOAT_PRINT(val) (double)(val)
#  endif
#endif
#define PRINT_COEF(val)                                                 \
  do {                                                                  \
    fprintf(stderr, "%s:%d " #val " " FLOAT_FMT "\n",                   \
            __FILE__, __LINE__, FLOAT_PRINT(val));                      \
  }while(0)

#define PRINT_REAL(val)                                                 \
  do {                                                                  \
    fprintf(stderr, "%s:%d " #val " %.20g \n",                          \
            __FILE__, __LINE__, (double)(val));                         \
  }while(0)
