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
   File: sopalin_option.c

   Implements a function tha will print PaStiX
   compilation option.
*/
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#include "common_pastix.h"
#include "out.h"
#include "sopalin_define.h"

#define print_onempi(fmt, ...) if( SOLV_PROCNUM == 0 )           fprintf(stdout, fmt, ##__VA_ARGS__)
#define print_one(fmt, ...)    if( me == 0 && SOLV_PROCNUM == 0) fprintf(stdout, fmt, ##__VA_ARGS__)
#define print_all(fmt, ...)    fprintf(stdout, fmt, ##__VA_ARGS__)
#define print_error(...)

/*********************************/
/*
  Function: sopalin_option

  Print PaStiX compile options.

  Parameters:

  Returns:
    void
 */
/*********************************/
void sopalin_option(void){

  int smp    = 0;
  int mpi    = 0;
  int stats  = 0;
  int napa   = 0;
  int irecv  = 0;
  int isend  = 0;
  int conso  = 0;
  int fanbl  = 0;
  int typint = 0;
  int dbl    = 0;
  int cplx   = 0;
  int metis  = 0;
  int scotch = 0;
  int bubb   = 0;
  int ooc    = 0;
  int dist   = 0;

#if (defined EXACT_THREAD)
  char *tag = "Exact Thread";
#elif (defined EXACT_TAG)
  char *tag = "Exact Tag";
#else
  char *tag = "Tag Fanin or Block";
#endif

#ifdef STATS_SOPALIN
  stats = 1;
#endif
#ifdef NAPA_SOPALIN
  napa  = 1;
#endif
#ifdef SMP_SOPALIN
  smp   = 1;
#endif
#ifdef TEST_IRECV
  irecv = 1;
#endif
#ifdef TEST_ISEND
  isend = 1;
#endif
#ifdef FORCE_CONSO
  conso = 1;
#endif
#ifdef RECV_FANIN_OR_BLOCK
  fanbl = 1;
#endif
#ifndef FORCE_NOMPI
  mpi   = 1;
#endif
#ifdef FORCE_LONG
  typint= 1;
#elif (defined FORCE_INT32)
  typint= 2;
#elif (defined FORCE_INT64)
  typint= 3;
#endif
#ifdef PREC_DOUBLE
  dbl   = 1;
#endif
#ifdef TYPE_COMPLEX
  cplx  = 1;
#endif
#ifdef METIS
  metis = 1;
#endif
#ifdef WITH_SCOTCH
  scotch = 1;
#endif
#ifdef PASTIX_DYNSCHED
  bubb  =  1;
#endif
#ifdef OOC
  ooc   =  1;
#endif
#ifdef DISTRIBUTED
  dist = 1;
#endif

  fprintf(stdout, OUT_OPT_HEAD1);
  fprintf(stdout, OUT_OPT_HEAD2);
  fprintf(stdout, OUT_OPT_HEAD3);
  //fprintf(stdout, OUT_OPT_VERS,  VERSION);
  fprintf(stdout, OUT_OPT_SMP,   smp?  "Defined":"Not defined");
  fprintf(stdout, OUT_OPT_MPI,   mpi?  "Defined":"Not defined");
  fprintf(stdout, OUT_OPT_DSCD,  bubb? "Defined":"Not defined");
  fprintf(stdout, OUT_OPT_STATS, stats?"Defined":"Not defined");
  fprintf(stdout, OUT_OPT_NAPA,  napa? "Defined":"Not defined");
  fprintf(stdout, OUT_OPT_IRECV, irecv?"Defined":"Not defined");
  fprintf(stdout, OUT_OPT_ISEND, isend?"Defined":"Not defined");
  fprintf(stdout, OUT_OPT_TAG,   tag);
  fprintf(stdout, OUT_OPT_FORCE, conso?"Defined":"Not defined");
  fprintf(stdout, OUT_OPT_RFOB,  fanbl?"Defined":"Not defined");
  fprintf(stdout, OUT_OPT_OOC,   ooc?  "Defined":"Not defined");
  fprintf(stdout, OUT_OPT_DIST,  dist? "Defined":"Not defined");
  fprintf(stdout, OUT_OPT_METIS, metis?"Defined":"Not defined");
  fprintf(stdout, OUT_OPT_SCOTCH, scotch?"Defined":"Not defined");
  fprintf(stdout, OUT_OPT_INT,   ((typint%2) ? (( typint == 1 ) ? "long":"int64_t")
                                             : (( typint == 0 ) ? "int":"int32_t")));
  fprintf(stdout, OUT_OPT_FLOAT, dbl?"double":"simple", cplx?"complex":"");
  fprintf(stdout, OUT_OPT_END);

  return;
}
