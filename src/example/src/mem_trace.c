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
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <pthread.h>
#ifndef FORCE_NOMPI
#include <mpi.h>
#endif

#include "mem_trace.h"
#ifdef __APPLE__
#include<mach/mach.h>
#endif

int trace_task = 0;

typedef struct {
  long size,resident,share,text,lib,data,dt;
} statm_t;

static
void read_off_memory_status(statm_t * result)
{
#ifdef __APPLE__

  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if ( KERN_SUCCESS != task_info( mach_task_self(),
                                  TASK_BASIC_INFO,
                                  (task_info_t)&t_info,
                                  &t_info_count))
  (*result).size     = t_info.virtual_size;
  (*result).resident = t_info.resident_size;
  (*result).share    = -1;
  (*result).text     = -1;
  (*result).lib      = -1;
  (*result).data     = -1;
  (*result).dt       = -1;
#else
  const char* statm_path = "/proc/self/statm";

  FILE *f = fopen(statm_path,"r");
  if(!f){
    perror(statm_path);
    abort();
  }
  if(7 != fscanf(f,"%ld %ld %ld %ld %ld %ld %ld",
                 &((*result).size),
                 &((*result).resident),
                 &((*result).share),
                 &((*result).text),
                 &((*result).lib),
                 &((*result).data),
                 &((*result).dt)))
  {
    perror(statm_path);
    abort();
  }
  fclose(f);
#endif /* __APPLE__ */
}

static int         flag_stop = 0;
static pthread_t   mem_thread;

static
void * follow_memory_usage(void * arg)
{
  int rank;
  char filename[256];
  FILE * f;
  statm_t mem;
  struct timeval tv;
  double t1, t2;

#ifndef FORCE_NOMPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = 0;
#endif

  sprintf(filename, "memory_usage_%d", rank);
  f = fopen(filename, "w");
  gettimeofday(&tv, NULL);
  t1 = ((double) tv.tv_sec + (double) tv.tv_usec * (double)1.0e-6L);
  fprintf(f, "%s %s %s\n", "time", "memSize", "resident");

  while (flag_stop == 0) {
    read_off_memory_status(&mem);
    gettimeofday(&tv, NULL);
    t2 = ((double) tv.tv_sec + (double) tv.tv_usec * (double)1.0e-6L);
    fprintf(f, "%e %ld %ld %d\n", t2-t1 + 10, mem.size, mem.resident,
            trace_task);
    usleep(10000);
  }
  fclose(f);
  return NULL;
}


void memtrace_start(void)
{
  pthread_create(&mem_thread, NULL, follow_memory_usage, NULL);

}

void memtrace_stop(void)
{
  flag_stop = 1;
  pthread_join(mem_thread, NULL);
}

