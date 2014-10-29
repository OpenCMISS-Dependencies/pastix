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
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <pthread.h>

#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#include "common_pastix.h"
#include "queue.h"
#include "bulles.h"
#include "sopalin_thread.h"

#if (defined X_ARCHalpha_compaq_osf1)
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/processor.h>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#define X_INCLUDE_CXML
#endif

#ifdef X_ARCHi686_pc_linux
#include <sched.h>
#ifdef X_ARCHi686_mac
#include <mach/thread_act.h>
#include <mach/mach_init.h>
#endif
#endif


#include "common_pastix.h"
#include "tools.h"
#ifdef PASTIX_EZTRACE
#include "pastix_eztrace.h"
#else
#include "trace.h"
#endif
#include "sopalin_define.h"
#include "symbol.h"
#include "ftgt.h"
#include "csc.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
#include "sopalin_thread.h"
#include "stack.h"
#include "sopalin3d.h"

#ifdef WITH_HWLOC
#include <hwloc.h>
#endif

#define print_one(fmt, ...) if(procnum == 0) fprintf(stdout, fmt, __VA_ARGS__)

/******************************************************************************/
#ifdef FORCE_NOSMP

/********************************************/
/*     Lancement du calcul non-threadé      */
/********************************************/
/*
 *  Function: sopalin_launch_thread:
 *
 *  Launch non-threaded computation.
 *
 *  Defined if *FORCE_NOSMP* is defined.
 *
 *  Parameters:
 *    procnum      - MPI processus number.
 *    procnbr      - Number of MPI processus.
 *    ptr          - .
 *    calc_thrdnbr -
 */

void sopalin_launch_thread(void * sopalin_data_ref,
                           PASTIX_INT procnum, PASTIX_INT procnbr, void *ptr, PASTIX_INT verbose,
                           PASTIX_INT calc_thrdnbr, void * (*calc_routine)(void *), void *calc_data,
                           PASTIX_INT comm_thrdnbr, void * (*comm_routine)(void *), void *comm_data,
                           PASTIX_INT ooc_thrdnbr,  void * (*ooc_routine)(void *),  void *ooc_data){
  sopthread_data_t  d_comm;
  sopthread_data_t  d_ooc;
  sopthread_data_t  d_calc;
  pthread_t         pthread_comm;
  pthread_t         pthread_ooc;
  pthread_attr_t    attr_comm;
  pthread_attr_t    attr_ooc;
  int               ret;
  Sopalin_Data_t   *sopalin_data = sopalin_data_ref;
  (void)procnbr; (void)ptr;

  /* Lancement d'un thread de chargement ooc si il est séparé */
  if (verbose > API_VERBOSE_NO)
    print_one("Launching %d threads"
              " (%d commputation, %d communication, %d out-of-core)",
              (int)(comm_thrdnbr+ooc_thrdnbr),
              0, (int)comm_thrdnbr, (int)ooc_thrdnbr);

  if (ooc_thrdnbr > 0)
    {
      pthread_attr_init(&attr_ooc);
      d_ooc.me   = 2;
      d_ooc.data = ooc_data;

      ret = pthread_create(&pthread_ooc, &attr_ooc, ooc_routine, (void *)&d_ooc);
      if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    }

  /* Lancement d'un thread de comm si il est séparé */
  if (comm_thrdnbr > 0)
    {

      pthread_attr_init(&attr_comm);
      d_comm.me   = 1;
      d_comm.data = comm_data;

      ret = pthread_create(&pthread_comm, &attr_comm, comm_routine, (void *)&d_comm);
      if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    }

  d_calc.me   = 0;
  d_calc.data = calc_data;

  calc_routine(&d_calc);

  /* Récupération du thread de comm */
  if (comm_thrdnbr > 0)
    {
      ret = pthread_join(pthread_comm,(void**)NULL);
      if (ret) {errorPrint("thread join."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    }

  /* Récupération du thread de chargement ooc */
  if (ooc_thrdnbr > 0)
    {
      ret = pthread_join(pthread_ooc,(void**)NULL);
      if (ret) {errorPrint("thread join."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    }
}

#else
/********************************************/
/*       Lancement des threads POSIX        */
/********************************************/
#ifndef PASTIX_BUBBLESCHED

void sopalin_launch_thread(void * sopalin_data_ref,
                           PASTIX_INT procnum, PASTIX_INT procnbr, void *ptr, PASTIX_INT verbose,
                           PASTIX_INT calc_thrdnbr, void * (*calc_routine)(void *), void *calc_data,
                           PASTIX_INT comm_thrdnbr, void * (*comm_routine)(void *), void *comm_data,
                           PASTIX_INT ooc_thrdnbr,  void * (*ooc_routine)(void *),  void *ooc_data)
{
  sopthread_data_t *d       = NULL;
  pthread_t        *calltab = NULL;
  PASTIX_INT               i;
  PASTIX_INT               ret;
  PASTIX_INT               thrdnbr;
  PASTIX_INT               thrdnbr_wo_ooc;
  Sopalin_Data_t   *sopalin_data = sopalin_data_ref;
  (void)procnbr; (void)ptr;

  thrdnbr = calc_thrdnbr + comm_thrdnbr + ooc_thrdnbr ;
  thrdnbr_wo_ooc = calc_thrdnbr + comm_thrdnbr;
  if (verbose > API_VERBOSE_NO)
    print_one("Launching %d threads"
              " (%d commputation, %d communication, %d out-of-core)\n",
              (int) thrdnbr, (int)calc_thrdnbr,
              (int)comm_thrdnbr, (int)ooc_thrdnbr);
  MALLOC_INTERN(calltab, thrdnbr, pthread_t);
  MALLOC_INTERN(d,       thrdnbr, sopthread_data_t);

  if (calc_thrdnbr > 1)
    {
      int comm_size;
      CHECK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &comm_size));
      if (comm_size > 1)
        CHECK_THREAD_LEVEL(sopalin_data->sopar->iparm[IPARM_THREAD_COMM_MODE]);
    }
  /* Lancement des threads de calcul */
  for (i=0;i<calc_thrdnbr;i++)
    {
      pthread_attr_t attr;
      pthread_attr_init(&attr);

#ifdef MARCEL2
      {
        int cpu = (i+procnum*thrdnbr)%sysconf(_SC_NPROCESSORS_ONLN);
        if (thrdnbr <= sysconf(_SC_NPROCESSORS_ONLN))
          marcel_attr_setvpset(&attr, MARCEL_VPSET_VP(cpu));
      }
#endif

      d[i].me   = i;
      d[i].data = calc_data;
      ret = pthread_create(&calltab[i],&attr,calc_routine,(void *)&d[i]);
      if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    }

  /* Lancement des threads de chargement ooc */
  for (i=thrdnbr_wo_ooc; i<thrdnbr; i++)
    {
      pthread_attr_t attr;
      pthread_attr_init(&attr);
#ifdef MARCEL2
      {
        int cpu = (i+procnum*thrdnbr)%sysconf(_SC_NPROCESSORS_ONLN);

        if (thrdnbr <= sysconf(_SC_NPROCESSORS_ONLN))
          marcel_attr_setvpset(&attr, MARCEL_VPSET_VP(cpu));
      }
#endif

      d[i].me   = i;
      d[i].data = ooc_data;
      ret = pthread_create(&calltab[i],&attr,ooc_routine,(void *)&d[i]);
      if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    }

  /* Lancement des threads de communication */
  if ((comm_thrdnbr > 0)  && (comm_routine != NULL))
    {
      /*       print_one("-- Options Communication --\n"); */
      /*       print_one("    - Type            : %d\n", sopar->type_comm); */
      /*       print_one("    - Nbthread        : %d\n", sopar->nbthrdcomm);   */

      if (comm_thrdnbr > 1)
        {
          for (i=calc_thrdnbr;i<thrdnbr_wo_ooc;i++)
            {
              pthread_attr_t attr;
              pthread_attr_init(&attr);

              d[i].me   = i;
              d[i].data = comm_data;
              ret = pthread_create(&calltab[i], &attr,
                                   comm_routine, (void *)&d[i]);
              if (ret)
                {
                  errorPrint("thread create.");
                  EXIT(MOD_SOPALIN,THREAD_ERR);
                }
            }
        }
      else
        {
          d[calc_thrdnbr].me   = calc_thrdnbr;
          d[calc_thrdnbr].data = comm_data;
          comm_routine((void*)&d[calc_thrdnbr]);
        }
    }

  /* Recuperation de tous les threads lancés */
  for (i=0;i<thrdnbr;i++)
    {
      /* On ne recupere pas le thread qd il a pas été lancé */
      if ((comm_thrdnbr == 1) && (i == calc_thrdnbr)) continue;

      ret = pthread_join(calltab[i],(void**)NULL);
      if (ret) {errorPrint("thread join."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    }

  memFree_null(calltab);
  memFree_null(d);
}

/**********************************************/
/* Lancement des threads MARCEL (avec bulles) */
/**********************************************/
#else
void sopalin_launch_thread(void * sopalin_data_ref,
                           PASTIX_INT procnum, PASTIX_INT procnbr, void *ptr, PASTIX_INT verbose,
                           PASTIX_INT calc_thrdnbr, void * (*calc_routine)(void *), void *calc_data,
                           PASTIX_INT comm_thrdnbr, void * (*comm_routine)(void *), void *comm_data,
                           PASTIX_INT ooc_thrdnbr,  void * (*ooc_routine)(void *),  void *ooc_data){
  char             *name;
  pthread_t        *calltab;
  marcel_bubble_t  *bubbletab;
  BubbleTree       *btree;
  sopthread_data_t *d;
  int               bubbnbr;       /* Nombre de bulles crées */
  int               thrdnbr;       /* Nombre total de thread à lancer */
  int               nodenbr;       /* Nombre de noeud dans l'arbre de bulles */
  int               nbbub2alloc;
  int               ret;
  int               me = 0;
  int               i;
  Sopalin_Data_t   *sopalin_data = sopalin_data_ref;
  (void)procnbr;

  btree        = (BubbleTree *)ptr;
  nodenbr      = btree->nodenbr;
  calc_thrdnbr = btree->leavesnbr;
  bubbnbr      = nodenbr - calc_thrdnbr;
  thrdnbr      = calc_thrdnbr + comm_thrdnbr;
  nbbub2alloc  = MAX(bubbnbr, 1);
  if (verbose > API_VERBOSE_NO)
    print_one("Launching %d threads"
              " (%d commputation, %d communication, %d out-of-core)",
              thrdnbr, calc_thrdnbr, comm_thrdnbr, 0);

  if (comm_thrdnbr > 0)
    nbbub2alloc++;

  MALLOC_INTERN(name,      30,          char);
  MALLOC_INTERN(calltab,   thrdnbr,     pthread_t);
  MALLOC_INTERN(d,         thrdnbr,     sopthread_data_t);
  MALLOC_INTERN(bubbletab, nbbub2alloc, marcel_bubble_t);

#ifdef PROFILE
  marcel_start_playing();
#endif

  /* Lancement dans le cas ou on a un arbre de bulles */
  if (bubbnbr > 0)
    {
      for (i=0; i<bubbnbr; i++)
  marcel_bubble_init(&bubbletab[i]);

      for (i=1; i<bubbnbr; i++)
  marcel_bubble_insertbubble(&bubbletab[ btree->nodetab[i+calc_thrdnbr].fathnum - calc_thrdnbr ],
           &bubbletab[i]);

      /* Creation des threads feuilles */
      for(i=0; i<calc_thrdnbr; i++)
  {
    marcel_attr_t attr;
    sprintf(name, "thread%03d", i);
    marcel_attr_init(&attr);
    marcel_attr_setinitbubble(&attr, &bubbletab[btree->nodetab[i].fathnum - calc_thrdnbr]);
    marcel_attr_setid(&attr,0);
    marcel_attr_setprio(&attr,MA_MAX_RT_PRIO);
    marcel_attr_setname(&attr,name);

    d[i].me   = i;
    d[i].data = calc_data;

    ret = pthread_create(&calltab[i], &attr, calc_routine, (void *)&d[i]);
    if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
  }

      /* Creation d'un thread par bulle */
/*       for (i=calc_thrdnbr; i<nodenbr; i++) */
/*  { */
/*    marcel_attr_t attr;  */
/*    int prio = MA_MAX_RT_PRIO+btree->nodetab[i].treelevel; /\* + [btsknbr[i][1]-btsknbr[i][0]+1);*\/ */
/*    sprintf(name, "threadBubble%03d", i); */
/*    if (prio > MA_DEF_PRIO) */
/*      prio = MA_DEF_PRIO; */
/*    if (prio < MA_MAX_RT_PRIO) */
/*      prio = MA_MAX_RT_PRIO; */
/*    marcel_attr_init(&attr); */
/*    marcel_attr_setinitbubble(&attr, &bubbletab[i-thrdnbr]); */
/*    marcel_attr_setid(&attr,1); */
/*    marcel_attr_setprio(&attr,prio); */
/*    marcel_attr_setname(&attr,name); */

/*    d[i].me   = i; */
/*    d[i].data = data; */

/*    ret = pthread_create(&calltab[i], &attr, calc_routine, (void *)&d[i]); */
/*    if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);} */
/*  } */

    }
  /* Cas ou on a qu'un seul thread par noeud */
  else
    {
      marcel_attr_t attr;

      marcel_bubble_init(&bubbletab[0]);

      marcel_attr_init(&attr);
      marcel_attr_setinitbubble(&attr, &bubbletab[0]);
      marcel_attr_setid(&attr,0);
      marcel_attr_setprio(&attr,MA_MAX_RT_PRIO);
      marcel_attr_setname(&attr,"thread");

      d[0].me   = 0;
      d[0].data = calc_data;

      ret = pthread_create(&calltab[0], &attr, calc_routine, (void *)&d[0]);
      if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    }

#ifdef BUBBLE_SCHED_NULL
  print_one("ODONNANCEUR : sans\n");
#elif (defined BUBBLE_SCHED_SPREAD)
  print_one("ODONNANCEUR : bubble_spread\n");
#elif (defined BUBBLE_SCHED_AFFINITY)
  print_one("ODONNANCEUR : bubble_affinity\n");
#elif (defined BUBBLE_SCHED_MEMAWARE)
  print_one("ODONNANCEUR : bubble_memaware\n");
#else
  print_one("ODONNANCEUR : inconnu\n");
#endif

  /* Lancement des threads de communication */
  if ((comm_thrdnbr > 0)  && (comm_routine != NULL))
    {
      if (comm_thrdnbr > 1)
  {
    marcel_bubble_init(&bubbletab[bubbnbr]);
    marcel_bubble_insertbubble(&bubbletab[0], &bubbletab[bubbnbr]);

    for (i=calc_thrdnbr;i<thrdnbr;i++)
      {
        marcel_attr_t attr;
        sprintf(name, "threadComm%03d", i);
        marcel_attr_init(&attr);
        marcel_attr_setinitbubble(&attr, &bubbletab[bubbnbr]);
        marcel_attr_setid(&attr,0);
        marcel_attr_setprio(&attr,MA_MAX_RT_PRIO);
        marcel_attr_setname(&attr,name);

        d[i].me   = i;
        d[i].data = comm_data;
        ret = pthread_create(&calltab[i], &attr, comm_routine, (void *)&d[i]);
        if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
      }
  }
      else
  {
    /* Reveil de la bulle principale */
    marcel_wake_up_bubble(&bubbletab[0]);
    d[calc_thrdnbr].me   = calc_thrdnbr;
    d[calc_thrdnbr].data = comm_data;
    thrdnbr--;
    comm_routine((void*)&d[calc_thrdnbr]);
  }
    }
  else
    {
      comm_thrdnbr = 0;
      thrdnbr = calc_thrdnbr;
    }

  /* Reveil de la bulle principale */
  marcel_wake_up_bubble(&bubbletab[0]);

  for (i=0;i<thrdnbr;i++)
    {
      ret = pthread_join(calltab[i],(void**)NULL);
      if (ret)
  {
    errorPrint("thread join bubbnum %d, ret %d", i, ret);
    perror("erreur join thread");
    EXIT(MOD_SOPALIN,THREAD_ERR);
  }
    }

  for(i=0; i<bubbnbr; i++)
    {
      marcel_bubble_destroy(&bubbletab[i]);
    }

  memFree_null(name);
  memFree_null(bubbletab);
  memFree_null(calltab);
  memFree_null(d);
}
#endif
#endif

/******************************************************************************/

/**********************************************/
/*        Bind des threads sur les procs      */
/**********************************************/
PASTIX_INT sopalin_bindthread(PASTIX_INT cpu)
{
#ifdef MARCEL

  {
    marcel_vpset_t vpset = MARCEL_VPSET_ZERO;
    marcel_vpset_vp(&vpset, cpu);
    marcel_apply_vpset(&vpset);
  }

#else /* Dans les autres cas on se preoccupe de l'archi */

#ifdef WITH_HWLOC
  {
    hwloc_topology_t topology; /* Topology object */
    hwloc_obj_t      obj;      /* Hwloc object    */
    hwloc_cpuset_t   cpuset;   /* HwLoc cpuset    */

    /* Allocate and initialize topology object.  */
    hwloc_topology_init(&topology);

    /* Perform the topology detection.  */
    hwloc_topology_load(topology);

    /* Get last one.  */
    obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, cpu);
    if (!obj)
      return 0;

    /* Get a copy of its cpuset that we may modify.  */
    /* Get only one logical processor (in case the core is SMT/hyperthreaded).  */
#if !defined(HWLOC_BITMAP_H)
    cpuset = hwloc_cpuset_dup(obj->cpuset);
    hwloc_cpuset_singlify(cpuset);
#else
    cpuset = hwloc_bitmap_dup(obj->cpuset);
    hwloc_bitmap_singlify(cpuset);
#endif

    /* And try to bind ourself there.  */
    if (hwloc_set_cpubind(topology, cpuset, HWLOC_CPUBIND_THREAD)) {
      char *str = NULL;
#if !defined(HWLOC_BITMAP_H)
      hwloc_cpuset_asprintf(&str, obj->cpuset);
#else
      hwloc_bitmap_asprintf(&str, obj->cpuset);
#endif
      printf("Couldn't bind to cpuset %s\n", str);
      free(str);
    }

    /* Get the number at Proc level */
    cpu = obj->children[0]->os_index;

    /* Free our cpuset copy */
#if !defined(HWLOC_BITMAP_H)
    hwloc_cpuset_free(cpuset);
#else
    hwloc_bitmap_free(cpuset);
#endif

    /* Destroy topology object.  */
    hwloc_topology_destroy(topology);
  }
#else /* WITH_HWLOC */
#ifdef X_ARCHpower_ibm_aix
  {
    tid_t self_ktid = thread_self ();

    bindprocessor(BINDTHREAD, self_ktid, cpu);
  }
#elif (defined X_ARCHalpha_compaq_osf1)
  {
    bind_to_cpu_id(getpid(), cpu, 0);
  }
#elif (defined X_ARCHi686_pc_linux)

#ifndef X_ARCHi686_mac
  {
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(cpu, &mask);

#ifdef HAVE_OLD_SCHED_SETAFFINITY
    if(sched_setaffinity(0,&mask) < 0)
#else /* HAVE_OLD_SCHED_SETAFFINITY */
    if(sched_setaffinity(0,sizeof(mask),&mask) < 0)
#endif /* HAVE_OLD_SCHED_SETAFFINITY */
      {
  perror("sched_setaffinity");
  EXIT(MOD_SOPALIN, INTERNAL_ERR);
      }
  }
#else /* X_ARCHi686_mac */
  {
    thread_affinity_policy_data_t ap;
    int                           ret;

    ap.affinity_tag = 1; /* non-null affinity tag */
    ret = thread_policy_set(
          mach_thread_self(),
          THREAD_AFFINITY_POLICY,
          (integer_t*) &ap,
          THREAD_AFFINITY_POLICY_COUNT
          );
    if(ret != 0)
      {
  perror("thread_policy_set");
  EXIT(MOD_SOPALIN, INTERNAL_ERR);
      }
  }
#endif /* X_ARCHi686_mac */
#endif /* X_ACHIxxx      */
#endif /* WITH_HWLOC     */
#endif /* MARCEL         */

  return cpu;
}

/******************************************************************************/

/**********************************************/
/*   Lancement des threads de communication   */
/**********************************************/
void sopalin_launch_comm(int nbthrdcomm, void * (*comm_routine)(void *), void *data)
{
  int ret, i;
  pthread_t        *calltab = NULL;
  sopthread_data_t *d       = NULL;

  MALLOC_INTERN(calltab, nbthrdcomm, pthread_t);
  MALLOC_INTERN(d,       nbthrdcomm, sopthread_data_t);

  for (i=0;i<nbthrdcomm;i++)
    {

      pthread_attr_t attr;
      pthread_attr_init(&attr);

      d[i].me   = -1-i;
      d[i].data = data;

      ret = pthread_create(&calltab[i],&attr,comm_routine,(void *)&d[i]);
      if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}

    }

  for (i=0; i<nbthrdcomm; i++)
    pthread_detach(calltab[i]);

  memFree_null(calltab);
  /* memFree_null(d);*/

}
