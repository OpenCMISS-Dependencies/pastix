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
/*+********************************************************+*/
/*+                                                        +*/
/*+   NAME       : simu.h                                  +*/
/*+                                                        +*/
/*+   AUTHORS    : Pascal HENON                            +*/
/*+                                                        +*/
/*+   FUNCTION   : Part of a parallel direct block solver. +*/
/*+                Structure used in the distribution      +*/
/*+                lead by simulation                      +*/
/*+   DATES      : # Version 0.0  : from : 28 sep 1998     +*/
/*+                                 to     05 oct 1998     +*/
/*+                                                        +*/
/*+********************************************************+*/


/*
**  The type and structure definitions.
*/
typedef struct SimuTimer_ {
  double s;
  double ms;
} SimuTimer;

typedef struct SimuCluster_ {
  PASTIX_INT                fprocnum;
  PASTIX_INT                lprocnum;
  ExtendVectorINT  * ftgtsend;         /*+ ftgt send by this proc (one vector per processor       +*/
  PASTIX_INT prionum;
} SimuCluster;

typedef struct SimuProc_ {
  SimuTimer          timer;            /*+ Simulated clock of the processor                       +*/
  Queue            * taskheap;         /*+ heap of cblk that all contrib have been receveive      +*/
  Queue            * taskheap2;        /*+ queue of cblk ordered by treelevel                     +*/

  PASTIX_INT                prionum;          /*+ Current priority to assign to a cblk mapp on this proc +*/
  ExtendVectorINT  * tasktab;
} SimuProc;

typedef struct SimuFtgt_ {
  FanInTarget           ftgt;             /*+ the ftgt structures info                              +*/
  PASTIX_INT                   clustnum;         /*+ the cluster that send the ftgt                        +*/
  SimuTimer             timerecv;         /*+ time the ftgt will be receive                         +*/
  double                costsend;         /*+ time for the ftgt go from procsrc to procdest         +*/
  double                costadd;          /*+ cost to add the ftgt to its cblk destination          +*/
} SimuFtgt;        

typedef struct SimuBlockTarget_ {
  PASTIX_INT           prionum;
  PASTIX_INT           taskid;
  PASTIX_INT           bloksrc;
  PASTIX_INT           cblkdest;
  PASTIX_INT           frownum;
  PASTIX_INT           lrownum;
} SimuBlockTarget;

typedef struct SimuBlok_ {
  PASTIX_INT                   tasknum;          /*+ Number of the task                                     +*/
  PASTIX_INT                   ftgtnum;          /*+ index of the first ftgt destinated to this cblk 
					    in the ftgttab. This index is also used to find
					    the first cblk timer (one per cand proc) in the
					    timetab                                                +*/
  PASTIX_INT                   ctrbcnt;          /*+ counter for contributions remaining                    +*/
  PASTIX_INT                   ctrbnbr;          /*+ OIMBE temporaire sert juste pour DRUNK                 +*/
    
} SimuBlok;


typedef struct SimuTask_ {
  PASTIX_INT                   taskid;           /*+ identification of the task type                        +*/
  PASTIX_INT                   prionum;          /*+ priority of the task                                   +*/
  PASTIX_INT                   cblknum;          /*+ Number of the cblknum the task deals with              +*/
  PASTIX_INT                   bloknum;          /*+ number of the block that the task deals with           +*/
  PASTIX_INT                   bloknum2;         
  PASTIX_INT                   facebloknum;      /*+ Number of the facing block for E2                      +*/  
  SimuTimer             time;             /*+ Time the task is ready if it doesn't need message      +*/   
  PASTIX_INT                   mesglen;         /*+ Time to send the block target                          +*/
  double                cost;             /*+ Cost of the task                                       +*/
  PASTIX_INT                   ctrbcnt;          /* nbr ftgt + le btgt (+ E1 pret si E2) */
  PASTIX_INT                   ftgtcnt;          /* nbr of contrib from fan-in target                       +*/
  PASTIX_INT                   tasknext;         /* chainage des E1 ou E2, si fin = -1 => liberer les btagptr */
} SimuTask;

/* task type allowed*/
#ifndef COMP_1D
#define COMP_1D    0
#define DIAG       1
#define E1         2
#define E2         3
#endif



typedef struct SimuCblk_ {
  PASTIX_INT                   ctrbcnt;          /*+ counter for contributions remaining                    +*/
} SimuCblk;


typedef struct SimuCtrl_ {
  PASTIX_INT                   cblknbr;          /*+ Number of cblk                                          +*/
  PASTIX_INT                   ftgtprio;         /*+ Priority to assign to current ftgts                     +*/
  PASTIX_INT                   tasknbr;          /*+ Number of tasks                                         +*/
  PASTIX_INT                   ftgtcnt;          /*+ Number of received communication                        +*/
  SimuTask         *    tasktab;          /*+ SimuTask vector                                         +*/
  SimuProc         *    proctab;          /*+ Virtual processor tab                                   +*/
  SimuCluster      *    clustab;          /*+ Virtual cluster tab                                     +*/
  PASTIX_INT              *    ownetab;          /*+ Vector containing the distribution of the diagonal blok +*/
  PASTIX_INT              *    blprtab;          /*+ Vector containing the distribution of the blok          +*/
  SimuCblk         *    cblktab;          /*+ SimuCblk vector                                         +*/
  SimuBlok         *    bloktab;          /*+ SimuBlok vector                                         +*/
  SimuFtgt         *    ftgttab;          /*+ Vector containing the fan in target                     +*/
  PASTIX_INT                   ftgtnbr;
  SimuTimer        *    tasktimetab;      /*+ Vector containing a timer for each cand proc on each task +*/
  SimuTimer        *    ftgttimetab;      /*+ Vector containing a timer for each cluster on each ftgt  +*/
} SimuCtrl;

/*
**  The function prototypes.
*/

#ifndef SIMU
#define static 
#endif

PASTIX_INT                 simuInit        (SimuCtrl *, SymbolMatrix *, PASTIX_INT, PASTIX_INT, PASTIX_INT, PASTIX_INT, Cand *);
PASTIX_INT                 simuRealloc     (SimuCtrl *, PASTIX_INT, PASTIX_INT);
void                simuExit        (SimuCtrl *, PASTIX_INT, PASTIX_INT, PASTIX_INT);
PASTIX_INT                 compTimer       (SimuTimer *, SimuTimer *);
void                timerAdd        (SimuTimer *, double);
double              timerVal        (SimuTimer *);
void                timerSet        (SimuTimer *, double);
#undef static



#define CLUST2INDEX(n,c) ((c) + simuctrl->bloktab[n].ftgtnum - ctrl->candtab[ctrl->egraph->ownetab[n]].fccandnum)
#define INDEX2CLUST(r,s) ((r) - simuctrl->bloktab[s].ftgtnum + ctrl->candtab[ctrl->egraph->ownetab[s]].fccandnum)
#define TIMER(pr)        (&(simuctrl->proctab[pr].timer))


