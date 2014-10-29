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
/************************************************************/
/**                                                        **/
/**   NAME       : param_blend.h                           **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the blend control structure.        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ The parameters structure definition +*/

typedef struct netperf_ {
  double startup;
  double bandwidth;
} netperf;


/*+ Structure containing the structure passed through the blend primitives +*/
typedef struct BlendCtrl_ {
  BlendParam        *option;
  netperf           *perfptr;

  PASTIX_INT                clustnum;      /*+ Local cluster ID                       +*/
  PASTIX_INT                clustnbr;      /*+ Number of MPI process                  +*/
  PASTIX_INT                procnbr;       /*+ Number total of processors             +*/
  PASTIX_INT                proclocnbr;    /*+ Number of processors for one clustnum  +*/
  PASTIX_INT                thrdnbr;       /*+ Number total of threads                +*/
  PASTIX_INT                thrdlocnbr;    /*+ Number of threads for one clustnum     +*/
  PASTIX_INT                cudanbr;       /*+ Number of cuda device for one clustnum +*/
  PASTIX_INT                bublnbr;       /*+ Number of threads for one clustnum     +*/
  PASTIX_INT               *proc2clust;    /*+ proc2clust[i] = cluster of proc i      +*/
  BubbleTree        *btree;         /*+ arbre de bulles +*/
  EliminGraph       *egraph;        /*+ the elimination graph (only in vertex) +*/
  EliminTree        *etree;         /*+ the elimination tree                   +*/
  CostMatrix        *costmtx;       /*+ the cost bounded to each cblk and blok +*/
  Cand              *candtab;       /*+ processor candidate tab                +*/
  Queue             *lheap;         /*+ Use to order leaves                    +*/
  ExtendVectorINT   *intvec;        /*+ vector of PASTIX_INT used by several routines.
                                     The aim of this variable is to avoid
                                     repetedly memAlloc and memFree call      +*/
  ExtendVectorINT   *intvec2;       /*+ Another one                            +*/
  FILE              *tracefile;
} BlendCtrl;



/* Calcul le numero du process MPI */
#define CLUSTNUM(proc) (proc/(ctrl->procnbr/ctrl->clustnbr))
/* Calcul le noeud SMP sur lequel se trouve le process MPI si plusieurs par noeud SMP */
#define SMPNUM(clust) ( (clust*(ctrl->procnbr/ctrl->clustnbr)) / ctrl->option->procnbr)


PASTIX_INT      blendCtrlInit (BlendCtrl *, PASTIX_INT, PASTIX_INT, PASTIX_INT, PASTIX_INT, BlendParam *);
void     blendCtrlExit (BlendCtrl *);
void     perfcluster2  (PASTIX_INT procsrc, PASTIX_INT procdst, PASTIX_INT sync_comm_nbr, 
			netperf *np, BlendCtrl *ctrl);

