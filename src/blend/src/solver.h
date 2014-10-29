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
/**   NAME       : solver.h                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the solver matrix.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     28 oct 1998     **/
/**                # Version 1.0  : from : 06 jun 2002     **/
/**                                 to     06 jun 2002     **/
/**                                                        **/
/************************************************************/

#ifndef SOLVER_H
#define SOLVER_H

/*
**  The type and structure definitions.
*/

#define COMP_1D                     0
#define DIAG                        1
#define E1                          2
#define E2                          3
#define DRUNK                       4

typedef struct Task_ {
  PASTIX_INT                       taskid;               /*+ COMP_1D DIAG E1 E2                                     +*/
  PASTIX_INT                       prionum;              /*+ Priority value for the factorization                   +*/
  PASTIX_INT                       prionum2;             /*+ Priority value for solve steps                         +*/
  PASTIX_INT                       cblknum;              /*+ Attached column block                                  +*/
  PASTIX_INT                       bloknum;              /*+ Attached block                                         +*/
  PASTIX_INT volatile              ftgtcnt;              /*+ Number of fan-in targets                               +*/
  PASTIX_INT volatile              ctrbcnt;              /*+ Total number of contributions                          +*/
  volatile BlockTarget *    btagptr;              /*+ si non local, pointeur sur la structure (NB reception) +*/
  PASTIX_INT                       indnum;               /*+ For E2 (COMP_1D), index of ftgt (>0) else if local = -taskdest
                                                      For DIAG and E1 , index of btag (>0) if there is a
						      local one it must be the first of the chain of local E1   +*/
  PASTIX_INT                       tasknext;             /*+ chainage des E1 ou E2, si fin = -1 => liberer les btagptr +*/
  PASTIX_INT                       taskmstr;             /*+ Master task for E1 or E2 tasks                         +*/
                                                  /*+ Index of DIAG (or first E1) task for E1                +*/
                                                  /*+ Index of E1 (or first E2) task for E2                  +*/
#if (defined PASTIX_DYNSCHED) || (defined TRACE_SOPALIN) || (defined BLEND_COMPUTE_THREADID)
  PASTIX_INT                       threadid;             /*+ Index of the bubble which contains the task +*/
  PASTIX_INT                       GPUid;
  PASTIX_INT                       cand;		  /*+ Thread candidate in static version          +*/
#endif
#ifdef TRACE_SOPALIN
  PASTIX_INT                       fcandnum;             /*+ First thread candidate                      +*/
  PASTIX_INT                       lcandnum;		  /*+ Last thread candidate                       +*/
  PASTIX_INT                       id;                   /*+ Global cblknum of the attached column block +*/
#endif
} Task;

/*+ Solver column block structure. +*/

typedef struct SolverCblk_  {
  PASTIX_INT                       fcolnum;              /*+ First column index                     +*/
  PASTIX_INT                       lcolnum;              /*+ Last column index (inclusive)          +*/
  PASTIX_INT                       bloknum;              /*+ First block in column (diagonal)       +*/
  PASTIX_INT                       stride;               /*+ Column block stride                    +*/
  PASTIX_INT			    color;		  /*+ Color of column block (PICL trace)     +*/
#ifdef STARPU_GET_TASK_CTX
  PASTIX_INT                       ctx;                  /*+ Context given to StarPU                +*/
#endif
  PASTIX_INT                       procdiag;             /*+ Processor owner of diagonal block      +*/
  PASTIX_INT                       cblkdiag;             /*+ Column block owner of diagonal block   +*/
  PASTIX_FLOAT * restrict          coeftab;              /*+ Coefficients access vector             +*/
  PASTIX_FLOAT * restrict          ucoeftab;             /*+ Coefficients access vector             +*/
} SolverCblk; 

/*+ Solver block structure. +*/

typedef struct SolverBlok_ {
  PASTIX_INT                       frownum;              /*+ First row index            +*/
  PASTIX_INT                       lrownum;              /*+ Last row index (inclusive) +*/
  PASTIX_INT                       cblknum;              /*+ Facing column block        +*/
  PASTIX_INT                       levfval;              /*+ Level-of-fill value        +*/
  PASTIX_INT                       coefind;              /*+ Index in coeftab           +*/
} SolverBlok;

/*+ Solver matrix structure. +*/

typedef struct SolverMatrix_ {
  PASTIX_INT                       baseval;              /*+ Base value for numberings                 +*/
  PASTIX_INT                       nodenbr;              /*+ Number of nodes in matrix                 +*/
  PASTIX_INT                       cblknbr;              /*+ Number of column blocks                   +*/
  PASTIX_INT                       bloknbr;              /*+ Number of blocks                          +*/
  SolverCblk * restrict     cblktab;              /*+ Array of solver column blocks             +*/
  SolverBlok * restrict     bloktab;              /*+ Array of solver blocks                    +*/
  PASTIX_INT                       coefnbr;              /*+ Number of coefficients                    +*/

  PASTIX_INT                       ftgtnbr;              /*+ Number of fanintargets                    +*/
  PASTIX_INT                       ftgtcnt;              /*+ Number of fanintargets to receive         +*/
  FanInTarget * restrict    ftgttab;              /*+ Fanintarget access vector                 +*/

  PASTIX_INT                       coefmax;              /*+ Working block max size (cblk coeff 1D)    +*/
  PASTIX_INT                       bpftmax;              /*+ Maximum of block size for btag to receive +*/
  PASTIX_INT                       cpftmax;              /*+ Maximum of block size for ftgt to receive +*/
  PASTIX_INT                       nbftmax;              /*+ Maximum block number in ftgt              +*/
  PASTIX_INT                       arftmax;              /*+ Maximum block area in ftgt                +*/

  PASTIX_INT                       clustnum;             /*+ current processor number                  +*/
  PASTIX_INT                       clustnbr;             /*+ number of processors                      +*/
  PASTIX_INT                       procnbr;              /*+ Number of physical processor used         +*/
  PASTIX_INT                       thrdnbr;              /*+ Number of local computation threads       +*/
  PASTIX_INT                       bublnbr;              /*+ Number of local computation threads       +*/
  BubbleTree  * restrict    btree;                /*+ Bubbles tree                              +*/

  BlockTarget * restrict    btagtab;              /*+ Blocktarget access vector                 +*/
  PASTIX_INT                       btagnbr;              /*+ Number of Blocktargets                    +*/
  PASTIX_INT                       btgsnbr;              /*+ Number of Blocktargets to send            +*/
  PASTIX_INT                       btgrnbr;              /*+ Number of Blocktargets to recv            +*/
  BlockCoeff  * restrict    bcoftab;              /*+ BlockCoeff access vector                  +*/
  PASTIX_INT                       bcofnbr;

  PASTIX_INT                       indnbr;
  PASTIX_INT * restrict            indtab;
  Task * restrict           tasktab;              /*+ Task access vector                        +*/
  PASTIX_INT                       tasknbr;              /*+ Number of Tasks                           +*/
  PASTIX_INT **                    ttsktab;              /*+ Task access vector by thread              +*/
  PASTIX_INT *                     ttsknbr;              /*+ Number of tasks by thread                 +*/

  PASTIX_INT *                     proc2clust;           /*+ proc -> cluster                           +*/
  PASTIX_INT                       gridldim;             /*+ Dimensions of the virtual processors      +*/
  PASTIX_INT                       gridcdim;             /*+ grid if dense end block                   +*/
  UpDownVector              updovct;              /*+ UpDown vector                             +*/
#ifdef STARPU_GET_TASK_CTX
  PASTIX_INT                       starpu_subtree_nbr;
#endif
} SolverMatrix;

#endif /* SOLVER_H */
