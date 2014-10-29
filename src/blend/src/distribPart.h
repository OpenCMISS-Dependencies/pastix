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
/**   NAME       : distribPart.h                           **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Distribution of the symbolic matrix     **/
/**                lead by simulation.                     **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 27 sep 1998     **/
/**                                 to     04 oct 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The function prototypes.
*/

#ifndef DISTRIB
#define static
#endif

void              distribPart              (SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_DIAG            (PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_E1              (PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_E2              (PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_COMP1D          (PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       computeTaskReceiveTime   (const PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static double      computeTaskOnProcReadyTime(PASTIX_INT, PASTIX_INT, SimuCtrl *, BlendCtrl *, Queue *);
static PASTIX_INT        getNextTaskNextProc      (SimuCtrl *, BlendCtrl *, PASTIX_INT *);
static PASTIX_INT        comp_int                 (const PASTIX_INT *, const PASTIX_INT *);
static PASTIX_INT        getNextProc              (SimuProc *, PASTIX_INT);
static PASTIX_INT        getTaskUnmapped          (Queue *, Queue *, SimuCtrl *);
static PASTIX_INT        getCostLowerTask         (SimuCtrl *, BlendCtrl *, PASTIX_INT *);
static PASTIX_INT        chooseCand               (PASTIX_INT, SimuCtrl *, BlendCtrl *);
static void       computeBlockCtrbNbr      (SimuCtrl *, SymbolMatrix *, BlendCtrl *);
static void       updateFtgtStruct         (PASTIX_INT, PASTIX_INT, PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *);
static void       queueReorder             (PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *ctrl);
static void       queueCostReorder         (PASTIX_INT, SimuCtrl *, BlendCtrl *ctrl);
static void       setClusterTime           (PASTIX_INT, SimuCtrl *);
static void       putInReadyQueue          (PASTIX_INT, PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *);
#undef static
