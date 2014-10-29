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
/**   NAME       : fanboth.h                               **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Compute a moderate amalgamation in the  **/
/**                fan-in target to save memory            **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 April 2000   **/
/**                                 to     20 April 2000   **/
/**                                                        **/
/************************************************************/
#ifndef FANBOTH
#define static
typedef struct {
  PASTIX_INT ctrbnbr; /*+ number of contribution of the partial ftgt +*/ 
  PASTIX_INT ctrbcnt; 
  PASTIX_INT prionum; /*+ Priority of the partial ftgt +*/
  PASTIX_INT indnum;  /*+ index where the ftgt must be insert in the initial ftgttab +*/
  PASTIX_INT ftgtnum; /*+ Index of the initial ftgt from which is issue the partial ftgt +*/
  PASTIX_INT ftgtnewnum; /*+ index of the ftgt in the final ftgttab +*/
  PASTIX_INT next;    /*+ Chain to the next partial ftgt of the initial ftgt +*/
} ExtraFtgt;
PASTIX_INT Malt2(SolverMatrix *, double);
PASTIX_INT getFtgtInd2(SolverMatrix *, PASTIX_INT *, Queue *, Queue *);
PASTIX_INT getFtgtNextAccess(PASTIX_INT ind, PASTIX_INT ftgtaccessnbr, PASTIX_INT *ftgtaccesstab);
#undef static
#endif
