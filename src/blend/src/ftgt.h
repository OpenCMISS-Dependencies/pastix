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
/**   NAME       : ftgt.h                                  **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the graph structure.                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 02 may 1998     **/
/**                                 to     16 oct 1998     **/
/**                                                        **/
/************************************************************/

#ifndef FTGT_H
#define FTGT_H

#ifdef CXREF_DOC
#include "common_pastix.h"
#endif /* CXREF_DOC */


/*
**  The type and structure definitions.
*/

/*
** WARNING : All structures should have a odd number of integer for memory alignment
*/

/*+ Fanintarget info type +*/

typedef enum {
  FTGT_CTRBNBR = 0,                               /*+ Number of contributions           +*/
  FTGT_CTRBCNT,                                   /*+ Number of contributions remaining +*/
  FTGT_PROCDST,                                   /*+ Destination for fanintarget       +*/
  FTGT_TASKDST,                                   /*+ Task  destination                 +*/
  FTGT_BLOKDST,                                   /*+ Block destination (->COMP_1D)     +*/
  FTGT_PRIONUM,                                   /*+ Fanintarget priority              +*/
  FTGT_FCOLNUM,                                   /*+ Fanintarget first column          +*/
  FTGT_LCOLNUM,                                   /*+ Fanintarget last column           +*/
  FTGT_FROWNUM,                                   /*+ Fanintarget first row             +*/
  FTGT_LROWNUM,                                   /*+ Fanintarget last row              +*/
#if (defined OOC) || (defined TRACE_SOPALIN)
  FTGT_GCBKDST,                                   /*+ Global Cblk destination(->COMP_1D)+*/
  FTGT_IDTRACE,                                   /*+ To have 12 integer in FanInTarget +*/
#endif
  MAXINFO
} FanInInfo;

typedef enum {
  BTAG_PRIONUM = 0,
  BTAG_TASKDST,
  BTAG_PROCDST,
  BTAG_TASKCNT,
#if (defined TRACE_SOPALIN)
  BTAG_NULL,                                      /*+ Global Cblk destination(->COMP_1D)+*/
  BTAG_IDTRACE,                                   /*+ To have 12 integer in FanInTarget +*/
#endif
  BTAGINFO
} BtagInfo;

typedef enum {
  BCOF_FROWNUM = 0,
  BCOF_LROWNUM,
  BCOF_FCOLNUM,
  BCOF_LCOLNUM,
  BCOFINFO
} BcofInfo;




/*+ Fanintarget structure +*/

typedef struct FanInTarget_ {
  PASTIX_INT                       infotab[MAXINFO];     /*+ Fanintarget descriptor (size MAXINFO) +*/
  PASTIX_FLOAT *                   coeftab;              /*+ Fanintarget vector access             +*/
} FanInTarget;

typedef struct BlockCoeff_ 
{
  PASTIX_INT              infotab[BCOFINFO];
  volatile PASTIX_FLOAT * coeftab; /* blocktarget coeff vector if != NULL envoi possible */
  PASTIX_INT              sendcnt; /* number of blocktarget send, if == 0 free coeftab */
} BlockCoeff;


/*+ BlockTarget structure +*/
typedef struct BlockTarget_ {
  PASTIX_INT           infotab[BTAGINFO];
  BlockCoeff*   bcofptr; /* index of the blockCoeff                        */
} BlockTarget;

#endif /* FTGT_H */
