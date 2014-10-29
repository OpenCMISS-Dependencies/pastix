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
/**   NAME       : symbol_faxi.h                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the incomplete symbolic             **/
/**                factorization routine.                  **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 03 jun 2002     **/
/**                                 to     03 jun 2002     **/
/**                                                        **/
/************************************************************/

#define SYMBOL_FAXI_H

#ifdef CXREF_DOC
#include "common_pastix.h"
#ifndef GRAPH_H
#include "graph.h"
#endif /* GRAPH_H */
#ifndef SYMBOL_H
#include "symbol.h"
#endif /* SYMBOL_H */
#ifndef ORDER_H
#include "order.h"
#endif /* ORDER_H */
#ifndef FAX_H
#include "fax.h"
#endif /* FAX_H */
#endif /* CXREF_DOC */

/*
**  The defines.
*/

/* Prime number for hashing vertex numbers. */

#define SYMBOL_FAXI_HASHPRIME       17            /*+ Prime number for hashing +*/

/*
**  The type and structure definitions.
*/

/*+ The chained column block structure. These
    blocks are chained in a single linked list
    for block merge with blocks of left columns. +*/

typedef struct SymbolFaxiTlok_ {
  PASTIX_INT                       frownum;              /*+ First row index            +*/
  PASTIX_INT                       lrownum;              /*+ Last row index (inclusive) +*/
  PASTIX_INT                       cblknum;              /*+ Facing column block        +*/
  PASTIX_INT                       levfval;              /*+ Level of fill of block     +*/
  PASTIX_INT                       nextnum;              /*+ Index of next block        +*/
} SymbolFaxiTlok;
