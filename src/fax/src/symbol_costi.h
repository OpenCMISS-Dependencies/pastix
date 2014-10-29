/* Copyright INRIA 2004
**
** This file is part of the Scotch distribution.
**
** The Scotch distribution is libre/free software; you can
** redistribute it and/or modify it under the terms of the
** GNU Lesser General Public License as published by the
** Free Software Foundation; either version 2.1 of the
** License, or (at your option) any later version.
**
** The Scotch distribution is distributed in the hope that
** it will be useful, but WITHOUT ANY WARRANTY; without even
** the implied warranty of MERCHANTABILITY or FITNESS FOR A
** PARTICULAR PURPOSE. See the GNU Lesser General Public
** License for more details.
**
** You should have received a copy of the GNU Lesser General
** Public License along with the Scotch distribution; if not,
** write to the Free Software Foundation, Inc.,
** 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
**
** $Id: symbol_costi.h 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_costi.h                          **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the incomplete symbolic matrix cost **/
/**                computing routine.                      **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 25 jun 2002     **/
/**                                 to     25 jun 2002     **/
/**                                                        **/
/************************************************************/

#define SYMBOL_COSTI_H

#ifdef CXREF_DOC
#include "common_pastix.h"
#ifndef GRAPH_H
#include "graph.h"
#endif /* GRAPH_H */
#ifndef DOF_H
#include "dof.h"
#endif /* DOF_H */
#ifndef SYMBOL_H
#include "symbol.h"
#endif /* SYMBOL_H */
#endif /* CXREF_DOC */

/*
**  The function prototypes.
*/

#ifndef SYMBOL_COSTI
#define static
#endif

static void                 symbolCosti2        (const SymbolCblk * const cblktax, const SymbolBlok * const bloktax, const Dof * const deofptr, const PASTIX_INT levfval, double * const nnzptr, double * const opcptr, const PASTIX_INT cblkmin, const PASTIX_INT cblknbr);

#undef static
