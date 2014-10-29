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
** $Id: symbol_base.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_base.c                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines hold the base changing      **/
/**                routine for symbol matrices.            **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 jun 1999     **/
/**                                 to     01 jun 1999     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL

#include "common_pastix.h"
#include "symbol.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/* This routine sets the base of the
** given symbol matrix structure to
** the given base value.
** It returns:
** - VOID  : in all cases.
*/

void
symbolBase (
SymbolMatrix * const        symbptr,              /*+ Symbol structure +*/
const PASTIX_INT                   baseval)              /*+ New base value   +*/
{
  PASTIX_INT               baseadj;                      /* Base adjust */
  PASTIX_INT               cblknum;
  PASTIX_INT               bloknum;

  baseadj = baseval - symbptr->baseval;           /* Set base adjust     */
  if (baseadj == 0)                               /* If base already set */
    return;

  symbptr->baseval = baseval;                     /* Set graph base */

  for (cblknum = 0; cblknum <= symbptr->cblknbr; cblknum ++) {
    symbptr->cblktab[cblknum].fcolnum += baseadj;
    symbptr->cblktab[cblknum].lcolnum += baseadj;
    symbptr->cblktab[cblknum].bloknum += baseadj;
  }

  for (bloknum = 0; bloknum < symbptr->bloknbr; bloknum ++) {
    symbptr->bloktab[bloknum].frownum += baseadj;
    symbptr->bloktab[bloknum].lrownum += baseadj;
    symbptr->bloktab[bloknum].cblknum += baseadj;
  }
}
