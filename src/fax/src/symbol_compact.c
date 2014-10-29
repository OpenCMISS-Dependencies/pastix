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
** $Id: symbol_compact.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/*
  File: symbol_compact.c 

  Part of a parallel direct block solver.
  This routine compacts all blocks which
  belong to the same facing column block,
  irrespective of their level of fill.

  Authors:
    - Francois PELLEGRINI

  Dates: 
    Version 1.3 - from 17 oct 2003 to 17 oct 2003
    Version 3.0 - from 29 feb 2004 to 29 feb 2004
*/

/*
**  The defines and includes.
*/

#include "common_pastix.h"
#include "symbol.h"
#include "order.h"
#include "fax.h"

/*
  Function: symbolCompact
  
  This routine merges the given symbolic matrix.

  Parameters:
    symbptr - The symbol matrix

  Returns:
    NO_ERR - on success.
    !0     - on error.
+*/

int
symbolCompact (SymbolMatrix * const        symbptr)
{
  SymbolBlok *              bloktax;
  PASTIX_INT                       cblknum;
  PASTIX_INT                       bloknew;

  bloktax = symbptr->bloktab - symbptr->baseval;

  for (cblknum = 0, bloknew = symbptr->baseval;
       cblknum < symbptr->cblknbr; cblknum ++) {
    PASTIX_INT                       bloknum;
    PASTIX_INT                       bloklst;

    bloknum = symbptr->cblktab[cblknum].bloknum;  /* Update block index in column block array */
    symbptr->cblktab[cblknum].bloknum = bloknew;

    bloktax[bloknew] = bloktax[bloknum];          /* Copy diagonal block to its new position */
    bloklst = bloknew;                            /* Set it as last block                    */

    for (bloknum ++, bloknew ++;
         bloknum < symbptr->cblktab[cblknum + 1].bloknum; bloknum ++) {
      if ((bloktax[bloknum].cblknum == bloktax[bloklst].cblknum) &&
          (bloktax[bloknum].frownum == (bloktax[bloklst].lrownum + 1))) {
        bloktax[bloklst].lrownum = bloktax[bloknum].lrownum;
        if (bloktax[bloknum].levfval < bloktax[bloklst].levfval)
          bloktax[bloklst].levfval = bloktax[bloknum].levfval;
      }
      else {
        bloktax[bloknew] = bloktax[bloknum];      /* Copy block to its new position  */
        bloklst = bloknew ++;                     /* One more block in symbol matrix */
      }
    }
  }
  symbptr->cblktab[cblknum].bloknum = bloknew;    /* Set end of column block array */

  symbptr->bloknbr = bloknew - symbptr->baseval;  /* Update number of blocks */

#ifdef FAX_DEBUG
  if (symbolCheck (symbptr) != 0) {               /* Should not happen */
    errorPrint ("symbolCompact: internal error");
    return (1);
  }
#endif /* FAX_DEBUG */

  return NO_ERR;
}
