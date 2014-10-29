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
** $Id: symbol_check.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_check.c                          **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module checks the consistency of   **/
/**                symbolic matrices.                      **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 29 sep 1998     **/
/**                                 to     07 oct 1998     **/
/**                # Version 1.0  : from : 03 jun 2002     **/
/**                                 to     03 jun 2002     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL_CHECK

#include "common_pastix.h"
#include "symbol.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/*+ This routine checks the consistency
*** of the given symbolic block matrix.
*** Because of incomplete factorization,
*** from version 1.0, no check is performed
*** regarding the existence of facing blocks
*** in facing columns.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
symbolCheck (
const SymbolMatrix * const  symbptr)
{
  PASTIX_INT                         baseval;            /* Base value                           */
  const SymbolCblk * restrict cblktax;            /* Based access to cblktab              */
  PASTIX_INT                         cblkmax;            /* Maximum column block index           */
  PASTIX_INT                         cblknum;            /* Based number of current column block */
  const SymbolBlok * restrict bloktax;            /* Based access to bloktab              */
  PASTIX_INT                         blokmax;            /* Maximum block index                  */
  PASTIX_INT                         bloknum;            /* Based number of current block        */
  PASTIX_INT                         nodemax;            /* Maximum node index                   */

  baseval = symbptr->baseval;
  cblktax = symbptr->cblktab - baseval;
  cblkmax = symbptr->cblknbr + (baseval - 1);
  bloktax = symbptr->bloktab - baseval;
  blokmax = symbptr->bloknbr + baseval;
  nodemax = symbptr->nodenbr + (baseval - 1);

  for (cblknum = bloknum = baseval;
       cblknum <= cblkmax; cblknum ++) {
    if ((cblktax[cblknum].fcolnum     <  baseval)                  ||
        (cblktax[cblknum].lcolnum     >  nodemax)                  ||
        (cblktax[cblknum].bloknum     >  blokmax)                  ||
        (cblktax[cblknum].fcolnum     >  cblktax[cblknum].lcolnum) ||
        (cblktax[cblknum + 1].fcolnum <= cblktax[cblknum].lcolnum) ||
        (cblktax[cblknum + 1].bloknum <= cblktax[cblknum].bloknum)) {
      errorPrint ("symbolCheck: invalid column block array");
      return     (1);
    }

    if ((bloktax[bloknum].frownum != cblktax[cblknum].fcolnum) ||
        (bloktax[bloknum].lrownum != cblktax[cblknum].lcolnum) ||
        (bloktax[bloknum].cblknum != cblknum)                  ||
        (bloktax[bloknum].levfval != 0)) {
      errorPrint ("symbolCheck: invalid diagonal block");
      return     (1);
    }

    for (bloknum ++; bloknum < cblktax[cblknum + 1].bloknum; bloknum ++) {
      if ((bloktax[bloknum].cblknum <  baseval)                      ||
          (bloktax[bloknum].cblknum >  cblkmax)                      ||
          (bloktax[bloknum].frownum <= bloktax[bloknum - 1].lrownum) ||
          (bloktax[bloknum].cblknum <  bloktax[bloknum - 1].cblknum)) {
        errorPrint ("symbolCheck: invalid block array");
        return     (1);
      }
    }
  }

  return (0);
}
