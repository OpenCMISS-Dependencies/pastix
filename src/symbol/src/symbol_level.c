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
** $Id: symbol_level.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_level.c                          **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module contains the incomplete     **/
/**                symbol matrix pruning routine.          **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 03 jun 2002     **/
/**                                 to     03 jun 2002     **/
/**                # Version 1.1  : from : 26 jun 2002     **/
/**                                 to     03 mar 2005     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL_LEVEL

#include "common_pastix.h"
#include "symbol.h"
#include "order.h"
#include "fax.h"

/*************************/
/*                       */
/* Level making routine. */
/*                       */
/*************************/

/*+ This routine computes an incomplete
*** block symbolic factorization structure
*** of given level by removing all blocks
*** the fill level of which is strictly
*** greater than the given cut-off value.
*** In order for the elimination tree to
*** have the same tree structure as the
*** complete factored matrix, first
*** extra-diagonal blocks of every column
*** are never removed.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
symbolLevel (
SymbolMatrix * const        dstsymbptr,           /*+ New symbolic block matrix [based] +*/
const SymbolMatrix * const  srcsymbptr,           /*+ Old symbolic block matrix [based] +*/
const PASTIX_INT                   levfval)              /*+ Cut-off level of fill             +*/
{
  PASTIX_INT                         baseval;            /* Base value                             */
  PASTIX_INT                         cblknum;            /* Based number of current column block   */
  const SymbolCblk * restrict srccblktax;         /* Based access to old column block array */
  SymbolCblk * restrict       dstcblktax;         /* Based access to new column block array */
  const SymbolBlok * restrict srcbloktax;         /* Based access to old block array        */
  SymbolBlok * restrict       dstbloktax;         /* Based access to new block array        */
  PASTIX_INT                         srcbloknum;         /* Based number of current old block      */
  PASTIX_INT                         dstbloknum;         /* Based number of current new block      */

  if (((dstsymbptr->cblktab = (SymbolCblk *) memAlloc ((srcsymbptr->cblknbr + 1) * sizeof (SymbolCblk))) == NULL) ||
      ((dstsymbptr->bloktab = (SymbolBlok *) memAlloc ( srcsymbptr->bloknbr      * sizeof (SymbolBlok))) == NULL)) {
    errorPrint ("symbolLevel: out of memory");
    if (dstsymbptr->cblktab != NULL) {
      memFree (dstsymbptr->cblktab);
      dstsymbptr->cblktab = NULL;
    }
    return (1);
  }
  baseval    = srcsymbptr->baseval;
  srccblktax = srcsymbptr->cblktab - baseval;     /* Set based accesses */
  dstcblktax = dstsymbptr->cblktab - baseval;
  srcbloktax = srcsymbptr->bloktab - baseval;
  dstbloktax = dstsymbptr->bloktab - baseval;

  for (cblknum = srcbloknum = dstbloknum = baseval; /* For all column blocks */
       cblknum < srcsymbptr->cblknbr + baseval; cblknum ++) {
    dstcblktax[cblknum].fcolnum = srccblktax[cblknum].fcolnum; /* Build new column block data */
    dstcblktax[cblknum].lcolnum = srccblktax[cblknum].lcolnum;
    dstcblktax[cblknum].bloknum = dstbloknum;

    dstbloktax[dstbloknum].frownum = dstcblktax[cblknum].fcolnum; /* Build diagonal block */
    dstbloktax[dstbloknum].lrownum = dstcblktax[cblknum].lcolnum;
    dstbloktax[dstbloknum].cblknum = cblknum;
    dstbloktax[dstbloknum].levfval = 0;

    if (srcbloknum >= srccblktax[cblknum + 1].bloknum) /* If no extra-diagonal block */
      continue;                                   /* Proceed with next column block  */

    dstbloktax[dstbloknum].frownum = srcbloktax[srcbloknum].frownum; /* Copy data of first column block as is */
    dstbloktax[dstbloknum].lrownum = srcbloktax[srcbloknum].lrownum;
    dstbloktax[dstbloknum].cblknum = cblknum;
    dstbloktax[dstbloknum].levfval = srcbloktax[srcbloknum].levfval;
    dstbloknum ++;
    srcbloknum ++;

    for ( ; srcbloknum < srccblktax[cblknum + 1].bloknum; srcbloknum ++) {
      if (srcbloktax[srcbloknum].levfval <= levfval) { /* If block should be kept */
        dstbloktax[dstbloknum].frownum = srcbloktax[srcbloknum].frownum;
        dstbloktax[dstbloknum].lrownum = srcbloktax[srcbloknum].lrownum;
        dstbloktax[dstbloknum].cblknum = cblknum;
        dstbloktax[dstbloknum].levfval = srcbloktax[srcbloknum].levfval;
        dstbloknum ++;
      }
    }
  }
  dstcblktax[cblknum].fcolnum = srccblktax[cblknum].fcolnum; /* Build terminating column block */
  dstcblktax[cblknum].lcolnum = srccblktax[cblknum].lcolnum;
  dstcblktax[cblknum].bloknum = dstbloknum;

  dstsymbptr->baseval = baseval;                  /* Fill in matrix fields */
  dstsymbptr->cblknbr = srcsymbptr->cblknbr;
  dstsymbptr->bloknbr = dstbloknum - baseval;
  dstsymbptr->nodenbr = srcsymbptr->nodenbr;

#ifdef SYMBOL_DEBUG
  if (symbolCheck (dstsymbptr) != 0) {
    errorPrint ("symbolLevel: internal error");
    symbolExit (dstsymbptr);
    return     (1);
  }
#endif /* SYMBOL_DEBUG */

  return (0);
}
