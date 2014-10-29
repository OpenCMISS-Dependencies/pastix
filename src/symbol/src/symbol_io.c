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
** $Id: symbol_io.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_io.c                             **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the input/output        **/
/**                routines for symbolic matrices.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 23 aug 1998     **/
/**                                 to     07 oct 1998     **/
/**                # Version 0.1  : from : 21 mar 2002     **/
/**                                 to     21 mar 2002     **/
/**                # Version 1.0  : from : 03 jun 2002     **/
/**                                 to     08 sep 2003     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL_IO

#include "common_pastix.h"
#include "symbol.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/*+ This routine loads the given
*** block matrix structure
*** from the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
symbolLoad (
SymbolMatrix * const        symbptr,
FILE * const                stream)
{
  PASTIX_INT                 versval;
  PASTIX_INT                 baseval;
  PASTIX_INT                 nodenbr;
  PASTIX_INT                 cblknbr;
  PASTIX_INT                 cblknum;
  PASTIX_INT                 bloknbr;
  PASTIX_INT                 bloknum;

  if ((intLoad (stream, &versval) +               /* Read header */
       intLoad (stream, &cblknbr) +
       intLoad (stream, &bloknbr) +
       intLoad (stream, &nodenbr) +
       intLoad (stream, &baseval) != 5) ||
      (versval < 0)                     ||        /* Version should be 0 or 1 */
      (versval > 1)                     ||
      (bloknbr < cblknbr)               ||
      (nodenbr < cblknbr)) {
    errorPrint ("symbolLoad: bad input (1)");
    return     (1);
  }

  if (((symbptr->cblktab = (SymbolCblk *) memAlloc ((cblknbr + 1) * sizeof (SymbolCblk))) == NULL) ||
      ((symbptr->bloktab = (SymbolBlok *) memAlloc ( bloknbr      * sizeof (SymbolBlok))) == NULL)) {
    errorPrint ("symbolLoad: out of memory");
    symbolExit (symbptr);
    symbolInit (symbptr);
    return     (1);
  }
  symbptr->baseval = baseval;
  symbptr->cblknbr = cblknbr;
  symbptr->bloknbr = bloknbr;
  symbptr->nodenbr = nodenbr;

  for (cblknum = 0; cblknum < cblknbr; cblknum ++) {
    if ((intLoad (stream, &symbptr->cblktab[cblknum].fcolnum) + /* Read column blocks */
         intLoad (stream, &symbptr->cblktab[cblknum].lcolnum) +
         intLoad (stream, &symbptr->cblktab[cblknum].bloknum) != 3) ||
        (symbptr->cblktab[cblknum].fcolnum > symbptr->cblktab[cblknum].lcolnum)) {
      errorPrint ("symbolLoad: bad input (2)");
      symbolExit (symbptr);
      symbolInit (symbptr);
      return     (1);
    }
  }
  symbptr->cblktab[cblknbr].fcolnum =             /* Set last column block */
  symbptr->cblktab[cblknbr].lcolnum = nodenbr + baseval;
  symbptr->cblktab[cblknbr].bloknum = bloknbr + baseval;

  for (bloknum = 0; bloknum < bloknbr; bloknum ++) {
    if ((intLoad (stream, &symbptr->bloktab[bloknum].frownum) + /* Read column blocks */
         intLoad (stream, &symbptr->bloktab[bloknum].lrownum) +
         intLoad (stream, &symbptr->bloktab[bloknum].cblknum) != 3) ||
        (symbptr->bloktab[bloknum].frownum > symbptr->bloktab[bloknum].lrownum)) {
      errorPrint ("symbolLoad: bad input (3)");
      symbolExit (symbptr);
      symbolInit (symbptr);
      return     (1);
    }

    symbptr->bloktab[bloknum].levfval = 0;        /* Assume version 0 */
    if ((versval > 0) &&
        ((intLoad (stream, &symbptr->bloktab[bloknum].levfval) != 1) ||
         (symbptr->bloktab[bloknum].levfval < 0))) {
      errorPrint ("symbolLoad: bad input (4)");
      symbolExit (symbptr);
      symbolInit (symbptr);
      return     (1);
    }
  }

  return (0);
}

/*+ This routine saves the given
*** block matrix structure
*** to the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
symbolSave (
const SymbolMatrix * const  symbptr,
FILE * const                stream)
{
  const SymbolCblk *  cblktnd;
  const SymbolCblk *  cblkptr;
  const SymbolBlok *  bloktnd;
  const SymbolBlok *  blokptr;
  int                 o;

  o = (fprintf (stream, "1\n%ld\t%ld\t%ld\t%ld\n", /* Write file header */
                (long) symbptr->cblknbr,
                (long) symbptr->bloknbr,
                (long) symbptr->nodenbr,
                (long) symbptr->baseval) == EOF);
  for (cblkptr = symbptr->cblktab, cblktnd = cblkptr + symbptr->cblknbr;
       (cblkptr < cblktnd) && (o == 0); cblkptr ++) {
    o = (fprintf (stream, "%ld\t%ld\t%ld\n",
                  (long) cblkptr->fcolnum,
                  (long) cblkptr->lcolnum,
                  (long) cblkptr->bloknum) == EOF);
  }
  for (blokptr = symbptr->bloktab, bloktnd = blokptr + symbptr->bloknbr;
       (blokptr < bloktnd) && (o == 0); blokptr ++) {
    o = (fprintf (stream, "%ld\t%ld\t%ld\t%ld\n",
                  (long) blokptr->frownum,
                  (long) blokptr->lrownum,
                  (long) blokptr->cblknum,
                  (long) blokptr->levfval) == EOF);
  }

  return (o);
}
