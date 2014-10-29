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
** $Id: symbol_levf.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_levf.c                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module computes some parameters of **/
/**                the incomplete elimination tree of      **/
/**                symbolic matrices.                      **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 24 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL_LEVF

#include "common_pastix.h"
#include "symbol.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/*+ This routine computes some parameters
*** of the incomplete elimination tree of
*** the given incomplete symbolic block matrix.
*** - 0   : on success.
*** - !0  : on error.
+*/

int
symbolLevf (
const SymbolMatrix * const  symbptr,              /*+ Symbolic matrix to evaluate                +*/
PASTIX_INT * const                 levfmax,              /*+ Maximum level of fill                      +*/
PASTIX_INT ** const                levfptr)              /*+ Array of number of blocks at level of fill +*/
{
  PASTIX_INT * restrict      levftab;                    /* Array of level of fill counts */
  PASTIX_INT                 levfmac;                    /* Current maximum level of fill */
  PASTIX_INT                 bloknum;                    /* Number of current block       */

  if (levfptr != NULL) {                          /* If level of fill distribution wanted */
    if ((levftab = (PASTIX_INT *) memAlloc (symbptr->cblknbr * sizeof (PASTIX_INT))) == NULL) {
      errorPrint ("symbolLevf: out of memory");
      return     (1);
    }

    memSet (levftab, 0, symbptr->cblknbr * sizeof (PASTIX_INT));

    levfmac = 0;
    for (bloknum = 0; bloknum < symbptr->bloknbr; bloknum ++) {
      PASTIX_INT                 levfval;

      levfval = symbptr->bloktab[bloknum].levfval;

#ifdef SYMBOL_DEBUG
      if (levfval >= symbptr->cblknbr) {
        errorPrint ("symbolLevf: internal error");
        memFree    (levftab);
        return     (1);
      }
#endif /* SYMBOL_DEBUG */

      levftab[levfval] ++;
      if (levfval > levfmac)
        levfmac = levfval;
    }

    *levfptr = (PASTIX_INT *) memRealloc (levftab, (levfmac + 1) * sizeof (PASTIX_INT));
  }
  else {                                          /* Only maximum level of fill wanted */
    levfmac = 0;
    for (bloknum = 0; bloknum < symbptr->bloknbr; bloknum ++) {
      PASTIX_INT                 levfval;

      levfval = symbptr->bloktab[bloknum].levfval;
      if (levfval > levfmac)
        levfmac = levfval;
    }
  }

  *levfmax = levfmac;

  return (0);
}
