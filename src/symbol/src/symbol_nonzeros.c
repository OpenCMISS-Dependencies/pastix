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
** $Id: symbol_nonzeros.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_nonzeros.c                       **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This test module writes the coordinates **/
/**                of all non-zeroes so as to check that   **/
/**                several factorization algorithms yield  **/
/**                the same results.                       **/
/**                                                        **/
/**   DATES      : # Version 1.2  : from : 28 aug 2002     **/
/**                                 to     28 aug 2002     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL_NONZEROS

#include "common_pastix.h"
#include "symbol.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/*+ This routine computes the factorization
*** and solving cost of the given symbolic
*** block matrix, whose nodes hold the number
*** of DOFs given by the proper DOF structure.
*** To ensure maximum accuracy and minimum loss
*** of precision, costs are summed-up recursively.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
symbolNonzeros (
const SymbolMatrix * const  symbptr,              /*+ Symbolic matrix to evaluate +*/
FILE * const                stream)               /*+ Output file                 +*/
{
  PASTIX_INT                   cblknum;                  /* Number of current column block */
  SymbolBlok * restrict bloktax;                  /* Based pointer to block array   */
  PASTIX_INT                   bloknum;                  /* Number of current block        */
  PASTIX_INT                   colnum;
  PASTIX_INT                   rownum;

  bloktax = symbptr->bloktab - symbptr->baseval;

  for (cblknum = 0; cblknum < symbptr->cblknbr; cblknum ++) {
    for (colnum  = symbptr->cblktab[cblknum].fcolnum;
         colnum <= symbptr->cblktab[cblknum].lcolnum; colnum ++) {
      for (bloknum = symbptr->cblktab[cblknum].bloknum;
           bloknum < symbptr->cblktab[cblknum + 1].bloknum; bloknum ++) {
        for (rownum  = bloktax[bloknum].frownum;
             rownum <= bloktax[bloknum].lrownum; rownum ++) {
          if (fprintf (stream, "%ld\t%ld\t%c\n",
                       (long) rownum,
                       (long) colnum,
                       (bloktax[bloknum].levfval > 0) ? 'n' : 'z') == EOF)
            return (1);
        }
      }
    }
  }

  return (0);
}
