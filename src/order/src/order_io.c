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
** $Id: order_io.c 2 2004-06-02 14:05:03Z ramet $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : order_io.c                              **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the input/output        **/
/**                routines for the ordering structure.    **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 04 jan 1999     **/
/**                                 to     05 jan 1999     **/
/**                # Version 2.0  : from : 28 feb 2004     **/
/**                                 to     28 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ORDER_IO

#include "common_pastix.h"
#include "order.h"

/***********************************/
/*                                 */
/* The ordering handling routines. */
/*                                 */
/***********************************/

/*+ This routine reads the given ordering
*** structure from the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
orderLoad (
Order * const               ordeptr,
FILE * const                stream)
{
  PASTIX_INT               versval;                      /* Version number */
  PASTIX_INT               cblknbr;
  PASTIX_INT               cblknum;
  PASTIX_INT               vertnbr;
  PASTIX_INT               vertnum;
  PASTIX_INT               vertnnd;
  PASTIX_INT *             permtax;
  PASTIX_INT *             peritax;
  int               i;

  if ((intLoad (stream, &versval) +
       intLoad (stream, &cblknbr) +
       intLoad (stream, &vertnbr) != 3) ||
      (versval != 0)                    ||
      (cblknbr > vertnbr)) {
    errorPrint ("orderLoad: bad input (1)");
    return     (1);
  }

  if (((ordeptr->rangtab = (PASTIX_INT *) memAlloc ((cblknbr + 1) * sizeof (PASTIX_INT))) == NULL) ||
      ((ordeptr->permtab = (PASTIX_INT *) memAlloc (vertnbr       * sizeof (PASTIX_INT))) == NULL) ||
      ((ordeptr->peritab = (PASTIX_INT *) memAlloc (vertnbr       * sizeof (PASTIX_INT))) == NULL)) {
    errorPrint ("orderLoad: out of memory");
    orderExit  (ordeptr);
    orderInit  (ordeptr);
    return     (1);
  }
  ordeptr->cblknbr = cblknbr;

  for (cblknum = 0, i = 1; (i == 1) && (cblknum <= cblknbr); cblknum ++) /* Read column-block data */
    i = intLoad (stream, &ordeptr->rangtab[cblknum]);

  for (vertnum = 0; (i == 1) && (vertnum < vertnbr); vertnum ++) /* Read direct permutation */
    i = intLoad (stream, &ordeptr->permtab[vertnum]);

  if (i != 1) {
    errorPrint ("orderLoad: bad input (2)");
    orderExit  (ordeptr);
    orderInit  (ordeptr);
    return     (1);
  }

  permtax = ordeptr->permtab - ordeptr->rangtab[0];
  peritax = ordeptr->peritab - ordeptr->rangtab[0];
  for (vertnum = ordeptr->rangtab[0], vertnnd = vertnum + vertnbr; /* Build inverse permutation */
       vertnum < vertnnd; vertnum ++)
    peritax[permtax[vertnum]] = vertnum;

  return (0);
}

/*+ This routine saves the given
*** ordering structure to the
*** given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
orderSave (
const Order * const         ordeptr,
FILE * const                stream)
{
  PASTIX_INT               vertnbr;
  PASTIX_INT               vertnum;
  PASTIX_INT               cblknum;
  int               o;

  if (ordeptr->rangtab == NULL) {
    errorPrint ("orderSave: cannot save ordering without column block data");
    return     (1);
  }
  if (ordeptr->permtab == NULL) {
    errorPrint ("orderSave: cannot save ordering without direct permutation data");
    return     (1);
  }

  vertnbr = ordeptr->rangtab[ordeptr->cblknbr] -  /* Get number of nodes */
            ordeptr->rangtab[0];

  if (fprintf (stream, "0\n%ld\t%ld\n",
               (long) ordeptr->cblknbr,
               (long) vertnbr) == EOF) {
    errorPrint ("orderSave: bad output (1)");
    return     (1);
  }

  for (cblknum = 0, o = 1; (o == 1) && (cblknum < ordeptr->cblknbr); cblknum ++) { /* Save column-block range array */
    o = intSave (stream, ordeptr->rangtab[cblknum]);
    putc (((cblknum & 7) == 7) ? '\n' : '\t', stream);
  }
  o = intSave (stream, ordeptr->rangtab[cblknum]);
  putc ('\n', stream);

  for (vertnum = 0; (o == 1) && (vertnum < (vertnbr - 1)); vertnum ++) { /* Save direct permutation */
    o = intSave (stream, ordeptr->permtab[vertnum]);
    putc (((vertnum & 7) == 7) ? '\n' : '\t', stream);
  }
  o = intSave (stream, ordeptr->permtab[vertnum]);
  putc ('\n', stream);

  if (o != 1)
    errorPrint ("orderSave: bad output (2)");

  return (1 - o);
}
