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
** $Id: order_base.c 2 2004-06-02 14:05:03Z ramet $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : order_base.c                            **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines hold the base changing      **/
/**                routine for the ordering structure.     **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 06 oct 2003     **/
/**                                 to     06 oct 2003     **/
/**                # Version 2.0  : from : 21 apr 2004     **/
/**                                 to     21 apr 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ORDER

#include "common_pastix.h"
#include "order.h"

/***********************************/
/*                                 */
/* The ordering handling routines. */
/*                                 */
/***********************************/

/* This routine sets the base of the
** given ordering structure to the given
** base value.
** It returns:
** - VOID  : in all cases.
*/

void
orderBase (
Order * restrict const      ordeptr,              /*+ Ordering structure +*/
const PASTIX_INT                   baseval)              /*+ New base value     +*/
{
  PASTIX_INT               baseadj;                      /* Base adjust */
  PASTIX_INT               cblknum;
  PASTIX_INT               vertnbr;
  PASTIX_INT               vertnum;

  if (ordeptr->rangtab == NULL)                   /* Cannot know old base if range array not provided */
    return;

  baseadj = baseval - ordeptr->rangtab[0];        /* Set base adjust     */
  if (baseadj == 0)                               /* If base already set */
    return;

  for (cblknum = 0; cblknum <= ordeptr->cblknbr; cblknum ++)
    ordeptr->rangtab[cblknum] += baseadj;

  vertnbr = ordeptr->rangtab[ordeptr->cblknbr] - ordeptr->rangtab[0];
  if (ordeptr->permtab != NULL) {
    for (vertnum = 0; vertnum < vertnbr; vertnum ++)
      ordeptr->permtab[vertnum] += baseadj;
  }
  if (ordeptr->peritab != NULL) {
    for (vertnum = 0; vertnum < vertnbr; vertnum ++)
      ordeptr->peritab[vertnum] += baseadj;
  }
}
