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
** $Id: order_check.c 2 2004-06-02 14:05:03Z ramet $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : order_check.c                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module checks the consistency of   **/
/**                orderings.                              **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 06 oct 1998     **/
/**                                 to     06 oct 1998     **/
/**                # Version 1.0  : from : 19 nov 2003     **/
/**                                 to     20 nov 2003     **/
/**                # Version 2.0  : from : 28 feb 2004     **/
/**                                 to     28 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ORDER_CHECK

#include "common_pastix.h"
#include "order.h"

/***********************************/
/*                                 */
/* The ordering handling routines. */
/*                                 */
/***********************************/

/*+ This routine checks the consistency
*** of the given ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
orderCheck (
const Order * restrict const  ordeptr)
{
  PASTIX_INT                   baseval;                  /* Node base value            */
  PASTIX_INT                   vnodmax;                  /* Maximum node value         */
  PASTIX_INT                   vnodnum;                  /* Number of current node     */
  PASTIX_INT                   rangnum;                  /* Current column block index */
  const PASTIX_INT * restrict  peritax;                  /* Based access to peritab    */
  const PASTIX_INT * restrict  permtax;                  /* Based access to permtab    */

  if (ordeptr->cblknbr < 0) {
    errorPrint ("orderCheck: invalid nunber of column blocks");
    return     (1);
  }

  baseval = ordeptr->rangtab[0];                  /* Get base value */
  if (baseval < 0) {
    errorPrint ("orderCheck: invalid vertex node base number");
    return     (1);
  }

  peritax = ordeptr->peritab - baseval;           /* Set based accesses */
  vnodmax = ordeptr->rangtab[ordeptr->cblknbr] - 1;

  for (rangnum = 0; rangnum < ordeptr->cblknbr; rangnum ++) {
    if ((ordeptr->rangtab[rangnum] <  baseval) ||
        (ordeptr->rangtab[rangnum] >  vnodmax) ||
        (ordeptr->rangtab[rangnum] >= ordeptr->rangtab[rangnum + 1])) {
      errorPrint ("orderCheck: invalid range array");
      return     (1);
    }
  }

  permtax = ordeptr->permtab - baseval;

  for (vnodnum = baseval;
       vnodnum <= vnodmax; vnodnum ++) {
    PASTIX_INT                   vnodold;

    vnodold = peritax[vnodnum];
    if ((vnodold < baseval) ||
        (vnodold > vnodmax) ||
        (permtax[vnodold] != vnodnum)) {
      errorPrint ("orderCheck: invalid permutation arrays");
      return     (1);
    }
  }

  return (0);
}
