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
** $Id: symbol_tree.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_tree.c                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module computes some parameters of **/
/**                the elimination tree of symbolic        **/
/**                matrices.                               **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 27 oct 1998     **/
/**                                 to     27 oct 1998     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL_TREE

#include "common_pastix.h"
#include "dof.h"
#include "symbol.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/*+ This routine computes some parameters
*** of the elimination tree of the given
*** symbolic block matrix, whose nodes hold
*** the number of DOFs given by the proper
*** DOF structure.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
symbolTree (
const SymbolMatrix * const  symbptr,              /*+ Symbolic matrix to evaluate              +*/
const Dof * const           deofptr,              /*+ DOF structure associated with the matrix +*/
PASTIX_INT * const                 leafnbrptr,           /*+ Number of leaves in elimination tree     +*/
PASTIX_INT * const                 heigminptr,           /*+ Minimum leaf height in elimination tree  +*/
PASTIX_INT * const                 heigmaxptr,           /*+ Maximum leaf height in elimination tree  +*/
double * const              heigavgptr,           /*+ Average leaf height in elimination tree  +*/
double * const              heigdltptr)           /*+ Deviation of heights                     +*/
{
  const SymbolCblk *  cblktnd;                    /* End of column block array */
  const SymbolCblk *  cblkptr;                    /* Pointer to current        */
  const SymbolBlok *  bloktax;                    /* Based access to bloktab   */
  PASTIX_INT                 nodenum;                    /* Number of current node    */
  unsigned char *              leaftab;                    /* Flag array                */
  unsigned char *              leaftax;                    /* Based access to leaftab   */
  unsigned char *              leafptr;                    /* Pointer to current        */
  PASTIX_INT                 leafnbr;                    /* Number of leaves          */
  PASTIX_INT *               heigtab;                    /* Array of node heights     */
  PASTIX_INT *               heigtax;                    /* Based access to heigtab   */
  PASTIX_INT *               heigtnd;                    /* End of height array       */
  PASTIX_INT *               heigptr;                    /* Pointer to current        */
  PASTIX_INT                 heigval;                    /* Current height value      */
  PASTIX_INT                 heigmin;
  PASTIX_INT                 heigmax;
  double              heigsum;
  double              heigavg;
  double              heigdlt;

  if (((heigtab = (PASTIX_INT *) memAlloc (symbptr->nodenbr * sizeof (PASTIX_INT) + symbptr->cblknbr * sizeof (unsigned char))) == NULL)) {
    errorPrint ("symbolTree: out of memory");
    return     (1);
  }
  leaftab = (unsigned char *) (heigtab + symbptr->nodenbr);

  memSet (leaftab, 0, symbptr->cblknbr * sizeof (unsigned char));
#ifdef SYMBOL_DEBUG
  memSet (heigtab, ~0, symbptr->nodenbr * sizeof (PASTIX_INT));
#endif /* SYMBOL_DEBUG */

  leaftax = leaftab          - symbptr->baseval;
  heigtax = heigtab          - symbptr->baseval;
  bloktax = symbptr->bloktab - symbptr->baseval;

  for (cblkptr = symbptr->cblktab + symbptr->cblknbr - 1,
       heigptr = heigtab + symbptr->nodenbr - 1,
       nodenum = symbptr->nodenbr - 1 + symbptr->baseval;
       cblkptr >= symbptr->cblktab;
       cblkptr --) {
    if (cblkptr[1].bloknum - cblkptr[0].bloknum > 1) { /* If extra-diagonal block present         */
      heigval = heigtax[bloktax[cblkptr->bloknum + 1].frownum] + 1; /* Get height of father, + 1  */
      leaftax[bloktax[cblkptr->bloknum + 1].cblknum] = 1; /* Flag father column block as not leaf */
    }
    else                                          /* Block is root or isolated */
      heigval = 1;

    for (heigtnd = heigtax + cblkptr->fcolnum;
         heigptr >= heigtnd; heigptr --, nodenum --) {
      *heigptr = heigval;                         /* Set node height                   */
      heigval += noddDlt (deofptr, nodenum);      /* Increase height by number of DOFs */
    }
  }

  leafnbr = 0;
  heigmin = INTVALMAX;
  heigmax = 0;
  heigsum = 0.0L;

  for (cblkptr = symbptr->cblktab, leafptr = leaftab, heigptr = heigtab,
       cblktnd = cblkptr + symbptr->cblknbr;
       cblkptr < cblktnd;
       cblkptr ++, leafptr ++) {
    if (*leafptr == 0) {                          /* If column block is leaf  */
      heigval = heigtax[cblkptr->fcolnum];        /* Get height of first node */
      *heigptr ++ = heigval;                      /* Compress height array    */

      leafnbr ++;                                 /* Update leaf data */
      heigsum += (double) heigval;
      if (heigval < heigmin)
        heigmin = heigval;
      if (heigval > heigmax)
        heigmax = heigval;
    }
  }

  heigavg = (double)((leafnbr == 0) ? 0.0L : heigsum / (double) leafnbr);
  heigdlt = 0.0L;
  for (heigptr = heigtab, heigtnd = heigptr + leafnbr;
       heigptr < heigtnd; heigptr ++)
    heigdlt += fabs ((double) *heigptr - heigavg);
  if (leafnbr > 0)
    heigdlt /= (heigavg * (double) leafnbr);

  *leafnbrptr = leafnbr;
  *heigminptr = heigmin;
  *heigmaxptr = heigmax;
  *heigavgptr = heigavg;
  *heigdltptr = heigdlt;

  memFree (heigtab);

  return (0);
}
