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
** $Id: symbol_costi.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/*
  File: symbol_costi.c                          
  
  Part of a parallel direct block solver.
  This module computes the block solving
  cost of incomplete symbolic matrices,
  using the cost functions of "Fonctions
  de Complexite pour la Resolution par
  Blocs de Systemes Lineaires Denses et
  Creux", by Pierre Ramet.

  Authors: 
    - Francois PELLEGRINI

  Dates: 
    Version 1.0 - from 24 jun 2002 to 26 jun 2002
    Version 3.0 - from 29 feb 2004 to 29 feb 2004

*/

/*
**  The defines and includes.
*/

#define SYMBOL_COSTI

#include "common_pastix.h"
#include "dof.h"
#include "symbol.h"
#include "symbol_costi.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/*
  Function: symbolCosti

  This routine computes the factorization
  and solving cost of the given symbolic
  block matrix, whose nodes hold the number
  of DOFs given by the proper DOF structure.

  To ensure maximum accuracy and minimum loss
  of precision, costs are summed-up recursively.

  In fact it calls <symbolCosti2>.

  Parameters:
    symbptr - Symbolic matrix to evaluate             
    deofptr - DOF structure associated with the matrix
    typeval - Type of cost computation                
    levfval - Level of fill                           
    nnzptr  - Size of the structure, to be filled     
    opcptr  - Operation count, to be filled

  Returns:
    0  - on success.
    !0 - on error.
*/

int
symbolCosti (const SymbolMatrix * const  symbptr,
	     const Dof * const           deofptr,
	     const SymbolCostType        typeval,
	     const PASTIX_INT                   levfval,
	     double * const              nnzptr,
	     double * const              opcptr)
{
  if (typeval != SYMBOLCOSTLDLT) {
    errorPrint ("symbolCosti: cost function not supported");
    return     (1);
  }

  *opcptr = 0.0L;
  *nnzptr = 0.0L;

  symbolCosti2 (symbptr->cblktab - symbptr->baseval, /* Perform recursion on column blocks */
                symbptr->bloktab - symbptr->baseval,
                deofptr, levfval, nnzptr, opcptr, symbptr->baseval, symbptr->cblknbr);

  return (0);
}

/*
  Function: symbolCosti2

  This routine computes the factorization
  and solving cost of the given symbolic
  block matrix, whose nodes hold the number
  of DOFs given by the proper DOF structure.

  To ensure maximum accuracy and minimum loss
  of precision, costs are summed-up recursively.

  Parameters:
    cblktax - Based access to cblktab
    bloktax - Based access to bloktab
    deofptr - DOF structure associated with the matrix
    levfval - Level of fill
    nnzptr  - Size of the structure, to be filled
    opcptr  - Operation count, to be filled
    cblkmin - Minimum column block index to consider
    cblknbr - Number of column blocks to consider
*/
static
void
symbolCosti2 (const SymbolCblk * restrict const cblktax,
	      const SymbolBlok * restrict const bloktax,
	      const Dof * restrict const        deofptr,
	      const PASTIX_INT                         levfval,
	      double * restrict const           nnzptr,
	      double * restrict const           opcptr,
	      const PASTIX_INT                         cblkmin,
	      const PASTIX_INT                         cblknbr)
{
  PASTIX_INT                 bloknum;                    /* Number of current extra-diagonal block             */
  PASTIX_INT                 cmednum;                    /* Median column block number                         */
  PASTIX_INT                 cfacnum;                    /* Number of facing column block                      */
  PASTIX_INT                 cdofnbr;                    /* Number of DOFs in column block (l_k)               */
  PASTIX_INT                 rdofsum;                    /* Number of DOFs in all row blocks (g_{ki} or g_{k}) */
  double              nnzval;                     /* Number of non-zeroes in subtree                    */
  double              opcval;                     /* Operation count in subtree                         */

  nnzval =                                        /* Initialize local values */
  opcval = 0.0L;

  if (cblknbr > 1) {                              /* If more than one column block, perform recursion */
    cmednum = cblknbr / 2;
    symbolCosti2 (cblktax, bloktax, deofptr, levfval, &nnzval, &opcval, cblkmin, cmednum);
    symbolCosti2 (cblktax, bloktax, deofptr, levfval, &nnzval, &opcval, cblkmin + cmednum, cblknbr - cmednum);

    *nnzptr += nnzval;                            /* Sum-up local values */
    *opcptr += opcval;
  }
  else {                                          /* Single column block                              */
    PASTIX_INT                 levffac;                  /* Minimum level of fill over facing block(s)       */
    PASTIX_INT                 rdounbr;                  /* Number of DOFs in undropped row blocks (h'_{ki}) */
    PASTIX_INT                 rdousum;                  /* Number of DOFs in undropped row blocks (h'_{ki}) */

    cdofnbr = noddVal (deofptr, cblktax[cblkmin].lcolnum + 1) -
              noddVal (deofptr, cblktax[cblkmin].fcolnum);

    bloknum = cblktax[cblkmin].bloknum + 1;       /* Get index of first extra-diagonal block */

    if (bloknum == cblktax[cblkmin + 1].bloknum)  /* If diagonal block only  */
      levffac = levfval;                          /* No contributions at all */
    else {
      levffac = bloktax[bloknum].levfval;         /* Get level of fill of first extra-diagonal block                */
      for (bloknum ++; (bloknum < cblktax[cblkmin + 1].bloknum) && /* For all first blocks facing same column block */
           (bloktax[bloknum].cblknum == bloktax[bloknum - 1].cblknum); bloknum ++) {
        if (bloktax[bloknum].levfval < levffac)   /* If facing block has smaller level of fill */
          levffac = bloktax[bloknum].levfval;     /* Keep smallest level of fill               */
      }
    }

    rdofsum =
    rdousum =
    rdounbr = 0;

    for (bloknum = cblktax[cblkmin + 1].bloknum - 1; /* Scan extra-diagonals, backwards */
         bloknum > cblktax[cblkmin].bloknum; ) {
      if (bloktax[bloknum].levfval > levfval) {   /* Skip dropped blocks */
        bloknum --;
        continue;
      }

      rdousum += rdounbr;
      rdounbr  = 0;

      cfacnum = bloktax[bloknum].cblknum;
      do {
        PASTIX_INT                 rdofblk;              /* Number of DOFs in local block */

        if (bloktax[bloknum].levfval > levfval)   /* Skip dropped blocks */
          continue;

        rdofblk  = noddVal (deofptr, bloktax[bloknum].lrownum + 1) -
                   noddVal (deofptr, bloktax[bloknum].frownum);
        rdofsum += rdofblk;                       /* Account for undropped blocks */

        if (MAX (bloktax[bloknum].levfval, levffac) < levfval) /* If contribution is not dropped */
          rdounbr += rdofblk;                     /* Add contribution to undropped set of blocks */
      } while (bloktax[-- bloknum].cblknum == cfacnum);

      opcval += ((double) (rdounbr)) *            /* Count C3'(k,i) + C3''(k,i) */
                ((double) (rdounbr + rdousum)) *
                ((double) (2 * cdofnbr + 1));
    }

    *nnzptr += ((double) (cdofnbr + rdofsum)) * ((double) cdofnbr); /* Sum-up stored coefficients */
    *opcptr += (double)(opcval +
			((double) cdofnbr) *               /* Count C1(k) + C2(k) */
			(((double) cdofnbr) * ((double) (2 * cdofnbr + 6 * rdofsum + 3)) + 1.0L) / 6.0L);
  }
}
