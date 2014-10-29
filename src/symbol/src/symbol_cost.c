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
** $Id: symbol_cost.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_cost.c                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module computes the block solving  **/
/**                cost of symbolic matrices, using the    **/
/**                cost functions of "Fonctions de Comple- **/
/**                xite pour la Resolution par Blocs de    **/
/**                Systemes Lineaires Denses et Creux", by **/
/**                Pierre Ramet.                           **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 14 oct 1998     **/
/**                                 to     15 oct 1998     **/
/**                # Version 1.0  : from : 24 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 1.2  : from : 02 sep 2002     **/
/**                                 to     02 sep 2002     **/
/**                # Version 3.0  : from : 29 sep 2004     **/
/**                                 to     29 sep 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL_COST

#include "common_pastix.h"
#include "dof.h"
#include "symbol.h"
#include "symbol_cost.h"

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
symbolCost (
const SymbolMatrix * const  symbptr,              /*+ Symbolic matrix to evaluate              +*/
const Dof * const           deofptr,              /*+ DOF structure associated with the matrix +*/
const SymbolCostType        typeval,              /*+ Type of cost computation                 +*/
double * const              nnzptr,               /*+ Size of the structure, to be filled      +*/
double * const              opcptr)               /*+ Operation count, to be filled            +*/
{
  if (typeval != SYMBOLCOSTLDLT) {
    errorPrint ("symbolCost: cost function not supported");
    return     (1);
  }

  *opcptr = 0.0L;
  *nnzptr = 0.0L;

  symbolCost2 (symbptr->cblktab - symbptr->baseval, /* Perform recursion on column blocks */
               symbptr->bloktab - symbptr->baseval,
               deofptr, nnzptr, opcptr, symbptr->baseval, symbptr->cblknbr);

  return (0);
}

static
void
symbolCost2 (
const SymbolCblk * restrict const cblktax,        /*+ Based access to cblktab                  +*/
const SymbolBlok * restrict const bloktax,        /*+ Based access to bloktab                  +*/
const Dof * restrict const        deofptr,        /*+ DOF structure associated with the matrix +*/
double * restrict const           nnzptr,         /*+ Size of the structure, to be filled      +*/
double * restrict const           opcptr,         /*+ Operation count, to be filled            +*/
const PASTIX_INT                         cblkmin,        /*+ Minimum column block index to consider   +*/
const PASTIX_INT                         cblknbr)        /*+ Number of column blocks to consider      +*/
{
  PASTIX_INT                 bloknum;                    /* Number of current extra-diagonal block             */
  PASTIX_INT                 cmednum;                    /* Median column block number                         */
  PASTIX_INT                 cfacnum;                    /* Number of facing column block                      */
  PASTIX_INT                 cdofnbr;                    /* Number of DOFs in column block (l_k)               */
  PASTIX_INT                 rdofnbr;                    /* Number of DOFs in row blocks (h_{ki})              */
  PASTIX_INT                 rdofsum;                    /* Number of DOFs in all row blocks (g_{ki} or g_{k}) */
  double              nnzval;                     /* Number of non-zeroes in subtree                    */
  double              opcval;                     /* Operation count in subtree                         */

  nnzval =                                        /* Initialize local values */
  opcval = 0.0L;

  if (cblknbr > 1) {                              /* If more than one column block, perform recursion */
    cmednum = cblknbr / 2;
    symbolCost2 (cblktax, bloktax, deofptr, &nnzval, &opcval, cblkmin, cmednum);
    symbolCost2 (cblktax, bloktax, deofptr, &nnzval, &opcval, cblkmin + cmednum, cblknbr - cmednum);

    *nnzptr += nnzval;                            /* Sum-up local values */
    *opcptr += opcval;
  }
  else {                                          /* Single column block */
    cdofnbr = noddVal (deofptr, cblktax[cblkmin].lcolnum + 1) -
              noddVal (deofptr, cblktax[cblkmin].fcolnum);
    rdofnbr =
    rdofsum = 0;

    for (bloknum = cblktax[cblkmin + 1].bloknum - 1; /* Scan extra-diagonals, backwards */
         bloknum > cblktax[cblkmin].bloknum; ) {
      rdofsum += rdofnbr;
      rdofnbr  = 0;

      cfacnum = bloktax[bloknum].cblknum;         /* Accumulate extra-diagonal blocks facing same column block */
      do {
        rdofnbr += noddVal (deofptr, bloktax[bloknum].lrownum + 1) -
                   noddVal (deofptr, bloktax[bloknum].frownum);
      } while (bloktax[-- bloknum].cblknum == cfacnum);

#ifndef DEAD_CODE
      opcval += (double)(((double) (rdofnbr)) *            /* Count Upsilon_{ki} */
			 ((double) (rdofnbr + 1 + 2 * rdofsum)) * 0.5L +
			 ((double) (rdofnbr + rdofsum)) *  /* Count Lambda_{k}(i) */
			 ((double) (cdofnbr)) * ((double) (2 * rdofnbr + 1)));
#else /* DEAD_CODE */
      opcval += (double)(((double) (rdofnbr)) *            /* Count C3'(k,i) + C3''(k,i) */
			 ((double) (rdofnbr + rdofsum)) *
			 ((double) (2 * cdofnbr + 1)));
#endif /* DEAD_CODE */
    }
    rdofsum += rdofnbr;                           /* Get overall sum */

    *nnzptr += ((double) (cdofnbr + rdofsum)) * ((double) cdofnbr); /* Sum-up stored coefficients */
    *opcptr += (double)(opcval +
			((double) cdofnbr) *               /* Count C1(k) + C2(k) */
			(((double) cdofnbr) * ((double) (2 * cdofnbr + 6 * rdofsum + 3)) + 1.0L) / 6.0L);
  }
}
