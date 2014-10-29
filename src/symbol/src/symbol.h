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
** $Id: symbol.h 316 2005-06-06 16:17:44Z ramet $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol.h                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the symbolic matrix.                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     07 oct 1998     **/
/**                # Version 0.1  : from : 21 mar 2002     **/
/**                                 to     21 mar 2002     **/
/**                # Version 1.0  : from : 03 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 1.3  : from : 10 apr 2003     **/
/**                                 to     10 jun 2003     **/
/**                # Version 3.0  : from : 28 feb 2004     **/
/**                                 to     03 mar 2005     **/
/**                                                        **/
/************************************************************/

#ifndef SYMBOL_H
#define SYMBOL_H
#define SYMBOL_VERSION              1

/*
**  The type and structure definitions.
*/

/*+ The column block structure. +*/

typedef struct SymbolCblk_ {
  PASTIX_INT                       fcolnum;              /*+ First column index               +*/
  PASTIX_INT                       lcolnum;              /*+ Last column index (inclusive)    +*/
  PASTIX_INT                       bloknum;              /*+ First block in column (diagonal) +*/
#ifdef STARPU_GET_TASK_CTX
  PASTIX_INT                       ctx;
#endif
} SymbolCblk;

/*+ The column block structure. +*/

typedef struct SymbolBlok_ {
  PASTIX_INT                       frownum;              /*+ First row index            +*/
  PASTIX_INT                       lrownum;              /*+ Last row index (inclusive) +*/
  PASTIX_INT                       cblknum;              /*+ Facing column block        +*/
  PASTIX_INT                       levfval;              /*+ Level-of-fill value        +*/
} SymbolBlok;

/*+ The symbolic block matrix. +*/

typedef struct SymbolMatrix_ {
  PASTIX_INT                       baseval;              /*+ Base value for numberings         +*/
  PASTIX_INT                       cblknbr;              /*+ Number of column blocks           +*/
  PASTIX_INT                       bloknbr;              /*+ Number of blocks                  +*/
  SymbolCblk * restrict     cblktab;              /*+ Array of column blocks [+1,based] +*/
  SymbolBlok * restrict     bloktab;              /*+ Array of blocks [based]           +*/
  PASTIX_INT                       nodenbr;              /*+ Number of nodes in matrix         +*/
#ifdef STARPU_GET_TASK_CTX
  PASTIX_INT                       starpu_subtree_nbr;
#endif
} SymbolMatrix;

/*+ The type of cost computations. +*/

typedef enum SymbolCostType_ {
  SYMBOLCOSTLDLT                                  /*+ Crout (i.e. LDLt) cost function +*/
} SymbolCostType;

/* Structure for keeping track of selected
   blocks in the matrix pattern. The values
   of the tables are the remaining values
   for the yet unselected blocks.           */

typedef struct SymbolKeepBlok_ {
  PASTIX_INT                       levfval;              /*+ Values for incomplete factorisation +*/
  PASTIX_INT                       nupdval;
  PASTIX_INT                       ctrival;
  PASTIX_INT                       ctroval;
  PASTIX_INT                       hghtval;
} SymbolKeepBlok;

typedef struct SymbolKeep_ {
  PASTIX_INT                       levfmax;              /*+ Maximum values for incomplete fax +*/
  PASTIX_INT                       nupdmax;
  PASTIX_INT                       ctrimax;
  PASTIX_INT                       ctromax;
  PASTIX_INT                       hghtmax;
  unsigned char * restrict           keeptab;              /*+ Flag array for kept blocks      +*/
  SymbolKeepBlok * restrict kblktab;              /*+ Block parameter array           +*/
  double * restrict         levftab;              /*+ Area arrays for selected blocks +*/
  double * restrict         nupdtab;
  double * restrict         ctritab;
  double * restrict         ctrotab;
  double * restrict         hghttab;
} SymbolKeep;

/*
**  The function prototypes.
*/

int  symbolInit          (SymbolMatrix * const symbptr);
void symbolExit          (SymbolMatrix * const symbptr);
void symbolBase          (SymbolMatrix * const symbptr, const PASTIX_INT baseval);
void symbolRealloc       (SymbolMatrix * const symbptr);
int  symbolLoad          (SymbolMatrix * const symbptr, FILE * const stream);
int  symbolSave          (const SymbolMatrix * const symbptr, FILE * const stream);
int  symbolCheck         (const SymbolMatrix * const symbptr);
int  symbolDraw          (const SymbolMatrix * const symbptr, FILE * const stream);
int  symbolDrawFunc      (const SymbolMatrix * const symbptr,
			  int (*) (const SymbolMatrix * const, const SymbolBlok * const, void * const, float * const),
			  int (*) (const SymbolMatrix * const, const SymbolBlok * const, void * const, float * const),
			  void * const, FILE * const stream);
void symbolDrawColor     (const PASTIX_INT labl, float color[]);
#ifdef DOF_H
int  symbolCost          (const SymbolMatrix * const symbptr, const Dof * const deofptr,
			  const SymbolCostType typeval, double * const nnzptr, double * const opcptr);
int  symbolCosti         (const SymbolMatrix * const symbptr, const Dof * const deofptr,
			  const SymbolCostType typeval, const PASTIX_INT levfval, double * const nnzptr, double * const opcptr);
int  symbolLevf          (const SymbolMatrix * const symbptr, PASTIX_INT * const levfmax, PASTIX_INT ** const levftab);
int  symbolTree          (const SymbolMatrix * const symbptr, const Dof * const deofptr,
			  PASTIX_INT * const leafnbr, PASTIX_INT * const heigmin, PASTIX_INT * const heigmax,
			  double * const heigavg, double * const heigdlt);
int  symbolNonzeros      (const SymbolMatrix * const symbptr, FILE * const stream);
#endif /* DOF_H */

int  symbolKeepInit      (SymbolKeep * restrict const keepptr, const SymbolMatrix * const symbptr);
void symbolKeepExit      (SymbolKeep * restrict const keepptr);
void symbolKeepAdd       (SymbolKeep * restrict const keepptr, const SymbolMatrix * const symbptr,
			  int (* funcptr) (const SymbolKeepBlok * const, void * const), void * dataptr);
void symbolKeepDel       (SymbolKeep * restrict const keepptr, const SymbolMatrix * const symbptr,
			  int (* funcptr) (const SymbolKeepBlok * const, void * const), void * dataptr);
int  symbolKeepCompute   (SymbolKeep * restrict const keepptr, const SymbolMatrix * const symbptr);
int  symbolKeepHisto     (SymbolKeep * const keepptr, const SymbolMatrix * const,
			  int (* funcptr) (const SymbolKeepBlok * const, void * const), void * dataptr);
int  symbolKeepPurge     (SymbolKeep * restrict const keepptr, SymbolMatrix * restrict const symbptr);
int  symbolKeepView      (const SymbolKeep * const keepptr, const double nnzlmax, const char * const nameptr);

#endif /* SYMBOL_H */
