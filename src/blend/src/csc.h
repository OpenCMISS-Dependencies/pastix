/* Copyright 2008 BORDEAUX I UNIVERSITY & INRIA 
**
** This file is part of the PaStiX parallel sparse matrix package.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
#ifndef CSC_H
#define CSC_H
#include "common_pastix.h"

/* Section: Macros */
/*
  Macro: CSC_FNBR

  Accessor to the number of column block.

  Parameters:
    a - Pointer to the CSC.

  Returns:
    Number of column block.
 */
#define CSC_FNBR(a)     (a)->cscfnbr /* cblk nbr */
/*
  Macro: CSC_FTAB

  Accessor to the array of column blocks.

  Parameters:
    a - Pointer to the CSC.

  Returns:
    Address of the array of column blocks.
*/
#define CSC_FTAB(a)     (a)->cscftab
/*
  Macro: CSC_COLNBR

  Accessor to the number of column in a block column.

  Parameters:
    a - Pointer to the CSC matrix.
    b - Column block index.

  Returns:
    The number of column in the block column
*/
#define CSC_COLNBR(a,b) (a)->cscftab[b].colnbr
/*
  Macro: CSC_COLTAB

  Accessor to the array of start for each column
  in the rows and values arrays.

  Parameters:
    a - Pointer to the CSC matrix.
    b - Column block index.

  Reurns:
    Address of the array of indexes of start for each column
    in the rows and values arrays.

*/
#define CSC_COLTAB(a,b) (a)->cscftab[b].coltab
/*
  Macro: CSC_COL

  Accessor to the index of first element of a column in rows
  and values.

  Parameters:
    a - Pointer to the CSC matrix.
    b - Column block index.
    c - Column index.
*/
#define CSC_COL(a,b,c)  (a)->cscftab[b].coltab[c]
/*
   Macro: CSC_ROWTAB

   Accessor to the array of rows.

   Parameters:
     a - Pointer to the CSC matrix.
*/
#define CSC_ROWTAB(a)   (a)->rowtab
/*
   Macro: CSC_ROW

   Accessor to a row in the CSC.

   Parameters:
     a - Pointer to the CSC matrix.
     b - Index of the row.
*/
#define CSC_ROW(a,b)    (a)->rowtab[b]
/*
   Macro: CSC_VALTAB

   Accessor to the array of values.

   Parameters:
     a - Pointer to the CSC matrix.
*/
#define CSC_VALTAB(a)   (a)->valtab
/*
   Macro: CSC_VAL

   Accessor to a value in the CSC.

   Parameters:
     a - Pointer to the CSC matrix.
     b - Index of the value.
*/
#define CSC_VAL(a,b)    (a)->valtab[b]
/*
   Macro: CSC_FROW

   Accessor to the first row of the column $c$ in
   the column block $b$.

   Parameters:
     a - Pointer to the CSC matrix.
     b - Column block index.
     c - Column index.
*/
#define CSC_FROW(a,b,c) (a)->rowtab[(a)->cscftab[b].coltab[c]]
/*
   Macro: CSC_FROW

   Accessor to the first value of the column $c$ in
   the column block $b$.

   Parameters:
     a - Pointer to the CSC matrix.
     b - Column block index.
     c - Column index.
*/
#define CSC_FVAL(a,b,c) (a)->valtab[(a)->cscftab[b].coltab[c]]
/*
  Macro: CSC_VALNBR

  Compute the Number of element on the matrix.

  Parameters:
    a - Pointer to the CSC matrix.
*/
#define CSC_VALNBR(a)   (a)->cscftab[(a)->cscfnbr\
				     -1].coltab[(a)->cscftab[(a)->cscfnbr \
							     -1].colnbr]

/* Section: Structures */
/*
   Structure: CscFormat_

   Internal block column structure.

   Contains:
     colnbr - Number of columns in the block column.
     coltab - Array of indexes of the start of each column in
	      the row and value arrays.
*/
struct CscFormat_ {
  PASTIX_INT   colnbr;
  PASTIX_INT * coltab;
};

/*
   Type: CscFormat

   See <CscFormat_> structure.
*/
typedef struct CscFormat_ CscFormat;

/*
  Structure: CscMatrix_

  Internal column block distributed CSC matrix.

  Contains:
    cscfnbr - Number of column block.
    cscftab - Array of Block column structures. (<CscFormat>)
    rowtab  - Array of rows in the matrix.
    valtab  - Array of values of the matrix.
    type    - 'S' for symmetric, 'H' for hermitian, U for unsymmetric.
*/
struct CscMatrix_ {
  PASTIX_INT         cscfnbr;
  CscFormat * cscftab;
  PASTIX_INT       * rowtab;
  PASTIX_FLOAT     * valtab;
  char         type;
};
/*
  Type: CscMatrix

  See <CscMatrix_> structure.
*/
typedef struct CscMatrix_ CscMatrix;
#endif /* CSC_H */
