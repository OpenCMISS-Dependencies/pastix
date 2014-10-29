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
/*
  File: bordi.h

  This is the test module for incomplete 
  block ordering strategies.            

  Authors:
    Xavier LACOSTE - lacoste@labri.fr

*/

#ifndef BORDI_H
#define BORDI_H

/*
  Function: orderSplit
  
  Subdivide column blocks in column clocks of size blocksize.

  WARNING: unused 

  Parameters: 
    ordeptr   - Ordering.
    blocksize - size of block wanted.
*/
void orderSplit (Order * const ordeptr, 
		 PASTIX_INT           blocksize);

#if (defined SCOTCH_SEQSCOTCH || defined SCOTCH_H || defined SCOTCH_PTSCOTCH || defined PTSCOTCH_H)
/* 
   Function: orderSplit2

   Subdivide block columns in blocks of minimum size *bsmin*.
   Size depends on *rho* parameter too.

   WARNING: unused

   Parameters:
     ordeptr - Ordering.
     grphptr - Graph corresponding to the matrix.
     rho     - Parameter to compute new bloc size.
     bsmin   - Minimal bloc size.

 */
void orderSplit2 (Order        * const ordeptr, 
		  SCOTCH_Graph * const grphptr, 
		  double               rho, 
		  PASTIX_INT                  bsmin);

/*
  Function: orderSplit3
  
  Splits column blocs such has there is one column block in front of each bloc.

  Parameters:
    ordeptr  - Ordering.
    grphptr  - Graph associated to the matrix.
    matrsymb - Symbol matrix.

 */
void orderSplit3 (Order        * const ordeptr, 
		  SCOTCH_Graph * const grphptr, 
		  SymbolMatrix * const matrsymb);

#endif /* SCOTCH */

/*
  Function: symbolSplit
  
  Splits the symbol matrix.

  Parameters:
    matrsymb - symbol matrix
 */
void symbolSplit (SymbolMatrix * matrsymb);

/*
  Function: symbolRustine

  DESCRIPTION TO FILL

  Parameters:
    matrsymb  - Symbol matrix
    matrsymb2 - Symbol matrix
 */
void
symbolRustine (SymbolMatrix *       matrsymb, 
	       SymbolMatrix * const matrsymb2);


/* *********************************************
   Functions: These are the cost functions.

   This routine computes the factorization
   and solving cost of the given symbolic
   block matrix, whose nodes hold the number
   of DOFs given by the proper DOF structure.
   To ensure maximum accuracy and minimum loss
   of precision, costs are summed-up recursively.
   It returns:
   - 0   : on success.
   - !0  : on error.
 * *********************************************/
/* 
   Function: symbolCostn 

   TO FILL 

   Parameters:
     symbptr - Symbolic matrix to evaluate              
     deofptr - DOF structure associated with the matrix 
     keeptab - Flag array for blocks to keep            
     typeval - Type of cost computation                 
     nnzptr  - Size of the structure, to be filled      
     opcptr  - Operation count, to be filled            


 */
int symbolCostn (const SymbolMatrix * const  symbptr,
		 const Dof * const           deofptr,
		 const unsigned char * const          keeptab,
		 const SymbolCostType        typeval,
		 double * const              nnzptr, 
		 double * const              opcptr);


#if (defined SCOTCH_SEQSCOTCH || defined SCOTCH_H || defined SCOTCH_PTSCOTCH || defined PTSCOTCH_H)
/* **************************************
   Functions: This is the main function. 
 ****************************************/
/*
  Function: bordi

  This is the main function
  TO FILL

  Parameters:
    alpha    -
    symbptr  -
    graphptr -
    orderptr -

 */
void bordi(int            alpha, 
	   SymbolMatrix * symbptr, 
	   SCOTCH_Graph * graphptr, 
	   Order        * orderptr);
#endif /* SCOTCH */
#endif /* BORDI_H */
