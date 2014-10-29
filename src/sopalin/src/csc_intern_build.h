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
  File: csc_intern_build.h
  
  Functions to build internal CSCd from user CSCd.

  Function to free internal CSCd.

*/
#ifndef CSC_INTERN_BUILD_H
#define CSC_INTERN_BUILD_H

/*
  Function: CscOrdistrib

  Fill in *thecsc* CSC matrix in column block representation.

  Parameters: 
  
  thecsc     - Matrix in block column CSC format to fill in.
  Type       - 3 charactères for matrix type : only Type[1] is used to check if matrix is Symetric(S) or not(U).
  transcsc   - Transpose of the CSC in non symetric mode.
  ord        - ordering
  Nrow       - Number of rows.
  Ncol       - Number of columns.
  Nnzero     - Number of non zeros in the matrix.
  colptr     - Index in *rowind* and *val* of the start of each column.
  rowind     - Index of the elements.
  val        - values of the elements.
  forcetrans - If matrix symetric, transcsc will be the copy of the CSC_VALTAB.
  symbmtx    - Solver matrix
  procnum    - MPI process number.
  dof        - Number of degree of freedom
*/
void CscOrdistrib(CscMatrix          *thecsc, 
		  char               *Type, 
		  PASTIX_FLOAT             **transcsc,
		  const Order        *ord, 
		  PASTIX_INT                 Nrow, 
		  PASTIX_INT                 Ncol,
		  PASTIX_INT                 Nnzero, 
		  PASTIX_INT                *colptr, 
		  PASTIX_INT                *rowind, 
		  PASTIX_FLOAT              *val, 
		  PASTIX_INT                 forcetrans,
		  const SolverMatrix *symbmtx, 
		  PASTIX_INT                 procnum, 
		  PASTIX_INT                 dof);

/*
  Function: CscdOrdistrib

  Fill in *thecsc* CSC matrix in column block representation.

  - Construct cachetab (sizeof(PASTIX_INT)*globalNbCol) which will contain
  the column block wich will own each column (internal numerotation), 
  or -1 if not local 

  - Build newcoltab (sizeof(PASTIX_INT)*globalNbCol) which will contain the 
  coltab corresponding to the local internal CSCd.
  This CSCd correspond to the given CSCd adding upper part in Symmetric matrix.
  Also count number of triples (i,j,v) to send to each other processors.

  - Send the information about how many triples will be sent
  
  - Fill-in the arrays containing triples to send and send them.

  - Receive those arrays and correct the newcoltab arrays with information 
  from others processors.

  - Build CSC_COLNBR from symbolic matrix informations and CSC_COL from newcoltab.

  - Construct transpose matrix, in symmetric mode, transcsc == CSC_VALTAB; in 
  unsymmetric mode, allocate trowtab (number of total local elements) , 
  and build trscltb which contains number of elements, 
  in each column of each column bloc.

  - fill-in internal CSC row and values from local given CSCd, 
  also fill-in trowtab and transcsc in unsymmetric mode.
  CSC_COL and trscltb are incremented for each element added. 

  - fill-in  internal CSC row and values from iniformation received,
  also fill in transposed CSCd in unsymmetric mode.
  CSC_COL and trscltb are incremented for each element added.

  - restore CSC_COL.
  
  - sort internal CSCd.

  - sort intranal transposed CSCd.

  Parameters: 
  
  thecsc     - Matrix in block column CSC format to fill in.
  Type       - 3 charactères for matrix type : only Type[1] is used to check if matrix is Symetric(S) or not(U).
  transcsc   - Transpose of the CSC in non symetric mode.
  ord        - ordering
  Ncol       - Number of columns.
  colptr     - Index in *rowind* and *val* of the start of each column.
  rowind     - Index of the elements.
  val        - values of the elements.
  l2g        - global numbers of local nodes.
  gNcol      - global number of columns.
  g2l        - local numbers of global nodes, if not local contains -owner
  forcetrans - If matrix symetric, transcsc will be the copy of the CSC_VALTAB.
  symbmtx    - Solver matrix
  procnum    - MPI process number.
  dof        - Number of degree of freedom
  comm       - MPI communicator.
*/
void CscdOrdistrib(CscMatrix          *thecsc, 
		   char               *Type, 
		   PASTIX_FLOAT             **transcsc,
		   const Order        *ord, 
		   PASTIX_INT                 Ncol,
		   PASTIX_INT                *colptr, 
		   PASTIX_INT                *rowind, 
		   PASTIX_FLOAT              *val, 
		   PASTIX_INT                *l2g,
		   PASTIX_INT                 gNcol,
		   PASTIX_INT                *g2l,
		   PASTIX_INT                 forcetrans,
		   const SolverMatrix *symbmtx, 
		   PASTIX_INT                 procnum,
		   PASTIX_INT                 dof,
		   MPI_Comm            comm);

/* 
   Function: CscExit
   
   Free the internal CSCd structure.
   
   Parameters:
     thecsc - Internal CSCd to free.
*/
void CscExit(CscMatrix *thecsc);
#endif /* CSC_INTERN_BUILD_H */
