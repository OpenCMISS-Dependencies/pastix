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
  File: csc_intern_solve.h

  Functions to copy internal CSCd data 
  onto solver matrix coeftab.
*/

#ifndef CSC_INTERN_SOLVE_H
#define CSC_INTERN_SOLVE_H

/*
  Function: Csc2solv_cblk

  Copy the part of the internal CSCd corresponding to
  the column bloc itercblk into the SolverMatrix structure 
  coeftab which will be used to compute the decomposition.

  Used in NUMA mode.
 
  Parameters:
    cscmtx   - The internal CSCd matrix.
    datacode - The SolverMatrix structure used during decomposition.
    trandcsc - The internal CSCd transpose used in LU decomposition.
    itercblk - Column bloc number in which we had the internal CSCd. 
  
  
*/
void Csc2solv_cblk(const CscMatrix *cscmtx, 
		   SolverMatrix    *solvmtx, 
		   PASTIX_FLOAT           *trandcsc, 
		   PASTIX_INT              itercblk);

#endif /* CSC_INTERN_SOLVE_H */
