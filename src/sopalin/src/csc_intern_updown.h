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
  File: csc_intern_updown.h

  Build UpDownVector from user right-hand-side and CSCd.
  Retrieve soltion from UpDownVector.
  Construct UpDownVector such as X[i] = 1, or X[i] = i.

*/
#ifndef CSC_INTERN_UPDOWN_H
#define CSC_INTERN_UPDOWN_H

/*
  Function: CscdUpdownRhs

  Fill-in UpDownVector structure from user right-hand-side member.

  Parameters:
    updovct - UpDownVector structure to fill-in.
    symbmtx - Solver matrix.
    rhs     - Right-hand-side member.
    perm    - reverse permutation tabular.
    dof      - Number of degree of freedom.
 */
void CscUpdownRhs(UpDownVector       *updovct,
		  const SolverMatrix *symbmtx, 
		  const PASTIX_FLOAT        *rhs, 
		  const PASTIX_INT          *perm,
		  int                 dof);

/*
  Function: CscdUpdownRhs

  Fill-in UpDownVector structure from user distributed right-hand-side member.

  Parameters:
    updovct - UpDownVector structure to fill-in.
    symbmtx - Solver matrix.
    rhs     - Right-hand-side member.
    invp    - reverse permutation tabular.
    g2l     - local numbers of global nodes, if not local contains -owner
    dof      - Number of degree of freedom.
 */
void CscdUpdownRhs(UpDownVector       *updovct,
		   const SolverMatrix *symbmtx, 
		   const PASTIX_FLOAT        *rhs, 
		   const PASTIX_INT          *invp,
		   const PASTIX_INT          *g2l,
		   const PASTIX_INT           ln,
		   int                 dof);

/*
  Function:CscdRhsUpdown

  Builds solution from UpDownVector structure

  Parameters:
    updovct  - UpDownVector structure containing the solution.
    symbmtx  - Solver matrix structure.
    rhs      - Solution to fill.
    ncol     - Number of columns in local matrix.
    dof      - Number of degree of freedom.
    comm     - MPI communicator.
  
 */
void CscRhsUpdown(const UpDownVector *updovct, 
		  const SolverMatrix *symbmtx, 
		  PASTIX_FLOAT              *rhs, 
		  const PASTIX_INT           ncol,
		  const PASTIX_INT          *invp,
		  const int           dof, 
		  const int           rhsmaking, 
		  MPI_Comm            comm);

/*
  Function:CscdRhsUpdown

  Builds distributed solution from
  UpDownVector structure

  Parameters:
    updovct  - UpDownVector structure containing the solution.
    symbmtx  - Solver matrix structure.
    x        - Solution to fill.
    ncol     - Number of columns in local matrix.
    g2l      - local numbers of global nodes, if not local contains -owner
    ord      - ordering
    dof      - Number of degree of freedom.
    comm     - MPI communicator.
  
 */
void CscdRhsUpdown(const UpDownVector *updovct, 
		   const SolverMatrix *symbmtx, 
		   PASTIX_FLOAT              *x,
		   const PASTIX_INT           ncol, 
		   const PASTIX_INT          *g2l,
		   const PASTIX_INT          *invp,
		   int                 dof,  
		   MPI_Comm            comm);

/*
  Function: Csc2updown

  Fill-in UpDownVector structure such as the solution of
  the system Ax=b is x_i=1 (API_RHS_1) or x_i=i (API_RHS_I). 

  Parameters:
    cscmtx   - internal CSCd matrix.
    updovct  - UpDownVector structure to fill-in.
    symbmtx  - Solver matrix.
    mode     - wanted solution API_RHS_1 or API_RHS_I.
    comm     - MPI communicator.
*/
void Csc2updown(const CscMatrix    *cscmtx, 
		UpDownVector       *updovct,
		const SolverMatrix *symbmtx, 
		int                 mode,
		MPI_Comm            comm);


/*
  Function: Csc2updown_X0

  Fill-in initial X0 for reffinement if we don't want to use
  Solve step.

  (iparm[IPARM_ONLY_RAFF] == API_YES)

  Parameters:
    updovct - UpDownVector structure were to copy B as the first X0 used for raffinement.
    symbmtx - Solver matrix.
    mode    - Rule to construct X0 (API_RHS_1 : X0[i] = 1, API_RHS_I : X0[i] = i).
    comm    - MPI_Communicator.
*/
void Csc2updown_X0(UpDownVector *updovct, 
		   /*const*/ SolverMatrix *symbmtx, 
		   int mode, 
		   MPI_Comm comm);

#endif /* CSC_INTERN_UPDOWN_H */
