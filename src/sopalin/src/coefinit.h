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
  File: coefinit.h

  Allocation and initialisation of the coeficient of the solver matrix.  
*/
#ifndef COEFINIT_H
#define COEFINIT_H

/* Section: Functions declarations*/

/*
  Function: CoefMatrix_Allocate
  
  Allocate matrix coefficients in coeftab and ucoeftab.

  Should be first called with me = -1 to allocated coeftab.
  Then, should be called with me set to thread ID 
  to allocate column blocks coefficients arrays.
  
  Parameters
 
     datacode  - solverMatrix 
     factotype - factorization type (LU, LLT ou LDLT)
     me        - thread number. (-1 for first call, 
                 from main thread. >=0 to allocate column blocks 
		 assigned to each thread.)
 
*/
void CoefMatrix_Allocate (SopalinParam    *sopar,
			  SolverMatrix    *datacode,
			  pthread_mutex_t *mutex,
			  PASTIX_INT              factotype, 
			  PASTIX_INT              me);

/*
  Function: CoefMatrix_Init

  Init coeftab and ucoeftab coefficients.

  Parameters:
     datacode     - solverMatrix 
     barrier      - Barrier used for thread synchronisation.
     me           - Thread ID 
     iparm        - Integer parameters array.
     transcsc     - vecteur transcsc
     sopalin_data - <Sopalin_Data_t> structure for NUMA version.
*/
void CoefMatrix_Init     (SolverMatrix         *datacode, 
			  sopthread_barrier_t  *barrier, 
			  PASTIX_INT                   me,
			  PASTIX_INT                  *iparm, 
			  PASTIX_FLOAT               **transcsc, 
			  Sopalin_Data_t       *sopalin_data);

/*
  Function: CoefMatrix_Free
  
  Free the solver matrix coefficient tabular : coeftab and ucoeftab.
  
  WARNING: Call it with one unnique thread. 

  Parameters:
    datacode   - solverMatrix 
    factotype  - factorisation type (<API_FACT>)
    
*/  
void CoefMatrix_Free     (SopalinParam *sopar,
			  SolverMatrix *datacode, 
			  PASTIX_INT           factotype);


#endif /* COEFINIT_H */
