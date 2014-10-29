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
  File: tools.h

  headers for some tools used in sopalin.

 */
#ifndef TOOLS_H
#define TOOLS_H

/* 
   Function: dim_dgeam


   Computes b = alph * a + b.
 
     - if transa and transb equals 'N'.
        > b := alpha * a + b 
     - if transa = 'T' and transb ='N'.
        > b := alpha * trans(a) + b 
     - if transa = 'N' and transb ='T'.
        > trans(b) := alpha * a + trans(b) 
     - if transa = 'T' and transb ='T'.
        > trans(b) := alpha * trans(a) + trans(b) 

   Parameters: 
     transa - indicates if a needs to be transposed.
     transb - indicates if b needs to be transposed.
     m      - number of row in a and b. 
     n      - number of colonnes in a and b.
     alpha  - scalar.
     a      - Matrix a.
     lda    - Stride between 2 columns of a.
     b      - Matrix b.
     ldb    - Stride between 2 columns of b.
*/
void         dim_dgeam(char *transa, char *transb, PASTIX_INT m, PASTIX_INT n, PASTIX_FLOAT alpha,
		       PASTIX_FLOAT *a, PASTIX_INT lda, PASTIX_FLOAT *b, PASTIX_INT ldb);

/* 
   Function: GetMpiType

   Construct a MPI type to store complex values.

   Parameters:
     none

   Return:
   
   the constructed MPI type for complex floating values.
 */
MPI_Datatype GetMpiType(void);

/* 
   Function: FreeMpiType

   Free the MPI type to store complex values.

   Parameters:
     none

 */
void FreeMpiType(void);

/*
  Function: mysum
  
  Computes the sum of *in* and *inout* vectors 
  and stores it in *inout*.
  
  
  Parameters: 
    in    - First vector.
    inout - Second vector wich will store the sum.
    len   - Size of each vector.
    dptr  - MPI datatype
 */
void         mysum(void *in, void *inout, int *len, MPI_Datatype *dptr);
/*
  Function: GetMpiSum
  
  Creates a new MPI sum operation.

  Parameters:
    none
    
  Return: the new operation created.
  
 */
MPI_Op       GetMpiSum(void);

/*
  Function: FreeMpiSum
  
  Free the MPI sum operation.

  Parameters:
    none
    
 */
void FreeMpiSum(void);

#endif
