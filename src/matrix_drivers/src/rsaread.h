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
  File rsaread.h

  Interface for the fortran driver writen in skitf.f
*/

#ifndef USE_NOFORTRAN
#if (defined X_ARCHpower_ibm_aix)
#define FORTRAN_CALL(nom) nom
#else
#define FORTRAN_CALL(nom) nom ## _
#endif
#else
#define FORTRAN_CALL(nom) 
#endif

/*
  Function: FORTRAN_CALL(wreadmtc)

  Declaration of the wreadmtc fortran function
  defined in skitf.f.

  Parameters:
     tmp1      - Maximum number of column                                    (INPUT)
     tmp2      - Maximum number of non zeros                                 (INPUT)
     tmp3      - job to be done (see skitf file)                             (INPUT)
     filename  - Path to the file to read from                               (INPUT)
     len       - length of *filname*                                         (INPUT)
     val       - Values of the elements of the matrix.			     (OUTPUT)
     row       - Rows of the elements of the matrix.			     (OUTPUT)
     col       - Index of first element of each column in *row* and *val*    (OUTPUT)
     crhs      - Right hand side(s).					     (OUTPUT)
     nrhs      - Number of right hand side(s).				     (OUTPUT)
     RhsType   - Right hand side type					     (OUTPUT)
     tmpNrow   - Number of rows.					     (OUTPUT)
     tmpNcol   - Number of columns.					     (OUTPUT)
     tmpNnzero - Number of non zeros.					     (OUTPUT)
     title     - name of the matrix.					     (OUTPUT)
     key       - key of the matrix (see skitf.f)			     (OUTPUT)
     Type      - Type of the matrix 					     (OUTPUT)
     ierr      - Error return value                                          (OUTPUT)

 */
void  FORTRAN_CALL(wreadmtc)(int        * tmp1,
			     int        * tmp2,
			     int        * tmp3,
			     const char * filename,
			     int        * len, 
			     double     * val, 
			     int        * row, 
			     int        * col, 
			     double     * crhs, 
			     int        * nrhs,
			     char       * RhsType, 
			     int        * tmpNrow, 
			     int        * tmpNcol, 
			     int        * tmpNnzero,
			     char       * title,
			     char       * key, 
			     char       * Type, 
			     int        * ierr);   
    
/*
  Function: rsaReadHeader

  Interface for <wreadmtc> to read header from file *filename*.

  Parameters:
    filename - Path to the file to read from
    Nrow     - Number of row
    Ncol     - Number of column
    Nnzero   - Number of non zeros
    Type     - Type of the matrix
    RhsType  - Type of the right hand side

*/
void rsaReadHeader(char const   *filename, 
		   pastix_int_t *Nrow, 
		   pastix_int_t *Ncol, 
		   pastix_int_t *Nnzero, 
		   char         *Type, 
		   char         *RhsType);

/*
  Function: rsaRead

  Read matrix using wreadmtc fortran driver defined in skitf.f.
  The return matrix is in CSC format, 
  Nrow is equal to Ncol or you should get an error.

  Parameters:
    filename - Path to the file to read from
    Nrow     - Number of rows in matrix
    Ncol     - Number of columns in matrix
    Nnzero   - Number of non zros in matrix
    col      - Index of first element of each column in *row* and *val*
    row      - Row of eah element
    val      - Value of each element
    Type     - Type of the matrix
    RhsType  - Type of the right hand side.
 */
void rsaRead(char const      *filename, 
	     pastix_int_t    *Nrow, 
	     pastix_int_t    *Ncol, 
	     pastix_int_t    *Nnzero, 
	     pastix_int_t   **col, 
	     pastix_int_t   **row, 
	     pastix_float_t **val, 
	     char           **Type, 
	     char           **RhsType);
