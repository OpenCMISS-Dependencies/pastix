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
  File: chbread.h

  Read a matrix in chb format.
 */



/*
  Function: chbReadHeader

  Reads header from a chb matrix file.

  header format is:
  > title 73-80 Key
  > 1-14 totcrd 15-28 ptrcrd 29-42 indcrd 43-56 valcrd 57-70 rhscrd
  > 1-3 mxtype 15-28 nrow 29-42 ncol 43-56 nnzero 57-70 neltvl 
  > 1-16 ptrfmt 17-32 indfmt 33-52 valfmt 53-72 rhsfmt
  > 1  2 rhstyp 3  15-28 nrhs 29-42 nrhsix

  Parameters
    infile  - File to read from
    Type    - Type of the matrix
    Nrow    - Number of row in the matrix
    Ncol    - Number of columns in the matrix
    Nnzero  - Number of non zeros in the matrix
    Nrhs    - Number of right-hand-side terms
    Ptrfmt  - 
    Indfmt  -
    Valfmt  -
    Rhsfmt  -
    Ptrcrd  -
    Indcrd  -
    Valcrd  -
    Rhscrd  -
    RhsType - Type of right-hand-side term(s)

 */
void chbReadHeader(FILE         *infile, 
		   char         *Type, 
		   pastix_int_t *Nrow, 
		   pastix_int_t *Ncol, 
		   pastix_int_t *Nnzero, 
		   pastix_int_t *Nrhs, 
		   char         *Ptrfmt, 
		   char         *Indfmt, 
		   char         *Valfmt, 
		   char         *Rhsfmt, 
		   pastix_int_t *Ptrcrd, 
		   pastix_int_t *Indcrd, 
		   pastix_int_t *Valcrd, 
		   pastix_int_t *Rhscrd,  
		   char         *RhsType);
/*
  Function: chbRead

  Reads a matrix in chb format.

  Header is described in <chbReadHeader>
  Formats are sicribed in <chbParseRfmt> and <chbParseIfmt>

  In our file we have
  header
  valuesFormat
  rowFormat
  columnFormat
  (rhsFormat)
  then the columns,
  the rows, 
  the values,
  (the rhs)
  

  Parameters:
    filename - Path to the file to read from
    Nrow     - Number of rows
    Ncol     - Number of columns
    Nnzero   - Number of non zeros
    col      - Index of first element of each column in *row* and *val*
    row      - Row of eah element				       	
    val      - Value of each element				       	
    Type     - Type of the matrix				       	
    RhsType  - Type of the right-hand-side terms.			       	
    rhs      - right-hand-side term(s)


 */
void chbRead(char const      *filename, 
	     pastix_int_t    *Nrow, 
	     pastix_int_t    *Ncol, 
	     pastix_int_t    *Nnzero, 
	     pastix_int_t   **col, 
	     pastix_int_t   **row,
	     pastix_float_t **val, 
	     char           **Type, 
	     char           **RhsType, 
	     pastix_float_t **rhs);


/*
  Function: hbParseRfmt

  CHB float format parser
  
  Format is like :
  > (3(1P,E25.16)) 
  or
  > (1P,3E25.16)
  or 
  > (1P3E25.16)
  or 
  > (3E25.16)
  or 
  > (3E25)
  for perline = 3, format = E, width = 25 and prec = 16


  Parameters:
    fmt      - format to parse
    perline  - number of element per line 
    width    - 
    prec     - Precision
    flag     - 
*/
void chbParseRfmt(char *fmt, pastix_int_t *perline, pastix_int_t *width, pastix_int_t *prec, char *flag);

/*
  Function:chbParseIfmt

  CHB integer format parser

  format is :
  > (perlineIwidth) or (X,perlineIwidth)
  Parameters:
    fmt      - format to parse
    perline  - number of element per line 
    width    - 
    flag     - 
*/
void chbParseIfmt(char *fmt, pastix_int_t *perline, pastix_int_t *width, pastix_int_t *flag);


