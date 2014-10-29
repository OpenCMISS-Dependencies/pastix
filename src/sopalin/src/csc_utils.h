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
#ifndef CSC_UTILS_H
#define CSC_UTILS_H

/*
  Function: csc_symgraph

  
  Modify the CSC to a symetric graph one.
  Don't use it on a lower symetric CSC 
  it would give you all the CSC upper + lower.
  
  External function

  Parameters: 
    n     - Number of columns/vertices
    ia	  - Starting index of each column in *ja* and *a*
    ja	  - Row index of each element
    a 	  - Value of each element,can be NULL    
    newn  - New number of column
    newia - Starting index of each column in *ja* and *a* 
    newja - Row index of each element
    newa  - Value of each element,can be NULL    

 */
int csc_symgraph(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, 
		 PASTIX_INT *newn, PASTIX_INT **newia, PASTIX_INT **newja, PASTIX_FLOAT **newa);



/*
  Function: csc_symgraph_int

  
  Modify the CSC to a symetric graph one.
  Don't use it on a lower symetric CSC 
  it would give you all the CSC upper + lower.
  
  Parameters: 
    n           - Number of columns/vertices
    ia	        - Starting index of each column in *ja* and *a*
    ja	        - Row index of each element
    a 	        - Value of each element,can be NULL    
    newn        - New number of column
    newia       - Starting index of each column in *ja* and *a* 
    newja       - Row index of each element
    newa        - Value of each element,can be NULL    
    malloc_flag - flag to indicate if function call is intern to pastix or extern.
 */
int csc_symgraph_int (PASTIX_INT n,     PASTIX_INT * ia,    PASTIX_INT * ja,    PASTIX_FLOAT * a, 
		      PASTIX_INT *newn, PASTIX_INT **newia, PASTIX_INT **newja, PASTIX_FLOAT **newa, 
		      int malloc_flag);



/** 
    Function: csc_noDiag
    
    Supress diagonal term.              
    After this call, *ja* can be reallocated to *ia[n] -1*.
    
    Parameters:
      n  - size of the matrix.
      ia - Index in *ja* and *a* of the first element of each column
      ja - row of each element
      a  - value of each element, can be set to NULL

    Returns:
      ia and ja tabulars modified.
*/
void csc_noDiag(PASTIX_INT baseval, PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a);

/*
  Function: csc_check_doubles
  
  Check if the csc contains doubles and if correct if asked

  Assumes that the CSC is sorted.

  Assumes that the CSC is Fortran numeroted (base 1)

  Parameters:
    n      - Size of the matrix.
    colptr - Index in *rows* and *values* of the first element of each column
    rows   - row of each element
    values - value of each element
    dof    - Number of degrees of freedom
    flag   - Indicate if user wants correction (<API_BOOLEAN>)
    flagalloc - indicate if allocation on CSC uses internal malloc. 

    
  Returns:
    API_YES - If the matrix contained no double or was successfully corrected.
    API_NO  - Otherwise.
*/
int csc_check_doubles(PASTIX_INT      n,
		      PASTIX_INT   *  colptr,
		      PASTIX_INT   ** rows,
		      PASTIX_FLOAT ** values, 
		      int      dof,
		      int      flag,
		      int      flagalloc);

/*
  Function: csc_checksym

    Check if the CSC graph is symetric.
    
    For all local column C, 
    
    For all row R in the column C,
    
    We look in column R if we have the row number C.
       
    If we can correct we had missing non zeros.
    
    Assumes that the CSC is Fortran numbered (1 based).
    
    Assumes that the matrix is sorted.

  Parameters:
    n        - Number of local columns
    colptr   - Starting index of each columns in *ja*
    rows     - Row of each element.
    values   - Value of each element.
    correct  - Flag indicating if we can correct the symmetry.
    alloc    - indicate if allocation on CSC uses internal malloc. 
    dof      - Number of degrees of freedom.
*/
int csc_checksym(PASTIX_INT      n, 
		 PASTIX_INT     *colptr, 
		 PASTIX_INT    **rows, 
		 PASTIX_FLOAT  **values, 
		 int      correct,
		 int      alloc,
		 int      dof);

void CSC_colPerm(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_INT *cperm);
void CSC_colScale(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_FLOAT *dcol);
void CSC_rowScale(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_FLOAT *drow);

void CSC_sort(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_INT ndof);
void CSC_Fnum2Cnum(PASTIX_INT *ja, PASTIX_INT *ia, PASTIX_INT n);
void CSC_Cnum2Fnum(PASTIX_INT *ja, PASTIX_INT *ia, PASTIX_INT n);

/*
  Function: CSC_buildZerosAndNonZerosGraphs
  
  Separate a graph in two graphs, following 
  wether the diagonal term of a column is null or not.

  Parameters:
    n, colptr, rows, values  - The initial CSC
    n_nz, colptr_nz, rows_nz - The graph of the non-null diagonal part.
    n_z, colptr_z, rows_z    - The graph of the null diagonal part.
    perm                     - Permutation to go from the first graph to 
                               the one composed of the two graph concatenated.
    revperm                  - Reverse permutation tabular.
    criteria                 - Value beside which a number is said null.
*/
int CSC_buildZerosAndNonZerosGraphs(PASTIX_INT     n,
				    PASTIX_INT    *colptr,
				    PASTIX_INT    *rows,
				    PASTIX_FLOAT  *values,
				    PASTIX_INT    *n_nz,
				    PASTIX_INT   **colptr_nz,
				    PASTIX_INT   **rows_nz,
				    PASTIX_INT    *n_z,
				    PASTIX_INT   **colptr_z,
				    PASTIX_INT   **rows_z,
				    PASTIX_INT    *perm, 
				    PASTIX_INT    *revperm,
				    double  criteria);

/*
  Function: CSC_isolate

  Isolate a list of unknowns at the end of the CSC.

  Parameters:
    n            - Number of columns.
    colptr       - Index of first element of each column in *ia*.
    rows         - Rows of each non zeros.	    
    n_isolate    - Number of unknow to isolate.
    isolate_list - List of unknown to isolate.
*/
int CSC_isolate(PASTIX_INT     n,
		PASTIX_INT    *colptr,
		PASTIX_INT    *rows,
		PASTIX_INT     n_isolate,
		PASTIX_INT    *isolate_list,
		PASTIX_INT    *perm,
		PASTIX_INT    *revperm);


/*
  Function: csc_save

  Save a csc on disk.

  Parameters:
    n       - number of columns
    colptr  - First cscd starting index of each column in *ja* and *a*
    rows    - Row of each element in first CSCD
    values  - value of each cscd in first CSCD (can be NULL)
    dof     - Number of degrees of freedom
    outfile - Output stream.

  Return:
    NO_ERR
  
*/
int csc_save(PASTIX_INT      n,
	     PASTIX_INT    * colptr,
	     PASTIX_INT    * rows,
	     PASTIX_FLOAT  * values,
	     int      dof,
	     FILE   * outfile);
/*
  Function: csc_load

  Load a csc from disk.

  Fill *n*, *colptr*, *rows*, *values* and *dof* from *infile*.

  Parameters:
    n       - number of columns
    colptr  - First cscd starting index of each column in *ja* and *a*
    rows    - Row of each element in first CSCD
    values  - value of each cscd in first CSCD (can be NULL)
    dof     - Number of degrees of freedom
    outfile - Output stream.

  Return:
    NO_ERR
  
*/
int csc_load(PASTIX_INT    *  n,
	     PASTIX_INT    ** colptr,
	     PASTIX_INT    ** rows,
	     PASTIX_FLOAT  ** values,
	     int    *  dof,
	     FILE   *  infile);

#endif /* CSC_UTILS_H */
