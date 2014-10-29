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
  File: rsaread.c

  Interface for the fortran driver writen in skitf.f
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#ifdef FORCE_NOMPI
#include "pastix_nompi.h"
#else
#include <mpi.h>
#endif


#ifdef TYPE_COMPLEX
#if (defined X_ARCHalpha_compaq_osf1)
#ifndef USE_CXX
#ifndef   _RWSTD_HEADER_REQUIRES_HPP
#include <complex>
#else  /* _RWSTD_HEADER_REQUIRES_HPP */
#include <complex.hpp>
#endif /* _RWSTD_HEADER_REQUIRES_HPP */
#endif /* USE_CXX */
#else  /* X_ARCHalpha_compaq_osf1 */
#include <complex.h>
#endif /* X_ARCHalpha_compaq_osf1 */
#endif /* TYPE_COMPLEX */

#ifdef X_ARCHsun
#include <inttypes.h>
#endif

#include "pastix.h"
#include "common_drivers.h"
#include "rsaread.h"

    
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
		   char         *RhsType) 
{
#ifdef USE_NOFORTRAN
  fprintf(stderr, "ERROR: rsaread is a fortran driver.\n\tPlease recompile without -DUSE_NOFORTRAN.\n");
  return;
#else
  int     tmp;
  int    *col=NULL;
  int    *row=NULL;
  char    title[72+1];
  char    key[8+1];
  int     nrhs;
  int     len;
  int     ierr;
  double *val=NULL;
  double *crhs=NULL;
  int     tmpNrow;
  int     tmpNcol;
  int     tmpNnzero;

  len=strlen(filename);
  tmp=0;
  FORTRAN_CALL(wreadmtc)
    (&tmp,&tmp,&tmp,filename,&len,val,row,col,crhs,&nrhs,
     RhsType,&tmpNrow,&tmpNcol,&tmpNnzero,title,key,Type,&ierr);
  if(ierr != 0) {
    fprintf(stderr, "cannot read matrix (job=0)\n");
  }
  *Nrow=(pastix_int_t)tmpNrow;
  *Ncol=(pastix_int_t)tmpNcol;
  *Nnzero=(pastix_int_t)tmpNnzero;
  Type[3]='\0';

  /*ASSERT(*Nrow==*Ncol,MOD_SI);*/
  if ((*Nrow==*Ncol) == 0)
    {
      fprintf(stderr,"ERROR : (*Nrow!=*Ncol)\n");
      exit(EXIT_FAILURE);
    }
#endif /* USE_NOFORTRAN */
}

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
	     char           **RhsType)
{
#ifdef USE_NOFORTRAN
  fprintf(stderr, "ERROR: rsaread is a fortran driver.\n\tPlease recompile without -DUSE_NOFORTRAN.\n");
  return;
#else
  int     i;
  int     tmp;
  char    title[72+1];
  char    key[8+1];
  int     nrhs;
  int     len;
  int     ierr;
  double *crhs=NULL;
  int     tmpNrow;
  int     tmpNcol;
  int     tmpNnzero;
  int    *tmpcol;
  int    *tmprow;
  double *tmpval;
  int     base;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));

  
#ifdef TYPE_COMPLEX 
  fprintf(stderr, "\nWARNING: This drivers reads non complex matrices, imaginary part will be 0\n\n");
#endif

  rsaReadHeader(filename, Nrow, Ncol, Nnzero, *Type, *RhsType);
  printf("RSA: Nrow=%ld Ncol=%ld Nnzero=%ld\n",(long)*Nrow,(long)*Ncol,(long)*Nnzero);
  tmpNrow=(int)*Nrow;
  tmpNcol=(int)*Ncol;
  tmpNnzero=(int)*Nnzero;

  len=strlen(filename);
  *col=(pastix_int_t*)malloc(((*Nrow)+1)*sizeof(pastix_int_t));
  ASSERT(*col!=NULL,MOD_SI);
  tmpcol=(int*)malloc((tmpNrow+1)*sizeof(int));
  ASSERT(tmpcol!=NULL,MOD_SI);
  *row=(pastix_int_t*)malloc((*Nnzero)*sizeof(pastix_int_t));

  ASSERT(*row!=NULL,MOD_SI);
  tmprow=(int*)malloc(tmpNnzero*sizeof(int));
  ASSERT(tmprow!=NULL,MOD_SI);

  *val=(pastix_float_t*)malloc((*Nnzero)*sizeof(pastix_float_t));
  tmpval=(double*)malloc((*Nnzero)*sizeof(double));

  ASSERT(*val!=NULL,MOD_SI);
  tmp=2;
  nrhs=0;

  FORTRAN_CALL(wreadmtc)
    (&tmpNrow,&tmpNnzero,&tmp,filename,&len,tmpval,tmprow,tmpcol,crhs,
     &nrhs,*RhsType,&tmpNrow,&tmpNcol,&tmpNnzero,title,key,*Type,&ierr);

  (*RhsType)[0]='\0';
  myupcase(*Type);
  if(ierr != 0) {
    fprintf(stderr, "cannot read matrix (job=2)\n");
  }

  if (tmpcol[0] == 0)
    base = 1;
  else
    base = 0;
    
  for (i=0;i<tmpNrow+1;i++) (*col)[i]=(pastix_int_t)(tmpcol[i] + base);
  for (i=0;i<tmpNnzero;i++) (*row)[i]=(pastix_int_t)(tmprow[i] + base);
#ifdef TYPE_COMPLEX
  for (i=0;i<tmpNnzero;i++) (*val)[i]=(pastix_float_t)(tmpval[i] + I*0);
#else
  for (i=0;i<tmpNnzero;i++) (*val)[i]=(pastix_float_t)(tmpval[i]);
#endif
  memFree_null(tmpcol);
  memFree_null(tmprow);
  memFree_null(tmpval);
  *Nrow=(pastix_int_t)tmpNrow;
  *Ncol=(pastix_int_t)tmpNcol;
  *Nnzero=(pastix_int_t)tmpNnzero;
#endif /* USE_NOFORTRAN */
}
