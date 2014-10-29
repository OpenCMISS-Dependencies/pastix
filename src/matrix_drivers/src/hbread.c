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
  File: hbread.c

  Interface to the Harwell-Boeing driver in C (iohb.c)
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
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
#include "hbread.h"
#include "iohb.h"
/*
  Function: HBRead
   
  Interface to the Harwell-Boeing driver in C (iohb.c)

  Parameters
    filename - Path to the file to read from
    Nrow     - Number of rows
    Ncol     - Number of columns
    Nnzero   - Number of non zeros
    col      - Index of first element of each column in *row* and *val*
    row      - Row of eah element				       
    val      - Value of each element				       
    Type     - Type of the matrix				       
    RhsType  - Type of the right hand side.			       
 */
void HBRead(char const      *filename, 
	    pastix_int_t    *Nrow, 
	    pastix_int_t    *Ncol, 
	    pastix_int_t    *Nnzero, 
	    pastix_int_t   **col, 
	    pastix_int_t   **row, 
	    pastix_float_t **val, 
	    char           **Type, 
	    char           **RhsType)
{
  int      i;
  int      nrhs;
  int      tmpNrow;
  int      tmpNcol;
  int      tmpNnzero;
  int     *tmpcol;
  int     *tmprow;
  int      Nrow2; 
  int      Ncol2;
  int      Nnzero2;
  double  *tmpval;
  int      ierr;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));

  readHB_info(filename, &Nrow2, &Ncol2, &Nnzero2, Type, &nrhs);

  *Nrow = Nrow2;
  *Ncol = Ncol2;
  *Nnzero = Nnzero2;

/*   fprintf(stderr,"Matrix in file %s is %ld x %ld, with %ld nonzeros with type %s;\n", */
/* 	  filename, (long)*Nrow, (long)*Ncol, (long)*Nnzero, *Type); */
/*   fprintf(stderr,"%d right-hand-side(s) available.\n",nrhs); */

/*   printf("RSA: Nrow=%ld Ncol=%ld Nnzero=%ld\n",(long)*Nrow,(long)*Ncol,(long)*Nnzero); */
#ifdef TYPE_COMPLEX
  fprintf(stderr,"Warning: HBRead is a real matrix driver, imaginary part will be 0\n");
  exit(EXIT_FAILURE);
#endif

  tmpNrow=(int)*Nrow;
  tmpNcol=(int)*Ncol;
  tmpNnzero=(int)*Nnzero;


  *col=(pastix_int_t*)malloc((*Nrow+1)*sizeof(pastix_int_t));
  ASSERT(*col!=NULL,MOD_SI);
  tmpcol=(int*)malloc((tmpNrow+1)*sizeof(int));
  ASSERT(tmpcol!=NULL,MOD_SI);
  *row=(pastix_int_t*)malloc(*Nnzero*sizeof(pastix_int_t));
  ASSERT(*row!=NULL,MOD_SI);
  tmprow=(int*)malloc(tmpNnzero*sizeof(int));
  ASSERT(tmprow!=NULL,MOD_SI);
  *val=(pastix_float_t*)malloc(*Nnzero*sizeof(pastix_float_t));
  ASSERT(*val!=NULL,MOD_SI);

  nrhs=0;
#if (defined PREC_DOUBLE && !defined TYPE_COMPLEX)
  tmpval = *val;
#else
  tmpval = (double*)malloc(*Nnzero*sizeof(double));
#endif

  ierr = readHB_mat_double(filename, tmpcol, tmprow, tmpval);
  if(ierr == 0) {
    fprintf(stderr, "cannot read matrix (job=2)\n");
  }

#if (!defined PREC_DOUBLE || defined TYPE_COMPLEX)
  for (i = 0; i < *Nnzero; i++)
    (*val)[i] = (pastix_float_t)tmpval[i];
#endif
  (*RhsType)[0]='\0';
  myupcase(*Type);
  for (i=0;i<tmpNrow+1;i++) (*col)[i]=(pastix_int_t)(tmpcol[i]);
  for (i=0;i<tmpNnzero;i++) (*row)[i]=(pastix_int_t)(tmprow[i]);
  memFree_null(tmpcol);
  memFree_null(tmprow);
  *Nrow=(pastix_int_t)tmpNrow;
  *Ncol=(pastix_int_t)tmpNcol;
  *Nnzero=(pastix_int_t)tmpNnzero;
}
