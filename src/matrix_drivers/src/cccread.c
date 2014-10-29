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
   File: cccread.c
   
   Reads file in ccc format.

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
#include "cccread.h"

/*
  Function: cccReadHeader

  Reads header from a file in ccc matrix format.
  Nrow is equal to Ncol.
  Type is "CSA" if SYMPART is defined (default) "CUA" otherwise.

  File format is like :
  > Ncol Nnzero

  
  Parameters:
    infile - File to read from.
    Nrow   - Number of rows
    Ncol   - Number of columns
    Nnzero - Number of non zeros
    Type   - Type of the matrix
 */
void cccReadHeader(FILE         *infile, 
		   pastix_int_t *Nrow, 
		   pastix_int_t *Ncol, 
		   pastix_int_t *Nnzero, 
		   char         *Type)
{
  long temp1,temp2;
  if (2 != fscanf(infile, "%ld %ld\n", &temp1, &temp2))
    {
      fprintf(stderr, "ERROR: Reading matrix header\n");
      exit(1);
    }
  *Nrow=(pastix_int_t)temp1;
  *Nnzero=(pastix_int_t)temp2;
  *Ncol = *Nrow;
  Type[0] = 'C';
  Type[1] = 'U';
  Type[2] = 'A';
  Type[3] = '\0';
#ifdef SYMPART
  Type[1] = 'S';
  *Nnzero=(*Nnzero-*Ncol)/2+*Ncol;
#endif
}


/* 
   Function: cccRead
   
   Reads Matrix in ccc format.

   Header format is described in <cccReadHeader>, 
   "filename"/hfile contains columns
   Enf of the matrix is in three files in CSC format :
   "filename"/ifile contains Ncol columns,
   "filename"/jfile contains Nnzeros rows,
   "filename"/afile contains Nnzeros values.

   Parameters:
     filename - Path to the directory containing hfile, ifile, jfile and afile
     Nrow     - Number of rows
     Ncol     - Number of columns
     Nnzero   - Number of non zeros
     col      - Index of first element of each column in *row* and *val*
     row      -	Row of eah element				       
     val      -	Value of each element				       
     Type     -	Type of the matrix				       
     RhsType  -	Type of the right hand side.			       
*/
void cccRead(char const      *filename, 
	     pastix_int_t    *Nrow, 
	     pastix_int_t    *Ncol, 
	     pastix_int_t    *Nnzero, 
	     pastix_int_t   **col, 
	     pastix_int_t   **row, 
	     pastix_float_t **val, 
	     char           **Type, 
	     char           **RhsType)
{
  FILE *infile,*infile1,*infile2;
  pastix_int_t iter,size,i=0,ii=0;
  double temp1,temp2;
  char   filename2[STR_SIZE];

  temp1=0;
  temp2=0;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));
  (*RhsType)[0] = '\0';

#ifndef TYPE_COMPLEX 
  fprintf(stderr, "\nWARNING: This drivers reads complex matrices, imaginary part will be dropped\n\n");
#endif

  sprintf(filename2, "%s/hfile", filename);
  infile = fopen(filename2, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "hfile");
      EXIT(MOD_SI,FILE_ERR);
    }
  cccReadHeader(infile, Nrow, Ncol, Nnzero, *Type);
  fclose(infile);

  printf("Nrow %ld Ncol %ld Nnzero %ld\n", (long)*Nrow, (long)*Ncol, (long)*Nnzero);
  
  (*col) = (pastix_int_t *) malloc((*Ncol+1)*sizeof(pastix_int_t));
  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
 
  if (((*col) == NULL) || ((*row) == NULL) || ((*val) == NULL))
    fprintf(stderr, "cccRead : Not enough memory for \n");

  sprintf(filename2, "%s/ifile", filename);
  infile = fopen(filename2, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "ifile");
      EXIT(MOD_SI,FILE_ERR);
    }
  for (iter=0; iter<(*Ncol+1); iter++)
    {
      long temp;
      if (1 != fscanf(infile, "%ld", &temp))
	{
	  fprintf(stderr, "ERROR: Reading Matrix\n");
	  exit(1);
	}
      (*col)[iter]=(pastix_int_t)temp;
    }
  fclose(infile);

#ifdef SYMPART
  size=2*(*Nnzero-*Ncol)+*Ncol;
#else
  size=*Nnzero;
#endif
  sprintf(filename2, "%s/jfile", filename);
  infile1 = fopen(filename2, "r");
  if (infile1==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "jfile");
      EXIT(MOD_SI,FILE_ERR);
    }
  sprintf(filename2, "%s/afile", filename);
  infile2 = fopen(filename2, "r");
  if (infile2==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "afile");
      EXIT(MOD_SI,FILE_ERR);
    }

  
  for (iter=0; iter<size; iter++)
    {
      long x;
      if (1 != fscanf(infile1, "%ld", &x))
	{
	  fprintf(stderr, "ERROR: Reading Matrix\n");
	  exit(1);
	}

#ifdef SYMPART
      if (iter+1>=(*col)[ii+1])
	{
	  (*col)[ii+1]=i+1;
	  ii++;
	}
      if ((pastix_int_t)x>=ii+1) 
	{
	  (*row)[i] = (pastix_int_t)x;
	}
#else
      (*row)[iter] = (pastix_int_t)x;
#endif
      if (2 != fscanf(infile2, "%lf %lf", &temp1, &temp2))
	{
	  fprintf(stderr, "ERROR: Reading Matrix\n");
	  exit(1);
	}

#ifdef SYMPART
      if ((pastix_int_t)x>=ii+1)
	{
#if (defined X_ARCHalpha_compaq_osf1)
#ifdef TYPE_COMPLEX
	  (*val)[i] = PASTIX_FLOAT(temp1,temp2);
#else
	  (*val)[i] =(pastix_float_t)temp1;
#endif
#else
	  (*val)[i] = (pastix_float_t)temp1;
	  (*val)[i] += ((pastix_float_t)temp2)*I;
#endif
	  i++;
	}
#else
      (*val)[iter] = (pastix_float_t)temp1;
      (*val)[iter] += ((pastix_float_t)temp2)*I;
#endif
    }
  fclose(infile1);
  fclose(infile2);

#ifdef SYMPART
  (*col)[*Ncol]=i+1;
  ASSERT(i==*Nnzero,MOD_SI);
#endif
}
