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
  File: olafread.c
  
  Driver for the olaf matrix format.
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
#include "olafread.h"

/*
  Function: olafReadHeader
  
  Reads header from file *infile*

  header format is :
  > Nrow
  > Nnzero

  Nrow is equal to Ncol

  Parameters:
    infile - File to read from
    Nrow   - Number of rows	
    Ncol   - Number of columns	
    Nnzero - Number of non zeros
    Type   - Type of the matrix 
 */
void olafReadHeader(FILE         *infile,
		    pastix_int_t *Nrow, 
		    pastix_int_t *Ncol, 
		    pastix_int_t *Nnzero, 
		    char         *Type)
{
  long temp1;
  if (1 != fscanf(infile, "%ld\n", &temp1))
    {
      fprintf(stderr, "ERROR: Reading matrix header\n");
      exit(1);
    }
  *Nrow = (pastix_int_t)temp1;
  if (1 != fscanf(infile, "%ld\n", &temp1))
    {
      fprintf(stderr, "ERROR: Reading matrix header\n");
      exit(1);
    }
  *Nnzero = (pastix_int_t)temp1;
  *Ncol = *Nrow;
  Type[0] = 'R';
  Type[1] = 'S';
  Type[2] = 'A';
  Type[3] = '\0';
}

/*
  Function: olafRead

  Reads a matrix in olaf format.

  Header format is described in <olafReadHeader>, 
  Olaf files contains :
  colptr, row and avals in the CSC format 
  > colptr[0]
  > colptr[1]
  >....
  > row[0]
  > row[1]
  > ...
  > avals[0]
  > avals[1]
  > ...
  

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
void olafRead(char const      *filename, 
	      pastix_int_t    *Nrow, 
	      pastix_int_t    *Ncol, 
	      pastix_int_t    *Nnzero, 
	      pastix_int_t   **col, 
	      pastix_int_t   **row, 
	      pastix_float_t **val, 
	      char           **Type, 
	      char           **RhsType, 
	      pastix_float_t **rhs)
{
  FILE *infile;
  pastix_int_t iter,size;
  long temp1;
  double temp2;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));
  (*RhsType)[0] = 'A';
  (*RhsType)[1] = 'A';
  (*RhsType)[2] = 'A';

#ifdef TYPE_COMPLEX 
  fprintf(stderr, "\nWARNING: This drivers reads non complex matrices, imaginary part will be 0\n\n");
#endif

  infile = fopen(filename, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "olafcsr");
      EXIT(MOD_SI,FILE_ERR);
    }
  olafReadHeader(infile, Nrow, Ncol, Nnzero, *Type);

  printf("Nrow %ld Ncol %ld Nnzero %ld\n", (long)*Nrow, (long)*Ncol, (long)*Nnzero);
  
  (*col) = (pastix_int_t *) malloc((*Ncol+1)*sizeof(pastix_int_t));
  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  (*rhs) = (pastix_float_t *) malloc((*Ncol)*sizeof(pastix_float_t));
 
  if (((*col) == NULL) || ((*row) == NULL) || ((*val) == NULL) || ((*rhs) == NULL))
    fprintf(stderr, "olafRead : Not enough memory for \n");

  for (iter=0; iter<(*Ncol+1); iter++)
    {
      if (1 != fscanf(infile, "%ld", &temp1))
	{
	  fprintf(stderr, "ERROR: Reading matrix header\n");
	  exit(1);
	}

      (*col)[iter] = (pastix_int_t)temp1;
    }

  size=*Nnzero;
  
  for (iter=0; iter<size; iter++)
    {
      if (1 != fscanf(infile, "%ld", &temp1))
	{
	  fprintf(stderr, "ERROR: Reading matrix header\n");
	  exit(1);
	}

      (*row)[iter] = (pastix_int_t)temp1;
    }

  for (iter=0; iter<size; iter++)
    {
      if (1 != fscanf(infile, "%lf", &temp2))
	{
	  fprintf(stderr, "ERROR: Reading matrix header\n");
	  exit(1);
	}

      (*val)[iter] = (pastix_float_t)temp2;
    }

  fclose(infile);

  infile = fopen("olafrhs", "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "olafrhs");
      EXIT(MOD_SI,FILE_ERR);
    }

  for (iter=0; iter<(*Ncol); iter++)
    {
      if (1 != fscanf(infile, "%lf", &temp2))
	{
	  fprintf(stderr, "ERROR: Reading matrix header\n");
	  exit(1);
	}

      (*rhs)[iter] = (pastix_float_t)temp2;
    }

  fclose(infile);
}
