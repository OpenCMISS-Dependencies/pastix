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
  File: threetilesread.c

  Reads matrix from three files in IJV separated format.

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
#include "threefilesread.h"

/* 
   Function: threeFilesReadHeader
   
   Read header from three file IJV format.

   Header contains:
   > Nrow Ncol Nnzero
   or 
   > Ncol
   > Nnzero

   Parameters:
     infile - file to read header from
     Nrow   - Number of rows
     Ncol   - Number of columns
     Nnzero - Number of non zeros
     Type   - Type of the matrix (always "RUA")
 */
void threeFilesReadHeader(FILE         *infile, 
			  pastix_int_t *Nrow, 
			  pastix_int_t *Ncol, 
			  pastix_int_t *Nnzero, 
			  char         *Type)
{
  char line[BUFSIZ];
  long temp1,temp2,temp3;
 
  Type[0] = 'R';
  Type[1] = 'U';
  Type[2] = 'A';
  Type[3] = '\0';

  /* ncol nrow nnzero */
  FGETS(line,BUFSIZ,infile);

  sscanf(line, "%ld %ld %ld", &temp1,&temp2,&temp3);
  if (temp1!=temp2)
    {
      temp2=temp1;
      FGETS(line,BUFSIZ,infile);
      sscanf(line, "%ld", &temp3);
    }

  *Nrow = (pastix_int_t)temp1;
  *Ncol = (pastix_int_t)temp2;
  *Nnzero = (pastix_int_t)temp3;
}

/*
  Function: threeFilesRead

  Read matrix from three files IJV

  header file is "filename"/header
  columns file is "filename"/ia_threeFiles
  rows file is "filename"/ja_threeFiles
  values file is "filename"/ra_threeFiles

  header is describde in <threeFilesReadHeader>
  each other file contain one element by line.

  Parameters:
    dirname - Path to the directory containing matrix 
    Ncol    - Number of columns					 
    Nrow    - Number of rows						 
    Nnzero  - Number of non zeros					 
    col     - Index of first element of each column in *row* and *val*  
    row     - Row of eah element				       	 
    val     - Value of each element				       	 
    Type    - Type of the matrix				       	 
    RhsType - Type of the right-hand-side.			         
 */
void threeFilesRead(char const      *dirname, 
		    pastix_int_t    *Ncol, 
		    pastix_int_t    *Nrow, 
		    pastix_int_t    *Nnzero, 
		    pastix_int_t   **col, 
		    pastix_int_t   **row, 
		    pastix_float_t **val, 
		    char           **Type, 
		    char           **RhsType)
{

  FILE * iaFile;
  FILE * jaFile;
  FILE * raFile;
  FILE * headerFile;
  char * filename;
  pastix_int_t * tempcol;
  pastix_int_t * temprow;
  pastix_float_t * tempval;
  pastix_int_t iter,baseval;
  pastix_int_t tmp,total,pos,limit; 
  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(1*sizeof(char));
  (*RhsType)[0] = '\0';

  filename = malloc(strlen(dirname)+10);

#ifdef TYPE_COMPLEX 
  fprintf(stderr, "\nWARNING: This drivers reads non complex matrices, imaginary part will be 0\n\n");
#endif
  
  sprintf(filename,"%s/header",dirname);
  headerFile = fopen (filename,"r");
  if (headerFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }
  threeFilesReadHeader(headerFile,Nrow,Ncol,Nnzero,*Type);
  fclose (headerFile);

  sprintf(filename,"%s/ia_threeFiles",dirname); 
  iaFile = fopen(filename,"r");  
  if (iaFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }

  sprintf(filename,"%s/ja_threeFiles",dirname);
  jaFile = fopen(filename,"r");
  if (jaFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }
  sprintf(filename,"%s/ra_threeFiles",dirname);
  raFile = fopen(filename,"r");
  if (raFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }

  /* Allocation memoire */
  tempcol = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  temprow = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  tempval = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  
  if ((tempcol==NULL) || (temprow == NULL) || (tempval == NULL))
    {
      fprintf(stderr, "threeFilesRead : Not enough memory (Nnzero %ld)\n",(long)*Nnzero);
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  
  /* Remplissage */
  for (iter=0; iter<(*Nnzero); iter++)
    {
      long temp1,temp2;
      double tempv;
      if ( 1 != fscanf(iaFile,"%ld\n", &temp1))
	{
	  fprintf(stderr, "ERROR: reading matrix\n");
	  exit(1);
	}
      temprow[iter]=(pastix_int_t)temp1;
      if (1 != fscanf(jaFile,"%ld\n", &temp2))
	{
	  fprintf(stderr, "ERROR: reading matrix\n");
	  exit(1);
	}
      tempcol[iter]=(pastix_int_t)temp2;
      if (1 != fscanf(raFile,"%le\n", &tempv))
	{
	  fprintf(stderr, "ERROR: reading matrix\n");
	  exit(1);
	}
      tempval[iter]= (pastix_float_t)tempv;
    }
  
  fclose (iaFile);
  fclose (jaFile);
  fclose (raFile);

  (*col) = (pastix_int_t *) malloc((*Nrow+1)*sizeof(pastix_int_t));  
  memset(*col,0,(*Nrow+1)*sizeof(pastix_int_t));  
  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  memset(*row,0,(*Nnzero)*sizeof(pastix_int_t));
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  if (((*col)==NULL) || ((*row) == NULL) || ((*val) == NULL))
    {
      fprintf(stderr, "threeFilesRead : Not enough memory (Nnzero %ld)\n",(long)*Nnzero);
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  
  for (iter = 0; iter < (*Nnzero); iter ++)      
    {
      (*col)[tempcol[iter]-1]++;
    }

  baseval=1; /* Attention on base a 1 */
  total = baseval;
  
  for (iter = 0; iter < (*Ncol)+1; iter ++)      
    {
      tmp = (*col)[iter];
      (*col)[iter]=total;
      total+=tmp;
    }

  for (iter = 0; iter < (*Nnzero); iter ++)      
    {
      
      pos = (*col)[tempcol[iter]-1]-1;
      limit = (*col)[tempcol[iter]]-1;
      while((*row)[pos] != 0 && pos < limit)
	{
	  pos++;
	}
      if (pos == limit)
	fprintf(stderr, "Erreur de lecture\n");
      
      (*row)[pos] = temprow[iter];
      (*val)[pos] = tempval[iter];
    }      
  
  memFree_null(tempval);
  memFree_null(temprow);
  memFree_null(tempcol);
}
