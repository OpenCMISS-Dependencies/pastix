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
  File: fdupread.c

  Interface for the fortran driver writen in driver_fdupros.f
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
#include "fdupread.h"
#ifdef FDUPROS
/*
  Function: driverFdupros

  C driver from Fabrice Dupros.
  Binary format.

*/
void driverFdupros(char const      *filename,
		   pastix_int_t    *Nrow,
		   pastix_int_t    *Ncol,
		   pastix_int_t    *Nnzero,
		   pastix_int_t   **col,
		   pastix_int_t   **row,
		   pastix_float_t **val,
		   pastix_float_t **rhs,
		   char           **Type,
		   char           **RhsType)
{

  double        *tmpval;
  pastix_int_t   i;
  pastix_int_t   index;
  pastix_int_t   index0;
  int            tmpint1=0;
  int            tmpint2=0;
  FILE          *fp;
  char          *filename2;
  int           *tmpcoord;
  pastix_int_t   onerow;
  pastix_int_t   onecol;
  pastix_float_t oneval;

  *Type    = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));

  (*Type)[0]='R';
  (*Type)[1]='S';
  (*Type)[2]='A';
  (*Type)[3]='\0';
    
  (*RhsType)[0]='\0';
  filename2 = (char *)malloc((strlen(filename)+20)*sizeof(char));
  sprintf(filename2, "%s/MATA01", filename);
  fp = fopen(filename2, "rb");
  fread(&tmpint1, sizeof(int), 1, fp);
  fread(&tmpint2, sizeof(int), 1, fp);
  printf("MATA - %d %d\n", tmpint1, tmpint2);
  tmpval   = (double*)malloc(tmpint2*sizeof(double));

  fread(tmpval, sizeof(double), tmpint2, fp);
  fclose(fp);

  *Ncol   = (pastix_int_t)tmpint1;
  *Nnzero = (pastix_int_t)tmpint2;
  *Nrow   = *Ncol;


  sprintf(filename2, "%s/COORD01", filename);
  fp = fopen(filename2, "rb");
  fread(&tmpint1, sizeof(int), 1, fp);
  fread(&tmpint2, sizeof(int), 1, fp);
  tmpcoord = (int*) malloc(2*(*Nnzero)*sizeof(int));
  printf("COORD - %d %d\n", tmpint1, tmpint2);
  
  fread(tmpcoord, sizeof(int), 2*(*Nnzero), fp);
  
  fclose(fp);

  if (NULL == ((*col) = (pastix_int_t *)malloc(sizeof(pastix_int_t)*((*Nrow)+1))))
    MALLOC_ERROR("col");;
  if (NULL == ((*row) = (pastix_int_t *)malloc(sizeof(pastix_int_t)*(*Nnzero))))
    MALLOC_ERROR("row");
  if (NULL == ((*val) = (pastix_float_t *)malloc(sizeof(pastix_float_t)*(*Nnzero))))
    MALLOC_ERROR("val");

  for ( i = 0; i < *Ncol+1; i++)
    {
      (*col)[i] = 0;
    }
  for ( i = 0; i < *Nnzero; i++)
    {
      (*col)[tmpcoord[i]-1]++;
    }

  index = 1;
  for ( i = 0; i < *Ncol+1; i++)
    {
      index0 = (*col)[i];
      (*col)[i] = index;
      index = index + index0;
    }
  for ( i = 0; i < *Nnzero; i++)
    {
      onecol = tmpcoord[i];
      onerow = tmpcoord[i+(*Nnzero)];
      oneval = (pastix_float_t)(tmpval[i]);
      index = (*col)[onecol-1]-1;
      if (index >= *Nnzero)
	{
	  fprintf(stdout, "onecol %ld Ncol %ld index %ld, *Nnzero %ld\n", onecol, *Ncol, index, *Nnzero);
	  exit(1);
	}
      (*val)[index] = oneval;
      (*row)[index] = onerow;
      (*col)[onecol-1] = index+1+1;
    }

  for ( i = *Ncol -1 ; i > 0 ; i--)
    {
      (*col)[i+1] = (*col)[i];
    }
  (*col)[0] = 1;

  free(tmpcoord);
  printf("FDuprosDriver: Nrow=%ld Ncol=%ld Nnzero=%ld\n",
	 (long)*Nrow,(long)*Ncol,(long)*Nnzero);

  if (NULL == ((*rhs) =
	       (pastix_float_t *)malloc(sizeof(pastix_float_t)*(*Nrow))))
    MALLOC_ERROR("rhs");
  if (*Nnzero < *Nrow)
    {
      free(tmpval);
      if (NULL == ((tmpval) =
		   (double *)malloc(sizeof(double)*(*Nrow))));
    }

  sprintf(filename2, "%s/ZRHS01", filename);
  fp = fopen(filename2, "rb");

  fread(&tmpint1, sizeof(int), 1, fp);

  fread(tmpval, sizeof(double), tmpint1, fp);

  for (i = 0; i < *Ncol; i++)
    (*rhs)[i] = (pastix_float_t)(tmpval[i]);

  free(tmpval);

}

/*
  Function: driverFdupros_dist

  C driver from Fabrice Dupros.
  Binary format, distributed.

*/
void driverFdupros_dist(char const      *filename,
			pastix_int_t    *Nrow,
			pastix_int_t    *Ncol,
			pastix_int_t    *Nnzero,
			pastix_int_t   **col,
			pastix_int_t   **row,
			pastix_int_t   **loc2glob,
			pastix_float_t **val,
			pastix_float_t **rhs,
			char           **Type,
			char           **RhsType,
			MPI_Comm         pastix_comm)
{

  double        *tmpval;
  pastix_int_t   i;
  pastix_int_t   index;
  int            tmpint1=0;
  int            tmpint2=0;
  int            tmpint3=0;
  int            tmpint4=0;
  FILE          *fp;
  char          *filename2;
  int           *tmpcoord;
  pastix_int_t   onerow;
  pastix_int_t   onecol;
  pastix_float_t oneval;
  int            rank;

  *Type    = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));

  (*Type)[0]='R';
  (*Type)[1]='S';
  (*Type)[2]='A';
  (*Type)[3]='\0';

  (*RhsType)[0]='\0';
  //	Lecture Matrice
  MPI_Comm_rank(pastix_comm, &rank);
  rank = rank +1;

  filename2 = (char *)malloc((strlen(filename)+20)*sizeof(char));
  sprintf(filename2, "%s/MATA%4.4d", filename, rank);
  fp = fopen(filename2, "rb");
  fread(&tmpint1, sizeof(int), 1, fp);
  fread(&tmpint2, sizeof(int), 1, fp);
  fread(&tmpint3, sizeof(int), 1, fp);
  fread(&tmpint4, sizeof(int), 1, fp);
  
  tmpval   = (double*)malloc(tmpint4*sizeof(double));

  fread(tmpval, sizeof(double), tmpint4, fp);
  fclose(fp);

  *Ncol   = (pastix_int_t)tmpint3;
  *Nnzero = (pastix_int_t)tmpint4;
  *Nrow   = *Ncol;


  sprintf(filename2, "%s/COORD%4.4d", filename, rank);
  fp = fopen(filename2, "rb");
  fread(&tmpint1, sizeof(int), 1, fp);
  fread(&tmpint2, sizeof(int), 1, fp);
  fread(&tmpint3, sizeof(int), 1, fp);
  fread(&tmpint4, sizeof(int), 1, fp);
  tmpcoord = (int*) malloc(2*(*Nnzero)*sizeof(int));

  fread(tmpcoord, sizeof(int), 2*(*Nnzero), fp);

  fclose(fp);

  if (NULL == ((*col) =
	       (pastix_int_t *)malloc(sizeof(pastix_int_t)*((*Nrow)+1))))
    MALLOC_ERROR("col");
  if (NULL == ((*row) = (pastix_int_t *)malloc(sizeof(pastix_int_t)*(*Nnzero))))
    MALLOC_ERROR("row");
  if (NULL == ((*loc2glob) =
	       (pastix_int_t *)malloc(sizeof(pastix_int_t)*(*Nrow))))
    MALLOC_ERROR("loc2glob");
  if (NULL == ((*val) =
	       (pastix_float_t *)malloc(sizeof(pastix_float_t)*(*Nnzero))))
    MALLOC_ERROR("val");

  for ( i = 0; i < *Ncol+1; i++)
    {
      (*col)[i] = 0;
    }

  index = 0;
  (*col)[index] = 1;
  (*loc2glob)[index] = tmpcoord[(*Nnzero)];

  tmpcoord[(*Nnzero)] = index+1;

  for ( i = 1; i < *Nnzero; i++)
    {
      if ((*loc2glob)[index] != tmpcoord[(*Nnzero)+i])
	{
	  index++;
	  (*col)[index] = i+1;
	  (*loc2glob)[index] = tmpcoord[(*Nnzero)+i];
	}
      tmpcoord[(*Nnzero)+i] = index+1;
    }
  (*col)[*Ncol] = *Nnzero;

  for ( i = 0; i < *Nnzero; i++)
    {
      onecol = tmpcoord[(*Nnzero)+i];
      onerow = tmpcoord[i];
      oneval = (pastix_float_t)(tmpval[i]);
      if (onecol > *Ncol)
	{
	  fprintf(stdout, "%ld: %ld %ld %ld\n", 
		  (long)rank, (long)i, (long)onecol, (long)(*Ncol));
	  exit(1);
	}
      index = (*col)[onecol-1]-1;
      if (index >= *Nnzero)
	{
	  fprintf(stdout, "%ld: onecol %ld index %ld, *Nnzero %ld\n", 
		  (long)rank, (long)onecol, (long)index, (long)(*Nnzero));
	  exit(1);
	}

      (*val)[index] = oneval;
      (*row)[index] = onerow;
      (*col)[onecol-1] ++;
    }
  for ( i = *Ncol -1 ; i >= 0 ; i--)
    {
      (*col)[i+1] = (*col)[i];
    }
  (*col)[0] = 1;


  free(tmpcoord);

  if (NULL == ((*rhs) =
	       (pastix_float_t *)malloc(sizeof(pastix_float_t)*(*Nrow))))
    MALLOC_ERROR("rhs");

  if (*Nnzero < *Nrow)
    {
      free(tmpval);
      if (NULL == ((tmpval) =
		   (double *)malloc(sizeof(double)*(*Nrow))));
    }
  
  sprintf(filename2, "%s/ZRHS%4.4d", filename, rank);
  fp = fopen(filename2, "rb");

  fread(&tmpint1, sizeof(int), 1, fp);

  fread(tmpval, sizeof(double), tmpint1, fp);

  for (i = 0; i < *Ncol; i++)
    (*rhs)[i] = (pastix_float_t)(tmpval[i]);

  free(tmpval);
  
}
#endif 
