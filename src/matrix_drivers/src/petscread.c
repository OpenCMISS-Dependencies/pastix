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
 * File: petscread.c
 *
 * Read Matrix written by PETSc in binary format.
 *
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

/* must be before complex definition */
#include "mmio.h"

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

/* must be after complex definition */
#include "pastix.h"
#include "common_drivers.h"

static inline
int swap_indians_int(int num)
{
  int byte0, byte1, byte2, byte3;

  byte0 = (num & 0x000000FF) >> 0 ;
  byte1 = (num & 0x0000FF00) >> 8 ;
  byte2 = (num & 0x00FF0000) >> 16 ;
  byte3 = (num & 0xFF000000) >> 24 ;

  return((byte0 << 24) | (byte1 << 16) | (byte2 << 8) | (byte3 << 0));
}

static inline
pastix_float_t swap_indians_pastix_float_t(pastix_float_t d)
{
    size_t i;
    union
    {
	pastix_float_t value;
	char bytes[sizeof(pastix_float_t)];
    } in, out;
    in.value = d;
    for (i = 0; i < sizeof(pastix_float_t); i++)
	out.bytes[i] = in.bytes[sizeof(pastix_float_t)-1-i];

    return out.value;
}

static inline
int swap_indians_int2(int d)
{
    size_t i;
    union
    {
	int value;
	char bytes[sizeof(int)];
    } in, out;
    in.value = d;
    for (i = 0; i < sizeof(int); i++)
	out.bytes[i] = in.bytes[sizeof(int)-1-i];

    return out.value;
}

#define SWAP_INDIANS(i)                                                 \
  ( sizeof(int) == sizeof(i) )?                                         \
  ( (need_convert)?swap_indians_int(i):(i)):                            \
  ( (need_convert)?swap_indians_pastix_float_t(i):(i))
/*
 *  Function: PETScRead
 *
 *  Reads a matrix in a binary PETScformat/
 *
 * Parameters:
 *   dirname - Path to the directory containing matrix
 *   Ncol    - Number of columns
 *   Nrow    - Number of rows
 *   Nnzero  - Number of non zeros
 *   col     - Index of first element of each column in *row* and *val*
 *   row     - Row of eah element
 *   val     - Value of each element
 *   Type    - Type of the matrix
 *   RhsType - Type of the right-hand-side.
 *
 */
void PETScRead(char const      *filename,
	       pastix_int_t    *Ncol,
	       pastix_int_t    *Nrow,
	       pastix_int_t    *Nnzero,
	       pastix_int_t   **col,
	       pastix_int_t   **row,
	       pastix_float_t **val,
	       char           **Type,
	       char           **RhsType)
{
  char           *  buffer;
  int            *  intbuff;
  pastix_float_t *  floatbuff;
  unsigned long     fileLen;
  int               i, j, iter, rowsize, baseval;
  pastix_int_t      total;
  pastix_int_t      tmp;
  pastix_int_t      pos;
  pastix_int_t      limit;
  int               need_convert;
  int             * tempcol;
  int             * temprow;
  pastix_float_t  * tempval;
  FILE            * file = fopen(filename, "rb");
  int rc;

  if (file==NULL)
  {
    fprintf(stderr,"cannot load %s\n", filename);
    EXIT(MOD_SI,FILE_ERR);
  }
  fseek(file, 0, SEEK_END);
  fileLen = ftell(file);
  fseek(file, 0, SEEK_SET);

  buffer = (char *)malloc(fileLen+1);
  if (!buffer)
  {
    fprintf(stderr, "Error in PETScRead : Not enough memory\n");
    EXIT(MOD_SI,OUTOFMEMORY_ERR);
  }

  rc = fread(buffer, fileLen, 1, file);
  if (rc != (int)fileLen ) {
      perror("Error while reading file");
  }
  fclose(file);
  intbuff = (int*)buffer;
  i = (int)(*(intbuff++));
  if (i == 1211216) need_convert = 0;
  else need_convert = 1;
  i = SWAP_INDIANS(i);
  if (i != 1211216)
  {
    fprintf(stderr, "Error in PETScRead : Incorrect file header\n");
    EXIT(MOD_SI,FILE_RR);
  }


  *Nrow = (pastix_int_t)(SWAP_INDIANS(*(intbuff++)));
  *Ncol = (pastix_int_t)(SWAP_INDIANS(*(intbuff++)));
  *Nnzero =  (pastix_int_t)(SWAP_INDIANS(*(intbuff++)));
  fprintf(stdout, "%ld %ld %ld %ld\n", (long)i, (long)(*Ncol), (long)(*Nrow), (long)(*Nnzero));
  if (fileLen*sizeof(char) != (4+(*Nrow)+(*Nnzero))*sizeof(int)+(*Nnzero)*sizeof(pastix_float_t))
  {
    fprintf(stderr, "Error in PETScRead : Incorrect size of file (%ld != %ld) \n",
	    (long)(fileLen*sizeof(char)),
	    (long)(1*sizeof(char)+(3+(*Nrow)+(*Nnzero))*sizeof(int)+(*Nnzero)*sizeof(pastix_float_t)));
    EXIT(MOD_SI,FILE_RR);
  }
  if (*Nnzero == -1)
  {
    fprintf(stderr, "Error in PETScRead : Full matrices not supported\n");
    EXIT(MOD_SI,FILE_ERR);
  }
  if (*Nrow != *Ncol)
  {
    fprintf(stderr, "Error in PETScRead : Non square matrices not supported\n");
    EXIT(MOD_SI,FILE_ERR);
  }

  tempcol = (int *) malloc((*Nnzero)*sizeof(int));
  temprow = (int *) malloc((*Nnzero)*sizeof(int));
  tempval = (pastix_float_t  *) malloc((*Nnzero)*sizeof(pastix_float_t));

  {
    int * ttrow = temprow;

    for (i = 0; i < *Nrow; i++)
    {
      rowsize = SWAP_INDIANS(*(intbuff++));
      for (j = 0; j < rowsize; j++)
	*(ttrow++) = i;
    }
  }
  for (i = 0; i < *Nnzero; i++)
  {
    tempcol[i] = SWAP_INDIANS(*(intbuff++));
  }
  floatbuff = (pastix_float_t*)(intbuff);
  for (i = 0; i < *Nnzero; i++)
  {
    tempval[i] = SWAP_INDIANS(*(floatbuff++));
  }
  free(buffer);

  (*col) = (pastix_int_t *) malloc((*Nrow+1)*sizeof(pastix_int_t));
  memset(*col,0,(*Nrow+1)*sizeof(pastix_int_t));
  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  memset(*row,0,(*Nnzero)*sizeof(pastix_int_t));
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  if (((*col)==NULL) || ((*row) == NULL) || ((*val) == NULL))
  {
    fprintf(stderr, "petscread : Not enough memory (Nnzero %ld)\n",(long)*Nnzero);
    EXIT(MOD_SI,OUTOFMEMORY_ERR);
  }

  baseval=1; /* Attention on base a 1 */
  for(iter=0; iter<(*Nnzero); iter++)
  {
    tempcol[iter]+=baseval;
    temprow[iter]+=baseval;
  }

  for (iter = 0; iter < (*Nnzero); iter ++)
  {
    (*col)[tempcol[iter]-1]++;
  }


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
  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(1*sizeof(char));
  (*RhsType)[0] = '\0';
  (*Type)[0] = 'R';
  (*Type)[1] = 'U';
  (*Type)[2] = 'A';
  (*Type)[3] = '\0';

  return;
}
