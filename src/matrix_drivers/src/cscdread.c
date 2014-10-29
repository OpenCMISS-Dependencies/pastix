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
   File: cscdread.c

   Read files in cscd format.
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
#include "cscdread.h"

/*
  Function: cscdRead

  Reads a matrix in cscd format
  Nrow is equal to Ncol.

  File format is like :
  ...


  Parameters:
    dirname     - Path to the directory containing matrix
    colptr      - Index of first element of each column in *row* and *val*
    row         - Row of eah element
    loc2glb     - Correspondance between local and global numbering
    avals       - Value of each element
    rhs         - Right Hand Side
    colnbr      - Number of columns
    nnz         - Number of non-zeros
    pastix_comm - MPI communicator
 */
void
cscdRead(char const      *dirname,
	 pastix_int_t   **colptr,
	 pastix_int_t   **row,
	 pastix_int_t   **loc2glb,
	 pastix_float_t **avals,
	 pastix_float_t **rhs,
	 pastix_int_t    *colnbr,
	 pastix_int_t    *nnz,
	 MPI_Comm         pastix_comm)
{
  const pastix_int_t nbreltperline = 4; /* nbr of elt per line */
  FILE              *infile;
  char               line[BUFSIZ], file[BUFSIZ];
  int                myrank, nbproc, tmpint;
  long               tempint1,   tempint2,   tempint3,   tempint4;
  long double        tempfloat1, tempfloat2, tempfloat3, tempfloat4;
  pastix_int_t       vertloc, edgeloc;
  pastix_int_t       iterelt;
  pastix_int_t      *vectsize    = NULL;
  pastix_int_t      *vectsizercv = NULL;
  pastix_int_t       offset      = 1;
  char              *filename;
  int                i;
  (void)pastix_comm;

#ifdef TYPE_COMPLEX
  fprintf(stderr, "\nWARNING: This drivers reads non complex matrices, imaginary part will be 0\n\n");
#endif

  /* verifier que le nombre de fichier = nb de proc MPI */
  MPI_Comm_rank( pastix_comm, &myrank);
  MPI_Comm_size( pastix_comm, &nbproc);

  filename = (char*)malloc(sizeof(char)*(strlen(dirname)+40));
  sprintf(filename,"%s/main",dirname);

  if (myrank == 0)
    {
      infile = fopen(filename, "r");
      if (infile==NULL)
	{
	  fprintf(stderr,"cannot load %s\n", filename);
	  exit(EXIT_FAILURE);
	}
      FGETS(line, BUFSIZ, infile);
      sscanf(line, "%d", &tmpint); /* Read number of filename */
      fprintf(stdout, "Nombre de fichier %d\n", tmpint);
      fclose(infile);

      if (nbproc != tmpint)
	{
	  /* pour l'instant rien. au choix : recreer un nouveau comm mpi, refusionner la csc et la redecouper */
	  if (myrank == 0)
	    fprintf(stderr, "Veuillez fournir un communicateur MPI de %d processus\nActuellement, le communicateur contient %d processus\n", tmpint, nbproc);
	  exit(EXIT_FAILURE);
	}
    }

  infile = fopen(filename, "r");
  if (infile == NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      exit(EXIT_FAILURE);
    }
  FGETS(line, BUFSIZ, infile);

  for (i=0; i<=myrank; i++)
    {
      FGETS(line, BUFSIZ, infile);
      sscanf(line, "%s", file);
    }
  fclose(infile);

  sprintf(filename,"%s/%s",dirname,file);
  infile = fopen(filename, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"[P%d] cannot load %s\n", myrank, filename);
      exit(EXIT_FAILURE);
    }

  FGETS(line, BUFSIZ, infile);
  sscanf(line, "%ld %ld", &tempint1,&tempint2);
  vertloc = tempint2;
/*   fprintf(stderr, "[P%d] rowlocal %ld\n",myrank, (long) vertloc); */
  FGETS(line, BUFSIZ, infile);
  sscanf(line, "%ld", &tempint1);
  edgeloc = tempint1;
  *colnbr = vertloc;
  *nnz = edgeloc;
/*   fprintf(stderr, "[P%d] nzlocal %ld\n", myrank, (long) edgeloc); */

  *colptr      = (pastix_int_t *)   malloc((vertloc+1) * sizeof(pastix_int_t));
  if (   (*colptr)   == NULL)
    fprintf(stderr, "[P%d] cscdRead : Not enough memory for *colptr\n",myrank);
  *loc2glb = (pastix_int_t *)   malloc(  vertloc   * sizeof(pastix_int_t));
  if ( (*loc2glb) == NULL)
    fprintf(stderr, "[P%d] cscdRead : Not enough memory for *loc2glb\n",myrank);
  *row      = (pastix_int_t *)   malloc(  edgeloc   * sizeof(pastix_int_t));
  if (  (*row)   == NULL)
    fprintf(stderr, "[P%d] cscdRead : Not enough memory for *row\n",myrank);
  *avals   = (pastix_float_t *) malloc(  edgeloc   * sizeof(pastix_float_t));
  if (  (*avals) == NULL)
    fprintf(stderr, "[P%d] cscdRead : Not enough memory for *aval\n",myrank);
  memset(*avals, 0, (edgeloc)*sizeof(pastix_float_t));
  *rhs     = (pastix_float_t *) malloc(  vertloc   * sizeof(pastix_float_t));
  if (   (*rhs)  == NULL)
    fprintf(stderr, "[P%d] cscdRead : Not enough memory for *rhs\n",myrank);
  memset(*rhs, 0, (vertloc)*sizeof(pastix_float_t));

  /* Recuperation de Loc2glb*/
  for (iterelt=0; iterelt<vertloc+1-nbreltperline;iterelt+=nbreltperline )
    {
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld %ld %ld",
	     &tempint1, &tempint2, &tempint3, &tempint4);
      (*loc2glb)[iterelt]   = (pastix_int_t)tempint1;
      (*loc2glb)[iterelt+1] = (pastix_int_t)tempint2;
      (*loc2glb)[iterelt+2] = (pastix_int_t)tempint3;
      (*loc2glb)[iterelt+3] = (pastix_int_t)tempint4;
    }
  switch (vertloc-iterelt)
    {
    case 1:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%ld",&tempint1);
      (*loc2glb)[iterelt] += (pastix_int_t)tempint1;
      iterelt++;
      break;
    case 2:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld", &tempint1, &tempint2);
      (*loc2glb)[iterelt]   = (pastix_int_t)tempint1;
      (*loc2glb)[iterelt+1] = (pastix_int_t)tempint2;
      iterelt+=2;
      break;
    case 3:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld %ld", &tempint1, &tempint2, &tempint3);
      (*loc2glb)[iterelt]   = (pastix_int_t)tempint1;
      (*loc2glb)[iterelt+1] = (pastix_int_t)tempint2;
      (*loc2glb)[iterelt+2] = (pastix_int_t)tempint3;
      iterelt+=3;
      break;
    default:
      break;
    }

  /* Recuperation de Colptr dans un style tres particulier... */
  for (iterelt=0; iterelt<vertloc+1+1-nbreltperline;iterelt+=nbreltperline )
    {
      FGETS(line,BUFSIZ,infile);
      if (4 != sscanf(line,"%ld %ld %ld %ld", &tempint1, &tempint2, &tempint3, &tempint4))
	{
	  fprintf(stderr, "ERROR: reading colptr\n");
	  exit(1);
	}
      (*colptr)[iterelt]   = (pastix_int_t)tempint1;
      (*colptr)[iterelt+1] = (pastix_int_t)tempint2;
      (*colptr)[iterelt+2] = (pastix_int_t)tempint3;
      (*colptr)[iterelt+3] = (pastix_int_t)tempint4;
    }

  switch (vertloc-iterelt+1)
    {
    case 1:
      FGETS(line,BUFSIZ,infile);
      if (1 != sscanf(line,"%ld",&tempint1))
	{
	  fprintf(stderr, "ERROR: reading colptr\n");
	  exit(1);
	}
      (*colptr)[iterelt] += (pastix_int_t)tempint1;
      iterelt++;
      break;
    case 2:
      FGETS(line,BUFSIZ,infile);
      if (2 != sscanf(line,"%ld %ld", &tempint1, &tempint2))
	{
	  fprintf(stderr, "ERROR: reading colptr\n");
	  exit(1);
	}
      (*colptr)[iterelt]   = (pastix_int_t)tempint1;
      (*colptr)[iterelt+1] = (pastix_int_t)tempint2;
      iterelt+=2;
      break;
    case 3:
      FGETS(line,BUFSIZ,infile);
      if (3 != sscanf(line,"%ld %ld %ld", &tempint1, &tempint2, &tempint3))
	{
	  fprintf(stderr, "ERROR: reading colptr\n");
	  exit(1);
	}
      (*colptr)[iterelt]   = (pastix_int_t)tempint1;
      (*colptr)[iterelt+1] = (pastix_int_t)tempint2;
      (*colptr)[iterelt+2] = (pastix_int_t)tempint3;
      iterelt+=3;
      break;
    default:
      break;
    }
  fprintf(stdout, "iterelt %ld, vertloc %ld\n", (long)iterelt, (long)vertloc);
  vectsize    = malloc(sizeof(pastix_int_t)*nbproc);
  vectsize    = memset(vectsize, 0, sizeof(pastix_int_t)*nbproc);
  vectsizercv = malloc(sizeof(pastix_int_t)*nbproc);
  vectsizercv = memset(vectsizercv, 0, sizeof(pastix_int_t)*nbproc);

  vectsize[myrank] = vertloc;

  if (vectsize    == NULL) fprintf(stderr, "[P%d] Erreur : alloc vectsize\n",    myrank);
  if (vectsizercv == NULL) fprintf(stderr, "[P%d] Erreur : alloc vectsizercv\n", myrank);

  MPI_Allreduce(vectsize, vectsizercv, nbproc, MPI_PASTIX_INT, MPI_SUM, pastix_comm);
  for (i=0; i<myrank; i++)
    offset += vectsizercv[i];

  /* Recuperation de ROW*/
  for (iterelt=0; iterelt<edgeloc+1-nbreltperline; iterelt+=4)
    {
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld %ld %ld", &tempint1,&tempint2,&tempint3,&tempint4);
      (*row)[iterelt]   = (pastix_int_t)tempint1;
      (*row)[iterelt+1] = (pastix_int_t)tempint2;
      (*row)[iterelt+2] = (pastix_int_t)tempint3;
      (*row)[iterelt+3] = (pastix_int_t)tempint4;

    }
  switch (edgeloc-iterelt)
    {
    case 1:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%ld",&tempint1);
      (*row)[iterelt] = (pastix_int_t)tempint1;
      iterelt++;
      break;
    case 2:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld", &tempint1, &tempint2);
      (*row)[iterelt]   = (pastix_int_t)tempint1;
      (*row)[iterelt+1] = (pastix_int_t)tempint2;
      iterelt+=2;
      break;
    case 3:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld %ld", &tempint1, &tempint2, &tempint3);
      (*row)[iterelt]   = (pastix_int_t)tempint1;
      (*row)[iterelt+1] = (pastix_int_t)tempint2;
      (*row)[iterelt+2] = (pastix_int_t)tempint3;
      iterelt+=3;
      break;
    }

  /* Recuperation de Aval */
  for (iterelt=0; iterelt<edgeloc+1-nbreltperline; iterelt+=4)
    {
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%Lf %Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3,&tempfloat4);
      (*avals)[iterelt]   = (pastix_float_t)tempfloat1;
      (*avals)[iterelt+1] = (pastix_float_t)tempfloat2;
      (*avals)[iterelt+2] = (pastix_float_t)tempfloat3;
      (*avals)[iterelt+3] = (pastix_float_t)tempfloat4;
    }

  switch (edgeloc-iterelt)
    {
    case 1:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%Lf",&tempfloat1);
      (*avals)[iterelt] = (pastix_float_t)tempfloat1;
      iterelt++;
      break;
    case 2:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%Lf %Lf",&tempfloat1,&tempfloat2);
      (*avals)[iterelt]   = (pastix_float_t)tempfloat1;
      (*avals)[iterelt+1] = (pastix_float_t)tempfloat2;
      iterelt+=2;
      break;
    case 3:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3);
      (*avals)[iterelt]   = (pastix_float_t)tempfloat1;
      (*avals)[iterelt+1] = (pastix_float_t)tempfloat2;
      (*avals)[iterelt+2] = (pastix_float_t)tempfloat3;
      iterelt+=3;
      break;
    }

  /* Recuperation de Rhs, second membre... */
  for (iterelt=0; iterelt<vertloc+1-nbreltperline; iterelt+=4)
    {
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%Lf %Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3,&tempfloat4);
      (*rhs)[iterelt]   = (pastix_float_t)tempfloat1;
      (*rhs)[iterelt+1] = (pastix_float_t)tempfloat2;
      (*rhs)[iterelt+2] = (pastix_float_t)tempfloat3;
      (*rhs)[iterelt+3] = (pastix_float_t)tempfloat4;
    }

  switch (vertloc-iterelt)
    {
    case 1:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%Lf",&tempfloat1);
      (*rhs)[iterelt] = (pastix_float_t)tempfloat1;
      iterelt++;
      break;
    case 2:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%Lf %Lf",&tempfloat1,&tempfloat2);
      (*rhs)[iterelt]   = (pastix_float_t)tempfloat1;
      (*rhs)[iterelt+1] = (pastix_float_t)tempfloat2;
      iterelt++;
      break;
    case 3:
      FGETS(line,BUFSIZ,infile);
      sscanf(line,"%Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3);
      (*rhs)[iterelt]   = (pastix_float_t)tempfloat1;
      (*rhs)[iterelt+1] = (pastix_float_t)tempfloat2;
      (*rhs)[iterelt+2] = (pastix_float_t)tempfloat3;
      iterelt++;
      break;
    default:
      break;
    }

  fclose(infile);
  free(filename);
  free(vectsize);
  free(vectsizercv);
}
