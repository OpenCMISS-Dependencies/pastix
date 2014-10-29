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
#include "common_pastix.h"
#include "errors.h"

#define NB_ELEMENTS 10000
int main(int argc, char ** argv)
{
  int i,j;
  PASTIX_INT   * index;
  PASTIX_FLOAT * vals;
  PASTIX_INT   * index2;
  PASTIX_FLOAT * vals2;
  Clock  clk;
  void  * ptr[2];

  if (NULL == (index  = memAlloc(sizeof(PASTIX_INT)*NB_ELEMENTS)))
    MALLOC_ERROR("index");
  if (NULL == (index2  = memAlloc(sizeof(PASTIX_INT)*NB_ELEMENTS)))
    MALLOC_ERROR("index2");
  if (NULL == (vals  = memAlloc(sizeof(PASTIX_FLOAT)*NB_ELEMENTS)))
    MALLOC_ERROR("vals");
  if (NULL == (vals2  = memAlloc(sizeof(PASTIX_FLOAT)*NB_ELEMENTS)))
    MALLOC_ERROR("vals2");

  memcpy(index2, index, sizeof(PASTIX_INT)*NB_ELEMENTS);
  memcpy(vals2, vals, sizeof(PASTIX_FLOAT)*NB_ELEMENTS);

  clockInit(&clk);
  clockStart(&clk);
  ptr[0] = index;
  ptr[1] = vals; 
  qsortIntFloatAsc(ptr,NB_ELEMENTS);

  clockStop (&clk);
  for (i = 0; i < NB_ELEMENTS -1; i++)
    if (index[i] > index[i+1])
      {
	errorPrint("Mauvais tri 1 : mal ordonné");
	return EXIT_FAILURE;
      }
  
  for (i = 0; i < NB_ELEMENTS; i++)
    {
      j = 0;
      while (index[j] != index2[i] && j < NB_ELEMENTS)
	j++;
      while (index[j] == index2[i] && vals[j] != vals2[i] && j < NB_ELEMENTS)
	j++;

      if (j == NB_ELEMENTS)
	{
	  errorPrint("Mauvais tri 1: élément disparu");
	  return EXIT_FAILURE;
	}
      if (index[j] != index2[i] || vals[j] != vals2[i])
	{
	  errorPrint("Mauvais tri 1 : mauvaise recopie de valeur associee a l'index");
	  return EXIT_FAILURE;
	}
    }
  fprintf(stdout,"tri de %ld elements aleatoires correct en %.3g s\n", (long)NB_ELEMENTS, (float)clockVal(&clk));


  /* Sort sorted elements */
  for (i = 0; i < NB_ELEMENTS; i++)
    {
      index[i] = i;
      vals[i] = (PASTIX_FLOAT)random();
    }

  memcpy(index2, index, sizeof(PASTIX_INT)*NB_ELEMENTS);
  memcpy(vals2, vals, sizeof(PASTIX_FLOAT)*NB_ELEMENTS);
  clockInit(&clk);
  clockStart(&clk);
  ptr[0] = index;
  ptr[1] = vals; 
  qsortIntFloatAsc(ptr,NB_ELEMENTS);
  clockStop (&clk);
  for (i = 0; i < NB_ELEMENTS -1; i++)
    if (index[i] > index[i+1])
      {
	errorPrint("Mauvais tri 2 : mal ordonné");
	return EXIT_FAILURE;
      }
  
  for (i = 0; i < NB_ELEMENTS; i++)
    {
      j = 0;
      while (index[j] != index2[i] && j < NB_ELEMENTS)
	j++;
      while (index[j] == index2[i] && vals[j] != vals2[i] && j < NB_ELEMENTS)
	j++;

      if (j == NB_ELEMENTS)
	{
	  errorPrint("Mauvais tri 2 : élément disparu");
	  return EXIT_FAILURE;
	}
      if (index[j] != index2[i] || vals[j] != vals2[i])
	{
	  errorPrint("Mauvais tri 2 : mauvaise recopie de valeur associee a l'index");
	  return EXIT_FAILURE;
	}
    }
  fprintf(stdout,"tri de %ld elements triés correct en %.3g s \n", (long)NB_ELEMENTS, (float)clockVal(&clk));

  memFree_null(index);
  memFree_null(index2);
  memFree_null(vals);
  memFree_null(vals2);

  return EXIT_SUCCESS;
}

