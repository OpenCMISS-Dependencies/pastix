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
#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#include "common_pastix.h"
#include "tools.h"
#include "cscd_utils.h"

int main(int argc, int argv)
{

  int i;

  PASTIX_INT    n1;
  PASTIX_INT   *colptr1 = NULL;
  PASTIX_INT   *rows1   = NULL;
  PASTIX_INT   *l2g1    = NULL;
  PASTIX_FLOAT *values1 = NULL;
  
  PASTIX_INT    n2;
  PASTIX_INT   *colptr2 = NULL;
  PASTIX_INT   *rows2   = NULL;
  PASTIX_INT   *l2g2    = NULL;
  PASTIX_FLOAT *values2 = NULL;

  PASTIX_INT    n3;
  PASTIX_INT   *colptr3 = NULL;
  PASTIX_INT   *rows3   = NULL;

  PASTIX_FLOAT *values3 = NULL;

  /*
    First matrix : 
      0 . 1 0 .
      1 . 0 0 .
      0 . 2 0 .
      3 . 7 1 .
      4 . 0 0 .
   */
  n1 = 3;
  l2g1    = malloc(n1*sizeof(PASTIX_INT));
  l2g1[0] = 1;
  l2g1[1] = 3;
  l2g1[2] = 4;

  colptr1    = malloc((n1+1)*sizeof(PASTIX_INT));
  colptr1[0] = 1;
  colptr1[1] = 4;
  colptr1[2] = 7;
  colptr1[3] = 8;

  rows1    = malloc((colptr1[n1] -1)*sizeof(PASTIX_INT));
  rows1[0] = 2;
  rows1[1] = 4;
  rows1[2] = 5;
  rows1[3] = 1;
  rows1[4] = 3;
  rows1[5] = 4;
  rows1[6] = 4;

  values1    = malloc((colptr1[n1] -1)*sizeof(PASTIX_FLOAT));
  values1[0] = 1;
  values1[1] = 3;
  values1[2] = 4;
  values1[3] = 1;
  values1[4] = 2;
  values1[5] = 7;
  values1[6] = 1;

  /* 
     Second matrix empty with local column 1 and 2
   */

  n2      = 2;
  l2g2    = malloc(n2*sizeof(PASTIX_INT));
  l2g2[0] = 1;
  l2g2[1] = 2;
  
  colptr2    = malloc((n2+1)*sizeof(PASTIX_INT));
  colptr2[0] = 1;
  colptr2[1] = 1;
  colptr2[2] = 1;
  
  rows2      = NULL;
  values2    = NULL;
  
  cscd_addlocal(n1,   colptr1,  rows1, NULL, l2g1,
		n2,   colptr2,  rows2, NULL, l2g2,
		&n3, &colptr3, &rows3, &values3, CSCD_ADD);

  /* CSCD 3 must be equal to CSCD1 */
  if (n1 != n3)
    {
      fprintf(stderr, "Error adding empty CSCD : n1 != n3\n");
      return EXIT_FAILURE;
    }
  for(i = 0; i < n2+1; i++)
    if (colptr1[i] != colptr3[i])
      {
	fprintf(stderr, "Error adding empty CSCD : colptr1[i] != colptr3[i]\n");
	return EXIT_FAILURE;
      }
  
  for(i = 0; i < colptr3[n2]-1; i++)
    if (rows1[i] != rows3[i])
      {
	fprintf(stderr, "Error adding empty CSCD : rows1[i] != rows3[i]\n");
	return EXIT_FAILURE;
      }
  if (values3 != NULL)
     {
	fprintf(stderr, "Error adding empty CSCD with no values : value3 != NULL\n");
	return EXIT_FAILURE;
      }

  free(colptr3);
  colptr3 = NULL;
  free(rows3);
  rows3   = NULL;

  cscd_addlocal(n1,   colptr1,  rows1, values1, l2g1,
		n2,   colptr2,  rows2, values2, l2g2,
		&n3, &colptr3, &rows3, &values3, CSCD_ADD);

  /* CSCD 3 must be equal to CSCD1 */
  if (n1 != n3)
    {
      fprintf(stderr, "Error adding empty CSCD : n1 != n3\n");
      return EXIT_FAILURE;
    }
  for(i = 0; i < n2+1; i++)
    if (colptr1[i] != colptr3[i])
      {
	fprintf(stderr, "Error adding empty CSCD : colptr1[i] != colptr3[i]\n");
	return EXIT_FAILURE;
      }
  
  for(i = 0; i < colptr3[n2]-1; i++)
    if (rows1[i] != rows3[i])
      {
	fprintf(stderr, "Error adding empty CSCD : rows1[i] != rows3[i]\n");
	return EXIT_FAILURE;
      }
 
  for(i = 0; i < colptr3[n2]-1; i++)
    if (values1[i] != values3[i])
      {
	fprintf(stderr, "Error adding empty CSCD : values1[i] != values3[i]\n");
	return EXIT_FAILURE;
      }
  free(colptr3);
  colptr3 = NULL;
  free(rows3);
  rows3   = NULL;
  free(values3);
  values3 = NULL;

  colptr2[1] = 4;
  colptr2[2] = 5;

  rows2    = malloc((colptr2[n2]-1)*sizeof(PASTIX_INT));
  rows2[0] = 1;
  rows2[1] = 2;
  rows2[2] = 3;
  rows2[3] = 4;

  values2    = malloc((colptr2[n2]-1)*sizeof(PASTIX_FLOAT));
  values2[0] = 1;
  values2[1] = 2;
  values2[2] = 3;
  values2[3] = 4;

  cscd_addlocal(n1,   colptr1,  rows1, values1, l2g1,
		n2,   colptr2,  rows2, values2, l2g2,
		&n3, &colptr3, &rows3, &values3, CSCD_ADD);

  if (values3[1] != 3)
    {
      fprintf(stderr, "Error adding  CSCD : values3[1] != 3\n");
      return EXIT_FAILURE;
    }

  free(colptr3);
  free(values3);
  free(rows3);
  
  cscd_addlocal(n1,   colptr1,  rows1, values1, l2g1,
		n2,   colptr2,  rows2, values2, l2g2,
		&n3, &colptr3, &rows3, &values3, CSCD_KEEP);

  if (values3[1] != 1)
    {
      fprintf(stderr, "Error adding  CSCD : values3[1] != 1\n");
      return EXIT_FAILURE;
    }
  free(colptr3);
  free(values3);
  free(rows3);
  
  cscd_addlocal(n1,   colptr1,  rows1, values1, l2g1,
		n2,   colptr2,  rows2, values2, l2g2,
		&n3, &colptr3, &rows3, &values3, CSCD_OVW);

  if (values3[1] != 2)
    {
      fprintf(stderr, "Error adding  CSCD : values3[1] != 2\n");
      return EXIT_FAILURE;
    }


  free(colptr1);
  free(colptr2);
  free(colptr3);

  free(rows1);
  free(rows2);
  free(rows3);

  free(values1);
  free(values2);
  free(values3);

  free(l2g1);
  free(l2g2);


  return EXIT_SUCCESS;
}
