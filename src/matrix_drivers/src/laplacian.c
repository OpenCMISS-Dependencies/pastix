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
 * File: laplacian.c
 *
 * Generate a laplacian
 *
 * Example :
 * >  2 -1  0  0
 * > -1  2 -1  0
 * >  0 -1  2 -1
 * >  0  0 -1  2
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
#include "laplacian.h"
/*
 * Function: genlaplacien
 *
 * Generate a laplacien of size *n*
 *
 * Parameters:
 *   n       - Size of the wanted matrix
 *   nnzeros - Number of non zeros in the produced matrice
 *   ia      - Index of first element of each column in *row* and *val*
 *   ja      - Row of eah element
 *   avals   - Value of each element
 *   rhs     - Right-hand-side member
 *   type    - Type of the matrix
 *   rhstype - Type of the right hand side.
 */
int genlaplacian(pastix_int_t     n,
                 pastix_int_t    *nnzeros,
                 pastix_int_t   **ia,
                 pastix_int_t   **ja,
                 pastix_float_t **avals,
                 pastix_float_t **rhs,
                 char           **type,
                 char           **rhstype)
{

  pastix_int_t i;
  pastix_int_t j;

  *nnzeros = 3*n - 2;
  *ia      = NULL;
  *ja      = NULL;
  *avals   = NULL;
  *rhs     = NULL;
  *type    = NULL;
  *rhstype = NULL;

  /* Allocating */
  if ((NULL == (*ia       = (pastix_int_t *)  malloc((n+1)     *sizeof(pastix_int_t))  )) ||
      (NULL == (*ja       = (pastix_int_t *)  malloc((*nnzeros)*sizeof(pastix_int_t))  )) ||
      (NULL == (*avals    = (pastix_float_t *)malloc((*nnzeros)*sizeof(pastix_float_t)))) ||
      (NULL == (*rhs      = (pastix_float_t *)malloc(n         *sizeof(pastix_float_t)))) ||
      (NULL == (*type     = (char *)          malloc(4         *sizeof(char))      )) ||
      (NULL == (*rhstype  = (char *)          malloc(4         *sizeof(char))          )))
  {
    fprintf(stderr, "Error in CSC allocation\n");
    if (*type != NULL)
    {
      free(*type);
      *type = NULL;
    }
    if (*rhs != NULL)
    {
      free(*rhs);
      *rhs = NULL;
    }
    if (*avals != NULL)
    {
      free(*avals);
      *avals = NULL;
    }
    if (*ja != NULL)
    {
      free(*ja);
      *ja = NULL;
    }
    if (*ia != NULL)
    {
      free(*ia);
      *ia = NULL;
    }
    return EXIT_FAILURE;
  }

  /* Building ia, ja and avals and rhs*/
  j=0;
  for (i = 0; i < n; i++)
    {
      (*ia)[i] = j+1;
      /* ONLY triangular inferior matrix */
      /*       if (i != 0) */
      /*  { */
      /*    (*ja)[j]    = i; */
      /*    (*avals)[j] = -1; */
      /*    j++; */
      /*  } */
      (*ja)[j]    = i+1;
      (*avals)[j] = 2;
      j++;
      if (i != n-1)
      {
        (*ja)[j]    = i+2;
#ifdef TYPE_COMPLEX
        (*avals)[j] = - 1 +  2* _Complex_I;
#else
        (*avals)[j] = -1;
#endif
        j++;
      }
      (*rhs)[i] = 0;
#if (defined TYPE_COMPLEX && defined SYMMETRIC_LAPLACIAN)
      (*rhs)[i] = -4 + 4 * _Complex_I;
#endif
    }
  (*ia)[i] = j+1;
#ifdef TYPE_COMPLEX
  (*rhs)[0] = 3 - _Complex_I;
#ifdef SYMMETRIC_LAPLACIAN
  (*rhs)[0] =  -1 + 3 *  _Complex_I;
#endif
  (*rhs)[n-1] =  -1 + 3 *  _Complex_I;
#else
  (*rhs)[0] = 1;
  (*rhs)[n-1] = 1;
#endif
  /* type and rhstype */
#if (defined TYPE_COMPLEX && !defined SYMMETRIC_LAPLACIAN)
  sprintf (*type, "RHA");
#else
  sprintf (*type, "RSA");
#endif
  sprintf (*rhstype,"???");
  return EXIT_SUCCESS;
}
