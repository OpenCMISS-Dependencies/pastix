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

/************************************************************/
/**                                                        **/
/**   NAME       : amalgamate.c                            **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 15/08/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "common_pastix.h"
/* #include "symbol.h" */

#include "sparRow.h"
#include "SF_level.h"



long SF_level(PASTIX_INT job, csptr A, PASTIX_INT level, csptr P)
{
  /***********************************************************************************************/
  /* This function computes the non zero pattern of the levelized incomplete factor              */
  /* for a sparse lower triangular                                                               */
  /* matrix in CSC format.  This pattern is exact iff the matrix has a SYMMETRIC non zero        */
  /* structure.                                                                                  */
  /* On entry:                                                                                   */
  /* job : = 0 alloc the matrix ; =1 fill the indices of the sparse pattern  =2 alloc (not cont) */
  /*       and fill the coefficients                                                             */
  /* A : pointer to the matrix                                                                   */
  /* level : level desired for the ilu(k) factorization                                          */
  /* P   : an empty csr matrix (initialized with dimension n)                                    */
  /* On return:                                                                                  */
  /*     P : a csr matrix containing the non zero pattern of the factorized matrix               */
  /*         the memory for the numerical values is allocated and initialized to zero            */
  /*     The total number of nnz in the lower part is returned                                   */
  /*                                                                                             */
  /* NOTE:                                                                                       */
  /*   1) This algorithm has been implemented according to the paper of David Hysom and          */
  /*     Alex Pothen : Level-based Incomplete LU factorization: Graph Model and Algorithm        */
  /***********************************************************************************************/
  PASTIX_INT *visited = NULL;
  PASTIX_INT *length  = NULL;
  PASTIX_INT *stack   = NULL;
  PASTIX_INT *adj     = NULL;
  PASTIX_INT *ja      = NULL;
  PASTIX_INT used;
  PASTIX_INT h, i,j,k, t;
  long nnz;

  if(A->n == 0)
    return 0;
  /** Allocated the working array **/
  MALLOC_INTERN(visited, A->n, PASTIX_INT);
  MALLOC_INTERN(length,  A->n, PASTIX_INT);
  MALLOC_INTERN(stack,   A->n, PASTIX_INT);
  MALLOC_INTERN(ja,      A->n, PASTIX_INT);
  nnz = 0;

  /** Initialized visited ***/
  for(j=0;j<A->n;j++)
    {
      visited[j] = -1;
      length[j] = 0;
    }


  /** Apply GS_Urow for each row **/
  for(i=0;i<A->n;i++)
    {
      used = 0; /** Reset the stack number of elements **/
      stack[0] = i;
      used++;
      length[i] = 0;
      visited[i] = i;

      ja[0] = i; /** Put the diagonal term **/
      k=1;
      /** BFS phase **/
      while(used > 0)
        {

          used--;
          h = stack[used];
          adj = A->ja[h];
          for(j=0;j<A->nnzrow[h];j++)
            {
              t = adj[j];
              if(visited[t] != i)
                {
                  visited[t] = i;
                  if(t<i && length[h] < level)
                    {

                      stack[used] = t;
                      used++;
                      length[t] = length[h]+1;
                    }
                  if(t>i)
                    {
                      ja[k++] = t;
                    }
                }

            }
        }

      /*** allocate the new row and fill it ***/
      switch(job)
        {
        case 0:
          P->nnzrow[i] = k;
          break;
        case 1:
          P->ja[i] = P->ja[0] + nnz;
          memcpy(P->ja[i], ja, sizeof(PASTIX_INT)*k);
          break;
        case 2:
          P->nnzrow[i] = k;
#ifdef DEBUG_KASS
          assert(k>0);
#endif
          MALLOC_INTERN(P->ja[i], k, PASTIX_INT);
          memcpy(P->ja[i], ja, sizeof(PASTIX_INT)*k);
          break;
        }
      nnz += k;
    }
  memFree_null(ja);
  memFree_null(visited);
  memFree_null(length);
  memFree_null(stack);

  if(job == 0)
    {
      /*fprintf(stderr, "NNZ in the triangular factor L = %ld \n", nnz);*/
      P->inarow = 1;
      MALLOC_INTERN(P->jatab, nnz, PASTIX_INT);
      P->ja[0] = P->jatab;

      P->matab = NULL;
      P->ma[0] = NULL;

      ASSERT(P->ja[0] != NULL,MOD_BLEND);
      /*ASSERT(P->ma[0] != NULL,MOD_BLEND);*/

    }

  return nnz;
}
