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
/**   NAME       : find_supernodes.c                       **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10/02/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "common_pastix.h"
#include "find_supernodes.h"

/** local function **/
void  post_order(PASTIX_INT n, PASTIX_INT *father, PASTIX_INT *T,  PASTIX_INT *perm, PASTIX_INT *iperm);
void compute_subtree_size(PASTIX_INT n, PASTIX_INT *father, PASTIX_INT *perm, PASTIX_INT *iperm, PASTIX_INT *T);
void compute_elimination_tree(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_INT *perm, PASTIX_INT *iperm, PASTIX_INT *father);

void  find_supernodes(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_INT *perm, PASTIX_INT *iperm, PASTIX_INT *snodenbr, PASTIX_INT *snodetab, PASTIX_INT *treetab)
{
  /********************************************************************************************/
  /* This function computes the supernodal partition of a reordered matrix.                   */
  /* The permutation of the matrix given on entry is modified to obtain a                     */
  /* postorder of the elimination tree: this does not affect the fill-in                      */
  /* properties of the initial ordering.                                                      */
  /* The matrix pattern is assumed to be symmetric                                            */
  /* NOTE :                                                                                   */
  /* This function can take on entry the lower triangular part or the whole matrix A:         */
  /* this does not change the results                                                         */
  /*                                                                                          */
  /* On entry:                                                                                */
  /* n, ia, ja : the initial matrix A in CSC format                                           */
  /*             in C numbering                                                               */
  /*  perm : the permutation vector                                                           */
  /* iperm : the inverse permutation vector                                                   */
  /* snodetab : allocated to contain at most n integers                                       */
  /* treetab  : allocated to contain at most n integers                                       */
  /*            it can be set to NULL in which case the return value is treetab = NULL        */
  /* On return :                                                                              */
  /* perm, iperm : postorder of the node in the elimination tree deduced from the initial     */
  /*               ordering                                                                   */
  /* snodenbr  : number of supernodes found                                                   */
  /* snodetab  : snodetab[i] is the beginning in the new ordering of the ith supernodes       */
  /* treetab   : treetab[s] is the number of the father of supernode i on the supernodal      */
  /*             elimination tree                                                             */
  /********************************************************************************************/

  PASTIX_INT *father     = NULL; /** father[i] is the father of node i in he elimination tree of A **/
  PASTIX_INT *T          = NULL; /** T[j] is the number of node in the subtree rooted in node j in
                              the elimination tree of A **/
  PASTIX_INT *S          = NULL; /** S[i] is the number of sons for node i in the elimination tree **/
  PASTIX_INT *isleaf     = NULL;
  PASTIX_INT *prev_rownz = NULL;
  PASTIX_INT i, j, k;
  PASTIX_INT pi, pj;
  PASTIX_INT dad;


  MALLOC_INTERN(T,          n, PASTIX_INT);
  MALLOC_INTERN(S,          n, PASTIX_INT);
  MALLOC_INTERN(father,     n, PASTIX_INT);
  MALLOC_INTERN(isleaf,     n, PASTIX_INT);
  MALLOC_INTERN(prev_rownz, n, PASTIX_INT);
#ifdef DEBUG_BLEND
  assert(ia[0] == 0);
#endif


#ifdef DEBUG_BLEND
  /** Check the permutation vector **/
  for(i=0;i<n;i++)
    {
      assert(perm[i] >= 0);
      assert(perm[i] < n);
    }

  bzero(S, sizeof(PASTIX_INT)*n);
  for(i=0;i<n;i++)
    S[perm[i]]++;

  k = 0;
  for(i=0;i<n;i++)
    if(S[i] != 1)
      k++;
  if(k>0)
    errorPrint("perm array is not valid, number of error =  %ld", (long)k);
  assert(k==0);
#endif


  /*** Compute the elimination tree of A ***/
  compute_elimination_tree(n, ia, ja, perm, iperm, father);


  /*** Compute the postorder of the elimination tree ***/
  /*** This operation modifies perm and iperm ***/
  post_order(n, father, T, perm, iperm);

  /*** Compute the number of descendant of each node i in the elimination tree ***/
  compute_subtree_size(n, father, perm, iperm, T);


  bzero(isleaf, sizeof(PASTIX_INT)*n);
  bzero(prev_rownz, sizeof(PASTIX_INT)*n);


  for(j=0;j<n;j++)
    {
      pj = iperm[j];
      for(i=ia[pj];i<ia[pj+1];i++)
        {
          pi = perm[ja[i]];
          if(pi > j)
            {
              k = prev_rownz[pi];
              if(k < j - T[pj]+1 )
                isleaf[j] = 1;

              prev_rownz[pi] = j;
            }
        }
    }


  /*** Compute the number of sons of each node in the elimination tree ***/
  bzero(S, sizeof(PASTIX_INT)*n);
  for(i=0;i<n;i++)
    if(father[i] != i)
      S[father[i]]++;

  for(i=0;i<n;i++)
    if(S[i] != 1)
      isleaf[perm[i]] = 1;

  (*snodenbr) = 0;
  for(i=0;i<n;i++)
    if(isleaf[i] == 1)
      {
        snodetab[(*snodenbr)] = i;
        (*snodenbr)++;
      }
  snodetab[(*snodenbr)] = n;


  if(treetab != NULL)
    {
      /*** Node to supernodet conversion vector ***/
      for(i=0;i<(*snodenbr);i++)
        for(j=snodetab[i];j<snodetab[i+1];j++)
          S[j] = i;

      /*** Fill the treetab info ***/
      for(i=0;i<(*snodenbr);i++)
        {
          k=(*snodenbr);
          for(j=snodetab[i];j<snodetab[i+1];j++)
            {
              dad = S[perm[father[iperm[j]]]];
              if( dad < k && dad > i)
                k = dad;
            }
          treetab[i] = k;
          if(k==(*snodenbr))
            {
              /*fprintf(stdout, "THIS IS A ROOT %d \n", i);*/
              treetab[i] = i; /** This is a root **/
            }
#ifdef DEBUG_BLEND
          assert(treetab[i] >= i);
#endif
        }


    }

  memFree(prev_rownz);
  memFree(isleaf);
  memFree(father);
  memFree(S);
  memFree(T);






}


void  post_order(PASTIX_INT n, PASTIX_INT *father, PASTIX_INT *T,  PASTIX_INT *perm, PASTIX_INT *iperm)
{
  /********************************************************************************/
  /* This function compute the post order of the elimination tree given on entry  */
  /* On entry:                                                                    */
  /* n : number of nodes                                                          */
  /* father : father[i] is the node number of the father of node i                */
  /*          if node i is a root then father[i] = i                              */
  /* T : a temporary vector of size n                                             */
  /* perm, iperm : ordering of the matrix (value optional: if set then the post   */
  /*                            ordering try to keep the initial ordering as much */
  /*                            as possible)                                      */
  /* On return :                                                                  */
  /* perm, iperm : permutation and inverse permutation vector of the postorder    */
  /********************************************************************************/
  PASTIX_INT i;
  PASTIX_INT j, k, t;


  /*** First compute the number of node in the subtree rooted in node i ***/
  compute_subtree_size(n, father, perm, iperm, T);

  /** When multiple roots are present we have to compute the start index of each root **/
  t=0;
  for(k=0;k<n;k++)
    {
      i = iperm[k];
      if(father[i] == i)
        {
          /** This is a root **/
          j = T[i];
          T[i] += t;
          t += j;
        }
    }

#ifdef DEBUG_BLEND
  for(i=0;i<n;i++)
    assert(T[i] <= T[father[i]]);
#endif



 for(k=n-1;k>=0;k--)
   {
     i = iperm[k];
     perm[i] = T[father[i]]; /** We MUST HAVE father[i] == i for a root ! **/
     T[father[i]]-= T[i];
     T[i] = perm[i]-1;
#ifdef DEBUG_BLEND
     assert(perm[father[i]] >= perm[i]);
#endif
   }



  /** We need to retrieve 1 for the C numbering compliance **/
  for(i=0;i<n;i++)
    perm[i]--;

#ifdef DEBUG_BLEND
  /** Check the permutation vector **/
  for(i=0;i<n;i++)
    {
      /*fprintf(stderr, "(%d = %d) ", i, perm[i]);*/
      assert(perm[i] >= 0);
      assert(perm[i] < n);
    }

  bzero(iperm, sizeof(PASTIX_INT)*n);
  for(i=0;i<n;i++)
    iperm[perm[i]]++;

  k = 0;
  for(i=0;i<n;i++)
    if(iperm[i] != 1)
      k++;
  if(k>0)
    errorPrint("Number of errors in perm vector in postorder %ld", (long)k);
  assert(k==0);
#endif


  /** Compute the iperm vector **/
  for(i=0;i<n;i++)
    iperm[perm[i]] = i;



}



void compute_subtree_size(PASTIX_INT n, PASTIX_INT *father, PASTIX_INT *perm, PASTIX_INT *iperm, PASTIX_INT *T)
{
  /********************************************/
  /*  Compute the size of each subtree given  */
  /*  the number of the father of each node   */
  /********************************************/
  PASTIX_INT k, i;
  (void)perm;

  /*** OIMBE pas la peine d'utiliser un tas; il suffit de parcourir iperm pour assurer
       de toujours traiter un fils avant son pere ***/

  bzero(T, sizeof(PASTIX_INT)*n);

  for(k=0;k<n;k++)
    {
      i = iperm[k];
      T[i]++;
      if(i!=father[i])
        T[father[i]] += T[i];
    }



#ifdef DEBUG_BLEND
 {
   PASTIX_INT sum;
   sum = 0;
   for(i=0;i<n;i++)
     if(father[i] == i)
       sum += T[i];

   if(sum != n)
     errorPrint("compute_subtree_size: sum of the subtree = %ld n = %ld", (long)sum, (long)n);
   assert(n == sum);
 }
#endif

}




void compute_elimination_tree(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_INT *perm, PASTIX_INT *iperm, PASTIX_INT *father)
{
  /******************************************************************************/
  /* Compute the elimination tree of a matrix A (without computing the symbolic */
  /* factorization) associated with a reordering of the matrix                  */
  /* On entry:                                                                  */
  /* n, ia, ja : the adjacency graph of the matrix (symmetrized)                */
  /* perm : a permutation vector of the matrix                                  */
  /* iperm : the inverse permutation vector                                     */
  /* On return:                                                                 */
  /* father : father[i] = father of node i on the eliminination tree            */
  /* If node i is a root then father[i] = i                                     */
  /* NOTE : father is allocated at a size of n interger on input                */
  /******************************************************************************/
  PASTIX_INT i, j, k;
  PASTIX_INT node;
  PASTIX_INT vroot;

  /** Optim **/
  PASTIX_INT flag, ind;
  PASTIX_INT *jrev = NULL;
  PASTIX_INT *jf   = NULL;


  MALLOC_INTERN(jrev, n, PASTIX_INT);
  for(i=0;i<n;i++)
    jrev[i] = -1;
  MALLOC_INTERN(jf, n, PASTIX_INT);
  bzero(jf, sizeof(PASTIX_INT)*n);

  for(i=0;i<n;i++)
    father[i] = -1;

  for(i=0;i<n;i++)
    {
      ind = 0;
      node = iperm[i];
      for(j=ia[node];j<ia[node+1];j++)
        {

          k = ja[j];
          if(perm[k] < perm[node])
            {
              flag = 1;
              vroot = k;
              while(father[vroot] != -1 && father[vroot] != node)
                {
                  if(jrev[vroot] >= 0)
                    {
                      flag = 0;
                      break;
                    }
                  jrev[vroot] = ind;
                  jf[ind] = vroot;
                  ind++;

                  vroot = father[vroot];
                }
              if(flag == 1)
                father[vroot] = node;
            }
        }
      /** reinit jrev **/
      for(j=0;j<ind;j++)
        jrev[jf[j]]=-1;
    }

  memFree_null(jrev);
  memFree_null(jf);

  for(i=0;i<n;i++)
    if(father[i] == -1)
      father[i]=i;

#ifdef DEBUG_KASS
  /*** Check to see if a father has a lower rank in the permutation array than one of its sons ***/
  for(i=0;i<n;i++)
    {
      if(perm[i] > perm[father[i]])
        {
          fprintf(stderr, "Node %ld perm=%ld Father %ld perm=%ld \n", (long)i, (long)perm[i], (long)father[i], (long)perm[father[i]]);
          assert(perm[i] <= perm[father[i]]);
        }
    }
#endif
}
