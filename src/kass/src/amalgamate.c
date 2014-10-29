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

#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#include "common_pastix.h"
/* #include "symbol.h" */
#include "queue.h"
#include "perf.h"
#include "sparRow.h"
/* #include "sort_row.h" */
#include "amalgamate.h"

/*#define BLAS_GAIN*/    /** Amalgamation use the best ratio time/nnzadd
                         to merge cblk **/
/*#define BLAS_SOLVE*/   /** Amalgamation seeks to optimize the triangular
                         solve **/

#ifdef BLAS_SOLVE
#define CBLKTIME  cblk_time_fact
#else
#define CBLKTIME  cblk_time_solve
#endif

/** 2 percents **/
#define RAT_CBLK 0.02

#define INFINI 10e6

#define print_one(fmt, ...)    if( procnum == 0) fprintf(stdout, fmt, ##__VA_ARGS__)

extern void UnionSet(PASTIX_INT *set1, PASTIX_INT n1, PASTIX_INT *set2, PASTIX_INT n2, PASTIX_INT *set, PASTIX_INT *n);

double cblk_time_fact (PASTIX_INT n, PASTIX_INT *ja, PASTIX_INT colnbr);
double cblk_time_solve(PASTIX_INT n, PASTIX_INT *ja, PASTIX_INT colnbr);
PASTIX_INT    merge_cost(PASTIX_INT a, PASTIX_INT b, csptr P, PASTIX_INT *colweight);
double merge_gain(PASTIX_INT a, PASTIX_INT b, csptr P, PASTIX_INT *colweight, PASTIX_INT *tmp);
void   merge_col (PASTIX_INT a, PASTIX_INT b, csptr P);
void   get_son(PASTIX_INT node, PASTIX_INT *sonindex, PASTIX_INT *sontab, PASTIX_INT *colweight, PASTIX_INT *ns, PASTIX_INT *list);



void amalgamate(double rat, csptr P, PASTIX_INT snodenbr, PASTIX_INT *snodetab, PASTIX_INT *treetab, PASTIX_INT *cblknbr, PASTIX_INT **rangtab, PASTIX_INT *nodetab, MPI_Comm pastix_comm)
{
  /********************************************************************/
  /* amagamate takes a supernode graph (P,colwgt) and performs some   */
  /* amalgation of the column until a fill tolerance (rat) is reached */
  /* It returns the partition of columns obtained                     */
  /* On return                                                        */
  /* P  is updated to form the column compressed graph (not the
     quotient graph !!) */
  /* treetab is also updated to fit with the new column compressed    */
  /*   graph
       */
  /* cblknbr, rangtab is the supernode partition of the initial
     matrix P */
  /********************************************************************/

  PASTIX_INT *nnzadd    = NULL;
  PASTIX_INT *sonindex  = NULL;

  PASTIX_INT *sontab    = NULL;
  PASTIX_INT *treetab2  = NULL; /** The treetab is updated as amalgamation comes
                             along **/
  PASTIX_INT *tmp       = NULL;
  PASTIX_INT *tmp2      = NULL;
  PASTIX_INT *newnum    = NULL;
  PASTIX_INT i,j,k,ind;
  PASTIX_INT father;
  PASTIX_INT cblknum;
  PASTIX_INT *colweight = NULL;
  double *gain   = NULL;
  PASTIX_INT toto;
  PASTIX_INT n, nn;
  double rat_cblk;
  int blas_gain = 0;

  double key;
  Queue heap;
  long fillmax, fill;
  long fillcblk, fillwhile;

  int procnum;
  (void)pastix_comm;

  MPI_Comm_rank(pastix_comm, &procnum);

  if(rat < 0)
    {
      rat = -rat;
      rat_cblk = MIN(rat, RAT_CBLK);
    }
  else
    rat_cblk = rat;


  blas_gain = 0; /** We always begin by amalgation leaded by the
                     reduction of the number of cblk (until rat_cblk)  **/

  n = P->n; /** Number of supernode **/
  if(snodetab != NULL)
    nn = snodetab[snodenbr]; /** Number of unknowns **/
  else
    nn = n;


#ifdef DEBUG_KASS

  if(snodetab != NULL)
    {
      ASSERT(n == snodenbr, MOD_KASS);
      for(i=0;i<n;i++)
        {
          ASSERT(P->nnzrow[i] >= snodetab[i+1]-snodetab[i], MOD_KASS);
          k = 0;
          for(j=snodetab[i]; j < snodetab[i+1];j++)
            ASSERT(P->ja[i][k++] == j, MOD_KASS);
        }
    }
#endif

  /*** Set the weight of each column ***/
  MALLOC_INTERN(colweight, n, PASTIX_INT);
  if(snodetab == NULL)
    {
      for(i=0;i<n;i++)
        colweight[i] = 1;
    }
  else
    {
      for(i=0;i<snodenbr;i++)
        colweight[i] = snodetab[i+1]-snodetab[i];
    }



#ifdef DEBUG_KASS
  for(i=0;i<n;i++)
    ASSERT(colweight[i] >= 0, MOD_KASS);
#endif


  /**********************************************/
  /*** Compute the maximum extra fill allowed ***/
  /**********************************************/
  /** fillmax is the limit of nnz added du to the whole amalgamation
      process **/
  /** fillcblk (<fillmax) is the limit under which the amalgamation
      tries to reduce the number of supernode **/
  /** between fillcblk and fillmax the amalgamation tries to reduce
      the time of the solve or factorization time **/

  if(snodetab == NULL)
    {
      fillmax = (long)(CSnnz(P)*rat);
      fillcblk = (long)(CSnnz(P)*rat_cblk);
    }
  else
    {
      fillmax = 0;
      for(i=0;i<P->n;i++)
        {
          fillmax += (colweight[i]*(colweight[i]+1))/2;
#ifdef DEBUG_KASS
          if(P->nnzrow[i] < colweight[i])
            fprintf(stderr, "nnzrow[%ld] = %ld colweight = %ld \n", (long)i, (long)P->nnzrow[i], (long)colweight[i]);
          ASSERT(P->nnzrow[i] >= colweight[i], MOD_KASS);
#endif
          fillmax += (P->nnzrow[i]-colweight[i])*colweight[i];
        }
#ifdef DEBUG_KASS
      fprintf(stderr, "(NNZL) = %ld \n", (long)fillmax);
#endif
      fillcblk = (long)(fillmax*rat_cblk);
      fillmax = (long)(fillmax*rat);
    }
#ifdef DEBUG_KASS
  print_one("fillcblk %ld fillmax %ld \n", (long)fillcblk, (long)fillmax);
#endif
  fillwhile = fillcblk;






  MALLOC_INTERN(treetab2, n, PASTIX_INT);
  MALLOC_INTERN(nnzadd,   n, PASTIX_INT);
  MALLOC_INTERN(gain,     n, double);
  memCpy(treetab2, treetab, sizeof(PASTIX_INT)*n);
  MALLOC_INTERN(tmp,  nn, PASTIX_INT);
  MALLOC_INTERN(tmp2, nn, PASTIX_INT);

#ifdef DEBUG_KASS
  key = 0.0;
  for(i=0;i<P->n;i++)
    key += CBLKTIME(P->nnzrow[i], P->ja[i], snodetab[i+1]-snodetab[i]);
  fprintf(stderr, "COST of the NON AMALGAMTED MATRIX = %g \n", key);
#endif



  /**********************************************/
  /*** Compute the son list of each supernode ***/
  /**********************************************/
  /*** Compute the number of sons of each cblk ***/
  bzero(tmp, sizeof(PASTIX_INT)*n);
  for(i=0;i<n-1;i++)
    if(treetab2[i] >= 0) /** IF THIS SNODE IS NOT A ROOT **/
      tmp[treetab2[i]]++;

  ind = 0;
  for(i=0;i<n;i++)
    ind += tmp[i];

  MALLOC_INTERN(sonindex, n+1, PASTIX_INT);
  MALLOC_INTERN(sontab,   ind, PASTIX_INT);
  ind = 0;
  for(i=0;i<n;i++)
    {
      sonindex[i] = ind;
      ind += tmp[i];
    }
  sonindex[n] = ind;

  bzero(tmp, sizeof(PASTIX_INT)*n);

  for(i=0;i<n-1;i++)
    {
      cblknum = treetab2[i];
      if(cblknum >= 0) /** IF THIS SNODE IS NOT A ROOT **/
        {
          sontab[sonindex[cblknum]+tmp[cblknum]] = i;
          tmp[cblknum]++;
        }
    }


  /***********************************************************/
  /* Compute the fill to merge a column of P with its father */
  /***********************************************************/
  for(i=0;i<n;i++)
    {
      father = treetab2[i];
      if(father == -1 || father == i)
        {
          nnzadd[i] = INFINI;
          gain[i] = INFINI;
          continue;
        }

      nnzadd[i] = merge_cost(i, father, P, colweight);
      /*#ifdef BLAS_GAIN*/

      /** @@@ OIMBE inutile : blas_gain vaut 0 au depart maintenant **/
      if(blas_gain == 1)
        gain[i] = merge_gain(i, father, P, colweight, tmp2)/nnzadd[i];
      /*#else*/
      else
        gain[i] = (double)(nnzadd[i]);
      /*#endif*/

      /*if(nnzadd[i] == 0)
        fprintf(stderr, "Cblk %ld \n", i);*/
    }

  /*exit(0);*/

  /*****************************************************/
  /** Merge all the columns so it doesn't add fill-in **/
  /*****************************************************/
  toto = 0;
  /*fprintf(stderr, "LOOP \n");*/
  for(i=0;i<n;i++)
    {

#ifdef DEBUG_KASS
      ASSERT(colweight[i] > 0, MOD_KASS);
#endif
      if(colweight[i] != 0 && nnzadd[i] == 0 && gain[i] <= 0)
        {
          father = treetab2[i];

          toto++;
#ifdef DEBUG_KASS
          ASSERT(father > 0 && father != i, MOD_KASS);
#endif
          /** We merge the snode i and its father **/
          merge_col(i, father, P);

          /*fprintf(stderr, "MERGE %ld (%ld) and %ld (%ld) \n", (long)i, (long)colweight[i],
            (long)father, (long)colweight[father]); */

          colweight[father] += colweight[i];
          colweight[i] = 0; /** mark this node as does not exist **/

          /**  we update nnzadd for the father **/
          k = treetab2[father];
          if(k != -1 && k != father)
            {
              nnzadd[father] = merge_cost(father, k, P, colweight);

              if(blas_gain == 1)
                gain[father] = merge_gain(father, k, P, colweight, tmp2)/nnzadd[father];
              else
                gain[father] = (double)(nnzadd[father]);
            }

          /** We update the sons of i now **/
          get_son(i, sonindex, sontab, colweight, &ind, tmp);
          for(j=0;j<ind;j++)
            {
              k = tmp[j];
              treetab2[k] = father;
              nnzadd[k] = merge_cost(k, father, P, colweight);


              /** @@@ OIMBE inutile : blas_gain vaut 0 au depart maintenant **/
              if(blas_gain == 1)
                gain[k] = merge_gain(k, father, P, colweight, tmp2)/nnzadd[k];
              else
                gain[k] = (double)(nnzadd[k]);

#ifdef DEBUG_KASS
              ASSERT(nnzadd[k] > 0, MOD_KASS);
#endif
            }

        }
    }
#ifdef DEBUG_KASS
  print_one("Apres amalgamation init cblk = %ld \n", (long)(n-toto));
#endif
  /*** Put in a sort heap the column sorted by their nnzadd ***/
 debut:
  queueInit(&heap, n);
  for(i=0;i<n;i++)
    if(colweight[i] > 0 && (treetab2[i] > 0 && treetab2[i] != i))
      {
#ifdef DEBUG_KASS
        if(blas_gain != 1)
          {
            if(nnzadd[i] <= 0)
              fprintf(stderr, "nnzadd[%ld] = %ld \n", (long)i, (long)nnzadd[i]);
            ASSERT(nnzadd[i]>0, MOD_KASS);
          }
#endif

        /*#ifdef BLAS_GAIN*/
        if(blas_gain == 1)
          {
            if(gain[i] <= 0)
              queueAdd(&heap, i, gain[i]);
          }
        /*#else*/
        else
          queueAdd(&heap, i, gain[i]);
        /*#endif*/
      }

  /*******************************************/
  /* Merge supernodes until we reach fillmax */
  /*******************************************/
  /*** Merge supernodes untill we reach the fillmax limit ****/
  fill = 0.0;


  while(queueSize(&heap)>0 && fill < fillwhile)
     {
       i = queueGet2(&heap, &key, NULL);


       /*if(nnzadd[i] != (PASTIX_INT)key || colweight[i] <= 0)*/
       if(gain[i] != key || colweight[i] <= 0)
         continue;

       if(fill + nnzadd[i] > fillmax)
         break;
       else
         fill += nnzadd[i];




       toto++;
       father = treetab2[i];
#ifdef DEBUG_KASS
       ASSERT(father > 0 && father != i, MOD_KASS);
       ASSERT(colweight[father]>0, MOD_KASS);
#endif

       /*fprintf(stderr, "Merge col %ld and %ld gain = %g \n", (long)i, (long)father, gain[i]);*/

       /** We merge the snode i and its father and
           we update treetab2, nnzadd, colweight **/
       merge_col(i, father, P);
       colweight[father] += colweight[i];
       colweight[i] = 0; /** mark this node as does not exist **/

       /**  we update nnzadd for the father **/
       k = treetab2[father];
       if(k != -1 && k != father)
         {
           nnzadd[father] = merge_cost(father, k, P, colweight);
           /*#ifdef BLAS_GAIN*/
           if(blas_gain == 1)
             gain[father] = merge_gain(father, k, P, colweight, tmp2)/nnzadd[father];
           /*#else*/
           else
             gain[father] = (double)(nnzadd[father]);
           /*#endif*/

           /*queueAdd(&heap, father, (double) nnzadd[father]);*/
           /*#ifdef BLAS_GAIN*/
           if(blas_gain == 1)
             {
               if(gain[father] <= 0)
                 queueAdd(&heap, father, gain[father]);
             }
           else
             /*#else*/
             queueAdd(&heap, father, gain[father]);
           /*#endif*/
         }

       /** We update the sons of i now **/
       get_son(i, sonindex, sontab, colweight, &ind, tmp);
       for(j=0;j<ind;j++)
         {
           k = tmp[j];
           treetab2[k] = father;
           nnzadd[k] = merge_cost(k, father, P, colweight);
           /*#ifdef BLAS_GAIN*/
           if(blas_gain == 1)
             {
               gain[k] = merge_gain(k, father, P, colweight, tmp2)/nnzadd[k];
               if(gain[k] <= 0)
                 queueAdd(&heap, k, gain[k]);
             }
           /*#else*/
           else
             {
               gain[k] = (double)(nnzadd[k]);
               queueAdd(&heap, k, gain[k]);
             }
           /*#endif*/
         }
     }
#ifdef DEBUG_KASS
  print_one("After amalg phase cblk = %ld fillwhile %ld fillmax %ld \n", (long)toto, (long)fillwhile, (long)fillmax);
#endif
  if(fillwhile < fillmax)
    {

      fillwhile = fillmax;

      /** Now the gain of amalgamation is based on the BLAS model **/
      queueExit(&heap);
#ifdef DEBUG_KASS
      ASSERT(blas_gain == 0, MOD_KASS);
#endif
      blas_gain = 1;

      /** Recompute the gain using BLAS model **/
      for(i=0;i<n;i++)
        {
          father = treetab2[i];
          if(father == -1 || father == i)
            {
              gain[i] = INFINI;
              continue;
            }
          gain[i] = merge_gain(i, father, P, colweight, tmp2)/nnzadd[i];
        }


      goto debut;
    }


  memFree(nnzadd);
  memFree(gain);
  queueExit(&heap);

  /*fprintf(stderr, "FINAL cblk = %ld \n", (long)(n-toto));*/



  /********************************/
  /* Compute the new partition    */
  /********************************/

  /** Count the number of supernodes **/

  /** tmp will be the newnum of node i in the rangtab **/
  newnum = tmp;
  bzero(newnum, sizeof(PASTIX_INT)*n);
  k = 0;
  for(i=0;i<n;i++)
    if(colweight[i] > 0)
      newnum[i] = k++;
  *cblknbr = k;
#ifdef DEBUG_KASS
  print_one("Number of cblk after amal = %ld \n", (long)*cblknbr);
#endif

  MALLOC_INTERN(*rangtab, k+1, PASTIX_INT);
  bzero(*rangtab, sizeof(PASTIX_INT)*(k+1));

  for(i=0;i<n;i++)
    if(colweight[i] > 0)
      (*rangtab)[newnum[i] + 1]+= colweight[i];

  for(i=1;i<= (*cblknbr);i++)
    (*rangtab)[i] += (*rangtab)[i-1];


  for(i=0;i<n;i++)
    {
      if(colweight[i] > 0)
        {
          if(snodetab != NULL)
            for(j=snodetab[i];j<snodetab[i+1];j++)
              nodetab[(*rangtab)[newnum[i]]++] = j;
          else
            nodetab[(*rangtab)[newnum[i]]++] = i;
        }
      else
        {
          /** find the cblk this node is in **/
          father = i;
          while(colweight[father] <= 0)
            {
              father = treetab2[father];
#ifdef DEBUG_KASS
              ASSERT(father > 0, MOD_KASS);
#endif
            }
          if(snodetab != NULL)
            for(j=snodetab[i];j<snodetab[i+1];j++)
              nodetab[(*rangtab)[newnum[father]]++] = j;
          else
            nodetab[(*rangtab)[newnum[father]]++] = i;
        }
    }

  /** reset rangtab to its real value **/
  for(i=*cblknbr; i>0; i--)
    (*rangtab)[i] = (*rangtab)[i-1];
  (*rangtab)[0] = 0;





#ifdef DEBUG_KASS
  /*for(i=0;i<*cblknbr+1;i++)
    fprintf(stderr, "rangtab[%ld] = %ld \n", (long)i, (long)(*rangtab)[i]);
    exit(0);*/
  for(i=0;i<n;i++)
    if(colweight[i] > 0)
      ASSERT(colweight[i] == (*rangtab)[newnum[i]+1]-(*rangtab)[newnum[i]], MOD_KASS);


  /** check the iperm vector (nodetab) **/
 {
   PASTIX_INT *flag;
   fprintf(stderr, "Cblknbr = %ld NN = %ld \n", *cblknbr, (long)nn);
   ASSERT( (*rangtab)[*cblknbr] == nn, MOD_KASS);

   MALLOC_INTERN(flag, nn, PASTIX_INT);


   bzero(flag, sizeof(PASTIX_INT)*nn);
   for(i=0;i<nn;i++)
     {
       ASSERT(nodetab[i] >= 0 && nodetab[i] < nn, MOD_KASS);
       flag[nodetab[i]]++;
     }
   for(i=0;i<nn;i++)
     {
       if(flag[nodetab[i]] != 1)
         fprintf(stderr, "(Nodetab[%ld] = %ld falg = %ld ) ", (long)i, (long)nodetab[i], (long)flag[nodetab[i]]);
       ASSERT(flag[nodetab[i]] == 1, MOD_KASS);
     }

   memFree(flag);
 }
#endif

  /** Compact P: each column from 1 to cblknbr will represent a cblk ***/
  ind = 0;
  while(ind < P->n && P->nnzrow[ind] >0 )
    ind++;

  for(i=ind;i<n;i++)
    if(P->nnzrow[i] > 0)
      {
#ifdef DEBUG_KASS
        ASSERT(colweight[i] > 0, MOD_KASS);
#endif
        P->nnzrow[ind] = P->nnzrow[i];
        P->ja[ind] = P->ja[i];
        P->nnzrow[i] = 0;
        P->ja[i] = NULL;
        ind++;
      }
  P->n = *cblknbr;

  /*** Apply the new permutation to P ****/
  /** tmp is the perm vector **/
  for(i=0;i<nn;i++)
    tmp[nodetab[i]] = i;

  for(i=0;i<P->n;i++)
    {
      PASTIX_INT *ja;
      ja = P->ja[i];
      for(j=0;j<P->nnzrow[i];j++)
        ja[j] = tmp[ja[j]];
    }


#ifdef DEBUG_KASS
  /** Check some things about P **/
  for(i=0;i<P->n;i++)
    {
      if(P->nnzrow[i] < (*rangtab)[i+1]-(*rangtab)[i])
        for(j=0;j< P->nnzrow[i];j++)
          fprintf(stderr, "ja = %ld, rang = %ld \n", (long)P->ja[i][j],(long)(*rangtab)[i]+j);

      ASSERT(P->nnzrow[i] >= (*rangtab)[i+1]-(*rangtab)[i], MOD_KASS);
      for(j=0;j< (*rangtab)[i+1]-(*rangtab)[i];j++)
        {
          if(P->ja[i][j] != (*rangtab)[i]+j)
            fprintf(stderr, "Cblk %ld j %ld ja %ld rangtab[%ld]=%ld rangtab[%ld] = %ld \n",
                    i, j, P->ja[i][j], i, (*rangtab)[i], i+1, (*rangtab)[i+1]);
          ASSERT(P->ja[i][j] == (*rangtab)[i]+j, MOD_KASS);
        }

      /** The matrix should be still sorted **/
      for(j=1;j<P->nnzrow[i];j++)
        ASSERT(P->ja[i][j] > P->ja[i][j-1], MOD_KASS);

    }

  /*for(i=0;i<nn;i++)
    fprintf(stderr, "%ld ", nodetab[i]);
    fprintf(stderr, "\n");*/

#endif

#ifdef DEBUG_KASS
  key = 0.0;
  for(i=0;i<P->n;i++)
    key += CBLKTIME(P->nnzrow[i], P->ja[i], (*rangtab)[i+1]-(*rangtab)[i]);
  fprintf(stderr, "COST of the AMALGAMATED MATRIX = %g \n", key);
#endif




  memFree(tmp);
  memFree(tmp2);
  memFree(sonindex);
  memFree(sontab);
  memFree(treetab2);
  memFree(colweight);

}



PASTIX_INT merge_cost(PASTIX_INT a, PASTIX_INT b, csptr P, PASTIX_INT *colweight)
{
  PASTIX_INT i1, n1, i2, n2;
  PASTIX_INT *ja1, *ja2;
  PASTIX_INT cost;

  ja1 = P->ja[a];
  ja2 = P->ja[b];


  n1 = P->nnzrow[a];
  n2 = P->nnzrow[b];

  i1 = i2 = 0;
  /** The diagonal elements of row a does not create fill-in **/
  while(i1 < n1 && ja1[i1] < ja2[0])
    i1++;


  /*fprintf(stderr, "MERGECOST %ld (%ld) + %ld (%ld)  i1 = %ld \n", a, n1, b, n2, i1);*/

  cost = 0;
  while(i1 < n1 && i2 < n2)
    {
      if(ja1[i1] < ja2[i2])
        {
          cost += colweight[b];
          /*fprintf(stderr, "TOTO cost1 %ld \n", cost);*/
          i1++;
          continue;
        }
      if(ja1[i1] > ja2[i2])
        {
          cost += colweight[a];
          /*fprintf(stderr, "TOTO cost2 %ld \n", cost);*/
          i2++;
          continue;
        }

      /** ja1[i1] == ja2[i2] **/
      i1++;
      i2++;
    }

  while(i1 < n1)
    {
      cost += colweight[b];
      /*fprintf(stderr, "TOTO costR1 %ld \n", cost);*/
      i1++;
      continue;
    }

  while(i2 < n2)
    {
      cost += colweight[a];
      /*fprintf(stderr, "TOTO costR2 %ld \n", cost);*/
      i2++;
      continue;
    }
#ifdef DEBUG_KASS
  ASSERT(cost >= 0, MOD_KASS);
#endif

  return cost;
}


void get_son(PASTIX_INT node, PASTIX_INT *sonindex, PASTIX_INT *sontab, PASTIX_INT *colweight, PASTIX_INT *ns, PASTIX_INT *list)
{
  PASTIX_INT i, s;
  PASTIX_INT nss;
  PASTIX_INT ind;
  ind = 0;
  for(i=sonindex[node];i<sonindex[node+1];i++)
    {
      s = sontab[i];
      if(colweight[s] <= 0)
        {
          get_son(s, sonindex, sontab, colweight, &nss, list+ind);
          ind += nss;
        }
      else
        list[ind++] = s;
    }
  *ns = ind;

}


void  merge_col(PASTIX_INT a, PASTIX_INT b, csptr P)
{
  PASTIX_INT i, i1, i2;
  PASTIX_INT n1, n2;
  PASTIX_INT *ja1, *ja2;
  PASTIX_INT *ja = NULL;


  ja1 = P->ja[a];
  ja2 = P->ja[b];
  n1 = P->nnzrow[a];
  n2 = P->nnzrow[b];


  MALLOC_INTERN(ja, n1+n2, PASTIX_INT);

  i1 = 0;
  i2 = 0;
  i = 0;

  /*fprintf(stderr, "MERGE %ld and %ld  \n", (long)a, (long)b);*/

  while(i1 < n1 && i2 < n2)
    {
      /*fprintf(stderr, "i1 %ld i2 %ld n1 %ld n2 %ld \n", (long)i1, (long)i2, (long)n1, (long)n2);*/
      if(ja1[i1] < ja2[i2])
        {
          ja[i] = ja1[i1];
          i1++;
          i++;
          continue;
        }

      if(ja1[i1] > ja2[i2])
        {
          ja[i] = ja2[i2];
          i2++;
          i++;
          continue;
        }

      ja[i] = ja1[i1];
      i++;
      i1++;
      i2++;
    }
  /*fprintf(stderr, "DONE LOPP \n");
    fprintf(stderr, "DONE i1 %ld i2 %ld n1 %ld n2 %ld \n", (long)i1, (long)i2, (long)n1, (long)n2);*/

#ifdef DEBUG_KASS
  assert(i1 == n1 || i2 == n2);
#endif

  for(;i1<n1;i1++)
    ja[i++] = ja1[i1];

  for(;i2<n2;i2++)
    ja[i++] = ja2[i2];

  /*fprintf(stderr, "E1 \n");*/

#ifdef DEBUG_KASS
  assert(i >= n1 && i >= n2);
  assert(i <= n1+n2);
#endif


  if(P->ja[a] != NULL)
    {
      P->nnzrow[a] = 0;
      memFree(P->ja[a]);
    }

  if(P->ja[b] != NULL)
    memFree(P->ja[b]);

  P->nnzrow[b] = i;
  /*fprintf(stderr, "E2: realloc i = %ld  \n", (long)i);*/


  /*P->ja[b] = (PASTIX_INT *)realloc(ja, sizeof(PASTIX_INT)*i);*/
  MALLOC_INTERN(P->ja[b], i, PASTIX_INT);
  memCpy(P->ja[b], ja, sizeof(PASTIX_INT)*i);
  memFree(ja);

  /*fprintf(stderr, "DONE2 i1 %ld i2 %ld n1 %ld n2 %ld \n", (long)i1, (long)i2, (long)n1, (long)n2);*/

}


double merge_gain(PASTIX_INT a, PASTIX_INT b, csptr P, PASTIX_INT *colweight, PASTIX_INT *tmp)
{
  double costa, costb, costm;
  PASTIX_INT nm;

  costa = CBLKTIME(P->nnzrow[a], P->ja[a], colweight[a]);
  costb = CBLKTIME(P->nnzrow[b], P->ja[b], colweight[b]);

  UnionSet(P->ja[a], P->nnzrow[a], P->ja[b], P->nnzrow[b], tmp, &nm);

  costm = CBLKTIME(nm, tmp, colweight[a] + colweight[b]);

  /*fprintf(stderr, "Cost(%ld) = %g cost(%ld) = %g costm = %g \n", (long)a, costa, (long)b, costb, costm);*/

#ifdef DEBUG_KASS
  /*if(costm > costa + costb)
    fprintf(stderr, "BLAS negative column %ld merge column %ld \n", (long)a, (long)b);*/
#endif

  return costm - costa - costb;

}


double cblk_time_fact(PASTIX_INT n, PASTIX_INT *ja, PASTIX_INT colnbr)
{
  /*******************************************/
  /* Compute the time to compute a cblk      */
  /* according to the BLAS modelization      */
  /*******************************************/
  double cost;
  PASTIX_INT i;
  PASTIX_INT L, G, H;

  /** The formula are based on the costfunc.c in blend **/
  /** @@@ OIMBE: il faudra faire les DOF_CONSTANT ***/

  /** Diagonal factorization and TRSM **/
  L = colnbr;
  G = n-L;
#define CHOLESKY
#ifndef CHOLESKY
  cost =(double)(L*PERF_COPY(L)+ PERF_PPF(L) + PERF_TRSM(L, G) + L*PERF_SCAL(G)
                  + L*PERF_COPY(G));
#else
  cost = (double)(PERF_POF(L) + PERF_TRSM(L, G)) ;
#endif

  /** Contributions **/
  i = colnbr;
  while(i<n)
    {
      H = 1;
      i++;
      while(i<n && ja[i] == ja[i-1]+1)
        {
          i++;
          H++;
        }

      cost += (double)(PERF_GEMM(G, H, L));
      G -= H;

    }

  return cost;
}


double cblk_time_solve(PASTIX_INT n, PASTIX_INT *ja, PASTIX_INT colnbr)
{
  double cost;
  PASTIX_INT L;
  (void)ja;

  L = colnbr;

  cost = (double)PERF_TRSV(L) + (double) PERF_GEMV(L, n-L);
  return cost;
}
