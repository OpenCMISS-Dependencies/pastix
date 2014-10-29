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
/**   NAME       : kass.c                                  **/
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

#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#include "common_pastix.h"
#include "dof.h"
#include "symbol.h"
/* #include "symbol_cost.h" */
#include "order.h"
#include "fax.h"
/* #include "ifax.h" */
#include "sparRow.h"
#include "sort_row.h"
#include "SF_level.h"
#include "SF_Direct.h"

#include "compact_graph.h"
#include "amalgamate.h"
#include "find_supernodes.h"
#include "KSupernodes.h"
#include "kass.h"

#define print_one(fmt, ...)    if( procnum == 0) fprintf(stdout, fmt, __VA_ARGS__)

/*#define KS*/ /** Decommenter ca pour remettre ancienne version de
                   KASS (amalgamation a l'interieur des supernoeuds ! **/


/*#define SCOTCH_SNODE*/ /*CA MARCHE PAS POUR ILUK  */


extern double nnz(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);
extern double recursive_sum(PASTIX_INT a, PASTIX_INT b,
                            double (*fval)(PASTIX_INT, const SymbolMatrix *, const Dof *),
                            const SymbolMatrix * symbmtx, const Dof * dofptr);

void kass_symbol(csptr mat, PASTIX_INT levelk, double rat,
                  PASTIX_INT *perm, PASTIX_INT *iperm, PASTIX_INT snodenbr, PASTIX_INT *snodetab, PASTIX_INT *streetab,
                  PASTIX_INT *cblknbr, PASTIX_INT **rangtab, SymbolMatrix *symbmtx, MPI_Comm pastix_comm);
void Build_SymbolMatrix(csptr P, PASTIX_INT cblknbr, PASTIX_INT *rangtab, SymbolMatrix *symbmtx);
void Patch_SymbolMatrix(SymbolMatrix *symbmtx);

void kass(int            levelk,
          int            rat,
          SymbolMatrix * symbptr,
          PASTIX_INT            baseval,
          PASTIX_INT            vertnbr,
          PASTIX_INT            edgenbr,
          PASTIX_INT          * verttab,
          PASTIX_INT          * edgetab,
          Order        * orderptr,
          MPI_Comm       pastix_comm)
{
  PASTIX_INT snodenbr;
  PASTIX_INT *snodetab   = NULL;
  PASTIX_INT *treetab    = NULL;
  PASTIX_INT *ia         = NULL;
  PASTIX_INT *ja         = NULL;
  PASTIX_INT i, j, n;
  PASTIX_INT ind;
  csptr mat;
  PASTIX_INT *tmpj       = NULL;
  PASTIX_INT *perm       = NULL;
  PASTIX_INT *iperm      = NULL;
  PASTIX_INT newcblknbr;
  PASTIX_INT *newrangtab = NULL;
  Dof dofstr;
  Clock timer1;
  double nnzS;
  int procnum;
  (void)edgenbr;

  MPI_Comm_rank(pastix_comm,&procnum);

#ifdef DEBUG_KASS
  print_one("--- kass begin ---\n");
#endif
/*   graphData (graphptr,  */
/*           (SCOTCH_Num * )&baseval,  */
/*           (SCOTCH_Num * )&vertnbr,  */
/*           (SCOTCH_Num **)&verttab,  */
/*           NULL, NULL, NULL,  */
/*           (SCOTCH_Num * )&edgenbr,  */
/*           (SCOTCH_Num **)&edgetab,  */
/*           NULL); */

  n = vertnbr;
  ia = verttab;
  ja = edgetab;
  perm = orderptr->permtab;
  iperm = orderptr->peritab;

  /*** Convert Fortran to C numbering ***/
  if(baseval == 1)
    {
      for(i=0;i<=n;i++)
          ia[i]--;
      for(i=0;i<n;i++)
        for(j=ia[i];j<ia[i+1];j++)
          ja[j]--;
      for(i=0;i<n;i++)
        orderptr->permtab[i]--;
      for(i=0;i<n;i++)
        orderptr->peritab[i]--;
    }

  MALLOC_INTERN(treetab, n, PASTIX_INT);
#ifndef SCOTCH_SNODE
  /*if(rat != -1 )*/
    {
      /***** FIND THE SUPERNODE PARTITION FROM SCRATCH ********/

      /*** Find the supernodes of the direct factorization  ***/
      MALLOC_INTERN(snodetab, n+1, PASTIX_INT);


      clockInit(&timer1);
      clockStart(&timer1);
      find_supernodes(n, ia, ja, perm, iperm, &snodenbr, snodetab, treetab);
      clockStop(&timer1);
      print_one("Time to find the supernode (direct) %.3g s \n", clockVal(&timer1));

      /*memfree(treetab);*/
      print_one("Number of supernode for direct factorization %ld \n", (long)snodenbr);
    }
#else
  /*else*/
    {
      /***** USE THE SUPERNODE PARTITION OF SCOTCH  ********/
      snodenbr = orderptr->cblknbr;
      MALLOC_INTERN(snodetab, n+1, PASTIX_INT);
      memCpy(snodetab, orderptr->rangtab, sizeof(PASTIX_INT)*(snodenbr+1));
      print_one("Number of column block found in scotch (direct) %ld \n", (long)snodenbr);

    }
#endif

  /****************************************/
  /*  Convert the graph                   */
  /****************************************/
    MALLOC_INTERN(mat, 1, struct SparRow);
  initCS(mat, n);
  MALLOC_INTERN(tmpj, n, PASTIX_INT);
  /**** Convert and permute the matrix in sparrow form  ****/
  /**** The diagonal is not present in the CSR matrix, we have to put it in the matrix ***/
  bzero(tmpj, sizeof(PASTIX_INT)*n);
  for(i=0;i<n;i++)
    {
      /*** THE GRAPH DOES NOT CONTAIN THE DIAGONAL WE ADD IT ***/
      tmpj[0] = i;
      ind = 1;
      for(j=ia[i];j<ia[i+1];j++)
        tmpj[ind++] = ja[j];

      mat->nnzrow[i] = ind;
      MALLOC_INTERN(mat->ja[i], ind, PASTIX_INT);
      memCpy(mat->ja[i], tmpj, sizeof(PASTIX_INT)*ind);
      mat->ma[i] = NULL;
    }
  CS_Perm(mat, perm);
  /*** Reorder the matrix ***/
  sort_row(mat);
  memFree(tmpj);


  /***** COMPUTE THE SYMBOL MATRIX OF ILU(K) WITH AMALGAMATION *****/
  kass_symbol(mat, levelk, (double)(rat)/100.0, perm,
              iperm, snodenbr, snodetab, treetab, &newcblknbr, &newrangtab,
              symbptr, pastix_comm);


  cleanCS(mat);
  memFree(mat);
  memFree(treetab);

  dofInit(&dofstr);
  dofConstant(&dofstr, 0, symbptr->nodenbr, 1);
  nnzS =  recursive_sum(0, symbptr->cblknbr-1, nnz, symbptr, &dofstr);
  print_one("Number of non zero in the non patched symbol matrix = %g, fillrate1 %.3g \n",
            nnzS+n, (nnzS+n)/(ia[n]/2.0 +n));
  dofExit(&dofstr);


  if(symbolCheck(symbptr) != 0)
    {
      errorPrint("SymbolCheck after kass_symbol.");
      ASSERT(0, MOD_KASS);
    }



  if(levelk != -1)
    {
      /********************************************************/
      /** ADD BLOCKS IN ORDER TO GET A REAL ELIMINATION TREE **/
      /********************************************************/
      Patch_SymbolMatrix(symbptr);
    }



  dofInit(&dofstr);
  dofConstant(&dofstr, 0, symbptr->nodenbr, 1);
  nnzS =  recursive_sum(0, symbptr->cblknbr-1, nnz, symbptr, &dofstr);

  dofExit(&dofstr);
  print_one("Number of block in final symbol matrix = %ld \n", (long)symbptr->bloknbr);
  print_one("Number of non zero in final symbol matrix = %g, fillrate2 %.3g \n",  nnzS+n, (nnzS+n)/(ia[n]/2.0 +n));
  if(symbolCheck(symbptr) != 0)
    {
      errorPrint("SymbolCheck after Patch_SymbolMatrix.");
      ASSERT(0, MOD_KASS);
    }
#ifdef DEBUG_KASS
  print_one("--- kass end ---\n");
#endif
  memFree(snodetab);
  orderptr->cblknbr = newcblknbr;
  memFree(orderptr->rangtab);
  orderptr->rangtab = newrangtab;

}

void kass_symbol(csptr mat, PASTIX_INT levelk, double rat, PASTIX_INT *perm, PASTIX_INT *iperm, PASTIX_INT snodenbr, PASTIX_INT *snodetab, PASTIX_INT *streetab, PASTIX_INT *cblknbr, PASTIX_INT **rangtab, SymbolMatrix *symbmtx, MPI_Comm pastix_comm)
{
  /**************************************************************************************/
  /* This function computes a symbolic factorization ILU(k) given a CSR matrix and an   */
  /* ordering. Then it computes a block partition of the factor to get BLAS3            */
  /* efficiency                                                                         */
  /* NOTE: the CSC matrix is given symmetrized and without the diagonal                 */
  /**************************************************************************************/

  PASTIX_INT i, j;
  PASTIX_INT nnzL;
  PASTIX_INT *iperm2  = NULL;
  PASTIX_INT *treetab = NULL;
  PASTIX_INT n;
  csptr P;
  Clock timer1;
  int procnum;

  MPI_Comm_rank(pastix_comm,&procnum);

  n = mat->n;
  MALLOC_INTERN(iperm2, n, PASTIX_INT);

  /*compact_graph(mat, NULL, NULL, NULL);*/

  /*** Compute the ILU(k) pattern of the quotient matrix ***/
  MALLOC_INTERN(P, 1, struct SparRow);
  initCS(P, n);
  print_one("Level of fill = %ld\nAmalgamation ratio = %d \n", (long)levelk, (int)(rat*100));
  clockInit(&timer1);
  clockStart(&timer1);

  if(levelk == -1)
    {

      /***** FACTORISATION DIRECT *******/
      /***** (Re)compute also the streetab (usefull when SCOTCH_SNODE
             is active) ***/
      SF_Direct(mat, snodenbr, snodetab, streetab, P);

      clockStop(&timer1);
      print_one("Time to compute scalar symbolic direct factorization  %.3g s \n", clockVal(&timer1));
#ifdef DEBUG_KASS
      print_one("non-zeros in P = %ld \n", (long)CSnnz(P));
#endif
      nnzL = 0;
      for(i=0;i<P->n;i++)
        {
          PASTIX_INT ncol;
          ncol = snodetab[i+1]-snodetab[i];
          nnzL += (ncol*(ncol+1))/2;
#ifdef DEBUG_KASS
          ASSERT(P->nnzrow[i] >= ncol, MOD_KASS);
          if(P->nnzrow[i] >= n)
            fprintf(stderr,"P->nnzrow[%ld] = %ld \n", (long)i, (long)P->nnzrow[i]);
          ASSERT(P->nnzrow[i] < n, MOD_KASS);
#endif
          nnzL += (P->nnzrow[i]-ncol)*ncol;
        }
#ifdef DEBUG_KASS
      print_one("NNZL = %ld \n", (long)nnzL);
#endif
    }
  else
    {
      /***** FACTORISATION INCOMPLETE *******/
      nnzL = SF_level(2, mat, levelk, P);

      clockStop(&timer1);
      print_one("Time to compute scalar symbolic factorization of ILU(%ld) %.3g s \n",
              (long)levelk, clockVal(&timer1));

    }
  print_one("Scalar nnza = %ld nnzlk = %ld, fillrate0 = %.3g \n",
            (long)( CSnnz(mat) + n)/2, (long)nnzL, (double)nnzL/(double)( (CSnnz(mat)+n)/2.0 ));



  /** Sort the rows of the symbolic matrix */
  sort_row(P);

  clockInit(&timer1);
  clockStart(&timer1);

  if(levelk != -1)
    {

      /********************************/
      /** Compute the "k-supernodes" **/
      /********************************/

#ifdef KS
      assert(levelk >= 0);
      KSupernodes(P, rat, snodenbr, snodetab, cblknbr, rangtab);
#else

#ifdef SCOTCH_SNODE
      if(rat == -1)
        assert(0); /** do not have treetab with this version of Scotch **/
#endif

      MALLOC_INTERN(treetab, P->n, PASTIX_INT);
      for(j=0;j<snodenbr;j++)
        {
          for(i=snodetab[j];i<snodetab[j+1]-1;i++)
            treetab[i] = i+1;

          /*** Version generale ****/
          if(streetab[j] == -1 || streetab[j] == j)
            treetab[i] = -1;
          else
            treetab[i]=snodetab[streetab[j]];
          /*** Version restricted inside the supernode (like KSupernodes) ***/
          /*treetab[snodetab[j+1]-1] = -1;*/  /** this should give the same results than
                                                  KSupernodes **/
        }

      /** NEW ILUK + DIRECT **/
      amalgamate(rat, P, -1, NULL, treetab, cblknbr, rangtab, iperm2, pastix_comm);

      memFree(treetab);
      for(i=0;i<n;i++)
        iperm2[i] = iperm[iperm2[i]];
      memcpy(iperm, iperm2, sizeof(PASTIX_INT)*n);
      for(i=0;i<n;i++)
        perm[iperm[i]] = i;
#endif
    }
  else{


    /*if(0)*/
      {
        amalgamate(rat, P, snodenbr, snodetab, streetab, cblknbr,
                   rangtab, iperm2, pastix_comm);

        /** iperm2 is the iperm vector of P **/
        for(i=0;i<n;i++)
          iperm2[i] = iperm[iperm2[i]];
        memcpy(iperm, iperm2, sizeof(PASTIX_INT)*n);
        for(i=0;i<n;i++)
          perm[iperm[i]] = i;
      }
      /*else
      {
        fprintf(stderr, "RAT = 0 SKIP amalgamation \n");
        *cblknbr = snodenbr;
        MALLOC_INTERN(*rangtab, snodenbr+1, PASTIX_INT);
        memcpy(*rangtab, snodetab, sizeof(PASTIX_INT)*(snodenbr+1));
        }*/
  }

  clockStop(&timer1);
  print_one("Time to compute the amalgamation of supernodes %.3g s\n", clockVal(&timer1));

  print_one("Number of cblk in the amalgamated symbol matrix = %ld \n", (long)*cblknbr);


  Build_SymbolMatrix(P, *cblknbr, *rangtab, symbmtx);


  print_one("Number of block in the non patched symbol matrix = %ld \n", (long)symbmtx->bloknbr);


  memFree(iperm2);
  cleanCS(P);
  memFree(P);

}





void  Build_SymbolMatrix(csptr P, PASTIX_INT cblknbr, PASTIX_INT *rangtab, SymbolMatrix *symbmtx)
{
  PASTIX_INT i, j, k, l;
  PASTIX_INT cblknum;
  PASTIX_INT ind;
  PASTIX_INT    *tmpj      = NULL;
  double *tmpa      = NULL;
  PASTIX_INT    *node2cblk = NULL;
  PASTIX_INT    *ja        = NULL;
  PASTIX_INT n;

  n = rangtab[cblknbr];

  /**** First we transform the P matrix to find the block ****/
  MALLOC_INTERN(tmpj, n, PASTIX_INT);
  MALLOC_INTERN(tmpa, n, double);
  MALLOC_INTERN(node2cblk, n, PASTIX_INT);


  for(k=0;k<cblknbr;k++)
    for(i=rangtab[k];i<rangtab[k+1];i++)
      node2cblk[i] = k;

  for(k=0;k<cblknbr;k++)
    {
      /*i = rangtab[k];*/ /*OLD VERSION QUAND P pas recompacte */
      i = k;
#ifdef DEBUG_KASS
      ASSERT(P->nnzrow[i] >= (rangtab[k+1]-rangtab[k]), MOD_KASS);

      for(l=0;l<rangtab[k+1]-rangtab[k];l++)
        {
          ASSERT(P->ja[i][l] == rangtab[k]+l, MOD_KASS);
          ASSERT(node2cblk[P->ja[i][l]] == i, MOD_KASS);
        }
#endif
      ja = P->ja[i];
      j = 0;
      ind = 0;
      while(j<P->nnzrow[i])
        {
          cblknum = node2cblk[ja[j]];
          l=j+1;
          while(l < P->nnzrow[i] && ja[l] == ja[l-1]+1 && node2cblk[ja[l]] == cblknum)
            l++;

          tmpj[ind] = ja[j];
          tmpa[ind] = (double)(l-j);
          j = l;
          ind++;
        }

      memFree(P->ja[i]);
      P->nnzrow[i] = ind;
      MALLOC_INTERN(P->ja[i], ind, PASTIX_INT);
      MALLOC_INTERN(P->ma[i], ind, double);
      memCpy(P->ja[i], tmpj, sizeof(PASTIX_INT)*ind);
      memCpy(P->ma[i], tmpa, sizeof(double)*ind);

    }


  memFree(tmpj);
  memFree(tmpa);

#ifdef DEBUG_KASS
  for(k=0;k<cblknbr;k++)
    {
      /*i = rangtab[k];*/
      i = k;
      assert(P->nnzrow[i] > 0);

      if(P->ma[i][0] != (double)(rangtab[k+1]-rangtab[k]))
        print_one("Cblk %ld ma %ld rg %ld \n", k, (PASTIX_INT)P->ma[i][0],rangtab[k+1]-rangtab[k]);

      assert(P->ma[i][0] == (double)(rangtab[k+1]-rangtab[k]));
    }
#endif

  /**********************************/
  /*** Compute the symbol matrix ****/
  /**********************************/
  symbmtx->baseval = 0;
  symbmtx->cblknbr = cblknbr;

  ind = 0;
  symbmtx->bloknbr = CSnnz(P);
  symbmtx->nodenbr = rangtab[cblknbr];

  MALLOC_INTERN(symbmtx->cblktab, cblknbr+1,        SymbolCblk);
  MALLOC_INTERN(symbmtx->bloktab, symbmtx->bloknbr, SymbolBlok);

  ind = 0;
  for(k=0;k<cblknbr;k++)
    {
      symbmtx->cblktab[k].fcolnum = rangtab[k];
      symbmtx->cblktab[k].lcolnum = rangtab[k+1]-1;
      symbmtx->cblktab[k].bloknum = ind;
      /*l = rangtab[k];*/ /** OLD VERSION **/
      l = k;
      for(i=0;i<P->nnzrow[l];i++)
        {
          j = P->ja[l][i];
          symbmtx->bloktab[ind].frownum = j;
          symbmtx->bloktab[ind].lrownum = j+(PASTIX_INT)(P->ma[l][i])-1;
          symbmtx->bloktab[ind].cblknum = node2cblk[j];
          symbmtx->bloktab[ind].levfval = 0;
          ind++;
        }
#ifdef DEBUG_KASS
      assert(symbmtx->bloktab[symbmtx->cblktab[k].bloknum].frownum == symbmtx->cblktab[k].fcolnum);
      assert(symbmtx->bloktab[symbmtx->cblktab[k].bloknum].lrownum == symbmtx->cblktab[k].lcolnum);
      assert(symbmtx->bloktab[symbmtx->cblktab[k].bloknum].cblknum == k);
#endif


    }
  /*  virtual cblk to avoid side effect in the loops on cblk bloks */
  symbmtx->cblktab[cblknbr].fcolnum = symbmtx->cblktab[cblknbr-1].lcolnum+1;
  symbmtx->cblktab[cblknbr].lcolnum = symbmtx->cblktab[cblknbr-1].lcolnum+1;
  symbmtx->cblktab[cblknbr].bloknum = ind;

#ifdef DEBUG_KASS
  if(ind != symbmtx->bloknbr)
    fprintf(stderr, "ind %ld bloknbr %ld \n", ind, symbmtx->bloknbr);
  assert(ind == symbmtx->bloknbr);
#endif


  memFree(node2cblk);
}


void Patch_SymbolMatrix(SymbolMatrix *symbmtx)
{
  PASTIX_INT i,j, k;
  PASTIX_INT vroot;
  PASTIX_INT        *father     = NULL; /** For the cblk of the symbol matrix **/
  SymbolBlok *newbloktab = NULL;
  SymbolCblk *cblktab    = NULL;
  SymbolBlok *bloktab    = NULL;
  csptr Q;


  cblktab = symbmtx->cblktab;
  bloktab = symbmtx->bloktab;

  MALLOC_INTERN(father, symbmtx->cblknbr, PASTIX_INT);
  MALLOC_INTERN(newbloktab, symbmtx->bloknbr + symbmtx->cblknbr, SymbolBlok);

  MALLOC_INTERN(Q, 1, struct SparRow);
  initCS(Q, symbmtx->cblknbr);

  /* Count how many extra-diagonal bloks are facing each diagonal blok
   */
  for(i=0;i<symbmtx->cblknbr;i++)
    for(j=cblktab[i].bloknum+1;j<cblktab[i+1].bloknum;j++)
      Q->nnzrow[bloktab[j].cblknum]++;

  /* Allocate nFacingBlok integer for each diagonal blok */
  for(i=0;i<symbmtx->cblknbr;i++)
    {
      MALLOC_INTERN(Q->ja[i], Q->nnzrow[i], PASTIX_INT);
      Q->ma[i] = NULL;
    }

  for(i=0;i<symbmtx->cblknbr;i++)
    Q->nnzrow[i] = 0;

  /* Q->ja[k] will contain, for each extra-diagonal facing blok
   * of the column blok k, its column blok.
   */
  for(i=0;i<symbmtx->cblknbr;i++)
    for(j=cblktab[i].bloknum+1;j<cblktab[i+1].bloknum;j++)
      {
        k = bloktab[j].cblknum;
        Q->ja[k][Q->nnzrow[k]++] = i;
      }

  for(i=0;i<Q->n;i++)
    father[i] = -1;

  for(i=0;i<Q->n;i++)
    {
      /* for each blok facing diagonal blok i,
       * belonging to column blok k.
       *
       */
      for(j=0;j<Q->nnzrow[i];j++)
        {
          k = Q->ja[i][j];
#ifdef DEBUG_KASS
          assert(k<i);
#endif
          vroot = k;
          while(father[vroot] != -1 && father[vroot] != i)
            vroot = father[vroot];
          father[vroot] = i;

        }
    }

  for(i=0;i<Q->n;i++)
    if(father[i] == -1)
      father[i]=i+1;

  cleanCS(Q);
  memFree(Q);



  k = 0;
  for(i=0;i<symbmtx->cblknbr-1;i++)
    {
      PASTIX_INT odb, fbloknum;

      fbloknum = cblktab[i].bloknum;
      memCpy(newbloktab+k, bloktab + fbloknum, sizeof(SymbolBlok));
      cblktab[i].bloknum = k;
      k++;
      odb = cblktab[i+1].bloknum-fbloknum;
      if(odb <= 1 || bloktab[fbloknum+1].cblknum != father[i])
        {
          /** Add a blok toward the father **/
          newbloktab[k].frownum = cblktab[ father[i] ].fcolnum;
          newbloktab[k].lrownum = cblktab[ father[i] ].fcolnum; /** OIMBE try lcolnum **/
          newbloktab[k].cblknum = father[i];
#ifdef DEBUG_KASS
          if(father[i] != i)
            assert(cblktab[father[i]].fcolnum > cblktab[i].lcolnum);
#endif

          newbloktab[k].levfval = 0;
          k++;
        }


      if( odb > 1)
        {
          memCpy(newbloktab +k, bloktab + fbloknum+1, sizeof(SymbolBlok)*(odb-1));
          k+=odb-1;
        }

    }
  /** Copy the last one **/
  memCpy(newbloktab+k, bloktab + symbmtx->cblktab[symbmtx->cblknbr-1].bloknum, sizeof(SymbolBlok));
  cblktab[symbmtx->cblknbr-1].bloknum = k;
  k++;
  /** Virtual cblk **/
  symbmtx->cblktab[symbmtx->cblknbr].bloknum = k;

#ifdef DEBUG_KASS
  assert(k >= symbmtx->bloknbr);
  assert(k < symbmtx->cblknbr+symbmtx->bloknbr);
#endif
  symbmtx->bloknbr = k;
  memFree(symbmtx->bloktab);
  MALLOC_INTERN(symbmtx->bloktab, k, SymbolBlok);
  memCpy( symbmtx->bloktab, newbloktab, sizeof(SymbolBlok)*symbmtx->bloknbr);
  /*  virtual cblk to avoid side effect in the loops on cblk bloks */
  cblktab[symbmtx->cblknbr].bloknum = k;

  memFree(father);
  memFree(newbloktab);
}
