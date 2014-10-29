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
  File: bordi.c

  This is the test module for incomplete
  block ordering strategies.

  Authors:
    Francois PELLEGRINI - .

  Dates:
    Version 0.0 - from : 20 oct 1998
                  to     27 oct 1998
    Version 0.1 - from : 01 jun 1999
                  to     01 jun 1999
    Version 1.0 - from : 05 jul 2002
                  to     05 jul 2002
    Version 1.2 - from : 28 aug 2002
                  to     15 sep 2003
*/

/*
**  The defines and includes.
*/

#define ALPHA_LEVFK

#include <assert.h>
#include "common_pastix.h"
#ifdef WITH_SCOTCH
#ifdef DISTRIBUTED
#include <mpi.h>
#include "ptscotch.h"
#else
#include "scotch.h"
#endif
#endif
#include "dof.h"
#include "symbol.h"
#include "order.h"
#include "fax.h"
#include "ifax.h"
#include "bordi.h"

/*
** Section: Defines
*/

/*
   Define: C_FLAGSMBOUT
   UNUSED define.

   Define: C_FLAGDRWOUT
   UNUSED define.

   Define: C_FLAGPLTOUT
   UNUSED define.

   Define: C_FLAGORDOUT
   UNUSED define.

   Define: C_FLAGSCALAR
   UNUSED define.
*/
#define C_FLAGSMBOUT                1
#define C_FLAGDRWOUT                2
#define C_FLAGPLTOUT                4
#define C_FLAGORDOUT                8
#define C_FLAGSCALAR                16
/* Define: C_FILENBR
   Number of files in list
*/
#define C_FILENBR                   6
/* Define: C_FILEARGNBR
   Number of files which can be arguments
*/
#define C_FILEARGNBR                2
/* Define: C_filenamesrcinp
   Source graph input file name.
   UNUSED
*/
#define C_filenamesrcinp            C_fileTab[0].name
/* Define: C_filenamelogout
   Log file name
   UNUSED
*/
#define C_filenamelogout            C_fileTab[1].name
/* Define: C_filenamesmbout
   Symbol matrix file name.
   UNUSED
*/
#define C_filenamesmbout            C_fileTab[2].name
/* Define: C_filenamedrwout
   Drawing file name.
   UNUSED
 */
#define C_filenamedrwout            C_fileTab[3].name
/* Define: C_filenamepltout
   Plot file name.
   UNUSED
*/
#define C_filenamepltout            C_fileTab[4].name
/* Define: C_filenameordout
   Ordering file name.
   UNUSED
*/
#define C_filenameordout            C_fileTab[5].name
/* Define: C_filepntrsrcinp
   Source graph input file.
   UNUSED
*/
#define C_filepntrsrcinp            C_fileTab[0].pntr
/* Define: C_filepntrlogout
   Log file.
   UNUSED
*/
#define C_filepntrlogout            C_fileTab[1].pntr
/* Define: C_filepntrsmbout
   Symbol matrix file.
   UNUSED
*/
#define C_filepntrsmbout            C_fileTab[2].pntr
/* Define: C_filepntrdrwout
   Drawing file.
   UNUSED
*/
#define C_filepntrdrwout            C_fileTab[3].pntr
/* Define: C_filepntrpltout
   Plot file.
   UNUSED
*/
#define C_filepntrpltout            C_fileTab[4].pntr
/* Define: C_filepntrordout
   Ordering file.
   UNUSED
*/
#define C_filepntrordout            C_fileTab[5].pntr

/*
**  The static and global definitions.
*/
/*
** Section: Functions
*/
/*
  Function: orderSplit

  Subdivide column blocks in column clocks of size blocksize.

  WARNING: unused

  Parameters:
    ordeptr   - Ordering.
    blocksize - size of block wanted.
*/
void orderSplit (Order * const ordeptr,
                 PASTIX_INT           blocksize)
{
  PASTIX_INT *rangtmp = NULL;
  PASTIX_INT i,j,k,cblktmp=0;

  for (i=0;i<ordeptr->cblknbr;i++)
    {
      cblktmp+=(ordeptr->rangtab[i+1]-ordeptr->rangtab[i])/blocksize;
      if ((ordeptr->rangtab[i+1]-ordeptr->rangtab[i])%blocksize)
        cblktmp++;
    }

  MALLOC_INTERN(rangtmp, cblktmp+1, PASTIX_INT);

  rangtmp[0]=ordeptr->rangtab[0];
  j=1;
  for (i=0;i<ordeptr->cblknbr;i++)
    {
      for (k=0;k<(ordeptr->rangtab[i+1]-ordeptr->rangtab[i])/blocksize;k++)
        rangtmp[j++]=ordeptr->rangtab[i]+(k+1)*blocksize;
      if ((ordeptr->rangtab[i+1]-ordeptr->rangtab[i])%blocksize)
        rangtmp[j++]=ordeptr->rangtab[i+1];
    }

  if (ordeptr->rangtab != NULL)
    memFree_null (ordeptr->rangtab);

  ordeptr->cblknbr=cblktmp;
  ordeptr->rangtab=rangtmp;

  /*
  printf("cblknbr=%ld\n",(long)ordeptr->cblknbr);
  for (i=0;i<=ordeptr->cblknbr;i++)
    printf("%ld ",(long)ordeptr->rangtab[i]);
  printf("\n");
  */
}

#if (defined SCOTCH_SEQSCOTCH || defined SCOTCH_H)
/*
   Function: orderSplit2

   Subdivide block columns in blocks of minimum size *bsmin*.
   Size depends on *rho* parameter too.

   WARNING: unused

   Parameters:
     ordeptr - Ordering.
     grphptr - Graph corresponding to the matrix.
     rho     - Parameter to compute new bloc size.
     bsmin   - Minimal bloc size.

 */
void orderSplit2 (Order        * const ordeptr,
                  SCOTCH_Graph * const grphptr,
                  double               rho,
                  PASTIX_INT                  bsmin)
{
  PASTIX_INT *rangtmp   = NULL;
  PASTIX_INT *blocksize = NULL;
  PASTIX_INT i,j,k,kk,cblktmp=0;
  double n;

  PASTIX_INT                   baseval;
  PASTIX_INT                   vertnbr;
  PASTIX_INT *                 verttab;
  PASTIX_INT                   edgenbr;
  PASTIX_INT *                 edgetab;

  SCOTCH_graphData (grphptr,
                    (SCOTCH_Num *) &baseval,
                    (SCOTCH_Num *) &vertnbr,
                    (SCOTCH_Num **)&verttab,
                    NULL, NULL, NULL,
                    (SCOTCH_Num *) &edgenbr,
                    (SCOTCH_Num **)&edgetab,
                    NULL);

  MALLOC_INTERN(blocksize, ordeptr->cblknbr, PASTIX_INT);

  for (i=0;i<ordeptr->cblknbr;i++)
    {
      /* find good blocksize */
      blocksize[i]=0;
      for (j=ordeptr->rangtab[i];j<ordeptr->rangtab[i+1];j++)
        {
          PASTIX_INT jj,max=0;
          jj=ordeptr->peritab[j];
          for (k=verttab[jj];k<verttab[jj+1];k++)
            {
              kk=ordeptr->permtab[edgetab[k]];
              if ((max<(j-kk)) && (kk>=ordeptr->rangtab[i]) && (kk<j))
                max=j-kk;
            }
          blocksize[i]+=max;
        }
      blocksize[i]++;
      n=(double)(ordeptr->rangtab[i+1]-ordeptr->rangtab[i]);
      blocksize[i]=(PASTIX_INT)(2.0*rho*((double)(blocksize[i]))/(1.0+sqrt((2.0*n-1.0)*(2.0*n-1.0)-8.0*((double)(blocksize[i]))))-1.0);
      if (blocksize[i]<bsmin)
        blocksize[i]=bsmin;
      if (blocksize[i]>ordeptr->rangtab[i+1]-ordeptr->rangtab[i])
        blocksize[i]=ordeptr->rangtab[i+1]-ordeptr->rangtab[i];

      printf("blocksize %ld [%ld]\n",
             (long)blocksize[i],
             (long)ordeptr->rangtab[i+1]-ordeptr->rangtab[i]);

      cblktmp+=(ordeptr->rangtab[i+1]-ordeptr->rangtab[i])/blocksize[i];
      if ((ordeptr->rangtab[i+1]-ordeptr->rangtab[i])%blocksize[i])
        cblktmp++;
    }

  MALLOC_INTERN(rangtmp, cblktmp+1, PASTIX_INT);

  rangtmp[0]=ordeptr->rangtab[0];
  j=1;
  for (i=0;i<ordeptr->cblknbr;i++)
    {
      for (k=0;k<(ordeptr->rangtab[i+1]-ordeptr->rangtab[i])/blocksize[i];k++)
        rangtmp[j++]=ordeptr->rangtab[i]+(k+1)*blocksize[i];
      if ((ordeptr->rangtab[i+1]-ordeptr->rangtab[i])%blocksize[i])
        rangtmp[j++]=ordeptr->rangtab[i+1];
    }

  if (ordeptr->rangtab != NULL)
    memFree_null (ordeptr->rangtab);

  ordeptr->cblknbr=cblktmp;
  ordeptr->rangtab=rangtmp;

  memFree_null(blocksize);

  /*
  printf("cblknbr=%ld\n",(long)ordeptr->cblknbr);
  for (i=0;i<=ordeptr->cblknbr;i++)
    printf("%ld ",(long)ordeptr->rangtab[i]);
  printf("\n");
  */
}

/*
  Function: orderSplit3

  Splits column blocs such has there is one column block in front of each bloc.

  Parameters:
    ordeptr  - Ordering.
    grphptr  - Graph associated to the matrix.
    matrsymb - Symbol matrix.

 */
void orderSplit3 (Order        * const ordeptr,
                  SCOTCH_Graph * const grphptr,
                  SymbolMatrix * const matrsymb)
{
  PASTIX_INT i,j,cblknum,bloknum,nodenbr;
  PASTIX_INT *mark1   = NULL;
  PASTIX_INT *rangtmp = NULL;
  PASTIX_INT cblktmp=0;

  PASTIX_INT                   baseval;
  PASTIX_INT                   vertnbr;
  PASTIX_INT *                 verttab;
  PASTIX_INT                   edgenbr;
  PASTIX_INT *                 edgetab;

  /*
  for (i=0;i<=ordeptr->cblknbr;i++)
    printf("%ld ",(long)ordeptr->rangtab[i]);
  printf("\n");
  */

  SCOTCH_graphData (grphptr,
                    (SCOTCH_Num *) &baseval,
                    (SCOTCH_Num *) &vertnbr,
                    (SCOTCH_Num **)&verttab,
                    NULL, NULL, NULL,
                    (SCOTCH_Num *) &edgenbr,
                    (SCOTCH_Num **)&edgetab,
                    NULL);

  nodenbr=ordeptr->rangtab[ordeptr->cblknbr];
  MALLOC_INTERN(mark1, nodenbr+1, PASTIX_INT);
  for (i=0;i<=nodenbr;i++)
    mark1[i]=0;

  for (cblknum=0;cblknum<matrsymb->cblknbr;cblknum++)
    for (bloknum=matrsymb->cblktab[cblknum].bloknum;bloknum<matrsymb->cblktab[cblknum+1].bloknum;bloknum++)
      {
        /*printf("%ld %ld %ld %ld\n",(long)cblknum,(long)bloknum,(long)matrsymb->bloktab[bloknum].frownum,(long)matrsymb->bloktab[bloknum].lrownum);*/
        /*if (matrsymb->bloktab[bloknum].levfval==0)*/
          {
            mark1[matrsymb->bloktab[bloknum].frownum]=1;
            mark1[matrsymb->bloktab[bloknum].lrownum+1]=1;
          }
      }

  /*
  for (i=0;i<nodenbr;i++)
    printf("%ld",(long)mark1[i]);
  printf("\n");
  */

  for (i=0;i<nodenbr;i++)
    if (mark1[i]==1) cblktmp++;
  MALLOC_INTERN(rangtmp, cblktmp+1, PASTIX_INT);
  j=0;
  for (i=0;i<nodenbr;i++)
    if (mark1[i]==1) rangtmp[j++]=i;
  rangtmp[j]=nodenbr;
  ASSERT(j==cblktmp,MOD_SOPALIN);

  if (ordeptr->rangtab != NULL)
    memFree_null (ordeptr->rangtab);

  ordeptr->cblknbr=cblktmp;
  ordeptr->rangtab=rangtmp;

  memFree_null(mark1);

  /*
  for (i=0;i<=ordeptr->cblknbr;i++)
    printf("%ld ",(long)ordeptr->rangtab[i]);
  printf("\n");
  */
}
#endif /* SCOTCH */

/*
  Function: symbolSplit

  Splits the symbol matrix.

  Parameters:
    matrsymb - symbol matrix
 */
void symbolSplit (SymbolMatrix * matrsymb)
{
  PASTIX_INT i,j,iter,add,cblknum,bloknum,num,frownum,lrownum;
  SymbolBlok *bloktmp   = NULL;
  SymbolCblk *cblktmp   = NULL;
  PASTIX_INT        *node2cblk = NULL;

  MALLOC_INTERN(node2cblk, matrsymb->nodenbr, PASTIX_INT);
  for (i=0;i<matrsymb->cblknbr+1;i++)
    for (j=matrsymb->cblktab[i].fcolnum;j<=matrsymb->cblktab[i].lcolnum;j++)
      node2cblk[j]=i;

  add=0;
  for (cblknum=0;cblknum<matrsymb->cblknbr-1;cblknum++)
    {
      for (bloknum=matrsymb->cblktab[cblknum].bloknum+1;bloknum<matrsymb->cblktab[cblknum+1].bloknum;bloknum++)
        {
          frownum=matrsymb->bloktab[bloknum].frownum;
          lrownum=matrsymb->bloktab[bloknum].lrownum;
          while (node2cblk[frownum]!=node2cblk[lrownum])
            {
              add++;
              num=frownum;
              while (node2cblk[num]==node2cblk[frownum]) num++;
              frownum=num;
            }
        }
    }

  MALLOC_INTERN(bloktmp, matrsymb->bloknbr+add, SymbolBlok);
  MALLOC_INTERN(cblktmp, matrsymb->cblknbr+1,   SymbolCblk);
  for (i=0;i<matrsymb->cblknbr+1;i++)
    {
      cblktmp[i].fcolnum=matrsymb->cblktab[i].fcolnum;
      cblktmp[i].lcolnum=matrsymb->cblktab[i].lcolnum;
      cblktmp[i].bloknum=matrsymb->cblktab[i].bloknum;
    }

  iter=0,add=0;
  for (cblknum=0;cblknum<matrsymb->cblknbr-1;cblknum++)
    {
      /* recopier le bloc diagonal */
      bloknum=matrsymb->cblktab[cblknum].bloknum;
      bloktmp[iter].frownum=matrsymb->bloktab[bloknum].frownum;
      bloktmp[iter].lrownum=matrsymb->bloktab[bloknum].lrownum;
      bloktmp[iter].cblknum=matrsymb->bloktab[bloknum].cblknum;
      bloktmp[iter].levfval=0;
      iter++;

      /* on recopie tous les blocs extra du bloc-colonne */
      for (bloknum=matrsymb->cblktab[cblknum].bloknum+1;bloknum<matrsymb->cblktab[cblknum+1].bloknum;bloknum++)
        {
          frownum=matrsymb->bloktab[bloknum].frownum;
          lrownum=matrsymb->bloktab[bloknum].lrownum;
          while (node2cblk[frownum]!=node2cblk[lrownum])
            {
              add++;
              num=frownum;
              while (node2cblk[num]==node2cblk[frownum]) num++;
              bloktmp[iter].frownum=frownum;
              bloktmp[iter].lrownum=num-1;
              bloktmp[iter].cblknum=node2cblk[frownum];
              bloktmp[iter].levfval=0;
              iter++;
              frownum=num;
            }
          bloktmp[iter].frownum=frownum;
          bloktmp[iter].lrownum=lrownum;
          bloktmp[iter].cblknum=node2cblk[frownum];
          bloktmp[iter].levfval=0;
          iter++;
        }

      cblktmp[cblknum+1].bloknum+=add;
    }

  bloktmp[iter].frownum=matrsymb->cblktab[matrsymb->cblknbr-1].fcolnum;
  bloktmp[iter].lrownum=matrsymb->cblktab[matrsymb->cblknbr-1].lcolnum;
  bloktmp[iter].cblknum=matrsymb->cblknbr-1;
  bloktmp[iter].levfval=0;
  cblktmp[matrsymb->cblknbr].bloknum+=add;

  memFree_null(matrsymb->bloktab);
  memFree_null(matrsymb->cblktab);
  memFree_null(node2cblk);
  matrsymb->bloktab=bloktmp;
  matrsymb->cblktab=cblktmp;
  matrsymb->bloknbr+=add;
}

/*
  Function: symbolRustine

  DESCRIPTION TO FILL

  Parameters:
    matrsymb  - Symbol matrix
    matrsymb2 - Symbol matrix
 */
void
symbolRustine (SymbolMatrix *       matrsymb,
               SymbolMatrix * const matrsymb2)
{
  PASTIX_INT i,iter,add,cblknum,bloknum,bloknum2;
  SymbolBlok *bloktmp = NULL;
  SymbolCblk *cblktmp = NULL;

  MALLOC_INTERN(bloktmp, matrsymb->bloknbr+matrsymb->cblknbr, SymbolBlok);
  MALLOC_INTERN(cblktmp, matrsymb->cblknbr+1,                 SymbolCblk);
  for (i=0;i<matrsymb->cblknbr+1;i++)
    {
      cblktmp[i].fcolnum=matrsymb->cblktab[i].fcolnum;
      cblktmp[i].lcolnum=matrsymb->cblktab[i].lcolnum;
      cblktmp[i].bloknum=matrsymb->cblktab[i].bloknum;
    }

  iter=0,add=0;
  for (cblknum=0;cblknum<matrsymb->cblknbr-1;cblknum++)
    {
      /* recopier le bloc diagonal */
      bloknum=matrsymb->cblktab[cblknum].bloknum;
      bloktmp[iter].frownum=matrsymb->bloktab[bloknum].frownum;
      bloktmp[iter].lrownum=matrsymb->bloktab[bloknum].lrownum;
      bloktmp[iter].cblknum=matrsymb->bloktab[bloknum].cblknum;
      bloktmp[iter].levfval=matrsymb->bloktab[bloknum].levfval;
      iter++;

      bloknum=matrsymb->cblktab[cblknum].bloknum+1;
      bloknum2=matrsymb2->cblktab[cblknum].bloknum+1;

      if (bloknum==matrsymb->cblktab[cblknum+1].bloknum)
        {
          /* pas d'extra diag */
          if (matrsymb==matrsymb2)
            {
              add++;
#ifdef RUSTIN_ADD_NEXT_CBLK
              bloktmp[iter].frownum=matrsymb->cblktab[cblknum+1].fcolnum;
              bloktmp[iter].lrownum=matrsymb->cblktab[cblknum+1].fcolnum;
              bloktmp[iter].cblknum=cblknum+1;
#else
              bloktmp[iter].frownum=matrsymb->cblktab[matrsymb->cblknbr-1].fcolnum;
              bloktmp[iter].lrownum=matrsymb->cblktab[matrsymb->cblknbr-1].fcolnum;
              bloktmp[iter].cblknum=matrsymb->cblknbr-1;
#endif
              bloktmp[iter].levfval=matrsymb->bloktab[bloknum].levfval;
              printf("add blok dans cblk %ld : %ld %ld %ld\n",(long)cblknum,
                     (long)bloktmp[iter].frownum,(long)bloktmp[iter].lrownum,
                     (long)bloktmp[iter].cblknum);
              iter++;
            }
          else
            {
              ASSERT(bloknum2!=matrsymb2->cblktab[cblknum+1].bloknum,
                     MOD_SOPALIN);
              add++;
              bloktmp[iter].frownum=matrsymb2->bloktab[bloknum2].frownum;
              bloktmp[iter].lrownum=matrsymb2->bloktab[bloknum2].frownum;
              bloktmp[iter].cblknum=matrsymb2->bloktab[bloknum2].cblknum;
              bloktmp[iter].levfval=matrsymb2->bloktab[bloknum2].levfval;
              iter++;
            }
        }
      else
        {
          if (matrsymb->bloktab[bloknum].cblknum!=
              matrsymb2->bloktab[bloknum2].cblknum)
            {
              /* le premier extra diag ne va pas */
              add++;
              bloktmp[iter].frownum=matrsymb2->bloktab[bloknum2].frownum;
              bloktmp[iter].lrownum=matrsymb2->bloktab[bloknum2].frownum;
              bloktmp[iter].cblknum=matrsymb2->bloktab[bloknum2].cblknum;
              bloktmp[iter].levfval=matrsymb2->bloktab[bloknum2].levfval;
              iter++;
            }

          /* on recopie tous les blocs extra du bloc-colonne */
          for (bloknum=matrsymb->cblktab[cblknum].bloknum+1;bloknum<matrsymb->cblktab[cblknum+1].bloknum;bloknum++)
            {
              bloktmp[iter].frownum=matrsymb->bloktab[bloknum].frownum;
              bloktmp[iter].lrownum=matrsymb->bloktab[bloknum].lrownum;
              bloktmp[iter].cblknum=matrsymb->bloktab[bloknum].cblknum;
              bloktmp[iter].levfval=matrsymb->bloktab[bloknum].levfval;
              iter++;
            }
        }

      cblktmp[cblknum+1].bloknum+=add;
    }

  bloktmp[iter].frownum=matrsymb->cblktab[matrsymb->cblknbr-1].fcolnum;
  bloktmp[iter].lrownum=matrsymb->cblktab[matrsymb->cblknbr-1].lcolnum;
  bloktmp[iter].cblknum=matrsymb->cblknbr-1;
  bloktmp[iter].levfval=0;
  cblktmp[matrsymb->cblknbr].bloknum+=add;

  memFree_null(matrsymb->bloktab);
  memFree_null(matrsymb->cblktab);
  matrsymb->bloktab=bloktmp;
  matrsymb->cblktab=cblktmp;
  matrsymb->bloknbr+=add;
  ASSERT(add<matrsymb->cblknbr,MOD_SOPALIN);
}

/*
  struct: KeepParam_

  TODO: Fill-in
*/
typedef struct KeepParam_ {
  SymbolKeep * restrict     keepptr;              /*+ Reference to keepdat for access   +*/
  PASTIX_INT                       levfmax;              /*+ Inclusive maximum level of fill   +*/
  PASTIX_INT                       ctrimin;              /*+ Minimum output contribution value +*/
} KeepParam;

/* ********************************************
   Functions: These are the selection functions

   This routine should return 1 if the
   given block is to be kept and 0 else.

   CAUTION: blocks with (levfval<1) should
   always be kept.
* *********************************************/


#if (defined SCOTCH_SEQSCOTCH || defined SCOTCH_H)
/*
   Function: keepFuncD1

   Test if kblkptr->levfval > kparptr->levfmax

   Parameters:
     kblkptr - Block to test
     dataptr - Data fields

   Returns:
     API_YES - if kblkptr->levfval > kparptr->levfmax
     API_NO  - otherwise
*/
static int keepFuncD1 (const SymbolKeepBlok * const  kblkptr,
                       void * const                  dataptr)
{
  KeepParam * restrict      kparptr;

  kparptr = (KeepParam *) dataptr;

  if (kblkptr->levfval > kparptr->levfmax)
    return (API_YES);

  return (API_NO);
}
#endif /* SCOTCH */

#ifdef DEADCODE
/*
   Function: keepFuncH1

   Test if kblkptr->levfval == kparptr->levfmax

   Parameters:
     kblkptr - Block to test
     dataptr - Data fields

   Returns:
     API_YES - if kblkptr->levfval == kparptr->levfmax
     API_NO  - otherwise
*/
static int keepFuncH1 (const SymbolKeepBlok * const  kblkptr,
                       void * const                  dataptr)
{
  KeepParam * restrict      kparptr;

  kparptr = (KeepParam *) dataptr;

  if (kblkptr->levfval == kparptr->levfmax)
    return (1);

  return (0);
}

/*
   Function: keepFuncD2

   Test if kblkptr->levfval == kparptr->levfmax and kblkptr->ctrival <= kparptr->ctrimin

   Parameters:
     kblkptr - Block to test
     dataptr - Data fields

   Returns:
     API_YES - if kblkptr->levfval == kparptr->levfmax and kblkptr->ctrival <= kparptr->ctrimin
     API_NO  - otherwise
*/
static int keepFuncD2 (const SymbolKeepBlok * const  kblkptr,            /*+ Block to test +*/
                       void * const                  dataptr)            /*+ Data fields   +*/
{
  KeepParam * restrict      kparptr;

  kparptr = (KeepParam *) dataptr;

  if ((kblkptr->levfval == kparptr->levfmax) && (kblkptr->ctrival <= kparptr->ctrimin))
    return (1);

  return (0);
}
/*
   Function: keepFuncH2

   Test if kparptr->keepptr->keeptab[kblkptr - kparptr->keepptr->kblktab] != 0

   Parameters:
     kblkptr - Block to test
     dataptr - Data fields

   Returns:
     API_YES - if kparptr->keepptr->keeptab[kblkptr - kparptr->keepptr->kblktab] != 0
     API_NO  - otherwise
*/
static int keepFuncH2 (const SymbolKeepBlok * const  kblkptr,            /*+ Block to test +*/
                       void * const                  dataptr)            /*+ Data fields   +*/
{
  KeepParam * restrict      kparptr;

  kparptr = (KeepParam *) dataptr;

  if (kparptr->keepptr->keeptab[kblkptr - kparptr->keepptr->kblktab] != 0) /* Consider only selected blocks */
    return (1);

  return (0);
}
#endif /* DEADCODE */

/***********************/
/*                     */
/* These are the color */
/* handling routines.  */
/*                     */
/***********************/
#if 0
static float ratcol[6][3] = { { 0.0, 0.0, 0.0 },
                              { 0.0, 0.0, 1.0 },
                              { 0.0, 1.0, 0.0 },
                              { 1.0, 0.8, 0.0 },
                              { 1.0, 0.3, 0.3 },
                              { 1.0, 0.3, 0.3 } };

#endif

#ifdef DEADCODE
static int blokColor (const SymbolMatrix * const  symbptr,              /*+ Symbol matrix +*/
                      const SymbolBlok * const    blokptr,              /*+ Block to test +*/
                      void * const                dataptr,              /*+ Data fields   +*/
                      float * const               colotab)
{
  const SymbolKeep * restrict     keepptr;
  const SymbolKeepBlok * restrict kblkptr;
  PASTIX_INT                             bloknum;

  keepptr = (SymbolKeep *) dataptr;
  bloknum = blokptr - symbptr->bloktab;
  kblkptr = keepptr->kblktab + bloknum;

#if 0
  rat = (float) kblkptr->levfval / (float) sydrptr->levfmax;
  idx = (int) (rat * 4.0);
  rat = (rat - ((float) idx * 0.25)) * 4.0;

  colotab[0] = ratcol[idx][0] * rat + ratcol[idx + 1][0] * (1.0 - rat);
  colotab[1] = ratcol[idx][1] * rat + ratcol[idx + 1][1] * (1.0 - rat);
  colotab[2] = ratcol[idx][2] * rat + ratcol[idx + 1][2] * (1.0 - rat);
#endif
#if 1
  if (kblkptr->levfval < 2) {
    colotab[0] = 1.0;
    colotab[1] =
    colotab[2] = 0.2;
  }
  else if (kblkptr->ctrival > ((3 * keepptr->ctrimax) / 5)) {
    colotab[0] = 0.0;
    colotab[1] = 1.0;
    colotab[2] = 0.2;
  }
  else {
    colotab[0] =
    colotab[1] =
    colotab[2] = 0.4;
  }
#endif
#if 0
  if (keepptr->keeptab[bloknum] != 0) {
    colotab[0] =
    colotab[1] =
    colotab[2] = 0.4;
  }
  else {
    colotab[0] = 1.0;
    colotab[1] =
    colotab[2] = 0.2;
  }
#endif
  return (1);
}
#endif /* DEADCODE */

/* *********************************************
   Functions: These are the cost functions.

   This routine computes the factorization
   and solving cost of the given symbolic
   block matrix, whose nodes hold the number
   of DOFs given by the proper DOF structure.
   To ensure maximum accuracy and minimum loss
   of precision, costs are summed-up recursively.
   It returns:
   - 0   : on success.
   - !0  : on error.
 * *********************************************/
/*
   Function: symbolCostn2

   Function declaration.
   TO FILL

   Parameters:
     cblktax - Based access to cblktab
     bloktax - Based access to bloktab
     deofptr - DOF structure associated with the matrix
     keeptab - Flag array for blocks to keep
     nnzptr  - Size of the structure, to be filled
     opcptr  - Operation count, to be filled
     cblkmin - Minimum column block index to consider
     cblknbr - Number of column blocks to consider
 */
static void symbolCostn2(const SymbolCblk * const cblktax,
                         const SymbolBlok * const bloktax,
                         const Dof        * const deofptr,
                         const unsigned char       * const keeptab,
                         double           * const nnzptr,
                         double           * const opcptr,
                         const PASTIX_INT                cblkmin,
                         const PASTIX_INT                cblknbr);

/*
   Function: symbolCostn

   TO FILL

   Parameters:
     symbptr - Symbolic matrix to evaluate
     deofptr - DOF structure associated with the matrix
     keeptab - Flag array for blocks to keep
     typeval - Type of cost computation
     nnzptr  - Size of the structure, to be filled
     opcptr  - Operation count, to be filled


 */
int symbolCostn (const SymbolMatrix * const  symbptr,
                 const Dof * const           deofptr,
                 const unsigned char * const          keeptab,
                 const SymbolCostType        typeval,
                 double * const              nnzptr,
                 double * const              opcptr)
{
  if (typeval != SYMBOLCOSTLDLT) {
    errorPrint ("symbolCostn: cost function not supported");
    return     (1);
  }

  *opcptr = 0.0L;
  *nnzptr = 0.0L;

  symbolCostn2 (symbptr->cblktab - symbptr->baseval, /* Perform recursion on column blocks */
                symbptr->bloktab - symbptr->baseval,
    deofptr, keeptab - symbptr->baseval, nnzptr, opcptr, symbptr->baseval, symbptr->cblknbr);

  return (0);
}

/*
   Function: symbolCostn2

   TO FILL

   Parameters:
     cblktax - Based access to cblktab
     bloktax - Based access to bloktab
     deofptr - DOF structure associated with the matrix
     keeptax - Flag array for blocks to keep
     nnzptr  - Size of the structure, to be filled
     opcptr  - Operation count, to be filled
     cblkmin - Minimum column block index to consider
     cblknbr - Number of column blocks to consider
*/
static void symbolCostn2 (
                          const SymbolCblk * restrict const cblktax,
                          const SymbolBlok * restrict const bloktax,
                          const Dof * restrict const        deofptr,
                          const unsigned char * const                keeptax,
                          double * restrict const           nnzptr,
                          double * restrict const           opcptr,
                          const PASTIX_INT                         cblkmin,
                          const PASTIX_INT                         cblknbr)
{
  PASTIX_INT                 bloknum;                    /* Number of current extra-diagonal block             */
  PASTIX_INT                 cmednum;                    /* Median column block number                         */
  PASTIX_INT                 cfacnum;                    /* Number of facing column block                      */
  PASTIX_INT                 cdofnbr;                    /* Number of DOFs in column block (l_k)               */
  PASTIX_INT                 rdofsum;                    /* Number of DOFs in all row blocks (g_{ki} or g_{k}) */
  double              nnzval;                     /* Number of non-zeroes in subtree                    */
  double              opcval;                     /* Operation count in subtree                         */

  nnzval =                                        /* Initialize local values */
  opcval = 0.0L;

  if (cblknbr > 1) {                              /* If more than one column block, perform recursion */
    cmednum = cblknbr / 2;
    symbolCostn2 (cblktax, bloktax, deofptr, keeptax, &nnzval, &opcval, cblkmin, cmednum);
    symbolCostn2 (cblktax, bloktax, deofptr, keeptax, &nnzval, &opcval, cblkmin + cmednum, cblknbr - cmednum);

    *nnzptr += nnzval;                            /* Sum-up local values */
    *opcptr += opcval;
  }
  else {                                          /* Single column block                              */
    PASTIX_INT                 rdounbr;                  /* Number of DOFs in undropped row blocks (h'_{ki}) */
    PASTIX_INT                 rdousum;                  /* Number of DOFs in undropped row blocks (h'_{ki}) */

    cdofnbr = noddVal (deofptr, cblktax[cblkmin].lcolnum + 1) -
              noddVal (deofptr, cblktax[cblkmin].fcolnum);

    bloknum = cblktax[cblkmin].bloknum + 1;       /* Get index of first extra-diagonal block */

    rdofsum =
    rdousum =
    rdounbr = 0;

    for (bloknum = cblktax[cblkmin + 1].bloknum - 1; /* Scan extra-diagonals, backwards */
         bloknum > cblktax[cblkmin].bloknum; ) {
      if (keeptax[bloknum] == 0) {                /* Skip dropped blocks */
        bloknum --;
        continue;
      }

      rdousum += rdounbr;
      rdounbr  = 0;

      cfacnum = bloktax[bloknum].cblknum;
      do {
        PASTIX_INT                 rdofblk;              /* Number of DOFs in local block */

        if (keeptax[bloknum] == 0)                /* Skip dropped blocks */
          continue;

        rdofblk  = noddVal (deofptr, bloktax[bloknum].lrownum + 1) -
                   noddVal (deofptr, bloktax[bloknum].frownum);
        rdofsum += rdofblk;                       /* Account for undropped blocks */

        rdounbr += rdofblk;                       /* Add contribution to undropped set of blocks */
      } while (bloktax[-- bloknum].cblknum == cfacnum);

      opcval += ((double) (rdounbr)) *            /* Count C3'(k,i) + C3''(k,i) */
                ((double) (rdounbr + rdousum)) *
                ((double) (2 * cdofnbr + 1));
    }

    *nnzptr += ((double) (cdofnbr + rdofsum)) * ((double) cdofnbr); /* Sum-up stored coefficients */
    *opcptr += opcval +
      ((double) cdofnbr) *               /* Count C1(k) + C2(k) */
      (double)(((double) cdofnbr) * ((double) (2 * cdofnbr + 6 * rdofsum + 3)) + 1.0L) /(double) 6.0L;
  }
}


#if (defined SCOTCH_SEQSCOTCH || defined SCOTCH_H)
/* **************************************
   Functions: This is the main function.
 ****************************************/
/*
  Function: bordi

  This is the main function
  TO FILL

  Parameters:
    alpha    -
    symbptr  -
    graphptr -
    orderptr -

 */
void bordi(int            alpha,
           SymbolMatrix * symbptr,
           SCOTCH_Graph * graphptr,
           Order        * orderptr)
{
  SymbolKeep          keepdat;
  PASTIX_INT                 baseval;                    /* Graph base value              */
  PASTIX_INT                 vertnbr;                    /* Number of vertices in graph   */
  PASTIX_INT                 edgenbr;                    /* Number of edges in graph      */
  PASTIX_INT *               verttab;
  PASTIX_INT *               edgetab;
  double              grafnnz;                    /* Number of symmetric non-zeros */
  Dof                 matrdeof;                   /* Matrix DOF structure          */
  double              matrnnz;                    /* Non-zeroes                    */
  double              matrnnzmax;                 /* Maximum non-zeroes            */
  double              matropc;                    /* Operation count               */
  double              matropcmax;                 /* Maximum operation count       */
  PASTIX_INT                 leafnbr;                    /* Number of leaves              */
  PASTIX_INT                 heigmin;                    /* Minimum height                */
  PASTIX_INT                 heigmax;                    /* Maximum height                */
  double              heigavg;                    /* Average height                */
  double              heigdlt;                    /* Deviation                     */
/*   PASTIX_INT                 levfnum;                    /\* Current level of fill         *\/ */
  PASTIX_INT                 bloksum;                    /* Accumulated number of blocks  */
  PASTIX_INT                 bloknum;
  Clock               runtime[3];                 /* Timing variables                       */
  double              fillrat;                    /* Fill ratio : 0.0 -> NNZA ; 1.0 -> NNZL */
/*   int                 i; */
  char                buftab1[1024];
  char                buftab2[1024];
/*   FILE *              stream; */

#ifdef ALPHA_LEVFK
  fillrat = 1.0;
#else
  fillrat = alpha/100;                                 /* Default fill ratio */
#endif

  buftab1[0] = 'Z';
  buftab1[1] = '\0';

  SCOTCH_graphData (graphptr, (SCOTCH_Num *)&baseval, (SCOTCH_Num *)&vertnbr, (SCOTCH_Num **)&verttab,
                    NULL, NULL, NULL, (SCOTCH_Num *)&edgenbr, (SCOTCH_Num **)&edgetab, NULL);

  {
    SymbolMatrix    symbptr2;

    clockInit  (&runtime[0]);
    clockStart (&runtime[0]);

    symbolInit      (&symbptr2);
    symbolFaxGraph  (&symbptr2, graphptr, orderptr);
    printf("Avant Split symb : cblkbr=%ld bloknbr=%ld\n",(long)symbptr2.cblknbr,(long)symbptr2.bloknbr);
    /*
      symbolFaxiGraph (&symbptr2, graphptr, orderptr);
      symbolFaxiGraphZero (&symbptr2, graphtr, orderptr);
    */
    symbolCheck     (&symbptr2);

    /*orderSplit      (orderptr, 1);*/
    /*orderSplit2     (orderptr, graphptr, 0.1, 4);*/
    orderSplit3     (orderptr, graphptr, &symbptr2);

    orderCheck      (orderptr);
    symbolExit      (&symbptr2);

    clockStop  (&runtime[0]);
  }

  dofInit     (&matrdeof);
  dofConstant (&matrdeof, baseval, vertnbr, 1);   /* One DOF per node */

  symbolInit      (symbptr);

  clockInit  (&runtime[1]);
  clockStart (&runtime[1]);

#ifdef ALPHA_LEVFK
  printf("symbolFaxiGraph with k=%ld\n",(long)alpha);
  /*symbolFaxiGraph (symbptr, graphptr, orderptr, alpha);*/
  ifax(orderptr->rangtab[orderptr->cblknbr],verttab,edgetab,alpha,orderptr->cblknbr,orderptr->rangtab,orderptr->permtab,orderptr->peritab,symbptr);
#else
  symbolFaxiGraph (symbptr, graphptr, orderptr);
#endif

  clockStop  (&runtime[1]);

  symbolCheck(symbptr);

  {
    SymbolMatrix    symbptr3;

    clockInit  (&runtime[2]);
    clockStart (&runtime[2]);

    symbolInit      (&symbptr3);
    symbolFaxGraph  (&symbptr3, graphptr, orderptr);
    symbolCheck     (&symbptr3);

    printf("Apres Split symb : cblkbr=%ld bloknbr=%ld\n",(long)symbptr3.cblknbr,(long)symbptr3.bloknbr);
    printf("Apres Split isymb : cblkbr=%ld bloknbr=%ld\n",(long)symbptr->cblknbr,(long)symbptr->bloknbr);
    printf("debut Rustine\n");
    symbolRustine   (symbptr, &symbptr3);
    printf("fin Rustine\n");
    printf("Apres Split isymb.rustine : cblkbr=%ld bloknbr=%ld\n",(long)symbptr->cblknbr,(long)symbptr->bloknbr);
    /*
    symbolSplit     (symbptr);
    printf("Apres Split isymb.split : cblkbr=%ld bloknbr=%ld\n",(long)symbptr->cblknbr,(long)symbptr->bloknbr);
    symbolCompact   (symbptr);
    printf("Apres Split isymb.compact : cblkbr=%ld bloknbr=%ld\n",(long)symbptr->cblknbr,(long)symbptr->bloknbr);
    */

    symbolExit      (&symbptr3);

    clockStop  (&runtime[2]);
  }

  symbolCheck (symbptr);

  symbolCost  (symbptr  , &matrdeof, SYMBOLCOSTLDLT, &matrnnzmax, &matropcmax);
  symbolTree  (symbptr, &matrdeof, &leafnbr,
               &heigmin, &heigmax, &heigavg, &heigdlt);

  grafnnz = (double) edgenbr / 2;                 /* Get number of symmetric non-zeros */

  printf ("A\tGrafnnz=%e\n",
          (double) grafnnz);
  printf ("B\tCblknbr=%6ld\nB\tBloknbr=%6ld\t(%5.3f)\nB\tNNZ=%e\t(%5.3f)\tStor=%e\nB\tOPC=%e\n",
          (long) symbptr->cblknbr, (long) symbptr->bloknbr,
          (double) ((double) symbptr->bloknbr / (double) symbptr->cblknbr),
          matrnnzmax,
          (double) ((double) matrnnzmax / (double) grafnnz),
          (double) (3 * symbptr->cblknbr + 2 * symbptr->bloknbr),
          matropcmax);
  printf ("B\tLeafnbr=%6ld\nB\tTree hmin=%ld\thmax=%ld\thavg=%6.2f\thdlt=%4.2f\n",
          (long) leafnbr, (long) heigmin, (long) heigmax, heigavg, (double) (heigdlt * 100.0));

  fprintf(stderr,"SF time (prepar) %lf\n",clockVal(&runtime[0]));
  fprintf(stderr,"SF time (QGifax) %lf\n",clockVal(&runtime[1]));
  fprintf(stderr,"SF time (bricol) %lf\n",clockVal(&runtime[2]));

#ifdef ALPHA_LEVFK
  return;
#endif

  /* Selection area. This is where one builds
  ** the list of actions to perform to select
  ** the blocks to keep and to remove.
  */

  {
    KeepParam           kpardat;
    PASTIX_INT                 levfnum;
    double              nnzsum;
    double              fillsiz;                  /* Resulting size to achieve */

    fillsiz = grafnnz + fillrat * (matrnnzmax - grafnnz); /* Compute requested fill size */

#ifndef PATCH_HENON_BLEND
    {
      PASTIX_INT                 cblknum;

      for (cblknum = 0; cblknum < symbptr->cblknbr; cblknum ++) {
        if (symbptr->cblktab[cblknum + 1].bloknum - symbptr->cblktab[cblknum].bloknum >= 2) /* If extra-diagonal present              */
          symbptr->bloktab[symbptr->cblktab[cblknum].bloknum + 1 - symbptr->baseval].levfval = 0; /* Always keep first extra-diagonal */
      }
    }
#endif /* PATCH_HENON_BLEND */

    kpardat.keepptr = &keepdat;                   /* Set reference to keep structure            */
    symbolKeepInit (&keepdat, symbptr);         /* All blocks selected at initialization time */

    {
      PASTIX_INT                 cblknum;
      long                bloknbr;
      long                blokkep;

      blokkep =
        bloknbr = symbptr->cblknbr;
      for (cblknum = 0; cblknum < symbptr->cblknbr; cblknum ++) {
        PASTIX_INT                 bloklst;

        bloklst = symbptr->cblktab[cblknum].bloknum - symbptr->baseval; /* Last is diagonal block */

        for (bloknum = symbptr->cblktab[cblknum].bloknum + 1;
             bloknum < symbptr->cblktab[cblknum + 1].bloknum; bloknum ++) {
          if (keepdat.keeptab[bloknum - symbptr->baseval] != 0) {
            bloknbr ++;
            if ((symbptr->bloktab[bloknum - symbptr->baseval].cblknum != symbptr->bloktab[bloklst].cblknum)     ||
                (symbptr->bloktab[bloknum - symbptr->baseval].frownum != symbptr->bloktab[bloklst].lrownum + 1) ||
                (keepdat.kblktab[bloknum - symbptr->baseval].levfval  != keepdat.kblktab[bloklst].levfval)      ||
                (keepdat.kblktab[bloknum - symbptr->baseval].ctrival  != keepdat.kblktab[bloklst].ctrival)) {
              blokkep ++;
              bloklst = bloknum - symbptr->baseval;
            }
          }
        }
      }
      printf ("K0\tIndepBlocks=%ld out of %ld\n",
               (long) blokkep, (long) bloknbr);
    }
    printf ("B\tLevfmax=%ld\tNupdmax=%ld\nB\tCtrimax=%ld\tCtromax=%ld\tHghtmax=%ld\n",
             (long) keepdat.levfmax, (long) keepdat.nupdmax, (long) keepdat.ctrimax, (long) keepdat.ctromax, (long) keepdat.hghtmax);

    symbolKeepHisto (&keepdat, symbptr, NULL, NULL); /* Compute histograms for whole matrix */
    sprintf (buftab2, "%s_init", buftab1);        /* Output initial histogram results         */
    /*symbolKeepView (&keepdat, grafnnz, buftab2);*/

    for (levfnum = 1, nnzsum = keepdat.levftab[0] + keepdat.levftab[1]; /* For all level-of-fill values, starting from 2 */
         levfnum < keepdat.levfmax; nnzsum += keepdat.levftab[++ levfnum]) {
      if (nnzsum >= fillsiz)                      /* Keep only blocks up to maximum fill size, and just above */
        break;
    }
#ifdef ALPHA_LEVFK
    kpardat.levfmax = alpha;
#else
    kpardat.levfmax = levfnum;                    /* Set levfmax for del and histo routines */
#endif

    printf ("S\tLevfnum=%ld\n",
            (long) kpardat.levfmax);

    symbolKeepDel (&keepdat, symbptr, keepFuncD1, &kpardat); /* Remove blocks such that (levfval > levfnum)   */
    symbolKeepCompute (&keepdat, symbptr);      /* Re-compute block parameters                                */

    {
      PASTIX_INT                 cblknum;
      long                bloknbr;
      long                blokkep;

      blokkep =
      bloknbr = symbptr->cblknbr;
      for (cblknum = 0; cblknum < symbptr->cblknbr; cblknum ++) {
        PASTIX_INT                 bloklst;

        bloklst = symbptr->cblktab[cblknum].bloknum - symbptr->baseval; /* Last is diagonal block */

        for (bloknum = symbptr->cblktab[cblknum].bloknum + 1;
             bloknum < symbptr->cblktab[cblknum + 1].bloknum; bloknum ++) {
          if (keepdat.keeptab[bloknum - symbptr->baseval] != 0) {
            bloknbr ++;
            if ((symbptr->bloktab[bloknum - symbptr->baseval].cblknum != symbptr->bloktab[bloklst].cblknum)     ||
                (symbptr->bloktab[bloknum - symbptr->baseval].frownum != symbptr->bloktab[bloklst].lrownum + 1) ||
                (keepdat.kblktab[bloknum - symbptr->baseval].levfval  != keepdat.kblktab[bloklst].levfval)      ||
                (keepdat.kblktab[bloknum - symbptr->baseval].ctrival  != keepdat.kblktab[bloklst].ctrival)) {
              blokkep ++;
              bloklst = bloknum - symbptr->baseval;
            }
          }
        }
      }
      printf ("K1\tIndepBlocks=%ld out of %ld\n",
              (long) blokkep, (long) bloknbr);

    }

#ifdef DEADCODE
    {
      PASTIX_INT                 ctrinum;

      symbolKeepHisto (&keepdat, symbptr, keepFuncH1, &kpardat); /* Compute histograms for (levfval == levfnum) */
      sprintf (buftab2, "%s_colv", buftab1);        /* Output histograms for these blocks                         */
      symbolKeepView (&keepdat, grafnnz, buftab2);

      for (ctrinum = keepdat.ctrimax, nnzsum -= keepdat.levftab[levfnum]; /* For all output contributions at (levfval == levfnum) */
           ctrinum >= 0; nnzsum += keepdat.ctritab[ctrinum], ctrinum --) {
        if ((nnzsum + keepdat.ctritab[ctrinum]) > fillsiz) /* Keep only blocks up to maximum fill size */
          break;
      }
      kpardat.ctrimin = ctrinum;                    /* Set minimum level of input contributions at (levfval == levfnum) to keep */

      printf ("S\tCtrinum=%ld\n", (long) kpardat.ctrimin);

      symbolKeepDel     (&keepdat, symbptr, keepFuncD2, &kpardat); /* Remove blocks according to funcD1 */
      symbolKeepCompute (&keepdat, symbptr);      /* Re-compute block parameters for remaining blocks   */

      {
        PASTIX_INT                 cblknum;
        PASTIX_INT                 bloknum;
        long                bloknbr;
        long                blokkep;

        blokkep =
          bloknbr = symbptr->cblknbr;
        for (cblknum = 0; cblknum < symbptr->cblknbr; cblknum ++) {
          PASTIX_INT                 bloklst;

          bloklst = symbptr->cblktab[cblknum].bloknum - symbptr->baseval; /* Last is diagonal block */

          for (bloknum = symbptr->cblktab[cblknum].bloknum + 1;
               bloknum < symbptr->cblktab[cblknum + 1].bloknum; bloknum ++) {
            if (keepdat.keeptab[bloknum - symbptr->baseval] != 0) {
              bloknbr ++;
              if ((symbptr->bloktab[bloknum - symbptr->baseval].cblknum != symbptr->bloktab[bloklst].cblknum)     ||
                  (symbptr->bloktab[bloknum - symbptr->baseval].frownum != symbptr->bloktab[bloklst].lrownum + 1) ||
                  (keepdat.kblktab[bloknum - symbptr->baseval].levfval  != keepdat.kblktab[bloklst].levfval)      ||
                  (keepdat.kblktab[bloknum - symbptr->baseval].ctrival  != keepdat.kblktab[bloklst].ctrival)) {
                blokkep ++;
                bloklst = bloknum - symbptr->baseval;
              }
            }
          }
        }
        printf ("K2\tIndepBlocks=%ld out of %ld\n",
                (long) blokkep, (long) bloknbr);
      }

      symbolKeepHisto (&keepdat, symbptr, keepFuncH2, &kpardat); /* Compute histograms for all selected blocks */
      sprintf (buftab2, "%s_finl", buftab1);        /* Output final results                                      */
      symbolKeepView (&keepdat, grafnnz, buftab2);
    }
#endif /* DEADCODE */

  }

/* End of selection area.
*/

  for (bloknum = 0, bloksum = 0; bloknum < symbptr->bloknbr; bloknum ++) {
    if (keepdat.keeptab[bloknum] != 0)
      bloksum ++;
  }
  symbolCostn (symbptr, &matrdeof, keepdat.keeptab, SYMBOLCOSTLDLT, &matrnnz, &matropc);
  printf ("F\tBloknbr=%.4ld\tNNZ=%e\t(%5.3lf)\n",
          (long) bloksum,
          (double) matrnnz,
          (double) (matrnnz / grafnnz));
  printf ("F\tOPC=%e\t(%e)\n",
          (double) matropc,
          (matropc * 100.0) / (double) matropcmax);
  printf ("F\tMemTot=%e\tNbiter=%5.3lf\n",
          (double) (2.0 * grafnnz + matrnnz + (double) (bloksum * 2) + (double) (symbptr->cblknbr * 3)),
          (((double) (matropcmax - matropc - grafnnz) / (double) (grafnnz + (4 * matrnnz - 3 * symbptr->nodenbr))) - 1.0) / 2.0);

  symbolKeepPurge (&keepdat, symbptr);        /* CAUTION: These two change the structure of the matrix ! */
  symbolCompact   (symbptr);

  dofExit(&matrdeof);
  symbolKeepExit(&keepdat);
}
#endif /* SCOTCH */
