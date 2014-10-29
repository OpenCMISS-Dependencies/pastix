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
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common_pastix.h"
#include "symbol.h"
#include "dof.h"
#include "perf.h"
#include "symbol_cost.h"

#define PlasmaLeft  141
#define PlasmaRight 142
#include "flops.h"

double flops_zgetrf(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);
double flops_dgetrf(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);
double flops_zpotrf(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);
double flops_dpotrf(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);


void symbCost(PASTIX_INT *iparm, double *dparm, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
  double flops = 0.;
  printf("SymbolCost: number of operations Cholesky %g \n",
          recursive_sum(0, symbmtx->cblknbr-1, cholesky, symbmtx, dofptr));
  printf("SymbolCost: number of operations Crout2t  %g \n",
          recursive_sum(0, symbmtx->cblknbr-1, crout_2t, symbmtx, dofptr));
  printf("SymbolCost: number of operations Crout3t  %g \n",
          recursive_sum(0, symbmtx->cblknbr-1, crout_3t, symbmtx, dofptr));
  printf("SymbolCost: number of operations CroutHyb %g \n",
          recursive_sum(0, symbmtx->cblknbr-1, crout_hyb, symbmtx, dofptr));
  printf("SymbolCost: number of operations CroutHyb blok %g \n",
          recursive_sum(0, symbmtx->cblknbr-1, crout_blok, symbmtx, dofptr));
  printf("SymbolCost: number of non-zero   %g \n",
          recursive_sum(0, symbmtx->cblknbr-1, nnz, symbmtx, dofptr));

  set_iparm(iparm, IPARM_NNZEROS,   (PASTIX_INT)recursive_sum(0, symbmtx->cblknbr-1, nnz,        symbmtx, dofptr));

  if ( iparm[IPARM_FACTORIZATION] == API_FACT_LU ) {
    if ( (iparm[IPARM_FLOAT] == API_COMPLEXDOUBLE) ||
         (iparm[IPARM_FLOAT] == API_COMPLEXSINGLE) ) {
      flops = recursive_sum(0, symbmtx->cblknbr-1, flops_zgetrf, symbmtx, dofptr);
    }
    else {
      flops = recursive_sum(0, symbmtx->cblknbr-1, flops_dgetrf, symbmtx, dofptr);
    }
  } else {
    if ( (iparm[IPARM_FLOAT] == API_COMPLEXDOUBLE) ||
         (iparm[IPARM_FLOAT] == API_COMPLEXSINGLE) ) {
      flops = recursive_sum(0, symbmtx->cblknbr-1, flops_zpotrf, symbmtx, dofptr);
    }
    else {
      flops = recursive_sum(0, symbmtx->cblknbr-1, flops_dpotrf, symbmtx, dofptr);
    }
  }
  set_dparm(dparm, DPARM_FACT_FLOPS, flops);
}



double recursive_sum(PASTIX_INT a, PASTIX_INT b, double (*fval)(PASTIX_INT, const SymbolMatrix *, const Dof *),
                     const SymbolMatrix * symbmtx, const Dof * dofptr)
{
  if(a != b)
    return recursive_sum(        a, (a+b)/2, fval, symbmtx, dofptr)
         + recursive_sum((a+b)/2+1,       b, fval, symbmtx, dofptr);

  return fval(a, symbmtx, dofptr);
}

double crout_2t(PASTIX_INT cblknum, const SymbolMatrix *symbmtx, const Dof * dofptr)
{ PASTIX_INT i;
  double gk = 0;
  double lk = 0;

  /* lk is the dimension of the diagonal blok */
  lk = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

  /* gk is the height of off-diag bloks */
  for(i=symbmtx->cblktab[cblknum].bloknum+1; i<symbmtx->cblktab[cblknum+1].bloknum; i++)
    gk += (double)(symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum +1);

  return( (2*lk*(dofptr->noddval)*lk*(dofptr->noddval)*lk*(dofptr->noddval) + (6*gk*(dofptr->noddval)+3)*lk*(dofptr->noddval)*lk*(dofptr->noddval) + (6*gk*(dofptr->noddval)*gk*(dofptr->noddval)+6*gk*(dofptr->noddval)-5)*lk*(dofptr->noddval))/6);

}


double crout_3t(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{ PASTIX_INT i;
  double gk = 0;
  double lk = 0;

  /* lk is the dimension of the diagonal blok */
  lk = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

  /* gk is the height of off-diag bloks */
  for(i=symbmtx->cblktab[cblknum].bloknum+1;i<symbmtx->cblktab[cblknum+1].bloknum;i++)
    gk += (double)(symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum +1);

  return( (lk*(dofptr->noddval)*lk*(dofptr->noddval)*lk*(dofptr->noddval) + (3*gk*(dofptr->noddval)+1)*lk*(dofptr->noddval)*lk*(dofptr->noddval) + (3*gk*(dofptr->noddval)*gk*(dofptr->noddval)-2+2*gk*(dofptr->noddval))*lk*(dofptr->noddval))/2 );

}


double crout_hyb(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{ PASTIX_INT i;
  double gk = 0;
  double lk = 0;

  /* lk is the dimension of the diagonal blok */
  lk =(double)( symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

  /* gk is the height of off-diag bloks */
  for(i=symbmtx->cblktab[cblknum].bloknum+1;i<symbmtx->cblktab[cblknum+1].bloknum;i++)
    gk += (double)(symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum +1);

  return( (lk*(dofptr->noddval)*lk*(dofptr->noddval)*lk*(dofptr->noddval) + 3*(gk*(dofptr->noddval)+1)*lk*(dofptr->noddval)*lk*(dofptr->noddval) + (3*gk*(dofptr->noddval)*gk*(dofptr->noddval) -4 + 6*gk*(dofptr->noddval))*lk*(dofptr->noddval))/3 );

}

double cholesky(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{ PASTIX_INT i;
  double gk = 0;
  double lk = 0;
#ifdef DOF_CONSTANT
  /* lk is the dimension of the diagonal blok */
  lk = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

  /* gk is the height of off-diag bloks */
  for(i=symbmtx->cblktab[cblknum].bloknum+1;i<symbmtx->cblktab[cblknum+1].bloknum;i++)
    gk += (double)(symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum +1);

  return( (2*lk*(dofptr->noddval)*lk*(dofptr->noddval)*lk*(dofptr->noddval) + (6*gk*(dofptr->noddval)-3)*lk*(dofptr->noddval)*lk*(dofptr->noddval) +(6*gk*(dofptr->noddval)*gk*(dofptr->noddval)+1-6*gk*(dofptr->noddval))*lk*(dofptr->noddval))/6 );
#endif
}

/*******************************************/
/* Number of non zero  extradiagonal terms */
/*******************************************/

double nnz(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{ PASTIX_INT i;
  double gk = 0;
  double lk = 0;
#ifdef DOF_CONSTANT
  /* lk is the dimension of the diagonal blok */
  lk = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

  /* gk is the height of off-diag bloks */
  for(i=symbmtx->cblktab[cblknum].bloknum+1;i<symbmtx->cblktab[cblknum+1].bloknum;i++)
    gk +=(double)( symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum +1);



  return( lk*(dofptr->noddval)*(lk*(dofptr->noddval)+1)/2 + gk*(dofptr->noddval)*lk*(dofptr->noddval) - lk*(dofptr->noddval));
#endif
}


double crout_blok(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
    double l, h, g;
    PASTIX_INT k;
    double nbops = 0;
    h=0;
    /** we need the height of cblk non empty lines  and the broadness
      of the cbl to compute the local compute cost **/
    l = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);
    g = 0;
    for(k=symbmtx->cblktab[cblknum].bloknum;k<symbmtx->cblktab[cblknum+1].bloknum;k++)
      g += (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);

    /** retrieve diag height so let g be the odb non empty lines height **/
    g -= l;

    /** compute the local compute cost **/
#ifdef DEBUG_BLEND
    ASSERT(l>0,MOD_BLEND);
#endif

#ifdef DOF_CONSTANT
    nbops = (double)(OPS_PPF(l*(dofptr->noddval)));
    if(g>0)
      nbops += (double)(OPS_TRSM(l*(dofptr->noddval),g*(dofptr->noddval))) + l*(double)(OPS_SCAL(g*(dofptr->noddval)));

    /** compute for each odb its contribution compute cost and add cost **/
    for(k=symbmtx->cblktab[cblknum].bloknum+1;k<symbmtx->cblktab[cblknum+1].bloknum;k++)
      {
        h = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
        /* g is the odb lines number above this odb (odb lines include)*/
        nbops += /*l*(double)(OPS_SCAL(g)) +*/(double)(OPS_GEMM(l*(dofptr->noddval),g*(dofptr->noddval),h*(dofptr->noddval))) + (double)(OPS_GEAM(g*(dofptr->noddval),h*(dofptr->noddval)));
#ifdef DEBUG_BLEND
        ASSERT(nbops>=0,MOD_BLEND);
#endif
        g -= h;
      }
#endif
    return nbops;
}

double flops_zgetrf(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
  double M, N, K;
  PASTIX_INT k;
  double nbops = 0.;

  /*
   * Size of the factorization kernel (square)
   */
  N = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

  /*
   * Height of the TRSM to which apply the TRSM
   */
  M = 0;
  for(k = symbmtx->cblktab[cblknum].bloknum+1;
      k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
      M += (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
    }

#ifdef DOF_CONSTANT
  N *= (double)dofptr->noddval;
  M *= (double)dofptr->noddval;
#endif

  nbops  = FLOPS_ZGETRF( N, N );
  nbops += 2. * FLOPS_ZTRSM( PlasmaRight, M, N );

  /*
   * Compute the cost of each GEMM
   */
  K = N;
  for(k = symbmtx->cblktab[cblknum].bloknum+1;
      k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
      N = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);

#ifdef DOF_CONSTANT
      N *= (double)dofptr->noddval;
#endif

      nbops += 2. * FLOPS_ZGEMM( M, N, K );

      M -= N;
    }

  return nbops;
}

double flops_dgetrf(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
  double M, N, K;
  PASTIX_INT k;
  double nbops = 0.;

  /*
   * Size of the factorization kernel (square)
   */
  N = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

  /*
   * Height of the TRSM to which apply the TRSM
   */
  M = 0;
  for(k = symbmtx->cblktab[cblknum].bloknum+1;
      k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
      M += (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
    }

#ifdef DOF_CONSTANT
  N *= (double)dofptr->noddval;
  M *= (double)dofptr->noddval;
#endif

  nbops  = FLOPS_DGETRF( N, N );
  nbops += 2. * FLOPS_DTRSM( PlasmaRight, M, N );

  /*
   * Compute the cost of each GEMM
   */
  K = N;
  for(k = symbmtx->cblktab[cblknum].bloknum+1;
      k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
      N = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);

#ifdef DOF_CONSTANT
      N *= (double)dofptr->noddval;
#endif

      nbops += 2. * FLOPS_DGEMM( M, N, K );

      M -= N;
    }

  return nbops;
}

double flops_zpotrf(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
  double M, N, K;
  PASTIX_INT k;
  double nbops = 0.;

  /*
   * Size of the factorization kernel (square)
   */
  N = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

  /*
   * Height of the TRSM to which apply the TRSM
   */
  M = 0;
  for(k = symbmtx->cblktab[cblknum].bloknum+1;
      k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
      M += (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
    }

#ifdef DOF_CONSTANT
  N *= (double)dofptr->noddval;
  M *= (double)dofptr->noddval;
#endif

  nbops  = FLOPS_ZPOTRF( N );
  nbops += FLOPS_ZTRSM( PlasmaRight, M, N );

  /*
   * Compute the cost of each GEMM
   */
  K = N;
  for(k = symbmtx->cblktab[cblknum].bloknum+1;
      k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
      N = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);

#ifdef DOF_CONSTANT
      N *= (double)dofptr->noddval;
#endif

      nbops += FLOPS_ZGEMM( M, N, K );

      M -= N;
    }

  return nbops;
}

double flops_dpotrf(PASTIX_INT cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
  double M, N, K;
  PASTIX_INT k;
  double nbops = 0.;

  /*
   * Size of the factorization kernel (square)
   */
  N = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

  /*
   * Height of the TRSM to which apply the TRSM
   */
  M = 0;
  for(k = symbmtx->cblktab[cblknum].bloknum+1;
      k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
      M += (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
    }

#ifdef DOF_CONSTANT
  N *= (double)dofptr->noddval;
  M *= (double)dofptr->noddval;
#endif

  nbops  = FLOPS_DPOTRF( N );
  nbops += FLOPS_DTRSM( PlasmaRight, M, N );

  /*
   * Compute the cost of each GEMM
   */
  K = N;
  for(k = symbmtx->cblktab[cblknum].bloknum+1;
      k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
      N = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);

#ifdef DOF_CONSTANT
      N *= (double)dofptr->noddval;
#endif

      nbops += FLOPS_DGEMM( M, N, K );

      M -= N;
    }

  return nbops;
}
