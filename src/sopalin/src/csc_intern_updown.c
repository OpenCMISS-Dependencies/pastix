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
  File: csc_intern_updown.c

  Build UpDownVector from user right-hand-side and CSCd.
  Retrieve soltion from UpDownVector.
  Construct UpDownVector such as X[i] = 1, or X[i] = i.

*/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#include "common_pastix.h"
#include "tools.h"
#include "order.h"
#include "csc.h"
#include "updown.h"


#include "ftgt.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"

#include "csc_intern_updown.h"

#ifdef DEBUG_RAFF
#define CSC_LOG
#endif

#ifdef CPLX
#define SMX_SOL (1.0+0.0*I)
#else
#define SMX_SOL 1.0
#endif

/*
  Function: CscdUpdownRhs

  Fill-in UpDownVector structure from user right-hand-side member.

  Parameters:
    updovct - UpDownVector structure to fill-in.
    solvmtx - Solver matrix.
    rhs     - Right-hand-side member.
    invp    - reverse permutation tabular.
    dof     - Number of degree of freedom.
 */
void CscUpdownRhs(UpDownVector       *updovct,
                  const SolverMatrix *solvmtx,
                  const PASTIX_FLOAT        *rhs,
                  const PASTIX_INT          *invp,
                  int                 dof)
{
  PASTIX_INT itercblk;
  PASTIX_INT itercol;
  PASTIX_INT itersm2x;
  PASTIX_INT indice;
  PASTIX_INT i;

  print_debug(DBG_CSC_LOG, "-> CscUpdownRhs \n");

  for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;
      for (itercol=solvmtx->cblktab[itercblk].fcolnum;
           itercol<solvmtx->cblktab[itercblk].lcolnum+1;
           itercol++)
        {
          indice = invp[(itercol-itercol%dof)/dof] *dof + itercol%dof;
          for (i=0; i<updovct->sm2xnbr; i++)
            {
              updovct->sm2xtab[itersm2x + i*updovct->sm2xsze] =
                rhs[indice + i*updovct->gnodenbr];
            }
          itersm2x++;
        }
    }

  print_debug(DBG_CSC_LOG, "<- CscUpdownRhs \n");
}

/*
  Function: CscdUpdownRhs

  Fill-in UpDownVector structure from user distributed right-hand-side member.

  Parameters:
    updovct - UpDownVector structure to fill-in.
    solvmtx - Solver matrix.
    rhs     - Right-hand-side member.
    invp    - reverse permutation tabular.
    g2l     - local numbers of global nodes, if not local contains -owner
    dof     - Number of degree of freedom.
 */
void CscdUpdownRhs(UpDownVector       *updovct,
                   const SolverMatrix *solvmtx,
                   const PASTIX_FLOAT        *rhs,
                   const PASTIX_INT          *invp,
                   const PASTIX_INT          *g2l,
                   const PASTIX_INT           ln,
                   int                 dof)
{
  PASTIX_INT itercblk;
  PASTIX_INT itercol;
  PASTIX_INT itersm2x;
  PASTIX_INT indice;
  PASTIX_INT i;

  print_debug(DBG_CSC_LOG, "-> CscdUpdownRhs \n");

  for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;
      for (itercol=solvmtx->cblktab[itercblk].fcolnum;
           itercol<solvmtx->cblktab[itercblk].lcolnum+1;
           itercol++)
        {
          /* No problem, Column is local */
          indice = (g2l[invp[(itercol - itercol%dof)/dof]]-1)*dof + itercol%dof;
          for (i=0; i<updovct->sm2xnbr; i++)
            {
              updovct->sm2xtab[itersm2x + i*updovct->sm2xsze] =
                rhs[indice + i*ln];
            }
          itersm2x++;
        }
    }

  print_debug(DBG_CSC_LOG, "<- CscdUpdownRhs \n");
}

/*
  Function:CscRhsUpdown

  Builds solution from UpDownVector structure

  Parameters:
    updovct  - UpDownVector structure containing the solution.
    solvmtx  - Solver matrix structure.
    rhs      - Solution to fill.
    ncol     - Number of columns in local matrix.
    dof      - Number of degree of freedom.
    comm     - MPI communicator.

 */
void CscRhsUpdown(const UpDownVector *updovct,
                  const SolverMatrix *solvmtx,
                  PASTIX_FLOAT              *rhs,
                  const PASTIX_INT           ncol,
                  const PASTIX_INT          *invp,
                  const int           dof,
                  const int           rhsmaking,
                  MPI_Comm            comm)
{
  PASTIX_INT    iter;
  PASTIX_INT    itercblk;
  PASTIX_INT    itersm2x;
  PASTIX_INT    itercol;
  PASTIX_INT    indice, i;
  PASTIX_INT    size = updovct->sm2xnbr*ncol*dof;
  PASTIX_FLOAT *rhs2 = NULL;
  (void)comm;

#ifdef INOUT_ALLREDUCE
  rhs2 = rhs;
#else
  MALLOC_INTERN(rhs2, size, PASTIX_FLOAT);
#endif

  print_debug(DBG_CSC_LOG, "-> CscRhsUpdown \n");

  for (iter=0; iter<size; iter++)
    {
      rhs2[iter] = 0.0;
#ifndef INOUT_ALLREDUCE
      rhs[iter]  = 0.0;
#endif
    }

  if (rhsmaking == API_RHS_B)
    {
      for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
        {
          itersm2x = updovct->cblktab[itercblk].sm2xind;
          for (itercol=solvmtx->cblktab[itercblk].fcolnum;
               itercol<solvmtx->cblktab[itercblk].lcolnum+1;
               itercol++)
            {
              indice = invp[(itercol-itercol%dof)/dof] *dof + itercol%dof;
              for (i=0; i<updovct->sm2xnbr; i++)
                {
                  rhs2[indice + i*updovct->gnodenbr] =
                    updovct->sm2xtab[itersm2x + i*updovct->sm2xsze];
                }
              itersm2x++;
            }
        }
    }
  else /* Rhs making */
    {
      for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
        {
          itersm2x = updovct->cblktab[itercblk].sm2xind;
          for (itercol=solvmtx->cblktab[itercblk].fcolnum;
               itercol<solvmtx->cblktab[itercblk].lcolnum+1;
               itercol++)
            {
              for (i=0; i<updovct->sm2xnbr; i++)
                {
                  rhs2[itercol + i*updovct->gnodenbr] =
                    updovct->sm2xtab[itersm2x + i*updovct->sm2xsze];
                }
              itersm2x++;
            }
        }
    }
  MPI_Allreduce((void *) rhs2, (void *) rhs, size,
                COMM_FLOAT, COMM_SUM, comm);

#ifndef INOUT_ALLREDUCE
  memFree_null(rhs2);
#endif

  print_debug(DBG_CSC_LOG, "<- CscRhsUpdown \n");
}

/*
  Function:CscdRhsUpdown

  Builds distributed solution from
  UpDownVector structure

  Parameters:
    updovct  - UpDownVector structure containing the solution.
    solvmtx  - Solver matrix structure.
    x        - Solution to fill.
    ncol     - Number of columns in local matrix.
    g2l      - local numbers of global nodes, if not local contains -owner
    ord      - ordering
    dof      - Number of degree of freedom.
    comm     - MPI communicator.

 */
void CscdRhsUpdown(const UpDownVector *updovct,
                   const SolverMatrix *solvmtx,
                   PASTIX_FLOAT              *x,
                   const PASTIX_INT           ncol,
                   const PASTIX_INT          *g2l,
                   const PASTIX_INT          *invp,
                   int                 dof,
                   MPI_Comm            comm)
{
  PASTIX_INT iter;
  PASTIX_INT itercblk;
  PASTIX_INT itersm2x;
  PASTIX_INT itercol;
  PASTIX_INT size = updovct->sm2xnbr*ncol*dof;
  (void)comm;

  print_debug(DBG_CSC_LOG, "-> CscdRhsUpdown \n");

  for (iter=0; iter<size; iter++)
    {
      x[iter]  = 0.0;
    }

  for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;
      for (itercol=solvmtx->cblktab[itercblk].fcolnum;
           itercol<solvmtx->cblktab[itercblk].lcolnum+1;
           itercol++)
        {
          PASTIX_INT i;
          for (i=0; i<updovct->sm2xnbr; i++)
            {
              x[(g2l[invp[(itercol-itercol%dof)/dof]] - 1)*dof +
                itercol%dof + i*ncol] = updovct->sm2xtab[itersm2x+i*
                                                         updovct->sm2xsze];
            }
          itersm2x++;
        }
    }

  print_debug(DBG_CSC_LOG, "<- CscdRhsUpdown \n");
}

/*
  Function: Csc2updown

  Fill-in UpDownVector structure such as the solution of
  the system Ax=b is x_i=1 (API_RHS_1) or x_i=i (API_RHS_I).

  Parameters:
    cscmtx   - internal CSCd matrix.
    updovct  - UpDownVector structure to fill-in.
    solvmtx  - Solver matrix.
    mode     - wanted solution API_RHS_1 or API_RHS_I.
    comm     - MPI communicator.
*/
void Csc2updown(const CscMatrix    *cscmtx,
                UpDownVector       *updovct,
                const SolverMatrix *solvmtx,
                int                 mode,
                MPI_Comm            comm)
{
  PASTIX_INT    itercblk;
  PASTIX_INT    itercol;
  PASTIX_INT    itertempy;
  PASTIX_INT    iterval;
  PASTIX_INT    itersmx;
  PASTIX_INT    cblknbr;
  PASTIX_FLOAT *smb   = NULL;
  PASTIX_FLOAT *tempy = NULL;
  (void)comm;

  print_debug(DBG_CSC_LOG, "-> Csc2updown \n");

  cblknbr = solvmtx->cblknbr;

  MALLOC_INTERN(tempy, updovct->gnodenbr, PASTIX_FLOAT);
#ifdef INOUT_ALLREDUCE
  smb = tempy;
#else
  MALLOC_INTERN(smb,   updovct->gnodenbr, PASTIX_FLOAT);
#endif

  print_debug(DBG_CSC_LOG, "nodenbr=%ld\n",(long)updovct->gnodenbr);

  for (itersmx=0; itersmx<updovct->sm2xnbr; itersmx++)
    {
      for (itertempy=0; itertempy<updovct->gnodenbr; itertempy++)
        {
          tempy[itertempy] = 0.0;
#ifndef INOUT_ALLREDUCE
          smb[itertempy]   = 0.0;
#endif
        }

      for (itercblk=0; itercblk < cblknbr; itercblk++)
        {
          PASTIX_INT colnbr   = solvmtx->cblktab[itercblk].lcolnum - solvmtx->cblktab[itercblk].fcolnum +1;

          for (itercol=0; itercol < colnbr; itercol++)
            {
              PASTIX_INT colvalidx  = CSC_COL(cscmtx,itercblk,itercol);
              PASTIX_INT ncolvalidx = CSC_COL(cscmtx,itercblk,itercol+1);

              for (iterval=colvalidx; iterval<ncolvalidx; iterval++)
                {
                  switch (mode)
                    {
                    case API_RHS_1:
                      tempy[CSC_ROW(cscmtx,iterval)] += (PASTIX_FLOAT)(itersmx+SMX_SOL)*CSC_VAL(cscmtx,iterval);
                      break;
                    case API_RHS_I:
                      tempy[CSC_ROW(cscmtx,iterval)] +=
                        ((PASTIX_FLOAT)(solvmtx->cblktab[itercblk].fcolnum+itercol))*
                        CSC_VAL(cscmtx,iterval);
                      break;
                    }
                }
            }
        }

      MPI_Allreduce((void *) tempy, (void *) smb, updovct->gnodenbr,
                    COMM_FLOAT, COMM_SUM, comm);

      for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
        {
          PASTIX_INT iterdval = updovct->cblktab[itercblk].sm2xind+itersmx*updovct->sm2xsze;

          for (iterval=0; iterval<CSC_COLNBR(cscmtx,itercblk); iterval++)
            {
              updovct->sm2xtab[iterdval+iterval] = smb[solvmtx->cblktab[itercblk].fcolnum+iterval];
            }
        }
    }

  memFree_null(tempy);
#ifndef INOUT_ALLREDUCE
  memFree_null(smb);
#endif

  print_debug(DBG_CSC_LOG, "<- Csc2updown \n");
}


/*
  Function: Csc2updown_X0

  Fill-in initial X0 for reffinement if we don't want to use
  Solve step.

  (iparm[IPARM_ONLY_RAFF] == API_YES)

  Parameters:
    updovct - UpDownVector structure were to copy B as the first X0 used for raffinement.
    solvmtx - Solver matrix.
    mode    - Rule to construct X0 (API_RHS_0 : X0[i] = 0, API_RHS_1 : X0[i] = 1, API_RHS_I : X0[i] = i).
    comm    - MPI_Communicator.
*/
void Csc2updown_X0(UpDownVector *updovct,
                   /*const*/ SolverMatrix *solvmtx,
                   int mode,
                   MPI_Comm comm)
{
  PASTIX_INT  itercblk;
  PASTIX_INT  iterval;
  PASTIX_INT  itersmx;
  PASTIX_INT  cblknbr = solvmtx->cblknbr;
  (void)comm;

  print_debug(DBG_CSC_LOG, "-> Csc2updown_X0 \n");
  print_debug(DBG_CSC_LOG, "nodenbr=%ld\n",(long)updovct->gnodenbr);

  for (itersmx=0; itersmx<updovct->sm2xnbr; itersmx++)
    {
      for (itercblk=0; itercblk < cblknbr; itercblk++)
        {
          PASTIX_INT colnbr   = solvmtx->cblktab[itercblk].lcolnum - solvmtx->cblktab[itercblk].fcolnum +1;
          PASTIX_INT iterdval = updovct->cblktab[itercblk].sm2xind+itersmx*updovct->sm2xsze;

          for (iterval=0; iterval < colnbr; iterval++)
            {
              switch (mode)
                {
                case API_RHS_0:
                  updovct->sm2xtab[iterdval+iterval] = (PASTIX_FLOAT)0.0;
                  break;
                case API_RHS_1:
                  updovct->sm2xtab[iterdval+iterval] = (PASTIX_FLOAT)SMX_SOL;
                  break;
                case API_RHS_I:
                  updovct->sm2xtab[iterdval+iterval] = ((PASTIX_FLOAT)(solvmtx->cblktab[itercblk].fcolnum+iterval));
                  break;

                }
            }
        }
    }

  print_debug(DBG_CSC_LOG, "<- Csc2updown_X0 \n");
}
