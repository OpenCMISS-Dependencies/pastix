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
/******************************************************************************
 * File: csc_intern_build.c                                                   *
 ******************************************************************************
 *                                                                            *
 * Functions to build internal CSCd from user CSCd.                           *
 *                                                                            *
 * Function to free internal CSCd.                                            *
 *                                                                            *
 ******************************************************************************/

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
#include "ftgt.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"

#include "sopalin_define.h"
#include "csc_intern_build.h"

/****                    Section: Macros and defines                       ****/

/******************************************************************************
 * Macro: COL_IS_LOCAL(col)                                                   *
 ******************************************************************************
 * Uses the result of glob2loc to determine if the                            *
 * column is local.                                                           *
 *                                                                            *
 * Parameters:                                                                *
 *   col - glob2loc value.                                                    *
 *                                                                            *
 ******************************************************************************/
#define COL_IS_LOCAL(col) (col > 0)


/******************************************************************************
 * Macro: SET_CSC_COL(newcoltab)                                              *
 ******************************************************************************
 * Allocate and fill CSC_COLTAB using newcoltab array.                        *
 *                                                                            *
 * CSC_COLNBR(thecsc,index) will contain the number of column of the          *
 * column block "index" of thecsc                                             *
 *                                                                            *
 * CSC_COLTAB(thecsc, index) is the tabular containing the starting index     *
 * of each column of the column bloc in CSC_ROWTAB and CSC_VALTAB             *
 *                                                                            *
 * CSC_COL(thecsc, index, iter) is the starting index for the column iter     *
 * of the block index.                                                        *
 *                                                                            *
 * Parameters:                                                                *
 *   newcoltab - array containing starting index of each column in rows and   *
 *               values, using new ordering.                                  *
 *                                                                            *
 ******************************************************************************/
#define SET_CSC_COL(newcoltab) do {                                     \
    /* Local coltab */                                                  \
    CSC_FNBR(thecsc) = solvmtx->cblknbr;                                \
    MALLOC_INTERN(CSC_FTAB(thecsc), CSC_FNBR(thecsc), CscFormat);       \
                                                                        \
    for (index=0; index<solvmtx->cblknbr; index++)                      \
    {                                                                   \
      PASTIX_INT fcolnum = solvmtx->cblktab[index].fcolnum;                    \
      PASTIX_INT lcolnum = solvmtx->cblktab[index].lcolnum;                    \
      CSC_COLNBR(thecsc,index) = (lcolnum - fcolnum+1);                 \
                                                                        \
      MALLOC_INTERN(CSC_COLTAB(thecsc,index),                           \
                    CSC_COLNBR(thecsc,index)+1, PASTIX_INT);                   \
                                                                        \
      if ((fcolnum)%dof != 0)                                           \
        errorPrint("dof doesn't divide fcolnum");                       \
                                                                        \
      colsize = 0;                                                      \
      for (iter=0; iter<(CSC_COLNBR(thecsc,index)+1); iter++)           \
      {                                                                 \
        /* fcolnum %dof = 0 */                                          \
        nodeidx = (fcolnum+(iter-iter%dof))/dof;                        \
        if (g2l != NULL &&                                              \
            iter != CSC_COLNBR(thecsc,index) &&                         \
            !COL_IS_LOCAL(g2l[ord->peritab[nodeidx]]))                  \
        {                                                               \
          errorPrint("Columns in internal CSCD must be in given CSCD"); \
        }                                                               \
                                                                        \
        CSC_COL(thecsc,index,iter) = colsize+ strdcol;                  \
        if (iter < CSC_COLNBR(thecsc,index))                            \
        {                                                               \
          colsize = (newcoltab[nodeidx+1] -                             \
                     newcoltab[nodeidx])*dof;                           \
        }                                                               \
        strdcol =  CSC_COL(thecsc,index,iter);                          \
      }                                                                 \
                                                                        \
      /* strdcol <- colptr[n] for the index_th column block */          \
      /*   ie : the number of element in the column block */            \
                                                                        \
      strdcol = CSC_COL(thecsc,index,CSC_COLNBR(thecsc,index));         \
    }                                                                   \
    memFree_null(newcoltab);                                            \
  } while (0)

/******************************************************************************
 * Macro: CSC_ALLOC                   *
 ******************************************************************************
 *                        *
 * Allocates CSC_ROWTAB, CSC_VALTAB and the transpose CSC arrays.             *
 *                        *
 ******************************************************************************/
#define CSC_ALLOC do {                                            \
                                                                  \
    /* Puting values */                                           \
    MALLOC_INTERN(CSC_ROWTAB(thecsc), strdcol, PASTIX_INT);              \
    MALLOC_INTERN(CSC_VALTAB(thecsc), strdcol, PASTIX_FLOAT);            \
                                                                  \
    if (transcsc != NULL)                                         \
    {                                                             \
      if (Type[1] == 'S' || Type[1] == 'H')                       \
      {                                                           \
        if (forcetrans == API_YES)                                \
        {                                                         \
          (*transcsc) = CSC_VALTAB(thecsc);                       \
        }                                                         \
      }                                                           \
      else                                                        \
      {                                                           \
        MALLOC_INTERN(*transcsc, strdcol, PASTIX_FLOAT);                 \
        MALLOC_INTERN(trowtab, strdcol, PASTIX_INT);                     \
        MALLOC_INTERN(trscltb, solvmtx->cblknbr, PASTIX_INT *);          \
                                                                  \
        for (index=0; index<solvmtx->cblknbr; index++)            \
        {                                                         \
          MALLOC_INTERN(trscltb[index],                           \
                        CSC_COLNBR(thecsc,index)+1, PASTIX_INT);         \
          for (iter=0; iter<(CSC_COLNBR(thecsc,index)+1); iter++) \
          {                                                       \
            trscltb[index][iter] = CSC_COL(thecsc,index,iter);    \
          }                                                       \
        }                                                         \
      }                                                           \
    }                                                             \
  } while (0)

/******************************************************************************
 * Macro: SET_CSC_ROW_VAL(itercblk, therow, thecol, val)          *
 ******************************************************************************
 *                        *
 * Fill next CSC_ROW and CSC_VAL.               *
 *                        *
 * Parameters:                      *
 *   itercblk - Column block index.               *
 *   therow   - Row index.                  *
 *   thecol   - Column index.                 *
 *   val      - Value array.                  *
 *                        *
 ******************************************************************************/
#define SET_CSC_ROW_VAL(itercblk, therow, thecol, val)        \
  _set_csc_row_val(solvmtx, thecsc, itercblk, therow, thecol, \
                   dof, iterdofcol, iterdofrow, iter,         \
                   colidx, strdcol, val)
static inline
void _set_csc_row_val(const SolverMatrix *solvmtx,
                      CscMatrix          *thecsc,
                      PASTIX_INT          itercblk,
                      PASTIX_INT          therow,
                      PASTIX_INT          thecol,
                      PASTIX_INT          dof,
                      PASTIX_INT          iterdofcol,
                      PASTIX_INT          iterdofrow,
                      PASTIX_INT          iter,
                      PASTIX_INT          colidx,
                      PASTIX_INT          strdcol,
                      PASTIX_FLOAT       *val) {
  PASTIX_INT fcolnum = solvmtx->cblktab[itercblk].fcolnum;
  PASTIX_INT validx  = (iter-1)*dof*dof + iterdofcol*dof + iterdofrow;

  colidx = CSC_COL(thecsc, itercblk, thecol - fcolnum);

  if (strdcol <= colidx)
    {
      errorPrint("%s:%d colidx %ld >= strdcol %ld",
                 __FILE__, __LINE__, (long)colidx, (long)strdcol);

      EXIT(MOD_SOPALIN, UNKNOWN_ERR);
    }

  CSC_ROW(thecsc, colidx) = therow;
  CSC_VAL(thecsc, colidx) = val[validx];
  CSC_COL(thecsc, itercblk, thecol-fcolnum)++;

}

#define SET_CSC_ROW_VAL_CONJ(itercblk, therow, thecol, val) do {    \
    PASTIX_INT fcolnum = solvmtx->cblktab[itercblk].fcolnum;               \
    PASTIX_INT validx  = (iter-1)*dof*dof + iterdofcol*dof + iterdofrow;   \
                                                                    \
    colidx = CSC_COL(thecsc, itercblk, thecol - fcolnum);           \
                                                                    \
    if (strdcol <= colidx)                                          \
    {                                                               \
      errorPrint("%s:%d colidx %ld >= strdcol %ld",                 \
                 __FILE__, __LINE__, (long)colidx, (long)strdcol);  \
                                                                    \
      EXIT(MOD_SOPALIN, UNKNOWN_ERR);                               \
    }                                                               \
                                                                    \
    CSC_ROW(thecsc, colidx) = therow;                               \
    CSC_VAL(thecsc, colidx) = CONJ_FLOAT(val[validx]);              \
    CSC_COL(thecsc, itercblk, thecol-fcolnum)++;                    \
  } while (0)

/******************************************************************************
 * Macro: SET_TRANS_ROW_VAL(itercblk, therow, thecol, val)          *
 ******************************************************************************
 *                        *
 * Fill next entries for transpose CSC.             *
 *                        *
 * Parameters:                      *
 *   itercblk - Column block index.               *
 *   therow   - Row index.                  *
 *   thecol   - Column index.                 *
 *   val      - Value array.                  *
 *                        *
 ******************************************************************************/
#define SET_TRANS_ROW_VAL(itercblk, therow, thecol, val) do {     \
    PASTIX_INT fcolnum = solvmtx->cblktab[itercblk].fcolnum;             \
    PASTIX_INT trsidx  = trscltb[itercblk][therow-fcolnum];              \
    PASTIX_INT validx  = (iter-1)*dof*dof + iterdofcol*dof + iterdofrow; \
                                                                  \
    (*transcsc)[trsidx] = val[validx];                            \
    trowtab[trsidx] = thecol;                                     \
                                                                  \
    ASSERTDBG(therow >= fcolnum, MOD_SOPALIN);                    \
    trscltb[itercblk][therow -fcolnum]++;                         \
  } while (0)
/******************************************************************************
 * Macro: CSC_SORT                    *
 ******************************************************************************
 *                        *
 * Sort the internal CSC                  *
 *                        *
 ******************************************************************************/
#define CSC_SORT do {                                                  \
    /* Sort */                                                         \
    for (index=0; index<solvmtx->cblknbr; index++) {                   \
      for (iter=0; iter<CSC_COLNBR(thecsc,index); iter++) {            \
        PASTIX_INT   *t = &(CSC_FROW(thecsc,index,iter));              \
        PASTIX_FLOAT *v = &(CSC_FVAL(thecsc,index,iter));              \
        PASTIX_INT    n = CSC_COL(thecsc,index,iter+1)-                \
          CSC_COL(thecsc,index,iter);                                  \
        PASTIX_INT    ndof2 = 1; /* internal CSC is with one DoF */    \
        void * sortptr[3];                                             \
        sortptr[0] = t;                                                \
        sortptr[1] = v;                                                \
        sortptr[2] = &ndof2;                                           \
        qsortIntFloatAsc(sortptr, n);                                  \
      }                                                                \
    }                                                                  \
    if (transcsc != NULL && (*transcsc) != NULL) {                     \
      if (!forcetrans) {                                               \
        for (index=0; index<solvmtx->cblknbr; index++) {               \
          for (iter=0; iter<CSC_COLNBR(thecsc,index); iter++) {        \
            PASTIX_INT   *t = &(trowtab[CSC_COL(thecsc,index,iter)]);  \
            PASTIX_FLOAT *v = &((*transcsc)[CSC_COL(thecsc,            \
                                                    index,iter)]);     \
            PASTIX_INT n;                                              \
            PASTIX_INT    ndof2 = 1; /* internal CSC is with 1 DoF */  \
            void * sortptr[3];                                         \
                                                                       \
            n = CSC_COL(thecsc,index,iter+1) -                         \
              CSC_COL(thecsc,index,iter);                              \
                                                                       \
            sortptr[0] = t;                                            \
            sortptr[1] = v;                                            \
            sortptr[2] = &ndof2;                                       \
            qsortIntFloatAsc(sortptr, n);                              \
          }                                                            \
        }                                                              \
        memFree_null(trowtab);                                         \
      }                                                                \
    }                                                                  \
  }while (0)

/****                        Section: Functions                            ****/

/******************************************************************************
 * Function: CscOrdistrib                                                     *
 ******************************************************************************
 *                                                                            *
 * Fill in *thecsc* CSC matrix in column block representation.                *
 *                                                                            *
 * Parameters:                                                                *
 *   thecsc     - Matrix in block column CSC format to fill in.               *
 *   Type       - 3 charactères for matrix type : only Type[1] is             *
 *                used to check if matrix is Symetric(S) or not(U).           *
 *   transcsc   - transpose of the CSC in non symetric mode.                  *
 *   ord        - ordering                                                    *
 *   Nrow       - Number of rows.                                             *
 *   Ncol       - Number of columns.                                          *
 *   Nnzero     - Number of non zeros in the matrix.                          *
 *   colptr     - Index in *rowind* and *val* of the start of each column.    *
 *   rowind     - Index of the elements.                                      *
 *   val        - values of the elements.                                     *
 *   forcetrans - If matrix symetric, transcsc will be the copy of the        *
 *                CSC_VALTAB.                                                 *
 *   solvmtx    - Solver matrix                                               *
 *   procnum    - MPI process number                                          *
 *   dof        - Number of degree of freedom.                                *
 *                                                                            *
 ******************************************************************************/
void CscOrdistrib(CscMatrix          *thecsc,
                  char               *Type,
                  PASTIX_FLOAT             **transcsc,
                  const Order        *ord,
                  PASTIX_INT                 Nrow,
                  PASTIX_INT                 Ncol,
                  PASTIX_INT                 Nnzero,
                  PASTIX_INT                *colptr,
                  PASTIX_INT                *rowind,
                  PASTIX_FLOAT              *val,
                  PASTIX_INT                 forcetrans,
                  const SolverMatrix *solvmtx,
                  PASTIX_INT                 procnum,
                  PASTIX_INT                 dof)
{
  PASTIX_INT   index, itercol, newcol, iter, rowp1,colidx;
  PASTIX_INT   itercblk;
  PASTIX_INT   itercblk2;
  PASTIX_INT  *globcoltab = NULL;
  PASTIX_INT   strdcol    = 0;
  PASTIX_INT **trscltb    = NULL;
  PASTIX_INT  *trowtab    = NULL;
  PASTIX_INT  *cachetab   = NULL;
  PASTIX_INT   therow;
  PASTIX_INT   iterdofrow;
  PASTIX_INT   iterdofcol;
  PASTIX_INT   nodeidx;
  PASTIX_INT   colsize;
  /* To use common macro with CscdOrdistrib */
  PASTIX_INT  *g2l        = NULL;
  (void)Nrow; (void)Nnzero; (void)procnum;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscOrdistrib \n");
#endif
  thecsc->type = Type[1];
  /* Global coltab */
  MALLOC_INTERN(globcoltab,Ncol+1, PASTIX_INT);

  /* globcoltab will contain the number of element in each column
   of the symetrized matrix in the new ordering */

#if (DBG_SOPALIN_TIME==1)
  Clock clk;
  clockInit(&clk);
  clockStart(&clk);
#endif

  for(index=0; index<(Ncol+1); index++)
    globcoltab[index] = 0;

  for (itercol=0; itercol<Ncol; itercol++)
  {
    newcol = ord->permtab[itercol];

    globcoltab[newcol] += colptr[itercol+1] - colptr[itercol];

    if (Type[1] == 'S' || Type[1] == 'H')
    {
      for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
      {
        if ((rowind[iter-1]-1) != itercol)
        {
          newcol = ord->permtab[rowind[iter-1]-1];
          (globcoltab[newcol])++;
        }
      }
    }
  }

#if (DBG_SOPALIN_TIME==1)
  clockStop(&(clk));
  fprintf(stdout, "CscOrdistrib step 1 : %.3g s\n",
          (double)clockVal(&clk));
  clockInit(&clk);
  clockStart(&clk);
#endif

  /* Now, globcoltab will contain starting index of each
   column of rows and values in new ordering */
  newcol = 0;
  for (index=0; index<(Ncol+1); index++)
  {
    colidx = globcoltab[index];
    globcoltab[index] = newcol;
    newcol += colidx;
  }

  SET_CSC_COL(globcoltab);

#if (DBG_SOPALIN_TIME==1)
  clockStop(&(clk));
  fprintf(stdout, "CscOrdistrib step 2 : %.3g s\n",
          (double)clockVal(&clk));
  clockInit(&clk);
  clockStart(&clk);
#endif

  CSC_ALLOC;

  /* Tab: contain the column block number or -1 if not local */
  MALLOC_INTERN(cachetab, (Ncol+1)*dof, PASTIX_INT);
  for (itercol=0; itercol<(Ncol+1)*dof; itercol++)
    cachetab[itercol] = -1;
  for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
  {
    for (itercol=solvmtx->cblktab[itercblk].fcolnum;
         itercol<solvmtx->cblktab[itercblk].lcolnum+1;
         itercol++)
    {
      cachetab[itercol] = itercblk;
    }
  }

  /* Filling in thecsc with values and rows*/
  for (itercol=0; itercol<Ncol; itercol++)
  {
    itercblk = cachetab[ord->permtab[itercol]*dof];

    /* ok put the value */
    for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
    {
      for (iterdofcol = 0; iterdofcol < dof; iterdofcol++)
      {
        for (iterdofrow = 0; iterdofrow < dof; iterdofrow++)
        {

          rowp1 = rowind[iter-1]-1;
          therow = ord->permtab[rowp1]*dof + iterdofrow;
          newcol = ord->permtab[itercol]*dof+iterdofcol;

          if (itercblk != -1)
          {
            SET_CSC_ROW_VAL(itercblk, therow, newcol, val);
          }

          itercblk2 = cachetab[therow];

          if (itercblk2 != -1)
          {
            switch (Type[1])
            {
            case 'S':
            {
              if (rowp1 != itercol)
              {
                /* newcol <-> therow */
                SET_CSC_ROW_VAL(itercblk2, newcol, therow,
                                val);
              }
              break;
            }
            case 'H':
            {
              if (rowp1 != itercol)
              {
                /* newcol <-> therow */
                SET_CSC_ROW_VAL_CONJ(itercblk2, newcol, therow,
                                     val);
              }
              break;
            }
            case 'U':
            {
              if (transcsc != NULL)
              {
                SET_TRANS_ROW_VAL(itercblk2, therow, newcol,
                                  val);
              }
              break;
            }
            default:
              errorPrint("Unknown Matrix Type");
              EXIT(MOD_SOPALIN, UNKNOWN_ERR);
            }
          }
        }
      }
    }
  }
  memFree_null(cachetab);

  /*
   memFree_null(colptr);
   memFree_null(rowind);
   memFree_null(val);
   */
  if (trscltb != NULL)
  {
    for (index=0; index<solvmtx->cblknbr; index++)
    {
      memFree_null(trscltb[index]);
    }
    memFree_null(trscltb);
  }

  /* 2nd membre */
  /* restore good coltab */
  colidx = 0;
  for (index=0; index<solvmtx->cblknbr; index++)
  {
    for(iter=0;iter<(CSC_COLNBR(thecsc,index)+1); iter++)
    {
      newcol = CSC_COL(thecsc,index,iter);
      CSC_COL(thecsc,index,iter) = colidx;
      colidx = newcol;
    }
  }
#if (DBG_SOPALIN_TIME==1)
  clockStop(&(clk));
  fprintf(stdout, "CscOrdistrib step 3 : %.3g s\n",
          (double)clockVal(&clk));
  clockInit(&clk);
  clockStart(&clk);
#endif
  CSC_SORT;

#if (DBG_SOPALIN_TIME==1)
  clockStop(&(clk));
  fprintf(stdout, "CscOrdistrib step 4 : %.3g s\n",
          (double)clockVal(&clk));
#endif
#ifdef CSC_LOG
  fprintf(stdout, "<- CscOrdistrib \n");
#endif
}

/******************************************************************************
 * Function: CscdOrdistrib                                                    *
 ******************************************************************************
 *                                                                            *
 * Fill in *thecsc* CSC matrix in column block representation.                *
 *                                                                            *
 * - Construct cachetab (sizeof(PASTIX_INT)*globalNbCol) which will contain        *
 *   the column block wich will own each column (internal numerotation),      *
 *   or -1 if not local                   *
 *                        *
 * - Build newcoltab (sizeof(PASTIX_INT)*globalNbCol) which will contain the         *
 *   coltab corresponding to the local internal CSCd.             *
 *   This CSCd correspond to the given CSCd adding upper part in        *
 *   Symmetric matrix.                    *
 *   Also count number of triples (i,j,v) to send to each other processors.   *
 *                        *
 * - Send the information about how many triples will be sent           *
 *                        *
 * - Fill-in the arrays containing triples to send and send them.         *
 *                        *
 * - Receive those arrays and correct the newcoltab arrays with information   *
 *   from others processors.                  *
 *                        *
 * - Build CSC_COLNBR from symbolic matrix informations and CSC_COL from      *
 *   newcoltab.                     *
 *                        *
 * - Construct transpose matrix, in symmetric mode, transcsc == CSC_VALTAB;   *
 *   in unsymmetric mode, allocate trowtab (number of total local elements),  *
 *   and build trscltb which contains number of elements,           *
 *   in each column of each column bloc.              *
 *                        *
 * - fill-in internal CSC row and values from local given CSCd,         *
 *   also fill-in trowtab and transcsc in unsymmetric mode.           *
 *   CSC_COL and trscltb are incremented for each element added.        *
 *                        *
 * - fill-in  internal CSC row and values from iniformation received,         *
 *   also fill in transposed CSCd in unsymmetric mode.            *
 *   CSC_COL and trscltb are incremented for each element added.        *
 *                        *
 * - restore CSC_COL.                     *
 *                        *
 * - sort internal CSCd.                  *
 *                        *
 * - sort intranal transposed CSCd.                 *
 *                        *
 * Parameters:                      *
 *                                                                            *
 *   thecsc     - Matrix in block column CSC format to fill in.         *
 *   Type       - 3 charactères for matrix Type : only Type[1] is used to     *
 *                check if matrix is Symetric(S) or not(U).           *
 *   transcsc   - Transpose of the CSC in non symetric mode.          *
 *   ord        - ordering                  *
 *   Ncol       - Number of columns.                *
 *   colptr     - Index in *rowind* and *val* of the start of each column.    *
 *   rowind     - Index of the elements.              *
 *   val        - values of the elements.               *
 *   l2g        - global numbers of local nodes.            *
 *   gNcol      - global number of columns.               *
 *   g2l        - local numbers of global nodes, if not local contains -owner *
 *   forcetrans - If matrix symetric, transcsc will be the copy of the        *
 *                CSC_VALTAB.                   *
 *   solvmtx    - Solver matrix                 *
 *   procnum    - MPI process number.                 *
 *   dof        - Number of degree of freedom.              *
 *   comm       - MPI communicator.                 *
 *                                                                            *
 ******************************************************************************/
void CscdOrdistrib(CscMatrix          *thecsc,
                   char               *Type,
                   PASTIX_FLOAT             **transcsc,
                   const Order        *ord,
                   PASTIX_INT                 Ncol,
                   PASTIX_INT                *colptr,
                   PASTIX_INT                *rowind,
                   PASTIX_FLOAT              *val,
                   PASTIX_INT                *l2g,
                   PASTIX_INT                 gNcol,
                   PASTIX_INT                *g2l,
                   PASTIX_INT                 forcetrans,
                   const SolverMatrix *solvmtx,
                   PASTIX_INT                 procnum,
                   PASTIX_INT                 dof,
                   MPI_Comm            comm)
{
  PASTIX_INT          index;
  PASTIX_INT          itercol;
  PASTIX_INT          newcol;
  PASTIX_INT          loccol;
  PASTIX_INT          iter;
  PASTIX_INT          rowp1;
  PASTIX_INT          colidx;
  PASTIX_INT          itercblk;
  PASTIX_INT          itercblk2;
  PASTIX_INT          strdcol     = 0;
  PASTIX_INT        **trscltb     = NULL;
  PASTIX_INT         *trowtab     = NULL;
  PASTIX_INT         *cachetab    = NULL;
  int          commSize;
  PASTIX_INT          proc;
  PASTIX_INT         *tosend      = NULL;
  MPI_Request *tosend_req  = NULL;
  MPI_Request *tosend_creq = NULL;
  MPI_Request *tosend_rreq = NULL;
  MPI_Request *tosend_vreq = NULL;
  PASTIX_INT         *torecv      = NULL;
  MPI_Request *torecv_req  = NULL;
  PASTIX_INT         *tosend_cnt  = NULL;
  PASTIX_INT        **tosend_col  = NULL;
  PASTIX_INT        **tosend_row  = NULL;
  PASTIX_FLOAT      **tosend_val  = NULL;
  PASTIX_INT        **torecv_col  = NULL;
  PASTIX_INT        **torecv_row  = NULL;
  PASTIX_FLOAT      **torecv_val  = NULL;
  PASTIX_INT         *newcoltab   = NULL;
  PASTIX_INT          owner;
  PASTIX_INT          therow;
#ifndef FORCE_NOMPI
  MPI_Status   status;
#endif
  PASTIX_INT          iterdofrow;
  PASTIX_INT          iterdofcol;
  PASTIX_INT          nodeidx;
  PASTIX_INT          colsize;
  (void)comm;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscdOrdistrib \n");
#endif
  thecsc->type = Type[1];
#if (DBG_SOPALIN_TIME==1)
  Clock clk;
  clockInit(&clk);
  clockStart(&clk);
#endif

  /* Initialize some MPI structures */
  CALL_MPI MPI_Comm_size(comm, &commSize);
  TEST_MPI("MPI_Comm_size");

  /* tosend will count the number of coefficient to send */
  MALLOC_INTERN(tosend, commSize, PASTIX_INT);
  for (proc = 0; proc < commSize; proc++)
    tosend[proc] = 0;
  MALLOC_INTERN(tosend_req,  commSize, MPI_Request);
  MALLOC_INTERN(tosend_creq, commSize, MPI_Request);
  MALLOC_INTERN(tosend_rreq, commSize, MPI_Request);
  MALLOC_INTERN(tosend_vreq, commSize, MPI_Request);

  /* cachetab: contain the column block or -1 if not local */
  MALLOC_INTERN(cachetab, (gNcol+1)*dof, PASTIX_INT);
  for (itercol=0; itercol< (gNcol+1)*dof; itercol++)
    cachetab[itercol] = -1;
  for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
  {
    for (itercol=solvmtx->cblktab[itercblk].fcolnum;
         itercol<solvmtx->cblktab[itercblk].lcolnum+1;
         itercol++)
    {
      cachetab[itercol] = itercblk;
    }
  }

  /* newcoltab will contain the number of element in each column
   of the symetrized matrix in the new ordering */
  MALLOC_INTERN(newcoltab, gNcol+1, PASTIX_INT);

  for(index=0; index<(gNcol+1); index++)
    newcoltab[index] = 0;

  for (itercol=0; itercol<Ncol; itercol++)
  {
    newcol = ord->permtab[l2g[itercol]-1];

    newcoltab[newcol] += colptr[itercol+1] - colptr[itercol];

    for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
    {

      loccol = g2l[rowind[iter-1]-1];
      if (COL_IS_LOCAL(loccol))
      {
        if (Type[1] == 'S' || Type[1] == 'H')
        {
          if (loccol -1 != itercol)
          {
            newcol = ord->permtab[rowind[iter-1]-1];
            newcoltab[newcol]++;
          }
        }
      }
      else
      {
        tosend[-loccol]++;
      }
    }
  }

  /* Will recv information about what will be received from other processors */
  MALLOC_INTERN(torecv, commSize, PASTIX_INT);
  for (proc = 0; proc < commSize; proc++)
    torecv[proc] = 0;

  MALLOC_INTERN(torecv_req, commSize, MPI_Request);

  for (proc = 0; proc < commSize; proc++)
  {
    if (proc != procnum)
    {
      CALL_MPI MPI_Isend(&tosend[proc], 1, COMM_INT, proc, 0,
                         comm, &tosend_req[proc]);
      TEST_MPI("Isend");
      CALL_MPI MPI_Irecv(&torecv[proc], 1, COMM_INT, proc, 0,
                         comm, &torecv_req[proc]);
      TEST_MPI("Irecv");
    }
  }

  /* Will contains values from other processors to add to the local
   internal CSCD */
  MALLOC_INTERN(tosend_col, commSize, PASTIX_INT*);
  MALLOC_INTERN(tosend_row, commSize, PASTIX_INT*);
  MALLOC_INTERN(tosend_val, commSize, PASTIX_FLOAT*);
  MALLOC_INTERN(tosend_cnt, commSize, PASTIX_INT);

  for (proc = 0; proc < commSize; proc++)
  {
    if (proc != procnum && tosend[proc] > 0)
    {
      MALLOC_INTERN(tosend_col[proc], tosend[proc], PASTIX_INT);
      MALLOC_INTERN(tosend_row[proc], tosend[proc], PASTIX_INT);
      MALLOC_INTERN(tosend_val[proc], tosend[proc]*dof*dof, PASTIX_FLOAT);
    }
    tosend_cnt[proc] = 0;
  }


  /* Filling in sending tabs with values and rows*/
  for (itercol=0; itercol<Ncol; itercol++)
  {
    itercblk = cachetab[(ord->permtab[l2g[itercol]-1])*dof];

    if (itercblk != -1)
    {
      /* ok put the value */
      for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
      {

        rowp1 = rowind[iter-1]-1;
        newcol = ord->permtab[l2g[itercol]-1];
        therow = ord->permtab[rowp1];

        itercblk2 = cachetab[therow*dof];

        if (itercblk2 != -1)
        {

        }
        else
        {
          /* Prepare to send row to the owner */
          owner = -g2l[ord->peritab[therow]];

          tosend_col[owner][tosend_cnt[owner]] = newcol;
          tosend_row[owner][tosend_cnt[owner]] = therow;

          memcpy(&(tosend_val[owner][tosend_cnt[owner]*dof*dof]),
                 &(val[(iter-1)*dof*dof]),
                 sizeof(PASTIX_FLOAT)*dof*dof);
          tosend_cnt[owner]++;
        }
      }

    }
    else
    {
      /* Impossible*/

    }
  }

  /* Sending values to other processors in IJV format. */
  for (proc = 0; proc < commSize; proc++)
  {
    if (proc != procnum)
    {
      CALL_MPI MPI_Wait(&tosend_req[proc], &status);
      TEST_MPI("MPI_Wait");
      if (tosend_cnt[proc] > 0)
      {
        CALL_MPI MPI_Isend(tosend_col[proc], (int)(tosend_cnt[proc]),
                           COMM_INT, proc, 1, comm, &tosend_creq[proc]);
        TEST_MPI("MPI_Isend");
        CALL_MPI MPI_Isend(tosend_row[proc], (int)(tosend_cnt[proc]),
                           COMM_INT, proc, 2, comm, &tosend_rreq[proc]);
        TEST_MPI("MPI_Isend");
        CALL_MPI MPI_Isend(tosend_val[proc],
                           (int)(tosend_cnt[proc]*dof*dof),
                           COMM_FLOAT, proc, 3, comm, &tosend_vreq[proc]);
        TEST_MPI("MPI_Isend");
      }
      CALL_MPI MPI_Wait(&torecv_req[proc], &status);
      TEST_MPI("MPI_Wait");
    }
  }

  memFree_null(tosend_req);
  memFree_null(torecv_req);

  /* Receiving information from other processors and updating newcoltab */
  MALLOC_INTERN(torecv_col, commSize, PASTIX_INT*);
  MALLOC_INTERN(torecv_row, commSize, PASTIX_INT*);
  MALLOC_INTERN(torecv_val, commSize, PASTIX_FLOAT*);
  for (proc = 0; proc < commSize; proc++)
  {
    if (proc != procnum && torecv[proc] > 0 )
    {
      MALLOC_INTERN(torecv_col[proc], torecv[proc], PASTIX_INT);
      MALLOC_INTERN(torecv_row[proc], torecv[proc], PASTIX_INT);
      MALLOC_INTERN(torecv_val[proc], torecv[proc]*dof*dof, PASTIX_FLOAT);
      CALL_MPI MPI_Recv(torecv_col[proc], torecv[proc], COMM_INT,
                        proc, 1, comm, &status );
      TEST_MPI("MPI_Recv");
      CALL_MPI MPI_Recv(torecv_row[proc], torecv[proc], COMM_INT,
                        proc, 2, comm, &status );
      TEST_MPI("MPI_Recv");
      CALL_MPI MPI_Recv(torecv_val[proc], torecv[proc]*dof*dof, COMM_FLOAT,
                        proc, 3, comm, &status );
      TEST_MPI("MPI_Recv");

      if (Type[1] == 'S' || Type[1] == 'H')
      {
        for (iter = 0; iter < torecv[proc]; iter++)
        {
          newcol= torecv_row[proc][iter];
          newcoltab[newcol] ++;
        }
      }
    }
  }


#if (DBG_SOPALIN_TIME==1)
  clockStop(&(clk));
  fprintf(stdout, "CscdOrdistrib step 1 : %.3g s\n",
          (double)clockVal(&clk));
  clockInit(&clk);
  clockStart(&clk);
#endif
  /* Finishing newcoltab construction :
   *
   * Now, newcoltab will contain starting index of each
   * column of rows and values in new ordering
   */
  newcol = 0;
  for (index=0; index<(gNcol+1); index++)
  {
    colidx = newcoltab[index];
    newcoltab[index] = newcol;
    newcol += colidx;
  }

  SET_CSC_COL(newcoltab);

#if (DBG_SOPALIN_TIME==1)
  clockStop(&(clk));
  fprintf(stdout, "CscdOrdistrib step 2 : %.3g s\n",
          (double)clockVal(&clk));
  clockInit(&clk);
  clockStart(&clk);
#endif

  CSC_ALLOC;

  /* Filling in thecsc with values and rows*/
  for (itercol=0; itercol<Ncol; itercol++)
  {
    itercblk = cachetab[(ord->permtab[l2g[itercol]-1])*dof];

    if (itercblk != -1)
    {
      /* ok put the value */
      for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
      {

        for (iterdofcol = 0; iterdofcol < dof; iterdofcol++)
        {
          for (iterdofrow = 0; iterdofrow < dof; iterdofrow++)
          {
            rowp1 = rowind[iter-1]-1;
            therow = ord->permtab[rowp1]*dof + iterdofrow;
            newcol = (ord->permtab[l2g[itercol]-1])*dof+iterdofcol;
            SET_CSC_ROW_VAL(itercblk, therow, newcol, val);

            itercblk2 = cachetab[therow];

            if (itercblk2 != -1)
            {
              switch (Type[1])
              {
              case 'S':
              {
                if (rowp1 != l2g[itercol]-1)
                {
                  /* same thing but newcol <-> therow */
                  SET_CSC_ROW_VAL(itercblk2, newcol, therow,
                                  val);

                }
                break;
              }
              case 'H':
              {
                if (rowp1 != l2g[itercol]-1)
                {
                  /* same thing but newcol <-> therow */
                  SET_CSC_ROW_VAL_CONJ(itercblk2, newcol, therow,
                                  val);

                }
                break;
              }
              case 'U':
              {
                if (transcsc != NULL)
                {
                  SET_TRANS_ROW_VAL(itercblk2, therow, newcol,
                                    val);
                }
                break;
              }
              default:
                errorPrint("Unknown Matrix Type");
                EXIT(MOD_SOPALIN, UNKNOWN_ERR);
              }

            }
          }
        }
      }

    }
    else
    {
      /* Impossible*/
      errorPrint("Error in CscdOrdistrib");
      EXIT(MOD_SOPALIN, UNKNOWN_ERR)
        }
  }

  memFree_null(tosend);

  for (proc = 0; proc < commSize; proc++)
  {
    if (proc != procnum)
    {
      for (iter = 0; iter < torecv[proc]; iter++)
      {
        for (iterdofcol = 0; iterdofcol < dof; iterdofcol++)
        {
          for (iterdofrow = 0; iterdofrow < dof; iterdofrow++)
          {
            switch (Type[1])
            {
            case 'S':
            {
              newcol  = torecv_col[proc][iter]*dof+iterdofcol;
              therow  = torecv_row[proc][iter]*dof+iterdofrow;
              itercblk2 = cachetab[therow];
              /* iter is 0 based here, not in SET_CSC_ROW_VAL */
              iter++;
              SET_CSC_ROW_VAL(itercblk2, newcol, therow,
                              torecv_val[proc]);
              iter--;
              break;
            }
            case 'H':
            {
              newcol  = torecv_col[proc][iter]*dof+iterdofcol;
              therow  = torecv_row[proc][iter]*dof+iterdofrow;
              itercblk2 = cachetab[therow];
              /* iter is 0 based here, not in SET_CSC_ROW_VAL */
              iter++;
              SET_CSC_ROW_VAL_CONJ(itercblk2, newcol, therow,
                                   torecv_val[proc]);
              iter--;
              break;
            }
            case 'U':
            {
              if (transcsc != NULL)
              {
                newcol = torecv_col[proc][iter]*dof+iterdofcol;
                therow = torecv_row[proc][iter]*dof+iterdofrow;
                itercblk2 = cachetab[therow];
                /* iter is 0 based here,
                 not in SET_TRANS_ROW_VAL */
                iter++;
                SET_TRANS_ROW_VAL(itercblk2, therow, newcol,
                                  torecv_val[proc]);
                iter--;
              }
              break;
            }
            default:
              errorPrint("Unknown Matrix Type");
              EXIT(MOD_SOPALIN, UNKNOWN_ERR);
            }
          }
        }
      }
      if (torecv[proc] > 0)
      {
        memFree_null(torecv_col[proc]);
        memFree_null(torecv_row[proc]);
        memFree_null(torecv_val[proc]);
      }
    }
  }
  memFree_null(torecv_col);
  memFree_null(torecv_row);
  memFree_null(torecv_val);
  memFree_null(cachetab);
  memFree_null(torecv);
  for (proc = 0; proc < commSize; proc++)
  {
    if (proc != procnum)
    {
      if (tosend_cnt[proc] > 0)
      {
        CALL_MPI MPI_Wait(&tosend_creq[proc], &status);
        TEST_MPI("MPI_Wait");
        memFree_null(tosend_col[proc]);
        CALL_MPI MPI_Wait(&tosend_rreq[proc], &status);
        TEST_MPI("MPI_Wait");
        memFree_null(tosend_row[proc]);
        CALL_MPI MPI_Wait(&tosend_vreq[proc], &status);
        TEST_MPI("MPI_Wait");
        memFree_null(tosend_val[proc]);
      }
    }
  }
  memFree_null(tosend_creq);
  memFree_null(tosend_rreq);
  memFree_null(tosend_vreq);
  memFree_null(tosend_cnt);
  memFree_null(tosend_col);
  memFree_null(tosend_row);
  memFree_null(tosend_val);

  if (trscltb != NULL)
  {
    for (index=0; index<solvmtx->cblknbr; index++)
    {
      memFree_null(trscltb[index]);
    }
    memFree_null(trscltb);
  }

  /* 2nd membre */
  /* restore good coltab */
  colidx = 0;
  for (index=0; index<solvmtx->cblknbr; index++)
  {
    for(iter=0;iter<(CSC_COLNBR(thecsc,index)+1); iter++)
    {
      newcol = CSC_COL(thecsc,index,iter);
      CSC_COL(thecsc,index,iter) = colidx;
      colidx = newcol;
    }
  }

#if (DBG_SOPALIN_TIME==1)
  clockStop(&(clk));
  fprintf(stdout, "CscdOrdistrib step 3 : %.3g s\n",
          (double)clockVal(&clk));
  clockInit(&clk);
  clockStart(&clk);
#endif

  CSC_SORT;

#if (DBG_SOPALIN_TIME==1)
  clockStop(&(clk));
  fprintf(stdout, "CscdOrdistrib step 4 : %.3g s\n",
          (double)clockVal(&clk));
#endif
#ifdef CSC_LOG
  fprintf(stdout, "<- CscdOrdistrib \n");
#endif
}


/******************************************************************************
 * Function: CscExit                                                          *
 ******************************************************************************
 *                                                                            *
 * Free the internal CSCd structure.                                          *
 *                                                                            *
 * Parameters:                                                                *
 *   thecsc - Internal CSCd to free.                                          *
 *                                                                            *
 ******************************************************************************/
void CscExit(CscMatrix *thecsc)
{
  PASTIX_INT itercscf;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscExit \n");
#endif

  if ( CSC_FTAB(thecsc) != NULL )
  {
    for (itercscf = 0;
         itercscf < CSC_FNBR(thecsc);
         itercscf++)
    {
      memFree_null(CSC_COLTAB(thecsc,itercscf));
    }

    memFree_null(CSC_FTAB(thecsc));

    /* Destruction coltab et valtab */
    memFree_null(CSC_ROWTAB(thecsc));
    memFree_null(CSC_VALTAB(thecsc));

    CSC_FTAB(thecsc)   = NULL;
    CSC_ROWTAB(thecsc) = NULL;
    CSC_VALTAB(thecsc) = NULL;
  }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscExit \n");
#endif
}
