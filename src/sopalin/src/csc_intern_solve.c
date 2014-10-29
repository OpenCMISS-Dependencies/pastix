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
#include <stdlib.h>

#include "common_pastix.h"
#include "csc.h"
#include "symbol.h"
#include "ftgt.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"

#include "sopalin_acces.h"
#include "csc_intern_solve.h"

#ifdef DEBUG_RAFF
#  define CSC_LOG
#endif

/*
 * Function: Csc2solv_cblk
 *
 * Copy the part of the internal CSCd corresponding to
 * the column bloc itercblk into the SolverMatrix structure
 * coeftab which will be used to compute the decomposition.
 *
 * Used in NUMA mode.
 *
 * Parameters:
 *   cscmtx   - The internal CSCd matrix.
 *   datacode - The SolverMatrix structure used during decomposition.
 *   trandcsc - The internal CSCd transpose used in LU decomposition.
 *   itercblk - Column bloc number in which we had the internal CSCd.
 */
void Csc2solv_cblk(const CscMatrix *cscmtx,
                   SolverMatrix    *datacode,
                   PASTIX_FLOAT           *trandcsc,
                   PASTIX_INT              itercblk)
{
  PASTIX_INT itercoltab;
  PASTIX_INT iterbloc;
  PASTIX_INT coefindx;
  PASTIX_INT iterval;

#ifdef CSC_LOG
  fprintf(stdout, "-> Csc2solv \n");
#endif

  if (itercblk < CSC_FNBR(cscmtx)){
    for (itercoltab=0;
         itercoltab < CSC_COLNBR(cscmtx,itercblk);
         itercoltab++)
      {
        for (iterval = CSC_COL(cscmtx,itercblk,itercoltab);
             iterval < CSC_COL(cscmtx,itercblk,itercoltab+1);
             iterval++)
          {
            if (CSC_ROW(cscmtx,iterval) >=
                SYMB_FCOLNUM(itercblk))
              {
                iterbloc = SYMB_BLOKNUM(itercblk);

                ASSERTDBG(iterbloc < SYMB_BLOKNBR, MOD_SOPALIN);
                while (( iterbloc < SYMB_BLOKNUM(itercblk+1)) &&
                       (( SYMB_LROWNUM(iterbloc) < CSC_ROW(cscmtx,iterval)) ||
                        ( SYMB_FROWNUM(iterbloc) > CSC_ROW(cscmtx,iterval))))
                  {
                    iterbloc++;
                  }

                if ( iterbloc < SYMB_BLOKNUM(itercblk+1) )
                  {
                    coefindx = SOLV_COEFIND(iterbloc);

                    coefindx += CSC_ROW(cscmtx,iterval) - SYMB_FROWNUM(iterbloc);

                    coefindx += SOLV_STRIDE(itercblk)*itercoltab;
                    SOLV_COEFTAB(itercblk)[coefindx] = CSC_VAL(cscmtx,iterval);

                    if (trandcsc != NULL && iterbloc != SYMB_BLOKNUM(itercblk))
                      {
                        if (cscmtx->type == 'H')
                          SOLV_UCOEFTAB(itercblk)[coefindx] = CONJ_FLOAT(trandcsc[iterval]);
                        else
                          SOLV_UCOEFTAB(itercblk)[coefindx] = trandcsc[iterval];
                      }
                  }
                else printf("ILU: csc2solv drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n",
                            (long)datacode->cblktab[itercblk].fcolnum+
                            (long)itercoltab,(long)itercoltab,
                            (long)CSC_ROW(cscmtx,iterval),(long)iterval,
                            (long)itercblk,
                            (long)datacode->cblktab[itercblk].fcolnum,
                            (long)datacode->cblktab[itercblk].lcolnum);
              }
          }
      }
  }
#ifdef CSC_LOG
  fprintf(stdout, "<- Csc2solv \n");
#endif
}
