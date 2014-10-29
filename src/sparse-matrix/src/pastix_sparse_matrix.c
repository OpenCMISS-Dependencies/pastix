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
#include "common_pastix.h"
#include "symbol.h"
#include "ftgt.h"
#include "csc.h"
#include "queue.h"
#include "bulles.h"
#include "updown.h"
#include "solver.h"
#include "sopalin_define.h"
#include "sopalin_thread.h"
#include "order.h"
#ifdef FORCE_NOMPI
#  include "nompi.h"
#else
#  include <mpi.h>
#endif
#include "sopalin3d.h"
#ifdef WITH_SCOTCH
#  ifdef    DISTRIBUTED
#    include <ptscotch.h>
#  else
#    include <scotch.h>
#  endif /* DISTRIBUTED */
#endif /* WITH_SCOTCH */
#ifdef WITH_SEM_BARRIER
#  include <semaphore.h>
#  include <fcntl.h>
#  include <sys/ipc.h>
#  include <sys/shm.h>
#endif
#include "pastixstr.h"
#include "pastix.h"
#include "pastix_sparse_matrix.h"
#include "pastix_sparse_matrix_intern.h"



int
symbol_to_sparse_matrix(SymbolMatrix             symbmtx,
                        pastix_sparse_matrix_t * sparsemtx,
                        pastix_coef_type_t       coeftype,
                        SolverCblk             * cblktab)
{
  PASTIX_INT cblknum, bloknum, cblksize;
  pastix_block_t * block;
  pastix_cblk_t  * cblk;
  sparsemtx->baseval  = symbmtx.baseval;
  sparsemtx->coeftype = coeftype;
  sparsemtx->coefnbr  = 0;
  sparsemtx->nodenbr  = symbmtx.nodenbr;
  sparsemtx->cblknbr  = symbmtx.cblknbr;
  sparsemtx->bloknbr  = symbmtx.bloknbr;
  sparsemtx->maxcblksize = 0;
  MALLOC_INTERN(sparsemtx->cblktab, sparsemtx->cblknbr, struct pastix_cblk_);
  MALLOC_INTERN(sparsemtx->bloktab, sparsemtx->bloknbr, struct pastix_block_);
  cblk = sparsemtx->cblktab;
  for (cblknum = 0; cblknum < sparsemtx->cblknbr; cblknum++, cblk++)
    {
      cblk->fblknbr = 0;
      cblk->fblktab = NULL;
      cblk->lvalues = NULL;
      cblk->uvalues = NULL;
      cblk->lschur  = NULL;
      cblk->uschur  = NULL;
    }

  cblk = sparsemtx->cblktab;
  for (cblknum = 0; cblknum < sparsemtx->cblknbr; cblknum++, cblk++)
    {
      cblk->fcolnum  = symbmtx.cblktab[cblknum].fcolnum;
      cblk->lcolnum  = symbmtx.cblktab[cblknum].lcolnum;
      cblk->bloknbr  =
        symbmtx.cblktab[cblknum+1].bloknum - symbmtx.cblktab[cblknum].bloknum;
      cblk->bloktab  = &(sparsemtx->bloktab[symbmtx.cblktab[cblknum].bloknum]);
      cblk->color    = cblktab[cblknum].color;
      cblk->procdiag = cblktab[cblknum].procdiag;
      cblk->cblkdiag = cblktab[cblknum].cblkdiag;

      cblk->stride = 0;

      block = cblk->bloktab;
      for (bloknum = symbmtx.cblktab[cblknum].bloknum;
           bloknum < symbmtx.cblktab[cblknum+1].bloknum;
           bloknum++, block++)
        {
          block->frownum = symbmtx.bloktab[bloknum].frownum;
          block->lrownum = symbmtx.bloktab[bloknum].lrownum;
          block->levfval = symbmtx.bloktab[bloknum].levfval;
          block->coefind = cblk->stride;
          block->cblk    = cblk;
          block->fcblk   =
            &(sparsemtx->cblktab[symbmtx.bloktab[bloknum].cblknum]);
          block->fcblk->fblknbr++;
          cblk->stride  += block->lrownum - block->frownum + 1;
        }
      cblksize      = (cblk->lcolnum - cblk->fcolnum + 1)*(cblk->stride);
      sparsemtx->maxcblksize = MAX(sparsemtx->maxcblksize, cblksize);
      sparsemtx->coefnbr += cblksize;
    }

  /* Allocate fblknbr */
  cblk = sparsemtx->cblktab;
  for (cblknum = 0; cblknum < sparsemtx->cblknbr; cblknum++, cblk++)
    {
      MALLOC_INTERN(cblk->fblktab,cblk->fblknbr, pastix_block_t*);
      cblk->fblknbr = 0;
    }

  /* Fill fblknbr */
  cblk = sparsemtx->cblktab;
  for (cblknum = 0; cblknum < sparsemtx->cblknbr; cblknum++, cblk++)
    {
      block = cblk->bloktab;
      for (bloknum = symbmtx.cblktab[cblknum].bloknum;
           bloknum < symbmtx.cblktab[cblknum+1].bloknum;
           bloknum++, block++)
        {
          pastix_cblk_t * fcblk =
            &(sparsemtx->cblktab[symbmtx.bloktab[bloknum].cblknum]);
          fcblk->fblktab[fcblk->fblknbr] = block;
          fcblk->fblknbr++;
        }
    }
  return NO_ERR;
}

int
pastix_get_sparse_matrix(pastix_data_t          * pastix_data,
                         pastix_coef_type_t       coeftype,
                         pastix_sparse_matrix_t * sparsemtx)
{
  return symbol_to_sparse_matrix(*(pastix_data->symbmtx),
                                 sparsemtx,
                                 coeftype,
                                 pastix_data->solvmatr.cblktab);
}


int
pastix_sparse_matrix_free(pastix_sparse_matrix_t * sparsemtx)
{
  PASTIX_INT cblknum;
  for (cblknum = 0; cblknum < sparsemtx->cblknbr; cblknum++)
    {
      if (NULL != sparsemtx->cblktab[cblknum].lvalues)
        memFree_null(sparsemtx->cblktab[cblknum].lvalues);
      if (NULL != sparsemtx->cblktab[cblknum].uvalues)
        memFree_null(sparsemtx->cblktab[cblknum].uvalues);
      if (NULL != sparsemtx->cblktab[cblknum].lschur)
        memFree_null(sparsemtx->cblktab[cblknum].lschur);
      if (NULL != sparsemtx->cblktab[cblknum].uschur)
        memFree_null(sparsemtx->cblktab[cblknum].uschur);
    }
  memFree_null(sparsemtx->cblktab);
  memFree_null(sparsemtx->bloktab);
  return NO_ERR;
}
