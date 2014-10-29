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
#include "dof.h"
#include "cost.h"
#include "symbol.h"
#include "elimin.h"
#include "extrastruct.h"
#include "extendVector.h"
#include "cand.h"
#include "param_blend.h"
#include "queue.h"
#include "bulles.h"
#include "blendctrl.h"
#include "ftgt.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "simu.h"
#include "costfunc.h"
#include "partbuild.h"
/* #include "eliminfunc.h" */
/* #include "splitpart.h" */

struct PASTIX_interval_s {
  PASTIX_INT start;
  PASTIX_INT end;
  PASTIX_INT nFacingBlocks;
};
typedef struct PASTIX_interval_s PASTIX_interval_t;

struct PASTIX_partition_s {
  PASTIX_interval_t * part;
  PASTIX_INT          size;
  PASTIX_INT          max_size;
};
typedef struct PASTIX_partition_s PASTIX_partition_t;

static
int cmp_interval(const void *p1, const void *p2)
{
  PASTIX_interval_t * i1 = (PASTIX_interval_t *)p1;
  PASTIX_interval_t * i2 = (PASTIX_interval_t *)p2;

  if (i1->start != i2->start)
    return i1->start - i2->start;

  return i1->end - i2->end;
}

static inline
int add_interval(PASTIX_partition_t * p,
                 PASTIX_INT s,
                 PASTIX_INT e,
                 PASTIX_INT n)
{
  ASSERTDBG(e  >= s, MOD_BLEND);
  p->part[p->size].start = s;
  p->part[p->size].end   = e;
  p->part[p->size].nFacingBlocks = n;
  p->size++;
  return 0;
}

static inline
int init_partition(PASTIX_INT size,
                   PASTIX_INT start,
                   PASTIX_INT end,
                   PASTIX_partition_t * p)
{
  p->size = 0;
  p->max_size = size;
  MALLOC_INTERN(p->part, size, PASTIX_interval_t);
  add_interval(p, start, end, 0);
  return 0;
}

static inline
int deinit_partition(PASTIX_partition_t * p)
{
  memFree_null(p->part);
  return 0;
}
#define ADD_TMP_INTERVAL(s, e, n) do {          \
    tmp--;                                      \
    ASSERTDBG(p->size-1 != tmp, MOD_BLEND);     \
    p->part[tmp].start = s;                     \
    p->part[tmp].end   = e;                     \
    p->part[tmp].nFacingBlocks = n;             \
  } while(0)

#define REMOVE_TMP_INTERVAL(j) do {             \
    if (j != tmp)                               \
      {                                         \
        p->part[j] = p->part[tmp];              \
      }                                         \
    tmp++;                                      \
  } while (0)

static inline
int insert_in_partition(PASTIX_INT start,
                        PASTIX_INT end,
                        PASTIX_partition_t *p)
{
  PASTIX_INT i,j, tmp;

  tmp = p->max_size-1;
  p->part[tmp].start = start;
  p->part[tmp].end   = end;

  for (i = 0; i < p->size; i++)
    {
      for (j = tmp;  j < p->max_size; j++)
        {
          PASTIX_interval_t * my_inter;
          my_inter = &(p->part[j]);

          if (p->part[i].end <  my_inter->start)
            /* the new intervale is after this one and does not intersect */
            continue;

          if (p->part[i].start > my_inter->end)
            /* the new intervale is before this one and does not intersect */
            continue;

          /* The new interval intersect this one */
          if (p->part[i].start > my_inter->start)
            {
              if (p->part[i].end > my_inter->end)
                {
                  /* my_int  i
                   * ----          ----
                   *  A    ----    ----
                   * ----       => ----
                   *         B
                   *       ----    ----
                   */
                  PASTIX_INT k = p->part[i].end;
                  /* split B in A inter B and B \ A */
                  /* B \ A */
                  add_interval(p, my_inter->end+1, k, p->part[i].nFacingBlocks);
                  ASSERTDBG(p->size-1 != tmp, MOD_BLEND);

                  /* A inter B */
                  p->part[i].end  = my_inter->end;
                  p->part[i].nFacingBlocks ++;

                  /* Now we search to intersect with A \ B */
                  my_inter->end   = p->part[i].start-1;
                }
              else
                {
                  /* p->part[i].end <= my_inter->end */
                  if (p->part[i].end == my_inter->end)
                    {
                      /* my_int  i
                       * ----          ----
                       *  A    ----    ----
                       *            =>
                       *         B
                       * ----  ----    ----
                       */
                      p->part[i].nFacingBlocks++;
                      /* Now we search to intersect with A \ B */
                      my_inter->end   = p->part[i].start-1;
                    }
                  else
                    {
                      /* p->part[i].end < my_inter->end */
                      /* my_int  i
                       * ----          ----
                       *  A    ----    ----
                       *         B  =>
                       *       ----    ----
                       * ----          ----
                       */
                      p->part[i].nFacingBlocks++;
                      /* A \ B makes 2 interrvals ! */
                      /* store first one in a new area */
                      ADD_TMP_INTERVAL(p->part[i].end+1, my_inter->end, 1);

                      /* second overwrite my_inter */
                      my_inter->end = p->part[i].start-1;
                    }
                }
            }
          else
            {
              /* p->part[i].start <= my_inter->start */
              if (p->part[i].start == my_inter->start)
                {
                  if (p->part[i].end > my_inter->end)
                    {
                      /* my_int  i
                       * ----  ----    ----
                       *  A
                       * ----       => ----
                       *         B
                       *       ----    ----
                       */

                      /* B \ A */
                      add_interval(p, my_inter->end+1, p->part[i].end,
                                   p->part[i].nFacingBlocks);
                      ASSERTDBG(p->size-1 != tmp, MOD_BLEND);

                      /* A inter B */
                      p->part[i].end = my_inter->end;
                      p->part[i].nFacingBlocks++;
                      REMOVE_TMP_INTERVAL(j);
                    }
                  else
                    {
                      /* p->part[i].end <= my_inter->end */
                      if ( p->part[i].end == my_inter->end)
                        {
                          /* my_int  i
                           * ----  ----    ----
                           *  A
                           *            =>
                           *         B
                           * ----  ----    ----
                           */

                          p->part[i].nFacingBlocks++;
                          /* we are done with my_inter */
                          REMOVE_TMP_INTERVAL(j);
                        }
                      else
                        {
                          /* p->part[i].end < my_inter->end */
                          /* my_int  i
                           * ----  ----    ----
                           *  A
                           *         B   =>
                           *       ----    ----
                           * ----          ----
                           */
                          /* B inter A */
                          p->part[i].nFacingBlocks++;
                          /* A \ B */
                          my_inter->start = p->part[i].end+1;
                        }
                    }
                }
              else
                {
                  /* p->part[i].start < my_inter->start */
                  if (p->part[i].end > my_inter->end)
                    {
                      /* my_int  i
                       *       ----    ----
                       * ----          ----
                       *  A
                       * ----       => ----
                       *         B
                       *       ----    ----
                       */

                      /* B \ A */
                      add_interval(p, p->part[i].start, my_inter->start-1,
                                   p->part[i].nFacingBlocks);
                      ASSERTDBG(p->size-1 != tmp, MOD_BLEND);

                      add_interval(p, my_inter->end+1, p->part[i].end,
                                   p->part[i].nFacingBlocks);
                      ASSERTDBG(p->size-1 != tmp, MOD_BLEND);

                      /* A inter B */
                      p->part[i].start = my_inter->start;
                      p->part[i].end   = my_inter->end;
                      p->part[i].nFacingBlocks++;
                      /* we are done with my_inter */
                      REMOVE_TMP_INTERVAL(j);
                    }
                  else
                    {
                      /* p->part[i].end <= my_inter->end */
                      if ( p->part[i].end == my_inter->end)
                        {
                          /* my_int  i
                           *       ----    ----
                           * ----          ----
                           *  A
                           *            =>
                           *         B
                           * ----  ----    ----
                           */
                          add_interval(p, p->part[i].start, my_inter->start-1,
                                       p->part[i].nFacingBlocks);
                          ASSERTDBG(p->size-1 != tmp, MOD_BLEND);

                          p->part[i].start = my_inter->start;
                          p->part[i].nFacingBlocks++;
                          /* we are done with my_inter */
                          REMOVE_TMP_INTERVAL(j);
                        }
                      else
                        {
                          /* p->part[i].end < my_inter->end */
                          /* my_int  i
                           *       ----    ----
                           * ----          ----
                           *  A
                           *         B   =>
                           *       ----    ----
                           * ----          ----
                           */
                          /* B \ A */
                          add_interval(p, p->part[i].start, my_inter->start-1,
                                       p->part[i].nFacingBlocks);
                          ASSERTDBG(p->size-1 != tmp, MOD_BLEND);

                          /* B inter A */
                          p->part[i].start = my_inter->start;
                          p->part[i].nFacingBlocks++;
                          /* A \ B */
                          my_inter->start = p->part[i].end+1;
                        }
                    }
                }
            }
          ASSERTDBG(p->part[i].end  >= p->part[i].start, MOD_BLEND);
        }
    }

  /* Those ones didn't match any other interval */
  for (i = tmp; i < p->max_size; i++)
    add_interval(p, p->part[i].start, p->part[i].end, 1);
  return NO_ERR;
}

int smart_cblk_split(BlendCtrl      * ctrl,
                     SymbolMatrix   * symbmtx,
                     PASTIX_INT       cblknum,
                     PASTIX_INT       procnbr,
                     PASTIX_INT       blas_min_col,
                     PASTIX_INT       blas_max_col,
                     PASTIX_INT     * nseq,
                     PASTIX_INT    ** seq)
{
  PASTIX_INT   i;
  PASTIX_partition_t part;
  PASTIX_INT   block_size;
  PASTIX_INT   block_size_min;
  PASTIX_INT   block_size_max;

  /* in the worst case we have 2*blocknbr+1 intervals */
  init_partition(2*ctrl->egraph->verttab[cblknum].innbr+1,
                 symbmtx->cblktab[cblknum].fcolnum,
                 symbmtx->cblktab[cblknum].lcolnum,
                 &part);

  for(i=0;i<ctrl->egraph->verttab[cblknum].innbr;i++)
    {
      PASTIX_INT bloknum;
      bloknum = ctrl->egraph->inbltab[ctrl->egraph->verttab[cblknum].innum+i];

      ASSERTDBG(symbmtx->bloktab[bloknum].cblknum == cblknum, MOD_BLEND);
      insert_in_partition(symbmtx->bloktab[bloknum].frownum,
                          symbmtx->bloktab[bloknum].lrownum, &part);
    }

  qsort(part.part,
        part.size,
        sizeof(PASTIX_interval_t),
        cmp_interval);

#ifdef FLAG_ASSERT
  for (i = 0; i < part.size; i++)
    {
      ASSERTDBG(part.part[i].start <= part.part[i].end, MOD_BLEND);
      ASSERTDBG(part.part[i].nFacingBlocks >=0 , MOD_BLEND);
      if (i < part.size - 1)
        ASSERTDBG(part.part[i].end < part.part[i+1].start, MOD_BLEND);
    }
#endif /* FLAG_ASSERT */

  if (procnbr == 1)
    {
      /* Make blocks of size blas_max_col */
      block_size = blas_max_col;
    }
  else
    {
      PASTIX_INT abs = ctrl->option->abs;
      if(procnbr > ctrl->option->ratiolimit)
        {
          abs *= 2; /* Increase abs for 2D */
        }

      /***  If option adaptative block size is set then compute the size
       ***  of a column block ***/
      if(abs > 0)
        {
          block_size = (symbmtx->cblktab[cblknum].lcolnum -
                        symbmtx->cblktab[cblknum].fcolnum + 1)/(abs * procnbr);

          block_size = MAX(block_size, blas_min_col);
          block_size = MIN(block_size, blas_max_col);

          }
        else
          block_size = blas_min_col;
    }

  block_size_max = 2*block_size;
  if (block_size_max == 0) block_size_max = 2;
  block_size_min = 0.5*block_size;
  if (block_size_min == 0) block_size_min++;

  {
    PASTIX_INT start, end;
    start = part.part[0].start;
    /* maximum number of entries */
    (*nseq) = (symbmtx->cblktab[cblknum].lcolnum -
               symbmtx->cblktab[cblknum].fcolnum + 1)/block_size_min;
    /* extend intvec to contain seq */
    extendint_ToSize(2*(*nseq), ctrl->intvec);
    *seq = ctrl->intvec->inttab;


    *nseq = 0;
    for (i = 0; i < part.size; i++)
      {
        end = part.part[i].end;


        if (part.part[i].nFacingBlocks == 0)
          {
            /* we cut all the interval with size ~block_size */
            PASTIX_INT blocknbr;
            PASTIX_INT average_size;
            blocknbr = (end-start)/block_size;
            average_size = (end-start)/blocknbr;

            if (start + average_size < part.part[i].start)
              {
                /* we don't won't tp cut in previous block,
                   else we would have done it earlier*/
                (*seq)[2*(*nseq)]   = start;
                (*seq)[2*(*nseq)+1] = part.part[i-1].end;
                (*nseq)++;
                start = part.part[i-1].end+1;
                blocknbr--;
                average_size = (end-start)/blocknbr;
              }

            while(end - start > average_size)
              {
                (*seq)[2*(*nseq)]   = start;
                (*seq)[2*(*nseq)+1] = start + average_size - 1;
                (*nseq)++;
                start = start + average_size;
              }

            if (end - start > block_size_min)
              {
                (*seq)[2*(*nseq)]   = start;

                (*seq)[2*(*nseq)+1] = end;
                (*nseq)++;
                start = end + 1;
              }

            if (end - start > 0)
              {
                (*seq)[2*(*nseq-1)+1] = end;
              }
          }
        else
          {
            /* If we have to cut and the number of blocks is increasing,
               then best is to cut at the end of one block first */
            if(start !=  part.part[i].start &&
               part.part[i].nFacingBlocks >= part.part[i-1].nFacingBlocks &&
               part.part[i].end + 1- start > block_size_max &&
               part.part[i].start - start > block_size_min)
              {
                (*seq)[2*(*nseq)]   = start;
                (*seq)[2*(*nseq)+1] = part.part[i].start - 1;
                (*nseq)++;
                start = part.part[i].start;
              }

            /* We cut the less we can : block_size_max */
            while(end - start > block_size_max)
              {
                (*seq)[2*(*nseq)]   = start;
                (*seq)[2*(*nseq)+1] = start + block_size_max - 1;
                (*nseq)++;
                start = start + block_size_max;
              }

          }
      }
  }

  if (symbmtx->cblktab[cblknum].lcolnum != (*seq)[2*(*nseq-1)+1])
    {
      (*seq)[2*(*nseq)]   = (*seq)[2*(*nseq-1)+1]+1;
      (*seq)[2*(*nseq)+1] = symbmtx->cblktab[cblknum].lcolnum;
      (*nseq)++;
    }
  deinit_partition(&part);

#ifdef FLAG_ASSERT
  for (i = 0; i < *nseq; i++)
    {
      ASSERTDBG((*seq)[2*i] < (*seq)[2*i+1], MOD_BLEND);
      if (i < *nseq -1)
        ASSERTDBG((*seq)[2*i+1] + 1 == (*seq)[2*(i+1)], MOD_BLEND);
    }
#endif /* FLAG_ASSERT */
  ASSERTDBG(symbmtx->cblktab[cblknum].fcolnum == (*seq)[0], MOD_BLEND);
  ASSERTDBG(symbmtx->cblktab[cblknum].lcolnum == (*seq)[2*(*nseq-1)+1],
            MOD_BLEND);
  return 0;
}
