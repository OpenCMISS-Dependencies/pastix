/* Copyright INRIA 2004
**
** This file is part of the Scotch distribution.
**
** The Scotch distribution is libre/free software; you can
** redistribute it and/or modify it under the terms of the
** GNU Lesser General Public License as published by the
** Free Software Foundation; either version 2.1 of the
** License, or (at your option) any later version.
**
** The Scotch distribution is distributed in the hope that
** it will be useful, but WITHOUT ANY WARRANTY; without even
** the implied warranty of MERCHANTABILITY or FITNESS FOR A
** PARTICULAR PURPOSE. See the GNU Lesser General Public
** License for more details.
**
** You should have received a copy of the GNU Lesser General
** Public License along with the Scotch distribution; if not,
** write to the Free Software Foundation, Inc.,
** 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
**
** $Id: symbol_keep.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_keep.c                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module contains the incomplete     **/
/**                symbol matrix carving routines.         **/
/**                                                        **/
/**   DATES      : # Version 1.3  : from : 30 jun 2003     **/
/**                                 to     16 jul 2003     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL_KEEP

#include "common_pastix.h"
#include "symbol.h"

/***************************/
/*                         */
/* Block keeping routines. */
/*                         */
/***************************/

/* This routine initializes the given
** graph structure.
** It returns:
** - 0  : in all cases.
*/

int
symbolKeepInit (
SymbolKeep * restrict const keepptr,              /*+ Keep structure to initialize  +*/
const SymbolMatrix * const  symbptr)              /*+ Symbol structure to work from +*/
{
  PASTIX_INT                 levfmax;                    /* Maximum level of fill     */
  PASTIX_INT                 nupdmax;                    /* Maximum number of updates */
  PASTIX_INT                 ctrimax;
  PASTIX_INT                 ctromax;
  PASTIX_INT                 hghtmax;
  PASTIX_INT                 bloknum;

  memSet (keepptr, 0, sizeof (SymbolKeep));

  if (memAllocGroup ((void **)
                     &keepptr->keeptab, (size_t) (symbptr->bloknbr * sizeof (unsigned char)),
                     &keepptr->kblktab, (size_t) (symbptr->bloknbr * sizeof (SymbolKeepBlok)), NULL) == NULL) {
    errorPrint ("symbolKeepInit: out of memory (1)");
    return     (1);
  }
  memSet (keepptr->keeptab, 1, symbptr->bloknbr); /* All blocks are kept at first */

  symbolKeepCompute (keepptr, symbptr);           /* Compute parameters of full factorization */

  levfmax = 0;                                    /* Extract maximums */
  nupdmax = 0;
  ctrimax = 0;
  ctromax = 0;
  hghtmax = 0;
  for (bloknum = 0; bloknum < symbptr->bloknbr; bloknum ++) {
    if (keepptr->kblktab[bloknum].levfval > levfmax)
      levfmax = keepptr->kblktab[bloknum].levfval;
    if (keepptr->kblktab[bloknum].nupdval > nupdmax)
      nupdmax = keepptr->kblktab[bloknum].nupdval;
    if (keepptr->kblktab[bloknum].ctrival > ctrimax)
      ctrimax = keepptr->kblktab[bloknum].ctrival;
    if (keepptr->kblktab[bloknum].ctroval > ctromax)
      ctromax = keepptr->kblktab[bloknum].ctroval;
    if (keepptr->kblktab[bloknum].hghtval > hghtmax)
      hghtmax = keepptr->kblktab[bloknum].hghtval;
  }
  keepptr->levfmax = levfmax;
  keepptr->nupdmax = nupdmax;
  keepptr->ctrimax = ctrimax;
  keepptr->ctromax = ctromax;
  keepptr->hghtmax = hghtmax;

  if (memAllocGroup ((void **)
                     &keepptr->levftab, (size_t) ((levfmax + 1)    * sizeof (double)),
                     &keepptr->nupdtab, (size_t) ((nupdmax + 1)    * sizeof (double)),
                     &keepptr->ctritab, (size_t) ((ctrimax + 1)    * sizeof (double)),
                     &keepptr->ctrotab, (size_t) ((ctromax + 1)    * sizeof (double)),
                     &keepptr->hghttab, (size_t) ((hghtmax + 1)    * sizeof (double)), NULL) == NULL) {
    errorPrint ("symbolKeepInit: out of memory (2)");
    memFree    (keepptr->keeptab);
    return     (1);
  }

  return (0);
}

/* This routine frees the contents
** of the given keep structure.
** It returns:
** - VOID  : in all cases.
*/

void
symbolKeepExit (
SymbolKeep * restrict const keepptr)              /*+ Keep structure to free +*/
{
  if (keepptr->keeptab != NULL)                   /* Free group leaders */
    memFree (keepptr->keeptab);
  if (keepptr->levftab != NULL)
    memFree (keepptr->levftab);
}

/* This routine adds the blocks for which
** the addition function returns !0 to the
** keep structure.
** It returns:
** - void  : in all cases.
*/

void
symbolKeepAdd (
SymbolKeep * restrict const keepptr,              /*+ Keep structure to update         +*/
const SymbolMatrix * const  symbptr,              /*+ Symbol structure to consider     +*/
int                      (* funcptr) (const SymbolKeepBlok * const, void * const),
void *                      dataptr)              /*+ Parameters for addition function +*/
{
  SymbolKeepBlok * restrict kblktax;              /* Based access to keep block array */
  unsigned char * restrict           keeptax;
  PASTIX_INT                       cblknum;
  PASTIX_INT                       bloknum;

  kblktax = keepptr->kblktab - symbptr->baseval;
  keeptax = keepptr->keeptab - symbptr->baseval;

  for (cblknum = 0; cblknum < symbptr->cblknbr; cblknum ++) {
    for (bloknum = symbptr->cblktab[cblknum].bloknum;
         bloknum < symbptr->cblktab[cblknum + 1].bloknum; bloknum ++) {
      if ((keeptax[bloknum] == 0) &&              /* If block can be added  */
          (funcptr (&kblktax[bloknum], dataptr) != 0)) /* And must be added */
        keeptax[bloknum] = 1;                     /* Add block              */
    }
  }
}

/* This routine deletes the blocks for which
** the deletion function returns !0 to the
** keep structure.
** It returns:
** - void  : in all cases.
*/

void
symbolKeepDel (
SymbolKeep * restrict const keepptr,              /*+ Keep structure to update         +*/
const SymbolMatrix * const  symbptr,              /*+ Symbol structure to consider     +*/
int                      (* funcptr) (const SymbolKeepBlok * const, void * const),
void *                      dataptr)              /*+ Parameters for addition function +*/
{
  SymbolKeepBlok * restrict kblktax;              /* Based access to keep block array */
  unsigned char * restrict           keeptax;
  PASTIX_INT                       cblknum;
  PASTIX_INT                       bloknum;

  kblktax = keepptr->kblktab - symbptr->baseval;
  keeptax = keepptr->keeptab - symbptr->baseval;

  for (cblknum = 0; cblknum < symbptr->cblknbr; cblknum ++) {
    for (bloknum = symbptr->cblktab[cblknum].bloknum;
         bloknum < symbptr->cblktab[cblknum + 1].bloknum; bloknum ++) {
      if ((keeptax[bloknum] != 0) &&              /* If block can be deleted  */
          (funcptr (&kblktax[bloknum], dataptr) != 0)) /* And must be deleted */
        keeptax[bloknum] = 0;                     /* Remove block             */
    }
  }
}

/* This routine recomputes the parameters
** of the blocks of the symbolic matrix
** with respect to the blocks which are
** flagged as kept.
*/

int
symbolKeepCompute (
SymbolKeep * restrict const keepptr,              /*+ Keep structure to update     +*/
const SymbolMatrix * const  symbptr)              /*+ Symbol structure to consider +*/
{
  SymbolKeepBlok * restrict kblktax;              /* Based access to keep block array              */
  SymbolCblk * restrict     cblktax;              /* Based access to column block array            */
  PASTIX_INT                       cblknum;              /* Based number of current column block          */
  SymbolBlok * restrict     bloktax;              /* Based access to incomplete block array        */
  PASTIX_INT                       bloknum;              /* Based number of current first free block slot */

  for (bloknum = 0; bloknum < symbptr->bloknbr; bloknum ++) { /* Reset parameter fields */
    keepptr->kblktab[bloknum].levfval = symbptr->bloktab[bloknum].levfval;
    keepptr->kblktab[bloknum].nupdval = 0;
    keepptr->kblktab[bloknum].ctrival = 0;
    keepptr->kblktab[bloknum].ctroval = 0;
  }

  cblktax = symbptr->cblktab - symbptr->baseval;
  bloktax = symbptr->bloktab - symbptr->baseval;
  kblktax = keepptr->kblktab - symbptr->baseval;

  for (cblknum = symbptr->baseval; cblknum < symbptr->baseval + symbptr->cblknbr - 1; cblknum ++) { /* Update pass (+ 3-loop pass, with forward contributions) */

    for (bloknum = cblktax[cblknum].bloknum + 1; ; bloknum ++) { /* For all existing off-diagonal blocks */
      PASTIX_INT                 levffac;
      PASTIX_INT                 blokctr;
      PASTIX_INT                 blokfac;
      PASTIX_INT                 nupdfac;

#ifdef SYMBOL_DEBUG
      if (bloknum >= cblktax[cblknum + 1].bloknum) { /* Should not happen */
        errorPrint ("symbolKeepCompute: internal error (1)");
        exit (1);
      }
#endif /* SYMBOL_DEBUG */

      kblktax[bloknum].ctroval ++;                /* Block sends contribution to facing diagonal block      */
      levffac = kblktax[bloknum].levfval;         /* Get level of fill of first block facing a column block */
      nupdfac = kblktax[bloknum].nupdval;
      while (((bloknum + 1) < cblktax[cblknum + 1].bloknum) &&
             (bloktax[bloknum + 1].cblknum == bloktax[bloknum].cblknum)) {
        bloknum ++;
        kblktax[bloknum].ctroval ++;              /* Block sends contribution to facing diagonal block            */
        if (kblktax[bloknum].levfval < levffac)   /* Get smallest fill value over blocks facing same column block */
          levffac = kblktax[bloknum].levfval;
        if (kblktax[bloknum].nupdval < nupdfac)   /* Get smallest fill value over blocks facing same column block */
          nupdfac = kblktax[bloknum].nupdval;
      }
      levffac ++;                                 /* Add +1 for computing levels of fill of contributions */

      blokctr = bloknum + 1;                      /* Search for first non-first-extra-diagonal block(s) */
      if (blokctr >= cblktax[cblknum + 1].bloknum) /* If none present, no contributions to account for  */
        break;                                    /* And end of processing for bloknum loop             */

      blokfac = cblktax[bloktax[bloknum].cblknum].bloknum + 1;
      do {
        kblktax[blokctr].ctroval ++;              /* All contributing extra-diagonals contribute to their facing block, which should exist */

        while (bloktax[blokfac].lrownum < bloktax[blokctr].frownum) { /* Skip non-facing blocks */
          blokfac ++;
          if (blokfac >= cblktax[bloktax[bloknum].cblknum + 1].bloknum) { /* Should not happen */
            errorPrint ("symbolKeepCompute: internal error (2)");
            exit (1);
          }
        }

        for ( ; ; ) {
          PASTIX_INT                 nupdval;            /* Number of updates of block(s) to be created */

          nupdval = 1 + ((kblktax[blokctr].nupdval > nupdfac) ? kblktax[blokctr].nupdval : nupdfac); /* Compute number of updates */

          kblktax[blokfac].ctrival ++;              /* Found facing block recieves contribution                        */
          if (kblktax[blokfac].levfval > (kblktax[blokctr].levfval + levffac)) /* Update level of fill of facing block */
            kblktax[blokfac].levfval = (kblktax[blokctr].levfval + levffac);
          kblktax[blokfac].nupdval = nupdval;

          if (bloktax[blokfac].lrownum >= bloktax[blokctr].lrownum) /* If facing block encompasses contributing block */
            break;                                  /* Skip to next contributing block                                */
          blokfac ++;                               /* Else find next facing block (which should exist)               */
          if (blokfac >= cblktax[bloktax[bloknum].cblknum + 1].bloknum) { /* Therefore, this should not happen        */
            errorPrint ("symbolKeepCompute: internal error (3)");
            exit (1);
          }
        }
      } while (++ blokctr < cblktax[cblknum + 1].bloknum);
    }
  }

  for (cblknum = symbptr->cblknbr - 1; cblknum >= 0; cblknum --) { /* Fill heights */
    if (symbptr->cblktab[cblknum + 1].bloknum - symbptr->cblktab[cblknum].bloknum > 1) { /* If extra-diagonal block present */
      PASTIX_INT                 hghtval;

      hghtval = kblktax[cblktax[bloktax[symbptr->cblktab[cblknum].bloknum + 1].cblknum].bloknum].hghtval + 1;

      for (bloknum = symbptr->cblktab[cblknum].bloknum; /* For all blocks of column block */
           bloknum < symbptr->cblktab[cblknum + 1].bloknum; bloknum ++)
        kblktax[bloknum].hghtval = hghtval;       /* Set height of blocks as height of column block */
    }
    else                                          /* Block is root or isolated */
      kblktax[symbptr->cblktab[cblknum].bloknum].hghtval = 1;
  }

  return (0);
}

/*+ This routine computes histogram
*** values for the given symbol matrix.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

static
void
symbolKeepHisto2 (
SymbolKeep * const          keepptr,              /*+ Keep structure to fill                 +*/
const SymbolMatrix * const  symbptr,              /*+ Symbol matrix to consider              +*/
int                      (* funcptr) (const SymbolKeepBlok * const, void * const),
void *                      dataptr,              /*+ Parameters for filtering function      +*/
const PASTIX_INT                   cblkmin,              /*+ Minimum column block index to consider +*/
const PASTIX_INT                   cblknbr)              /*+ Number of column blocks to consider    +*/
{
  if (cblknbr > 1) {                              /* If more than one column block, perform recursion */
    PASTIX_INT                 cmednum;
    double * restrict   levftmp;
    double * restrict   nupdtmp;
    double * restrict   ctritmp;
    double * restrict   ctrotmp;
    double * restrict   hghttmp;
    PASTIX_INT                 i;

    if (memAllocGroup ((void **)
        &levftmp, (size_t) ((keepptr->levfmax + 1) * sizeof (double)),
        &nupdtmp, (size_t) ((keepptr->nupdmax + 1) * sizeof (double)),
        &ctritmp, (size_t) ((keepptr->ctrimax + 1) * sizeof (double)),
        &ctrotmp, (size_t) ((keepptr->ctromax + 1) * sizeof (double)),
        &hghttmp, (size_t) ((keepptr->hghtmax + 1) * sizeof (double)), NULL) == NULL)
      return;
    memSet (levftmp, 0, (keepptr->levfmax + 1) * sizeof (double));
    memSet (nupdtmp, 0, (keepptr->nupdmax + 1) * sizeof (double));
    memSet (ctritmp, 0, (keepptr->ctrimax + 1) * sizeof (double));
    memSet (ctrotmp, 0, (keepptr->ctromax + 1) * sizeof (double));
    memSet (hghttmp, 0, (keepptr->hghtmax + 1) * sizeof (double));

    cmednum = cblknbr / 2;
    symbolKeepHisto2 (keepptr, symbptr, funcptr, dataptr, cblkmin, cmednum);
    symbolKeepHisto2 (keepptr, symbptr, funcptr, dataptr, cblkmin + cmednum, cblknbr - cmednum);

    for (i = 0; i <= keepptr->levfmax; i ++)
      keepptr->levftab[i] += levftmp[i];
    for (i = 0; i <= keepptr->nupdmax; i ++)
      keepptr->nupdtab[i] += nupdtmp[i];
    for (i = 0; i <= keepptr->ctrimax; i ++)
      keepptr->ctritab[i] += ctritmp[i];
    for (i = 0; i <= keepptr->ctromax; i ++)
      keepptr->ctrotab[i] += ctrotmp[i];
    for (i = 0; i <= keepptr->hghtmax; i ++)
      keepptr->hghttab[i] += hghttmp[i];
  }
  else {                                          /* Single column block */
    const SymbolKeepBlok * restrict kblktax;
    const SymbolBlok * restrict     bloktax;
    PASTIX_INT                             bloknum;
    double                          cblksiz;

    cblksiz  = (double) (symbptr->cblktab[cblkmin - symbptr->baseval].lcolnum - symbptr->cblktab[cblkmin - symbptr->baseval].fcolnum + 1);

    kblktax = keepptr->kblktab - symbptr->baseval;
    bloktax = symbptr->bloktab - symbptr->baseval;
    for (bloknum = symbptr->cblktab[cblkmin - symbptr->baseval].bloknum;
         bloknum < symbptr->cblktab[cblkmin - symbptr->baseval + 1].bloknum; bloknum ++) {
      if ((funcptr == NULL) || (funcptr (&kblktax[bloknum], dataptr) != 0)) {
        double                bloksiz;

        bloksiz = cblksiz * (double) (bloktax[bloknum].lrownum - bloktax[bloknum].frownum + 1);

        keepptr->levftab[kblktax[bloknum].levfval] += bloksiz;
        keepptr->nupdtab[kblktax[bloknum].nupdval] += bloksiz;
        keepptr->ctritab[kblktax[bloknum].ctrival] += bloksiz;
        keepptr->ctrotab[kblktax[bloknum].ctroval] += bloksiz;
        keepptr->hghttab[kblktax[bloknum].hghtval] += bloksiz;
      }
    }
  }
}

int
symbolKeepHisto (
SymbolKeep * const          keepptr,              /*+ Keep structure to fill    +*/
const SymbolMatrix * const  symbptr,              /*+ Symbol matrix to consider +*/
int                      (* funcptr) (const SymbolKeepBlok * const, void * const),
void *                      dataptr)              /*+ Parameters for addition function +*/
{
  memSet (keepptr->levftab, 0, (keepptr->levfmax + 1) * sizeof (double));
  memSet (keepptr->nupdtab, 0, (keepptr->nupdmax + 1) * sizeof (double));
  memSet (keepptr->ctritab, 0, (keepptr->ctrimax + 1) * sizeof (double));
  memSet (keepptr->ctrotab, 0, (keepptr->ctromax + 1) * sizeof (double));
  memSet (keepptr->hghttab, 0, (keepptr->hghtmax + 1) * sizeof (double));

  symbolKeepHisto2 (keepptr, symbptr, funcptr, dataptr, symbptr->baseval, symbptr->cblknbr);

  return (0);
}

/* This routine purges the symbol matrix
** of its removed blocks and compacts
** the block structure.
** It returns:
** - void  : in all cases.
*/

int
symbolKeepPurge (
SymbolKeep * restrict const   keepptr,            /*+ Keep structure for kept blocks   +*/
SymbolMatrix * restrict const symbptr)            /*+ Symbol structure to purge        +*/
{
  unsigned char * restrict           keeptax;              /* Based access to keep block flag array */
  SymbolBlok *              bloktax;
  PASTIX_INT                       cblknum;
  PASTIX_INT                       bloknew;

  keeptax = keepptr->keeptab - symbptr->baseval;
  bloktax = symbptr->bloktab - symbptr->baseval;

  for (cblknum = 0, bloknew = symbptr->baseval;
       cblknum < symbptr->cblknbr; cblknum ++) {
    PASTIX_INT                       bloknum;

    for (bloknum = symbptr->cblktab[cblknum].bloknum,
         symbptr->cblktab[cblknum].bloknum = bloknew;
         bloknum < symbptr->cblktab[cblknum + 1].bloknum; bloknum ++) {
      if (keeptax[bloknum] != 0) {                /* If block must be kept           */
        bloktax[bloknew] = bloktax[bloknum];      /* Copy block to its new position  */
        keeptax[bloknew] = 1;                     /* Keep block in new numbering     */
        bloknew ++;                               /* One more block in symbol matrix */
      }
    }
  }
  symbptr->cblktab[cblknum].bloknum = bloknew;    /* Set end of column block array */

  symbptr->bloknbr = bloknew - symbptr->baseval;  /* Update number of blocks */

#ifdef SYMBOL_DEBUG
  if (symbolCheck (symbptr) != 0) {               /* Should not happen */
    errorPrint ("symbolKeepPurge: internal error");
    return (1);
  }
#endif /* SYMBOL_DEBUG */

  return (0);
}

/* This routine creates the set of GnuPlot
** histogram files for the kept blocks.
** It returns:
** - 0   : in case of success.
** - !0  : in case of error.
*/

static
int
symbolKeepView2 (
const double * const        dataptr,              /*+ Data field to write       +*/
const PASTIX_INT                   datanbr,              /*+ Number of data fields     +*/
const double                datamax,              /*+ Sum of all of block data  +*/
const char * const          nameptr,              /*+ Base name for output file +*/
const char * const          commptr)              /*+ Comment string            +*/
{
  FILE *                      fileptr;
  double                      datasum;
  PASTIX_INT                         i;

  if ((fileptr = fopen (nameptr, "w")) == NULL) {
    errorPrint ("symbolKeepView2: cannot open output file");
    return     (1);
  }

  fprintf (fileptr, "%s", commptr);                     /* Write comment */

  for (i = 0, datasum = 0.0L; i < datanbr; i ++) { /* Accumulate and write data */
    datasum += dataptr[i];
    fprintf (fileptr, "%ld\t%le\n", (long) i, (double) datasum / (double) datamax);
  }

  fclose (fileptr);

  return (0);
}

int
symbolKeepView (
const SymbolKeep * const    keepptr,              /*+ Symbol matrix to           +*/
const double                nnzlmax,              /*+ Sum of all of block data   +*/
const char * const          nameptr)              /*+ Base name for output files +*/
{
  char                        buftab1[1024];
  char                        buftab2[1024];

  sprintf (buftab1, "%s_levf.dat", nameptr);
  sprintf (buftab2, "# Plot of levf for %s\n", nameptr);
  symbolKeepView2 (keepptr->levftab, keepptr->levfmax + 1, nnzlmax, buftab1, buftab2);
  sprintf (buftab1, "%s_nupd.dat", nameptr);
  sprintf (buftab2, "# Plot of nupd for %s\n", nameptr);
  symbolKeepView2 (keepptr->nupdtab, keepptr->nupdmax + 1, nnzlmax, buftab1, buftab2);
  sprintf (buftab1, "%s_ctri.dat", nameptr);
  sprintf (buftab2, "# Plot of ctri for %s\n", nameptr);
  symbolKeepView2 (keepptr->ctritab, keepptr->ctrimax + 1, nnzlmax, buftab1, buftab2);
  sprintf (buftab1, "%s_ctro.dat", nameptr);
  sprintf (buftab2, "# Plot of ctro for %s\n", nameptr);
  symbolKeepView2 (keepptr->ctrotab, keepptr->ctromax + 1, nnzlmax, buftab1, buftab2);
  sprintf (buftab1, "%s_hght.dat", nameptr);
  sprintf (buftab2, "# Plot of hght for %s\n", nameptr);
  symbolKeepView2 (keepptr->hghttab, keepptr->hghtmax + 1, nnzlmax, buftab1, buftab2);

  return (0);
}
