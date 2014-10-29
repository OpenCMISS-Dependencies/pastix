/* Copyright 2007 ENSEIRB, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
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
/**   NAME       : ocv.c                                   **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of an ordering file converter.     **/
/**                This module contains the main function. **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 16 dec 2007     **/
/**                                 to     16 dec 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define OCV

#include "common.h"
#include "order.h"
#include "ocv.h"

/*
**  The static and global variables
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
File                        C_fileTab[C_FILENBR] = { /* The file array; public  */
                              { "-", NULL, "r" },
                              { "-", NULL, "r" },
                              { "-", NULL, "w" } };

static const char *         C_usageList[] = {
  "ocv [<input ordering file> [<inpput mapping file> [<output block ordering file>]]]",
  "  -h  : Display this help",
  "  -V  : Print program version and copyright",
  NULL };

/*****************************/
/*                           */
/* This is the main function */
/*                           */
/*****************************/

int
main (
int                         argc,
char *                      argv[])
{
  Order               ordedat;
  int                 i;

  errorProg ("ocv");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  for (i = 0; i < C_FILENBR; i ++)                /* Set default stream pointers */
    C_fileTab[i].pntr = (C_fileTab[i].mode[0] == 'r') ? stdin : stdout;
  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes */
    if ((argv[i][0] != '+') &&                    /* If found a file name      */
        ((argv[i][0] != '-') || (argv[i][1] == '\0'))) {
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        C_fileTab[C_fileNum ++].name = argv[i];
      else {
        errorPrint ("main: too many file names given");
        return     (1);
      }
    }
    else {                                       /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                               /* Give help */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'V' :
          fprintf (stderr, "ocv - F. Pellegrini\n");
          fprintf (stderr, "Copyright 2007 ENSEIRB, INRIA & CNRS, France\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        default :
          errorPrint ("main: unprocessed option (\"%s\")", argv[i]);
          return     (1);
      }
    }
  }

  for (i = 0; i < C_FILENBR; i ++) {             /* For all file names     */
    if ((C_fileTab[i].name[0] != '-') ||         /* If not standard stream */
        (C_fileTab[i].name[1] != '\0')) {
      if ((C_fileTab[i].pntr = fopen (C_fileTab[i].name, C_fileTab[i].mode)) == NULL) { /* Open the file */
        errorPrint ("main: cannot open file (%d)", i);
        return     (1);
      }
    }
  }

  orderInit (&ordedat);

  if (orderDataLoad (&ordedat, C_filepntrordinp, C_filepntrmapinp) != 0) {
    errorPrint ("main: cannot load Scotch files");
    return     (1);
  }    
  if (orderCheck (&ordedat) != 0) {
    errorPrint ("main: invalid ordering");
    return     (1);
  }

  if (orderSave (&ordedat, C_filepntrordout) != 0) {
    errorPrint ("main: cannot save block ordering");
    return     (1);
  }

  for (i = 0; i < C_FILENBR; i ++) {             /* For all file names     */
    if ((C_fileTab[i].name[0] != '-') ||         /* If not standard stream */
        (C_fileTab[i].name[1] != '\0')) {
      fclose (C_fileTab[i].pntr);                /* Close the stream */
    }
  }

  return (0);
}

int
orderDataLoad (
Order * restrict const      ordeptr,
FILE * restrict const       permstream,
FILE * restrict const       mappstream)
{
  PASTIX_INT * restrict      peritax;
  PASTIX_INT * restrict      permtax;
  PASTIX_INT * restrict      mapptax;
  PASTIX_INT                 mappval;
  PASTIX_INT                 vertnbr;
  PASTIX_INT                 vertnnd;
  PASTIX_INT                 vertmin;
  PASTIX_INT                 vertmax;
  PASTIX_INT                 vertnum;
  PASTIX_INT                 cblknum;
  PASTIX_INT                 baseval;

  if (intLoad (permstream, &vertnbr) != 1) {
    errorPrint ("orderDataLoad: bad input (1)");
    return     (1);
  }

  if (memAllocGroup ((void **) (void *)
                     &ordeptr->rangtab, (size_t) ((vertnbr + 1) * sizeof (PASTIX_INT)),
                     &ordeptr->permtab, (size_t) ( vertnbr      * sizeof (PASTIX_INT)),
                     &ordeptr->peritab, (size_t) ( vertnbr      * sizeof (PASTIX_INT)), NULL) == NULL) {
    errorPrint ("orderDataLoad: out of memory (1)");
    return     (1);
  }
  if ((mapptax = memAlloc (vertnbr * sizeof (PASTIX_INT))) == NULL) {
    errorPrint ("orderDataLoad: out of memory (2)");
    memFree    (ordeptr->rangtab);
    return     (1);
  }

  vertmin = vertnbr + 1;
  vertmax = -1;
  for (vertnum = 0; vertnum < vertnbr; vertnum ++) {
    PASTIX_INT                 vertval;

    if (intLoad (permstream, &vertval) +          /* Read item data */
        intLoad (permstream, &ordeptr->permtab[vertnum]) != 2) {
      errorPrint ("orderDataLoad: bad input (2)");
      memFree    (mapptax);
      memFree    (ordeptr->rangtab);
      return     (1);
    }
    if (vertval > vertmax)
      vertmax = vertval;
    if (vertval < vertmin)
      vertmin = vertval;
  }
  if ((vertmin < 0) || (vertmax != (vertmin + vertnbr - 1))) {
    errorPrint ("orderDataLoad: bad input (3)");
    memFree    (mapptax);
    memFree    (ordeptr->rangtab);
    return     (1);
  }
  baseval = vertmin;

  peritax = ordeptr->peritab - baseval;
  permtax = ordeptr->permtab - baseval;
  for (vertnum = baseval, vertnnd = vertnbr + baseval; vertnum < vertnnd; vertnum ++) /* Compute inverse permutation */
    peritax[permtax[vertnum]] = vertnum;

  if ((intLoad (mappstream, &mappval) != 1) ||
      (mappval != vertnbr)) {
    errorPrint ("orderDataLoad: bad input (4)");
    memFree    (mapptax);
    memFree    (ordeptr->rangtab);
    return     (1);
  }

  vertmin = vertnbr + 1;
  vertmax = -1;
  for (vertnum = 0; vertnum < vertnbr; vertnum ++) {
    PASTIX_INT                 vertval;

    if (intLoad (mappstream, &vertval) +          /* Read item data */
        intLoad (mappstream, &mapptax[vertnum]) != 2) {
      errorPrint ("orderDataLoad: bad input (5)");
      memFree    (mapptax);
      memFree    (ordeptr->rangtab);
      return     (1);
    }
    if (vertval > vertmax)
      vertmax = vertval;
    if (vertval < vertmin)
      vertmin = vertval;
  }
  if ((vertmin < baseval) || (vertmax != (vertmin + vertnbr - 1))) {
    errorPrint ("orderDataLoad: bad input (6)");
    memFree    (mapptax);
    memFree    (ordeptr->rangtab);
    return     (1);
  }

  mapptax -= baseval;
  for (vertnum = cblknum = 0, mappval = -1; vertnum < vertnbr; vertnum ++) {
    PASTIX_INT                 mapptmp;

    mapptmp = mapptax[ordeptr->peritab[vertnum]];
    if (mapptmp != mappval) {
      mappval = mapptmp;
      ordeptr->rangtab[cblknum ++] = vertnum + baseval;
    }
  }
  ordeptr->rangtab[cblknum] = vertnum + baseval;
  ordeptr->cblknbr = cblknum;

  memFree (mapptax + baseval);

  return (0);
}
