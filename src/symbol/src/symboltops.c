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
** $Id: symboltops.c 316 2005-06-06 16:17:44Z ramet $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symboltops.c                            **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This program produces a PostScript(tm)  **/
/**                file from a symbol matrix file.         **/
/**                                                        **/
/**   DATES      : # Version 0.1  : from : 21 mar 2002     **/
/**                                 to     21 mar 2002     **/
/**                # Version 3.0  : from : 29 sep 2004     **/
/**                                 to     29 sep 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "common_pastix.h"
#include "symbol.h"

/******************************/
/*                            */
/* This is the main function. */
/*                            */
/******************************/

int
main (int argc, char *argv[])
{
  FILE *              symbfile;                   /* Symbol input stream      */
  SymbolMatrix        symbdat;                    /* Mesh graph               */
  FILE *              postfile;                   /* PostScript output stream */
  int                 errval;
  char *              errstr;

  errorProg ("symboltops");

  if (argc > 3) {
    errorPrint ("usage: symboltops symbol_file postscript_tm_file");
    return     (1);
  }

  errval = 0;
  errstr = NULL;

  symbfile = stdin;                               /* Assume standard streams */
  postfile = stdout;

  symbolInit (&symbdat);                          /* Initialize symbol structure */

  if (argc > 1) {                                 /* If input file name provided */
    if ((strcmp (argv[1], "-") != 0) &&
        ((symbfile = fopen (argv[1], "r")) == NULL)) {
      errval = 1;
      errstr = "cannot open symbol file";
    }
    else if (argc > 2) {                          /* If output file name provided */
      if ((strcmp (argv[2], "-") != 0) &&
          ((postfile = fopen (argv[2], "w+")) == NULL)) {
        errval = 1;
        errstr = "cannot open PostScript(tm) file";
      }
    }
  }

  if (errval == 0) {
    if (symbolLoad (&symbdat, symbfile) != 0) {   /* Read symbol data */
      errval = 1;
      errstr = "cannot read symbol file";
    }
    else if (symbolDraw (&symbdat, postfile) != 0) {
      errval = 1;
      errstr = "cannot write PostScript(tm) file";
    }
  }

  if (errstr != NULL)
    errorPrint (errstr);

  symbolExit (&symbdat);

  return (errval);
}
