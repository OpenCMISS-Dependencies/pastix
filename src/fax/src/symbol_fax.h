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
** $Id: symbol_fax.h 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_fax.h                            **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Jean ROMAN (v0.0)                       **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the symbolic factorization routine. **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     24 sep 1998     **/
/**                # Version 0.1  : from : 04 apr 1999     **/
/**                                 to     01 may 1999     **/
/**                # Version 1.0  : from : 01 jun 2002     **/
/**                                 to     05 jun 2002     **/
/**                # Version 1.1  : from : 26 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 3.0  : from : 03 mar 2004     **/
/**                                 to     03 mar 2004     **/
/**                                                        **/
/************************************************************/

#define SYMBOL_FAX_H

/*
**  The defines.
*/

/* Prime number for hashing vertex numbers. */

#define SYMBOL_FAX_HASHPRIME        17            /*+ Prime number for hashing +*/

/*
**  The type and structure definitions.
*/

/*+ The chained column block structure. These
    blocks are chained in a single linked list
    for block merge with blocks of left columns. +*/

typedef struct SymbolFaxTlok_ {
  PASTIX_INT                       frownum;              /*+ First row index            +*/
  PASTIX_INT                       lrownum;              /*+ Last row index (inclusive) +*/
  PASTIX_INT                       cblknum;              /*+ Facing column block        +*/
  PASTIX_INT                       nextnum;              /*+ Index of next block        +*/
} SymbolFaxTlok;
