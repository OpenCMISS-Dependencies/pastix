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
** $Id: symbol_draw.h 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_draw.h                           **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the symbolic matrix drawing         **/
/**                routine.                                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 29 sep 1998     **/
/**                                 to     30 sep 1998     **/
/**                # Version 1.0  : from : 26 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 1.3  : from : 10 jun 2003     **/
/**                                 to     10 jun 2003     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

#define SYMBOL_DRAW_H

#ifdef CXREF_DOC
#include "common_pastix.h"
#ifndef SYMBOL_H
#include "symbol.h"
#endif /* SYMBOL_H */
#endif /* CXREF_DOC */

/*
**  The definitions.
*/

/*+ Generic PostScript (tm) output definitions. +*/

#define SYMBOL_PSDPI                72            /*+ PostScript dots-per-inch            +*/
#define SYMBOL_PSPICTSIZE           6.6           /*+ PostScript picture size (in inches) +*/
/*
**  The function prototypes.
*/

#ifndef SYMBOL_DRAW
#define static
#endif

#undef static
