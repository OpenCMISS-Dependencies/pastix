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
** $Id: dof_io.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : dof_io.c                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the input/output        **/
/**                routines for the DOF structure.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 07 oct 1998     **/
/**                                 to     07 oct 1998     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/
#define DOF_IO

#include "common_pastix.h"
#include "dof.h"

/****************************************/
/*                                      */
/* The DOF structure handling routines. */
/*                                      */
/****************************************/

/*+ This routine saves the given DOF
*** structure to the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
dofSave (
const Dof * const           deofptr,
FILE * const                stream)
{
  const PASTIX_INT *         noddtnd;
  const PASTIX_INT *         noddptr;
  PASTIX_INT                 noddnum;
  int                 o;

  o = (fprintf (stream, "0\n%ld\t%ld\t%ld\n\n",   /* Write file header */
                (long) deofptr->nodenbr,
                (long) deofptr->noddval,
                (long) deofptr->baseval) == EOF);
  if (deofptr->noddtab != NULL) {
    for (noddptr = deofptr->noddtab, noddtnd = noddptr + deofptr->nodenbr, noddnum = 1;
         (noddptr < noddtnd) && (o == 0); noddptr ++, noddnum ++) {
      o = (fprintf (stream, "%ld%c",
                    (long) *noddptr,
                    ((noddnum % 8) == 0) ? '\n' : '\t') == EOF);
    }
    o |= (fprintf (stream, "%ld\n",
                   (long) *noddptr) == EOF);
  }

  return (o);
}
