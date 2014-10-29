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
** $Id: dof.h 316 2005-06-06 16:17:44Z ramet $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : dof.h                                   **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the DOF handling structure.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 07 oct 1998     **/
/**                                 to     16 oct 1998     **/
/**                # Version 1.0  : from : 06 jun 2002     **/
/**                                 to     06 jun 2002     **/
/**                # Version 3.0  : from : 28 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

#ifndef DOF_H
#define DOF_H

/*
**  The type and structure definitions.
*/

/*+ The DOF structure. This structure is
    always associated to a Graph structure,
    which holds the base value.             +*/

typedef struct Dof_ {
  PASTIX_INT                       baseval;              /*+ Base value for indexing                                       +*/
  PASTIX_INT                       nodenbr;              /*+ Number of nodes in DOF array                                  +*/
  PASTIX_INT                       noddval;              /*+ DOF value for every node (if noddtab == NULL, 0 else)         +*/
  PASTIX_INT * restrict            noddtab;              /*+ Array of node->first DOF indexes (if noddval == 0) [+1,based] +*/
} Dof;

/*
**  The function prototypes.
*/

#ifndef DOF
#define static
#endif

int                         dofInit             (Dof * const deofptr);
void                        dofExit             (Dof * const deofptr);
int                         dofLoad             (Dof * const deofptr, FILE * const stream);
int                         dofSave             (const Dof * const deofptr, FILE * const stream);
void                        dofConstant         (Dof * const deofptr, const PASTIX_INT baseval, const PASTIX_INT nodenbr, const PASTIX_INT noddval);
#ifdef GRAPH_H
int                         dofGraph            (Dof * const deofptr, const Graph * const grafptr, const PASTIX_INT, const PASTIX_INT * const peritab);
#endif /* GRAPH_H */

#undef static

/*
**  The macro definitions.
*/

#ifdef DOF_CONSTANT
#define noddVal(deofptr,nodenum)    ((deofptr)->baseval + (deofptr)->noddval * ((nodenum) - (deofptr)->baseval))
#define noddDlt(deofptr,nodenum)    ((deofptr)->noddval)
#else /* DOF_CONSTANT */
#define noddVal(deofptr,nodenum)    (((deofptr)->noddtab != NULL) ? (deofptr)->noddtab[(deofptr)->baseval + (nodenum)] : ((deofptr)->baseval + (deofptr)->noddval * ((nodenum) - (deofptr)->baseval)))
#define noddDlt(deofptr,nodenum)    (((deofptr)->noddtab != NULL) ? ((deofptr)->noddtab[(deofptr)->baseval + (nodenum) + 1] - (deofptr)->noddtab[(deofptr)->baseval + (nodenum)]) : (deofptr)->noddval)
#endif /* DOF_CONSTANT */
#endif /* DOF_H */
