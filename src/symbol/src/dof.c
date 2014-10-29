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
** $Id: dof.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : dof.c                                   **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the general purpose     **/
/**                routines for the DOF structure.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 07 oct 1998     **/
/**                                 to     14 oct 1998     **/
/**                # Version 3.0  : from : 28 feb 2004     **/
/**                                 to     28 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DOF

#include "common_pastix.h"
#ifdef WITH_SCOTCH
#ifdef DISTRIBUTED
#include <mpi.h>
#include "ptscotch.h"
#else
#include "scotch.h"
#endif
#endif
#include "dof.h"

/******************************/
/*                            */
/* The DOF handling routines. */
/*                            */
/******************************/

/*+ This routine initializes
*** the given DOF structure.
*** It returns:
*** - 0  : in all cases.
+*/

int
dofInit (
Dof * const                 deofptr)
{
  deofptr->baseval = 0;
  deofptr->nodenbr = 0;
  deofptr->noddval = 1;                           /* Set constant, non zero, number of DOFs */
  deofptr->noddtab = NULL;

  return (0);
}

/*+ This routine frees the contents
*** of the given DOF structure.
*** It returns:
*** - VOID  : in all cases.
+*/

void
dofExit (
Dof * const                 deofptr)
{
  if (deofptr->noddtab != NULL)
    memFree (deofptr->noddtab);

#ifdef DOF_DEBUG
  dofInit (deofptr);
#endif /* DOF_DEBUG */
}

/*+ This routine sets the number of DOFs
*** per node to a constant value.
*** It returns:
*** - VOID  : in all cases.
+*/

void
dofConstant (
Dof * const                 deofptr,
const PASTIX_INT                   baseval,
const PASTIX_INT                   nodenbr,
const PASTIX_INT                   noddval)
{
  deofptr->baseval = baseval;
  deofptr->nodenbr = nodenbr;
  if (deofptr->noddtab != NULL) {                 /* If DOF array already allocated */
    memFree (deofptr->noddtab);                   /* It is no longer of use         */
    deofptr->noddtab = NULL;
  }
  deofptr->noddval = noddval;
}

/*+ This routine builds the DOF index
*** array from the graph vertex array.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

#ifdef WITH_SCOTCH
int
dofGraph (
Dof * const                 deofptr,              /*+ DOF index array to build [based]            +*/
const SCOTCH_Graph * const         grafptr,              /*+ Matrix adjacency structure [based]          +*/
const PASTIX_INT                   deofval,              /*+ DOFs per node if no graph vertex load array +*/
const PASTIX_INT * const           peritab)              /*+ Inverse vertex->node permutation array      +*/
{
  PASTIX_INT                 baseval;
  PASTIX_INT                 vertnbr;
  PASTIX_INT *               velotab;
  PASTIX_INT                 edgenbr;
  (void)peritab;

  SCOTCH_graphData (grafptr,
                    (SCOTCH_Num *) &baseval,
                    (SCOTCH_Num *) &vertnbr,
                    NULL, NULL,
                    (SCOTCH_Num **)&velotab,
                    NULL,
                    (SCOTCH_Num *) &edgenbr,
                    NULL,
                    NULL);

  deofptr->baseval = baseval;
  deofptr->nodenbr = vertnbr;
  if (velotab == NULL) {                          /* If no vertex weight (i.e. DOF) array */
    deofptr->noddtab = NULL;                      /* No DOF array                         */
    deofptr->noddval = deofval;                   /* Get node DOF value                   */
  }
  else {                                          /* Vertex load array present */
#ifdef DOF_CONSTANT
    deofptr->noddtab = NULL;                      /* No DOF array */
    deofptr->noddval = deofval;
#else /* DOF_CONSTANT */
    const PASTIX_INT * restrict  velotax;                /* Based access to grafptr->velotab  */
    PASTIX_INT                   nodenum;                /* Number of current node            */
    PASTIX_INT *                 noddtnd;                /* Pointer to end of DOF index array */
    PASTIX_INT *                 noddptr;                /* Pointer to current DOF index      */
    const PASTIX_INT *           periptr;

    deofptr->noddval = 0;                         /* DOF values are not constant */
    if ((deofptr->noddtab = (PASTIX_INT *) memAlloc ((vertnbr + 1) * sizeof (PASTIX_INT))) == NULL) {
      errorPrint ("dofGraph: out of memory");
      return     (1);
    }
    for (noddptr = deofptr->noddtab, noddtnd = noddptr + vertnbr,
         periptr = peritab, nodenum = baseval,
         velotax = grafptr->velotab - baseval;
         noddptr < noddtnd; noddptr ++, periptr ++) {
      *noddptr = nodenum;                         /* Set index to DOF array        */
      nodenum += velotax[*periptr];               /* Add number of DOFs for vertex */
    }
    *noddptr = nodenum;                           /* Set end of DOF array */
#endif /* DOF_CONSTANT */
  }

  return (0);
}
#endif
