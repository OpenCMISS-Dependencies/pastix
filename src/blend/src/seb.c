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
#include "error.h"
#include "elimin.h"
#include "cand.h"
#include "cost.h"
#ifdef WITH_SCOTCH
#ifdef DISTRIBUTED
#include <mpi.h>
#include "ptscotch.h"
#else
#include "scotch.h"
#endif
#endif

void tree2graph(const EliminTree * etree, const CostMatrix *costmtx,  SCOTCH_Graph *graph)
{
  PASTIX_INT i, j;
  PASTIX_INT *sonsnbr;
  PASTIX_INT *verttab;
  PASTIX_INT *verttabtmp;
  PASTIX_INT *velotab;
  PASTIX_INT *edgetab;
  PASTIX_INT vertnbr = 0;
  PASTIX_INT edgenbr = 0;

  /* MALLOC_INTERN(graph, 1, SCOTCH_Graph); */

  /* Calcul du nombre de fils par noeud */
  MALLOC_INTERN(sonsnbr, etree->nodenbr, PASTIX_INT);
  for(i=0; i<etree->nodenbr; i++)
    {
      sonsnbr[i] = etree->nodetab[i].sonsnbr+1;
      edgenbr   += etree->nodetab[i].sonsnbr+1;
    }
 
  /* Suppression de l'arete de la racine vers -1 */
  edgenbr--;
  sonsnbr[ROOT(etree)]--;

  /* Ajout de la symmetrie */
/*   for(i=0; i<etree->nodenbr; i++) */
/*     { */
/*       PASTIX_INT deb = etree->nodetab[i].fsonnum; */
/*       PASTIX_INT fin = etree->nodetab[i].sonsnbr+deb; */
      
/*       for(j=deb; j<fin; j++) */
/* 	sonsnbr[j]++; */
/*     } */

  vertnbr    = etree->nodenbr;
  MALLOC_INTERN(verttab,    vertnbr+1, PASTIX_INT);
  MALLOC_INTERN(verttabtmp, vertnbr,   PASTIX_INT);
  MALLOC_INTERN(velotab,    vertnbr,   PASTIX_INT);
  MALLOC_INTERN(edgetab,    edgenbr,   PASTIX_INT);  

  /* Calcul de verttab (indice du premier fils dans edgetab) */
  verttab[0]    = 0;
  verttabtmp[0] = 0;
  velotab[0]    = MAX(1, (PASTIX_INT)(costmtx->cblktab[0].total*1000000));
  for(i=1; i<vertnbr; i++)
    {
      verttab[i]    = verttab[i-1] + sonsnbr[i-1];
      verttabtmp[i] = verttab[i-1] + sonsnbr[i-1];
      velotab[i]    = MAX(1,(PASTIX_INT)(costmtx->cblktab[i].total*1000000));
    }
  verttab[vertnbr] = verttab[vertnbr-1] + sonsnbr[vertnbr-1];

  /* Remplissage de edgetab */
  for(i=0; i<vertnbr; i++)
    {
      PASTIX_INT deb = etree->nodetab[i].fsonnum;
      PASTIX_INT fin = etree->nodetab[i].sonsnbr+deb;
      
      for(j=deb; j<fin; j++)
	{
	  edgetab[verttabtmp[i]] = etree->sonstab[j];
	  verttabtmp[i]++;
	  edgetab[verttabtmp[etree->sonstab[j]]] = i;
	  verttabtmp[etree->sonstab[j]]++;
	}
    }

  SCOTCH_graphBuild(graph,
		    0,
		    etree->nodenbr,
		    verttab,
		    verttab+1,
		    velotab,
		    NULL,
		    edgenbr,
		    edgetab,
		    NULL);
  
  if (SCOTCH_graphCheck(graph))
    {
      EXIT(MOD_BLEND, INTERNAL_ERR);
    }
}

void SetCandtab(Cand *candtab, const PASTIX_INT *parttab, const PASTIX_INT vertnbr)
{
  PASTIX_INT i;

  for (i=0; i<vertnbr; i++)
    {
      candtab[i].fcandnum = parttab[i];
      candtab[i].lcandnum = parttab[i];
    }
}
