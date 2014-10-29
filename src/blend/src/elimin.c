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
#include <assert.h>

#include "common_pastix.h"
#include "elimin.h"


PASTIX_INT egraphInit(EliminGraph *egraph)
{
  egraph->baseval = 0;
  egraph->vertnbr = 0;
  egraph->verttab = NULL;
  egraph->inbltab = NULL;
  egraph->ownetab = NULL;
  return 1;
}

void egraphExit(EliminGraph *egraph)
{
  memFree_null(egraph->verttab);
  memFree_null(egraph->inbltab);
  memFree_null(egraph->ownetab);
  memFree_null(egraph);
}


PASTIX_INT treeInit(EliminTree *etree)
{
  etree->baseval = 0;
  etree->nodenbr = 0;
  etree->nodetab = NULL;
  etree->sonstab = NULL;
  return 1;
}

void treeExit(EliminTree *etree)
{
    memFree_null(etree->nodetab);
    memFree_null(etree->sonstab);
    memFree_null(etree);
}


void treePlot(EliminTree *etree, FILE *out)
{
  PASTIX_INT i;

  fprintf(out,
	  "digraph G {\n"
	  "\tcolor=white\n"
	  "rankdir=BT;\n");

  for (i=0;  i < etree->nodenbr; i++)
    {
      if ((etree->nodetab[i]).fathnum == -1)
	continue;
      fprintf(out, "\t\"%ld\"->\"%ld\"\n", (long)i, (long)((etree->nodetab[i]).fathnum));
    }

  fprintf(out, "}\n");
}
	    
