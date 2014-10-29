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
#ifndef BULLES_H
#define BULLES_H

/*+ The node structure. +*/
typedef struct BubbleTreeNode_ {
  int                       fathnum;      /*+ index of the father node               +*/
  int                       sonsnbr;
  int                       fsonnum;
  int                       fcandnum;     /*+ first bubble proc                      +*/
  int                       lcandnum;     /*+ last bubble proc                       +*/
  double                    costlevel;    /*+ cost of way from the root to this node +*/
  int                       treelevel;    /*+ cost of way from the root to this node +*/
  PASTIX_INT                       priomin;      /*+ Minimal priority of tasks owned by the bubble +*/
  PASTIX_INT                       priomax;      /*+ Maximal priority of tasks owned by the bubble +*/
  Queue *                   taskheap;     /*+ Liste de taches de la bulle            +*/
} BubbleTreeNode;


/*+ The bubble tree. +*/
typedef struct BubbleTree_ {
  int                       leavesnbr;    /*+ Number of leaves in tree  +*/
  int                       nodenbr;      /*+ Number of nodes           +*/
  int                       nodemax;      /*+ Number max of nodes       +*/
  int                      *sonstab;      /*+ Array of node             +*/
  BubbleTreeNode           *nodetab;      /*+ Array of node             +*/
} BubbleTree;


#define BFATHER(btree, r)     (btree)->nodetab[r].fathnum
#define BNBSON(btree, r)      (btree)->nodetab[r].sonsnbr
#define BROOT(btree)          (btree)->nodetab[(btree)->leavesnbr]
#define BSON(btree, n, r)     (btree)->sonstab[(btree)->nodetab[n].fsonnum + r ]

void  Bubble_InitTree  (BubbleTree *, int);
void  Bubble_Free      (BubbleTree *);
int   Bubble_Add       (BubbleTree *, PASTIX_INT, PASTIX_INT, double, PASTIX_INT);
void  Bubble_BuildTree (const BubbleTree *);
void  Bubble_Print     (const BubbleTree *, const double *, double, FILE*);

#endif /* BULLES_H */






