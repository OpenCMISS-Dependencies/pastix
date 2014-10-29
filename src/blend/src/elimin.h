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
/************************************************************/
/**                                                        **/
/**   NAME       : elimin.h                                **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the elimination tree.               **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     27 jul 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ The node structure. +*/
typedef struct TreeNode_ {
  PASTIX_INT                       sonsnbr;              /*+ Number of sons                          +*/
  PASTIX_INT                       fathnum;              /*+ index of the father node                +*/
  PASTIX_INT                       fsonnum;              /*+ index of first son                      +*/
} TreeNode;

/*+ The elimination tree. +*/

typedef struct EliminTree_ {
  PASTIX_INT                       baseval;              /*+ Base value for numberings         +*/
  PASTIX_INT                       nodenbr;              /*+ Number of nodes                   +*/
  TreeNode   *              nodetab;              /*+ Array of node          [+1,based] +*/
  PASTIX_INT        *              sonstab;              /*+ Sons index of nodes               +*/  
} EliminTree;


/*+ The elimination graph. +*/
/*+ we only need the in-edges between graph vertex
    Out-edges can be found with the symbol matrix data +*/
/* OIMBE innbr ne sert pas necessairement !*/
typedef struct EliminVertex_ {
  PASTIX_INT                       innum;                /*+ index of first in-bloc            +*/
  PASTIX_INT                       innbr;                /*+ number of in-blocs                +*/   
} EliminVertex;

typedef struct EliminGraph_ {
  PASTIX_INT                       baseval;              /*+ Base value for numberings         +*/
  PASTIX_INT                       vertnbr;              /*+ number of vertex in the graph     +*/
  EliminVertex   *          verttab;              /*+ Array of vertex                   +*/           
  PASTIX_INT            *          inbltab;              /*+ Array of in-blocs index           +*/
  PASTIX_INT            *          ownetab;              /*+ Array of cbl owner bloc           +*/           
} EliminGraph;

/*
**  The function prototypes.
*/

#ifndef STRUCT_ELIMINTREE
#define static
#endif
PASTIX_INT                         egraphInit        (EliminGraph *);
void                        egraphExit        (EliminGraph *);
PASTIX_INT                         egraphLoad        (EliminGraph *, FILE *);
PASTIX_INT                         egraphSave        (EliminGraph *, FILE *);
PASTIX_INT                         treeInit          (EliminTree *);
void                        treeExit          (EliminTree *);
PASTIX_INT                         treeLoad          (EliminTree *, FILE *);
PASTIX_INT                         treeSave          (EliminTree *, FILE *);
void                        treePlot          (EliminTree *, FILE *);

#undef static

#define TFATHER(treenode, r)     treenode->nodetab[treenode->nodetab[r].fathnum]
/* return the i_th sons num of the node n */ 
#define TSON(treenode, n, r)     treenode->sonstab[ treenode->nodetab[n].fsonnum + r ]
#define ROOT(treenode)           treenode->nodenbr-1
