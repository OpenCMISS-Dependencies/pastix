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
/**   NAME       : splitpart.h                             **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                repartition and make processor          **/
/**                candidate groups                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     09 sep 1998     **/
/**                                                        **/
/************************************************************/
#ifndef SPLITPART_H

#define CLUSTER 1
#define NOCLUSTER 0


#define static

#endif

void  splitPart     (SymbolMatrix *, BlendCtrl *, const Dof *);
PASTIX_INT   check_candidat(SymbolMatrix *, BlendCtrl *);

void        setTreeLevel         (Cand *, const EliminTree *);
void        setTreeCostLevel     (Cand *, const EliminTree *, const CostMatrix *);
static void setSubtreeLevel      (PASTIX_INT, Cand *, const EliminTree *);
static void setSubtreeCostLevel  (PASTIX_INT, Cand *, const EliminTree *, const CostMatrix *);
static void setDistribType       (const PASTIX_INT, SymbolMatrix *, Cand *, const PASTIX_INT);
static void setSubtreeDistribType(const SymbolMatrix *, const CostMatrix *, PASTIX_INT , const BlendCtrl *, PASTIX_INT);

static void splitOnProcs    (SymbolMatrix *, ExtraSymbolMatrix *, ExtraCostMatrix *, BlendCtrl *, 
			     const Dof *, PASTIX_INT, PASTIX_INT);
static void splitCblk       (SymbolMatrix *, ExtraSymbolMatrix *, ExtraCostMatrix *, BlendCtrl *, 
			     const Dof *, PASTIX_INT, PASTIX_INT, PASTIX_INT *);

static void printTree          (FILE*, const EliminTree *, PASTIX_INT);
static void propMappTree       (SymbolMatrix *, ExtraSymbolMatrix *, ExtraCostMatrix *, BlendCtrl *, const Dof *);
static void propMappSubtree    (SymbolMatrix *, ExtraSymbolMatrix *, ExtraCostMatrix *, BlendCtrl *, const Dof *,
				PASTIX_INT, PASTIX_INT, PASTIX_INT, PASTIX_INT, double *);
static void propMappSubtreeNC  (SymbolMatrix *, ExtraSymbolMatrix *, ExtraCostMatrix *, BlendCtrl *, const Dof *,
				PASTIX_INT, PASTIX_INT, PASTIX_INT, PASTIX_INT, double *);
static void propMappSubtreeOn1P(SymbolMatrix *, ExtraSymbolMatrix *, ExtraCostMatrix *, BlendCtrl *, const Dof *,
				PASTIX_INT, PASTIX_INT, PASTIX_INT, PASTIX_INT);

static void propMappTreeNoSplit    (SymbolMatrix *, BlendCtrl *, const Dof *);
static void propMappSubtreeNoSplit (SymbolMatrix *, BlendCtrl *, const Dof *, PASTIX_INT, PASTIX_INT, PASTIX_INT, double *);


static double maxProcCost     (double *, PASTIX_INT);
static void   subtreeSetCand  (PASTIX_INT, PASTIX_INT, BlendCtrl *, double);
static double blokUpdateCost  (PASTIX_INT, PASTIX_INT, CostMatrix *, ExtraCostMatrix *, const SymbolMatrix *, 
			       const ExtraSymbolMatrix *, BlendCtrl *, const Dof *);

static PASTIX_INT    countBlok            (PASTIX_INT, SymbolMatrix *, PASTIX_INT);
static PASTIX_INT    setSubtreeBlokNbr    (PASTIX_INT, const EliminTree *, SymbolMatrix *, ExtraSymbolMatrix *, PASTIX_INT);
static void   clusterCandCorrect   (PASTIX_INT, Cand *, const EliminTree *, BlendCtrl *);
static void   setClusterCand       (PASTIX_INT, Cand *, const EliminTree *, PASTIX_INT, PASTIX_INT);

#undef static

