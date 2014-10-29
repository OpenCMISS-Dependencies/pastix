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
/*+********************************************************+*/
/*+                                                        +*/
/*+   NAME       : symbol_cost.h                           +*/
/*+                                                        +*/
/*+   AUTHORS    : Pascal HENON                            +*/
/*+                                                        +*/
/*+   FUNCTION   : Part of a parallel direct block solver. +*/
/*+                Cost compute functions                  +*/
/*+                                                        +*/
/*+   DATES      : # Version 0.0  : from : 01 jan 1999     +*/
/*+                                 to     20 jan 1999     +*/
/*+                                                        +*/
/*+********************************************************+*/

#ifndef SYMBOL_COST_H
#define static
#endif

void           symbCost      (PASTIX_INT *, double *, const SymbolMatrix *, const Dof *);
double         recursive_sum (PASTIX_INT, PASTIX_INT, double (*fval)(PASTIX_INT, const SymbolMatrix *, const Dof *), 
			      const SymbolMatrix *, const Dof *);
double         cholesky      (PASTIX_INT, const SymbolMatrix *, const Dof *);
double         crout_hyb     (PASTIX_INT, const SymbolMatrix *, const Dof *);
double         crout_2t      (PASTIX_INT, const SymbolMatrix *, const Dof *);
double         crout_3t      (PASTIX_INT, const SymbolMatrix *, const Dof *);
double         nnz           (PASTIX_INT, const SymbolMatrix *, const Dof *);
double         crout_blok    (PASTIX_INT, const SymbolMatrix *, const Dof *);

#undef static
