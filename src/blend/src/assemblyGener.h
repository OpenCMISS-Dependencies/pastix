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
/**   NAME       : assemblyGener.h                         **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Generation of the assembly strucuture   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     09 sep 1998     **/
/**                                                        **/
/************************************************************/
#ifndef ASSEMBLY
#define static
#endif

#ifdef OLD_ASSEMBLY
void       assemblyGener(Assembly1D *, PASTIX_INT, const SymbolMatrix *, const PASTIX_INT *);
#else
void assemblyGener(PASTIX_INT clustnum, Assembly1D *assemb1D, Assembly2D *assemb2D, PASTIX_INT clustnbr, const SymbolMatrix *symbmtx, const PASTIX_INT *blprtab, BlendCtrl *ctrl, const Dof * const dofptr);
static void symbolGener(PASTIX_INT clustnum, const PASTIX_INT *cblklocalnum1D, PASTIX_INT bloknbr1D, PASTIX_INT cblknbr1D, const PASTIX_INT *cbprtab, const SymbolMatrix *symbmtx, SymbolMatrix *symb1D, const Dof * const dofptr);
#endif

#undef static
