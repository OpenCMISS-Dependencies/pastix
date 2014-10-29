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
/**   NAME       : extrastruct.h                           **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Structures used to add extra cblk       **/
/**                when a cblk is splitted.                **/
/**                The goal is to avoid a big amount of    **/
/**                dynamic memory allocations and          **/
/**                memory copies                           **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/
typedef struct ExtraCostMatrix_ {
  CostCblk   *              cblktab;
  CostBlok   *              bloktab;
} ExtraCostMatrix;



typedef struct ExtraSymbolMatrix_ {
  PASTIX_INT                       baseval;              /*+ Base value for numberings                    +*/
  PASTIX_INT                       addcblk;              /*+ Number of cblk created                       +*/
  PASTIX_INT                       addblok;              /*+ Number of blok created                       +*/
  PASTIX_INT        *              sptcblk;              /*+ Index for splitted cblk in the cblktab       +*/
  PASTIX_INT        *              sptcbnb;              /*+ Number of splitted cblk for a cblk           +*/
  PASTIX_INT        *              sptblok;              /*+ Index for splitted blok in the bloktab       +*/
  PASTIX_INT        *              sptblnb;              /*+ Number of splitted blok for a blok           +*/
  PASTIX_INT        *              subtreeblnbr;         /*+ Number of blok in the subtree                +*/
  PASTIX_INT                       curcblk;              /*+ Cursor for cblktab                           +*/
  PASTIX_INT                       sizcblk;              /*+ Size of allocated cblktab                    +*/
  PASTIX_INT                       curblok;              /*+ Cursor for bloktab                           +*/
  PASTIX_INT                       sizblok;              /*+ Size of allocated bloktab                    +*/
  SymbolCblk *              cblktab;              /*+ Array of column blocks [+1,based]            +*/
  SymbolBlok *              bloktab;              /*+ Array of blocks [based]                      +*/
} ExtraSymbolMatrix;

/*
**  The function prototypes.
*/

#ifndef EXTRASYMBOL
#define static
#endif
void                        extra_inc_cblk           (ExtraSymbolMatrix *, ExtraCostMatrix *);
void                        extra_inc_blok           (ExtraSymbolMatrix *, ExtraCostMatrix *);
PASTIX_INT                         extrasymbolInit          (ExtraSymbolMatrix *);
void                        extrasymbolExit          (ExtraSymbolMatrix *);
PASTIX_INT                         extrasymbolLoad          (ExtraSymbolMatrix *, FILE *);
PASTIX_INT                         extrasymbolSave          (ExtraSymbolMatrix *, FILE *);
PASTIX_INT                         extracostInit            (ExtraCostMatrix *);
void                        extracostExit            (ExtraCostMatrix *);
PASTIX_INT                         extracostLoad            (ExtraCostMatrix *, FILE *);
PASTIX_INT                         extracostSave            (ExtraCostMatrix *, FILE *);

#undef static


