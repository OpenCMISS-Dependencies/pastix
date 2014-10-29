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
/**   NAME       : assembly.h                              **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the assembly process.               **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

#define ASSEMBLY_H

#ifdef CXREF_DOC
#include "common_pastix.h"
#endif /* CXREF_DOC */

/*
**  The type and structure definitions.
*/

/*+ 1D Assembly structure based on the 1D distribution for assembly. +*/

typedef struct Assembly1D_ {
  PASTIX_INT *                     blprtab;              /*+ Array of block to processor mapping        +*/
  PASTIX_INT *                     nocbtab;              /*+ Array of node to column block (1D) mapping +*/
  PASTIX_INT *                     rnumtab;              /*+ Absolute number to relative number tab     +*/
} Assembly1D;


/*+ 2D Assembly structure. +*/
typedef struct Assembly2D_ {
  SymbolMatrix *            symbmtx;              /*+ Symbol matrix in 1D distribution needed for assembling the matrix   +*/ 
  PASTIX_INT *                     blok2proc_tab;        /*+ local block i in 1D --> processor owner in 2D distribution +*/
  PASTIX_INT *                     blok2cblk_tab;        /*+ local block i in 2D --> local cblk on the same processor in the 2D distribution  +*/
  PASTIX_INT *                     blok2blok_tab;        /*+ local block i in 1D --> local block i in the 2D distribution +*/
} Assembly2D;
