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
/**   NAME       : extendVector.h                          **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Vector that can extend its size         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 27 jul 1998     **/
/**                                 to     03 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

typedef struct ExtendVectorINT_ {
    PASTIX_INT          vecsize;
    PASTIX_INT          eltnbr;          /*+ number of elements +*/
    PASTIX_INT    *     inttab;          /*+ array of PASTIX_INT       +*/
} ExtendVectorINT;


/*
**  The function prototypes.
*/

#ifndef EXTENDVECTOR
#define static
#endif

PASTIX_INT                     *extendint_Init    (ExtendVectorINT *, PASTIX_INT);
void                     extendint_Exit    (ExtendVectorINT *);
void                     extendint_Add     (ExtendVectorINT *, PASTIX_INT);
PASTIX_INT                      extendint_Size    (ExtendVectorINT *);
PASTIX_INT                      extendint_Read    (ExtendVectorINT *, PASTIX_INT);
void                     extendint_Clear   (ExtendVectorINT *);
void                     extendint_ToSize  (PASTIX_INT, ExtendVectorINT *);
void                     extendint_incr    (ExtendVectorINT *);
#undef static
