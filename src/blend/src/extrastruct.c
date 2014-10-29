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
#include <stdlib.h>
#include <stdio.h>

#include "common_pastix.h"
#include "cost.h"
#include "symbol.h"
#include "extrastruct.h"




PASTIX_INT extrasymbolInit(ExtraSymbolMatrix *extrasymb)
{
  extrasymb->baseval = 0;
  extrasymb->addcblk = 0;
  extrasymb->addblok = 0;
  extrasymb->sptcblk = NULL;
  extrasymb->sptblok = NULL;
  extrasymb->sptblnb = NULL;
  extrasymb->sptcbnb = NULL;
  extrasymb->subtreeblnbr = NULL;
  extrasymb->curcblk = 0;
  extrasymb->sizcblk = 0;
  extrasymb->curblok = 0;
  extrasymb->sizblok = 0;
  extrasymb->cblktab = NULL;
  extrasymb->bloktab = NULL;
  return 1;
}

void extrasymbolExit(ExtraSymbolMatrix *extrasymb)
{
  memFree_null(extrasymb->sptcblk);
  memFree_null(extrasymb->sptcbnb);
  memFree_null(extrasymb->sptblok);
  memFree_null(extrasymb->sptblnb);
  memFree_null(extrasymb->subtreeblnbr);
  if(extrasymb->sizcblk > 0)
    memFree_null(extrasymb->cblktab);
  if(extrasymb->sizblok > 0)
    memFree_null(extrasymb->bloktab);
  memFree_null(extrasymb);
}

PASTIX_INT extracostInit(ExtraCostMatrix *extracost)
{
    extracost->cblktab = NULL;
    extracost->bloktab = NULL;
    return 1;
}

void extracostExit(ExtraCostMatrix *extracost)
{
    if(extracost->cblktab != NULL)
	memFree_null(extracost->cblktab);
    if(extracost->bloktab != NULL)
	memFree_null(extracost->bloktab);
    memFree_null(extracost);
}


void extra_inc_blok(ExtraSymbolMatrix *extrasymb, ExtraCostMatrix *extracost)
{
    extrasymb->curblok++;
    
    /** if extra-blokktab is not big enough, make it bigger !! **/
    if(extrasymb->curblok + 1 >= extrasymb->sizblok)
	{
	    SymbolBlok *tmp;
	    CostBlok   *tmp2;
	    tmp2 = extracost->bloktab;
	    tmp  = extrasymb->bloktab;
	    /* add memory space to extra symbol matrix */ 
	    MALLOC_INTERN(extrasymb->bloktab, 
			  extrasymb->sizblok + extrasymb->sizblok/2 + 1, 
			  SymbolBlok);
	    memcpy(extrasymb->bloktab, tmp, sizeof(SymbolBlok)*extrasymb->curblok);
	    /* add memory space to extra cost matrix */
	    MALLOC_INTERN(extracost->bloktab, 
			  extrasymb->sizblok + extrasymb->sizblok/2 + 1, 
			  CostBlok);
	    /*fprintf(stderr, "Size %ld curblok %ld NewSize %ld \n", extrasymb->sizblok, extrasymb->curblok, (extrasymb->sizblok + extrasymb->sizblok/2 +1));
	    ASSERT( extracost->bloktab != NULL,MOD_BLEND);*/
	    memCpy(extracost->bloktab, tmp2, sizeof(CostBlok)*extrasymb->curblok);

	    extrasymb->sizblok = extrasymb->sizblok + extrasymb->sizblok/2 + 1;
	    memFree_null(tmp);
	    memFree_null(tmp2);
	}
}

void extra_inc_cblk(ExtraSymbolMatrix *extrasymb, ExtraCostMatrix *extracost)
{
    extrasymb->curcblk++;
    /** if extra-cblktab is not big enough, make it bigger !! **/
    if(extrasymb->curcblk + 1 >= extrasymb->sizcblk)
	{
	    SymbolCblk *tmp;
	    CostCblk   *tmp2;
	    tmp  = extrasymb->cblktab;
	    tmp2 = extracost->cblktab;
	    /* add memory space to extra symbol matrix */ 
	    MALLOC_INTERN(extrasymb->cblktab,
			  extrasymb->sizcblk + extrasymb->sizcblk/5 + 1,
			  SymbolCblk);
	    memcpy(extrasymb->cblktab, tmp, sizeof(SymbolCblk)*extrasymb->curcblk);
	    /* add memory space to extra cost matrix */
	    MALLOC_INTERN(extracost->cblktab, 
			  extrasymb->sizcblk + extrasymb->sizcblk/5 + 1,
			  CostCblk);
	    memcpy(extracost->cblktab, tmp2, sizeof(CostCblk)*extrasymb->curcblk);

	    extrasymb->sizcblk = extrasymb->sizcblk + extrasymb->sizcblk/5 + 1;
	    memFree_null(tmp);
	    memFree_null(tmp2);
	}
}
