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
#include "extendVector.h"

PASTIX_INT * extendint_Init(ExtendVectorINT *vec, PASTIX_INT size)
{
    vec->vecsize = size;
    vec->eltnbr  = 0;
    vec->inttab  = NULL;
    MALLOC_INTERN(vec->inttab, size, PASTIX_INT);
    return vec->inttab;
}

    
void extendint_Exit(ExtendVectorINT *vec)
{
    if(vec->inttab != NULL)
	memFree_null(vec->inttab);
    /*memFree_null(vec);*/
}

void extendint_Add(ExtendVectorINT *vec, PASTIX_INT elt)
{
    vec->inttab[vec->eltnbr] = elt;
    extendint_incr(vec);
}

PASTIX_INT extendint_Size(ExtendVectorINT *vec)
{
  return vec->eltnbr;
}

PASTIX_INT extendint_Read(ExtendVectorINT *vec, PASTIX_INT eltnum)
{
  ASSERT(eltnum <= vec->eltnbr,MOD_BLEND);
  return vec->inttab[eltnum];
}



void extendint_ToSize(PASTIX_INT size, ExtendVectorINT *vec)
{
    extendint_Clear(vec);

    if(size <= vec->vecsize)  /* there 's enough space */
	return;
    

    if(vec->inttab != NULL)   
	memFree_null(vec->inttab);

    MALLOC_INTERN(vec->inttab, size, PASTIX_INT);
    vec->vecsize = size;
}
    
void extendint_incr(ExtendVectorINT *vec)
{
    vec->eltnbr++;
    /** if the vector is not big enough, make it bigger !! **/
    if(!(vec->eltnbr < vec->vecsize))
	{
	    PASTIX_INT *tmp;
	    tmp = vec->inttab;
	    /* add memory space */
	    MALLOC_INTERN(vec->inttab, vec->vecsize + vec->vecsize/2 +1, PASTIX_INT);
	    memcpy(vec->inttab, tmp, sizeof(PASTIX_INT)*vec->eltnbr);
	    vec->vecsize = vec->vecsize + vec->vecsize/2 +1;
	    memFree_null(tmp);
	}
}

void extendint_Clear(ExtendVectorINT *vec)
{
    vec->eltnbr = 0;
}
