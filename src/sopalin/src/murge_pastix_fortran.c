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
#include "murge_pastix.h"

FORTRAN_NAME(MURGE_ANALYZE,
             murge_analyze,
             (INTS *id,  INTS *ierror),
             (id,  ierror))
{
  *ierror = MURGE_Analyze(*id);
}


FORTRAN_NAME(MURGE_FACTORIZE,
             murge_factorize,
             (INTS *id,  INTS *ierror),
             (id,  ierror))
{
  *ierror = MURGE_Factorize(*id);
}

FORTRAN_NAME(MURGE_SETORDERING,
             murge_setordering,
             (INTS *id, INTS * permutation, INTS * ierror),
             (id, permutation, ierror))
{
  *ierror = MURGE_SetOrdering(*id, permutation);
}

FORTRAN_NAME(MURGE_PRODUCTSETLOCALNODENBR,
             murge_productsetlocalnodenbr,
             (INTS *id, INTS *n, INTS * ierror),
             (id, n, ierror))
{
  *ierror =  MURGE_ProductSetLocalNodeNbr(*id, *n);
}

FORTRAN_NAME(MURGE_PRODUCTSETGLOBALNODENBR,
             murge_productsetglobalnodenbr,
             (INTS *id, INTS *n, INTS * ierror),
             (id, n, ierror))
{
  *ierror =  MURGE_ProductSetGlobalNodeNbr(*id, *n);
}

FORTRAN_NAME(MURGE_PRODUCTSETLOCALNODELIST,
             murge_productsetlocalnodelist,
             (INTS *id, INTS *l2g, INTS * ierror),
             (id, l2g, ierror))
{
  *ierror =  MURGE_ProductSetLocalNodeList(*id, l2g);
}
FORTRAN_NAME(MURGE_GETLOCALPRODUCT,
             murge_getlocalproduct,
             (INTS *id,  COEF * x, INTS * ierror),
             (id, x, ierror))
{
  *ierror = MURGE_GetLocalProduct(*id, x);
}

FORTRAN_NAME(MURGE_GETGLOBALPRODUCT,
             murge_getglobalproduct,
             (INTS *id,  COEF * x, INTS *root, INTS * ierror),
             (id, x, root, ierror))
{
  *ierror = MURGE_GetGlobalProduct(*id, x, *root);
}

FORTRAN_NAME(MURGE_FORCENOFACTO,
             murge_forcenofacto,
             (INTS *id, INTS * ierror),
             (id, ierror))
{
  *ierror =  MURGE_ForceNoFacto(*id);
}

FORTRAN_NAME(MURGE_SETLOCALNODELIST,
             murge_setlocalnodelist,
             (INTS *id, INTS * n, INTS * list, INTS * ierror),
             (id, n, list, ierror))
{
  *ierror =  MURGE_SetLocalNodeList(*id, *n, list);
}


FORTRAN_NAME(MURGE_ASSEMBLYSETSEQUENCE,
             murge_assemblysetsequence,
             (INTS *id, INTL *coefnbr, INTS * ROWs, INTS * COLs,
              INTS *op, INTS *op2, INTS *mode, INTS *nodes,
              INTS * id_seq, INTS *ierror),
             (id, coefnbr, ROWs, COLs, op, op2, mode, nodes, id_seq, ierror))
{
  *ierror = MURGE_AssemblySetSequence(*id, *coefnbr, ROWs, COLs,
                                      *op, *op2, *mode, *nodes, id_seq);
}

FORTRAN_NAME(MURGE_ASSEMBLYUSESEQUENCE,
             murge_assemblyusesequence,
             (INTS *id, INTS *id_seq, COEF * values, INTS * ierror),
             (id, id_seq, values, ierror))
{
  *ierror = MURGE_AssemblyUseSequence(*id, *id_seq, values);
}

FORTRAN_NAME(MURGE_ASSEMBLYDELETESEQUENCE,
             murge_assemblydeletesequence,
             (INTS *id, INTS *id_seq, INTS *ierror),
             (id, id_seq, ierror))
{
  *ierror = MURGE_AssemblyDeleteSequence(*id, *id_seq);
}

FORTRAN_NAME(MURGE_SETDROPNODES,
	     murge_setdropnodes,
	     (INTS * id, INTS * nodenbr, INTS * dropmask, INTS * ierror),
	     (id, nodenbr, dropmask, ierror))
{
  *ierror = MURGE_SetDropNodes(*id, *nodenbr, dropmask);
}

FORTRAN_NAME(MURGE_SETDROPCOLS,
	     murge_setdropcols,
	     (INTS * id, INTS * nodenbr, INTS * dropmask, INTS * ierror),
	     (id, nodenbr, dropmask, ierror))
{
  *ierror = MURGE_SetDropCols(*id, *nodenbr, dropmask);
}

FORTRAN_NAME(MURGE_SETDROPROWS,
	     murge_setdroprows,
	     (INTS * id, INTS * nodenbr, INTS * dropmask, INTS * ierror),
	     (id, nodenbr, dropmask, ierror))
{
  *ierror = MURGE_SetDropRows(*id, *nodenbr, dropmask);
}

FORTRAN_NAME(MURGE_COLGETNONZEROSNBR,
	     murge_colgetnonzerosnbr,
	     (INTS * id, INTS * COL, INTS * nnzNbr,INTS * ierror),
	     (id, COL, nnzNbr, ierror))
{
  *ierror = MURGE_ColGetNonZerosNbr(*id, *COL, nnzNbr);
}


FORTRAN_NAME(MURGE_COLGETNONZEROSIDX,
	     murge_colgetnonzerosidx,
	     (INTS * id, INTS * COL, INTS * nnzIdx,INTS * ierror),
	     (id, COL, nnzIdx, ierror))
{
  *ierror = MURGE_ColGetNonZerosIdx(*id, *COL, nnzIdx);
}
