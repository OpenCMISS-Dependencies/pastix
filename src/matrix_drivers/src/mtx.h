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

#define MTX_ISSYM(a) (a)[1]=='S'
#define MTX_ISHER(a) (a)[1]=='H'
#define MTX_ISCOM(a) (a)[0]=='C'
#define MTX_ISRHX(a) (a)[2]=='X'
#define MTX_ISRHS(a) (a)[0]!='\0'


void MatrixMarketRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType);
void mumpsReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type);
void  mumpsRead(char const *dirname, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType);
void mtxSym(pastix_int_t Nrow, pastix_int_t Ncol, pastix_int_t Nnzero, pastix_int_t *col, pastix_int_t *row, pastix_float_t *val);
void mtxCheck(pastix_int_t Nrow, pastix_int_t Ncol, pastix_int_t Nnzero, pastix_int_t *col, pastix_int_t *row, pastix_float_t *val);
void mtxReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type);
void mtxRead(char *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType, pastix_int_t flagsort);
void ijvReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type);
void ijvRead(char *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType, pastix_int_t flagsort);
void rsaReadHeader(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type, char *RhsType);
void rsaRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType);
#ifdef PREC_DOUBLE
void HBRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType);
#endif
void chbRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType, pastix_float_t **rhs);
void cccReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type);
void cccRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType);
void olafReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type);
void olafRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType, pastix_float_t **rhs);
void peerRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col,pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType, pastix_float_t **rhs);

void diag_dominance(pastix_int_t baseval, pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a);
void diag_unite(pastix_int_t baseval, pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a);
void no_diag(pastix_int_t baseval, pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja);
void symmetrize_pattern(pastix_int_t n, pastix_int_t **iaptr, pastix_int_t **japtr, pastix_float_t **aaptr);
void dimsym(pastix_int_t n, pastix_int_t **iaptr, pastix_int_t **japtr);
void checkStrucSym(pastix_int_t n, pastix_int_t *nz, pastix_int_t **colptr, pastix_int_t **row, pastix_float_t **avals);

void driverFdupros(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val,
		   pastix_float_t ** rhs, char **Type, char **RhsType);

#include "csparse.h"
