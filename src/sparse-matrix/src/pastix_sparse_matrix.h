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
#ifndef PASTIX_SPARSE_MATRIX_H
#define PASTIX_SPARSE_MATRIX_H

/*+ Enumerate the possible type of floating point values +*/
enum pastix_coef_type_ {
  PASTIX_SIMPLE,          /*+ Float value          +*/
  PASTIX_DOUBLE,          /*+ Double value         +*/
  PASTIX_SIMPLE_COMPLEX,  /*+ Complex value        +*/
  PASTIX_DOUBLE_COMPLEX   /*+ Double complex value +*/
};

typedef enum   pastix_coef_type_     pastix_coef_type_t;
typedef struct pastix_cblk_          pastix_cblk_t;
typedef struct pastix_block_         pastix_block_t;
typedef struct pastix_sparse_matrix_ pastix_sparse_matrix_t;

/*+
 Free a pastix_sparse_matrix_t.
 +*/
int
pastix_sparse_matrix_free(pastix_sparse_matrix_t * sparsemtx);

/*+
 Get a pastix_sparse_matrix_t
 +*/
int
pastix_get_sparse_matrix(pastix_data_t          * pastix_data,
                         pastix_coef_type_t       coeftype,
                         pastix_sparse_matrix_t * sparsemtx);

/*+ The block structure. +*/
struct pastix_block_ {
  PASTIX_INT                   frownum;   /*+ First row index                  +*/
  PASTIX_INT                   lrownum;   /*+ Last row index (inclusive)       +*/
  PASTIX_INT                   levfval;   /*+ Level-of-fill value              +*/
  PASTIX_INT                   coefind;   /*+ Index of the first coefficient of
                                     the bloc in the data array.        +*/
  pastix_cblk_t       * cblk;      /*+ Pointer to the column block the
                                     block belongs to.                  +*/
  pastix_cblk_t       * fcblk;     /*+ Pointer to the column block the
                                     block if facing.                   +*/
};

/*+ The column block structure. +*/
struct pastix_cblk_ {
  PASTIX_INT                      fcolnum; /*+ First column index                   +*/
  PASTIX_INT                      lcolnum; /*+ Last column index (inclusive)        +*/
  PASTIX_INT                      bloknbr; /*+ Number of blocks in the cblk.        +*/
  pastix_block_t *         bloktab; /*+ List of blocks of the cblk.          +*/
  PASTIX_INT                      fblknbr; /*+ Number of blocks contributing to
                                      the column block.                      +*/
  pastix_block_t **        fblktab; /*+ List of pointers to the blocks
                                      contributing to the column block.      +*/
  void                   * lvalues; /*+ Values in the L part of the cblk.    +*/
  void                   * uvalues; /*+ Values in the U part of the cblk.    +*/
  void                   * lschur;  /*+ Values in the L Schur part of the
                                      cblk.                                  +*/
  void                   * uschur;  /*+ Values in the U Schur part of the
                                      cblk.                                  +*/
  PASTIX_INT                      stride;  /*+ Number of rows in th cblk            +*/
  PASTIX_INT                      color;   /*+ Color of column block (PICL trace)   +*/
  PASTIX_INT                      procdiag;/*+ Processor owner of diagonal block    +*/
  PASTIX_INT                      cblkdiag;/*+ Column block owner of diagonal block +*/
};


/*+ The sparse column block matrix +*/
struct pastix_sparse_matrix_ {
  PASTIX_INT                   baseval;      /*+ Base value for numberings         +*/
  pastix_coef_type_t    coeftype;     /*+ Type of floats                    +*/
  PASTIX_INT                   coefnbr;      /*+ Total number of entries.          +*/
  PASTIX_INT                   nodenbr;      /*+ Rank of the matrix.               +*/
  PASTIX_INT                   cblknbr;      /*+ Number of column blocks           +*/
  PASTIX_INT                   bloknbr;      /*+ Number of blocks                  +*/
  pastix_cblk_t       * cblktab;      /*+ Array of column blocks [+1,based] +*/
  pastix_block_t      * bloktab;      /*+ Array of blocks [based]           +*/
  PASTIX_INT                   maxcblksize;  /*+ Maximum size of a column block.   +*/
};

#endif /* PASTIX_SPARSE_MATRIX_H */
