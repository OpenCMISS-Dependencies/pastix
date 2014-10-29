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
/*
 * File: petscread.h
 *
 * Read Matrix written by PETSc in binary format.
 *
 */

/*
 *  Function: PETScRead
 *
 *  Reads a matrix in a binary PETScformat/
 *
 * Parameters:
 *   dirname - Path to the directory containing matrix
 *   Ncol    - Number of columns
 *   Nrow    - Number of rows
 *   Nnzero  - Number of non zeros
 *   col     - Index of first element of each column in *row* and *val*
 *   row     - Row of eah element
 *   val     - Value of each element
 *   Type    - Type of the matrix
 *   RhsType - Type of the right-hand-side.
 *
 */
void PETScRead(char const      *filename,
               pastix_int_t    *Ncol,
               pastix_int_t    *Nrow,
               pastix_int_t    *Nnzero,
               pastix_int_t   **col,
               pastix_int_t   **row,
               pastix_float_t **val,
               char           **Type,
               char           **RhsType);
