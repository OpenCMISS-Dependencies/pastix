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
  File: debug_dump.h
  
  Functions to dump informations on disk.
*/

#ifndef DEBUG_DUMP_H
#define DEBUG_DUMP_H

/*
  Function: dump1
  
  Dumps ord->permtab on disk.

  Format: 
    > i -> permtab[i]

  parameters:
    ord    - Order structure to print permtab from.
    stream - FILE *, opened in write mode, in which permtab will be writen.
    colnbr - Number of elements in permtab.
 */
void dump1(Order *ord,
	   FILE  *stream, 
	   PASTIX_INT    colnbr);


/*
  Function: dump2

  Prints internal CSCd, in (i,j,v) format, in a file.

  Parameters:
    datacode - SolverMatrix.
    stream   - FILE * opened in write mode.
  
 */
void dump2(const SolverMatrix * datacode,
           CscMatrix          * cscmtx,
	   PASTIX_FLOAT              * trandcsc,
	   FILE               *stream);


/*
  Function: dump3

  Prints solver matrix informations, in (i,j,v) format, in a file.

  Parameters:
    datacode - SolverMatrix.
    stream   - FILE * opened in write mode.
*/
void dump3(const SolverMatrix *datacode, 
	   FILE               *stream);

/*
  Function: dump3_LU

  Prints solver matrix informations, in (i,j,v) format, in a file, 
  for LU decomposition.

  Parameters:
    datacode - SolverMatrix.
    streamL  - FILE * opened in write mode.
    streamU  - FILE * opened in write mode.
*/
void dump3_LU(const SolverMatrix * datacode, 
	      FILE               * streamL, 
	      FILE               * streamU);


/*
  Function: dump4
  
  Writes column blocks and blocs dimension in a file.
  
  Parameters:
    datacode - SolverMatrix containing informations about blocs
    stream   - FILE * opened in write mode.
*/
void dump4(const SolverMatrix *datacode, 
	   FILE               *stream);


/*
  Function: dump5

  Writes right-hand-side memeber in a file.

  Parameters:
    datacode - SolverMatrix containing right-hand-side member.
    stream   - FILE * opened in write mode.
*/
void dump5(const SolverMatrix *datacode, 
	   FILE               *stream);


/*
  Function: dump6

  Prints diagonal blocks in the folowing format :
  > ** block diag <cblknbr> **
  > <line1> [<value1> <value2> ... ]
  > <line2> [...                   ]
  
  Prints one file dor L and one for U.

  Parameters:
    datacode - SolverMatrix.
    streamL  - FILE * into which L diagonal blocs will be writen.
    streamU  - FILE * into which U diagonal blocs will be writen.
*/
void dump6(const SolverMatrix *datacode, 
	   FILE               *streamL, 
	   FILE               *streamU);


/*
  Function: dump7
  
  Writes a vector in the folowing format :
  > <line1> <value[line1]>
  > <line2> <value[line1]>

  Parameters: 
    v      - vector to write.
    stream - FILE * opened in write mode.
    nbr    - Size of the vector v.
*/
void dump7(PASTIX_FLOAT *v, 
	   FILE  *stream, 
	   PASTIX_INT    colnbr);

#endif /* DEBUG_DUMP_H */
