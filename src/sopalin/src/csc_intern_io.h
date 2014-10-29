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
  File: csc_intern_io.h

  Functions to save or load internal CSC in binary or ascii mode.
  
*/

#ifndef CSC_INTERN_IO_H
#define CSC_INTERN_IO_H
/*
  Function: CscSave

  Writes on disk an internal CSCd in text format.

  Format is :
  
  > CSC_FNBR(cscptr)
  > CSC_COLNBR(cscptr,iter)    ! iter = 0 to CSC_FNBR(cscptr) - 1
  > CSC_COL(cscptr,iter,iter2) ! iter2 = 0 to CSC_COLNBR(cscptr,iter)
  > ...
  > CSC_ROW(cscptr,iter) ! For all rows and values (iter)
  > CSC_VAL(cscptr,iter)

  Parameters :
    cscprt - the internal CSCd structure to save.
    stream - the FILE to write into, open in write mode.
*/
PASTIX_INT CscSave(const CscMatrix * const cscptr, 
	    FILE            * const stream);

/*
  Function: CscBSave

  Writes on disk an internal CSCd in binary format.

  Parameters :
    cscprt - the internal CSCd structure to save.
    stream - the FILE to write into, open in write mode.
*/
PASTIX_INT CscBSave(const CscMatrix * const cscptr, 
	     FILE            * const stream);

/* 
   Function: CscLoad

   Reads an internal CSCd from disk.

   Format is :
   
   > CSC_FNBR(cscptr)
   > CSC_COLNBR(cscptr,iter)    ! iter = 0 to CSC_FNBR(cscptr) - 1
   > CSC_COL(cscptr,iter,iter2) ! iter2 = 0 to CSC_COLNBR(cscptr,iter)
   > ...
   > CSC_ROW(cscptr,iter) ! For all rows and values (iter)
   > CSC_VAL(cscptr,iter)

   Parameters :
     cscprt - the internal CSCd structure to load.
     stream - the FILE to write into, open in read mode. 
*/
PASTIX_INT CscLoad(CscMatrix * cscptr, 
	    FILE      * stream);

/*
  Function: CscBLoad
  
  Loads an internal CSCd from a file saved in binary mode.

  Parameters :
    cscprt - the internal CSCd structure to load.
    stream - the FILE to write into, open in read mode. 
*/
PASTIX_INT CscBLoad(CscMatrix * cscptr, 
	     FILE      * stream);

#endif /* CSC_INTERN_IO_H */
