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
  File: tools.c

  Some tools used in sopalin.

 */
#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#include "common_pastix.h"
#include "sopalin_compute.h"
#include "tools.h"

/*
   Function: dim_dgeam


   Computes b = alpha * a + b.

     - if transa and transb equals 'N'.
        > b := alpha * a + b
     - if transa = 'T' and transb ='N'.
        > b := alpha * trans(a) + b
     - if transa = 'N' and transb ='T'.
        > trans(b) := alpha * a + trans(b)
     - if transa = 'T' and transb ='T'.
        > trans(b) := alpha * trans(a) + trans(b)

   Parameters:
     transa - indicates if a needs to be transposed.
     transb - indicates if b needs to be transposed.
     m      - number of row in a and b.
     n      - number of colonnes in a and b.
     alpha  - scalar.
     a      - Matrix a.
     lda    - Stride between 2 columns of a.
     b      - Matrix b.
     ldb    - Stride between 2 columns of b.
*/
void dim_dgeam(char *transa,char *transb,PASTIX_INT m,PASTIX_INT n,PASTIX_FLOAT alpha,PASTIX_FLOAT *a,
               PASTIX_INT lda,PASTIX_FLOAT *b,PASTIX_INT ldb)
{
  PASTIX_INT i;

  if (*transa=='N')
    {
      if (*transb=='N')
        for (i=0;i<n;i++)
          {
#ifdef CPLX
          {SOPALIN_AXPY(m,alpha,(PASTIX_FLOAT *) &(a[i*lda]),1, (PASTIX_FLOAT *) &(b[i*ldb]),1);}
#else
          {SOPALIN_AXPY(m,alpha,&(a[i*lda]),1,&(b[i*ldb]),1);}
#endif
          }
      else
        for (i=0;i<n;i++)
          {
#ifdef CPLX
          {SOPALIN_AXPY(m,alpha,(PASTIX_FLOAT *) &(a[i*lda]),1, (PASTIX_FLOAT *) &(b[i]),ldb);}
#else
          {SOPALIN_AXPY(m,alpha,&(a[i*lda]),1,&(b[i]),ldb);}
#endif
          }
    }
  else
    {
      if (*transb=='N')
        for (i=0;i<n;i++)
          {
#ifdef CPLX
          {SOPALIN_AXPY(m,alpha,(PASTIX_FLOAT *) &(a[i]),lda, (PASTIX_FLOAT *) &(b[i*ldb]),1);}
#else
          {SOPALIN_AXPY(m,alpha,&(a[i]),lda,&(b[i*ldb]),1);}
#endif
          }
      else
        for (i=0;i<n;i++)
          {
#ifdef CPLX
          {SOPALIN_AXPY(m,alpha,(PASTIX_FLOAT *) &(a[i]),lda, (PASTIX_FLOAT *) &(b[i]),ldb);}
#else
          {SOPALIN_AXPY(m,alpha,&(a[i]),lda,&(b[i]),ldb);}
#endif
          }
    }
}

/*
   Function: GetMpiType

   Construct a MPI type to store complex values.

   Parameters:
     none

   Return:

   the constructed MPI type for complex floating values.
 */
MPI_Datatype GetMpiType(void)
{
  static int first=1;
  static MPI_Datatype floattype;

  if (first==1)
    {
#ifdef FORCE_NOMPI
#ifdef PREC_DOUBLE
      floattype=MPI_DOUBLE_COMPLEX;
#else
      floattype=MPI_COMPLEX;
#endif
#else
#ifdef PREC_DOUBLE
#ifdef MPI_DOUBLE_COMPLEX
      floattype=MPI_DOUBLE_COMPLEX;
#else
      MPI_Type_contiguous(2,MPI_DOUBLE,&floattype);
      MPI_Type_commit(&floattype);
#endif
#else
#ifdef MPI_COMPLEX
      floattype=MPI_COMPLEX;
#else
      MPI_Type_contiguous(2,MPI_FLOAT,&floattype);
      MPI_Type_commit(&floattype);
#endif
#endif
#endif
      first=0;
    }

  return floattype;
}

/*
   Function: FreeMpiType

   Free the MPI type to store complex values.

   Parameters:
     none

 */
void FreeMpiType(void)
{

#ifndef FORCE_NOMPI
#if ((defined PREC_DOUBLE && (!defined MPI_DOUBLE_COMPLEX)) || (!defined MPI_COMPLEX))
  MPI_Datatype floattype = GetMpiType();
  MPI_Type_free(&(floattype));
#endif
#endif
}

/*
  Function: mysum

  Computes the sum of *in* and *inout* vectors
  and stores it in *inout*.


  Parameters:
    in    - First vector.
    inout - Second vector wich will store the sum.
    len   - Size of each vector.
    dptr  - MPI datatype
 */
void mysum(void *in, void *inout, int *len, MPI_Datatype *dptr)
{
  PASTIX_INT i;
  PASTIX_FLOAT *a = (PASTIX_FLOAT *) in;
  PASTIX_FLOAT *b = (PASTIX_FLOAT *) inout;
  (void)dptr;

  for (i=0; i<*len; ++i)
    {
      b[i] += a[i];
    }
}

/*
  Function: GetMpiSum

  Creates a new MPI sum operation.

  Parameters:
    none

  Return: the new operation created.

 */
MPI_Op GetMpiSum(void)
{
  static MPI_Op myop;
#if ((defined PREC_DOUBLE && (!defined MPI_DOUBLE_COMPLEX)) || (!defined MPI_COMPLEX))
  static int first=1;
  if (first==1)
    {
      MPI_Op_create(&mysum,1,&myop);
      first=0;
    }
#else
  myop = MPI_SUM;
#endif
  return myop;
}

/*
  Function: FreeMpiSum

  Free the MPI sum operation.

  Parameters:
    none

 */
void FreeMpiSum(void)
{
#if ((defined PREC_DOUBLE && (!defined MPI_DOUBLE_COMPLEX)) || (!defined MPI_COMPLEX))
  MPI_Op myop = GetMpiSum();
  MPI_Op_free(&(myop));
#endif
}
