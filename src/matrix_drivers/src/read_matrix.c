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
 * File: read_matrix.c
 *
 * Definition of a global function to read all type of matrices.
 *
 */
#include <stdio.h>
#include <math.h>

#ifdef FORCE_NOMPI
#include "pastix_nompi.h"
#else
#include <mpi.h>
#endif
#include <stdlib.h>
#include <sys/types.h>
#include <stdint.h>

#ifdef TYPE_COMPLEX
#if (defined X_ARCHalpha_compaq_osf1)
#ifndef USE_CXX
#ifndef   _RWSTD_HEADER_REQUIRES_HPP
#include <complex>
#else  /* _RWSTD_HEADER_REQUIRES_HPP */
#include <complex.hpp>
#endif /* _RWSTD_HEADER_REQUIRES_HPP */
#endif /* USE_CXX */

#else  /* X_ARCHalpha_compaq_osf1 */

#include <complex.h>
#endif /* X_ARCHalpha_compaq_osf1 */
#endif /* TYPE_COMPLEX */

#ifdef X_ARCHsun
#include <inttypes.h>
#endif

#include "pastix.h"
#include "cscd_utils.h"
#include "string.h"

#include "common_drivers.h"
#include "read_matrix.h"
#include "rsaread.h"
#include "hbread.h"
#include "mmread.h"
#include "mmdread.h"
#include "cccread.h"
#include "olafread.h"
#include "chbread.h"
#include "cscdread.h"
#include "peerread.h"
#include "threefilesread.h"
#include "laplacian.h"
#include "petscread.h"
#ifdef FDUPROS
#include "fdupread.h"
#endif
#ifdef __INTEL_COMPILER
/* remark #1419: external declaration in primary source file */
#pragma warning(disable:1419)
#endif

/*
 * Function: read_matrix_common
 *
 * Common part to <read_matrix> and <dread_matrix>.
 *
 * Reads a matrix from a file in the format given by
 * driver_type (see <driver_type_enum>).
 *
 * Parameters:
 *   filename    - Name of the file to read from.
 *   ncol        - Number of column in the matrix (output).
 *   colptr      - Indexes in rows and avals of first element of
 *                 each column of the matrix.(output)
 *   rows        - Row of each element of the matrix.
 *   loc2glob    - Local to global column number indirection.
 *   values      - Values of each element of the matrix.
 *   type        - type of the matrix.
 *   rhstype     - type of the right and side.
 *   driver_type - driver to use to read the matrix.
 *   pastix_comm - MPI communicator containing all processes wich
 *                 call read_matrix.
 */
int read_matrix_common(char            *filename,    pastix_int_t    *ncol,
                       pastix_int_t   **colptr,      pastix_int_t   **rows,
                       pastix_int_t   **loc2glob,
                       pastix_float_t **values,      pastix_float_t **rhs,
                       char           **type,        char           **rhstype,
                       driver_type_t    driver_type, MPI_Comm         pastix_comm)
{
  pastix_int_t   nrows;
  pastix_int_t   nnz;
  int mpid;

  *rhs = NULL;

  MPI_Comm_rank(pastix_comm,&mpid);

  if ( mpid ==0 ||
       driver_type == LAPLACIAN ||
       driver_type == CSCD ||
       driver_type == FDUP_DIST ||
       driver_type == MMD)
    {
      switch(driver_type)
        {
#ifndef USE_NOFORTRAN
        case RSA:
          printf("driver: RSA file: %s\n", filename);
          rsaRead(filename,
                  ncol, &nrows, &nnz,
                  colptr, rows, values,
                  type, rhstype);
          break;
#endif /* USE_NOFORTRAN */
        case CHB:
          printf("driver: CHB file: %s\n", filename);
          chbRead(filename,
                  ncol, &nrows, &nnz,
                  colptr, rows, values,
                  type, rhstype, rhs);
          break;
        case CCC:
          printf("driver: CCC file: %s\n", filename);
          cccRead(filename,
                  ncol, &nrows, &nnz,
                  colptr, rows, values,
                  type, rhstype);
          break;
        case RCC:
          printf("driver: RCC file: %s\n", filename);
          cccRead(filename,
                  ncol, &nrows, &nnz,
                  colptr, rows, values,
                  type, rhstype);
          (*type)[0]='R';
          break;
        case OLAF:
          printf("driver: OLAF file: %s\n", filename);
          olafRead(filename,
                   ncol, &nrows, &nnz,
                   colptr, rows, values,
                   type, rhstype,
                   rhs);
          break;
        case PEER:
          printf("driver: PEER file: %s\n", filename);
          peerRead(filename,
                   ncol, &nrows, &nnz,
                   colptr, rows, values,
                   type, rhstype,
                   rhs);
      /* peerRead2("rsaname", ncol, &nrows, &nnz, colptr, rows, values, type, rhstype, rhs); */
          break;
        case HB:
          printf("driver: HB file: %s\n", filename);
          HBRead(filename,
                 ncol, &nrows, &nnz,
                 colptr, rows, values,
                 type, rhstype);
          break;
        case THREEFILES:
          printf("driver: 3files file: %s\n", filename);
          threeFilesRead(filename,
                         ncol, &nrows, &nnz,
                         colptr, rows, values,
                         type, rhstype);
          break;
        case MM:
          printf("driver: MatrixMarket file: %s\n", filename);
          MatrixMarketRead(filename,
                           ncol, &nrows, &nnz,
                           colptr, rows, values,
                           type, rhstype);
          break;
        case MMD:
          printf("driver: DistributedMatrixMarket file: %s\n", filename);
          DistributedMatrixMarketRead(filename,
                                      ncol, &nrows, &nnz,
                                      colptr, rows, values, loc2glob,
                                      type, rhstype);
          break;
        case PETSCS:
        case PETSCU:
        case PETSCH:
          printf("driver: PETSc file: %s\n", filename);
          PETScRead(filename,
                    ncol, &nrows, &nnz,
                    colptr, rows, values,
                    type, rhstype);
          if (driver_type == PETSCS) (*type)[1] = 'S';
          if (driver_type == PETSCH) (*type)[1] = 'H';
          break;
        case CSCD:
          printf("driver CSCdt file: %s\n", filename);
          cscdRead(filename,
                   colptr, rows, loc2glob, values,
                   rhs, ncol, &nnz,
                   pastix_comm);

          *type = (char *) malloc(4*sizeof(char));
          sprintf(*type,"RSA");
          break;
        case LAPLACIAN:
          if (mpid == 0)
            printf("driver Laplacian\n");
          genlaplacian(*ncol, &nnz, colptr, rows, values, rhs, type, rhstype);
          return EXIT_SUCCESS;
          break;
#ifdef FDUPROS
        case FDUP:
          printf("driver: FDupros file: %s\n",filename);
          driverFdupros(filename,
                        ncol, &nrows, &nnz,
                        colptr, rows, values,
                        rhs,
                        type, rhstype);
          break;
        case FDUP_DIST:
          printf("driver: FDupros file: %s\n",filename);
          driverFdupros_dist(filename,
                             ncol, &nrows, &nnz,
                             colptr, rows, loc2glob, values,
                             rhs,
                             type, rhstype,
                             pastix_comm);
          free(*rhs);
          *rhs = NULL;
          break;
#endif
        default:
          printf("default driver\n");
          rsaRead("rsaname", ncol, &nrows, &nnz, colptr, rows, values, type, rhstype);
          break;
          /*
           printf("matrix driver unknown\n");
           EXIT(MOD_SOPALIN,NOTIMPLEMENTED_ERR);
           */
        }
    }
#ifndef TYPE_COMPLEX
  if (*type)
    if ((*type)[1] == 'H')
      (*type)[1] = 'S';
#endif

    /* read RHS file */
  if (mpid == 0)
    {
      FILE *file;
      char filename2[256];
      pastix_int_t i;
      double re;

      sprintf(filename2,"%s.rhs",filename);
      fprintf(stderr,"open RHS file : %s\n",filename2);
      file = fopen(filename2,"r");
      if (file==NULL)
        {
          fprintf(stderr,"cannot load %s\n", filename2);
        }
      else
        {
          *rhs = (pastix_float_t *) malloc((*ncol)*sizeof(pastix_float_t));
          for (i = 0; i < *ncol; i++)
            {
              (*rhs)[i] = 0.0;
              if (1 != fscanf(file,"%lg\n", &re))
                {
                  fprintf(stderr, "ERROR: reading rhs(%ld)\n", (long int)i);
                  exit(1);
                }
              (*rhs)[i] = (pastix_float_t)re;
            }
          fclose(file);
        }
    }

  if (*rhs == NULL  && ( mpid == 0 ||
                         driver_type == LAPLACIAN) &&
      driver_type != CSCD && driver_type != FDUP_DIST && driver_type != MMD)
    {
      *rhs = (pastix_float_t *) malloc((*ncol)*sizeof(pastix_float_t));
      pastix_int_t i,j;
      for (i = 0; i < *ncol; i++)
        (*rhs)[i] = 0.0;

#ifdef TYPE_COMPLEX
      fprintf(stdout, "Setting right-hand-side member such as X[i] = i + i*I\n");
#else
      fprintf(stdout, "Setting right-hand-side member such as X[i] = i\n");
#endif
      for (i = 0; i < *ncol; i++)
        {
          for (j = (*colptr)[i] -1; j < (*colptr)[i+1]-1; j++)
            {
#ifdef TYPE_COMPLEX
              (*rhs)[(*rows)[j]-1] += (i+1 + I*(i+1))*(*values)[j];
#else
              (*rhs)[(*rows)[j]-1] += (i+1)*(*values)[j];
#endif
              if (MTX_ISSYM((*type)) && i != (*rows)[j]-1)
                {
#ifdef TYPE_COMPLEX
                  (*rhs)[i] += ((*rows)[j] + (*rows)[j] * I)*(*values)[j];
#else
                  (*rhs)[i] += ((*rows)[j])*(*values)[j];
#endif
                }
            }
        }
    }
  if (driver_type == CSCD || driver_type == FDUP_DIST || driver_type == MMD)
    {
      pastix_int_t i,j;
      pastix_int_t N;
      pastix_int_t send[2], recv[2];
      int comm_size;
      MPI_Comm_size(pastix_comm,&comm_size);
      send[0] = *ncol;
      send[1] = 0;
      if (*ncol != 0 && *rhs == NULL)
        send[1] = 1;
      MPI_Allreduce(send, recv, 2, MPI_PASTIX_INT, MPI_SUM, pastix_comm);
      N = recv[0];
      if (recv[1] > 0)
        {
          pastix_float_t *RHS;
          pastix_float_t *RHS_recv;
          RHS  = (pastix_float_t *) malloc((N)*sizeof(pastix_float_t));

          for (i = 0; i < N; i++)
            (RHS)[i] = 0.0;
          if (mpid == 0)
            fprintf(stdout, "Setting RHS such as X[i] = i\n");
          for (i = 0; i < *ncol; i++)
            {
              for (j = (*colptr)[i] -1; j < (*colptr)[i+1]-1; j++)
                {
                  (RHS)[(*rows)[j]-1] += ((*loc2glob)[i])*(*values)[j];
                  if (MTX_ISSYM((*type)) && i != (*rows)[j]-1)
                    (RHS)[(*loc2glob)[i]-1] += ((*rows)[j])*(*values)[j];
                }
            }
          RHS_recv  = (pastix_float_t *) malloc((N)*sizeof(pastix_float_t));
          MPI_Allreduce(RHS, RHS_recv, N, MPI_PASTIX_FLOAT, MPI_SUM, pastix_comm);
          free(RHS);
          *rhs = (pastix_float_t *) malloc((*ncol)*sizeof(pastix_float_t));

          for (i = 0; i < *ncol; i++)
            (*rhs)[i] = RHS_recv[(*loc2glob)[i]-1];
        }
    }

  if ( driver_type != LAPLACIAN && driver_type != CSCD &&
       driver_type != FDUP_DIST && driver_type != MMD)
    {
      MPI_Bcast(ncol,1,MPI_PASTIX_INT,0,pastix_comm);
      MPI_Bcast(&nnz,1,MPI_PASTIX_INT,0,pastix_comm);

      if (mpid!=0)
        {
          *colptr = (pastix_int_t *)   malloc((*ncol+1)*sizeof(pastix_int_t));
          *rows   = (pastix_int_t *)   malloc(nnz*sizeof(pastix_int_t));
          *values = (pastix_float_t *) malloc(nnz*sizeof(pastix_float_t));
          *rhs    = (pastix_float_t *) malloc((*ncol)*sizeof(pastix_float_t));
          *type   = (char *)  malloc(4*sizeof(char));
        }

      MPI_Bcast(*colptr, *ncol+1, MPI_PASTIX_INT,   0, pastix_comm);
      MPI_Bcast(*rows,    nnz,    MPI_PASTIX_INT,   0, pastix_comm);
      MPI_Bcast(*values,  nnz,    MPI_PASTIX_FLOAT, 0, pastix_comm);
      MPI_Bcast(*rhs,    *ncol,   MPI_PASTIX_FLOAT, 0, pastix_comm);
      MPI_Bcast(*type,    4,      MPI_CHAR,         0, pastix_comm);
  }

  return EXIT_SUCCESS;
}

/*
 * Function: read_matrix
 *
 * Reads a matrix from a file in the format given by
 * driver_type (see <driver_type_enum>).
 *
 * Parameters:
 *   filename    - Name of the file to read from.
 *   ncol        - Number of column in the matrix (output).
 *   colptr      - Indexes in rows and avals of first element of
 *                 each column of the matrix.(output)
 *   rows        - Row of each element of the matrix.
 *   values      - Values of each element of the matrix.
 *   type        - type of the matrix.
 *   rhstype     - type of the right and side.
 *   driver_type - driver to use to read the matrix.
 *   pastix_comm - MPI communicator containing all processes wich
 *                 call read_matrix.
 */
int read_matrix(char            *filename,    pastix_int_t    *ncol,
                pastix_int_t   **colptr,      pastix_int_t   **rows,
                pastix_float_t **values,      pastix_float_t **rhs,
                char           **type,        char           **rhstype,
                driver_type_t    driver_type, MPI_Comm         pastix_comm)
{
  pastix_int_t  *loc2glob;
  int mpid;

  MPI_Comm_rank(pastix_comm,&mpid);
  {
    int ret;
    if (EXIT_SUCCESS != ( ret = read_matrix_common(filename,
                                                   ncol, colptr, rows,
                                                   &loc2glob, values, rhs,
                                                   type, rhstype, driver_type,
                                                   pastix_comm)))
      return ret;
  }

  if (driver_type == CSCD || driver_type == FDUP_DIST || driver_type == MMD)
    {
      pastix_int_t    gn;
      pastix_int_t   *gcolptr;
      pastix_int_t   *grow;
      pastix_float_t *gavals;
      pastix_float_t *grhs;

      PASTIX_EXTERN(cscd2csc)(*ncol,*colptr,  *rows, *values, *rhs,  NULL, NULL,
                              &gn,  &gcolptr, &grow, &gavals, &grhs, NULL, NULL,
                              loc2glob, pastix_comm, 1);

      free(*colptr);
      free(*rows);
      free(*values);
      free(*rhs);
      free(loc2glob);
      *ncol   = gn;
      *colptr = gcolptr;
      *rows   = grow;
      *values = gavals;
      *rhs    = grhs;

    }

  return EXIT_SUCCESS;
}


/*
 * Function: dread_matrix
 *
 * Reads a matrix from a file in the format given by driver_type
 * (see <driver_type_enum>) and distribute the matrix on the MPI processors.
 *
 * Parameters:
 *   filename    - Name of the file to read from.
 *   ncol        - Number of column in the matrix (output).
 *   colptr      - Indexes in rows and avals of first element of each column of
 *                 the matrix.(output)
 *   rows        - Row of each element of the matrix.
 *   loc2glob    - Local to global column number indirection.
 *   values      - Values of each element of the matrix.
 *   type        - type of the matrix.
 *   rhstype     - type of the right and side.
 *   driver_type - driver to use to read the matrix.
 *   pastix_comm - MPI communicator containing all processes wich call
 *                 <dread_matrix>.
 */
int dread_matrix(char            *filename,    pastix_int_t    *ncol,
                 pastix_int_t   **colptr,      pastix_int_t   **rows,
                 pastix_int_t   **loc2glob,
                 pastix_float_t **values,      pastix_float_t **rhs,
                 char           **type,        char           **rhstype,
                 driver_type_t    driver_type, MPI_Comm         pastix_comm)
{
  int mpid;

  MPI_Comm_rank(pastix_comm,&mpid);
  {
    int ret;
    if (EXIT_SUCCESS != ( ret = read_matrix_common(filename,
                                                   ncol, colptr, rows,
                                                   loc2glob, values, rhs,
                                                   type, rhstype, driver_type,
                                                   pastix_comm)))
        return ret;
  }

  if (driver_type != CSCD && driver_type != FDUP_DIST && driver_type != MMD)
    {

      pastix_int_t    lN;
      pastix_int_t   *lcolptr;
      pastix_int_t   *lrow;
      pastix_float_t *lavals;
      pastix_float_t *lrhs;

      PASTIX_EXTERN(csc_dispatch)(*ncol, *colptr,  *rows, *values, *rhs,  NULL, NULL,
                                  &lN,   &lcolptr, &lrow, &lavals, &lrhs, NULL,
                                  loc2glob, CSC_DISP_SIMPLE, MPI_COMM_WORLD);

      memFree_null(*colptr);
      memFree_null(*rows);
      memFree_null(*values);
      memFree_null(*rhs);

      *ncol   = lN;
      *colptr = lcolptr;
      *rows   = lrow;
      *values = lavals;
      *rhs    = lrhs;
    }

  return EXIT_SUCCESS;
}


/*
 * structure: couple_
 *
 * Structure used to sort couples.
 */
typedef struct couple_
{
  pastix_int_t i,j;
} couple;

/*
 * Function: comparcouple
 *
 * Compares 2 couples (i,j).
 *
 * Parameters:
 *   a - one couple
 *   b - one other couple
 *
 * Returns:
 *  -1 - if a_i < b_i
 *   1 - if a_i > b_i
 *  -1 - if a_i = b_i and a_j < b_j
 *   1 - else
 */
int comparcouple(const void *a, const void *b)
{
  couple const * const aa=(const couple*) a;
  couple const * const bb=(const couple*) b;

  if (aa->i < bb->i)
  {
    return -1;
  }
  else
  {
    if (aa->i > bb->i)
    {
      return 1;
    }
    else
    {
      if (aa->j < bb->j)
      {
        return -1;
      }
      else
      {
        return 1;
      }
    }
  }
}

/*
 * Function: checkStrucSym
*
 * Check if unsymmetric matrix has a symmetric pattern, if not correct it.
*
 * Parameters:
 *   n      - size of the matrix
 *   nz     - number of elements
 *   colptr - Index of first element of each column in *row* and *avals*
 *   row    - row of each element
 *   avals  - value of each element
 */
void checkStrucSym(pastix_int_t     n,
                   pastix_int_t    *nz,
                   pastix_int_t   **colptr,
                   pastix_int_t   **row,
                   pastix_float_t **avals)
{
  pastix_int_t    itercol;
  pastix_int_t    iterval;
  pastix_int_t    iterval2;
  pastix_int_t    nbrlostelt=0; /* Number of element to have a symmetric pattern */
  double          maxvalue=0;
  couple         *tablostelt=NULL;
  pastix_int_t    tempnz;
  pastix_int_t   *tempcol;
  pastix_int_t   *temprow;
  pastix_float_t *tempval;
  pastix_int_t    iterlost;

  for (itercol=0; itercol<n;itercol++)
  {
    for (iterval=(*colptr)[itercol]-1; iterval<(*colptr)[itercol+1]-1; iterval++)
    {
      if ((*row)[iterval]-1 != itercol)
      {
        const pastix_int_t rowidx=(*row)[iterval]-1;
        pastix_int_t find=0;

        for (iterval2=(*colptr)[rowidx]-1; iterval2<(*colptr)[rowidx+1]-1;iterval2++)
        {
          if ((*row)[iterval2]-1 == itercol)
          {
            find=1;
            break;
          }
        }
        if (find == 0)
        {
          if (ABS_FLOAT((*avals)[iterval]) > maxvalue)
            maxvalue = ABS_FLOAT((*avals)[iterval]);

          nbrlostelt++;
        }

      }
    }
  }

  fprintf(stderr, "nbrlostelt %ld maxvalue %e \n", (long) nbrlostelt, maxvalue);

  if (nbrlostelt!=0)
  {
    /* repair csc to have symmetric pattern */
    if (!(tablostelt = (couple *) malloc(nbrlostelt*sizeof(couple))))
      MALLOC_ERROR("tablostelt");

    nbrlostelt=0;

    for (itercol=0; itercol<n; itercol++)
    {
      for (iterval=(*colptr)[itercol]-1; iterval<(*colptr)[itercol+1]-1; iterval++)
      {
        if ((*row)[iterval]-1 != itercol)
        {
          const pastix_int_t rowidx = (*row)[iterval]-1;
          pastix_int_t find=0;

          for (iterval2=(*colptr)[rowidx]-1; iterval2<(*colptr)[rowidx+1]-1; iterval2++)
          {
            if ((*row)[iterval2]-1 == itercol)
            {
              find=1;
              break;
            }
          }
          if (find == 0)
          {
            tablostelt[nbrlostelt].i = (*row)[iterval]-1;
            tablostelt[nbrlostelt].j = itercol;
            nbrlostelt++;
          }
        }
      }
    }
    /* sort tablostelt */
    qsort(tablostelt,nbrlostelt,sizeof(couple),comparcouple);

    /* rebuild good format */
    tempnz = (*nz)+nbrlostelt;
    if (!(tempcol = (pastix_int_t *) malloc(n*sizeof(pastix_int_t))))
      MALLOC_ERROR("tempcol");
    if (!(temprow = (pastix_int_t *) malloc(tempnz*sizeof(pastix_int_t))))
      MALLOC_ERROR("temprow");
    if (!(tempval = (pastix_float_t *) malloc(tempnz*sizeof(pastix_float_t))))
      MALLOC_ERROR("tempval");

    iterlost = 0;
    tempcol[0] = (*colptr)[0];
    for (itercol=0; itercol<n; itercol++)
    {
      for (iterval=(*colptr)[itercol]-1; iterval<(*colptr)[itercol+1]-1; iterval++)
      {
        while ((tablostelt[iterlost].i == itercol) && (tablostelt[iterlost].j < (*row)[iterval]-1))
        {
          /* put elt */
          temprow[iterval+iterlost] = tablostelt[iterlost].j+1;
          tempval[iterval+iterlost] = 0;
          iterlost++;
        }
        temprow[iterval+iterlost] = (*row)[iterval];
        tempval[iterval+iterlost] = (*avals)[iterval];
      }

      while (tablostelt[iterlost].i == itercol)
      {
        /* put elt */
        temprow[iterval+iterlost] = tablostelt[iterlost].j+1;
        tempval[iterval+iterlost] = 0;
        iterlost++;
      }
      tempcol[itercol+1] = (*colptr)[itercol+1]+iterlost;
    }
    *nz=tempnz;
    memFree_null((*colptr));
    memFree_null((*row));
    memFree_null((*avals));
    *colptr=tempcol;
    *row=temprow;
    *avals=tempval;
  }
}
