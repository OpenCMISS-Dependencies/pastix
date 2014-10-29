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
 * File: read_matrix.h
 *
 * Definition of a global function to read all type of matrices.
 *
 */

#ifndef _GLIBCXX_HAVE_COMPLEX_H
#  define _GLIBCXX_HAVE_COMPLEX_H 0
#endif

#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || defined __COMPLEX_H__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX || defined _STLP_template_complex  || defined _LIBCPP_COMPLEX )
#  define PASTIX_HAS_COMPLEX
#endif

#ifdef PASTIX_HAS_COMPLEX
#  ifdef   __cplusplus
#    define  COMPLEX  std::complex<float>
#    define  DCOMPLEX std::complex<double>
#  else /* not __cplusplus */
#    define  COMPLEX float complex
#    define  DCOMPLEX double complex
#  endif /* not __cplusplus */
#endif

/*
 * typedef: driver_type_enum
 *
 * RSA        - rsa file format.
 * CHB        - chb file format.
 * CCC        - ccc file format.
 * RCC        - rcc file format.
 * OLAF       - olaf file format.
 * PEER       - peer file format.
 * HB         - Harwell-boeing file format.
 * THREEFILES - IJV file format with one file for each element of the triple.
 * MM         - Matrix Market file format.
 * PETSCS     - PETSc binary file format symmetric.
 * PETSCH     - PETSc binary file format hermitian.
 * PETSCU     - PETSc binary file format unsymmetric.
 * CSCD       - CSCd file format.
 * LAPLACIAN  - Generating laplacian.
 * FDUP       - Binary matrix driver from Fabrice dupros.
 */
#ifndef READ_MATRIX_H
#define READ_MATRIX_H
enum driver_type_enum {
  RSA,
  CHB,
  CCC,
  RCC,
  OLAF,
  PEER,
  HB,
  THREEFILES,
  MM,
  MMD,
  PETSCS,
  PETSCU,
  PETSCH,
  CSCD,
  LAPLACIAN,
  FDUP,
  FDUP_DIST
};
/*
  typedef: driver_type_t
 */
typedef enum driver_type_enum driver_type_t;

#ifndef MTX_ISSYM
#  define MTX_ISSYM(a) ((a)[1]=='S')
#  define MTX_ISHER(a) ((a)[1]=='H')
#  define MTX_ISCOM(a) ((a)[0]=='C')
#  define MTX_ISRHX(a) ((a)[2]=='X')
#  define MTX_ISRHS(a) ((a)[0]!='\0')
#endif

#ifdef PASTIX_INT_T_AND_SO_ON
/*
 * Function: read_matrix
 *
 * Reads a matrix from a file in the format given by driver_type
 *  (see <driver_type_enum>).
 *
 * Parameters:
 *   filename    - Name of the file to read from.
 *   ncol        - Number of column in the matrix (output).
 *   colptr      - Indexes in row and avals of first element of each column
 *                 of the matrix.(output)
 *   row         - Row of each element of the matrix.
 *   values      - Values of each element of the matrix.
 *   type        - type of the matrix.
 *   rhstype     - type of the right and side.
 *   driver_type - driver to use to read the matrix.
 *   pastix_comm - MPI communicator containing all processes wich call
 *                 <read_matrix>.
 */
int read_matrix(char            *filename,
                pastix_int_t    *ncol,
                pastix_int_t   **colptr,
                pastix_int_t   **row,
                pastix_float_t **values,
                pastix_float_t **rhs,
                char           **type,
                char           **rhstype,
                driver_type_t    driver_type,
                MPI_Comm         pastix_comm);
#endif
int s_read_matrix(char            *filename,
                  pastix_int_t    *ncol,
                  pastix_int_t   **colptr,
                  pastix_int_t   **row,
                  float **values,
                  float **rhs,
                  char           **type,
                  char           **rhstype,
                  driver_type_t    driver_type,
                  MPI_Comm         pastix_comm);
int d_read_matrix(char            *filename,
                  pastix_int_t    *ncol,
                  pastix_int_t   **colptr,
                  pastix_int_t   **row,
                  double **values,
                  double **rhs,
                  char           **type,
                  char           **rhstype,
                  driver_type_t    driver_type,
                  MPI_Comm         pastix_comm);
#ifdef PASTIX_HAS_COMPLEX
int c_read_matrix(char            *filename,
                  pastix_int_t    *ncol,
                  pastix_int_t   **colptr,
                  pastix_int_t   **row,
                  COMPLEX **values,
                  COMPLEX **rhs,
                  char           **type,
                  char           **rhstype,
                  driver_type_t    driver_type,
                  MPI_Comm         pastix_comm);
int z_read_matrix(char            *filename,
                  pastix_int_t    *ncol,
                  pastix_int_t   **colptr,
                  pastix_int_t   **row,
                  DCOMPLEX **values,
                  DCOMPLEX **rhs,
                  char           **type,
                  char           **rhstype,
                  driver_type_t    driver_type,
                  MPI_Comm         pastix_comm);
#endif

#ifdef PASTIX_INT_T_AND_SO_ON
/*
 * Function: dread_matrix
 *
 * Reads a matrix from a file in the format given by driver_type
 * (see <driver_type_enum>) and distribute the matrix on the MPI processors.
 *
 * Parameters:
 *   filename    - Name of the file to read from.
 *   ncol        - Number of column in the matrix (output).
 *   colptr      - Indexes in rows and avals of first element of each column
 *                 of the matrix.(output)
 *   rows        - Row of each element of the matrix.
 *   values      - Values of each element of the matrix.
 *   type        - type of the matrix.
 *   rhstype     - type of the right and side.
 *   driver_type - driver to use to read the matrix.
 *   pastix_comm - MPI communicator containing all processes wich call
 *                 <read_matrix>.
 */
int dread_matrix(char            *filename,
                 pastix_int_t    *ncol,
                 pastix_int_t   **colptr,
                 pastix_int_t   **rows,
                 pastix_int_t   **loc2glob,
                 pastix_float_t **values,
                 pastix_float_t **rhs,
                 char           **type,
                 char           **rhstype,
                 driver_type_t    driver_type,
                 MPI_Comm         pastix_comm);
#endif
int s_dread_matrix(char            *filename,
                   pastix_int_t    *ncol,
                   pastix_int_t   **colptr,
                   pastix_int_t   **rows,
                   pastix_int_t   **loc2glob,
                   float **values,
                   float **rhs,
                   char           **type,
                   char           **rhstype,
                   driver_type_t    driver_type,
                   MPI_Comm         pastix_comm);
int d_dread_matrix(char            *filename,
                   pastix_int_t    *ncol,
                   pastix_int_t   **colptr,
                   pastix_int_t   **rows,
                   pastix_int_t   **loc2glob,
                   double **values,
                   double **rhs,
                   char           **type,
                   char           **rhstype,
                   driver_type_t    driver_type,
                   MPI_Comm         pastix_comm);
#ifdef PASTIX_HAS_COMPLEX
int c_dread_matrix(char            *filename,
                   pastix_int_t    *ncol,
                   pastix_int_t   **colptr,
                   pastix_int_t   **rows,
                   pastix_int_t   **loc2glob,
                   COMPLEX **values,
                   COMPLEX **rhs,
                   char           **type,
                   char           **rhstype,
                   driver_type_t    driver_type,
                   MPI_Comm         pastix_comm);
int z_dread_matrix(char            *filename,
                   pastix_int_t    *ncol,
                   pastix_int_t   **colptr,
                   pastix_int_t   **rows,
                   pastix_int_t   **loc2glob,
                   DCOMPLEX **values,
                   DCOMPLEX **rhs,
                   char           **type,
                   char           **rhstype,
                   driver_type_t    driver_type,
                   MPI_Comm         pastix_comm);

#endif
/*
 *  Function: comparcouple
 *
 *  Compares 2 couples (i,j).
 *
 *  Parameters:
 *    a - one couple
 *    b - one other couple
 *
 *  Returns:
 *    -1 - if a_i < b_i
 *    1  - if a_i > b_i
 *    -1 - if a_i = b_i and a_j < b_j
 *    1  - else
 */
int comparcouple(const void *a, const void *b);
int s_comparcouple(const void *a, const void *b);
int d_comparcouple(const void *a, const void *b);
#ifdef PASTIX_HAS_COMPLEX
int c_comparcouple(const void *a, const void *b);
int z_comparcouple(const void *a, const void *b);
#endif

#ifdef PASTIX_INT_T_AND_SO_ON
/*
 *  Function: checkStrucSym
 *
 *  Check if unsymmetric matrix has a symmetric pattern, if not correct it.
 *
 *  Parameters:
 *    n      - size of the matrix
 *    nz     - number of elements
 *    colptr - Index of first element of each column in *row* and *avals*
 *    row    - row of each element
 *    avals  - value of each element
 */
void checkStrucSym(pastix_int_t     n,
                   pastix_int_t    *nz,
                   pastix_int_t   **colptr,
                   pastix_int_t   **row,
                   pastix_float_t **avals);
#endif
void s_checkStrucSym(pastix_int_t     n,
                     pastix_int_t    *nz,
                     pastix_int_t   **colptr,
                     pastix_int_t   **row,
                     float **avals);
void d_checkStrucSym(pastix_int_t     n,
                     pastix_int_t    *nz,
                     pastix_int_t   **colptr,
                     pastix_int_t   **row,
                     double **avals);
#ifdef PASTIX_HAS_COMPLEX
void c_checkStrucSym(pastix_int_t     n,
                     pastix_int_t    *nz,
                     pastix_int_t   **colptr,
                     pastix_int_t   **row,
                     COMPLEX **avals);
void z_checkStrucSym(pastix_int_t     n,
                     pastix_int_t    *nz,
                     pastix_int_t   **colptr,
                     pastix_int_t   **row,
                     DCOMPLEX **avals);
#endif
#undef PASTIX_HAS_COMPLEX
#undef COMPLEX
#undef DCOMPLEX
#endif /* not READ_MATRIX_H */
