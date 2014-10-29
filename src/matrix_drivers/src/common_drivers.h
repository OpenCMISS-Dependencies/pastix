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
 *  File: common_drivers.h
 *
 *  Definition common to all drivers.
 */
#ifndef COMMON_DRIVER_H
#define COMMON_DRIVER_H
#define STR_SIZE 256

#define MTX_ISSYM(a) ((a)[1]=='S')
#define MTX_ISHER(a) ((a)[1]=='H')
#define MTX_ISCOM(a) ((a)[0]=='C')
#define MTX_ISRHX(a) ((a)[2]=='X')
#define MTX_ISRHS(a) ((a)[0]!='\0')

#define memFree_null(x) {if (x ==NULL) {fprintf(stdout,"%s:%d freeing NULL\n",__FILE__,__LINE__);} free(x); x=NULL;}

#ifdef TYPE_COMPLEX
#if (defined X_ARCHalpha_compaq_osf1)
#ifndef USE_CXX
#define ABS_FLOAT(x) abs(x)
#endif /* USE_CXX */
#else /*X_ARCHalpha_compaq_osf1*/
#ifdef    PREC_DOUBLE
#define ABS_FLOAT(x) cabs(x)
#else  /* PREC_DOUBLE */
#define ABS_FLOAT(x) cabsf(x)
#endif /* PREC_DOUBLE */
#endif /* X_ARCHalpha_compaq_osf1 */
#else /* TYPE_COMPLEX */
#ifdef PREC_DOUBLE
#define ABS_FLOAT(x) fabs(x)
#else /* PREC_DOUBLE */
#define ABS_FLOAT(x) fabsf(x)
#endif /* FORCEDOUBLE */

#endif /* TYPE_COMPLEX */

#define FGETS(line, BUFSIZ, infile) {					\
  if (NULL == fgets(line, BUFSIZ, infile))				\
    {									\
      fprintf(stderr, "ERROR: %s:%d fgets\n", __FILE__, __LINE__);	\
      exit(1);								\
    }									\
  }
#define ASSERT(expr,module){						\
    if((expr) == 0){							\
      fprintf(stderr,"error in assert (line=%d,file=%s)\n",		\
        __LINE__,__FILE__);					\
      exit(EXIT_FAILURE);}}

#define EXIT(module,error) {			\
    exit(EXIT_FAILURE);}

#define MALLOC_ERROR(x) {						\
    fprintf(stderr,"\nERROR %s allocation (line=%d,file=%s)\n\n",	\
      (x),__LINE__,__FILE__);					\
    EXIT(MOD_UNKNOWN,ALLOC_ERR);}

#ifdef PASTIX_RENAME

#ifdef FORGET_TYPE
#ifdef PREC_DOUBLE
#ifdef TYPE_COMPLEX
#define PASTIX_PREFIX(x) _PASTIX_Z_ ## x
#define PASTIX_EXTERN(x) z_ ## x
#else
#define PASTIX_PREFIX(x) _PASTIX_D_ ## x
#define PASTIX_EXTERN(x) d_ ## x
#endif
#else
#ifdef TYPE_COMPLEX
#define PASTIX_PREFIX(x) _PASTIX_C_ ## x
#define PASTIX_EXTERN(x) c_ ## x
#else
#define PASTIX_PREFIX(x) _PASTIX_S_ ## x
#define PASTIX_EXTERN(x) s_ ## x
#endif
#endif

#else  /* FORGET_TYPE*/

#define PASTIX_PREFIX(x) _PASTIX_ ## x
#define PASTIX_EXTERN(x) x

#endif

#else /* PASTIX_RENAME */

#ifdef FORGET_TYPE

#ifdef PREC_DOUBLE
#ifdef TYPE_COMPLEX
#define PASTIX_PREFIX(x) Z_ ## x
#define PASTIX_EXTERN(x) z_ ## x
#else
#define PASTIX_PREFIX(x) D_ ## x
#define PASTIX_EXTERN(x) d_ ## x
#endif
#else
#ifdef TYPE_COMPLEX
#define PASTIX_PREFIX(x) C_ ## x
#define PASTIX_EXTERN(x) c_ ## x
#else
#define PASTIX_PREFIX(x) S_ ## x
#define PASTIX_EXTERN(x) s_ ## x
#endif
#endif

#else  /* FORGET_TYPE*/

#define PASTIX_PREFIX(x) x
#define PASTIX_EXTERN(x) x

#endif
#endif /* PASTIX_RENAME */

#ifdef FORGET_TYPE
#undef pastix_float_t
#undef MPI_PASTIX_FLOAT

#ifdef PREC_DOUBLE
#ifdef TYPE_COMPLEX
#define pastix_float_t complex double
#define MPI_PASTIX_FLOAT MPI_DOUBLE_COMPLEX
#else
#define pastix_float_t double
#define MPI_PASTIX_FLOAT MPI_DOUBLE
#endif
#else /* PREC_DOUBLE */
#ifdef TYPE_COMPLEX
#define pastix_float_t complex float
#define MPI_PASTIX_FLOAT MPI_COMPLEX
#else
#define pastix_float_t float
#define MPI_PASTIX_FLOAT MPI_FLOAT
#endif
#endif

#endif

#define IOHBTerminate        PASTIX_PREFIX(IOHBTerminate)
#define ParseIfmt            PASTIX_PREFIX(ParseIfmt)
#define ParseRfmt            PASTIX_PREFIX(ParseRfmt)
#define substr               PASTIX_PREFIX(substr)
#define upcase               PASTIX_PREFIX(upcase)
#define writeHB_mat_char     PASTIX_PREFIX(writeHB_mat_char)
#define writeHB_mat_double   PASTIX_PREFIX(writeHB_mat_double)
#define checkStrucSym        PASTIX_PREFIX(checkStrucSym)
#define comparcouple         PASTIX_PREFIX(comparcouple)
#define dread_matrix         PASTIX_EXTERN(dread_matrix)
#define read_matrix          PASTIX_EXTERN(read_matrix)
#define read_matrix_common   PASTIX_EXTERN(read_matrix_common)
#define rsaRead              PASTIX_PREFIX(rsaRead)
#define rsaReadHeader        PASTIX_PREFIX(rsaReadHeader)
#define HBRead               PASTIX_PREFIX(HBRead)
#define MatrixMarketRead     PASTIX_PREFIX(MatrixMarketRead)
#define cccRead              PASTIX_PREFIX(cccRead)
#define cccReadHeader        PASTIX_PREFIX(cccReadHeader)
#define olafRead             PASTIX_PREFIX(olafRead)
#define olafReadHeader       PASTIX_PREFIX(olafReadHeader)
#define chbParseIfmt         PASTIX_PREFIX(chbParseIfmt)
#define chbParseRfmt         PASTIX_PREFIX(chbParseRfmt)
#define chbRead              PASTIX_PREFIX(chbRead)
#define chbReadHeader        PASTIX_PREFIX(chbReadHeader)
#define cscdRead             PASTIX_PREFIX(cscdRead)
#define peerRead             PASTIX_PREFIX(peerRead)
#define peerRead2            PASTIX_PREFIX(peerRead2)
#define threeFilesRead       PASTIX_PREFIX(threeFilesRead)
#define threeFilesReadHeader PASTIX_PREFIX(threeFilesReadHeader)
#define genlaplacian         PASTIX_PREFIX(genlaplacian)
#define getfilename          PASTIX_PREFIX(getfilename)
#define getordering          PASTIX_PREFIX(getordering)
#define global_usage         PASTIX_PREFIX(global_usage)
#define str_tolower          PASTIX_PREFIX(str_tolower)
#define driverFdupros        PASTIX_PREFIX(driverFdupros)
#define driverFdupros_dist   PASTIX_PREFIX(driverFdupros_dist)
#define PETScRead            PASTIX_PREFIX(PETScRead)
/*
  Function: myupcase

  Rewrites *s* to upper case.

  Parameters:
    s - string to rexwrite in upcase.
*/
void myupcase(char *S);

/*
  Function:  mysubstr

  Copy len element, from *S[pos]* in *s*.

  Parameters:
    s   - destination
    S   - Source
    pos - sarting position
    len - Size of the string to copy.
*/
void mysubstr(char *s, const char *S, const pastix_int_t pos, const pastix_int_t len);

/*
  Function:  mysubstr2

  Copy the number placed between a and b in fmt.

  Parameters:
    fmt - String in which there is a and b
    a   - first element
    b   - last element
    val - the integer between a and b

*/
void mysubstr2(const char *fmt, const char a, const char b, pastix_int_t *val);
#endif /* not COMMON_DRIVER_H */
