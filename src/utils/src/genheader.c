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
/****************************************
 * genheader.c
 *
 * generate an executable to generate
 * headers for pastix users.
 * Used by makefile
 ****************************************/

#include "common_pastix.h"
#include <errno.h>

#include <string.h>

/*
  Program: genheader

  Generate headers (for fortran and C ) setting PASTIX_INT and
  MPI_PASTIX_INT maccros to the correct value.

  Parameters:
  argc - must be 3
  argv - contains the names wanted for the C header and Fortran header

  Returns:
  EXIT_SUCCESS

*/
int main (int argc, char ** argv)
{
  char * header_c  = NULL;
  char * header_f  = NULL;
  char * header_m  = NULL;
  char * header_mp = NULL;
  char * murge_h   = NULL;
  char * genfort   = NULL;
  char * insertinc = NULL;
  char * cmd       = NULL;
  char * argument  = NULL;
  FILE * file      = NULL;
  int    cplx      = 0;
  int    realsize  = (int)(sizeof(PASTIX_FLOAT)/sizeof(unsigned char));
  int    intsize   = (int)(sizeof(PASTIX_INT)/sizeof(unsigned char));
#ifdef INTSSIZE64
  int    sintsize  = (int)(sizeof(int64_t)/sizeof(unsigned char));
#else
  int    sintsize  = (int)(sizeof(int32_t)/sizeof(unsigned char));
#endif
#ifdef TYPE_COMPLEX
  cplx     = 1;
  realsize = (int)(sizeof(PASTIX_FLOAT)/sizeof(unsigned char)/2);
#endif

  if (argc != 9)
    {
      fprintf(stderr, "usage : %s headerNameC.h headerNameFortran.h"
              " headerMurge.inc headerMurgePastix.inc murge.h"
              " genfort.pl insert-fortran-inc.sh [C|Fortran|Murge]\n",argv[0]);
      return EXIT_FAILURE;
    }

  header_c  = argv[1];
  header_f  = argv[2];
  header_m  = argv[3];
  header_mp = argv[4];
  murge_h   = argv[5];
  genfort   = argv[6];
  insertinc = argv[7];
  argument  = argv[8];

  if (!strcmp(argument, "C")) {
    file = fopen (header_c,"w");
    if (file == NULL) {
      fprintf(stdout, "Error opening %s: %s\n", header_c, strerror(errno));
      abort();
    }

    /* #ifdef TYPE_COMPLEX */
    /* #if (defined X_ARCHalpha_compaq_osf1) */
    /* #ifndef   _RWSTD_HEADER_REQUIRES_HPP */
    /*   fprintf(file, "#include <complex>\n\n"); */
    /* #else  /\* _RWSTD_HEADER_REQUIRES_HPP *\/ */
    /*   fprintf(file, "#include <complex.hpp>\n\n"); */
    /* #endif /\* _RWSTD_HEADER_REQUIRES_HPP *\/ */
    /* #else /\*X_ARCHalpha_compaq_osf1*\/ */
    /*   fprintf(file, "#include <complex.h>\n\n"); */
    /* #endif */
    /* #endif */

    fprintf(file,
            "#ifndef PASTIX_INT_T_AND_SO_ON\n"
            "#define PASTIX_INT_T_AND_SO_ON\n");
#ifdef FORCE_LONG
    fprintf(file, "typedef long          pastix_int_t;\n");
    fprintf(file, "typedef unsigned long pastix_uint_t;\n");
    fprintf(file, "#  define MPI_PASTIX_INT   MPI_LONG\n");
#else
#ifdef FORCE_INT32
    fprintf(file, "typedef int32_t       pastix_int_t;\n");
    fprintf(file, "typedef uint32_t      pastix_uint_t;\n");
    fprintf(file, "#  define MPI_PASTIX_INT   MPI_INTEGER4\n");
#else
#ifdef FORCE_INT64
    fprintf(file, "typedef int64_t       pastix_int_t;\n");
    fprintf(file, "typedef uint64_t      pastix_uint_t;\n");
    fprintf(file, "#  define MPI_PASTIX_INT   MPI_INTEGER8\n");
#else
    fprintf(file, "typedef int           pastix_int_t;\n");
    fprintf(file, "typedef unsigned int  pastix_uint_t;\n");
    fprintf(file, "#  define MPI_PASTIX_INT   MPI_INT\n");
#endif
#endif
#endif
#ifdef PREC_DOUBLE
#ifdef TYPE_COMPLEX
    fprintf(file, "#  if (defined _COMPLEX_H || defined _H_COMPLEX"
            " || defined __COMPLEX__ || defined __COMPLEX_H__ ||"
            " (defined _GLIBCXX_HAVE_COMPLEX_H && _GLIBCXX_HAVE_COMPLEX_H == 1)"
            " || defined __STD_COMPLEX || defined _STLP_template_complex || defined _LIBCPP_COMPLEX )\n");
    fprintf(file,
      "#    ifdef __cplusplus\n"
      "#      ifndef COMPLEXDOUBLE_\n"
      "#        define COMPLEXDOUBLE_\n"
      "#      endif\n"
      "       typedef std::complex<double>  pastix_float_t;\n"
      "#    else\n"
      "       typedef double complex pastix_float_t;\n"
      "#    endif\n");
    fprintf(file, "#    define MPI_PASTIX_FLOAT MPI_DOUBLE_COMPLEX\n");
    fprintf(file, "#    define PASTIX_FLOAT            pastix_float_t\n");
    fprintf(file, "#  endif\n");
#else
    fprintf(file, "typedef double pastix_float_t;\n");
    fprintf(file, "#  define MPI_PASTIX_FLOAT MPI_DOUBLE\n");
    fprintf(file, "#  define PASTIX_FLOAT            pastix_float_t\n");
#endif
#else
#ifdef TYPE_COMPLEX
    fprintf(file, "#  if (defined _COMPLEX_H || defined _H_COMPLEX "
            " || defined __COMPLEX__ || __COMPLEX_H__ ||"
            " ( defined _GLIBCXX_HAVE_COMPLEX_H && _GLIBCXX_HAVE_COMPLEX_H == 1)"
            " || defined __STD_COMPLEX || defined _STLP_template_complex || defined _LIBCPP_COMPLEX )\n");
    fprintf(file,
      "#    ifdef __cplusplus\n"
      "#      ifndef COMPLEXFLOAT_\n"
      "#        define COMPLEXFLOAT_\n"
      "#      endif\n"
      "      typedef std::complex<float> pastix_float_t;\n"
      "#    else\n"
      "      typedef float complex pastix_float_t;\n"
      "#    endif\n");
    fprintf(file, "#    define MPI_PASTIX_FLOAT MPI_COMPLEX\n");
    fprintf(file, "#    define PASTIX_FLOAT            pastix_float_t\n");
    fprintf(file, "#  endif\n");
#else
    fprintf(file, "typedef float pastix_float_t;\n");
    fprintf(file, "#  define MPI_PASTIX_FLOAT MPI_FLOAT\n");
    fprintf(file, "#  define PASTIX_FLOAT            pastix_float_t\n");
#endif
#endif
    fprintf(file, "#  define PASTIX_INT              pastix_int_t\n");
    fprintf(file, "#  define PASTIX_UINT             pastix_uint_t\n");
#ifdef FORCE_NOMPI
    fprintf(file, "#  ifdef  MPI_Comm\n");
    fprintf(file, "#    undef  MPI_Comm\n");
    fprintf(file, "#  endif\n");
    fprintf(file, "#  define MPI_Comm int\n");
#endif
    fprintf(file, "#endif /* PASTIX_INT_T_AND_SO_ON */\n");

    fclose(file);
  }
  else if (!strcmp(argument, "Fortran")) {

    /* Fortran header */
    file = fopen (header_f,"w");

    fprintf(file, "#define PASTIX_INT_KIND    %d\n",
            (int)(sizeof(PASTIX_INT)/sizeof(unsigned char)));
    fprintf(file, "#define pastix_int_t       INTEGER(kind=%d)\n",
            (int)(sizeof(PASTIX_INT)/sizeof(unsigned char)));
    fprintf(file, "#define pastix_uint_t      unsigned INTEGER(kind=%d)\n",
            (int)(sizeof(PASTIX_INT)/sizeof(unsigned char)));
    fprintf(file, "#define pastix_data_ptr_t  INTEGER(kind=%d)\n",
            (int)(sizeof(PASTIX_INT*)/sizeof(unsigned char)));
    fprintf(file, "#define MPI_PASTIX_INT     MPI_INTEGER%d\n",
            (int)(sizeof(PASTIX_INT)/sizeof(unsigned char)));

#ifdef TYPE_COMPLEX
    fprintf(file, "#define pastix_float_t     COMPLEX(kind=%d)\n",
            (int)(sizeof(PASTIX_FLOAT)/sizeof(unsigned char)/2));
    fprintf(file, "#define MPI_PASTIX_FLOAT   MPI_COMPLEX%d\n",
            (int)(sizeof(PASTIX_FLOAT)/sizeof(unsigned char)));
#else
    fprintf(file, "#define pastix_float_t     REAL(kind=%d)\n",
            (int)(sizeof(PASTIX_FLOAT)/sizeof(unsigned char)));
    fprintf(file, "#define MPI_PASTIX_FLOAT   MPI_REAL%d\n",
            (int)(sizeof(PASTIX_FLOAT)/sizeof(unsigned char)));
#endif

    fclose(file);
  }
  else if (! strcmp(argument, "Murge")) {

    /* Murge Fortran header */
    if (NULL == (cmd = (char*)malloc((strlen(genfort)+
                                      strlen(murge_h)+
                                      strlen(header_m)+128)*sizeof(char))))
      return EXIT_FAILURE;
    sprintf(cmd, "perl %s -f %s -c %d -r %d -s %d -l %d > %s;",
            genfort, murge_h, cplx,
            realsize, sintsize, intsize, header_m);
    if (-1 == system(cmd))
      return EXIT_FAILURE;
    free(cmd);
    if (NULL == (cmd = (char*)malloc((strlen(header_mp)*3+
                                      strlen(insertinc)+
                                      strlen(header_m)+512)*sizeof(char))))
      return EXIT_FAILURE;
    sprintf(cmd, "sed -e 's/INTS/INTEGER(KIND=%d)/g'", sintsize);
    sprintf(cmd, "%s  -e 's/INTL/INTEGER(KIND=%d)/g'", cmd, intsize);
#ifdef TYPE_COMPLEX
    sprintf(cmd, "%s  -e 's/COEF/COMPLEX(KIND=%d)/g'", cmd, realsize);
#else
    sprintf(cmd, "%s  -e 's/COEF/REAL(KIND=%d)/g'", cmd, realsize);
#endif
    sprintf(cmd, "%s %s > %s.tmp; %s %s %s.tmp;",cmd,
            header_mp, header_mp, insertinc, header_m, header_mp);
    //fprintf(stdout, "GENHEADER DEBUG: %s\n", cmd);
    if (0 != system(cmd))
      return 1;
    free(cmd);
  }
  return 0;
}
