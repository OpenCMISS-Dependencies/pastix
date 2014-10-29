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
  File: get_options.c

  Definition of a global function to get exemple parameters.

*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/types.h>
#include <stdint.h>

#ifdef FORCE_NOMPI
#include "pastix_nompi.h"
#else
#include <mpi.h>
#endif


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
#include "common_drivers.h"
#include "read_matrix.h"
#include "get_options.h"
#include "api_str_to_int.h"
#include <string.h>

int api_iparmreader(char * filename, pastix_int_t *iparmtab);
int api_dparmreader(char * filename, double *dparmtab);

/*
  Function: str_tolower

  Rewrites *string* in lower case.

  Parameters:
    string - string to rewrite in lower case.
*/
int str_tolower(char * string)
{
  int j = 0;
  while (string[j] != '\0')
    {
      string[j] = (char)tolower(string[j]);
      j++;
    }
  return EXIT_SUCCESS;
}

/*
  Function: getfilename

  Sets filename to source if source doesn't starts with '-'.
  Otherwise, filename is set to defaultname.

  Parameters:
    filename    - string to set to correct filename.
    source      - possible source for filename.
    defaultname - default filename.

  Returns:
    0 if set to default.
    1 if set to source.
*/
int getfilename(char ** filename, char * source, char * defaultname)
{
  if (source == NULL || source[0] == '-')
    {
      *filename = (char *) malloc((strlen(defaultname)+1)*sizeof(char));
      strcpy(*filename,defaultname);
      return 0;
    }
  *filename = (char *) malloc((strlen(source)+1)*sizeof(char));
  strcpy(*filename,source);
  return 1;
}

/*
  Function: getordering

  Sets *ordering* from source.

  Parameters:
    ordering    - integer to set to correct ordering.
    source      - source for ordering name.

  Returns:
    EXIT_SUCCESS if ordering exists.
    EXIT_FAILURE if ordering doesn't exists.
*/
int getordering(int  * ordering,
                char * source)
{
  if (strcmp(source, "scotch") == 0)
    {
      *ordering = API_ORDER_SCOTCH;
      return EXIT_SUCCESS;
    }
  if (strcmp(source, "metis") == 0)
    {
      *ordering = API_ORDER_METIS;
      return EXIT_SUCCESS;
    }
  return EXIT_FAILURE;
}

/*
  Function: global_usage

  Print usage corresponding to all pastix exemples.

  Parameters:
    mpi_comm - MPI communicator.
    argv     - program argument

*/
void global_usage(MPI_Comm mpi_comm, char ** argv)
{
  int rank;

  MPI_Comm_rank(mpi_comm, &rank);
  if (rank == 0)
    {
      fprintf(stdout, "Usage : %s [option] \n",argv[0]);
      fprintf(stdout, "\toptions : \n");
      fprintf(stdout, "\t\t -rsa     [filename]          driver RSA (use Fortran) \n");
      fprintf(stdout, "\t\t -chb     [filename]          driver CHB\n");
      fprintf(stdout, "\t\t -ccc     [filename]          driver CCC\n");
      fprintf(stdout, "\t\t -rcc     [filename]          driver RCC\n");
      fprintf(stdout, "\t\t -olaf    [filename]          driver OLAF\n");
      fprintf(stdout, "\t\t -peer    [filename]          driver PEER\n");
      fprintf(stdout, "\t\t -petsc_s [filename]          driver PETSc symmetric\n");
      fprintf(stdout, "\t\t -petsc_h [filename]          driver PETSc hermitian\n");
      fprintf(stdout, "\t\t -petsc_u [filename]          driver PETSc unsymmetric\n");
      fprintf(stdout, "\t\t -hb      [filename]          driver HB (double)\n");
      fprintf(stdout, "\t\t -3files  [filename]          driver IJV 3files \n");
      fprintf(stdout, "\t\t -mm      [filename]          driver Matrix Market\n");
      fprintf(stdout, "\t\t -dmm     [filename]          driver Matrix Market (distributed)\n");
#ifdef FDUPROS
      fprintf(stdout, "\t\t -fdup    [filename]          driver from Fabrice Dupros\n");
      fprintf(stdout, "\t\t -fdupd   [filename]          driver from Fabrice Dupros, distributed\n");
#endif
      fprintf(stdout, "\t\t -ord     <scotch|metis>      select ordering library\n");
      fprintf(stdout, "\t\t -lap     <integer>           generate a laplacian of size <integer>\n");
      fprintf(stdout, "\t\t -incomp  <integer> <integer> incomplete factorization, with the given level of fill [1-5],\n");
      fprintf(stdout, "\t\t                              and amalgamation [10-70]\n");
      fprintf(stdout, "\t\t -ooc     <integer>           Memory limit in Mo/percent depending on compilation options\n");
      fprintf(stdout, "\t\t -kass    <integer>           kass, with the given amalgamation\n");
      fprintf(stdout, "\t\t -t       <integer>           define thread number\n");
      fprintf(stdout, "\t\t -v       <integer>           define verbose level (1,2 or 3)\n");
      fprintf(stdout, "\t\t -iparm   <IPARM_ID> <value>  set an integer parameter\n");
      fprintf(stdout, "\t\t -dparm   <DPARM_ID> <value>  set a floating parameter\n");

      /*       fprintf(stdout, "\t\t b         driver \"Fabrice Dupros\"\n"); */
      fprintf(stdout, "\t\t -h                          print this help\n");
    }
}
/*
  Function: get_options

  Get options from argv.

  Parameters:
  argc          - number of arguments.
  argv          - argument tabular.
  driver_type   - type of driver (output, -1 if not set).
  filename      - Matrix filename (output).
  nbmatrices    - number of matrices in arguments.
  nbthread      - number of thread (output, 1 if not set).
  verbose       - verbose level 1,2 or 3
  ordering      - ordering to choose (see <API_ORDER>).
  incomplete    - indicate if -incomp is present
  level_of_fill - Level of fill for incomplete factorization.
  amalgamation  - Amalgamation for kass.
  ooc           - Out-of-core limite (Mo or percent depending on compilation option)
  size          - Size of the matrix (generated matrix only)
*/
int get_options(int              argc,
                char           **argv,
                driver_type_t  **driver_type,
                char          ***filename,
                int             *nbmatrices,
                int             *nbthread,
                int             *verbose,
                int             *ordering,
                int             *incomplete,
                int             *level_of_fill,
                int             *amalgamation,
                int             *ooc,
                pastix_int_t    *size)
{

  int i = 1;
  int maxmatrices = 10;

  (*driver_type) = (driver_type_t*)malloc(maxmatrices*sizeof(driver_type_t));
  (*filename)    = (char **       )malloc(maxmatrices*sizeof(char*));
  *nbmatrices    = 0;
  *nbthread      = 1;
  *verbose       = 1;
  *size          = 0;
  *ordering      = API_ORDER_SCOTCH;
  *incomplete    = API_NO;
  *level_of_fill = 0;
  *amalgamation  = 5;
  *ooc           = 2000;

  if (argc == 1)
    goto usage;
  while(i < argc)
    {
      if (argv[i][0] == '-')
        {

          switch (argv[i][1]) {

          case 'c':
          case 'C':
            str_tolower(argv[i]);
            if (strcmp(argv[i], "-chb") == 0)
              {
                (*driver_type)[(*nbmatrices)] = CHB;
                i+=getfilename(&(*filename)[(*nbmatrices)],
                               (i+1<argc)?argv[i+1]:NULL, "rsaname");
                (*nbmatrices)++;
              }
            else
              {
                if (strcmp(argv[i], "-ccc") == 0)
                  {
                    (*driver_type)[(*nbmatrices)] = CCC;
                    i+= getfilename(&(*filename)[(*nbmatrices)],
                                    (i+1<argc)?argv[i+1]:NULL, "dirname");
                    (*nbmatrices)++;
                  }
                else
                  {
                    if (strcmp(argv[i], "-cscd") == 0)
                      {
                        (*driver_type)[(*nbmatrices)] = CSCD;
                        i+= getfilename(&(*filename)[(*nbmatrices)],
                                        (i+1<argc)?argv[i+1]:NULL, "dirname");
                        (*nbmatrices)++;
                      }
                    else
                      goto unknown_option;
                  }
              }
            break;
          case 'd':
          case 'D':
            str_tolower(argv[i]);
            if (strcmp(argv[i], "-dparmfile") == 0)
              {
                i++;
              }
            else if (strcmp(argv[i], "-dparm") == 0)
              {
                i+=2;
              }
            else if (strcmp(argv[i],"-dmm") == 0 ||
                     strcmp(argv[i],"distributedmatrixmarket") == 0)
              {
                (*driver_type)[(*nbmatrices)] = MMD;
                i += getfilename(&(*filename)[(*nbmatrices)],
                                 (i+1<argc)?argv[i+1]:NULL, "mmname");
                (*nbmatrices)++;
              }
            else
              goto unknown_option;
            break;
#ifdef FDUPROS
          case 'f':
          case 'F':
            {
              str_tolower(argv[i]);
              if (strcmp(argv[i], "-fdup") == 0)
                {
                  (*driver_type)[(*nbmatrices)] = FDUP;
                  i+=getfilename(&(*filename)[(*nbmatrices)],
                                 (i+1<argc)?argv[i+1]:NULL, "dirname");
                  (*nbmatrices)++;
                }
              else
                {
                  if (strcmp(argv[i], "-fdupd") == 0)
                    {
                      (*driver_type)[(*nbmatrices)] = FDUP_DIST;
                      i+=getfilename(&(*filename)[(*nbmatrices)],
                                     (i+1<argc)?argv[i+1]:NULL, "dirname");
                      (*nbmatrices)++;
                    }
                  else
                    {
                      goto unknown_option;
                    }
                }
            }
            break;
#endif

          case 'h':
          case 'H':
            str_tolower(argv[i]);
            if (strcmp(argv[i],"-h") ==0 || strcmp(argv[i],"-help") ==0)
              goto usage;
            else
              if (strcmp(argv[i],"-hb") == 0 ||
                  strcmp(argv[i],"-harwell-boeing") == 0 ||
                  strcmp(argv[i],"-harwellboeing") == 0)
                {
                  (*driver_type)[(*nbmatrices)] = HB;
                  i+=getfilename(&(*filename)[(*nbmatrices)],
                                 (i+1<argc)?argv[i+1]:NULL, "rsaname");
                  (*nbmatrices)++;
                }
            break;
          case 'i':
          case 'I':
            str_tolower(argv[i]);
            if (strcmp(argv[i], "-incomp") == 0)
              {
                *incomplete    = API_YES;
                *level_of_fill = atoi(argv[i+1]);
                i++;
                *amalgamation  = atoi(argv[i+1]);
                i++;
              }
            else if (strcmp(argv[i], "-iparmfile") == 0)
              {
                i++;
              }
            else if (strcmp(argv[i], "-iparm") == 0)
              {
                i+=2;
              }
            else
              goto unknown_option;
            break;

          case 'k':
          case 'K':
            str_tolower(argv[i]);
            if (strcmp(argv[i], "-kass") == 0)
              {
                *level_of_fill = -1;
                *amalgamation  = atoi(argv[i+1]);
                i++;
              }
            else
              goto unknown_option;
            break;
          case 'l':
          case 'L':
            str_tolower(argv[i]);
            if (strcmp(argv[i], "-lap") == 0)
              {
                (*driver_type)[(*nbmatrices)] = LAPLACIAN;
                *size = atoi(argv[i+1]);
                (*filename)[(*nbmatrices)] = NULL;
                if (0 == *size)
                  goto unknown_option;
                i++;
                (*nbmatrices)++;
              }
            else
              goto unknown_option;
            break;

          case 'm':
          case 'M':
            str_tolower(argv[i]);
            if (strcmp(argv[i],"-mm") == 0 ||
                strcmp(argv[i],"matrixmarket") == 0)
              {
                (*driver_type)[(*nbmatrices)] = MM;
                i += getfilename(&(*filename)[(*nbmatrices)],
                                 (i+1<argc)?argv[i+1]:NULL, "mmname");
                (*nbmatrices)++;
              }
            else
              goto unknown_option;
            break;

          case 'o':
          case 'O':
            str_tolower(argv[i]);
            if (strcmp(argv[i],"-olaf") == 0)
              {
                (*driver_type)[(*nbmatrices)] = OLAF;
                i+= getfilename(&(*filename)[(*nbmatrices)],
                                (i+1<argc)?argv[i+1]:NULL, "olafcsr");
                (*nbmatrices)++;
              }
            else
              {
                if (strcmp(argv[i],"-ord") == 0)
                  {
                    if (EXIT_FAILURE ==
                        getordering(ordering,(i+1<argc)?argv[i+1]:NULL))
                      goto usage;
                    else
                      i++;
                  }
                else
                  {
                    if (strcmp(argv[i], "-ooc") == 0)
                      {
                        *ooc = atoi(argv[i+1]);
                        if (0 == *ooc)
                          goto unknown_option;
                        i++;
                      }
                    else
                      goto unknown_option;
                  }
              }
            break;

          case 'p':
          case 'P':
            str_tolower(argv[i]);
            if (strcmp(argv[i],"-peer") == 0)
              {
                (*driver_type)[(*nbmatrices)] = PEER;
                i+= getfilename(&(*filename)[(*nbmatrices)],
                                (i+1<argc)?argv[i+1]:NULL, "rsaname");
                (*nbmatrices)++;
              }
            else
              {
                if ( strcmp(argv[i], "-petsc_s") == 0 )
                  {
                    (*driver_type)[(*nbmatrices)] = PETSCS;
                    i+= getfilename(&(*filename)[(*nbmatrices)],
                                    (i+1<argc)?argv[i+1]:NULL, "PETSCFILE");
                    (*nbmatrices)++;
                  }
                else
                  {
                    if ( strcmp(argv[i], "-petsc_u") == 0 )
                      {
                        (*driver_type)[(*nbmatrices)] = PETSCU;
                        i+= getfilename(&(*filename)[(*nbmatrices)],
                                        (i+1<argc)?argv[i+1]:NULL, "PETSCFILE");
                        (*nbmatrices)++;
                      }
                  else
                    {
                      if ( strcmp(argv[i], "-petsc_h") == 0 )
                        {
                          (*driver_type)[(*nbmatrices)] = PETSCH;
                          i+= getfilename(&(*filename)[(*nbmatrices)],
                                          (i+1<argc)?argv[i+1]:NULL, "PETSCFILE");
                          (*nbmatrices)++;
                        }
                      else {
                        goto unknown_option;
                      }
                    }
                  }
              }
            break;

          case 'r':
          case 'R':
            str_tolower(argv[i]);
            if (strcmp(argv[i],"-rsa") == 0)
              {
                (*driver_type)[(*nbmatrices)] = RSA;
                i+= getfilename(&(*filename)[(*nbmatrices)],
                                (i+1<argc)?argv[i+1]:NULL, "rsaname");
                (*nbmatrices)++;
              }
            else
              {
                if (strcmp(argv[i],"-rcc") == 0)
                  {
                    (*driver_type)[(*nbmatrices)] = RCC;
                    i+= getfilename(&(*filename)[(*nbmatrices)],
                                    (i+1<argc)?argv[i+1]:NULL, "dirname");
                    (*nbmatrices)++;
                  }
                else
                  goto unknown_option;
              }
            break;

          case 't':
          case 'T':
            str_tolower(argv[i]);
            if (strcmp(argv[i], "-t") == 0)
              {
                *nbthread = atoi(argv[i+1]);
                if (0 == *nbthread)
                  goto unknown_option;
                i++;
              }
            else
              goto unknown_option;
            break;

          case 'v':
            str_tolower(argv[i]);
            if (strcmp(argv[i], "-v") == 0)
              {
                *verbose = atoi(argv[i+1]);
                if (0 == *verbose)
                  goto unknown_option;
                i++;
                (*verbose)--;
              }
            else
              goto unknown_option;
            break;
          case '3':
            str_tolower(argv[i]);
            if (strcmp(argv[i],"-3files") == 0)
              {
                (*driver_type)[(*nbmatrices)] = THREEFILES;
                i+= getfilename(&(*filename)[(*nbmatrices)],
                                (i+1<argc)?argv[i+1]:NULL, "dirname");
                (*nbmatrices)++;
              }
            else
              goto unknown_option;
            break;

          default:
          unknown_option:
            fprintf(stderr,
                    "ERROR: main: unprocessed option (\"%s\")\n", argv[i]);
          usage:
            global_usage(MPI_COMM_WORLD,argv);
            free(*driver_type);
            free(*filename);
            MPI_Finalize();
            return EXIT_FAILURE;

          }
        }
      if (maxmatrices == (*nbmatrices))
        {
          maxmatrices *=2;
          (*driver_type) = (driver_type_t*)realloc((*driver_type),
                                                   maxmatrices *
                                                   sizeof(driver_type_t));
          (*filename)    = (char **       )realloc((*filename),
                                                   maxmatrices*sizeof(char*));
        }
      i++;
    }

  /* default driver */
  if ((*nbmatrices) == 0)
    {
      (*driver_type)[0] = RSA;
      (*filename)[0] = malloc((strlen("rsaname")+1)*sizeof(char));
      strcpy((*filename)[0],"rsaname");
      (*nbmatrices) ++;
    }
  return EXIT_SUCCESS;
}

/*
  Function: get_idparm

  Get options from argv.

  Parameters:
  argc          - number of arguments.
  argv          - argument tabular.
  iparm         - type of driver (output, -1 if not set).
  dparm         - type of driver (output, -1 if not set).
*/
int get_idparm(int            argc,
               char         **argv,
               pastix_int_t  *iparm,
               double        *dparm)
{
  int i             = 1;

  while(i < argc)
    {
      if (argv[i][0] == '-')
        {
          switch (argv[i][1]) {

          case 'd':
          case 'D':
            str_tolower(argv[i]);
            if (strcmp(argv[i], "-dparmfile") == 0)
              {
                i++;
                api_dparmreader(argv[i], dparm);
              }
            else if (strcmp(argv[i], "-dparm") == 0)
              {
                int    dparm_idx;
                double value;
                char * endptr;
                i++;
                dparm_idx = (int)strtol(argv[i], &endptr, 10);
                if (endptr == argv[i])
                  {
                    if( 1 == api_str_to_int(argv[i], &dparm_idx))
                      goto unknown_option;
                  }
                i++;
                value = (double)strtod(argv[i], &endptr);
                if (endptr == argv[i])
                  goto unknown_option;
                dparm[dparm_idx] = value;

              }
            break;
          case 'i':
          case 'I':
            str_tolower(argv[i]);
            if (strcmp(argv[i], "-iparmfile") == 0)
              {
                i++;
                api_iparmreader(argv[i], iparm);
              }
            else if (strcmp(argv[i], "-iparm") == 0)
              {
                int iparm_idx, value;
                char * endptr;
                i++;
                iparm_idx = (int)strtol(argv[i], &endptr, 10);
                if (endptr == argv[i])
                  {
                    if( 1 == api_str_to_int(argv[i], &iparm_idx))
                      goto unknown_option;
                  }
                i++;
                value = (int)strtol(argv[i], &endptr, 10);
                if (endptr == argv[i])
                  if( 1 == api_str_to_int(argv[i], &value))
                    goto unknown_option;
                iparm[iparm_idx] = value;
              }
            break;
          }
        }
      i++;
    }

  return EXIT_SUCCESS;
 unknown_option:
  fprintf(stderr,
          "ERROR: main: unprocessed option (\"%s\")\n", argv[i]);
  global_usage(MPI_COMM_WORLD,argv);
  MPI_Finalize();
  return EXIT_FAILURE;

}
