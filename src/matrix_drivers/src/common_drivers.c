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
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>

#ifdef FORCE_NOMPI
#include "pastix_nompi.h"
#else
#include <mpi.h>
#endif
#include <stdlib.h>

#ifdef X_ARCHsun
#include <inttypes.h>
#endif

#include "pastix.h"

#include "common_drivers.h"

#ifdef __INTEL_COMPILER
/* Ignore icc remark : "external declaration in primary source file" */
#pragma warning(disable:1419)
#endif

/* Trouver une solution plus propre pour cette 
   déclarer fonction interne de la libpastix */

#include<ctype.h>

/*
  Function: myupcase

  Rewrites *s* to upper case.
  
  Parameters: 
    s - string to rexwrite in upcase.
*/
void myupcase(char *S)
{
  pastix_int_t iter=0;

  while (S[iter] != '\0')
    {
      S[iter] = (char)toupper(S[iter]);
      iter++;
    }
}


/*
  Function:  mysubstr

  Copy len element, from *S[pos]* in *s*.
  
  Parameters:
    s   - destination
    S   - Source
    pos - sarting position
    len - Size of the string to copy.
*/
void mysubstr(char *s, const char *S, const pastix_int_t pos, const pastix_int_t len)
{
  pastix_int_t iter;
  for (iter=0; iter<len;iter++)
    {
      s[iter] = S[pos+iter];
    }
  s[len] = '\0';
}

/*
  Function:  mysubstr2

  Copy the number placed between a and b in fmt.
  
  Parameters:
    fmt - String in which there is a and b
    a   - first element
    b   - last element
    val - the integer between a and b

*/
void mysubstr2(const char *fmt, const char a, const char b, pastix_int_t *val)
{
  char * posb = strchr(fmt,b);
  char * posa = strchr(fmt,a);
  int len = posb - posa - 1;
  char *tmp = (char *) malloc(len+1);
  if (tmp == NULL)
    {
      fprintf(stderr, "mysubstr2 : Not enough memory for tmp\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  mysubstr(tmp, fmt, strchr(fmt, a) - fmt +1, len);
  *val = atoi(tmp);
  memFree_null(tmp);
}
