/* Copyright INRIA 2004
**
** This file is part of the Scotch distribution.
**
** The Scotch distribution is libre/free software; you can
** redistribute it and/or modify it under the terms of the
** GNU Lesser General Public License as published by the
** Free Software Foundation; either version 2.1 of the
** License, or (at your option) any later version.
**
** The Scotch distribution is distributed in the hope that
** it will be useful, but WITHOUT ANY WARRANTY; without even
** the implied warranty of MERCHANTABILITY or FITNESS FOR A
** PARTICULAR PURPOSE. See the GNU Lesser General Public
** License for more details.
**
** You should have received a copy of the GNU Lesser General
** Public License along with the Scotch distribution; if not,
** write to the Free Software Foundation, Inc.,
** 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
**
** $Id: common_error.c 2 2004-06-02 14:05:03Z ramet $
*/
/*
  File: common_error
  
  Part of a parallel direct block solver.
  
  This module handles errors.

  Authors: 
    Mathieu Faverge    - faverge@labri.fr     
    David    GOUDIN     - .
    Pascal   HENON      - .
    Xavier   LACOSTE    - lacoste@labri.fr
    Francois PELLEGRINI - .
    Pierre   RAMET      - ramet@labri.fr

  
  Dates: 
    Version 0.0  - from 08 may 1998
                   to   02 oct 1998
*/

/*
**  The defines and includes.
*/

#define COMMON_ERROR

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "common_pastix.h" 

/* 
   Group: Variables

   string: errorProgName
     Name of the programm running
 */
static char                 errorProgName[32] = "";

/*
  Group: Functions
  
  The error handling routines.

*/


/* 
   Function: errorProg
   
   This routine sets the program name for
   error reporting.
   
   Parameters:
     progstr - Name of the program, will be truncated to 29 characters.

   Returns:
     VOID - in all cases.
*/
void errorProg (const char * const progstr)
{
  strncpy (errorProgName, progstr, 29);
  errorProgName[29] = '\0';
  strcat  (errorProgName, ": ");
}

/*
  Function: errorPrint

  This routine prints an error message with
  a variable number of arguments, as printf ()
  does, and exits.

  Parameters:
    errstr - Format for the error to string.
    ...    - arguments depending on the format.
               printf-like variable argument list.

  Returns:
    VOID - in all cases.
*/

void errorPrint (const char * const errstr,             
		 ...)
{
  va_list             errlist;                    /* Argument list of the call */

  fprintf  (stderr, "\n%sERROR: ", errorProgName);
  va_start (errlist, errstr);
  vfprintf (stderr, errstr, errlist);             /* Print arguments */
  va_end   (errlist);
  fprintf  (stderr, "\n\n");
  fflush   (stderr);                              /* In case it has been set to buffered mode */
}

/* 
   Funciton: errorPrintW

   This routine prints a warning message with
   a variable number of arguments, as printf ()
   does.
   
   Parameters:
     warnstr - Format for the warning to print.
     ...     - arguments depending on the format, 
               printf-like variable argument list.

   Returns:
     VOID - in all cases.
*/

void errorPrintW (const char * const warnstr,            
		  ...)
{
  va_list             errlist;                    /* Argument list of the call */

  fprintf  (stderr, "\n%sWARNING: ", errorProgName);
  va_start (errlist, warnstr);
  vfprintf (stderr, warnstr, errlist);             /* Print arguments */
  va_end   (errlist);
  fprintf  (stderr, "\n\n");
  fflush   (stderr);                              /* In case it has been set to buffered mode */
}
