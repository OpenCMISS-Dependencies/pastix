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
  File: chbread.c

  Read a matrix in chb format.

  Matrix can be complex or not.
  
  If the matrix is complex and TYPE_COMPLEX is not defined, 
  imaginary part will be dropped.

  If the matrix is real and TYPE_COMPLEX is defined, the 
  imaginary part will be 0.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
#include "chbread.h"


/*
  Function: chbReadHeader

  Reads header from a chb matrix file.

  header format is:
  > title 73-80 Key
  > 1-14 totcrd 15-28 ptrcrd 29-42 indcrd 43-56 valcrd 57-70 rhscrd
  > 1-3 mxtype 15-28 nrow 29-42 ncol 43-56 nnzero 57-70 neltvl 
  > 1-16 ptrfmt 17-32 indfmt 33-52 valfmt 53-72 rhsfmt
  > 1  2 rhstyp 3  15-28 nrhs 29-42 nrhsix

  Parameters
    infile  - File to read from
    Type    - Type of the matrix
    Nrow    - Number of row in the matrix
    Ncol    - Number of columns in the matrix
    Nnzero  - Number of non zeros in the matrix
    Nrhs    - Number of right-hand-side terms
    Ptrfmt  - 
    Indfmt  -
    Valfmt  -
    Rhsfmt  -
    Ptrcrd  -
    Indcrd  -
    Valcrd  -
    Rhscrd  -
    RhsType - Type of right-hand-side term(s)

 */
void chbReadHeader(FILE         *infile, 
		   char         *Type, 
		   pastix_int_t *Nrow, 
		   pastix_int_t *Ncol, 
		   pastix_int_t *Nnzero, 
		   pastix_int_t *Nrhs, 
		   char         *Ptrfmt, 
		   char         *Indfmt, 
		   char         *Valfmt, 
		   char         *Rhsfmt, 
		   pastix_int_t *Ptrcrd, 
		   pastix_int_t *Indcrd, 
		   pastix_int_t *Valcrd, 
		   pastix_int_t *Rhscrd,  
		   char         *RhsType)
{
  char line[BUFSIZ];
  long totcrd;
  long neltvl;
  long temp1,temp2,temp3,temp4;

  /* first line */
  /* 1-72 title 73-80 Key */
  FGETS(line, BUFSIZ, infile);
  
  /* Seconde line */
  /* 1-14 totcrd 15-28 ptrcrd 29-42 indcrd 43-56 valcrd 57-70 rhscrd */
  FGETS(line, BUFSIZ ,infile);
  sscanf(line, "%14ld%14ld%14ld%14ld%14ld", &totcrd, &temp1,&temp2,&temp3,&temp4);
  *Ptrcrd = temp1;
  *Indcrd = temp2;
  *Valcrd = temp3;
  *Rhscrd = temp4;

  /* Third line*/
  /* 1-3 mxtype 15-28 nrow 29-42 ncol 43-56 nnzero 57-70 neltvl */
  FGETS(line, BUFSIZ, infile);
  sscanf(line, "%3c%ld%ld%ld%ld", Type,&temp1,&temp2,&temp3,&neltvl);
  *Nrow = temp1;
  *Ncol = temp2;
  *Nnzero = temp3;
  Type[3] = '\0';
  myupcase(Type);

  /* fourth line */
  /* 1-16 ptrfmt 17-32 indfmt 33-52 valfmt 53-72 rhsfmt */
  FGETS(line, BUFSIZ, infile);

  sscanf(line, "%16c%16c%20c%20c", Ptrfmt, Indfmt, Valfmt, Rhsfmt);
  Ptrfmt[16] = '\0';
  Indfmt[16] = '\0';
  Valfmt[20] = '\0';
  Rhsfmt[20] = '\0';

  /* fifth line optional */
  /* 1  2 rhstyp 3  15-28 nrhs 29-42 nrhsix */
  if (*Rhscrd != 0)
    {
      FGETS(line, BUFSIZ, infile);
      sscanf(line,"%3c%ld", RhsType, &temp1);
      *Nrhs = temp1;
      RhsType[3] = '\0';
    }
  else
    {
      RhsType[0] = '\0';
    }
}

/*
  Function: chbRead

  Reads a matrix in chb format.

  Header is described in <chbReadHeader>
  Formats are sicribed in <chbParseRfmt> and <chbParseIfmt>

  In our file we have
  header
  valuesFormat
  rowFormat
  columnFormat
  (rhsFormat)
  then the columns,
  the rows, 
  the values,
  (the rhs)
  

  Parameters:
    filename - Path to the file to read from
    Nrow     - Number of rows
    Ncol     - Number of columns
    Nnzero   - Number of non zeros
    col      - Index of first element of each column in *row* and *val*
    row      - Row of eah element				       	
    val      - Value of each element				       	
    Type     - Type of the matrix				       	
    RhsType  - Type of the right-hand-side terms.			       	
    rhs      - right-hand-side term(s)


 */
void chbRead(char const      *filename, 
	     pastix_int_t    *Nrow, 
	     pastix_int_t    *Ncol, 
	     pastix_int_t    *Nnzero, 
	     pastix_int_t   **col, 
	     pastix_int_t   **row,
	     pastix_float_t **val, 
	     char           **Type, 
	     char           **RhsType, 
	     pastix_float_t **rhs)
{
  FILE *infile;
  pastix_int_t Nrhs;
  char Ptrfmt[17];
  char Indfmt[17];
  char Valfmt[21];
  char Rhsfmt[21];
  pastix_int_t Ptrcrd, Indcrd, Valcrd, Rhscrd=0;
  pastix_int_t Valperline, Valwidth, Valprec;
  char Valflag;
  pastix_int_t Rhsperline, Rhswidth, Rhsprec;
  char Rhsflag;
  pastix_int_t Indperline, Indwidth, Indflag;
  pastix_int_t Ptrperline, Ptrwidth, Ptrflag;
  char *element;
  pastix_int_t count;
  pastix_int_t iter,item,colcur;
  char line[BUFSIZ];
  pastix_float_t * rhs2 = NULL;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));
  
  infile = fopen(filename, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }
  
  chbReadHeader(infile, *Type, Nrow, Ncol, Nnzero, &Nrhs, Ptrfmt, Indfmt, Valfmt, Rhsfmt,
		&Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, *RhsType);
  printf("CHB: Nrow=%ld Ncol=%ld Nnzero=%ld\n",(long)*Nrow,(long)*Ncol,(long)*Nnzero);

  chbParseRfmt(Valfmt, &Valperline, &Valwidth, &Valprec, &Valflag);
  chbParseIfmt(Indfmt, &Indperline, &Indwidth, &Indflag);
  chbParseIfmt(Ptrfmt, &Ptrperline, &Ptrwidth, &Ptrflag);

  if (Rhscrd != 0)
    {
      chbParseRfmt(Rhsfmt, &Rhsperline, &Rhswidth, &Rhsprec, &Rhsflag);
    }

  element = (char *) malloc(Ptrwidth+1);
  if (element == NULL)
    {
      fprintf(stderr, "chbRead : Not enough memory for element\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  element[Ptrwidth] = '\0';

  (*col) = (pastix_int_t *) malloc((*Ncol+1)*sizeof(pastix_int_t));
  if ((*col) == NULL)
    {
      fprintf(stderr, "chbRead : Not enough memory for *col\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }

  count = 0;
  for (iter=0; iter<Ptrcrd; iter++)
    {
      FGETS(line, BUFSIZ, infile);

      colcur = Ptrflag;

      for (item=0; item<Ptrperline; item++)
	{
	  if (count > (*Ncol))
	    break;

	  strncpy(element, line+colcur, Ptrwidth);
	  (*col)[count] = atoi(element);

	  /*
	    if ((iter==0) || (iter==1))
	    {
	    fprintf(stderr, "count %ld element %s col %ld\n",
	    count,element, (*col)[count]);
	    }
	  */
	  count++;
	  colcur += Ptrwidth;

	}
    }
  memFree_null(element);

  element = (char *) malloc(Indwidth+1);
  if (element == NULL)
    {
      fprintf(stderr, "chbRead : Not enough memory for element\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  element[Indwidth] = '\0';

  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  if ((*row) == NULL)
    {
      fprintf(stderr, "chbRead : Not enough memory for *row\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }

  count = 0;
  for (iter=0; iter<Indcrd; iter++)
    {
      FGETS(line, BUFSIZ, infile);

      colcur = Indflag;
      for (item=0; item<Indperline; item++)
	{
	  if (count == (*Nnzero))
	    break;

	  strncpy(element, line+colcur, Indwidth);
	  (*row)[count] = atoi(element);

	  /*
	    if ((iter==0) || (iter==1))
	    {
	    fprintf(stderr, "count %ld element %s row %ld\n",
	    count, element, (*row)[count]);
	    }
	  */
	  count++;
	  colcur += Indwidth;
	}
    }
  memFree_null(element);


  element = (char *) malloc(Valwidth+1);
  if (element == NULL)
    {
      fprintf(stderr, "ChbRead : Not enough memory for element\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  element[Valwidth] = '\0';
  
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  if ((*val) == NULL)
    {
      fprintf(stderr, "chbRead : Not enough memory for *val\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }

  count = 0;
  for (iter=0; iter<Valcrd; iter++)
    {
      FGETS(line, BUFSIZ, infile);

      colcur = 0;
      for (item=0; item<Valperline; item++)
	{
	  if (MTX_ISCOM(*Type))
	    {
	      if (count == 2*(*Nnzero))
		break;
	    }
	  else
	    {
	      if (count == (*Nnzero))
		break;
	    }

	  strncpy(element, line+colcur, Valwidth);
	  if (Valflag == 'D')
	    if (strchr(element, 'D') != NULL)
	      *(strchr(element, 'D')) = 'E';



	  if (MTX_ISCOM(*Type))
	    {	

	      if (count %2)
		{
#ifdef TYPE_COMPLEX
#if (defined X_ARCHalpha_compaq_osf1)
		  (*val)[(count-1)/2] += PASTIX_FLOAT(0.0, atof(element));
#else
		  (*val)[(count-1)/2] += (pastix_float_t)atof(element)*I;
#endif
#endif
		}
	      else
		{
		  (*val)[count/2] = (pastix_float_t)atof(element);
		}
	    }
	  else
	    {
	      (*val)[count] = (pastix_float_t)atof(element);
	    }

	  count++;
	  colcur += Valwidth;
	}
    }
  memFree_null(element);
  fprintf(stderr, "Bufsiz %d\n", BUFSIZ);

  if (Rhscrd != 0)
    {
      pastix_int_t nbline = Rhscrd;

      if (MTX_ISRHX(*RhsType))
	nbline = Rhscrd/2;

      element = (char *) malloc(Rhswidth+1);
      if (element == NULL)
	{
	  fprintf(stderr, "ChbRead : Not enough memory for element\n");
	  EXIT(MOD_SI,OUTOFMEMORY_ERR);
	}
      element[Rhswidth] = '\0';

      if (MTX_ISRHX(*RhsType))
	(*rhs) = (pastix_float_t *) malloc(2*(*Ncol)*sizeof(pastix_float_t));
      else
	(*rhs) = (pastix_float_t *) malloc((*Ncol)*sizeof(pastix_float_t));
      if ((*rhs) == NULL)
	{
	  fprintf(stderr, "chbRead : Not enough memory for *rhs\n");
          EXIT(MOD_SI,OUTOFMEMORY_ERR);
	}

      count = 0;
      for (iter=0; iter<nbline; iter++)
	{
	  FGETS(line, BUFSIZ, infile);

	  colcur=0;
	  for (item=0; item<Rhsperline; item++)
	    {
	      if (MTX_ISCOM(*Type))
		{
		  if (count == 2*(*Ncol))
		    break;
		}
	      else
		{
		  if (count == (*Ncol))
		    break;
		}

	      strncpy(element, line+colcur, Rhswidth);


	      if (MTX_ISCOM(*Type))
		{
		  if (count % 2)
		    {
#ifdef TYPE_COMPLEX
#if (defined X_ARCHalpha_compaq_osf1)
		      (*rhs)[(count-1)/2] += PASTIX_FLOAT(0.0, atof(element));
#else
		      (*rhs)[(count-1)/2] += (pastix_float_t)atof(element)*I;
#endif
#endif
		    }
		  else
		    {
		      (*rhs)[count/2] = (pastix_float_t)atof(element);
		    }
		}
	      else
		{
		  (*rhs)[count] = (pastix_float_t)atof(element);
		}

	      count++;
	      colcur += Rhswidth;
	    }
	}

      if (MTX_ISRHX(*RhsType))
	{
	  rhs2 = &(*rhs)[*Ncol];
	  if (rhs2 == NULL)
	    {
	      fprintf(stderr, "chbRead : Not enough memory for *rhs2\n");
	      EXIT(MOD_SI,OUTOFMEMORY_ERR);
	    }

	  count = 0;
	  for (iter=0; iter<nbline; iter++)
	    {
	      FGETS(line, BUFSIZ, infile);

	      colcur = 0;
	      for (item=0; item<Rhsperline; item++)
		{
		  if (MTX_ISCOM(*Type))
		    {
		      if (count == 2*(*Ncol))
			break;
		    }
		  else
		    {
		      if (count == (*Ncol))
			break;
		    }

		  strncpy(element, line+colcur, Rhswidth);

		  if (MTX_ISCOM(*Type))
		    {
		      if (count % 2)
			{
#ifdef TYPE_COMPLEX
#if (defined X_ARCHalpha_compaq_osf1)
			  rhs2[(count-1)/2] += PASTIX_FLOAT(0.0, atof(element));
#else
			  rhs2[(count-1)/2] += (pastix_float_t)atof(element)*I;
#endif
#endif
			}
		      else
			{
			  rhs2[count/2] = (pastix_float_t)atof(element);
			}
		    }
		  else
		    {
		      rhs2[count] = (pastix_float_t)atof(element);
		    }

		  count++;
		  colcur += Rhswidth;
		}
	    }
	}
      else
	{
	  rhs2 = NULL;
	}
      memFree_null(element);
    }
  else
    {
      (*rhs) = NULL;
      rhs2 = NULL;
    }
  
  fclose(infile);
}







/*
  Function: hbParseRfmt

  CHB float format parser
  
  Format is like :
  > (3(1P,E25.16)) 
  or
  > (1P,3E25.16)
  or 
  > (1P3E25.16)
  or 
  > (3E25.16)
  or 
  > (3E25)
  for perline = 3, format = E, width = 25 and prec = 16


  Parameters:
    fmt      - format to parse
    perline  - number of element per line 
    width    - 
    prec     - Precision
    flag     - 
*/
void chbParseRfmt(char *fmt, pastix_int_t *perline, pastix_int_t *width, pastix_int_t *prec, char *flag)
{
  myupcase(fmt);
  
  if (strchr(fmt,'E') != NULL)
    {
      *flag = 'E';
    }
  else if (strchr(fmt,'D') != NULL)
    { 
      *flag = 'D';
    }
  else if (strchr(fmt,'F') != NULL)
    { 
      *flag = 'F';
    }

  if (strchr(fmt,'P') != NULL)
    {
      if (strchr(fmt,'P')[1] == ',')
	{
	  if (strchr(fmt,'P')[2] == *flag)
	    {
	      /* (3(1P,E24.16)) */
/* 	      mysubstr2(fmt, '(', *flag, perline); */
 	      mysubstr2(fmt, '(', '(', perline); /* XL : A mon avis c'est plutôt ça... */
	    }
	  else
	    {
	      /* (1P,3E25.16E3) */
	      mysubstr2(fmt, ',', *flag, perline);
	    }
	}

      else
	mysubstr2(fmt, 'P', *flag, perline);
    }
  else
    mysubstr2(fmt, '(', *flag, perline);
  
  if ( strchr(fmt,'.') )
    {
      mysubstr2(fmt, '.', ')', prec);
      mysubstr2(fmt, *flag, '.', width);
    }
  else
    {
      mysubstr2(fmt, *flag, ')', width);
    }
}

/*
  Function:chbParseIfmt

  CHB integer format parser

  format is :
  > (perlineIwidth) or (X,perlineIwidth)
  Parameters:
    fmt      - format to parse
    perline  - number of element per line 
    width    - 
    flag     - 
*/
void chbParseIfmt(char *fmt, pastix_int_t *perline, pastix_int_t *width, pastix_int_t *flag)
{
  myupcase(fmt);

  *flag=0;

  if ( strchr(fmt, ',') )
    {
      mysubstr2(fmt, ',', 'I', perline);
      mysubstr2(fmt, 'I', ')', width);
      if (strchr(fmt, 'X'))
	{
	  *flag=1;
	}
    }
  else
    {
      mysubstr2(fmt, '(', 'I', perline);
      mysubstr2(fmt, 'I', ')', width);
      if (strchr(fmt, 'X'))
	{
	  *flag=1;
	}
    }
}
